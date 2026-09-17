"""The reference backend fixture: ``ampere.backends.reference``, the shipped package.

Until W2.1 this module *was* the reference backend — ``ampere.core`` plus a
handful of concrete models and transformations, written here because
``ampere/backends/reference/`` did not exist yet. It does now, so this file has
become what it was always meant to be: a thin adapter turning the battery's
neutral declarations into the shipped classes.

What the battery therefore exercises for real, on every row:

* :class:`~ampere.backends.reference.PowerLaw` — the ``POWER_LAW`` kind;
* :class:`~ampere.backends.reference.CalibrationScale` — ``SCALE``;
* :class:`~ampere.backends.reference.Resample` — ``REBIN``;
* :class:`~ampere.backends.reference.FractionalModelNoise` and
  :class:`~ampere.backends.reference.FractionalModelGPNoise` — the X-1 rows in
  ``test_likelihoods.py``.

Two pieces stay local, and deliberately:

``LinearModel``
    ``ModelKind.LINEAR`` is a *test* construct — the battery needs a closed
    form whose parameters both lower to ``Identity`` — and a straight line in
    wavelength is not on ``DEVELOPMENT_PLAN.md`` §5's list of models the
    reference backend ships. Shipping one only to satisfy a fixture would put
    test scaffolding in the installed package.
``Photometry``
    ``TransformationKind.PHOTOMETRY`` exists, per ``protocol.py``, "so a row
    can build a chain whose kinds do not compose": it is a *kind*-changing
    fixture, declared by bare pivots and names. The shipped
    :class:`~ampere.backends.reference.SyntheticPhotometry` is the physics — it
    takes tabulated response curves, publishes ``points=`` at its own
    tabulation and reads the container through ``Axis.locate`` — which the
    battery's declaration cannot supply. It is covered by ``tests/backends/``
    instead, including Gap 1's two-instrument scenario.

The counting wrapper is likewise test-only: :class:`CountingModel` is
``protocol.py``'s requirement, not something a shipped model should carry, so
it is mixed in here rather than in the package.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.backends.reference import (
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    CalibrationScale,
    ClosurePhase,
    EpochSample,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    PSFConvolution,
    PowerLaw,
    ReflexOrbit,
    Resample,
    TimeSmearing,
    UniformDisc,
    UniformDiscVisibilities,
)
from ampere.backends.reference.models import _SpectralModel
from ampere.core import (
    SHO,
    DenseGP,
    HilbertSpaceGP,
    GaussianProcessNoise,
    JointGaussianProcessNoise,
    GPSolver,
    HierarchicalPrior,
    IndependentNoise,
    Kernel,
    Matern12,
    Matern32,
    Matern52,
    NoiseModel,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Plate,
    QuasisepGP,
    RotationTerm,
    Spectrum,
    SquaredExponential,
    Transformation,
    propagate_mask,
)

from ._kernels import build_kernel
from ..protocol import (
    AstrometryPieces,
    ImagePieces,
    BackendCapabilities,
    CovarianceSpec,
    InterferometryPieces,
    KernelFamily,
    ModelKind,
    ModelSpec,
    SolverKind,
    TransformationKind,
    TransformationSpec,
)

__all__ = [
    "REFERENCE_ASTROMETRY",
    "REFERENCE_INTERFEROMETRY",
    "LinearModel",
    "Photometry",
    "PowerLawModel",
    "ReferenceBackend",
    "ReferenceParameterSpace",
    "influence_matrix",
    "offset_plate",
]

FLUX_UNIT = u.Jy
COORDINATE_UNIT = u.micron


# ---------------------------------------------------------------------------
# Shared numerics
# ---------------------------------------------------------------------------


def offset_plate(size: int) -> Plate:
    """The plate :attr:`ModelSpec.plated` declares: one offset per channel.

    ``mu`` and ``sigma`` are the population hyperparameters and ``offsets`` is
    the array-valued member drawn from them — the smallest declaration that
    gives ``emit`` an array parameter with a *named* dimension to write.
    """
    return Plate(
        "objects",
        size=size,
        hyperparameters=[
            Parameter("mu", st.norm(0.0, 1.0)),
            Parameter("sigma", st.halfnorm(0.0, 1.0)),
        ],
        members=[Parameter("offsets", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
    )


def influence_matrix(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """An ``(n_out, n_in)`` resampling matrix from *source* onto *target*.

    Every input sample is assigned to its nearest target coordinate and the
    members of each target bin are averaged; a target that claims no input
    falls back to the single nearest one, so no row is empty and no output
    sample is silently zero. Deliberately mixing — a diagonal matrix would
    make the ANY mask rule vacuous.
    """
    distance = np.abs(target[:, None] - source[None, :])
    nearest = np.argmin(distance, axis=0)
    weights = np.zeros((target.size, source.size), dtype=float)
    for out in range(target.size):
        members = np.flatnonzero(nearest == out)
        if members.size == 0:
            members = np.array([int(np.argmin(distance[out]))])
        weights[out, members] = 1.0 / members.size
    return weights


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class _CountingModel(_SpectralModel):
    """The battery's counter and its opt-out from ``compile_for``.

    ``protocol.py``'s :class:`~tests.conformance.protocol.CountingModel` asks
    for an ``evaluations`` counter so ``test_inference.py`` can prove an
    out-of-support theta is refused *without* running the model. That is a test
    obligation, not something a shipped model should carry, so it is mixed in
    here.

    ``ModelSpec.compiled=False`` also has to be reproducible, and the shipped
    models always honour ``compile_for``. Refusing it here is a one-line
    override rather than a flag on the installed class.
    """

    #: W2.12 item 4: every part a fixture composes declares the fixture's own
    #: ``name`` as its backend, so the row asserting ``problem.backend ==
    #: backend.name`` is a real assertion. Stated here rather than inherited
    #: from ``Model`` so that the tie to :data:`BACKEND` is visible.
    BACKEND: ClassVar[str] = "reference"

    def __init__(self, spec: ModelSpec, **kwargs: Any) -> None:
        super().__init__(spec.coordinates, channels=spec.channels, **kwargs)
        self.spec = spec
        self.evaluations = 0
        if spec.plated:
            for parameter in offset_plate(len(spec.channels)).expand():
                self.register_parameter(parameter)

    def reset_evaluations(self) -> None:
        self.evaluations = 0

    def compile_for(self, requirements: Mapping[str, Any]) -> Model:
        if not self.spec.compiled:
            return Model.compile_for(self, requirements)
        return super().compile_for(requirements)

    def evaluate(self, **values: Any) -> ModelResult:
        """The shipped evaluation, plus the plate's per-channel offset.

        ``ModelSpec.plated`` adds ``offsets[i]`` to channel *i* so that
        ``results.md`` §14's row has an array-valued parameter with a named
        dimension to write. That is a declaration the battery needs and the
        shipped models have no reason to carry, so it is applied here.
        """
        self.evaluations += 1
        context = self.context(values)
        offsets = context.get("objects.offsets")
        emitted: dict[str, Spectrum] = {}
        for index, channel in enumerate(self.channels):
            grid = self._grid(channel, context)
            flux = self._flux(grid, context)
            if offsets is not None:
                flux = flux + np.asarray(offsets, dtype=float)[index]
            emitted[channel] = self._emit(channel, grid, flux)
        return ModelResult(emitted)


class LinearModel(_CountingModel):
    """``f(x) = offset + slope * x`` — the closed form, evaluated directly.

    Local to the battery: a straight line in wavelength is a test construct,
    not one of the models ``DEVELOPMENT_PLAN.md`` §5 has the reference backend
    ship. It exists because the battery needs a form whose parameters are both
    unbounded, so both lower to ``Identity``.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return context["offset"] + context["slope"] * grid


class PowerLawModel(_CountingModel, PowerLaw):
    """The battery's ``POWER_LAW`` kind, on the **shipped** power law.

    The maths, the buffer and the parameter declaration all come from
    :class:`~ampere.backends.reference.PowerLaw`; only the counter and the
    ``compiled=False`` opt-out are added. The battery's priors are passed in,
    because they are the declaration ``ModelKind.POWER_LAW`` fixes and the
    shipped class rightly takes whatever the user asks for.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(
            spec,
            norm=st.loguniform(0.1, 10.0),
            index=st.norm(-1.0, 0.5),
            reference_wavelength=spec.reference_coordinate,
        )


class Photometry(Transformation):
    """``Spectrum -> PhotometricPoints``: the battery's kind-changing step.

    Local to the battery. ``protocol.py`` says this kind "exists so a row can
    build a chain whose kinds do not compose", and declares it by bare pivots
    and filter names — there are no response curves to integrate. The shipped
    :class:`~ampere.backends.reference.SyntheticPhotometry` is the real step,
    and is covered by ``tests/backends/``.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = PhotometricPoints
    #: W2.12 item 4, as on :class:`_CountingModel`.
    BACKEND: ClassVar[str] = "reference"

    def __init__(self, target: Any, filters: Any, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.target = np.asarray(target, dtype=float)
        self.filters = tuple(str(name) for name in filters)

    def influence(self, source: np.ndarray) -> np.ndarray:
        return influence_matrix(np.asarray(source, dtype=float), self.target)

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = self.influence(samples.spectral_axis.values)
        return PhotometricPoints(
            self.filters,
            self.target * COORDINATE_UNIT,
            weights @ samples.values,
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )


# ---------------------------------------------------------------------------
# Parameter space
# ---------------------------------------------------------------------------


class ReferenceParameterSpace:
    """``ParameterSet`` itself, behind the :class:`ParameterSpace` protocol."""

    def __init__(self, declaration: ParameterSet) -> None:
        self._declaration = declaration

    @property
    def declaration(self) -> ParameterSet:
        return self._declaration

    @property
    def free_size(self) -> int:
        return self._declaration.free_size

    def free_labels(self) -> tuple[str, ...]:
        return self._declaration.free_labels()

    def pack(self, values: Mapping[str, Any]) -> np.ndarray:
        return self._declaration.pack(values)

    def unpack(self, theta: Any) -> dict[str, Any]:
        return self._declaration.unpack(theta)

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        return self._declaration.prior_transform(unit_cube)

    def lnprior(self, values: Any) -> float:
        return self._declaration.lnprior(values)

    def constrain(self, unconstrained: Any) -> np.ndarray:
        return self._declaration.constrain(unconstrained)

    def unconstrain(self, values: Any) -> np.ndarray:
        return self._declaration.unconstrain(values)

    def lnprior_unconstrained(self, unconstrained: Any) -> float:
        return self._declaration.lnprior_unconstrained(unconstrained)


# ---------------------------------------------------------------------------
# The backend
# ---------------------------------------------------------------------------


_MODELS = {ModelKind.LINEAR: LinearModel, ModelKind.POWER_LAW: PowerLawModel}
_KERNELS: dict[KernelFamily, type[Kernel]] = {
    KernelFamily.MATERN12: Matern12,
    KernelFamily.MATERN32: Matern32,
    KernelFamily.MATERN52: Matern52,
    KernelFamily.SHO: SHO,
    KernelFamily.ROTATION: RotationTerm,
    KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential,
}


class ReferenceBackend:
    """The first fixture: the shipped ``ampere.backends.reference`` package.

    Its ``name`` is the shipped backend's own ``BACKEND`` string, and W2.12
    fixed that these are the same thing: one name per backend, everywhere.
    """

    name = "reference"
    capabilities = BackendCapabilities(
        differentiable=False,
        batchable=False,
        device="cpu",
        float64=True,
        # W2.3 filled the quasiseparable slot in with celerite2 (an exact
        # rank-2 representation of Matern-3/2, not celerite2's eps
        # approximation), so the DenseGP<->QuasisepGP agreement rows run.
        solvers=frozenset({SolverKind.DENSE, SolverKind.QUASISEP, SolverKind.HILBERT}),
        # W4.1: this is the only backend with an interferometric vocabulary
        # until W4.3 writes the native twins, so every interferometric row
        # skips elsewhere with a reason naming that.
        interferometry=True,
        # W4.9: likewise the only backend with an astrometric vocabulary
        # until the torch and jax fixtures declare their own twins below.
        astrometry=True,
        joint_noise=True,  # W5.9
        # W5.5: likewise the only backend with a gridded-image vocabulary
        # until the torch and jax fixtures declare their own twins below.
        image=True,
    )

    def model(self, spec: ModelSpec) -> Model:
        return _MODELS[spec.kind](spec)

    def transformation(self, spec: TransformationSpec) -> Transformation:
        if spec.kind is TransformationKind.SCALE:
            return CalibrationScale(st.lognorm(0.2), label=spec.label)
        if spec.kind is TransformationKind.REBIN:
            return Resample(spec.target, label=spec.label)
        return Photometry(spec.target, spec.filters, label=spec.label)

    def kernel(self, spec: CovarianceSpec) -> Kernel:
        return build_kernel(spec, _KERNELS)

    def gp_solver(
        self,
        kind: SolverKind,
        *,
        basis_size: int | Sequence[int] = 32,
        boundary_factor: float = 2.0,
    ) -> GPSolver:
        if kind is SolverKind.DENSE:
            return DenseGP()
        if kind is SolverKind.HILBERT:
            return HilbertSpaceGP(basis_size=basis_size, boundary_factor=boundary_factor)
        return QuasisepGP()

    def independent_noise(self) -> NoiseModel:
        # The core classes, unmodified: since W2.13 they declare
        # ``BACKEND = "reference"`` (the ABC's default), which is the truth
        # here and the reason this fixture alone keeps them.
        return IndependentNoise()

    def gp_noise(self, kernel: Kernel, solver: GPSolver, *, jitter: Any = None) -> NoiseModel:
        return GaussianProcessNoise(kernel, solver, jitter=jitter)

    def joint_gp_noise(
        self,
        kernel: Kernel,
        solver: GPSolver,
        *,
        datasets: Sequence[str],
        coupling: Any,
    ) -> NoiseModel:
        # W5.9 -- appended.
        return JointGaussianProcessNoise(kernel, solver, datasets=datasets, coupling=coupling)

    def parameter_space(self, declaration: ParameterSet) -> ReferenceParameterSpace:
        return ReferenceParameterSpace(declaration)

    def interferometry(self) -> InterferometryPieces:
        # The shipped classes, unmodified: they declare BACKEND = "reference",
        # which is the truth here and the reason this fixture alone keeps them.
        return REFERENCE_INTERFEROMETRY

    def image(self) -> ImagePieces:
        return REFERENCE_IMAGE

    def astrometry(self) -> AstrometryPieces:
        # The shipped classes, unmodified, for the same reason as above.
        return REFERENCE_ASTROMETRY

    def to_numpy(self, values: Any) -> np.ndarray:
        return np.asarray(values)


#: The shipped astrometric vocabulary, as one record (W4.9). Module-level so
#: that :class:`MirrorBackend` can say which of it is its own.
REFERENCE_IMAGE = ImagePieces(
    psf_convolution=PSFConvolution,
    gaussian_source=GaussianSource,
    binary=Binary,
)

#: The shipped gridded-image vocabulary, as one record (W5.5).
REFERENCE_ASTROMETRY = AstrometryPieces(
    epoch_sample=EpochSample,
    reflex_orbit=ReflexOrbit,
)

#: The shipped interferometric vocabulary, as one record (W4.1). Module-level
#: so that :class:`MirrorBackend` can say which of it is its own.
REFERENCE_INTERFEROMETRY = InterferometryPieces(
    fourier_sample=FourierSample,
    closure_phase=ClosurePhase,
    amplitude=Amplitude,
    bandwidth_smearing=BandwidthSmearing,
    time_smearing=TimeSmearing,
    uniform_disc=UniformDisc,
    gaussian_source=GaussianSource,
    binary=Binary,
    uniform_disc_visibilities=UniformDiscVisibilities,
    gaussian_source_visibilities=GaussianSourceVisibilities,
    binary_visibilities=BinaryVisibilities,
)
