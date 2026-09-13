"""The jax backend fixture: ``ampere.backends.jax``, the shipped package.

The adapter turning the battery's neutral declarations into this backend's
concrete classes. Nothing here is a second implementation of anything: the
models, instrument steps, kernels and solver all come from
``ampere.backends.jax``, and what stays local is exactly what stayed local in
``reference.py`` and for the same reasons —

``LinearModel``
    ``ModelKind.LINEAR`` is a *test* construct (the battery needs a closed form
    whose parameters both lower to ``Identity``), and a straight line in
    wavelength is not on ``DEVELOPMENT_PLAN.md`` §5's list of models a backend
    ships.
``Photometry``
    ``TransformationKind.PHOTOMETRY`` exists "so a row can build a chain whose
    kinds do not compose": a kind-changing fixture declared by bare pivots and
    names. The shipped :class:`~ampere.backends.jax.SyntheticPhotometry` is the
    physics — tabulated response curves, ``points=`` at its own tabulation,
    ``Axis.locate`` — which the battery's declaration cannot supply.

and the ``evaluations`` counter, which is ``protocol.py``'s requirement rather
than something a shipped model should carry.

x64
---
:class:`JaxBackend` calls ``configure_x64()`` in its **constructor**. That is
the point of ``lowering.md`` §10.2: ampere never flips the process-global flag
as an import side effect, and the application — here, the test suite — turns it
on explicitly before any jax work. Constructing the fixture is the first jax
work this suite does, so the call belongs there. Every class in
``ampere.backends.jax`` raises if it is skipped, which is conformance §11's
row 11 and is asserted directly in ``tests/backends/test_jax.py``.

Solvers
-------
``DENSE`` and ``QUASISEP``, both this backend's own since W2.5 slice 2. The
quasiseparable one is :class:`ampere.backends.jax.QuasisepGP` — celerite2's
*jax* interface over ampere's own exact rank-2 Matérn-3/2 representation, not
``ampere.core.QuasisepGP``, which is celerite2's numpy solver and would
satisfy every agreement row while proving nothing (which is exactly what
``TestSolverAgreement``'s "a declared quasiseparable strategy must be a
different strategy from the dense one" row exists to catch, and what its
backend flag now catches too).
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import jax.numpy as jnp
import numpy as np
import scipy.stats as st

from ampere.backends.jax import (
    SHO,
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    CalibrationScale,
    ClosurePhase,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    TimeSmearing,
    UniformDisc,
    UniformDiscVisibilities,
    DenseGP,
    GaussianProcessNoise,
    IndependentNoise,
    Matern12,
    Matern32,
    Matern52,
    PowerLaw,
    QuasisepGP,
    Resample,
    RotationTerm,
    SquaredExponential,
    configure_x64,
)
from ampere.backends.jax._config import require_x64
from ampere.backends.jax._device import DEVICE
from ampere.backends.jax.models import _SpectralModel
from ampere.backends.jax.parameters import LoweredParameterSet
from ampere.core import (
    GPSolver,
    Kernel,
    Model,
    NoiseModel,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Spectrum,
    Transformation,
    VisibilitySet,
    propagate_mask,
)

from ._kernels import build_kernel
from ..protocol import (
    BackendCapabilities,
    CovarianceSpec,
    InterferometryPieces,
    KernelFamily,
    ModelKind,
    ModelSpec,
    SolverKind,
    TransformationKind,
    TransformationSpec,
    complex_axes,
)
from .reference import influence_matrix, offset_plate

#: This backend's one name everywhere (W2.12): ``ConformanceBackend.name``, the
#: ``BACKEND`` every part it composes declares, ``FittingProblem.backend``,
#: ``ampere_backend`` in a run's provenance, and the key
#: ``ampere.core.lowering``'s registry is consulted with.
BACKEND = "jax"

FLUX_UNIT = u.Jy
COORDINATE_UNIT = u.micron

__all__ = [
    "BACKEND",
    "JAX_INTERFEROMETRY",
    "JaxBackend",
    "JaxLinearModel",
    "JaxPhotometry",
    "JaxPointSourceModel",
    "JaxPowerLawModel",
]


class _CountingModel(_SpectralModel):
    """The battery's counter, its ``compile_for`` opt-out, and the plate offsets.

    ``protocol.py``'s ``CountingModel`` asks for an ``evaluations`` counter so
    ``test_inference.py`` can prove an out-of-support θ is refused *without*
    running the model. That is a test obligation, not something a shipped model
    should carry, so it is mixed in here — as it is on the reference fixture.
    """

    BACKEND: ClassVar[str] = BACKEND

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

    def _offset(self, context: Mapping[str, Any], channel: str) -> Any:
        """``offsets[i]`` for *channel*, or zero when the spec declares no plate.

        ``ModelSpec.plated`` adds one offset per channel so that ``results.md``
        §14's row has an array-valued parameter with a *named* dimension to
        write. Indexed in jax, so the native path keeps its gradient with
        respect to the plate member.
        """
        offsets = context.get("objects.offsets")
        if offsets is None:
            return jnp.asarray(0.0, dtype=jnp.float64)
        return jnp.asarray(offsets, dtype=jnp.float64)[self.channels.index(channel)]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> Any:
        context = self.context(values)
        return self._flux(self.grid(channel), context) + self._offset(context, channel)

    def evaluate(self, **values: Any) -> ModelResult:
        self.evaluations += 1
        context = self.context(values)
        emitted = {
            channel: self._emit(
                channel,
                self._flux(self.grids[channel], context) + self._offset(context, channel),
            )
            for channel in self.channels
        }
        return ModelResult(emitted)


class JaxLinearModel(_CountingModel):
    """``f(x) = offset + slope * x`` — the closed form, in ``jax.numpy``.

    Local to the battery: a straight line in wavelength is a test construct,
    not one of the models ``DEVELOPMENT_PLAN.md`` §5 has a backend ship. It
    exists because the battery needs a form whose parameters are both
    unbounded, so both lower to ``Identity``.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: Any, context: Mapping[str, Any]) -> Any:
        return jnp.asarray(context["offset"]) + jnp.asarray(context["slope"]) * grid


class JaxPowerLawModel(_CountingModel, PowerLaw):
    """The battery's ``POWER_LAW`` kind, on the **shipped** jax power law.

    The maths, the buffer and the parameter declaration all come from
    :class:`ampere.backends.jax.PowerLaw`; only the counter, the plate offsets
    and the ``compiled=False`` opt-out are added. The battery's priors are
    passed in, because they are the declaration ``ModelKind.POWER_LAW`` fixes
    and the shipped class rightly takes whatever the user asks for.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(
            spec,
            norm=st.loguniform(0.1, 10.0),
            index=st.norm(-1.0, 0.5),
            reference_wavelength=spec.reference_coordinate,
        )


class JaxPointSourceModel(Model):
    """``V(x) = norm * exp(i * index * x)`` — the battery's ``COMPLEX`` kind, in jax.

    The jax twin of ``torch_backend.PointSourceModel``, added at **W4.3** so that
    W2.4 slice 3's ``complex_gaussian`` rows and W4.2's circular-complex-GP rows
    run on this backend instead of skipping: until now ``ModelKind.COMPLEX``
    existed on the torch fixture alone, and a closed form checked on one backend
    is a closed form checked once.

    Local to the battery, and complex all the way through: the buffers are real
    (the ``u`` and ``v`` axes), the flux is a ``complex128`` array, and
    :meth:`evaluate` emits a :class:`~ampere.core.VisibilitySet`, the only kind
    ``results_schema.md`` §16 allows complex values in.

    It does not subclass ``ampere.backends.jax.models._SpectralModel``, and that
    is the point rather than an omission: that class is spectrum shaped — a
    micron axis, a Jy ``Spectrum`` per channel, a real ``_flux`` — and a
    visibility model shares none of it but the parameter plumbing. What it does
    share is the **native surface** a realisation walks (``grid`` and ``flux``),
    written out here so that a reviewer can see the whole of what
    :mod:`ampere.backends.jax.problem` requires of a model in one screen. The
    shipped visibility models are :mod:`ampere.backends.jax.interferometry`'s;
    this one is the battery's, for the same reason its straight line is.

    The second (``v``) and third (``spectral_axis``) axes come from
    :func:`~tests.conformance.protocol.complex_axes`, the same rule the observed
    container uses, so the predicted and observed axes are equal by construction
    rather than by coincidence.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, spec: ModelSpec) -> None:
        require_x64("a jax JaxPointSourceModel")
        u_axis, v_axis, wavelength = complex_axes(spec.coordinates)
        self._wavelength = wavelength
        self.spec = spec
        self.evaluations = 0
        self.channels = tuple(spec.channels)
        self.register_buffer("u", u_axis)
        self.register_buffer("v", v_axis)
        self._u = jnp.asarray(u_axis, dtype=jnp.float64)
        self._v = jnp.asarray(v_axis, dtype=jnp.float64)
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))

    def reset_evaluations(self) -> None:
        self.evaluations = 0

    # -- the native surface a realisation composes ---------------------------

    def grid(self, channel: str) -> Any:
        """The ``u`` axis. One channel, so *channel* is only validated."""
        if channel not in self.channels:
            raise KeyError(channel)
        return self._u

    def flux(self, channel: str, values: Any = None) -> Any:
        """``norm * exp(i * index * u)`` as a complex array, gradient intact.

        Built from a real modulus and a real angle rather than through
        ``exp(1j * phase)``, so both parameters stay real leaves of the graph and
        the derivative is the one a real-parameter model should have — the same
        choice ``torch.polar`` makes on the other backend.
        """
        context = self.context({} if values is None else values)
        norm = jnp.asarray(context["norm"], dtype=jnp.float64)
        angle = jnp.asarray(context["index"], dtype=jnp.float64) * self._u
        return norm * (jnp.cos(angle) + 1j * jnp.sin(angle))

    # -- the contract surface ------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        self.evaluations += 1
        emitted = {
            channel: VisibilitySet(
                np.asarray(self._u),
                np.asarray(self._v),
                self._wavelength,
                np.asarray(self.flux(channel, values)).astype(np.complex128, copy=False),
                unit=FLUX_UNIT,
            )
            for channel in self.channels
        }
        return ModelResult(emitted)


class JaxPhotometry(Transformation):
    """``Spectrum -> PhotometricPoints``: the battery's kind-changing step, in jax.

    Local to the battery, for the reason ``protocol.py`` gives: this kind
    "exists so a row can build a chain whose kinds do not compose", and is
    declared by bare pivots and filter names — there are no response curves to
    integrate. The shipped :class:`~ampere.backends.jax.SyntheticPhotometry` is
    the real step, and is covered by ``tests/backends/``.

    The influence matrix is the same nearest-neighbour assignment the reference
    fixture uses (:func:`~tests.conformance.backends.reference.influence_matrix`),
    imported rather than rewritten: it is the *declaration* of what this test
    kind does, and two fixtures that computed different matrices would make the
    cross-backend rows compare two different steps.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = PhotometricPoints
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, target: Any, filters: Any, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.target = np.asarray(target, dtype=float)
        self.filters = tuple(str(name) for name in filters)

    def influence(self, source: Any) -> np.ndarray:
        return influence_matrix(np.asarray(source, dtype=float), self.target)

    def apply_flux(self, flux: Any, grid: Any, values: Any) -> tuple[Any, Any]:
        weights = jnp.asarray(self.influence(np.asarray(grid)), dtype=jnp.float64)
        return weights @ flux, jnp.asarray(self.target, dtype=jnp.float64)

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = self.influence(samples.spectral_axis.values)
        integrated = jnp.asarray(weights, dtype=jnp.float64) @ jnp.asarray(
            samples.values, dtype=jnp.float64
        )
        return PhotometricPoints(
            self.filters,
            self.target * COORDINATE_UNIT,
            np.asarray(integrated),
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )


_MODELS = {
    ModelKind.LINEAR: JaxLinearModel,
    ModelKind.POWER_LAW: JaxPowerLawModel,
    ModelKind.COMPLEX: JaxPointSourceModel,
}
_KERNELS: dict[KernelFamily, type[Kernel]] = {
    KernelFamily.MATERN12: Matern12,
    KernelFamily.MATERN32: Matern32,
    KernelFamily.MATERN52: Matern52,
    KernelFamily.SHO: SHO,
    KernelFamily.ROTATION: RotationTerm,
    KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential,
}


#: This backend's interferometric vocabulary, as one record (*W4.3*). The
#: shipped classes, unmodified: they declare ``BACKEND = "jax"``, which is what
#: makes the composed problems in ``test_interferometry.py`` single-backend
#: ones. A twin associates with its reference by class name in the backend's own
#: module, the native surface it exposes and this declaration — there is no step
#: registry (``phase4_placement_memo.md`` §1.2) — so this record is the only
#: place the battery needs to learn them.
JAX_INTERFEROMETRY = InterferometryPieces(
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


class JaxBackend:
    """The jax fixture: the shipped ``ampere.backends.jax`` package."""

    name = BACKEND
    capabilities = BackendCapabilities(
        # The capability the whole backend exists for. It is the flag every
        # part declares, and `ampere.backends.jax.lower_problem` is what cashes
        # it: a pure jax log-density a NUTS kernel can differentiate.
        differentiable=True,
        # Slice 2's `vmap` work, landed: every model, step, kernel, noise
        # model and the dense solver is whole-array `jax.numpy`, so
        # `jax.vmap` over the realised density maps as it maps any pure
        # function -- and `LoweredProblem.log_prob_unconstrained_batched` is
        # the surface. `QuasisepGP` is the one part that still says False
        # (celerite2 registers no batching rule), which is why this claim is
        # about `SINGLE` -- the shape the flags row composes -- and why a
        # quasiseparable problem honestly reports `batchable=False`.
        batchable=True,
        device="cpu",
        # `configure_x64()` below, and every class raises without it.
        float64=True,
        # W3.1: a jax array carries a `jaxlib` `Device` handle, which is
        # process-local and has no pickle reduction, so a jax-composed problem
        # cannot be sent to a worker process — and jax warns, loudly and
        # correctly, that forking a jax process is likely to deadlock. Not a
        # gap: throughput on a device backend is slice 2's per-chunk `vmap`,
        # and `simulate_many` refuses the process pool here by name.
        picklable=False,
        # **W4.3**: ``JaxPointSourceModel`` realises ModelKind.COMPLEX, so the
        # ``complex_gaussian`` rows (W2.4 slice 3) and the circular complex GP
        # rows (W4.2) run here rather than skipping. Until now the kind existed
        # on the torch fixture alone, and a closed form checked on one backend
        # is a closed form checked once.
        complex_models=True,
        # **W4.3**: the native interferometric twins, so every row in
        # ``test_interferometry.py`` runs on this column instead of skipping
        # with a reason naming what W4.3 owed.
        interferometry=True,
        # Both, since W2.5 slice 2 chose celerite2.jax for the O(N) solve.
        solvers=frozenset({SolverKind.DENSE, SolverKind.QUASISEP}),
    )

    def __init__(self) -> None:
        # `lowering.md` §10.2: the *application* turns x64 on, explicitly,
        # before any jax work — never ampere as an import side effect. For this
        # suite the application is the fixture, and constructing it is the
        # first jax work the process does.
        configure_x64()

    def model(self, spec: ModelSpec) -> Model:
        return _MODELS[spec.kind](spec)

    def transformation(self, spec: TransformationSpec) -> Transformation:
        if spec.kind is TransformationKind.SCALE:
            return CalibrationScale(st.lognorm(0.2), label=spec.label)
        if spec.kind is TransformationKind.REBIN:
            return Resample(spec.target, label=spec.label)
        return JaxPhotometry(spec.target, spec.filters, label=spec.label)

    def kernel(self, spec: CovarianceSpec) -> Kernel:
        return build_kernel(spec, _KERNELS)

    def gp_solver(self, kind: SolverKind) -> GPSolver:
        if kind is SolverKind.DENSE:
            return DenseGP()
        if kind is SolverKind.QUASISEP:
            return QuasisepGP()
        raise NotImplementedError(f"the {BACKEND!r} backend declares no {kind.value} solver.")

    def independent_noise(self) -> NoiseModel:
        # This backend's own, since W2.13: a noise model is a capability part
        # now, so composing ampere.core's would declare two backends.
        return IndependentNoise()

    def gp_noise(self, kernel: Kernel, solver: GPSolver, *, jitter: Any = None) -> NoiseModel:
        return GaussianProcessNoise(kernel, solver, jitter=jitter)

    def interferometry(self) -> InterferometryPieces:
        return JAX_INTERFEROMETRY

    def parameter_space(self, declaration: ParameterSet) -> LoweredParameterSet:
        return LoweredParameterSet(declaration)

    def to_numpy(self, values: Any) -> np.ndarray:
        return np.asarray(values)
