"""The torch backend fixture: ``ampere.backends.torch``, the shipped package.

Registered through ``_optional`` in ``tests/conformance/backends/__init__.py``,
so an environment without the ``torch`` extra simply does not have this column
— the rows are backend-agnostic, and a missing backend is reported by its
absence from the fixture ids rather than by a skipped row
(``architecture.md`` §4 rule 2 forbids making ``import ampere`` require torch).

What the battery therefore exercises for real, on every row:

* :class:`~ampere.backends.torch.PowerLaw` — the ``POWER_LAW`` kind;
* :class:`~ampere.backends.torch.CalibrationScale` — ``SCALE``;
* :class:`~ampere.backends.torch.Resample` — ``REBIN``;
* :class:`~ampere.backends.torch.TorchParameterSpace` — every row in
  ``test_parameters.py``, against ``scipy``, and every row in
  ``test_cross_backend.py``, against the reference backend's own answer;
* :class:`~ampere.backends.torch.DenseGP` — the GP rows, through torch's
  Cholesky rather than scipy's.

Three pieces stay local, for the same reasons they do in the reference fixture.
``LinearModel`` realises ``ModelKind.LINEAR``, which is a *test* construct — a
straight line in wavelength is not on ``DEVELOPMENT_PLAN.md`` §5's list of
models a backend ships, and shipping one only to satisfy a fixture would put
test scaffolding in the installed package. ``Photometry`` realises
``TransformationKind.PHOTOMETRY``, which exists "so a row can build a chain
whose kinds do not compose" and is declared by bare pivots and names; the
shipped :class:`~ampere.backends.torch.SyntheticPhotometry` is the physics,
takes tabulated response curves and publishes ``points=``, which the battery's
declaration cannot supply. ``PointSourceModel`` realises ``ModelKind.COMPLEX``
(W2.4 slice 3): a flat complex visibility with a winding phase, which is what
the ``complex_gaussian`` rows need a channel of. A real visibility model
belongs with the modality it serves — plan §5 puts that in Phase 4 — so the
battery's is the battery's, exactly as its straight line is. All three are
written in torch here, like everything else this fixture returns, so no row in
the battery runs numpy arithmetic under the ``torch`` id.

The counting wrapper is likewise test-only: ``protocol.py``'s ``CountingModel``
is the battery's requirement, not something a shipped model should carry.

``SolverKind.QUASISEP`` is declared since W2.4 slice 2
------------------------------------------------------
It was not in slice 1, and the reason it was not is worth keeping: declaring
``QUASISEP`` and handing back a dense solver would satisfy every agreement row
and prove nothing, which is what ``tests/conformance/README.md`` §6 says the
no-fallback row exists to catch. Slice 2 made the measurement
``DEVELOPMENT_PLAN.md`` §6 asked for — GPyTorch's structured solvers against
celerite2's compiled kernels — and :class:`~ampere.backends.torch.QuasisepGP`
is the result: celerite2's semiseparable factorisation under
``torch.autograd``, exact and differentiable in the hyperparameters. So the
five rows that skipped here now run, against a solver that really is a
different recursion from the dense one.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st
import torch

from ampere.backends.torch import (
    BACKEND,
    CalibrationScale,
    DenseGP,
    GaussianProcessNoise,
    IndependentNoise,
    Matern32,
    PowerLaw,
    QuasisepGP,
    Resample,
    SquaredExponential,
    TorchParameterSpace,
    TorchSpectralModel,
    as_tensor,
    to_numpy,
)
from ampere.core import (
    GPSolver,
    HierarchicalPrior,
    Kernel,
    Model,
    ModelResult,
    NoiseModel,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Plate,
    Spectrum,
    Transformation,
    VisibilitySet,
    propagate_mask,
)

from ..protocol import (
    BackendCapabilities,
    CovarianceSpec,
    KernelFamily,
    ModelKind,
    ModelSpec,
    SolverKind,
    TransformationKind,
    TransformationSpec,
    complex_axes,
)

__all__ = [
    "BACKEND",
    "LinearModel",
    "Photometry",
    "PointSourceModel",
    "PowerLawModel",
    "TorchBackend",
]

FLUX_UNIT = u.Jy
COORDINATE_UNIT = u.micron


def offset_plate(size: int) -> Plate:
    """The plate :attr:`ModelSpec.plated` declares: one offset per channel."""
    return Plate(
        "objects",
        size=size,
        hyperparameters=[
            Parameter("mu", st.norm(0.0, 1.0)),
            Parameter("sigma", st.halfnorm(0.0, 1.0)),
        ],
        members=[Parameter("offsets", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
    )


def influence_matrix(source: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """An ``(n_out, n_in)`` resampling matrix from *source* onto *target*.

    Every input sample is assigned to its nearest target coordinate and the
    members of each target bin are averaged; a target that claims no input
    falls back to the single nearest one, so no row is empty and no output
    sample is silently zero. Deliberately mixing — a diagonal matrix would make
    the ANY mask rule vacuous.
    """
    distance = torch.abs(target[:, None] - source[None, :])
    nearest = torch.argmin(distance, dim=0)
    weights = torch.zeros(
        (int(target.numel()), int(source.numel())), dtype=source.dtype, device=source.device
    )
    for out in range(int(target.numel())):
        members = torch.nonzero(nearest == out, as_tuple=True)[0]
        if members.numel() == 0:
            members = torch.argmin(distance[out]).reshape(1)
        weights[out, members] = 1.0 / float(members.numel())
    return weights


class _CountingModel(TorchSpectralModel):
    """The battery's counter and its opt-out from ``compile_for``.

    ``protocol.py``'s ``CountingModel`` asks for an ``evaluations`` counter so
    ``test_inference.py`` can prove an out-of-support theta is refused
    *without* running the model. That is a test obligation, not something a
    shipped model should carry, so it is mixed in here.
    """

    #: W2.12 item 4: every part a fixture composes declares the fixture's own
    #: ``name`` as its backend, so the row asserting ``problem.backend ==
    #: backend.name`` is a real assertion.
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

    def evaluate_tensor(self, **values: Any) -> dict[str, torch.Tensor]:
        """The shipped evaluation, plus the plate's per-channel offset.

        ``ModelSpec.plated`` adds ``offsets[i]`` to channel *i* so that
        ``results.md`` §14's row has an array-valued parameter with a named
        dimension to write. That is a declaration the battery needs and the
        shipped models have no reason to carry, so it is applied here.
        """
        self.evaluations += 1
        context = self._context_tensors(values)
        offsets = context.get("objects.offsets")
        emitted: dict[str, torch.Tensor] = {}
        for index, channel in enumerate(self.channels):
            grid = self.grid_tensor(channel)
            flux = self._flux(grid, context)
            if offsets is not None:
                flux = flux + offsets.reshape(-1)[index]
            emitted[channel] = flux
        return emitted


class LinearModel(_CountingModel):
    """``f(x) = offset + slope * x`` — the closed form, evaluated directly.

    Local to the battery: a straight line in wavelength is a test construct,
    not one of the models ``DEVELOPMENT_PLAN.md`` §5 has a backend ship. It
    exists because the battery needs a form whose parameters are both
    unbounded, so both lower to ``Identity``.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: torch.Tensor, context: Mapping[str, torch.Tensor]) -> torch.Tensor:
        return context["offset"] + context["slope"] * grid


class PowerLawModel(_CountingModel, PowerLaw):
    """The battery's ``POWER_LAW`` kind, on the **shipped** torch power law.

    The maths, the buffer and the parameter declaration all come from
    :class:`~ampere.backends.torch.PowerLaw`; only the counter and the
    ``compiled=False`` opt-out are added.
    """

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(
            spec,
            norm=st.loguniform(0.1, 10.0),
            index=st.norm(-1.0, 0.5),
            reference_wavelength=spec.reference_coordinate,
        )


class PointSourceModel(Model):
    """``V(x) = norm * exp(i * index * x)`` — the battery's ``COMPLEX`` kind.

    Local to the battery, and complex all the way through: the buffer is real
    (the ``u`` axis), the flux is a ``torch.complex128`` tensor, and
    :meth:`evaluate` emits a :class:`~ampere.core.VisibilitySet`, which is the
    only container kind ``results_schema.md`` §16 allows complex values in.

    It does not subclass :class:`~ampere.backends.torch.TorchSpectralModel`,
    and that is the point rather than an omission: that class is spectrum
    shaped — a micron axis, a Jy ``Spectrum`` per channel, a real ``_flux`` —
    and a visibility model shares none of it but the parameter plumbing. What
    it does share is the **native surface** a realisation walks (``grid`` and
    ``flux``), which is written out here so that a reviewer can see the whole
    of what :mod:`ampere.backends.torch.problem` requires of a model in one
    screen.

    The second (``v``) axis comes from
    :func:`~tests.conformance.protocol.complex_axes`, the same rule the
    observed container uses; see that function for why it is not written twice.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, spec: ModelSpec) -> None:
        u_axis, v_axis = complex_axes(spec.coordinates)
        self.spec = spec
        self.evaluations = 0
        self.channels = tuple(spec.channels)
        self.dtype = torch.complex128
        self.device = torch.device("cpu")
        self.register_buffer("u", u_axis)
        self.register_buffer("v", v_axis)
        self._u = as_tensor(u_axis)
        self._v = as_tensor(v_axis)
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))

    def reset_evaluations(self) -> None:
        self.evaluations = 0

    # -- the native surface a realisation composes ---------------------------

    def grid(self, channel: str) -> torch.Tensor:
        """The ``u`` axis. One channel, so *channel* is only validated."""
        if channel not in self.channels:
            raise KeyError(channel)
        return self._u

    def flux(self, channel: str, values: Any = None) -> torch.Tensor:
        """``norm * exp(i * index * u)`` as a complex tensor, gradient intact.

        ``torch.polar`` rather than ``exp(1j * phase)``: it takes a real
        modulus and a real angle and builds the complex number from them, so
        both parameters stay real leaves of the graph and autograd never has to
        differentiate through a complex-valued exponential. The value is the
        same; the derivative is the one a real-parameter model should have.
        """
        context = self.context({} if values is None else values)
        norm = as_tensor(context["norm"])
        index = as_tensor(context["index"])
        return torch.polar(norm.expand(self._u.shape), index * self._u)

    # -- the contract surface -------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        self.evaluations += 1
        emitted = {
            channel: VisibilitySet(
                to_numpy(self._u),
                to_numpy(self._v),
                to_numpy(self.flux(channel, values)).astype(np.complex128, copy=False),
                unit=FLUX_UNIT,
            )
            for channel in self.channels
        }
        return ModelResult(emitted)


class Photometry(Transformation):
    """``Spectrum -> PhotometricPoints``: the battery's kind-changing step.

    Local to the battery. ``protocol.py`` says this kind "exists so a row can
    build a chain whose kinds do not compose", and declares it by bare pivots
    and filter names — there are no response curves to integrate. The shipped
    :class:`~ampere.backends.torch.SyntheticPhotometry` is the real step.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = PhotometricPoints
    DIFFERENTIABLE: ClassVar[bool] = True
    #: As the shipped steps, since W2.4 slice 2; see the fixture's
    #: ``capabilities``.
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    #: W2.12 item 4, as on :class:`_CountingModel`.
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, target: Any, filters: Any, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.target = as_tensor(np.asarray(target, dtype=float))
        self.filters = tuple(str(name) for name in filters)

    def influence(self, source: Any) -> np.ndarray:
        return to_numpy(influence_matrix(as_tensor(source), self.target))

    def apply_flux(self, flux: Any, grid: Any, values: Any) -> tuple[Any, Any]:
        """The native surface: the same matrix as :meth:`apply`, on bare tensors.

        Added at W2.5 slice 2, because the battery gained a shape that needs it
        (``PHOTOMETRIC`` in ``test_inference.py``). Until then no
        ``ProblemSpec`` ended in a kind-changing step, so this fixture's step
        was never composed into a problem and never reached a realisation —
        which is exactly the hole that shape closes, and why the jax backend's
        shipped photometry could be broken for a whole slice without anything
        noticing.

        The coordinates come back as this step's own target rather than the
        incoming grid, because photometry changes the axis as well as the kind.
        """
        return influence_matrix(as_tensor(grid), self.target) @ flux, self.target

    def apply(self, samples: Any, values: Any) -> PhotometricPoints:
        weights = influence_matrix(as_tensor(samples.spectral_axis.values), self.target)
        return PhotometricPoints(
            self.filters,
            to_numpy(self.target) * COORDINATE_UNIT,
            to_numpy(weights @ as_tensor(samples.values)),
            unit=samples.unit,
            mask=propagate_mask(samples, to_numpy(weights)),
        )


_MODELS = {
    ModelKind.LINEAR: LinearModel,
    ModelKind.POWER_LAW: PowerLawModel,
    ModelKind.COMPLEX: PointSourceModel,
}
_KERNELS = {KernelFamily.MATERN32: Matern32, KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential}


class TorchBackend:
    """The torch fixture: the shipped ``ampere.backends.torch`` package.

    Its ``name`` is the package's own ``BACKEND`` string, imported rather than
    retyped, because W2.12 fixed that these are the same thing: one name per
    backend, everywhere.
    """

    name = BACKEND
    capabilities = BackendCapabilities(
        # The pieces are differentiable — every model and step this fixture
        # returns has a tensor-valued entry point whose gradient
        # ``tests/backends/`` takes. The end-to-end gradient through a composed
        # problem additionally needs a tensor-valued evaluation path that
        # ``ampere.core``'s numpy containers do not admit; see the package
        # docstring, and W2.4's report.
        differentiable=True,
        # Honest, and it changed at W2.4 slice 2: every piece this fixture
        # returns declares BATCHABLE, and the realised density is evaluated
        # over a stack in one call through ``torch.func.vmap``
        # (``LoweredProblem.log_prob_unconstrained_batched``). The claim is
        # conjunctive, so a problem carrying ``QuasisepGP`` — a compiled
        # extension vmap cannot see through — still reports False, which the
        # row asserting problem.batchable == this value would catch.
        batchable=True,
        device="cpu",
        # architecture.md §5's policy, not a preference: every tensor is built
        # float64 and torch's global default dtype is never touched.
        float64=True,
        # W2.4 slice 3: ``PointSourceModel`` realises ModelKind.COMPLEX, so the
        # ``complex_gaussian`` rows run here rather than skipping.
        complex_models=True,
        # Both, since W2.4 slice 2: see the module docstring.
        solvers=frozenset({SolverKind.DENSE, SolverKind.QUASISEP}),
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
        # This backend's own since W2.13. The *declaration* is still
        # ``ampere.core``'s -- these subclass it, so same FAMILY, same
        # HYPERPARAMETERS, same QUASISEPARABLE flag and therefore the same spec
        # hash -- but the covariance is built in torch, which is what makes a
        # GP hyperparameter differentiable. With the core kernels the Cholesky
        # had a gradient and the amplitude did not: W2.4's carried finding.
        return _KERNELS[spec.family](spec.amplitude, spec.length_scale)

    def gp_solver(self, kind: SolverKind) -> GPSolver:
        return DenseGP() if kind is SolverKind.DENSE else QuasisepGP()

    def independent_noise(self) -> NoiseModel:
        # This backend's own, since W2.13: a noise model is a capability part
        # now, so composing ampere.core's would declare two backends.
        return IndependentNoise()

    def gp_noise(self, kernel: Kernel, solver: GPSolver, *, jitter: Any = None) -> NoiseModel:
        return GaussianProcessNoise(kernel, solver, jitter=jitter)

    def parameter_space(self, declaration: ParameterSet) -> TorchParameterSpace:
        return TorchParameterSpace(declaration)

    def to_numpy(self, values: Any) -> np.ndarray:
        return to_numpy(values)
