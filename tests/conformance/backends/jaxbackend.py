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
``DENSE`` only. ``ampere.core.QuasisepGP`` is celerite2's *numpy* solver, and
handing it back here would satisfy every agreement row while proving nothing —
which is exactly what ``TestSolverAgreement``'s "a declared quasiseparable
strategy must be a different strategy from the dense one" row exists to catch.
Slice 2 of W2.5 chooses between ``tinygp``'s ``QuasisepSolver`` and
``celerite2.jax`` by measuring both against this battery; until then the
quasiseparable rows skip, with the suite's own reason naming what is owed.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import jax.numpy as jnp
import numpy as np
import scipy.stats as st

from ampere.backends.jax import (
    CalibrationScale,
    DenseGP,
    Matern32,
    PowerLaw,
    Resample,
    SquaredExponential,
    configure_x64,
)
from ampere.backends.jax.models import _SpectralModel
from ampere.backends.jax.parameters import LoweredParameterSet
from ampere.core import (
    GPSolver,
    Kernel,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Spectrum,
    Transformation,
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
    "JaxBackend",
    "JaxLinearModel",
    "JaxPhotometry",
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


_MODELS = {ModelKind.LINEAR: JaxLinearModel, ModelKind.POWER_LAW: JaxPowerLawModel}
_KERNELS = {KernelFamily.MATERN32: Matern32, KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential}


class JaxBackend:
    """The jax fixture: the shipped ``ampere.backends.jax`` package."""

    name = BACKEND
    capabilities = BackendCapabilities(
        # The capability the whole backend exists for. It is the flag every
        # part declares, and `ampere.backends.jax.lower_problem` is what cashes
        # it: a pure jax log-density a NUTS kernel can differentiate.
        differentiable=True,
        # Slice 2's `vmap` work; claiming it today would be a promise unkept.
        batchable=False,
        device="cpu",
        # `configure_x64()` below, and every class raises without it.
        float64=True,
        # See this module's docstring: DENSE only until slice 2 has chosen
        # between tinygp's QuasisepSolver and celerite2.jax.
        solvers=frozenset({SolverKind.DENSE}),
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
        return _KERNELS[spec.family](spec.amplitude, spec.length_scale)

    def gp_solver(self, kind: SolverKind) -> GPSolver:
        if kind is not SolverKind.DENSE:
            raise NotImplementedError(
                f"the {BACKEND!r} backend declares no {kind.value} solver; see this module's "
                f"docstring for what slice 2 owes."
            )
        return DenseGP()

    def parameter_space(self, declaration: ParameterSet) -> LoweredParameterSet:
        return LoweredParameterSet(declaration)

    def to_numpy(self, values: Any) -> np.ndarray:
        return np.asarray(values)
