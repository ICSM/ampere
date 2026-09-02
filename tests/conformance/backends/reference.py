"""The reference backend fixture: pure numpy ``ampere.core``, as it stands.

``architecture.md`` §1 rung 1 makes the reference backend "pure numpy/scipy",
and §2 settles that it is a real backend rather than a promoted legacy path —
"closer in size to a large contract-conformance test fixture than to a
production backend". Until ``ampere/backends/reference/`` lands in Phase 2, the
reference path *is* ``ampere.core`` plus a handful of concrete models and
transformations, and those live here.

Everything in this module is ordinary numpy. The classes are deliberately the
smallest thing that satisfies :mod:`tests.conformance.protocol`, so that a
Phase-2 author reading it sees the shape of the obligation and not a second
implementation to keep in step.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    AxisRequirement,
    ChannelRequirements,
    DenseGP,
    GPSolver,
    HierarchicalPrior,
    Kernel,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Plate,
    QuasisepGP,
    Spectrum,
    SquaredExponential,
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

__all__ = [
    "LinearModel",
    "Photometry",
    "PowerLawModel",
    "Rebin",
    "ReferenceBackend",
    "ReferenceParameterSpace",
    "Scale",
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


class _CountingModel(Model):
    """Shared plumbing: the evaluation counter and the optional template cache.

    The counter is what lets ``test_inference.py`` assert that an
    out-of-support θ is refused *without* running the model (``inference.md``
    §18). The template cache is the ``compile_for`` behaviour
    ``transformations.md`` §14 asks a compiled model for.
    """

    def __init__(self, spec: ModelSpec) -> None:
        self.spec = spec
        self.evaluations = 0
        self.templates: dict[str, Spectrum] = {}
        self.register_buffer(
            "wavelength", np.asarray(spec.coordinates, dtype=float), unit=COORDINATE_UNIT
        )
        if spec.plated:
            for parameter in offset_plate(len(spec.channels)).expand():
                self.register_parameter(parameter)

    def reset_evaluations(self) -> None:
        self.evaluations = 0

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        if not self.spec.compiled:
            return super().compile_for(requirements)
        for channel in self.spec.channels:
            asked = requirements.get(channel)
            if asked is None or "spectral_axis" not in asked:
                continue
            grid = _to_micron(asked["spectral_axis"].coordinates())
            self.templates[channel] = Spectrum(
                grid * COORDINATE_UNIT, np.zeros(grid.size), unit=FLUX_UNIT
            )
        return self

    def _grid(self, channel: str, context: Mapping[str, Any]) -> np.ndarray:
        template = self.templates.get(channel)
        if template is None:
            return np.asarray(context["wavelength"], dtype=float)
        return template.spectral_axis.values

    def _emit(self, channel: str, grid: np.ndarray, flux: np.ndarray) -> Spectrum:
        template = self.templates.get(channel)
        if template is None:
            return Spectrum(grid * COORDINATE_UNIT, flux, unit=FLUX_UNIT)
        return template.with_values(flux)

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        raise NotImplementedError

    def evaluate(self, **values: Any) -> ModelResult:
        self.evaluations += 1
        context = self.context(values)
        offsets = context.get("objects.offsets")
        channels: dict[str, Spectrum] = {}
        for index, channel in enumerate(self.spec.channels):
            grid = self._grid(channel, context)
            flux = self._flux(grid, context)
            if offsets is not None:
                flux = flux + np.asarray(offsets, dtype=float)[index]
            channels[channel] = self._emit(channel, grid, flux)
        return ModelResult(channels)


class LinearModel(_CountingModel):
    """``f(x) = offset + slope * x`` — the closed form, evaluated directly."""

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return context["offset"] + context["slope"] * grid


class PowerLawModel(_CountingModel):
    """``f(x) = norm * (x / x_ref) ** index`` — evaluated as written."""

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_buffer("reference", float(spec.reference_coordinate))
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return context["norm"] * (grid / context["reference"]) ** context["index"]


# ---------------------------------------------------------------------------
# Transformations
# ---------------------------------------------------------------------------


class Scale(Transformation):
    """A multiplicative calibration with one free parameter."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.register_parameter(Parameter("scale", st.lognorm(0.2)))

    def apply(self, samples: Any, values: Any) -> Spectrum:
        return samples.with_values(samples.values * self.context(values)["scale"])


class Rebin(Transformation):
    """Resample onto a coarser grid, publishing a requirement and a mask."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(self, target: Any, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.target = np.asarray(target, dtype=float)

    def requirements(self) -> tuple[AxisRequirement, ...]:
        return (
            AxisRequirement(
                "spectral_axis",
                unit=COORDINATE_UNIT,
                intervals=(float(self.target[0]), float(self.target[-1])),
                max_step=float(np.diff(self.target).min()) / 2.0,
            ),
        )

    def influence(self, source: np.ndarray) -> np.ndarray:
        return influence_matrix(np.asarray(source, dtype=float), self.target)

    def apply(self, samples: Any, values: Any) -> Spectrum:
        weights = self.influence(samples.spectral_axis.values)
        return Spectrum(
            self.target * COORDINATE_UNIT,
            weights @ samples.values,
            unit=samples.unit,
            mask=propagate_mask(samples, weights),
        )


class Photometry(Transformation):
    """``Spectrum → PhotometricPoints``: the kind-changing step."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = PhotometricPoints

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


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=float)
    return np.asarray(coordinates, dtype=float)


_MODELS = {ModelKind.LINEAR: LinearModel, ModelKind.POWER_LAW: PowerLawModel}
_KERNELS = {KernelFamily.MATERN32: Matern32, KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential}


class ReferenceBackend:
    """The first fixture: ``ampere.core``'s numpy path."""

    name = "reference"
    capabilities = BackendCapabilities(
        differentiable=False,
        batchable=False,
        device="cpu",
        float64=True,
        # QuasisepGP is a declared strategy slot with no implementation
        # (likelihoods.md §5); Phase 2 fills it in and adds it here.
        solvers=frozenset({SolverKind.DENSE}),
    )

    def model(self, spec: ModelSpec) -> Model:
        return _MODELS[spec.kind](spec)

    def transformation(self, spec: TransformationSpec) -> Transformation:
        if spec.kind is TransformationKind.SCALE:
            return Scale(label=spec.label)
        if spec.kind is TransformationKind.REBIN:
            return Rebin(spec.target, label=spec.label)
        return Photometry(spec.target, spec.filters, label=spec.label)

    def kernel(self, spec: CovarianceSpec) -> Kernel:
        return _KERNELS[spec.family](spec.amplitude, spec.length_scale)

    def gp_solver(self, kind: SolverKind) -> GPSolver:
        return DenseGP() if kind is SolverKind.DENSE else QuasisepGP()

    def parameter_space(self, declaration: ParameterSet) -> ReferenceParameterSpace:
        return ReferenceParameterSpace(declaration)

    def to_numpy(self, values: Any) -> np.ndarray:
        return np.asarray(values)
