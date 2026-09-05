"""A second in-repo backend, so the battery's parametrisation is not notional.

W1.10's acceptance criterion is that "adding a backend requires only a
fixture". A suite with one fixture cannot demonstrate that: the parametrisation
would be untested machinery, and the first Phase-2 author would find out the
hard way which rows had quietly baked in an assumption about the reference
path. This module is the cheapest honest demonstration — a second backend that
runs every row with **zero test-body changes**.

It is *not* a backend in ``architecture.md`` §1's sense, and it must not be
mistaken for one. It shares the reference backend's transformations, kernels
and solvers, because it has nothing different to offer there. What it does
supply of its own is exactly the two things a real backend owns:

* **models** — the same two closed forms declared under the same parameter
  names, computed by a deliberately different route (Horner's rule; ``exp`` of
  a logarithm) in classes of their own. Two backends that agreed only because
  they shared an implementation would prove nothing. (Being honest about how
  far that goes: the power law really does land on different bits — around
  ``3e-17`` on these grids — while Horner's rule on ``offset + slope * x`` is
  the same sequence of operations and agrees exactly. The mirror's value is
  that it is a *separate class the suite never names*, not that every one of
  its numbers differs in the last place.)
* **a parameter space** — one that round-trips the declaration through
  ``to_spec``/``from_spec`` before delegating, so every parameter row is also
  a serialisation row on this fixture.

Delete this module the day a real second backend registers, or keep it: it
costs one fixture and it is the only thing standing between the battery and
silently degenerating to a single-backend suite again.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import Model, Parameter, ParameterSet

from ..protocol import BackendCapabilities, ModelKind, ModelSpec
from .reference import ReferenceBackend, _CountingModel

__all__ = [
    "MirrorBackend",
    "MirrorLinearModel",
    "MirrorParameterSpace",
    "MirrorPowerLawModel",
]


class MirrorLinearModel(_CountingModel):
    """``offset + slope * x`` by Horner's rule rather than as written."""

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return np.polyval([context["slope"], context["offset"]], grid)


class MirrorPowerLawModel(_CountingModel):
    """``norm * (x / x_ref) ** index`` through the logarithm."""

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        # The unit is part of the declaration, and W2.1's neutral model
        # identity compares buffers: a reference wavelength declared unitless
        # here and in micron there is a real disagreement, not a detail.
        self.register_buffer(
            "reference_wavelength", float(spec.reference_coordinate), unit=u.micron
        )
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        ratio = np.log(grid / context["reference_wavelength"])
        return np.exp(context["index"] * ratio + np.log(context["norm"]))


class MirrorParameterSpace:
    """A parameter space reached through the serialisable declaration.

    ``ParameterSet.to_spec()`` is what ``parameters.md`` §13 calls the
    serialisable declaration and what ``results.md`` hashes into the
    provenance attrs. Realising the space from *that* rather than from the
    object it was built from means every row in ``test_parameters.py`` also
    asserts, on this fixture, that the declaration survives the trip.
    """

    def __init__(self, declaration: ParameterSet) -> None:
        self._declaration = ParameterSet.from_spec(declaration.to_spec())

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


_MODELS = {ModelKind.LINEAR: MirrorLinearModel, ModelKind.POWER_LAW: MirrorPowerLawModel}


class MirrorBackend(ReferenceBackend):
    """The second fixture. Same contracts, different arithmetic."""

    name = "mirror"
    capabilities = BackendCapabilities(
        differentiable=False,
        batchable=False,
        device="cpu",
        float64=True,
        solvers=ReferenceBackend.capabilities.solvers,
        tolerances=ReferenceBackend.capabilities.tolerances,
    )

    def model(self, spec: ModelSpec) -> Model:
        return _MODELS[spec.kind](spec)

    def parameter_space(self, declaration: ParameterSet) -> MirrorParameterSpace:
        return MirrorParameterSpace(declaration)
