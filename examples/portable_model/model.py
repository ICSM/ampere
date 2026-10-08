"""The quickstart's ``LinearModel``, written once against the portable base."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import scipy.stats as st

from ampere.core import Parameter, PortableModel

__all__ = ["LinearModel"]


class LinearModel(PortableModel):
    """``F_nu = slope * lambda + intercept``, in Jy, on the channel ``"sed"``.

    The quickstart's declaration -- the same two parameters with the same flat
    priors on [-10, 10] -- with its arithmetic moved into :meth:`_flux`. There
    is no ``np.asarray`` on the grid and no ``float(...)`` on a parameter:
    ``grid`` arrives as the backend's array and every parameter value is one
    too, so this one method runs on numpy, jax and torch.
    """

    def __init__(self, wavelength: Any, *, channels: str | Sequence[str] = "sed") -> None:
        super().__init__(wavelength, channels=channels)
        self.register_parameter(Parameter("slope", st.uniform(-10, 20)))
        self.register_parameter(Parameter("intercept", st.uniform(-10, 20)))

    def _flux(self, grid: Any, context: Mapping[str, Any]) -> Any:
        return context["slope"] * grid + context["intercept"]
