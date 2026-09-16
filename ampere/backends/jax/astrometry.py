"""Astrometric time series on the jax path.

The native twin of :mod:`ampere.backends.reference.astrometry` (W4.9), as
:mod:`ampere.backends.jax.interferometry` is on the other modern backend: the
same reflex-orbit model and the same epoch-sampling step, with the arithmetic
in ``jax.numpy``.

**The inheriting pattern**, for the reason
:mod:`ampere.backends.torch.astrometry`'s module docstring gives in full:
:class:`EpochSample` has no arithmetic at all (``apply`` is the identity), and
:class:`ReflexOrbit`'s declaration — the two fixed channels, the per-channel
negotiated epoch grid, ``compile_for`` — is exactly what both backends must
agree on bit-for-bit, while its arithmetic is four lines. Both classes here
derive from their reference counterpart and override only the capability
flags, this backend's device placement, and — for ``ReflexOrbit`` alone — the
four-line expression a gradient passes through.

``ReflexOrbit``'s native surface is spelled ``native_grid``/``native_flux``,
which since *W5.20* is the canonical spelling on every backend and every
modality — not a special case the interferometry twins were driven into by a
parameter of their own called ``flux``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np

from ampere.backends.reference.astrometry import DAYS_PER_YEAR
from ampere.backends.reference.astrometry import EpochSample as _ReferenceEpochSample
from ampere.backends.reference.astrometry import ReflexOrbit as _ReferenceReflexOrbit
from ampere.core import ChannelRequirements, TransformationError

from ._config import BACKEND, require_x64
from ._device import DEVICE, device_flag, place_on, resolve_device

__all__ = ["EpochSample", "ReflexOrbit"]


class _JaxAstrometryStep:
    """The two capability flags this step needs, and jax's placement.

    A minimal mixin, for the reason
    :mod:`ampere.backends.torch.astrometry`'s ``TorchAstrometryStep`` gives:
    :class:`EpochSample` has no array this backend's arithmetic ever touches.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        self._device = resolved
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))


class EpochSample(_JaxAstrometryStep, _ReferenceEpochSample):
    """Pin a model channel onto the observed epochs, in jax — still no arithmetic."""

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """The identity: the coordinates and the values are already correct."""
        return flux, grid


class ReflexOrbit(_ReferenceReflexOrbit):
    """A linear proper motion plus a periodic reflex wobble, in jax.

    Every declarative part is the reference class's, unchanged; what is
    forked is :meth:`_offset`, the four-line expression a gradient with
    respect to ``pmra``, ``pmdec``, ``period``, ``phase``, ``amp_ra`` or
    ``amp_dec`` passes through.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        self._grids: dict[str, jax.Array] = {}
        self._own_grid = self._place(self._data("time"))

    def _place(self, values: Any) -> jax.Array:
        return place_on(
            jnp.asarray(np.asarray(values, dtype=float)), getattr(self, "_device", None)
        )

    # -- negotiation -----------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> ReflexOrbit:
        """Adopt the negotiated epochs, then cache them as jax arrays."""
        super().compile_for(requirements)
        for channel, template in self.templates.items():
            self._grids[channel] = self._place(template.time.values)
        return self

    # -- the native surface a realisation composes ---------------------------

    def grid_tensor(self, channel: str) -> jax.Array:
        """The time array this channel evaluates on: negotiated, or the model's own."""
        found = self._grids.get(channel)
        if found is not None:
            return found
        return self._own_grid

    def native_grid(self, channel: str) -> jax.Array:
        """:meth:`grid_tensor` under the name the realisation looks for."""
        return self.grid_tensor(channel)

    def native_flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's offset on *channel*, mas, as a differentiable array."""
        return self._offset_native(channel, self.grid_tensor(channel), self.context(values))

    # -- the contract surface ------------------------------------------------

    def _offset(self, channel: str, t: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        """The jax offset, back on the numpy side of the container boundary."""
        return np.asarray(self._offset_native(channel, self._place(t), context))

    def _offset_native(self, channel: str, t: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        cycle = 2.0 * math.pi * t / context["period"] + context["phase"]
        if channel == "ra":
            return context["pmra"] * t / DAYS_PER_YEAR + context["amp_ra"] * jnp.sin(cycle)
        return context["pmdec"] * t / DAYS_PER_YEAR + context["amp_dec"] * jnp.cos(cycle)
