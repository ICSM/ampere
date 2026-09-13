"""Astrometric time series on the torch path.

The native twin of :mod:`ampere.backends.reference.astrometry` (W4.9): the
same reflex-orbit model and the same epoch-sampling step, with every number a
gradient passes through computed in ``torch``, so a fit to injected astrometry
runs under NUTS.

**The inheriting pattern**, the same choice W4.3 made for interferometry and
for the same shape of reason, restated for this modality: :class:`EpochSample`
here is the identity on its container (``apply`` returns what it is given —
see the reference module), so there is *no* arithmetic to fork at all, and
:class:`ReflexOrbit`'s declaration (the two-channel template, the
per-channel negotiated grid, ``compile_for``) is exactly what both backends
must agree on bit-for-bit, while its arithmetic is four lines of ``sin``/
``cos``. Both are cases where re-declaring risks the two backends drifting on
something that matters (which epochs a channel's template actually holds)
while buying nothing (there is nothing to re-declare that is not shared), so
both classes here derive from their reference counterpart and override only
the four capability flags, the dtype/device plumbing, and — for
``ReflexOrbit`` alone — the four lines that are actually differentiable.

:class:`ReflexOrbit`'s native surface is spelled ``grid``/``flux`` rather than
``native_grid``/``native_flux`` (the interferometry twins' spelling,
:mod:`ampere.backends.torch.interferometry`): that spelling exists only
because an interferometric source model already has a parameter called
``flux``, and ``Parameterised._check_free_name`` refuses a method that would
shadow it. Nothing this model declares is called ``grid`` or ``flux``, so the
plain names are used, exactly as :class:`~ampere.backends.torch.models.
TorchSpectralModel` uses them for the other one-axis, ``Layout.POINTS`` kind.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.reference.astrometry import DAYS_PER_YEAR
from ampere.backends.reference.astrometry import EpochSample as _ReferenceEpochSample
from ampere.backends.reference.astrometry import ReflexOrbit as _ReferenceReflexOrbit
from ampere.core import ChannelRequirements

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, move, place, to_numpy
from .parameters import LoweredParameters

__all__ = ["EpochSample", "ReflexOrbit"]


class TorchAstrometryStep:
    """The two capability flags this step needs, and torch's placement.

    A minimal mixin rather than :mod:`ampere.backends.torch.interferometry`'s
    ``TorchInterferometryStep``: :class:`EpochSample` registers no buffer this
    backend's arithmetic ever touches (its one buffer, ``epochs``, is read by
    :meth:`~ampere.backends.reference.astrometry.EpochSample.requirements`
    alone, in numpy, at composition time), so there is nothing here for a
    tensor cache to hold.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        place(self, dtype, device)

    def to(self, *, dtype: Any = None, device: Any = None) -> TorchAstrometryStep:
        """Re-declare where this step lives. There is nothing to move."""
        return move(self, dtype=dtype, device=device)


class EpochSample(TorchAstrometryStep, _ReferenceEpochSample):
    """Pin a model channel onto the observed epochs, in torch — still no arithmetic.

    ``apply`` is inherited unchanged (it returns the container it is given,
    on either backend); the native surface is the identity for the same
    reason.
    """

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
        """The identity: the coordinates and the values are already correct."""
        return flux, grid


class ReflexOrbit(_ReferenceReflexOrbit):
    """A linear proper motion plus a periodic reflex wobble, in torch.

    Every declarative part — the two fixed channels, ``compile_for``'s
    per-channel adoption of the negotiated epochs — is the reference class's,
    unchanged; what is forked is :meth:`_offset`, the four-line expression a
    gradient with respect to ``pmra``, ``pmdec``, ``period``, ``phase``,
    ``amp_ra`` or ``amp_dec`` passes through.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        place(self, dtype, device)
        #: The torch-side home for this model's constant arrays
        #: (``lowering.md`` §7): the model's own default time grid.
        self.tensors = LoweredParameters()
        self.tensors.register_buffer(
            "time",
            as_tensor(self._data("time"), dtype=self.dtype, device=self.device),
            persistent=True,
        )
        self._grids: dict[str, torch.Tensor] = {}

    # -- negotiation -----------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> ReflexOrbit:
        """Adopt the negotiated epochs, then cache them as tensors."""
        super().compile_for(requirements)
        for channel, template in self.templates.items():
            self._grids[channel] = as_tensor(
                template.time.values, dtype=self.dtype, device=self.device
            )
        return self

    def to(self, *, dtype: Any = None, device: Any = None) -> ReflexOrbit:
        """Move this model's grids, and re-declare where they live."""
        move(self, dtype=dtype, device=device)
        self._grids = {
            channel: grid.to(dtype=self.dtype, device=self.device)
            for channel, grid in self._grids.items()
        }
        return self

    # -- the native surface a realisation composes ---------------------------

    def grid_tensor(self, channel: str) -> torch.Tensor:
        """The time tensor this channel evaluates on: negotiated, or the model's own."""
        found = self._grids.get(channel)
        if found is not None:
            return found
        return self.tensors.get_buffer("time")

    def grid(self, channel: str) -> torch.Tensor:
        """:meth:`grid_tensor` under the name the realisation looks for."""
        return self.grid_tensor(channel)

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> torch.Tensor:
        """This model's offset on *channel*, mas, as a differentiable tensor."""
        return self._offset_native(
            channel, self.grid_tensor(channel), self._context_tensors(self.context(values))
        )

    def _context_tensors(self, context: Mapping[str, Any]) -> dict[str, torch.Tensor]:
        """This model's vocabulary as tensors, keeping any graph it arrived with."""
        return {
            name: as_tensor(value, dtype=self.dtype, device=self.device)
            for name, value in context.items()
        }

    # -- the contract surface ------------------------------------------------

    def _offset(self, channel: str, t: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        """The torch offset, back on the numpy side of the container boundary.

        The reference class's ``evaluate`` calls this, so the contract path —
        and therefore every conformance row that compares this model against
        the closed-form ephemeris — runs torch's arithmetic rather than
        numpy's.
        """
        return to_numpy(
            self._offset_native(
                channel,
                as_tensor(t, dtype=self.dtype, device=self.device),
                self._context_tensors(context),
            )
        )

    def _offset_native(
        self, channel: str, t: torch.Tensor, context: Mapping[str, torch.Tensor]
    ) -> torch.Tensor:
        cycle = 2.0 * math.pi * t / context["period"] + context["phase"]
        if channel == "ra":
            return context["pmra"] * t / DAYS_PER_YEAR + context["amp_ra"] * torch.sin(cycle)
        return context["pmdec"] * t / DAYS_PER_YEAR + context["amp_dec"] * torch.cos(cycle)
