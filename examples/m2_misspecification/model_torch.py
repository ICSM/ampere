"""The toy model in ``torch``: the same physics, differentiably.

Imported only when the study runs on the torch backend — ``import
examples.m2_misspecification`` does not reach this module, because
``architecture.md`` §4 rule 2's discipline ("``import ampere`` must never
require torch") is worth honouring in an example as well as in the package.

What makes this class differentiable is not a flag but two methods.
:mod:`ampere.backends.torch.problem` composes a problem's log-density out of
its parts' **native** surfaces, and for a model those are ``grid(channel)`` and
``flux(channel, values)``; it finds them by ``hasattr``, so this hand-written
class lowers exactly as ``ampere.backends.torch.PowerLaw`` does, without
subclassing it or registering anything.
:meth:`~ampere.core.Model.evaluate` is still the contract surface and still
returns numpy containers — that boundary is where a gradient stops, which is
precisely why the native surface exists.

The declaration — the parameter names, their priors, the buffers, the fixed
line centres — is imported from :mod:`.model` rather than repeated, so the
three variants cannot drift anywhere except in their arithmetic.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.torch import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, to_numpy
from ampere.core import DTYPE, ChannelRequirements, Model, ModelResult, Spectrum

from .model import (
    FLUX_UNIT,
    LINE1_CENTRE,
    LINE1_WIDTH,
    LINE2_CENTRE,
    LINE2_WIDTH,
    REFERENCE_WAVELENGTH,
    WAVELENGTH_UNIT,
    _as_parameter,
    _channels,
    _declared,
    _to_micron,
)

__all__ = ["AbsorptionLines", "flux_tensor"]


def flux_tensor(
    wavelength: torch.Tensor,
    A: torch.Tensor,
    B: torch.Tensor,
    d1: torch.Tensor,
    d2: torch.Tensor,
) -> torch.Tensor:
    """:func:`examples.m2_misspecification.model.flux_at` in ``torch``.

    A free function, so ``torch.autograd.grad`` can be taken through it in a
    test without composing a whole problem first. Every operation is
    whole-tensor: nothing branches on a value and nothing indexes by one, which
    is what makes ``BATCHABLE = True`` below a fact rather than a hope.
    """
    continuum = A + B * (wavelength - REFERENCE_WAVELENGTH)
    first = torch.exp(-0.5 * ((wavelength - LINE1_CENTRE) / LINE1_WIDTH) ** 2)
    second = torch.exp(-0.5 * ((wavelength - LINE2_CENTRE) / LINE2_WIDTH) ** 2)
    return continuum * (1.0 - d1 * first - d2 * second)


class AbsorptionLines(Model):
    """The toy model on the torch path.

    Parameters
    ----------
    wavelength
        The model's own grid, micron.
    A, B, d1, d2
        As :class:`examples.m2_misspecification.model.AbsorptionLines`.
    channels
        Name of the emitted channel.
    dtype, device
        Threaded explicitly (``lowering.md`` §10.1): float64 on the CPU, and
        never taken from ``torch.get_default_dtype()``.
    """

    AXIS: ClassVar[str] = "spectral_axis"
    #: The four flags, declared rather than inherited. ``DIFFERENTIABLE`` is
    #: the claim :meth:`flux` makes good on; ``BATCHABLE`` is true because
    #: every operation in :func:`flux_tensor` is whole-tensor, so
    #: ``torch.func.vmap`` maps it as it maps any pure function. (A problem
    #: built from this model and ``QuasisepGP`` still reports
    #: ``batchable=False``: the flags aggregate conjunctively and that solver
    #: says ``False`` for its own reasons.)
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        wavelength: Any,
        *,
        A: Any = None,
        B: Any = None,
        d1: Any = None,
        d2: Any = None,
        channels: str | Sequence[str] = "default",
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        grid = _to_micron(wavelength)
        self.channels = _channels(channels)
        self.dtype = dtype
        self.device = device
        self.register_buffer("wavelength", grid, unit=WAVELENGTH_UNIT)
        for name, given in _declared(A=A, B=B, d1=d1, d2=d2).items():
            self.register_parameter(_as_parameter(name, given))
        self.templates: dict[str, Spectrum] = {}
        self.grids: dict[str, torch.Tensor] = {
            channel: as_tensor(grid, dtype=dtype, device=device) for channel in self.channels
        }

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each channel, once."""
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * WAVELENGTH_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
            self.grids[channel] = as_tensor(grid, dtype=self.dtype, device=self.device)
        return self

    # -- the native surface the realisation composes -------------------------

    def grid(self, channel: str) -> torch.Tensor:
        """The tensor coordinates *channel* is evaluated on, micron."""
        return self.grids[channel]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> torch.Tensor:
        """This model's flux on *channel*, Jy, as a differentiable tensor."""
        context = self._tensors(values)
        return flux_tensor(
            self.grid(channel), context["A"], context["B"], context["d1"], context["d2"]
        )

    def _tensors(self, values: Mapping[str, Any] | None) -> dict[str, torch.Tensor]:
        """This model's vocabulary as tensors, keeping any autograd graph intact."""
        return {
            name: as_tensor(value, dtype=self.dtype, device=self.device)
            for name, value in self.context(values).items()
        }

    # -- the contract surface ------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        """``ampere.core``'s contract: one container per channel, in numpy.

        The tensors come back to numpy here because they must: a
        :class:`~ampere.core.Spectrum` validates its axes in numpy and cannot
        hold a tensor carrying a graph. :meth:`flux` is the same computation
        without that boundary.
        """
        emitted = {}
        for channel in self.channels:
            flux = to_numpy(self.flux(channel, values)).astype(DTYPE, copy=False)
            template = self.templates.get(channel)
            emitted[channel] = (
                Spectrum(to_numpy(self.grid(channel)) * WAVELENGTH_UNIT, flux, unit=FLUX_UNIT)
                if template is None
                else template.with_values(flux)
            )
        return ModelResult(emitted)
