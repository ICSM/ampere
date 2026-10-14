"""The torch twin of :class:`.model.EmissionLadder`: one line.

The twin base comes first, so its namespace and capability flags win.
Imported only when the study asks for torch, so ``import examples.line_fluxes``
needs neither extra.
"""

from __future__ import annotations

import ampere.backends.torch

from .model import EmissionLadder

__all__ = ["TorchEmissionLadder"]


class TorchEmissionLadder(ampere.backends.torch.PortableModel, EmissionLadder):
    """:class:`examples.line_fluxes.model.EmissionLadder` on torch."""
