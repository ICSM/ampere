"""The torch twin of :class:`.model.LinearModel`: one line.

The twin base comes first, so its namespace
(:class:`ampere.backends.torch.TorchOps`, float64 on the CPU) and capability
flags win.
"""

from __future__ import annotations

import ampere.backends.torch

from .model import LinearModel

__all__ = ["TorchLinearModel"]


class TorchLinearModel(ampere.backends.torch.PortableModel, LinearModel):
    """:class:`examples.portable_model.model.LinearModel` on torch."""
