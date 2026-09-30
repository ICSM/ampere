"""The PHOENIX emulator and star in ``torch``: the same arithmetic, differentiably.

Imported only when a torch problem is built (:func:`.phoenix_star.star_class`),
so ``import examples.phoenix_star`` never needs torch -- the discipline
:mod:`examples.m2_misspecification.model_torch` follows. As there, what makes
these classes differentiable is the native surface ``grid(channel)`` /
``flux(channel, values)``, which :mod:`ampere.backends.torch` finds by
``hasattr``; :meth:`~ampere.core.Model.evaluate` still returns numpy
containers.

Nothing is re-declared: the parameters, buffers, channels and the committed
arrays come from :mod:`.emulator` and :mod:`.phoenix_star`, and the arithmetic
(:func:`.emulator.log_shape`, :func:`.emulator.star_log_flux`) is the very same
function run on :data:`TORCH_OPS` -- ``torch.exp``, ``torch.sum``,
``torch.stack`` and ``@`` in place of numpy's -- so the three backends cannot
drift anywhere but in floating-point rounding.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.torch import BACKEND, as_tensor

from . import emulator as _reference
from . import phoenix_star as _star
from .emulator import Ops

__all__ = ["TORCH_OPS", "PhoenixEmulator", "PhoenixStar"]


def _asindex(value: Any) -> torch.Tensor:
    return torch.as_tensor(np.asarray(value, dtype=np.int64))


#: float64 on the CPU, never ``torch.get_default_dtype()`` (``lowering.md`` 10.1).
TORCH_OPS = Ops(torch, as_tensor, _asindex)


class PhoenixEmulator(_reference.PhoenixEmulator):
    """:class:`examples.phoenix_star.emulator.PhoenixEmulator` on torch."""

    OPS: ClassVar[Ops] = TORCH_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND


class PhoenixStar(_star.PhoenixStar):
    """:class:`examples.phoenix_star.phoenix_star.PhoenixStar` on torch."""

    OPS: ClassVar[Ops] = TORCH_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND
