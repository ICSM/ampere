"""The PHOENIX emulator and star in ``jax``: the mirror of :mod:`.emulator_torch`.

The same arithmetic on :data:`JAX_OPS` (``jax.numpy``), float64 required: the
constructors call :func:`ampere.backends.jax.require_x64`, and
:func:`.phoenix_star.backend_module` calls
:func:`ampere.backends.jax.configure_x64` before any jax work, as
:mod:`examples.sed_composition` does.
"""

from __future__ import annotations

from typing import Any, ClassVar

import jax.numpy as jnp

from ampere.backends.jax import BACKEND, require_x64

from . import emulator as _reference
from . import phoenix_star as _star
from .emulator import Ops

__all__ = ["JAX_OPS", "PhoenixEmulator", "PhoenixStar"]


def _asarray(value: Any) -> Any:
    return jnp.asarray(value, dtype=jnp.float64)


def _asindex(value: Any) -> Any:
    return jnp.asarray(value, dtype=jnp.int64)


JAX_OPS = Ops(jnp, _asarray, _asindex)


class PhoenixEmulator(_reference.PhoenixEmulator):
    """:class:`examples.phoenix_star.emulator.PhoenixEmulator` on jax."""

    OPS: ClassVar[Ops] = JAX_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(**kwargs)


class PhoenixStar(_star.PhoenixStar):
    """:class:`examples.phoenix_star.phoenix_star.PhoenixStar` on jax."""

    OPS: ClassVar[Ops] = JAX_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(*args, **kwargs)
