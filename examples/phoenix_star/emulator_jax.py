"""The PHOENIX emulator and star in ``jax``: the mirror of :mod:`.emulator_torch`.

The same arithmetic on :data:`JAX_OPS` (:class:`ampere.backends.jax.JaxOps`), float64 required: the
constructors call :func:`ampere.backends.jax.require_x64`, and
:func:`.phoenix_star.backend_module` calls
:func:`ampere.backends.jax.configure_x64` before any jax work, as
:mod:`examples.sed_composition` does.
"""

from __future__ import annotations

from typing import Any, ClassVar

from ampere.backends.jax import BACKEND, JaxOps, require_x64
from ampere.core import ArrayOps

from . import emulator as _reference
from . import phoenix_star as _star

__all__ = ["JAX_OPS", "PhoenixEmulator", "PhoenixStar"]


#: The public :class:`~ampere.backends.jax.JaxOps` (W7.13): float64, the
#: default device.
JAX_OPS: ArrayOps = JaxOps()


class PhoenixEmulator(_reference.PhoenixEmulator):
    """:class:`examples.phoenix_star.emulator.PhoenixEmulator` on jax."""

    OPS: ClassVar[ArrayOps] = JAX_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(**kwargs)


class PhoenixStar(_star.PhoenixStar):
    """:class:`examples.phoenix_star.phoenix_star.PhoenixStar` on jax."""

    OPS: ClassVar[ArrayOps] = JAX_OPS
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(*args, **kwargs)
