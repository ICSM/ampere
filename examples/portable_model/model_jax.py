"""The jax twin of :class:`.model.LinearModel`: one line.

The twin base comes first, so its namespace (:class:`ampere.backends.jax.JaxOps`)
and capability flags win, and the x64 check runs before the model is built;
call :func:`ampere.backends.jax.configure_x64` before any jax work.
"""

from __future__ import annotations

import ampere.backends.jax

from .model import LinearModel

__all__ = ["JaxLinearModel"]


class JaxLinearModel(ampere.backends.jax.PortableModel, LinearModel):
    """:class:`examples.portable_model.model.LinearModel` on jax."""
