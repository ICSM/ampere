"""The jax twin of :class:`.model.EmissionLadder`: one line.

The twin base comes first, so its namespace and capability flags win and its
x64 check runs before the model is built; call
:func:`ampere.backends.jax.configure_x64` before any jax work (the study's
``model_for`` does).
"""

from __future__ import annotations

import ampere.backends.jax

from .model import EmissionLadder

__all__ = ["JaxEmissionLadder"]


class JaxEmissionLadder(ampere.backends.jax.PortableModel, EmissionLadder):
    """:class:`examples.line_fluxes.model.EmissionLadder` on jax."""
