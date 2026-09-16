"""``SourceWithBackground``, bound to the jax backend.

See :mod:`.model`'s docstring for why one shared class rather than a
per-backend reimplementation. Imports jax only when called — and calls
``configure_x64()`` first, because ``lowering.md`` §10.2 makes turning x64 on
the *application's* job and never an import side effect of ``ampere``; for this
package the application is this function.
"""

from __future__ import annotations

from typing import Any

from .model import SourceWithBackground, source_with_background

__all__ = ["SourceWithBackground", "build", "source_with_background"]


def _itf() -> Any:
    from ampere.backends.jax import configure_x64, interferometry

    configure_x64()
    return interferometry


def build(x: Any, y: Any, **kwargs: Any) -> SourceWithBackground:
    """``SourceWithBackground`` on the jax backend's own ``GaussianSource``."""
    return source_with_background(_itf(), x, y, **kwargs)
