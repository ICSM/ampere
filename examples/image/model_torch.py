"""``SourceWithBackground``, bound to the torch backend.

See :mod:`.model`'s docstring for why one shared class rather than a
per-backend reimplementation: this class has no arithmetic of its own to fork,
only a grid negotiation and an addition. Imports torch only when called, so
``import examples.image`` still needs neither extra.
"""

from __future__ import annotations

from typing import Any

from .model import SourceWithBackground, source_with_background

__all__ = ["SourceWithBackground", "build", "source_with_background"]


def _itf() -> Any:
    from ampere.backends.torch import interferometry

    return interferometry


def build(x: Any, y: Any, **kwargs: Any) -> SourceWithBackground:
    """``SourceWithBackground`` on the torch backend's own ``GaussianSource``."""
    return source_with_background(_itf(), x, y, **kwargs)
