"""``BinaryWithDisc``, bound to the torch backend.

See :mod:`.model`'s docstring for why one shared class, rather than a
per-backend reimplementation as :mod:`examples.m2_misspecification`'s toy
model has: this class has no arithmetic of its own to fork. Imports torch
only when called, so ``import examples.interferometry`` still needs neither
extra (:mod:`.study`'s ``model_for`` is the one place that reaches here).
"""

from __future__ import annotations

from typing import Any

from .model import BinaryWithDisc, binary_with_disc

__all__ = ["BinaryWithDisc", "binary_with_disc"]


def _itf() -> Any:
    from ampere.backends.torch import interferometry

    return interferometry


def build(x: Any, y: Any, *, disc_flux: float, disc_fwhm: float, **binary_kwargs: Any) -> BinaryWithDisc:
    """``BinaryWithDisc`` on the torch backend's own ``Binary``/``GaussianSource``."""
    return binary_with_disc(_itf(), x, y, disc_flux=disc_flux, disc_fwhm=disc_fwhm, **binary_kwargs)
