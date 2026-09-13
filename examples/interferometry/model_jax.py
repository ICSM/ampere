"""``BinaryWithDisc``, bound to the jax backend.

See :mod:`.model`'s docstring for why one shared class, rather than a
per-backend reimplementation as :mod:`examples.m2_misspecification`'s toy
model has: this class has no arithmetic of its own to fork. Imports jax only
when called, so ``import examples.interferometry`` still needs neither extra
(:mod:`.study`'s ``model_for`` is the one place that reaches here); the
caller is responsible for calling ``ampere.backends.jax.configure_x64()``
first, exactly as :mod:`examples.m2_misspecification.study` requires of its
own jax path.
"""

from __future__ import annotations

from typing import Any

from .model import BinaryWithDisc, binary_with_disc

__all__ = ["BinaryWithDisc", "binary_with_disc"]


def _itf() -> Any:
    from ampere.backends.jax import interferometry

    return interferometry


def build(x: Any, y: Any, *, disc_flux: float, disc_fwhm: float, **binary_kwargs: Any) -> BinaryWithDisc:
    """``BinaryWithDisc`` on the jax backend's own ``Binary``/``GaussianSource``."""
    return binary_with_disc(_itf(), x, y, disc_flux=disc_flux, disc_fwhm=disc_fwhm, **binary_kwargs)
