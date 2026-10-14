"""Spectroscopy in ampere: the observable's front door (W7.3).

One import for everything a fit to spectra and photometry needs on the
reference (numpy) path: the two container kinds, the likelihood families they
are scored with, the reference backend's instrument steps, and the reader that
fills a ``Spectrum`` from a JWST pipeline product::

    from ampere.spectroscopy import (
        read_jwst, Spectrum, LSFConvolution, Resample, GaussianFamily,
    )

**What this package is.** The front-door pattern of ``ampere.interferometry``
(W6.12, the placement memo's option D, ``docs/design/phase4_placement_memo.md``
§2.2 to 2.4) applied to spectroscopy: a user-facing home *layered over* the
placement that the contracts fix, not a place code moves to. The kinds and the
families stay in :mod:`ampere.core`, where every backend reads them; the steps
stay in each backend, under the backend rule (one name per backend, a native
problem composed entirely from one backend's pieces). What is defined here and
nowhere else is the reader, :func:`read_jwst` (:mod:`ampere.spectroscopy.jwst`).

**What is re-exported, and what is not.** The kinds :class:`~ampere.core.Spectrum`
and :class:`~ampere.core.PhotometricPoints`; the four families a flux is scored
with (:class:`~ampere.core.GaussianFamily`, :class:`~ampere.core.StudentTFamily`,
:class:`~ampere.core.CauchyFamily`, :class:`~ampere.core.PoissonFamily` — the
complex, Rice and von Mises families are interferometric and live behind
``ampere.interferometry``); the reference steps that act on a spectrum
(:class:`~ampere.backends.reference.Resample`,
:class:`~ampere.backends.reference.LSFConvolution`,
:class:`~ampere.backends.reference.CalibrationScale`,
:class:`~ampere.backends.reference.SyntheticPhotometry`); and the reader with
its error. The noise models (``IndependentNoise``, ``GaussianProcessNoise``)
and the kernels are not, following the interferometry door's line: they belong
to every observable alike, and a fit imports them from :mod:`ampere.core`. Nor
are the source models: a spectral model is the user's choice of physics, not
part of the observable.

**Torch and jax users** compose from their own backend's pieces, as before.
Nothing native is re-exported here, because a reference step in a torch
problem is refused at composition, and a front door that offered one would be
offering a refusal. The ``Spectrum`` :func:`read_jwst` returns is
backend-neutral and feeds any backend's steps unchanged.

``import ampere`` imports none of this; the package is imported on its own.
"""

from __future__ import annotations

from ampere.backends.reference import (
    CalibrationScale,
    LSFConvolution,
    Resample,
    SyntheticPhotometry,
)
from ampere.core import (
    CauchyFamily,
    GaussianFamily,
    PhotometricPoints,
    PoissonFamily,
    Spectrum,
    StudentTFamily,
)

from .jwst import JWSTError, read_jwst

__all__ = [
    "CalibrationScale",
    "CauchyFamily",
    "GaussianFamily",
    "JWSTError",
    "LSFConvolution",
    "PhotometricPoints",
    "PoissonFamily",
    "Resample",
    "Spectrum",
    "StudentTFamily",
    "SyntheticPhotometry",
    "read_jwst",
]
