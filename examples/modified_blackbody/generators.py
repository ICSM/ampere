"""Synthetic data for ``examples.modified_blackbody`` (W6.13 (3)).

The legacy script (``examples/examples_paper/modifiedblackbody.py``) makes
one dataset -- ten-band AKARI/Herschel photometry, 10 % noise -- from its
four-parameter model at one truth. This module is that generation, against
:mod:`ampere.core`'s negotiate-then-``compile_for`` dance
(:mod:`examples.sed_composition.generators` documents it at length): one
instrument, :meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`,
binds the model's one channel, so the model is compiled onto the negotiated
tabulation grid before it is evaluated.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from ampere.core import Instrument, Model, PhotometricPoints, negotiate

__all__ = [
    "FRACTIONAL_NOISE",
    "SEED",
    "TRUTH",
    "synthetic_data",
]

#: Reproducible everywhere this example is run. The legacy script seeds
#: nothing (``np.random.randn`` unseeded); this is a fixed seed in the v2
#: convention :mod:`examples.sed_composition` and :mod:`examples.linear_sed`
#: also use.
SEED = 20260928

#: The legacy script's own truth (``examples/examples_paper/
#: modifiedblackbody.py`` lines 68-71: ``t_true = 30.``, ``logm_true = 1``,
#: ``beta_true = -2``, ``d_true = 0.11``).
TRUTH: dict[str, float] = {"temperature": 30.0, "logmass": 1.0, "beta": -2.0, "distance": 0.11}

#: ``input_noise_phot`` in the legacy script.
FRACTIONAL_NOISE = 0.1


def synthetic_data(model: Model, instrument: Instrument, *, seed: int = SEED) -> PhotometricPoints:
    """A noisy ten-band photometric catalogue at :data:`TRUTH`.

    Negotiates *instrument*'s requirement, compiles *model* onto it,
    evaluates once at the truth, and perturbs by :data:`FRACTIONAL_NOISE`.
    """
    rng = np.random.default_rng(seed)
    requirements = negotiate([instrument])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)
    photometry_truth = instrument(truth)

    sigma = FRACTIONAL_NOISE * np.abs(photometry_truth.values)
    values = photometry_truth.values + rng.normal(0.0, sigma)
    return PhotometricPoints(
        photometry_truth.filters,
        photometry_truth.spectral_axis.values * u.um,
        values * u.Jy,
        uncertainty=sigma * u.Jy,
    )
