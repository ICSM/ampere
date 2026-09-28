"""Synthetic data for the photometry-and-spectra example: push the truth through
the *negotiated* grid, not the model's own one.

Extends :mod:`examples.sed_composition.generators`'s dance (read that
module's docstring for why a synthetic-data generator cannot just call the
model on its declared grid) to three instruments instead of two:
:func:`~ampere.core.negotiate`, then
:meth:`~ampere.core.transform.Model.compile_for`, *then* evaluate — the same
:class:`~ampere.core.dataset.FittingProblem` does in its own constructor.

Two things this module's truth carries that the sibling's does not, both the
point of :mod:`.photometry_spectra`'s tutorial:

**Two different calibration factors.** :data:`CALIBRATION_TRUTH` gives the
short-wavelength spectrograph (``"sl"``) and the long-wavelength one
(``"ll"``) different miscalibrations, 0.92 and 1.08, so a fit that *ties*
the two factors together has a real truth to fail to recover — the whole
point of showing what tying does.

**A smooth residual the model cannot express.** :func:`synthetic_data`
injects a Gaussian bump into the ``"ll"`` spectrum only — 10 % amplitude,
4 micron width, centred at 25 micron, the classic "the dust is not one
blackbody" feature. No modified blackbody can fit a bump like that; the
tutorial page's point, measured rather than asserted, is what absorbs it
instead under :class:`~ampere.core.IndependentNoise` (the calibration factor
and the temperature, biased) versus under
:class:`~ampere.core.GaussianProcessNoise` (the physical parameters and the
calibration factors recover, and the GP's posterior mean localises the
bump).
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate

__all__ = [
    "BUMP_CENTRE",
    "BUMP_FRACTION",
    "BUMP_WIDTH",
    "CALIBRATION_TRUTH",
    "PHOTOMETRY_FRACTIONAL_NOISE",
    "SEED",
    "SPECTRUM_FRACTIONAL_NOISE",
    "TRUTH",
    "synthetic_data",
]

#: Reproducible everywhere this example is run: the noise draw, and (through
#: ``FittingProblem(..., seed=SEED)``) the sampler's own streams.
SEED = 20260928

#: The injected physical truth — the same cool, optically thin dust greybody
#: :data:`examples.sed_composition.generators.TRUTH` uses, since this
#: example's point is the third observation and the calibration story, not a
#: different source.
TRUTH: dict[str, float] = {"temperature": 180.0, "beta": 1.6, "scale": 1.0e-15}

#: The two spectrographs' calibration nuisance parameters
#: (``sl.instrument.calibration_scale.scale`` and
#: ``ll.instrument.calibration_scale.scale`` once bound into a problem) —
#: deliberately different from each other (2026-09-28 ruling on W6.2), so an
#: untied fit has two truths to recover and a tied one has a genuine conflict
#: to split the difference on.
CALIBRATION_TRUTH: dict[str, float] = {"sl": 0.92, "ll": 1.08}

#: The injected bump's centre, micron -- inside the "ll" spectrograph's 14-38
#: micron coverage.
BUMP_CENTRE = 25.0

#: The injected bump's width, micron -- taken as the Gaussian's standard
#: deviation (the ruling says "width" without specifying which; the
#: standard-deviation reading is the one that keeps the bump's 1-sigma
#: extent, 21-29 micron, comfortably inside the "ll" spectrograph's coverage).
BUMP_WIDTH = 4.0

#: The injected bump's peak amplitude, as a fraction of the noiseless "ll"
#: flux at each wavelength (so the bump scales with the local continuum
#: rather than being a fixed flux).
BUMP_FRACTION = 0.10

#: Fractional 1-sigma noise, applied to each instrument's own truth values --
#: the same fractions :mod:`examples.sed_composition.generators` uses.
SPECTRUM_FRACTIONAL_NOISE = 0.03
PHOTOMETRY_FRACTIONAL_NOISE = 0.05


def synthetic_data(
    model: Model,
    sl: Instrument,
    ll: Instrument,
    camera: Instrument,
    *,
    seed: int = SEED,
) -> tuple[Spectrum, Spectrum, PhotometricPoints]:
    """Noisy observed containers for *sl*, *ll* and *camera*, from ``TRUTH``.

    Negotiates the three instruments' requirements, compiles *model* onto the
    union grid, evaluates it once at ``TRUTH``, and pushes the result through
    each instrument with its own calibration truth -- exactly what
    :class:`~ampere.core.dataset.FittingProblem` does when it is built, run
    here by hand because there is no fitting problem yet: this *produces* the
    data one goes into. ``model`` is mutated in place
    (:meth:`~ampere.core.transform.Model.compile_for` adopts the negotiated
    grid as its template) and handed back to the caller unchanged in every
    other respect, so it can be passed on to
    :class:`~ampere.core.dataset.FittingProblem` as-is.

    Parameters
    ----------
    model
        The (uncompiled) model. Compiled onto the union grid as a side effect.
    sl, ll, camera
        The three instruments, already bound to the model's channel with
        distinct labels.
    seed
        Seeds the noise draw *and* the bump's placement (the bump itself is
        deterministic given :data:`BUMP_CENTRE`/:data:`BUMP_WIDTH`/
        :data:`BUMP_FRACTION`; only the noise realisation is random). The
        fitting problem built from the result gets its own seed
        independently -- this is only about the data.

    Returns
    -------
    tuple
        ``(observed_sl, observed_ll, observed_photometry)``, ready to bind
        into :class:`~ampere.core.dataset.Dataset`\\ s.
    """
    rng = np.random.default_rng(seed)
    requirements = negotiate([sl, ll, camera])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)

    # calibration_scale.scale is each spectrograph's own step's name
    # (relative to the instrument, not yet qualified by a dataset label --
    # that qualification only exists once the instrument is bound into a
    # Dataset -- see examples.sed_composition.generators for the same point).
    sl_truth = sl(truth, {"calibration_scale.scale": CALIBRATION_TRUTH["sl"]})
    ll_truth = ll(truth, {"calibration_scale.scale": CALIBRATION_TRUTH["ll"]})
    photometry_truth = camera(truth)

    # The complement: a smooth feature no modified blackbody can express,
    # injected into "ll" only, after its own calibration factor has already
    # been applied.
    wavelength = ll_truth.spectral_axis.values
    bump = BUMP_FRACTION * np.exp(-0.5 * ((wavelength - BUMP_CENTRE) / BUMP_WIDTH) ** 2)
    ll_truth = ll_truth.with_values(ll_truth.values * (1.0 + bump))

    sl_sigma = SPECTRUM_FRACTIONAL_NOISE * sl_truth.values
    ll_sigma = SPECTRUM_FRACTIONAL_NOISE * ll_truth.values
    photometry_sigma = PHOTOMETRY_FRACTIONAL_NOISE * photometry_truth.values

    observed_sl = Spectrum(
        sl_truth.spectral_axis.values * u.um,
        (sl_truth.values + rng.normal(0.0, sl_sigma)) * u.Jy,
        uncertainty=sl_sigma * u.Jy,
    )
    observed_ll = Spectrum(
        ll_truth.spectral_axis.values * u.um,
        (ll_truth.values + rng.normal(0.0, ll_sigma)) * u.Jy,
        uncertainty=ll_sigma * u.Jy,
    )
    observed_photometry = PhotometricPoints(
        photometry_truth.filters,
        photometry_truth.spectral_axis.values * u.um,
        (photometry_truth.values + rng.normal(0.0, photometry_sigma)) * u.Jy,
        uncertainty=photometry_sigma * u.Jy,
    )
    return observed_sl, observed_ll, observed_photometry
