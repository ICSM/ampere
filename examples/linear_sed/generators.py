"""Synthetic data for ``examples.linear_sed``, and the tracked Spitzer IRS grid.

Two independent jobs live here. :func:`irs_wavelength_grids` reproduces, with
:mod:`astropy.io.fits` and :mod:`astropy.table` alone, the SL/LL wavelength
split ``ampere.legacy.data.Spectrum.fromFile(..., format="SPITZER-YAAAR")``
reads out of the tracked ``examples/test_data/cassis_yaaar_spcfw_14191360t.fits``
(``ampere/legacy/data/spectrum.py`` lines 660-694): the file is a primary
image HDU of shape ``(387, 16)`` whose fifteen named columns
(``COL01DEF``-``COL15DEF``) plus one dummy column become an
:class:`~astropy.table.Table`; rows with a non-finite value in any column are
dropped, the table is sorted by wavelength, and ``module`` splits it into SL
(0 or 1) and LL (2 or 3). :mod:`tests.examples.test_linear_sed` checks the two
grids this produces against the legacy reader's own, byte for byte.

:func:`deduplicated_grids` is the grid a v2 :class:`~ampere.core.Spectrum` can
actually hold. The CASSIS reduction repeats a handful of wavelength samples
where adjacent spectral orders overlap (17 of SL's 200 points, 10 of LL's
187) — legacy's :class:`~ampere.legacy.data.Spectrum` never validated
monotonicity and stored them as given, but v2's ``Spectrum`` "must be
strictly increasing, which is validated here rather than assumed. Legacy
ampere assumed it silently" (``ampere/core/results_schema.py``). Since every
flux value this module attaches to the grid is synthetic (the FITS file
supplies only the *sampling*, never real flux), the fix costs nothing beyond
the few repeated points themselves: :func:`numpy.unique` on an already-sorted
array both drops the duplicates and keeps the order, so ``deduplicated_grids``
is what :func:`synthetic_data` and :mod:`.linear_sed` actually build
containers on.

:func:`synthetic_data` is the same negotiate-then-``compile_for`` dance
:mod:`examples.sed_composition.generators` documents at length: three
instruments (a photometric catalogue and the two IRS chunks) bind the model's
one channel, so the model must be compiled onto their negotiated union before
it is evaluated, and each instrument then reads its own slice of the result.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate

__all__ = [
    "CALIBRATION_TRUTH",
    "IRS_FILE",
    "PHOTOMETRY_FRACTIONAL_NOISE",
    "SEED",
    "SPECTRUM_FRACTIONAL_NOISE",
    "TRUTH",
    "deduplicated_grids",
    "irs_wavelength_grids",
    "synthetic_data",
]

#: Reproducible everywhere this example is run: the noise draw, and (through
#: ``FittingProblem(..., seed=SEED)``) the sampler's own streams. The legacy
#: script (``examples/minimal_working_example.py``) seeds nothing -- it calls
#: ``np.random.randn`` unseeded -- so there is no legacy value to match; this
#: is simply a fixed seed in the v2 convention :mod:`examples.sed_composition`
#: and :mod:`examples.astrometry` also use.
SEED = 20260928

#: The legacy script's own truth (``examples/minimal_working_example.py``
#: lines 96-98: ``slope = 1.0``, ``intercept = 1.0``).
TRUTH: dict[str, float] = {"slope": 1.0, "intercept": 1.0}

#: Legacy injected no deliberate miscalibration between the photometry and
#: the spectrum (unlike :data:`examples.sed_composition.generators.
#: CALIBRATION_TRUTH`'s 2 %) -- both chunks' calibration truth is exactly 1.0.
CALIBRATION_TRUTH = 1.0

#: Fractional 1-sigma noise -- ``input_noise_phot``/``input_noise_spec`` in
#: the legacy script, both 0.1.
PHOTOMETRY_FRACTIONAL_NOISE = 0.1
SPECTRUM_FRACTIONAL_NOISE = 0.1

#: The tracked CASSIS file the legacy script reads through
#: ``Spectrum.fromFile(..., format="SPITZER-YAAAR")``.
IRS_FILE = Path(__file__).resolve().parents[1] / "test_data" / "cassis_yaaar_spcfw_14191360t.fits"


def irs_wavelength_grids(path: Path | str = IRS_FILE) -> tuple[np.ndarray, np.ndarray]:
    """The SL and LL wavelength grids, exactly as the legacy reader produces them.

    Byte-for-byte equal to ``ampere.legacy.data.Spectrum.fromFile(path,
    format="SPITZER-YAAAR")``'s two chunks -- checked once, in this module's
    own smoke test -- duplicated wavelength samples from overlapping IRS
    orders included. See :func:`deduplicated_grids` for the grid a v2
    ``Spectrum`` can hold.
    """
    with fits.open(path) as hdul:
        header = hdul[0].header
        names = [header[f"COL{index:02d}DEF"] for index in range(1, 16)] + ["DUMMY"]
        table = Table(hdul[0].data, names=names)
    table.rename_column("error (RMS+SYS)", "uncertainty")
    finite = np.logical_and.reduce([np.isfinite(column) for column in table.columns.values()])
    table = table[finite]
    table.sort(keys="wavelength")
    sl = np.logical_or(table["module"] == 0.0, table["module"] == 1.0)
    ll = np.logical_or(table["module"] == 2.0, table["module"] == 3.0)
    sl_wavelength = np.asarray(table["wavelength"][sl], dtype=float)
    ll_wavelength = np.asarray(table["wavelength"][ll], dtype=float)
    return sl_wavelength, ll_wavelength


def deduplicated_grids(path: Path | str = IRS_FILE) -> tuple[np.ndarray, np.ndarray]:
    """:func:`irs_wavelength_grids`, with repeated wavelength samples dropped.

    Both grids are already sorted, so :func:`numpy.unique` both removes the
    duplicates and leaves the order untouched.
    """
    sl_wavelength, ll_wavelength = irs_wavelength_grids(path)
    return np.unique(sl_wavelength), np.unique(ll_wavelength)


def synthetic_data(
    model: Model,
    catalogue: Instrument,
    sl_instrument: Instrument,
    ll_instrument: Instrument,
    *,
    seed: int = SEED,
) -> tuple[Spectrum, Spectrum, PhotometricPoints]:
    """Noisy observed containers for the catalogue and the two IRS chunks.

    Negotiates the three instruments' requirements, compiles *model* onto the
    union grid, evaluates it once at :data:`TRUTH`, and pushes the result
    through each instrument -- the same dance
    :func:`examples.sed_composition.generators.synthetic_data` documents at
    length, run here by hand because there is no fitting problem yet.

    Parameters
    ----------
    model
        The (uncompiled) model. Compiled onto the negotiated grid as a side
        effect, then handed back to the caller unchanged in every other
        respect -- :meth:`~ampere.core.transform.FittingProblem`'s own
        constructor calls ``compile_for`` again, on the same requirements,
        idempotently.
    catalogue, sl_instrument, ll_instrument
        The three instruments, already bound to the model's channel with
        distinct labels.
    seed
        Seeds the noise draw. The fitting problem built from the result gets
        its own seed independently.

    Returns
    -------
    tuple
        ``(observed_sl, observed_ll, observed_photometry)``.
    """
    rng = np.random.default_rng(seed)
    requirements = negotiate([catalogue, sl_instrument, ll_instrument])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)

    calibration: dict[str, Any] = {"calibration_scale.scale": CALIBRATION_TRUTH}
    sl_truth = sl_instrument(truth, calibration)
    ll_truth = ll_instrument(truth, calibration)
    photometry_truth = catalogue(truth)

    def _perturb(values: np.ndarray, fraction: float) -> tuple[np.ndarray, np.ndarray]:
        sigma = fraction * np.abs(values)
        return values + rng.normal(0.0, sigma), sigma

    sl_values, sl_sigma = _perturb(sl_truth.values, SPECTRUM_FRACTIONAL_NOISE)
    ll_values, ll_sigma = _perturb(ll_truth.values, SPECTRUM_FRACTIONAL_NOISE)
    phot_values, phot_sigma = _perturb(photometry_truth.values, PHOTOMETRY_FRACTIONAL_NOISE)

    observed_sl = Spectrum(
        sl_truth.spectral_axis.values * u.um, sl_values * u.Jy, uncertainty=sl_sigma * u.Jy
    )
    observed_ll = Spectrum(
        ll_truth.spectral_axis.values * u.um, ll_values * u.Jy, uncertainty=ll_sigma * u.Jy
    )
    observed_photometry = PhotometricPoints(
        photometry_truth.filters,
        photometry_truth.spectral_axis.values * u.um,
        phot_values * u.Jy,
        uncertainty=phot_sigma * u.Jy,
    )
    return observed_sl, observed_ll, observed_photometry
