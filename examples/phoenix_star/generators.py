"""Synthetic data for ``examples.phoenix_star`` -- the legacy truth, one seed.

The legacy script (``examples/examples_paper/phoenixstar.py``) is synthetic by
construction: it evaluates its own model at ``theta_true = [log10(4.68), 6750,
4.37, 1.0, 3.2]`` (lines 169-171), integrates it through eight filters with
5 % noise, and resamples it onto a Gaia-RVS grid with 1 % noise (lines
217-249). :func:`synthetic_data` does the same with the twin's own model and
instruments: negotiate the two instruments, compile the model onto the
result, evaluate once at :data:`TRUTH`, push the result through each
instrument, and perturb with fractional Gaussian noise at a fixed seed (the
legacy seeds nothing).
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate

__all__ = [
    "CALIBRATION_TRUTH",
    "PHOTOMETRY_FRACTIONAL_NOISE",
    "SEED",
    "SPECTRUM_FRACTIONAL_NOISE",
    "TRUTH",
    "synthetic_data",
]

SEED = 20260930

#: ``theta_true`` (``phoenixstar.py`` lines 169-171), by the twin's names.
TRUTH: dict[str, float] = {
    "log_luminosity": float(np.log10(4.68)),
    "teff": 6750.0,
    "logg": 4.37,
    "a_v": 1.0,
    "r_v": 3.2,
}

#: The legacy injects no miscalibration.
CALIBRATION_TRUTH = 1.0

#: ``input_noise_phot = 0.05`` and ``input_noise_spec = 0.01``.
PHOTOMETRY_FRACTIONAL_NOISE = 0.05
SPECTRUM_FRACTIONAL_NOISE = 0.01


def synthetic_data(
    model: Model,
    catalogue: Instrument,
    rvs: Instrument,
    *,
    truth: dict[str, float] | None = None,
    seed: int = SEED,
    photometry_noise: float = PHOTOMETRY_FRACTIONAL_NOISE,
    spectrum_noise: float = SPECTRUM_FRACTIONAL_NOISE,
) -> tuple[PhotometricPoints, Spectrum]:
    """``(observed_photometry, observed_rvs)`` at *truth* (default :data:`TRUTH`).

    *model* is compiled onto the two instruments' negotiated requirements as a
    side effect (idempotently -- the fitting problem does it again).
    """
    rng = np.random.default_rng(seed)
    compiled = model.compile_for(negotiate([catalogue, rvs]))
    produced = compiled(**(TRUTH if truth is None else truth))
    calibration: dict[str, Any] = {"calibration_scale.scale": CALIBRATION_TRUTH}
    photometry = catalogue(produced)
    spectrum = rvs(produced, calibration)

    def perturb(values: np.ndarray, fraction: float) -> tuple[np.ndarray, np.ndarray]:
        sigma = fraction * np.abs(values)
        return values + rng.normal(0.0, sigma), sigma

    phot_values, phot_sigma = perturb(photometry.values, photometry_noise)
    spec_values, spec_sigma = perturb(spectrum.values, spectrum_noise)
    observed_photometry = PhotometricPoints(
        photometry.filters,
        photometry.spectral_axis.values * u.um,
        phot_values * u.Jy,
        uncertainty=phot_sigma * u.Jy,
    )
    observed_rvs = Spectrum(
        spectrum.spectral_axis.values * u.um, spec_values * u.Jy, uncertainty=spec_sigma * u.Jy
    )
    return observed_photometry, observed_rvs
