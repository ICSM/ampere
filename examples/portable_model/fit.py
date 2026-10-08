"""The guide's fit: the quickstart's line, one spectrograph and the photometry, on any backend.

The composition is the quickstart's (``docs/source/notebooks/quickstart.ipynb``)
reduced to two datasets -- the WISE W1 / MIPS 70 catalogue and the IRS
short-low chunk -- so that the guide runs in the documentation build. What
differs between backends is **one module**: :func:`build_problem` takes its
model twin from :func:`model_class` and its instrument steps and noise model
from :func:`backend_module`, and nothing else changes. The data are made once,
on the reference path, so every backend fits the same numbers.

:func:`fit_emcee` and :func:`fit_nuts` are the guide's two engine cells, and
``tests/examples/test_portable_model.py`` runs exactly these functions, so the
NUTS result the guide quotes is the one CI's jax leg checks.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Instrument,
    Likelihood,
    PhotometricPoints,
    Spectrum,
    negotiate,
)
from examples.linear_sed.generators import deduplicated_grids

from .model import LinearModel

__all__ = [
    "BACKENDS",
    "GRID",
    "NUTS_BUDGET",
    "SEED",
    "TRUTH",
    "backend_module",
    "build_problem",
    "fit_emcee",
    "fit_nuts",
    "model_class",
    "synthetic_data",
]

#: The three backends the guide names.
BACKENDS: tuple[str, ...] = ("reference", "jax", "torch")
#: The quickstart's model grid, micron.
GRID = 10 ** np.linspace(0.0, 1.9, 2000)
#: The quickstart's truth: slope 1, intercept 1, the calibration factor 1.
TRUTH: dict[str, float] = {"slope": 1.0, "intercept": 1.0}
SEED = 20260928
FILTERS = ("WISE_RSR_W1", "SPITZER_MIPS_70")
PHOTOMETRY_TABULATION = np.geomspace(1.0, 79.0, 500)
CALIBRATION_PRIOR = st.lognorm(0.0025, scale=1.0)
#: The guide's NUTS budget, also the one its test runs.
NUTS_BUDGET: dict[str, int] = {"draws": 500, "warmup": 500, "chains": 2}


def backend_module(name: str) -> ModuleType:
    """``ampere.backends.<name>``, with jax's float64 switched on first."""
    module = importlib.import_module(f"ampere.backends.{name}")
    if name == "jax":
        module.configure_x64()
    return module


def model_class(name: str) -> type[LinearModel]:
    """The guide's model on *name*: the one source, or its one-line native twin."""
    if name == "reference":
        return LinearModel
    if name == "jax":
        from .model_jax import JaxLinearModel

        return JaxLinearModel
    if name == "torch":
        from .model_torch import TorchLinearModel

        return TorchLinearModel
    raise ValueError(f"unknown backend {name!r}; the guide names {BACKENDS}.")


def _instruments(pieces: ModuleType, sl_grid: np.ndarray) -> tuple[Instrument, Instrument]:
    catalogue = Instrument(
        [pieces.SyntheticPhotometry.from_library(list(FILTERS), PHOTOMETRY_TABULATION)],
        channel="sed",
        label="catalogue",
    )
    sl = Instrument(
        [pieces.Resample(sl_grid), pieces.CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed",
        label="sl",
    )
    return catalogue, sl


def synthetic_data(seed: int = SEED) -> tuple[PhotometricPoints, Spectrum]:
    """The catalogue and the SL spectrum at :data:`TRUTH`, with 10 per cent noise, made on numpy."""
    sl_grid, _ = deduplicated_grids()
    catalogue, sl = _instruments(backend_module("reference"), sl_grid)
    model = LinearModel(GRID)
    truth = model.compile_for(negotiate([catalogue, sl]))(**TRUTH)
    rng = np.random.default_rng(seed)
    clean_photometry = catalogue(truth)
    clean_sl = sl(truth, {"calibration_scale.scale": 1.0})

    def noisy(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        sigma = 0.1 * np.abs(values)
        return values + rng.normal(0.0, sigma), sigma

    values, sigma = noisy(np.asarray(clean_photometry.values))
    photometry = PhotometricPoints(
        clean_photometry.filters,
        clean_photometry.spectral_axis.values * u.um,
        values * u.Jy,
        uncertainty=sigma * u.Jy,
    )
    values, sigma = noisy(np.asarray(clean_sl.values))
    spectrum = Spectrum(
        clean_sl.spectral_axis.values * u.um, values * u.Jy, uncertainty=sigma * u.Jy
    )
    return photometry, spectrum


def build_problem(backend: str = "reference", *, seed: int = SEED) -> FittingProblem:
    """The guide's problem on *backend*: its model twin, its instrument steps, its noise model."""
    pieces = backend_module(backend)
    photometry, spectrum = synthetic_data(seed)
    sl_grid, _ = deduplicated_grids()
    catalogue, sl = _instruments(pieces, sl_grid)
    # The reference path's noise models are ampere.core's own; each native
    # backend ships its own IndependentNoise.
    noise = getattr(pieces, "IndependentNoise", IndependentNoise)
    likelihood = Likelihood(GaussianFamily(), noise())
    datasets = DatasetCollection(
        {
            "catalogue": Dataset(photometry, catalogue, likelihood=likelihood),
            "sl": Dataset(spectrum, sl, likelihood=likelihood),
        }
    )
    return FittingProblem(model_class(backend)(GRID), datasets, seed=seed)


def fit_emcee(
    problem: FittingProblem, *, walkers: int = 24, steps: int = 1000, burn_in: int = 500
) -> Any:
    """The guide's emcee cell: the quickstart's 24 walkers, 1000 steps, the optimiser's start."""
    from ampere.inference import EmceeEngine

    return EmceeEngine(problem, walkers=walkers).run(steps, burn_in=burn_in)


def fit_nuts(problem: FittingProblem, **budget: int) -> Any:
    """The guide's NUTS cell, at :data:`NUTS_BUDGET` unless *budget* overrides it."""
    from ampere.inference import NUTSEngine

    settings = {**NUTS_BUDGET, **budget}
    return NUTSEngine(problem).run(
        settings["draws"], warmup=settings["warmup"], chains=settings["chains"]
    )
