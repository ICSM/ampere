"""Shared fits for the gross-case figures of ``reading_the_diagnostics.rst``.

Not a figure script (the leading underscore keeps it out of the page): the two
gross-case figures import it, so the IRS spectrum is fitted once per build.

It re-expresses the calls of ``docs/design/walkthroughs/persona_b2.py`` (the
user-journeys memo's Appendix B.2) at the docs budget: the independent fit
and the wide-kernel-prior GP, both started from ``optimise`` (the memo's
``independent/optimiser`` and ``gp-wide/optimiser`` rows; from the prior start
neither converges at a docs budget), each with
``EmceeEngine(p, walkers=24).run(DOC_STEPS, burn_in=DOC_BURN)`` instead of
``run(1500, burn_in=500)``, and ``gp_localisation`` over ``THIN``-thinned
draws instead of every draw (the memo's 130 s per run).
"""

from __future__ import annotations

import functools

import astropy.units as u
import numpy as np
import scipy.stats as st
from astropy.io import fits
from astropy.table import Table

from ampere.backends.reference import PowerLaw, Resample
from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Matern32,
    QuasisepGP,
    Spectrum,
)
from ampere.inference import EmceeEngine, optimise
from ampere.results import add_residuals, gp_localisation

DOC_STEPS, DOC_BURN, DOC_WALKERS, THIN = 300, 150, 24, 20
WIDE = (0.05, 5.0)  # persona_b2.PRIORS["wide"]: amplitude scale (Jy), length scale (um)
FILE = "examples/test_data/cassis_yaaar_spcfw_14191360t.fits"


@functools.cache
def spectrum() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The IRS spectrum as persona_b2.py reads it: finite, positive error, sorted, unique."""
    with fits.open(FILE) as handle:
        header = handle[0].header
        names = [header[f"COL{i:02d}DEF"] for i in range(1, 16)] + ["DUMMY"]
        table = Table(handle[0].data, names=names)
    ok = np.isfinite(table["flux"]) & np.isfinite(table["error (RMS+SYS)"])
    table = table[ok & (table["error (RMS+SYS)"] > 0)]
    table.sort("wavelength")
    _, keep = np.unique(table["wavelength"], return_index=True)
    table = table[keep]
    return (
        np.asarray(table["wavelength"], float),
        np.asarray(table["flux"], float),
        np.asarray(table["error (RMS+SYS)"], float),
    )


def _problem(likelihood: Likelihood) -> FittingProblem:
    wl, fl, er = spectrum()
    model = PowerLaw(
        np.geomspace(5.0, 40.0, 800), norm=st.loguniform(1e-3, 10.0), index=st.uniform(-3.0, 6.0)
    )
    data = Dataset(
        Spectrum(wl * u.um, fl * u.Jy, uncertainty=er * u.Jy),
        Instrument([Resample(wl)], channel="default", label="irs"),
        likelihood=likelihood,
    )
    return FittingProblem(model, {"irs": data}, seed=1)


def _kernel() -> Matern32:
    amplitude, length = WIDE
    return Matern32(
        st.halfnorm(scale=amplitude),
        st.halfnorm(scale=length),
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
        axes=("spectral_axis",),
    )


def _sample(problem: FittingProblem):
    """persona_b2.py's ``start == "optimiser"`` arm: a ball around ``optimise(problem)``."""
    engine = EmceeEngine(problem, walkers=DOC_WALKERS)
    start = engine.initial_positions(DOC_WALKERS, around=optimise(problem))
    return engine.run(DOC_STEPS, burn_in=DOC_BURN, initial=start)


@functools.cache
def independent():
    """The power law under independent noise: ``(tree, problem)`` with residuals added."""
    problem = _problem(Likelihood(GaussianFamily(), IndependentNoise()))
    tree = _sample(problem)
    return add_residuals(tree, problem, thin=THIN), problem


@functools.cache
def gp_wide():
    """The power law under a GP with the wide kernel prior, started from the optimiser."""
    likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(_kernel(), QuasisepGP()))
    problem = _problem(likelihood)
    tree = _sample(problem)
    return gp_localisation(tree, problem, thin=THIN), tree, problem
