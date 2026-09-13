"""The synthetic reflex orbit both native backends' suites are held to.

One module rather than two copies, exactly
:mod:`tests.backends.interferometry_fixtures`'s own reason: the claim the
torch and jax astrometry suites make is the same claim — "this backend's
arithmetic is the reference backend's, and its gradient is real" — and a
claim written twice is a claim that eventually differs. Nothing here imports
a backend: every function takes the backend's ``astrometry`` module as its
first argument, exactly as :mod:`tests.conformance` takes a fixture, so the
same body runs under ``-e torch`` and under ``-e jax``.

The epochs and the truth are ``tests/conformance/test_astrometry.py``'s,
deliberately: a row here that disagrees with a row there is a disagreement
about arithmetic rather than about which orbit was observed.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    TimeSeries,
    negotiate,
)

#: The observation epochs, days — the conformance battery's own.
EPOCHS = np.array([0.0, 45.0, 130.0, 210.0, 305.0, 390.0, 460.0, 540.0, 615.0, 700.0, 780.0])

#: The synthetic orbit, at its true values.
ORBIT: dict[str, float] = {
    "pmra": 1.2,
    "pmdec": -0.6,
    "period": 400.0,
    "phase": 0.7,
    "amp_ra": 0.6,
    "amp_dec": 0.35,
}

#: Measurement uncertainty, mas.
SIGMA = 0.03

#: The free parameters of a fitted problem, and their true values.
TRUTH: dict[str, float] = {
    "model.pmra": ORBIT["pmra"],
    "model.pmdec": ORBIT["pmdec"],
}

#: The seed every synthetic dataset here is drawn at.
SEED = 20260913


def observed_channel(values: np.ndarray | None = None) -> TimeSeries:
    """A ``TimeSeries`` on :data:`EPOCHS`, with a uniform sigma."""
    filled = np.zeros(EPOCHS.size) if values is None else values
    return TimeSeries(
        EPOCHS * u.day, filled * u.mas, uncertainty=np.full(EPOCHS.size, SIGMA) * u.mas
    )


def truth_model(astro: Any) -> Any:
    """The reflex orbit at its true values, on its own default epoch grid."""
    return astro.ReflexOrbit(EPOCHS * u.day, **ORBIT)


def chain(astro: Any, channel: str, label: str) -> Instrument:
    """The one-step epoch-sampling instrument for one coordinate."""
    return Instrument([astro.EpochSample(EPOCHS)], channel=channel, label=label)


def synthetic(astro: Any, *, seed: int = SEED) -> tuple[Any, TimeSeries, TimeSeries]:
    """Negotiate, compile, evaluate and observe: the data a fit is given.

    Returned with the compiled truth model, because its negotiated epochs are
    what the *fitted* model must be built on.
    """
    ra_instrument = chain(astro, "ra", "astrom_ra")
    dec_instrument = chain(astro, "dec", "astrom_dec")
    compiled = truth_model(astro).compile_for(negotiate([ra_instrument, dec_instrument]))
    result = compiled.evaluate()
    clean_ra = np.asarray(ra_instrument(result).values)
    clean_dec = np.asarray(dec_instrument(result).values)
    rng = np.random.default_rng(seed)
    noisy_ra = clean_ra + SIGMA * rng.standard_normal(clean_ra.size)
    noisy_dec = clean_dec + SIGMA * rng.standard_normal(clean_dec.size)
    return compiled, observed_channel(noisy_ra), observed_channel(noisy_dec)


def fitted_model(astro: Any, compiled: Any) -> Any:
    """The reflex orbit with every orbital parameter free, on the negotiated epochs."""
    return astro.ReflexOrbit(
        compiled.buffers["time"].value * u.day,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        period=st.norm(400.0, 30.0),
        phase=st.uniform(0.0, 2.0 * np.pi),
        amp_ra=st.uniform(0.0, 2.0),
        amp_dec=st.uniform(0.0, 2.0),
    )


def two_channel_problem(
    backend: Any, astro: Any, *, gp: bool = False, seed: int = SEED
) -> FittingProblem:
    """RA and Dec from one orbit, on one backend.

    Parameters
    ----------
    backend
        The backend package (``ampere.backends.torch`` or ``.jax``), for its
        noise model, kernel and solver.
    astro
        That backend's ``astrometry`` module.
    gp
        Whether the two channels carry the flexible likelihood
        (``GaussianProcessNoise(Matern32(axes=("time",)), QuasisepGP())``) or
        independent Gaussian noise.
    seed
        The noise realisation.
    """
    compiled, observed_ra, observed_dec = synthetic(astro, seed=seed)
    ra_instrument = chain(astro, "ra", "astrom_ra")
    dec_instrument = chain(astro, "dec", "astrom_dec")
    if gp:
        kernel = backend.Matern32(0.03, 120.0, axes=("time",))
        ra_noise: Any = backend.GaussianProcessNoise(kernel, backend.QuasisepGP())
        dec_noise: Any = backend.GaussianProcessNoise(kernel, backend.QuasisepGP())
    else:
        ra_noise = backend.IndependentNoise()
        dec_noise = backend.IndependentNoise()
    datasets = DatasetCollection(
        {
            "ra": Dataset(
                observed_ra, ra_instrument, likelihood=Likelihood(GaussianFamily(), ra_noise), label="ra"
            ),
            "dec": Dataset(
                observed_dec,
                dec_instrument,
                likelihood=Likelihood(GaussianFamily(), dec_noise),
                label="dec",
            ),
        }
    )
    return FittingProblem(fitted_model(astro, compiled), datasets, seed=SEED)
