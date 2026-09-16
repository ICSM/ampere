"""Synthetic data for the astrometry example: negotiate, then compile, then evaluate.

Exactly :mod:`examples.sed_composition.generators`'s own dance, cited there in
full: an instrument bound to a model that has not adopted the negotiated grid
is refused, so a synthetic-data generator does what
:class:`~ampere.core.dataset.FittingProblem` does in its own constructor
(negotiate, ``compile_for``, evaluate) rather than calling the model on its
own fallback epoch grid.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from ampere.core import Instrument, Model, TimeSeries, negotiate

__all__ = [
    "EPOCHS",
    "SEED",
    "SIGMA",
    "TRUTH",
    "synthetic_data",
]

#: Reproducible everywhere this example is run.
SEED = 20260913

#: Observation epochs, days — irregular (real astrometric scheduling is),
#: spanning about three periods (the truth's ``period`` is 400 days) so the
#: wobble and the linear drift are both constrained and the period search is
#: not left free to alias against a baseline shorter than the true cycle.
EPOCHS = np.array(
    [
        0.0,
        35.0,
        80.0,
        130.0,
        175.0,
        210.0,
        255.0,
        305.0,
        350.0,
        390.0,
        430.0,
        460.0,
        505.0,
        540.0,
        580.0,
        615.0,
        660.0,
        700.0,
        740.0,
        780.0,
        830.0,
        860.0,
        900.0,
        940.0,
        985.0,
        1020.0,
        1060.0,
        1100.0,
    ]
)

#: The injected truth: a modest proper motion and a reflex wobble from an
#: unseen companion, at a period and semi-amplitude a decade of epochs can
#: actually constrain.
TRUTH: dict[str, float] = {
    "pmra": 1.2,
    "pmdec": -0.6,
    "period": 400.0,
    "phase": 0.7,
    "amp_ra": 0.6,
    "amp_dec": 0.35,
}

#: Per-epoch measurement uncertainty, mas — a modest ground-based-astrometry
#: number, small enough that the wobble (a few tenths of a mas) is detectable
#: over a decade of epochs.
SIGMA = 0.03


def synthetic_data(
    model: Model,
    ra_instrument: Instrument,
    dec_instrument: Instrument,
    *,
    seed: int = SEED,
) -> tuple[TimeSeries, TimeSeries]:
    """Noisy observed ``(ra, dec)`` containers, from :data:`TRUTH`.

    Negotiates the two instruments' requirements, compiles *model* onto the
    union epochs (independently per channel — see
    :meth:`~ampere.backends.reference.astrometry.ReflexOrbit.compile_for`),
    evaluates it once at :data:`TRUTH`, and pushes the result through each
    instrument. ``model`` is mutated in place (``compile_for`` adopts the
    negotiated grid) and handed back to the caller unchanged in every other
    respect, so it can be passed on to
    :class:`~ampere.core.dataset.FittingProblem` as-is.
    """
    rng = np.random.default_rng(seed)
    requirements = negotiate([ra_instrument, dec_instrument])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)

    ra_truth = ra_instrument(truth)
    dec_truth = dec_instrument(truth)

    ra_noisy = np.asarray(ra_truth.values) + rng.normal(0.0, SIGMA, ra_truth.values.shape)
    dec_noisy = np.asarray(dec_truth.values) + rng.normal(0.0, SIGMA, dec_truth.values.shape)

    observed_ra = TimeSeries(
        ra_truth.time.values * u.day,
        ra_noisy * u.mas,
        uncertainty=np.full(ra_noisy.shape, SIGMA) * u.mas,
    )
    observed_dec = TimeSeries(
        dec_truth.time.values * u.day,
        dec_noisy * u.mas,
        uncertainty=np.full(dec_noisy.shape, SIGMA) * u.mas,
    )
    return observed_ra, observed_dec
