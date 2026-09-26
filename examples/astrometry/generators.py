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

from ampere.core import Instrument, Matern32, Model, RotationCoupling, TimeSeries, negotiate

__all__ = [
    "EPOCHS",
    "JOINT_LENGTH_SCALE",
    "JOINT_TRUTH",
    "SEED",
    "SIGMA",
    "SIGMA_SPREAD",
    "TRUTH",
    "channel_sigmas",
    "coupling_matrix",
    "marginal_amplitudes",
    "synthetic_data",
    "synthetic_joint_data",
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


# ---------------------------------------------------------------------------
# The injected correlated error (W5.9)
# ---------------------------------------------------------------------------

#: Correlation length of the injected centroiding systematic, days. Long
#: compared with the epoch spacing, so the systematic is *smooth* -- which is
#: what makes it a systematic rather than extra white noise, and what an
#: independent GP on each axis can partly absorb.
JOINT_LENGTH_SCALE = 150.0

#: The injected coupling ``B``, in the physical parameterisation
#: :class:`~ampere.core.RotationCoupling` declares: a centroiding error
#: elongated along a position angle of 0.7 rad, with semi-axis standard
#: deviations of 0.030 mas and 0.008 mas.
#:
#: The correlation it implies between the ``ra`` and ``dec`` residuals at one
#: epoch is about 0.8, which is the number the whole study turns on: two
#: independent GPs can reproduce each axis's *marginal* scatter exactly and
#: can say nothing at all about that 0.8, so they treat two strongly dependent
#: measurements as two independent ones and report intervals that are too
#: narrow.
JOINT_TRUTH: dict[str, float] = {
    "angle": 0.7,
    "log_variance_0": float(np.log(0.030**2)),
    "log_variance_1": float(np.log(0.008**2)),
}


#: How far each channel's per-epoch sigma may wander from :data:`SIGMA` in the
#: heteroscedastic arm (**W5.24**): a log-uniform factor of up to two either
#: way. Real astrometric solutions do this --- seeing, airmass and the number of
#: reference stars change epoch by epoch, and the two sky axes of one centroid
#: rarely share an error bar --- so the heteroscedastic arm is the realistic
#: one, and the equal-sigma arm W5.9 pinned is the special case.
SIGMA_SPREAD = float(np.log(2.0))


def channel_sigmas(seed: int = SEED) -> np.ndarray:
    """Each channel's own per-epoch sigma, ``(2, N)``: row 0 ``ra``, row 1 ``dec``.

    Drawn from its own stream (``seed + 2``), so switching the heteroscedastic
    arm on changes the error bars and the white noise scaled by them and
    nothing else: the injected systematic is the same draw either way.
    """
    rng = np.random.default_rng(seed + 2)
    return SIGMA * np.exp(rng.uniform(-SIGMA_SPREAD, SIGMA_SPREAD, size=(2, EPOCHS.size)))


def coupling_matrix(truth: dict[str, float] | None = None) -> np.ndarray:
    """``B`` at *truth* (:data:`JOINT_TRUTH` by default), as a 2x2 array."""
    return np.asarray(
        RotationCoupling(0.0, 0.0, 0.0).matrix(JOINT_TRUTH if truth is None else truth)
    )


def marginal_amplitudes(truth: dict[str, float] | None = None) -> dict[str, float]:
    """Each channel's own standard deviation under ``B (x) K_x``, in mas.

    ``sqrt(B_tt)``: the marginal scatter the injected systematic gives channel
    ``t``, with the cross-covariance ``B_01`` dropped. This is exactly what two
    *independent* GPs can reproduce and exactly where they stop --- which is why
    the comparison arm of the calibration study is given these numbers rather
    than a prior over them. Handing the independent model the right marginals
    leaves the missing cross-covariance as the only difference between the two
    arms, which is the difference the study is about.
    """
    matrix = coupling_matrix(truth)
    return {"ra": float(np.sqrt(matrix[0, 0])), "dec": float(np.sqrt(matrix[1, 1]))}


def synthetic_joint_data(
    model: Model,
    ra_instrument: Instrument,
    dec_instrument: Instrument,
    *,
    seed: int = SEED,
    heteroscedastic: bool = False,
) -> tuple[TimeSeries, TimeSeries]:
    """:func:`synthetic_data` with a **correlated** centroiding systematic injected.

    The misspecification W5.9's study is about. On top of the orbit and the
    per-epoch white noise, both channels are perturbed by one draw from
    ``B (x) K_x`` -- a smooth error shared by the two sky axes, of the kind a
    centroiding solution with a preferred direction produces. The draw is made
    here, in numpy, from the materialised Kronecker covariance rather than
    through the noise model, so that the data this study fits are generated
    independently of the code that scores them.

    ``heteroscedastic=True`` is **W5.24**'s arm: each channel carries its own
    sigma per epoch (:func:`channel_sigmas`), the white noise is drawn at it,
    and the containers say so --- which is what sends a
    :class:`~ampere.core.JointGaussianProcessNoise` fit of these data off
    W5.9's rotated path and onto the dense or reduced-rank route. Left
    ``False``, every draw is W5.9's.
    """
    # Two streams, and the split is deliberate: the **white** noise comes from
    # the same generator, in the same order, as :func:`synthetic_data`'s, so
    # the difference between the two data sets is *exactly* the injected
    # systematic and nothing else. That makes the injection inspectable --- a
    # reader (and ``tests/examples``) can subtract one from the other and look
    # at what was added --- rather than something only the code knows.
    systematic_rng = np.random.default_rng(seed + 1)
    rng = np.random.default_rng(seed)
    requirements = negotiate([ra_instrument, dec_instrument])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)
    ra_truth = np.asarray(ra_instrument(truth).values).ravel()
    dec_truth = np.asarray(dec_instrument(truth).values).ravel()

    kernel = Matern32(1.0, JOINT_LENGTH_SCALE, axes=("time",))
    covariance = kernel.matrix(EPOCHS[:, None], EPOCHS[:, None], kernel.resolve({}))
    stacked = np.kron(coupling_matrix(), covariance)
    # A relative floor on the diagonal: a Matern-3/2 covariance over twenty-odd
    # epochs at a 150-day length scale is near-singular in float64, and this is
    # a data generator rather than a likelihood, so the stabiliser is the
    # honest thing rather than a hidden one.
    stacked = stacked + np.eye(stacked.shape[0]) * 1e-10 * float(np.max(np.diag(stacked)))
    systematic = systematic_rng.multivariate_normal(np.zeros(stacked.shape[0]), stacked).reshape(
        2, -1
    )

    if heteroscedastic:
        sigmas = channel_sigmas(seed)
        ra_noisy = ra_truth + systematic[0] + rng.normal(0.0, sigmas[0], ra_truth.shape)
        dec_noisy = dec_truth + systematic[1] + rng.normal(0.0, sigmas[1], dec_truth.shape)
    else:
        sigmas = np.full((2, EPOCHS.size), SIGMA)
        ra_noisy = ra_truth + systematic[0] + rng.normal(0.0, SIGMA, ra_truth.shape)
        dec_noisy = dec_truth + systematic[1] + rng.normal(0.0, SIGMA, dec_truth.shape)
    observed_ra = TimeSeries(
        EPOCHS * u.day,
        ra_noisy * u.mas,
        uncertainty=sigmas[0] * u.mas,
    )
    observed_dec = TimeSeries(
        EPOCHS * u.day,
        dec_noisy * u.mas,
        uncertainty=sigmas[1] * u.mas,
    )
    return observed_ra, observed_dec
