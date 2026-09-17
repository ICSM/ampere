"""One :class:`~ampere.backends.reference.ReflexOrbit`, two ``TimeSeries`` datasets:
right ascension and declination, sharing one orbital-parameter evaluation —
the second modality W4.4's template was written to prove (W4.9). See
:doc:`the tutorial page </astrometry>` for the walk-through; this module is
the code it walks through.

The physics
-----------
A source drifting under proper motion and wobbling under a companion's reflex
motion (a circular, node-aligned orbit projected onto the sky — see
:mod:`ampere.backends.reference.astrometry`'s module docstring for why this
is the design sketch's own choice rather than an omission) emits two
channels, ``"ra"`` and ``"dec"``, from one model evaluation: ``period`` and
``phase`` are shared at the language level, needing no ``Tie``.

The two instruments
--------------------
**astrom_ra**, **astrom_dec** — each one ``EpochSample`` step, pinning the
model onto the observed epochs; no parameters, no arithmetic. Both bind
different channels (``"ra"``, ``"dec"``), so — unlike
:mod:`examples.sed_composition` and :mod:`examples.interferometry`, whose two
instruments share one channel — there is no label collision to walk through
here; distinct labels are given anyway, for symmetry with the other pages.

The likelihood
--------------
Three arms. ``--gp`` and the default are described below; ``--joint``
(**W5.9**) is one correlated process over *both* channels,
:class:`~ampere.core.JointGaussianProcessNoise` with ``K = B (x) K_x``,
fitted to data carrying an injected centroiding systematic shared by the two
sky axes. ``--sbc`` runs the calibration study that compares the arms: the
joint fit's proper-motion intervals cover at their nominal rate under that
systematic and two independent GPs' do not, because two independent GPs can
reproduce each axis's marginal scatter exactly and can say nothing at all
about the correlation between them.

The other two arms: independent Gaussian noise (the rigid
comparison), or :class:`~ampere.core.GaussianProcessNoise` with
:class:`~ampere.core.Matern32` on the ``QuasisepGP`` solver — the O(N) path a
:class:`~ampere.core.TimeSeries`'s one ordered ``time`` axis is exactly the
right shape for. Independent noise is the default because the recovery this
item is accepted on is about the orbital parameters, not about the flexible
likelihood's own calibration claim (that is M2's question, asked of a
misspecified sky model, not of this modality's simple truth).

Backends
--------
:func:`build_problem` takes a backend name and builds every piece from that
backend's own namespace (the composition's one-backend rule). On
``"reference"`` the model has no gradient, so :func:`fit` samples with
:class:`~ampere.inference.EmceeEngine`; on ``"torch"`` or ``"jax"`` it samples
with :class:`~ampere.inference.NUTSEngine`.

Run it::

    python -m examples.astrometry                      # reference, emcee
    python -m examples.astrometry --backend torch       # NUTS
    python -m examples.astrometry --gp                  # the flexible likelihood
    python -m examples.astrometry --joint               # the joint channel noise (W5.9)
    python -m examples.astrometry --sbc joint           # the calibration study

:mod:`tests.examples.test_astrometry_example` is this module's own coverage:
a fast, always-on suite at a tiny budget. The full-budget recovery this item
is accepted on is verified by running this module directly and by
``tests/astrometry``'s own pinned rows at a per-PR budget larger than the
smoke test's but far short of the full one.
"""

from __future__ import annotations

import argparse
import sys
import time
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    RotationCoupling,
)

from . import generators

__all__ = [
    "BACKENDS",
    "DEFAULT_BURN_IN",
    "DEFAULT_CHAINS",
    "DEFAULT_DRAWS",
    "DEFAULT_SBC_BURN_IN",
    "DEFAULT_SBC_COUNT",
    "DEFAULT_SBC_DRAWS",
    "DEFAULT_SBC_STEPS",
    "DEFAULT_SBC_WALKERS",
    "DEFAULT_STEPS",
    "DEFAULT_WALKERS",
    "DEFAULT_WARMUP",
    "JOINT_LOG_VARIANCE_PRIOR",
    "QUALIFIED_TRUTH",
    "SBC_PARAMETERS",
    "SBC_PHASE_WIDTH",
    "SBC_SHARED_PARAMETER",
    "backend_module",
    "build_instruments",
    "build_model",
    "build_problem",
    "calibrate",
    "coverage_at",
    "fit",
    "joint_noise",
    "main",
    "noise_module",
    "recovers_truth",
    "report",
    "sbc_problem",
]

BACKENDS = ("reference", "torch", "jax")

#: The truth, qualified by the names a built ``FittingProblem`` actually
#: samples.
QUALIFIED_TRUTH: dict[str, float] = {
    "model.pmra": generators.TRUTH["pmra"],
    "model.pmdec": generators.TRUTH["pmdec"],
    "model.period": generators.TRUTH["period"],
    "model.phase": generators.TRUTH["phase"],
    "model.amp_ra": generators.TRUTH["amp_ra"],
    "model.amp_dec": generators.TRUTH["amp_dec"],
}

# Budgets. Twelve epochs and six free parameters is a small problem, so a
# reference/emcee fit is fast even at a generous walker/step count; torch/jax
# NUTS needs far fewer gradient-based draws for the same coverage.
DEFAULT_WALKERS = 24
DEFAULT_STEPS = 1500
DEFAULT_BURN_IN = 500
DEFAULT_DRAWS = 500
DEFAULT_WARMUP = 500
DEFAULT_CHAINS = 2

#: The flexible likelihood's own hyperparameters, held fixed here (the point
#: of this study is recovering the orbit, not fitting the noise process).
GP_AMPLITUDE = 0.03
GP_LENGTH_SCALE = 120.0

# -- the joint arm (W5.9) ---------------------------------------------------

#: Prior on each of ``B``'s log-eigenvalues. It contains both of
#: :data:`~examples.astrometry.generators.JOINT_TRUTH`'s (-7.0 and -9.7, at one
#: and 1.7 standard deviations) and it is *informative*, deliberately: a
#: centroiding systematic is known a priori to be of order tens of microarcsec,
#: and a coupling free to reach a milliarcsecond would simply absorb the reflex
#: wobble instead --- the flexible likelihood's standing hazard, which M2 studies
#: at length and which a prior on the noise *scale* is the ordinary answer to.
#: On the **log** scale for the reason ``RotationCoupling`` asks for
#: log-variances at all: a variance spans orders of magnitude and NUTS wants the
#: parameterisation that makes it a location.
JOINT_LOG_VARIANCE_PRIOR = st.norm(-8.0, 1.0)

#: The comparison arm's GP amplitudes are **not** fitted: they are held at
#: :func:`~examples.astrometry.generators.marginal_amplitudes`, the standard
#: deviation the injected systematic actually gives each channel. That is the
#: whole design of the comparison. Two independent GPs can reproduce each
#: axis's marginal scatter exactly, and giving them a prior to discover it
#: instead would let them *over*-inflate and hide the effect under a fitted
#: nuisance --- which is what an earlier version of this study measured, and
#: why it measured nothing. Handing them the right marginals leaves the missing
#: cross-covariance as the only difference between the arms.

#: What the SBC study ranks. The orbital parameters only: the claim under test
#: is that the *physical* posterior is calibrated, and the noise model's own
#: parameters mean different things in the two arms (there is no coupling in an
#: independent-GP fit at all), so ranking them would compare two different
#: questions.
#:
#: ``model.phase`` is the one that carries the result, and the reason is worth
#: stating because it is the whole mechanism. ``pmra`` enters only the ``ra``
#: channel and ``pmdec`` only ``dec``, so **each is constrained by one channel**
#: --- and ignoring a correlation between two channels that constrain different
#: parameters throws information away, which makes an interval too *wide*, not
#: too narrow. The orbital phase is shared by both channels, so both constrain
#: it, and there ignoring the correlation is **double-counting**: two error-laden
#: estimates that move together are combined as though they moved independently,
#: and the combined interval is narrower than the truth's own scatter. That is
#: undercoverage, and it is what a cross-channel systematic does to every
#: parameter a multi-channel instrument measures jointly.
SBC_PARAMETERS: tuple[str, ...] = ("model.pmra", "model.pmdec", "model.phase")

#: The parameter the coverage comparison is pinned on --- see above.
SBC_SHARED_PARAMETER = "model.phase"

#: Prior width on the shared orbital phase in the calibration study, radians.
#: Informative for :func:`build_model`'s own reason for the period: a phase
#: search over the whole circle at this sampling is a different (and harder)
#: problem, and the realistic case is a phase already known to a tenth of a
#: radian from a previous epoch. Ten times wider than the likelihood's own
#: width, so the posterior is the data's and not the prior's.
SBC_PHASE_WIDTH = 0.15

#: SBC budgets. Deliberately small enough to run inside a test: 32 refits of a
#: two-parameter problem over 28 epochs. `python -m examples.astrometry --sbc`
#: with `--sbc-count` raised is the real study.
#: How much of ``B`` the calibration study fits: **none of it**. Both arms are
#: given the noise process they are entitled to know --- the joint arm the whole
#: of ``B``, the comparison arm each channel's correct marginal --- so that the
#: one thing that differs between them is the cross-covariance, which is the one
#: thing the study is about. Fitting ``B`` as well is what ``--joint`` does and
#: what the NUTS rows in ``tests/inference`` exercise; doing it *here* would mix
#: a question about the model with a question about a sampler's ability to
#: explore a bimodal noise sector (the coupling's own relabelling symmetry, see
#: ``RotationCoupling``) at a budget a test can afford.
COUPLING_FIT = "none"

DEFAULT_SBC_COUNT = 32
DEFAULT_SBC_DRAWS = 200
DEFAULT_SBC_WALKERS = 16
DEFAULT_SBC_STEPS = 900
DEFAULT_SBC_BURN_IN = 400


def backend_module(backend: str) -> Any:
    """The namespace *backend*'s model and instrument pieces live in."""
    if backend == "reference":
        from ampere.backends import reference as module

        return module
    if backend == "torch":
        from ampere.backends import torch as module  # type: ignore[assignment]

        return module
    if backend == "jax":
        from ampere.backends import jax as module  # type: ignore[assignment]

        module.configure_x64()
        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {BACKENDS}.")


def noise_module(backend: str) -> Any:
    """Where ``IndependentNoise`` comes from for *backend*.

    ``ampere.core.IndependentNoise`` already declares ``BACKEND =
    "reference"``, so on the reference backend it is imported from
    ``ampere.core`` directly rather than through the backend package, exactly
    as :mod:`examples.sed_composition.sed_composition` does.
    """
    if backend == "reference":
        import ampere.core as module

        return module
    return backend_module(backend)


def build_model(backend: str) -> Any:
    """The one :class:`~ampere.backends.reference.ReflexOrbit`, on *backend*."""
    module = backend_module(backend)
    return module.ReflexOrbit(
        generators.EPOCHS,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        # A period search over many decades is a genuinely multi-modal
        # problem at this modality's sparse, irregular sampling -- unlike
        # every model the interferometry template fits, a reflex orbit is
        # periodic, and a period prior wide enough to reach past the epochs'
        # own baseline lets a sampler alias onto a spurious cycle indefinitely
        # (see docs/source/astrometry.rst's closing section). A period search
        # informed to within a few tens of days -- the realistic case once a
        # periodogram or a previous epoch has suggested roughly where to look
        # -- is what this study fits, rather than a global period search,
        # which is a different (and harder) problem this item does not claim
        # to solve.
        period=st.norm(400.0, 30.0),
        phase=st.uniform(0.0, 2.0 * np.pi),
        amp_ra=st.uniform(0.0, 2.0),
        amp_dec=st.uniform(0.0, 2.0),
    )


def build_instruments(backend: str) -> tuple[Instrument, Instrument]:
    """The two epoch-sampling instruments, ``"astrom_ra"`` and ``"astrom_dec"``."""
    module = backend_module(backend)
    ra_instrument = Instrument(
        [module.EpochSample(generators.EPOCHS)], channel="ra", label="astrom_ra"
    )
    dec_instrument = Instrument(
        [module.EpochSample(generators.EPOCHS)], channel="dec", label="astrom_dec"
    )
    return ra_instrument, dec_instrument


def joint_noise(backend: str, *, fit: str = "full") -> Any:
    """The joint noise model over the two channels (**W5.9**).

    ``K = B (x) K_x``: one correlated process over ``ra`` and ``dec``, with
    ``B`` in :class:`~ampere.core.RotationCoupling`'s physical parameterisation
    (a position angle and two log-variances) and ``K_x`` a Matern-3/2 along the
    shared epoch grid.

    ``fit`` chooses how much of ``B`` is fitted: ``"full"`` frees the angle and
    both log-variances, ``"variances"`` holds the angle at the injected value
    (which is what the calibration study does --- see :func:`sbc_problem` for
    why), and ``"none"`` holds all three.

    Two things are held fixed and both are deliberate. The **kernel's
    amplitude** must be: ``B (x) (a^2 K) = (a^2 B) (x) K``, so a free amplitude
    beside a free ``B`` is one degree of freedom written twice, and
    ``JointGaussianProcessNoise`` refuses the pair by name. The **length
    scale** is fixed at the injected one because this study is about the
    cross-channel structure rather than about the correlation length --- the
    flexible likelihood's own calibration claim is M2's question, asked of a
    misspecified physical model, not of this modality's known systematic.
    """
    noise = noise_module(backend)
    kernel = noise.Matern32(1.0, generators.JOINT_LENGTH_SCALE, axes=("time",))
    angle: Any = st.uniform(0.0, np.pi) if fit == "full" else generators.JOINT_TRUTH["angle"]
    variances: list[Any] = (
        [JOINT_LOG_VARIANCE_PRIOR, JOINT_LOG_VARIANCE_PRIOR]
        if fit in ("full", "variances")
        else [generators.JOINT_TRUTH["log_variance_0"], generators.JOINT_TRUTH["log_variance_1"]]
    )
    coupling = RotationCoupling(angle, *variances)
    return noise.JointGaussianProcessNoise(
        kernel, noise.QuasisepGP(), datasets=("ra", "dec"), coupling=coupling
    )


def build_problem(
    backend: str = "reference",
    *,
    gp: bool = False,
    joint: bool = False,
    injected: bool | None = None,
    seed: int = generators.SEED,
) -> FittingProblem:
    """The composed problem: one model, two channels, distinct labels.

    ``gp=True`` swaps independent noise for the flexible likelihood
    (``GaussianProcessNoise(Matern32(...), QuasisepGP())``) on both channels.

    ``joint=True`` is **W5.9**'s arm: one correlated process over both
    channels (:func:`joint_noise`), declared on the ``DatasetCollection``
    rather than on either dataset's likelihood, because it is one covariance
    over two datasets. It also switches the data to the ones carrying an
    injected correlated centroiding systematic
    (:func:`~examples.astrometry.generators.synthetic_joint_data`) --- there is
    no point fitting a cross-channel model to data with no cross-channel
    structure. ``injected=`` overrides that pairing, which is how the
    comparison arm is built: ``build_problem(gp=True, injected=True)`` fits the
    *same* systematic-bearing data with two independent GPs, and is the fit
    whose coverage the study shows is not nominal.
    """
    model = build_model(backend)
    ra_instrument, dec_instrument = build_instruments(backend)
    generate = (
        generators.synthetic_joint_data
        if (joint if injected is None else injected)
        else generators.synthetic_data
    )
    observed_ra, observed_dec = generate(model, ra_instrument, dec_instrument, seed=seed)
    noise = noise_module(backend)

    def likelihood() -> Likelihood:
        if joint or not gp:
            # A channel of a joint group carries the family and nothing else:
            # the group owns the whole covariance, diagonal included, and a
            # per-channel scale or jitter would give the channels different
            # diagonals, which is the one thing the rotation may not have.
            return Likelihood(GaussianFamily(), noise.IndependentNoise())
        kernel = noise.Matern32(GP_AMPLITUDE, GP_LENGTH_SCALE, axes=("time",))
        return Likelihood(GaussianFamily(), noise.GaussianProcessNoise(kernel, noise.QuasisepGP()))

    datasets = DatasetCollection(
        {
            "ra": Dataset(observed_ra, ra_instrument, likelihood=likelihood(), label="ra"),
            "dec": Dataset(observed_dec, dec_instrument, likelihood=likelihood(), label="dec"),
        },
        joint={"astrom": joint_noise(backend)} if joint else None,
    )
    return FittingProblem(model, datasets, seed=seed)


def fit(
    problem: FittingProblem,
    *,
    backend: str,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    draws: int | None = None,
    warmup: int | None = None,
    chains: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the engine its backend can feed."""
    if backend == "reference":
        from ampere.inference import EmceeEngine

        engine = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return engine.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    from ampere.inference import NUTSEngine

    engine = NUTSEngine(problem)
    return engine.run(
        DEFAULT_DRAWS if draws is None else draws,
        warmup=DEFAULT_WARMUP if warmup is None else warmup,
        chains=DEFAULT_CHAINS if chains is None else chains,
        progress=progress,
    )


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each parameter's central *level* interval contains its truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any) -> str:
    """A human-readable posterior summary, truth in brackets, 95 % coverage flagged."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    covered = recovers_truth(run)
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        "  posterior (truth in brackets, 95 % coverage flagged):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [2.5, 97.5])
        truth = QUALIFIED_TRUTH.get(name)
        flag = "" if truth is None else ("  ok" if covered[name] else "  MISS")
        bracket = "" if truth is None else f"  (truth {truth:+.6g})"
        lines.append(
            f"    {name:20s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# The calibration study: does the joint model buy coverage back? (W5.9)
# ---------------------------------------------------------------------------


def sbc_problem(
    backend: str = "reference",
    *,
    arm: str = "joint",
    observed: tuple[Any, Any] | None = None,
    seed: int = generators.SEED,
) -> FittingProblem:
    """The **two-parameter** problem the calibration study simulates and fits.

    Not :func:`build_problem`. Simulation-based calibration refits once per
    simulation, so the problem has to be small; and it draws the truth from the
    prior, so the problem has to be *unimodal*, which a reflex orbit with a free
    period is emphatically not at this modality's sparse sampling (see the
    closing section of :doc:`the tutorial page </astrometry>`). So the orbit's
    period, phase and semi-amplitudes are held at the truth and the two proper
    motions are fitted --- which is the pair the coverage claim is about
    anyway, because they are the physical parameters a shared centroiding
    systematic actually biases.

    Three arms:

    Three parameters are fitted: the two proper motions, one per channel, and
    the orbital **phase**, which both channels measure --- see
    :data:`SBC_PARAMETERS` for why the shared one is where the result lives.

    ``"joint"``
        ``JointGaussianProcessNoise`` over both channels: the model that
        matches the generating process.
    ``"independent"``
        Two ordinary ``GaussianProcessNoise`` likelihoods, each with the
        **correct marginal amplitude** for its own axis (see
        :data:`JOINT_LOG_VARIANCE_PRIOR`'s neighbour above). Each reproduces
        its own axis's scatter exactly and neither can say anything about the
        correlation between them, so this arm treats two strongly dependent
        measurements as independent --- and that, and nothing else, is what
        separates it from the joint arm. ``likelihoods.md`` §15's "two scalar
        GPs with tied hyperparameters is the nearest approximation and is a
        different model", made into an experiment.
    ``"rigid"``
        Independent white noise, the arm with no flexible likelihood at all.
        Kept because it is the comparison a reader reaches for first.

    *observed* replaces the two containers, which is how the comparison arm
    refits the simulating arm's own data.
    """
    module = backend_module(backend)
    noise = noise_module(backend)
    model = module.ReflexOrbit(
        generators.EPOCHS,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        period=generators.TRUTH["period"],
        phase=st.norm(generators.TRUTH["phase"], SBC_PHASE_WIDTH),
        amp_ra=generators.TRUTH["amp_ra"],
        amp_dec=generators.TRUTH["amp_dec"],
    )
    ra_instrument, dec_instrument = build_instruments(backend)
    if observed is None:
        observed = generators.synthetic_joint_data(
            module.ReflexOrbit(generators.EPOCHS, **generators.TRUTH),
            ra_instrument,
            dec_instrument,
            seed=seed,
        )
    observed_ra, observed_dec = observed

    amplitudes = generators.marginal_amplitudes()

    def likelihood(channel: str) -> Likelihood:
        if arm == "independent":
            kernel = noise.Matern32(
                amplitudes[channel], generators.JOINT_LENGTH_SCALE, axes=("time",)
            )
            return Likelihood(
                GaussianFamily(), noise.GaussianProcessNoise(kernel, noise.QuasisepGP())
            )
        return Likelihood(GaussianFamily(), noise.IndependentNoise())

    datasets = DatasetCollection(
        {
            "ra": Dataset(observed_ra, ra_instrument, likelihood=likelihood("ra"), label="ra"),
            "dec": Dataset(observed_dec, dec_instrument, likelihood=likelihood("dec"), label="dec"),
        },
        joint={"astrom": joint_noise(backend, fit=COUPLING_FIT)} if arm == "joint" else None,
    )
    return FittingProblem(model, datasets, seed=seed)


def calibrate(
    backend: str = "reference",
    *,
    arm: str = "joint",
    count: int = DEFAULT_SBC_COUNT,
    draws: int = DEFAULT_SBC_DRAWS,
    walkers: int = DEFAULT_SBC_WALKERS,
    steps: int = DEFAULT_SBC_STEPS,
    burn_in: int = DEFAULT_SBC_BURN_IN,
    seed: int = generators.SEED,
) -> Any:
    """Simulation-based calibration of one arm, against the **joint** generator.

    The simulating problem is always the joint one, so every replicate's data
    carry a correlated centroiding systematic drawn from ``B (x) K_x`` at a
    ``B`` drawn from its own prior --- the correlated draw
    ``DatasetCollection.draw_group`` makes, not two marginal ones. What changes
    between arms is only what is *fitted*: the same data, scored by a model
    that knows about the cross-channel structure or by one that does not.

    Returns the ``calibration`` group :func:`ampere.results.sbc` produces:
    ranks, the coverage curve and the uniformity p-value, for
    :data:`SBC_PARAMETERS`.
    """
    from ampere.inference import EmceeEngine
    from ampere.results import sbc

    simulating = sbc_problem(backend, arm="joint", seed=seed)

    def factory(replica: FittingProblem) -> Any:
        if arm == "joint":
            fitted = replica
        else:
            fitted = sbc_problem(
                backend,
                arm=arm,
                observed=(replica.datasets["ra"].observed, replica.datasets["dec"].observed),
                seed=seed,
            )
        return EmceeEngine(fitted, walkers=walkers)

    return sbc(
        simulating,
        factory,
        count=count,
        draws=draws,
        run_options={"steps": steps, "burn_in": burn_in},
        parameters=list(SBC_PARAMETERS),
        seed=seed,
        label=f"astrometry {arm} arm",
    )


def coverage_at(calibration: Any, level: float = 0.9, parameter: str | None = None) -> float:
    """Empirical coverage at a nominal *level*, for one parameter or averaged.

    One number, because the claim is one claim: "the central *level* interval
    contains the truth *level* of the time". *parameter* names one of
    :data:`SBC_PARAMETERS` --- pass :data:`SBC_SHARED_PARAMETER` for the one the
    comparison turns on --- and ``None`` averages over all of them. The curve is
    on the returned dataset for anyone who wants the rest of it.
    """
    coverage = calibration["coverage"]
    nearest = int(np.argmin(np.abs(np.asarray(coverage["level"], dtype=float) - float(level))))
    row = coverage.isel(level=nearest)
    if parameter is not None:
        names = [str(name) for name in np.asarray(coverage["parameter"])]
        return float(np.asarray(row, dtype=float)[names.index(parameter)])
    return float(np.mean(np.asarray(row, dtype=float)))


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", default="reference", choices=list(BACKENDS))
    parser.add_argument("--gp", action="store_true", help="the flexible likelihood")
    parser.add_argument(
        "--joint",
        action="store_true",
        help="one correlated process over both channels, on injected correlated data (W5.9)",
    )
    parser.add_argument(
        "--sbc",
        choices=("joint", "independent", "rigid"),
        default=None,
        help="run the calibration study for one arm instead of a single fit",
    )
    parser.add_argument("--sbc-count", type=int, default=DEFAULT_SBC_COUNT)
    parser.add_argument("--sbc-draws", type=int, default=DEFAULT_SBC_DRAWS)
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None, help="emcee only")
    parser.add_argument("--steps", type=int, default=None, help="emcee only")
    parser.add_argument("--burn-in", type=int, default=None, help="emcee only")
    parser.add_argument("--draws", type=int, default=None, help="NUTS only")
    parser.add_argument("--warmup", type=int, default=None, help="NUTS only")
    parser.add_argument("--chains", type=int, default=None, help="NUTS only")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    if args.sbc is not None:
        started = time.perf_counter()
        calibration = calibrate(
            args.backend, arm=args.sbc, count=args.sbc_count, draws=args.sbc_draws, seed=args.seed
        )
        elapsed = time.perf_counter() - started
        print(f"SBC, {args.sbc} arm: {args.sbc_count} simulation(s), {args.sbc_draws} draw(s)")
        for level in (0.5, 0.9, 0.95):
            shared = coverage_at(calibration, level, SBC_SHARED_PARAMETER)
            print(
                f"  coverage at {level:.2f}: {coverage_at(calibration, level):.3f} "
                f"(mean), {shared:.3f} ({SBC_SHARED_PARAMETER})"
            )
        pvalues = np.asarray(calibration["ks_pvalue"], dtype=float)
        for name, pvalue in zip(SBC_PARAMETERS, pvalues, strict=True):
            print(f"  rank uniformity p({name}) = {pvalue:.3f}")
        print(f"  {elapsed:.1f} s wall clock")
        return 0

    problem = build_problem(args.backend, gp=args.gp, joint=args.joint, seed=args.seed)
    print(f"negotiated channels: {list(problem.requirements['model'])}")
    for channel_name in problem.requirements["model"]:
        req = problem.requirements["model"][channel_name]
        print(f"sources asking of channel {channel_name!r}: {req.sources}")
    print(f"free parameters: {problem.parameters.free_names}")

    started = time.perf_counter()
    run = fit(
        problem,
        backend=args.backend,
        walkers=args.walkers,
        steps=args.steps,
        burn_in=args.burn_in,
        draws=args.draws,
        warmup=args.warmup,
        chains=args.chains,
    )
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
