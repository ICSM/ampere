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

``--joint --heteroscedastic`` (**W5.24**) is the same arm on data whose two
channels carry *different* per-epoch error bars, as real astrometric
solutions do. The rotation that makes W5.9's arm ``T·O(N)`` is exact only
when the channels share one sigma vector, so this arm binds ``DenseGP`` and
the group takes the dense route --- ``B (x) K_x + blockdiag(diag sigma_t^2)``
factorised directly. At this study's 28 epochs that is the cheapest route as
well as the exact one (measured per likelihood call on one core: dense
0.38 ms, W5.9's rotated ``QuasisepGP`` path 0.46 ms, the reduced-rank
``HilbertSpaceGP`` route 0.64 ms at ``m = 32`` with its error still 6e-3);
the reduced-rank route earns its keep at larger ``N`` and under NUTS, where
its ``T·m`` whitened block is what ``tests/inference`` exercises.

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
    python -m examples.astrometry --joint --heteroscedastic   # unequal channel sigmas (W5.24)
    python -m examples.astrometry --sbc joint           # the calibration study
    python -m examples.astrometry --wide-prior          # the aliasing hazard, its remedy (W5.15)
    python -m examples.astrometry --wide-prior --engine nautilus  # or ultranest, in `-e nested`

The period-aliasing arm (W5.15)
--------------------------------
``--wide-prior`` swaps the informed period prior (``norm(400, 30)``) for
``loguniform(50, 2000)`` — W4.9's exact measured prior, on the same
twenty-eight epochs — which reproduces, rather than merely describes, the
hazard :doc:`the tutorial page </astrometry>` §6 used to only narrate: an
emcee ensemble or a multi-chain NUTS run both alias onto one spurious period
and stay there. The remedy this arm exists to demonstrate is nested
sampling, which does not need to choose a mode at all: ``--engine dynesty``
(the default once ``--wide-prior`` is given), ``--engine nautilus`` or
``--engine ultranest`` return every alias as a separately-weighed mode,
which :func:`period_modes` extracts from the equal-weight draws and
``main`` prints as a table, mass fraction and local evidence
(``ln Z_k = ln Z + ln f_k``) included. §6 and §9 of the tutorial page carry
the measured table; this module's docstring and code are the thing it
measures.

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
    "ENGINES",
    "HETEROSCEDASTIC_SOLVER",
    "JOINT_LOG_VARIANCE_PRIOR",
    "QUALIFIED_TRUTH",
    "SBC_DIRECTION",
    "SBC_PARAMETERS",
    "WIDE_PERIOD_PRIOR",
    "backend_module",
    "build_instruments",
    "build_model",
    "build_problem",
    "calibrate",
    "coverage_at",
    "direction_of",
    "fit",
    "joint_noise",
    "main",
    "noise_module",
    "period_modes",
    "recovers_truth",
    "report",
    "sbc_problem",
]

BACKENDS = ("reference", "torch", "jax")

#: The engines :func:`fit` can be asked for by name. ``None`` (the CLI's own
#: default) keeps the pre-W5.15 behaviour: emcee on the reference backend,
#: NUTS on torch/jax. The three nested samplers are the standard remedy for
#: the multi-modal posterior :data:`WIDE_PERIOD_PRIOR` produces (W5.15).
ENGINES = ("emcee", "nuts", "dynesty", "nautilus", "ultranest")

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

#: W4.9's exact measured period prior (**W5.15**): wide enough to reach past
#: the epochs' own baseline and genuinely multi-modal as a result. Reproduced
#: verbatim by ``build_model(wide_prior=True)`` rather than re-measured, so
#: the aliasing this module now demonstrates is the same finding
#: ``docs/source/astrometry.rst`` §6 always described, not a new one.
WIDE_PERIOD_PRIOR = st.loguniform(50.0, 2000.0)

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

#: The route the heteroscedastic arm takes (**W5.24**): ``"dense"``, chosen by
#: measured cost at this study's 28 epochs (see the module docstring). The
#: other value :func:`joint_noise` accepts is ``"hilbert"``, the reduced-rank
#: route; ``"quasisep"`` is W5.9's rotated path and refuses unequal channels.
HETEROSCEDASTIC_SOLVER = "dense"

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
SBC_PARAMETERS: tuple[str, ...] = ("model.pmra", "model.pmdec")

#: The derived quantity the comparison is pinned on, and the reason it is a
#: *derived* one rather than a parameter.
#:
#: A cross-channel systematic does not bias a single channel's own parameter ---
#: it correlates the **errors** of the two channels' parameters. ``pmra`` enters
#: only ``ra`` and ``pmdec`` only ``dec``, and both are measured with the same
#: weight along the epoch grid (a proper motion is a linear trend), so a
#: centroiding error shared by the two axes makes their two errors move
#: together. Each *marginal* posterior is then still about the right width ---
#: which is why both arms' marginal coverage comes out nominal, and this study
#: reports that --- while the **joint** posterior has the wrong shape: the
#: independent-GP fit reports two errors as uncorrelated when they are
#: correlated at about 0.86.
#:
#: The projection that sees it is the one along which the two errors add:
#: ``(pmra + pmdec)/sqrt(2)``, the **diagonal** of the proper-motion plane. Its
#: true variance is ``v(1 + rho)``; a model that believes the errors
#: independent reports ``v``, understating the interval by ``sqrt(1 + rho)`` ---
#: about 1.36 here. That is undercoverage, and it is what a shared centroiding
#: systematic does to the *direction* of a measured proper motion, which is a
#: quantity astronomers publish.
SBC_DIRECTION = "model.direction"

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


def build_model(backend: str, *, wide_prior: bool = False) -> Any:
    """The one :class:`~ampere.backends.reference.ReflexOrbit`, on *backend*.

    ``wide_prior=True`` swaps the informed period prior for
    :data:`WIDE_PERIOD_PRIOR` (**W5.15**), everything else unchanged.
    """
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
        # (see docs/source/astrometry.rst §6). A period informed to within a
        # few tens of days -- the realistic case once a periodogram or a
        # previous epoch has suggested roughly where to look -- is what this
        # study fits by default, rather than a global period search over
        # decades. **W5.15** runs that global search anyway, deliberately,
        # under ``wide_prior=True``: not because it is now easy, but because
        # nested sampling turns "genuinely multi-modal" from a hazard an
        # ensemble or a gradient sampler cannot escape into a posterior whose
        # modes are separately weighed and reported, evidence included (see
        # the module docstring and :func:`period_modes`).
        period=WIDE_PERIOD_PRIOR if wide_prior else st.norm(400.0, 30.0),
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


def joint_noise(backend: str, *, fit: str = "full", solver: str = "quasisep") -> Any:
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

    ``solver`` (**W5.24**) picks the bound strategy, and with it the route
    unequal channel sigmas take: ``"quasisep"`` (W5.9's rotated O(N) path,
    equal sigmas only), ``"dense"`` (exact, any sigmas) or ``"hilbert"`` (the
    reduced-rank route, a ``HilbertSpaceGP`` whose box is three data
    half-extents wide --- a 150-day length scale against a 1100-day baseline
    needs the room).
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
    if solver == "quasisep":
        strategy = noise.QuasisepGP()
    elif solver == "dense":
        strategy = noise.DenseGP()
    elif solver == "hilbert":
        strategy = noise.HilbertSpaceGP(basis_size=32, boundary_factor=3.0)
    else:
        raise ValueError(f"unknown solver {solver!r}; the three are quasisep, dense, hilbert.")
    return noise.JointGaussianProcessNoise(
        kernel, strategy, datasets=("ra", "dec"), coupling=coupling
    )


def build_problem(
    backend: str = "reference",
    *,
    gp: bool = False,
    joint: bool = False,
    injected: bool | None = None,
    heteroscedastic: bool = False,
    wide_prior: bool = False,
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

    ``heteroscedastic=True`` (**W5.24**) draws the injected data with a sigma
    per channel per epoch and binds the joint arm's
    :data:`HETEROSCEDASTIC_SOLVER`, so the group takes the dense route.

    ``wide_prior=True`` (**W5.15**) is :func:`build_model`'s own flag,
    threaded through: the period prior becomes :data:`WIDE_PERIOD_PRIOR`,
    everything else --- data, instruments, noise, the other five priors ---
    unchanged. It is the arm :func:`fit` and the CLI's ``--wide-prior`` exist
    to demonstrate the second remedy on.
    """
    model = build_model(backend, wide_prior=wide_prior)
    ra_instrument, dec_instrument = build_instruments(backend)
    if joint if injected is None else injected:
        observed_ra, observed_dec = generators.synthetic_joint_data(
            model, ra_instrument, dec_instrument, seed=seed, heteroscedastic=heteroscedastic
        )
    else:
        observed_ra, observed_dec = generators.synthetic_data(
            model, ra_instrument, dec_instrument, seed=seed
        )
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
        joint=(
            {
                "astrom": joint_noise(
                    backend, solver=HETEROSCEDASTIC_SOLVER if heteroscedastic else "quasisep"
                )
            }
            if joint
            else None
        ),
    )
    return FittingProblem(model, datasets, seed=seed)


def fit(
    problem: FittingProblem,
    *,
    backend: str,
    engine: str | None = None,
    live_points: int | None = None,
    dlogz: float | None = None,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    draws: int | None = None,
    warmup: int | None = None,
    chains: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the engine its backend can feed, or the one asked for.

    ``engine=None`` (the default) is unchanged from before **W5.15**: emcee
    on the reference backend, NUTS on torch/jax. Set ``engine`` to
    ``"dynesty"``, ``"nautilus"`` or ``"ultranest"`` to run one of the three
    nested samplers instead --- the standard second remedy for a genuinely
    multi-modal posterior, such as :func:`build_problem`'s ``wide_prior=True``
    arm produces (``docs/source/astrometry.rst`` §6).

    ``live_points`` is every nested sampler's own live-set size (``None``
    defaults to :func:`~ampere.inference.engine.default_live_points`, shared
    by all three). ``dlogz`` is dynesty's and ultranest's stopping criterion
    on the remaining evidence, passed through when given; nautilus stops on
    its own ``f_live``/``n_eff`` instead and does not take a ``dlogz``, so
    this argument is silently unused when ``engine="nautilus"`` rather than
    raising --- the same "extra keywords the chosen engine does not use are
    ignored" contract :func:`fit` already has for ``walkers``/``draws`` and
    the rest.
    """
    if engine is None:
        engine = "emcee" if backend == "reference" else "nuts"
    if engine not in ENGINES:
        raise ValueError(f"unknown engine {engine!r}; choices are {ENGINES}.")

    if engine == "emcee":
        from ampere.inference import EmceeEngine

        sampler = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return sampler.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "nuts":
        from ampere.inference import NUTSEngine

        sampler = NUTSEngine(problem)
        return sampler.run(
            DEFAULT_DRAWS if draws is None else draws,
            warmup=DEFAULT_WARMUP if warmup is None else warmup,
            chains=DEFAULT_CHAINS if chains is None else chains,
            progress=progress,
        )
    if engine == "dynesty":
        from ampere.inference import DynestyEngine

        sampler = DynestyEngine(problem, live_points=live_points)
        run_options = {} if dlogz is None else {"dlogz": dlogz}
        return sampler.run(progress=progress, **run_options)
    if engine == "nautilus":
        from ampere.inference import NautilusEngine

        sampler = NautilusEngine(problem, live_points=live_points)
        return sampler.run(progress=progress)
    from ampere.inference import UltranestEngine

    sampler = UltranestEngine(problem, live_points=live_points)
    run_options = {} if dlogz is None else {"dlogz": dlogz}
    return sampler.run(progress=progress, **run_options)


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
    """A human-readable posterior summary, truth in brackets, 95 % coverage flagged.

    Appends the engine-neutral evidence triple (**W5.15**) when *run* carries
    one --- every nested-sampling run does, no ensemble or gradient run does
    (``results.md`` §9).
    """
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
    if "ampere_log_evidence" in attrs:
        evidence = float(attrs["ampere_log_evidence"])
        err = attrs.get("ampere_log_evidence_err")
        method = attrs.get("ampere_evidence_method", "?")
        err_bit = "" if err is None else f" +- {float(err):.3f}"
        lines.append(f"  ln Z = {evidence:+.3f}{err_bit}  ({method})")
    return "\n".join(lines)


def period_modes(
    run: Any,
    *,
    parameter: str = "model.period",
    gap: float = 0.1,
    minimum_mass: float = 0.02,
) -> list[dict[str, float]]:
    """Split *run*'s equal-weight *parameter* draws into aliasing modes (**W5.15**).

    A period search under :data:`WIDE_PERIOD_PRIOR` is genuinely multi-modal
    (see the module docstring); a nested-sampling run's equal-weight draws
    already carry every mode the sampler found, unlike an ensemble or a
    gradient chain's single-mode draws, so recovering them is a matter of
    *finding the modes in the draws already there*, not of re-sampling.

    The rule is deliberately simple, because the aliases at this study's
    baseline are not subtle: sort the draws in log space (a period search's
    natural scale --- aliases sit at roughly one cycle apart in ``1/P``, which
    is a large, roughly constant separation in ``ln P``, not in ``P`` itself)
    and split wherever a consecutive gap exceeds *gap*. A cluster below
    *minimum_mass* of the total draws --- sampling noise in the sampler's own
    boundary rather than a real mode --- is dropped.

    ``gap=0.1`` is measured, not guessed (**W5.15**): the true period and its
    nearest aliases, as an emcee ensemble at this study's default budget
    actually splits across (``docs/source/astrometry.rst`` §6), sit about
    ``ln(762.9 / 399.8) ~= 0.65`` apart in log space, while a single
    well-resolved mode's own 16/84 % spread is under ``0.01`` -- ``0.1`` sits
    comfortably in the two-order-of-magnitude gap between them.

    Each surviving mode is returned as a dict with:

    ``median``, ``lower``, ``upper``
        The draws' median and central 16/84 % interval, on *parameter*'s own
        scale (not logged back out of convenience; the split, not the
        report, is what needs the log scale).
    ``mass_fraction``
        ``f_k``, the fraction of *all* the run's draws in this mode.
    ``count``
        The number of draws in this mode, for a sanity check against
        ``mass_fraction`` and the run's total.
    ``log_evidence``
        ``ln Z_k = ln Z + ln f_k`` --- present only when *run* carries
        ``ampere_log_evidence`` (nested-sampling runs; an ensemble or
        gradient run has no evidence to apportion and this key is omitted
        rather than written as ``None``, so a caller's ``"log_evidence" in
        mode`` is the one check it needs).

    Modes are returned ordered by increasing *parameter*, not by mass ---
    call sorted on ``mass_fraction`` for "which mode is biggest".
    """
    posterior = run["posterior"].dataset
    draws = np.sort(np.asarray(posterior[parameter], dtype=float).ravel())
    total = draws.size
    log_draws = np.log(draws)
    split_at = np.flatnonzero(np.diff(log_draws) > gap) + 1
    clusters = np.split(draws, split_at)

    log_evidence = run.attrs.get("ampere_log_evidence")

    modes: list[dict[str, float]] = []
    for cluster in clusters:
        fraction = cluster.size / total
        if fraction < minimum_mass:
            continue
        lower, upper = np.percentile(cluster, [16.0, 84.0])
        mode: dict[str, float] = {
            "median": float(np.median(cluster)),
            "lower": float(lower),
            "upper": float(upper),
            "mass_fraction": float(fraction),
            "count": float(cluster.size),
        }
        if log_evidence is not None:
            mode["log_evidence"] = float(log_evidence) + float(np.log(fraction))
        modes.append(mode)
    return modes


# ---------------------------------------------------------------------------
# The calibration study: does the joint model buy coverage back? (W5.9)
# ---------------------------------------------------------------------------


def sbc_problem(
    backend: str = "reference",
    *,
    arm: str = "joint",
    observed: tuple[Any, Any] | None = None,
    heteroscedastic: bool = False,
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

    Two parameters are fitted, one per channel: the proper motions. Everything
    else about the orbit is held at the truth, so that a study which refits once
    per simulation stays affordable and stays unimodal --- a reflex orbit with a
    free period is neither (see the closing section of :doc:`the tutorial page
    </astrometry>`).

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

    ``heteroscedastic=True`` (**W5.24**) draws the data with a sigma per channel
    per epoch and binds :data:`HETEROSCEDASTIC_SOLVER` on the joint arm; the
    comparison arms need nothing new, because each of their likelihoods reads
    its own container's sigmas already.
    """
    module = backend_module(backend)
    noise = noise_module(backend)
    model = module.ReflexOrbit(
        generators.EPOCHS,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        period=generators.TRUTH["period"],
        phase=generators.TRUTH["phase"],
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
            heteroscedastic=heteroscedastic,
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
        joint=(
            {
                "astrom": joint_noise(
                    backend,
                    fit=COUPLING_FIT,
                    solver=HETEROSCEDASTIC_SOLVER if heteroscedastic else "quasisep",
                )
            }
            if arm == "joint"
            else None
        ),
    )
    return FittingProblem(model, datasets, seed=seed)


def direction_of(pmra: Any, pmdec: Any) -> Any:
    """``(pmra + pmdec)/sqrt(2)``: the diagonal of the proper-motion plane.

    The projection along which the two channels' errors *add*, and therefore
    the one a cross-channel systematic corrupts. See :data:`SBC_DIRECTION`.
    """
    return (np.asarray(pmra, dtype=float) + np.asarray(pmdec, dtype=float)) / np.sqrt(2.0)


def _thinned(values: Any, draws: int) -> np.ndarray:
    """*draws* posterior samples, taken evenly out of a run's chains.

    Evenly rather than from the head, which is what breaks the autocorrelation
    an MCMC run's neighbouring draws carry --- ``ampere.results.sbc``'s own
    rule, applied here because this study ranks a *derived* quantity and so
    walks the runs itself.
    """
    flat = np.asarray(values, dtype=float).reshape(-1)
    if flat.size < draws:
        raise ValueError(
            f"a rank against {draws} draw(s) needs at least that many, got {flat.size}."
        )
    return flat[np.linspace(0, flat.size - 1, draws).astype(int)]


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
    heteroscedastic: bool = False,
) -> Any:
    """Simulation-based calibration of one arm, against the **joint** generator.

    Talts et al. (2018) by refitting, in ``ampere.results.sbc``'s own shape and
    returning its own ``calibration`` group (built with the library's
    :func:`~ampere.results.calibration_dataset`, so the schema is one schema).
    The loop is written out here rather than delegated for one reason: what
    this study ranks is a **derived** quantity, ``(pmra + pmdec)/sqrt(2)``, and
    ``sbc`` ranks posterior *variables*. See :data:`SBC_DIRECTION` for why the
    derived one is where the answer is.

    The simulating problem is always the joint one, so every replicate's data
    carry a correlated centroiding systematic drawn from ``B (x) K_x`` --- the
    correlated draw ``DatasetCollection.draw_group`` makes, not two marginal
    ones. What changes between arms is only what is *fitted*: the same data,
    scored by a model that knows about the cross-channel structure or by one
    that does not.

    ``heteroscedastic=True`` (**W5.24**) runs the same study on channels with
    their own per-epoch sigmas: the simulating problem's containers carry them,
    so ``draw_group`` draws each replicate's white noise at them and the joint
    arm scores it on the dense route.
    """
    from ampere.inference import EmceeEngine
    from ampere.results import REFIT_ROUTE, calibration_dataset, replace_observations

    simulating = sbc_problem(backend, arm="joint", seed=seed, heteroscedastic=heteroscedastic)
    rng = np.random.default_rng(seed)
    rows: list[list[int]] = []
    for index in range(int(count)):
        theta = simulating.sample_prior(rng)
        simulation = simulating.simulate(theta, observe=True)
        if simulation.failed or simulation.observations is None:
            continue
        replica = replace_observations(simulating, simulation.observations, seed=seed + index)
        fitted = (
            replica
            if arm == "joint"
            else sbc_problem(
                backend,
                arm=arm,
                observed=(replica.datasets["ra"].observed, replica.datasets["dec"].observed),
                seed=seed + index,
            )
        )
        run = EmceeEngine(fitted, walkers=walkers).run(steps, burn_in=burn_in, progress=False)
        posterior = run["posterior"].dataset
        row = [
            int(np.sum(_thinned(posterior[name], draws) < float(theta[name])))
            for name in SBC_PARAMETERS
        ]
        drawn = direction_of(
            _thinned(posterior["model.pmra"], draws), _thinned(posterior["model.pmdec"], draws)
        )
        row.append(
            int(np.sum(drawn < float(direction_of(theta["model.pmra"], theta["model.pmdec"]))))
        )
        rows.append(row)
    if not rows:
        raise RuntimeError("no simulation produced a usable fit.")
    return calibration_dataset(
        np.asarray(rows, dtype=int),
        [*SBC_PARAMETERS, SBC_DIRECTION],
        posterior_draws=int(draws),
        route=REFIT_ROUTE,
        attrs={
            "ampere_calibration_label": (
                f"astrometry {arm} arm" + (", heteroscedastic" if heteroscedastic else "")
            )
        },
    )


def coverage_at(calibration: Any, level: float = 0.9, parameter: str | None = None) -> float:
    """Empirical coverage at a nominal *level*, for one parameter or averaged.

    One number, because the claim is one claim: "the central *level* interval
    contains the truth *level* of the time". *parameter* names one of
    :data:`SBC_PARAMETERS` or :data:`SBC_DIRECTION` --- pass the latter for the
    one the comparison turns on --- and ``None`` averages over all of them. The curve is
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
        "--heteroscedastic",
        action="store_true",
        help="a sigma per channel per epoch; the joint arm takes the dense route (W5.24)",
    )
    parser.add_argument(
        "--wide-prior",
        action="store_true",
        help=(
            "the period prior W4.9 measured as multi-modal, loguniform(50, 2000); "
            "reaches for dynesty when --engine is not also given (W5.15)"
        ),
    )
    parser.add_argument(
        "--engine",
        choices=list(ENGINES),
        default=None,
        help="default: emcee on reference, NUTS on torch/jax; the three nested samplers (W5.15)",
    )
    parser.add_argument("--live-points", type=int, default=None, help="nested samplers only")
    parser.add_argument(
        "--dlogz", type=float, default=None, help="dynesty/ultranest's stopping criterion"
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
            args.backend,
            arm=args.sbc,
            count=args.sbc_count,
            draws=args.sbc_draws,
            seed=args.seed,
            heteroscedastic=args.heteroscedastic,
        )
        elapsed = time.perf_counter() - started
        print(f"SBC, {args.sbc} arm: {args.sbc_count} simulation(s), {args.sbc_draws} draw(s)")
        for level in (0.5, 0.9, 0.95):
            direction = coverage_at(calibration, level, SBC_DIRECTION)
            print(
                f"  coverage at {level:.2f}: {coverage_at(calibration, level):.3f} "
                f"(mean), {direction:.3f} ({SBC_DIRECTION})"
            )
        pvalues = np.asarray(calibration["ks_pvalue"], dtype=float)
        names = [*SBC_PARAMETERS, SBC_DIRECTION]
        for name, pvalue in zip(names, pvalues, strict=True):
            print(f"  rank uniformity p({name}) = {pvalue:.3f}")
        print(f"  {elapsed:.1f} s wall clock")
        return 0

    problem = build_problem(
        args.backend,
        gp=args.gp,
        joint=args.joint,
        heteroscedastic=args.heteroscedastic,
        wide_prior=args.wide_prior,
        seed=args.seed,
    )
    print(f"negotiated channels: {list(problem.requirements['model'])}")
    for channel_name in problem.requirements["model"]:
        req = problem.requirements["model"][channel_name]
        print(f"sources asking of channel {channel_name!r}: {req.sources}")
    print(f"free parameters: {problem.parameters.free_names}")

    engine = args.engine
    if args.wide_prior and engine is None:
        engine = "dynesty"
        print(
            "--wide-prior with no --engine: reaching for dynesty -- nested sampling is the "
            "second remedy for a period prior wide enough to alias (W5.15)"
        )

    started = time.perf_counter()
    run = fit(
        problem,
        backend=args.backend,
        engine=engine,
        live_points=args.live_points,
        dlogz=args.dlogz,
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

    if args.wide_prior:
        modes = period_modes(run)
        print(f"  period modes ({len(modes)} found):")
        for mode in sorted(modes, key=lambda mode: mode["mass_fraction"], reverse=True):
            evidence_bit = f"  ln Z_k {mode['log_evidence']:+.3f}" if "log_evidence" in mode else ""
            print(
                f"    period {mode['median']:8.3f}  16/84%[{mode['lower']:8.3f}, "
                f"{mode['upper']:8.3f}]  f_k={mode['mass_fraction']:.3f}{evidence_bit}"
            )
    return 0
