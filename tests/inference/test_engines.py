"""The engine drivers, end to end: emcee, dynesty and zeus over §4.5's surface.

Five things are under test here, and the third is the one this item exists for.

1. **Every engine runs the same problem.** The two-dataset joint fit of
   ``inference.md`` §15 — one model, two channels, an instrument chain each and
   a calibration nuisance tied across them — is sampled by all three, and the
   rows asserted of the result are identical for all three. A driver that
   needed the problem shaped its way would fail here.
2. **Every run emits the run.** The ArviZ ``DataTree``, with the per-draw
   ``log_prior``/``log_likelihood`` split, the per-dataset decomposition, the
   observed data and the provenance attrs, and a netCDF round trip that keeps
   the hashes intact by value.
3. **The three agree on a known posterior.** The agreement problem is
   deliberately conjugate — a power law with its index held fixed is *linear*
   in ``norm``, so a Gaussian prior and Gaussian noise give a Gaussian
   posterior whose mean and variance are written down in :func:`analytic`
   rather than sampled. Each engine is held to that, and then to each other.
   Comparing samplers only against each other would pass three samplers that
   are wrong in the same way.
4. **The declared failure path is consumed as declared.** A model that raises
   its declared ``simulator_failures`` type over part of the prior produces
   recorded failures, a completed run and a surfaced summary — and the *same*
   model without the declaration kills the run, which is the documented sharp
   edge asserted rather than described.
5. **The drivers are backend-neutral.** No module under ``ampere.inference``
   imports ``ampere.backends`` — checked by parsing the import graph, not by
   grepping the prose — and importing the namespace does not pull a backend in.

Budgets, seeds and the trade-off
--------------------------------
Every run here is small and seeded: the module is a couple of minutes, which is
what lets it sit in the per-PR gate beside the other new-namespace suites
rather than on a nightly. That costs statistical power, so the agreement
tolerances are stated in units of the analytic posterior's own standard
deviation and are set several times the Monte Carlo standard error at these
budgets — loose enough that a correct sampler passes essentially always, and
far tighter than the errors a driver bug produces, since mis-transposed draws
or an undiscarded burn-in move a summary by whole standard deviations rather
than fractions of one.

Seeds are fixed throughout (``FittingProblem(seed=...)`` drives every stream
the drivers use), so a failure here is reproducible rather than a flake. The
runs are module- or class-scoped fixtures, so each happens once.

The expensive rows use the one-dimensional conjugate problem and the cheap
structural rows use the joint one, which is the opposite of the obvious
arrangement and is deliberate: structure needs a *complicated* problem and a
handful of draws, while agreement needs a *simple* problem and many.
"""

from __future__ import annotations

import ast
import json
import math
import subprocess
import sys
import warnings
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any, NamedTuple

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

import ampere.inference
from ampere.backends.reference import CalibrationScale, PowerLaw, Resample
from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Matern32,
    Model,
    Parameter,
    PoissonFamily,
    Spectrum,
    Tie,
)
from ampere.core.exceptions import LikelihoodError, OptionalDependencyError
from ampere.inference import (
    DynestyEngine,
    EmceeEngine,
    EngineError,
    SamplingFailureWarning,
    ZeusEngine,
)
from ampere.results import from_netcdf, to_netcdf

REFERENCE_WAVELENGTH = 1.0
TRUTH = {"norm": 2.0, "index": -1.2, "calibration": 1.0}
SEED = 20260905

FINE = np.geomspace(1.0, 20.0, 40)
COARSE = np.geomspace(1.5, 15.0, 12)


# ---------------------------------------------------------------------------
# The problems
# ---------------------------------------------------------------------------


def power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def noisy(grid: np.ndarray, values: np.ndarray, sigma: float, seed: int) -> Spectrum:
    """One observed spectrum: the truth plus Gaussian noise from a fixed stream."""
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, sigma, values.size)) * u.Jy,
        uncertainty=np.full(values.size, sigma) * u.Jy,
    )


BLUE_DATA = noisy(FINE, power_law(FINE, TRUTH["norm"], TRUTH["index"]), 0.08, seed=11)
RED_DATA = noisy(COARSE, power_law(COARSE, TRUTH["norm"], TRUTH["index"]), 0.05, seed=13)


def joint_problem(seed: int | None = SEED) -> FittingProblem:
    """``inference.md`` §15's shape, on the shipped reference backend.

    One power law published on two channels and observed twice: directly on the
    model's own grid, and resampled onto a coarser one — genuinely different
    instrument chains rather than the same fit written down twice. The two
    calibration factors are tied, so the problem has three free dimensions and
    not four, which is what makes it a joint fit rather than two fits sharing a
    figure.

    Every prior has the support its parameter has (a lognormal on both the
    normalisation and the calibration factor), which is the documented remedy
    for the trap in ``ampere.inference``'s docstring rather than an accident.
    """
    model = PowerLaw(
        FINE,
        norm=st.lognorm(0.4, scale=2.0),
        index=st.norm(-1.2, 0.3),
        reference_wavelength=REFERENCE_WAVELENGTH,
        channels=("blue", "red"),
    )
    return FittingProblem(
        model,
        DatasetCollection(
            {
                "blue": Dataset(
                    BLUE_DATA,
                    Instrument(
                        [CalibrationScale(st.lognorm(0.05), label="calibration")],
                        channel="blue",
                    ),
                ),
                "red": Dataset(
                    RED_DATA,
                    Instrument(
                        [
                            Resample(COARSE),
                            CalibrationScale(st.lognorm(0.05), label="calibration"),
                        ],
                        channel="red",
                    ),
                ),
            }
        ),
        ties=[
            Tie(
                "calibration",
                (
                    "blue.instrument.calibration.scale",
                    "red.instrument.calibration.scale",
                ),
            )
        ],
        seed=seed,
    )


# -- the conjugate problem the engines are held to --------------------------

AGREEMENT_GRID = np.geomspace(1.0, 10.0, 20)
AGREEMENT_SIGMA = 0.2
AGREEMENT_INDEX = -1.0
AGREEMENT_PRIOR = (2.0, 0.5)  # (mean, sd) of the Gaussian prior on `norm`
AGREEMENT_DATA = noisy(
    AGREEMENT_GRID, power_law(AGREEMENT_GRID, 2.0, AGREEMENT_INDEX), AGREEMENT_SIGMA, seed=7
)


class Analytic(NamedTuple):
    mean: float
    sd: float


def analytic() -> Analytic:
    """The posterior on ``norm``, in closed form.

    With ``index`` fixed the power law is *linear* in ``norm`` — the prediction
    is ``norm * x`` for a known design vector ``x`` — so a Gaussian prior and
    known Gaussian noise conjugate exactly:

    ``1/s^2 = 1/tau^2 + sum(x_i^2 / sigma_i^2)`` and
    ``m = s^2 (mu / tau^2 + sum(x_i y_i / sigma_i^2))``.

    Written out here rather than estimated from one long reference run, because
    an oracle a sampler produced is not an oracle.
    """
    x = (AGREEMENT_GRID / REFERENCE_WAVELENGTH) ** AGREEMENT_INDEX
    y = np.asarray(AGREEMENT_DATA.values)
    mu, tau = AGREEMENT_PRIOR
    precision = 1.0 / tau**2 + float(np.sum(x**2)) / AGREEMENT_SIGMA**2
    mean = (mu / tau**2 + float(np.sum(x * y)) / AGREEMENT_SIGMA**2) / precision
    return Analytic(mean=mean, sd=1.0 / math.sqrt(precision))


def agreement_problem(seed: int | None = SEED) -> FittingProblem:
    """One dataset, one free parameter, an exactly known posterior."""
    return FittingProblem(
        PowerLaw(
            AGREEMENT_GRID,
            norm=st.norm(*AGREEMENT_PRIOR),
            index=AGREEMENT_INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(AGREEMENT_DATA)],
        seed=seed,
    )


# ---------------------------------------------------------------------------
# Driving the three engines through one interface
# ---------------------------------------------------------------------------

ENGINES = ("emcee", "dynesty", "zeus")

#: Constructors, for the rows that check refusals rather than runs.
FACTORIES: Mapping[str, Callable[..., Any]] = {
    "emcee": EmceeEngine,
    "dynesty": DynestyEngine,
    "zeus": ZeusEngine,
}

#: Structural budget: enough draws to have a shape and a netCDF file, not
#: enough to have a posterior. Run on the (expensive) joint problem.
SMOKE: Mapping[str, Callable[[FittingProblem], Any]] = {
    "emcee": lambda problem: EmceeEngine(problem, walkers=8).run(steps=50, burn_in=10),
    "dynesty": lambda problem: DynestyEngine(problem, live_points=40).run(maxcall=1200),
    "zeus": lambda problem: ZeusEngine(problem, walkers=8).run(steps=20, burn_in=5),
}

#: Agreement budget: enough draws for a mean and a standard deviation to mean
#: something. Run on the (cheap) one-dimensional conjugate problem.
AGREEMENT: Mapping[str, Callable[[FittingProblem], Any]] = {
    "emcee": lambda problem: EmceeEngine(problem, walkers=16).run(steps=350, burn_in=100),
    "dynesty": lambda problem: DynestyEngine(problem, live_points=120).run(),
    "zeus": lambda problem: ZeusEngine(problem, walkers=12).run(steps=80, burn_in=25),
}

#: Reproducibility budget: the smallest run that still exercises every stream a
#: driver derives from the problem's seed.
TINY: Mapping[str, Callable[[FittingProblem], Any]] = {
    "emcee": lambda problem: EmceeEngine(problem, walkers=8).run(steps=20, burn_in=5),
    "dynesty": lambda problem: DynestyEngine(problem, live_points=30).run(maxcall=600),
    "zeus": lambda problem: ZeusEngine(problem, walkers=8).run(steps=10, burn_in=2),
}


@pytest.fixture(scope="module")
def joint_runs() -> dict[str, Any]:
    """One short run per engine on the joint problem, shared by the rows below."""
    return {name: SMOKE[name](joint_problem()) for name in ENGINES}


@pytest.fixture(scope="module")
def agreement_runs() -> dict[str, Any]:
    """One long-enough run per engine on the conjugate problem."""
    return {name: AGREEMENT[name](agreement_problem()) for name in ENGINES}


def run_stats(run: Any) -> dict[str, np.ndarray]:
    """``sample_stats`` as plain arrays, which is all these rows need."""
    dataset = run["sample_stats"].dataset
    return {name: np.asarray(dataset[name]) for name in ("lp", "log_prior", "log_likelihood")}


# ---------------------------------------------------------------------------
# 1. The toy two-dataset joint problem samples end to end on all three
# ---------------------------------------------------------------------------


class TestTheJointProblemSamples:
    """Accept criterion 1, and the emission obligations that ride with it."""

    def test_the_problem_is_the_one_the_contract_describes(self) -> None:
        """Three free dimensions, not four: the tie costs one (``inference.md`` §15)."""
        problem = joint_problem()
        assert problem.free_labels() == ("model.norm", "model.index", "calibration")
        assert problem.tied_names == ("calibration",)
        assert sorted(problem.datasets) == ["blue", "red"]
        assert len(problem.sites()["calibration"]) == 2

    @pytest.mark.parametrize("engine", ENGINES)
    def test_it_emits_the_arviz_groups(self, engine: str, joint_runs: dict[str, Any]) -> None:
        run = joint_runs[engine]
        assert {"posterior", "sample_stats", "log_likelihood", "observed_data"} <= set(run.children)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == {"model.norm", "model.index", "calibration"}
        assert posterior.sizes["draw"] > 0

    @pytest.mark.parametrize("engine", ENGINES)
    def test_every_draw_carries_the_split(self, engine: str, joint_runs: dict[str, Any]) -> None:
        """``DEVELOPMENT_PLAN.md`` §4.6: per-draw ``log_likelihood`` and ``log_prior``."""
        stats = run_stats(joint_runs[engine])
        assert np.all(np.isfinite(stats["log_prior"]))
        assert np.all(np.isfinite(stats["log_likelihood"]))
        assert stats["lp"] == pytest.approx(stats["log_prior"] + stats["log_likelihood"])

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_log_likelihood_group_decomposes_per_dataset(
        self, engine: str, joint_runs: dict[str, Any]
    ) -> None:
        """``results.md`` §6 / ``inference.md`` §18(a): one term per dataset, and they sum."""
        run = joint_runs[engine]
        group = run["log_likelihood"].dataset
        assert set(group.data_vars) == {"blue", "red"}
        total = np.asarray(group["blue"]) + np.asarray(group["red"])
        assert total == pytest.approx(run_stats(run)["log_likelihood"])
        assert run["log_likelihood"].attrs["ampere_decomposition"] == "per_dataset"

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_stored_draws_and_their_evaluations_belong_together(
        self, engine: str, joint_runs: dict[str, Any]
    ) -> None:
        """Re-score stored draws: the split must be the one filed against them.

        This is the row that catches the mistake an ensemble driver is most
        likely to make — ``(step, walker)`` stored as ``(walker, step)`` leaves
        every shape plausible and pairs each draw with somebody else's
        log-probability.
        """
        run = joint_runs[engine]
        problem = joint_problem()
        posterior = run["posterior"].dataset
        stats = run_stats(run)
        chains, draws = stats["lp"].shape
        rng = np.random.default_rng(3)
        for _ in range(5):
            c, d = int(rng.integers(chains)), int(rng.integers(draws))
            theta = np.array(
                [float(posterior[name][c, d]) for name in problem.parameters.free_names]
            )
            evaluation = problem.evaluate(theta)
            assert evaluation.log_prob == pytest.approx(stats["lp"][c, d], rel=1e-10)
            assert evaluation.log_prior == pytest.approx(stats["log_prior"][c, d], rel=1e-10)

    @pytest.mark.parametrize("engine", ENGINES)
    def test_nothing_was_recomputed(self, engine: str, joint_runs: dict[str, Any]) -> None:
        """The evaluations came from the sampler's own calls, not a second pass.

        Not merely an efficiency claim: for a *stochastic* model a recomputed
        log-likelihood is not the one the sampler accepted on, so the cache
        hitting is what makes the stored ``lp`` the number the run was actually
        driven by. At these budgets the cache cannot overflow, so any
        recomputation here means a driver looked up a θ it never scored.
        """
        assert joint_runs[engine].attrs["ampere_engine_draws_recomputed"] == 0

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_observed_data_travel_with_the_run(
        self, engine: str, joint_runs: dict[str, Any]
    ) -> None:
        observed = joint_runs[engine]["observed_data"].dataset
        assert set(observed.data_vars) == {"blue", "red"}
        assert observed["blue"].shape == (FINE.size,)
        assert observed["red"].shape == (COARSE.size,)


class TestProvenance:
    """The attrs every run carries (``DEVELOPMENT_PLAN.md`` §4.6, ``results.md`` §9)."""

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_run_says_what_produced_it(self, engine: str, joint_runs: dict[str, Any]) -> None:
        attrs = joint_runs[engine].attrs
        assert attrs["ampere_engine"] == engine
        # W2.12: derived from the problem's own pieces, not declared by the
        # driver -- ampere_backend is a fact about the run.
        assert attrs["ampere_backend"] == "reference"
        assert attrs["ampere_seed"] == SEED
        assert attrs["ampere_seed_source"] == "explicit"
        assert json.loads(attrs["ampere_free_names"]) == [
            "model.norm",
            "model.index",
            "calibration",
        ]
        assert "numpy" in json.loads(attrs["ampere_library_versions"])
        # The engine's own settings ride in the same attrs, so an archived run
        # says how it was configured and not only by what.
        assert any(key.startswith(f"ampere_{engine}_") for key in attrs)

    @pytest.mark.parametrize("engine", ENGINES)
    def test_it_carries_the_hashes_that_identify_the_composition(
        self, engine: str, joint_runs: dict[str, Any]
    ) -> None:
        attrs = joint_runs[engine].attrs
        assert len(attrs["ampere_spec_hash"]) == 32
        assert len(attrs["ampere_data_hash"]) == 32
        assert len(attrs["ampere_problem_hash"]) == 32
        identities = json.loads(attrs["ampere_model_identity_hashes"])
        assert sorted(identities) == ["model"]
        assert len(identities["model"]) == 32

    def test_the_three_engines_agree_about_what_was_fitted(
        self, joint_runs: dict[str, Any]
    ) -> None:
        """The spec and data hashes are properties of the problem, not the run.

        A run's identity has two halves and they must not be confused: *what
        was fitted* is the same whichever engine fitted it, while *which
        engine* is exactly what ``ampere_engine`` records.
        """
        hashes = {
            name: (
                joint_runs[name].attrs["ampere_spec_hash"],
                joint_runs[name].attrs["ampere_data_hash"],
                joint_runs[name].attrs["ampere_problem_hash"],
            )
            for name in ENGINES
        }
        assert len(set(hashes.values())) == 1


class TestBackendIdentity:
    """W2.12: the driver reads the backend off the problem and cannot assert one.

    Before W2.12 ``Engine(problem, backend=...)`` declared it, so
    ``ampere_backend`` recorded whatever the caller typed. The argument is gone
    (pre-release, so no shim), and the flag is aggregated from what the models
    and instrument steps declare.
    """

    @pytest.mark.parametrize("engine", ENGINES)
    def test_a_driver_refuses_a_declared_backend(self, engine: str) -> None:
        # emcee and dynesty simply have no such keyword, so Python refuses it.
        # zeus forwards **sampler_settings to its sampler, which would have
        # swallowed the name silently, so it refuses explicitly and says why.
        with pytest.raises((TypeError, EngineError), match="backend"):
            FACTORIES[engine](joint_problem(), backend="torch")

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_driver_reports_what_the_problem_declares(self, engine: str) -> None:
        problem = joint_problem()
        assert FACTORIES[engine](problem).backend == problem.backend == "reference"

    def test_the_emitted_attrs_follow_the_problem_not_the_driver(self) -> None:
        # A problem whose pieces declare a different backend emits that one,
        # with nothing said anywhere in the driver: the derivation is the whole
        # mechanism.
        class NativeCalibrationScale(CalibrationScale):
            BACKEND = "mirror"

        class NativePowerLaw(PowerLaw):
            BACKEND = "mirror"

        class NativeNoise(IndependentNoise):
            # W2.13: the noise model is a capability part too, so the
            # dataset's default (which declares "reference") would make this
            # a two-backend problem rather than a "mirror" one.
            BACKEND = "mirror"

        model = NativePowerLaw(
            FINE,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        )
        problem = FittingProblem(
            model,
            [
                Dataset(
                    BLUE_DATA,
                    Instrument([NativeCalibrationScale(st.lognorm(0.05), label="calibration")]),
                    Likelihood(GaussianFamily(), NativeNoise()),
                )
            ],
            seed=SEED,
        )
        assert problem.backend == "mirror"
        run = EmceeEngine(problem, walkers=8).run(steps=20, burn_in=5)
        assert run.attrs["ampere_backend"] == "mirror"


class TestNetcdfRoundTrip:
    """Accept criterion 3: the stored run round-trips with hashes intact."""

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_hashes_survive_by_value(
        self, engine: str, joint_runs: dict[str, Any], tmp_path: Path
    ) -> None:
        run = joint_runs[engine]
        path = tmp_path / f"{engine}.nc"
        to_netcdf(run, path)
        back = from_netcdf(path)
        for key in (
            "ampere_spec_hash",
            "ampere_data_hash",
            "ampere_problem_hash",
            "ampere_model_identity_hashes",
            "ampere_component_spec_hashes",
            "ampere_data_hashes",
            "ampere_engine",
            "ampere_seed",
        ):
            assert back.attrs[key] == run.attrs[key], key

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_draws_and_the_split_survive(
        self, engine: str, joint_runs: dict[str, Any], tmp_path: Path
    ) -> None:
        run = joint_runs[engine]
        path = tmp_path / f"{engine}-values.nc"
        to_netcdf(run, path)
        back = from_netcdf(path)
        for name in ("model.norm", "model.index", "calibration"):
            assert np.asarray(back["posterior"][name]) == pytest.approx(
                np.asarray(run["posterior"][name])
            )
        assert np.asarray(back["sample_stats"]["lp"]) == pytest.approx(run_stats(run)["lp"])
        assert np.asarray(back["log_likelihood"]["blue"]) == pytest.approx(
            np.asarray(run["log_likelihood"]["blue"])
        )


# ---------------------------------------------------------------------------
# 2. Posterior summaries within tolerance — of the truth, and of each other
# ---------------------------------------------------------------------------

#: How far a run's posterior mean may sit from the analytic one, in units of
#: the analytic standard deviation. At these budgets the Monte Carlo standard
#: error of the mean is comfortably under 0.1 sd for all three engines, so this
#: is a wide margin against a flake and a narrow one against a driver bug: the
#: smallest real mistakes here — a burn-in not discarded, draws paired with the
#: wrong evaluations — move the mean by whole standard deviations.
MEAN_TOLERANCE = 0.25

#: The same for the standard deviation, as a fraction. Same width, and it buys
#: less, because a spread estimated from correlated draws is the noisier of the
#: two summaries.
SD_TOLERANCE = 0.25


def posterior_summary(run: Any) -> tuple[float, float]:
    values = np.asarray(run["posterior"]["model.norm"]).ravel()
    return float(values.mean()), float(values.std(ddof=1))


class TestPosteriorAgreement:
    """Accept criterion 2, against an oracle rather than only against each other."""

    def test_the_analytic_posterior_is_the_one_the_problem_declares(self) -> None:
        """Guard the oracle: ``log_prob`` must be the Gaussian :func:`analytic` claims.

        Without this row the agreement tests would still pass if the formula
        and the problem had drifted apart — three samplers agreeing with each
        other, and all three compared against the wrong number.
        """
        problem = agreement_problem()
        truth = analytic()
        grid = np.array([truth.mean - truth.sd, truth.mean, truth.mean + truth.sd])
        computed = np.array([problem.log_prob([value]) for value in grid])
        expected = st.norm(truth.mean, truth.sd).logpdf(grid)
        # Equal up to the normalising constant, which log_prob does not carry.
        offset = computed[1] - expected[1]
        assert (computed - expected) == pytest.approx(np.full(3, offset), abs=1e-8)

    @pytest.mark.parametrize("engine", ENGINES)
    def test_each_engine_recovers_the_analytic_posterior(
        self, engine: str, agreement_runs: dict[str, Any]
    ) -> None:
        truth = analytic()
        mean, sd = posterior_summary(agreement_runs[engine])
        assert abs(mean - truth.mean) < MEAN_TOLERANCE * truth.sd
        assert abs(sd / truth.sd - 1.0) < SD_TOLERANCE

    @pytest.mark.parametrize(("left", "right"), [("emcee", "dynesty"), ("emcee", "zeus")])
    def test_the_engines_agree_with_each_other(
        self, left: str, right: str, agreement_runs: dict[str, Any]
    ) -> None:
        """Two independent errors could cancel against the oracle; they will not here."""
        truth = analytic()
        left_mean, left_sd = posterior_summary(agreement_runs[left])
        right_mean, right_sd = posterior_summary(agreement_runs[right])
        assert abs(left_mean - right_mean) < 2.0 * MEAN_TOLERANCE * truth.sd
        assert abs(left_sd / right_sd - 1.0) < 2.0 * SD_TOLERANCE

    def test_dynesty_also_reports_the_evidence(self, agreement_runs: dict[str, Any]) -> None:
        """The one summary the ensemble engines cannot give, and a finite one at that."""
        attrs = agreement_runs["dynesty"].attrs
        assert math.isfinite(attrs["ampere_dynesty_logz"])
        assert attrs["ampere_dynesty_logzerr"] > 0.0
        assert attrs["ampere_dynesty_dead_points"] > 0


class TestReproducibility:
    """``inference.md`` §12: one seed per run, named sub-streams derived from it."""

    @pytest.mark.parametrize("engine", ENGINES)
    def test_the_same_seed_gives_the_same_run(self, engine: str) -> None:
        one = TINY[engine](agreement_problem(seed=4242))
        two = TINY[engine](agreement_problem(seed=4242))
        assert np.asarray(one["posterior"]["model.norm"]) == pytest.approx(
            np.asarray(two["posterior"]["model.norm"])
        )

    @pytest.mark.parametrize("engine", ENGINES)
    def test_different_seeds_give_different_runs(self, engine: str) -> None:
        """The other half: a seed that is recorded but does nothing is worse than none.

        Compared shape-first, because a nested sampler's run length is itself
        seed-dependent — the number of dead points before the evidence
        criterion is met differs between seeds — so two dynesty runs can differ
        by being different *sizes* before they differ in any value.
        """
        one = np.asarray(TINY[engine](agreement_problem(seed=4242))["posterior"]["model.norm"])
        two = np.asarray(TINY[engine](agreement_problem(seed=99))["posterior"]["model.norm"])
        assert one.shape != two.shape or not np.allclose(one, two)

    def test_the_streams_are_named_per_engine_and_concern(self) -> None:
        """Initialisation must not share a stream with sampling (``inference.md`` §12)."""
        problem = agreement_problem()
        engine = EmceeEngine(problem, walkers=8)
        assert engine.stream("initialisation") is not engine.stream("sampler")
        assert engine.stream("initialisation") is problem.rng("emcee.initialisation")

    def test_an_unseeded_problem_still_runs(self) -> None:
        """``seed=None`` means nothing is reproducible, not that nothing works."""
        run = TINY["emcee"](agreement_problem(seed=None))
        assert run.attrs["ampere_seed_source"] == "entropy"
        assert "ampere_seed" not in run.attrs


# ---------------------------------------------------------------------------
# 3. Refusals: check_engine, and the settings an engine cannot honour
# ---------------------------------------------------------------------------


def latent_gp_problem() -> FittingProblem:
    """A Poisson count model with a GP noise term, needing a latent block.

    ``inference.md`` §10's own ``check_engine`` example: the marginalisation
    this composition needs is not one a gradient-free engine can deliver, so
    the refusal is the point.
    """
    grid = np.array([1.0, 2.0, 3.0])
    counts = Spectrum(grid * u.micron, np.array([4.0, 7.0, 2.0]))

    class Rate(Model):
        def __init__(self) -> None:
            self.register_buffer("grid", grid, unit=u.micron)
            self.register_parameter(Parameter("rate", st.loguniform(0.5, 50.0)))

        def evaluate(self, **values: Any) -> Any:
            context = self.context(values)
            return Spectrum(
                context["grid"] * u.micron,
                np.full(np.shape(context["grid"]), context["rate"]),
            )

    return FittingProblem(
        Rate(),
        [
            Dataset(
                counts,
                likelihood=Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
                label="counts",
            )
        ],
    )


class TestRefusals:
    @pytest.mark.parametrize("engine", ENGINES)
    def test_a_likelihood_no_gradient_free_engine_can_run_is_refused_at_construction(
        self, engine: str
    ) -> None:
        """``check_engine`` is called before any sampling, with ``observed=``.

        The drivers pass ``differentiable=False`` explicitly rather than
        letting it default to the problem's own capability: the question is
        "can this engine run this likelihood?", and its answer must not change
        because the problem happens to sit on a differentiable backend.
        """
        with pytest.raises(LikelihoodError, match=engine):
            FACTORIES[engine](latent_gp_problem())

    @pytest.mark.parametrize("engine", ENGINES)
    def test_a_problem_with_nothing_to_sample_is_refused(self, engine: str) -> None:
        problem = FittingProblem(
            PowerLaw(AGREEMENT_GRID, norm=2.0, index=AGREEMENT_INDEX),
            [Dataset(AGREEMENT_DATA)],
        )
        with pytest.raises(EngineError, match="nothing to sample"):
            FACTORIES[engine](problem)

    @pytest.mark.parametrize("engine", ["emcee", "zeus"])
    @pytest.mark.parametrize(
        ("walkers", "message"), [(3, "even"), (2, "2 x n_dim")], ids=["odd", "too-few"]
    )
    def test_the_walker_count_is_checked(self, engine: str, walkers: int, message: str) -> None:
        with pytest.raises(EngineError, match=message):
            FACTORIES[engine](joint_problem(), walkers=walkers)

    @pytest.mark.parametrize("engine", ["emcee", "zeus"])
    def test_a_burn_in_that_keeps_nothing_is_refused(self, engine: str) -> None:
        engine_object = FACTORIES[engine](agreement_problem(), walkers=8)
        with pytest.raises(EngineError, match="nothing to emit"):
            engine_object.run(steps=10, burn_in=10)

    def test_bad_start_positions_are_refused(self) -> None:
        with pytest.raises(EngineError, match="walkers, n_dim"):
            EmceeEngine(agreement_problem(), walkers=8).run(steps=5, initial=np.zeros((8, 3)))

    def test_dynesty_needs_live_points(self) -> None:
        with pytest.raises(EngineError, match="live points"):
            DynestyEngine(agreement_problem(), live_points=1)


# ---------------------------------------------------------------------------
# 4. Failure signalling, consumed as declared
# ---------------------------------------------------------------------------


class SimulatorCrash(RuntimeError):
    """Stands in for a wrapped external code exiting non-zero."""


def flaky_problem(*, declare: bool = True) -> FittingProblem:
    """A model that refuses part of its own prior's support.

    Deliberately the shape of the trap ``ampere.inference``'s documentation
    warns about — a prior whose support is wider than the model's — with the
    declared-exception remedy applied or not. ``simulator_failures=`` is the
    whole difference between a run that records failures and a run that dies,
    so the same problem is built both ways and both outcomes are asserted.
    """
    grid = np.array([1.0, 2.0, 3.0])
    observed = Spectrum(grid * u.micron, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.3, 0.3, 0.3] * u.Jy)

    class Fragile(Model):
        def __init__(self) -> None:
            self.register_buffer("grid", grid, unit=u.micron)
            self.register_parameter(Parameter("slope", st.norm(1.0, 2.0)))

        def evaluate(self, **values: Any) -> Any:
            context = self.context(values)
            slope = float(context["slope"])
            if slope <= 0.0:
                raise SimulatorCrash("the external code refuses a non-positive slope")
            return Spectrum(context["grid"] * u.micron, slope * context["grid"] * u.Jy)

    return FittingProblem(
        Fragile(),
        [Dataset(observed, label="d")],
        seed=SEED,
        simulator_failures=(SimulatorCrash,) if declare else (),
    )


class TestFailureSignalling:
    """``inference.md`` §11, consumed rather than reimplemented."""

    @pytest.fixture(scope="class")
    def flaky_run(self) -> tuple[Any, list[warnings.WarningMessage]]:
        engine = EmceeEngine(flaky_problem(), walkers=8)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            run = engine.run(steps=150, burn_in=50)
        return run, [w for w in caught if issubclass(w.category, SamplingFailureWarning)]

    def test_the_run_completes_rather_than_raising(
        self, flaky_run: tuple[Any, list[warnings.WarningMessage]]
    ) -> None:
        run, _ = flaky_run
        assert run["posterior"].dataset.sizes["draw"] == 100

    def test_the_failures_are_counted_and_named(
        self, flaky_run: tuple[Any, list[warnings.WarningMessage]]
    ) -> None:
        run, _ = flaky_run
        counts = json.loads(run.attrs["ampere_failure_counts"])
        assert counts.get("model_failed", 0) > 0
        recorded = json.loads(run.attrs["ampere_failures"])
        assert {entry["exception_type"] for entry in recorded} == {"SimulatorCrash"}
        assert {entry["where"] for entry in recorded} == {"model"}

    def test_the_summary_is_surfaced_three_ways(
        self, flaky_run: tuple[Any, list[warnings.WarningMessage]]
    ) -> None:
        """A warning for whoever is watching, an attr for the archive, an attribute for a script."""
        run, caught = flaky_run
        assert len(caught) == 1
        assert "model_failed" in str(caught[0].message)
        assert "strict=True" in str(caught[0].message)
        assert "model_failed" in run.attrs["ampere_failure_summary"]

    def test_one_warning_per_run_not_one_per_failed_draw(
        self, flaky_run: tuple[Any, list[warnings.WarningMessage]]
    ) -> None:
        """``inference.md`` §11: aggregated, so a high failure rate does not spam."""
        run, caught = flaky_run
        failures = sum(json.loads(run.attrs["ampere_failure_counts"]).values())
        assert failures > len(caught) == 1

    def test_every_stored_draw_is_still_scoreable(
        self, flaky_run: tuple[Any, list[warnings.WarningMessage]]
    ) -> None:
        """A stored draw is an accepted position, so it is never a failed one.

        The failures live in the proposals emcee rejected. What must hold of
        the run is that no stored draw carries a stale reason and every one has
        a finite log-probability.
        """
        run, _ = flaky_run
        stats = run["sample_stats"].dataset
        assert not np.any(np.asarray(stats["failed"]))
        assert np.all(np.isfinite(np.asarray(stats["lp"])))

    def test_a_clean_run_warns_about_nothing(self, joint_runs: dict[str, Any]) -> None:
        assert "ampere_failure_summary" not in joint_runs["emcee"].attrs
        assert json.loads(joint_runs["emcee"].attrs["ampere_failure_counts"]) == {}

    def test_the_undeclared_version_of_the_same_model_kills_the_run(self) -> None:
        """The other half of the documented trap, asserted rather than described.

        Without ``simulator_failures=``, the very same model stops the run —
        because the default catch set is narrow on purpose, and a bug turned
        into ``-inf`` is a fit that runs, converges and is wrong. This row is
        what makes the prior-support warning in ``ampere.inference``'s
        documentation a checked statement rather than folklore.
        """
        with pytest.raises(SimulatorCrash):
            EmceeEngine(flaky_problem(declare=False), walkers=8).run(steps=150)


# ---------------------------------------------------------------------------
# 5. The zeus extra, and the rule that keeps all three portable
# ---------------------------------------------------------------------------


class TestTheZeusExtra:
    def test_the_refusal_names_the_extra(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """zeus stays an extra (``architecture.md`` §3), so its absence must say so."""
        monkeypatch.setitem(sys.modules, "zeus", None)
        with pytest.raises(OptionalDependencyError) as raised:
            ZeusEngine(agreement_problem())
        assert raised.value.package == "zeus"
        assert raised.value.extra == "zeus"
        assert 'pip install "ampere[zeus]"' in str(raised.value)

    def test_the_refusal_comes_before_any_work(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """At construction, not mid-run: the user has not pressed go yet."""
        monkeypatch.setitem(sys.modules, "zeus", None)
        with pytest.raises(OptionalDependencyError):
            ZeusEngine(agreement_problem(), walkers=8)


class TestBackendNeutrality:
    """``inference.md`` §10's claim, enforced mechanically rather than by review.

    The engines are written "against §4.5's surface and nothing else", so they
    "must work unchanged with every backend". The check is that this namespace
    cannot *import* a backend: not a style preference, but what stops a driver
    acquiring a reference-backend-shaped assumption that the torch and jax
    tracks then have to work around.

    Parsed rather than grepped, deliberately. The documentation does name a
    backend — the prior-support trap is best explained with the concrete
    ``BlackBody`` that raises on a negative temperature — and a textual search
    would either forbid that or be defeated by a string-built import. The
    import graph is what actually constrains the code.
    """

    @pytest.mark.parametrize("module", sorted(Path(ampere.inference.__file__).parent.glob("*.py")))
    def test_no_module_imports_a_backend(self, module: Path) -> None:
        tree = ast.parse(module.read_text(encoding="utf-8"), filename=str(module))
        imported: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module is not None:
                imported.append(node.module)
        offenders = [name for name in imported if name.split(".")[:2] == ["ampere", "backends"]]
        assert offenders == [], f"{module.name} imports {offenders}"

    def test_the_only_ampere_imports_are_core_and_results(self) -> None:
        """The positive half of the same claim: two dependencies, both contracts."""
        root = Path(ampere.inference.__file__).parent
        namespaces = set()
        for path in sorted(root.glob("*.py")):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            for node in ast.walk(tree):
                names = []
                if isinstance(node, ast.Import):
                    names = [alias.name for alias in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    names = [node.module]
                namespaces.update(
                    ".".join(name.split(".")[:2]) for name in names if name.startswith("ampere.")
                )
        assert namespaces == {"ampere.core", "ampere.results"}

    def test_importing_it_does_not_import_a_backend(self) -> None:
        probe = subprocess.run(
            [
                sys.executable,
                "-c",
                (
                    "import sys, ampere.inference; "
                    "print(any(name.startswith('ampere.backends') for name in sys.modules))"
                ),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        assert probe.stdout.strip() == "False"

    def test_a_hand_written_model_needs_no_backend_at_all(self) -> None:
        """The portability claim at its smallest: no ``ampere.backends`` in sight.

        The docstring examples make the same point; this row makes it a
        first-class assertion rather than a consequence of documentation.
        """
        run = TINY["emcee"](flaky_problem())
        assert run["posterior"].dataset.sizes["draw"] > 0
        assert run.attrs["ampere_engine"] == "emcee"
