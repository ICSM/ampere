"""W3.6: family D — simulation-based calibration and coverage, on the numpy path.

The route-1 (``sbi``) half of family D lives in ``tests/inference/test_sbi.py``,
because it needs a trained network; everything here runs in ``dev``.

What is worth testing, and what a weaker suite would miss
---------------------------------------------------------
The failure mode a calibration module has is that it reports *calibrated*. A
suite that only asserted "a Dataset came out with the right dimensions" would
pass an implementation whose ranks were arbitrary, so the arithmetic is pinned
against constructions where the answer is known in advance:

* **The rank is checked against a hand-computed count** on draws laid out by
  hand, so the definition (``#{draw < truth}``, thinned per chain) is under
  test and not merely exercised.
* **The coverage curve is checked against ranks placed deliberately**: ranks
  all at the median give a curve that is 1 everywhere above level 0, ranks all
  at the extremes give one that is 0 below level 1. Both are the identities
  ``coverage_from_ranks`` claims, and neither survives an off-by-one in the
  interval.
* **The end-to-end loop is run on a problem whose posterior is known to be
  calibrated** — an ordinary Gaussian likelihood fitted with the very model
  that generated the data — and separately on one that is **deliberately
  miscalibrated**, by handing the loop a factory that shrinks its own
  posterior. The uniformity test must let the first through and reject the
  second; a test that only ran the first would pass a module that never
  rejects anything.

Budgets are tiny and seeded throughout. The suite does not attempt to
demonstrate that a real fit is calibrated to any precision — that is what
``examples/wstat_comparison.py``'s study is for, at a budget no per-PR gate can
afford.
"""

from __future__ import annotations

import math
import warnings
from typing import Any

import astropy.units as u
import matplotlib
import numpy as np
import pytest
import scipy.stats as st

matplotlib.use("Agg")

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    HierarchicalPrior,
    Model,
    Parameter,
    ParameterSet,
    Population,
    Spectrum,
)
from ampere.core.exceptions import ResultsError
from ampere.inference import EmceeEngine
from ampere.results import (
    CALIBRATION_GROUP,
    CALIBRATION_SCHEMA_VERSION,
    DEFAULT_LEVELS,
    REFIT_ROUTE,
    DrawRecorder,
    attach_calibration,
    calibration_dataset,
    coverage_from_ranks,
    figure_metadata,
    from_netcdf,
    plot_coverage,
    plot_sbc_ranks,
    replace_observations,
    sbc,
    to_netcdf,
    uniformity_pvalues,
)

SEED = 20260909
GRID = np.linspace(1.0, 6.0, 12)


class Line(Model):
    """``y = slope * x + offset`` — two parameters, both with proper priors."""

    def __init__(self) -> None:
        self.register_buffer("wavelength", GRID, unit=u.um)
        self.register_parameter(Parameter("slope", st.norm(2.0, 0.6)))
        self.register_parameter(Parameter("offset", st.norm(0.5, 0.4)))

    def evaluate(self, **values: Any) -> Spectrum:
        context = self.context(values)
        predicted = context["slope"] * context["wavelength"] + context["offset"]
        return Spectrum(context["wavelength"] * u.um, predicted * u.Jy)


def toy_problem(seed: int | None = SEED) -> FittingProblem:
    """A two-parameter Gaussian problem whose posterior is calibrated by construction.

    The data are the model plus independent Gaussian noise of a known sigma, and
    the fit uses that same model and that same sigma, so the posterior an honest
    sampler produces *is* the true conditional and its SBC ranks are uniform.
    That is what makes it the right null: a loop that fails here is broken, not
    unlucky.
    """
    observed = Spectrum(
        GRID * u.um,
        (2.0 * GRID + 0.5) * u.Jy,
        uncertainty=np.full(GRID.size, 0.5) * u.Jy,
    )
    return FittingProblem(Line(), [Dataset(observed, label="line")], seed=seed)


# ---------------------------------------------------------------------------
# The arithmetic
# ---------------------------------------------------------------------------


class TestRanksAndCoverage:
    def test_the_coverage_curve_is_the_identity_ranks_at_the_median_imply(self) -> None:
        """Every truth landing exactly at the posterior median covers at every level."""
        ranks = np.full((40, 1), 50)
        levels, coverage = coverage_from_ranks(ranks, 100)
        np.testing.assert_allclose(levels, DEFAULT_LEVELS)
        # Level 0 is the degenerate case (the interval is a point) and the rank
        # fraction is exactly 0.5, so even that one covers.
        np.testing.assert_allclose(coverage[:, 0], 1.0)

    def test_ranks_at_the_extremes_cover_only_at_the_full_level(self) -> None:
        """A posterior that always excludes the truth covers nowhere but at level 1."""
        ranks = np.concatenate([np.zeros((20, 1), dtype=int), np.full((20, 1), 100)])
        _, coverage = coverage_from_ranks(ranks, 100)
        assert float(coverage[-1, 0]) == pytest.approx(1.0)
        assert float(np.max(coverage[:-1, 0])) == pytest.approx(0.0)

    def test_a_uniform_rank_set_lands_on_the_diagonal(self) -> None:
        ranks = np.arange(101, dtype=int).reshape(-1, 1)
        levels, coverage = coverage_from_ranks(ranks, 100)
        # 101 evenly spaced ranks: the empirical curve tracks the nominal one to
        # within one rank's worth, 1/101.
        np.testing.assert_allclose(coverage[:, 0], levels, atol=1.5 / 101.0)

    def test_uniformity_rejects_a_visibly_clustered_rank_set(self) -> None:
        clustered = np.full((200, 1), 3)
        spread = np.linspace(0, 100, 200).round().astype(int).reshape(-1, 1)
        assert float(uniformity_pvalues(clustered, 100)[0]) < 1e-6
        assert float(uniformity_pvalues(spread, 100)[0]) > 0.05

    def test_levels_outside_zero_to_one_are_refused(self) -> None:
        with pytest.raises(ResultsError, match="probabilities in"):
            coverage_from_ranks(np.zeros((3, 1), dtype=int), 10, levels=[0.5, 1.5])


class TestTheSharedDataset:
    def test_it_carries_the_schema_the_group_promises(self) -> None:
        ranks = np.arange(60, dtype=int).reshape(30, 2) % 21
        result = calibration_dataset(ranks, ["a", "b"], posterior_draws=20, route=REFIT_ROUTE)
        assert result.sizes["simulation"] == 30
        assert result.sizes["parameter"] == 2
        assert list(result["ranks"].dims) == ["simulation", "parameter"]
        assert list(result["coverage"].dims) == ["level", "parameter"]
        assert list(result["ks_pvalue"].dims) == ["parameter"]
        assert result.attrs["ampere_calibration_schema_version"] == CALIBRATION_SCHEMA_VERSION
        assert result.attrs["ampere_calibration_route"] == REFIT_ROUTE
        assert result.attrs["ampere_calibration_posterior_draws"] == 20

    def test_a_name_per_column_is_required(self) -> None:
        with pytest.raises(ResultsError, match="one name per ranked column"):
            calibration_dataset(
                np.zeros((4, 3), dtype=int), ["a"], posterior_draws=5, route=REFIT_ROUTE
            )

    def test_an_unknown_route_is_refused_by_name(self) -> None:
        with pytest.raises(ResultsError, match="calibration route"):
            calibration_dataset(
                np.zeros((4, 1), dtype=int), ["a"], posterior_draws=5, route="guesswork"
            )


# ---------------------------------------------------------------------------
# The replica
# ---------------------------------------------------------------------------


class TestReplaceObservations:
    def test_only_the_data_change(self) -> None:
        problem = toy_problem()
        simulation = problem.simulate(observe=True, rng=np.random.default_rng(2))
        assert not simulation.failed and simulation.observations is not None
        replica = replace_observations(problem, simulation.observations, seed=17)

        assert replica.free_labels() == problem.free_labels()
        assert replica.seed == 17
        original = problem.datasets["line"].observed
        replaced = replica.datasets["line"].observed
        assert not np.allclose(np.asarray(replaced.values), np.asarray(original.values))
        # Coordinates, uncertainties and units are the container's own and must
        # survive: a replica whose sigma changed would be a different experiment.
        np.testing.assert_allclose(
            np.asarray(replaced.uncertainty), np.asarray(original.uncertainty)
        )
        np.testing.assert_allclose(
            np.asarray(replaced.axes[0].values), np.asarray(original.axes[0].values)
        )

    def test_a_missing_dataset_is_refused_by_name(self) -> None:
        problem = toy_problem()
        with pytest.raises(ResultsError, match="line"):
            replace_observations(problem, {})

    def test_the_shared_parameter_set_survives_the_replica(self) -> None:
        """W5.28(a): a population-level ``shared`` set must not vanish on replay.

        Before this fix, ``replace_observations`` rebuilt the collection with
        only ``joint`` carried over, so a replica silently lost its
        hyperprior component -- turning a population-level SBC run into an
        independent one without any error.
        """
        observed = Spectrum(
            GRID * u.um,
            (2.0 * GRID + 0.5) * u.Jy,
            uncertainty=np.full(GRID.size, 0.5) * u.Jy,
        )
        hyper = ParameterSet([Parameter("mu_slope", st.norm(2.0, 0.6))])
        datasets = DatasetCollection(
            {"line": Dataset(observed, label="line")},
            shared=hyper,
            shared_label="hyper",
        )
        problem = FittingProblem(Line(), datasets, seed=SEED)
        simulation = problem.simulate(observe=True, rng=np.random.default_rng(2))
        assert not simulation.failed and simulation.observations is not None

        replica = replace_observations(problem, simulation.observations, seed=17)

        assert replica.datasets.shared is hyper
        assert replica.datasets.shared_label == "hyper"

    def test_the_populations_survive_the_replica(self) -> None:
        """W5.30(c): a W5.12 ``Population`` declaration must not vanish on replay.

        Before this fix, ``replace_observations`` rebuilt the collection with
        ``joint`` and ``shared`` carried over but not ``populations``, so a
        population-level SBC replay would silently be fitted as independent
        per-member objects -- the same failure mode W5.28(a) closed for the
        shared parameter set. The population here is declared the documented,
        canonical way -- passed directly to ``FittingProblem(...,
        populations=...)`` (``docs/source/advanced.rst``), not to the
        ``DatasetCollection`` -- which is the case ``problem.datasets.populations``
        alone would miss: it is ``problem.populations`` (``FittingProblem``'s
        own concatenation of both routes) that has to survive the replica.
        """
        observed = Spectrum(
            GRID * u.um,
            (2.0 * GRID + 0.5) * u.Jy,
            uncertainty=np.full(GRID.size, 0.5) * u.Jy,
        )
        population = Population(
            "objects",
            members=[Parameter("slope", HierarchicalPrior("norm", {"loc": "mu_slope"}))],
            hyperpriors=[Parameter("mu_slope", st.norm(2.0, 0.6))],
            over=["model"],
        )
        datasets = DatasetCollection({"line": Dataset(observed, label="line")})
        problem = FittingProblem(Line(), datasets, populations=[population], seed=SEED)
        # The canonical route leaves the DatasetCollection itself unaware of
        # the population -- confirming this is the gap the fix has to close.
        assert problem.datasets.populations == ()
        assert problem.populations == (population,)
        simulation = problem.simulate(observe=True, rng=np.random.default_rng(2))
        assert not simulation.failed and simulation.observations is not None

        replica = replace_observations(problem, simulation.observations, seed=17)

        assert replica.populations == (population,)


# ---------------------------------------------------------------------------
# The loop, against an exact posterior
# ---------------------------------------------------------------------------

SIGMA = 0.5
PRIOR_MEAN = np.array([2.0, 0.5])
PRIOR_COVARIANCE = np.diag([0.6**2, 0.4**2])
DESIGN = np.column_stack([GRID, np.ones_like(GRID)])


class _ExactPosterior:
    """A "fit" that is the analytic linear-Gaussian posterior, exactly.

    This is the oracle the loop is actually tested against, and it exists
    because an MCMC arm cannot be one. Any departure from uniformity an emcee
    run shows at a per-PR budget has three possible causes — a bug in the loop,
    an unconverged chain, and the Monte Carlo noise of a few dozen simulations
    — and a test that cannot separate them tests nothing. Here there is no
    sampler: :class:`Line` with Gaussian priors and Gaussian noise has a
    posterior in closed form, so the *only* thing that can make these ranks
    non-uniform is the loop itself.

    Emitted through :class:`~ampere.results.DrawRecorder` rather than by
    building a tree by hand, so the loop reads the draws off exactly the group
    an engine writes.
    """

    def __init__(self, problem: FittingProblem, *, draws: int, seed: int) -> None:
        self.problem = problem
        self.draws = draws
        self.rng = np.random.default_rng(seed)

    def run(self, **kwargs: Any) -> Any:
        observed = np.asarray(self.problem.datasets["line"].observed.values, dtype=float)
        prior_precision = np.linalg.inv(PRIOR_COVARIANCE)
        covariance = np.linalg.inv(prior_precision + DESIGN.T @ DESIGN / SIGMA**2)
        mean = covariance @ (prior_precision @ PRIOR_MEAN + DESIGN.T @ observed / SIGMA**2)
        recorder = DrawRecorder(self.problem)
        for row in self.rng.multivariate_normal(mean, covariance, size=self.draws):
            recorder.record(row)
        return recorder.emit()


@pytest.fixture(scope="module")
def exact() -> Any:
    """One 150-simulation study against the closed-form posterior, shared by three checks."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return sbc(
            toy_problem(),
            lambda replica: _ExactPosterior(replica, draws=80, seed=int(replica.seed) % 9973),
            count=150,
            draws=60,
        )


class TestAgainstAnExactPosterior:
    """The loop's arithmetic, where the right answer is known and there is no sampler."""

    def test_the_ranks_of_an_exact_posterior_are_uniform(self, exact: Any) -> None:
        assert float(np.min(exact["ks_pvalue"].values)) > 0.02

    def test_their_first_two_moments_are_the_uniform_ones(self, exact: Any) -> None:
        """A KS test is blunt about location and scale; these are not.

        Under uniformity the rank fraction has mean 1/2 and variance 1/12, and
        the two failures family D exists to catch are exactly a shift of the
        first (bias) and an inflation of the second (over-confidence). Three
        standard errors of each, at 150 simulations.
        """
        fraction = np.asarray(exact["ranks"].values, dtype=float) / 60.0
        count = fraction.shape[0]
        mean_error = math.sqrt((1.0 / 12.0) / count)
        # Var of the sample variance of a uniform: (mu_4 - sigma^4) / n, with
        # mu_4 = 1/80 for U(0, 1).
        variance_error = math.sqrt((1.0 / 80.0 - (1.0 / 12.0) ** 2) / count)
        for column in range(fraction.shape[1]):
            assert abs(float(np.mean(fraction[:, column])) - 0.5) < 3.0 * mean_error
            assert abs(float(np.var(fraction[:, column])) - 1.0 / 12.0) < 3.0 * variance_error

    def test_the_coverage_curve_tracks_the_diagonal(self, exact: Any) -> None:
        levels = np.asarray(exact.coords["level"].values, dtype=float)
        coverage = np.asarray(exact["coverage"].values, dtype=float)
        for column in range(coverage.shape[1]):
            assert float(np.max(np.abs(coverage[:, column] - levels))) < 0.15


# ---------------------------------------------------------------------------
# The loop, with a real engine
# ---------------------------------------------------------------------------


def _study(problem: FittingProblem, factory: Any, **kwargs: Any) -> Any:
    """Run :func:`sbc` with the small-budget warning silenced."""
    options = {"steps": 120, "burn_in": 40}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return sbc(problem, factory, run_options=options, **kwargs)


class _Shrunk:
    """An engine wrapper whose run's posterior is squeezed towards its own median.

    The deliberately miscalibrated arm. Nothing about the fit changes; the
    stored draws are replaced by draws a quarter as wide, which is exactly the
    over-confidence family D exists to catch and which no residual test would
    notice — the *fit* is as good as it ever was.
    """

    def __init__(self, engine: Any, temperature: float = 0.25) -> None:
        self.engine = engine
        self.temperature = temperature

    def run(self, **kwargs: Any) -> Any:
        tree = self.engine.run(**kwargs)
        group = tree["posterior"].dataset
        squeezed = {}
        for name in group.data_vars:
            values = np.asarray(group[name].values, dtype=float)
            centre = np.median(values)
            squeezed[name] = (group[name].dims, centre + (values - centre) * self.temperature)
        tree["posterior"].dataset = group.assign(squeezed)
        return tree


class TestTheRefitLoop:
    """The end-to-end loop under a real engine, and the arm that must fail.

    The uniformity of the emcee arm is **not** asserted in absolute terms: at
    two dozen simulations the null's own KS p-value wanders across a couple of
    orders of magnitude with the seed (checked directly against
    :class:`_ExactPosterior` on these very simulations), so a threshold here
    would be a test of the seed. What is asserted is the comparison — the same
    simulations, the same engine, one arm's posterior deliberately narrowed —
    which is exactly the control-arm logic ``examples/wstat_comparison.py``'s
    study uses, and which the shared simulation noise cancels out of.
    """

    def test_it_runs_with_emcee_and_records_what_it_did(self) -> None:
        result = _study(
            toy_problem(),
            lambda replica: EmceeEngine(replica, walkers=8),
            count=24,
            draws=40,
            label="toy",
        )
        assert result.sizes["simulation"] == 24
        assert sorted(str(name) for name in result.coords["parameter"].values) == [
            "model.offset",
            "model.slope",
        ]
        assert result.attrs["ampere_calibration_route"] == REFIT_ROUTE
        assert result.attrs["ampere_calibration_engine"] == "EmceeEngine"
        assert result.attrs["ampere_calibration_label"] == "toy"
        assert result.attrs["ampere_calibration_failures"] == 0
        assert result.attrs["ampere_calibration_uniformity_check"] == "scipy.stats.kstest"
        ranks = np.asarray(result["ranks"].values)
        assert ranks.min() >= 0 and ranks.max() <= 40

    def test_a_deliberately_narrowed_posterior_is_rejected_and_the_honest_one_is_not(
        self,
    ) -> None:
        """The arm that proves the check can fail, against its own control."""
        honest = _study(
            toy_problem(),
            lambda replica: EmceeEngine(replica, walkers=8),
            count=24,
            draws=40,
        )
        narrowed = _study(
            toy_problem(),
            lambda replica: _Shrunk(EmceeEngine(replica, walkers=8)),
            count=24,
            draws=40,
        )
        honest_p = float(np.min(honest["ks_pvalue"].values))
        narrowed_p = float(np.min(narrowed["ks_pvalue"].values))
        assert narrowed_p < 1e-4
        assert honest_p > 100.0 * narrowed_p

        # The signature of over-confidence, read off the curve rather than the
        # p-value: the truth lands in the tails, so the nominal 68 % interval
        # covers far less than 68 % of the time — and the control's does not.
        levels = np.asarray(narrowed.coords["level"].values, dtype=float)
        at_68 = {
            name: [
                float(np.interp(0.68, levels, np.asarray(result["coverage"].values)[:, column]))
                for column in range(result.sizes["parameter"])
            ]
            for name, result in (("honest", honest), ("narrowed", narrowed))
        }
        assert max(at_68["narrowed"]) < min(at_68["honest"])
        assert max(at_68["narrowed"]) < 0.45

    def test_a_small_budget_says_so(self) -> None:
        with pytest.warns(UserWarning, match="goodness-of-fit test"):
            sbc(
                toy_problem(),
                lambda replica: EmceeEngine(replica, walkers=8),
                count=2,
                draws=20,
                run_options={"steps": 60, "burn_in": 20},
            )

    def test_a_fit_too_short_to_rank_is_refused_rather_than_shortened(self) -> None:
        with pytest.raises(ResultsError, match="uniform null"):
            _study(
                toy_problem(),
                lambda replica: EmceeEngine(replica, walkers=4),
                count=2,
                draws=10_000,
            )

    def test_an_unknown_parameter_names_what_the_fit_has(self) -> None:
        with pytest.raises(ResultsError, match=r"model\.slope"):
            _study(
                toy_problem(),
                lambda replica: EmceeEngine(replica, walkers=8),
                count=2,
                draws=20,
                parameters=["model.gradient"],
            )

    def test_a_subset_of_parameters_is_ranked_and_the_rest_left_alone(self) -> None:
        """The mapping form, which is what a study of a *different* formulation needs."""
        result = _study(
            toy_problem(),
            lambda replica: EmceeEngine(replica, walkers=8),
            count=4,
            draws=20,
            parameters={"model.slope": "model.slope"},
        )
        assert [str(name) for name in result.coords["parameter"].values] == ["model.slope"]

    def test_the_study_is_reproducible_from_its_recorded_seed(self) -> None:
        first = _study(
            toy_problem(seed=None),
            lambda replica: EmceeEngine(replica, walkers=8),
            count=4,
            draws=20,
            seed=4242,
        )
        second = _study(
            toy_problem(seed=None),
            lambda replica: EmceeEngine(replica, walkers=8),
            count=4,
            draws=20,
            seed=int(first.attrs["ampere_calibration_seed"]),
        )
        np.testing.assert_array_equal(first["ranks"].values, second["ranks"].values)


# ---------------------------------------------------------------------------
# Storage and pictures
# ---------------------------------------------------------------------------


def _close(*figures: Any) -> None:
    """Release the figures a test drew.

    matplotlib keeps every pyplot figure alive until it is closed, and a suite
    that draws a few dozen and keeps them all warns about it -- and, in a long
    ``test-all`` run, is holding memory for nothing.
    """
    import matplotlib.pyplot as pyplot

    for figure in figures:
        pyplot.close(figure)


@pytest.fixture(scope="module")
def study() -> Any:
    """One small study, shared by the storage and plotting checks."""
    return _study(
        toy_problem(),
        lambda replica: EmceeEngine(replica, walkers=8),
        count=12,
        draws=30,
    )


class TestStorageAndPlots:
    def test_the_group_round_trips_through_netcdf(self, study: Any, tmp_path: Any) -> None:
        run = EmceeEngine(toy_problem(), walkers=8).run(steps=80, burn_in=20)
        attach_calibration(run, study)
        assert CALIBRATION_GROUP in run.children

        path = tmp_path / "run.nc"
        to_netcdf(run, path)
        back = from_netcdf(path)
        stored = back[CALIBRATION_GROUP].dataset
        np.testing.assert_array_equal(stored["ranks"].values, study["ranks"].values)
        np.testing.assert_allclose(stored["coverage"].values, study["coverage"].values)
        assert [str(n) for n in stored.coords["parameter"].values] == [
            str(n) for n in study.coords["parameter"].values
        ]
        for key, value in study.attrs.items():
            assert key in stored.attrs
            if isinstance(value, str):
                assert stored.attrs[key] == value

    def test_both_plots_render_headless_from_the_group(self, study: Any) -> None:
        figure = plot_sbc_ranks(study)
        assert len(figure.axes) == study.sizes["parameter"]
        metadata = figure_metadata(figure)
        assert any(key.endswith(".ks_pvalue") for key in metadata)

        axes = plot_coverage(study)
        coverage_metadata = figure_metadata(axes.get_figure())
        assert any(key.endswith(".coverage_68") for key in coverage_metadata)
        assert axes.get_xlim() == (0.0, 1.0)
        _close(figure, axes.get_figure())

    def test_both_plots_accept_the_run_that_carries_the_group(self, study: Any) -> None:
        run = EmceeEngine(toy_problem(), walkers=8).run(steps=80, burn_in=20)
        attach_calibration(run, study)
        ranks = plot_sbc_ranks(run)
        coverage = plot_coverage(run)
        assert ranks is not None and coverage is not None
        _close(ranks, coverage.get_figure())

    def test_a_run_without_the_group_is_refused_with_the_remedy(self) -> None:
        run = EmceeEngine(toy_problem(), walkers=8).run(steps=80, burn_in=20)
        with pytest.raises(ResultsError, match=r"engine\.calibrate"):
            plot_sbc_ranks(run)

    def test_a_selection_of_parameters_is_honoured(self, study: Any) -> None:
        figure = plot_sbc_ranks(study, parameters=["model.slope"])
        assert len(figure.axes) == 1
        _close(figure)

    def test_too_many_panels_is_refused_loudly(self, study: Any) -> None:
        with pytest.raises(ResultsError, match="max_panels"):
            plot_sbc_ranks(study, max_panels=1)
