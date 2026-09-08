"""W2.7: ``diagnostics.md``'s post-fit families, end to end on a toy problem.

The acceptance criterion is that both post-fit families produce output on a toy
problem, that every :class:`~ampere.core.AnomalyScore` renders with its
provenance shown, and that the whiteness statistic actually *works* — which is
the row worth arguing about, because a test that only asserts "a number came
out" would pass an implementation that returned a constant.

So the toy problem comes in two versions built from the same model and the same
noise budget, differing only in whether the data carry structure the model
cannot fit:

* :func:`white_problem` — the data *are* the model plus independent Gaussian
  noise. The whiteness test must **not** reject: this is the null, and a test
  that fires here is a test that would send every user to the flexible
  likelihood whether they needed it or not.
* :func:`structured_problem` — the same, plus a smooth sinusoidal deviation
  with a correlation length several samples wide. The whiteness test must
  reject: this is exactly the "leftover structure a smooth model cannot
  explain" that ``DEVELOPMENT_PLAN.md`` §4.8 wants to replace an act of faith
  with.

Both are sampled on an **irregular** coordinate axis, because the whole point
of the separation-binned statistic is that it never assumes otherwise
(``results_schema.md`` §4.2's "no regular-grid assumption anywhere"), and the
regular case is checked separately against the classical answer it must reduce
to.

Figures are rendered on the Agg backend into memory and never written: no
binary artefact belongs in this repository, and a diagnostic plot that needs a
display to be tested is a diagnostic plot nobody will run in CI.
"""

from __future__ import annotations

import warnings
from typing import Any

import astropy.units as u
import matplotlib
import numpy as np
import pytest
import scipy.stats as st

matplotlib.use("Agg")

from ampere.core import (
    AnomalyScore,
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    PoissonFamily,
    Spectrum,
)
from ampere.core.exceptions import ResultsError
from ampere.results import (
    GP_LOCALISATION_CAVEAT,
    GP_LOCALISATION_GROUP,
    RESIDUALS_GROUP,
    DrawRecorder,
    add_residuals,
    chi_square_pvalue,
    figure_metadata,
    gp_localisation,
    gp_localisation_score,
    plot_anomaly_score,
    plot_gp_localisation,
    plot_residuals,
    residual_whiteness,
    separation_binned_autocorrelation,
)

pytest.importorskip("arviz", reason="ampere.results needs arviz")

pyplot = pytest.importorskip("matplotlib.pyplot")


# ---------------------------------------------------------------------------
# The toy problem
# ---------------------------------------------------------------------------

#: An irregular axis: a linear grid jittered by a fixed, reproducible amount,
#: so no two spacings are equal and nothing can accidentally take a regular
#: fast path.
GRID = np.sort(np.linspace(1.0, 9.0, 90) + np.random.default_rng(20260907).uniform(-0.02, 0.02, 90))

#: The noise budget, the same in both versions of the problem.
SIGMA = 0.02

#: Amplitude and length of the structure the misspecified data carry. The
#: length is several samples wide on purpose: structure narrower than the
#: sampling is indistinguishable from noise, and a diagnostic that claimed
#: otherwise would be lying.
BUMP_AMPLITUDE = 0.05
BUMP_LENGTH = 1.2


class Powerlaw(Model):
    """``norm * x ** index`` on a fixed grid — the smooth model that gets fitted."""

    def __init__(self, grid: np.ndarray) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.3)))
        self.register_parameter(Parameter("norm", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            Spectrum(ctx["grid"] * u.micron, ctx["norm"] * ctx["grid"] ** ctx["index"])
        )


def _observed(structured: bool) -> Spectrum:
    """The data: the truth, plus noise, plus (optionally) unmodelled structure."""
    rng = np.random.default_rng(20260908)
    truth = 1.0 * GRID**-1.0
    values = truth + rng.normal(0.0, SIGMA, GRID.size)
    if structured:
        values = values + BUMP_AMPLITUDE * np.sin(2.0 * np.pi * GRID / BUMP_LENGTH)
    return Spectrum(GRID * u.micron, values, uncertainty=np.full(GRID.size, SIGMA))


def white_problem(**kwargs: Any) -> FittingProblem:
    """The well-specified fit: independent noise, and the data really are white."""
    return FittingProblem(
        Powerlaw(GRID),
        DatasetCollection({"sed": Dataset(_observed(structured=False), label="sed")}),
        seed=20260908,
        **kwargs,
    )


def structured_problem(**kwargs: Any) -> FittingProblem:
    """The misspecified fit: the same smooth model against structured data."""
    return FittingProblem(
        Powerlaw(GRID),
        DatasetCollection({"sed": Dataset(_observed(structured=True), label="sed")}),
        seed=20260908,
        **kwargs,
    )


def gp_problem(**kwargs: Any) -> FittingProblem:
    """The flexible fit: the same structured data, with the GP switched on.

    The GP hyperparameters are declared as ordinary parameters on the noise
    model, which is what ``parameters.md`` §13 anticipated and what family C's
    input is (``diagnostics.md`` §4.2).
    """
    return FittingProblem(
        Powerlaw(GRID),
        DatasetCollection(
            {
                "sed": Dataset(
                    _observed(structured=True),
                    label="sed",
                    likelihood=Likelihood(
                        GaussianFamily(),
                        GaussianProcessNoise(
                            Matern32(
                                length_scale=st.loguniform(0.2, 5.0),
                                amplitude=st.loguniform(0.01, 1.0),
                            ),
                            solver=DenseGP(),
                        ),
                    ),
                )
            }
        ),
        seed=20260908,
        **kwargs,
    )


def near_truth(problem: FittingProblem, chains: int = 1, draws: int = 6) -> Any:
    """A small genuine run concentrated near the truth.

    Emcee would do the same job and take a hundred times as long: what these
    tests need is draws that are *good fits*, so that a residual is a residual
    and not a statement about a sampler's burn-in. Every draw is a real
    ``FittingProblem.evaluate``, so the emitted run is a real run — the same
    ``DrawRecorder`` path ``ampere.inference``'s drivers use.
    """
    recorder = DrawRecorder(problem, chains=chains)
    rng = np.random.default_rng(4242)
    for chain in range(chains):
        for _ in range(draws):
            values = {"model.index": -1.0 + rng.normal(0.0, 0.004), "model.norm": 1.0}
            for name in problem.parameters.free_names:
                if name.endswith("length_scale"):
                    values[name] = 0.4
                elif name.endswith("amplitude"):
                    values[name] = 0.05
            recorder.record(values, chain=chain)
    return recorder.emit(engine="fixture")


# ---------------------------------------------------------------------------
# The statistic itself
# ---------------------------------------------------------------------------


class TestSeparationBinnedAutocorrelation:
    """The estimator, checked against answers that are known without it."""

    def test_it_reduces_to_the_classical_lag_autocorrelation(self) -> None:
        # diagnostics.md §3.2's whole worry was that generalising to irregular
        # coordinates would give up the classical statistic. On an evenly spaced
        # axis with the default window it gives up nothing: bin b IS lag b + 1.
        x = np.arange(200.0)
        rng = np.random.default_rng(7)
        r = rng.normal(size=200)
        r[1:] += 0.8 * r[:-1]  # an AR(1) process: a known autocorrelation shape
        _, rho, counts = separation_binned_autocorrelation(x, r, bins=3)
        expected = [
            float(np.mean(r[:-lag] * r[lag:]) / np.mean(r**2)) * (r.size - lag) / (r.size - lag)
            for lag in (1, 2, 3)
        ]
        assert counts.tolist() == [199, 198, 197]
        assert np.allclose(rho, expected, atol=1e-12)

    def test_white_noise_has_no_binned_correlation(self) -> None:
        x = np.sort(np.random.default_rng(3).uniform(0.0, 50.0, 400))
        r = np.random.default_rng(4).normal(size=400)
        _, rho, _ = separation_binned_autocorrelation(x, r, bins=4)
        assert np.all(np.abs(rho) < 0.2)

    def test_smooth_structure_shows_as_positive_short_separation_correlation(self) -> None:
        x = np.sort(np.random.default_rng(5).uniform(0.0, 50.0, 400))
        r = np.sin(2.0 * np.pi * x / 20.0)
        _, rho, _ = separation_binned_autocorrelation(x, r, bins=4)
        assert rho[0] > 0.9

    def test_masked_samples_take_no_part(self) -> None:
        x = np.arange(100.0)
        r = np.where(np.arange(100) % 2 == 0, 1.0, -1.0)
        r[10:20] = np.nan
        _, rho, counts = separation_binned_autocorrelation(x, r, bins=1)
        assert np.isclose(rho[0], -1.0)
        assert int(counts.sum()) < 99

    def test_it_refuses_a_length_mismatch(self) -> None:
        with pytest.raises(ResultsError, match="one residual per coordinate"):
            separation_binned_autocorrelation(np.arange(5.0), np.zeros(4))

    def test_it_refuses_too_few_samples(self) -> None:
        with pytest.raises(ResultsError, match="at least three"):
            separation_binned_autocorrelation(np.arange(2.0), np.zeros(2))


# ---------------------------------------------------------------------------
# Family B, end to end
# ---------------------------------------------------------------------------


class TestSignedResiduals:
    def test_they_are_not_stored_by_default_and_are_derived_on_demand(self) -> None:
        problem = white_problem()
        tree = near_truth(problem)
        assert RESIDUALS_GROUP not in tree.children
        add_residuals(tree, problem)
        assert RESIDUALS_GROUP in tree.children
        residuals = tree[RESIDUALS_GROUP].dataset
        assert residuals["sed"].dims == ("chain", "draw", "sed_spectral_axis")
        # results.md §7's derivation, checked against itself: a standardised
        # residual of a good fit to sigma-level noise is order 1.
        values = np.asarray(residuals["sed"].values)
        assert np.all(np.isfinite(values))
        assert 0.5 < float(np.std(values)) < 2.0

    def test_the_sign_survives(self) -> None:
        # diagnostics.md §3.3's precise gap: |residual| is recoverable from the
        # log-likelihood and the sign is not, which is why this group exists.
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        values = np.asarray(tree[RESIDUALS_GROUP]["sed"].values)
        assert np.any(values > 0.0) and np.any(values < 0.0)

    def test_unstandardised_residuals_are_in_the_data_units(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem), problem, standardised=False)
        values = np.asarray(tree[RESIDUALS_GROUP]["sed"].values)
        assert float(np.std(values)) < 0.5 * SIGMA * 3.0
        assert int(tree[RESIDUALS_GROUP].attrs["ampere_standardised"]) == 0

    def test_masked_samples_are_nan_and_not_dropped(self) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[10:15] = True
        observed = Spectrum(
            GRID * u.micron,
            _observed(structured=False).values,
            uncertainty=np.full(GRID.size, SIGMA),
            mask=mask,
        )
        problem = FittingProblem(
            Powerlaw(GRID),
            DatasetCollection({"sed": Dataset(observed, label="sed")}),
            seed=20260908,
        )
        tree = add_residuals(near_truth(problem), problem)
        values = np.asarray(tree[RESIDUALS_GROUP]["sed"].values)
        assert values.shape[-1] == GRID.size, "the coordinate axis stays the container's own"
        assert np.all(np.isnan(values[..., 10:15]))
        assert np.all(np.isfinite(values[..., :10]))

    def test_thinning_keeps_the_draw_correspondence(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem, draws=6), problem, thin=3)
        assert tree[RESIDUALS_GROUP]["draw"].values.tolist() == [0, 3]

    def test_it_refuses_a_different_problem(self) -> None:
        tree = near_truth(white_problem())
        with pytest.raises(ResultsError, match="different problem"):
            add_residuals(tree, structured_problem())

    def test_it_refuses_to_standardise_without_uncertainties(self) -> None:
        # A Poisson dataset models its own dispersion, so it legitimately has
        # no per-sample sigma — and a "standardised" residual then has no
        # denominator. Better to say so than to invent one.
        counts = Spectrum(GRID * u.micron, np.round(60.0 * GRID**-1.0 + 3.0))

        class Rate(Model):
            def __init__(self, grid: np.ndarray) -> None:
                self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
                self.register_parameter(Parameter("norm", st.loguniform(10.0, 200.0)))

            def evaluate(self, **values: Any) -> ModelResult:
                ctx = self.context(values)
                return ModelResult(
                    Spectrum(ctx["grid"] * u.micron, ctx["norm"] * ctx["grid"] ** -1.0 + 3.0)
                )

        problem = FittingProblem(
            Rate(GRID),
            DatasetCollection(
                {
                    "sed": Dataset(
                        counts,
                        label="sed",
                        likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
                    )
                }
            ),
            seed=20260908,
        )
        recorder = DrawRecorder(problem, chains=1)
        recorder.record({"model.norm": 60.0})
        tree = recorder.emit(engine="fixture")
        with pytest.raises(ResultsError, match="no uncertainties"):
            add_residuals(tree, problem)
        # ... and says what to do instead, which works.
        add_residuals(tree, problem, standardised=False)
        assert RESIDUALS_GROUP in tree.children


class TestResidualWhiteness:
    """The acceptance row: the test fires on structure and not on noise."""

    def test_white_residuals_are_not_rejected(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem), problem)
        test = residual_whiteness(tree, n_permutations=199)
        assert test.dataset == "sed"
        assert test.p_value > 0.1, f"the null was rejected at p={test.p_value}"
        assert test.n_pairs > 0 and test.n_draws == 6

    def test_structured_residuals_are_rejected(self) -> None:
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        test = residual_whiteness(tree, n_permutations=199)
        assert test.p_value <= 0.01, f"structure was not detected (p={test.p_value})"
        # And it is detected as *positive* correlation at short separation,
        # which is the shape a Matern kernel is built to absorb.
        assert test.autocorrelation[0] > 0.5

    def test_the_statistic_is_larger_where_the_structure_is(self) -> None:
        white = residual_whiteness(
            add_residuals(near_truth(white_problem()), white_problem()), n_permutations=99
        )
        structured = residual_whiteness(
            add_residuals(near_truth(structured_problem()), structured_problem()),
            n_permutations=99,
        )
        assert structured.statistic > 10.0 * white.statistic

    def test_it_is_reproducible_from_the_stored_run_alone(self) -> None:
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        first = residual_whiteness(tree, n_permutations=99)
        second = residual_whiteness(tree, n_permutations=99)
        assert first.p_value == second.p_value
        assert np.array_equal(first.statistics, second.statistics)

    def test_a_different_seed_gives_a_different_null_draw(self) -> None:
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        one = residual_whiteness(tree, n_permutations=39, seed=1)
        two = residual_whiteness(tree, n_permutations=39, seed=2)
        assert one.statistic == pytest.approx(two.statistic), "the statistic is not random"
        assert one.seed != two.seed

    def test_the_p_value_cannot_beat_its_own_resolution(self) -> None:
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        test = residual_whiteness(tree, n_permutations=19)
        assert test.p_value >= test.resolution
        assert test.resolution == pytest.approx(1.0 / 20.0)

    def test_it_names_add_residuals_when_the_group_is_missing(self) -> None:
        with pytest.raises(ResultsError, match="add_residuals"):
            residual_whiteness(near_truth(white_problem()))

    def test_it_refuses_at_least_one_permutation(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem), problem)
        with pytest.raises(ResultsError, match="at least one"):
            residual_whiteness(tree, n_permutations=0)


class TestChiSquarePosteriorPredictive:
    """The cheap half of family B: a Bayesian p-value from the stored numbers."""

    def test_a_good_fit_is_not_rejected(self) -> None:
        check = chi_square_pvalue(near_truth(white_problem()))
        assert 0.01 < check.p_value < 0.99
        assert check.degrees_of_freedom == GRID.size

    def test_a_misspecified_fit_is_rejected(self) -> None:
        check = chi_square_pvalue(near_truth(structured_problem()))
        assert check.p_value < 1e-3

    def test_it_agrees_with_the_residual_group_it_did_not_use(self) -> None:
        # The point of the check: the log-likelihood route and the residual
        # route are the same number, so "cheap" is not "approximate".
        problem = structured_problem()
        tree = near_truth(problem, draws=3)
        check = chi_square_pvalue(tree)
        add_residuals(tree, problem)
        residuals = np.asarray(tree[RESIDUALS_GROUP]["sed"].values).reshape(-1, GRID.size)
        from_residuals = np.sort(np.nansum(residuals**2, axis=1))
        assert np.allclose(np.sort(check.statistics), from_residuals, rtol=1e-8)

    def test_it_refuses_a_gp_fit_by_name(self) -> None:
        with pytest.raises(ResultsError, match="does not factorise"):
            chi_square_pvalue(near_truth(gp_problem()))

    def test_it_refuses_a_fitted_error_scale_by_name(self) -> None:
        problem = FittingProblem(
            Powerlaw(GRID),
            DatasetCollection(
                {
                    "sed": Dataset(
                        _observed(structured=False),
                        label="sed",
                        likelihood=Likelihood(
                            GaussianFamily(), IndependentNoise(scale=st.loguniform(0.5, 2.0))
                        ),
                    )
                }
            ),
            seed=20260908,
        )
        recorder = DrawRecorder(problem, chains=1)
        recorder.record({"model.index": -1.0, "model.norm": 1.0, "sed.likelihood.scale": 1.0})
        with pytest.raises(ResultsError, match="sigma moves from draw to draw"):
            chi_square_pvalue(recorder.emit(engine="fixture"))


# ---------------------------------------------------------------------------
# Family C, end to end
# ---------------------------------------------------------------------------


class TestGpLocalisation:
    def test_the_conditioned_mean_localises_the_structure(self) -> None:
        problem = gp_problem()
        tree = gp_localisation(near_truth(problem, draws=2), problem)
        assert GP_LOCALISATION_GROUP in tree.children
        group = tree[GP_LOCALISATION_GROUP].dataset
        mean = np.asarray(group["sed_mean"].values)[0, 0]
        assert mean.shape == (GRID.size,)
        assert np.all(np.asarray(group["sed_variance"].values) >= 0.0)
        # It is signed, and it tracks the deviation the model could not fit.
        truth = BUMP_AMPLITUDE * np.sin(2.0 * np.pi * GRID / BUMP_LENGTH)
        assert float(np.corrcoef(mean, truth)[0, 1]) > 0.8

    def test_it_covers_the_full_axis_masked_samples_included(self) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[40:50] = True
        observed = Spectrum(
            GRID * u.micron,
            _observed(structured=True).values,
            uncertainty=np.full(GRID.size, SIGMA),
            mask=mask,
        )
        problem = FittingProblem(
            Powerlaw(GRID),
            DatasetCollection(
                {
                    "sed": Dataset(
                        observed,
                        label="sed",
                        likelihood=Likelihood(
                            GaussianFamily(),
                            GaussianProcessNoise(Matern32(0.4, 0.05), solver=DenseGP()),
                        ),
                    )
                }
            ),
            seed=20260908,
        )
        tree = gp_localisation(near_truth(problem, draws=1), problem)
        mean = np.asarray(tree[GP_LOCALISATION_GROUP]["sed_mean"].values)[0, 0]
        assert mean.shape == (GRID.size,), "the excluded region is on the axis"
        assert np.all(np.isfinite(mean)), "'what would the GP have said here?' is answered"

    def test_an_explicit_grid_is_honoured(self) -> None:
        problem = gp_problem()
        grid = np.linspace(1.0, 9.0, 31)
        tree = gp_localisation(near_truth(problem, draws=1), problem, at=grid)
        group = tree[GP_LOCALISATION_GROUP].dataset
        assert np.asarray(group["sed_mean"].values).shape[-1] == 31
        assert group.attrs["ampere_evaluation_grid"] == "explicit"

    def test_it_refuses_a_problem_with_no_gp(self) -> None:
        problem = white_problem()
        with pytest.raises(ResultsError, match="no dataset in this problem"):
            gp_localisation(near_truth(problem), problem)

    def test_a_named_non_gp_dataset_propagates_conditional_refusal(self) -> None:
        from ampere.core.exceptions import LikelihoodError

        problem = white_problem()
        with pytest.raises(LikelihoodError, match="needs a GaussianProcessNoise"):
            gp_localisation(near_truth(problem), problem, datasets=["sed"])


class TestAnomalyScore:
    def test_the_score_carries_its_provenance_and_the_caveat(self) -> None:
        problem = gp_problem()
        tree = gp_localisation(near_truth(problem, draws=2), problem)
        score = gp_localisation_score(tree)
        assert isinstance(score, AnomalyScore)
        assert score.provenance == "gp_localisation_postfit"
        # diagnostics.md §4.3: the caveat must reach a user who never looks at
        # a figure. It is in the container, not only in the plot.
        assert GP_LOCALISATION_CAVEAT in score.interpretation_notes
        assert score.n_samples == GRID.size
        assert np.all(score.values >= 0.0), "higher = worse, on a documented range"

    def test_the_score_is_largest_where_the_deviation_is(self) -> None:
        problem = gp_problem()
        score = gp_localisation_score(gp_localisation(near_truth(problem, draws=2), problem))
        truth = np.abs(BUMP_AMPLITUDE * np.sin(2.0 * np.pi * GRID / BUMP_LENGTH))
        assert float(np.corrcoef(score.values, truth)[0, 1]) > 0.6

    def test_it_names_gp_localisation_when_the_group_is_missing(self) -> None:
        with pytest.raises(ResultsError, match="gp_localisation"):
            gp_localisation_score(near_truth(gp_problem()))


# ---------------------------------------------------------------------------
# The renderers
# ---------------------------------------------------------------------------


class TestPlotResiduals:
    def test_it_draws_both_panels_and_reports_the_test(self) -> None:
        problem = structured_problem()
        tree = add_residuals(near_truth(problem), problem)
        figure = plot_residuals(tree, n_permutations=39)
        try:
            assert len(figure.axes) == 2
            metadata = figure_metadata(figure)
            assert float(metadata["sed.whiteness_p_value"]) <= 0.05
            assert float(metadata["sed.whiteness_statistic"]) > 0.0
        finally:
            pyplot.close(figure)

    def test_whiteness_false_draws_only_the_residual_panel(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem), problem)
        figure = plot_residuals(tree, whiteness=False)
        try:
            assert len(figure.axes) == 1
        finally:
            pyplot.close(figure)

    def test_it_warns_on_a_gp_fit_and_points_at_family_c(self) -> None:
        # diagnostics.md §3.1: a GP-augmented fit's residuals are whitened by
        # construction, so testing them for whiteness is close to circular.
        problem = gp_problem()
        tree = add_residuals(near_truth(problem, draws=2), problem)
        with pytest.warns(UserWarning, match="plot_gp_localisation"):
            figure = plot_residuals(tree, whiteness=False)
        pyplot.close(figure)

    def test_a_standard_fit_draws_without_a_warning(self) -> None:
        problem = white_problem()
        tree = add_residuals(near_truth(problem), problem)
        with warnings.catch_warnings():
            warnings.simplefilter("error", UserWarning)
            figure = plot_residuals(tree, whiteness=False)
        pyplot.close(figure)

    def test_it_names_add_residuals_rather_than_a_missing_group(self) -> None:
        with pytest.raises(ResultsError, match="add_residuals"):
            plot_residuals(near_truth(white_problem()))


class TestPlotGpLocalisation:
    def test_it_draws_and_carries_the_caveat_in_metadata(self) -> None:
        problem = gp_problem()
        tree = gp_localisation(near_truth(problem, draws=2), problem)
        figure = plot_gp_localisation(tree)
        try:
            assert figure_metadata(figure)["gp_localisation_caveat"] == GP_LOCALISATION_CAVEAT
            drawn = " ".join(text.get_text() for text in figure.texts)
            assert "does not say why" in " ".join(drawn.split())
        finally:
            pyplot.close(figure)

    def test_show_caveat_false_suppresses_only_the_annotation(self) -> None:
        # results.md §8 is explicit: setting it False "must still leave the
        # caveat on the returned figure's metadata... and only suppresses the
        # drawn annotation".
        problem = gp_problem()
        tree = gp_localisation(near_truth(problem, draws=2), problem)
        figure = plot_gp_localisation(tree, show_caveat=False)
        try:
            assert figure_metadata(figure)["gp_localisation_caveat"] == GP_LOCALISATION_CAVEAT
            drawn = " ".join(text.get_text() for text in figure.texts)
            assert "does not say why" not in drawn
        finally:
            pyplot.close(figure)

    def test_it_refuses_an_impossible_band(self) -> None:
        with pytest.raises(ResultsError, match="credible-interval mass"):
            plot_gp_localisation(near_truth(gp_problem()), band=1.5)

    def test_it_names_gp_localisation_rather_than_a_missing_group(self) -> None:
        with pytest.raises(ResultsError, match=r"ampere\.results\.gp_localisation"):
            plot_gp_localisation(near_truth(gp_problem()))


class TestPlotAnomalyScore:
    @staticmethod
    def _score(provenance: str) -> AnomalyScore:
        return AnomalyScore(
            coordinates=GRID,
            values=np.abs(np.sin(GRID)),
            provenance=provenance,
            interpretation_notes=f"a {provenance} score, for testing",
        )

    def test_one_renderer_takes_either_provenance(self) -> None:
        for provenance in ("rhmf_prefit", "gp_localisation_postfit"):
            axes = plot_anomaly_score(self._score(provenance))
            try:
                assert axes.get_legend() is not None
                assert provenance in axes.get_legend_handles_labels()[1]
            finally:
                pyplot.close(axes.get_figure())

    def test_provenance_reaches_a_programmatic_consumer(self) -> None:
        score = self._score("rhmf_prefit")
        axes = plot_anomaly_score(score, show_provenance=False)
        try:
            metadata = figure_metadata(axes.get_figure())
            assert metadata["anomaly_score.rhmf_prefit"] == score.interpretation_notes
        finally:
            pyplot.close(axes.get_figure())

    def test_show_provenance_false_is_honoured_for_one_score(self) -> None:
        axes = plot_anomaly_score(self._score("rhmf_prefit"), show_provenance=False)
        try:
            assert axes.get_legend() is None
        finally:
            pyplot.close(axes.get_figure())

    def test_two_provenances_on_one_axes_force_provenance_back_on(self) -> None:
        # results.md §8: "show_provenance=False may not suppress them when two
        # scores of different provenance are drawn together". Together means on
        # one axes, and only the axes knows.
        axes = plot_anomaly_score(self._score("rhmf_prefit"), show_provenance=False)
        try:
            assert axes.get_legend() is None
            with pytest.warns(UserWarning, match="show_provenance=False was ignored"):
                plot_anomaly_score(
                    self._score("gp_localisation_postfit"), ax=axes, show_provenance=False
                )
            labels = axes.get_legend_handles_labels()[1]
            assert sorted(labels) == ["gp_localisation_postfit", "rhmf_prefit"]
        finally:
            pyplot.close(axes.get_figure())

    def test_two_scores_of_the_same_provenance_are_not_forced(self) -> None:
        axes = plot_anomaly_score(self._score("rhmf_prefit"), show_provenance=False)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", UserWarning)
                plot_anomaly_score(self._score("rhmf_prefit"), ax=axes, show_provenance=False)
            assert axes.get_legend() is None
        finally:
            pyplot.close(axes.get_figure())

    def test_the_end_to_end_family_c_score_renders_with_its_provenance(self) -> None:
        problem = gp_problem()
        score = gp_localisation_score(gp_localisation(near_truth(problem, draws=2), problem))
        axes = plot_anomaly_score(score)
        try:
            assert "gp_localisation_postfit" in axes.get_legend_handles_labels()[1]
            assert GP_LOCALISATION_CAVEAT in figure_metadata(axes.get_figure()).get(
                "anomaly_score.gp_localisation_postfit", ""
            )
        finally:
            pyplot.close(axes.get_figure())

    def test_it_refuses_something_that_is_not_a_score(self) -> None:
        with pytest.raises(ResultsError, match="AnomalyScoreLike"):
            plot_anomaly_score(object())

    def test_it_refuses_a_multi_axis_score(self) -> None:
        score = AnomalyScore(
            coordinates=np.column_stack([GRID, GRID]),
            values=np.abs(np.sin(GRID)),
            provenance="rhmf_prefit",
            interpretation_notes="two coordinates",
        )
        with pytest.raises(ResultsError, match="1-D deficiency map"):
            plot_anomaly_score(score)


class TestImportPolicy:
    def test_importing_ampere_results_does_not_import_matplotlib_or_arviz(self) -> None:
        import subprocess
        import sys

        code = (
            "import sys, ampere.results; "
            "print(any(n == 'arviz' or n.startswith('arviz.') for n in sys.modules), "
            "any(n == 'matplotlib.pyplot' for n in sys.modules))"
        )
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        )
        assert result.stdout.strip() == "False False"
