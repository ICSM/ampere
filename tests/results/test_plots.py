"""W2.8: the general-purpose plots and the per-observation group, end to end.

``results.md`` §8's acceptance criterion for this item is that "each plot
renders from a stored run with merged names as labels", and §6's that "a run
with the pointwise group survives netCDF and ``arviz.loo`` consumes it". Both
are asserted here against a real emitted run rather than a hand-built tree,
because a plotting function that took the emitted run and nothing else is the
whole of §8 and a fixture that skipped emission would not be testing it.

The toy problem is deliberately the smallest one that exercises what the plots
have to get right: two free parameters (so a corner plot has a grid), one
array-valued block behind a separate fixture (so the refusal threshold has
something to refuse), a prior-rejected draw (so a trace has a gap to draw), and
one dataset with uncertainties (so a predictive check has a scale).

Figures are rendered on the Agg backend into memory and never written: no
binary artefact belongs in this repository (``AGENTS.md`` ground rule 7).
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Any

import astropy.units as u
import matplotlib
import numpy as np
import pytest
import scipy.stats as st

matplotlib.use("Agg")

from ampere.core import (
    AxisSpec,
    ClosurePhases,
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    FunctionSamples,
    GaussianFamily,
    GaussianProcessNoise,
    Layout,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    QuasisepGP,
    Spectrum,
    VisibilitySet,
    VonMisesFamily,
)
from ampere.core.exceptions import ResultsError
from ampere.results import (
    CONDITIONAL_LOO_DECOMPOSITION,
    FACTORISED_DECOMPOSITION,
    LOG_LIKELIHOOD_GROUP,
    POINTWISE_LOG_LIKELIHOOD_GROUP,
    POSTERIOR_PREDICTIVE_GROUP,
    DrawRecorder,
    ResultsWarning,
    add_pointwise_log_likelihood,
    add_posterior_predictive,
    add_residuals,
    figure_metadata,
    from_netcdf,
    gp_localisation,
    gp_localisation_score,
    plot_anomaly_score,
    plot_corner,
    plot_gp_localisation,
    plot_posterior_predictive,
    plot_residuals,
    plot_trace,
    pointwise_as_log_likelihood,
    to_netcdf,
)

pytest.importorskip("arviz", reason="ampere.results needs arviz")

pyplot = pytest.importorskip("matplotlib.pyplot")

GRID = np.linspace(1.0, 9.0, 24)
SIGMA = 0.05
TRUTH = {"model.index": -1.0, "model.norm": 1.0}


class Powerlaw(Model):
    """``norm * x ** index`` on a fixed grid."""

    def __init__(self, grid: np.ndarray) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.3)))
        self.register_parameter(Parameter("norm", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            Spectrum(
                ctx["grid"] * u.micron,
                ctx["norm"] * ctx["grid"] ** ctx["index"] * u.Jy,
            )
        )


def observed(mask: np.ndarray | None = None) -> Spectrum:
    rng = np.random.default_rng(20260908)
    return Spectrum(
        GRID * u.micron,
        (1.0 * GRID**-1.0 + rng.normal(0.0, SIGMA, GRID.size)) * u.Jy,
        uncertainty=np.full(GRID.size, SIGMA) * u.Jy,
        mask=mask,
    )


def toy(likelihood: Likelihood | None = None, mask: np.ndarray | None = None) -> FittingProblem:
    return FittingProblem(
        Powerlaw(GRID),
        DatasetCollection({"sed": Dataset(observed(mask), label="sed", likelihood=likelihood)}),
        seed=20260908,
    )


def gp_toy(solver: Any) -> FittingProblem:
    return toy(
        Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(
                Matern32(length_scale=st.loguniform(0.2, 5.0), amplitude=st.loguniform(0.01, 1.0)),
                solver=solver,
            ),
        )
    )


class PointSource(Model):
    """A flat complex visibility — the smallest complex-valued channel there is."""

    def __init__(self, u_axis: np.ndarray, v_axis: np.ndarray, wave: np.ndarray) -> None:
        self.register_buffer("u", np.asarray(u_axis, dtype=float))
        self.register_buffer("v", np.asarray(v_axis, dtype=float))
        self.register_buffer("wave", np.asarray(wave, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("flux", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            VisibilitySet(
                ctx["u"],
                ctx["v"],
                ctx["wave"] * u.micron,
                ctx["flux"] * np.ones(ctx["u"].size, dtype=complex),
            )
        )


def visibility_problem() -> FittingProblem:
    axes = (np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([2.2, 2.2]))
    return FittingProblem(
        PointSource(*axes),
        DatasetCollection(
            {
                "vis": Dataset(
                    VisibilitySet(
                        axes[0],
                        axes[1],
                        axes[2] * u.micron,
                        np.array([1 + 0j, 1 + 0j]),
                        uncertainty=np.array([0.1, 0.1]),
                    ),
                    label="vis",
                    likelihood=Likelihood(ComplexGaussianFamily()),
                )
            }
        ),
        seed=1,
    )


def gp_visibility_problem() -> FittingProblem:
    """:func:`visibility_problem`, GP-fitted over ``(u, v)`` — W5.3's family C row."""
    axes = (np.array([1.0, 2.0, -3.0]), np.array([3.0, -4.0, 1.5]), np.array([2.2, 2.2, 2.2]))
    noise = GaussianProcessNoise(Matern32(0.3, 1.0, axes=("u", "v")), solver=DenseGP())
    return FittingProblem(
        PointSource(*axes),
        DatasetCollection(
            {
                "vis": Dataset(
                    VisibilitySet(
                        axes[0],
                        axes[1],
                        axes[2] * u.micron,
                        np.array([1 + 0j, 1 + 0j, 1 + 0j]),
                        uncertainty=np.array([0.1, 0.1, 0.1]),
                    ),
                    label="vis",
                    likelihood=Likelihood(ComplexGaussianFamily(), noise),
                )
            }
        ),
        seed=1,
    )


class ClosureSource(Model):
    """A flat closure phase — the smallest :class:`~ampere.core.ClosurePhases` channel."""

    def __init__(
        self,
        u1: np.ndarray,
        v1: np.ndarray,
        u2: np.ndarray,
        v2: np.ndarray,
        wave: np.ndarray,
    ) -> None:
        self.register_buffer("u1", np.asarray(u1, dtype=float))
        self.register_buffer("v1", np.asarray(v1, dtype=float))
        self.register_buffer("u2", np.asarray(u2, dtype=float))
        self.register_buffer("v2", np.asarray(v2, dtype=float))
        self.register_buffer("wave", np.asarray(wave, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("phase", st.uniform(-1.0, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            ClosurePhases(
                ctx["u1"],
                ctx["v1"],
                ctx["u2"],
                ctx["v2"],
                ctx["wave"] * u.micron,
                ctx["phase"] * np.ones(ctx["u1"].size) * u.rad,
            )
        )


def closure_problem() -> FittingProblem:
    axes = (
        np.array([1.0, 2.0]),
        np.array([3.0, -4.0]),
        np.array([-1.0, 0.5]),
        np.array([2.0, 1.5]),
        np.array([2.2, 2.2]),
    )
    return FittingProblem(
        ClosureSource(*axes),
        DatasetCollection(
            {
                "t3": Dataset(
                    ClosurePhases(
                        *axes[:4],
                        axes[4] * u.micron,
                        np.array([0.1, -0.2]) * u.rad,
                        uncertainty=np.array([0.05, 0.05]) * u.rad,
                    ),
                    label="t3",
                    likelihood=Likelihood(VonMisesFamily()),
                )
            }
        ),
        seed=2,
    )


class TwoAxisPoint(FunctionSamples):
    """An out-of-tree point kind with two axes and no default coordinate.

    ``tests/core/thirdparty_polarimeter.py`` is the pattern for declaring a
    kind outside ``ampere`` itself; this one is deliberately local to this
    test module rather than a shared fixture, since its only job is to prove
    :func:`~ampere.results._plotting.coordinate_of` refuses a multi-axis kind
    that names no default and is given no ``coordinate=`` — DEVELOPMENT_PLAN.md
    §4.3's "a kind is class attributes" applies exactly the same to a kind
    that declares nothing beyond ``AXES``/``LAYOUT``/``ALLOW_COMPLEX``.
    """

    AXES = (
        AxisSpec("p", physical_types=("dimensionless",)),
        AxisSpec("q", physical_types=("dimensionless",)),
    )
    LAYOUT = Layout.POINTS
    ALLOW_COMPLEX = False

    def __init__(self, p: np.ndarray, q: np.ndarray, values: np.ndarray, **kwargs: Any) -> None:
        super().__init__({"p": p, "q": q}, values, **kwargs)


class PairSource(Model):
    """A flat level on a ``(p, q)`` point set — the smallest two-axis channel."""

    def __init__(self, p: np.ndarray, q: np.ndarray) -> None:
        self.register_buffer("p", np.asarray(p, dtype=float))
        self.register_buffer("q", np.asarray(q, dtype=float))
        self.register_parameter(Parameter("level", st.uniform(-1.0, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(TwoAxisPoint(ctx["p"], ctx["q"], ctx["level"] * np.ones(ctx["p"].size)))


def pair_problem() -> FittingProblem:
    p, q = np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])
    return FittingProblem(
        PairSource(p, q),
        DatasetCollection(
            {
                "pair": Dataset(
                    TwoAxisPoint(p, q, np.array([0.1, 0.2, 0.3]), uncertainty=np.full(3, 0.05)),
                    label="pair",
                    likelihood=Likelihood(GaussianFamily()),
                )
            }
        ),
        seed=3,
    )


def _tiny_run(problem: FittingProblem, *, values: dict[str, float], draws: int = 3) -> Any:
    """The smallest genuine run for a fixture whose only free names are *values*'."""
    recorder = DrawRecorder(problem, chains=1)
    for _ in range(draws):
        recorder.record(dict(values))
    return recorder.emit(engine="fixture")


def run(problem: FittingProblem, *, chains: int = 2, draws: int = 30, reject: bool = False) -> Any:
    """A small genuine run near the truth, through the driver-facing recorder."""
    recorder = DrawRecorder(problem, chains=chains)
    rng = np.random.default_rng(4242)
    for chain in range(chains):
        for draw in range(draws):
            values = {
                "model.index": -1.0 + rng.normal(0.0, 0.02),
                "model.norm": 1.0 + rng.normal(0.0, 0.02),
            }
            for name in problem.parameters.free_names:
                if name.endswith("length_scale"):
                    values[name] = 0.4
                elif name.endswith("amplitude"):
                    values[name] = 0.05
            if reject and chain == 0 and draw == 1:
                # Outside the loguniform support: -inf prior, NaN likelihood,
                # no failure (results.md §5). This is the draw a trace must
                # show as a gap.
                values["model.norm"] = 1e9
            recorder.record(values, chain=chain)
    return recorder.emit(engine="fixture")


def _with_extra_scalars(tree: Any, count: int, *, prefix: str = "extra", seed: int = 7) -> Any:
    """*tree*, with *count* extra scalar posterior variables assigned.

    W3.10's paging tests want more merged names than any real toy problem
    declares; assigning bare scalar variables straight onto the stored
    posterior — as the existing array-block tests already do for one
    variable — is the fixture-free way to get there.
    """
    posterior = tree["posterior"].dataset
    rng = np.random.default_rng(seed)
    shape = (posterior.sizes["chain"], posterior.sizes["draw"])
    extra = {f"{prefix}.p{i}": (("chain", "draw"), rng.normal(size=shape)) for i in range(count)}
    tree["posterior"] = posterior.assign(extra)
    return tree


# ---------------------------------------------------------------------------
# plot_corner
# ---------------------------------------------------------------------------


class TestCorner:
    def test_it_renders_from_a_stored_run_with_merged_names_as_labels(self) -> None:
        figure = plot_corner(run(toy()))
        assert figure_metadata(figure)["corner.variables"] == "model.index, model.norm"

    def test_a_subset_is_selected_by_merged_name(self) -> None:
        figure = plot_corner(run(toy()), var_names=["model.index"])
        assert figure_metadata(figure)["corner.variables"] == "model.index"

    def test_an_unknown_name_is_refused_with_the_list_beside_it(self) -> None:
        with pytest.raises(ResultsError, match="merged parameter name"):
            plot_corner(run(toy()), var_names=["index"])

    def test_prior_rejected_draws_are_excluded(self) -> None:
        # They are stored (results.md §5) and they are not posterior samples.
        full = plot_corner(run(toy(), reject=False))
        rejected = plot_corner(run(toy(), reject=True))
        assert int(figure_metadata(full)["corner.draws"]) == 60
        assert int(figure_metadata(rejected)["corner.draws"]) == 59

    def test_labels_must_be_one_per_column(self) -> None:
        with pytest.raises(ResultsError, match="per column, not per variable"):
            plot_corner(run(toy()), labels=["only one"])

    def test_truths_may_be_a_mapping_of_merged_names(self) -> None:
        # A user has theta as a mapping — it is what simulate() took and what
        # Simulation.parameters gives back — not as a positional list.
        assert plot_corner(run(toy()), truths=dict(TRUTH)) is not None

    def test_a_truth_for_a_column_that_is_not_drawn_is_refused(self) -> None:
        with pytest.raises(ResultsError, match="flatters the fit"):
            plot_corner(run(toy()), var_names=["model.index"], truths=dict(TRUTH))

    def test_paginate_false_refuses_an_oversized_block_exactly_as_before_w3_10(self) -> None:
        # results.md §8 pre-W3.10: "a corner plot of a 10^5-element latent
        # block must be refused loudly rather than attempted". W3.10 turned
        # the *default* behaviour into paging (see TestCornerPaging below),
        # but paginate=False must restore this refusal text unchanged, for a
        # caller who needs a single figure or a hard failure.
        tree = run(toy())
        block = np.zeros((2, 30, 500))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"latent.z": (("chain", "draw", "latent.z_dim_0"), block)}
        )
        with pytest.raises(ResultsError, match="refused loudly rather than attempted"):
            plot_corner(tree, var_names=["latent.z"], paginate=False)
        with pytest.raises(ResultsError, match="max_variables"):
            plot_corner(tree, var_names=["latent.z"], paginate=False)

    def test_the_threshold_is_an_override_and_not_a_wall(self) -> None:
        tree = run(toy())
        block = np.random.default_rng(1).normal(size=(2, 30, 3))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"block": (("chain", "draw", "block_dim_0"), block)}
        )
        figure = plot_corner(tree, var_names=["block"], max_variables=3)
        assert figure_metadata(figure)["corner.variables"] == "block[0], block[1], block[2]"

    def test_it_refuses_a_tree_that_is_not_a_run(self) -> None:
        with pytest.raises(ResultsError, match="not an emitted run"):
            plot_corner(None)


class TestCornerPaging:
    """W3.10: paging above MAX_CORNER_VARIABLES, with a loud ResultsWarning."""

    def test_25_scalar_parameters_give_two_pages_and_one_warning(self) -> None:
        # results.md §8's accept criterion, literally.
        tree = _with_extra_scalars(run(toy()), 23)  # 2 (toy's own) + 23 = 25
        with pytest.warns(ResultsWarning, match="25 variables"):
            figures = plot_corner(tree)
        assert isinstance(figures, list)
        assert [
            len(figure_metadata(figure)["corner.variables"].split(", ")) for figure in figures
        ] == [20, 5]
        assert [figure_metadata(figure)["page"] for figure in figures] == ["1 of 2", "2 of 2"]

    def test_a_200_element_plate_gives_ten_pages(self) -> None:
        # results.md §8's other accept criterion, literally.
        tree = run(toy())
        block = np.random.default_rng(3).normal(size=(2, 30, 200))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"latent.z": (("chain", "draw", "latent.z_dim_0"), block)}
        )
        with pytest.warns(ResultsWarning, match="10 pages"):
            figures = plot_corner(tree, var_names=["latent.z"])
        assert len(figures) == 10
        assert all(
            len(figure_metadata(figure)["corner.variables"].split(", ")) == 20 for figure in figures
        )
        assert [figure_metadata(figure)["page"] for figure in figures] == [
            f"{index} of 10" for index in range(1, 11)
        ]

    def test_a_call_within_the_cap_still_returns_one_figure(self) -> None:
        # The return-type contract: paginate defaults to True, and that must
        # not change the return of a call that never pages.
        figure = plot_corner(run(toy()))
        assert not isinstance(figure, list)
        assert "page" not in figure_metadata(figure)

    def test_an_array_block_moves_whole_to_a_fresh_page_rather_than_splitting(self) -> None:
        # results.md §8: "array blocks kept whole where they fit". Fifteen
        # scalars leave five slots on page one; a six-element block does not
        # fit in those five, but does fit a fresh page of its own, so the
        # whole block moves there rather than splitting five-and-one.
        tree = _with_extra_scalars(run(toy()), 13)  # 2 + 13 = 15 scalars
        block = np.random.default_rng(4).normal(size=(2, 30, 6))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"block": (("chain", "draw", "block_dim_0"), block)}
        )
        with pytest.warns(ResultsWarning):
            figures = plot_corner(tree)
        assert len(figures) == 2
        assert len(figure_metadata(figures[0])["corner.variables"].split(", ")) == 15
        assert figure_metadata(figures[1])["corner.variables"] == (
            "block[0], block[1], block[2], block[3], block[4], block[5]"
        )

    def test_paginate_false_refuses_many_scalars_with_the_pre_w3_10_text(self) -> None:
        tree = _with_extra_scalars(run(toy()), 23)
        with pytest.raises(ResultsError, match="the limit is 20"):
            plot_corner(tree, paginate=False)

    def test_the_warning_names_the_cap_and_the_var_names_route(self) -> None:
        tree = _with_extra_scalars(run(toy()), 23)
        with pytest.warns(ResultsWarning, match="max_variables=20") as caught:
            plot_corner(tree)
        assert any("var_names=" in str(warning.message) for warning in caught)

    def test_labels_and_truths_are_resolved_across_every_page(self) -> None:
        tree = _with_extra_scalars(run(toy()), 23)
        names = [str(name) for name in tree["posterior"].dataset.data_vars]
        labels = [f"L{index}" for index in range(len(names))]
        truths = dict(zip(names, range(len(names)), strict=True))
        with pytest.warns(ResultsWarning):
            figures = plot_corner(tree, labels=labels, truths=truths)
        assert len(figures) == 2
        assert (
            figure_metadata(figures[0])["corner.draws"]
            == figure_metadata(figures[1])["corner.draws"]
        )

    def test_a_mismatched_label_count_is_refused_before_any_page_is_drawn(self) -> None:
        tree = _with_extra_scalars(run(toy()), 23)
        with pytest.raises(ResultsError, match="per column, not per variable"):
            plot_corner(tree, labels=["only one"])


# ---------------------------------------------------------------------------
# plot_trace
# ---------------------------------------------------------------------------


class TestTrace:
    def test_it_renders_and_counts_the_gaps(self) -> None:
        figure = plot_trace(run(toy(), reject=True))
        assert figure_metadata(figure)["trace.rejected_draws"] == "1"

    def test_a_prior_rejected_draw_is_a_gap_and_not_a_zero(self) -> None:
        # inference.md §18(c)'s distinction, made visible: the trace line must
        # break at the rejected draw rather than dropping to zero, which is a
        # value the parameter never took.
        figure = plot_trace(run(toy(), reject=True))
        trace = figure.axes[1]
        drawn = np.concatenate([line.get_ydata() for line in trace.get_lines()])
        assert np.isnan(drawn).sum() == 1
        assert not np.any(drawn == 0.0)

    def test_it_reads_sample_stats_and_draws_lp_beside_the_parameters(self) -> None:
        figure = plot_trace(run(toy()))
        # Two parameters plus lp, two panels each.
        assert len(figure.axes) == 6

    def test_combined_pools_the_chains(self) -> None:
        separate = plot_trace(run(toy()))
        combined = plot_trace(run(toy()), combined=True)
        assert len(separate.axes[1].get_lines()) == 2
        assert len(combined.axes[1].get_lines()) == 1

    def test_it_refuses_a_tree_that_is_not_a_run(self) -> None:
        with pytest.raises(ResultsError, match="not an emitted run"):
            plot_trace(None)

    def test_it_warns_once_on_an_approximate_run(self) -> None:
        """W5.0: ``ampere_approximation`` not ``"none"`` is one loud warning."""
        tree = run(toy())
        tree.attrs["ampere_approximation"] = "mean_field"
        with pytest.warns(ResultsWarning, match="mean_field"):
            plot_trace(tree)

    def test_it_says_nothing_for_an_exact_runs_default(self) -> None:
        with warnings.catch_warnings():
            warnings.simplefilter("error", ResultsWarning)
            plot_trace(run(toy()))


class TestTracePaging:
    """W3.10: plot_trace pages above MAX_TRACE_VARIABLES exactly as plot_corner does."""

    def test_it_pages_above_the_cap_with_a_warning_and_lp_on_every_page(self) -> None:
        tree = run(toy())
        block = np.random.default_rng(5).normal(size=(2, 30, 90))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"latent.z": (("chain", "draw", "latent.z_dim_0"), block)}
        )
        with pytest.warns(ResultsWarning, match="3 pages"):
            figures = plot_trace(tree, var_names=["latent.z"])
        assert len(figures) == 3
        # 40 + 40 + 10 rows of latent.z, plus one lp row on every page.
        assert [len(figure.axes) // 2 for figure in figures] == [41, 41, 11]
        assert [figure_metadata(figure)["page"] for figure in figures] == [
            "1 of 3",
            "2 of 3",
            "3 of 3",
        ]

    def test_a_call_within_the_cap_still_returns_one_figure(self) -> None:
        figure = plot_trace(run(toy()))
        assert not isinstance(figure, list)
        assert "page" not in figure_metadata(figure)

    def test_paginate_false_refuses_above_the_cap_exactly_as_before_w3_10(self) -> None:
        tree = run(toy())
        block = np.random.default_rng(5).normal(size=(2, 30, 90))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"latent.z": (("chain", "draw", "latent.z_dim_0"), block)}
        )
        with pytest.raises(ResultsError, match="refused loudly rather than attempted"):
            plot_trace(tree, var_names=["latent.z"], paginate=False)

    def test_rejected_draws_are_counted_on_every_page(self) -> None:
        tree = run(toy(), reject=True)
        block = np.random.default_rng(6).normal(size=(2, 30, 45))
        tree["posterior"] = tree["posterior"].dataset.assign(
            {"latent.z": (("chain", "draw", "latent.z_dim_0"), block)}
        )
        with pytest.warns(ResultsWarning):
            figures = plot_trace(tree, var_names=["latent.z"])
        assert len(figures) == 2
        assert all(figure_metadata(figure)["trace.rejected_draws"] == "1" for figure in figures)


# ---------------------------------------------------------------------------
# add_posterior_predictive and plot_posterior_predictive
# ---------------------------------------------------------------------------


class TestPosteriorPredictive:
    def test_replicates_are_derived_on_demand_and_not_stored_by_default(self) -> None:
        problem = toy()
        tree = run(problem, draws=8)
        assert POSTERIOR_PREDICTIVE_GROUP not in tree.children
        add_posterior_predictive(tree, problem)
        replicates = tree[POSTERIOR_PREDICTIVE_GROUP].dataset
        assert replicates["sed"].dims == ("chain", "draw", "sed_spectral_axis")
        values = np.asarray(replicates["sed"].values)
        assert np.all(np.isfinite(values))
        # They are *drawn*, not the observations copied: the scatter about the
        # observations is of order sigma.
        spread = float(np.std(values - np.asarray(tree["observed_data"]["sed"].values)))
        assert 0.3 * SIGMA < spread < 3.0 * SIGMA

    def test_the_named_sub_stream_is_recorded_and_is_not_simulate(self) -> None:
        # lowering.md §9.2: adding a predictive check must not change an SBI
        # budget's draws, which is why the stream is separate by default.
        problem = toy()
        tree = add_posterior_predictive(run(problem, draws=4), problem)
        attrs = tree[POSTERIOR_PREDICTIVE_GROUP].attrs
        assert attrs["ampere_seed_stream"] == "posterior_predictive"

    def test_masked_samples_are_nan_rather_than_the_observed_value(self) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[5:8] = True
        problem = toy(mask=mask)
        tree = add_posterior_predictive(run(problem, draws=4), problem)
        values = np.asarray(tree[POSTERIOR_PREDICTIVE_GROUP]["sed"].values)
        assert values.shape[-1] == GRID.size
        assert np.all(np.isnan(values[..., 5:8]))
        assert np.all(np.isfinite(values[..., :5]))

    def test_thinning_keeps_the_draw_correspondence(self) -> None:
        problem = toy()
        tree = add_posterior_predictive(run(problem, draws=6), problem, thin=3)
        assert tree[POSTERIOR_PREDICTIVE_GROUP]["draw"].values.tolist() == [0, 3]

    def test_the_plot_names_its_precondition_rather_than_a_missing_group(self) -> None:
        # results.md §8, literally: the refusal must say "you have not computed
        # this yet", naming the function.
        with pytest.raises(ResultsError, match="add_posterior_predictive"):
            plot_posterior_predictive(run(toy(), draws=4))

    def test_it_renders_and_reports_a_p_value(self) -> None:
        problem = toy()
        tree = add_posterior_predictive(run(problem, draws=20), problem)
        figure = plot_posterior_predictive(tree)
        metadata = figure_metadata(figure)
        assert 0.0 <= float(metadata["sed.posterior_predictive_p_value"]) <= 1.0
        assert float(metadata["sed.observed_statistic"]) > 0.0

    def test_a_caller_may_supply_the_discrepancy(self) -> None:
        problem = toy()
        tree = add_posterior_predictive(run(problem, draws=8), problem)
        figure = plot_posterior_predictive(
            tree, statistic=lambda values, sigma: float(np.nanmax(np.abs(values)))
        )
        assert "sed.posterior_predictive_p_value" in figure_metadata(figure)

    def test_a_complex_dataset_is_refused_by_name(self) -> None:
        # The derived groups do not mirror results.md §4's real/imag split yet,
        # and inventing half a replicate would be worse than saying so.
        problem = visibility_problem()
        recorder = DrawRecorder(problem, chains=1)
        for _ in range(2):
            recorder.record({"model.flux": 1.0})
        tree = recorder.emit(engine="fixture")
        with pytest.raises(ResultsError, match="complex-valued"):
            add_posterior_predictive(tree, problem)


# ---------------------------------------------------------------------------
# The per-observation log-likelihood group (results.md §6)
# ---------------------------------------------------------------------------


class TestPointwiseLogLikelihood:
    def test_the_factorised_terms_sum_to_the_stored_joint(self) -> None:
        # results.md §6: for independent noise the decomposition is exact.
        problem = toy()
        tree = add_pointwise_log_likelihood(run(problem, draws=6), problem)
        group = tree[POINTWISE_LOG_LIKELIHOOD_GROUP]
        assert group.attrs["ampere_decomposition"] == FACTORISED_DECOMPOSITION
        assert group["sed"].attrs["ampere_decomposition"] == FACTORISED_DECOMPOSITION
        terms = np.asarray(group["sed"].values)
        joint = np.asarray(tree[LOG_LIKELIHOOD_GROUP]["sed"].values)
        assert np.allclose(np.sum(terms, axis=-1), joint)

    def test_a_gp_declares_the_conditional_loo_decomposition(self) -> None:
        # And it deliberately does NOT sum to the joint: a GP likelihood has no
        # per-observation factorisation (likelihoods.md §16).
        problem = gp_toy(DenseGP())
        tree = add_pointwise_log_likelihood(run(problem, draws=4), problem)
        group = tree[POINTWISE_LOG_LIKELIHOOD_GROUP]
        assert group["sed"].attrs["ampere_decomposition"] == CONDITIONAL_LOO_DECOMPOSITION
        terms = np.asarray(group["sed"].values)
        joint = np.asarray(tree[LOG_LIKELIHOOD_GROUP]["sed"].values)
        assert not np.allclose(np.sum(terms, axis=-1), joint)

    def test_a_solver_that_refuses_conditional_loo_refuses_by_name(self) -> None:
        # W2.3 deferred QuasisepGP's O(N) recursion. The refusal must carry the
        # solver's own sentence and must NOT fall back to a dense solve.
        problem = gp_toy(QuasisepGP())
        tree = run(problem, draws=2)
        with pytest.raises(ResultsError, match="QuasisepGP"):
            add_pointwise_log_likelihood(tree, problem)
        assert POINTWISE_LOG_LIKELIHOOD_GROUP not in tree.children

    def test_masked_samples_are_nan_on_the_container_axis(self) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[2:4] = True
        problem = toy(mask=mask)
        tree = add_pointwise_log_likelihood(run(problem, draws=4), problem)
        terms = np.asarray(tree[POINTWISE_LOG_LIKELIHOOD_GROUP]["sed"].values)
        assert terms.shape[-1] == GRID.size
        assert np.all(np.isnan(terms[..., 2:4]))
        assert np.all(np.isfinite(terms[..., :2]))

    def test_it_survives_netcdf_with_its_hashes(self, tmp_path: Path) -> None:
        problem = toy()
        tree = add_pointwise_log_likelihood(run(problem, draws=6), problem)
        written = to_netcdf(tree, tmp_path / "run.nc")
        back = from_netcdf(written)
        assert POINTWISE_LOG_LIKELIHOOD_GROUP in back.children
        assert back.attrs["ampere_spec_hash"] == tree.attrs["ampere_spec_hash"]
        assert back.attrs["ampere_problem_hash"] == tree.attrs["ampere_problem_hash"]
        assert (
            back[POINTWISE_LOG_LIKELIHOOD_GROUP].attrs["ampere_decomposition"]
            == FACTORISED_DECOMPOSITION
        )
        assert np.allclose(
            np.asarray(back[POINTWISE_LOG_LIKELIHOOD_GROUP]["sed"].values),
            np.asarray(tree[POINTWISE_LOG_LIKELIHOOD_GROUP]["sed"].values),
            equal_nan=True,
        )

    def test_arviz_loo_consumes_it(self) -> None:
        # results.md §6: "these are what arviz.loo consumes". The bridge is
        # explicit because the run's own log_likelihood group is per dataset
        # (§15 R6) and calling loo on it computes leave-one-DATASET-out.
        import arviz

        problem = toy()
        tree = add_pointwise_log_likelihood(run(problem, chains=2, draws=40), problem)
        view = pointwise_as_log_likelihood(tree)
        result = arviz.loo(view, var_name="sed")
        assert np.isfinite(float(result.elpd))
        # The original run is untouched: its log_likelihood is still per dataset.
        assert tree[LOG_LIKELIHOOD_GROUP]["sed"].dims == ("chain", "draw")

    def test_the_bridge_names_its_precondition(self) -> None:
        with pytest.raises(ResultsError, match="add_pointwise_log_likelihood"):
            pointwise_as_log_likelihood(run(toy(), draws=2))


# ---------------------------------------------------------------------------
# W5.3: the four plots on a point kind with several axes
# ---------------------------------------------------------------------------


class TestMultiAxisCoordinate:
    """``VisibilitySet``/``ClosurePhases`` (several axes) against a plain
    single-axis kind's existing rows, which stay exactly as they were.
    """

    def test_a_single_axis_kind_still_resolves_exactly_as_before(self) -> None:
        # The byte-identical guarantee, pinned: passing coordinate=None (the
        # new default) to a Spectrum-backed run is the pre-W5.3 call.
        problem = toy()
        tree = add_posterior_predictive(run(problem, draws=6), problem)
        figure = plot_posterior_predictive(tree)
        assert figure is not None
        pyplot.close(figure)

    def test_visibility_set_renders_with_its_default_baseline_length_coordinate(
        self,
    ) -> None:
        problem = visibility_problem()
        tree = _tiny_run(problem, values={"model.flux": 1.0})
        tree = add_posterior_predictive(tree, problem, datasets=["vis"], component="abs")
        figure = plot_posterior_predictive(tree, datasets=["vis"], component="abs")
        try:
            assert figure is not None
            axis = figure.axes[0]
            assert "baseline length" in axis.get_xlabel()
        finally:
            pyplot.close(figure)

    def test_visibility_set_residuals_render_with_the_default_coordinate(self) -> None:
        problem = visibility_problem()
        tree = _tiny_run(problem, values={"model.flux": 1.0})
        tree = add_residuals(tree, problem, datasets=["vis"], component="abs")
        figure = plot_residuals(tree, datasets=["vis"], component="abs", whiteness=False)
        try:
            assert "baseline length" in figure.axes[0].get_xlabel()
        finally:
            pyplot.close(figure)

    def test_closure_phases_renders_with_its_default_longest_baseline_coordinate(
        self,
    ) -> None:
        problem = closure_problem()
        tree = _tiny_run(problem, values={"model.phase": 0.1})
        tree = add_posterior_predictive(tree, problem, datasets=["t3"])
        figure = plot_posterior_predictive(tree, datasets=["t3"])
        try:
            assert "longest baseline" in figure.axes[0].get_xlabel()
        finally:
            pyplot.close(figure)

    def test_closure_phases_residuals_render_with_the_default_coordinate(self) -> None:
        problem = closure_problem()
        tree = _tiny_run(problem, values={"model.phase": 0.1})
        tree = add_residuals(tree, problem, datasets=["t3"])
        figure = plot_residuals(tree, datasets=["t3"], whiteness=False)
        try:
            assert "longest baseline" in figure.axes[0].get_xlabel()
        finally:
            pyplot.close(figure)

    def test_gp_localisation_and_anomaly_score_render_for_the_gp_fitted_visibility_set(
        self,
    ) -> None:
        problem = gp_visibility_problem()
        tree = _tiny_run(problem, values={"model.flux": 1.0})
        tree = gp_localisation(tree, problem, datasets=["vis"], component="abs")
        figure = plot_gp_localisation(tree, datasets=["vis"], component="abs")
        try:
            assert "baseline length" in figure.axes[0].get_xlabel()
        finally:
            pyplot.close(figure)
        score = gp_localisation_score(tree, dataset="vis", component="abs")
        axes = plot_anomaly_score(score)
        try:
            assert axes.get_legend() is not None
        finally:
            pyplot.close(axes.get_figure())

    def test_it_refuses_a_multi_axis_kind_with_no_default_and_no_coordinate(self) -> None:
        # tests/core/thirdparty_polarimeter.py's pattern, applied to a kind
        # that declares nothing beyond AXES/LAYOUT/ALLOW_COMPLEX (TwoAxisPoint
        # above): results.md §13 item 14's refusal, lifted everywhere else,
        # still stands where a kind opts into nothing.
        problem = pair_problem()
        tree = _tiny_run(problem, values={"model.level": 0.15})
        tree = add_posterior_predictive(tree, problem, datasets=["pair"])
        with pytest.raises(ResultsError, match="no default plotted coordinate"):
            plot_posterior_predictive(tree, datasets=["pair"])

    def test_a_coordinate_argument_names_one_of_the_kinds_own_axes(self) -> None:
        problem = pair_problem()
        tree = _tiny_run(problem, values={"model.level": 0.15})
        tree = add_posterior_predictive(tree, problem, datasets=["pair"])
        figure = plot_posterior_predictive(tree, datasets=["pair"], coordinate="p")
        try:
            assert figure.axes[0].get_xlabel() == "p"
        finally:
            pyplot.close(figure)

    def test_a_coordinate_callable_computes_one_from_the_kinds_axes(self) -> None:
        def total(axes: Any) -> tuple[np.ndarray, str]:
            return axes["p"].values + axes["q"].values, "p + q"

        problem = pair_problem()
        tree = _tiny_run(problem, values={"model.level": 0.15})
        tree = add_posterior_predictive(tree, problem, datasets=["pair"])
        figure = plot_posterior_predictive(tree, datasets=["pair"], coordinate=total)
        try:
            assert figure.axes[0].get_xlabel() == "p + q"
        finally:
            pyplot.close(figure)

    def test_a_complex_dataset_with_no_component_names_the_four_choices(self) -> None:
        problem = visibility_problem()
        tree = _tiny_run(problem, values={"model.flux": 1.0})
        with pytest.raises(ResultsError, match='"real", "imag", "abs" or "phase"'):
            add_posterior_predictive(tree, problem, datasets=["vis"])

    def test_component_on_a_real_dataset_is_refused(self) -> None:
        problem = closure_problem()
        tree = _tiny_run(problem, values={"model.phase": 0.1})
        with pytest.raises(ResultsError, match="real-valued"):
            add_posterior_predictive(tree, problem, datasets=["t3"], component="abs")

    def test_an_unknown_component_is_refused_by_name(self) -> None:
        problem = visibility_problem()
        tree = _tiny_run(problem, values={"model.flux": 1.0})
        with pytest.raises(ResultsError, match="component must be one of"):
            add_posterior_predictive(tree, problem, datasets=["vis"], component="modulus")
