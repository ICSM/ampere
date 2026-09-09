"""``examples/wstat_comparison.py`` runs end to end (W2.9's accept criterion).

The docs build cannot be trusted as the sole proof that this example runs:
``docs/source/conf.py`` sets ``nbsphinx_allow_errors = True`` (so a failed
notebook does not fail the Sphinx build), and this example is a literal
script included via ``:download:`` rather than an executed notebook in any
case (see ``docs/source/wstat_comparison.rst`` for why: this repository
commits no run outputs or other binary artefacts, and an executed notebook
with embedded figures would be exactly that). This module is the reliable,
environment-independent gate instead — it imports the example as a module
and exercises every code path :func:`wstat_comparison.main` uses, plus the
two specific correctness claims the work item's binding decisions make:
masking safety (``NoiseParams.retain`` keeps the background buffer aligned)
and ``sample()``'s refusal.

Budgets here are deliberately smaller than the example's own "doc" budget
(:data:`wstat_comparison.DOC_STEPS` and friends, tuned to run in well under a
minute standalone) — this suite only needs to prove the machinery runs and
produces sane output, not to reproduce the doc page's own posterior.
"""

from __future__ import annotations

import warnings
from types import ModuleType

import numpy as np
import pytest

from ampere.core import Likelihood
from ampere.core.exceptions import LikelihoodError


@pytest.fixture(scope="module")
def synthetic_data(wstat_example: ModuleType) -> tuple[np.ndarray, np.ndarray]:
    return wstat_example.synthetic_xray_counts()


class TestSyntheticData:
    def test_shapes_and_nonnegativity(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        assert source_counts.shape == (wstat_example.N_BINS,)
        assert background_counts.shape == (wstat_example.N_BINS,)
        assert np.all(source_counts >= 0.0)
        assert np.all(background_counts >= 0.0)
        assert np.all(source_counts == np.round(source_counts))
        assert np.all(background_counts == np.round(background_counts))

    def test_deterministic(self, wstat_example: ModuleType) -> None:
        first = wstat_example.synthetic_xray_counts()
        second = wstat_example.synthetic_xray_counts()
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[1], second[1])


class TestProfiledCashWithBackgroundFamily:
    """The binding decisions: masking safety, ``sample()``'s refusal, preconditions."""

    def test_masking_excises_the_background_buffer_the_same_way_as_dropping_channels(
        self, wstat_example: ModuleType
    ) -> None:
        """The awkward_instrument.md §5 demonstration, automated for this family.

        A masked channel must be scored identically to one that was never
        there at all — the same claim ``likelihoods.md`` §8 makes about
        ``weights()``/excision in general, checked here for a family that
        carries its *own* aligned data (gap X-2) rather than only data
        ``Likelihood`` already excises on its behalf.
        """
        import astropy.units as u

        from ampere.core import IndependentNoise, Spectrum

        energy = np.array([1.0, 2.0, 3.0, 4.0])
        background = np.array([2.0, 3.0, 1.0, 4.0])
        ratio = 0.5

        family_full = wstat_example.ProfiledCashWithBackground(background, ratio)
        likelihood = Likelihood(family_full, IndependentNoise())
        observed = Spectrum(
            energy * u.keV, [5.0, 6.0, 2.0, 7.0], mask=np.array([False, True, False, False])
        )
        predicted = observed.with_values([4.0, 5.0, 3.0, 6.0])
        masked_log_prob = likelihood.log_prob(predicted, observed)

        family_excised = wstat_example.ProfiledCashWithBackground(background[[0, 2, 3]], ratio)
        excised_likelihood = Likelihood(family_excised, IndependentNoise())
        excised_observed = Spectrum(np.array([1.0, 3.0, 4.0]) * u.keV, [5.0, 2.0, 7.0])
        excised_predicted = excised_observed.with_values([4.0, 3.0, 6.0])
        excised_log_prob = excised_likelihood.log_prob(excised_predicted, excised_observed)

        assert masked_log_prob == pytest.approx(excised_log_prob, abs=1e-12)

    def test_sample_refuses_with_a_reason_specific_to_the_profiling(
        self, wstat_example: ModuleType
    ) -> None:
        family = wstat_example.ProfiledCashWithBackground(np.array([1.0, 2.0]), 0.4)
        with pytest.raises(LikelihoodError, match="profiled background estimate"):
            family.sample(np.array([1.0, 2.0]), None, np.random.default_rng(0))

    def test_check_observed_rejects_non_integer_source_counts(
        self, wstat_example: ModuleType
    ) -> None:
        import astropy.units as u

        from ampere.core import IndependentNoise, Spectrum

        family = wstat_example.ProfiledCashWithBackground(np.array([1.0, 2.0]), 0.4)
        likelihood = Likelihood(family, IndependentNoise())
        counts = Spectrum([1.0, 2.0] * u.keV, [4.5, 7.0])
        with pytest.raises(LikelihoodError, match="non-negative integer source counts"):
            likelihood.check_alignment(counts.with_values([4.0, 7.0]), counts)

    def test_construction_refuses_a_negative_ratio(self, wstat_example: ModuleType) -> None:
        with pytest.raises(LikelihoodError, match="positive exposure/area ratio"):
            wstat_example.ProfiledCashWithBackground(np.array([1.0, 2.0]), -0.1)

    def test_log_prob_matches_a_brute_force_profile_search(self, wstat_example: ModuleType) -> None:
        """The closed-form profile MLE against numerical optimisation, several draws."""
        import scipy.optimize as opt
        import scipy.stats as st

        rng = np.random.default_rng(2026090801)
        for _ in range(25):
            source = float(rng.integers(0, 20))
            background = float(rng.integers(0, 20))
            model = rng.uniform(0.05, 20.0)
            ratio = rng.uniform(0.05, 2.0)

            family = wstat_example.ProfiledCashWithBackground(np.array([background]), ratio)
            from ampere.core import NoiseParams

            noise = NoiseParams(sigma=None, values={}, retain=np.array([True]))
            analytic = family.log_prob(np.array([model]), np.array([source]), noise)

            def negative_log_likelihood(
                b: float, *, source=source, background=background, model=model, ratio=ratio
            ) -> float:
                b = max(b, 1e-300)
                return -(
                    st.poisson.logpmf(source, model + ratio * b) + st.poisson.logpmf(background, b)
                )

            numeric = -opt.minimize_scalar(
                negative_log_likelihood, bounds=(1e-12, 1e4), method="bounded"
            ).fun
            assert analytic == pytest.approx(numeric, abs=1e-6)


class TestBothProblemsBuildAndEvaluate:
    def test_wstat_problem_evaluates_finitely_at_truth(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        problem = wstat_example.build_wstat_problem(source_counts, background_counts)
        assert problem.parameters.free_names == ("model.norm", "model.index")
        theta = {
            "model.norm": wstat_example.TRUTH["src_norm"],
            "model.index": wstat_example.TRUTH["src_index"],
        }
        value = problem.log_prob(np.array([theta[n] for n in problem.parameters.free_names]))
        assert np.isfinite(value)

    def test_joint_problem_evaluates_finitely_at_truth(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        problem = wstat_example.build_joint_problem(source_counts, background_counts)
        assert problem.parameters.free_names == (
            "model.src_norm",
            "model.src_index",
            "model.bkg_norm",
        )
        theta = {
            "model.src_norm": wstat_example.TRUTH["src_norm"],
            "model.src_index": wstat_example.TRUTH["src_index"],
            "model.bkg_norm": wstat_example.TRUTH["bkg_norm"],
        }
        value = problem.log_prob(np.array([theta[n] for n in problem.parameters.free_names]))
        assert np.isfinite(value)


class TestEndToEnd:
    """The work item's own escape hatch: exercise the example end to end here.

    A much smaller budget than the doc page's own — enough draws for a
    ``(chain, draw)`` shape and finite log-likelihoods, not enough for a
    converged posterior, which the doc page's own run (§ ``run_doc_budget``
    below) is what claims that.
    """

    def test_wstat_route_runs_and_produces_a_summary(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        problem = wstat_example.build_wstat_problem(source_counts, background_counts)
        run = wstat_example.run_engine(problem, walkers=8, steps=30, burn_in=10)
        assert run["posterior"].dataset.sizes["chain"] == 8
        assert run["posterior"].dataset.sizes["draw"] == 20
        summary = wstat_example.summarise(run, ("model.norm", "model.index"))
        for median, low, high in summary.values():
            assert np.isfinite(median) and np.isfinite(low) and np.isfinite(high)
            assert low <= median <= high

    def test_joint_route_runs_and_produces_a_summary(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        problem = wstat_example.build_joint_problem(source_counts, background_counts)
        run = wstat_example.run_engine(problem, walkers=8, steps=30, burn_in=10)
        assert run["posterior"].dataset.sizes["chain"] == 8
        assert run["posterior"].dataset.sizes["draw"] == 20
        summary = wstat_example.summarise(
            run, ("model.src_norm", "model.src_index", "model.bkg_norm")
        )
        for median, low, high in summary.values():
            assert np.isfinite(median) and np.isfinite(low) and np.isfinite(high)
            assert low <= median <= high

    def test_compare_states_the_recommendation_not_false_balance(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        source_counts, background_counts = synthetic_data
        wstat_problem = wstat_example.build_wstat_problem(source_counts, background_counts)
        wstat_run = wstat_example.run_engine(wstat_problem, walkers=8, steps=30, burn_in=10)
        wstat_summary = wstat_example.summarise(wstat_run, ("model.norm", "model.index"))

        joint_problem = wstat_example.build_joint_problem(source_counts, background_counts)
        joint_run = wstat_example.run_engine(joint_problem, walkers=8, steps=30, burn_in=10)
        joint_summary = wstat_example.summarise(
            joint_run, ("model.src_norm", "model.src_index", "model.bkg_norm")
        )

        text = wstat_example.compare(wstat_summary, joint_summary)
        assert "Recommendation: prefer the two-dataset Bayesian joint fit" in text
        assert "Trade-off" in text

    def test_main_runs_end_to_end_at_the_doc_budget(
        self, wstat_example: ModuleType, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """The exact entry point ``python examples/wstat_comparison.py`` uses.

        This is the closest thing to "the example ran, as shipped" a pytest
        gate can offer: it calls the real ``main()`` at the real doc budget
        (well under a minute — see ``wstat_comparison.py``'s own module
        docstring), not a scaled-down stand-in.
        """
        wstat_example.main()
        printed = capsys.readouterr().out
        assert "Recommendation: prefer the two-dataset Bayesian joint fit" in printed


@pytest.fixture(scope="module")
def reduced_study(wstat_example: ModuleType) -> dict[str, object]:
    """One reduced coverage study, shared by the checks below.

    Six simulations, short chains: enough to prove the machinery runs and that
    both routes are ranked against the same datasets, and nowhere near enough
    to say anything about coverage. Module-scoped because it is the only
    expensive thing in this file.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        return wstat_example.coverage_study(count=6, draws=40, walkers=8, steps=120, burn_in=40)


class TestTheCoverageStudy:
    """W3.6: the repeated-trial coverage study the single run cannot substitute for.

    ``examples/wstat_comparison.py``'s own paragraph on the profiled nuisance
    says that "a repeated-trial coverage study would show it reliably; a single
    seeded run, by construction, only shows one draw from that distribution".
    W3.6 makes that study cheap enough to ship, and this class is the reduced
    version of it — the same code path at a budget a per-PR gate can afford, so
    that the docs page's numbers come from machinery this suite has actually
    run.

    What is *not* asserted here is the direction of the effect. At eight
    simulations the uniformity test has no power at all, and asserting that
    WStat's coverage comes out low would be asserting a coin flip; the full
    budget (``python examples/wstat_comparison.py --coverage --full``) is what
    the docs page quotes. What is asserted is that the machinery runs, that
    both routes are ranked against the same simulations, and that the
    generative subclass really does make the joint problem simulable — which is
    the one piece that could break silently, because ``PoissonFamily.sample()``
    refuses by default.
    """

    def test_the_generative_subclass_makes_the_joint_problem_simulable(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """``PoissonFamily.sample()`` refuses; :class:`CountingPoisson`'s does not."""
        source_counts, background_counts = synthetic_data
        plain = wstat_example.build_joint_problem(source_counts, background_counts)
        with pytest.raises(Exception, match="does not implement sample"):
            plain.simulate(observe=True, rng=np.random.default_rng(0), stream="probe")

        generative = wstat_example.build_joint_problem(
            source_counts, background_counts, generative=True
        )
        simulation = generative.simulate(observe=True, rng=np.random.default_rng(0))
        assert not simulation.failed
        assert simulation.observations is not None
        for label in ("src", "bkg"):
            drawn = np.asarray(simulation.observations[label].values, dtype=float)
            assert drawn.shape == (wstat_example.N_BINS,)
            assert np.all(drawn >= 0.0)
            assert np.all(drawn == np.round(drawn))

    def test_the_generative_subclass_scores_identically_to_the_base_family(
        self, wstat_example: ModuleType, synthetic_data: tuple[np.ndarray, np.ndarray]
    ) -> None:
        """Adding ``sample()`` must not change the density, or the study is of another model."""
        source_counts, background_counts = synthetic_data
        plain = wstat_example.build_joint_problem(source_counts, background_counts)
        generative = wstat_example.build_joint_problem(
            source_counts, background_counts, generative=True
        )
        values = {
            "model.src_norm": wstat_example.TRUTH["src_norm"],
            "model.src_index": wstat_example.TRUTH["src_index"],
            "model.bkg_norm": wstat_example.TRUTH["bkg_norm"],
        }
        assert generative.log_prob(values) == pytest.approx(plain.log_prob(values))

    def test_the_reduced_study_runs_and_ranks_both_routes(self, wstat_example: ModuleType) -> None:
        """The whole study, at a budget the gate can afford."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            study = wstat_example.coverage_study(
                count=6, draws=40, walkers=8, steps=120, burn_in=40
            )
        assert sorted(study) == ["joint", "wstat"]
        for route, expected in (
            ("wstat", ["model.index", "model.norm"]),
            ("joint", ["model.src_index", "model.src_norm"]),
        ):
            result = study[route]
            assert result.sizes["simulation"] == 6
            assert sorted(str(n) for n in result.coords["parameter"].values) == expected
            ranks = np.asarray(result["ranks"].values)
            assert ranks.min() >= 0 and ranks.max() <= 40
            assert result.attrs["ampere_calibration_route"] == "refit"
            assert result.attrs["ampere_calibration_engine"] == "EmceeEngine"
        # Both routes were ranked against the *same* simulations, which is what
        # makes the joint route a control rather than a second experiment.
        assert (
            study["wstat"].attrs["ampere_calibration_seed"]
            == study["joint"].attrs["ampere_calibration_seed"]
        )

    def test_the_summary_states_how_to_read_it_and_quotes_the_numbers(
        self, wstat_example: ModuleType, reduced_study: dict[str, object]
    ) -> None:
        text = wstat_example.coverage_summary(reduced_study)
        assert "The joint route is the control" in text
        assert "mean rank" in text and "68% cov." in text
        # The prior caution is not decoration: a study read without it answers
        # a question about the prior it was run under, not about the statistic.
        assert "averages over" in text
        # The text rank histogram is the page's figure, this repository
        # committing no binary artefacts.
        assert "Rank histograms" in text
        for name in ("model.norm", "model.index", "model.src_norm", "model.src_index"):
            assert name in text

    def test_the_figures_render_headless(
        self, wstat_example: ModuleType, reduced_study: dict[str, object]
    ) -> None:
        import matplotlib

        matplotlib.use("Agg")
        from ampere.results import figure_metadata

        figures = wstat_example.coverage_figure(reduced_study)
        try:
            assert sorted(figures) == [
                "joint_coverage",
                "joint_ranks",
                "wstat_coverage",
                "wstat_ranks",
            ]
            assert any(
                key.endswith(".ks_pvalue") for key in figure_metadata(figures["wstat_ranks"])
            )
        finally:
            # matplotlib keeps every pyplot figure alive until it is closed, and
            # four of them per run of this file is enough to warn about in a
            # whole-suite gate.
            import matplotlib.pyplot as pyplot

            for figure in figures.values():
                pyplot.close(figure)

    def test_the_coverage_flag_is_off_by_default(
        self, wstat_example: ModuleType, capsys: pytest.CaptureFixture[str]
    ) -> None:
        """``main()`` must stay the fast worked example it was.

        The study is minutes, not seconds, so the default entry point still runs
        only the single-run comparison — and says, in the output, that the flag
        exists, because a study nobody knows about is a study nobody runs.
        """
        wstat_example.main([])
        printed = capsys.readouterr().out
        assert "Repeated-trial coverage" not in printed
        assert "--coverage" in printed
