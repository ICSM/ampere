"""``examples/image`` builds, negotiates and fits end to end (W5.5).

``tests/examples/test_interferometry_study.py``'s shape, applied to the
gridded modality: fast, structural, and a tiny-budget fit that runs end to end
on the reference backend, plus the CLI. The study's own longer runs — the
documentation budget, the full SBC count, and the 128x128 benchmark cell — are
behind the ``image_full`` marker for ``tests/m2/conftest.py``'s reason: they
are tens of minutes and do not belong in a per-PR gate.

The pinned claim here is the **direction**, not a number: the misspecified arm
is biased in units of its own posterior width, the flexible arm is less so, and
twelve SBC replicas are a smoke budget in Talts et al.'s sense — evidence of a
gross effect. Anything finer needs the marked run.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import Image, Layout

from examples.image import generators, model, study

#: Small enough to run in a few seconds; large enough that every rule the study
#: depends on is exercised (a padded grid, a crop, a 2-D GP).
TINY = study.EmceeBudget(walkers=8, steps=40, burn_in=10)


def _prior_median(problem: object) -> dict[str, float]:
    """Each free parameter at the median of its own prior.

    ``prior_transform`` maps the unit cube, so a vector of halves is the
    per-margin median — a legal, finite point on every arm without any test
    having to know what the flexible arm's kernel parameters are called.
    """
    values = problem.prior_transform(np.full(problem.free_size, 0.5))  # type: ignore[attr-defined]
    return dict(zip(problem.parameters.free_names, (float(v) for v in values), strict=True))  # type: ignore[attr-defined]


#: The image the fast rows fit. Smaller than the study's own ``SMALL_PIXELS``:
#: a 16x16 image is a 256-square dense GP solve, which a 40-step ensemble can
#: afford in a gate.
TINY_PIXELS = 16


class TestTheTruthComposes:
    def test_synthetic_produces_a_gridded_image_with_uncertainties(self) -> None:
        import ampere.backends.reference as ref
        import ampere.backends.reference.interferometry as itf

        _, observed = generators.synthetic(ref, itf, TINY_PIXELS)
        assert isinstance(observed, Image)
        assert observed.LAYOUT is Layout.GRID
        assert observed.values.shape == (TINY_PIXELS, TINY_PIXELS)
        assert observed.uncertainty is not None
        assert np.all(np.asarray(observed.uncertainty) > 0.0)

    def test_the_background_is_fainter_per_pixel_and_much_broader(self) -> None:
        """The misspecification has to be *smooth*, or no kernel could absorb it."""
        assert generators.BACKGROUND["fwhm"] > 4.0 * generators.SOURCE["fwhm"]
        assert generators.BACKGROUND["flux"] < generators.SOURCE["flux"]

    def test_the_psf_is_wider_than_the_pixel_scale_at_every_benchmark_size(self) -> None:
        """Otherwise the convolution is a near-identity and the study is vacuous."""
        for pixels in study.BENCHMARK_SIZES:
            assert generators.PSF_FWHM > generators.pixel_scale(pixels)


class TestTheStepNegotiates:
    def test_the_negotiated_grid_is_padded_and_the_prediction_is_not(self) -> None:
        import ampere.backends.reference as ref
        import ampere.backends.reference.interferometry as itf

        from ampere.core import negotiate

        _, observed = generators.synthetic(ref, itf, TINY_PIXELS)
        instrument = generators.camera(ref, observed)
        needed = negotiate([instrument])
        padded = np.asarray(needed["sky"]["x"].coordinates().to_value("mas"))
        assert padded.size > TINY_PIXELS
        compiled = generators.truth_model(itf, TINY_PIXELS).compile_for(needed)
        model_image = compiled.evaluate()["sky"]
        # The model adopted the padded grid, not the observed one: the native
        # ``native_grid`` surface is the modern backends' and does not exist on
        # the reference path, so the container's own axis is what to read here.
        assert model_image.x.values.size == padded.size
        predicted = instrument(compiled.evaluate())
        assert predicted.values.shape == (TINY_PIXELS, TINY_PIXELS)
        assert np.array_equal(predicted.x.values, observed.x.values)


class TestBuildProblem:
    """One truth, three arms, each a one-dataset problem on one ``sky`` channel."""

    @pytest.mark.parametrize("arm", study.ARMS)
    def test_each_arm_composes_a_one_dataset_problem(self, arm: str) -> None:
        problem = study.build_problem("reference", arm, pixels=TINY_PIXELS)
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"image"}
        # ``fwhm`` before ``flux``: GaussianSource registers them in that order,
        # and a parameter's position is its declaration order (parameters.md §13).
        assert problem.parameters.free_names[:2] == ("model.fwhm", "model.flux")

    def test_the_flexible_arm_has_two_more_free_parameters(self) -> None:
        incomplete = study.build_problem("reference", "incomplete", pixels=TINY_PIXELS)
        flexible = study.build_problem("reference", "flexible", pixels=TINY_PIXELS)
        assert flexible.free_size == incomplete.free_size + 2

    def test_only_the_correct_arm_carries_the_background(self) -> None:
        """And the misspecified arms do not have it at zero: they do not have it."""
        correct = study.model_for("reference", "correct")
        incomplete = study.model_for("reference", "incomplete")
        assert correct._background is not None
        assert incomplete._background is None

    def test_an_unknown_arm_is_refused_by_name(self) -> None:
        with pytest.raises(ValueError, match="unknown arm"):
            study.model_for("reference", "wrong")

    def test_the_log_prob_is_finite_on_every_arm(self) -> None:
        for arm in study.ARMS:
            problem = study.build_problem("reference", arm, pixels=TINY_PIXELS)
            assert np.isfinite(problem.log_prob(_prior_median(problem)))


class TestATinyFitRunsEndToEnd:
    @pytest.mark.parametrize("arm", ["correct", "flexible"])
    def test_the_arm_recovers_something_finite(self, arm: str) -> None:
        problem = study.build_problem("reference", arm, pixels=TINY_PIXELS, seed=1)
        run = study.run(problem, TINY)
        summaries = study.summarise(run, names=list(generators.TRUTH))
        for summary in summaries.values():
            assert np.isfinite(summary.median)
            assert summary.width > 0.0

    def test_simulate_draws_an_image(self) -> None:
        """The ``Layout.GRID`` draw path, through a whole problem rather than a step."""
        problem = study.build_problem("reference", "correct", pixels=TINY_PIXELS, seed=2)
        simulation = problem.simulate(_prior_median(problem), observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        drawn = simulation.observations["image"]
        assert isinstance(drawn, Image)
        assert drawn.values.shape == (TINY_PIXELS, TINY_PIXELS)


class TestTheBenchmark:
    """The table's *shape* is pinned; its numbers are informational."""

    def test_the_benchmark_reports_one_cell_per_size_and_solver(self) -> None:
        costs = study.benchmark_solvers(sizes=(12, 16))
        assert len(costs) == 4
        assert {cost.solver for cost in costs} == {
            "DenseGP",
            f"HSGP({study.BASIS_PER_AXIS}x{study.BASIS_PER_AXIS})",
        }
        for cost in costs:
            assert cost.samples == cost.pixels**2
            assert cost.seconds > 0.0
            assert cost.peak_mib > 0.0
            assert np.isfinite(cost.log_prob)

    def test_the_table_renders_every_cell(self) -> None:
        costs = study.benchmark_solvers(sizes=(12,))
        table = study.benchmark_table(costs)
        assert "solver" in table
        assert table.count("\n") == len(costs) + 1

    @pytest.mark.image_full
    def test_the_full_three_size_benchmark(self) -> None:
        """The documentation budget: 128x128 is a 16,384-square dense Cholesky."""
        costs = study.benchmark_solvers()
        assert len(costs) == 2 * len(study.BENCHMARK_SIZES)


class TestCalibration:
    """The coverage claim, at a smoke budget."""

    def test_the_flexible_arm_covers_at_a_smoke_budget(self) -> None:
        calibration = study.run_calibration(
            "flexible",
            pixels=TINY_PIXELS,
            count=6,
            draws=60,
            budget=study.EmceeBudget(walkers=8, steps=60, burn_in=20),
        )
        coverage = study.coverage_at(calibration, 0.9)
        assert coverage.shape == (len(generators.TRUTH),)
        assert np.all(coverage >= 0.0)
        assert np.all(coverage <= 1.0)
        # Six replicas cannot resolve 0.9 from 0.75; what they can say is that
        # the interval is not grossly under-covering, which is the direction
        # this arm exists to demonstrate.
        assert float(np.mean(coverage)) >= 0.5

    @pytest.mark.image_full
    def test_the_full_calibration_row(self) -> None:
        """The pinned coverage comparison, at the study's own budget."""
        covered = {
            arm: float(np.mean(study.coverage_at(study.run_calibration(arm), 0.9)))
            for arm in ("incomplete", "flexible")
        }
        assert covered["flexible"] >= covered["incomplete"]


class TestTheCli:
    def test_main_runs_a_tiny_study_and_prints_a_table(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        code = main_with(["--arms", "correct", "incomplete", "--pixels", str(TINY_PIXELS)])
        assert code == 0
        captured = capsys.readouterr().out
        assert "model.flux" in captured
        assert "correct" in captured
        assert "incomplete" in captured

    def test_main_benchmark_flag_prints_the_table(self, capsys: pytest.CaptureFixture[str]) -> None:
        code = main_with(["--benchmark", "--sizes", "12"])
        assert code == 0
        captured = capsys.readouterr().out
        assert "DenseGP" in captured
        assert "peak MiB" in captured

    def test_no_dense_leaves_the_exact_solver_out(self, capsys: pytest.CaptureFixture[str]) -> None:
        """The escape from the largest cell's ~10 GiB, reachable from the CLI.

        ``benchmark_solvers`` has had ``include_dense`` since the study landed;
        without this flag the documented entry point could not reach it, so the
        only way to measure the approximate solver at 128x128 was to allocate
        five copies of a 16,384-square covariance first.
        """
        code = main_with(["--benchmark", "--sizes", "12", "--no-dense"])
        assert code == 0
        captured = capsys.readouterr().out
        assert "HSGP" in captured
        assert "DenseGP" not in captured


def main_with(argv: list[str]) -> int:
    """The CLI, imported lazily so this module's collection stays cheap."""
    from examples.image.__main__ import main

    return main(argv)


class TestSourceWithBackgroundOwnsExactlyTheSourcesParameters:
    """A regression row for :mod:`examples.image.model`'s own claim."""

    def test_the_background_never_becomes_a_free_parameter(self) -> None:
        import ampere.backends.reference.interferometry as itf

        grid = generators.seed_grid()
        composite = model.source_with_background(
            itf,
            grid,
            grid,
            flux=generators.SOURCE["flux"],
            fwhm=generators.SOURCE["fwhm"],
            background_flux=generators.BACKGROUND["flux"],
            background_fwhm=generators.BACKGROUND["fwhm"],
        )
        assert composite.parameters.free_names == ()
