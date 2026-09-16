"""``examples/interferometry`` builds, negotiates and fits end to end (W4.4).

The calibration study's pinned science lives in
``tests/interferometry/test_calibration.py`` (it needs the shared,
session-scoped ``calibrations`` fixture that suite's ``conftest.py``
provides). This module is this item's ``tests/examples`` coverage in
``tests/examples/test_sed_composition.py``'s own shape: fast, structural, and
a tiny-budget single fit that runs end to end on the reference backend —
plus the CLI, ``python -m examples.interferometry``.
"""

from __future__ import annotations

import numpy as np
import pytest

from examples.interferometry import generators, model, study
from examples.interferometry.__main__ import main

TINY = study.EmceeBudget(walkers=8, steps=40, burn_in=10)


class TestTheTruthComposes:
    def test_synthetic_produces_two_observables_on_one_geometry(self) -> None:
        import ampere.backends.reference.interferometry as itf

        _, vis, t3 = generators.synthetic(itf)
        assert vis.n_samples == 24
        assert t3.n_samples == 16

    def test_the_disc_is_fainter_and_more_extended_than_the_binary(self) -> None:
        assert generators.DISC_FLUX < generators.BINARY["flux"]
        assert generators.DISC_FWHM > generators.COMPONENT_FWHM


class TestBuildProblem:
    """One truth, three arms, each a two-dataset problem on one ``sky`` channel."""

    @pytest.mark.parametrize("arm", study.ARMS)
    def test_each_arm_composes_a_two_dataset_problem(self, arm: str) -> None:
        problem = study.build_problem("reference", arm)
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"vis", "t3"}
        assert problem.parameters.free_names[:2] == ("model.separation", "model.flux_ratio")

    def test_the_flexible_arm_has_two_more_free_parameters(self) -> None:
        correct = study.build_problem("reference", "correct")
        flexible = study.build_problem("reference", "flexible")
        assert flexible.free_size == correct.free_size + 2

    def test_an_unknown_arm_is_refused_by_name(self) -> None:
        with pytest.raises(ValueError, match="unknown arm"):
            study.model_for("reference", "wrong", generators.seed_grid(), generators.seed_grid())


class TestATinyFitRunsEndToEnd:
    def test_the_correct_arm_recovers_something_finite(self) -> None:
        problem = study.build_problem("reference", "correct", seed=1)
        run = study.run(problem, TINY)
        summaries = study.summarise(run, names=list(generators.TRUTH))
        for name, summary in summaries.items():
            assert np.isfinite(summary.median)
            assert summary.width > 0.0


class TestTheChromaticArm:
    def test_the_three_kernels_are_distinguishable_objects(self) -> None:
        kernels = {
            kind: study._chromatic_kernel("reference", kind) for kind in study.CHROMATIC_KERNELS
        }
        assert kernels["spatial"] is not kernels["spectral"]
        assert type(kernels["product"]).__name__ == "Product"

    def test_the_chromatic_problem_composes(self) -> None:
        problem = study.chromatic_problem("reference", "product", seed=1)
        assert set(problem.datasets) == {"vis"}


class TestTheCli:
    def test_main_runs_a_tiny_study_and_prints_a_table(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        code = main(["--arms", "correct", "incomplete"])
        assert code == 0
        captured = capsys.readouterr().out
        assert "model.separation" in captured
        assert "correct" in captured
        assert "incomplete" in captured

    def test_main_chromatic_flag_runs_the_three_kernels(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        code = main(["--chromatic"])
        assert code == 0
        captured = capsys.readouterr().out
        for kind in study.CHROMATIC_KERNELS:
            assert kind in captured


class TestBinaryWithDiscOwnsExactlyTheBinarysParameters:
    """A regression row for :mod:`examples.interferometry.model`'s own claim."""

    def test_binary_with_disc_never_frees_the_discs_parameters(self) -> None:
        import ampere.backends.reference.interferometry as itf

        grid = generators.seed_grid()
        composite = model.binary_with_disc(
            itf,
            grid,
            grid,
            disc_flux=generators.DISC_FLUX,
            disc_fwhm=generators.DISC_FWHM,
            **generators.BINARY,
            component_fwhm=generators.COMPONENT_FWHM,
        )
        assert composite.parameters.free_names == ()
