"""``examples/photometry_spectra`` builds, negotiates, ties and fits end to end (W6.2).

Full-budget recovery -- the numbers :doc:`the tutorial page
</photometry_spectra>`'s table quotes, for each of the four ``tie``/``gp``
combinations -- is **not** part of this suite, for the same reason
:mod:`tests.examples.test_sed_composition` gives: the script's default
budget is designed to be informative rather than fast, and milestone M2's
own full ladder sets the precedent for keeping that kind of run out of the
per-PR gate. Instead, there is no pytest coverage of the full budget: it is
verified by running ``python -m examples.photometry_spectra`` directly, once
per combination, and the branch report quotes the recovered intervals and
the wall-clock time from doing exactly that.

What *is* here, and fast: the composition builds three datasets under both
tie modes, and tying removes exactly one free parameter; the generator's two
spectrographs carry their own, different calibration truths (checked
against the noiseless prediction, not asserted); a tiny-budget fit runs in
every combination of ``tie``/``gp`` and returns exactly the free parameters
the problem itself declares; ``main`` prints the negotiated requirements and
a report, with and without ``--tie``/``--gp``.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import negotiate
from examples.photometry_spectra import generators
from examples.photometry_spectra.photometry_spectra import (
    FILTERS,
    LL_OBSERVED_WAVELENGTH,
    SL_OBSERVED_WAVELENGTH,
    build_instruments,
    build_model,
    build_problem,
    fit,
    main,
    qualified_truth,
    recovers_truth,
    report,
)

#: 20 walkers clears emcee's ``2 x n_dim`` floor for every combination this
#: module builds (at most nine free parameters, under ``tie=False, gp=True``).
TINY = {"walkers": 20, "steps": 20, "burn_in": 5}


class TestTheCompositionBuilds:
    """One model, three datasets, distinct labels, one negotiated channel."""

    @pytest.mark.parametrize("tie", [False, True])
    def test_build_problem_composes_three_datasets(self, tie: bool) -> None:
        problem = build_problem(tie=tie)
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"sl", "ll", "catalogue"}
        assert set(problem.parameters.free_names) == set(qualified_truth(tie=tie))

    def test_tying_removes_exactly_one_free_parameter(self) -> None:
        untied = build_problem(tie=False)
        tied = build_problem(tie=True)
        assert tied.free_size == untied.free_size - 1
        assert "calibration" in tied.parameters.free_names
        assert "sl.instrument.calibration_scale.scale" not in tied.parameters.free_names
        assert "ll.instrument.calibration_scale.scale" not in tied.parameters.free_names

    def test_gp_adds_four_free_parameters_per_tie_mode(self) -> None:
        for tie in (False, True):
            without_gp = build_problem(tie=tie, gp=False)
            with_gp = build_problem(tie=tie, gp=True)
            assert with_gp.free_size == without_gp.free_size + 4

    def test_build_model_and_instruments_are_reusable_alone(self) -> None:
        model = build_model()
        sl, ll, camera = build_instruments()
        assert sl.label == "sl"
        assert ll.label == "ll"
        assert camera.label == "catalogue"
        assert sl.channel == ll.channel == camera.channel == "sed"
        observed_sl, observed_ll, observed_photometry = generators.synthetic_data(
            model, sl, ll, camera, seed=generators.SEED
        )
        assert observed_sl.values.size == SL_OBSERVED_WAVELENGTH.size
        assert observed_ll.values.size == LL_OBSERVED_WAVELENGTH.size
        assert tuple(observed_photometry.filters) == FILTERS


class TestTheCalibrationFactorsAreInjected:
    """The two spectrographs' generated spectra carry their own, different truths."""

    def test_sl_near_10_micron_matches_its_own_calibration_truth(self) -> None:
        model = build_model()
        sl, ll, camera = build_instruments()
        observed_sl, _, _ = generators.synthetic_data(model, sl, ll, camera, seed=generators.SEED)

        # The same negotiate/compile_for dance generators.synthetic_data does
        # internally, on a fresh model and instruments, uncalibrated -- the
        # baseline the noisy, calibrated spectrum is checked against.
        baseline_model = build_model()
        baseline_sl, baseline_ll, baseline_camera = build_instruments()
        requirements = negotiate([baseline_sl, baseline_ll, baseline_camera])
        compiled = baseline_model.compile_for(requirements)
        truth = compiled(**generators.TRUTH)
        raw_sl = baseline_sl(truth, {"calibration_scale.scale": 1.0})

        wavelength = observed_sl.spectral_axis.values
        index = int(np.argmin(np.abs(wavelength - 10.0)))
        ratio = float(observed_sl.values[index] / raw_sl.values[index])
        assert ratio == pytest.approx(generators.CALIBRATION_TRUTH["sl"], abs=0.1)

    def test_the_two_truths_differ(self) -> None:
        assert generators.CALIBRATION_TRUTH["sl"] != generators.CALIBRATION_TRUTH["ll"]


class TestTheFitRuns:
    """A tiny-budget fit, end to end, returns exactly the declared free parameters."""

    @pytest.mark.parametrize("gp", [False, True])
    @pytest.mark.parametrize("tie", [False, True])
    def test_fit_returns_the_expected_parameters(self, tie: bool, gp: bool) -> None:
        problem = build_problem(tie=tie, gp=gp)
        run = fit(problem, **TINY)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(problem.parameters.free_names)
        assert posterior.sizes["chain"] == TINY["walkers"]
        assert posterior.sizes["draw"] == TINY["steps"] - TINY["burn_in"]

    def test_report_names_every_parameter(self) -> None:
        problem = build_problem()
        run = fit(problem, **TINY)
        text = report(run, tie=False)
        assert "emcee on reference" in text
        for name in problem.parameters.free_names:
            assert name in text

    def test_recovers_truth_reports_one_bool_per_qualified_parameter(self) -> None:
        problem = build_problem(tie=True)
        run = fit(problem, **TINY)
        covered = recovers_truth(run, tie=True)
        assert set(covered) == set(qualified_truth(tie=True))
        assert all(isinstance(value, bool) for value in covered.values())


class TestMain:
    """The CLI entry point, at the tiny budget."""

    def test_main_prints_requirements_and_a_report(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        assert (
            main(
                [
                    "--walkers",
                    str(TINY["walkers"]),
                    "--steps",
                    str(TINY["steps"]),
                    "--burn-in",
                    str(TINY["burn_in"]),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "negotiated channels:" in out
        assert "sources asking of channel 'sed':" in out
        assert "free parameters:" in out
        assert "emcee on reference" in out

    def test_main_with_tie_and_gp(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert (
            main(
                [
                    "--tie",
                    "--gp",
                    "--walkers",
                    str(TINY["walkers"]),
                    "--steps",
                    str(TINY["steps"]),
                    "--burn-in",
                    str(TINY["burn_in"]),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "calibration" in out
        assert "likelihood.amplitude" in out
