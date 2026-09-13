"""``examples/astrometry`` builds, negotiates, and fits end to end (W4.9).

Full-budget recovery is **not** part of this suite, for
``tests/examples/test_sed_composition.py``'s own reason: it is one to two
orders of magnitude slower than everything else here. It is verified by
running ``python -m examples.astrometry`` directly (the branch report quotes
the recovered intervals) and by ``tests/astrometry``'s own pinned rows at a
budget larger than this smoke test's.

What *is* here, and fast: the composition builds and negotiates two
channels from one model; a tiny-budget fit runs end to end on the reference
backend and returns the six expected parameters; ``main`` prints the
negotiated requirements and a report.
"""

from __future__ import annotations

import numpy as np
import pytest

from examples.astrometry import generators
from examples.astrometry.astrometry import (
    QUALIFIED_TRUTH,
    build_instruments,
    build_model,
    build_problem,
    fit,
    main,
    recovers_truth,
    report,
)

#: emcee needs at least 2 x n_dim walkers to span this six-parameter problem.
TINY = {"walkers": 16, "steps": 20, "burn_in": 5}


class TestTheCompositionBuilds:
    """One model, two channels, one negotiated epoch grid per channel."""

    def test_build_problem_composes_two_channels(self) -> None:
        problem = build_problem("reference")
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"ra", "dec"}
        assert set(problem.parameters.free_names) == set(QUALIFIED_TRUTH)

    def test_each_channel_is_negotiated_from_its_own_instrument(self) -> None:
        problem = build_problem("reference")
        for channel in ("ra", "dec"):
            requirements = problem.requirements["model"][channel]
            assert "time" in requirements.axes
            axis = requirements.axes["time"]
            assert axis.points is not None
            assert axis.points.size == generators.EPOCHS.size

    def test_build_model_and_instruments_are_reusable_alone(self) -> None:
        model = build_model("reference")
        ra_instrument, dec_instrument = build_instruments("reference")
        assert ra_instrument.label == "astrom_ra"
        assert dec_instrument.label == "astrom_dec"
        assert ra_instrument.channel == "ra"
        assert dec_instrument.channel == "dec"
        observed_ra, observed_dec = generators.synthetic_data(
            model, ra_instrument, dec_instrument, seed=generators.SEED
        )
        assert observed_ra.values.size == generators.EPOCHS.size
        assert observed_dec.values.size == generators.EPOCHS.size

    def test_the_gp_arm_composes_too(self) -> None:
        problem = build_problem("reference", gp=True)
        assert set(problem.parameters.free_names) >= set(QUALIFIED_TRUTH)
        assert np.isfinite(
            problem.log_prob(
                {
                    **{name: value for name, value in QUALIFIED_TRUTH.items()},
                    "ra.likelihood.amplitude": 0.02,
                    "ra.likelihood.length_scale": 100.0,
                    "dec.likelihood.amplitude": 0.02,
                    "dec.likelihood.length_scale": 100.0,
                }
            )
        )


class TestTheFitRuns:
    """A tiny-budget fit, end to end, returns the six expected parameters."""

    def test_fit_returns_the_expected_parameters(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == set(QUALIFIED_TRUTH)
        assert posterior.sizes["chain"] == TINY["walkers"]
        assert posterior.sizes["draw"] == TINY["steps"] - TINY["burn_in"]

    def test_report_names_every_parameter(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        text = report(run)
        assert "emcee on reference" in text
        for name in QUALIFIED_TRUTH:
            assert name in text

    def test_recovers_truth_reports_one_bool_per_parameter(self) -> None:
        problem = build_problem("reference")
        run = fit(problem, backend="reference", **TINY)
        covered = recovers_truth(run)
        assert set(covered) == set(QUALIFIED_TRUTH)
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
        assert "sources asking of channel 'ra':" in out
        assert "sources asking of channel 'dec':" in out
        assert "free parameters:" in out
        assert "emcee on reference" in out
        assert "wall clock" in out
