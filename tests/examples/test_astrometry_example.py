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
    SBC_PARAMETERS,
    SBC_SHARED_PARAMETER,
    build_instruments,
    build_model,
    build_problem,
    calibrate,
    coverage_at,
    fit,
    main,
    recovers_truth,
    report,
    sbc_problem,
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


# ---------------------------------------------------------------------------
# The joint arm (W5.9)
# ---------------------------------------------------------------------------

#: The tiny budget the always-on calibration rows use. Twelve refits is a
#: smoke test of the *machinery*, not a calibration result, and the row that
#: uses it says so: `ampere.results.sbc` itself warns below a hundred.
TINY_SBC = {"count": 12, "draws": 120, "walkers": 12, "steps": 250, "burn_in": 100}

#: The budget the pinned claim is measured at. Two arms x 48 refits of a
#: five-parameter problem is well over ten minutes on one core, which is why
#: this row carries the ``astrometry_full`` marker and the row above does not.
FULL_SBC = {"count": 48, "draws": 200}

#: Pinned margins, measured on the reference backend at
#: ``examples.astrometry.generators.SEED`` (the run the branch report quotes):
#: the joint arm covers 1.000 and the independent arm 0.646 at the 0.90 level
#: on ``model.phase``. The margin below nominal is generous because 48
#: simulations give a standard error of about 0.043 on a coverage of 0.9, and
#: the *gap* is what the claim is about --- it is 0.35 as measured and the row
#: asks for a fifth of that.
JOINT_COVERAGE_FLOOR = 0.80
COVERAGE_GAP = 0.07


class TestTheJointArm:
    """One correlated process over both channels, on injected correlated data."""

    def test_the_joint_arm_composes_and_scores(self) -> None:
        problem = build_problem("reference", joint=True)
        assert set(problem.datasets) == {"ra", "dec"}
        # The group is one component of the parameter space...
        assert set(problem.parameters.free_names) >= {
            "astrom.angle",
            "astrom.log_variance_0",
            "astrom.log_variance_1",
        }
        # ... and one contribution, in place of its members'.
        assert problem.datasets.contribution_labels() == ("astrom",)
        assert problem.datasets.group_of("ra") == "astrom"
        theta = {
            **QUALIFIED_TRUTH,
            "astrom.angle": generators.JOINT_TRUTH["angle"],
            "astrom.log_variance_0": generators.JOINT_TRUTH["log_variance_0"],
            "astrom.log_variance_1": generators.JOINT_TRUTH["log_variance_1"],
        }
        evaluation = problem.evaluate(theta)
        assert set(evaluation.contributions) == {"astrom"}
        assert np.isfinite(evaluation.log_likelihood)

    def test_the_injected_systematic_is_correlated_between_the_channels(self) -> None:
        """The data the study fits really do carry a cross-channel error."""
        model = build_model("reference")
        ra_instrument, dec_instrument = build_instruments("reference")
        plain = generators.synthetic_data(model, ra_instrument, dec_instrument)
        injected = generators.synthetic_joint_data(model, ra_instrument, dec_instrument)
        differences = [
            np.asarray(injected[index].values, dtype=float)
            - np.asarray(plain[index].values, dtype=float)
            for index in (0, 1)
        ]
        # The difference between the two generators is exactly the injected
        # systematic (they share a white-noise stream), so it is a draw from
        # B (x) K_x and nothing else.
        assert np.max(np.abs(differences[0])) > generators.SIGMA
        correlation = float(np.corrcoef(differences[0], differences[1])[0, 1])
        implied = generators.coupling_matrix()
        expected = implied[0, 1] / np.sqrt(implied[0, 0] * implied[1, 1])
        assert correlation == pytest.approx(expected, abs=0.25)
        assert correlation > 0.5

    def test_the_joint_fit_recovers_the_orbit_at_a_small_budget(self) -> None:
        """A short chain on the joint arm runs and returns the six parameters."""
        problem = build_problem("reference", joint=True)
        # Nine dimensions here, not six: emcee needs 2 x n_dim walkers, and the
        # group's three parameters are three of them.
        run = fit(problem, backend="reference", **{**TINY, "walkers": 24})
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) >= set(QUALIFIED_TRUTH)
        assert "astrom.angle" in posterior.data_vars

    def test_simulate_draws_the_channels_together(self) -> None:
        problem = build_problem("reference", joint=True)
        theta = {
            **QUALIFIED_TRUTH,
            "astrom.angle": generators.JOINT_TRUTH["angle"],
            "astrom.log_variance_0": generators.JOINT_TRUTH["log_variance_0"],
            "astrom.log_variance_1": generators.JOINT_TRUTH["log_variance_1"],
        }
        simulation = problem.simulate(theta, observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        assert set(simulation.observations) == {"ra", "dec"}


class TestTheCalibrationStudy:
    """SBC on each arm, against the joint generator (*W5.9*)."""

    def test_the_study_runs_at_a_smoke_budget(self) -> None:
        """The machinery, not the claim: twelve refits is not a calibration result."""
        calibration = calibrate("reference", arm="joint", seed=generators.SEED, **TINY_SBC)
        assert calibration["ranks"].shape == (TINY_SBC["count"], len(SBC_PARAMETERS))
        coverage = coverage_at(calibration, 0.9, SBC_SHARED_PARAMETER)
        assert 0.0 <= coverage <= 1.0

    def test_the_independent_arm_refits_the_same_data(self) -> None:
        """The comparison arm is a different *model* on the simulating arm's data."""
        joint = sbc_problem("reference", arm="joint")
        independent = sbc_problem(
            "reference",
            arm="independent",
            observed=(joint.datasets["ra"].observed, joint.datasets["dec"].observed),
        )
        assert independent.datasets.contribution_labels() == ("ra", "dec")
        assert independent.datasets.group_of("ra") is None
        assert set(independent.parameters.free_names) >= {
            "ra.likelihood.amplitude",
            "dec.likelihood.amplitude",
        }
        for label in ("ra", "dec"):
            assert np.array_equal(
                np.asarray(independent.datasets[label].observed.values),
                np.asarray(joint.datasets[label].observed.values),
            )

    @pytest.mark.astrometry_full
    def test_the_joint_arm_is_calibrated_where_the_independent_one_is_not(self) -> None:
        """The pinned claim, at the count where SBC has power.

        Under an injected centroiding systematic shared by the two sky axes,
        the joint fit's interval on the **shared** orbital phase covers at its
        nominal rate and the independent-GP fit's does not. Two independent
        GPs reproduce each axis's marginal scatter exactly and can say nothing
        about the correlation between them, so they combine two error-laden
        estimates of the phase as though the errors were independent --- and
        the combined interval comes out narrower than the truth's own scatter.
        """
        joint = calibrate("reference", arm="joint", seed=generators.SEED, **FULL_SBC)
        independent = calibrate("reference", arm="independent", seed=generators.SEED, **FULL_SBC)
        joint_coverage = coverage_at(joint, 0.9, SBC_SHARED_PARAMETER)
        independent_coverage = coverage_at(independent, 0.9, SBC_SHARED_PARAMETER)
        assert joint_coverage >= JOINT_COVERAGE_FLOOR, (
            f"the joint arm covered {joint_coverage:.3f} at the 0.90 level, below the pinned "
            f"floor of {JOINT_COVERAGE_FLOOR}."
        )
        assert independent_coverage <= joint_coverage - COVERAGE_GAP, (
            f"the independent arm covered {independent_coverage:.3f} against the joint arm's "
            f"{joint_coverage:.3f}: a gap of {joint_coverage - independent_coverage:.3f}, below "
            f"the pinned {COVERAGE_GAP}."
        )
