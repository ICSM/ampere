"""The item's pinned claim: coverage under "correct" and "flexible", failure under "incomplete".

The route is simulation-based calibration
(:func:`~ampere.results.calibration.sbc`), in
``tests/m2/test_visibility_calibration.py``'s own shape and for its own
reason: a single noisy draw's posterior either happens to cover the truth or
does not, and that is a statement about one seed, not about a likelihood. SBC
draws a fresh binary from the same prior the fitted models use, simulates a
fresh two-dataset observation from the **correct** (disc-included) truth, and
refits it under one arm's formulation, :data:`~examples.interferometry.study.SIMULATIONS`
times — so "arm X covers" is a statement about coverage over many draws.

The budget (:data:`~examples.interferometry.study.SIMULATIONS` = 12) is
small — Talts et al. (2018) want a hundred or more for the uniformity test to
have power — so, as that module's docstring states for its own row, the
numbers here are evidence of a gross effect and nothing finer; ``sbc`` warns
about exactly that (``tests/interferometry/conftest.py`` allows the warning
through). What is asserted is the **direction** plus a margin, not an exact
figure:

- measured at the pinned seed (``examples.interferometry.study.gen.SEED``):
  ``correct`` coverage@0.9 = ``[1.00, 1.00]``, ``incomplete`` = ``[0.25, 0.00]``,
  ``flexible`` = ``[0.92, 0.92]``, on ``(separation, flux_ratio)``.
- the thresholds below leave several simulations of slack in every
  direction (one simulation is worth 1/12 = 0.083 of coverage) before the
  row would fail.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from examples.interferometry import study

NOMINAL = 0.9

#: "correct" and "flexible" must cover at least this often at the nominal
#: level — well below what was measured (1.00 and 0.92), several
#: simulations of slack.
MIN_COVERAGE = 0.55
#: "incomplete" must cover no more often than this — well above what was
#: measured (0.25 and 0.00).
MAX_COVERAGE = 0.55
#: The gap between "flexible" and "incomplete" must be at least this wide —
#: measured 0.67 and 0.92.
COVERAGE_GAP = 0.3


class TestTheModalityComposes:
    """Structural: cheap, and it has to hold before the coverage numbers mean anything."""

    def test_the_two_datasets_share_one_negotiation(self) -> None:
        problem = study.build_problem("reference", "correct")
        assert set(problem.datasets) == {"vis", "t3"}
        assert problem.requirements["model"]["sky"].sources == ("vis", "t3")

    def test_the_flexible_arm_puts_the_gp_on_visibilities_only(self) -> None:
        from ampere.core import Marginalisation

        problem = study.build_problem("reference", "flexible")
        vis = problem.datasets["vis"]
        t3 = problem.datasets["t3"]
        assert vis.likelihood.marginalisation is Marginalisation.ANALYTIC
        assert vis.latent is None
        assert type(t3.likelihood.noise).__name__ == "IndependentNoise"

    def test_the_incomplete_and_correct_arms_differ_only_by_the_disc(self) -> None:
        correct = study.build_problem("reference", "correct")
        incomplete = study.build_problem("reference", "incomplete")
        assert (
            correct.parameters.free_names
            == incomplete.parameters.free_names
            == (
                "model.separation",
                "model.flux_ratio",
            )
        )


class TestCalibration:
    """The pinned coverage claim, over :data:`~examples.interferometry.study.SIMULATIONS` refits."""

    def test_the_correct_arm_is_calibrated(self, calibrations: dict[str, Any]) -> None:
        coverage = study.coverage_at(calibrations["correct"], NOMINAL)
        print(f"\ncorrect coverage at {NOMINAL}: {coverage}")
        assert np.all(coverage >= MIN_COVERAGE)

    def test_the_incomplete_arm_fails(self, calibrations: dict[str, Any]) -> None:
        coverage = study.coverage_at(calibrations["incomplete"], NOMINAL)
        print(f"\nincomplete coverage at {NOMINAL}: {coverage}")
        assert np.all(coverage <= MAX_COVERAGE)

    def test_the_flexible_arm_recovers(self, calibrations: dict[str, Any]) -> None:
        coverage = study.coverage_at(calibrations["flexible"], NOMINAL)
        print(f"\nflexible coverage at {NOMINAL}: {coverage}")
        assert np.all(coverage >= MIN_COVERAGE)

    def test_the_gap_between_flexible_and_incomplete_is_wide(
        self, calibrations: dict[str, Any]
    ) -> None:
        flexible = study.coverage_at(calibrations["flexible"], NOMINAL)
        incomplete = study.coverage_at(calibrations["incomplete"], NOMINAL)
        gap = float(np.min(flexible - incomplete))
        print(f"\ncoverage gap (flexible - incomplete) at {NOMINAL}: {gap:.3f}")
        assert gap >= COVERAGE_GAP

    def test_every_simulation_produced_a_usable_fit(self, calibrations: dict[str, Any]) -> None:
        for arm, calibration in calibrations.items():
            failures = int(calibration.attrs["ampere_calibration_failures"])
            assert failures == 0, f"arm {arm!r} had {failures} failed simulation(s)."
            assert calibration["ranks"].shape == (study.SIMULATIONS, 2)
