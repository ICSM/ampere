"""W5.1's pinned claim: a latent GP on closure phases survives a smooth phase error.

The route is W4.4's (:mod:`tests.interferometry.test_calibration`), with the
engine changed and the misspecification moved. The simulating problem is
W4.4's two-dataset binary without the disc — visibilities under independent
circular complex noise, closure phases under independent von Mises noise;
every replica then has a fixed, smooth per-triangle phase error of 0.6 rad
added to its closure phases
(:func:`~examples.interferometry.study.with_phase_error`) and is refitted
twice, :data:`~examples.interferometry.study.SIMULATIONS` times each:

``"rigid"``
    independent von Mises noise — the error is invisible to it, so it can only
    be absorbed by moving the binary;
``"latent"``
    W5.1's latent GP over the triangle axes
    (``Matern32(axes=("u1", "v1", "u2", "v2"))`` on ``DenseGP``), under NUTS
    on jax — the composition is fitted on the native path only.

The visibilities carry independent noise in both arms and the same engine
runs both, so the two differ in the closure phases' likelihood alone. The
budget (twelve simulations) is W4.4's: evidence of a gross effect, not a
uniformity test — ``sbc`` warns about exactly that, and the warning is
allowed through. What is asserted is the direction plus a margin:

- measured at seed 20260925 (``run_phase_calibration``'s default):
  ``latent`` coverage@0.9 = ``[1.00, 1.00]``, ``rigid`` = ``[0.42, 0.08]``,
  on ``(separation, flux_ratio)``; at seed 7, ``[0.92, 0.92]`` and
  ``[0.25, 0.33]``;
- the thresholds are :mod:`tests.interferometry.test_calibration`'s own
  floor, ceiling and gap, unchanged, so "covers" means the same thing for
  both modalities' flexible likelihoods. The tightest margin is the rigid
  separation's: 0.42 against a 0.55 ceiling, between one and two
  simulations.

**Budget.** The two calibrations take 206 s on the jax backend (measured:
``pytest -m interferometry_full`` on this file, 211 s wall clock), so the
rows are ``interferometry_full``-gated like W4.4's full-budget rows; the
default run keeps the structural rows only, which take well under a second.
"""

from __future__ import annotations

import importlib.util
import warnings
from typing import Any

import numpy as np
import pytest

from examples.interferometry import study

from .test_calibration import COVERAGE_GAP, MAX_COVERAGE, MIN_COVERAGE, NOMINAL

requires_jax = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None,
    reason="the closure-phase latent GP is fitted on the native path only; this arm runs on jax",
)


@pytest.fixture(scope="module")
def phase_calibrations() -> dict[str, Any]:
    """Both closure-phase arms' SBC datasets, computed once per module."""
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*goodness-of-fit test.*")
        return {arm: study.run_phase_calibration(arm) for arm in study.PHASE_ARMS}


class TestThePhaseArmComposes:
    """Structural, on the reference backend: cheap, and true before any coverage means anything."""

    def test_the_injected_error_is_smooth_and_several_sigma(self) -> None:
        import interferometry_fixtures as fixtures

        error = study.phase_error(fixtures.closure_phases())
        assert error.shape == (fixtures.closure_phases().n_samples,)
        assert np.max(np.abs(error)) > 3.0 * study.gen.SIGMA_CLOSURE
        assert np.max(np.abs(error)) <= study.PHASE_ERROR_AMPLITUDE

    def test_the_injection_wraps(self) -> None:
        import interferometry_fixtures as fixtures

        near_cut = fixtures.closure_phases(np.full(fixtures.closure_phases().n_samples, 3.1))
        shifted = np.asarray(study.with_phase_error(near_cut).values)
        assert np.all(shifted > -np.pi) and np.all(shifted <= np.pi)

    def test_the_latent_arm_declares_a_latent_block_and_the_rigid_arm_none(self) -> None:
        import interferometry_fixtures as fixtures

        from ampere.core import Marginalisation

        observed = fixtures.closure_phases()
        latent = study.phase_problem("reference", "latent", observed).datasets["t3"]
        rigid = study.phase_problem("reference", "rigid", observed).datasets["t3"]
        assert latent.likelihood.marginalisation is Marginalisation.LATENT
        assert latent.latent is not None and latent.latent.size == observed.n_samples
        assert rigid.latent is None

    def test_no_gradient_free_engine_takes_the_latent_arm(self) -> None:
        """Native path only: emcee is refused by name on the reference backend."""
        import interferometry_fixtures as fixtures

        from ampere.core.exceptions import ContractError
        from ampere.inference import EmceeEngine

        problem = study.phase_problem("reference", "latent", fixtures.closure_phases())
        with pytest.raises(ContractError, match="emcee cannot run this likelihood"):
            EmceeEngine(problem, walkers=8)


@requires_jax
@pytest.mark.interferometry_full
class TestPhaseCalibration:
    """The pinned coverage claim, over :data:`~examples.interferometry.study.SIMULATIONS` refits."""

    def test_the_latent_arm_covers(self, phase_calibrations: dict[str, Any]) -> None:
        coverage = study.coverage_at(phase_calibrations["latent"], NOMINAL)
        print(f"\nlatent coverage at {NOMINAL}: {coverage}")
        assert np.all(coverage >= MIN_COVERAGE)

    def test_the_rigid_arm_does_not(self, phase_calibrations: dict[str, Any]) -> None:
        coverage = study.coverage_at(phase_calibrations["rigid"], NOMINAL)
        print(f"\nrigid coverage at {NOMINAL}: {coverage}")
        assert np.all(coverage <= MAX_COVERAGE)

    def test_the_gap_between_latent_and_rigid_is_wide(
        self, phase_calibrations: dict[str, Any]
    ) -> None:
        latent = study.coverage_at(phase_calibrations["latent"], NOMINAL)
        rigid = study.coverage_at(phase_calibrations["rigid"], NOMINAL)
        gap = float(np.min(latent - rigid))
        print(f"\ncoverage gap (latent - rigid) at {NOMINAL}: {gap:.3f}")
        assert gap >= COVERAGE_GAP

    def test_every_simulation_produced_a_usable_fit(
        self, phase_calibrations: dict[str, Any]
    ) -> None:
        for arm, calibration in phase_calibrations.items():
            failures = int(calibration.attrs["ampere_calibration_failures"])
            assert failures == 0, f"arm {arm!r} had {failures} failed simulation(s)."
            assert calibration["ranks"].shape == (study.SIMULATIONS, 2)
