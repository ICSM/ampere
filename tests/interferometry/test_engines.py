"""The study "under every engine": NUTS on torch, NPE through :class:`~ampere.inference.SBIEngine`.

W4.3 already proved the *modality* runs under every engine
(``tests/inference/test_interferometry.py``); this file's job is narrower —
that the study's own two-dataset, binary-plus-disc problem (not the plain
binary those tests use) composes and fits under each, at a smoke budget. The
pinned science stays on the reference backend
(``tests/interferometry/test_calibration.py``), for
``tests/m2/conftest.py``'s own reason: the claim is about the likelihood, not
about which array library ran it.

Both classes below are skipped by name when their extra is not installed, so
this file collects (and mostly skips) under plain ``-e dev``.
"""

from __future__ import annotations

import importlib.util
from typing import Any

import numpy as np
import pytest

from examples.interferometry import study

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_sbi = pytest.mark.skipif(importlib.util.find_spec("sbi") is None, reason="needs ampere[sbi]")

#: NUTS's smoke budget: a handful of draws, enough to say the realised density
#: is finite and the chain moves, not enough to say anything about coverage.
NUTS_DRAWS = 60
NUTS_WARMUP = 60


@needs_torch
class TestNutsOnTorch:
    """The "flexible" arm (the study's most demanding composition) under NUTS."""

    @pytest.fixture(scope="class")
    def run(self) -> Any:
        from ampere.inference import NUTSEngine

        problem = study.build_problem("torch", "flexible")
        return NUTSEngine(problem).run(NUTS_DRAWS, warmup=NUTS_WARMUP, chains=1, progress=False)

    def test_the_composite_and_the_plain_binary_both_lower_natively(self, run: Any) -> None:
        posterior = run["posterior"]
        for name in study.gen.TRUTH:
            values = np.asarray(posterior[name].values, dtype=float)
            assert values.size == NUTS_DRAWS
            assert np.all(np.isfinite(values))

    def test_the_correct_arms_composite_also_lowers_natively(self) -> None:
        """``BinaryWithDisc``'s ``native_grid``/``native_flux`` surface, exercised."""
        from ampere.inference import NUTSEngine

        problem = study.build_problem("torch", "correct")
        run = NUTSEngine(problem).run(NUTS_DRAWS, warmup=NUTS_WARMUP, chains=1, progress=False)
        posterior = run["posterior"]
        for name in study.gen.TRUTH:
            values = np.asarray(posterior[name].values, dtype=float)
            assert np.all(np.isfinite(values))


#: The smoke budget ``tests/inference/test_interferometry.py`` uses for this
#: same modality's own NPE row, unchanged: large enough that the density
#: estimator is trained rather than initialised, small enough for a per-PR gate.
SBI_BUDGET = 1200
SBI_EPOCHS = 300
SBI_DRAWS = 200
CALIBRATION_COUNT = 100
CALIBRATION_DRAWS = 100


@needs_sbi
class TestNpeOnTheFlexibleArm:
    """W4.2's guidance for this row: TARP and the coverage curve, not the marginal KS p-value."""

    @pytest.fixture(scope="class")
    def engine(self) -> Any:
        from ampere.inference import SBIEngine

        problem = study.build_problem("reference", "flexible")
        engine = SBIEngine(problem, method="npe", budget=SBI_BUDGET, layout="set", embedding="set")
        engine.run(draws=SBI_DRAWS, training={"max_num_epochs": SBI_EPOCHS})
        return engine

    @pytest.fixture(scope="class")
    def calibration(self, engine: Any) -> Any:
        return engine.calibrate(count=CALIBRATION_COUNT, posterior_draws=CALIBRATION_DRAWS)

    def test_the_fit_runs_under_the_set_layout(self, engine: Any) -> None:
        assert engine.encoding.kind == "set"
        assert set(engine.encoding.labels) == {"vis", "t3"}
        assert engine.encoding.complex_columns is True

    def test_the_trained_posterior_passes_the_tarp_coverage_check(self, calibration: Any) -> None:
        assert abs(float(calibration.attrs["ampere_calibration_tarp_atc"])) < 0.15
        levels = np.asarray(calibration.coords["level"].values)
        curve = np.asarray(calibration["coverage"].values)
        for index in range(curve.shape[1]):
            assert abs(float(np.interp(0.68, levels, curve[:, index])) - 0.68) < 0.2
