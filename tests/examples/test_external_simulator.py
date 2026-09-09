"""``examples/sbi/external_simulator.py`` runs, under the process pool (W3.1).

The example is the *canonical* SBI simulator — compiled code behind a Python
call — so it is run here rather than merely imported: a subprocess wrapper that
does not survive being pickled to a worker, or whose crashes are not flagged, is
a wrapper that would fail on the first real budget.

The example directory goes on ``sys.path`` and the module is imported by name
rather than loaded from its path. That is not tidiness: a worker process must be
able to *re-import* the module to unpickle the model, and a module loaded under
a synthetic name from a file path cannot be re-imported anywhere.
"""

from __future__ import annotations

import pathlib
import sys
from collections.abc import Iterator
from types import ModuleType

import numpy as np
import pytest

from ampere.core import FailureReason, ProcessExecutor

EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / "examples" / "sbi"


@pytest.fixture(scope="module")
def example() -> Iterator[ModuleType]:
    """The imported example, importable by name in a worker for the length of the module."""
    path = str(EXAMPLES)
    inserted = path not in sys.path
    if inserted:
        sys.path.insert(0, path)
    try:
        import external_simulator

        yield external_simulator
    finally:
        if inserted and path in sys.path:
            sys.path.remove(path)


class TestTheExampleRuns:
    def test_a_budget_runs_under_the_process_pool(self, example: ModuleType) -> None:
        problem, batch = example.run(8, 2)
        assert len(batch) == 8
        assert len(batch.usable) >= 1
        assert batch.theta.shape == (8, 2)
        assert batch.observations is not None
        assert batch.observations["sed"].shape == (8, example.GRID.size)
        assert int(batch.failed.sum()) == sum(problem.failure_counts.values())

    def test_the_pooled_budget_is_the_serial_one(self, example: ModuleType) -> None:
        """The executor is an implementation detail of *how*, never of *what*."""
        serial = example.build_problem()
        pooled = example.build_problem()
        one = serial.simulate_many(6, observe=True)
        other = pooled.simulate_many(6, observe=True, executor=ProcessExecutor(2), chunk_size=4)
        assert np.array_equal(one.theta, other.theta)
        for left, right in zip(one, other, strict=True):
            assert left.failed == right.failed
            if left.failed:
                continue
            assert np.array_equal(
                left.observations["sed"].values, right.observations["sed"].values
            )

    def test_a_simulator_crash_is_flagged_with_its_own_diagnostic(
        self, example: ModuleType
    ) -> None:
        """A non-zero exit reaches the failure record as the simulator's own stderr."""
        problem = example.build_problem()
        table = np.array([[-1.0, 2.0], [-1.0, 9.9], [-1.0, 3.0]])
        batch = problem.simulate_many(3, values=table, executor=ProcessExecutor(2))
        assert batch.failed.tolist() == [False, True, False]
        failure = batch.failures[1]
        assert failure is not None
        assert failure.reason is FailureReason.MODEL_FAILED
        assert failure.exception_type == "SimulatorFailed"
        assert "failed to converge" in failure.message
        assert "exited 2" in failure.message
        assert len(batch.usable) == 2

    def test_a_draw_that_overruns_is_killed_and_counted_separately(
        self, example: ModuleType
    ) -> None:
        problem = example.build_problem(sleep=30.0)
        table = np.array([[-1.0, 2.0], [-2.4, 2.0], [-1.0, 3.0]])
        batch = problem.simulate_many(
            3, values=table, executor=ProcessExecutor(2, timeout=3.0)
        )
        assert batch.failed.tolist() == [False, True, False]
        failure = batch.failures[1]
        assert failure is not None
        assert failure.reason is FailureReason.EXECUTION_FAILED
        assert problem.failure_counts[FailureReason.EXECUTION_FAILED] == 1
        assert len(batch.usable) == 2

    def test_the_prediction_is_the_power_law_the_simulator_computed(
        self, example: ModuleType
    ) -> None:
        problem = example.build_problem()
        drawn = problem.simulate({"model.index": -1.0, "model.norm": 2.0})
        assert not drawn.failed
        assert np.allclose(
            drawn.predicted["sed"].values, 2.0 * example.GRID**-1.0, rtol=1e-12
        )

    def test_the_working_directory_is_per_process_and_made_once(
        self, example: ModuleType
    ) -> None:
        first = example.working_directory()
        assert first.is_dir()
        assert example.working_directory() == first

    def test_the_report_names_the_usable_pairs_and_the_reasons(
        self, example: ModuleType
    ) -> None:
        problem, batch = example.run(6, 2)
        text = example.report(problem, batch)
        assert "usable" in text and "theta stacked as (6, 2)" in text

    def test_main_runs_the_whole_example_end_to_end(
        self, example: ModuleType, capsys: pytest.CaptureFixture[str]
    ) -> None:
        assert example.main(["6", "2"]) == 0
        assert "draw(s) simulated" in capsys.readouterr().out
