"""``examples/sbi/npe_native.py`` runs, and fits a native problem in-process (W4.0 (5)).

Carried from W3.13's finding: no ``examples/sbi`` script did a bare NPE fit of
a native (torch/jax) problem for its own sake — every one of the five either
wrapped a subprocess simulator or stayed on the reference backend. This is
that script's coverage, matching ``tests/examples/test_cached_fit.py``'s
pattern: the directory goes on ``sys.path`` and the module is imported by
name, because a module loaded from a file path under a synthetic name cannot
be cleanly re-imported (``tests/examples/test_external_simulator.py``'s
reasoning, restated there in full).
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys
from collections.abc import Iterator
from types import ModuleType

import pytest

HAS_SBI = importlib.util.find_spec("sbi") is not None

needs_sbi = pytest.mark.skipif(
    not HAS_SBI,
    reason="needs the 'sbi' extra (pixi run -e sbi ...), which brings in torch too",
)

EXAMPLES = pathlib.Path(__file__).resolve().parents[2] / "examples" / "sbi"


@pytest.fixture(scope="module")
def example() -> Iterator[ModuleType]:
    """The imported example module."""
    path = str(EXAMPLES)
    inserted = path not in sys.path
    if inserted:
        sys.path.insert(0, path)
    try:
        import npe_native

        yield npe_native
    finally:
        if inserted and path in sys.path:
            sys.path.remove(path)


@needs_sbi
class TestTheExampleFitsANativeProblem:
    """The problem really is native, and the engine fits it with nothing else going on."""

    def test_the_problem_is_torch_and_batchable(self, example: ModuleType) -> None:
        problem = example.build_problem()
        assert problem.backend == "torch"
        assert problem.batchable is True

    def test_it_fits_and_reports(self, example: ModuleType) -> None:
        problem = example.build_problem()
        run = example.fit(problem, budget=40, draws=20)
        report = example.report(run)
        assert "npe on torch" in report
        assert "posterior (truth in brackets):" in report
        assert "model.norm" in report and "model.index" in report
        assert run.attrs["ampere_backend"] == "torch"
        assert run.attrs["ampere_sbi_simulations"] == 40
        assert run["posterior"].dataset.sizes["draw"] == 20

    def test_main_reports(self, example: ModuleType, capsys: pytest.CaptureFixture[str]) -> None:
        assert example.main(["--budget", "40", "--draws", "20"]) == 0
        out = capsys.readouterr().out
        assert "npe on torch" in out
        assert "wall clock" in out
