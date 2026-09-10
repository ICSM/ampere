"""``examples/sbi/cached_fit.py`` runs, and a second identical run is a hit (W3.13).

The example demonstrates :class:`~ampere.results.ArtefactStore` /
:func:`~ampere.results.artefact_key` directly, one level below where
``SBIEngine(cache=...)`` calls the same seam (it predates W3.3's encoding
layout, which is why it builds its own small training routine rather than
using the engine — see its own docstring). It had no test coverage before
this item; the docs tutorial (``docs/source/sbi.rst``) quotes real output
from running it, so this is also what keeps that quoted output honest.

The directory goes on ``sys.path`` and the module is imported by name rather
than by path, matching ``tests/examples/test_external_simulator.py``'s
reasoning: a module loaded from a file path under a synthetic name cannot be
cleanly re-imported, and importing by name is what every other example test
in this directory does.
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pytest

HAS_SBI = importlib.util.find_spec("sbi") is not None

needs_sbi = pytest.mark.skipif(
    not HAS_SBI,
    reason="needs the 'sbi' extra (pixi run -e sbi ...)",
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
        import cached_fit

        yield cached_fit
    finally:
        if inserted and path in sys.path:
            sys.path.remove(path)


@needs_sbi
class TestTheExampleCaches:
    """A first run trains; an identical second run is a hit; a changed setting is a named miss."""

    def test_a_first_run_trains_and_writes(self, example: ModuleType, tmp_path: Path) -> None:
        posterior, hit, diff = example.fit_cached(budget=40, epochs=5, cache_dir=tmp_path)
        assert hit is False
        assert diff == {}
        draws = posterior.sample((3,), show_progress_bars=False)
        assert draws.reshape(-1).numel() == 3

    def test_a_second_identical_run_is_a_hit(self, example: ModuleType, tmp_path: Path) -> None:
        example.fit_cached(budget=40, epochs=5, cache_dir=tmp_path)
        posterior, hit, diff = example.fit_cached(budget=40, epochs=5, cache_dir=tmp_path)
        assert hit is True
        assert diff == {}
        # The round trip actually happened: the cached posterior still samples.
        draws = posterior.sample((3,), show_progress_bars=False)
        assert draws.reshape(-1).numel() == 3

    def test_a_changed_budget_is_a_miss_that_names_it(
        self, example: ModuleType, tmp_path: Path
    ) -> None:
        example.fit_cached(budget=40, epochs=5, cache_dir=tmp_path)
        _, hit, diff = example.fit_cached(budget=80, epochs=5, cache_dir=tmp_path)
        assert hit is False
        assert diff == {"budget": (40, 80)}

    def test_main_reports_hit_or_miss(
        self, example: ModuleType, tmp_path: Path, capsys: pytest.CaptureFixture[str]
    ) -> None:
        args = ["--budget", "40", "--epochs", "5", "--cache-dir", str(tmp_path)]
        assert example.main(args) == 0
        assert "cache miss" in capsys.readouterr().out
        assert example.main(args) == 0
        assert "cache hit" in capsys.readouterr().out
