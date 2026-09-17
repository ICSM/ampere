"""Shared collection hooks for ``tests/results``."""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``population_full`` rows unless ``-m population_full`` asked for them (W5.22).

    ``tests/m2/conftest.py``'s hook, restated for this directory:
    ``tests/results/test_population.py``'s 200-object archive (fitting and
    timing 200 tiny ``emcee`` runs) is minutes of wall clock a per-PR gate
    should not pay for by default, so it is opt-in. A local hook rather than
    an ``addopts`` deselection in ``pyproject.toml``, so it says why it
    skipped in the report and leaves every other suite's invocation alone.
    """
    selected = config.getoption("-m", default="") or ""
    if "population_full" in selected:
        return
    skip = pytest.mark.skip(
        reason="the population study's full 200-object archive: run with `pytest -m population_full`"
    )
    for item in items:
        if "population_full" in item.keywords:
            item.add_marker(skip)
