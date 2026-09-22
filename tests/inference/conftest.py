"""Shared collection rules for ``tests/inference``.

One job today: keep the engine battery's full-budget calibration rows
(W5.14, ``tests/inference/test_nested.py``) out of the per-PR gate unless
they are asked for by name.
"""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``engines_full`` rows unless ``-m engines_full`` asked for them.

    The full budget is 100 refits *per engine* — Talts et al.'s uniformity
    test only has power at that count, and at that count it is minutes per
    engine rather than seconds. The reduced sibling that always runs checks
    that the machinery produces well-shaped ranks with no failed simulation,
    which is what a per-PR gate can afford and what catches a driver that
    cannot be refitted.

    Skipping here rather than with a global ``addopts = "-m 'not
    engines_full'"`` in ``pyproject.toml`` is deliberate, for the three
    reasons ``tests/m2/conftest.py``'s identical hook gives: an ``addopts``
    deselection is invisible from the command line, applies to every suite in
    the repository, and silently conflicts with any other ``-m`` a caller
    passes.
    """
    selected = config.getoption("-m", default="") or ""
    if "engines_full" in selected:
        return
    skip = pytest.mark.skip(
        reason="the engine battery's full SBC budget: run with `pytest -m engines_full`"
    )
    for item in items:
        if "engines_full" in item.keywords:
            item.add_marker(skip)
