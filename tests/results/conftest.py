"""Shared collection hooks for ``tests/results``."""

from __future__ import annotations

import importlib.util

import pytest


def _skip_study_rows_outside_dev(items: list[pytest.Item]) -> None:
    """Skip ``study``-marked rows wherever a modern backend is installed (W5.26).

    ``tests/m2/conftest.py``'s hook, restated for this directory:
    ``tests/results/test_population.py``'s reweighting rows sample a
    reference-backend ``emcee`` archive (``ConstantModel``/``TwoParameterModel``,
    neither backend-specific) and assert claims about
    :func:`~ampere.results.fit_population`, not about the array library, so
    running them again in ``torch``, ``jax`` or ``sbi`` buys no additional
    evidence. ``TestTheApproximateEngineRow`` is deliberately left unmarked:
    it is the one row in this module that *needs* a variational backend
    (``requires_vi``) and already skips itself in ``dev`` for that reason, so
    marking it ``study`` too would make it skip everywhere and never run.
    """
    reasons = [name for name in ("torch", "jax") if importlib.util.find_spec(name) is not None]
    if not reasons:
        return
    skip = pytest.mark.skip(
        reason=(
            "study row: numpy-path-only, redundant with the backend-agreement rows "
            f"(run on `dev` instead); this environment has {' and '.join(reasons)} installed"
        )
    )
    for item in items:
        if "study" in item.keywords:
            item.add_marker(skip)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``population_full`` rows unless ``-m population_full`` asked for them
    (W5.22), and the ``study`` rows outside ``dev`` (W5.26).

    ``tests/m2/conftest.py``'s hook, restated for this directory:
    ``tests/results/test_population.py``'s 200-object archive (fitting and
    timing 200 tiny ``emcee`` runs) is minutes of wall clock a per-PR gate
    should not pay for by default, so it is opt-in. A local hook rather than
    an ``addopts`` deselection in ``pyproject.toml``, so it says why it
    skipped in the report and leaves every other suite's invocation alone.
    """
    selected = config.getoption("-m", default="") or ""
    if "population_full" not in selected:
        skip = pytest.mark.skip(
            reason="the population study's full 200-object archive: run with "
            "`pytest -m population_full`"
        )
        for item in items:
            if "population_full" in item.keywords:
                item.add_marker(skip)
    _skip_study_rows_outside_dev(items)
