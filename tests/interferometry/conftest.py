"""Shared fixtures for ``tests/interferometry`` — W4.4's CI-runnable assertions.

The study lives in ``examples/interferometry`` because it is an example: a
reader is meant to run it, change an arm and see what happens.
``tests/m2/conftest.py``'s own two jobs, done here in its own shape.

**Importing it.** ``examples/`` is not an installed package, so the
repository root must reach ``sys.path`` some way that does not depend on
how ``ampere`` itself was installed. ``tests/conftest.py`` (**W5.28(k)**)
does this once for the whole suite now, rather than a copy of the
insertion living here as well.

**Paying for the calibration study once.** The three SBC runs (one per arm,
:data:`~examples.interferometry.study.SIMULATIONS` simulations each) that the
science assertions and the calibration figures both read cost under three
minutes on the reference backend combined — see the module's own timing note
in the branch report. That is the whole budget for this suite's pinned
assertions, so the fixture is **session**-scoped.
"""

from __future__ import annotations

from typing import Any

import pytest

from examples.interferometry import study


@pytest.fixture(scope="session")
def calibrations() -> dict[str, Any]:
    """The three arms' SBC calibration datasets, computed once per session."""
    import warnings

    with warnings.catch_warnings():
        # SIMULATIONS is a per-PR budget, below sbc()'s own goodness-of-fit
        # power floor — it warns about exactly that, and the warning is
        # allowed through, as tests/m2/test_visibility_calibration.py's own
        # row does.
        warnings.filterwarnings("ignore", message=r".*goodness-of-fit test.*")
        return {arm: study.run_calibration(arm) for arm in study.ARMS}


@pytest.fixture(scope="session")
def agg_backend():
    """Select matplotlib's non-interactive backend before anything imports pyplot."""
    import matplotlib

    previous = matplotlib.get_backend()
    matplotlib.use("Agg", force=True)
    try:
        yield
    finally:
        matplotlib.use(previous, force=True)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``interferometry_full`` rows unless ``-m interferometry_full`` asked for them.

    ``tests/m2/conftest.py``'s own reasoning, unchanged: a local hook rather
    than an ``addopts`` deselection, so it says why it skipped and does not
    reach into any other suite's invocation.
    """
    selected = config.getoption("-m", default="") or ""
    if "interferometry_full" in selected:
        return
    skip = pytest.mark.skip(
        reason="the interferometry study's full budget: run with `pytest -m interferometry_full`"
    )
    for item in items:
        if "interferometry_full" in item.keywords:
            item.add_marker(skip)
