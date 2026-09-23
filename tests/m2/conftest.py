"""Shared fixtures for ``tests/m2`` — milestone M2's CI-runnable assertions.

The study lives in ``examples/m2_misspecification`` because it is an example: a
reader is meant to run it, change a scenario and see what happens. That leaves
this suite two jobs, and both are done here rather than repeated per module.

**Importing it.** ``examples/`` is not an installed package — the wheel
declares ``include = ["ampere*"]`` — so the repository root must reach
``sys.path`` some way that does not depend on how ``ampere`` itself was
installed. ``tests/conftest.py`` (**W5.28(k)**) now does this once for the
whole suite, rather than a copy of the insertion living here as well.

**Paying for the chains once.** The eight runs the science assertions and the
figures both read (four scenarios x two likelihoods, at
``study.TEST_EMCEE``) cost about eighty seconds on the reference backend.
That is the whole budget for this suite, so the fixture is **session**-scoped
and every module reads the same eight runs. The consequence to be aware of is
that :func:`~examples.m2_misspecification.study.prepare` mutates the entries
in place — attaching the derived ``residuals`` / ``gp_localisation`` groups —
so a test must not assume it received a freshly emitted run.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Any

import pytest

from examples.m2_misspecification import study

#: The size every CI assertion is made at: the paper study's own, and the one
#: rung of the ladder whose chains fit in a per-PR budget. The other two rungs
#: are ``pytest -m m2_full`` (``tests/m2/test_ladder.py``).
CI_SIZE = 200


@pytest.fixture(scope="session")
def _study_session() -> tuple[
    dict[tuple[str, str], dict[str, Any]], dict[tuple[str, str], study.Diagnosis]
]:
    """The eight reference-backend runs and their diagnostics, computed once.

    Reference-backend only, deliberately, and in every environment: the science
    claim is a claim about the *likelihood*, not about the array library, and
    making it once — on the backend that is always present — keeps the torch and
    jax jobs spending their time on what only they can test, which is
    ``tests/m2/test_backend_agreement.py``.

    The runs and the diagnostics come back together rather than from two
    fixtures, because :func:`~examples.m2_misspecification.study.prepare` does
    both jobs in one pass: it derives each run's ``residuals`` or
    ``gp_localisation`` group *and* returns the statistic computed from it.
    Splitting them would either derive the groups twice or leave the second
    fixture reading numbers the first had already thrown away.
    """
    results = study.run_study(size=CI_SIZE, budget=study.TEST_EMCEE)
    return results, study.prepare(results, thin=20)


@pytest.fixture(scope="session")
def study_results(
    _study_session: tuple[dict[tuple[str, str], dict[str, Any]], Any],
) -> dict[tuple[str, str], dict[str, Any]]:
    """The eight runs, each carrying the derived group its diagnostic needed."""
    return _study_session[0]


@pytest.fixture(scope="session")
def diagnoses(
    _study_session: tuple[Any, dict[tuple[str, str], study.Diagnosis]],
) -> dict[tuple[str, str], study.Diagnosis]:
    """The whiteness / localisation answers for those same eight runs."""
    return _study_session[1]


@pytest.fixture(scope="session")
def agg_backend() -> Iterator[None]:
    """Select matplotlib's non-interactive backend before anything imports pyplot."""
    import matplotlib

    previous = matplotlib.get_backend()
    matplotlib.use("Agg", force=True)
    try:
        yield
    finally:
        matplotlib.use(previous, force=True)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``m2_full`` rows unless ``-m m2_full`` asked for them.

    The full ladder samples 2 000- and 20 000-point spectra at the milestone
    budget: tens of minutes, on three backends. That is evidence worth having
    and recording, and it is not something a per-PR gate can pay for, so it is
    opt-in.

    Skipping here rather than with a global ``addopts = "-m 'not m2_full'"`` in
    ``pyproject.toml`` is deliberate: an ``addopts`` deselection is invisible
    from the command line, applies to every suite in the repository, and
    silently conflicts with any other ``-m`` a caller passes. This hook is
    local to ``tests/m2``, says why it skipped in the report, and leaves every
    other suite's invocation exactly as it was.
    """
    selected = config.getoption("-m", default="") or ""
    if "m2_full" in selected:
        return
    skip = pytest.mark.skip(reason="milestone M2's full ladder: run with `pytest -m m2_full`")
    for item in items:
        if "m2_full" in item.keywords:
            item.add_marker(skip)
