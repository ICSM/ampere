"""Repository-wide pytest setup: the one place the root goes on ``sys.path``.

**W5.28(k), deciding and recording the question the wave-3/4 reviews asked
of the sys.path block this file replaces.** ``tests/`` is not, and stays
not, an importable package: no ``__init__.py`` anywhere under it, and none
added here. Making it one (so that pytest's own ``rootdir``-walking import
machinery would plant the repository root on ``sys.path`` for free) was
considered and rejected -- not because it would not work today, but because
it would be re-deriving, implicitly and fragile-ly, exactly the one fact
this file now states explicitly: several suites (``tests/m2``,
``tests/interferometry``, ``tests/benchmarks``, the many-lines rows of
``tests/inference/test_nuts.py``) import :mod:`examples` — a directory the
built wheel deliberately does not ship (``pyproject.toml``'s ``include =
["ampere*"]``) — so the repository root must reach ``sys.path`` by a route
that does not depend on how ``ampere`` itself was installed. Two of those
suites' own docstrings had already reasoned this through and landed on the
same answer *without* being able to point at each other's file, which is
what "the conftest sys.path pattern is used four times" was: not a defect
to design away, but one decision, taken independently, four times over,
because there was nowhere to write it down once. This file is that once.

**Why explicit, rather than pytest's own rootdir insertion.** Two concrete
ways the implicit route was seen to fail, both already lived through:
``tests/m2/conftest.py``'s original wording -- "an editable install's path
hook reaching this far would break the first time somebody ran it against a
wheel" (pytest's own insertion is keyed to *where a test file lives*, which
answers a different question from "does :mod:`examples` exist on this
machine at all") -- and ``tests/benchmarks/conftest.py``'s -- "``pixi run
bench`` invokes the bare ``pytest`` entry point and that inserts the test
file's directory and not the project's" (a package-chain walk changes
*which* directory pytest inserts, but the entry point still decides whether
it inserts one for :mod:`examples` at all). A conftest as a plain file
sidesteps both: pytest finds and runs every ``conftest.py`` from
``rootdir`` down to a collected test's own directory purely by walking the
filesystem, never by asking Python's import system where a package's root
is -- which is a guarantee ``__init__.py`` files would not add anything to,
and unconditional collection of *this* file already gives for free.

**What this replaces.** The identical block --
``_ROOT = pathlib.Path(__file__).resolve().parents[N]; if str(_ROOT) not in
sys.path: sys.path.insert(0, str(_ROOT))`` -- used to be repeated in
``tests/m2/conftest.py``, ``tests/interferometry/conftest.py``,
``tests/benchmarks/conftest.py``, ``tests/benchmarks/test_m2_sampling.py``,
``tests/benchmarks/test_m2_likelihood.py`` and
``tests/inference/test_nuts.py``, each paying to rediscover the same
reasoning. pytest imports this file before any test or conftest below it in
the tree, so one insertion here reaches every one of them; each has been
trimmed to rely on it rather than repeat it.

**What this does not replace.** ``tests/inference/test_interferometry.py``
puts a *different* directory (``tests/backends``, for a sibling fixture
module) on ``sys.path``, for a different reason -- a same-tree neighbour,
not the repository root -- and keeps doing so unchanged; the ``examples/``
single-script loaders (``tests/examples/test_cached_fit.py`` and its
siblings) insert and then remove a path around one dynamic import by name,
which is a different problem (loading one standalone script as a module)
solved a different way (temporarily, per import, because two such scripts
could collide by basename) and has nothing to gain from a permanent,
process-wide insertion.

**The ``study`` marker (W5.32 (i)).** ``_skip_study_rows_outside_dev`` used
to be defined twice -- identically, down to the docstring's reasoning --
in ``tests/m2/conftest.py`` and ``tests/results/conftest.py``, which meant
the marker only did anything for items collected under those two
directories; ``tests/astrometry/test_recovery.py`` carried its own local
``dev_only`` skipif for exactly the same purpose because there was no
project-wide hook it could lean on. ``pytest_collection_modifyitems`` is
not a "first result wins" hook -- pytest calls every conftest that defines
it, each with the *whole* session's item list, not just the items under
its own directory -- so one copy here, rather than one per suite, is
enough to cover every ``study``-marked row anywhere in the tree. Each
directory that still has its own opt-in ``*_full`` marker to skip
(``tests/m2/conftest.py``'s ``m2_full``, ``tests/results/conftest.py``'s
``population_full``) keeps that local hook for its own marker and no
longer needs to call this one itself.
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys

import pytest

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def _skip_study_rows_outside_dev(items: list[pytest.Item]) -> None:
    """Skip ``study``-marked rows wherever a modern backend is installed (W5.26).

    ``study`` marks a row that samples on the reference backend only and
    asserts a claim about a *likelihood*, about ``ampere.results`` machinery,
    or about a study's own recovery -- never about the array library. Once a
    per-backend agreement test (or the modern-backend NUTS/nested-sampling
    rows a suite already runs) has proven the backends agree with the
    reference implementation, re-running a ``study`` row's own numpy-path
    sampling in ``torch``, ``jax`` or ``sbi`` buys no additional evidence,
    only their wall clock -- so it runs in ``dev`` only.

    Detected by whether ``torch``/``jax`` import here, rather than by
    environment name -- a ``sbi``-environment run has torch installed and
    should skip these rows for the same reason a ``torch`` run does.
    """
    reasons = [name for name in ("torch", "jax") if importlib.util.find_spec(name) is not None]
    if not reasons:
        return
    skip = pytest.mark.skip(
        reason=(
            "study row: numpy-path-only, run on `dev` only (see pyproject.toml's `study` "
            f"marker); this environment has {' and '.join(reasons)} installed"
        )
    )
    for item in items:
        if "study" in item.keywords:
            item.add_marker(skip)


# W6.10 (the xdist audit, docs/development.md "The merged gate"): suites that
# are correct under ``-n`` but slower for it, kept to ONE worker each by
# ``--dist loadgroup`` (the task lines pass it; without it, or without
# pytest-xdist, the marker is inert and nothing here applies).
# ``interferometry``: three class-scoped calibration fits (~100 s of setup
# each) are rebuilt in every worker that draws one of their rows; ``astrometry``:
# two nested-sampling rows (~270 s and ~160 s) bound its wall clock whatever
# the worker count.
_SERIAL_UNDER_XDIST = ("interferometry", "astrometry")


def _group_serial_suites_for_xdist(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Pin each of :data:`_SERIAL_UNDER_XDIST` to one xdist worker (W6.10)."""
    if not config.pluginmanager.hasplugin("xdist"):
        return
    here = pathlib.Path(__file__).resolve().parent
    for item in items:
        try:
            suite = item.path.resolve().relative_to(here).parts[0]
        except (ValueError, IndexError):
            continue
        if suite in _SERIAL_UNDER_XDIST:
            item.add_marker(pytest.mark.xdist_group(f"serial-{suite}"))


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``study`` rows outside ``dev`` (W5.32 (i)); group the serial suites (W6.10)."""
    _skip_study_rows_outside_dev(items)
    _group_serial_suites_for_xdist(config, items)
