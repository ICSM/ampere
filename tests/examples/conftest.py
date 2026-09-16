"""Shared fixtures for ``tests/examples``.

``examples/wstat_comparison.py`` is a real, standalone user example
(``python examples/wstat_comparison.py`` is how it is meant to be run), so
it registers its ``ProfiledCashWithBackground`` family with
``@register_family`` at *module* scope, exactly the way a genuine user
module would. Loading it more than once by class identity would then trip
``register_family``'s own duplicate-name guard (a freshly executed module
gets a fresh class object, and the guard rejects a second, distinct class
under the same name) — so it is loaded exactly once per test module, here,
rather than once per test.

That registration is still the same module-global ``_FAMILIES`` mutation
``tests/core/conftest.py`` documents at length: left alone, it would stay
visible to every test that runs after it in the same pytest process,
including ``tests/conformance/test_likelihoods.py``'s check that the
marginalisation table covers exactly the registered family set. This
directory is not part of any gate that shares a process with
``tests/conformance`` today (``test-all`` names ``tests/core tests/results
tests/conformance tests/backends tests/inference`` explicitly), but the
fixture below snapshots and restores ``_FAMILIES`` around the one import
regardless, so that remains true if that ever changes.
"""

from __future__ import annotations

import importlib.util
import pathlib
from collections.abc import Iterator
from types import ModuleType

import pytest

_EXAMPLE_PATH = pathlib.Path(__file__).resolve().parents[2] / "examples" / "wstat_comparison.py"


@pytest.fixture(scope="module")
def wstat_example() -> Iterator[ModuleType]:
    """The loaded ``examples/wstat_comparison.py`` module, family registry restored after."""
    from ampere.core import likelihood as _likelihood_module

    snapshot = dict(_likelihood_module._FAMILIES)
    spec = importlib.util.spec_from_file_location("wstat_comparison_example", _EXAMPLE_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    try:
        yield module
    finally:
        _likelihood_module._FAMILIES.clear()
        _likelihood_module._FAMILIES.update(snapshot)


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``image_full`` rows unless ``-m image_full`` asked for them (W5.5).

    ``tests/interferometry/conftest.py``'s hook, under this directory's roof
    because that is where ``tests/examples/test_image_study.py`` lives: a local
    hook rather than an ``addopts`` deselection, so it says why it skipped and
    does not reach into any other suite's invocation.
    """
    selected = config.getoption("-m", default="") or ""
    if "image_full" in selected:
        return
    skip = pytest.mark.skip(reason="the image study's full budget: run with `pytest -m image_full`")
    for item in items:
        if "image_full" in item.keywords:
            item.add_marker(skip)
