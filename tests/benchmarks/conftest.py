"""Shared setup for ``tests/benchmarks``.

``test_solver_bakeoff.py`` imports :mod:`examples.image`, which is a
repository-root package rather than an installed one, so the root must reach
``sys.path`` some way that does not depend on how ``ampere`` itself was
installed -- ``tests/conftest.py`` (**W5.28(k)**) does this once for the
whole suite now, including for ``pixi run bench``'s bare ``pytest`` entry
point, which inserts only the *test file's* own directory and never the
project's on its own.

The ``image_full`` hook is ``tests/examples/conftest.py``'s, under this
directory's roof for the same reason it is under that one: ``pixi run bench``
is a per-PR job and the 128x128 and 256x256 cells of the bake-off are not
per-PR work.
"""

from __future__ import annotations

import pytest


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``image_full`` rows unless ``-m image_full`` asked for them."""
    selected = config.getoption("-m", default="") or ""
    if "image_full" in selected:
        return
    skip = pytest.mark.skip(reason="the image study's full budget: run with `pytest -m image_full`")
    for item in items:
        if "image_full" in item.keywords:
            item.add_marker(skip)
