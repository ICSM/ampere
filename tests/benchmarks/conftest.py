"""Shared setup for ``tests/benchmarks``.

``test_solver_bakeoff.py`` imports :mod:`examples.image`, which is a
repository-root package rather than an installed one, so the root goes on
``sys.path`` explicitly — ``tests/m2/conftest.py``'s pattern and its reasoning:
explicitly, rather than by relying on pytest's rootdir insertion, because
``pixi run bench`` invokes the bare ``pytest`` entry point and that inserts the
*test file's* directory and not the project's.

The ``image_full`` hook is ``tests/examples/conftest.py``'s, under this
directory's roof for the same reason it is under that one: ``pixi run bench``
is a per-PR job and the 128x128 and 256x256 cells of the bake-off are not
per-PR work.
"""

from __future__ import annotations

import pathlib
import sys

import pytest

_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    """Skip the ``image_full`` rows unless ``-m image_full`` asked for them."""
    selected = config.getoption("-m", default="") or ""
    if "image_full" in selected:
        return
    skip = pytest.mark.skip(reason="the image study's full budget: run with `pytest -m image_full`")
    for item in items:
        if "image_full" in item.keywords:
            item.add_marker(skip)
