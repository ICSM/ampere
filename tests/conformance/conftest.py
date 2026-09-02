"""Fixtures for the conformance battery.

Three fixtures and one helper, and that is the whole of the machinery:

``backend``
    One :class:`~tests.conformance.protocol.ConformanceBackend` per registered
    backend, so every test function in this directory runs once per backend
    with its name as the pytest id. Registration lives in
    ``tests/conformance/backends/__init__.py``.
``tolerances``
    That backend's :class:`~tests.conformance.protocol.Tolerances` table.
``backends``
    *All* registered backends at once, for the cross-backend rows — the ones
    that compare two implementations rather than holding one to an oracle.
    Those rows skip when only one backend is installed.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from .backends import available_backends, backend_ids
from .protocol import ConformanceBackend, Tolerances

BACKENDS = available_backends()


@pytest.fixture(params=BACKENDS, ids=backend_ids(BACKENDS))
def backend(request: pytest.FixtureRequest) -> ConformanceBackend:
    """The backend under test. Every row in this directory takes it."""
    return request.param


@pytest.fixture
def tolerances(backend: ConformanceBackend) -> Tolerances:
    """The per-comparison tolerance table *backend* declares."""
    return backend.capabilities.tolerances


@pytest.fixture(scope="session")
def backends() -> Iterator[tuple[ConformanceBackend, ...]]:
    """Every registered backend, for the rows that compare implementations."""
    yield BACKENDS


def cross_backend_pairs() -> list[tuple[ConformanceBackend, ConformanceBackend]]:
    """Every ordered pair of distinct backends, for parametrising agreement rows."""
    return [
        (first, second) for index, first in enumerate(BACKENDS) for second in BACKENDS[index + 1 :]
    ]


def pair_ids() -> list[str]:
    """Ids for :func:`cross_backend_pairs`."""
    return [f"{first.name}-vs-{second.name}" for first, second in cross_backend_pairs()]
