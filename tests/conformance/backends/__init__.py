"""The backend registry: the one place a new backend is named.

W1.10's acceptance criterion is that adding a backend to the conformance
battery "requires only a fixture". This module *is* that fixture. To add one:

1. write a class satisfying :class:`tests.conformance.protocol.ConformanceBackend`
   (``tests/conformance/README.md`` is the specification);
2. add one entry to :data:`_REGISTRY` below, guarded by
   :func:`_optional` if it needs a dependency the base install does not have.

Nothing else changes. No test body names a backend, and no test body may.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence

from ..protocol import ConformanceBackend
from .mirror import MirrorBackend
from .reference import ReferenceBackend

__all__ = ["available_backends", "backend_ids"]


def _optional(factory: Callable[[], ConformanceBackend]) -> ConformanceBackend | None:
    """Build a backend, or return ``None`` if its dependency is absent.

    Phase 2's torch and jax fixtures go through here: ``architecture.md`` §4
    rule 2 forbids making ``import ampere`` require either, so the battery has
    to degrade to the backends actually installed rather than erroring at
    collection. A missing backend is reported by its absence from the fixture
    ids, not by a skipped row — the rows themselves are backend-agnostic.
    """
    try:
        return factory()
    except ImportError:  # pragma: no cover - no optional backend exists yet
        return None


_REGISTRY: tuple[ConformanceBackend | None, ...] = (
    ReferenceBackend(),
    MirrorBackend(),
    # Phase 2:
    #   _optional(lambda: TorchBackend()),
    #   _optional(lambda: JaxBackend()),
)


def available_backends() -> tuple[ConformanceBackend, ...]:
    """Every backend installed in this environment, in registration order."""
    return tuple(backend for backend in _REGISTRY if backend is not None)


def backend_ids(backends: Sequence[ConformanceBackend]) -> list[str]:
    """Fixture ids for *backends* — used for both parametrisation and reports."""
    return [backend.name for backend in backends]
