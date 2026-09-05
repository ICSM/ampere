"""Concrete implementations of the ``ampere.core`` contracts, one per array library.

``architecture.md`` §1 makes this an escalating capability ladder rather than a
set of peer backends: ``reference`` is pure numpy/scipy (rung 1, correctness as
its only goal), and Phase 2's ``torch`` and ``jax`` subpackages add
differentiability (rung 2). A backend is precisely "an array library for
writing models and transformations, a set of GP solver implementations, and the
extra inference engines its differentiability unlocks" — everything else is
shared, and lives in ``ampere.core``.

**This module imports nothing.** ``architecture.md`` §4 rule 2 names
``ampere/backends/__init__.py`` specifically: it is reachable without opting
into any particular backend, so it must not pull in a subpackage's
dependencies. Import the backend you want by name::

    from ampere.backends import reference

``reference`` needs only the base install; ``torch`` and ``jax`` will need
their extras.
"""

from __future__ import annotations

__all__: list[str] = []
