"""Shared fixtures for ``tests/core``.

``ampere.core.likelihood`` keeps its family registry (``_FAMILIES``) as
module-global state, and ``@register_family`` has no public deregistration
route (W1.6: registering your own family is meant to be a one-line, no-import
affair, not one that also asks you to unregister it). Several tests in this
directory register throwaway families to exercise that "a user can add their
own without touching ampere" contract -- as does the worked example in
``docs/design/contracts/likelihoods.md`` (a "laplace" family), run via
``test_spec_doctests.py``. Left alone, a family registered by one test stays
visible to every test that runs after it in the same process, including
``tests/conformance/test_likelihoods.py``'s check that the marginalisation
table covers exactly the registered set -- W0.10 finding (c): running
``tests/core`` and ``tests/conformance`` together in one pytest invocation
failed that check even though each suite passed alone.

``ampere.core.lowering`` (W2.6) has the identical shape: ``_REGISTRY`` is
module-global mutable state, ``register_lowering``/``register_bijection_lowering``
have no deregistration route (the same one-line, no-import contract), and
both ``tests/core/test_lowering.py`` and ``ampere.core.lowering``'s own
module docstring (run as a doctest by ``test_spec_doctests.py``) register
throwaway rows on invented backend names. Left alone, those rows would leak
across tests -- and, per Finding 1 of the W2.6 review, across the whole
``test-all`` process -- the same way an unrestored ``_FAMILIES`` did.

Snapshotting and restoring each registry's underlying mapping from the test
side, autouse around every test in this directory, closes that gap without
touching the (frozen, merged) Phase-1 contract module or the W2.6 registry
module itself.
"""

from __future__ import annotations

from collections.abc import Iterator

import pytest


@pytest.fixture(autouse=True)
def _restore_likelihood_family_registry() -> Iterator[None]:
    """Snapshot ``ampere.core.likelihood._FAMILIES`` and restore it after each test."""
    from ampere.core import likelihood as _likelihood_module

    snapshot = dict(_likelihood_module._FAMILIES)
    try:
        yield
    finally:
        _likelihood_module._FAMILIES.clear()
        _likelihood_module._FAMILIES.update(snapshot)


@pytest.fixture(autouse=True)
def _restore_lowering_registry() -> Iterator[None]:
    """Snapshot ``ampere.core.lowering._REGISTRY`` and restore it after each test."""
    from ampere.core import lowering as _lowering_module

    snapshot = dict(_lowering_module._REGISTRY)
    try:
        yield
    finally:
        _lowering_module._REGISTRY.clear()
        _lowering_module._REGISTRY.update(snapshot)


@pytest.fixture(autouse=True)
def _restore_realisation_registry() -> Iterator[None]:
    """Snapshot ``ampere.core.realisation._REALISATIONS`` and restore it after each test.

    The same leak class as the lowering registry above (W0.10 finding (c)):
    a test that registers a realisation for a fake backend must not leave it
    behind for the single-process ``test-all`` run.
    """
    from ampere.core import realisation as _realisation_module

    snapshot = dict(_realisation_module._REALISATIONS)
    try:
        yield
    finally:
        _realisation_module._REALISATIONS.clear()
        _realisation_module._REALISATIONS.update(snapshot)
