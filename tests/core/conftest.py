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

Snapshotting and restoring the registry's underlying mapping from the test
side, autouse around every test in this directory, closes that gap without
touching the (frozen, merged) Phase-1 contract module itself.
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
