"""Make ampere's frozen mappings picklable — the one thing that stood in the way.

Every container, parameter record, result and simulation in ``ampere.core``
freezes its mappings with :class:`types.MappingProxyType`, which is how a
``ModelResult``'s channels or a ``Failure``'s recorded values are handed out
without handing out something a caller can edit. That decision is right and
predates this module. Its one consequence is that **nothing holding such a
mapping could be pickled**: ``mappingproxy`` has no pickle reduction of its
own, so ``pickle.dumps(problem)`` failed with ``cannot pickle 'mappingproxy'
object`` — as did a :class:`~ampere.core.dataset.Simulation`, a
:class:`~ampere.core.results_schema.ModelResult` and every ``FunctionSamples``
subclass.

W3.1 slice 1 needs exactly that: ``simulate_many``'s process-pool executor
sends a :class:`~ampere.core.dataset.FittingProblem` to each worker once and
receives ``Simulation`` records back, which is how a slow external simulator
(``DEVELOPMENT_PLAN.md`` §2's "Batched simulation" row, note (2)) is run in
parallel at all.

:mod:`copyreg` is the mechanism the standard library provides for exactly this
situation — a reduction for a type you do not own — so the fix is one
registration rather than a ``__getstate__``/``__setstate__`` pair on each of
the dozen classes that freeze a mapping, which would have been a dozen chances
to miss one. The registration is *additive*: ``mappingproxy`` could not be
pickled before, so nothing that used to work changes behaviour, and a proxy
round-trips as a proxy rather than degrading to a plain ``dict``.

Examples
--------
>>> import pickle, types
>>> proxy = types.MappingProxyType({"a": 1})
>>> restored = pickle.loads(pickle.dumps(proxy))
>>> restored == proxy, type(restored) is type(proxy)
(True, True)
"""

from __future__ import annotations

import copyreg
import types
from collections.abc import Mapping
from typing import Any

__all__ = ["rebuild_mappingproxy", "register"]

#: The ``mappingproxy`` type. It is not exposed under a public name anywhere in
#: the standard library, so it is obtained the way everyone obtains it.
_MAPPING_PROXY = type(types.MappingProxyType({}))


def rebuild_mappingproxy(mapping: Mapping[Any, Any]) -> types.MappingProxyType[Any, Any]:
    """Reconstruct a frozen mapping from the plain ``dict`` that was pickled.

    A module-level function rather than :class:`types.MappingProxyType` itself,
    because pickle stores the *callable* by qualified name and ``mappingproxy``
    is not reachable as an attribute of any module.
    """
    return types.MappingProxyType(dict(mapping))


def _reduce_mappingproxy(
    proxy: types.MappingProxyType[Any, Any],
) -> tuple[Any, tuple[dict[Any, Any]]]:
    return rebuild_mappingproxy, (dict(proxy),)


def register() -> None:
    """Install the reduction. Idempotent; called once when ``ampere.core`` loads."""
    copyreg.pickle(_MAPPING_PROXY, _reduce_mappingproxy)
