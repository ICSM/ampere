"""Realisation: a composed problem as one backend-native, differentiable density.

**The contract this module implements is ``inference.md`` §10a** (added by
W2.13, ruled 2026-09-07; the decision-log row "Realisation surface (W2.13)" in
``DEVELOPMENT_PLAN.md`` §2 is the record). It closes the gap both Phase 2
backend tracks found independently: ``DEVELOPMENT_PLAN.md`` §4.5's
engine-facing surface —
``FittingProblem.log_prob_unconstrained`` and friends — is **not traceable**
on any backend, because every ``ampere.core`` container coerces its values
with ``numpy.asarray`` (``results_schema.py``) and
``ParameterSet.lnprior_unconstrained`` short-circuits on ``math.isfinite``. A
torch tensor that requires grad cannot be put into a ``Spectrum`` at all; a jax
tracer raises ``TracerArrayConversionError`` at the same line. So a
gradient-based engine cannot be written against §4.5 alone, and
``ampere.inference`` — which may import nothing but ``ampere.core`` and
``ampere.results`` — cannot reach a backend's native path by import either.

``lowering.md`` §0 already says lowering "happens once, when a
``FittingProblem`` is *realised* on a backend — never per evaluation", and
names nothing further. This module names it. A **realisation** is the
backend's one-way translation of a composed problem into a native object
exposing the log density as a pure function of the unconstrained free vector,
built from the backend's own model, step, noise and solver evaluations — never
from core containers in the hot loop. A backend registers a factory here, at
import, keyed on its one name (the ``BACKEND`` capability flag, W2.12); an
engine calls :func:`realise` through ``ampere.core`` and dispatches on
``problem.backend``, so the import-graph rule stays literally true.

The registry has the shape of :mod:`ampere.core.lowering`'s, for the same
reasons: module-global, keyed on the backend string, no silent overwrites,
built-in rows distinguished from user rows. It is deliberately smaller —
one slot, one key — because a backend realises a *problem* in one way.

**The numpy contract path remains the oracle.** :func:`realise` checks the
cheapest necessary condition at construction: the realised density agrees with
``problem.log_prob_unconstrained`` at the problem's reference values. That is
not a proof — the conformance suite makes the full comparison — but a
realisation that had drifted from its contract path would otherwise sample a
healthy-looking posterior of a different problem.

What a realisation is **not**: it is not a replacement for §4.5. The
gradient-free engines keep consuming the contract path unchanged; a backend
that registers nothing (the reference backend, which has no differentiable
path) simply cannot be realised, and :func:`realise` says so by name with the
remedy.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

import numpy as np

from .exceptions import LoweringError, ParameterError

if TYPE_CHECKING:
    from .dataset import FittingProblem

__all__ = [
    "Realisation",
    "RealisationFactory",
    "log_likelihood_terms_of",
    "realise",
    "register_realisation",
    "registered_realisations",
]


@runtime_checkable
class Realisation(Protocol):
    """What a backend hands back from its realisation factory.

    ``inference.md`` §10a, "The surface". The three members below are
    **mandatory** and are the whole of what a gradient-based engine needs;
    ``theta`` is the unconstrained free vector in the backend's own array type
    (a numpy array is always acceptable input — the backend converts), and the
    return is a native scalar the backend's autodiff can differentiate.

    Deliberately minimal, because every member is a member three backends must
    implement in lockstep: this Protocol is the whole of the obligation a new
    backend takes on to unlock NUTS, and §10a fixes it at the point where
    ``ampere.inference`` stops needing anything more.

    One member is **optional** and is looked for with ``getattr`` rather than
    declared here, so that a realisation supplying only the mandatory three
    still satisfies :func:`isinstance` against this runtime-checkable
    Protocol: see :func:`log_likelihood_terms_of`.

    Failure signalling narrows ``inference.md`` §11 on this path, and §10a
    says so explicitly: inside a trace nothing can raise and no ``Failure``
    record can be built, so :meth:`log_prob_unconstrained` returns a **bare**
    ``-inf`` (computed with the backend's ``where``) and the reasons §11
    promises are recovered afterwards, by re-evaluating stored draws on the
    numpy path. ``strict=True`` on a realised problem means the *factory*
    refuses at construction, by name; a runtime evaluation failure is always
    ``-inf``.
    """

    @property
    def backend(self) -> str:
        """The backend's one name; must equal the realised problem's ``backend``."""

    @property
    def free_size(self) -> int:
        """Length of the unconstrained free vector; must equal the problem's."""

    def log_prob_unconstrained(self, theta: Any) -> Any:
        """``log p(θ) + log p(data | θ)`` plus the change of variables, natively."""


def log_likelihood_terms_of(realised: Realisation) -> Callable[[Any], Mapping[str, Any]] | None:
    """*realised*'s optional per-dataset decomposition, or ``None``.

    ``inference.md`` §10a's optional member::

        log_likelihood_terms(theta) -> mapping of dataset label to native scalar

    A driver **uses it when it is there**, to emit
    ``results.md``'s per-dataset decomposition natively and for every draw the
    sampler took; **absent, the driver recomputes the split on the numpy path
    for stored draws only** and records that it did, in
    ``engine_draws_recomputed``. Both halves of that sentence are the ruling
    (sub-decision 1), and this function is the one place the "is it there?"
    question is asked, so a driver cannot spell the member's name differently
    from the backend that supplies it.

    It is fetched rather than declared on :class:`Realisation` because
    ``runtime_checkable`` Protocols check *presence*: declaring it would make
    every minimal realisation fail ``isinstance``, which is exactly the check
    :func:`realise` uses to give a backend a legible error.
    """
    terms = getattr(realised, "log_likelihood_terms", None)
    return terms if callable(terms) else None


#: ``FittingProblem -> Realisation``. Must refuse, by name, at construction —
#: never inside the traced density — anything it cannot lower
#: (``likelihoods.md`` §17 Q1's trace-purity ruling).
RealisationFactory = Callable[["FittingProblem"], Realisation]


class _Row:
    __slots__ = ("backend", "builtin", "factory")

    def __init__(self, backend: str, factory: RealisationFactory, builtin: bool) -> None:
        self.backend = backend
        self.factory = factory
        self.builtin = builtin


#: Module-global, like ``ampere.core.lowering._REGISTRY``; ``tests/core/conftest.py``
#: snapshots and restores it per test for the same reason.
_REALISATIONS: dict[str, _Row] = {}


def register_realisation(
    backend: str,
    factory: RealisationFactory,
    *,
    override: bool = False,
    builtin: bool = False,
) -> None:
    """Register *factory* as the way problems on *backend* are realised.

    ``inference.md`` §10a, "The registry": one slot per backend, keyed on the
    backend's one name, no silent overwrites, built-in rows distinguished. A
    backend registers **at import**, so that importing
    ``ampere.backends.jax`` is both the user's opt-in to jax and the moment
    the jax realisation becomes reachable.

    Parameters
    ----------
    backend
        The backend's one name (its ``BACKEND`` capability flag, the lowering
        registry's key, its conformance fixture's ``name``).
    factory
        ``FittingProblem -> Realisation``. Called once per :func:`realise`;
        must refuse at construction anything it cannot lower.
    override
        Required to replace an existing row; otherwise this raises naming it.
    builtin
        Reserved for the backends ampere itself ships; a built-in row still
        needs ``override=True`` to replace.

    Raises
    ------
    ParameterError
        On a malformed key or a non-callable factory.
    LoweringError
        If a row already exists for *backend* and ``override`` is not ``True``.
    """
    if not isinstance(backend, str) or not backend:
        raise ParameterError(f"a backend name must be a non-empty string, got {backend!r}")
    if not callable(factory):
        raise ParameterError(f"a realisation factory must be callable, got {factory!r}")
    existing = _REALISATIONS.get(backend)
    if existing is not None and not override:
        kind = "ampere's own" if existing.builtin else "a user-registered"
        raise LoweringError(
            "realisation",
            backend=backend,
            detail=(
                f"a realisation is already registered for backend {backend!r} ({kind} row). "
                f"Pass override=True to replace it deliberately; silent replacement is refused "
                f"because whichever factory won would decide what every gradient-based engine "
                f"samples."
            ),
        )
    _REALISATIONS[backend] = _Row(backend, factory, builtin)


def registered_realisations() -> Mapping[str, bool]:
    """``{backend: builtin}`` for every registered row.

    ``inference.md`` §10a lists this as the registry's third member, "for
    provenance and tests": a run's attrs record whether it sampled through a
    realisation, and a conformance row is generated per registered backend.
    """
    return {name: row.builtin for name, row in _REALISATIONS.items()}


def realise(problem: FittingProblem) -> Realisation:
    """The backend-native, differentiable form of *problem*.

    ``inference.md`` §10a, "The registry": dispatches on ``problem.backend``
    (W2.12's derived flag) to the registered factory, then **checks the result
    on use**, the way an engine would otherwise have to — it names the same
    backend, has the problem's ``free_size``, and agrees with the contract path
    at the problem's reference values.

    The one-point agreement check is a **guard, not the proof** (sub-decision
    4): the numpy path stays the oracle, and ``tests/conformance`` compares
    the two at many points, including near a support boundary, once per
    registered realisation. What one point buys is that a driver handed a
    realisation which had drifted from its own contract path fails here rather
    than sampling a healthy-looking posterior of a different problem.

    Raises
    ------
    LoweringError
        If no realisation is registered for the problem's backend, naming the
        backend and the remedy; if the factory's result does not satisfy
        :class:`Realisation`; or if it disagrees with
        ``problem.log_prob_unconstrained`` at the reference point.
    """
    backend = problem.backend
    row = _REALISATIONS.get(backend)
    if row is None:
        known = ", ".join(sorted(_REALISATIONS)) or "(none in this environment)"
        raise LoweringError(
            "realisation",
            backend=backend,
            detail=(
                f"no realisation is registered for backend {backend!r}: this problem has no "
                f"differentiable native form. Backends with one registered here: {known}. "
                f"Import the backend package that supplies it (importing ampere.backends.jax "
                f"registers the jax realisation), build the problem from that backend's "
                f"pieces, or use a gradient-free engine (emcee, dynesty, zeus), which run the "
                f"contract path on every backend."
            ),
        )
    realised = row.factory(problem)
    if not isinstance(realised, Realisation):
        raise LoweringError(
            "realisation",
            backend=backend,
            detail=(
                f"the realisation factory for backend {backend!r} returned {realised!r}, which "
                f"does not provide backend, free_size and log_prob_unconstrained."
            ),
        )
    if realised.backend != backend:
        raise LoweringError(
            "realisation",
            backend=backend,
            detail=(
                f"the realisation factory for backend {backend!r} returned an object naming "
                f"backend {realised.backend!r}."
            ),
        )
    if int(realised.free_size) != problem.free_size:
        raise LoweringError(
            "realisation",
            backend=backend,
            detail=(
                f"the realisation has free_size {int(realised.free_size)} but the problem has "
                f"{problem.free_size}; the two must lower the same declaration."
            ),
        )
    _check_agrees_at_reference(problem, realised)
    return realised


def _check_agrees_at_reference(problem: FittingProblem, realised: Realisation) -> None:
    """The cheapest necessary condition that the realisation is *this* problem's."""
    reference = problem.unconstrain(problem.reference_values)
    expected = float(problem.log_prob_unconstrained(reference))
    try:
        got = float(np.asarray(realised.log_prob_unconstrained(reference)))
    except Exception as error:
        raise LoweringError(
            "realisation",
            backend=problem.backend,
            detail=(
                f"the realised density could not be evaluated at the problem's reference values "
                f"(a vector of length {problem.free_size}): {error}"
            ),
        ) from error
    if not math.isclose(got, expected, rel_tol=1e-6, abs_tol=1e-6):
        raise LoweringError(
            "realisation",
            backend=problem.backend,
            detail=(
                f"the realised density disagrees with the contract path at the problem's "
                f"reference values: realised {got!r}, FittingProblem.log_prob_unconstrained "
                f"{expected!r}. They are one quantity computed twice; sampling the first while "
                f"recording the second would describe a different posterior."
            ),
        )
