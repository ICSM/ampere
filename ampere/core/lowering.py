"""The lowering registry: ``register_lowering`` and the bijection slot (W2.6).

``docs/design/lowering.md`` §12.8 (ruled 2026-09-03, ``DEVELOPMENT_PLAN.md``
§2's "User-registered lowerings" row) accepts, in principle and hardened, a
registration hook that lets a user supply a backend-native equivalent for a
prior family or a custom :class:`~ampere.core.parameter.Bijection` that
ampere's own built-in tables (``lowering.md`` §3.2, §4) do not cover. Without
it, a family torch or numpyro does not implement — or a user-defined
:class:`Bijection` — is confined to the reference path forever
(``lowering.md`` §3.4 rule 4, §4's "unsupported — raises at lowering" row);
with it, a registrant can run their own lowering natively, at the cost of a
self-certification step (§12.8's *opt-in* battery, :func:`run_registrant_battery`
below) that stands in for the shared conformance suite's guarantee.

This module is the plumbing, not a backend. It knows nothing about torch,
jax, numpyro or paramax — a *backend* is just a string key here — and it is
consulted by backends (W2.4/W2.5, ``tests/conformance``'s in-repo fixtures)
the same way a third party's own code would consult it. Two slots, one
mechanism underneath (see "One store, two slots" below):

* :func:`register_lowering` / :func:`lookup_lowering` — prior families, keyed
  on the neutral name :func:`~ampere.core.parameter.describe_prior` reads off
  a frozen ``scipy.stats`` distribution (``PriorSpec.family``). A constructor
  takes the :class:`~ampere.core.parameter.PriorSpec` and returns a
  backend-native distribution object.
* :func:`register_bijection_lowering` / :func:`lookup_bijection_lowering` —
  custom :class:`~ampere.core.parameter.Bijection` classes, keyed on the class
  (``LoweringError``'s docstring already calls this "a bijection class name").
  A constructor takes the :class:`Bijection` *instance* and returns a
  backend-native transform object.

Three hardenings, each load-bearing
------------------------------------
**No silent overwrite.** Registering over a built-in row or an existing
registration without ``override=True`` raises :class:`LoweringError`, naming
what it refused to overwrite (built-in or user, and by which constructor).

**Trace purity is a rule on constructors, not a runtime check.** The registry
is consulted once, at lowering time, before any tracing
(``lowering.md`` §0: "[lowering] happens once, when a ``FittingProblem`` is
realised on a backend — never per evaluation"). That is *why* the mechanism
itself cannot interact with ``jax.jit``/``vmap``/gradients: a dict lookup and
a function call, both outside any traced region, trace nothing. What can
still go wrong is a constructor that closes over mutable state, branches on
untraced Python values that vary per call, or otherwise returns an object
that is not a pure function of its argument — exactly the failure mode the
built-in table rows (``lowering.md`` §3.2, §4) already avoid by construction.
**A registered constructor must be a pure function of its ``PriorSpec`` or
``Bijection`` argument, with no side effects and no closed-over mutable
state**, so that the *object it returns* is as trace-safe as a built-in row's.

**Self-certification, not silent trust.** :func:`run_registrant_battery`
compares the registered constructor's native ``log_prob`` against the
scipy-backed reference at a handful of points spanning the prior's support —
``lowering.md`` §11 row 3's comparison, run by the registrant against their
own row rather than by the shared suite.

One store, two slots
---------------------
Both slots share one underlying mapping, keyed by ``(kind, name, backend)``
with ``kind`` in ``{"prior", "bijection"}`` — the overwrite refusal,
provenance-entry shape and lookup-failure message are therefore identical
prose for both, which is why they live in one module and one table rather
than two independent registries with two sets of rules to keep in sync. This
is an implementation choice the spec does not pin down (flagged in the W2.6
report); nothing about the *public* API depends on it, and splitting the
store later is a non-breaking change if review prefers it.

Provenance stamping
--------------------
The ruling asks that "user-registered rows [be] stamped in provenance."
:func:`provenance_entries` renders every non-built-in
:class:`LoweringResolution` a run consulted into a small, netCDF-safe
mapping. ``ampere.results.provenance.provenance_attrs`` already accepts
arbitrary ``extra=`` entries for exactly this kind of engine-supplied fact
(pass ``extra={"registered_lowerings": provenance_entries(resolutions)}``),
so stamping needs no change to the provenance schema or
``PROVENANCE_SCHEMA_VERSION`` — see the W2.6 report for the alternative (a
first-class ``ampere_registered_lowerings`` key baked into
``provenance_attrs`` itself) flagged there for review.

Worked example
---------------
>>> import numpy as np, scipy.stats as st
>>> from ampere.core.lowering import (
...     register_lowering, lookup_lowering, run_registrant_battery,
...     provenance_entries,
... )

A prior family scipy has but a hypothetical "stub" backend does not.
Register a constructor that builds the backend's own object from the
canonical spec — here, standing in for a real backend, a tiny class computing
the exact normal log-density from ``loc``/``scale``:

>>> class _StubNormal:
...     def __init__(self, loc, scale):
...         self.loc, self.scale = loc, scale
...     def log_prob(self, x):
...         z = (np.asarray(x) - self.loc) / self.scale
...         return -0.5 * z**2 - np.log(self.scale) - 0.5 * np.log(2 * np.pi)
>>> def _build_normal(spec):
...     return _StubNormal(spec.kwds["loc"], spec.kwds["scale"])
>>> register_lowering("norm", "stub", _build_normal)  # doctest: +ELLIPSIS
LoweringResolution(kind='prior', name='norm', backend='stub', ...)

Registering the same ``(family, backend)`` again without ``override=True``
is refused, naming what it refused to overwrite:

>>> register_lowering("norm", "stub", _build_normal)  # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
ampere.core.exceptions.LoweringError: prior family 'norm' cannot be lowered...
>>> _ = register_lowering("norm", "stub", _build_normal, override=True)  # explicit, so it proceeds

A registrant self-certifies before trusting the row for sampling —
reference (scipy) versus native, at points spanning the support:

>>> report = run_registrant_battery(st.norm(2.0, 0.5), "stub")
>>> report.passed
True

A wrong registration fails loudly, naming the point and both values:

>>> def _wrong_scale(spec):
...     return _StubNormal(spec.kwds["loc"], spec.kwds["scale"] * 2.0)
>>> _ = register_lowering("norm", "wrong-stub", _wrong_scale)
>>> run_registrant_battery(st.norm(2.0, 0.5), "wrong-stub")  # doctest: +ELLIPSIS
Traceback (most recent call last):
    ...
ampere.core.exceptions.LoweringError: prior family 'norm' cannot be lowered...

The resolution the registry handed back is exactly what a caller stamps
into a run's provenance record
(``ampere.results.provenance.provenance_attrs``'s ``extra=``); built-in rows
are omitted, since only genuinely user-registered lowerings need flagging:

>>> resolution = lookup_lowering("norm", "stub")
>>> resolution.builtin
False
>>> entries = provenance_entries([resolution])
>>> entries[0]["name"], entries[0]["backend"], entries[0]["builtin"]
('norm', 'stub', False)
"""

from __future__ import annotations

import dataclasses
from collections.abc import Callable, Iterable, Sequence
from typing import Any

import numpy as np

from .exceptions import LoweringError, ParameterError
from .parameter import AnyPrior, Bijection, PriorSpec, describe_prior, log_density

__all__ = [
    "BatteryReport",
    "LoweringResolution",
    "lookup_bijection_lowering",
    "lookup_lowering",
    "provenance_entries",
    "register_bijection_lowering",
    "register_lowering",
    "registered_lowerings",
    "run_registrant_battery",
]

_Kind = str  # "prior" | "bijection" -- not an enum: kept as the plain strings
# LoweringResolution.kind and every message already use, so nothing has to
# convert between an enum and its value at the message-building boundary.

_PRIOR_KIND: _Kind = "prior"
_BIJECTION_KIND: _Kind = "bijection"


@dataclasses.dataclass(frozen=True)
class LoweringResolution:
    """One resolved ``(kind, name, backend)`` row: a constructor plus its provenance.

    Returned by :func:`register_lowering`/:func:`register_bijection_lowering`
    (what was just stored) and by :func:`lookup_lowering`/
    :func:`lookup_bijection_lowering` (what was found). ``constructor`` is
    called with a :class:`~ampere.core.parameter.PriorSpec` for a ``"prior"``
    row or a :class:`~ampere.core.parameter.Bijection` instance for a
    ``"bijection"`` row.
    """

    kind: _Kind
    name: str
    backend: str
    constructor: Callable[[Any], Any]
    builtin: bool
    constructor_module: str
    constructor_qualname: str

    def to_provenance_entry(self) -> dict[str, Any]:
        """A JSON-plain, netCDF-safe record of this row for a run's provenance."""
        return {
            "kind": self.kind,
            "name": self.name,
            "backend": self.backend,
            "builtin": self.builtin,
            "constructor": f"{self.constructor_module}.{self.constructor_qualname}",
        }


#: The one store both slots share (see the module docstring, "One store, two
#: slots"). Keyed on ``(kind, name, backend)``; values are the resolutions
#: themselves, since a resolution already carries everything a lookup needs.
_REGISTRY: dict[tuple[_Kind, str, str], LoweringResolution] = {}


def _register(
    kind: _Kind,
    name: str,
    backend: str,
    constructor: Callable[[Any], Any],
    *,
    override: bool,
    builtin: bool,
    what: str,
) -> LoweringResolution:
    if not isinstance(name, str) or not name:
        raise ParameterError(f"a {what} name must be a non-empty string, got {name!r}")
    if not isinstance(backend, str) or not backend:
        raise ParameterError(f"a backend name must be a non-empty string, got {backend!r}")
    if not callable(constructor):
        raise ParameterError(f"a lowering constructor must be callable, got {constructor!r}")
    key = (kind, name, backend)
    existing = _REGISTRY.get(key)
    if existing is not None and not override:
        origin = "a built-in" if existing.builtin else "a user"
        raise LoweringError(
            name,
            backend=backend,
            detail=(
                f"{origin} lowering is already registered for {what} {name!r} on backend "
                f"{backend!r} ({existing.constructor_module}.{existing.constructor_qualname}); "
                f"pass override=True to replace it deliberately."
            ),
        )
    resolution = LoweringResolution(
        kind=kind,
        name=name,
        backend=backend,
        constructor=constructor,
        builtin=builtin,
        constructor_module=getattr(constructor, "__module__", "<unknown>"),
        constructor_qualname=getattr(constructor, "__qualname__", repr(constructor)),
    )
    _REGISTRY[key] = resolution
    return resolution


def _lookup(kind: _Kind, name: str, backend: str, *, what: str) -> LoweringResolution:
    try:
        return _REGISTRY[(kind, name, backend)]
    except KeyError:
        raise LoweringError(
            name,
            backend=backend,
            detail=f"no lowering is registered for {what} {name!r} on backend {backend!r}",
        ) from None


def register_lowering(
    family: str,
    backend: str,
    constructor: Callable[[PriorSpec], Any],
    *,
    override: bool = False,
    builtin: bool = False,
) -> LoweringResolution:
    """Register a backend-native constructor for a prior family.

    ``family`` is the neutral name :func:`~ampere.core.parameter.describe_prior`
    reads off a frozen ``scipy.stats`` distribution
    (:attr:`~ampere.core.parameter.PriorSpec.family`) — the same name
    ``lowering.md`` §3.2's built-in table keys on, so a registration for
    ``"truncnorm"`` on ``"torch"`` is exactly what closes that table's gap.
    ``constructor`` receives the :class:`~ampere.core.parameter.PriorSpec`
    (never the raw scipy object — ``lowering.md`` §1.5, a backend "works from
    ``PriorSpec`` throughout") and must return the backend-native distribution
    object; see the module docstring's trace-purity rule.

    Parameters
    ----------
    family, backend
        The registry key.
    constructor
        ``PriorSpec -> native distribution``. Must be trace-pure (module
        docstring).
    override
        Required to replace a row that already exists (built-in or
        user-registered); otherwise this raises, naming what it refused to
        overwrite.
    builtin
        Reserved for the backends ampere itself ships (W2.4/W2.5's own
        registrations of ``lowering.md`` §3.2's table): a built-in row still
        requires ``override=True`` to replace, and is never counted as
        "user-registered" by :func:`provenance_entries`. Ordinary callers
        should not pass this.

    Raises
    ------
    LoweringError
        If a row already exists at this key and ``override`` is not ``True``.
    """
    return _register(
        _PRIOR_KIND,
        family,
        backend,
        constructor,
        override=override,
        builtin=builtin,
        what="prior family",
    )


def lookup_lowering(family: str, backend: str) -> LoweringResolution:
    """The registered resolution for ``family`` on ``backend``.

    Raises
    ------
    LoweringError
        Naming the family and the backend, with the standard three remedies
        (change the prior, register a lowering, or run on a backend that has
        it — see :class:`~ampere.core.exceptions.LoweringError`).
    """
    return _lookup(_PRIOR_KIND, family, backend, what="prior family")


def _looks_like_bijection(cls: type) -> bool:
    return all(hasattr(cls, attr) for attr in ("constrain", "unconstrain", "log_abs_det_jacobian"))


def register_bijection_lowering(
    bijection: type[Bijection] | Bijection,
    backend: str,
    constructor: Callable[[Bijection], Any],
    *,
    override: bool = False,
    builtin: bool = False,
) -> LoweringResolution:
    """Register a backend-native constructor for a custom :class:`Bijection`.

    The analogous per-backend slot ``lowering.md`` §12.8 asks for alongside
    :func:`register_lowering`: §4's table row for "user-supplied custom
    Bijection" is "unsupported — raises at lowering" on torch and numpyro;
    this closes it for a specific class. ``bijection`` may be the class
    itself or an instance (an instance is a convenience — only its type is
    used as the key, since a :class:`Bijection`'s neutral identity is its
    class, the same convention
    :class:`~ampere.core.exceptions.LoweringError`'s docstring already uses
    ("a bijection class name")). ``constructor`` receives the *instance*
    being lowered (its ``lower``/``upper`` or other fields may matter) and
    must return the backend-native transform object.

    Parameters, raises
    -------------------
    As :func:`register_lowering`, with ``bijection`` in place of ``family``.
    """
    cls = bijection if isinstance(bijection, type) else type(bijection)
    if not isinstance(cls, type) or not _looks_like_bijection(cls):
        raise ParameterError(
            f"register_bijection_lowering expects a Bijection class or instance, got {bijection!r}"
        )
    return _register(
        _BIJECTION_KIND,
        cls.__name__,
        backend,
        constructor,
        override=override,
        builtin=builtin,
        what="bijection class",
    )


def lookup_bijection_lowering(
    bijection: type[Bijection] | Bijection, backend: str
) -> LoweringResolution:
    """The registered resolution for a custom :class:`Bijection` class on ``backend``.

    Raises
    ------
    LoweringError
        Naming the bijection's class and the backend.
    """
    cls = bijection if isinstance(bijection, type) else type(bijection)
    return _lookup(_BIJECTION_KIND, cls.__name__, backend, what="bijection class")


def registered_lowerings(
    *, kind: _Kind | None = None, backend: str | None = None
) -> tuple[LoweringResolution, ...]:
    """Every registered row, optionally filtered by kind and/or backend.

    Introspection only — not part of the lowering path itself. Useful for a
    backend enumerating what it has, or a test asserting a registration
    landed.
    """
    return tuple(
        resolution
        for resolution in _REGISTRY.values()
        if (kind is None or resolution.kind == kind)
        and (backend is None or resolution.backend == backend)
    )


def provenance_entries(resolutions: Iterable[LoweringResolution]) -> list[dict[str, Any]]:
    """Canonical, netCDF-safe records of the *user-registered* rows a run consulted.

    ``lowering.md`` §12.8's hardening: "every registered row is stamped
    user-registered in provenance." Built-in rows are omitted — they are not
    what the ruling means by "user-registered," and stamping them would bury
    the signal a reviewer actually wants (did this run depend on something
    outside ampere's own conformance guarantees?) under every ordinary row.

    Pass the result straight through to
    ``ampere.results.provenance.provenance_attrs``'s ``extra=`` (see the
    module docstring's "Provenance stamping" section for why that needs no
    change to the provenance schema)::

        provenance_attrs(problem, extra={"registered_lowerings": provenance_entries(used)})
    """
    return [
        resolution.to_provenance_entry() for resolution in resolutions if not resolution.builtin
    ]


@dataclasses.dataclass(frozen=True)
class BatteryReport:
    """The result of one :func:`run_registrant_battery` self-certification run."""

    family: str
    backend: str
    points: tuple[float, ...]
    max_abs_error: float
    tolerance: float

    @property
    def passed(self) -> bool:
        """Whether every compared point agreed within ``tolerance``."""
        return self.max_abs_error <= self.tolerance


def run_registrant_battery(
    prior: AnyPrior,
    backend: str,
    *,
    points: Sequence[float] | None = None,
    tolerance: float = 1e-9,
    to_numpy: Callable[[Any], np.ndarray] = np.asarray,
    log_prob: Callable[[Any, np.ndarray], Any] | None = None,
) -> BatteryReport:
    """Self-certify a registered prior-family lowering: reference-vs-native agreement.

    ``lowering.md`` §12.8's hardening turns "a user-registered lowering
    bypasses the conformance suite's guarantees" into "can self-certify": a
    registrant runs this against their own row before trusting it for
    sampling. It is deliberately **not** wired into ``tests/conformance`` —
    that suite runs once per *ampere-shipped* backend fixture
    (``tests/conformance/backends/__init__.py``'s registry); a third party's
    custom family has no fixture there, and forcing one would misrepresent an
    opt-in self-check as ampere's own guarantee.

    Compares :func:`~ampere.core.parameter.log_density` (the reference,
    scipy-backed density — always right, because the declaration *is* the
    implementation, ``lowering.md`` §0) against the registered constructor's
    own density at a handful of points spanning the prior's support, the same
    comparison ``lowering.md`` §11 row 3 makes for the built-in table. The
    default points are the prior's own 10th/25th/50th/75th/90th percentiles
    (``prior.ppf``), which stay inside the support for any family.

    The lowered object is expected to expose ``log_prob(x)`` returning
    something ``to_numpy`` can convert — the ``torch.distributions`` /
    ``numpyro.distributions`` convention both real backends will use — but
    ``log_prob=`` overrides this for a native object with a different API.

    Parameters
    ----------
    prior
        The frozen prior to check (passed to
        :func:`~ampere.core.parameter.describe_prior` and to ``prior.ppf``
        directly).
    backend
        Which registered row to check —
        :func:`lookup_lowering` (``describe_prior(prior).family``, ``backend``)
        supplies the constructor.
    points
        Evaluation points; defaults to five quantiles spanning the support.
    tolerance
        Maximum tolerated absolute log-density difference at any point.
    to_numpy
        Converts the native ``log_prob`` output to a plain numpy array —
        override for a backend whose native type ``np.asarray`` cannot read
        directly (a torch/jax array typically can be, via ``__array__``).
    log_prob
        ``(native, points) -> array-like``, overriding the ``.log_prob``
        convention above.

    Returns
    -------
    BatteryReport
        On agreement within tolerance.

    Raises
    ------
    LoweringError
        If any point disagrees beyond ``tolerance``, naming the point and
        both values — or if the native object has no ``log_prob`` and
        ``log_prob=`` was not supplied.
    """
    spec = describe_prior(prior)
    resolution = lookup_lowering(spec.family, backend)
    native = resolution.constructor(spec)
    if points is None:
        quantiles = np.array([0.10, 0.25, 0.50, 0.75, 0.90])
        xs = np.asarray(prior.ppf(quantiles), dtype=float)
    else:
        xs = np.asarray(points, dtype=float)
    reference = np.asarray(log_density(prior, xs), dtype=float)
    if log_prob is not None:
        native_values = np.asarray(to_numpy(log_prob(native, xs)), dtype=float)
    else:
        native_log_prob = getattr(native, "log_prob", None)
        if native_log_prob is None:
            raise LoweringError(
                spec.family,
                backend=backend,
                detail=(
                    "registrant battery failed: the registered constructor's output has no "
                    ".log_prob(x) method (the torch.distributions/numpyro.distributions "
                    "convention this battery assumes by default); pass log_prob= explicitly "
                    "for a different API."
                ),
            )
        native_values = np.asarray(to_numpy(native_log_prob(xs)), dtype=float)
    errors = np.abs(native_values - reference)
    max_error = float(np.max(errors))
    report = BatteryReport(
        family=spec.family,
        backend=backend,
        points=tuple(float(x) for x in xs),
        max_abs_error=max_error,
        tolerance=tolerance,
    )
    if not report.passed:
        worst = int(np.argmax(errors))
        raise LoweringError(
            spec.family,
            backend=backend,
            detail=(
                f"registrant battery failed: at x={xs[worst]!r} the registered lowering's "
                f"log_prob gives {native_values[worst]!r} but the reference (scipy) gives "
                f"{reference[worst]!r} (|diff|={errors[worst]!r} > tolerance={tolerance!r}). "
                f"Fix the constructor registered for backend {backend!r}, or re-register it."
            ),
        )
    return report
