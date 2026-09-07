"""``lowering.md`` §4's bijection table, on numpyro.

Three ampere bijections, three numpyro transforms, and one direction convention
that is the named silent-bug class these contracts exist to prevent.

The direction, written out (§2)
-------------------------------
``biject_to(constraint)`` returns a transform whose **forward** direction maps
the unconstrained real line onto the constrained support. That is the same
direction as ampere's ``constrain``::

    ampere:  x = bijection.constrain(y)
    numpyro: x = biject_to(constraint)(y)

and the inverse direction is ``bijection.unconstrain(x)`` ↔ ``t.inv(x)``.

The log-determinant is where the bug lives. numpyro's signature is
``log_abs_det_jacobian(x, y)`` where ``x`` is the *forward transform's input* —
the unconstrained point — and ``y`` its output. Ampere's convention, following
the statistics literature and ``parameters.md`` §6, is the opposite: ``x`` is
constrained and ``y`` unconstrained. **The two libraries' ``x`` is ampere's
``y``.** So the correct call, and the only one that matches ampere's reference
semantics, is::

    ladj = t.log_abs_det_jacobian(y_unconstrained, x_constrained)

Passing the arguments the other way round does not crash, does not produce
NaNs, and does not change any shape or dtype. It produces a subtly wrong
posterior that looks plausible. :func:`log_abs_det_jacobian` below is the one
place in this package that makes the call, so there is one place to get it
right and one place to check.

``biject_to``, never ``transform_to`` (§2 b')
---------------------------------------------
Both registries map the real line to a constraint in the same direction, and
for simple constraints they often return the same object — so a substitution
passes casual inspection. They differ in guarantees: ``biject_to`` is
guaranteed bijective and guaranteed to implement ``log_abs_det_jacobian``;
``transform_to`` guarantees neither and exists for unconstrained optimisation.
Ampere needs the Jacobian on every path reaching ``lnprior_unconstrained``, so
the answer is always ``biject_to``. The name ``transform_to`` does not appear
in this package.

Why the constraint and not the distribution
-------------------------------------------
§4 recommends ``biject_to(lowered_distribution.support)`` over rebuilding the
compose chain by hand, and that recommendation is followed — but the constraint
is built from **ampere's own** :class:`~ampere.core.parameter.Bijection`, not
read off the lowered distribution. Two reasons, and the second is a fact about
numpyro rather than a preference:

1. a user may *declare* a bijection that is not the one inferred from the
   prior's support (``Identity`` over a bounded prior is a legitimate
   fixed-support reparametrisation, and is a conformance row), and reading the
   distribution's support would silently ignore the declaration;
2. a distribution built by §3.3's rule 1 is a ``TransformedDistribution``, and
   its ``support`` is only honest if the affine transform was given an explicit
   ``domain=``. Driving the bijection off ampere's declaration keeps this
   module independent of that.

Either way the library stays the authority on its own constraint registry —
the compose chains are numpyro's, not ampere's — and §4's table then documents
what comes back and serves as the cross-check. **Do not pattern-match on
transform type**: ``constraints.positive`` and ``constraints.greater_than(0)``
are different constraint *types* registered separately, so two priors with the
same support can hand back different objects that behave identically.
"""

from __future__ import annotations

from typing import Any

import jax.numpy as jnp
from numpyro.distributions import constraints
from numpyro.distributions.transforms import Transform, biject_to

from ampere.core.lowering import (
    lookup_bijection_lowering,
    register_bijection_lowering,
    registered_lowerings,
)
from ampere.core.parameter import Bijection, Identity, Log, Logit

from ._config import BACKEND

__all__ = [
    "BUILTIN_BIJECTIONS",
    "log_abs_det_jacobian",
    "lower_bijection",
    "register_builtin_bijection_lowerings",
]


def _identity(_: Identity) -> Transform:
    """``x = y`` → ``biject_to(constraints.real)``."""
    return biject_to(constraints.real)


def _log(bijection: Log) -> Transform:
    """``x = l + exp(y)`` → a bare ``ExpTransform`` at ``l = 0``, composed otherwise.

    ``constraints.positive`` and ``constraints.greater_than(0.0)`` are
    mathematically the same support and numpyro registers them separately, so
    asking for the first at ``l = 0`` is what produces §4's documented bare
    ``ExpTransform`` rather than a one-element compose chain. Behaviourally the
    two agree; the object differs, which is the note §4 attaches to this row.
    """
    lower = float(bijection.lower)
    if lower == 0.0:
        return biject_to(constraints.positive)
    return biject_to(constraints.greater_than(lower))


def _logit(bijection: Logit) -> Transform:
    """``x = l + (h - l)*sigma(y)``: a bare ``SigmoidTransform`` on the unit
    interval, composed with an affine map otherwise."""
    lower, upper = float(bijection.lower), float(bijection.upper)
    if (lower, upper) == (0.0, 1.0):
        return biject_to(constraints.unit_interval)
    return biject_to(constraints.interval(lower, upper))


#: ``lowering.md`` §4's table: ampere bijection class → numpyro constructor.
BUILTIN_BIJECTIONS: dict[type, Any] = {
    Identity: _identity,
    Log: _log,
    Logit: _logit,
}


def register_builtin_bijection_lowerings(*, override: bool = True) -> None:
    """Register §4's three rows in ``ampere.core.lowering``'s bijection slot.

    The same discipline as the prior table: the built-in rows go through the
    public registry, so a user registering a lowering for a custom
    :class:`~ampere.core.parameter.Bijection` reaches this backend by exactly
    the route ampere's own rows do.
    """
    for cls, constructor in BUILTIN_BIJECTIONS.items():
        register_bijection_lowering(cls, BACKEND, constructor, builtin=True, override=override)


def _ensure_registered() -> None:
    present = {row.name for row in registered_lowerings(kind="bijection", backend=BACKEND)}
    expected = {f"{cls.__module__}.{cls.__qualname__}" for cls in BUILTIN_BIJECTIONS}
    if not present.issuperset(expected):
        register_builtin_bijection_lowerings()


def lower_bijection(bijection: Bijection) -> Transform:
    """The numpyro transform *bijection* declares.

    A custom :class:`~ampere.core.parameter.Bijection` — one satisfying
    ampere's protocol with numpy operations a jax tracer will not accept — is a
    reference-path-only feature (§4's last row) unless its author has registered
    a lowering for it. Either way the answer comes from the registry, so an
    unregistered custom class raises
    :class:`~ampere.core.exceptions.LoweringError` naming the class and this
    backend rather than failing later inside a trace.
    """
    _ensure_registered()
    return lookup_bijection_lowering(type(bijection), BACKEND).constructor(bijection)


def log_abs_det_jacobian(transform: Transform, unconstrained: Any, constrained: Any) -> Any:
    """``log|d constrain / dy|`` at the **unconstrained** point.

    The one call site in this package, so the argument order of §2(b) is
    written down once. *unconstrained* is the transform's input and
    *constrained* its output; numpyro's parameter names for the two are ``x``
    and ``y``, which are ampere's ``y`` and ``x``.
    """
    return transform.log_abs_det_jacobian(jnp.asarray(unconstrained), jnp.asarray(constrained))


register_builtin_bijection_lowerings()
