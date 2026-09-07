"""``lowering.md`` §3's distribution table, on numpyro.

Every row of §3.2 — both tiers — plus §3.3's parametrisation conversions and
§3.4's fallback rule, registered into ``ampere.core.lowering``'s registry as
built-in rows keyed ``("prior", family, "jax")``. Nothing here is looked up by
name from anywhere else in this package: :func:`lower_prior` goes through
:func:`~ampere.core.lowering.lookup_lowering`, exactly as a third party's own
registration would be, so the built-in and the user-supplied paths cannot
diverge.

Why the conversions matter, restated once
-----------------------------------------
``scipy`` gives every continuous family a ``loc``/``scale`` shift-and-stretch
that numpyro's equivalents do not have, and dropping one silently moves a
prior's support. ``lowering.md`` §3.3's rule is therefore that a non-default
``loc`` the target family cannot express is either **composed away exactly**
(rule 1) or **raises** (rule 2) — never dropped, never approximated. The
conversions that change *numbers* rather than merely move them are worth naming
here because each one passes casual inspection:

``uniform``
    scipy's second argument is a **width**: ``high = loc + scale``.
``lognorm``
    scipy parametrises by the shape ``s`` and ``scale = exp(mu)``; numpyro's
    ``LogNormal`` takes ``loc = mu`` directly, so ``loc_target = log(scale)``.
    Passing ``scale`` through is wrong by an exponential.
``truncnorm``
    scipy's ``a``/``b`` are in **standardised** units — the support is
    ``[loc + a·scale, loc + b·scale]``. numpyro's ``low``/``high`` are not, and
    they are keyword-only. Passing ``a``/``b`` straight through gives a prior
    truncated around the wrong place, which reads as a badly behaved sampler
    rather than as a translation bug.
``halfnorm`` with ``loc ≠ 0``
    Lowered to ``TruncatedNormal(loc, scale, low=loc)`` rather than to a
    ``TransformedDistribution``: a normal centred at ``l`` and truncated below
    at ``l`` *is* a half-normal shifted to ``l``, exactly (the truncation
    renormalises by two, which is the half-normal's factor), and the result's
    ``support`` is ``greater_than(l)`` so numpyro's own ``biject_to`` yields
    the right bijection without anything being reasoned about by hand.

Composed shifts carry their domain
----------------------------------
Where rule 1 does apply (a shifted ``expon``, ``gamma`` or ``lognorm``; a
``beta`` on an interval other than the unit one) the affine transform is
constructed with an explicit ``domain=``. This is not decoration: numpyro's
``AffineTransform`` defaults its domain to ``real``, so
``TransformedDistribution(Exponential(...), AffineTransform(1.5, 1.0))``
reports ``support = Real()`` — a lie about a distribution supported on
``[1.5, ∞)``, and one that would give the wrong ``biject_to``. Verified against
numpyro 0.21.0: with ``domain=constraints.positive`` the same object reports
``GreaterThan(lower_bound=1.5)``.

Discrete families (§3.5)
------------------------
``poisson`` lowers cleanly as a *distribution* and is registered here. It does
not lower onto a gradient path, and it does not need to be refused here:
``ampere.core.default_bijection_for`` already raises ``CapabilityError`` for a
discrete family, so a discrete parameter cannot acquire an unconstraining
bijection at all. Lowering-as-a-distribution — which nested sampling and an SBI
simulator both use — stays open, which is exactly the door that ruling
(2026-09-03) kept ajar.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import jax.numpy as jnp
import numpyro.distributions as npd
from numpyro.distributions import constraints
from numpyro.distributions.transforms import AffineTransform

from ampere.core.exceptions import LoweringError
from ampere.core.lowering import register_lowering, registered_lowerings
from ampere.core.parameter import PriorSpec, describe_prior, prior_from_spec

from ._config import BACKEND

__all__ = [
    "FAMILIES",
    "NATIVE_ICDF",
    "has_native_icdf",
    "lower_prior",
    "register_builtin_lowerings",
]

#: Families for which numpyro's ``icdf`` is expected to work (``lowering.md``
#: §3.6, as amended by what this backend actually observed — see
#: :func:`has_native_icdf`). Advisory only: the decision is taken by probing the
#: constructed object, because a declared method that raises is not an ``icdf``.
NATIVE_ICDF: frozenset[str] = frozenset(
    {"norm", "uniform", "halfnorm", "loguniform", "lognorm", "expon", "truncnorm"}
)


# ---------------------------------------------------------------------------
# Argument plumbing
# ---------------------------------------------------------------------------


def _canonical(spec: PriorSpec) -> PriorSpec:
    """*spec* with every positional argument mapped onto its keyword name.

    ``describe_prior`` already emits the canonical keyword-only form
    (``lowering.md`` §3.1), so this is a no-op for a described prior and a
    normalisation for a hand-written :class:`~ampere.core.parameter.PriorSpec`
    built positionally. Round-tripping through ``scipy`` is the cheapest way to
    get scipy's own argument-name ordering without copying it.
    """
    if not spec.args:
        return spec
    return describe_prior(prior_from_spec(spec))


def _argument(values: Mapping[str, Any], name: str, default: float) -> Any:
    """One distribution argument, defaulted. May be an array (hierarchical)."""
    found = values.get(name)
    return default if found is None else found


def _structural(value: Any, argument: str, family: str) -> float:
    """A distribution argument whose *value* decides the shape of the lowering.

    ``halfnorm``'s ``loc``, ``expon``'s ``loc``, ``beta``'s ``loc``/``scale``:
    whether the lowering is the bare family or a composed one depends on
    whether these are at their defaults, so they have to be *known numbers* at
    lowering time. A :class:`~ampere.core.parameter.HierarchicalPrior` that
    referenced one would be asking for a branch on a value that does not exist
    until the sampler proposes it — the lowering would have to be re-decided
    every evaluation, which is the thing ``lowering.md`` §8 says a hierarchical
    prior may not be. It is refused by name rather than guessed at.
    """
    try:
        return float(value)
    except Exception:  # a tracer, an array, anything not a scalar
        raise LoweringError(
            family,
            backend=BACKEND,
            detail=(
                f"the {argument!r} argument decides the shape of this family's lowering (whether "
                f"it is the bare numpyro family or a composed one), so it must be a known number "
                f"at lowering time; it was given as {value!r}, which is not. A HierarchicalPrior "
                f"may reference the arguments that are pure numbers ('loc' and 'scale' of "
                f"'norm', say), not the ones that change which distribution object is built."
            ),
        ) from None


def _shifted(base: npd.Distribution, loc: Any, domain: Any) -> npd.Distribution:
    """*base* translated by *loc* — ``lowering.md`` §3.3's rule 1, exactly.

    ``domain=`` is mandatory rather than optional: without it numpyro's
    ``AffineTransform`` reports a real-line codomain and the composed
    distribution's ``support`` becomes a lie (see the module docstring).
    """
    return npd.TransformedDistribution(base, AffineTransform(loc=loc, scale=1.0, domain=domain))


# ---------------------------------------------------------------------------
# The table (lowering.md §3.2)
# ---------------------------------------------------------------------------


def _norm(spec: PriorSpec) -> npd.Distribution:
    """``norm(loc, scale)`` → ``Normal(loc, scale)`` — the identity row."""
    kwds = spec.kwds
    return npd.Normal(_argument(kwds, "loc", 0.0), _argument(kwds, "scale", 1.0))


def _uniform(spec: PriorSpec) -> npd.Distribution:
    """``uniform(loc, scale)`` → ``Uniform(low=loc, high=loc + scale)``.

    scipy's second argument is a **width**. The row most likely to be got right
    by accident and wrong under maintenance (§3.3).
    """
    kwds = spec.kwds
    low = _argument(kwds, "loc", 0.0)
    width = _argument(kwds, "scale", 1.0)
    return npd.Uniform(low, low + width)


def _halfnorm(spec: PriorSpec) -> npd.Distribution:
    """``halfnorm(loc, scale)`` → ``HalfNormal(scale)``, or the exact shift."""
    kwds = spec.kwds
    scale = _argument(kwds, "scale", 1.0)
    loc = _structural(_argument(kwds, "loc", 0.0), "loc", "halfnorm")
    if loc == 0.0:
        return npd.HalfNormal(scale)
    # A normal centred at `loc` and truncated below at `loc` is exactly a
    # half-normal shifted to `loc` (§3.3). Preferred over the affine route
    # because its support is greater_than(loc), so biject_to is right for free.
    return npd.TruncatedNormal(loc=loc, scale=scale, low=loc)


def _loguniform(spec: PriorSpec) -> npd.Distribution:
    """``loguniform(a, b)`` → ``LogUniform(low=a, high=b)`` — native on numpyro."""
    kwds = spec.kwds
    if _structural(_argument(kwds, "loc", 0.0), "loc", "loguniform") != 0.0:
        raise LoweringError(
            "loguniform",
            backend=BACKEND,
            detail="numpyro's LogUniform has no location parameter, and a shifted log-uniform "
            "is not an affine transform of an unshifted one, so no exact construction exists.",
        )
    return npd.LogUniform(_argument(kwds, "a", 1.0), _argument(kwds, "b", 2.0))


def _poisson(spec: PriorSpec) -> npd.Distribution:
    """``poisson(mu)`` → ``Poisson(rate=mu)``. A distribution, not a sampler path (§3.5)."""
    kwds = spec.kwds
    if _structural(_argument(kwds, "loc", 0.0), "loc", "poisson") != 0.0:
        raise LoweringError(
            "poisson",
            backend=BACKEND,
            detail="numpyro's Poisson has no location parameter; a shifted Poisson is a "
            "different distribution on the integers, not an affine transform of this one.",
        )
    return npd.Poisson(_argument(kwds, "mu", 1.0))


def _truncnorm(spec: PriorSpec) -> npd.Distribution:
    """``truncnorm(a, b, loc, scale)`` → ``TruncatedNormal``, **unstandardised**.

    scipy's ``a``/``b`` are in units of ``scale`` about ``loc``; numpyro's
    ``low``/``high`` are in the data's own units and are keyword-only. An
    infinite bound is passed as ``None``, which is how numpyro spells
    "untruncated on that side".
    """
    kwds = spec.kwds
    loc = _argument(kwds, "loc", 0.0)
    scale = _argument(kwds, "scale", 1.0)
    a = _structural(_argument(kwds, "a", -math.inf), "a", "truncnorm")
    b = _structural(_argument(kwds, "b", math.inf), "b", "truncnorm")
    low = None if a == -math.inf else loc + a * scale
    high = None if b == math.inf else loc + b * scale
    return npd.TruncatedNormal(loc=loc, scale=scale, low=low, high=high)


def _lognorm(spec: PriorSpec) -> npd.Distribution:
    """``lognorm(s, loc, scale)`` → ``LogNormal(loc=log(scale), scale=s)``.

    Wrong by an exponential if ``scale`` is passed through unchanged. A
    non-zero ``loc`` is a pure translation of the support, so §3.3's rule 1
    composes it away exactly.
    """
    kwds = spec.kwds
    shape = _argument(kwds, "s", 1.0)
    scale = _argument(kwds, "scale", 1.0)
    base = npd.LogNormal(loc=jnp.log(jnp.asarray(scale)), scale=shape)
    loc = _structural(_argument(kwds, "loc", 0.0), "loc", "lognorm")
    return base if loc == 0.0 else _shifted(base, loc, constraints.positive)


def _expon(spec: PriorSpec) -> npd.Distribution:
    """``expon(loc, scale)`` → ``Exponential(rate=1/scale)``, shifted if need be."""
    kwds = spec.kwds
    scale = _argument(kwds, "scale", 1.0)
    base = npd.Exponential(rate=1.0 / jnp.asarray(scale))
    loc = _structural(_argument(kwds, "loc", 0.0), "loc", "expon")
    return base if loc == 0.0 else _shifted(base, loc, constraints.positive)


def _gamma(spec: PriorSpec) -> npd.Distribution:
    """``gamma(a, loc, scale)`` → ``Gamma(concentration=a, rate=1/scale)``."""
    kwds = spec.kwds
    concentration = _argument(kwds, "a", 1.0)
    scale = _argument(kwds, "scale", 1.0)
    base = npd.Gamma(concentration=concentration, rate=1.0 / jnp.asarray(scale))
    loc = _structural(_argument(kwds, "loc", 0.0), "loc", "gamma")
    return base if loc == 0.0 else _shifted(base, loc, constraints.positive)


def _beta(spec: PriorSpec) -> npd.Distribution:
    """``beta(a, b, loc, scale)`` → ``Beta``, affinely mapped onto ``[loc, loc+scale]``.

    §3.2's row says ``loc = 0`` and ``scale = 1``; §3.3's rule 1 says compose an
    exactly-equivalent affine map away rather than raise, and for a Beta the
    ``loc``/``scale`` pair *is* exactly an affine map of the unit interval. The
    rule beats the row here — the row records what numpyro's ``Beta`` itself
    takes — and the result is a distribution whose support is
    ``interval(loc, loc + scale)``, which is what scipy's is.
    """
    kwds = spec.kwds
    base = npd.Beta(_argument(kwds, "a", 1.0), _argument(kwds, "b", 1.0))
    loc = _structural(_argument(kwds, "loc", 0.0), "loc", "beta")
    scale = _structural(_argument(kwds, "scale", 1.0), "scale", "beta")
    if (loc, scale) == (0.0, 1.0):
        return base
    return npd.TransformedDistribution(
        base, AffineTransform(loc=loc, scale=scale, domain=constraints.unit_interval)
    )


#: The §3.2 table itself: neutral family name → constructor. Registered into
#: ``ampere.core.lowering`` as built-in rows by
#: :func:`register_builtin_lowerings`, which this module calls at import.
FAMILIES: Mapping[str, Callable[[PriorSpec], npd.Distribution]] = {
    # Required tier -- the families the Phase 1 specs actually use.
    "norm": _norm,
    "uniform": _uniform,
    "halfnorm": _halfnorm,
    "loguniform": _loguniform,
    "poisson": _poisson,
    "truncnorm": _truncnorm,
    # Second tier -- not used by the Phase 1 specs, but common enough that
    # leaving them to be rediscovered per backend invites divergence.
    "lognorm": _lognorm,
    "expon": _expon,
    "gamma": _gamma,
    "beta": _beta,
}


def register_builtin_lowerings(*, override: bool = True) -> None:
    """Put every :data:`FAMILIES` row into ``ampere.core.lowering``'s registry.

    Called once at import, and again (cheaply) by :func:`lower_prior` if a row
    it needs has gone missing. That second call is not defensive
    over-engineering: ``tests/core/conftest.py`` snapshots and restores the
    registry around every test in that directory precisely because it is
    module-global mutable state, and a package whose built-in rows can be
    removed by someone else's teardown should be able to put them back rather
    than fail a later run with "no lowering is registered".

    ``builtin=True`` on every row, so :func:`~ampere.core.lowering.provenance_entries`
    omits them: they are ampere's own table, covered by the conformance suite,
    and stamping them would bury the signal a reviewer wants (did this run
    depend on something *outside* that guarantee?) under ten ordinary rows.
    """
    for family, constructor in FAMILIES.items():
        register_lowering(family, BACKEND, constructor, builtin=True, override=override)


def _ensure_registered() -> None:
    """Re-register the built-in rows if any has been removed. See above."""
    present = {row.name for row in registered_lowerings(kind="prior", backend=BACKEND)}
    if not present.issuperset(FAMILIES):
        register_builtin_lowerings()


def lower_prior(spec: PriorSpec, *, parameter: str | None = None) -> npd.Distribution:
    """The numpyro distribution *spec* declares.

    The whole of §3, applied: the registry is consulted (so a user-registered
    row for a family ampere has no table entry for is found the same way a
    built-in one is), the constructor is run, and a family with no row raises
    :class:`~ampere.core.exceptions.LoweringError` naming the family, the
    parameter and this backend — **never** a silent fallback to the reference
    path (§3.4 rule 3).

    Parameters
    ----------
    spec
        The neutral prior description. Positional arguments are canonicalised
        onto their keyword names first, so a hand-written spec lowers the same
        way a described one does.
    parameter
        The merged parameter name, for the error message. A refusal that does
        not say *which* prior it was about is a refusal the user has to go
        hunting with.
    """
    _ensure_registered()
    canonical = _canonical(spec)
    from ampere.core.lowering import lookup_lowering

    try:
        resolution = lookup_lowering(canonical.family, BACKEND)
    except LoweringError as error:
        raise LoweringError(
            canonical.family,
            backend=BACKEND,
            parameter=parameter,
            detail=(
                f"{error.detail} — ampere's built-in table for this backend covers "
                f"{', '.join(sorted(FAMILIES))}."
            ),
        ) from None
    try:
        return resolution.constructor(canonical)
    except LoweringError as error:
        # Re-raise with the parameter attached when the constructor did not
        # know it (the registry hands constructors a spec, not a name).
        if error.parameter is not None or parameter is None:
            raise
        raise LoweringError(
            error.family, backend=error.backend, parameter=parameter, detail=error.detail
        ) from None


def has_native_icdf(distribution: npd.Distribution) -> bool:
    """Whether *distribution* can actually evaluate an inverse CDF.

    **Probed, not declared** (``lowering.md`` §3.6). numpyro declares ``icdf``
    on the base :class:`~numpyro.distributions.Distribution` and every family
    inherits the declaration, so ``hasattr`` answers ``True`` for families that
    cannot do it. Worse, on numpyro 0.21.0 ``Gamma.icdf`` and ``Beta.icdf``
    *are* implemented but delegate to TensorFlow Probability and raise
    ``ImportError`` when it is absent — so even "the method is overridden" is
    not the question. Whether the call works in **this** environment is.

    Called once per family at lowering time, never inside a traced function.
    """
    try:
        distribution.icdf(jnp.asarray(0.5))
    except Exception:  # NotImplementedError, ImportError, whatever numpyro raises
        return False
    return True


register_builtin_lowerings()
