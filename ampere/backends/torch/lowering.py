"""``lowering.md`` §3 and §4 for torch: prior families and bijections.

Lowering is the one-way translation of ``ampere.core``'s neutral declarations
into a specific backend's objects, and it happens **once**, when a problem is
realised on a backend, never per evaluation (``lowering.md`` §0). This module
is the torch half of that translation for the two things a declaration is made
of: the prior family (§3) and the unconstraining bijection (§4).

Every row goes through ``ampere.core.lowering``'s registry
-----------------------------------------------------------
W2.6 landed ``register_lowering``/``lookup_lowering`` and the analogous
bijection slot, and ``lowering.md`` §12.8 keys them on the backend's one name
— ``"torch"`` here, the same string this package's models declare as their
``BACKEND``. The tables below are therefore *registered*, at import, with
``builtin=True``, rather than being a private dict this module consults
directly. Two things follow, and both are the point:

* a user who registers a lowering for a family torch does not have (with
  ``override=True`` where a built-in row already exists) is served by exactly
  the same lookup this module uses, so their row is not a second-class path;
* ``provenance_entries`` can report which non-built-in rows a run consulted,
  because every consultation goes through one place —
  :func:`consulted_resolutions` is what a driver stamps.

Refusals are typed and loud, never silent
------------------------------------------
``lowering.md`` §3.4 rule 3 forbids a "close enough" substitution outright: a
family with no exact construction from torch's primitives raises
:class:`~ampere.core.exceptions.LoweringError` naming the family, the
parameter and the backend. ``truncnorm`` is the row that exercises this — it is
absent from ``torch.distributions`` *and* from pyro — and it is §11 row 5 of
the conformance obligations.

Where a non-default ``loc``/``scale`` can be composed away exactly, it is
(§3.3's rule 1, "compose it away", preferred over raising): scipy gives every
continuous family a shift-and-stretch that torch's equivalents do not have, so
``halfnorm(loc=3, scale=2)`` becomes
``TransformedDistribution(HalfNormal(2), AffineTransform(loc=3, scale=1))``
rather than a ``HalfNormal(2)`` whose support has silently moved by three
units.

The trap this module exists to contain: torch reports a *codomain*, not a support
--------------------------------------------------------------------------------
``lowering.md`` §4 recommends asking the library for
``biject_to(lowered_distribution.support)`` rather than rebuilding the
transform by hand, so that the library stays the authority on its own
constraint registry. That advice is right and this module follows it — but it
rests on the lowered distribution reporting its true support, and for
``torch.distributions.TransformedDistribution`` **it does not**:
``TransformedDistribution.support`` is the *codomain of the last transform*,
not the base's support pushed forward. Measured against torch 2.13:

* ``TransformedDistribution(Uniform(log a, log b), ExpTransform())`` — the §3.2
  row for ``loguniform`` — reports ``GreaterThan(0.0)``, not
  ``Interval(a, b)``. ``biject_to`` would then hand back an exponential where
  ampere infers a logit, and the sampler would explore values above ``b``
  where the prior is ``-inf``;
* ``TransformedDistribution(HalfNormal(s), AffineTransform(loc=l, scale=1))`` —
  §3.3's shifted half-normal — reports ``Real()``, because an affine map's
  codomain is the whole line. ``biject_to`` would then hand back the identity,
  and the sampler would spend half its proposals below ``l``.

Neither is a torch bug: the codomain is what a `TransformedDistribution` can
know cheaply. It does mean a lowering that trusts it is wrong in a way that
still runs. So every composed distribution this module builds declares its
support explicitly (:class:`_SupportedTransformed`), taken from the scipy
family it is a translation of, and ``biject_to`` is then asked about *that*.
The agreement between the two routes is pinned by a test rather than assumed
(``tests/backends/test_torch_lowering.py``).

Examples
--------
>>> import scipy.stats as st
>>> from ampere.core import describe_prior
>>> from ampere.backends.torch.lowering import lower_prior
>>> prior = st.norm(1.0, 2.0)
>>> lowered = lower_prior(describe_prior(prior))
>>> native = float(lowered.log_prob(torch.tensor(1.0, dtype=torch.float64)))
>>> abs(native - float(prior.logpdf(1.0))) < 1e-12
True

A family torch cannot express is refused by name, never approximated:

>>> lower_prior(describe_prior(st.truncnorm(-1, 2, 5, 3)), parameter="temperature")
Traceback (most recent call last):
    ...
ampere.core.exceptions.LoweringError: prior family 'truncnorm' (parameter 'temperature')...
"""

from __future__ import annotations

import dataclasses
import math
import warnings
from typing import Any

import torch
import torch.distributions as dist
from torch.distributions import biject_to, constraints, transforms

from ampere.core.exceptions import LoweringError, LoweringFallbackWarning
from ampere.core.lowering import (
    LoweringResolution,
    lookup_bijection_lowering,
    lookup_lowering,
    register_bijection_lowering,
    register_lowering,
)
from ampere.core.parameter import (
    Bijection,
    HierarchicalPrior,
    Identity,
    Log,
    Logit,
    PriorSpec,
)

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor

__all__ = [
    "IcdfFallbackWarning",
    "LoweredPrior",
    "LoweringFallbackWarning",
    "consulted_resolutions",
    "lower_bijection",
    "lower_hierarchical",
    "lower_prior",
    "warn_icdf_fallback",
]


#: The shared warning, under the name W2.4 gave it.
#:
#: **W2.13 replaced this backend-local class** (fold-in 8, ruled 2026-09-07):
#: it lives in :mod:`ampere.core.exceptions` as
#: :class:`~ampere.core.exceptions.LoweringFallbackWarning`, so ``pytest.warns``
#: on one type catches whichever backend raised it — which is what both slice-1
#: docstrings asked for. The old name is kept as an alias because it is the one
#: a torch user has seen, and because the two backends' warnings were never
#: different things.
IcdfFallbackWarning = LoweringFallbackWarning


# ---------------------------------------------------------------------------
# The lowered object
# ---------------------------------------------------------------------------


class _SupportedTransformed(dist.TransformedDistribution):
    """A ``TransformedDistribution`` that reports the support it actually has.

    See the module docstring: torch's own ``support`` is the last transform's
    codomain, which for the two composed rows ``lowering.md`` §3.2/§3.3 need
    (``loguniform``, and any shifted family) is strictly larger than the
    distribution's real support. Every lowering that composes therefore states
    the support it is translating from, and this class is where that statement
    lives.

    Overriding ``support`` is enough to fix both consumers that matter:
    ``biject_to(lowered.support)`` (``lowering.md`` §4's recommended route to
    the bijection) and the out-of-support masking in
    :meth:`LoweredPrior.log_prob`.
    """

    def __init__(
        self,
        base_distribution: dist.Distribution,
        transform_list: list[transforms.Transform],
        support: constraints.Constraint,
    ) -> None:
        super().__init__(base_distribution, transform_list, validate_args=False)
        self._ampere_support = support

    @property
    def support(self) -> constraints.Constraint:  # type: ignore[override]
        return self._ampere_support


@dataclasses.dataclass(frozen=True)
class LoweredPrior:
    """One prior family, lowered: the torch object plus what a caller needs about it.

    ``distribution`` is the ``torch.distributions`` object §3.2's table names.
    The other three fields exist because a caller that had to rediscover them
    would rediscover them differently in each backend:

    ``interior``
        A point strictly inside the support, taken from the scipy family's own
        median. :meth:`log_prob` substitutes it for out-of-support inputs
        before evaluating, so that a ``nan`` from (say) ``Gamma.log_prob`` at a
        negative argument never reaches the ``torch.where`` that masks it —
        ``where`` propagates ``nan`` gradients from the branch it discards,
        which is the standard way a masked density silently poisons an HMC
        trajectory.
    ``has_icdf``
        Whether the lowered object implements ``icdf``, established by calling
        it once at lowering time rather than by consulting a table that can go
        stale. ``lowering.md`` §3.6's gap — ``gamma``, ``beta`` and ``poisson``
        have none on torch — is what the fallback contract is about.
    ``spec``
        The declaration this was lowered from, kept so the fallback can rebuild
        the scipy prior without the caller carrying it separately.
    """

    spec: PriorSpec
    distribution: dist.Distribution
    interior: torch.Tensor
    has_icdf: bool

    @property
    def support(self) -> constraints.Constraint:
        """The constraint ``biject_to`` should be asked about (``lowering.md`` §4)."""
        return self.distribution.support  # type: ignore[return-value]

    def log_prob(self, value: torch.Tensor) -> torch.Tensor:
        """``log p(x)``, ``-inf`` outside the support — never ``nan``.

        ``ampere.core``'s reference density (``log_density``) returns ``-inf``
        outside the support, so a lowering that raised (torch's
        ``validate_args=True`` behaviour) or returned ``nan`` (its
        ``validate_args=False`` behaviour for several families) would not be
        the same function. The substitution of :attr:`interior` before
        evaluating is what keeps the discarded branch finite; see the class
        docstring.
        """
        inside = self.distribution.support.check(value)
        safe = torch.where(inside, value, self.interior.to(value.dtype))
        density = self.distribution.log_prob(safe)
        return torch.where(inside, density, torch.full_like(density, -math.inf))

    def icdf(self, quantile: torch.Tensor) -> torch.Tensor:
        """The native inverse CDF. Only call it when :attr:`has_icdf`."""
        return self.distribution.icdf(quantile)


# ---------------------------------------------------------------------------
# §3.2's table, family by family
# ---------------------------------------------------------------------------


def _kwd(spec: PriorSpec, name: str, default: float) -> float:
    """One canonical keyword of a ``PriorSpec``, with scipy's own default.

    ``describe_prior`` omits an argument left at its scipy default, so
    ``describe_prior(st.beta(2, 3)).kwds`` has no ``loc`` or ``scale`` at all
    (``lowering.md`` §3.1: the emitted spec is canonical and keyword-only, but
    it is not exhaustive).
    """
    return float(spec.kwds.get(name, default))


def _shifted(
    base: dist.Distribution,
    loc: float,
    lower: float,
    *,
    dtype: torch.dtype,
    device: torch.device,
    scale: float = 1.0,
) -> dist.Distribution:
    """``lowering.md`` §3.3 rule 1: compose an affine ``loc``/``scale`` away, exactly.

    Preferred over raising for every family whose scipy ``loc`` is a pure
    translation of the support (``expon``, ``halfnorm``, ``gamma``,
    ``lognorm``), because it keeps the prior a first-class distribution object
    on the backend and is exact rather than an approximation. The support is
    declared explicitly — see :class:`_SupportedTransformed`.
    """
    return _SupportedTransformed(
        base,
        [
            transforms.AffineTransform(
                loc=as_tensor(loc, dtype=dtype, device=device),
                scale=as_tensor(scale, dtype=dtype, device=device),
            )
        ],
        constraints.greater_than(as_tensor(lower, dtype=dtype, device=device)),
    )


def _build_norm(spec: PriorSpec, *, dtype: torch.dtype, device: torch.device) -> dist.Distribution:
    return dist.Normal(
        as_tensor(_kwd(spec, "loc", 0.0), dtype=dtype, device=device),
        as_tensor(_kwd(spec, "scale", 1.0), dtype=dtype, device=device),
        validate_args=False,
    )


def _build_uniform(
    spec: PriorSpec, *, dtype: torch.dtype, device: torch.device
) -> dist.Distribution:
    """scipy's second argument is a **width**, not an upper bound (§3.3)."""
    low = _kwd(spec, "loc", 0.0)
    width = _kwd(spec, "scale", 1.0)
    return dist.Uniform(
        as_tensor(low, dtype=dtype, device=device),
        as_tensor(low + width, dtype=dtype, device=device),
        validate_args=False,
    )


def _build_halfnorm(
    spec: PriorSpec, *, dtype: torch.dtype, device: torch.device
) -> dist.Distribution:
    """``HalfNormal(scale)``, shifted by ``loc`` when there is one (§3.3)."""
    loc = _kwd(spec, "loc", 0.0)
    scale = _kwd(spec, "scale", 1.0)
    base = dist.HalfNormal(as_tensor(scale, dtype=dtype, device=device), validate_args=False)
    if loc == 0.0:
        return base
    return _shifted(base, loc, loc, dtype=dtype, device=device)


def _build_loguniform(
    spec: PriorSpec, *, dtype: torch.dtype, device: torch.device
) -> dist.Distribution:
    """No torch ``LogUniform``: §3.4 rule 1's exact construction from primitives.

    ``exp`` of a uniform on ``[log a, log b]`` is a log-uniform on ``[a, b]``
    exactly — the same construction torch's own ``LogNormal`` uses over
    ``Normal``. The support has to be declared, because the composed object
    would otherwise report ``GreaterThan(0)``; see the module docstring.
    """
    low = _kwd(spec, "a", 1.0)
    high = _kwd(spec, "b", 2.0)
    if not (low > 0.0 and high > low):
        raise LoweringError(
            spec.family,
            backend=BACKEND,
            detail=f"loguniform needs 0 < a < b, got a={low!r}, b={high!r}",
        )
    base = dist.Uniform(
        as_tensor(math.log(low), dtype=dtype, device=device),
        as_tensor(math.log(high), dtype=dtype, device=device),
        validate_args=False,
    )
    return _SupportedTransformed(
        base,
        [transforms.ExpTransform()],
        constraints.interval(
            as_tensor(low, dtype=dtype, device=device),
            as_tensor(high, dtype=dtype, device=device),
        ),
    )


def _build_poisson(
    spec: PriorSpec, *, dtype: torch.dtype, device: torch.device
) -> dist.Distribution:
    """Lowers cleanly as a *distribution*; never onto a gradient path (§3.5).

    ``default_bijection_for`` refuses a discrete family with a typed
    ``CapabilityError``, so a poisson parameter never reaches
    :func:`lower_bijection` — the refusal lives upstream, at the source, and
    this row does not need to repeat it.
    """
    return dist.Poisson(
        as_tensor(_kwd(spec, "mu", 1.0), dtype=dtype, device=device), validate_args=False
    )


def _build_lognorm(
    spec: PriorSpec, *, dtype: torch.dtype, device: torch.device
) -> dist.Distribution:
    """``loc_torch = log(scale_scipy)`` — passing ``scale`` through is wrong by an exp.

    scipy parametrises by the shape ``s`` (the log-scale standard deviation)
    and ``scale = exp(mu)``; torch's ``LogNormal`` takes ``mu`` directly
    (§3.3).
    """
    shape = _kwd(spec, "s", 1.0)
    loc = _kwd(spec, "loc", 0.0)
    scale = _kwd(spec, "scale", 1.0)
    if scale <= 0.0:
        raise LoweringError(
            spec.family, backend=BACKEND, detail=f"lognorm needs scale > 0, got {scale!r}"
        )
    base = dist.LogNormal(
        as_tensor(math.log(scale), dtype=dtype, device=device),
        as_tensor(shape, dtype=dtype, device=device),
        validate_args=False,
    )
    if loc == 0.0:
        return base
    return _shifted(base, loc, loc, dtype=dtype, device=device)


def _build_expon(spec: PriorSpec, *, dtype: torch.dtype, device: torch.device) -> dist.Distribution:
    """``rate = 1 / scale``, verified against ``torch.distributions.Exponential``."""
    loc = _kwd(spec, "loc", 0.0)
    scale = _kwd(spec, "scale", 1.0)
    if scale <= 0.0:
        raise LoweringError(
            spec.family, backend=BACKEND, detail=f"expon needs scale > 0, got {scale!r}"
        )
    base = dist.Exponential(as_tensor(1.0 / scale, dtype=dtype, device=device), validate_args=False)
    if loc == 0.0:
        return base
    return _shifted(base, loc, loc, dtype=dtype, device=device)


def _build_gamma(spec: PriorSpec, *, dtype: torch.dtype, device: torch.device) -> dist.Distribution:
    """``Gamma(concentration=a, rate=1/scale)``, with ``loc`` composed away."""
    shape = _kwd(spec, "a", 1.0)
    loc = _kwd(spec, "loc", 0.0)
    scale = _kwd(spec, "scale", 1.0)
    if scale <= 0.0:
        raise LoweringError(
            spec.family, backend=BACKEND, detail=f"gamma needs scale > 0, got {scale!r}"
        )
    base = dist.Gamma(
        as_tensor(shape, dtype=dtype, device=device),
        as_tensor(1.0 / scale, dtype=dtype, device=device),
        validate_args=False,
    )
    if loc == 0.0:
        return base
    return _shifted(base, loc, loc, dtype=dtype, device=device)


def _build_beta(spec: PriorSpec, *, dtype: torch.dtype, device: torch.device) -> dist.Distribution:
    """``Beta(concentration1=a, concentration0=b)`` on ``[loc, loc+scale]``.

    §3.2's row says ``loc`` must be 0 and ``scale`` 1, which is the constraint
    on the *naive* mapping; §3.3's rule 1 ranks composing an exact affine map
    above raising, and for ``beta`` the map ``x = loc + scale * u`` is exactly
    that. So a shifted, stretched beta lowers rather than being refused, and
    its support is declared as the interval it really is.
    """
    a = _kwd(spec, "a", 1.0)
    b = _kwd(spec, "b", 1.0)
    loc = _kwd(spec, "loc", 0.0)
    scale = _kwd(spec, "scale", 1.0)
    if scale <= 0.0:
        raise LoweringError(
            spec.family, backend=BACKEND, detail=f"beta needs scale > 0, got {scale!r}"
        )
    base = dist.Beta(
        as_tensor(a, dtype=dtype, device=device),
        as_tensor(b, dtype=dtype, device=device),
        validate_args=False,
    )
    if loc == 0.0 and scale == 1.0:
        return base
    return _SupportedTransformed(
        base,
        [
            transforms.AffineTransform(
                loc=as_tensor(loc, dtype=dtype, device=device),
                scale=as_tensor(scale, dtype=dtype, device=device),
            )
        ],
        constraints.interval(
            as_tensor(loc, dtype=dtype, device=device),
            as_tensor(loc + scale, dtype=dtype, device=device),
        ),
    )


#: Families ``lowering.md`` §3.2 records as having no exact torch construction,
#: with the reason the refusal message carries. ``truncnorm`` is the §11 row 5
#: case: absent from ``torch.distributions`` *and* from pyro, and not
#: expressible by truncating a torch primitive, because torch has no truncation
#: wrapper at all.
_NO_EXACT_CONSTRUCTION: dict[str, str] = {
    "truncnorm": (
        "torch.distributions has no TruncatedNormal (nor any truncation wrapper to build one "
        "from), and neither does pyro — so there is no exact construction from torch's own "
        "primitives, and lowering.md §3.4 rule 3 forbids approximating it with an untruncated "
        "Normal"
    ),
}


# The §3.2 table, as registry rows. One thin wrapper per family so that the
# registered constructor has a real ``__module__``/``__qualname__`` for
# provenance (``LoweringResolution.constructor``), and so that dtype and device
# are passed explicitly at every construction site rather than inherited from
# ``torch.get_default_dtype()`` (``lowering.md`` §10.1).


def _norm(spec: PriorSpec) -> dist.Distribution:
    return _build_norm(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _uniform(spec: PriorSpec) -> dist.Distribution:
    return _build_uniform(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _halfnorm(spec: PriorSpec) -> dist.Distribution:
    return _build_halfnorm(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _loguniform(spec: PriorSpec) -> dist.Distribution:
    return _build_loguniform(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _poisson(spec: PriorSpec) -> dist.Distribution:
    return _build_poisson(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _lognorm(spec: PriorSpec) -> dist.Distribution:
    return _build_lognorm(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _expon(spec: PriorSpec) -> dist.Distribution:
    return _build_expon(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _gamma(spec: PriorSpec) -> dist.Distribution:
    return _build_gamma(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


def _beta(spec: PriorSpec) -> dist.Distribution:
    return _build_beta(spec, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


#: The dtype/device-aware builders behind the registered wrappers. A built-in
#: row is dispatched through here rather than through the zero-argument wrapper
#: so that ``lower_prior``'s ``dtype=``/``device=`` reach the construction —
#: the registry's constructor signature is ``PriorSpec -> object`` and has
#: nowhere to carry them. A user-registered row is called through the registry
#: unchanged, and is responsible for its own precision.
_BUILTIN_BUILDERS = {
    "norm": _build_norm,
    "uniform": _build_uniform,
    "halfnorm": _build_halfnorm,
    "loguniform": _build_loguniform,
    "poisson": _build_poisson,
    "lognorm": _build_lognorm,
    "expon": _build_expon,
    "gamma": _build_gamma,
    "beta": _build_beta,
}


_BUILTIN_PRIORS = {
    "norm": _norm,
    "uniform": _uniform,
    "halfnorm": _halfnorm,
    "loguniform": _loguniform,
    "poisson": _poisson,
    "lognorm": _lognorm,
    "expon": _expon,
    "gamma": _gamma,
    "beta": _beta,
}


# ---------------------------------------------------------------------------
# §4's table: bijections
# ---------------------------------------------------------------------------


def _identity_transform(_bijection: Bijection) -> transforms.Transform:
    """``x = y``. What ``biject_to(constraints.real)`` returns, by construction."""
    return biject_to(constraints.real)


def _log_transform(bijection: Bijection) -> transforms.Transform:
    """``x = lower + exp(y)``.

    ``ComposeTransform([f, g])`` applies ``f`` first in both target libraries,
    so the exponential comes first and the shift second (§4's compose-order
    note). The other order gives ``exp(lower + y)`` — a different, and still
    perfectly finite, prior.

    §4's note that "numpyro's composed forms carry an explicit ``domain=`` on
    the inner ``AffineTransform``" is about numpyro and **only** numpyro:
    ``torch.distributions.transforms.AffineTransform`` has no such parameter,
    and passing one is a ``TypeError``. The distinction is easy to lose because
    §4's table puts the two libraries in adjacent columns.
    """
    lower = float(getattr(bijection, "lower", 0.0))
    exponential = transforms.ExpTransform()
    if lower == 0.0:
        return exponential
    return transforms.ComposeTransform(
        [
            exponential,
            transforms.AffineTransform(loc=as_tensor(lower), scale=as_tensor(1.0)),
        ]
    )


def _logit_transform(bijection: Bijection) -> transforms.Transform:
    """``x = lower + (upper - lower) * sigmoid(y)``."""
    lower = float(getattr(bijection, "lower", 0.0))
    upper = float(getattr(bijection, "upper", 1.0))
    sigmoid = transforms.SigmoidTransform()
    if lower == 0.0 and upper == 1.0:
        return sigmoid
    return transforms.ComposeTransform(
        [
            sigmoid,
            transforms.AffineTransform(loc=as_tensor(lower), scale=as_tensor(upper - lower)),
        ]
    )


_BUILTIN_BIJECTIONS: dict[type, Any] = {
    Identity: _identity_transform,
    Log: _log_transform,
    Logit: _logit_transform,
}


# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------


def _register_builtins() -> None:
    """Put §3.2's and §4's tables in ``ampere.core.lowering``, keyed on ``"torch"``.

    Idempotent: importing this module twice (or importing it after a user has
    deliberately overridden a row) must not raise, and must not clobber a
    user's ``override=True`` registration either. So a row is registered only
    if the key is empty.
    """
    for family, constructor in _BUILTIN_PRIORS.items():
        try:
            lookup_lowering(family, BACKEND)
        except LoweringError:
            register_lowering(family, BACKEND, constructor, builtin=True)
    for bijection_class, transform_constructor in _BUILTIN_BIJECTIONS.items():
        try:
            lookup_bijection_lowering(bijection_class, BACKEND)
        except LoweringError:
            register_bijection_lowering(
                bijection_class, BACKEND, transform_constructor, builtin=True
            )


_register_builtins()


# ---------------------------------------------------------------------------
# The lowering entry points
# ---------------------------------------------------------------------------


def _interior_point(
    support: constraints.Constraint, *, dtype: torch.dtype, device: torch.device
) -> torch.Tensor:
    """A point strictly inside *support*, read off the constraint's own bounds.

    Only ever used as the substitute value :meth:`LoweredPrior.log_prob` feeds
    to the density in place of an out-of-support input, so all that is asked of
    it is that the density be finite there. Derived from the constraint rather
    than from the scipy family's median because a hierarchical prior is lowered
    afresh on every evaluation (``lowering.md`` §8) and a scipy round trip per
    evaluation would put the reference path back in the hot loop, which is
    exactly what a native backend exists to avoid.
    """
    lower = getattr(support, "lower_bound", None)
    upper = getattr(support, "upper_bound", None)
    if lower is not None and upper is not None:
        point = 0.5 * (float(lower) + float(upper))
    elif lower is not None:
        point = float(lower) + 1.0
    elif upper is not None:
        point = float(upper) - 1.0
    else:
        point = 0.0
    return as_tensor(point, dtype=dtype, device=device)


def _has_native_icdf(distribution: dist.Distribution, interior: torch.Tensor) -> bool:
    """Whether the lowered object implements ``icdf`` — asked, not looked up.

    ``torch.distributions.Distribution.icdf`` raises ``NotImplementedError``
    unless a family overrides it, so the honest check is to call it once, at
    lowering time. ``lowering.md`` §3.6's table (``gamma``, ``beta`` and
    ``poisson`` have none) is then documentation rather than the mechanism, and
    a torch release that fills one of those gaps is picked up without an
    ampere change.
    """
    try:
        distribution.icdf(torch.full_like(interior, 0.5))
    except NotImplementedError:
        return False
    return True


def lower_prior(
    spec: PriorSpec,
    *,
    parameter: str | None = None,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> LoweredPrior:
    """``lowering.md`` §3.2's table, for one declared prior.

    The registry is the authority on which constructor applies (see the module
    docstring), so a user-registered row for a family torch lacks — or an
    ``override=True`` replacement of a built-in one — is served here without
    this module knowing about it.

    Parameters
    ----------
    spec
        The canonical, keyword-only description ``describe_prior`` emits.
    parameter
        The merged parameter name, named in the error if lowering fails. A
        refusal that does not say *which* parameter caused it makes the user
        bisect their own declaration.
    dtype, device
        Threaded explicitly into every tensor the construction creates
        (``lowering.md`` §10.1). ``torch.set_default_dtype`` is never called.

    Raises
    ------
    LoweringError
        If no row exists for this family on this backend — §3.4 rule 2, never
        rule 3's silent approximation.
    """
    try:
        resolution = lookup_lowering(spec.family, BACKEND)
    except LoweringError:
        detail = _NO_EXACT_CONSTRUCTION.get(spec.family)
        raise LoweringError(
            spec.family,
            backend=BACKEND,
            parameter=parameter,
            detail=detail,
        ) from None
    builder = _BUILTIN_PRIORS.get(spec.family)
    if resolution.builtin and builder is not None and builder is resolution.constructor:
        distribution = _BUILTIN_BUILDERS[spec.family](spec, dtype=dtype, device=device)
    else:
        distribution = resolution.constructor(spec)
    interior = _interior_point(distribution.support, dtype=dtype, device=device)
    return LoweredPrior(
        spec=spec,
        distribution=distribution,
        interior=interior,
        has_icdf=_has_native_icdf(distribution, interior),
    )


def _hierarchical_norm(arguments: dict[str, torch.Tensor]) -> dist.Distribution:
    return dist.Normal(arguments["loc"], arguments["scale"], validate_args=False)


def _hierarchical_uniform(arguments: dict[str, torch.Tensor]) -> dist.Distribution:
    low = arguments["loc"]
    return dist.Uniform(low, low + arguments["scale"], validate_args=False)


def _hierarchical_halfnorm(arguments: dict[str, torch.Tensor]) -> dist.Distribution:
    base = dist.HalfNormal(arguments["scale"], validate_args=False)
    loc = arguments["loc"]
    if bool(torch.all(loc == 0.0)):
        return base
    return _SupportedTransformed(
        base,
        [transforms.AffineTransform(loc=loc, scale=torch.ones_like(loc))],
        constraints.greater_than(loc),
    )


def _hierarchical_expon(arguments: dict[str, torch.Tensor]) -> dist.Distribution:
    base = dist.Exponential(1.0 / arguments["scale"], validate_args=False)
    loc = arguments["loc"]
    if bool(torch.all(loc == 0.0)):
        return base
    return _SupportedTransformed(
        base,
        [transforms.AffineTransform(loc=loc, scale=torch.ones_like(loc))],
        constraints.greater_than(loc),
    )


#: Hierarchical families this backend can build from *tensor* arguments.
#:
#: A ``HierarchicalPrior`` "cannot be lowered once and cached" (``lowering.md``
#: §8): its arguments are other parameters' values, which change every
#: evaluation, so the distribution object is constructed per evaluation from
#: the current tensors. That construction has to keep the autograd graph — the
#: point of a hierarchical prior on a differentiable backend is that the
#: hyperparameters get a gradient — which rules out the ``float()`` coercion
#: the §3.2 table rows use, and with it every family whose lowering needs a
#: Python-level branch on an argument's value.
#:
#: What is left is close to what ``ampere.core`` can infer a bijection for, and
#: that is not a coincidence: a hierarchical family taking shape arguments is
#: already refused upstream by ``_default_bijection_for_hierarchical``, because
#: its support is not determined by ``loc``/``scale`` alone.
_HIERARCHICAL_BUILDERS: dict[str, Any] = {
    "norm": _hierarchical_norm,
    "uniform": _hierarchical_uniform,
    "halfnorm": _hierarchical_halfnorm,
    "expon": _hierarchical_expon,
}

#: scipy's defaults for the two location-scale arguments every family above
#: takes, used when the declaration supplies neither a constant nor a reference.
_HIERARCHICAL_DEFAULTS = {"loc": 0.0, "scale": 1.0}

#: scipy's positional order for a frozen location-scale family, so a
#: ``HierarchicalPrior`` declared with ``args`` rather than ``kwds`` lowers the
#: same way ``HierarchicalPrior.bind`` would freeze it.
_HIERARCHICAL_POSITIONAL = ("loc", "scale")


def lower_hierarchical(
    prior: HierarchicalPrior,
    resolved: dict[str, torch.Tensor],
    *,
    parameter: str | None = None,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> LoweredPrior:
    """``lowering.md`` §8: a hierarchical prior, built afresh from current tensors.

    "A ``HierarchicalPrior`` cannot be lowered once and cached." Its arguments
    are other parameters' *values*, so hoisting the construction out of the
    evaluation loop as an optimisation would freeze the prior at its initial
    hyperparameters and break into a plausible-looking wrong answer — the
    failure §8 names explicitly, and the one ``lowering.md`` §11 row 7 checks
    for.

    Parameters
    ----------
    prior
        The declaration: its family, its constant arguments, and its
        references by name.
    resolved
        Current values of every parameter this prior references, as tensors,
        so the constructed distribution keeps the autograd graph back to them.
    parameter
        The merged name, for the error message.

    Raises
    ------
    LoweringError
        For a family outside :data:`_HIERARCHICAL_BUILDERS` — see its comment
        for why that set is what it is — or for a reference that has not been
        resolved yet, which means the evaluation order is wrong.
    """
    builder = _HIERARCHICAL_BUILDERS.get(prior.family)
    if builder is None:
        raise LoweringError(
            prior.family,
            backend=BACKEND,
            parameter=parameter,
            detail=(
                f"hierarchical priors lower on the {BACKEND!r} backend only for the "
                f"location-scale families {sorted(_HIERARCHICAL_BUILDERS)}, whose arguments can "
                f"be passed as tensors so the hyperparameters keep a gradient; family "
                f"{prior.family!r} is not one of them"
            ),
        )
    arguments: dict[str, torch.Tensor] = {
        keyword: as_tensor(value, dtype=dtype, device=device)
        for keyword, value in _HIERARCHICAL_DEFAULTS.items()
    }
    for index, positional in enumerate(prior.args):
        if index < len(_HIERARCHICAL_POSITIONAL):
            arguments[_HIERARCHICAL_POSITIONAL[index]] = as_tensor(
                positional, dtype=dtype, device=device
            )
    for keyword, constant in prior.kwds.items():
        arguments[keyword] = as_tensor(constant, dtype=dtype, device=device)
    for keyword, reference in prior.hyperparameters.items():
        try:
            arguments[keyword] = resolved[reference]
        except KeyError:
            raise LoweringError(
                prior.family,
                backend=BACKEND,
                parameter=parameter,
                detail=(
                    f"hierarchical prior references parameter {reference!r}, which is not "
                    f"resolved yet: a hyperparameter must be lowered before anything "
                    f"referencing it (lowering.md §8's topological emission order)"
                ),
            ) from None
    distribution = builder(arguments)
    interior = _interior_point(distribution.support, dtype=dtype, device=device)
    return LoweredPrior(
        spec=PriorSpec.of(prior.family),
        distribution=distribution,
        interior=interior,
        has_icdf=_has_native_icdf(distribution, interior),
    )


def lower_bijection(
    bijection: Bijection | None,
    lowered: LoweredPrior,
    *,
    parameter: str | None = None,
) -> transforms.Transform:
    """``lowering.md`` §4's table, for one parameter's unconstraining bijection.

    Two routes, and which one applies is fixed by §4 rather than by taste:

    * *bijection is None* — the parameter took the **inferred default**, read
      off its prior's support. §4's rule is then to ask the target library:
      ``biject_to(lowered.support)``, so the library stays the authority on its
      own constraint registry and the table serves as the cross-check. The two
      agree because ampere adopted numpyro's rule, and that agreement is
      pinned by a test rather than assumed.
    * *bijection declared* — the user asked for something specific, so
      ``biject_to`` would give the wrong answer. The declared class is looked
      up in the registry (``Identity``/``Log``/``Logit`` are registered as
      built-in rows; anything else needs
      ``register_bijection_lowering(..., "torch", ...)``) and built by hand
      from §4's table.

    Raises
    ------
    LoweringError
        For a custom :class:`~ampere.core.parameter.Bijection` with no
        registered torch row — §4's "user-supplied custom ``Bijection``:
        unsupported — raises at lowering". A custom bijection is a
        reference-path-only feature until someone registers it, exactly as a
        duck-typed prior is.
    """
    if bijection is None:
        return biject_to(lowered.support)
    try:
        resolution = lookup_bijection_lowering(bijection, BACKEND)
    except LoweringError:
        raise LoweringError(
            type(bijection).__name__,
            backend=BACKEND,
            parameter=parameter,
            detail=(
                f"{type(bijection).__name__} is a user-supplied Bijection with no torch "
                f"lowering: it satisfies ampere's protocol with numpy operations a torch "
                f"tensor will not accept. Register one with "
                f"register_bijection_lowering({type(bijection).__name__}, 'torch', ...)"
            ),
        ) from None
    return resolution.constructor(bijection)


def warn_icdf_fallback(families: set[str], *, strict: bool, where: str) -> None:
    """``lowering.md`` §3.6's contract, in one place: warn once, or raise under ``strict``.

    Called once per lowering — never per call — because the decision is made
    once, before any sampling, and repeating it per proposal would bury the
    run's output.

    Parameters
    ----------
    families
        The families with no native ``icdf``, so the message names what to
        change rather than saying only that something was missing.
    strict
        The problem's own ``strict`` flag (``FittingProblem(strict=...)``), not
        a new per-lowering knob: one flag, one meaning.
    where
        What is taking the fallback, for the message ("prior_transform").
    """
    if not families:
        return
    listed = ", ".join(sorted(families))
    message = (
        f"{where} for prior families [{listed}] has no native icdf on the {BACKEND!r} backend, "
        f"so it is evaluated on the reference (scipy) path. The value is identical — "
        f"prior_transform has one mathematical definition — but this run's prior transform is "
        f"not computed by {BACKEND}. Pass FittingProblem(strict=True) to refuse the fallback "
        f"instead, or choose priors from the families torch implements icdf for."
    )
    if strict:
        raise LoweringError(
            listed,
            backend=BACKEND,
            detail=(
                f"strict=True refuses the reference fallback: {where} needs a native icdf for "
                f"prior families [{listed}], and torch implements none for them. Either choose "
                f"a family torch can invert (norm, uniform, halfnorm, lognorm, loguniform, "
                f"expon), run on a backend that has it, or drop strict=True to accept the "
                f"reference (scipy) prior transform with a warning."
            ),
        )
    warnings.warn(message, IcdfFallbackWarning, stacklevel=3)


def consulted_resolutions(families: set[str]) -> tuple[LoweringResolution, ...]:
    """Every registry row a lowering of *families* consulted, for provenance.

    Pass through ``ampere.core.lowering.provenance_entries`` — which keeps only
    the non-built-in rows — and on to ``provenance_attrs(extra=...)``. That is
    ``lowering.md`` §12.8's "every registered row is stamped user-registered in
    provenance", discharged by the backend that did the consulting rather than
    by every driver separately.
    """
    found: list[LoweringResolution] = []
    for family in sorted(families):
        try:
            found.append(lookup_lowering(family, BACKEND))
        except LoweringError:  # pragma: no cover - a lowered family always has a row
            continue
    return tuple(found)
