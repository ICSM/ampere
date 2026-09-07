"""The likelihood families, lowered into jax for the realised density.

What this module is, and what it is not
---------------------------------------
``ampere.core``'s :class:`~ampere.core.LikelihoodFamily` classes are the
declaration *and* the reference implementation: they say what ``student_t``
means and they compute it, in numpy and ``scipy.stats``. Those computations
are not traceable — ``scipy.stats`` is numpy — so a differentiable backend has
to compute the same quantity in its own array library.

This module is that transcription and nothing else. Every function here is
``ampere.core.likelihood``'s own closed form written in ``jax.numpy``, with
the core's line beside it, so a reader can check the two against each other in
one place instead of hunting through a composition. The numpy path stays the
oracle: ``tests/conformance``'s ``TestTheRealisation`` compares the two at 24
points per declared shape, and the agreement is the reason this file may exist
at all.

Which families are here
-----------------------
Exactly the ones ``ampere.core`` **implements** — ``gaussian``,
``student_t``, ``cauchy``, ``complex_gaussian`` and ``poisson``. ``rice`` and
``von_mises`` are declared-but-unimplemented slots on the reference path
(their ``log_prob`` raises), so there is nothing here to agree with: they are
refused by name at lowering time, which is the same answer the numpy path
gives, arrived at earlier.

Censoring
---------
``likelihoods.md``'s Tobit construction, which ``ampere.core``'s
``_location_scale_log_prob`` writes once for every location-scale family: a
detection contributes ``log f(z) - log sigma``, an upper limit ``log F(z)``,
a lower limit ``log (1 - F(z))``. It is written once here too, in
:func:`_tobit`, over three callables — and the only reason it needed thought
is that jax has no ``logcdf`` for Student's t. It has ``betainc``, which is
what ``scipy``'s own t-distribution CDF is built from, so the closed form is
written out (:func:`_student_t_logcdf`) and agrees with ``scipy.stats.t`` to
1.1e-15 over the range the conformance rows exercise. Nothing here
approximates a censored likelihood; a family whose CDF this module could not
compute would be refused, not estimated.

The masked-sample rule needs no code at all: masking is **excision**
(``likelihoods.md`` §8), applied by
:mod:`ampere.backends.jax.problem` before anything here is called, so a family
never sees a mask and cannot forget to honour one.

Trace purity
------------
No function here raises. Every refusal this module makes is made by
:func:`lower_family` at **lowering** time, on a concrete declaration; inside
the traced density a failure is ``-inf`` through ``jnp.where``
(``likelihoods.md`` §17 Q1). The Poisson positivity guard is the worked
example: the reference path raises ``LikelihoodError`` when the model predicts
a non-positive expected count, this path returns ``-inf``, and both are §4.5's
answer — one recorded with a reason, one not, exactly as ``inference.md`` §10a
sub-decision 2 says a realised density behaves.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import jax
import jax.numpy as jnp
import jax.scipy.stats as jst
from jax.scipy.special import betainc, gammaln

from ampere.core import LikelihoodFamily, LimitKind
from ampere.core.exceptions import LoweringError

from ._config import BACKEND

__all__ = ["SUPPORTED_FAMILIES", "lower_family"]

_LOG_2PI = math.log(2.0 * math.pi)

#: One term per sample, so the Tobit sum can be written as a single ``where``
#: rather than three boolean-indexed gathers — which a trace cannot do, since
#: the number of detections is a value.
_Terms = Callable[[jax.Array], jax.Array]


def _tobit(
    logpdf: _Terms,
    logcdf: _Terms,
    logsf: _Terms,
    standardised: jax.Array,
    sigma: jax.Array,
    limits: jax.Array | None,
) -> jax.Array:
    """``ampere.core.likelihood._location_scale_log_prob``, in jax.

    The whole of the censored-likelihood mathematics for a location-scale
    family, and the reason ``likelihoods.md`` can say upper limits are
    "fairly straightforward" for uncorrelated noise: a detection contributes
    ``log f(z) - log sigma``, an upper limit ``log F(z)`` — the probability the
    truth lies below the recorded value — and a lower limit ``log (1 - F(z))``.

    The core selects with boolean indexing and sums three pieces. That cannot
    be traced: which samples are detections is a *value*, and the shapes would
    depend on it. So all three terms are computed for every sample and one
    ``jnp.where`` selects, which is the same arithmetic with a fixed shape.
    Both are O(N) and the difference is three cheap evaluations per sample —
    a price worth paying to keep the whole density one traceable expression.
    """
    detected = logpdf(standardised) - jnp.log(sigma)
    if limits is None:
        return jnp.sum(detected)
    codes = jnp.asarray(limits)
    terms = jnp.where(
        codes == int(LimitKind.UPPER_LIMIT),
        logcdf(standardised),
        jnp.where(codes == int(LimitKind.LOWER_LIMIT), logsf(standardised), detected),
    )
    return jnp.sum(terms)


def _student_t_tail(standardised: jax.Array, nu: jax.Array) -> jax.Array:
    """``I_x(nu/2, 1/2)`` with ``x = nu / (nu + z²)`` — twice the smaller tail.

    Student's t has no ``logcdf`` in ``jax.scipy.stats``, and this is the
    regularised incomplete beta function ``scipy``'s own ``t.cdf`` is built
    from, so writing it out is transcription rather than approximation. The
    identity is ``F(z) = I_x(nu/2, 1/2) / 2`` for ``z <= 0`` and ``1 - I_x / 2`` for
    ``z > 0``; measured against ``scipy.stats.t.logcdf`` it agrees to 1.1e-15.
    """
    z = jnp.asarray(standardised, dtype=jnp.float64)
    return betainc(0.5 * nu, 0.5, nu / (nu + z * z))


def _student_t_logcdf(standardised: jax.Array, nu: jax.Array) -> jax.Array:
    tail = _student_t_tail(standardised, nu)
    return jnp.where(standardised <= 0.0, jnp.log(0.5 * tail), jnp.log1p(-0.5 * tail))


def _student_t_logsf(standardised: jax.Array, nu: jax.Array) -> jax.Array:
    tail = _student_t_tail(standardised, nu)
    return jnp.where(standardised >= 0.0, jnp.log(0.5 * tail), jnp.log1p(-0.5 * tail))


def _own_values(family: LikelihoodFamily, values: Mapping[str, Any]) -> Mapping[str, Any]:
    """The family's own parameters out of the likelihood-wide resolved values.

    The same filtering ``ampere.core``'s families do
    (``{k: v for k, v in noise.values.items() if k in self.parameters}``), and
    for the same reason: ``Likelihood`` passes one flat mapping holding the
    family's parameters *and* the noise model's, and completing a declaration
    against names it does not own is an error.
    """
    return family.context({key: value for key, value in values.items() if key in family.parameters})


# ---------------------------------------------------------------------------
# The families
# ---------------------------------------------------------------------------


def _gaussian(
    predicted: jax.Array,
    observed: jax.Array,
    sigma: jax.Array | None,
    family: LikelihoodFamily,
    values: Mapping[str, Any],
    limits: jax.Array | None,
    latent: jax.Array | None,
) -> jax.Array:
    """``ampere.core.GaussianFamily.log_prob``, uncorrelated branch, in jax."""
    assert sigma is not None  # lower_family has checked the declaration
    standardised = (observed - predicted) / sigma
    if limits is None:
        # The core writes the uncensored case out rather than delegating,
        # because it is the hot loop of every ordinary fit; so does this.
        return jnp.sum(-0.5 * (standardised**2 + _LOG_2PI) - jnp.log(sigma))
    return _tobit(
        jst.norm.logpdf,
        jst.norm.logcdf,
        lambda z: jst.norm.logcdf(-z),
        standardised,
        sigma,
        limits,
    )


def _student_t(
    predicted: jax.Array,
    observed: jax.Array,
    sigma: jax.Array | None,
    family: LikelihoodFamily,
    values: Mapping[str, Any],
    limits: jax.Array | None,
    latent: jax.Array | None,
) -> jax.Array:
    """``ampere.core.StudentTFamily.log_prob`` in jax. ``nu`` keeps its gradient."""
    assert sigma is not None
    nu = jnp.asarray(_own_values(family, values)["nu"], dtype=jnp.float64)
    standardised = (observed - predicted) / sigma
    return _tobit(
        lambda z: jst.t.logpdf(z, nu),
        lambda z: _student_t_logcdf(z, nu),
        lambda z: _student_t_logsf(z, nu),
        standardised,
        sigma,
        limits,
    )


def _cauchy(
    predicted: jax.Array,
    observed: jax.Array,
    sigma: jax.Array | None,
    family: LikelihoodFamily,
    values: Mapping[str, Any],
    limits: jax.Array | None,
    latent: jax.Array | None,
) -> jax.Array:
    """``ampere.core.CauchyFamily.log_prob`` in jax.

    The Cauchy CDF is elementary — ``½ + arctan(z)/π`` — so unlike Student's t
    it needs no incomplete beta, and the survival function is its reflection.
    """
    assert sigma is not None
    standardised = (observed - predicted) / sigma
    return _tobit(
        jst.cauchy.logpdf,
        lambda z: jnp.log(0.5 + jnp.arctan(z) / jnp.pi),
        lambda z: jnp.log(0.5 - jnp.arctan(z) / jnp.pi),
        standardised,
        sigma,
        limits,
    )


def _complex_gaussian(
    predicted: jax.Array,
    observed: jax.Array,
    sigma: jax.Array | None,
    family: LikelihoodFamily,
    values: Mapping[str, Any],
    limits: jax.Array | None,
    latent: jax.Array | None,
) -> jax.Array:
    """``ampere.core.ComplexGaussianFamily.log_prob``, uncorrelated branch, in jax.

    The circular complex Gaussian: real and imaginary parts independent
    ``Normal(0, sigma**2)``, so ``log p = sum[-log(2 pi sigma**2) - |y - mu|**2 / (2 sigma**2)]``.
    The GP branch is refused on the reference path too
    (``GP_ANALYTIC_IMPLEMENTED = False``, Phase 4's), so this never sees one.
    """
    assert sigma is not None
    residual = jnp.abs(observed - predicted)
    variance = sigma**2
    return jnp.sum(-(residual**2) / (2.0 * variance) - _LOG_2PI - jnp.log(variance))


def _poisson(
    predicted: jax.Array,
    observed: jax.Array,
    sigma: jax.Array | None,
    family: LikelihoodFamily,
    values: Mapping[str, Any],
    limits: jax.Array | None,
    latent: jax.Array | None,
) -> jax.Array:
    """``ampere.core.PoissonFamily.log_prob`` in jax, latent branch included.

    ``log p = Σ [k log λ - λ - log k!]`` with ``log k!`` from ``gammaln``,
    which is what ``scipy.stats.poisson.logpmf`` computes (agreement measured
    at 1.8e-15). With a latent GP the rate is ``λ = rate · exp(f)``, which is
    the case ``DEVELOPMENT_PLAN.md`` §4.4 singles out and the reason this
    family exists on a differentiable backend at all: there is no covariance
    to add to a Poisson, so robustness needs the latent formulation, and
    ``Likelihood.check_engine`` refuses it to every gradient-free engine.

    The reference path *raises* for a non-positive expected count. Here it is
    ``-inf``: a Python ``if`` on a traced value is an error rather than a
    branch, and ``-inf`` is the same §4.5 answer without the recorded reason
    (``inference.md`` §10a sub-decision 2).
    """
    counts = jnp.asarray(observed, dtype=jnp.float64)
    rate = jnp.asarray(predicted, dtype=jnp.float64)
    if latent is not None:
        rate = rate * jnp.exp(latent)
    value = jnp.sum(counts * jnp.log(rate) - rate - gammaln(counts + 1.0))
    return jnp.where(jnp.all(rate > 0.0), value, -jnp.inf)


#: Family ``NAME`` -> its lowered evaluation. The keys are exactly the
#: families ``ampere.core`` implements; see this module's docstring for why
#: ``rice`` and ``von_mises`` are absent rather than stubbed.
_FAMILIES: dict[str, Any] = {
    "gaussian": _gaussian,
    "student_t": _student_t,
    "cauchy": _cauchy,
    "complex_gaussian": _complex_gaussian,
    "poisson": _poisson,
}

#: The family names this backend can compose into a differentiable density.
SUPPORTED_FAMILIES: frozenset[str] = frozenset(_FAMILIES)


def lower_family(family: LikelihoodFamily, label: str) -> Any:
    """The jax evaluation of *family*, or a :class:`LoweringError` naming it.

    Called once, when a problem is realised — never per evaluation, and never
    inside a trace (``lowering.md`` §0, ``likelihoods.md`` §17 Q1).

    Parameters
    ----------
    family
        The declared family, as composed into the dataset's ``Likelihood``.
    label
        The dataset's label, so the refusal says *which* dataset is the
        problem rather than only which family.
    """
    lowered = _FAMILIES.get(family.NAME)
    if lowered is not None:
        return lowered
    known = ", ".join(sorted(SUPPORTED_FAMILIES))
    if not family.IMPLEMENTED:
        raise LoweringError(
            family.NAME,
            backend=BACKEND,
            detail=(
                f"dataset {label!r} declares the {family.NAME!r} likelihood family, which "
                f"ampere.core declares but does not implement either — its own log_prob raises "
                f"(Phase 4, with the visibility modality). There is nothing for this backend to "
                f"agree with yet. Families this backend lowers: {known}."
            ),
        )
    raise LoweringError(
        family.NAME,
        backend=BACKEND,
        detail=(
            f"dataset {label!r} declares the {family.NAME!r} likelihood family, which this "
            f"backend does not lower into a differentiable density. Families it does: {known}. "
            f"The gradient-free engines run every implemented family, on this problem as "
            f"declared."
        ),
    )
