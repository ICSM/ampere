"""The likelihood families, transcribed into torch for the realised path.

Why a transcription and not a call
----------------------------------
``ampere.core``'s families are the definition of these quantities and this
module changes none of them: every function here is one of
``ampere.core.likelihood``'s ``log_prob`` bodies written in torch, line for
line, and the conformance suite compares the two at many points
(``inference.md`` §10a, "What the conformance suite owes").

They cannot simply be *called*, for the reason §10a gives: a
``LikelihoodFamily.log_prob`` returns ``float(...)`` of a numpy expression
built from scipy distributions, so a tensor entering one leaves as a number
and the graph is cut. The same is true of ``NoiseModel.sigma``, of
``Kernel.matrix`` and of the containers; each is transcribed where the
gradient has to survive, and nowhere else.

What is here and what is refused
--------------------------------
W2.13 fixed the coverage floor at the Gaussian family alone. **W2.4 slice 2
widens it to every family ``ampere.core`` implements for real data**:

====================  ====================================================
``gaussian``          uncensored and censored (the Tobit form)
``student_t``         uncensored; censored is refused, see below
``cauchy``            uncensored and censored
``poisson``           counts, and the **latent** form ``rate * exp(f)``
``complex_gaussian``  circular complex, uncorrelated (W2.4 slice 3)
====================  ====================================================

``complex_gaussian`` and the circular complex GP
------------------------------------------------
W2.4 slice 2 refused this family outright, on the grounds that the realised
path carried real tensors from the model's ``flux`` to the observed values.
Slice 3 makes both ends complex — :mod:`ampere.backends.torch.problem` keeps a
complex observed tensor where the container is complex, and a model whose
``flux`` returns a complex tensor composes — so the family is a transcription
like any other:

``log p = sum[ -|y - mu|**2 / (2 sigma**2) - log(2 pi) - log(sigma**2) ]``,

``ampere.core.ComplexGaussianFamily.log_prob``'s closed form line for line,
with ``torch.abs`` doing the work ``np.abs`` does there. ``sigma`` stays
**real**: it is the per-component standard deviation of a circular complex
Gaussian, which is what ``results_schema.md`` §16 says a
:class:`~ampere.core.VisibilitySet`'s real-valued uncertainty encodes. The
normalisation looks as though it is missing a half and is not — two
independent real Gaussians contribute ``-log(2 pi sigma**2)/2`` each, which is
the ``-log(2 pi) - log(sigma**2)`` written above.

**The correlated case no longer reaches this module at all (W4.2).**
``likelihoods.md`` §4 and the family's own
:attr:`~ampere.core.ComplexGaussianFamily.ANALYTIC_WITH_GP` declare that a
circular complex GP *does* marginalise in closed form — one real kernel
applied independently to the real and imaginary parts, equal component
covariances and zero pseudo-covariance, so the marginal likelihood is the sum
of two real Gaussian marginals over the same ``K + diag(σ²)``. Until W4.2 the
family declared ``GP_ANALYTIC_IMPLEMENTED = False`` and :func:`refuse_family`
refused the combination here, because ``ampere.core`` refused it first and
``inference.md`` §10a makes the numpy path the thing a realisation is checked
against: a native path with no counterpart is a backend inventing a likelihood.

W4.2 implemented the closed form in ``ampere.core``, and it landed where the
Gaussian family's GP already lives — on the **solver**, not here. A family whose
``ANALYTIC_WITH_GP`` is true hands its whole covariance to
``DenseGP.log_marginal_likelihood_native``
(:mod:`ampere.backends.torch.problem`'s ``gp_marginal`` branch), and what makes
the complex case work is that the branch stacks the complex residual into the
two real columns the solver now accepts. So the density below stays the
uncorrelated one, and :func:`refuse_family`'s ``GP_ANALYTIC_IMPLEMENTED`` guard
stays as the **generic** staging check for whatever family next declares a
closed form before writing it.

Two families stay refused, by name, at construction:

``rice``, ``von_mises``
    Declared slots with no implementation anywhere yet, including the
    reference path. A backend that implemented one first would be inventing
    the definition rather than transcribing it.

**Censored Student-t is refused specifically**, and the reason is worth
stating because it is a library gap rather than a decision: the Tobit form
needs the family's log-CDF, and torch has no incomplete beta function
(``torch.special`` offers no ``betainc``), so the Student-t CDF cannot be
written here without implementing a continued-fraction expansion whose
accuracy would then be ampere's problem. ``scipy`` has it, so the *contract*
path computes censored Student-t perfectly well — a gradient-free engine runs
that problem today. The refusal names the remedy.

Censoring, and where the branch lives
-------------------------------------
A censoring declaration is fixed at composition: which samples are limits is
a fact about the data, not about θ. So the three groups are selected once,
with boolean index masks built in :mod:`ampere.backends.torch.problem`, and
each group's contribution is summed separately. That is a Python branch on a
**constant**, which trace purity permits and jax's ``lax.cond`` would not
need either; the rule §10a states is that nothing branches on a *value*.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping
from typing import Any

import torch

from ampere.core import LimitKind
from ampere.core.exceptions import LoweringError

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor

__all__ = ["NATIVE_FAMILIES", "FamilyInputs", "native_log_prob", "refuse_family"]

_LOG_2PI = math.log(2.0 * math.pi)
_LOG_PI = math.log(math.pi)


@dataclasses.dataclass(frozen=True)
class FamilyInputs:
    """Everything a native family body reads, assembled once per evaluation.

    A struct rather than eight positional arguments because the four bodies
    take *different* subsets of it — Poisson never looks at ``sigma``, only
    Poisson looks at ``latent`` — and a signature that lists them all is a
    signature whose argument order can silently rotate.
    """

    #: The instrument-transformed prediction, retained samples only.
    predicted: torch.Tensor
    #: The observed values, retained samples only, as a constant tensor.
    observed: torch.Tensor
    #: Effective per-sample sigma, or ``None`` where the family declares
    #: ``REQUIRES_UNCERTAINTY = False`` and the data carry none.
    sigma: torch.Tensor | None
    #: The likelihood-wide resolved values (family and noise parameters,
    #: flat), exactly as ``Likelihood.log_prob`` hands them to a family.
    values: Mapping[str, Any]
    #: Index masks for the three limit kinds, or ``None`` when the dataset
    #: declares no censoring. Constants: see the module docstring.
    detection: torch.Tensor | None = None
    upper: torch.Tensor | None = None
    lower: torch.Tensor | None = None
    #: The latent GP values ``f``, for a family that consumes them.
    latent: torch.Tensor | None = None
    #: The family object, for its own parameters.
    family: Any = None
    #: Where a scalar this body needs to build (a Student-t ``nu``, a ``-inf``)
    #: is built. The lowering's own, so a problem composed on a GPU never
    #: silently materialises a CPU constant in the middle of its density
    #: (W2.4 slice 3).
    dtype: torch.dtype = DEFAULT_DTYPE
    device: torch.device = DEFAULT_DEVICE


def _own_values(inputs: FamilyInputs) -> Mapping[str, Any]:
    """The family's own parameters, resolved, without coercing to float.

    The filtering is ``ampere.core``'s established idiom — ``Likelihood``
    passes the flat, likelihood-wide values and each part completes its
    declaration against the names it owns — and the *absence* of a
    ``float()`` is what makes ``nu`` differentiable here where the contract
    path's ``_positive`` check makes it a number.
    """
    family = inputs.family
    own = {key: value for key, value in inputs.values.items() if key in family.parameters}
    return family.context(own)


def _scalar(inputs: FamilyInputs, value: Any) -> torch.Tensor:
    """*value* as a tensor where this evaluation's arithmetic is happening."""
    return as_tensor(value, dtype=inputs.dtype, device=inputs.device)


def _censored_total(
    inputs: FamilyInputs,
    log_density: torch.Tensor,
    log_cdf: torch.Tensor,
    log_sf: torch.Tensor,
) -> torch.Tensor:
    """``ampere.core.likelihood._location_scale_log_prob``, in torch.

    "A detection contributes ``log f(z) - log sigma``; an upper limit
    contributes ``log F(z)``, the probability that the truth lies below the
    recorded value; a lower limit contributes ``log (1 - F(z))``." The Tobit
    construction, transcribed — and note that *log_density* already carries
    the ``- log sigma``, exactly as the core's ``logpdf(z) - log(sigma)``
    does, while the two tail terms do not: a limit's contribution is a
    probability, not a density, and has no Jacobian.
    """
    if inputs.detection is None:
        return log_density.sum()
    total = log_density[inputs.detection].sum()
    if inputs.upper is not None:
        total = total + log_cdf[inputs.upper].sum()
    if inputs.lower is not None:
        total = total + log_sf[inputs.lower].sum()
    return total


def _standardise(inputs: FamilyInputs) -> tuple[torch.Tensor, torch.Tensor]:
    """``(observed - predicted) / sigma`` and ``log sigma``."""
    sigma = inputs.sigma
    if sigma is None:  # pragma: no cover - composition refuses this first
        raise LoweringError(
            "uncertainty",
            backend=BACKEND,
            detail="this family needs per-sample uncertainties and none were supplied.",
        )
    return (inputs.observed - inputs.predicted) / sigma, torch.log(sigma)


def _gaussian(inputs: FamilyInputs) -> torch.Tensor:
    """``GaussianFamily.log_prob``'s uncorrelated branch, censoring included.

    The uncensored case is written out rather than assembled from a
    distribution object for the reason the core gives for writing it out in
    numpy: it is the hot loop of every ordinary fit ampere runs.
    """
    z, log_sigma = _standardise(inputs)
    density = -0.5 * (z * z + _LOG_2PI) - log_sigma
    if inputs.detection is None:
        return density.sum()
    # log_ndtr rather than log(Phi): the whole point of a limit is that it
    # lives in a tail, where Phi(z) underflows float64 around z = -38 and the
    # log of the underflow is -inf rather than the -725 it should be.
    return _censored_total(inputs, density, torch.special.log_ndtr(z), torch.special.log_ndtr(-z))


def _student_t(inputs: FamilyInputs) -> torch.Tensor:
    """``StudentTFamily.log_prob``, uncensored (see the module docstring)."""
    z, log_sigma = _standardise(inputs)
    nu = _scalar(inputs, _own_values(inputs)["nu"])
    density = (
        torch.lgamma(0.5 * (nu + 1.0))
        - torch.lgamma(0.5 * nu)
        - 0.5 * torch.log(nu * math.pi)
        - 0.5 * (nu + 1.0) * torch.log1p(z * z / nu)
        - log_sigma
    )
    return density.sum()


def _cauchy(inputs: FamilyInputs) -> torch.Tensor:
    """``CauchyFamily.log_prob``: Student-t at ``nu = 1``, with a closed CDF.

    ``atan2(1, -z) / pi`` is the standard Cauchy CDF written so that it stays
    accurate in both tails: the naive ``0.5 + atan(z)/pi`` loses every
    significant digit of the lower tail to cancellation, which is precisely
    where an upper limit lives.
    """
    z, log_sigma = _standardise(inputs)
    density = -_LOG_PI - torch.log1p(z * z) - log_sigma
    if inputs.detection is None:
        return density.sum()
    ones = torch.ones_like(z)
    return _censored_total(
        inputs,
        density,
        torch.log(torch.atan2(ones, -z)) - _LOG_PI,
        torch.log(torch.atan2(ones, z)) - _LOG_PI,
    )


def _poisson(inputs: FamilyInputs) -> torch.Tensor:
    """``PoissonFamily.log_prob``, including the latent-GP form.

    ``counts ~ Poisson(rate * exp(f))`` is the case ``DEVELOPMENT_PLAN.md``
    §4.4 singles out as "effectively a modern-backend capability": there is no
    covariance matrix to marginalise, so the GP enters as N latent values that
    a gradient-based engine samples. This is the family that consumes them.

    A non-positive rate is a ``-inf`` rather than a raise
    (:attr:`inference.md` §10a), computed through :func:`torch.where` so the
    result stays attached to the graph.
    """
    rate = inputs.predicted
    if inputs.latent is not None:
        rate = rate * torch.exp(inputs.latent)
    positive = rate > 0.0
    safe = torch.where(positive, rate, torch.ones_like(rate))
    counts = inputs.observed
    terms = counts * torch.log(safe) - safe - torch.lgamma(counts + 1.0)
    return torch.where(positive.all(), terms.sum(), torch.full_like(terms.sum(), -math.inf))


def _complex_gaussian(inputs: FamilyInputs) -> torch.Tensor:
    """``ComplexGaussianFamily.log_prob``, transcribed. The circular complex model.

    The residual is complex and its modulus is what enters the exponent;
    ``sigma`` is real, and is the standard deviation of *each* component. Both
    facts are ``likelihoods.md`` §4's, and the density is the core's expression
    with ``torch.abs`` in place of ``np.abs``:

    ``-|y - mu|**2 / (2 sigma**2) - log(2 pi) - log(sigma**2)``.

    Censoring never reaches here: a limit on a complex value is not defined, so
    the family declares ``SUPPORTS_CENSORING = False`` and ``ampere.core``
    refuses the declaration when the ``Likelihood`` is built. Unlike the real
    Gaussian, then, there is no Tobit branch and no ``detection`` mask to
    consult.
    """
    sigma = inputs.sigma
    if sigma is None:  # pragma: no cover - composition refuses this first
        raise LoweringError(
            "uncertainty",
            backend=BACKEND,
            detail="the complex_gaussian family needs per-sample uncertainties.",
        )
    residual = torch.abs(inputs.observed - inputs.predicted)
    variance = sigma * sigma
    return torch.sum(-(residual * residual) / (2.0 * variance) - _LOG_2PI - torch.log(variance))


#: Family neutral name -> its torch body. The table a lowering consults, and
#: the one place a new family joins the realised path: add a body, add a row,
#: and the conformance suite compares it against ``ampere.core``'s at every
#: point it checks the realisation.
NATIVE_FAMILIES: dict[str, Any] = {
    "gaussian": _gaussian,
    "student_t": _student_t,
    "cauchy": _cauchy,
    "poisson": _poisson,
    "complex_gaussian": _complex_gaussian,
}

#: Families whose torch body cannot consume a censoring declaration, and the
#: reason, phrased for the user who hits it. Keyed by neutral name.
_CENSORING_REFUSALS: dict[str, str] = {
    "student_t": (
        "the Tobit form needs the Student-t log-CDF, which is a regularised incomplete beta "
        "function; torch.special has none, and writing a continued-fraction expansion here "
        "would make its accuracy ampere's problem rather than a library's. scipy has it, so "
        "the contract path computes this likelihood exactly — run it on a gradient-free "
        "engine (emcee, dynesty, zeus), or use the cauchy family, whose CDF is closed-form "
        "and which is the nu = 1 member of the same family"
    ),
}


def refuse_family(
    family: Any, *, censored: bool, latent: bool, correlated: bool = False
) -> LoweringError | None:
    """Whether this backend can lower *family*, and why not when it cannot.

    Returns the refusal rather than raising it, so the caller —
    :class:`ampere.backends.torch.problem._LoweredDataset` — raises every one
    of its refusals from one place, at construction, in the order that gives
    the most useful first message.

    *correlated* says whether the dataset's noise model induces correlations. It
    was added for one family — a ``complex_gaussian`` under a
    :class:`~ampere.core.GaussianProcessNoise`, declared analytic and unwritten
    — and **W4.2 wrote that closed form**, so the guard below now has no shipped
    instance and is kept as the generic one: any family that declares
    ``ANALYTIC_WITH_GP`` without ``GP_ANALYTIC_IMPLEMENTED`` is refused here,
    because a realisation is checked against the numpy path and there would be
    nothing to check it against. See the module docstring.
    """
    name = family.NAME
    if correlated and not bool(getattr(family, "GP_ANALYTIC_IMPLEMENTED", True)):
        return LoweringError(
            name,
            backend=BACKEND,
            detail=(
                f"the {name!r} family declares that it marginalises a Gaussian process "
                f"analytically (ANALYTIC_WITH_GP) but that the closed form is not implemented "
                f"yet (GP_ANALYTIC_IMPLEMENTED is False), so ampere.core has no numpy path for "
                f"this combination. A realisation is checked against that path (inference.md "
                f"§10a), so lowering it here would be this backend inventing a likelihood "
                f"nothing could check. Implement log_prob's correlated branch in ampere.core "
                f"and set GP_ANALYTIC_IMPLEMENTED = True, or use IndependentNoise with this "
                f"family."
            ),
        )
    if name not in NATIVE_FAMILIES:
        known = ", ".join(sorted(NATIVE_FAMILIES))
        return LoweringError(
            name,
            backend=BACKEND,
            detail=(
                f"this backend has no native torch body for the {name!r} likelihood family, so "
                f"it cannot be composed into a differentiable log-density. Families it does "
                f"lower: {known}. The gradient-free engines run every family, on this problem "
                f"as declared."
            ),
        )
    if censored and name in _CENSORING_REFUSALS:
        return LoweringError(
            name,
            backend=BACKEND,
            detail=(
                f"the {name!r} family is lowered here, but not together with a censoring "
                f"declaration: {_CENSORING_REFUSALS[name]}."
            ),
        )
    if latent and not family.CONSUMES_LATENT_GP:
        return LoweringError(
            name,
            backend=BACKEND,
            detail=(
                f"a latent GP was declared, but the {name!r} family does not consume latent "
                f"values (LikelihoodFamily.CONSUMES_LATENT_GP is False). ampere.core refuses "
                f"this combination too; the flag exists because a family that ignored the "
                f"latent values would return the uncorrelated likelihood and every value the "
                f"engine sampled would leave the log-probability untouched."
            ),
        )
    return None


def native_log_prob(inputs: FamilyInputs) -> torch.Tensor:
    """Dispatch to *inputs*'s family body. The refusals happened at construction."""
    return NATIVE_FAMILIES[inputs.family.NAME](inputs)


def limit_masks(
    kinds: Any, *, device: torch.device = DEFAULT_DEVICE
) -> tuple[torch.Tensor | None, torch.Tensor | None, torch.Tensor | None]:
    """The three constant index masks a censoring declaration reduces to.

    ``None`` for a group with no members, so a dataset whose limits were all
    masked away costs nothing and takes the uncensored path — ``likelihoods.md``
    §9's "masking beats censoring", arriving here as a shape rather than as a
    branch.

    *device* is the lowering's, because these masks index tensors that live
    there: a CPU mask against a CUDA density is an error torch raises once per
    evaluation rather than once at composition.
    """
    codes = torch.as_tensor(kinds, dtype=torch.int64, device=device)

    def group(kind: LimitKind) -> torch.Tensor | None:
        mask = codes == int(kind)
        return mask if bool(mask.any()) else None

    detection = group(LimitKind.DETECTION)
    upper = group(LimitKind.UPPER_LIMIT)
    lower = group(LimitKind.LOWER_LIMIT)
    if upper is None and lower is None:
        return None, None, None
    if detection is None:
        detection = torch.zeros_like(codes, dtype=torch.bool)
    return detection, upper, lower
