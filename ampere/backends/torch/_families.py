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

===================  ====================================================
``gaussian``         uncensored and censored (the Tobit form)
``student_t``        uncensored; censored is refused, see below
``cauchy``           uncensored and censored
``poisson``          counts, and the **latent** form ``rate * exp(f)``
===================  ====================================================

Two families stay refused, by name, at construction:

``complex_gaussian``
    Its data are complex, and the realised path carries real tensors from the
    model's ``flux`` to the observed values. Making it complex end to end is
    the visibility modality's work (plan §5, Phase 4), not a line in this
    module; ``ampere.core`` itself declares ``GP_ANALYTIC_IMPLEMENTED =
    False`` for the same combination and for the same reason.
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


def _scalar(value: Any) -> torch.Tensor:
    return as_tensor(value, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


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
    nu = _scalar(_own_values(inputs)["nu"])
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


#: Family neutral name -> its torch body. The table a lowering consults, and
#: the one place a new family joins the realised path: add a body, add a row,
#: and the conformance suite compares it against ``ampere.core``'s at every
#: point it checks the realisation.
NATIVE_FAMILIES: dict[str, Any] = {
    "gaussian": _gaussian,
    "student_t": _student_t,
    "cauchy": _cauchy,
    "poisson": _poisson,
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


def refuse_family(family: Any, *, censored: bool, latent: bool) -> LoweringError | None:
    """Whether this backend can lower *family*, and why not when it cannot.

    Returns the refusal rather than raising it, so the caller —
    :class:`ampere.backends.torch.problem._LoweredDataset` — raises every one
    of its refusals from one place, at construction, in the order that gives
    the most useful first message.
    """
    name = family.NAME
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


def limit_masks(kinds: Any) -> tuple[torch.Tensor | None, torch.Tensor | None, torch.Tensor | None]:
    """The three constant index masks a censoring declaration reduces to.

    ``None`` for a group with no members, so a dataset whose limits were all
    masked away costs nothing and takes the uncensored path — ``likelihoods.md``
    §9's "masking beats censoring", arriving here as a shape rather than as a
    branch.
    """
    codes = torch.as_tensor(kinds, dtype=torch.int64, device=DEFAULT_DEVICE)

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
