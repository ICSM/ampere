"""The plotting surface, declared once against the single results format.

``DEVELOPMENT_PLAN.md`` §4.6: "corner/trace/posterior-predictive plotting is
written once against it in ``ampere.results``... This is how the ``mixins.py``
monolith is retired and how sampler plotting parity stops regressing." The
monolith's defect was that each sampler driver carried its own plotting, so a
fix to one never reached the others. The cure is not a tidier monolith; it is
that **every function here takes the emitted run and nothing else** — no
sampler, no problem-specific state, no per-engine variant. A function that
needed to know which sampler produced its input would be reintroducing the bug.

Every function below is a declared signature with no implementation. Phase 2
lands them; W1.8's job is to fix the surface, its arguments and its obligations
first, so two backend tracks and the diagnostics families all write against one
already-agreed API. ``diagnostics.md`` §6 places families B (residual whiteness,
posterior-predictive checks) and C (GP localisation) in this namespace, and
their entry points are here beside the general-purpose ones.

Two obligations that are contract, not style
--------------------------------------------
1. **Every GP-localisation plot carries the degeneracy caveat by construction**
   (``diagnostics.md`` §4.3): "a plotting-function requirement for W1.8, not a
   'please remember to mention this' note". :data:`GP_LOCALISATION_CAVEAT` is
   that text, :func:`plot_gp_localisation` is required to render it, and
   :func:`gp_localisation_caveat` exposes it to a caller who is extracting
   numbers rather than looking at a figure.
2. **An anomaly score is never rendered without its provenance**
   (``diagnostics.md`` §5). One renderer serves both families, precisely so the
   two are read in one visual grammar — and the comparability risk that shared
   grammar creates is policed by showing where each score came from, which is
   why ``provenance`` is a required field of the input rather than an optional
   label.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any, Protocol, runtime_checkable

import numpy as np

__all__ = [
    "GP_LOCALISATION_CAVEAT",
    "AnomalyScoreLike",
    "gp_localisation_caveat",
    "plot_anomaly_score",
    "plot_corner",
    "plot_gp_localisation",
    "plot_posterior_predictive",
    "plot_residuals",
    "plot_trace",
]

#: The mandatory interpretation caveat on every GP-localisation output.
#:
#: ``prior_art.md`` §4 lesson S4, carried into ``diagnostics.md`` §4.3 by name:
#: Starfish found its global kernel amplitude and its explicit local kernels
#: traded off against one another on real data. Ampere has no separate local
#: components to trade off, but the ambiguity survives in another form, and it
#: is the difference between a diagnostic and a conclusion.
GP_LOCALISATION_CAVEAT = (
    "A large fitted GP amplitude localises where the model is deficient; it does not say why. "
    "Model error, an underestimated noise budget and a genuinely correlated astrophysical "
    "process are all consistent with the same posterior shape, and a large amplitude at a short "
    "length-scale may equally mean that the kernel's smooth global component is under-amplitude "
    "and compensating locally."
)

_PHASE_2 = (
    "Its signature and obligations are fixed by W1.8 "
    "(docs/design/contracts/results.md §8); the drawing lands in Phase 2."
)


def gp_localisation_caveat() -> str:
    """:data:`GP_LOCALISATION_CAVEAT`, for a programmatic consumer.

    ``diagnostics.md`` §4.3 requires the caveat to reach a user who extracts the
    score rather than viewing the plot — "a caption a user can silently crop out
    of a screenshot is not durable protection against over-interpretation".

    >>> gp_localisation_caveat().startswith("A large fitted GP amplitude localises")
    True
    """
    return GP_LOCALISATION_CAVEAT


@runtime_checkable
class AnomalyScoreLike(Protocol):
    """The coordinate-indexed deficiency map both diagnostic families produce.

    ``diagnostics.md`` §5 proposes an ``AnomalyScore`` container in
    ``ampere.core`` so that ``ampere.diagnostics`` (family A, pre-fit RHMF) and
    ``ampere.results`` (family C, post-fit GP localisation) can each produce one
    without either namespace depending on the other. That class does not exist —
    W1.4 landed before the proposal — so the renderer here is typed against the
    *shape* instead. Nothing changes at the call site when the class lands;
    ``docs/design/contracts/results.md`` §13 asks W1.13 to land it and this
    protocol then describes it exactly.

    ``provenance`` and ``interpretation_notes`` are **not** optional: they are
    what keeps a shared visual grammar from implying a comparability the two
    statistics do not have.
    """

    coordinates: Any
    values: np.ndarray
    mask: np.ndarray | None
    provenance: str
    interpretation_notes: str


# ---------------------------------------------------------------------------
# General-purpose posterior plots
# ---------------------------------------------------------------------------


def plot_corner(
    tree: Any,
    *,
    var_names: Sequence[str] | None = None,
    group: str = "posterior",
    labels: Sequence[str] | None = None,
    truths: Any = None,
    **kwargs: Any,
) -> Any:
    """The pairwise marginal grid.

    ``var_names`` selects merged parameter names (``model.index``,
    ``blue.instrument.calibrate.scale``); the merged name *is* the axis label,
    which is one of the three reasons ``inference.md`` §4.5 gives for the nested
    merge topology — ``sed.instrument.calibrate.scale`` says where to look and a
    flattened name does not.

    An array-valued parameter (a plate member, a latent GP block) is one
    variable with a named dimension, not ``N`` scalars, so selecting it selects
    the whole block; pass an ArviZ selection to narrow it. A corner plot of a
    10⁵-element latent block is a mistake this function should refuse loudly
    rather than attempt.
    """
    raise NotImplementedError(f"plot_corner is not implemented yet. {_PHASE_2}")


def plot_trace(
    tree: Any,
    *,
    var_names: Sequence[str] | None = None,
    group: str = "posterior",
    combined: bool = False,
    **kwargs: Any,
) -> Any:
    """Per-chain traces and marginals, the convergence eyeball.

    Reads ``sample_stats`` alongside the posterior, so a run whose draws include
    prior-rejected points shows them: those have ``lp = -inf`` and a **NaN**
    ``log_likelihood``, and rendering NaN as a gap rather than as zero is the
    visible half of ``inference.md`` §18(c)'s distinction between "not
    evaluated" and "impossible".
    """
    raise NotImplementedError(f"plot_trace is not implemented yet. {_PHASE_2}")


# ---------------------------------------------------------------------------
# Family B — post-fit residual whiteness and posterior-predictive checks
# ---------------------------------------------------------------------------


def plot_posterior_predictive(
    tree: Any,
    *,
    datasets: Sequence[str] | None = None,
    statistic: Any = None,
    **kwargs: Any,
) -> Any:
    """Replicate data against the observations — ``diagnostics.md`` family B.

    Consumes the ``posterior_predictive`` group. That group is **not** stored by
    default (``ampere.results.derived`` says why, and holds the on-demand
    builder), so this function's precondition is that
    :func:`~ampere.results.derived.add_posterior_predictive` has been called, and
    its refusal must say so rather than reporting a missing group.

    ``statistic`` is the discrepancy measure ``T(y, θ)`` whose posterior
    predictive p-value is reported; the default is a standardised-residual sum,
    and the Ljung-Box statistic of :func:`plot_residuals` is a legitimate
    alternative, which is the tie between the two halves of family B.
    """
    raise NotImplementedError(f"plot_posterior_predictive is not implemented yet. {_PHASE_2}")


def plot_residuals(
    tree: Any,
    *,
    datasets: Sequence[str] | None = None,
    whiteness: bool = True,
    **kwargs: Any,
) -> Any:
    """Signed residuals against coordinate, plus the whiteness panel.

    Consumes the ``residuals`` group (see
    :func:`~ampere.results.derived.add_residuals`). Signed, because
    autocorrelation is sign-sensitive and per-observation log-likelihood is not
    — ``diagnostics.md`` §3.3's precise gap.

    ``diagnostics.md`` §3.1 scopes family B to **standard-likelihood** fits: a
    GP-augmented fit's residuals are whitened by construction, so testing them
    for whiteness is close to circular. This function is therefore required to
    warn when the run's likelihood provenance says a GP was fitted, and to point
    at :func:`plot_gp_localisation` instead.
    """
    raise NotImplementedError(f"plot_residuals is not implemented yet. {_PHASE_2}")


# ---------------------------------------------------------------------------
# Family C — GP localisation
# ---------------------------------------------------------------------------


def plot_gp_localisation(
    tree: Any,
    *,
    datasets: Sequence[str] | None = None,
    band: float = 0.68,
    show_caveat: bool = True,
    **kwargs: Any,
) -> Any:
    """Where the flexible likelihood says the model is deficient.

    The signed conditioned GP mean with its posterior band, on the data's own
    coordinate axis — **masked samples included**, because
    :meth:`ampere.core.likelihood.Likelihood.conditional` defaults to the full
    axis for exactly this plot: "what would the GP have said here?" is the
    question a user asks about a region they excluded.

    ``show_caveat`` is a keyword rather than a fact of the caller's choosing:
    ``diagnostics.md`` §4.3 makes the caveat mandatory, so setting it ``False``
    must still leave the caveat on the returned figure's metadata and in this
    function's own docstring, and only suppresses the drawn annotation.

    Caveat rendered on every such plot
    ----------------------------------
    A large fitted GP amplitude localises where the model is deficient; it does
    not say why. Model error, an underestimated noise budget and a genuinely
    correlated astrophysical process are all consistent with the same posterior
    shape, and a large amplitude at a short length-scale may equally mean that
    the kernel's smooth global component is under-amplitude and compensating
    locally. See :data:`GP_LOCALISATION_CAVEAT`.
    """
    raise NotImplementedError(f"plot_gp_localisation is not implemented yet. {_PHASE_2}")


def plot_anomaly_score(
    score: AnomalyScoreLike,
    *,
    ax: Any = None,
    show_provenance: bool = True,
    **kwargs: Any,
) -> Any:
    """One renderer for both families' deficiency maps.

    ``diagnostics.md`` §5's resolved Tension 5: pre-fit RHMF flags and post-fit
    GP localisation share a visual grammar — a sequential scale, a documented
    range, "higher = worse" in both — so that moving from "should I screen this
    collection?" to "did my fit's residuals agree?" costs no re-reading. What
    the shared grammar must not do is imply the two statistics are
    interchangeable, and the guard against that is metadata rather than style:
    ``score.provenance`` and ``score.interpretation_notes`` are displayed, and
    ``show_provenance=False`` is not permitted to suppress them when two scores
    of different provenance are drawn together.

    This function deliberately does not import :mod:`ampere.diagnostics`: family
    A's namespace carries a JAX dependency, and rendering a score somebody hands
    over must not drag it in.
    """
    raise NotImplementedError(f"plot_anomaly_score is not implemented yet. {_PHASE_2}")
