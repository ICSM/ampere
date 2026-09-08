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

import warnings
from collections.abc import Sequence
from typing import Any, Protocol, runtime_checkable

import numpy as np
import scipy.stats as st

from ampere.core.exceptions import ResultsError

from . import _plotting as _p
from .derived import GP_LOCALISATION_GROUP, RESIDUALS_GROUP
from .diagnostics import residual_whiteness
from .provenance import ATTR_PREFIX

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

    :class:`ampere.core.AnomalyScore` is the class — ``diagnostics.md`` §5's
    proposal, landed at the freeze (ruled 2026-09-03, ``results.md`` §15 R4)
    in ``ampere.core`` so that ``ampere.diagnostics`` (family A, pre-fit
    RHMF) and ``ampere.results`` (family C, post-fit GP localisation) can
    each produce one without either namespace depending on the other. This
    protocol describes it exactly, and the renderer here stays typed against
    the *shape*: a caller may hand it the real class or anything matching.

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

    Each dataset gets a residual panel — the across-draw median with a 68 %
    band, on the data's own coordinate axis, masked samples left as gaps — and,
    unless ``whiteness=False``, a second panel holding the separation-binned
    autocorrelation with its permutation envelope and the ``(Q, p)`` the test
    reports. The statistic is defined in :mod:`ampere.results.diagnostics`;
    ``**kwargs`` are its knobs (``bins``, ``n_permutations``, ``max_separation``,
    ``max_draws``, ``seed``), so the figure and the number can never disagree
    about what was computed.

    Returns
    -------
    matplotlib.figure.Figure
        With :func:`~ampere.results.figure_metadata` carrying each dataset's
        ``(statistic, p_value)``, so a caller reading numbers off the figure
        does not have to re-run the test to get them.
    """
    group = _p.require_group(
        tree,
        RESIDUALS_GROUP,
        remedy="call ampere.results.add_residuals(tree, problem) first, which derives the "
        "signed standardised residuals from the stored draws and the problem "
        "(results.md §7).",
    )
    available = [str(name) for name in group.data_vars]
    labels = _p.select_datasets(available, datasets, what=RESIDUALS_GROUP)
    fitted_with_gp = [label for label in labels if label in _p.gp_datasets(tree)]
    if fitted_with_gp:
        warnings.warn(
            f"dataset(s) {sorted(fitted_with_gp)} were fitted with a Gaussian-process noise "
            f"model, and diagnostics.md §3.1 scopes the residual-whiteness family to "
            f"standard-likelihood fits: a GP-augmented fit's residuals are whitened by "
            f"construction, so testing them for whiteness is close to circular. The question "
            f"'where is this model deficient?' is answered for a GP fit by "
            f"ampere.results.plot_gp_localisation instead.",
            UserWarning,
            stacklevel=2,
        )
    panels = 2 if whiteness else 1
    figure, axes = _p.new_axes(nrows=panels * len(labels))
    for index, label in enumerate(labels):
        axis_name, coordinates = _p.coordinate_of(group, label)
        coordinate_unit, value_unit = _p.dataset_units(tree, label, axis_name)
        values = np.asarray(group[label].values, dtype=float)
        flat = values.reshape(-1, values.shape[-1])
        panel = axes[panels * index]
        with warnings.catch_warnings():
            # An all-masked sample is NaN in every draw by construction; that is
            # a gap in the plot, not a numerical problem worth a warning.
            warnings.simplefilter("ignore", RuntimeWarning)
            centre = np.nanmedian(flat, axis=0)
            lower = np.nanpercentile(flat, 16.0, axis=0)
            upper = np.nanpercentile(flat, 84.0, axis=0)
        panel.axhline(0.0, color="0.6", linewidth=0.8)
        panel.fill_between(coordinates, lower, upper, alpha=0.3, label="68 % of draws")
        panel.plot(coordinates, centre, marker=".", linewidth=1.0, label="posterior median")
        panel.set_xlabel(_p.axis_label(axis_name, coordinate_unit))
        panel.set_ylabel(
            "standardised residual"
            if int(group.attrs.get(f"{ATTR_PREFIX}standardised", 1))
            else _p.axis_label("residual", value_unit)
        )
        panel.set_title(f"{label}: signed residuals")
        panel.legend(loc="best", fontsize="small")
        if not whiteness:
            continue
        test = residual_whiteness(tree, dataset=label, **kwargs)
        _draw_whiteness(axes[panels * index + 1], test, coordinate_unit)
        _p.attach_metadata(figure, f"{label}.whiteness_statistic", f"{test.statistic:.6g}")
        _p.attach_metadata(figure, f"{label}.whiteness_p_value", f"{test.p_value:.6g}")
    figure.tight_layout()
    return figure


def _draw_whiteness(panel: Any, test: Any, coordinate_unit: str) -> None:
    """The autocorrelation panel: rho_b, its pair counts, and the verdict.

    The permutation resolution is drawn on the caption rather than implied: a
    p-value of ``1 / (1 + n_permutations)`` is the *floor*, not evidence of
    extraordinary structure, and a reader who cannot see how many permutations
    were run cannot tell the two apart.
    """
    panel.axhline(0.0, color="0.6", linewidth=0.8)
    width = float(np.min(np.diff(test.separations))) * 0.6 if test.separations.size > 1 else None
    panel.bar(test.separations, test.autocorrelation, width=width, alpha=0.7)
    panel.set_xlabel(_p.axis_label("separation", coordinate_unit))
    panel.set_ylabel("binned autocorrelation")
    verdict = "structure detected" if test.p_value < 0.05 else "consistent with white"
    panel.set_title(
        f"{test.dataset}: Q = {test.statistic:.3g}, p = {test.p_value:.3g} — {verdict} "
        f"({test.n_permutations} permutations, floor p = {test.resolution:.3g}; "
        f"{test.n_pairs} pairs, {test.n_draws} draw(s))",
        fontsize="small",
    )


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

    The band combines both sources of uncertainty by the law of total variance
    — the mean of the conditional variances plus the variance of the
    conditional means across draws — because either alone understates what the
    posterior actually says about the GP's amplitude here.

    Caveat rendered on every such plot
    ----------------------------------
    A large fitted GP amplitude localises where the model is deficient; it does
    not say why. Model error, an underestimated noise budget and a genuinely
    correlated astrophysical process are all consistent with the same posterior
    shape, and a large amplitude at a short length-scale may equally mean that
    the kernel's smooth global component is under-amplitude and compensating
    locally. See :data:`GP_LOCALISATION_CAVEAT`.
    """
    if not 0.0 < band < 1.0:
        raise ResultsError(f"band is a credible-interval mass in (0, 1), got {band!r}.")
    group = _p.require_group(
        tree,
        GP_LOCALISATION_GROUP,
        remedy="call ampere.results.gp_localisation(tree, problem) first, which evaluates "
        "Likelihood.conditional across the stored draws (results.md §7).",
    )
    available = sorted({str(name).rsplit("_", 1)[0] for name in group.data_vars})
    labels = _p.select_datasets(available, datasets, what=GP_LOCALISATION_GROUP)
    figure, axes = _p.new_axes(nrows=len(labels))
    z = float(st.norm.ppf(0.5 + band / 2.0))
    for index, label in enumerate(labels):
        axis_name, coordinates = _p.coordinate_of(group, f"{label}_mean")
        coordinate_unit, value_unit = _p.dataset_units(tree, label, axis_name)
        means = np.asarray(group[f"{label}_mean"].values, dtype=float)
        variances = np.asarray(group[f"{label}_variance"].values, dtype=float)
        flat_mean = means.reshape(-1, means.shape[-1])
        flat_variance = variances.reshape(-1, variances.shape[-1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            centre = np.nanmedian(flat_mean, axis=0)
            total = np.nanmean(flat_variance, axis=0) + np.nanvar(flat_mean, axis=0)
        spread = z * np.sqrt(np.clip(total, 0.0, None))
        panel = axes[index]
        panel.axhline(0.0, color="0.6", linewidth=0.8)
        panel.fill_between(
            coordinates,
            centre - spread,
            centre + spread,
            alpha=0.3,
            label=f"{band:.0%} credible band",
        )
        panel.plot(coordinates, centre, linewidth=1.2, label="conditioned GP mean")
        panel.set_xlabel(_p.axis_label(axis_name, coordinate_unit))
        panel.set_ylabel(_p.axis_label("GP mean", value_unit))
        panel.set_title(f"{label}: GP localisation")
        panel.legend(loc="best", fontsize="small")
    # The caveat travels with the figure whatever the caller asked for. Only the
    # drawn annotation is optional — diagnostics.md §4.3's "a caption a user can
    # silently crop out of a screenshot is not durable protection".
    _p.attach_metadata(figure, "gp_localisation_caveat", GP_LOCALISATION_CAVEAT)
    figure.tight_layout()
    if show_caveat:
        figure.subplots_adjust(bottom=0.28)
        figure.text(0.01, 0.01, _p.wrap(GP_LOCALISATION_CAVEAT), fontsize="x-small", va="bottom")
    return figure


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

    "Drawn together" means "on one axes", and only the axes can know: the
    second call cannot see the first call's arguments. So the axes keeps the
    record, and the moment a second provenance lands on it the labels come back
    on — **for both scores**, including the one already drawn under
    ``show_provenance=False``.

    This function deliberately does not import :mod:`ampere.diagnostics`: family
    A's namespace carries a JAX dependency, and rendering a score somebody hands
    over must not drag it in.

    Parameters
    ----------
    score
        Anything matching :class:`AnomalyScoreLike` —
        :class:`ampere.core.AnomalyScore` or a caller's own shape.
    ax
        Draw on this axes rather than a new figure. Passing the same axes twice
        is how two scores are compared, and is what triggers the rule above.
    show_provenance
        Draw the provenance label and the interpretation notes. Ignored — with
        a warning — as soon as the axes holds two provenances.
    **kwargs
        Forwarded to the line artist.

    Returns
    -------
    matplotlib.axes.Axes
        The axes drawn on, with :func:`~ampere.results.figure_metadata` on its
        figure carrying every drawn score's provenance and interpretation
        notes, whatever ``show_provenance`` was.
    """
    if not isinstance(score, AnomalyScoreLike):
        raise ResultsError(
            f"plot_anomaly_score needs an AnomalyScoreLike — coordinates, values, mask, "
            f"provenance and interpretation_notes — and got a {type(score).__name__}. "
            f"ampere.core.AnomalyScore is the class; anything with those attributes will do."
        )
    coordinates = np.asarray(score.coordinates, dtype=float)
    if coordinates.ndim == 2 and coordinates.shape[1] != 1:
        raise ResultsError(
            f"this renderer draws a 1-D deficiency map, and the score is indexed by "
            f"{coordinates.shape[1]} coordinates. Multi-axis and gridded scores are the Phase 5 "
            f"staging (diagnostics.md §9 limitation 1)."
        )
    coordinates = coordinates.reshape(coordinates.shape[0])
    values = np.asarray(score.values, dtype=float)
    mask = None if score.mask is None else np.asarray(score.mask, dtype=bool)
    drawn = np.where(mask, np.nan, values) if mask is not None else values

    figure, axes = _p.new_axes(ax)
    panel = axes[0]
    line = panel.plot(coordinates, drawn, linewidth=1.2, **kwargs)[0]
    fill = panel.fill_between(coordinates, 0.0, np.nan_to_num(drawn), alpha=0.25)
    on_axes = _p.register_score(panel, score.provenance, [line, fill])
    provenances = _p.distinct_provenances(on_axes)
    forced = len(provenances) > 1
    if forced and not show_provenance:
        warnings.warn(
            "show_provenance=False was ignored: this axes now holds anomaly scores of "
            f"{len(provenances)} different provenances {list(provenances)}, and diagnostics.md "
            "§5 does not permit two differently-computed scores to be drawn together without "
            "their provenance shown — the shared visual grammar would otherwise imply a "
            "comparability the statistics do not have.",
            UserWarning,
            stacklevel=2,
        )
    if show_provenance or forced:
        for provenance, artists in on_axes:
            artists[0].set_label(provenance)
        panel.legend(loc="best", fontsize="small", title="provenance")
    panel.set_ylabel("anomaly score (higher = worse)")
    panel.set_xlabel("coordinate")
    if float(np.nanmin(values)) >= 0.0:
        panel.set_ylim(bottom=0.0)
    # Whatever was drawn, the notes reach a programmatic consumer.
    _p.attach_metadata(figure, f"anomaly_score.{score.provenance}", score.interpretation_notes)
    if show_provenance or forced:
        notes = "\n\n".join(
            f"{provenance}: {score.interpretation_notes}"
            if provenance == score.provenance
            else provenance
            for provenance in provenances
        )
        panel.set_title(_p.wrap(notes, width=110), fontsize="xx-small", loc="left")
    return panel
