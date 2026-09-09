"""The plotting surface, declared once against the single results format.

``DEVELOPMENT_PLAN.md`` §4.6: "corner/trace/posterior-predictive plotting is
written once against it in ``ampere.results``... This is how the ``mixins.py``
monolith is retired and how sampler plotting parity stops regressing." The
monolith's defect was that each sampler driver carried its own plotting, so a
fix to one never reached the others. The cure is not a tidier monolith; it is
that **every function here takes the emitted run and nothing else** — no
sampler, no problem-specific state, no per-engine variant. A function that
needed to know which sampler produced its input would be reintroducing the bug.

W1.8 declared every function below as a signature with no implementation, so
that two backend tracks and the diagnostics families could write against one
already-agreed API before anything drew a line. Phase 2 filled them in: **W2.7**
landed ``diagnostics.md``'s families B and C (:func:`plot_residuals`,
:func:`plot_gp_localisation`, :func:`plot_anomaly_score`) and **W2.8** the three
general-purpose ones (:func:`plot_corner`, :func:`plot_trace`,
:func:`plot_posterior_predictive`). Nothing in the surface moved to make that
possible, which is the whole point of having fixed it first. ``diagnostics.md``
§6 places families B (residual whiteness, posterior-predictive checks) and C
(GP localisation) in this namespace, and their entry points are here beside the
general-purpose ones. **W3.10** turned :func:`plot_corner` and
:func:`plot_trace`'s refusal above their caps into automatic paging with a
loud :class:`ResultsWarning` (``results.md`` §8, *Amended W3.10*), so that a
run with more variables than fit on one figure gets several instead of
nothing; ``paginate=False`` keeps the original one-figure-or-refuse
behaviour for a caller who wants it.

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

import math
import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import Any, Protocol, runtime_checkable

import numpy as np
import scipy.stats as st

from ampere.core.exceptions import ResultsError

from . import _plotting as _p
from .calibration import (
    CALIBRATION_GROUP,
    LEVEL_DIM,
    PARAMETER_DIM,
    SIMULATION_DIM,
    TARP_LEVEL_DIM,
)
from .derived import GP_LOCALISATION_GROUP, POSTERIOR_PREDICTIVE_GROUP, RESIDUALS_GROUP
from .diagnostics import residual_whiteness
from .provenance import ATTR_PREFIX

__all__ = [
    "GP_LOCALISATION_CAVEAT",
    "MAX_CORNER_VARIABLES",
    "MAX_RANK_PANELS",
    "MAX_TRACE_VARIABLES",
    "AnomalyScoreLike",
    "ResultsWarning",
    "gp_localisation_caveat",
    "plot_anomaly_score",
    "plot_corner",
    "plot_coverage",
    "plot_gp_localisation",
    "plot_posterior_predictive",
    "plot_residuals",
    "plot_sbc_ranks",
    "plot_trace",
]

#: How many scalar columns :func:`plot_corner` will draw on one figure before
#: paging.
#:
#: ``results.md`` §8 required "a corner plot of a 10⁵-element latent block" to
#: be refused loudly rather than attempted, and left the threshold to the
#: implementation; **W3.10** (ruled 2026-09-08, on W2.8's confirmed caps)
#: changed the refusal itself into automatic paging with a loud warning, so
#: what this constant now bounds is a *page*, not the whole figure. Twenty is
#: where a pairwise grid stops being readable — 400 panels — long before it
#: stops being *drawable*, and :func:`plot_corner` names ``max_variables=`` in
#: its warning so a user who genuinely wants thirty per page gets thirty by
#: saying so, or ``paginate=False`` for the one-figure-or-refuse behaviour this
#: constant originally described. There is no defensible value here that is
#: not a judgement; what is not a judgement is that the default must be small
#: enough to catch the mistake and overridable enough not to be an obstacle.
MAX_CORNER_VARIABLES = 20

#: The same guard for :func:`plot_trace`, where the cost is linear rather than
#: quadratic in the column count — so the threshold is higher, and it is still
#: a threshold: a figure of 200 stacked panels is not a diagnostic. Bounds a
#: page since **W3.10**, exactly as :data:`MAX_CORNER_VARIABLES` does.
MAX_TRACE_VARIABLES = 40


class ResultsWarning(UserWarning):
    """A plotting call was honoured, but not exactly as asked.

    **Landed at W3.10**, for the one case that surface has so far: paging
    above :data:`MAX_CORNER_VARIABLES` or :data:`MAX_TRACE_VARIABLES` instead
    of the refusal ``results.md`` §8 used to require (*Amended W3.10*, on
    W2.8's confirmed caps — Peter's ruling 2026-09-08 was explicit that
    "automatic paging with a loud warning is wanted later"). Splitting a
    corner or trace figure into several pages changes what a caller gets back
    — a :class:`list` of figures rather than one — silently enough that a
    script written against the single-figure return could misbehave without
    ever raising, which is exactly the shape of bug a warning exists to catch
    before it does. Deliberately a :class:`UserWarning` subclass and not an
    :class:`~ampere.core.exceptions.AmpereError`: nothing failed, and a
    caller who wants the old failure back for one call still has it via
    ``paginate=False``.
    """


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

#: Passed to :func:`~ampere.results._plotting.parameter_columns` in place of a
#: caller's ``max_variables`` when paging, so that function's own "refuse
#: above this many columns" checks never fire — pagination, not refusal, is
#: what W3.10 wants above the real cap. Not ``math.inf``: the parameter is
#: typed ``int`` and a run with more columns than this is not a thing that
#: happens.
_UNPAGED = 2**30


def _paginate_columns(
    columns: Sequence[tuple[str, np.ndarray]], *, limit: int
) -> list[list[tuple[str, np.ndarray]]]:
    """Split *columns* into pages of at most *limit*, in the given order.

    ``results.md`` §8 (*Amended W3.10*): pages are cut in merged-name order,
    "array blocks kept whole where they fit". A block here is a run of
    consecutive columns sharing one parameter name — an array-valued
    parameter's expanded elements, or a lone scalar's single column — so
    blocks are recovered from the column labels :func:`~ampere.results.
    _plotting.parameter_columns` already produced (``"name"`` for a scalar,
    ``"name[element]"`` per element of a block) rather than threaded through
    as separate state.

    A block that fits in a fresh page but not in the page currently being
    filled starts a new page rather than splitting — that is "kept whole
    where they fit". A block bigger than *limit* cannot fit any page whole,
    so it alone is split into consecutive full pages of exactly *limit*
    columns, in element order; nothing else the run wants next is added to
    those pages, which keeps the split block from picking up a stray
    unrelated column as a neighbour.
    """
    blocks: list[list[tuple[str, np.ndarray]]] = []
    for column in columns:
        name = column[0].partition("[")[0]
        if blocks and blocks[-1][0][0].partition("[")[0] == name:
            blocks[-1].append(column)
        else:
            blocks.append([column])
    pages: list[list[tuple[str, np.ndarray]]] = []
    current: list[tuple[str, np.ndarray]] = []
    for block in blocks:
        if len(block) > limit:
            if current:
                pages.append(current)
                current = []
            for start in range(0, len(block), limit):
                pages.append(block[start : start + limit])
            continue
        if current and len(current) + len(block) > limit:
            pages.append(current)
            current = []
        current.extend(block)
    if current:
        pages.append(current)
    return pages or [[]]


def _warn_paged(*, what: str, total: int, limit: int, n_pages: int) -> None:
    """The loud warning ``results.md`` §8 requires when paging fires.

    Names the page count, the cap and the ``var_names=`` route to a smaller
    figure — the three things the item text asks the warning to name — plus
    ``paginate=False`` for a caller who would rather have today's refusal
    back for this one call.
    """
    warnings.warn(
        f"{what} of {total} variables exceeds max_variables={limit}, so it was split into "
        f"{n_pages} pages instead of refused (results.md §8, amended W3.10); each page's "
        f"figure_metadata records which. Narrow the run with var_names=, raise max_variables= "
        f"to change how many pages are cut, or pass paginate=False for the old refusal.",
        ResultsWarning,
        stacklevel=3,
    )


def plot_corner(
    tree: Any,
    *,
    var_names: Sequence[str] | None = None,
    group: str = "posterior",
    labels: Sequence[str] | None = None,
    truths: Any = None,
    max_variables: int = MAX_CORNER_VARIABLES,
    paginate: bool = True,
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
    the whole block; its elements are labelled from the dimension's own
    coordinate where it has one — a plate's coordinate is the dataset labels
    (``hierarchical_population.md`` §10.2) — and by integer index where it does
    not.

    **Paging (W3.10).** More than ``max_variables`` columns used to be refused
    outright; now they are **paged** instead — split into consecutive figures
    each within the cap, in merged-name order, with an array-valued block kept
    whole on one page where it fits on one at all. A block bigger than
    ``max_variables`` on its own cannot fit any page whole, so it alone is
    split across full pages of exactly ``max_variables`` columns, in element
    order. Paging fires a :class:`ResultsWarning` naming the page count, the
    cap and the ``var_names=`` route to a smaller figure instead, and every
    page's :func:`~ampere.results.figure_metadata` carries ``"page"`` as
    ``"i of n"``. ``paginate=False`` restores the original behaviour exactly:
    one figure or a :class:`~ampere.core.exceptions.ResultsError` naming the
    variable, its size and ``max_variables``, for a caller who needs a single
    figure and would rather fail than receive a list.

    Prior-rejected draws are **excluded**, not plotted. Their θ is a perfectly
    good number and stored as one (``results.md`` §5), but the point carries
    zero prior mass and is not a posterior sample; leaving it in would put mass
    where the posterior has none, and silently.

    Parameters
    ----------
    tree
        The run, and nothing else (``results.md`` §8).
    var_names
        Merged parameter names to draw; all of the group's by default.
    group
        Which group to read. ``"posterior"`` normally; ``"prior"`` for a run
        that stored one.
    labels
        Axis labels, one per **column** after array-valued blocks have been
        expanded, across *every* page. The merged names are the default and
        are usually the right answer; a mismatched length is refused rather
        than silently truncated.
    truths
        Reference values: a mapping of merged name to value (arrays allowed,
        matching the block), or a sequence in column order across every page.
        ``None`` draws no reference lines.
    max_variables
        The per-page cap above. Raising it is a deliberate act.
    paginate
        Page above the cap (the default) rather than refuse. ``False``
        restores the pre-W3.10 refusal, unchanged, for a caller who needs
        exactly one figure or a hard failure.
    **kwargs
        Forwarded to :func:`corner.corner`.

    Returns
    -------
    matplotlib.figure.Figure | list[matplotlib.figure.Figure]
        A single :class:`~matplotlib.figure.Figure` — exactly today's return
        — for every call whose columns fit within ``max_variables``, paginated
        or not. Only a call that actually pages returns a :class:`list` of
        figures, in page order. Each figure carries
        :func:`~ampere.results.figure_metadata` with the column labels drawn
        on it and the number of draws actually used, so a caller can tell a
        thinned or prior-rejected-heavy run from a full one without
        re-reading it; a paged figure's metadata also carries ``"page"``.
    """
    corner = _p.require_corner()
    dataset = _p.require_sampling_group(tree, group)
    names = _p.select_names(
        [str(name) for name in dataset.data_vars], var_names, what=f"{group} variable"
    )
    keep = _p.scored_draws(tree)
    if not paginate:
        columns = _p.parameter_columns(dataset, names, limit=max_variables, what="a corner plot")
        drawn_labels = _resolve_labels(columns, labels)
        return _render_corner_page(
            corner,
            columns,
            keep=keep,
            drawn_labels=drawn_labels,
            truths=_corner_truths(truths, columns),
            page=None,
            **kwargs,
        )
    columns = _p.parameter_columns(dataset, names, limit=_UNPAGED, what="a corner plot")
    drawn_labels = _resolve_labels(columns, labels)
    resolved_truths = _corner_truths(truths, columns)
    pages = _paginate_columns(columns, limit=max_variables)
    if len(pages) == 1:
        return _render_corner_page(
            corner,
            pages[0],
            keep=keep,
            drawn_labels=drawn_labels,
            truths=resolved_truths,
            page=None,
            **kwargs,
        )
    _warn_paged(what="a corner plot", total=len(columns), limit=max_variables, n_pages=len(pages))
    figures = []
    offset = 0
    for index, page_columns in enumerate(pages, start=1):
        span = len(page_columns)
        figures.append(
            _render_corner_page(
                corner,
                page_columns,
                keep=keep,
                drawn_labels=drawn_labels[offset : offset + span],
                truths=None if resolved_truths is None else resolved_truths[offset : offset + span],
                page=(index, len(pages)),
                **kwargs,
            )
        )
        offset += span
    return figures


def _resolve_labels(
    columns: Sequence[tuple[str, np.ndarray]], labels: Sequence[str] | None
) -> list[str]:
    """The axis labels to draw, one per column across every page.

    Shared by :func:`plot_corner`'s paginated and unpaginated paths so a
    caller's ``labels=`` is validated once, against the *total* column count,
    rather than once per page — a mismatch is reported the same way whether
    the run pages or not.
    """
    drawn_labels = [label for label, _ in columns]
    if labels is None:
        return drawn_labels
    if len(labels) != len(drawn_labels):
        raise ResultsError(
            f"{len(labels)} label(s) were given for {len(drawn_labels)} column(s) "
            f"{drawn_labels}. An array-valued parameter is one variable and several "
            f"columns (results.md §4), so the labels are per column, not per variable."
        )
    return [str(label) for label in labels]


def _render_corner_page(
    corner: Any,
    columns: Sequence[tuple[str, np.ndarray]],
    *,
    keep: np.ndarray | None,
    drawn_labels: Sequence[str],
    truths: Any,
    page: tuple[int, int] | None,
    **kwargs: Any,
) -> Any:
    """One figure's worth of :func:`plot_corner`, unpaged or one page of many.

    Split out of :func:`plot_corner` so that the unpaged path (one call, one
    figure) and the paged path (one call per page) draw identically — a
    figure this function returns cannot tell which path produced it, which is
    the point: paging must not change what one page looks like, only how many
    there are. ``drawn_labels`` and ``truths`` arrive already resolved and
    sliced to *columns*; this function only draws.
    """
    samples = np.column_stack(
        [(values if keep is None else values[keep]).ravel() for _, values in columns]
    )
    if samples.shape[0] == 0:
        raise ResultsError(
            "every stored draw was rejected by the prior (lp = -inf), so there is no posterior "
            "to draw. That is a statement about the run, not about this plot: check the priors "
            "and the sampler's initialisation."
        )
    figure = corner.corner(samples, labels=list(drawn_labels), truths=truths, **kwargs)
    _p.attach_metadata(figure, "corner.variables", ", ".join(label for label, _ in columns))
    _p.attach_metadata(figure, "corner.draws", str(int(samples.shape[0])))
    if page is not None:
        index, total_pages = page
        _p.attach_metadata(figure, "page", f"{index} of {total_pages}")
    return figure


def _corner_truths(truths: Any, columns: Sequence[tuple[str, np.ndarray]]) -> Any:
    """Reference values in column order, from a mapping or a sequence.

    A mapping keyed by merged name is the form a user has to hand — it is what
    ``FittingProblem.simulate`` took and what ``Simulation.parameters`` gives
    back — and turning it into corner's positional list here is the difference
    between "plot the truth" and "count the columns of an expanded plate by
    hand". A name the run does not hold is refused: a truth silently dropped is
    a recovery plot that looks better than the fit was.
    """
    if truths is None or not isinstance(truths, Mapping):
        return truths
    values: list[float | None] = []
    for label, _ in columns:
        name, _, element = label.partition("[")
        if name not in truths:
            values.append(None)
            continue
        entry = np.asarray(truths[name], dtype=float)
        if not element:
            values.append(float(entry.reshape(())))
            continue
        flat = entry.ravel()
        index = _element_position(label, element.rstrip("]"), flat.size)
        values.append(float(flat[index]))
    unknown = sorted(set(truths) - {label.partition("[")[0] for label, _ in columns})
    if unknown:
        raise ResultsError(
            f"truths were given for {unknown}, which are not columns of this plot; the columns "
            f"are {[label for label, _ in columns]}. A truth quietly dropped is a recovery plot "
            f"that flatters the fit."
        )
    return values


def _element_position(label: str, element: str, size: int) -> int:
    """Which element of a block a column label refers to.

    The label is the coordinate value where the dimension had one, so an
    integer index is tried first and a positional fallback is not attempted:
    a plate coordinate of dataset labels has no numeric reading, and guessing
    one would attach a truth to the wrong object.
    """
    try:
        index = int(element)
    except ValueError:
        raise ResultsError(
            f"column {label!r} is element {element!r} of an array-valued block whose dimension "
            f"carries named coordinates, so a truth for it cannot be positioned by name. Pass "
            f"truths as a sequence in column order instead."
        ) from None
    if not 0 <= index < size:
        raise ResultsError(
            f"the truth given for column {label!r} has {size} element(s), which does not reach "
            f"index {index}."
        )
    return index


def plot_trace(
    tree: Any,
    *,
    var_names: Sequence[str] | None = None,
    group: str = "posterior",
    combined: bool = False,
    max_variables: int = MAX_TRACE_VARIABLES,
    paginate: bool = True,
    **kwargs: Any,
) -> Any:
    """Per-chain traces and marginals, the convergence eyeball.

    Reads ``sample_stats`` alongside the posterior, so a run whose draws include
    prior-rejected points shows them: those have ``lp = -inf`` and a **NaN**
    ``log_likelihood``, and rendering NaN as a gap rather than as zero is the
    visible half of ``inference.md`` §18(c)'s distinction between "not
    evaluated" and "impossible". A zero would be a value the parameter took;
    a gap is the truth, which is that nothing was evaluated there.

    Each variable gets a row of two panels — the marginal on the left, the
    trace against draw index on the right — and ``lp`` from ``sample_stats``
    gets a row of its own at the bottom of **every page**, because a trace of
    the parameters without the log-density beside it hides the commonest
    failure a trace plot exists to catch, whichever page is being read.

    **Paging (W3.10).** More than ``max_variables`` parameter rows used to be
    refused outright; now they are **paged** instead — split into consecutive
    figures each within the cap, in merged-name order, with an array-valued
    block kept whole on one page where it fits on one at all. A block bigger
    than ``max_variables`` on its own cannot fit any page whole, so it alone
    is split across full pages of exactly ``max_variables`` rows, in element
    order. Paging fires a :class:`ResultsWarning` naming the page count, the
    cap and the ``var_names=`` route to a smaller figure instead, and every
    page's :func:`~ampere.results.figure_metadata` carries ``"page"`` as
    ``"i of n"``. ``paginate=False`` restores the original behaviour exactly:
    one figure or a :class:`~ampere.core.exceptions.ResultsError` naming the
    count and ``max_variables``, for a caller who needs a single figure and
    would rather fail than receive a list.

    Parameters
    ----------
    tree
        The run, and nothing else.
    var_names
        Merged parameter names to draw; all of the group's by default.
    group
        Which group to read.
    combined
        Pool the chains into one trace and one marginal rather than drawing
        each chain separately. The default is per chain, because chains that
        disagree are exactly what this plot is looked at for.
    max_variables
        The per-page cap above. Raising it is a deliberate act.
    paginate
        Page above the cap (the default) rather than refuse. ``False``
        restores the pre-W3.10 refusal, unchanged, for a caller who needs
        exactly one figure or a hard failure.
    **kwargs
        Forwarded to the trace line artist.

    Returns
    -------
    matplotlib.figure.Figure | list[matplotlib.figure.Figure]
        A single :class:`~matplotlib.figure.Figure` — exactly today's return
        — for every call whose rows fit within ``max_variables``, paginated
        or not. Only a call that actually pages returns a :class:`list` of
        figures, in page order. Each figure carries
        :func:`~ampere.results.figure_metadata` with the number of
        prior-rejected draws, so the gaps are countable as well as visible; a
        paged figure's metadata also carries ``"page"``.
    """
    dataset = _p.require_sampling_group(tree, group)
    names = _p.select_names(
        [str(name) for name in dataset.data_vars], var_names, what=f"{group} variable"
    )
    keep = _p.scored_draws(tree)
    stats = tree["sample_stats"].dataset if keep is not None else None
    if not paginate:
        columns = _p.parameter_columns(dataset, names, limit=max_variables, what="a trace plot")
        return _render_trace_page(
            columns, keep=keep, stats=stats, combined=combined, page=None, **kwargs
        )
    columns = _p.parameter_columns(dataset, names, limit=_UNPAGED, what="a trace plot")
    pages = _paginate_columns(columns, limit=max_variables)
    if len(pages) == 1:
        return _render_trace_page(
            pages[0], keep=keep, stats=stats, combined=combined, page=None, **kwargs
        )
    _warn_paged(what="a trace plot", total=len(columns), limit=max_variables, n_pages=len(pages))
    return [
        _render_trace_page(
            page_columns,
            keep=keep,
            stats=stats,
            combined=combined,
            page=(index, len(pages)),
            **kwargs,
        )
        for index, page_columns in enumerate(pages, start=1)
    ]


def _render_trace_page(
    columns: Sequence[tuple[str, np.ndarray]],
    *,
    keep: np.ndarray | None,
    stats: Any,
    combined: bool,
    page: tuple[int, int] | None,
    **kwargs: Any,
) -> Any:
    """One figure's worth of :func:`plot_trace`, unpaged or one page of many.

    ``lp`` is appended to *every* page's rows, not only the last: the trace
    plot's own reason for drawing it beside the parameters — catching the
    commonest convergence failure — applies to whichever page is on screen,
    and it costs one row, not the whole cap.
    """
    rows: list[tuple[str, np.ndarray]] = list(columns)
    if stats is not None and "lp" in stats.variables:
        rows.append(("lp", np.asarray(stats["lp"].values, dtype=float)))
    rejected = 0 if keep is None else int(np.count_nonzero(~keep))
    figure, axes = _p.grid_axes(len(rows), 2)
    for row, (label, values) in enumerate(rows):
        # A prior-rejected draw is a gap, never a zero and never a level line
        # joining the points either side of it: NaN is what matplotlib breaks
        # a line at, and it is also what keeps the draw out of the marginal.
        if keep is None:
            gapped = np.where(np.isfinite(values), values, np.nan)
        else:
            gapped = np.where(keep, values, np.nan)
        _draw_marginal(axes[row, 0], gapped, combined=combined)
        _draw_trace(axes[row, 1], gapped, combined=combined, **kwargs)
        axes[row, 0].set_ylabel(label)
        axes[row, 0].set_xlabel(label)
        axes[row, 1].set_xlabel("draw")
        if row == 0:
            axes[row, 0].set_title("marginal", fontsize="small")
            axes[row, 1].set_title(
                "trace (gaps are prior-rejected draws)" if rejected else "trace", fontsize="small"
            )
    _p.attach_metadata(figure, "trace.rejected_draws", str(rejected))
    if page is not None:
        index, total_pages = page
        _p.attach_metadata(figure, "page", f"{index} of {total_pages}")
    figure.tight_layout()
    return figure


def _draw_marginal(panel: Any, values: np.ndarray, *, combined: bool) -> None:
    """A histogram per chain, or one over the pool.

    A histogram rather than a kernel density: a KDE of a handful of draws
    invents a shape, and a run small enough for that to matter is exactly the
    run a user is squinting at to decide whether to run more.
    """
    series = [values.ravel()] if combined else [values[chain] for chain in range(values.shape[0])]
    for index, chain in enumerate(series):
        finite = chain[np.isfinite(chain)]
        if finite.size == 0:
            continue
        panel.hist(
            finite,
            bins=min(30, max(5, finite.size // 4)),
            histtype="step",
            density=True,
            label=None if combined else f"chain {index}",
        )
    if not combined and values.shape[0] > 1:
        panel.legend(loc="best", fontsize="xx-small")
    panel.set_ylabel("density")


def _draw_trace(panel: Any, values: np.ndarray, *, combined: bool, **kwargs: Any) -> None:
    """The trace itself: one line per chain, NaN left as a break."""
    style = {"linewidth": 0.8, **kwargs}
    if combined:
        panel.plot(np.arange(values.size), values.ravel(), **style)
        return
    draws = np.arange(values.shape[1])
    for chain in range(values.shape[0]):
        panel.plot(draws, values[chain], **style)


# ---------------------------------------------------------------------------
# Family B — post-fit residual whiteness and posterior-predictive checks
# ---------------------------------------------------------------------------


def plot_posterior_predictive(
    tree: Any,
    *,
    datasets: Sequence[str] | None = None,
    statistic: Any = None,
    band: float = 0.68,
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

    Each dataset gets two panels: the observations with their uncertainties
    against the replicate median and its credible band, and the replicates'
    distribution of ``T`` with the observed value marked and the p-value
    reported.

    What the default statistic is, precisely
    ----------------------------------------
    ``T(y) = sum_i ((y_i - m_i) / sigma_i)^2``, with ``m_i`` the replicate mean
    at sample *i* and ``sigma_i`` the observed container's own uncertainty (or
    1 where there is none, which the axis label says). It is a function of
    ``y`` alone rather than of ``(y, θ)``, because the group stores replicate
    observations and not the per-draw prediction that produced them — the
    posterior-predictive group's shape is one variable per dataset
    (``results.md`` §7), and it is not this plot's business to widen it. A
    caller who wants a genuinely θ-dependent discrepancy passes one:
    ``statistic(values, sigma) -> float``, applied to the observations and to
    every replicate alike.

    Parameters
    ----------
    tree
        The run, with the ``posterior_predictive`` group attached.
    datasets
        Which datasets to draw; all of the group's by default.
    statistic
        ``statistic(values, sigma) -> float``. ``sigma`` is ``None`` where the
        dataset carries no uncertainties.
    band
        Credible-interval mass of the replicate band, in ``(0, 1)``.
    **kwargs
        Forwarded to the replicate-median line artist.

    Returns
    -------
    matplotlib.figure.Figure
        With :func:`~ampere.results.figure_metadata` carrying each dataset's
        p-value and observed statistic, so the number on the figure is
        reachable without reading it off the title.
    """
    if not 0.0 < band < 1.0:
        raise ResultsError(f"band is a credible-interval mass in (0, 1), got {band!r}.")
    group = _p.require_group(
        tree,
        POSTERIOR_PREDICTIVE_GROUP,
        remedy="call ampere.results.add_posterior_predictive(tree, problem) first, which draws "
        "y_rep through FittingProblem.simulate(observe=True) at the stored draws "
        "(results.md §7).",
    )
    observed_group = _p.require_group(
        tree,
        "observed_data",
        remedy="emit the run with observed=True (the default), which stores each dataset's "
        "observations; a predictive check has nothing to check against without them.",
    )
    available = [str(name) for name in group.data_vars]
    labels = _p.select_datasets(available, datasets, what=POSTERIOR_PREDICTIVE_GROUP)
    measure = _standardised_sum_of_squares if statistic is None else statistic
    figure, axes = _p.grid_axes(len(labels), 2, size=(6.0, 2.6))
    lower, upper = 50.0 * (1.0 - band), 50.0 * (1.0 + band)
    for index, label in enumerate(labels):
        if label not in observed_group.variables:
            raise ResultsError(
                f"this run holds replicates for dataset {label!r} but no observations of it, so "
                f"there is nothing to compare them against."
            )
        axis_name, coordinates = _p.coordinate_of(group, label)
        coordinate_unit, value_unit = _p.dataset_units(tree, label, axis_name)
        replicates = np.asarray(group[label].values, dtype=float)
        flat = replicates.reshape(-1, replicates.shape[-1])
        observations = np.asarray(observed_group[label].values, dtype=float).ravel()
        sigma = _stored_uncertainty(tree, label)
        _draw_replicates(
            axes[index, 0],
            coordinates,
            observations,
            sigma,
            flat,
            band=(lower, upper, band),
            labels=(_p.axis_label(axis_name, coordinate_unit), _p.axis_label(label, value_unit)),
            **kwargs,
        )
        check = _predictive_pvalue(measure, observations, sigma, flat, label)
        _draw_discrepancy(axes[index, 1], check, label)
        _p.attach_metadata(figure, f"{label}.posterior_predictive_p_value", f"{check[0]:.6g}")
        _p.attach_metadata(figure, f"{label}.observed_statistic", f"{check[1]:.6g}")
    figure.tight_layout()
    return figure


def _stored_uncertainty(tree: Any, label: str) -> np.ndarray | None:
    """One dataset's uncertainties from ``constant_data``, or ``None``.

    Read off the run rather than asked of a problem, because ``results.md`` §8
    is that every function here takes the emitted run and nothing else. Masked
    samples are left alone: the replicate group already carries NaN there, and
    NaN times anything stays NaN.
    """
    children = getattr(tree, "children", {})
    if "constant_data" not in children:
        return None
    constant = tree["constant_data"].dataset
    name = f"{label}_uncertainty"
    if name not in constant.variables:
        return None
    return np.asarray(constant[name].values, dtype=float).ravel()


def _standardised_sum_of_squares(values: np.ndarray, sigma: np.ndarray | None) -> float:
    """``sum_i (y_i / sigma_i)^2`` over the samples that are not gaps.

    Centring is the caller's — :func:`_predictive_pvalue` subtracts the
    replicate mean before calling, so that this is the same function of the
    observations and of every replicate. NaN samples (masked, or a draw the
    forward model could not complete) take no part, which is the same rule the
    whiteness statistic follows.
    """
    finite = np.isfinite(values)
    if sigma is not None:
        scaled = np.where(finite, values, 0.0) / np.where(np.isfinite(sigma), sigma, np.inf)
    else:
        scaled = np.where(finite, values, 0.0)
    return float(np.sum(scaled[finite] ** 2))


def _predictive_pvalue(
    measure: Callable[[np.ndarray, np.ndarray | None], float],
    observations: np.ndarray,
    sigma: np.ndarray | None,
    replicates: np.ndarray,
    label: str,
) -> tuple[float, float, np.ndarray]:
    """``(p, T(y), T(y_rep))`` — the Bayesian p-value and what it came from.

    ``p = mean_k [T(y_rep_k) >= T(y)]``: the probability that a replicate is at
    least as discrepant as the data. Near 0 means the fit is worse than the
    model can explain; near 1 means it is *better*, which usually means the
    uncertainties are overstated — the same reading
    :class:`~ampere.results.diagnostics.ChiSquareCheck` documents, and
    deliberately the same convention, so the cheap check and the replicate one
    cannot be read in opposite directions.

    The standard caveat applies and is not a defect of this implementation: the
    replicates are drawn from a posterior conditioned on the same data the
    statistic is evaluated at, so a posterior-predictive p-value is
    *conservative* — it is pulled towards 0.5 relative to a p-value from
    held-out data, and a value near 0.5 is therefore weaker evidence of a good
    fit than it looks. It is a screening device, not a test.
    """
    usable = np.array([row for row in replicates if np.any(np.isfinite(row))])
    if usable.size == 0:
        raise ResultsError(
            f"every replicate of dataset {label!r} is NaN, so there is no predictive "
            f"distribution to compare against. That happens when no stored draw could be "
            f"simulated — check the run's sample_stats for failures."
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        centre = np.nanmean(usable, axis=0)
    observed = float(measure(observations - centre, sigma))
    drawn = np.array([float(measure(row - centre, sigma)) for row in usable])
    return float(np.mean(drawn >= observed)), observed, drawn


def _draw_replicates(
    panel: Any,
    coordinates: np.ndarray,
    observations: np.ndarray,
    sigma: np.ndarray | None,
    replicates: np.ndarray,
    *,
    band: tuple[float, float, float],
    labels: tuple[str, str],
    **kwargs: Any,
) -> None:
    """Observations with their errors, against the replicate median and band."""
    lower, upper, mass = band
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        centre = np.nanmedian(replicates, axis=0)
        low = np.nanpercentile(replicates, lower, axis=0)
        high = np.nanpercentile(replicates, upper, axis=0)
    panel.fill_between(coordinates, low, high, alpha=0.3, label=f"{mass:.0%} of replicates")
    panel.plot(coordinates, centre, linewidth=1.2, label="replicate median", **kwargs)
    panel.errorbar(
        coordinates,
        observations,
        yerr=None if sigma is None else sigma,
        fmt=".",
        color="k",
        markersize=4,
        linewidth=0.8,
        label="observed",
    )
    panel.set_xlabel(labels[0])
    panel.set_ylabel(labels[1])
    panel.legend(loc="best", fontsize="small")


def _draw_discrepancy(panel: Any, check: tuple[float, float, np.ndarray], label: str) -> None:
    """The replicates' distribution of ``T`` with the observed value marked.

    The p-value's own resolution is the draw count, and it is on the title for
    the same reason the permutation floor is on the whiteness panel: a p of
    ``0`` computed from eight replicates is not evidence of anything, and a
    reader who cannot see the count cannot tell that from a real rejection.
    """
    p_value, observed, drawn = check
    panel.hist(drawn, bins=min(25, max(5, drawn.size // 3)), histtype="step", density=True)
    panel.axvline(observed, color="k", linewidth=1.2, label="observed T")
    verdict = (
        "replicates disagree with the data" if p_value < 0.05 or p_value > 0.95 else "consistent"
    )
    panel.set_title(
        f"{label}: T = {observed:.4g}, p = {p_value:.3g} — {verdict} "
        f"({drawn.size} replicate(s); resolution 1/{drawn.size})",
        fontsize="small",
    )
    panel.set_xlabel("discrepancy T")
    panel.set_ylabel("density")
    panel.legend(loc="best", fontsize="small")


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


# ---------------------------------------------------------------------------
# Family D — posterior calibration
# ---------------------------------------------------------------------------


#: How many ranked columns :func:`plot_sbc_ranks` will draw before refusing.
#: The same judgement :data:`MAX_TRACE_VARIABLES` records, for the same reason:
#: a rank histogram per column of a 10⁵-element latent block is a mistake, and
#: refusing it with the count named is the difference between a refusal and a
#: hung process.
MAX_RANK_PANELS = 40


def _calibration_of(result: Any) -> Any:
    """A calibration group, from either the group itself or the run holding it.

    Both are worth accepting. A caller who has just run
    :func:`~ampere.results.calibration.sbc` holds the dataset; a caller reading
    an archived run holds the tree, and having to reach into it by group name
    would be exactly the "reporting a missing group" ``results.md`` §8 says a
    plot must not do.
    """
    if hasattr(result, "data_vars") and "ranks" in getattr(result, "data_vars", {}):
        return result
    children = getattr(result, "children", None)
    if children is not None and CALIBRATION_GROUP in children:
        return result[CALIBRATION_GROUP].dataset
    raise ResultsError(
        f"this plot draws a {CALIBRATION_GROUP!r} group and was given a "
        f"{type(result).__name__} that is neither one nor a run carrying one. Compute it with "
        f"ampere.results.sbc(problem, engine_factory, ...) or, for an SBIEngine's own trained "
        f"posterior, engine.calibrate(count=...), and attach it with "
        f"ampere.results.attach_calibration(run, calibration)."
    )


def _selected_columns(
    group: Any, parameters: Sequence[str] | None, limit: int, *, what: str, override: str
) -> list[int]:
    """Which ranked columns to draw, by position, refusing an unreadable figure."""
    available = [str(name) for name in np.asarray(group.coords[PARAMETER_DIM].values).ravel()]
    chosen = _p.select_names(available, parameters, what="calibrated parameter")
    if len(chosen) > limit:
        raise ResultsError(
            f"{len(chosen)} {what} were asked for and the limit is {limit}. A figure that "
            f"crowded is not readable and is usually a selection mistake; narrow it with "
            f"parameters={override}."
        )
    return [available.index(name) for name in chosen]


def plot_sbc_ranks(
    result: Any,
    *,
    parameters: Sequence[str] | None = None,
    bins: int | None = None,
    max_panels: int = MAX_RANK_PANELS,
    **kwargs: Any,
) -> Any:
    """The SBC rank histogram, one panel per calibrated parameter.

    ``diagnostics.md`` §11's first output. Under a calibrated posterior the
    rank of the true value among ``L`` posterior draws is uniform on
    ``{0, ..., L}``, so a flat histogram is the pass and the *shape* of a
    failure says what kind it is: **U-shaped** means the posteriors are too
    narrow (the truth keeps landing in the tails), **inverted-U-shaped** means they are
    too wide, and a **slope** means they are biased. That reading is drawn on
    each panel rather than left to the caption, because a histogram without it
    is a picture a reader has to already know how to interpret.

    The grey band is the 99 % interval of the binomial null — the spread a
    *correctly* calibrated posterior produces at this many simulations — so a
    small study's noise is visibly noise rather than a finding. The panel title
    carries the Kolmogorov—Smirnov p-value the group stores.

    Parameters
    ----------
    result
        A ``calibration`` group (:func:`~ampere.results.calibration.sbc` or
        :meth:`ampere.inference.SBIEngine.calibrate`) or a run carrying one.
    parameters
        Which calibrated parameters to draw; all of them by default.
    bins
        Histogram bins. The default is at most 20, and never more than there
        are simulations divided by 5 — a rank histogram with fewer than a
        handful of counts per bin shows sampling noise and nothing else.
    max_panels
        The refusal threshold; see :data:`MAX_RANK_PANELS`.
    **kwargs
        Forwarded to the bar artist.

    Returns
    -------
    matplotlib.figure.Figure
        With :func:`~ampere.results.figure_metadata` carrying each parameter's
        KS p-value, so a caller reading numbers off the figure does not have to
        open the group to get them.
    """
    group = _calibration_of(result)
    positions = _selected_columns(
        group, parameters, max_panels, what="rank histograms", override=", or raise max_panels="
    )
    labels = [str(name) for name in np.asarray(group.coords[PARAMETER_DIM].values).ravel()]
    ranks = np.asarray(group["ranks"].values, dtype=float)
    draws = int(group.attrs.get(f"{ATTR_PREFIX}calibration_posterior_draws", ranks.max() or 1))
    simulations = ranks.shape[0]
    count = _rank_bins(bins, simulations)
    pvalues = np.asarray(group["ks_pvalue"].values, dtype=float) if "ks_pvalue" in group else None

    figure, axes = _p.new_axes(nrows=len(positions), size=(8.0, 2.6))
    edges = np.linspace(0.0, float(draws), count + 1)
    expected = simulations / count
    # The 99 % binomial interval of a single bin's count under uniformity. Drawn
    # rather than described: at 20 simulations the null is wide enough that an
    # eye-catching histogram is usually nothing at all.
    low, high = st.binom(simulations, 1.0 / count).ppf([0.005, 0.995])
    for panel_index, position in enumerate(positions):
        panel = axes[panel_index]
        panel.axhspan(
            float(low), float(high), color="0.85", zorder=0, label="99 % of a uniform null"
        )
        panel.axhline(expected, color="0.5", linewidth=0.8, zorder=1)
        panel.hist(ranks[:, position], bins=edges, zorder=2, **kwargs)
        panel.set_xlabel(f"rank of the truth among {draws} posterior draws")
        panel.set_ylabel("simulations")
        verdict = _rank_verdict(ranks[:, position], draws)
        title = f"{labels[position]}: {verdict}"
        if pvalues is not None:
            title += f" (KS p = {float(pvalues[position]):.3g}, {simulations} simulation(s))"
            _p.attach_metadata(
                figure, f"{labels[position]}.ks_pvalue", f"{float(pvalues[position]):.6g}"
            )
        panel.set_title(title, fontsize="small")
        panel.legend(loc="best", fontsize="xx-small")
    figure.tight_layout()
    return figure


def _rank_bins(bins: int | None, simulations: int) -> int:
    """How many bins a rank histogram gets, capped by the number of simulations."""
    if bins is not None:
        chosen = int(bins)
        if chosen < 1:
            raise ResultsError(f"a rank histogram needs at least one bin, got {bins}.")
        return chosen
    return max(1, min(20, simulations // 5))


def _rank_verdict(ranks: np.ndarray, draws: int) -> str:
    """The shape of the histogram, named: U, inverted-U, sloped, or flat.

    Read off two moments of the rank fraction rather than off the drawn bars,
    so the words and the picture are two views of one number: a variance above
    the uniform's says the truth keeps landing in the tails (too narrow), one
    below says it keeps landing in the middle (too wide), and a mean away from
    a half says the posterior is displaced.
    """
    fraction = np.asarray(ranks, dtype=float) / max(float(draws), 1.0)
    if fraction.size == 0:
        return "no simulations"
    mean = float(np.mean(fraction))
    spread = float(np.var(fraction))
    uniform_spread = 1.0 / 12.0
    error = math.sqrt(uniform_spread / max(fraction.size, 1))
    if abs(mean - 0.5) > 3.0 * error:
        return "sloped — the posterior is biased " + ("low" if mean > 0.5 else "high")
    if spread > uniform_spread * 1.4:
        return "U-shaped — the posteriors are too narrow"
    if spread < uniform_spread * 0.6:
        return "inverted-U-shaped — the posteriors are too wide"
    return "consistent with uniform"


def plot_coverage(
    result: Any,
    *,
    parameters: Sequence[str] | None = None,
    tarp: bool = True,
    ax: Any = None,
    **kwargs: Any,
) -> Any:
    """The coverage curve: empirical against nominal credible level.

    ``diagnostics.md`` §11's second output, and the one a reader who does not
    think in ranks can act on directly: at nominal level ``alpha``, what fraction
    of the simulated truths actually fell inside the posterior's central ``alpha``
    credible interval? A calibrated posterior puts the curve on the diagonal;
    a curve **below** it is overconfident (the intervals are too small for
    their label) and one **above** it is conservative.

    One line per parameter, from the marginal ranks. Where the group carries
    TARP's own curve (:meth:`ampere.inference.SBIEngine.calibrate` with
    ``tarp=True``), it is drawn beside them in black, and the distinction is
    stated in the legend rather than left implicit: TARP measures the *joint*
    coverage, and a posterior can be perfectly calibrated in every
    one-dimensional margin while being wrong about the correlations between
    them (Lemos et al. 2023).

    Parameters
    ----------
    result
        A ``calibration`` group or a run carrying one.
    parameters
        Which calibrated parameters to draw; all of them by default.
    tarp
        Draw TARP's joint curve when the group has one.
    ax
        Draw on this axes rather than a new figure.
    **kwargs
        Forwarded to the line artist.

    Returns
    -------
    matplotlib.axes.Axes
        The axes drawn on, with
        :func:`~ampere.results.figure_metadata` carrying each parameter's
        coverage at the 68 % and 95 % levels — the two a reader quotes — and
        TARP's area-to-curve where there is one.
    """
    group = _calibration_of(result)
    positions = _selected_columns(
        group, parameters, MAX_RANK_PANELS, what="coverage curves", override=""
    )
    labels = [str(name) for name in np.asarray(group.coords[PARAMETER_DIM].values).ravel()]
    levels = np.asarray(group.coords[LEVEL_DIM].values, dtype=float)
    coverage = np.asarray(group["coverage"].values, dtype=float)
    simulations = int(
        group.attrs.get(f"{ATTR_PREFIX}calibration_simulations", group.sizes[SIMULATION_DIM])
    )

    figure, axes = _p.new_axes(ax, size=(5.5, 5.0))
    panel = axes[0]
    panel.plot(
        [0.0, 1.0], [0.0, 1.0], color="0.5", linestyle="--", linewidth=0.9, label="calibrated"
    )
    for position in positions:
        panel.plot(levels, coverage[:, position], marker=".", label=labels[position], **kwargs)
        for quoted in (0.68, 0.95):
            value = float(np.interp(quoted, levels, coverage[:, position]))
            _p.attach_metadata(
                figure, f"{labels[position]}.coverage_{round(quoted * 100)}", f"{value:.6g}"
            )
    if tarp and "tarp_coverage" in group:
        panel.plot(
            np.asarray(group.coords[TARP_LEVEL_DIM].values, dtype=float),
            np.asarray(group["tarp_coverage"].values, dtype=float),
            color="black",
            linewidth=1.4,
            label="TARP (joint)",
        )
        area = group.attrs.get(f"{ATTR_PREFIX}calibration_tarp_atc")
        if area is not None:
            _p.attach_metadata(figure, "tarp.area_to_curve", f"{float(area):.6g}")
    panel.set_xlabel("nominal credible level")
    panel.set_ylabel("empirical coverage")
    panel.set_xlim(0.0, 1.0)
    panel.set_ylim(0.0, 1.0)
    panel.set_aspect("equal", adjustable="box")
    panel.set_title(
        f"coverage over {simulations} simulation(s) — below the diagonal is overconfident",
        fontsize="small",
    )
    panel.legend(loc="best", fontsize="x-small")
    return panel
