"""Rendering and run-reading helpers shared by the diagnostic plots.

Private on purpose, and separate from :mod:`ampere.results.plots` on purpose.
``plots.py`` is the frozen surface ``results.md`` §8 fixes, and **two** work
items write bodies into it — W2.7's three diagnostic renderers
(:func:`~ampere.results.plots.plot_residuals`,
:func:`~ampere.results.plots.plot_gp_localisation`,
:func:`~ampere.results.plots.plot_anomaly_score`) and W2.8's three
general-purpose ones. Anything both would otherwise define at the top of that
one module lives here instead, so that neither item has to move the other's
code to land.

Nothing here is part of the contract. The two pieces that *are* contract, and
that this module only supplies the mechanism for, are:

* the GP-localisation caveat must survive ``show_caveat=False`` on the returned
  figure (``results.md`` §8), which :func:`attach_metadata` /
  :func:`figure_metadata` are how; and
* two anomaly scores of different provenance drawn on one axes may never have
  their provenance suppressed (``diagnostics.md`` §5), which
  :func:`register_score` is how — it is the axes, not the call, that knows a
  second score has arrived.
"""

from __future__ import annotations

import json
import weakref
from collections.abc import Callable, Sequence
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core.exceptions import OptionalDependencyError, ResultsError
from ampere.core.results_schema import Axis

from .provenance import ATTR_PREFIX

__all__ = [
    "attach_metadata",
    "axis_label",
    "coordinate_of",
    "dataset_units",
    "distinct_provenances",
    "figure_metadata",
    "gp_datasets",
    "grid_axes",
    "likelihood_specs",
    "new_axes",
    "parameter_columns",
    "register_score",
    "require_corner",
    "require_group",
    "require_matplotlib",
    "require_sampling_group",
    "run_seed",
    "scored_draws",
    "select_datasets",
    "select_names",
    "wrap",
]

#: Where :func:`attach_metadata` hangs its dictionary on a figure.
_METADATA_ATTRIBUTE = "_ampere_metadata"

#: Per-axes record of which provenances have been drawn on it, and by which
#: artists, so that a second score of a *different* provenance can force the
#: labels back on — including on the first score, which was drawn before
#: anybody knew a second one was coming. Weak, so holding a figure open is the
#: caller's business and closing one frees this.
_SCORES_ON_AXES: weakref.WeakKeyDictionary[Any, list[tuple[str, list[Any]]]] = (
    weakref.WeakKeyDictionary()
)


def require_matplotlib() -> Any:
    """Import matplotlib's pyplot on use, never on import.

    matplotlib is a **base** dependency of ampere, so ``extra=None``: an
    environment without it is incomplete rather than missing an extra. The
    import stays lazy for the same two reasons arviz's does
    (:mod:`ampere.results`): it is expensive, and it is on the path of anyone
    who merely wanted :func:`~ampere.results.provenance.hash_container`.
    """
    try:
        import matplotlib.pyplot as pyplot
    except ImportError as error:  # pragma: no cover - exercised by a minimal install
        raise OptionalDependencyError(
            "matplotlib",
            context="drawing a diagnostic plot (ampere.results is where all plotting lives, "
            "DEVELOPMENT_PLAN.md §4.6). matplotlib is a base dependency of ampere, so this "
            "environment is incomplete rather than merely missing an extra",
        ) from error
    return pyplot


def require_corner() -> Any:
    """Import ``corner`` on use, never on import.

    ``corner`` is a **base** dependency (``pyproject.toml``), so ``extra=None``
    exactly as for matplotlib: an environment without it is incomplete rather
    than missing an extra. It is what draws the pairwise grid — the library
    every astronomer already reads — rather than a grid hand-rolled here, which
    is the same argument ``DEVELOPMENT_PLAN.md`` §4.6 makes for ArviZ being the
    single results format.
    """
    try:
        import corner
    except ImportError as error:  # pragma: no cover - exercised by a minimal install
        raise OptionalDependencyError(
            "corner",
            context="drawing the pairwise marginal grid (ampere.results.plot_corner). corner is "
            "a base dependency of ampere, so this environment is incomplete rather than merely "
            "missing an extra",
        ) from error
    return corner


def new_axes(ax: Any = None, *, nrows: int = 1, size: tuple[float, float] = (8.0, 3.0)) -> Any:
    """``(figure, axes_array)``, either fresh or wrapped around a caller's axes.

    A caller who passes ``ax`` gets it back as a one-element array and owns the
    figure; otherwise a figure of *nrows* stacked panels is created here.
    Passing ``ax`` with ``nrows > 1`` is refused rather than silently drawing
    two panels' worth of content into one.
    """
    pyplot = require_matplotlib()
    if ax is not None:
        if nrows != 1:
            raise ResultsError(
                f"this plot draws {nrows} panels, so it cannot be given a single ax=; let it "
                f"create its own figure, or ask for the panels you want individually."
            )
        return ax.get_figure(), np.array([ax], dtype=object)
    figure, axes = pyplot.subplots(nrows=nrows, figsize=(size[0], size[1] * nrows), squeeze=False)
    return figure, axes[:, 0]


def grid_axes(nrows: int, ncols: int, size: tuple[float, float] = (5.0, 2.2)) -> Any:
    """``(figure, axes)`` for a fresh ``nrows x ncols`` grid, never squeezed.

    :func:`new_axes` is the one-column form the diagnostic panels want;
    :func:`~ampere.results.plots.plot_trace` wants two columns (a marginal and
    a trace) per variable, and squeezing a 1-row grid into a bare axes is the
    kind of shape surprise that makes a loop over panels wrong only sometimes.
    """
    pyplot = require_matplotlib()
    figure, axes = pyplot.subplots(
        nrows=nrows, ncols=ncols, figsize=(size[0] * ncols, size[1] * nrows), squeeze=False
    )
    return figure, axes


def attach_metadata(figure: Any, key: str, value: str) -> None:
    """Record *value* on *figure* under *key*, whatever was drawn on it.

    ``results.md`` §8 requires ``show_caveat=False`` to suppress only the drawn
    annotation and "still leave the caveat on the returned figure's metadata" —
    a caption a user can crop out of a screenshot is not durable protection
    against over-interpretation (``diagnostics.md`` §4.3), and neither is one a
    keyword argument can delete.
    """
    metadata = getattr(figure, _METADATA_ATTRIBUTE, None)
    if metadata is None:
        metadata = {}
        setattr(figure, _METADATA_ATTRIBUTE, metadata)
    metadata[key] = value


def figure_metadata(figure: Any) -> dict[str, str]:
    """What :func:`attach_metadata` recorded on *figure*, as a plain dictionary.

    Empty for a figure ampere did not draw. This is the programmatic half of
    the caveat requirement: a user extracting numbers rather than reading a
    caption still gets the warning, from the same object the numbers came on.
    """
    return dict(getattr(figure, _METADATA_ATTRIBUTE, None) or {})


def require_group(tree: Any, group: str, *, remedy: str) -> Any:
    """The named group's dataset, or a refusal that names the remedy.

    ``results.md`` §8 is explicit that a missing derived group must be reported
    as "you have not computed this yet", naming the function that computes it,
    "rather than reporting a missing group" — the second is a fact about
    xarray and the first is an answer.
    """
    children = getattr(tree, "children", None)
    if children is None or group not in children:
        raise ResultsError(
            f"this run has no {group!r} group. It is not stored by default — the cost is "
            f"N_draws x N_obs per dataset (results.md §7) — so it is computed on demand: "
            f"{remedy}"
        )
    return tree[group].dataset


def require_sampling_group(tree: Any, group: str) -> Any:
    """A group every emitted run carries, or a refusal that says what is wrong.

    Distinct from :func:`require_group`, whose message explains that a
    *derived* group is not stored by default and names the function that
    computes it. ``posterior`` and ``sample_stats`` are not derived: a tree
    without them was not produced by :func:`~ampere.results.emit`, and telling
    a user to "compute it on demand" would send them the wrong way.
    """
    children = getattr(tree, "children", None)
    if children is None or group not in children:
        raise ResultsError(
            f"this run has no {group!r} group. Every run ampere.results.emit produces carries "
            f"it, so this is not an emitted run — the plotting surface takes the emitted run "
            f"and nothing else (results.md §8)."
        )
    return tree[group].dataset


def select_names(
    available: Sequence[str], requested: Sequence[str] | None, *, what: str
) -> tuple[str, ...]:
    """Which named variables to draw: all of *available*, or the caller's subset.

    The sibling of :func:`select_datasets` for things that are not datasets —
    posterior variables, whose names are the **merged** parameter names, which
    is what makes an unknown one worth reporting with the full list beside it:
    ``model.index`` versus ``index`` is the commonest way to get this wrong,
    and the merged name is exactly what ``inference.md`` §4.5 says a user
    types.
    """
    if requested is None:
        chosen = tuple(available)
    else:
        unknown = [name for name in requested if name not in available]
        if unknown:
            raise ResultsError(
                f"this run has no {what} named {sorted(unknown)}; it has {sorted(available)}. "
                f"Variables are keyed by merged parameter name (results.md §4)."
            )
        chosen = tuple(requested)
    if not chosen:
        raise ResultsError(f"no {what} was selected, so there is nothing to draw.")
    return chosen


def parameter_columns(
    group: Any, names: Sequence[str], *, limit: int, what: str
) -> list[tuple[str, np.ndarray]]:
    """``(label, values)`` per scalar column, expanding array-valued blocks.

    ``results.md`` §4: an array-valued parameter — a plate member, a latent GP
    block — is **one variable with a named dimension**, never ``N`` scalar
    names. A pairwise grid or a trace stack is nevertheless per scalar, so the
    block is expanded here, and its element labels come from the dimension's
    own coordinate when it has one (a plate's coordinate is the dataset labels,
    ``hierarchical_population.md`` §10.2) and from the integer index when it
    does not.

    *limit* is where the expansion is **refused rather than attempted**:
    ``likelihoods.md`` §16(a)'s 10⁵-element latent block is exactly the input a
    corner plot must decline, and declining it with the count and the offending
    variable named is the difference between a refusal and a hung process.
    ``values`` are ``(chain, draw)``.
    """
    columns: list[tuple[str, np.ndarray]] = []
    for name in names:
        data = group[name]
        extra = [dim for dim in data.dims if dim not in ("chain", "draw")]
        values = np.asarray(data.values)
        if not extra:
            columns.append((str(name), values.astype(float)))
            continue
        block = values.reshape(values.shape[0], values.shape[1], -1)
        size = block.shape[-1]
        if size > limit:
            raise ResultsError(
                f"parameter {name!r} is one array-valued block of {size} elements (dimension(s) "
                f"{extra}), and {what} of it would be {size} panels. results.md §8 requires this "
                f"to be refused loudly rather than attempted. Select the variables you want with "
                f"var_names=, or raise max_variables= if you really mean it."
            )
        labels = _element_labels(data, extra, size)
        for index in range(size):
            columns.append((f"{name}[{labels[index]}]", block[:, :, index].astype(float)))
    if len(columns) > limit:
        raise ResultsError(
            f"{what} of {len(columns)} variables was asked for, and the limit is {limit}. "
            f"A grid that large is not readable and is usually a selection mistake; narrow it "
            f"with var_names=, or raise max_variables= deliberately."
        )
    return columns


def _element_labels(data: Any, extra: Sequence[str], size: int) -> list[str]:
    """One label per element of an array-valued block, coordinates preferred."""
    if len(extra) == 1 and extra[0] in data.coords:
        return [str(value) for value in np.asarray(data.coords[extra[0]].values).ravel()]
    shape = tuple(int(data.sizes[dim]) for dim in extra)
    if len(shape) == 1:
        return [str(index) for index in range(size)]
    return [", ".join(str(part) for part in index) for index in np.ndindex(*shape)]


def scored_draws(tree: Any) -> np.ndarray | None:
    """``(chain, draw)`` boolean: which draws the prior did not reject.

    ``inference.md`` §18(c) and ``results.md`` §5: a prior-rejected draw stores
    ``lp = -inf`` and a **NaN** log-likelihood, and the two are different
    statements that must not be flattened. What that means for a plot is that
    such a draw is not a posterior sample: it belongs in a trace as a *gap*,
    and it does not belong in a marginal at all.

    ``None`` when the run carries no ``sample_stats`` — the honest answer for a
    tree that was not emitted by ampere, and the caller then treats every draw
    as scored.
    """
    children = getattr(tree, "children", None)
    if children is None or "sample_stats" not in children:
        return None
    stats = tree["sample_stats"].dataset
    if "lp" not in stats.variables:  # pragma: no cover - a hand-built run
        return None
    return np.isfinite(np.asarray(stats["lp"].values, dtype=float))


def likelihood_specs(tree: Any) -> dict[str, dict[str, Any]]:
    """Each dataset's likelihood declaration, read back off the run's own attrs.

    ``ampere_likelihoods`` is written by
    :func:`~ampere.results.provenance.provenance_attrs` from
    :meth:`ampere.core.likelihood.Likelihood.to_spec`, so a stored run knows
    what noise model it was fitted with without the problem being at hand.
    That is what lets :func:`~ampere.results.plots.plot_residuals` obey
    ``diagnostics.md`` §3.1 from the run alone.
    """
    raw = getattr(tree, "attrs", {}).get(f"{ATTR_PREFIX}likelihoods")
    if not isinstance(raw, str):
        return {}
    try:
        decoded = json.loads(raw)
    except json.JSONDecodeError:  # pragma: no cover - a corrupted run
        return {}
    if not isinstance(decoded, dict):  # pragma: no cover - a corrupted run
        return {}
    return {str(label): spec for label, spec in decoded.items() if isinstance(spec, dict)}


def gp_datasets(tree: Any) -> tuple[str, ...]:
    """Labels whose likelihood declared a :class:`GaussianProcessNoise` model.

    The test family B is scoped away from (``diagnostics.md`` §3.1) and family
    C is scoped to (§4.1) — one predicate, so the two cannot drift apart.
    """
    return tuple(
        label
        for label, spec in sorted(likelihood_specs(tree).items())
        if spec.get("noise") == "GaussianProcessNoise"
    )


def select_datasets(
    available: Sequence[str], datasets: Sequence[str] | None, *, what: str
) -> tuple[str, ...]:
    """Which labels to act on: all of *available*, or the caller's subset.

    An unknown label is an error rather than an empty panel — asking for a
    dataset a run does not hold is a typo far more often than it is a
    deliberate no-op.
    """
    if datasets is None:
        chosen = tuple(available)
    else:
        unknown = [label for label in datasets if label not in available]
        if unknown:
            raise ResultsError(
                f"this run has no {what} for dataset(s) {sorted(unknown)}; it has "
                f"{sorted(available)}."
            )
        chosen = tuple(datasets)
    if not chosen:
        raise ResultsError(f"this run holds no {what} to plot.")
    return chosen


#: **W5.3.** Attribute on a stored per-axis variable naming which of a
#: multi-axis kind's :attr:`~ampere.core.results_schema.FunctionSamples.AXES`
#: it is (``"u"``, ``"v"``, ...) — the marker that tells one of the kind's
#: own raw axes apart from anything else sharing its joint sample dimension
#: (a role-suffixed value variable such as ``<label>_variance``, in
#: particular). Written by ``ampere.results.derived`` at the point a group is
#: built, from the *live* container, where the real axis units are.
PLOT_AXIS_ATTR = f"{ATTR_PREFIX}plot_axis"

#: Suffix of the variable a kind's own
#: :attr:`~ampere.core.results_schema.FunctionSamples.PLOT_COORDINATE`
#: default is precomputed into. Precomputed at group-build time rather than
#: here, for the same reason: only there is the live container available.
PLOT_COORDINATE_VARIABLE_SUFFIX = "__plot_coordinate"

#: Group attribute prefix carrying that precomputed default's axis label
#: (``"baseline length [m]"``), keyed by dataset label.
PLOT_COORDINATE_LABEL_ATTR_PREFIX = f"{ATTR_PREFIX}plot_coordinate_label__"


def _stored_axes(dataset: Any, label: str, index_dim: str) -> dict[str, Axis]:
    """A multi-axis kind's own raw axes, as real :class:`Axis` objects.

    Found by :data:`PLOT_AXIS_ATTR` rather than by name pattern, so a
    role-suffixed value variable that happens to share the joint sample
    dimension (``<label>_variance``, for instance) is never mistaken for one
    of the kind's axes.
    """
    prefix = f"{label}_"
    found: dict[str, Axis] = {}
    for name, variable in dataset.variables.items():
        axis_name = variable.attrs.get(PLOT_AXIS_ATTR)
        if axis_name is None or tuple(variable.dims) != (index_dim,):
            continue
        if not str(name).startswith(prefix):
            continue
        unit_text = str(variable.attrs.get("units", ""))
        unit = u.Unit(unit_text) if unit_text else None
        found[str(axis_name)] = Axis.build(
            str(axis_name), np.asarray(variable.values, dtype=float), unit
        )
    return found


def coordinate_of(
    dataset: Any,
    variable: str,
    *,
    coordinate: str | Callable[[dict[str, Axis]], tuple[np.ndarray, str]] | None = None,
) -> tuple[str, np.ndarray]:
    """The coordinate axis a dataset variable is plotted against.

    A point kind with **one** axis carries it as ``<label>_<axis>`` and
    ``coordinate`` plays no part — that single axis is the answer, exactly as
    before this argument existed (**W5.3**), which is what keeps every
    existing plot byte-identical.

    A point kind with **several** axes (``results.md`` §4's ``<label>_index``
    joint sample dimension) has no coordinate of its own, and this is the one
    place the rule for finding one lives, in order:

    1. ``coordinate``, if given: an axis name (one of the kind's own,
       :class:`str`), or a callable taking a mapping of axis name to
       :class:`~ampere.core.results_schema.Axis` and returning
       ``(coordinates, label)``.
    2. The kind's own default, :attr:`~ampere.core.results_schema.
       FunctionSamples.PLOT_COORDINATE`, precomputed by ``ampere.results.derived``
       at the point the group was built and stored on it.
    3. Refused by name, listing the kind's own axes — never a silent 1-D
       projection of data with no natural order.
    """
    data = dataset[variable]
    spatial = [name for name in data.dims if name not in ("chain", "draw")]
    if len(spatial) != 1:
        raise ResultsError(
            f"dataset {variable!r} is indexed by {spatial}, and the 1-D diagnostics of "
            f"diagnostics.md families B and C need exactly one ordered coordinate axis. "
            f"Gridded kinds are out of scope here (DEVELOPMENT_PLAN.md §4.4)."
        )
    name = str(spatial[0])
    if not name.endswith("_index"):
        if coordinate is not None:
            raise ResultsError(
                f"dataset {variable!r} already has one ordered coordinate axis ({name!r}); "
                f"coordinate= is for a point kind with several axes (results.md §4), which "
                f"this one is not."
            )
        if name not in data.coords:
            raise ResultsError(
                f"dataset {variable!r} has no coordinate values on its {name!r} axis, so there "
                f"is nothing to measure separations against."
            )
        return name, np.asarray(data.coords[name].values, dtype=float)
    label = name[: -len("_index")]
    axes = _stored_axes(dataset, label, name)
    if coordinate is None:
        default_variable = f"{label}{PLOT_COORDINATE_VARIABLE_SUFFIX}"
        if default_variable in dataset.variables:
            label_text = dataset.attrs.get(
                f"{PLOT_COORDINATE_LABEL_ATTR_PREFIX}{label}", default_variable
            )
            return str(label_text), np.asarray(dataset[default_variable].values, dtype=float)
        raise ResultsError(
            f"dataset {variable!r} is a point kind with several axes {sorted(axes)} "
            f"(results.md §4) and this kind declares no default plotted coordinate. Pass "
            f"coordinate= naming one of {sorted(axes)}, or a callable(axes) -> "
            f"(coordinates, label) computing one."
        )
    if isinstance(coordinate, str):
        if coordinate not in axes:
            raise ResultsError(
                f"dataset {variable!r} has no axis named {coordinate!r}; its axes are "
                f"{sorted(axes)}."
            )
        axis = axes[coordinate]
        return axis_label(axis.name, "" if axis.unit is None else str(axis.unit)), np.asarray(
            axis.values, dtype=float
        )
    if callable(coordinate):
        coordinates, label_text = coordinate(axes)
        return str(label_text), np.asarray(coordinates, dtype=float)
    raise ResultsError(f"coordinate= must be a string or a callable, got {type(coordinate)!r}.")


def dataset_units(tree: Any, label: str, axis: str) -> tuple[str, str]:
    """``(coordinate unit, value unit)`` for one dataset, from the run's attrs.

    Empty strings where a container carried no unit, which is legitimate: the
    axis label then simply says the name.
    """
    coordinate_unit = ""
    value_unit = ""
    children = getattr(tree, "children", {})
    if "observed_data" in children:
        observed = tree["observed_data"].dataset
        if label in observed.variables:
            value_unit = str(observed[label].attrs.get("units", ""))
        if axis in observed.coords:
            coordinate_unit = str(observed.coords[axis].attrs.get("units", ""))
    return coordinate_unit, value_unit


def axis_label(name: str, unit: str) -> str:
    """``"wavelength [um]"``, or just the name where there is no unit."""
    return f"{name} [{unit}]" if unit else name


def run_seed(tree: Any) -> int | None:
    """The run's recorded seed, or ``None`` for an entropy-seeded run.

    A diagnostic that draws its own randomness (family B's permutation
    calibration) derives a named sub-stream from this, so re-running the check
    on a stored run reproduces the same p-value.
    """
    seed = getattr(tree, "attrs", {}).get(f"{ATTR_PREFIX}seed")
    return None if seed is None else int(seed)


def register_score(
    axes: Any, provenance: str, artists: Sequence[Any]
) -> list[tuple[str, list[Any]]]:
    """Record that *artists* on *axes* carry *provenance*; return everything on it.

    ``diagnostics.md`` §5 polices the comparability risk of a shared visual
    grammar with metadata rather than with visual distinctness, and
    ``results.md`` §8 turns that into a rule about this function's caller:
    ``show_provenance=False`` "may not suppress them when two scores of
    different provenance are drawn together". Two scores are drawn *together*
    when they share an axes, and only the axes knows that — the second call
    cannot see the first's arguments. So the record lives here, and the caller
    labels every artist on the axes, not only its own, as soon as the returned
    list holds more than one distinct provenance.
    """
    drawn = _SCORES_ON_AXES.setdefault(axes, [])
    drawn.append((provenance, list(artists)))
    return drawn


def distinct_provenances(drawn: Sequence[tuple[str, Sequence[Any]]]) -> tuple[str, ...]:
    """The distinct provenances in a :func:`register_score` record, in order."""
    seen: dict[str, None] = {}
    for provenance, _ in drawn:
        seen.setdefault(provenance, None)
    return tuple(seen)


def wrap(text: str, width: int = 96) -> str:
    """Hard-wrap a caption so a long caveat does not run off the figure."""
    words = text.split()
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        if current and sum(len(piece) + 1 for piece in current) + len(word) > width:
            lines.append(" ".join(current))
            current = [word]
        else:
            current.append(word)
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines)
