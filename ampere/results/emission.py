"""Turning a run into the one results format, and back off disk again.

``DEVELOPMENT_PLAN.md`` §4.6: "ArviZ ``InferenceData`` is the single results
format. Every engine emits it; corner/trace/posterior-predictive plotting is
written once against it in ``ampere.results``; serialisation is netCDF." This
module is the emitting half. :func:`emit` takes the two things an engine driver
already has — the draws, and one
:class:`~ampere.core.dataset.Evaluation` per draw — and produces the object;
:func:`to_netcdf` and :func:`from_netcdf` put it on disk and get it back.

What every run stores, and why each piece is not optional
---------------------------------------------------------
``DEVELOPMENT_PLAN.md`` §4.6 asks for per-sample ``log_likelihood`` and
``log_prior``; ``inference.md`` §18 sharpens that into three obligations, all
discharged here.

* The **scalar joint** ``log_likelihood`` and ``log_prior`` go in
  ``sample_stats`` (with ``lp`` for ``log_prob``, ArviZ's own name for the log
  posterior density). This is the reading design horizon (b) needs: population
  importance reweighting wants one number per draw per archived object.
* The **per-dataset decomposition** goes in the ``log_likelihood`` group, one
  variable per dataset label. ``likelihoods.md`` §16 records that "per-sample
  log_likelihood" has two legitimate readings and that the per-*observation* one
  is **not well defined for a GP likelihood**, whose samples do not factorise.
  The per-dataset one always is, so it is what ampere emits, and the group says
  so in its own ``ampere_decomposition`` attribute rather than leaving a
  consumer to assume ArviZ's usual per-observation convention.
  ``docs/design/contracts/results.md`` §6 names the per-observation
  decomposition ampere would compute if asked, and why it is not computed by
  default.
* ``log_likelihood`` is **NaN**, never ``-inf``, for a draw the prior rejected.
  ``inference.md`` §10 makes that distinction at the source and §18(c) asks this
  contract not to flatten it: "not evaluated" and "impossible" are different
  statements, and reweighting must be able to tell them apart.

A note on the object's identity
-------------------------------
ArviZ 1.0 replaced the ``InferenceData`` *class* with :class:`xarray.DataTree`.
The **format** is unchanged — the same named groups, the same netCDF layout, and
``arviz.from_netcdf`` still reads files written by ArviZ 0.x — so
``DEVELOPMENT_PLAN.md`` §4.6's decision stands untouched; only the Python type
moved. Everything here returns and accepts a ``DataTree``.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from ampere.core.dataset import Evaluation, FittingProblem
from ampere.core.exceptions import OptionalDependencyError
from ampere.core.parameter import Binding, Parameter
from ampere.core.results_schema import FunctionSamples, Layout

from .exceptions import ResultsError
from .provenance import ATTR_PREFIX, provenance_attrs

__all__ = [
    "CHAIN_DIM",
    "CONSTANT_DATA_GROUP",
    "DRAW_DIM",
    "LOG_LIKELIHOOD_DECOMPOSITION",
    "LOG_LIKELIHOOD_GROUP",
    "OBSERVED_DATA_GROUP",
    "POSTERIOR_GROUP",
    "SAMPLE_STATS_GROUP",
    "DrawRecorder",
    "emit",
    "from_netcdf",
    "to_netcdf",
]

#: ArviZ's two sampling dimensions, in its own order.
CHAIN_DIM = "chain"
DRAW_DIM = "draw"

POSTERIOR_GROUP = "posterior"
SAMPLE_STATS_GROUP = "sample_stats"
LOG_LIKELIHOOD_GROUP = "log_likelihood"
OBSERVED_DATA_GROUP = "observed_data"
CONSTANT_DATA_GROUP = "constant_data"

#: What the ``log_likelihood`` group's variables are a decomposition *into*.
#: Recorded on the group so no consumer has to assume ArviZ's per-observation
#: convention applies (``likelihoods.md`` §16; results.md §6).
LOG_LIKELIHOOD_DECOMPOSITION = "per_dataset"


def _require_arviz() -> Any:
    """Import arviz on use, never on import (``architecture.md`` §4, rules 2 and 3)."""
    try:
        import arviz
    except ImportError as error:  # pragma: no cover - exercised by the minimal-install job
        raise OptionalDependencyError(
            "arviz",
            extra="arviz",
            context="emitting a run's results (ampere.results uses ArviZ's DataTree as the single "
            "results format, DEVELOPMENT_PLAN.md §4.6)",
        ) from error
    return arviz


# ---------------------------------------------------------------------------
# Shaping the draws
# ---------------------------------------------------------------------------


def _as_draw_array(theta: Any, free_size: int) -> np.ndarray:
    array = np.asarray(theta, dtype=float)
    if array.ndim == 2:
        array = array[np.newaxis, ...]
    if array.ndim != 3:
        raise ResultsError(
            f"draws must be shaped (draw, free_size) or (chain, draw, free_size); got an array "
            f"with {array.ndim} dimension(s) and shape {array.shape}."
        )
    if array.shape[-1] != free_size:
        raise ResultsError(
            f"draws have {array.shape[-1]} column(s) but this problem has {free_size} free "
            f"dimension(s). An engine must emit the flat vector in FittingProblem.parameters "
            f"order — the order free_labels() names."
        )
    if array.shape[0] == 0 or array.shape[1] == 0:
        # An empty run is refused rather than emitted. It carries no
        # information, and what it would produce is an InferenceData with
        # zero-length sampling dimensions that every consumer downstream then
        # has to guard against -- the same reason DrawRecorder.emit refuses to
        # emit before anything has been recorded.
        raise ResultsError(
            f"a run needs at least one chain and one draw; the draws have shape {array.shape}."
        )
    return array


def _as_evaluation_grid(
    evaluations: Sequence[Any], chains: int, draws: int
) -> list[list[Evaluation]]:
    listed = list(evaluations)
    if listed and isinstance(listed[0], Evaluation):
        if chains != 1:
            raise ResultsError(
                f"a flat sequence of Evaluations was given for {chains} chains; nest it as one "
                f"sequence per chain."
            )
        grid = [listed]
    else:
        grid = []
        for chain in listed:
            if isinstance(chain, str) or not isinstance(chain, Sequence):
                raise ResultsError(
                    f"every recorded draw must be a FittingProblem.evaluate() result, got "
                    f"{type(chain).__name__}. Call problem.evaluate(theta) rather than "
                    f"log_prob(theta), which cannot carry the split or the failure."
                )
            grid.append(list(chain))
    if len(grid) != chains or any(len(chain) != draws for chain in grid):
        shape = [len(chain) for chain in grid]
        raise ResultsError(
            f"there are {chains} chain(s) of {draws} draw(s) but {len(grid)} sequence(s) of "
            f"{shape} evaluation(s). One Evaluation per draw, in draw order."
        )
    for chain in grid:
        for evaluation in chain:
            if not isinstance(evaluation, Evaluation):
                raise ResultsError(
                    f"every recorded draw must be a FittingProblem.evaluate() result, got "
                    f"{type(evaluation).__name__}. Call problem.evaluate(theta) rather than "
                    f"log_prob(theta), which cannot carry the split or the failure."
                )
    return grid


# ---------------------------------------------------------------------------
# Posterior group
# ---------------------------------------------------------------------------


def _dimension_names(parameter: Parameter) -> tuple[str, ...]:
    """Dimension names for an array-valued parameter.

    A plate's member takes the **plate's own name** as its dimension, so every
    member of one plate shares one dimension and one coordinate — which is what
    makes ``hierarchical_population.md`` §10.5's "``objects.theta`` should be one
    dimension with the dataset labels as its coordinate" expressible at all.
    Anything else array-valued (a latent GP block, a per-channel offset) gets
    ``<name>_dim_<i>``, ArviZ's own default shape of name.
    """
    if parameter.plate is not None and len(parameter.shape) == 1:
        return (parameter.plate,)
    return tuple(f"{parameter.name}_dim_{axis}" for axis in range(len(parameter.shape)))


def _index_coordinate(bindings: Sequence[Binding], parameter: Parameter) -> list[str] | None:
    """The component each element of a plate's array-valued parameter feeds.

    Read off the mapping's own ``Binding.index`` entries rather than guessed
    from the dataset ordering, so the coordinate states a fact the composition
    already recorded: element *i* of this block is consumed *there*.
    ``hierarchical_population.md`` §10.2 requires this — "W1.8 must record the
    dataset labels as the plate's coordinate, not an integer range" — and reading
    the bindings is how it stays true when the wiring is not in dataset order.

    Returns ``None`` when the block is not routed element by element, in which
    case ArviZ's integer range is the honest label.
    """
    size = parameter.shape[0] if parameter.shape else 0
    found: dict[int, str] = {}
    for binding in bindings:
        if binding.global_name != parameter.name or binding.index is None:
            continue
        index = binding.index
        if not isinstance(index, int):
            return None
        found[index] = binding.component
    if len(found) != size or set(found) != set(range(size)):
        return None
    return [found[index] for index in range(size)]


def _posterior_variables(
    problem: FittingProblem, draws: np.ndarray, coords: Mapping[str, Sequence[Any]] | None
) -> tuple[dict[str, np.ndarray], dict[str, list[str]], dict[str, Sequence[Any]]]:
    parameters = problem.parameters
    bindings = problem.mapping.bindings
    chains, count = draws.shape[0], draws.shape[1]
    variables: dict[str, np.ndarray] = {}
    dims: dict[str, list[str]] = {}
    resolved: dict[str, Sequence[Any]] = dict(coords or {})
    for parameter in parameters:
        if parameter.is_fixed:
            continue
        block = draws[:, :, parameters.free_slice(parameter.name)]
        if not parameter.shape:
            variables[parameter.name] = block[:, :, 0]
            continue
        # One variable with a named dimension, never `parameter.size` scalar
        # names: `likelihoods.md` §16 is explicit that 10^5 scalar labels is the
        # wrong ArviZ representation for a latent block.
        variables[parameter.name] = block.reshape(chains, count, *parameter.shape)
        names = _dimension_names(parameter)
        dims[parameter.name] = list(names)
        if names[0] not in resolved:
            labels = _index_coordinate(bindings, parameter)
            if labels is not None:
                resolved[names[0]] = labels
    return variables, dims, resolved


# ---------------------------------------------------------------------------
# sample_stats and log_likelihood groups
# ---------------------------------------------------------------------------


def _sample_stats(grid: Sequence[Sequence[Evaluation]]) -> dict[str, np.ndarray]:
    shape = (len(grid), len(grid[0]))
    log_prob = np.empty(shape, dtype=float)
    log_prior = np.empty(shape, dtype=float)
    log_likelihood = np.empty(shape, dtype=float)
    failed = np.zeros(shape, dtype=bool)
    reason = np.full(shape, "", dtype=object)
    where = np.full(shape, "", dtype=object)
    for c, chain in enumerate(grid):
        for d, evaluation in enumerate(chain):
            log_prob[c, d] = evaluation.log_prob
            log_prior[c, d] = evaluation.log_prior
            # Deliberately not coerced: NaN means the prior rejected the point
            # before the likelihood was ever called (`inference.md` §18c).
            log_likelihood[c, d] = evaluation.log_likelihood
            if evaluation.failure is not None:
                failed[c, d] = True
                reason[c, d] = str(evaluation.failure.reason)
                where[c, d] = evaluation.failure.where
    return {
        "lp": log_prob,
        "log_prior": log_prior,
        "log_likelihood": log_likelihood,
        "failed": failed,
        "failure_reason": reason.astype(str),
        "failure_where": where.astype(str),
    }


def _log_likelihood_variables(
    problem: FittingProblem, grid: Sequence[Sequence[Evaluation]]
) -> dict[str, np.ndarray]:
    shape = (len(grid), len(grid[0]))
    variables = {label: np.full(shape, np.nan) for label in problem.datasets}
    for c, chain in enumerate(grid):
        for d, evaluation in enumerate(chain):
            for label, value in evaluation.contributions.items():
                if label in variables:
                    variables[label][c, d] = value
    return variables


# ---------------------------------------------------------------------------
# observed_data and constant_data groups
# ---------------------------------------------------------------------------


def _container_dims(label: str, container: FunctionSamples) -> tuple[str, ...]:
    """Dimension names for one observed container.

    Scoped by dataset label, because two datasets may both have a
    ``spectral_axis`` of different length and one xarray group cannot hold two
    dimensions of the same name. A gridded kind takes one dimension per axis; a
    point kind with a single axis takes that axis (so a spectrum plots against
    wavelength without further work); a point kind with several axes — a
    visibility set's *u* and *v*, which index samples jointly rather than
    separately — takes one sample dimension, and the axes become ordinary
    variables on it.
    """
    if container.LAYOUT is Layout.GRID:
        return tuple(f"{label}_{axis.name}" for axis in container.axes)
    if len(container.axes) == 1:
        return (f"{label}_{container.axes[0].name}",)
    return (f"{label}_index",)


#: ``(groups, dims, coords, units)`` — everything ``arviz.from_dict`` needs for
#: the two data groups, plus the per-variable units applied afterwards.
_DataGroups = tuple[
    dict[str, dict[str, np.ndarray]],
    dict[str, list[str]],
    dict[str, Sequence[Any]],
    dict[str, str],
]


def _data_groups(problem: FittingProblem) -> _DataGroups:
    observed: dict[str, np.ndarray] = {}
    constant: dict[str, np.ndarray] = {}
    dims: dict[str, list[str]] = {}
    coords: dict[str, Sequence[Any]] = {}
    units: dict[str, str] = {}
    taken: set[str] = set()

    def claim(name: str, label: str, what: str) -> str:
        """Reserve a variable name, refusing a collision rather than overwriting.

        Data-group names are built by joining a dataset label to an axis or
        role name, so two datasets can in principle produce the same one — a
        dataset ``a`` with an extra coordinate ``b_c`` against a dataset
        ``a_b`` with one called ``c``. That is the flattening collision
        ``inference.md`` §4.4 names, and here it would silently drop one
        dataset's data on the floor: a stored run that quietly lost a
        dataset is worse than one that refuses to be written.
        """
        if name in taken:
            raise ResultsError(
                f"dataset {label!r}'s {what} would be stored as {name!r}, which another dataset "
                f"has already claimed. Data-group names join the dataset label to an axis or "
                f"role name, so labels that differ only by where an underscore falls can "
                f"collide. Rename one of the datasets."
            )
        taken.add(name)
        return name

    for label in problem.datasets:
        container = problem.datasets[label].observed
        names = _container_dims(label, container)
        as_list = list(names)
        value_unit = None if container.unit is None else str(container.unit.to_string())
        if container.values.dtype.kind == "c":
            # netCDF has no complex type; the two parts are stored separately
            # and named so that nothing mistakes one for the whole.
            for part, values in (
                ("real", container.values.real),
                ("imag", container.values.imag),
            ):
                name = claim(f"{label}_{part}", label, f"{part} part")
                observed[name] = np.asarray(values)
                dims[name] = as_list
                if value_unit is not None:
                    units[name] = value_unit
        else:
            name = claim(label, label, "values")
            observed[name] = np.asarray(container.values)
            dims[name] = as_list
            if value_unit is not None:
                units[name] = value_unit
        if len(names) == len(container.axes):
            for axis, dim in zip(container.axes, names, strict=True):
                coords[dim] = axis.values.tolist()
                if axis.unit is not None:
                    units[dim] = str(axis.unit.to_string())
        else:
            for axis in container.axes:
                name = claim(f"{label}_{axis.name}", label, f"{axis.name} axis")
                constant[name] = np.asarray(axis.values)
                dims[name] = as_list
                if axis.unit is not None:
                    units[name] = str(axis.unit.to_string())
        if container.uncertainty is not None:
            name = claim(f"{label}_uncertainty", label, "uncertainties")
            constant[name] = np.asarray(container.uncertainty)
            dims[name] = as_list
            if value_unit is not None:
                units[name] = value_unit
        if container.mask is not None:
            # int8, because netCDF has no boolean type; the sense is the
            # container's own (True excludes the sample), stated in the attr.
            name = claim(f"{label}_mask", label, "mask")
            constant[name] = np.asarray(container.mask, dtype=np.int8)
            dims[name] = as_list
        for extra, values in container.extra_coords.items():
            name = claim(f"{label}_{extra}", label, f"extra coordinate {extra!r}")
            constant[name] = np.asarray(values)
            dims[name] = as_list
    groups: dict[str, dict[str, np.ndarray]] = {}
    if observed:
        groups[OBSERVED_DATA_GROUP] = observed
    if constant:
        groups[CONSTANT_DATA_GROUP] = constant
    return groups, dims, coords, units


# ---------------------------------------------------------------------------
# emit
# ---------------------------------------------------------------------------


def emit(
    problem: FittingProblem,
    draws: Any,
    evaluations: Sequence[Any],
    *,
    engine: str | None = None,
    backend: str = "reference",
    coords: Mapping[str, Sequence[Any]] | None = None,
    observed: bool = True,
    extra_attrs: Mapping[str, object] | None = None,
) -> Any:
    """Build the run's :class:`xarray.DataTree` from its draws and evaluations.

    Parameters
    ----------
    problem
        The composed problem the draws are over. Everything structural — the
        parameter names, the dataset labels, the observed containers, the
        provenance — is read from it, so an engine driver supplies only what it
        alone knows.
    draws
        ``(chain, draw, free_size)`` or ``(draw, free_size)``, in
        ``problem.parameters`` order.
    evaluations
        One :class:`~ampere.core.dataset.Evaluation` per draw, nested per chain
        when there is more than one. ``FittingProblem.evaluate`` produces
        everything needed in one pass, which is why it, and not ``log_prob``,
        is what a driver should call.
    engine, backend
        Recorded in the provenance attrs.
    coords
        Coordinate values for named dimensions, e.g. ``{"objects": labels}`` for
        a plate. A plate routed element by element gets its coordinate from the
        bindings automatically; this overrides that and supplies one where the
        bindings cannot.
    observed
        Whether to emit the ``observed_data``/``constant_data`` groups. On by
        default: a stored run that cannot be plotted against its own data is
        half a record.
    extra_attrs
        Engine-specific provenance (step size, live points, walker count).

    Raises
    ------
    ampere.core.exceptions.OptionalDependencyError
        If arviz is not installed.
    ResultsError
        If the draws and evaluations do not describe the same run.
    """
    arviz = _require_arviz()
    array = _as_draw_array(draws, problem.free_size)
    chains, count = int(array.shape[0]), int(array.shape[1])
    grid = _as_evaluation_grid(evaluations, chains, count)

    posterior, dims, resolved = _posterior_variables(problem, array, coords)
    sampling: dict[str, dict[str, np.ndarray]] = {
        POSTERIOR_GROUP: posterior,
        SAMPLE_STATS_GROUP: _sample_stats(grid),
        LOG_LIKELIHOOD_GROUP: _log_likelihood_variables(problem, grid),
    }
    tree = arviz.from_dict(
        sampling,
        dims=dims,
        coords=resolved,
        attrs={
            LOG_LIKELIHOOD_GROUP: {
                f"{ATTR_PREFIX}decomposition": LOG_LIKELIHOOD_DECOMPOSITION,
                f"{ATTR_PREFIX}decomposition_note": (
                    "one term per dataset, not per observation: a GP likelihood does not "
                    "factorise over observations (likelihoods.md §16)."
                ),
            }
        },
    )
    if observed:
        groups, data_dims, data_coords, units = _data_groups(problem)
        if groups:
            data = arviz.from_dict(groups, dims=data_dims, coords=data_coords)
            for name, child in data.children.items():
                tree[name] = child
            _apply_units(tree, groups, units)
    tree.attrs.update(provenance_attrs(problem, engine=engine, backend=backend, extra=extra_attrs))
    tree.attrs[f"{ATTR_PREFIX}chains"] = chains
    tree.attrs[f"{ATTR_PREFIX}draws"] = count
    return tree


def _apply_units(
    tree: Any, groups: Mapping[str, Mapping[str, np.ndarray]], units: Mapping[str, str]
) -> None:
    """Record each data variable's and coordinate's unit as a CF-style attr.

    A mask is stored as ``int8`` because netCDF has no boolean type, so the group
    holding one also records which way round the ones read: the container's own
    convention, where ``True`` **excludes** a sample.
    """
    for group, variables in groups.items():
        dataset = tree[group].dataset
        for name, unit in units.items():
            if name in dataset.variables:
                tree[group][name].attrs["units"] = unit
        if any(name.endswith("_mask") for name in variables):
            tree[group].attrs[f"{ATTR_PREFIX}mask_convention"] = (
                "1 marks an excluded sample, following numpy/astropy (results_schema.md §7)."
            )


# ---------------------------------------------------------------------------
# A driver-facing accumulator
# ---------------------------------------------------------------------------


class DrawRecorder:
    """Accumulate draws and evaluations, then :meth:`emit` them.

    The shape an engine driver actually wants: it has one θ at a time and does
    not know how many it will end up with. Nothing here is required — an engine
    that already holds its chain as an array should call :func:`emit` directly —
    but it keeps the "record everything, every run" policy of
    ``DEVELOPMENT_PLAN.md`` §4.6 to one line at the call site rather than a
    bookkeeping structure per driver.
    """

    def __init__(self, problem: FittingProblem, *, chains: int = 1) -> None:
        if chains < 1:
            raise ResultsError(f"a run has at least one chain, got {chains}.")
        self.problem = problem
        self.chains = int(chains)
        self._draws: list[list[np.ndarray]] = [[] for _ in range(self.chains)]
        self._evaluations: list[list[Evaluation]] = [[] for _ in range(self.chains)]

    def record(
        self, values: Any = None, evaluation: Evaluation | None = None, *, chain: int = 0
    ) -> Evaluation:
        """Record one draw. Evaluates the problem when no evaluation is given.

        ``values`` may be a mapping of merged names or the flat free-parameter
        vector, and ``None`` means the reference θ — the same three things
        :meth:`~ampere.core.dataset.FittingProblem.evaluate` accepts, resolved
        here **before** the draw is packed. Resolving it in one place is the
        point: ``evaluate(None)`` scores the reference θ, so packing ``None``
        separately would store a row of NaN beside a perfectly good
        log-probability, and the two halves of a draw would disagree.
        """
        if not 0 <= chain < self.chains:
            raise ResultsError(f"chain {chain} is out of range for {self.chains} chain(s).")
        if values is None:
            values = self.problem.reference_values
        if evaluation is None:
            evaluation = self.problem.evaluate(values)
        theta = self.problem.parameters.pack(
            values if isinstance(values, Mapping) else self.problem.parameters.unpack(values)
        )
        self._draws[chain].append(theta)
        self._evaluations[chain].append(evaluation)
        return evaluation

    def __len__(self) -> int:
        return sum(len(chain) for chain in self._draws)

    def emit(self, **kwargs: Any) -> Any:
        """Build the :class:`xarray.DataTree`; keyword arguments go to :func:`emit`."""
        counts = {len(chain) for chain in self._draws}
        if counts == {0}:
            raise ResultsError("nothing has been recorded yet, so there is no run to emit.")
        if len(counts) != 1:
            raise ResultsError(
                f"the chains hold different numbers of draws ({sorted(counts)}); ArviZ's "
                f"(chain, draw) layout is rectangular. Trim them to a common length first."
            )
        array = np.asarray([np.asarray(chain) for chain in self._draws], dtype=float)
        return emit(self.problem, array, self._evaluations, **kwargs)


# ---------------------------------------------------------------------------
# netCDF
# ---------------------------------------------------------------------------


def to_netcdf(tree: Any, path: str | Path, *, engine: str | None = None) -> str:
    """Write a run to netCDF, the format ``DEVELOPMENT_PLAN.md`` §4.6 names.

    ``engine`` selects the netCDF backend (``"h5netcdf"`` or ``"netcdf4"``);
    ``None`` lets xarray choose whichever is installed.
    """
    _require_arviz()
    target = str(path)
    try:
        tree.to_netcdf(target, engine=engine)
    except (ImportError, ValueError) as error:
        raise _netcdf_dependency_error(error) from error
    return target


def from_netcdf(path: str | Path, *, engine: str | None = None) -> Any:
    """Read a run back. Accepts anything ``arviz.from_netcdf`` accepts."""
    arviz = _require_arviz()
    try:
        return arviz.from_netcdf(str(path), engine=engine)
    except (ImportError, ValueError) as error:
        raise _netcdf_dependency_error(error) from error


def _netcdf_dependency_error(error: Exception) -> Exception:
    """Turn xarray's "no backend" complaint into ampere's own, with the remedy.

    xarray reports a missing netCDF backend as a bare ``ValueError`` naming
    engines, which is not obviously an *installation* problem. arviz does not
    require an engine of its own, so ``pip install ampere[arviz]`` alone cannot
    write a file — that is worth saying plainly rather than leaving a user to
    decode a backend list.
    """
    text = str(error).lower()
    if "engine" in text or "backend" in text or "h5netcdf" in text or "netcdf4" in text:
        return OptionalDependencyError(
            "h5netcdf",
            extra="arviz",
            context="reading or writing a run as netCDF (any of h5netcdf+h5py or netCDF4 will "
            "do; arviz itself does not require one)",
        )
    return error
