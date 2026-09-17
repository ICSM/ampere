"""The on-disk training-set format: netCDF, the same as everything else.

``results.md`` §11 layer 2 fixed this format and deliberately did not implement
it — "the writer lands with Phase 2's SBI and emulator work, against the format
fixed here" (§13.9). This module is that writer, and the format is §11's table
unchanged:

===================== ======================================================
Group                 Contents
===================== ======================================================
root attrs            :func:`~ampere.results.provenance.provenance_attrs` —
                      the spec, problem, data and *(W3.12)* model hashes, the
                      versions, the seed; §9's recipe, not a second one
``theta``             one variable per merged parameter name,
                      ``(sample,)`` plus the parameter's own dimensions
``<model>.<channel>`` one group per model channel: ``values``, and
                      ``uncertainty``/``mask`` where present,
                      ``(sample,)`` + the channel's coordinate dimensions
``coordinates``       each channel's coordinate arrays, units as attributes
``sample_stats``      ``failed`` and the failure record, per draw
``observations``      *(added W2.8)* the noisy draws, one subgroup per
                      dataset label, when the budget drew any
===================== ======================================================

Three properties are why the format is netCDF and not JSON or a pickle, and
each of them is a decision the writer has to keep true rather than a fact about
the container format. **NaN is native**, so a masked or crashed sample needs no
sentinel. **The coordinate arrays are stored once for the whole set** rather
than per pair — which is what makes a coordinate-conditioned (neural-operator
style) emulator's training set the same size as a fixed-grid one, and which
this writer therefore *enforces*: a batch whose channel coordinates move from
sample to sample is refused rather than silently written with the last one. And
**the spec hash — and, since W3.12, the model hash — sit in the attributes**,
so ``DEVELOPMENT_PLAN.md`` §7's "spec-hash invalidation of trained artefacts"
is a comparison of strings — which :func:`append_training_set` makes on every
append, because a budget extended after the model was edited is exactly the
poisoned cache §7 warns about. The spec hash alone is only the *parameter*
declaration, though, so a likelihood family, noise model, solver or kernel
swap that leaves every parameter's name and prior unchanged used to pass it
unnoticed: :func:`~ampere.results.provenance.model_hash` (``ampere_model_hash``,
``PROVENANCE_SCHEMA_VERSION`` 6) closes that gap, and a file written before it
existed is refused by name rather than treated as an agreement it cannot make.

What this deliberately does not store is the model. An emulator is trained on
``(θ, ModelResult)`` pairs and validated against the spec hash; reconstructing
the simulator that produced them is out of scope for a training set and is what
the composed problem is for.

Two losses this closes
----------------------
``serialisation_review.md`` §4 named exactly two things the writer owed, "cheap
while nothing is archived". Both are closed here: θ keeps its dtype (through
:func:`~ampere.results.serialisation.training_pair_to_dict`, whose plain-data
form now carries one — ``CONTAINER_SCHEMA_VERSION`` 2), and
``Simulation.observations`` and the ``Failure`` *detail* — the reason, the
message, where, and the exception type, not merely the flag — are carried.

The append story, and its cost
------------------------------
:func:`append_training_set` is a real operation on an existing file: the
``sample`` dimension grows, and the pairs already written are still there
afterwards. It is implemented as read-concatenate-rewrite, which is ``O(existing
+ new)`` per call rather than ``O(new)``. That is a deliberate trade rather than
an oversight: an in-place HDF5 resize needs the whole file written with
unlimited dimensions and every variable resized in every group by hand, through
h5netcdf rather than xarray, which is a second serialisation path to maintain
for a saving that matters only once a budget outgrows memory. The extension
point is exactly that, and it is named here so nobody has to rediscover the
trade. What append buys today is that a budget is written in batches — the
simulations themselves, which are far larger than their arrays, need never be
held all at once.
"""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Iterable, Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any

import numpy as np

from ampere.core.dataset import FittingProblem, Simulation
from ampere.core.encoding import axis_identity_complaint
from ampere.core.exceptions import OptionalDependencyError, ResultsError
from ampere.core.results_schema import FunctionSamples, ModelResult
from ampere.core.simulate import SimulationBatch

from .emission import SAMPLE_STATS_GROUP, _container_dims
from .provenance import ATTR_PREFIX, PROVENANCE_SCHEMA_VERSION, provenance_attrs
from .serialisation import CONTAINER_SCHEMA_VERSION, container_from_dict

__all__ = [
    "COORDINATES_GROUP",
    "OBSERVATIONS_GROUP",
    "SAMPLE_DIM",
    "SAMPLE_STATS_GROUP",
    "THETA_GROUP",
    "TRAINING_SET_SCHEMA_VERSION",
    "TrainingSet",
    "append_training_set",
    "read_training_set",
    "write_training_set",
]

#: The dimension every per-draw variable is indexed by.
SAMPLE_DIM = "sample"

#: One variable per merged parameter name.
THETA_GROUP = "theta"

#: Every channel's coordinate arrays, stored once for the whole set.
COORDINATES_GROUP = "coordinates"

# SAMPLE_STATS_GROUP ("failed" and the failure record, per draw) is defined
# once in ampere.results.emission and imported above, so the two run-record
# formats cannot drift apart on the group name.

#: The parent of one subgroup per dataset label, for a budget drawn with
#: ``observe=True``. Added W2.8 (``serialisation_review.md`` §4's second loss);
#: absent from a file whose budget drew no observations.
OBSERVATIONS_GROUP = "observations"

#: Bumped whenever the meaning of a group or attribute above changes.
TRAINING_SET_SCHEMA_VERSION = 1


def _require_xarray() -> Any:
    """Import xarray on use, never on import.

    xarray arrives with arviz, which is a **base** dependency as of W2.2, so
    ``extra=None`` exactly as for arviz itself: an environment without it is
    incomplete rather than missing an extra. The training set is a plain
    ``DataTree`` rather than an ArviZ run — it holds no chains and no draws —
    so this module reaches for xarray directly instead of going through arviz
    to get at the same class.
    """
    try:
        import xarray
    except ImportError as error:  # pragma: no cover - exercised by a minimal install
        raise OptionalDependencyError(
            "xarray",
            context="reading or writing a training set (results.md §11's layer 2 is netCDF, "
            "through xarray's DataTree). xarray arrives with arviz, which is a base dependency "
            "of ampere, so this environment is incomplete rather than merely missing an extra",
        ) from error
    return xarray


# ---------------------------------------------------------------------------
# The shape of one channel, read off the batch
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _Slot:
    """Everything about one channel that is the same for every sample.

    A slot is a model channel (``"model.blue"``) or an observed dataset
    (``"observations/sed"``). What is *per sample* is the values, and the
    uncertainties, mask and extra coordinates where a container carries them;
    everything here — the kind, the axes, the coordinate arrays, the unit — is
    a property of the set, stored once. That division is ``results.md`` §11's
    "the coordinate arrays are stored once for the whole set rather than per
    pair", and it is checked rather than assumed (:func:`_check_slot`).
    """

    path: str
    kind: str
    axes: tuple[str, ...]
    dims: tuple[str, ...]
    shape: tuple[int, ...]
    unit: str | None
    fidelity: str | None
    coordinates: tuple[tuple[str, np.ndarray, str | None], ...]
    axis_variables: tuple[tuple[str, np.ndarray, str | None], ...]
    dtype: str
    complex: bool
    uncertainty: bool
    mask: bool
    extra: tuple[str, ...]
    meta: str | None


def _slot_of(path: str, container: FunctionSamples) -> _Slot:
    """The slot a container implies, with :func:`_container_dims`' naming.

    Reusing ``emission``'s dimension naming rather than inventing a second one
    is not tidiness: it means a training set's channel dimension is spelled the
    same way a stored run's observed dimension is, so a coordinate array read
    off one lines up with the other without a translation table.
    """
    # The dimension prefix is the group path with its separator flattened:
    # xarray refuses a variable or dimension name containing "/", since that is
    # what addresses a node. "observations/sed" is a path and
    # "observations.sed_spectral_axis" is its dimension.
    dims = _container_dims(path.replace("/", "."), container)
    coordinates: list[tuple[str, np.ndarray, str | None]] = []
    axis_variables: list[tuple[str, np.ndarray, str | None]] = []
    if len(dims) == len(container.axes):
        for axis, dim in zip(container.axes, dims, strict=True):
            coordinates.append(
                (dim, np.asarray(axis.values), None if axis.unit is None else str(axis.unit))
            )
    else:
        # A point kind with several axes indexes its samples jointly, so the
        # axes are ordinary arrays on one sample dimension rather than
        # coordinates — the same rule results.md §4 states for observed_data.
        for axis in container.axes:
            axis_variables.append(
                (axis.name, np.asarray(axis.values), None if axis.unit is None else str(axis.unit))
            )
        coordinates.append((dims[0], np.arange(container.n_samples), None))
    values = np.asarray(container.values)
    return _Slot(
        path=path,
        kind=type(container).__name__,
        axes=tuple(axis.name for axis in container.axes),
        dims=tuple(dims),
        shape=tuple(int(size) for size in values.shape),
        unit=None if container.unit is None else str(container.unit.to_string()),
        fidelity=container.fidelity,
        coordinates=tuple(coordinates),
        axis_variables=tuple(axis_variables),
        dtype=str(values.dtype.newbyteorder("<")),
        complex=values.dtype.kind == "c",
        uncertainty=container.uncertainty is not None,
        mask=container.mask is not None,
        extra=tuple(sorted(container.extra_coords)),
        meta=json.dumps(dict(container.meta), sort_keys=True) if container.meta else None,
    )


def _check_slot(slot: _Slot, container: FunctionSamples, index: int) -> None:
    """Refuse a container that does not fit the slot the set was opened with.

    Each of these is a way a training set would otherwise become quietly
    unreadable: a channel that changed kind, a coordinate grid that moved, an
    uncertainty that appeared halfway through. §11's "stored once for the whole
    set" is only a saving if it is also true.
    """
    other = _slot_of(slot.path, container)
    if other.kind != slot.kind or other.shape != slot.shape or other.dims != slot.dims:
        raise ResultsError(
            f"sample {index} of channel {slot.path!r} is a {other.kind} of shape {other.shape}; "
            f"the set holds {slot.kind}s of shape {slot.shape}. A training set is one channel "
            f"shape throughout — that is what lets its coordinates be stored once."
        )
    for (name, values, _), (_, theirs, _) in zip(slot.coordinates, other.coordinates, strict=False):
        if not np.array_equal(values, theirs):
            raise ResultsError(
                f"sample {index} of channel {slot.path!r} is on a different {name!r} grid from "
                f"the first. results.md §11 stores each channel's coordinates once for the whole "
                f"set, which is what makes a coordinate-conditioned emulator's training set the "
                f"size of a fixed-grid one; a set whose grid moves per sample cannot use that "
                f"format. Write one set per grid."
            )
    if other.uncertainty != slot.uncertainty or other.mask != slot.mask:
        raise ResultsError(
            f"sample {index} of channel {slot.path!r} "
            f"{'has' if other.uncertainty else 'has no'} uncertainties and "
            f"{'has' if other.mask else 'has no'} mask, where the set "
            f"{'has' if slot.uncertainty else 'has no'} and "
            f"{'has' if slot.mask else 'has no'} respectively. Presence is a property of the "
            f"set, since the variable either exists in the file or does not."
        )
    if other.extra != slot.extra:
        raise ResultsError(
            f"sample {index} of channel {slot.path!r} carries extra coordinates {other.extra}, "
            f"where the set carries {slot.extra}."
        )


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


#: What the writers accept: loose simulations, one batch, or an iterator of
#: chunks. ``simulate_many(..., as_chunks=True)`` produces the third, which is
#: the form a budget larger than memory arrives in.
Budget = Iterable[Simulation] | SimulationBatch | Iterable[SimulationBatch]


def _chunks_of(simulations: Budget) -> Iterator[Sequence[Simulation]]:
    """Normalise a budget into chunks, holding one chunk at a time.

    Three shapes reach the writers and all three mean the same thing —
    "these draws, in this order". A bare
    :class:`~ampere.core.simulate.SimulationBatch` is one chunk; an iterator of
    them (what ``simulate_many(as_chunks=True)`` yields) is one chunk each,
    pulled only when the previous one has been written, which is the whole
    point: nothing before the current chunk is still in memory. Anything else
    is a plain iterable of :class:`~ampere.core.dataset.Simulation`, and is one
    chunk, because a caller holding a list already holds it all.
    """
    if isinstance(simulations, SimulationBatch):
        yield simulations
        return
    iterator = iter(simulations)
    first = next(iterator, None)
    if first is None:
        return
    if isinstance(first, SimulationBatch):
        yield first
        for chunk in iterator:
            if not isinstance(chunk, SimulationBatch):
                raise ResultsError(
                    f"this budget mixes SimulationBatch chunks with "
                    f"{type(chunk).__name__}; pass either the chunks "
                    f"simulate_many(as_chunks=True) yields or a flat sequence of Simulations, "
                    f"not both."
                )
            yield chunk
        return
    first_draw: Simulation = first
    yield [first_draw, *iterator]


def write_training_set(
    path: str | Path,
    simulations: Budget,
    problem: FittingProblem,
    *,
    engine: str | None = None,
) -> str:
    """Write a simulation budget as ``results.md`` §11's netCDF training set.

    Parameters
    ----------
    path
        Where to write. An existing file is replaced; use
        :func:`append_training_set` to grow one.
    simulations
        What :meth:`~ampere.core.dataset.FittingProblem.simulate` produced, in
        order. Failed draws are **written**, not dropped: ``inference.md`` §13's
        "reject-and-record rather than train on garbage" is only useful if the
        record reaches the file, and a budget's failure rate is a property of
        the prior worth measuring.

        Since W3.1 this also takes a
        :class:`~ampere.core.simulate.SimulationBatch`, or the **iterator of
        chunks** ``simulate_many(..., as_chunks=True)`` yields. The chunked form
        is written a chunk at a time — the first chunk written, the rest
        appended — so a budget larger than memory reaches the file without ever
        being held. That makes §13's limitation 9 append the cost of a very
        large budget rather than the peak memory: each chunk after the first
        pays ``O(existing + new)``, which is quadratic in the number of chunks.
        Prefer few large chunks over many small ones.
    problem
        The problem the budget was over. Its hashes go in the root attributes,
        and its ``ampere_spec_hash`` is what invalidates a stale artefact.
    engine
        The netCDF backend; ``None`` lets xarray choose.

    Returns
    -------
    str
        The path written.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the budget is empty, if every draw failed (so no channel shape is
        known and there is nothing to train on), or if the channels are not one
        shape on one grid throughout.
    """
    written: str | None = None
    for chunk in _chunks_of(simulations):
        if not chunk:
            continue
        if written is None:
            written = _write_first(path, chunk, problem, engine=engine)
        else:
            written = append_training_set(path, chunk, problem, engine=engine)
    if written is None:
        raise ResultsError("a training set needs at least one simulation; none was given.")
    return written


def _write_first(
    path: str | Path,
    batch: Sequence[Simulation],
    problem: FittingProblem,
    *,
    engine: str | None,
) -> str:
    """The first (or only) chunk: the file is created, replacing anything there."""
    slots = _slots_from(batch)
    if not slots:
        raise ResultsError(
            "every simulation in this budget failed, so there is no channel to write and "
            "nothing to train on. The failure counts are on the problem "
            "(FittingProblem.failure_counts) and are worth reading before spending another "
            "budget."
        )
    tree = _tree_from(batch, slots, problem, offset=0)
    return _write(tree, path, engine=engine)


def append_training_set(
    path: str | Path,
    simulations: Budget,
    problem: FittingProblem,
    *,
    engine: str | None = None,
) -> str:
    """Grow an existing training set by one batch. The ``sample`` dimension grows.

    Three checks run first, and any of them can **refuse**. The spec hash: appending
    draws from an edited model to a budget written before the edit produces a
    file whose two halves came from different simulators and whose
    provenance says they did not — precisely ``DEVELOPMENT_PLAN.md`` §7's
    poisoned cache, and a string comparison because §11 put the hash in the
    attributes so that it could be. And, since **W3.12**, the model hash
    (``ampere_model_hash``, schema 6): the spec hash alone is only the
    *parameter* declaration, so a likelihood family, noise model, solver or
    kernel swap that leaves every parameter's name and prior unchanged used
    to pass the spec-hash check unnoticed — plan §7's trap arriving by the
    other door, and the real gap W2.8's original check left open. A file
    written before the model hash existed (schema < 6: no
    ``ampere_model_hash`` attribute at all) is refused too, by name, rather
    than treated as an agreement it cannot actually make.

    And, since **W5.11**, the recorded **encoding layout**
    (:func:`_check_encoding_layout`): a budget packed before the
    axis-identity columns existed has narrower rows than one packed now, so
    appending would put two packings in one file under one hash. That
    refusal names the axis identity rather than reporting a hash mismatch,
    because the reason is knowable and the remedy follows from it.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the file's spec hash differs from the problem's, if the file's
        model hash differs from the problem's, if the file predates the
        model hash, if the file records an encoding layout from before
        W5.11, if the file is not a training set, or if the batch's channels
        do not match the file's.
    """
    chunks = [chunk for chunk in _chunks_of(simulations) if chunk]
    if not chunks:
        raise ResultsError("nothing to append; the batch is empty.")
    if len(chunks) > 1:
        written = path
        for chunk in chunks:
            written = append_training_set(path, chunk, problem, engine=engine)
        return str(written)
    batch = chunks[0]
    xarray = _require_xarray()
    existing = _open(path)
    current = provenance_attrs(problem)
    stored_spec = existing.attrs.get(f"{ATTR_PREFIX}spec_hash")
    if stored_spec != current.get(f"{ATTR_PREFIX}spec_hash"):
        raise ResultsError(
            f"this training set was written from a different declaration: it records "
            f"spec_hash {stored_spec!r} and the problem given hashes to "
            f"{current.get(f'{ATTR_PREFIX}spec_hash')!r}. Appending would leave one file whose "
            f"halves came from two simulators, which is the stale-artefact trap the hash exists "
            f"to catch (DEVELOPMENT_PLAN.md §7). Write a new set."
        )
    model_hash_key = f"{ATTR_PREFIX}model_hash"
    if model_hash_key not in existing.attrs:
        raise ResultsError(
            f"this training set predates ampere_model_hash: it was written under "
            f"PROVENANCE_SCHEMA_VERSION < {PROVENANCE_SCHEMA_VERSION} (W3.12), so its spec hash "
            f"agrees with the problem given but there is no recorded model hash to compare — a "
            f"likelihood family, noise model, solver or kernel swap that left every parameter's "
            f"name and prior unchanged could have been appended onto this file unnoticed. Write "
            f"a new set to record the model hash from the start."
        )
    stored_model = existing.attrs[model_hash_key]
    if stored_model != current.get(model_hash_key):
        raise ResultsError(
            f"this training set was written from a different model: it records "
            f"model_hash {stored_model!r} and the problem given hashes to "
            f"{current.get(model_hash_key)!r}. The parameter declaration agrees (the spec hash "
            f"matches) but the likelihood family, noise model, solver or kernel does not — "
            f"appending would leave one file whose halves came from two simulators, which is the "
            f"stale-artefact trap the hash exists to catch (DEVELOPMENT_PLAN.md §7). Write a new "
            f"set."
        )
    _check_encoding_layout(existing)
    offset = int(existing.attrs.get(f"{ATTR_PREFIX}samples", 0))
    slots = _slots_from(batch)
    if slots:
        _check_against_file(existing, slots)
    else:
        # A batch in which every draw failed carries no container to read a
        # shape off, and is still worth appending: a budget's failure rate is
        # data. The file already knows the shapes, so they are read back from
        # it and the batch is written as NaN rows with their failure records.
        slots = _slots_from_file(existing)
    addition = _tree_from(batch, slots, problem, offset=offset)
    merged = _concatenate(xarray, existing, addition)
    merged.attrs[f"{ATTR_PREFIX}samples"] = offset + len(batch)
    return _write(merged, path, engine=engine)


def _check_encoding_layout(existing: Any) -> None:
    """Refuse a training set whose recorded packing predates **W5.11**.

    A budget written under an older encoding was packed without the
    axis-identity columns, so its rows are narrower than this ampere's and
    every layout hash it recorded moved when W5.11 landed. Appending would
    leave one file whose halves were encoded under two packings — the
    poisoned cache of ``DEVELOPMENT_PLAN.md`` §7 again, arriving by the
    encoding's door rather than the model's — and a reader comparing the two
    hashes would see only that they differ.

    The refusal therefore names the reason:
    :func:`~ampere.core.encoding.axis_identity_complaint` writes it once and
    this raises it as the results-layer refusal every other check here
    raises. A file that records no layout at all (a budget written without an
    SBI run's encoding provenance) is not refused: there is nothing to compare,
    and inventing a disagreement would be as dishonest as missing one.
    """
    stored = existing.attrs.get(f"{ATTR_PREFIX}encoding_layout")
    if not stored:
        return
    try:
        record = json.loads(stored) if isinstance(stored, str) else stored
    except (TypeError, ValueError):
        return
    complaint = axis_identity_complaint(record)
    if complaint is not None:
        raise ResultsError(
            f"this training set was written under a different encoding: {complaint} The rows "
            f"already in this file are narrower than the ones being appended, so the two halves "
            f"would be packed differently under one hash (DEVELOPMENT_PLAN.md §7). Write a new "
            f"set."
        )


def _slots_from(batch: Sequence[Simulation]) -> dict[str, _Slot]:
    """The slots a batch implies, checked for consistency across its samples."""
    slots: dict[str, _Slot] = {}
    for index, simulation in enumerate(batch):
        for path, container in _containers_of(simulation):
            if path not in slots:
                slots[path] = _slot_of(path, container)
            else:
                _check_slot(slots[path], container, index)
    return dict(sorted(slots.items()))


def _containers_of(simulation: Simulation) -> list[tuple[str, FunctionSamples]]:
    """``(group path, container)`` for one simulation's channels and observations."""
    found: list[tuple[str, FunctionSamples]] = []
    for model, result in simulation.results.items():
        for channel in result:
            found.append((f"{model}.{channel}", result[channel]))
    if simulation.observations is not None:
        for label, container in simulation.observations.items():
            found.append((f"{OBSERVATIONS_GROUP}/{label}", container))
    return found


def _tree_from(
    batch: Sequence[Simulation],
    slots: Mapping[str, _Slot],
    problem: FittingProblem,
    *,
    offset: int,
) -> Any:
    """Build the ``DataTree`` for one batch, in §11's layout."""
    xarray = _require_xarray()
    count = len(batch)
    groups: dict[str, Any] = {
        THETA_GROUP: xarray.Dataset(_theta_variables(batch, problem)),
        SAMPLE_STATS_GROUP: xarray.Dataset(_sample_stats(batch)),
    }
    coordinates: dict[str, Any] = {}
    for path, slot in slots.items():
        groups[path] = _slot_dataset(xarray, slot, batch, count)
        for dim, values, unit in slot.coordinates:
            array = xarray.DataArray(values, dims=(dim,))
            if unit is not None:
                array.attrs["units"] = unit
            coordinates[dim] = array
    if coordinates:
        groups[COORDINATES_GROUP] = xarray.Dataset(coordinates)
    tree = xarray.DataTree.from_dict(groups)
    tree.attrs.update(
        provenance_attrs(
            problem,
            extra={
                "training_set_version": TRAINING_SET_SCHEMA_VERSION,
                "container_schema": CONTAINER_SCHEMA_VERSION,
                "samples": offset + count,
                "channels": sorted(slots),
                **dict(getattr(batch, "provenance", None) or {}),
            },
        )
    )
    return tree


def _theta_variables(batch: Sequence[Simulation], problem: FittingProblem) -> dict[str, Any]:
    """One variable per merged parameter name, dtype and all.

    ``Simulation.parameters`` is keyed by merged name and is what §11's table
    asks for. The dtype is taken from the first sample and the rest are cast to
    it, so an integer plate index stays integral through the file — the loss
    ``serialisation_review.md`` §4 named, closed here as well as in the
    plain-data form.

    The **free** parameters only. A fixed value is part of the declaration, not
    an input an emulator varies, and it is already in ``ampere_spec_hash``,
    which is what a training set is validated against; writing a column of one
    repeated number per fixed parameter would grow every file for nothing.
    """
    variables: dict[str, Any] = {}
    for parameter in problem.parameters:
        if parameter.is_fixed:
            continue
        name = parameter.name
        column = [np.asarray(simulation.parameters[name]) for simulation in batch]
        stacked = np.stack([entry.astype(column[0].dtype, copy=False) for entry in column])
        dims = (SAMPLE_DIM, *(f"{name}_dim_{axis}" for axis in range(stacked.ndim - 1)))
        variables[name] = (dims, stacked)
    return variables


def _sample_stats(batch: Sequence[Simulation]) -> dict[str, Any]:
    """``failed`` plus the whole :class:`~ampere.core.dataset.Failure` record.

    ``serialisation_review.md`` §4's second named loss was that "``Failure``
    detail [is] not carried (only the ``failed`` flag)". The flag says a draw
    was rejected; the reason, the message, where it happened and the exception's
    type say whether the budget is worth re-running with a wider prior or the
    simulator is broken, which is the question a 2 % crash rate actually poses.
    ``int8`` because netCDF has no boolean type.
    """
    failed = np.array([int(bool(simulation.failed)) for simulation in batch], dtype=np.int8)
    records = [
        {} if simulation.failure is None else simulation.failure.to_dict() for simulation in batch
    ]
    text = {
        "failure_reason": "reason",
        "failure_message": "message",
        "failure_where": "where",
        "failure_exception_type": "exception_type",
    }
    variables: dict[str, Any] = {"failed": (SAMPLE_DIM, failed)}
    for variable, key in text.items():
        variables[variable] = (
            SAMPLE_DIM,
            np.array([str(record.get(key, "")) for record in records], dtype=object).astype(str),
        )
    variables["failure_values"] = (
        SAMPLE_DIM,
        np.array(
            [json.dumps(record.get("values", {}), sort_keys=True) for record in records],
            dtype=object,
        ).astype(str),
    )
    return variables


def _slot_dataset(xarray: Any, slot: _Slot, batch: Sequence[Simulation], count: int) -> Any:
    """One channel's per-sample arrays, with the slot's constants as attributes."""
    dims = (SAMPLE_DIM, *slot.dims)
    values = _filled(slot.dtype, (count, *slot.shape))
    uncertainty = _filled("float64", (count, *slot.shape)) if slot.uncertainty else None
    mask = np.zeros((count, *slot.shape), dtype=np.int8) if slot.mask else None
    extra = {name: _filled("float64", (count, *slot.shape)) for name in slot.extra}
    for index, simulation in enumerate(batch):
        container = _container_at(simulation, slot.path)
        if container is None:
            continue
        values[index] = np.asarray(container.values)
        if uncertainty is not None and container.uncertainty is not None:
            uncertainty[index] = np.asarray(container.uncertainty)
        if mask is not None and container.mask is not None:
            mask[index] = np.asarray(container.mask, dtype=np.int8)
        for name in slot.extra:
            extra[name][index] = np.asarray(container.extra_coords[name])
    variables: dict[str, Any] = {}
    if slot.complex:
        # netCDF has no complex type; the parts are named so nothing mistakes
        # one for the whole, exactly as results.md §4 does for observed data.
        variables["values_real"] = (dims, values.real)
        variables["values_imag"] = (dims, values.imag)
    else:
        variables["values"] = (dims, values)
    if uncertainty is not None:
        variables["uncertainty"] = (dims, uncertainty)
    if mask is not None:
        variables["mask"] = (dims, mask)
    for name, array in extra.items():
        variables[f"extra_{name}"] = (dims, array)
    for name, axis_values, _ in slot.axis_variables:
        variables[f"axis_{name}"] = (slot.dims, axis_values)
    dataset = xarray.Dataset(variables)
    dataset.attrs.update(
        {
            f"{ATTR_PREFIX}kind": slot.kind,
            f"{ATTR_PREFIX}axes": json.dumps(list(slot.axes)),
            f"{ATTR_PREFIX}unit": "" if slot.unit is None else slot.unit,
            f"{ATTR_PREFIX}fidelity": "" if slot.fidelity is None else slot.fidelity,
            f"{ATTR_PREFIX}value_dtype": slot.dtype,
            f"{ATTR_PREFIX}complex": int(slot.complex),
            f"{ATTR_PREFIX}extra_coords": json.dumps(list(slot.extra)),
            f"{ATTR_PREFIX}axis_units": json.dumps(
                {name: unit for name, _, unit in slot.axis_variables if unit is not None}
            ),
            f"{ATTR_PREFIX}coordinate_units": json.dumps(
                {dim: unit for dim, _, unit in slot.coordinates if unit is not None}
            ),
        }
    )
    if slot.meta is not None:
        dataset.attrs[f"{ATTR_PREFIX}meta"] = slot.meta
    return dataset


def _filled(dtype: str, shape: tuple[int, ...]) -> np.ndarray:
    """An array whose "not written" value says so where the dtype can.

    NaN for anything floating or complex, which is §11's whole reason for
    netCDF: "NaN is native, so a masked or crashed sample needs no sentinel".
    An integer channel has no NaN, so it is zero-filled and ``sample_stats``'
    ``failed`` flag is the authority for that sample — recorded here rather
    than left as a surprise.
    """
    kind = np.dtype(dtype).kind
    if kind in "fc":
        return np.full(shape, np.nan, dtype=dtype)
    return np.zeros(shape, dtype=dtype)


def _container_at(simulation: Simulation, path: str) -> FunctionSamples | None:
    """The container one simulation holds at a slot path, or ``None``."""
    if path.startswith(f"{OBSERVATIONS_GROUP}/"):
        if simulation.observations is None:
            return None
        return simulation.observations.get(path.split("/", 1)[1])
    model, _, channel = path.partition(".")
    result = simulation.results.get(model)
    if result is None or channel not in result:
        return None
    return result[channel]


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class TrainingSet:
    """A training set read back off disk, by value.

    The pairs are rebuilt lazily rather than all at once: an emulator trainer
    iterates a budget it deliberately did not hold in memory while writing, and
    handing it every :class:`~ampere.core.results_schema.ModelResult` at once
    would undo that.

    Attributes
    ----------
    attrs
        The root attributes, ``ampere_spec_hash`` among them — the string a
        cached artefact is validated against.
    theta
        One array per merged parameter name, ``(sample,)`` plus the parameter's
        own dimensions, in the dtype it was written with.
    failed
        Boolean per sample.
    channels
        Model-channel paths, in file order (``"model.blue"``).
    observed
        Dataset labels whose observations were carried, or an empty tuple.
    """

    attrs: Mapping[str, Any]
    theta: Mapping[str, np.ndarray]
    failed: np.ndarray
    failures: tuple[Mapping[str, Any] | None, ...]
    channels: tuple[str, ...]
    observed: tuple[str, ...]
    _groups: Mapping[str, Any] = dataclasses.field(repr=False, default_factory=dict)
    _coordinates: Any = dataclasses.field(repr=False, default=None)

    def __len__(self) -> int:
        return int(self.failed.size)

    @property
    def spec_hash(self) -> str:
        """``ampere_spec_hash`` — what invalidates an artefact trained on this."""
        return str(self.attrs.get(f"{ATTR_PREFIX}spec_hash", ""))

    def parameters(self, index: int) -> dict[str, Any]:
        """θ at one sample, by merged name."""
        return {name: values[index] for name, values in self.theta.items()}

    def result(self, index: int) -> ModelResult:
        """The :class:`~ampere.core.results_schema.ModelResult` at one sample.

        Refuses a failed sample by name rather than returning a container of
        NaN: "this draw was rejected" is a different statement from "this draw
        predicted nothing", and a trainer that silently fitted the second would
        be training on exactly the garbage ``inference.md`` §13 rejects.
        """
        if bool(self.failed[index]):
            record = self.failures[index] or {}
            raise ResultsError(
                f"sample {index} was a failed simulation "
                f"({record.get('reason', 'reason not recorded')}: "
                f"{record.get('message', '')}), so it has no ModelResult. Filter on "
                f"TrainingSet.failed before training."
            )
        channels = {path.partition(".")[2]: self._container(path, index) for path in self.channels}
        return ModelResult(channels, parameters=self.parameters(index))

    def observations(self, index: int) -> dict[str, FunctionSamples]:
        """The noisy draws at one sample, by dataset label. Empty when none."""
        return {
            label: self._container(f"{OBSERVATIONS_GROUP}/{label}", index)
            for label in self.observed
        }

    def pair(self, index: int) -> dict[str, Any]:
        """One training pair in :func:`training_pair_from_dict`'s form.

        The same mapping a caller gets from the in-memory round trip, so code
        written against one works against the other — which is the point of
        having both.
        """
        observations = self.observations(index)
        pair: dict[str, Any] = {
            "theta": self.parameters(index),
            "failed": bool(self.failed[index]),
        }
        if not pair["failed"]:
            pair["result"] = self.result(index)
        if observations:
            pair["observations"] = observations
        record = self.failures[index]
        if record is not None:
            pair["failure"] = dict(record)
        return pair

    def _container(self, path: str, index: int) -> FunctionSamples:
        """Rebuild one container by handing its plain-data form to the reader.

        Going through :func:`~ampere.results.serialisation.container_from_dict`
        rather than constructing the class here means the whole base contract is
        re-checked on the way back in — axis names and physical types, shapes,
        non-negative uncertainties, a strictly boolean mask — and that an
        out-of-tree kind registered with ``register_kind`` reads back for free.
        """
        return container_from_dict(_record_at(self._groups[path], self._coordinates, index))


def read_training_set(path: str | Path, *, engine: str | None = None) -> TrainingSet:
    """Read a training set written by :func:`write_training_set`.

    Everything is loaded eagerly and the file is closed, because the alternative
    — a lazy tree holding an open handle — makes
    :func:`append_training_set`'s rewrite of the same path fail on Windows and
    succeed confusingly everywhere else.

    ``engine`` selects the netCDF backend; ``None`` lets xarray choose, which is
    what the writer does too — and reading a file through a different HDF5
    binding from the one that wrote it is a deadlock rather than an error (see
    :func:`_open`).
    """
    tree = _open(path, engine=engine)
    groups = {
        str(name): tree[name].dataset
        for name in tree.children
        if name not in (THETA_GROUP, SAMPLE_STATS_GROUP, COORDINATES_GROUP, OBSERVATIONS_GROUP)
    }
    observed: list[str] = []
    if OBSERVATIONS_GROUP in tree.children:
        for label in tree[OBSERVATIONS_GROUP].children:
            groups[f"{OBSERVATIONS_GROUP}/{label}"] = tree[f"{OBSERVATIONS_GROUP}/{label}"].dataset
            observed.append(str(label))
    coordinates = tree[COORDINATES_GROUP].dataset if COORDINATES_GROUP in tree.children else None
    stats = tree[SAMPLE_STATS_GROUP].dataset
    failed = np.asarray(stats["failed"].values).astype(bool)
    return TrainingSet(
        attrs=dict(tree.attrs),
        theta={
            str(name): np.asarray(tree[THETA_GROUP][name].values)
            for name in tree[THETA_GROUP].dataset.data_vars
        },
        failed=failed,
        failures=tuple(_failure_records(stats, failed)),
        channels=tuple(sorted(name for name in groups if "/" not in name)),
        observed=tuple(sorted(observed)),
        _groups=groups,
        _coordinates=coordinates,
    )


def _failure_records(stats: Any, failed: np.ndarray) -> list[Mapping[str, Any] | None]:
    """One :meth:`~ampere.core.dataset.Failure.to_dict` mapping per failed sample."""
    records: list[Mapping[str, Any] | None] = []
    for index in range(int(failed.size)):
        if not bool(failed[index]):
            records.append(None)
            continue
        record = {
            key: str(stats[f"failure_{key}"].values[index])
            for key in ("reason", "message", "where", "exception_type")
            if f"failure_{key}" in stats.variables
        }
        if "failure_values" in stats.variables:
            raw = str(stats["failure_values"].values[index])
            record["values"] = json.loads(raw) if raw else {}
        records.append(record)
    return records


def _record_at(dataset: Any, coordinates: Any, index: int) -> dict[str, Any]:
    """One sample of one channel, as :func:`container_to_dict`'s plain form."""
    attrs = dataset.attrs
    axes = json.loads(attrs[f"{ATTR_PREFIX}axes"])
    axis_units = json.loads(attrs.get(f"{ATTR_PREFIX}axis_units", "{}"))
    coordinate_units = json.loads(attrs.get(f"{ATTR_PREFIX}coordinate_units", "{}"))
    complex_valued = bool(int(attrs.get(f"{ATTR_PREFIX}complex", 0)))
    dims = list(dataset["values_imag" if complex_valued else "values"].dims)[1:]
    encoded_coordinates: dict[str, Any] = {}
    if len(axes) == len(dims):
        for axis, dim in zip(axes, dims, strict=True):
            values = np.asarray(coordinates[dim].values)
            encoded_coordinates[axis] = {
                "values": values.tolist(),
                "unit": coordinate_units.get(dim),
            }
    else:
        for axis in axes:
            values = np.asarray(dataset[f"axis_{axis}"].values)
            encoded_coordinates[axis] = {
                "values": values.tolist(),
                "unit": axis_units.get(axis),
            }
    dtype = str(attrs[f"{ATTR_PREFIX}value_dtype"])
    if complex_valued:
        values_record = {
            "dtype": dtype,
            "real": np.asarray(dataset["values_real"].values[index]).tolist(),
            "imag": np.asarray(dataset["values_imag"].values[index]).tolist(),
        }
    else:
        values_record = {
            "dtype": dtype,
            "data": np.asarray(dataset["values"].values[index]).tolist(),
        }
    unit = str(attrs.get(f"{ATTR_PREFIX}unit", ""))
    fidelity = str(attrs.get(f"{ATTR_PREFIX}fidelity", ""))
    record: dict[str, Any] = {
        "version": CONTAINER_SCHEMA_VERSION,
        "kind": str(attrs[f"{ATTR_PREFIX}kind"]),
        "coordinates": encoded_coordinates,
        "values": values_record,
        "unit": unit or None,
        "uncertainty": (
            {
                "dtype": "float64",
                "data": np.asarray(dataset["uncertainty"].values[index]).tolist(),
            }
            if "uncertainty" in dataset.variables
            else None
        ),
        "mask": (
            np.asarray(dataset["mask"].values[index]).astype(bool).tolist()
            if "mask" in dataset.variables
            else None
        ),
        "extra_coords": {
            name: {
                "dtype": "float64",
                "data": np.asarray(dataset[f"extra_{name}"].values[index]).tolist(),
            }
            for name in json.loads(attrs.get(f"{ATTR_PREFIX}extra_coords", "[]"))
        },
        "fidelity": fidelity or None,
    }
    meta = attrs.get(f"{ATTR_PREFIX}meta")
    if meta:
        record["meta"] = json.loads(meta)
    return record


# ---------------------------------------------------------------------------
# File handling
# ---------------------------------------------------------------------------


def _write(tree: Any, path: str | Path, *, engine: str | None) -> str:
    """Write the tree, turning a missing netCDF backend into ampere's own error."""
    target = str(path)
    try:
        tree.to_netcdf(target, engine=engine)
    except (ImportError, TypeError, ValueError) as error:
        text = str(error).lower()
        if "engine" in text or "backend" in text or "h5netcdf" in text or "netcdf4" in text:
            raise OptionalDependencyError(
                "h5netcdf",
                context="writing a training set as netCDF (results.md §11). h5netcdf is a base "
                "dependency of ampere, so this environment is incomplete",
            ) from error
        raise
    return target


def _open(path: str | Path, *, engine: str | None = None) -> Any:
    """Read a training set into memory and close the file behind us.

    ``engine=None`` — xarray's own choice, which is what
    :func:`~ampere.results.emission.from_netcdf` does too, and it matters here
    for a reason beyond consistency. Forcing ``"h5netcdf"`` on the read while
    the write took xarray's default meant one file written through netCDF4's
    HDF5 and reopened through h5py's, and HDF5's file locking will *deadlock*
    that pair rather than fail it — a hang with no traceback, seen once in a
    full-suite run and not reproducible on the suite alone. One engine per
    file, chosen the same way at both ends, is the fix; the parameter is still
    there for a caller who has a reason.

    Everything is loaded eagerly and the file closed, because a lazy tree
    holding an open handle is what :func:`append_training_set` would then have
    to rewrite underneath itself.
    """
    xarray = _require_xarray()
    target = Path(path)
    if not target.exists():
        raise ResultsError(f"no training set at {str(target)!r}.")
    opened = xarray.open_datatree(str(target), engine=engine)
    try:
        tree = opened.load()
    finally:
        opened.close()
    if f"{ATTR_PREFIX}training_set_version" not in tree.attrs:
        raise ResultsError(
            f"{str(target)!r} is not an ampere training set: it carries no "
            f"{ATTR_PREFIX}training_set_version attribute. A stored *run* is read with "
            f"ampere.results.from_netcdf instead."
        )
    version = int(tree.attrs[f"{ATTR_PREFIX}training_set_version"])
    if version != TRAINING_SET_SCHEMA_VERSION:
        raise ResultsError(
            f"unsupported training-set version {version!r}; this ampere writes and reads "
            f"version {TRAINING_SET_SCHEMA_VERSION}."
        )
    return tree


def _slots_from_file(tree: Any) -> dict[str, _Slot]:
    """Rebuild the slots a file already holds, from its own groups.

    Everything a slot needs is in the file, because everything a slot holds is
    what the writer stored once: the kind, the axes, the unit, the dtype, the
    coordinate arrays. Reading it back is how an all-failed batch can still be
    appended — it has no container to take a shape from and the file does.
    """
    coordinates = tree[COORDINATES_GROUP].dataset if COORDINATES_GROUP in tree.children else None
    slots: dict[str, _Slot] = {}
    for path in _group_paths(tree):
        if path in (THETA_GROUP, SAMPLE_STATS_GROUP, COORDINATES_GROUP):
            continue
        dataset = tree[path].dataset
        attrs = dataset.attrs
        if f"{ATTR_PREFIX}kind" not in attrs:  # pragma: no cover - a foreign group
            continue
        values = dataset["values_real" if "values_real" in dataset.variables else "values"]
        dims = tuple(str(dim) for dim in values.dims[1:])
        axes = tuple(json.loads(attrs[f"{ATTR_PREFIX}axes"]))
        axis_units = json.loads(attrs.get(f"{ATTR_PREFIX}axis_units", "{}"))
        coordinate_units = json.loads(attrs.get(f"{ATTR_PREFIX}coordinate_units", "{}"))
        unit = str(attrs.get(f"{ATTR_PREFIX}unit", ""))
        fidelity = str(attrs.get(f"{ATTR_PREFIX}fidelity", ""))
        slots[path] = _Slot(
            path=path,
            kind=str(attrs[f"{ATTR_PREFIX}kind"]),
            axes=axes,
            dims=dims,
            shape=tuple(int(size) for size in values.shape[1:]),
            unit=unit or None,
            fidelity=fidelity or None,
            coordinates=tuple(
                (
                    dim,
                    np.asarray(coordinates[dim].values),
                    coordinate_units.get(dim),
                )
                for dim in dims
                if coordinates is not None and dim in coordinates.variables
            ),
            axis_variables=tuple(
                (name, np.asarray(dataset[f"axis_{name}"].values), axis_units.get(name))
                for name in axes
                if f"axis_{name}" in dataset.variables
            ),
            dtype=str(attrs[f"{ATTR_PREFIX}value_dtype"]),
            complex=bool(int(attrs.get(f"{ATTR_PREFIX}complex", 0))),
            uncertainty="uncertainty" in dataset.variables,
            mask="mask" in dataset.variables,
            extra=tuple(json.loads(attrs.get(f"{ATTR_PREFIX}extra_coords", "[]"))),
            meta=attrs.get(f"{ATTR_PREFIX}meta"),
        )
    return dict(sorted(slots.items()))


def _check_against_file(tree: Any, slots: Mapping[str, _Slot]) -> None:
    """Refuse a batch whose channels are not the file's."""
    stored = sorted(json.loads(tree.attrs.get(f"{ATTR_PREFIX}channels", "[]")))
    if stored != sorted(slots):
        raise ResultsError(
            f"this batch holds channels {sorted(slots)} and the training set holds {stored}. "
            f"A set is one channel layout throughout; appending a different one would give it "
            f"two shapes and no way to say which sample has which."
        )


def _concatenate(xarray: Any, existing: Any, addition: Any) -> Any:
    """Grow the ``sample`` dimension, leaving the shared coordinates alone."""
    groups: dict[str, Any] = {}
    for name in sorted(_group_paths(existing)):
        old = existing[name].dataset
        if name == COORDINATES_GROUP:
            _check_coordinates(old, addition[name].dataset)
            groups[name] = old
            continue
        new = addition[name].dataset
        groups[name] = xarray.concat([old, new], dim=SAMPLE_DIM, data_vars="minimal")
    merged = xarray.DataTree.from_dict(groups)
    merged.attrs.update(dict(existing.attrs))
    return merged


def _group_paths(tree: Any) -> list[str]:
    """Every leaf group path in a training set, ``observations/<label>`` included."""
    paths: list[str] = []
    for name in tree.children:
        node = tree[name]
        if node.children:
            paths.extend(f"{name}/{child}" for child in node.children)
        else:
            paths.append(str(name))
    return paths


def _check_coordinates(old: Any, new: Any) -> None:
    """Refuse a batch on a different grid from the file's."""
    for name in old.variables:
        if name not in new.variables:
            continue
        if not np.array_equal(np.asarray(old[name].values), np.asarray(new[name].values)):
            raise ResultsError(
                f"this batch is on a different {str(name)!r} grid from the training set's. "
                f"results.md §11 stores each channel's coordinates once for the whole set, so a "
                f"set cannot hold two grids."
            )
