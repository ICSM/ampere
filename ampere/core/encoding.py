"""The coordinate--value--mask encoding: one tensor per observation (W3.3).

``docs/design/contracts/encoding.md`` is this module's contract, and this module
is what binds it. The problem it solves is stated there in one sentence and is
worth repeating, because every decision below follows from it: **an embedding
network in the SBI layer receives one tensor** -- that is the whole of what the
``sbi`` package's interface allows -- so everything a network could condition on
must travel *inside* it. Where each sample sits, what was measured, how well,
whether it counts, and which dataset it came from are therefore **columns**, and
this module fixes the packing: the rows, the column groups, their order, their
standardisation, and the frozen description (the :class:`EncodingLayout`) whose
hash says whether two tensors mean the same thing.

It deliberately fixes **nothing about the network**. :func:`unpack` is the one
way a network is meant to read the tensor, and every embedding -- ``sbi``'s set
and transformer nets today, neural-process or operator encoders later -- is a
module over the unpacked view. That is what makes the packing the contract and
the networks a free choice.

Two consequences are the point of the whole exercise:

* a set-based network over these rows is invariant to how many samples there are
  and where they sit, which is what amortisation across differently-sampled
  datasets needs;
* uncertainties are columns, so a network can condition on the error bars of the
  observation at hand, which is what amortisation across noise realisations
  needs -- provided the simulator varied them (the reserved observation context;
  ``docs/design/horizon_notes.md`` §4).

Backend-neutral by construction: numpy in, numpy out, no torch and no jax. The
torch-side wrappers live in :mod:`ampere.inference._sbi`, and they consume
:func:`unpack`, which is written in plain slicing so that it works unchanged on a
``torch.Tensor``.

Examples
--------
>>> import astropy.units as u
>>> import numpy as np
>>> from ampere.core import Dataset, Spectrum
>>> from ampere.core.encoding import EncodingLayout, encode_observations, unpack
>>> observed = Spectrum(
...     [1.0, 2.0, 3.0] * u.um,
...     [1.0, 2.0, 3.0] * u.Jy,
...     uncertainty=[0.1, 0.1, 0.1] * u.Jy,
...     mask=[False, True, False],
... )
>>> datasets = {"sed": Dataset(observed, label="sed")}
>>> layout = EncodingLayout.from_datasets(datasets)
>>> layout.kind, layout.row_cap, len(layout.hash)
('set', 3, 32)
>>> encoded = encode_observations(datasets, layout=layout)
>>> encoded.values.shape == (1, layout.row_cap, layout.columns_total)
True
>>> unpack(encoded.values, layout).valid.tolist()
[[True, False, True]]
"""

from __future__ import annotations

import dataclasses
import math
import sys
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .exceptions import ContractError
from .results_schema import Layout as ContainerLayout

__all__ = [
    "COLUMN_GROUPS",
    "DEFAULT_FOURIER_BANDS",
    "ENCODING_VERSION",
    "FLAT_KIND",
    "SET_KIND",
    "ColumnGroup",
    "DatasetLayout",
    "DatasetView",
    "Decoded",
    "Encoded",
    "EncodingError",
    "EncodingLayout",
    "Unpacked",
    "decode",
    "encode",
    "encode_observations",
    "unpack",
]

#: Bumped when a column group is added, removed, reordered or redefined. The
#: reserved zero-width ``context`` group exists precisely so that the first real
#: observation context is a *width* change under version 1 rather than a version
#: bump (``encoding.md`` §8).
ENCODING_VERSION = 1

#: ``B`` in ``encoding.md`` §3 group 3: how many NeRF-style Fourier bands each
#: standardised coordinate is expanded into. ``0`` removes the group.
DEFAULT_FOURIER_BANDS = 4

#: The general packing this module exists for: one row per sample, the column
#: groups of §3, and a set-based network over them.
SET_KIND = "set"

#: W3.2's fixed-size summary — one row, every dataset's unmasked observed
#: values concatenated in ``datasets`` order. Kept as a degenerate layout under
#: the same object so that its hash covers it too, and it remains the right
#: default for a single fitting problem whose data layout is fixed anyway.
FLAT_KIND = "flat"

#: The column groups of a ``"set"`` layout, **in order**. Names are contract.
COLUMN_GROUPS: tuple[str, ...] = (
    "dataset",
    "coordinate",
    "coordinate_features",
    "value",
    "value_asinh",
    "log_sigma",
    "mask",
    "set_features",
    "context",
)

#: Width of the per-set feature group: ``log(N_valid)``, ``has_sigma``,
#: ``is_complex``.
_SET_FEATURE_WIDTH = 3


class EncodingError(ContractError):
    """The encoding refused: a layout mismatch, a row cap, or unencodable data.

    Declared here rather than in :mod:`ampere.core.exceptions` because the
    encoding is a **post-freeze** §4 addition (``encoding.md``) and its
    vocabulary should arrive and, if it ever had to, leave with it. It is a
    :class:`~ampere.core.exceptions.ContractError` like every other contract
    refusal, so a caller catching that catches this.
    """


# ---------------------------------------------------------------------------
# The layout
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class ColumnGroup:
    """One named, contiguous block of columns with one meaning."""

    name: str
    offset: int
    width: int

    @property
    def stop(self) -> int:
        """One past the group's last column."""
        return self.offset + self.width

    def to_dict(self) -> dict[str, Any]:
        """The JSON-safe record that goes into the layout's hash."""
        return {"name": self.name, "offset": self.offset, "width": self.width}


@dataclasses.dataclass(frozen=True)
class DatasetLayout:
    """Everything the packing needs to know about one dataset's observation.

    Every field is computed **once, from the observed container**, at layout
    construction, and never from a simulation (``encoding.md`` §4). That is what
    makes a layout a function of the *problem*: fixed before any draw, identical
    at training and at inference, independent of the budget.

    Attributes
    ----------
    label
        The dataset's label, in ``datasets`` order.
    kind
        The container class name (``"Spectrum"``, ``"VisibilitySet"``, ...).
    container_layout
        ``"points"`` or ``"grid"`` -- :class:`ampere.core.results_schema.Layout`
        by value, so the record stays plain data.
    axes
        The axis names, in the container kind's declared order.
    grid_shape
        The container's value shape for a grid kind, else ``None``. It is what
        lets :attr:`Unpacked.per_dataset` reshape a slice back to its axes
        exactly, and it is why rows are stored in C-order.
    rows
        Total samples, **masked ones included**: a masked sample is a row with
        ``mask = 0``, not an absent row.
    valid
        Samples the effective mask retains.
    excluded
        Flat indices (C-order) the effective mask excludes. Stored rather than
        re-read at every encode, so that the observation and every simulated
        draw carry byte-identical mask columns even though
        :attr:`~ampere.core.dataset.Dataset.effective_mask` is resolved lazily.
    is_complex, has_sigma
        Whether the values are complex, and whether the container carries
        uncertainties. Both are recorded rather than inferred per call.
    coordinate_ranges
        ``(min, max)`` of each axis's observed values, the map to ``[-1, 1]``.
    value_scale
        ``s``: the median of ``|y|`` over valid samples, falling back to the
        median sigma and then to 1.0. Everything dimensionless is divided by it.
    """

    label: str
    kind: str
    container_layout: str
    axes: tuple[str, ...]
    grid_shape: tuple[int, ...] | None
    rows: int
    valid: int
    excluded: tuple[int, ...]
    is_complex: bool
    has_sigma: bool
    coordinate_ranges: tuple[tuple[float, float], ...]
    value_scale: float

    @property
    def axis_count(self) -> int:
        """How many coordinate axes this dataset's container declares."""
        return len(self.axes)

    def mask_column(self) -> np.ndarray:
        """``1.0`` where a row counts, ``0.0`` where the effective mask excludes it."""
        column = np.ones(self.rows, dtype=float)
        if self.excluded:
            column[np.asarray(self.excluded, dtype=int)] = 0.0
        return column

    def to_dict(self) -> dict[str, Any]:
        """The JSON-safe record that goes into the layout's hash."""
        return {
            "label": self.label,
            "kind": self.kind,
            "container_layout": self.container_layout,
            "axes": list(self.axes),
            "grid_shape": None if self.grid_shape is None else list(self.grid_shape),
            "rows": self.rows,
            "valid": self.valid,
            "excluded": list(self.excluded),
            "is_complex": self.is_complex,
            "has_sigma": self.has_sigma,
            "coordinate_ranges": [list(pair) for pair in self.coordinate_ranges],
            "value_scale": self.value_scale,
        }


@dataclasses.dataclass(frozen=True)
class EncodingLayout:
    """The frozen, hashable description of one problem's packing.

    Two tensors with the same :attr:`hash` mean the same thing. A network
    trained under one hash **refuses** an observation encoded under another, by
    name, saying which field differs (:meth:`check_against`) -- which is the
    whole reason this object exists rather than a handful of loose constants.

    Build one with :meth:`from_datasets`; the constructor is for round-tripping
    a stored record.
    """

    kind: str
    version: int
    datasets: tuple[DatasetLayout, ...]
    coordinates: int
    fourier_bands: int
    row_cap: int
    columns: tuple[ColumnGroup, ...]
    complex_columns: bool
    hash: str = dataclasses.field(default="", compare=False, repr=False)

    def __post_init__(self) -> None:
        if self.kind not in (SET_KIND, FLAT_KIND):
            raise EncodingError(
                f"an encoding layout is {SET_KIND!r} or {FLAT_KIND!r}, not {self.kind!r}."
            )
        object.__setattr__(self, "datasets", tuple(self.datasets))
        object.__setattr__(self, "columns", tuple(self.columns))
        object.__setattr__(self, "hash", _hash_of(self.to_dict()))

    # -- construction ---------------------------------------------------------

    @classmethod
    def from_datasets(
        cls,
        datasets: Mapping[str, Any],
        *,
        kind: str = SET_KIND,
        fourier_bands: int = DEFAULT_FOURIER_BANDS,
        row_cap: int | None = None,
    ) -> EncodingLayout:
        """The layout of *datasets*' observed containers.

        Parameters
        ----------
        datasets
            The problem's own collection (a
            :class:`~ampere.core.dataset.DatasetCollection`, or any mapping of
            label to :class:`~ampere.core.dataset.Dataset`). **The order is read
            from here**, never from a mapping built elsewhere: a dict that
            iterated differently would silently permute the network's columns.
        kind
            ``"set"`` (the general packing) or ``"flat"`` (W3.2's summary).
        fourier_bands
            ``B``; ``0`` removes the coordinate-feature group.
        row_cap
            ``R``. Defaults to the observation's own total row count, so the
            default layout fits the observation exactly and a differently
            sampled observation needs a wider cap **chosen on purpose**. A cap
            below the observation's own rows is refused here rather than at the
            first encode.
        """
        if not datasets:
            raise EncodingError(
                "there is nothing to encode: this problem declares no datasets, so the "
                "encoding has no rows and no columns."
            )
        bands = int(fourier_bands)
        if bands < 0:
            raise EncodingError(f"fourier_bands counts Fourier bands and cannot be {bands}.")
        records = tuple(_dataset_layout(label, dataset) for label, dataset in _ordered(datasets))
        if not any(record.valid for record in records):
            raise EncodingError(
                "every sample of every dataset is masked out, so the encoding would be "
                "empty: there is nothing for a density estimator to condition on."
            )
        observed_rows = sum(record.rows for record in records)
        cap = observed_rows if row_cap is None else int(row_cap)
        if cap < observed_rows:
            raise EncodingError(
                f"this observation has {observed_rows} row(s) but the layout's row cap is "
                f"{cap}. Rows are samples, masked ones included, so the cap has to be at "
                f"least the sample count; raise row_cap, or mask what should not be fitted."
            )
        axis_count = max((record.axis_count for record in records), default=0)
        complex_columns = any(record.is_complex for record in records)
        if kind == FLAT_KIND:
            for record in records:
                if record.is_complex:
                    raise EncodingError(
                        f"dataset {record.label!r} holds complex values, and the {FLAT_KIND!r} "
                        f"summary layout has no encoding for one -- flattening a complex array "
                        f"into a real feature vector is a choice (real and imaginary parts as "
                        f"two columns) that the {SET_KIND!r} layout makes and this one does "
                        f"not. Use layout='set', or fit a complex dataset with a "
                        f"likelihood-based engine."
                    )
            width = sum(record.valid for record in records)
            columns: tuple[ColumnGroup, ...] = (ColumnGroup("flat", 0, width),)
            return cls(
                kind=FLAT_KIND,
                version=ENCODING_VERSION,
                datasets=records,
                coordinates=axis_count,
                fourier_bands=0,
                row_cap=1,
                columns=columns,
                complex_columns=complex_columns,
            )
        return cls(
            kind=SET_KIND,
            version=ENCODING_VERSION,
            datasets=records,
            coordinates=axis_count,
            fourier_bands=bands,
            row_cap=cap,
            columns=_column_groups(
                axis_count=axis_count, bands=bands, complex_columns=complex_columns
            ),
            complex_columns=complex_columns,
        )

    # -- introspection --------------------------------------------------------

    @property
    def name(self) -> str:
        """The layout's short name, as a run and a training set record it."""
        return self.kind

    @property
    def columns_total(self) -> int:
        """Total width of the encoded tensor's last axis."""
        return sum(group.width for group in self.columns)

    @property
    def has_sigma(self) -> bool:
        """Whether **every** dataset carries uncertainties.

        The sigma column is present whatever this says (``encoding.md`` §3 group
        6): what changes is whether it means anything, which is recorded per
        dataset in :attr:`DatasetLayout.has_sigma` and broadcast onto every row
        as the second per-set feature.
        """
        return all(record.has_sigma for record in self.datasets)

    @property
    def labels(self) -> tuple[str, ...]:
        """The dataset labels, in order."""
        return tuple(record.label for record in self.datasets)

    @property
    def bounds(self) -> tuple[tuple[int, int], ...]:
        """``(start, stop)`` row bounds of each dataset's contiguous block."""
        found: list[tuple[int, int]] = []
        offset = 0
        for record in self.datasets:
            found.append((offset, offset + record.rows))
            offset += record.rows
        return tuple(found)

    def group(self, name: str) -> ColumnGroup:
        """The named column group."""
        for candidate in self.columns:
            if candidate.name == name:
                return candidate
        known = ", ".join(group.name for group in self.columns)
        raise EncodingError(f"this layout has no column group {name!r}; it has {known}.")

    def to_dict(self) -> dict[str, Any]:
        """The JSON-safe record the hash is taken over, and the attrs carry."""
        return {
            "kind": self.kind,
            "version": self.version,
            "datasets": [record.to_dict() for record in self.datasets],
            "coordinates": self.coordinates,
            "fourier_bands": self.fourier_bands,
            "row_cap": self.row_cap,
            "columns": [group.to_dict() for group in self.columns],
            "complex_columns": self.complex_columns,
        }

    @classmethod
    def from_dict(cls, record: Mapping[str, Any]) -> EncodingLayout:
        """Rebuild a layout from :meth:`to_dict` -- a stored run's attrs, say."""
        return cls(
            kind=str(record["kind"]),
            version=int(record["version"]),
            datasets=tuple(
                DatasetLayout(
                    label=str(one["label"]),
                    kind=str(one["kind"]),
                    container_layout=str(one["container_layout"]),
                    axes=tuple(str(name) for name in one["axes"]),
                    grid_shape=(
                        None
                        if one["grid_shape"] is None
                        else tuple(int(size) for size in one["grid_shape"])
                    ),
                    rows=int(one["rows"]),
                    valid=int(one["valid"]),
                    excluded=tuple(int(index) for index in one["excluded"]),
                    is_complex=bool(one["is_complex"]),
                    has_sigma=bool(one["has_sigma"]),
                    coordinate_ranges=tuple(
                        (float(pair[0]), float(pair[1])) for pair in one["coordinate_ranges"]
                    ),
                    value_scale=float(one["value_scale"]),
                )
                for one in record["datasets"]
            ),
            coordinates=int(record["coordinates"]),
            fourier_bands=int(record["fourier_bands"]),
            row_cap=int(record["row_cap"]),
            columns=tuple(
                ColumnGroup(str(one["name"]), int(one["offset"]), int(one["width"]))
                for one in record["columns"]
            ),
            complex_columns=bool(record["complex_columns"]),
        )

    # -- refusals -------------------------------------------------------------

    def differences(self, datasets: Mapping[str, Any]) -> list[str]:
        """Which fields of *datasets* this layout does not describe.

        A **different observation of the same shape** produces no differences --
        that is the amortisation the whole encoding exists for. A different
        *shape* produces one line per field that differs, which is what
        :meth:`check_against` puts in its refusal.
        """
        observed = EncodingLayout.from_datasets(
            datasets,
            kind=self.kind,
            fourier_bands=self.fourier_bands,
            row_cap=max(self.row_cap, sum(_rows_of(dataset) for _, dataset in _ordered(datasets))),
        )
        return self.compare(observed)

    def compare(self, other: EncodingLayout) -> list[str]:
        """Which fields differ between this layout and *other*, in words."""
        found: list[str] = []
        if self.kind != other.kind:
            found.append(f"kind: {self.kind!r} here, {other.kind!r} there")
        if self.version != other.version:
            found.append(f"version: {self.version} here, {other.version} there")
        if self.labels != other.labels:
            found.append(f"dataset labels: {list(self.labels)} here, {list(other.labels)} there")
            return found
        for mine, theirs in zip(self.datasets, other.datasets, strict=True):
            if mine.excluded != theirs.excluded:
                found.append(
                    f"dataset {mine.label!r} masked samples: {len(mine.excluded)} here, "
                    f"{len(theirs.excluded)} there (the mask is frozen into the layout, so a "
                    f"differently masked observation is a different layout)"
                )
            for field in ("kind", "axes", "rows", "is_complex", "has_sigma", "grid_shape"):
                left = getattr(mine, field)
                right = getattr(theirs, field)
                if left != right:
                    found.append(f"dataset {mine.label!r} {field}: {left!r} here, {right!r} there")
        if self.coordinates != other.coordinates:
            found.append(f"coordinates: {self.coordinates} here, {other.coordinates} there")
        if self.fourier_bands != other.fourier_bands:
            found.append(f"fourier_bands: {self.fourier_bands} here, {other.fourier_bands} there")
        return found

    def check_against(self, datasets: Mapping[str, Any]) -> None:
        """Refuse *datasets* by name if this layout does not describe them."""
        differences = self.differences(datasets)
        if differences:
            joined = "; ".join(differences)
            raise EncodingError(
                f"the encoding layout {self.name!r} (hash {self.hash}) does not describe this "
                f"problem: {joined}. A different observation of the same shape is fine -- that "
                f"is what amortisation means -- but a different shape needs its own layout, "
                f"and therefore its own trained network."
            )


# ---------------------------------------------------------------------------
# The encoded tensor
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Encoded:
    """One or more observations packed under one layout.

    ``values`` is always three-dimensional -- ``(count, rows, columns)``, with
    ``count = 1`` for a single observation -- so that the observation the
    posterior is conditioned on and the rows a network trained on have the same
    rank and the same code path. A ``"flat"`` layout is the degenerate case:
    ``rows == 1``, and :attr:`matrix` is the ``(count, features)`` summary W3.2
    trains on.
    """

    values: np.ndarray
    layout: EncodingLayout

    @property
    def count(self) -> int:
        """How many observations this holds."""
        return int(self.values.shape[0])

    @property
    def matrix(self) -> np.ndarray:
        """``(count, rows * columns)``: the tensor flattened per observation."""
        return np.asarray(self.values).reshape(self.count, -1)

    def __len__(self) -> int:
        return self.count


@dataclasses.dataclass(frozen=True)
class Decoded:
    """One dataset read back out of an encoded tensor (:func:`decode`).

    The round trip is exact for the coordinates, the mask and -- to floating
    point -- the values and, where they are strictly positive, the
    uncertainties. It is *not* a way to recover what the encoding does not
    carry: a padded row has no dataset, and a sample whose sigma is zero has no
    ``log sigma``.
    """

    label: str
    coordinates: dict[str, np.ndarray]
    values: np.ndarray
    uncertainty: np.ndarray | None
    mask: np.ndarray | None

    def as_container(self, template: Any) -> Any:
        """Rebuild a container from *template*'s axes and this dataset's arrays.

        :meth:`~ampere.core.results_schema.FunctionSamples.with_values` is the
        constructor used, so the axes, the unit, the extra coordinates and the
        metadata are the template's -- none of which the encoding carries, and
        none of which a decode could honestly invent.
        """
        return template.with_values(self.values, uncertainty=self.uncertainty, mask=self.mask)


# ---------------------------------------------------------------------------
# unpack: the one way a network reads the tensor
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class DatasetView:
    """One dataset's rows of an unpacked tensor (:attr:`Unpacked.per_dataset`).

    Every view is a **slice**, so this costs no copy on numpy and no graph node
    on torch. ``grid_shape`` is what lets a first-stage encoder reshape an
    ``Image`` or ``Cube`` slice back to its axes exactly -- rows are stored in
    C-order for precisely this reason -- so a CNN or ConvCNP encoder over a grid
    kind is a wrapper rather than a packing change.
    """

    label: str
    index: int
    start: int
    stop: int
    grid_shape: tuple[int, ...] | None
    dataset: Any
    coordinate: Any
    coordinate_features: Any
    value: Any
    value_asinh: Any
    log_sigma: Any
    mask: Any
    set_features: Any
    context: Any

    @property
    def features(self) -> Any:
        """Every non-mask group of this dataset's rows, concatenated."""
        return _concatenate(
            [
                self.dataset,
                self.coordinate,
                self.coordinate_features,
                self.value,
                self.value_asinh,
                self.log_sigma,
                self.set_features,
                self.context,
            ]
        )

    @property
    def valid(self) -> Any:
        """``(..., rows)`` boolean: which of this dataset's rows count."""
        return self.mask[..., 0] > 0.5


@dataclasses.dataclass(frozen=True)
class Unpacked:
    """Named views of an encoded tensor's column groups.

    **This is the only way a network is meant to read the tensor.** It is
    written in plain slicing and one concatenation, so the same function serves
    a numpy array in a test and a ``torch.Tensor`` inside an embedding net's
    ``forward``.
    """

    layout: EncodingLayout
    dataset: Any
    coordinate: Any
    coordinate_features: Any
    value: Any
    value_asinh: Any
    log_sigma: Any
    mask: Any
    set_features: Any
    context: Any
    per_dataset: tuple[DatasetView, ...]

    @property
    def features(self) -> Any:
        """Every non-mask group, concatenated: what a per-row network consumes."""
        return _concatenate(
            [
                self.dataset,
                self.coordinate,
                self.coordinate_features,
                self.value,
                self.value_asinh,
                self.log_sigma,
                self.set_features,
                self.context,
            ]
        )

    @property
    def valid(self) -> Any:
        """``(..., rows)`` boolean: ``True`` where a row counts."""
        return self.mask[..., 0] > 0.5


def unpack(x: Any, layout: EncodingLayout) -> Unpacked:
    """Named views of *x*'s column groups, plus the per-dataset slices.

    Parameters
    ----------
    x
        ``(..., rows, columns)`` under *layout*. A ``(..., rows * columns)``
        tensor is accepted and reshaped, because ``sbi`` flattens ``x`` on some
        paths and a wrapper must not have to care which one it is on.
    layout
        The layout *x* was encoded under. A tensor whose width is neither
        ``columns`` nor ``rows * columns`` is refused by name.
    """
    if layout.kind != SET_KIND:
        raise EncodingError(
            f"unpack() reads the {SET_KIND!r} packing's column groups, and this layout is "
            f"{layout.kind!r}, which has none: it is one row of concatenated values. Build a "
            f"{SET_KIND!r} layout if a network needs coordinates, sigma and masks."
        )
    columns = layout.columns_total
    rows = layout.row_cap
    width = int(x.shape[-1])
    if width == rows * columns and width != columns:
        x = x.reshape(*tuple(x.shape[:-1]), rows, columns)
    elif width != columns:
        raise EncodingError(
            f"this tensor is {width} wide, and the layout {layout.name!r} (hash "
            f"{layout.hash}) packs {columns} column(s) over {rows} row(s). It was encoded "
            f"under a different layout; a network trained on one layout must not be believed "
            f"about another."
        )
    if int(x.shape[-2]) != rows:
        raise EncodingError(
            f"this tensor has {int(x.shape[-2])} row(s) and the layout {layout.name!r} (hash "
            f"{layout.hash}) packs {rows}."
        )
    views = {group.name: x[..., group.offset : group.stop] for group in layout.columns}
    per_dataset: list[DatasetView] = []
    for index, record in enumerate(layout.datasets):
        first, last = layout.bounds[index]
        per_dataset.append(
            DatasetView(
                label=record.label,
                index=index,
                start=first,
                stop=last,
                grid_shape=record.grid_shape,
                **{
                    group.name: x[..., first:last, group.offset : group.stop]
                    for group in layout.columns
                },
            )
        )
    return Unpacked(layout=layout, per_dataset=tuple(per_dataset), **views)


# ---------------------------------------------------------------------------
# encode
# ---------------------------------------------------------------------------


def encode_observations(datasets: Mapping[str, Any], *, layout: EncodingLayout) -> Encoded:
    """Encode a collection's **observed** containers under *layout*.

    The observation the trained posterior is conditioned on. Equivalent to
    :func:`encode` over ``{label: dataset.observed}``, and spelled separately
    because that is the call every driver makes first.
    """
    return encode(
        {label: dataset.observed for label, dataset in _ordered(datasets)},
        layout=layout,
        batched=False,
    )


def encode(
    observations: Mapping[str, Any] | None,
    *,
    layout: EncodingLayout,
    batched: bool = False,
) -> Encoded:
    """Pack *observations* into one tensor under *layout*.

    Parameters
    ----------
    observations
        Label to container. With *batched* true these are
        :class:`~ampere.core.simulate.ContainerBatch`\\ es, whose leading axis is
        the sample axis; with it false they are single
        :class:`~ampere.core.results_schema.FunctionSamples`. A label the layout
        knows and this mapping lacks is refused by name -- a simulation that
        produced no observation is not an empty one.
    layout
        The frozen packing. Its statistics, its masks and its row bounds are
        used **as they stand**: nothing here is recomputed from *observations*,
        which is what makes training and inference identical by construction.
    batched
        Whether *observations* carries a leading sample axis.
    """
    sources = {} if observations is None else observations
    pieces = [
        _dataset_source(record, sources.get(record.label), batched=batched)
        for record in layout.datasets
    ]
    count = _one_count([piece[0] for piece in pieces], labels=layout.labels)
    if layout.kind == FLAT_KIND:
        return Encoded(values=_flat_tensor(layout, pieces, count), layout=layout)
    return Encoded(values=_set_tensor(layout, pieces, count), layout=layout)


def decode(encoded: Encoded, *, index: int = 0) -> dict[str, Decoded]:
    """Read one observation back out of an encoded tensor.

    The inverse of :func:`encode` for everything the packing carries: the
    coordinates come back through the layout's ranges, the values through
    ``s * sinh(value_asinh)`` -- the ``asinh`` column rather than the whitened
    one, because it is defined whether or not the dataset has uncertainties --
    and the uncertainties through ``s * exp(log_sigma)``. It exists for the
    round trip a test asserts and for a reader inspecting a stored tensor; it is
    not part of any network's path.
    """
    layout = encoded.layout
    if layout.kind != SET_KIND:
        raise EncodingError(
            f"decode() inverts the {SET_KIND!r} packing; the {FLAT_KIND!r} layout drops the "
            f"coordinates, the uncertainties and the masked samples, so there is nothing to "
            f"invert."
        )
    tensor = np.asarray(encoded.values, dtype=float)[int(index)]
    value_group = layout.group("value_asinh")
    sigma_group = layout.group("log_sigma")
    coordinate_group = layout.group("coordinate")
    mask_group = layout.group("mask")
    out: dict[str, Decoded] = {}
    for record, (start, stop) in zip(layout.datasets, layout.bounds, strict=True):
        block = tensor[start:stop]
        scale = record.value_scale
        real = scale * np.sinh(block[:, value_group.offset])
        if record.is_complex:
            imaginary = scale * np.sinh(block[:, value_group.offset + 1])
            values: np.ndarray = real + 1j * imaginary
        else:
            values = real
        uncertainty = scale * np.exp(block[:, sigma_group.offset]) if record.has_sigma else None
        keep = block[:, mask_group.offset] > 0.5
        mask = None if bool(np.all(keep)) else ~keep
        shape = record.grid_shape if record.grid_shape is not None else (record.rows,)
        coordinates: dict[str, np.ndarray] = {}
        for axis_index, name in enumerate(record.axes):
            column = block[:, coordinate_group.offset + axis_index]
            low, high = record.coordinate_ranges[axis_index]
            span = high - low
            raw = low + 0.5 * (column + 1.0) * span if span != 0.0 else np.full(column.shape, low)
            if record.container_layout == ContainerLayout.GRID.value:
                grid = raw.reshape(shape)
                picker: list[Any] = [0] * len(shape)
                picker[axis_index] = slice(None)
                coordinates[name] = np.asarray(grid[tuple(picker)])
            else:
                coordinates[name] = raw
        out[record.label] = Decoded(
            label=record.label,
            coordinates=coordinates,
            values=values.reshape(shape),
            uncertainty=None if uncertainty is None else uncertainty.reshape(shape),
            mask=None if mask is None else mask.reshape(shape),
        )
    return out


# ---------------------------------------------------------------------------
# The packing itself
# ---------------------------------------------------------------------------


def _set_tensor(
    layout: EncodingLayout,
    pieces: Sequence[tuple[int, np.ndarray, np.ndarray | None, Any]],
    count: int,
) -> np.ndarray:
    """``(count, row_cap, columns)``, built one dataset's block at a time."""
    tensor = np.zeros((count, layout.row_cap, layout.columns_total), dtype=float)
    groups = {group.name: group for group in layout.columns}
    for index, (record, (start, stop)) in enumerate(
        zip(layout.datasets, layout.bounds, strict=True)
    ):
        _, values, sigma, axes = pieces[index]
        block = tensor[:, start:stop, :]
        keep = record.mask_column()

        block[:, :, groups["dataset"].offset] = float(index)
        block[:, :, groups["mask"].offset] = keep

        coordinates = _standardised_coordinates(record, axes)
        group = groups["coordinate"]
        block[:, :, group.offset : group.offset + record.axis_count] = coordinates
        if layout.fourier_bands:
            group = groups["coordinate_features"]
            features = _fourier_features(coordinates, layout.fourier_bands)
            block[:, :, group.offset : group.offset + features.shape[1]] = features

        group = groups["set_features"]
        block[:, :, group.offset] = math.log(max(record.valid, 1))
        block[:, :, group.offset + 1] = 1.0 if record.has_sigma else 0.0
        block[:, :, group.offset + 2] = 1.0 if record.is_complex else 0.0

        scale = record.value_scale
        parts = (values.real, values.imag) if record.is_complex else (values,)
        stacked = np.stack(parts, axis=-1)
        # Masked rows are allowed to hold anything, including a NaN or a zero
        # sigma -- masking is what a user does *about* those -- so the arithmetic
        # below is expected to produce non-finite entries on those rows and
        # :func:`_sanitised` zeroes them. Silencing numpy here rather than
        # letting a warning escape is deliberate: the retained rows are checked
        # explicitly, by name, a few lines down.
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            asinh = np.arcsinh(stacked / scale)
            if record.has_sigma and sigma is not None:
                _check_sigma(record, sigma, keep)
                whitened = stacked / sigma[..., None]
                log_sigma = np.log(sigma) - math.log(scale)
            else:
                whitened = asinh
                log_sigma = np.zeros(values.shape, dtype=float)
        _check_finite(record, asinh, keep)

        group = groups["value"]
        block[:, :, group.offset : group.offset + asinh.shape[-1]] = _sanitised(whitened)
        group = groups["value_asinh"]
        block[:, :, group.offset : group.offset + asinh.shape[-1]] = _sanitised(asinh)
        block[:, :, groups["log_sigma"].offset] = _sanitised(log_sigma[..., None])[..., 0]
    return tensor


def _flat_tensor(
    layout: EncodingLayout,
    pieces: Sequence[tuple[int, np.ndarray, np.ndarray | None, Any]],
    count: int,
) -> np.ndarray:
    """``(count, 1, features)``: W3.2's summary, as a degenerate layout.

    Each dataset's observed values, flattened, with the samples its effective
    mask excludes **dropped** -- the same columns absent from every draw and
    from the observation the posterior is conditioned on, which is the whole of
    what "fixed layout" means.
    """
    columns: list[np.ndarray] = []
    for index, record in enumerate(layout.datasets):
        _, values, _, _ = pieces[index]
        keep = record.mask_column() > 0.5
        columns.append(np.asarray(values, dtype=float)[:, keep])
    summary = np.concatenate(columns, axis=1) if columns else np.zeros((count, 0))
    if summary.shape[1] == 0:
        raise EncodingError(
            "the summary vector is empty: every sample of every dataset is masked out, so "
            "there is nothing for the density estimator to condition on."
        )
    return summary.reshape(count, 1, summary.shape[1])


def _standardised_coordinates(record: DatasetLayout, axes: Any) -> np.ndarray:
    """``(rows, axis_count)`` in ``[-1, 1]``, per dataset and per axis."""
    raw = _coordinate_matrix(record, axes)
    out = np.zeros_like(raw)
    for index in range(record.axis_count):
        low, high = record.coordinate_ranges[index]
        span = high - low
        if span == 0.0:
            continue
        out[:, index] = 2.0 * (raw[:, index] - low) / span - 1.0
    return out


def _coordinate_matrix(record: DatasetLayout, axes: Any) -> np.ndarray:
    """``(rows, axis_count)`` of raw coordinates, in the container's C-order."""
    if record.axis_count == 0:
        return np.zeros((record.rows, 0), dtype=float)
    if record.container_layout == ContainerLayout.GRID.value:
        shape = record.grid_shape or ()
        indices = np.indices(shape)
        return np.column_stack(
            [
                np.asarray(axis.values, dtype=float)[indices[position]].reshape(-1)
                for position, axis in enumerate(axes)
            ]
        )
    return np.column_stack([np.asarray(axis.values, dtype=float).reshape(-1) for axis in axes])


def _fourier_features(coordinates: np.ndarray, bands: int) -> np.ndarray:
    """``(rows, axis_count * 2 * bands)``: ``sin`` and ``cos`` at ``pi * 2**k``.

    Axis-major, then band, then ``(sin, cos)``. An axis the layout pads (a
    dataset with fewer axes than the widest) contributes zeros rather than
    ``cos(0) = 1``: an axis that does not exist should not look like an axis
    sitting at the centre of its range.
    """
    rows, axis_count = coordinates.shape
    out = np.zeros((rows, axis_count * 2 * bands), dtype=float)
    for axis in range(axis_count):
        for band in range(bands):
            angle = math.pi * (2.0**band) * coordinates[:, axis]
            base = axis * 2 * bands + 2 * band
            out[:, base] = np.sin(angle)
            out[:, base + 1] = np.cos(angle)
    return out


def _sanitised(block: np.ndarray) -> np.ndarray:
    """Zero every non-finite entry, and every entry of an excluded row.

    ``encoding.md`` §3 leaves a masked sample's value as the container holds it
    and lets the mask remove it downstream; **amended by W3.3**, because a
    masked sample may legitimately be a NaN or an infinity in the observed
    container (masking it is what a user does *about* that), and a non-finite
    number in the stored tensor poisons the training loss, the netCDF round trip
    and ``sbi``'s own shape validation even where the mask says the row is
    absent. Finite values on excluded rows are kept, so the round trip is exact
    there too.
    """
    out = np.nan_to_num(np.asarray(block, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    return out


def _check_sigma(record: DatasetLayout, sigma: np.ndarray, keep: np.ndarray) -> None:
    """A valid row's sigma must be strictly positive and finite."""
    valid = np.asarray(keep, dtype=bool)
    offending = sigma[:, valid]
    bad = int(np.count_nonzero(~(np.isfinite(offending) & (offending > 0.0))))
    if bad:
        raise EncodingError(
            f"dataset {record.label!r}: {bad} retained sample(s) have an uncertainty that is "
            f"not strictly positive and finite, and the encoding's whitened value (y/sigma) "
            f"and log-sigma columns cannot represent one. Mask those samples, or give the "
            f"dataset no uncertainties at all -- the layout records that as has_sigma=False "
            f"and the asinh-scaled value stands in."
        )


def _check_finite(record: DatasetLayout, asinh: np.ndarray, keep: np.ndarray) -> None:
    """A valid row's value must be finite; a masked one may be anything."""
    valid = np.asarray(keep, dtype=bool)
    bad = int(np.count_nonzero(~np.isfinite(asinh[:, valid])))
    if bad:
        raise EncodingError(
            f"dataset {record.label!r}: {bad} retained sample(s) hold a non-finite value, "
            f"which no column of the encoding can carry. Mask them -- a masked sample's value "
            f"may be anything, because the mask column removes it."
        )


# ---------------------------------------------------------------------------
# Reading the containers
# ---------------------------------------------------------------------------


def _ordered(datasets: Mapping[str, Any]) -> list[tuple[str, Any]]:
    """``(label, dataset)`` in the collection's own order."""
    return [(str(label), datasets[label]) for label in datasets]


def _rows_of(dataset: Any) -> int:
    """The observed container's total sample count."""
    return int(np.asarray(dataset.observed.values).size)


def _dataset_layout(label: str, dataset: Any) -> DatasetLayout:
    """One dataset's frozen record, from its **observed** container alone."""
    observed = dataset.observed
    values = np.asarray(observed.values)
    excluded = _excluded_of(dataset, values.size)
    keep = np.ones(values.size, dtype=bool)
    keep[excluded] = False
    flat = values.reshape(-1)
    magnitude = np.abs(flat[keep])
    sigma = None if observed.uncertainty is None else np.asarray(observed.uncertainty).reshape(-1)
    scale = _value_scale(magnitude, None if sigma is None else sigma[keep])
    container_layout = getattr(type(observed), "LAYOUT", ContainerLayout.POINTS)
    is_grid = container_layout is ContainerLayout.GRID
    return DatasetLayout(
        label=label,
        kind=type(observed).__name__,
        container_layout=container_layout.value,
        axes=tuple(axis.name for axis in observed.axes),
        grid_shape=tuple(int(size) for size in values.shape) if is_grid else None,
        rows=int(values.size),
        valid=int(np.count_nonzero(keep)),
        excluded=tuple(int(index) for index in np.flatnonzero(~keep)),
        is_complex=bool(np.iscomplexobj(values)),
        has_sigma=observed.uncertainty is not None,
        coordinate_ranges=tuple(
            (float(np.min(axis.values)), float(np.max(axis.values))) for axis in observed.axes
        ),
        value_scale=scale,
    )


def _excluded_of(dataset: Any, size: int) -> np.ndarray:
    """Flat indices the effective mask excludes.

    :attr:`~ampere.core.dataset.Dataset.effective_mask` is the union of the
    observed and predicted containers' masks and is the right answer, but it is
    resolved lazily, at the first prediction. Where it has not been resolved yet
    the observed container's own mask is the honest fallback -- it is what the
    union is built from -- and either way the answer is frozen into the layout,
    so the observation and every simulated draw carry the same mask column.
    """
    mask = getattr(dataset, "effective_mask", None)
    if mask is None:
        observed = dataset.observed
        mask = None if observed.mask is None else np.asarray(observed.mask)
    if mask is None:
        return np.zeros(0, dtype=int)
    return np.flatnonzero(np.asarray(mask, dtype=bool).reshape(-1)[:size])


def _value_scale(magnitude: np.ndarray, sigma: np.ndarray | None) -> float:
    """``s``: the median of ``|y|``, then the median sigma, then 1.0."""
    for candidate in (magnitude, sigma):
        if candidate is None or candidate.size == 0:
            continue
        finite = candidate[np.isfinite(candidate)]
        if finite.size == 0:
            continue
        scale = float(np.median(np.abs(finite)))
        if math.isfinite(scale) and scale > 0.0:
            return scale
    return 1.0


def _dataset_source(
    record: DatasetLayout, source: Any, *, batched: bool
) -> tuple[int, np.ndarray, np.ndarray | None, Any]:
    """``(count, values (count, rows), sigma (count, rows) | None, axes)``."""
    if source is None:
        raise EncodingError(
            f"the encoding has no observation for dataset {record.label!r}: the simulation "
            f"produced no observation for it. simulate_many(..., observe=True) is what draws "
            f"them, and a dataset whose likelihood family cannot sample refuses by name "
            f"(inference.md §13)."
        )
    template = source.template if batched else source
    values = np.asarray(source.values)
    uncertainty = source.uncertainty
    if batched:
        count = int(values.shape[0])
        flat = values.reshape(count, -1)
        sigma = (
            None if uncertainty is None else np.asarray(uncertainty, dtype=float).reshape(count, -1)
        )
    else:
        count = 1
        flat = values.reshape(1, -1)
        sigma = None if uncertainty is None else np.asarray(uncertainty, dtype=float).reshape(1, -1)
    if flat.shape[1] != record.rows:
        raise EncodingError(
            f"dataset {record.label!r}: the layout packs {record.rows} row(s) and this "
            f"container has {flat.shape[1]}. The row count is part of the layout, so a "
            f"differently sampled observation needs its own layout and its own trained "
            f"network."
        )
    if record.has_sigma and sigma is None:
        raise EncodingError(
            f"dataset {record.label!r}: the layout records has_sigma=True, so the encoding's "
            f"whitened value and log-sigma columns expect uncertainties, and this container "
            f"carries none."
        )
    return count, flat, sigma, template.axes


def _one_count(counts: Sequence[int], *, labels: Sequence[str]) -> int:
    """Every dataset must contribute the same number of draws."""
    distinct = sorted(set(counts))
    if len(distinct) > 1:
        pairs = ", ".join(f"{label}={count}" for label, count in zip(labels, counts, strict=True))
        raise EncodingError(
            f"the datasets of one batch must carry the same number of draws, and these carry "
            f"{pairs}."
        )
    return distinct[0] if distinct else 0


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------


def _column_groups(
    *, axis_count: int, bands: int, complex_columns: bool
) -> tuple[ColumnGroup, ...]:
    """The nine groups of ``encoding.md`` §3, in order, with their widths."""
    value_width = 2 if complex_columns else 1
    widths = {
        "dataset": 1,
        "coordinate": axis_count,
        "coordinate_features": axis_count * 2 * bands,
        "value": value_width,
        "value_asinh": value_width,
        "log_sigma": 1,
        "mask": 1,
        "set_features": _SET_FEATURE_WIDTH,
        "context": 0,
    }
    groups: list[ColumnGroup] = []
    offset = 0
    for name in COLUMN_GROUPS:
        groups.append(ColumnGroup(name, offset, widths[name]))
        offset += widths[name]
    return tuple(groups)


def _concatenate(parts: Sequence[Any]) -> Any:
    """Join along the last axis, on whichever array library *parts* came from.

    The one place this module has to know that a tensor might not be numpy.
    Nothing is imported: a torch tensor can only exist if torch is already in
    ``sys.modules``, so the module is looked up rather than imported, and
    :mod:`ampere.core` keeps importing no backend.
    """
    kept = [part for part in parts if int(part.shape[-1]) > 0]
    if len(kept) == 1:
        return kept[0]
    library = type(kept[0]).__module__.partition(".")[0]
    if library == "torch":
        return sys.modules["torch"].cat(kept, dim=-1)
    if library in ("jax", "jaxlib"):
        return sys.modules["jax"].numpy.concatenate(kept, axis=-1)
    return np.concatenate(kept, axis=-1)


def _hash_of(record: Mapping[str, Any]) -> str:
    """:func:`ampere.results.provenance.hash_of`, imported on use.

    The import is function-local and deliberately so: :mod:`ampere.core` imports
    no other ampere namespace at module scope, and the layout's hash has to be
    *the same recipe* every other ampere fingerprint uses -- canonical JSON, a
    BLAKE2b digest with ampere's personalisation string -- so that a layout hash
    sits beside a spec hash in a run's attrs and means the same kind of thing.
    """
    from ampere.results.provenance import hash_of

    return hash_of(record)
