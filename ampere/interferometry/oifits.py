"""The OIFITS reader: one file into the shipped interferometric containers (W6.12).

OIFITS is the exchange format of optical and infrared interferometry —
revision 1 (Pauls et al. 2005, PASP 117, 1255) and revision 2 (Duvert, Young
& Hummel 2017, A&A 597, A8). This module reads the three data tables a fit
uses into ampere's containers and nothing else: ``OI_VIS2`` (squared
visibilities), ``OI_T3`` (closure phases) and ``OI_VIS`` (complex
visibilities), each with the ``OI_WAVELENGTH`` table of its instrument, the
``OI_ARRAY`` table that names its stations and the ``OI_TARGET`` table that
names its source. ``OI_FLUX`` and ``OI_INSPOL`` (revision 2), and ``T3AMP``,
are not read; they are recorded here so that their absence is a decision.

The reader is :func:`read_oifits`; it uses :mod:`astropy.io.fits` and no
other dependency. Its docstring states the conventions it enforces, which are
the containers' (``ampere/core/results_schema.py``) and the reference
``FourierSample``/``ClosurePhase`` steps' — the reader's whole job is to
translate the file's conventions into those, once, at composition time.
"""

from __future__ import annotations

import dataclasses
import os
import warnings
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits

from ampere.core import ClosurePhases, VisibilitySet

__all__ = ["OIFITSData", "OIFITSError", "read_oifits"]

#: The three data tables this reader turns into containers.
DATA_TABLES = ("OI_VIS", "OI_VIS2", "OI_T3")

#: The unit OIFITS mandates for ``EFF_WAVE``, ``EFF_BAND`` and the ``(u, v)``
#: columns, used when a file omits ``TUNITn``.
_METRES = u.m

#: The unit the containers' spectral axis is declared in: micron, the unit every
#: shipped ``FourierSample`` (reference, torch, jax: their ``SPECTRAL_UNIT``)
#: stores its wavelength buffer in and emits its predicted container's axis in.
#: The likelihood's alignment check compares the predicted and observed axes
#: **exactly**, unit included, and a container keeps the unit it is given, so
#: an axis in metres — however correct — would never align with a prediction.
SPECTRAL_UNIT = u.micron


class OIFITSError(ValueError):
    """An OIFITS file this reader cannot turn into containers, named with the reason."""


@dataclasses.dataclass(frozen=True)
class OIFITSData:
    """What :func:`read_oifits` returns: one container per observable the file holds.

    A table the file does not have is ``None``, never an empty container.

    Attributes
    ----------
    visibilities
        From ``OI_VIS``: complex visibilities ``VISAMP * exp(i VISPHI)``.
    squared_visibilities
        From ``OI_VIS2``: ``VIS2DATA`` as measured, real and unitless; fitted
        with :class:`~ampere.backends.reference.SquaredAmplitude` on the
        prediction.
    closure_phases
        From ``OI_T3``: ``T3PHI`` in radians, in the canonical triangle order.
    target, instrument, array
        The ``OI_TARGET`` name, the ``INSNAME`` and the ``ARRNAME`` (``None``
        when the file names no array) the containers were selected by.
    revision
        The OIFITS revision (``OI_REVN``) of the data tables read.
    wavelengths, bandwidths
        The instrument's channels: ``EFF_WAVE`` and ``EFF_BAND``, in metres.
    """

    visibilities: VisibilitySet | None
    squared_visibilities: VisibilitySet | None
    closure_phases: ClosurePhases | None
    target: str
    instrument: str
    array: str | None
    revision: int
    wavelengths: u.Quantity
    bandwidths: u.Quantity


def read_oifits(
    path: str | os.PathLike[str],
    *,
    target: str | None = None,
    insname: str | None = None,
    arrname: str | None = None,
) -> OIFITSData:
    """Read an OIFITS file into ampere's ``VisibilitySet`` and ``ClosurePhases`` containers.

    **Selection.** A container holds one target seen by one instrument on one
    array. A file with one of each needs no keywords; a file with several
    targets, instruments (``INSNAME``) or arrays (``ARRNAME``) among its data
    tables is **refused** unless the corresponding keyword says which, and the
    refusal lists the choices — the reader never picks the first one for you.
    Several data tables of one kind that match the selection (two nights in
    two ``OI_VIS2`` tables, say) are concatenated, in file order.

    **Samples.** One container sample per table row and spectral channel,
    row-major (all channels of row 0, then row 1, ...). The coordinates are
    ``u = UCOORD / EFF_WAVE[c]`` and ``v = VCOORD / EFF_WAVE[c]``,
    dimensionless (baseline in wavelengths), and the spectral axis is
    ``EFF_WAVE[c]`` in **micron** (:data:`SPECTRAL_UNIT`): the unit the shipped
    ``FourierSample`` steps emit their prediction's axis in, which the
    likelihood's exact axis comparison requires the data to share.
    ``extra_coords`` carry, per sample, the ``channel`` index, the ``mjd``
    (days), the channel's ``eff_band`` (metres; ``extra_coords`` are unitless
    labels, ``results_schema.md`` §15.4), and the station labels: a
    ``baseline`` such as ``"S1-S2"`` on a visibility, and on a closure phase
    the ``triangle`` (``"S1-S2-E1"``) with the three baselines it was formed
    from as ``baseline_ij``, ``baseline_jk`` and ``baseline_ki``. Stations are
    named from ``OI_ARRAY``'s ``STA_NAME``, or by their ``STA_INDEX`` as text
    when the file has no ``OI_ARRAY`` for the array.

    **Masks.** ``FLAG`` becomes the mask directly — OIFITS ``True`` is bad,
    and ``results_schema.md`` §7's ``True`` is excluded. A sample whose value
    or uncertainty is not finite, or whose uncertainty is not positive (an
    infinite weight), is masked too, and its stored value and uncertainty
    replaced by the placeholders ``0`` and ``1`` (the container refuses
    negative uncertainties, and a masked sample is not data); each container's
    ``meta`` counts ``n_flagged``, ``n_nonfinite`` and ``n_nonpositive_error``.

    **Baselines: the canonical order.** Every baseline is stored as ``b_ij =
    r_j - r_i`` with ``i < j`` by ``STA_INDEX`` — the convention of
    :class:`~ampere.core.ClosurePhases` and of OIFITS itself, whose
    ``(UCOORD, VCOORD)`` for ``STA_INDEX = (a, b)`` is ``r_b - r_a``. A row
    listed as ``(b, a)`` with ``b > a`` has its ``(u, v)`` negated (and, on
    ``OI_VIS``, its phase, since ``V(-b) = conj V(b)``); on ``OI_VIS2`` the
    flip changes no value.

    **Closure phases: the canonical triangle.** ``OI_T3`` row ``STA_INDEX =
    (a, b, c)`` stores ``(U1, V1) = b_ab`` and ``(U2, V2) = b_bc``, the third
    baseline ``b_ca = -(b_ab + b_bc)``, and ``T3PHI = arg(V_ab V_bc V_ca)``
    in degrees. The container wants stations ``i < j < k`` with ``(u1, v1) =
    b_ij``, ``(u2, v2) = b_jk`` and the value ``arg(V_ij V_jk V_ki)``, which
    is exactly what :class:`~ampere.backends.reference.ClosurePhase` computes
    from the three baselines :meth:`~ampere.backends.reference.FourierSample.from_observed`
    lays out — the convention this reader must match. So:

    * ``b_ij`` and ``b_jk`` are taken from the three file baselines, each
      negated when the file traverses it the other way (``b_xy = -b_yx``);
    * the phase: the loop ``a -> b -> c -> a`` is the loop ``i -> j -> k ->
      i`` traversed either the same way or backwards. A **cyclic** shift of
      ``(a, b, c)`` is the same loop, so the bispectrum and its phase are
      unchanged; a **transposition** reverses the loop, which conjugates the
      bispectrum (``V_yx = conj V_xy`` for a real sky) and negates the phase.
      The phase is therefore negated exactly when the permutation ``(a, b, c)
      -> (i, j, k)`` is **odd**, kept when it is even, and then wrapped into
      ``(-180, 180]`` degrees before conversion to radians, i.e. ``(-pi, pi]``;
    * ``T3PHIERR`` is converted from degrees to radians the same way.

    **What is not converted.** ``VIS2DATA`` stays a squared visibility: the
    container holds what was measured, with ``meta["observable"] =
    "squared_visibility"``, and the prediction is squared to meet it (the
    ``SquaredAmplitude`` step, with ``normalisation`` set to the model's fixed
    total flux, since ``VIS2DATA`` is normalised to the zero spacing). Rooting
    the data instead would need an error propagation that is not symmetric
    near zero.

    **OI_VIS.** OIFITS gives an amplitude error ``VISAMPERR`` and a phase error
    ``VISPHIERR`` separately; the complex container has one circular
    per-component sigma, and the reader takes ``VISAMPERR`` — the symmetric
    approximation, good where ``VISAMP * VISPHIERR`` (radians) is close to
    ``VISAMPERR``. A fit that needs the full model fits the amplitudes with
    ``Amplitude`` and :class:`~ampere.core.RiceFamily` instead. A revision-2
    table whose ``AMPTYP`` or ``PHITYP`` is ``differential`` is read with a
    loud warning: a differential phase is not the absolute phase the
    container's complex value asserts.

    Parameters
    ----------
    path
        The OIFITS file.
    target
        The ``OI_TARGET`` ``TARGET`` name to read, required when the selected
        tables hold more than one.
    insname
        The ``INSNAME`` to read, required when the data tables name more than
        one instrument.
    arrname
        The ``ARRNAME`` to read, required when the selected data tables name
        more than one array.

    Returns
    -------
    OIFITSData

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    OIFITSError
        When the file has no ``OI_VIS``, ``OI_VIS2`` or ``OI_T3`` table, lacks
        a table the data need (``OI_WAVELENGTH`` for the instrument,
        ``OI_TARGET``), is ambiguous without a keyword, or names a selection
        that is not there.
    """
    source = Path(os.fspath(path))
    if not source.is_file():
        raise FileNotFoundError(
            f"read_oifits: no OIFITS file at {str(source)!r} (resolved to "
            f"{str(source.resolve())!r})."
        )
    with fits.open(source, memmap=False) as hdus:
        tables = [hdu for hdu in hdus if isinstance(hdu, fits.BinTableHDU)]
        positions = {id(hdu): number for number, hdu in enumerate(hdus)}
        return _read(source, tables, positions, target=target, insname=insname, arrname=arrname)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def _name(value: Any) -> str:
    """A FITS string cell or header value as a stripped ``str``."""
    if isinstance(value, bytes):
        value = value.decode("ascii", errors="replace")
    return str(value).strip()


def _header_name(hdu: fits.BinTableHDU, key: str) -> str | None:
    value = hdu.header.get(key)
    return None if value is None else _name(value)


def _choose(
    kind: str, keyword: str, choices: Sequence[str], wanted: str | None, file: Path
) -> str | None:
    """Pick one of *choices* by *wanted*, refusing ambiguity and absence by name."""
    unique = sorted(set(choices))
    if wanted is not None:
        if wanted not in unique:
            raise OIFITSError(
                f"read_oifits: {file.name} has no {kind} {wanted!r}; the data tables name "
                f"{unique or 'none'}. Pass {keyword}= one of those."
            )
        return wanted
    if len(unique) > 1:
        raise OIFITSError(
            f"read_oifits: {file.name} holds data for {len(unique)} {kind}s, {unique}, and a "
            f"container holds one. Pass {keyword}= to choose; the reader does not pick the "
            f"first for you."
        )
    return unique[0] if unique else None


def _read(
    file: Path,
    tables: Sequence[fits.BinTableHDU],
    positions: Mapping[int, int],
    *,
    target: str | None,
    insname: str | None,
    arrname: str | None,
) -> OIFITSData:
    data = [hdu for hdu in tables if hdu.name in DATA_TABLES]
    if not data:
        names = sorted({hdu.name for hdu in tables}) or ["no binary tables"]
        raise OIFITSError(
            f"read_oifits: {file.name} has no OI_VIS, OI_VIS2 or OI_T3 table (it has {names}), "
            f"so there is nothing to put in a container. Is it an OIFITS file?"
        )

    instrument = _choose(
        "instrument", "insname", [_header_name(h, "INSNAME") or "" for h in data], insname, file
    )
    assert instrument is not None
    data = [h for h in data if (_header_name(h, "INSNAME") or "") == instrument]
    array_names = [n for n in (_header_name(h, "ARRNAME") for h in data) if n]
    array = _choose("array", "arrname", array_names, arrname, file)
    if array is not None:
        data = [h for h in data if _header_name(h, "ARRNAME") in (array, None)]

    wavelength_tables = [
        h for h in tables if h.name == "OI_WAVELENGTH" and _header_name(h, "INSNAME") == instrument
    ]
    if not wavelength_tables:
        raise OIFITSError(
            f"read_oifits: {file.name} has data for instrument {instrument!r} but no "
            f"OI_WAVELENGTH table with INSNAME = {instrument!r}, so its channels have no "
            f"wavelength."
        )
    channels = wavelength_tables[0]
    waves = _column(channels, "EFF_WAVE", _METRES)
    bands = _column(channels, "EFF_BAND", _METRES)

    target_tables = [h for h in tables if h.name == "OI_TARGET"]
    if not target_tables:
        raise OIFITSError(
            f"read_oifits: {file.name} has no OI_TARGET table, so its rows name no source. "
            f"OIFITS requires one."
        )
    targets = {
        int(row_id): _name(row_name)
        for row_id, row_name in zip(
            target_tables[0].data["TARGET_ID"], target_tables[0].data["TARGET"], strict=True
        )
    }
    present = sorted({int(i) for h in data for i in np.atleast_1d(h.data["TARGET_ID"])})
    unknown = [i for i in present if i not in targets]
    if unknown:
        raise OIFITSError(
            f"read_oifits: {file.name}'s data rows use TARGET_ID {unknown}, which OI_TARGET "
            f"does not list."
        )
    source_name = _choose("target", "target", [targets[i] for i in present], target, file)
    assert source_name is not None
    target_ids = {i for i in present if targets[i] == source_name}

    stations = _station_names(tables, array)
    revisions = sorted({int(h.header.get("OI_REVN", 1)) for h in data})
    context = _Context(
        file=file,
        waves=waves,
        bands=bands,
        stations=stations,
        target_ids=target_ids,
        positions=positions,
        common={
            "file": file.name,
            "target": source_name,
            "instrument": instrument,
            "array": array,
            "oi_revn": revisions[-1],
        },
    )
    # A table none of whose rows belong to the target contributes nothing, and
    # a kind with no such table at all is absent: ``None``, never empty.
    by_kind = {
        kind: [h for h in data if h.name == kind and _rows(h, context).any()]
        for kind in DATA_TABLES
    }
    return OIFITSData(
        visibilities=_visibilities(by_kind["OI_VIS"], context) if by_kind["OI_VIS"] else None,
        squared_visibilities=(
            _squared_visibilities(by_kind["OI_VIS2"], context) if by_kind["OI_VIS2"] else None
        ),
        closure_phases=_closure_phases(by_kind["OI_T3"], context) if by_kind["OI_T3"] else None,
        target=source_name,
        instrument=instrument,
        array=array,
        revision=revisions[-1],
        wavelengths=waves * _METRES,
        bandwidths=bands * _METRES,
    )


def _station_names(tables: Sequence[fits.BinTableHDU], array: str | None) -> dict[int, str]:
    """``STA_INDEX -> STA_NAME`` from the array's ``OI_ARRAY``, or empty when there is none."""
    arrays = [h for h in tables if h.name == "OI_ARRAY"]
    if array is not None:
        arrays = [h for h in arrays if _header_name(h, "ARRNAME") == array]
    if len(arrays) != 1:
        return {}
    rows = arrays[0].data
    return {
        int(index): _name(name)
        for index, name in zip(rows["STA_INDEX"], rows["STA_NAME"], strict=True)
    }


# ---------------------------------------------------------------------------
# Columns
# ---------------------------------------------------------------------------


def _column_unit(hdu: fits.BinTableHDU, name: str) -> u.UnitBase | None:
    position = hdu.columns.names.index(name)
    unit = hdu.columns[position].unit
    if not unit or not str(unit).strip():
        return None
    try:
        return u.Unit(str(unit).strip())
    except ValueError:
        return None


def _column(hdu: fits.BinTableHDU, name: str, unit: u.UnitBase | None = None) -> np.ndarray:
    """A column as float64, converted to *unit* from its ``TUNIT`` when it declares one."""
    if name not in hdu.columns.names:
        raise OIFITSError(f"read_oifits: the {hdu.name} table has no {name} column.")
    values = np.asarray(hdu.data[name], dtype=np.float64)
    declared = _column_unit(hdu, name)
    if unit is not None and declared is not None:
        try:
            values = (values * declared).to_value(unit)
        except u.UnitConversionError as exc:
            raise OIFITSError(
                f"read_oifits: {hdu.name}.{name} is in {declared}, which is not convertible to "
                f"{unit}. ({exc})"
            ) from exc
    return values


def _per_channel(hdu: fits.BinTableHDU, name: str, n_channels: int, unit: Any = None) -> Any:
    """A per-channel column as ``(n_rows, n_channels)``, checked against the instrument."""
    values = _column(hdu, name, unit) if name != "FLAG" else np.asarray(hdu.data[name], dtype=bool)
    values = values.reshape(len(hdu.data), -1)
    if values.shape[1] != n_channels:
        raise OIFITSError(
            f"read_oifits: {hdu.name}.{name} has {values.shape[1]} channels per row but the "
            f"instrument's OI_WAVELENGTH has {n_channels}."
        )
    return values


@dataclasses.dataclass(frozen=True)
class _Context:
    file: Path
    waves: np.ndarray
    bands: np.ndarray
    stations: dict[int, str]
    target_ids: set[int]
    positions: Mapping[int, int]
    common: Mapping[str, Any]

    def station(self, index: int) -> str:
        return self.stations.get(index, str(index))

    def label(self, *indices: int) -> str:
        return "-".join(self.station(i) for i in indices)


@dataclasses.dataclass
class _Block:
    """The flattened samples of one kind, accumulated over the matching tables."""

    columns: dict[str, list[np.ndarray]] = dataclasses.field(default_factory=dict)
    extensions: list[str] = dataclasses.field(default_factory=list)
    dates: list[str] = dataclasses.field(default_factory=list)

    def add(self, **arrays: np.ndarray) -> None:
        for key, value in arrays.items():
            self.columns.setdefault(key, []).append(np.asarray(value).reshape(-1))

    def get(self, key: str) -> np.ndarray:
        return np.concatenate(self.columns[key])


def _rows(hdu: fits.BinTableHDU, context: _Context) -> np.ndarray:
    """The rows of *hdu* that belong to the selected target."""
    ids = np.atleast_1d(np.asarray(hdu.data["TARGET_ID"], dtype=int))
    return np.isin(ids, sorted(context.target_ids))


def _grid(per_row: np.ndarray, n_channels: int) -> np.ndarray:
    """A per-row value repeated across the channels, ``(n_rows, n_channels)``."""
    return np.repeat(np.asarray(per_row).reshape(-1, 1), n_channels, axis=1)


def _clean(
    values: np.ndarray, errors: np.ndarray, flags: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Mask flagged, non-finite and non-positive-error samples; placeholders where masked."""
    nonfinite = ~(np.isfinite(values) & np.isfinite(errors))
    nonpositive = np.isfinite(errors) & (errors <= 0.0)
    mask = flags | nonfinite | nonpositive
    counts = {
        "n_flagged": int(flags.sum()),
        "n_nonfinite": int(nonfinite.sum()),
        "n_nonpositive_error": int(nonpositive.sum()),
        "n_masked": int(mask.sum()),
    }
    return np.where(mask, 0.0, values), np.where(mask, 1.0, errors), mask, counts


def _meta(context: _Context, block: _Block, observable: str, counts: dict[str, int]) -> dict:
    meta: dict[str, Any] = {
        **context.common,
        "observable": observable,
        "tables": list(block.extensions),
        **counts,
    }
    dates = sorted(set(block.dates))
    if dates:
        meta["date_obs"] = dates[0] if len(dates) == 1 else dates
    return meta


def _start(hdu: fits.BinTableHDU, context: _Context, block: _Block) -> None:
    """Record the table a block's rows came from, by its HDU number in the file."""
    block.extensions.append(f"{hdu.name} (HDU {context.positions.get(id(hdu), -1)})")
    date = hdu.header.get("DATE-OBS")
    if date is not None and _name(date):
        block.dates.append(_name(date))


def _canonical_pair(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Sorted station pairs, and ``-1`` where the file listed the pair reversed."""
    flip = np.where(a > b, -1.0, 1.0)
    return np.minimum(a, b), np.maximum(a, b), flip


# ---------------------------------------------------------------------------
# The three tables
# ---------------------------------------------------------------------------


def _baseline_rows(
    hdu: fits.BinTableHDU, context: _Context, block: _Block
) -> tuple[np.ndarray, np.ndarray, int]:
    """The per-row geometry of an ``OI_VIS``/``OI_VIS2`` table into *block*; the flip and size."""
    keep = _rows(hdu, context)
    n_channels = context.waves.size
    stations = np.asarray(hdu.data["STA_INDEX"], dtype=int).reshape(len(hdu.data), 2)[keep]
    first, second, flip = _canonical_pair(stations[:, 0], stations[:, 1])
    ucoord = _column(hdu, "UCOORD", _METRES)[keep] * flip
    vcoord = _column(hdu, "VCOORD", _METRES)[keep] * flip
    labels = np.array([context.label(int(i), int(j)) for i, j in zip(first, second, strict=True)])
    block.add(
        u=ucoord.reshape(-1, 1) / context.waves,
        v=vcoord.reshape(-1, 1) / context.waves,
        wave=_grid(np.ones(keep.sum()), n_channels) * context.waves,
        band=_grid(np.ones(keep.sum()), n_channels) * context.bands,
        channel=_grid(np.ones(keep.sum(), dtype=int), n_channels) * np.arange(n_channels),
        mjd=_grid(_column(hdu, "MJD")[keep], n_channels),
        baseline=_grid(labels, n_channels),
    )
    _start(hdu, context, block)
    return keep, flip, n_channels


def _squared_visibilities(tables: Sequence[fits.BinTableHDU], context: _Context) -> VisibilitySet:
    block = _Block()
    for hdu in tables:
        keep, _, n_channels = _baseline_rows(hdu, context, block)
        block.add(
            value=_per_channel(hdu, "VIS2DATA", n_channels)[keep],
            error=_per_channel(hdu, "VIS2ERR", n_channels)[keep],
            flag=_per_channel(hdu, "FLAG", n_channels)[keep],
        )
    values, errors, mask, counts = _clean(
        block.get("value"), block.get("error"), block.get("flag").astype(bool)
    )
    return VisibilitySet(
        block.get("u"),
        block.get("v"),
        (block.get("wave") * _METRES).to(SPECTRAL_UNIT),
        values,
        uncertainty=errors,
        mask=mask,
        extra_coords=_extra(block, ("baseline",)),
        meta=_meta(context, block, "squared_visibility", counts),
    )


def _visibilities(tables: Sequence[fits.BinTableHDU], context: _Context) -> VisibilitySet:
    block = _Block()
    unit: u.UnitBase | None = None
    for index, hdu in enumerate(tables):
        for keyword in ("AMPTYP", "PHITYP"):
            kind = _header_name(hdu, keyword)
            if kind is not None and kind.lower().startswith("differential"):
                warnings.warn(
                    f"ampere.interferometry.read_oifits: {context.file.name} {hdu.name} has "
                    f"{keyword} = {kind!r}; the complex VisibilitySet treats VISAMP and VISPHI as "
                    f"an absolute amplitude and phase, which a differential one is not. Fit the "
                    f"amplitudes alone, or do not fit this table.",
                    UserWarning,
                    stacklevel=4,
                )
        keep, flip, n_channels = _baseline_rows(hdu, context, block)
        declared = _column_unit(hdu, "VISAMP")
        if declared is not None and declared.physical_type == "dimensionless":
            declared = None
        if index and declared != unit:
            raise OIFITSError(
                f"read_oifits: {context.file.name}'s OI_VIS tables disagree on VISAMP's unit "
                f"({unit} and {declared}); read them separately with insname=/arrname=."
            )
        unit = declared
        amplitude = _per_channel(hdu, "VISAMP", n_channels)[keep]
        phase = np.deg2rad(_per_channel(hdu, "VISPHI", n_channels, u.deg)[keep]) * flip.reshape(
            -1, 1
        )
        block.add(
            value=amplitude * np.exp(1j * phase),
            error=_per_channel(hdu, "VISAMPERR", n_channels)[keep],
            flag=_per_channel(hdu, "FLAG", n_channels)[keep],
        )
    raw = block.get("value")
    # ``_clean`` judges finiteness on a real array; the modulus is non-finite
    # exactly when either component is.
    _, errors, mask, counts = _clean(
        np.abs(raw), block.get("error"), block.get("flag").astype(bool)
    )
    complex_values = np.where(mask, 0.0 + 0.0j, raw)
    scale: Any = 1.0 if unit is None else unit
    return VisibilitySet(
        block.get("u"),
        block.get("v"),
        (block.get("wave") * _METRES).to(SPECTRAL_UNIT),
        complex_values * scale,
        uncertainty=errors * scale,
        mask=mask,
        extra_coords=_extra(block, ("baseline",)),
        meta=_meta(context, block, "complex_visibility", counts),
    )


def _permutation_is_odd(triple: np.ndarray) -> np.ndarray:
    """Per row, whether sorting the three station indices is an odd permutation."""
    a, b, c = triple[:, 0], triple[:, 1], triple[:, 2]
    inversions = (a > b).astype(int) + (a > c).astype(int) + (b > c).astype(int)
    return inversions % 2 == 1


def _wrap_degrees(phase: np.ndarray) -> np.ndarray:
    """Into ``(-180, 180]``: ``180`` stays ``180`` and ``-180`` becomes ``180``."""
    return 180.0 - np.mod(180.0 - phase, 360.0)


def _closure_phases(tables: Sequence[fits.BinTableHDU], context: _Context) -> ClosurePhases:
    block = _Block()
    for hdu in tables:
        keep = _rows(hdu, context)
        n_channels = context.waves.size
        triples = np.asarray(hdu.data["STA_INDEX"], dtype=int).reshape(len(hdu.data), 3)[keep]
        if np.any(
            (triples[:, 0] == triples[:, 1])
            | (triples[:, 1] == triples[:, 2])
            | (triples[:, 0] == triples[:, 2])
        ):
            raise OIFITSError(
                f"read_oifits: {context.file.name} {hdu.name} has a row whose STA_INDEX repeats "
                f"a station, which is not a triangle."
            )
        u1 = _column(hdu, "U1COORD", _METRES)[keep]
        v1 = _column(hdu, "V1COORD", _METRES)[keep]
        u2 = _column(hdu, "U2COORD", _METRES)[keep]
        v2 = _column(hdu, "V2COORD", _METRES)[keep]
        rows = []
        for (a, b, c), ab_u, ab_v, bc_u, bc_v in zip(triples, u1, v1, u2, v2, strict=True):
            # The three oriented baselines the file states, and their reverses.
            edges = {
                (a, b): (ab_u, ab_v),
                (b, c): (bc_u, bc_v),
                (c, a): (-(ab_u + bc_u), -(ab_v + bc_v)),
            }
            edges.update({(y, x): (-p, -q) for (x, y), (p, q) in list(edges.items())})
            i, j, k = sorted((int(a), int(b), int(c)))
            rows.append(
                (
                    *edges[(i, j)],
                    *edges[(j, k)],
                    context.label(i, j, k),
                    context.label(i, j),
                    context.label(j, k),
                    context.label(k, i),
                )
            )
        geometry = np.array([row[:4] for row in rows], dtype=np.float64).reshape(-1, 4)
        labels = np.array([row[4:] for row in rows], dtype=str).reshape(-1, 4)
        sign = np.where(_permutation_is_odd(triples), -1.0, 1.0).reshape(-1, 1)
        phase = _per_channel(hdu, "T3PHI", n_channels, u.deg)[keep]
        rows_n = int(keep.sum())
        block.add(
            u1=geometry[:, 0:1] / context.waves,
            v1=geometry[:, 1:2] / context.waves,
            u2=geometry[:, 2:3] / context.waves,
            v2=geometry[:, 3:4] / context.waves,
            wave=_grid(np.ones(rows_n), n_channels) * context.waves,
            band=_grid(np.ones(rows_n), n_channels) * context.bands,
            channel=_grid(np.ones(rows_n, dtype=int), n_channels) * np.arange(n_channels),
            mjd=_grid(_column(hdu, "MJD")[keep], n_channels),
            triangle=_grid(labels[:, 0], n_channels),
            baseline_ij=_grid(labels[:, 1], n_channels),
            baseline_jk=_grid(labels[:, 2], n_channels),
            baseline_ki=_grid(labels[:, 3], n_channels),
            value=_wrap_degrees(sign * phase),
            error=_per_channel(hdu, "T3PHIERR", n_channels, u.deg)[keep],
            flag=_per_channel(hdu, "FLAG", n_channels)[keep],
        )
        _start(hdu, context, block)
    values, errors, mask, counts = _clean(
        block.get("value"), block.get("error"), block.get("flag").astype(bool)
    )
    return ClosurePhases(
        block.get("u1"),
        block.get("v1"),
        block.get("u2"),
        block.get("v2"),
        (block.get("wave") * _METRES).to(SPECTRAL_UNIT),
        (values * u.deg).to(u.rad),
        uncertainty=(errors * u.deg).to(u.rad),
        mask=mask,
        extra_coords=_extra(block, ("triangle", "baseline_ij", "baseline_jk", "baseline_ki")),
        meta=_meta(context, block, "closure_phase", counts),
    )


def _extra(block: _Block, labels: Sequence[str]) -> dict[str, np.ndarray]:
    """The per-sample ``extra_coords``: the station labels, ``channel``, ``mjd``, ``eff_band``."""
    extra: dict[str, np.ndarray] = {name: block.get(name).astype(str) for name in labels}
    extra["channel"] = block.get("channel").astype(int)
    extra["mjd"] = block.get("mjd").astype(np.float64)
    extra["eff_band"] = block.get("band").astype(np.float64)
    return extra
