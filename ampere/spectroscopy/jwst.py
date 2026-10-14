"""The JWST reader: one pipeline one-dimensional spectrum into a ``Spectrum`` (W7.3).

The JWST calibration pipeline writes its extracted one-dimensional spectra as
FITS files of the ``MultiSpecModel`` data model: an ``x1d`` product holds one
``EXTRACT1D`` binary table per extracted source, slit or spectral order
(numbered by ``EXTVER``), and a ``c1d`` product one ``COMBINE1D`` table per
combined source (``jwst-pipeline.readthedocs.io``, "Data products"). Each table
carries ``WAVELENGTH``, ``FLUX`` and its uncertainty (``FLUX_ERROR`` in an
``x1d``, ``ERROR`` in a ``c1d``), the surface-brightness pair
``SURF_BRIGHT``/``SB_ERROR``, the variance components, ``DQ``, the background
columns and ``NPIXELS``, every unit in its ``TUNITn``; the primary header names
the instrument, the optical elements, the exposure and the pipeline version.
Only ``x1d`` and ``c1d`` are read (D6): cubes (``s3d``), rectified 2-D spectra
(``s2d``) and the per-integration ``x1dints`` of time-series modes are not.

The reader is :func:`read_jwst`; it uses :mod:`astropy.io.fits` and no other
dependency (no ``jwst``, no ``stdatamodels``). Its docstring states the
conventions it enforces, which are :class:`~ampere.core.Spectrum`'s
(``ampere/core/results_schema.py``) and the reference steps': the reader
translates the file's conventions into those, once, at composition time.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits

from ampere.core import Spectrum

__all__ = ["DO_NOT_USE", "JWSTError", "read_jwst"]

#: The pipeline's ``DO_NOT_USE`` data-quality flag: bit 0, value 1, of the
#: ``DQ`` column (``stdatamodels``' ``dqflags.pixel["DO_NOT_USE"]``). The only
#: bit the reader masks on; every other bit is kept, raw, in
#: ``extra_coords["dq"]`` for a user who wants a stricter mask.
DO_NOT_USE = 1

#: The binary tables that hold a spectrum, and the product each one marks.
SPECTRUM_TABLES = {"EXTRACT1D": "x1d", "COMBINE1D": "c1d"}

#: The unit the container's spectral axis is declared in: micron, the unit the
#: reference steps (``Resample``, ``LSFConvolution``) and models take their
#: grids in, and the OIFITS reader's precedent. The likelihood's alignment
#: check compares the predicted and observed axes exactly, unit included.
SPECTRAL_UNIT = u.micron

#: Primary-header keywords copied into ``meta``, under the names given.
_PRIMARY_META = (
    ("instrument", "INSTRUME"),
    ("detector", "DETECTOR"),
    ("filter", "FILTER"),
    ("grating", "GRATING"),
    ("exp_type", "EXP_TYPE"),
    ("program", "PROGRAM"),
    ("observation", "OBSERVTN"),
    ("visit", "VISIT"),
    ("exposure", "EXPOSURE"),
    ("date_obs", "DATE-OBS"),
    ("pipeline_version", "CAL_VER"),
    ("crds_context", "CRDS_CTX"),
)


class JWSTError(ValueError):
    """A JWST file this reader cannot turn into a ``Spectrum``, named with the reason."""


def read_jwst(
    path: str | os.PathLike[str],
    *,
    source: str | int | None = None,
    extension: int | None = None,
    surface_brightness: bool = False,
) -> Spectrum:
    """Read a JWST pipeline ``x1d`` or ``c1d`` file into a :class:`~ampere.core.Spectrum`.

    **Selection.** A ``Spectrum`` holds one spectrum, so the reader reads one
    ``EXTRACT1D`` (``x1d``) or ``COMBINE1D`` (``c1d``) table per call. A file
    with one such table needs no keyword. A file with several (a NIRSpec
    fixed-slit exposure with one table per slit, a multi-object product with
    one per source) is **refused** unless a keyword says which, and the
    refusal lists the choices; the reader never picks the first one for you.

    * ``source=`` is matched against each table's header, in this order:
      ``SLTNAME`` (the slit or aperture name, ``"S200A1"``), then
      ``SOURCEID`` (compared as text, so ``41`` and ``"41"`` both match), then
      ``SRCNAME``. The first of the three that matches any table decides; a
      value matching several tables (one source in several spectral orders)
      is refused, naming them, so that ``extension=`` can say which.
    * ``extension=`` is the table's ``EXTVER``: the pipeline's own numbering
      of its spectra, the ``Ver`` column of ``astropy.io.fits``'s
      ``info()``, so that ``hdus["EXTRACT1D", 2]`` and ``extension=2`` name
      the same table. Given with ``source=``, both must name the same table.

    A selection that names nothing in the file is refused, listing what is
    there.

    **Columns and units.** The spectral axis is ``WAVELENGTH``, converted from
    its ``TUNIT`` to **micron** (:data:`SPECTRAL_UNIT`). The values and
    uncertainties are ``FLUX``/``FLUX_ERROR`` (``ERROR`` in a ``c1d``;
    ``surface_brightness=False``, the default) or ``SURF_BRIGHT``/``SB_ERROR``
    (``True``), in the unit the file declares for the value column (typically
    Jy and MJy/sr): the unit is
    read from ``TUNIT``, never assumed, and a column with no unit, or one
    astropy cannot parse, is refused by name. The uncertainty is converted to
    the value's unit.

    **Masks.** ``mask = (DQ & DO_NOT_USE) | non-finite value or uncertainty |
    non-positive uncertainty`` — ``results_schema.md`` §7's ``True`` is
    excluded. Only the pipeline's ``DO_NOT_USE`` bit masks; the other ``DQ``
    bits are informational and stay in ``extra_coords["dq"]``. A masked
    sample's value and uncertainty are replaced by the placeholders ``0`` and
    ``1`` (the container refuses negative uncertainties, and a masked sample is
    not data). ``meta`` counts ``n_do_not_use``, ``n_nonfinite``,
    ``n_nonpositive_error`` and ``n_masked`` over the container's samples.

    **Ordering.** ``Spectrum`` wants a strictly increasing axis. A table whose
    wavelengths decrease (MIRI's low-resolution spectrometer writes its
    spectra from long to short wavelength) or are otherwise out of order is
    sorted, with every column following, and ``meta["reordered"]`` says
    ``"reversed"`` or ``"sorted"`` (``None`` when the file was already in
    order). A row whose wavelength is not finite cannot be placed on the axis
    and is dropped (``meta["n_nonfinite_wavelength"]``). Rows sharing one
    wavelength (an order overlap, or a padding row repeating the last
    wavelength, which MIRI's LRS products carry) cannot all be kept, since the
    axis refuses duplicated coordinates: the first unmasked row of each group,
    in file order, is kept — or the first row when all are masked — and the
    rest dropped (``meta["n_duplicate"]``), with a warning when a dropped row
    was unmasked data. ``meta["n_rows"]`` is the table's row count.

    **What is kept, and what is not.** ``extra_coords`` carry, per sample, the
    raw ``dq`` bits, ``npixels`` (the extraction aperture's width in pixels)
    and, in a ``c1d``, ``n_input`` (the number of spectra combined) where the
    table has them. ``meta`` carries the primary header's
    ``instrument``, ``detector``, ``filter``, ``grating`` (``None`` when the
    instrument has none), ``exp_type``, ``program``, ``observation``,
    ``visit``, ``exposure``, ``target`` (``TARGNAME``, else ``TARGPROP``),
    ``date_obs``, ``pipeline_version`` (``CAL_VER``) and ``crds_context``; the
    table's ``product`` (``"x1d"`` or ``"c1d"``), ``extension`` (``EXTVER``),
    ``hdu`` (its index in the file), ``slit`` (``SLTNAME``), ``source_id``
    (``SOURCEID``), ``srctype`` (``POINT`` or ``EXTENDED``) and ``sporder``;
    and the ``file`` name and ``columns`` read. The background, the variance
    components and MIRI MRS's residual-fringe-corrected ``RF_*`` columns are
    not read. **The line-spread function** is not in the file:
    the pipeline does not tabulate a resolving power per mode in its products,
    so the reader carries none, and the user supplies one to
    :class:`~ampere.backends.reference.LSFConvolution`.

    **Extended sources.** A table with ``SRCTYPE = EXTENDED`` still carries a
    ``FLUX`` column, which the pipeline fills by integrating the surface
    brightness over the extraction aperture; ``surface_brightness=True`` reads
    the surface brightness itself. A table lacking the requested pair is
    refused, naming the keyword that reads the other.

    Parameters
    ----------
    path
        The ``x1d`` or ``c1d`` FITS file.
    source
        The slit name, source identifier or source name of the table to read,
        required when the file holds more than one.
    extension
        The ``EXTVER`` of the table to read; the alternative to *source*.
    surface_brightness
        ``True`` reads ``SURF_BRIGHT``/``SB_ERROR`` instead of
        ``FLUX``/``FLUX_ERROR``.

    Returns
    -------
    Spectrum

    Raises
    ------
    FileNotFoundError
        When *path* does not exist.
    JWSTError
        When the file has no ``EXTRACT1D`` or ``COMBINE1D`` table, holds
        several without a keyword saying which, names a selection that is not
        there, lacks a column the read needs, declares a unit the reader
        cannot use, or holds a per-row (vector) spectrum.
    """
    file = Path(os.fspath(path))
    if not file.is_file():
        raise FileNotFoundError(
            f"read_jwst: no JWST file at {str(file)!r} (resolved to {str(file.resolve())!r})."
        )
    with fits.open(file, memmap=False) as hdus:
        primary = hdus[0].header
        tables = [
            (number, hdu)
            for number, hdu in enumerate(hdus)
            if isinstance(hdu, fits.BinTableHDU) and hdu.name in SPECTRUM_TABLES
        ]
        if not tables:
            names = sorted({hdu.name for hdu in hdus[1:]}) or ["no extensions"]
            raise JWSTError(
                f"read_jwst: {file.name} has no EXTRACT1D or COMBINE1D table (it has {names}), "
                f"so it holds no JWST one-dimensional spectrum. Is it an x1d or c1d product?"
            )
        number, table = _select(file, tables, source=source, extension=extension)
        return _read(file, primary, number, table, surface_brightness=surface_brightness)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def _text(value: Any) -> str | None:
    """A header value as stripped text, ``None`` when absent or blank."""
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("ascii", errors="replace")
    text = str(value).strip()
    return text or None


def _label(table: fits.BinTableHDU) -> str:
    """How a refusal lists one table: its EXTVER and whatever names it carries."""
    parts = [f"{table.name} EXTVER {table.ver}"]
    for key in ("SLTNAME", "SOURCEID", "SRCNAME", "SRCTYPE", "SPORDER"):
        value = _text(table.header.get(key))
        if value is not None:
            parts.append(f"{key} {value!r}")
    return "(" + ", ".join(parts) + ")"


def _select(
    file: Path,
    tables: Sequence[tuple[int, fits.BinTableHDU]],
    *,
    source: str | int | None,
    extension: int | None,
) -> tuple[int, fits.BinTableHDU]:
    """Pick one spectrum table, refusing ambiguity and absence by name."""
    choices = ", ".join(_label(table) for _, table in tables)
    chosen = list(tables)
    if extension is not None:
        chosen = [(n, t) for n, t in chosen if t.ver == extension]
        if not chosen:
            raise JWSTError(
                f"read_jwst: {file.name} has no spectrum table with EXTVER {extension!r}; it "
                f"has {choices}. Pass extension= one of those EXTVERs."
            )
    if source is not None:
        wanted = str(source).strip()
        for key in ("SLTNAME", "SOURCEID", "SRCNAME"):
            matched = [(n, t) for n, t in chosen if _text(t.header.get(key)) == wanted]
            if matched:
                chosen = matched
                break
        else:
            raise JWSTError(
                f"read_jwst: {file.name} has no spectrum table with SLTNAME, SOURCEID or "
                f"SRCNAME {wanted!r}"
                + (f" and EXTVER {extension!r}" if extension is not None else "")
                + f"; it has {choices}. Pass source= one of those names."
            )
    if len(chosen) > 1:
        listed = ", ".join(_label(table) for _, table in chosen)
        hint = (
            "Pass extension= (the EXTVER) to choose among them"
            if source is not None
            else "Pass source= (a slit name or source id) or extension= (the EXTVER) to choose"
        )
        matching = f" matching {source!r}" if source is not None else ""
        raise JWSTError(
            f"read_jwst: {file.name} holds {len(chosen)} spectra{matching}, "
            f"{listed}, and a Spectrum holds one. {hint}; the reader does not pick the first "
            f"for you."
        )
    return chosen[0]


# ---------------------------------------------------------------------------
# Columns
# ---------------------------------------------------------------------------


def _unit(file: Path, table: fits.BinTableHDU, name: str) -> u.UnitBase:
    """A column's ``TUNIT`` as an astropy unit, refused by name when absent or unparseable."""
    declared = table.columns[table.columns.names.index(name)].unit
    text = _text(declared)
    if text is None:
        raise JWSTError(
            f"read_jwst: {file.name}'s {table.name} column {name} declares no unit (TUNIT); "
            f"the reader reads units from the file and does not assume one."
        )
    try:
        return u.Unit(text)
    except ValueError as exc:
        raise JWSTError(
            f"read_jwst: {file.name}'s {table.name} column {name} declares the unit {text!r}, "
            f"which astropy cannot parse. ({exc})"
        ) from exc


def _column(file: Path, table: fits.BinTableHDU, name: str) -> np.ndarray:
    """A one-value-per-row column as a float64 vector; a vector column is refused."""
    values = np.asarray(table.data[name], dtype=np.float64)
    if values.ndim != 1:
        raise JWSTError(
            f"read_jwst: {file.name}'s {table.name} column {name} holds {values.shape[1:]} "
            f"values per row; the reader reads one spectrum per table, one sample per row "
            f"(the x1d/c1d layout), not a per-row or per-integration (x1dints) layout."
        )
    return values


def _value_columns(
    file: Path, table: fits.BinTableHDU, surface_brightness: bool
) -> tuple[str, str]:
    """The value and uncertainty columns to read, refused by name when absent.

    The flux's uncertainty is ``FLUX_ERROR`` in an ``x1d`` table and ``ERROR``
    in the ``c1d`` layout the pipeline documents; either is accepted.
    """
    names = table.columns.names
    flux_error = "FLUX_ERROR" if "FLUX_ERROR" in names or "ERROR" not in names else "ERROR"
    flux, brightness = ("FLUX", flux_error), ("SURF_BRIGHT", "SB_ERROR")
    wanted, other, keyword = (
        (brightness, flux, "surface_brightness=False")
        if surface_brightness
        else (flux, brightness, "surface_brightness=True")
    )
    missing = [name for name in ("WAVELENGTH", *wanted) if name not in names]
    if missing:
        alternative = (
            f" It has {other[0]}/{other[1]}: pass {keyword} to read those."
            if all(name in names for name in other)
            else ""
        )
        raise JWSTError(
            f"read_jwst: {file.name}'s {table.name} table (EXTVER {table.ver}) has no "
            f"{', '.join(missing)} column (it has {list(names)}).{alternative}"
        )
    return wanted


# ---------------------------------------------------------------------------
# The read
# ---------------------------------------------------------------------------


def _keep(wavelength: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, int, int]:
    """Indices of the rows that go on the axis; the counts of the two kinds dropped.

    Rows with a non-finite wavelength are dropped; of rows sharing a wavelength,
    the first unmasked one in file order is kept (the first row when all are
    masked).
    """
    finite = np.flatnonzero(np.isfinite(wavelength))
    n_nonfinite_wavelength = wavelength.size - finite.size
    kept: dict[float, int] = {}
    for row in finite:
        value = float(wavelength[row])
        if value not in kept or (mask[kept[value]] and not mask[row]):
            kept[value] = int(row)
    rows = np.array(sorted(kept.values()), dtype=np.intp)
    n_duplicate = finite.size - rows.size
    dropped = np.setdiff1d(finite, rows)
    if np.any(~mask[dropped]):
        warnings.warn(
            f"ampere.spectroscopy.read_jwst: {int(np.sum(~mask[dropped]))} unmasked row(s) "
            f"share a wavelength with an earlier unmasked row and were dropped, since a "
            f"Spectrum's axis cannot repeat a coordinate; the first in file order was kept. "
            f"Read the table with astropy.io.fits to treat them otherwise.",
            UserWarning,
            stacklevel=4,
        )
    return rows, n_nonfinite_wavelength, n_duplicate


def _read(
    file: Path,
    primary: fits.Header,
    number: int,
    table: fits.BinTableHDU,
    *,
    surface_brightness: bool,
) -> Spectrum:
    value_name, error_name = _value_columns(file, table, surface_brightness)
    names = table.columns.names

    wave_unit = _unit(file, table, "WAVELENGTH")
    value_unit = _unit(file, table, value_name)
    error_unit = _unit(file, table, error_name)
    wavelength = _column(file, table, "WAVELENGTH")
    values = _column(file, table, value_name)
    errors = _column(file, table, error_name)
    try:
        wavelength = (wavelength * wave_unit).to_value(SPECTRAL_UNIT, equivalencies=u.spectral())
        errors = (errors * error_unit).to_value(value_unit)
    except u.UnitConversionError as exc:
        raise JWSTError(
            f"read_jwst: {file.name}'s {table.name} table declares units the reader cannot "
            f"reconcile (WAVELENGTH in {wave_unit}, {value_name} in {value_unit}, {error_name} "
            f"in {error_unit}). ({exc})"
        ) from exc

    if "DQ" in names:
        dq = np.asarray(table.data["DQ"]).astype(np.int64).reshape(-1)
    else:
        dq = np.zeros(values.size, dtype=np.int64)
    flagged = (dq & DO_NOT_USE) != 0
    nonfinite = ~(np.isfinite(values) & np.isfinite(errors))
    nonpositive = np.isfinite(errors) & (errors <= 0.0)
    mask = flagged | nonfinite | nonpositive

    rows, n_nonfinite_wavelength, n_duplicate = _keep(wavelength, mask)
    order = np.argsort(wavelength[rows], kind="stable")
    if np.all(order == np.arange(rows.size)):
        reordered = None
    elif np.all(order == np.arange(rows.size)[::-1]):
        reordered = "reversed"
    else:
        reordered = "sorted"
    rows = rows[order]

    flagged, nonfinite, nonpositive, mask = (
        a[rows] for a in (flagged, nonfinite, nonpositive, mask)
    )
    extra: dict[str, np.ndarray] = {}
    if "DQ" in names:
        extra["dq"] = dq[rows]
    if "NPIXELS" in names:
        extra["npixels"] = _column(file, table, "NPIXELS")[rows]
    if "N_INPUT" in names:
        extra["n_input"] = _column(file, table, "N_INPUT")[rows]

    meta: dict[str, Any] = {key: _text(primary.get(card)) for key, card in _PRIMARY_META}
    meta.update(
        target=_text(primary.get("TARGNAME")) or _text(primary.get("TARGPROP")),
        product=SPECTRUM_TABLES[table.name],
        extension=int(table.ver),
        hdu=number,
        slit=_text(table.header.get("SLTNAME")),
        source_id=_text(table.header.get("SOURCEID")),
        srctype=_text(table.header.get("SRCTYPE")) or _text(primary.get("SRCTYPE")),
        sporder=_text(table.header.get("SPORDER")),
        file=file.name,
        columns=[value_name, error_name],
        n_rows=int(values.size),
        n_nonfinite_wavelength=int(n_nonfinite_wavelength),
        n_duplicate=int(n_duplicate),
        reordered=reordered,
        n_do_not_use=int(flagged.sum()),
        n_nonfinite=int(nonfinite.sum()),
        n_nonpositive_error=int(nonpositive.sum()),
        n_masked=int(mask.sum()),
    )
    return Spectrum(
        wavelength[rows] * SPECTRAL_UNIT,
        np.where(mask, 0.0, values[rows]) * value_unit,
        uncertainty=np.where(mask, 1.0, errors[rows]) * value_unit,
        mask=mask,
        extra_coords=extra,
        meta=meta,
    )
