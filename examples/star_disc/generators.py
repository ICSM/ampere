"""The data for ``examples.star_disc``: HD 105's photometry, the RVS spectrum, IRS.

* :func:`read_votable` reads the tracked ``HD105_SED.vot`` (nineteen rows) --
  the file the legacy ``examples/star_disc.py`` reads through
  ``Photometry.fromFile(..., format="votable")``.
* :func:`read_limits` reads the two rows of ``HD105_SED.csv`` flagged as upper
  limits (500 and 880 micron), which the votable does not carry, under the
  twin's filter names.
* :func:`rvs_grid` is the legacy RVS window (``star_disc.py`` lines 110-116):
  0.847 micron multiplied by ``1 + 1/11000`` until past 0.871.
* :func:`read_irs` reads a CASSIS ``SPITZER-YAAAR`` file the user has fetched
  (``--irs PATH``): the wavelength split is
  :func:`examples.linear_sed.generators.irs_wavelength_grids`'s, imported, and
  the repeated samples are dropped as
  :func:`~examples.linear_sed.generators.deduplicated_grids` drops them.
* :func:`observe` pushes a model through the instruments at a truth and draws
  noise; :func:`observed_data` assembles the observation for either mode.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate
from examples.linear_sed.generators import irs_wavelength_grids

__all__ = [
    "CSV",
    "DATA_DIR",
    "LIMIT_NAMES",
    "MARSHALL",
    "RVS_FRACTIONAL_NOISE",
    "SEED",
    "SYNTHETIC_TRUTH",
    "VOTABLE",
    "observe",
    "read_irs",
    "read_limits",
    "read_votable",
    "rvs_grid",
]

SEED = 20260930

#: The tracked data sit in this package directory (both files byte-identical).
DATA_DIR = Path(__file__).resolve().parent
VOTABLE = DATA_DIR / "HD105_SED.vot"
CSV = DATA_DIR / "HD105_SED.csv"

#: The CSV's filter codes (the legacy's short names) for its two flagged rows,
#: by the twin's names: SPIRE 500 micron is the library's ``HERSCHEL_SPIRE_PLW``
#: (the votable names PSW and PMW the same way); LABOCA is not in the library
#: and is a top-hat (:data:`.star_disc.TOPHATS`).
LIMIT_NAMES: dict[str, str] = {
    "HRSL_SPRE": "HERSCHEL_SPIRE_PLW",
    "APEX_LBCA": "APEX/LABOCA.870",
}

#: Marshall et al. (2018)'s HD 105 (``star_disc.py`` lines 96-101), by the
#: twin's names: the synthetic RVS spectrum is drawn here in both modes.
MARSHALL: dict[str, float] = {"luminosity": 1.216, "teff": 6034.0, "logg": 4.478, "feh": 0.02}

#: ``--synthetic``'s truth: the Marshall star plus a dust belt inside the box.
SYNTHETIC_TRUTH: dict[str, float] = {
    **MARSHALL,
    "log_area": 0.5,
    "t_dust": 60.0,
    "lambda_0": 150.0,
    "beta": 1.0,
}

#: ``input_noise_spec = 0.05`` (S/N 20).
RVS_FRACTIONAL_NOISE = 0.05


def rvs_grid() -> np.ndarray:
    """The legacy RVS window exactly (``star_disc.py`` lines 110-116), micron."""
    waves = [0.847]
    while waves[-1] < 0.871:
        waves.append(waves[-1] * (1 + (1 / 11000)))
    return np.array(waves)


def read_votable(path: Path | str = VOTABLE) -> tuple[list[str], np.ndarray, np.ndarray]:
    """``(filter names, flux Jy, uncertainty Jy)`` from the votable, in file order."""
    table = parse_single_table(str(path)).to_table()
    names = [str(name).strip() for name in table["sed_filter"]]
    flux = np.asarray(table["sed_flux"].to(u.Jy).value, dtype=float)
    error = np.asarray(table["sed_eflux"].to(u.Jy).value, dtype=float)
    return names, flux, error


def read_limits(path: Path | str = CSV) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """``(names, wavelength um, flux Jy, uncertainty Jy)`` of the CSV's upper-limit rows.

    Only the rows with ``upper_limit == 1`` are read, in file order, with the
    legacy filter codes mapped through :data:`LIMIT_NAMES`. The flux and
    uncertainty are the CSV's exactly; the wavelength is the catalogue's (the
    photometry step's own pivot replaces it in the fit).
    """
    table = Table.read(str(path), format="ascii.csv")
    flagged = table[np.asarray(table["upper_limit"]) == 1]
    names = [LIMIT_NAMES[str(code).strip()] for code in flagged["filter"]]
    return (
        names,
        np.asarray(flagged["wave_um"], dtype=float),
        np.asarray(flagged["flux_Jy"], dtype=float),
        np.asarray(flagged["error_Jy"], dtype=float),
    )


def read_irs(path: Path | str) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """SL and LL of a CASSIS file as ``(wavelength, flux, uncertainty)``, deduplicated.

    The wavelength split and order are
    :func:`examples.linear_sed.generators.irs_wavelength_grids`'s (asserted);
    at a repeated wavelength the first sample is kept.
    """
    with fits.open(path) as hdul:
        header = hdul[0].header
        names = [header[f"COL{index:02d}DEF"] for index in range(1, 16)] + ["DUMMY"]
        table = Table(hdul[0].data, names=names)
    finite = np.logical_and.reduce([np.isfinite(column) for column in table.columns.values()])
    table = table[finite]
    table.sort(keys="wavelength")
    chunks = []
    grids = irs_wavelength_grids(path)
    for modules, grid in zip(((0.0, 1.0), (2.0, 3.0)), grids, strict=True):
        rows = table[np.isin(table["module"], modules)]
        wavelength = np.asarray(rows["wavelength"], dtype=float)
        assert np.array_equal(wavelength, grid)
        unique, first = np.unique(wavelength, return_index=True)
        flux = np.asarray(rows["flux"], dtype=float)[first]
        error = np.asarray(rows["error (RMS+SYS)"], dtype=float)[first]
        chunks.append((unique, flux, error))
    return chunks


def observe(
    model: Model,
    instruments: dict[str, Instrument],
    truth: dict[str, float],
    calibration: dict[str, dict[str, float]] | None = None,
) -> dict[str, Any]:
    """Each instrument's noiseless output of *model* at *truth* (model compiled first)."""
    compiled = model.compile_for(negotiate(list(instruments.values())))
    produced = compiled(**truth)
    calibration = calibration or {}
    return {
        label: instrument(produced, calibration[label])
        if label in calibration
        else instrument(produced)
        for label, instrument in instruments.items()
    }


def with_values(container: Any, values: np.ndarray, sigma: np.ndarray) -> Any:
    """A copy of an observed-shape *container* carrying *values* +- *sigma*, Jy."""
    if isinstance(container, PhotometricPoints):
        return PhotometricPoints(
            container.filters,
            container.spectral_axis.values * u.um,
            values * u.Jy,
            uncertainty=sigma * u.Jy,
        )
    return Spectrum(container.spectral_axis.values * u.um, values * u.Jy, uncertainty=sigma * u.Jy)
