"""The data for ``examples.cstar``: the tracked photometry, IRS spectrum and photosphere.

Everything here reads ``examples/cstar_data/`` (tracked, read-only for this
package) with :mod:`astropy` and :mod:`numpy` alone -- no legacy import.

* :func:`read_votable` reads ``Observed_SED.vot``: eleven points, ``MCPS_B``
  to ``SPITZER_MIPS_24``, every one a filter in ampere's bundled library. The
  file's flux columns are declared in **mJy** (the legacy reader passed them
  through); they are converted to Jy here, which is what the model emits.
* :func:`read_irs` reads ``SPEC_OGLE_CAGB_IRS.csv`` (a ``#``-prefixed header
  line, 364 rows) into **one chunk per value of its ``chunk`` column** -- the
  legacy reader's split (``ampere/legacy/data/spectrum.py`` lines 866-873,
  ``for sp in s`` in the script): chunk 1 (5.13-14.29 micron, 192 points) and
  chunk 2 (14.24-36.85 micron, 172 points). The legacy stored each chunk in
  file order, which interleaves the IRS orders and is **not monotonic** (the
  SL2/SL1 and LL2/LL1 orders are listed one after the other). A v2
  :class:`~ampere.core.Spectrum` must be strictly increasing, so each chunk
  is sorted and then deduplicated with :func:`numpy.unique` -- the pattern of
  :func:`examples.linear_sed.generators.deduplicated_grids`, **copied rather
  than imported**, since that function reads a CASSIS FITS file and this is a
  CSV. Neither chunk actually repeats a wavelength (the two chunks overlap in
  range, 14.24-14.29 micron, but they are separate spectra), so the
  deduplication is a guard, not a change.
* :func:`read_photosphere` reads ``photosphere_interpolated.csv`` (``nu`` Hz
  descending, ``F_nu``, no header) **exactly as the legacy does**:
  ``np.loadtxt(..., delimiter=",", skiprows=1)``, which drops the file's
  first data row (the legacy ``skiprows=1`` assumes a header the file does
  not have). The twin keeps that, so both hand Hyperion the same 300 rows.
* :func:`synthetic_observations` is the ``--synthetic`` mode: one simulation
  at :data:`SYNTHETIC_TRUTH`, through every instrument, plus Gaussian noise
  at the data's own uncertainties.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
from astropy.io.votable import parse_single_table
from astropy.table import Table

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate

__all__ = [
    "DATA_DIR",
    "IRS_FILE",
    "PHOTOSPHERE_FILE",
    "SEED",
    "SYNTHETIC_TRUTH",
    "VOTABLE",
    "read_irs",
    "read_photosphere",
    "read_votable",
    "synthetic_observations",
    "with_values",
]

SEED = 20260930

#: The legacy's tracked data (never written by this package).
DATA_DIR = Path(__file__).resolve().parents[1] / "cstar_data"
VOTABLE = DATA_DIR / "Observed_SED.vot"
IRS_FILE = DATA_DIR / "SPEC_OGLE_CAGB_IRS.csv"
PHOTOSPHERE_FILE = DATA_DIR / "photosphere_interpolated.csv"

#: ``--synthetic``'s truth: the legacy ``__init__`` defaults, by the twin's
#: names. ``envelope_mass = log10(6.985718e-6)``; ``envelope_rin`` and
#: ``envelope_r0 = log10(4.4859)`` stellar radii; ``envelope_rout = 3.65``,
#: placed inside the prior box (the legacy default, 1000 *inner* radii, is
#: log10(4486) = 3.65 stellar radii, so this is that default); the star's
#: legacy mass and luminosity; ``sic_fraction = 0.1``.
SYNTHETIC_TRUTH: dict[str, float] = {
    "envelope_mass": float(np.log10(6.985718e-6)),
    "envelope_rin": float(np.log10(4.4859)),
    "envelope_rout": 3.65,
    "envelope_r0": float(np.log10(4.4859)),
    "stellar_mass": 2.0,
    "stellar_luminosity": 6165.95,
    "sic_fraction": 0.1,
}


def read_votable(path: Path | str = VOTABLE) -> tuple[list[str], np.ndarray, np.ndarray]:
    """``(filter names, flux Jy, uncertainty Jy)`` from the votable, in file order."""
    table = parse_single_table(str(path)).to_table()
    names = [str(name).strip() for name in table["sed_filter"]]
    flux = np.asarray(table["sed_flux"].to(u.Jy).value, dtype=float)
    error = np.asarray(table["sed_eflux"].to(u.Jy).value, dtype=float)
    return names, flux, error


def read_irs(path: Path | str = IRS_FILE) -> dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """``{chunk: (wavelength micron, flux Jy, uncertainty Jy)}``, sorted and deduplicated.

    At a repeated wavelength the first sample in sorted order is kept.
    """
    table = Table.read(
        path,
        format="ascii.csv",
        comment="#",
        names=["globalsourceid", "aorkey", "wavelength", "fluxdensity", "uncertainty", "chunk"],
    )
    chunks: dict[int, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for chunk in np.unique(np.asarray(table["chunk"], dtype=int)):
        rows = table[np.asarray(table["chunk"], dtype=int) == chunk]
        wavelength = np.asarray(rows["wavelength"], dtype=float)
        order = np.argsort(wavelength, kind="stable")
        unique, first = np.unique(wavelength[order], return_index=True)
        flux = np.asarray(rows["fluxdensity"], dtype=float)[order][first]
        error = np.asarray(rows["uncertainty"], dtype=float)[order][first]
        chunks[int(chunk)] = (unique, flux, error)
    return chunks


def read_photosphere(path: Path | str = PHOTOSPHERE_FILE) -> tuple[np.ndarray, np.ndarray]:
    """``(nu Hz, F_nu)`` as the legacy reads it -- ``skiprows=1`` drops the first row."""
    nu, fnu = np.loadtxt(path, delimiter=",", skiprows=1, unpack=True)
    return nu, fnu


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


def synthetic_observations(
    model: Model,
    instruments: dict[str, Instrument],
    sigmas: dict[str, np.ndarray],
    *,
    truth: dict[str, float] | None = None,
    seed: int = SEED,
) -> dict[str, Any]:
    """One simulation at *truth* through each instrument, plus noise at *sigmas*.

    Calibration factors are held at 1.0 (no injected miscalibration). *model*
    runs Hyperion once here, at whatever photon preset it was built with.
    """
    rng = np.random.default_rng(seed)
    compiled = model.compile_for(negotiate(list(instruments.values())))
    produced = compiled(**(SYNTHETIC_TRUTH if truth is None else truth))
    unity = {"calibration_scale.scale": 1.0}
    observed: dict[str, Any] = {}
    for label, instrument in instruments.items():
        shape = instrument(produced, unity) if label != "photometry" else instrument(produced)
        sigma = np.asarray(sigmas[label], dtype=float)
        observed[label] = with_values(shape, shape.values + rng.normal(0.0, sigma), sigma)
    return observed
