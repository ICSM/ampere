"""The carbon star's dust without ``bhmie``: Mie opacities from ``miepython``.

The legacy ``HyperionCStarRTModel.__call__`` (``ampere/legacy/models/Hyperion.py``
lines 515-575) wrote a size table per species, a ``bhmie`` parameter file, and
ran the unpackaged ``bhmie`` Fortran program on it, then read the result back
with Hyperion's ``BHDust``. Ruled 2026-09-28 on
``docs/design/example_dependencies_memo.md`` §3, that step is replaced here by
``miepython`` (conda-forge, in the pixi ``hyperion`` feature). What is kept
from the legacy exactly:

* **The size distribution** (:func:`size_distribution`): ``sizes =
  linspace(amin, amax, na)`` with ``amin = 0.01``, ``amax = 1.0`` micron,
  ``na = 101``, and ``numbers = sizes**q * exp(-sizes / abrk)`` with ``q =
  -3.5`` and ``abrk = 1.0`` micron -- the legacy's own ``.size`` recipe, whose
  tracked output (``examples/cstar_data/*.size``) the smoke tests compare
  against to 1e-6.
* **The materials and densities**: amorphous carbon (Rouleau & Martin 1991,
  ``rouleau91_ac.optc``, 1.80 g cm^-3) and silicon carbide (Pégourié 1988,
  ``SiC_Pegourie1988.optc``, 3.22 g cm^-3). ``carbon="zubko96"`` swaps in the
  tracked Zubko et al. (1996) ACAR constants at the same density.
* **The opacity wavelength grid**: ``lmin, lmax, nl = 0.1, 200, 101``, as
  :data:`OPACITY_WAVELENGTH` (log-spaced, as ``bhmie``'s own grid is).
* **The Hyperion post-processing**: ``optical_properties.extrapolate_wav(0.05,
  1200)`` and ``set_lte_emissivities(n_temp=101, temp_min=2.7,
  temp_max=2000)``.

What changes
------------
**The per-species tables.** For each species of grain density rho, the tabulated
``n`` and ``k`` are interpolated in log λ onto :data:`OPACITY_WAVELENGTH`
(extrapolated flat beyond a table's ends, which never happens on this grid:
all three tables span 0.1-200 micron), ``miepython.efficiencies_mx(m, x)``
gives ``Q_ext``, ``Q_sca`` and ``g`` for every (λ, a) with ``m = n - ik``
(miepython's convention; see below) and ``x = 2πa/λ``, and the size
distribution is integrated with the trapezoid rule::

    κ_ext(λ) = ∫ πa² Q_ext n(a) da / ∫ (4/3)πa³ rho n(a) da      [cm² g⁻¹]
    κ_sca(λ) = ∫ πa² Q_sca n(a) da / ∫ (4/3)πa³ rho n(a) da
    g(λ)     = ∫ πa² Q_sca g n(a) da / ∫ πa² Q_sca n(a) da

**The mixture** at ``sic_fraction = f`` is mass-weighted, which is what the
legacy's ``bhmie`` mass fractions meant::

    κ      = (1 - f) κ_ext,C + f κ_ext,SiC
    albedo = ((1 - f) κ_sca,C + f κ_sca,SiC) / κ
    g      = ((1 - f) κ_sca,C g_C + f κ_sca,SiC g_SiC) / (κ albedo)

It is linear in the per-species tables, so the Mie computation -- the only
expensive part -- runs **once per process** at first use and is cached on
the module (:func:`species_tables`, a :func:`functools.cache` keyed on the
species); each simulation only mixes.

**The phase function.** ``bhmie`` tabulated the full Mie scattering matrix,
which ``BHDust`` hands to Hyperion. miepython gives the asymmetry parameter
``g`` per wavelength, and the array-form counterpart of that is
:class:`hyperion.dust.HenyeyGreensteinDust` (with ``p_lin_max = 0``, no linear
polarisation) -- so the twin scatters with a Henyey-Greenstein phase function
of the right ``g``, not ``IsotropicDust``, which would throw ``g`` away.
The difference matters only for the scattered-light share of the optical
SED; the infrared, which the photometry and the IRS spectrum constrain, is
thermal emission set by ``κ``.

**The refractive-index sign.** miepython's documented convention is ``m = n
- ik`` for an absorber. (miepython 3.3 in fact returns the same efficiencies
for ``n + ik``, so the sign cannot silently make an absorber transparent in
this version; the smoke test's SiC-feature and carbon-slope rows are the
guard should a later version stop being forgiving.)

**Hyperion 0.9.11 and NumPy 2.** Hyperion's Python front end writes every
HDF5 string attribute with ``np.string_``, which NumPy 2.0 removed (it was an
exact alias of ``np.bytes_``); conda-forge's 0.9.11 build nevertheless
declares ``numpy <3``, so in the ``hyperion`` environment (NumPy 2) the first
``model.write`` fails. :func:`hyperion_numpy_alias` restores that one alias,
``np.string_ = np.bytes_``, before this package imports Hyperion -- the only
removed NumPy name the 0.9.11 front end uses -- and does nothing on a NumPy
that still has it. It is scoped to the process that imports Hyperion through
this package, and it is to be dropped when upstream's 2026 release (which
the pin in ``pyproject.toml`` waits for) lands.
"""

from __future__ import annotations

import functools
from pathlib import Path
from typing import Any

import numpy as np

__all__ = [
    "CARBON",
    "DATA_DIR",
    "OPACITY_WAVELENGTH",
    "SIC",
    "SIZES",
    "SIZE_Q",
    "SIZE_SCALE",
    "Species",
    "build_dust",
    "hyperion_numpy_alias",
    "mixture",
    "read_optical_constants",
    "size_distribution",
    "species_opacity",
    "species_tables",
]

#: The legacy's tracked data directory (read, never written).
DATA_DIR = Path(__file__).resolve().parents[1] / "cstar_data"

#: ``amin=0.01, amax=1.0, na=101`` micron (legacy ``__init__``).
SIZES = np.linspace(0.01, 1.0, 101)
#: ``q=-3.5``.
SIZE_Q = -3.5
#: ``abrk=1.0`` micron: the exponential cut-off scale.
SIZE_SCALE = 1.0

#: ``lmin, lmax, nl = 0.1, 200.0, 101`` micron.
OPACITY_WAVELENGTH = np.geomspace(0.1, 200.0, 101)

#: The Hyperion post-processing, legacy values.
EXTRAPOLATE_WAV = (0.05, 1200.0)
LTE_TEMPERATURES = {"n_temp": 101, "temp_min": 2.7, "temp_max": 2000.0}


class Species:
    """A grain material: its optical-constants file and its bulk density."""

    def __init__(self, name: str, filename: str, density: float) -> None:
        self.name = name
        self.path = DATA_DIR / filename
        self.density = float(density)

    def __repr__(self) -> str:
        return f"Species({self.name!r}, {self.path.name!r}, density={self.density})"


#: The amorphous-carbon choices (``--carbon``), the legacy's Rouleau default first.
CARBON: dict[str, Species] = {
    "rouleau91": Species("rouleau91", "rouleau91_ac.optc", 1.80),
    "zubko96": Species("zubko96", "zubko96_ac_acar.optc", 1.80),
}
#: Silicon carbide (Pégourié 1988).
SIC = Species("sic", "SiC_Pegourie1988.optc", 3.22)


def hyperion_numpy_alias() -> None:
    """Restore ``np.string_`` (= ``np.bytes_``) for Hyperion 0.9.11 on NumPy 2; see above."""
    if "string_" not in np.__dict__:
        np.string_ = np.bytes_  # type: ignore[attr-defined]


def size_distribution(
    sizes: np.ndarray = SIZES, q: float = SIZE_Q, scale: float = SIZE_SCALE
) -> tuple[np.ndarray, np.ndarray]:
    """``(sizes, numbers)``, micron: the legacy ``.size`` recipe, in memory."""
    sizes = np.asarray(sizes, dtype=float)
    return sizes, sizes**q * np.exp(-sizes / scale)


def read_optical_constants(path: Path | str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(wavelength micron, n, k)`` from a three-column ``.optc`` table, ascending."""
    table = np.loadtxt(path)
    order = np.argsort(table[:, 0])
    return table[order, 0], table[order, 1], table[order, 2]


def species_opacity(
    species: Species, wavelength: np.ndarray = OPACITY_WAVELENGTH
) -> dict[str, np.ndarray]:
    """Mass opacities ``kext``, ``ksca`` (cm² g⁻¹) and asymmetry ``g`` on *wavelength*.

    The Mie sum for one species over the legacy size distribution; see the
    module docstring for the integrals. Needs ``miepython``.
    """
    import miepython

    wavelength = np.asarray(wavelength, dtype=float)
    table_wav, table_n, table_k = read_optical_constants(species.path)
    # Log-lambda interpolation; np.interp holds the end values flat beyond the table.
    log_wav = np.log(wavelength)
    n = np.interp(log_wav, np.log(table_wav), table_n)
    k = np.interp(log_wav, np.log(table_wav), table_k)
    sizes, numbers = size_distribution()

    m = np.repeat((n - 1j * k)[:, None], sizes.size, axis=1)
    x = 2.0 * np.pi * sizes[None, :] / wavelength[:, None]
    qext, qsca, _qback, g = miepython.efficiencies_mx(m.ravel(), x.ravel())
    shape = (wavelength.size, sizes.size)
    qext = np.asarray(qext, dtype=float).reshape(shape)
    qsca = np.asarray(qsca, dtype=float).reshape(shape)
    g = np.asarray(g, dtype=float).reshape(shape)

    a_cm = sizes * 1.0e-4
    area = np.pi * a_cm**2 * numbers
    mass = np.trapezoid(4.0 / 3.0 * np.pi * a_cm**3 * species.density * numbers, sizes)
    kext = np.trapezoid(area[None, :] * qext, sizes, axis=1) / mass
    sca = np.trapezoid(area[None, :] * qsca, sizes, axis=1)
    ksca = sca / mass
    asym = np.trapezoid(area[None, :] * qsca * g, sizes, axis=1) / sca
    return {"wavelength": wavelength, "kext": kext, "ksca": ksca, "g": asym}


@functools.cache
def species_tables(
    carbon: str = "rouleau91",
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """``(carbon, sic)`` opacity tables, computed once per process and cached."""
    if carbon not in CARBON:
        raise ValueError(f"unknown carbon {carbon!r}; choose one of {sorted(CARBON)}.")
    return species_opacity(CARBON[carbon]), species_opacity(SIC)


def mixture(sic_fraction: float, carbon: str = "rouleau91") -> dict[str, np.ndarray]:
    """The mass-weighted mixture: ``wavelength``, ``chi`` (cm² g⁻¹), ``albedo``, ``g``."""
    f = float(sic_fraction)
    if not 0.0 <= f <= 1.0:
        raise ValueError(f"sic_fraction must lie in [0, 1], got {f}.")
    c, s = species_tables(carbon)
    chi = (1.0 - f) * c["kext"] + f * s["kext"]
    scattering = (1.0 - f) * c["ksca"] + f * s["ksca"]
    albedo = scattering / chi
    g = ((1.0 - f) * c["ksca"] * c["g"] + f * s["ksca"] * s["g"]) / scattering
    return {"wavelength": c["wavelength"], "chi": chi, "albedo": albedo, "g": g}


def build_dust(sic_fraction: float, carbon: str = "rouleau91") -> Any:
    """A Hyperion ``HenyeyGreensteinDust`` for the mixture, legacy post-processing applied.

    Needs ``hyperion``. Hyperion wants frequency ascending, so the arrays are
    reversed from the wavelength-ascending tables.
    """
    hyperion_numpy_alias()
    from hyperion.dust import HenyeyGreensteinDust

    mixed = mixture(sic_fraction, carbon)
    nu = 2.99792458e14 / mixed["wavelength"]  # micron -> Hz
    order = np.argsort(nu)
    dust = HenyeyGreensteinDust(
        nu[order], mixed["albedo"][order], mixed["chi"][order], mixed["g"][order], 0.0
    )
    dust.optical_properties.extrapolate_wav(*EXTRAPOLATE_WAV)
    dust.set_lte_emissivities(**LTE_TEMPERATURES)
    return dust
