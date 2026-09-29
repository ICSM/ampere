"""The Kemper et al. (2002) NGC 6302 two-shell dust model -- W6.13 (2).

The v2 twin of ``examples/NGC6302.py`` (emcee) and ``examples/NGC6302_zeus.py``
(zeus); the legacy files stay exactly as they are (this module does **not**
touch them, nor the tracked data under ``examples/NGC6302/`` and
``examples/NGC6302-opacities.txt``). This first commit lands the model
itself, :class:`KemperTwoShell`, and its equality row against the legacy
``SpectrumNGC6302``; the data, the problem, the engines and the CLI follow
in the next commit.

The model
---------
:class:`KemperTwoShell` transcribes ``SpectrumNGC6302.__call__``,
``ckmodbb`` and ``shbb`` exactly (equations 6 and 7 of Kemper et al. 2002,
A&A 394, 679) -- checked bit-for-bit (``rtol=1e-10``) against
``examples.NGC6302.SpectrumNGC6302`` in the smoke test, vectorised over the
fourteen temperature integration steps (``ckmodbb``'s own loop runs
``range(steps - 1)`` with ``steps=15``, i.e. fourteen terms -- kept exactly,
not rounded up to fifteen) rather than looped in Python. Eight opacity
tables (:data:`~examples.ngc6302.generators.SPECIES`, the tracked files
named in ``examples/NGC6302-opacities.txt``) are loaded once, per-species,
as sixteen buffers (a native wavelength grid and an opacity array per
species) resolved by walking up from ``__file__``
(:data:`~examples.ngc6302.generators.OPACITY_DIRECTORY`) -- see
:mod:`.generators`'s module docstring for why not ``importlib.resources``.
The two unit conversions legacy applies (``examples/NGC6302.py`` lines
74-82) are folded into the buffers at load time, before interpolation:
enstatite and diopside are multiplied by ``1e-4`` (Q/a -> Q), calcite and
dolomite by ``density * 4/3 * 1e-4`` (M.A.C. -> Q, densities 2.71 and 2.87 g
cm\\ :sup:`-3`). Applying the conversion before rather than after the
log-space interpolation legacy performs is exact (a multiplicative factor is
an additive shift in log space, and interpolation commutes with an additive
shift), not an approximation -- confirmed by the equality test.

Interpolation onto whichever grid is actually being evaluated (the
constructor's own fallback grid before negotiation, or the negotiated grid
:meth:`~KemperTwoShell.compile_for` adopts) happens in
:func:`scipy.interpolate.interp1d` on ``log10`` of each species' tabulated
opacity, ``fill_value="extrapolate"`` -- legacy's own recipe, because the
tables do not all cover the full 2.36-196.6 micron range. ``compile_for``
interpolates once and caches the result, so a fit's hot loop pays for this
once, not per posterior draw.

The fifteen parameters and their priors
----------------------------------------
Legacy's own ``lims`` array (``examples/NGC6302.py`` lines 99-110, confirmed
by printing it: ``[-6, 0]`` for the eleven log-abundances, ``[10, 80]`` for
``Tcold0``/``Tcold1``, ``[80, 180]`` for ``Twarm0``/``Twarm1``) is
reproduced as eleven ``st.uniform(-6.0, 6.0)`` priors (seven cold --
species 0-4, 6, 7 -- and four warm -- species 1, 2, 5, 7) plus four
``st.uniform(10.0, 70.0)``/``st.uniform(80.0, 100.0)`` priors.

**Temperature ordering.** Legacy's own ``lnprior`` (``examples/NGC6302.py``
lines 257-274) does enforce ``Tcold1 > Tcold0`` and ``Twarm1 > Twarm0`` --
but as a hard rejection *on top of* the identical box for both parameters,
giving a flat joint density over the ordered triangle
``{10 <= Tcold0 < Tcold1 <= 80}`` (and likewise for the warm pair). v2's
:class:`~ampere.core.Parameter` priors are independent per parameter, with
one conditioning mechanism (:class:`~ampere.core.parameter.HierarchicalPrior`)
for a prior whose distribution parameters are themselves other parameters'
*values*; there is no way to express "loc a, scale (80 - a)" through it (the
mapping is a name reference, not an arithmetic expression), so a
hierarchical ``Tcold1 | Tcold0 ~ Uniform(Tcold0, 80)`` is not directly
constructible, and the nearest one that is (a fixed ``scale``) would break
the "same box" the ruling asks to keep. Per ruling 1's own escape hatch for
exactly this situation, this twin gives ``Tcold0``/``Tcold1`` (and
``Twarm0``/``Twarm1``) the **same, independent** box legacy's ``lims``
literally declares for both, named the same way (``Tcold0`` the outer/cooler
radius, ``Tcold1`` the inner/hotter one -- legacy's own convention). The
prior itself does not structurally forbid ``Tcold0 > Tcold1``; the forward
model does not silently reorder a disordered draw either -- ``Tcold1`` is
always passed as ``ckmodbb``'s ``tin`` and ``Tcold0`` as its ``tout``,
exactly as legacy calls it, so a disordered draw is scored by the same
physics, not given a free pass -- and the likelihood is what does the
discriminating (the two shells' data-implied temperatures are far enough
apart, and the model's ``tin``-referenced power law is not symmetric under
exchange, that the posterior is not expected to reward disorder). The
coverage run (a later commit) is the check that this in fact works.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st
from scipy.interpolate import interp1d

from ampere.core import Model, ModelResult, Parameter, Spectrum

from . import generators

__all__ = [
    "DEFAULT_GRID",
    "KemperTwoShell",
    "build_model",
]

# Legacy's own wavelength grid (``examples/NGC6302.py`` lines 295-297): a
# dense synthetic grid, not the observed data's own sampling -- the model's
# fallback grid before an instrument negotiates the real one.
_WAVE1 = np.linspace(2.3603, 35.0603, 327)
_WAVE2 = np.linspace(1.0 / 196.6261, 1.0 / 35.1, 117)
DEFAULT_GRID = np.concatenate((_WAVE1, 1.0 / _WAVE2[::-1]))

#: Fixed in ``ckmodbb`` (``examples/NGC6302.py`` line 203), never fitted.
_INDEX = 0.5
_RADIUS0_CM = 1e15
_DISTANCE_PC = 910.0
_GRAIN_RADIUS_UM = 0.1
_STEPS = 15


def _as_parameter(name: str, spec: Any, *, unit: u.UnitBase | None = None) -> Parameter:
    """A prior (a frozen scipy distribution) or a bare number, as a Parameter."""
    if isinstance(spec, Parameter):
        return spec
    if hasattr(spec, "logpdf") or hasattr(spec, "rvs"):
        return Parameter(name, spec, unit=unit)
    return Parameter(name, None, value=float(spec), fixed=True, unit=unit)


def _shbb(grid: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """``shbb`` (``examples/NGC6302.py`` lines 235-248), ``pinda=0`` always.

    Vectorised over *temperature* (shape ``(steps,)``) against *grid* (shape
    ``(n_wave,)``); returns shape ``(steps, n_wave)``. With ``pinda=0`` the
    ``mbb[1, :] = bbflux * wl ** pinda`` line is just ``bbflux`` (``wl ** 0
    == 1``), so it is not separately computed.
    """
    a1 = 3.97296e19
    a2 = 1.43875e4
    wavelength = grid[None, :]
    temp = temperature[:, None]
    return a1 / (wavelength**3) / (np.exp(a2 / (wavelength * temp)) - 1.0)


def _ckmodbb(
    opacity: np.ndarray,
    *,
    tin: float,
    tout: float,
    n0: float,
    grid: np.ndarray,
) -> np.ndarray:
    """``ckmodbb`` (``examples/NGC6302.py`` lines 202-233), equations 6-7 of
    Kemper et al. (2002), vectorised over the fourteen temperature steps.
    """
    distance_cm = _DISTANCE_PC * 3.0857e18
    grain_radius_cm = _GRAIN_RADIUS_UM * 1e-4
    pindex = qindex = _INDEX

    step_index = np.arange(_STEPS - 1, dtype=float)
    temperature = tin - step_index * (tin - tout) / _STEPS
    power = (temperature / tin) ** (-(3.0 - pindex) / qindex)
    weight = power * ((tin - tout) / _STEPS)
    blackbody = _shbb(grid, temperature)
    fnu = np.sum(opacity[None, :] * blackbody * weight[:, None], axis=0)

    factor = (4.0 * math.pi * grain_radius_cm**2 * _RADIUS0_CM**3 * n0) / (
        (3.0 - pindex) * distance_cm**2
    )
    return fnu * factor


class KemperTwoShell(Model):
    """The Kemper et al. (2002) two-shell dust model -- see the module docstring.

    Parameters
    ----------
    wavelength
        The model's own fallback grid, micron (negotiation replaces it the
        moment an instrument is bound; see :meth:`compile_for`).
    logacold0, logacold1, logacold2, logacold3, logacold4, logacold6, logacold7
        Cold-component log10 abundance, species 0-4, 6, 7 (calcite,
        enstatite, forsterite, diopside, dolomite, ice, olivine). A prior to
        fit each, or a number to hold it fixed.
    logawarm1, logawarm2, logawarm5, logawarm7
        Warm-component log10 abundance, species 1, 2, 5, 7 (enstatite,
        forsterite, iron, olivine).
    Tcold0, Tcold1, Twarm0, Twarm1
        Shell temperatures, kelvin -- see the module docstring's note on the
        ordering these are *not* structurally constrained to keep.
    channel
        Name of the channel the emitted :class:`~ampere.core.Spectrum`
        appears under.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        logacold0: Any,
        logacold1: Any,
        logacold2: Any,
        logacold3: Any,
        logacold4: Any,
        logacold6: Any,
        logacold7: Any,
        logawarm1: Any,
        logawarm2: Any,
        logawarm5: Any,
        logawarm7: Any,
        Tcold0: Any,
        Tcold1: Any,
        Twarm0: Any,
        Twarm1: Any,
        channel: str = "sed",
    ) -> None:
        grid = np.asarray(wavelength, dtype=float)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(
                f"KemperTwoShell needs a 1-D, non-empty grid, got shape {grid.shape}."
            )
        self.channel = str(channel)
        self.register_buffer("wavelength", grid, unit=u.um)

        self._opacity_tables: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        filenames = np.loadtxt(generators.OPACITY_FILE_LIST, dtype=str)
        for name, filename in zip(generators.SPECIES, filenames, strict=True):
            table_wavelength, table_opacity = np.loadtxt(
                generators.OPACITY_DIRECTORY / str(filename), comments="#", unpack=True
            )
            # The two unit conversions, applied once, before interpolation
            # (exact -- see the module docstring).
            if name == "calcite":
                table_opacity = table_opacity * 2.71 * (4.0 / 3.0) * 1e-4
            elif name in ("enstatite", "diopside"):
                table_opacity = table_opacity * 1e-4
            elif name == "dolomite":
                table_opacity = table_opacity * 2.87 * (4.0 / 3.0) * 1e-4
            wavelength_buffer = self.register_buffer(f"{name}_wavelength", table_wavelength, unit=u.um)
            opacity_buffer = self.register_buffer(
                f"{name}_opacity", table_opacity, unit=u.dimensionless_unscaled
            )
            self._opacity_tables[name] = (wavelength_buffer, opacity_buffer)

        self.register_parameter(_as_parameter("logacold0", logacold0))
        self.register_parameter(_as_parameter("logacold1", logacold1))
        self.register_parameter(_as_parameter("logacold2", logacold2))
        self.register_parameter(_as_parameter("logacold3", logacold3))
        self.register_parameter(_as_parameter("logacold4", logacold4))
        self.register_parameter(_as_parameter("logacold6", logacold6))
        self.register_parameter(_as_parameter("logacold7", logacold7))
        self.register_parameter(_as_parameter("logawarm1", logawarm1))
        self.register_parameter(_as_parameter("logawarm2", logawarm2))
        self.register_parameter(_as_parameter("logawarm5", logawarm5))
        self.register_parameter(_as_parameter("logawarm7", logawarm7))
        self.register_parameter(_as_parameter("Tcold0", Tcold0, unit=u.K))
        self.register_parameter(_as_parameter("Tcold1", Tcold1, unit=u.K))
        self.register_parameter(_as_parameter("Twarm0", Twarm0, unit=u.K))
        self.register_parameter(_as_parameter("Twarm1", Twarm1, unit=u.K))

        self._template: Spectrum | None = None
        self._opacity_on_grid: np.ndarray | None = None

    def _interpolate_opacity(self, grid: np.ndarray) -> np.ndarray:
        """Each species' opacity onto *grid*, log-space, legacy's own recipe.

        ``examples/NGC6302.py`` lines 52-56: ``interp1d`` on ``log10`` of the
        tabulated opacity, ``fill_value="extrapolate"``, because the tables
        do not all cover 2.36-196.6 micron. Returns shape ``(grid.size, 8)``.
        """
        columns = []
        for name in generators.SPECIES:
            table_wavelength, table_opacity = self._opacity_tables[name]
            spline = interp1d(
                table_wavelength,
                np.log10(table_opacity),
                assume_sorted=False,
                fill_value="extrapolate",
            )
            columns.append(10.0 ** spline(grid))
        return np.stack(columns, axis=1)

    def compile_for(self, requirements: Any) -> KemperTwoShell:
        """Adopt the negotiated grid for :attr:`channel`, once (W2.12's contract).

        Also re-interpolates the opacity tables onto that grid and caches
        the result, so the fit's hot loop (:meth:`evaluate`) pays for the
        interpolation once, not per posterior draw.
        """
        asked = requirements.get(self.channel)
        if asked is not None and "spectral_axis" in asked:
            grid = np.asarray(asked["spectral_axis"].coordinates(), dtype=float)
            self._template = Spectrum(grid * u.um, np.zeros(grid.size), unit=u.Jy)
            self._opacity_on_grid = self._interpolate_opacity(grid)
        return self

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        if self._template is not None:
            grid = self._template.spectral_axis.values
            opacity = self._opacity_on_grid
        else:
            grid = np.asarray(ctx["wavelength"], dtype=float)
            opacity = self._interpolate_opacity(grid)

        # acold/awarm -- examples/NGC6302.py lines 158-178. Species 5 (iron)
        # is never fitted cold; species 0, 3, 4, 6 (calcite, diopside,
        # dolomite, ice) are never fitted warm. Every one of the eight
        # species is still summed (n0 = 0 for the excluded ones), exactly as
        # legacy's own loop runs over all eight unconditionally.
        acold = [
            10.0 ** float(ctx["logacold0"]),
            10.0 ** float(ctx["logacold1"]),
            10.0 ** float(ctx["logacold2"]),
            10.0 ** float(ctx["logacold3"]),
            10.0 ** float(ctx["logacold4"]),
            0.0,
            10.0 ** float(ctx["logacold6"]),
            10.0 ** float(ctx["logacold7"]),
        ]
        awarm = [
            0.0,
            10.0 ** float(ctx["logawarm1"]),
            10.0 ** float(ctx["logawarm2"]),
            0.0,
            0.0,
            10.0 ** float(ctx["logawarm5"]),
            0.0,
            10.0 ** float(ctx["logawarm7"]),
        ]

        cold = np.zeros(grid.size)
        for index, n0 in enumerate(acold):
            cold = cold + _ckmodbb(
                opacity[:, index],
                tin=float(ctx["Tcold1"]),
                tout=float(ctx["Tcold0"]),
                n0=n0,
                grid=grid,
            )
        warm = np.zeros(grid.size)
        for index, n0 in enumerate(awarm):
            warm = warm + _ckmodbb(
                opacity[:, index],
                tin=float(ctx["Twarm1"]),
                tout=float(ctx["Twarm0"]),
                n0=n0,
                grid=grid,
            )
        flux = cold + warm

        spectrum = (
            self._template.with_values(flux)
            if self._template is not None
            else Spectrum(grid * u.um, flux, unit=u.Jy)
        )
        return ModelResult({self.channel: spectrum})


def build_model() -> KemperTwoShell:
    """The one :class:`KemperTwoShell`, with legacy's own box priors."""
    abundance_prior = st.uniform(-6.0, 6.0)
    return KemperTwoShell(
        DEFAULT_GRID,
        logacold0=abundance_prior,
        logacold1=abundance_prior,
        logacold2=abundance_prior,
        logacold3=abundance_prior,
        logacold4=abundance_prior,
        logacold6=abundance_prior,
        logacold7=abundance_prior,
        logawarm1=abundance_prior,
        logawarm2=abundance_prior,
        logawarm5=abundance_prior,
        logawarm7=abundance_prior,
        Tcold0=st.uniform(10.0, 70.0),
        Tcold1=st.uniform(10.0, 70.0),
        Twarm0=st.uniform(80.0, 100.0),
        Twarm1=st.uniform(80.0, 100.0),
    )
