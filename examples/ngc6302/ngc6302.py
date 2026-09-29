"""The Kemper et al. (2002) NGC 6302 two-shell dust model -- W6.13 (2).

The v2 twin of ``examples/NGC6302.py`` (emcee) and ``examples/NGC6302_zeus.py``
(zeus); the legacy files stay exactly as they are (this module does **not**
touch them, nor the tracked data under ``examples/NGC6302/`` and
``examples/NGC6302-opacities.txt``). One model, the same data, the same
question, written once against :mod:`ampere.core` and the reference backend,
with the sampler chosen by ``--engine`` rather than by which file you ran.

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
species 0-4, 6, 7 -- and four warm -- species 1, 2, 5, 7) plus the four
temperature-family parameters below.

**Temperature ordering.** Legacy's own ``lnprior`` (``examples/NGC6302.py``
lines 257-274) enforces ``Tcold1 > Tcold0`` and ``Twarm1 > Twarm0`` as a
hard rejection *on top of* the identical box for both parameters, giving a
flat joint density over the ordered triangle
``{10 <= Tcold0 < Tcold1 <= 80}`` (and likewise for the warm pair, over
``{80 <= Twarm0 < Twarm1 <= 180}``). This twin reproduces that triangle
**exactly**, without a rejection step, via one derived quantity per pair
(Peter's ruling, 2026-09-29): ``Tcold0 ~ st.triang(c=0, loc=10.0,
scale=70.0)`` -- density proportional to ``80 - Tcold0`` on ``[10, 80]``,
confirmed 2026-09-29 -- and ``Tcold_fraction ~ st.uniform(0.0, 1.0)``,
independent, with ``Tcold1 = Tcold0 + Tcold_fraction * (80 - Tcold0)``
computed in :meth:`KemperTwoShell.evaluate` (and, identically,
``Twarm0 ~ st.triang(c=0, loc=80.0, scale=100.0)``,
``Twarm_fraction ~ st.uniform(0.0, 1.0)``,
``Twarm1 = Twarm0 + Twarm_fraction * (180 - Twarm0)``).

**Why it is exact.** The map ``(T0, f) -> (T0, T1)`` (fixing ``T0``,
``T1 = T0 + f * (80 - T0)``) has Jacobian ``dT1/df = 80 - T0``, so the joint
density of the pair transforms as
``p(T0, T1) = p(T0) p(f) / |dT1/df| = [k (80 - T0)] * 1 / (80 - T0) = k`` --
constant on the triangle ``{10 <= T0 < T1 <= 80}``, precisely legacy's own
flat ordered prior, reproduced through two independent priors and one
derived quantity rather than a rejection step.
:class:`~ampere.core.parameter.HierarchicalPrior` still cannot express this
directly (it binds a raw value to a keyword, not an arithmetic expression
like "loc a, scale 80 - a"), which is why the derived quantity lives in
:meth:`~KemperTwoShell.evaluate` instead. A row in
:mod:`tests.examples.test_ngc6302` checks the exactness numerically: at
1 000 random points of the triangle, ``log p(T0) + log p(f) - log(80 - T0)``
is constant to ``1e-10`` (and likewise for the warm pair, against ``180``).

The legacy names ``Tcold1``/``Twarm1`` survive as **derived quantities**
(:func:`derived_temperatures`), not fitted parameters: :func:`report` prints
all four physical temperatures from the posterior draws beside the declared
ones, and :func:`recovers_truth` scores the four declared temperature-family
parameters plus the two derived ones.

History: the first draft of this twin gave ``Tcold0``/``Tcold1`` (and
``Twarm0``/``Twarm1``) the same, independent, *unordered* box legacy's
``lims`` literally declares for both, reasoning that v2's per-parameter
priors have no direct way to express "loc a, scale (80 - a)"; the 900-step
coverage run below (see "History" in that section) showed that choice
lets the ensemble settle into disordered, data-compatible modes that never
recombine, which the exact reparameterisation above forecloses by
construction.

The data
--------
One instrument, ``"iso"`` (the real spectrum is ISO SWS/LWS) --
:class:`~ampere.backends.reference.Resample` onto the observed 25-120 micron
selection (:func:`~examples.ngc6302.generators.load_observed_spectrum`) then
:class:`~ampere.backends.reference.CalibrationScale`. Legacy's own
``calUnc=1e-10`` declared "no calibration freedom"; this twin gives it a
five per cent log-normal instead, as :mod:`examples.linear_sed` documents
doing for its own ``calUnc``. :class:`~ampere.core.GaussianProcessNoise`
with a :class:`~ampere.core.Matern32` kernel in wavelength is on by default
(``--no-gp`` swaps in :class:`~ampere.core.IndependentNoise`); the
length-scale prior reproduces legacy's own ``scalelengthPrior=0.1`` (micron,
half-normal -- ``examples/NGC6302.py`` line ~367) directly, and the
amplitude prior is weakly informative, scaled to *this* spectrum's own flux
magnitude (a few hundred Jy across the 25-120 micron window) rather than
:mod:`examples.linear_sed`'s O(1) Jy synthetic scale, which would be badly
mismatched here.

Engines
-------
``--engine emcee|zeus`` (:class:`~ampere.inference.EmceeEngine`,
:class:`~ampere.inference.ZeusEngine`), both at legacy's own 50 walkers and
its emcee script's own 50 000 steps / 40 000 burn-in as the shared full
budget (zeus's own legacy script uses a smaller 5 000 / 3 500 budget for a
quicker look; this twin keeps one full budget across both engines, as
:mod:`examples.linear_sed` and :mod:`examples.modified_blackbody` do, so a
caller comparing engines compares them at one cost). ``--quick`` uses 2 000
steps / 1 000 burn-in instead, for a fast look, not a coverage claim.
``NGC6302_zeus.py`` also drops ``logawarm1`` (a fourteen-parameter model);
this twin keeps all fifteen parameters on both engines and does not
reproduce that reduction -- one model, two samplers, as the item text asks.

Synthetic mode and coverage
----------------------------
``--synthetic`` builds :func:`~examples.ngc6302.generators.synthetic_data`
instead of the real spectrum: the model evaluated at the "2002 solution"
(:data:`examples.ngc6302.generators.TRUTH` -- see that module's docstring
for its source and the species correspondence), on the observed wavelength
selection, plus Gaussian noise at the data's own five-per-cent uncertainty
rule. :func:`recovers_truth` checks the fifteen model parameters plus the
calibration factor (``iso.instrument.calibration_scale.scale``, truth 1.0),
**and** the two derived temperatures ``Tcold1``/``Twarm1`` (ruling 2) -- the
GP's own hyperparameters have a prior but no injected truth, as in
:mod:`examples.linear_sed`.

::

    python -m examples.ngc6302 --synthetic --engine emcee --quick   # a look
    python -m examples.ngc6302 --synthetic --engine emcee           # the coverage run's engine

Coverage run (Accept criterion) -- every parameter covers truth, but the
chain still does not converge; reported as a finding (ruling 3, tranche B)
-----------------------------------------------------------------------------
``python -m examples.ngc6302 --synthetic --engine emcee --seed 20260928
--walkers 50 --steps 10000 --burn-in 5000`` (GP on, ``QuasisepGP``), run
2026-09-29, wall clock 4748.3 s (79.1 min). ``ampere.results.summary`` on
the result, all eighteen free parameters (R-hat, ESS (bulk))::

    iso.instrument.calibration_scale.scale  r_hat 1.70  ess_bulk 77
    iso.likelihood.amplitude                r_hat 2.13  ess_bulk 64
    iso.likelihood.length_scale             r_hat 1.69  ess_bulk 78
    model.Tcold0                            r_hat 2.41  ess_bulk 60
    model.Tcold_fraction                    r_hat 2.51  ess_bulk 59
    model.Twarm0                            r_hat 2.66  ess_bulk 58
    model.Twarm_fraction                    r_hat 2.79  ess_bulk 57
    model.logacold0                         r_hat 2.49  ess_bulk 60
    model.logacold1                         r_hat 2.45  ess_bulk 60
    model.logacold2                         r_hat 2.73  ess_bulk 58
    model.logacold3                         r_hat 2.46  ess_bulk 60
    model.logacold4                         r_hat 2.43  ess_bulk 60
    model.logacold6                         r_hat 2.50  ess_bulk 60
    model.logacold7                         r_hat 2.90  ess_bulk 57
    model.logawarm1                         r_hat 2.73  ess_bulk 58
    model.logawarm2                         r_hat 2.30  ess_bulk 62
    model.logawarm5                         r_hat 2.48  ess_bulk 60
    model.logawarm7                         r_hat 2.34  ess_bulk 61

    max R-hat: 2.90   min ESS (bulk): 57.2

Per-parameter posterior (mean +- std, central 95 %, truth in brackets;
``Tcold1``/``Twarm1`` are the two derived temperatures)::

    iso.instrument.calibration_scale.scale  +0.978 +- 0.038   95%[+0.904, +1.062]  (truth +1)
    iso.likelihood.amplitude                +36.7  +- 60.9    95%[+0.19,  +170.5]  (no truth)
    iso.likelihood.length_scale             +0.167 +- 0.070   95%[+0.016, +0.289]  (no truth)
    Tcold0                                  +32.0  +- 3.4     95%[+23.0,  +41.1]   (truth +36.13)
    Tcold_fraction                          +0.546 +- 0.158   95%[+0.130, +0.674]  (truth +0.4764)
    Tcold1 (derived)                        +58.2  +- 7.9     95%[+39.0,  +64.0]   (truth +57.03)
    Twarm0                                  +136.8 +- 24.0    95%[+81.5,  +167.3]  (truth +105.20)
    Twarm_fraction                          +0.057 +- 0.052   95%[+0.0002,+0.252]  (truth +0.2340)
    Twarm1 (derived)                        +139.7 +- 21.3    95%[+89.9,  +167.9]  (truth +122.70)
    logacold0..7, logawarm1/2/5/7                                                  (all within 95%)

**All eighteen qualified parameters plus both derived temperatures cover
truth at 95 %** -- a real improvement on the first draft's four misses,
confirming the exact ordered prior does what it was meant to (the "unordered
but data-compatible" failure mode ruling 2 fixed cannot occur any more).
But **R-hat (1.69-2.90) and ESS (bulk, 57-78) are barely changed from the
900-step run below** despite an eleven-fold larger budget: the walkers still
have not mixed. Neither diagnostic is within a factor of two of ruling 3's
1.05/400 thresholds (2.90 is 2.76x 1.05; 57.2 is 7.0x short of 400), so the
one permitted extension to 20 000/10 000 was not run -- it would not be
expected to close a gap this large, and ruling 3 asks for it only when the
run is close. This is reported as a finding, not reseeded.

**Diagnosis.** With the temperature-ordering degeneracy eliminated by
construction, the coverage run isolates the *other* cause the original
diagnosis flagged: the eleven independent log-abundance parameters give the
eight dust species room to trade off against each other (and against the
two shells' temperatures -- ``Twarm0``'s posterior mean, +136.8, sits nowhere
near its truth, +105.2, though its wide 95 % interval still covers it) in
more than one way that fits the data comparably well, and emcee's default
stretch move lets an ensemble split across those modes and never recombine.
That this persists essentially unchanged at 10 000 steps (R-hat/ESS are the
same order of magnitude as the 900-step run's -- and an ESS near the walker
count, fifty, is what arviz reports when each walker occupies its own
region) is consistent with two readings the run cannot tell apart:
structurally distinct modes, or a stretch-move autocorrelation time in
eighteen dimensions longer than the run itself. Either way, more steps at
the same move are the wrong lever; a different move or sampler is. A
different move set (the legacy script's own
``DEMove``/``DESnookerMove`` comment, ``examples/NGC6302.py`` lines
387-394), ``zeus`` (this twin's other engine), or a reduced/marginalised
abundance parameterisation are the natural next things to try, and are out
of this item's scope.

History: the first draft's four independent, unordered temperature boxes
were run at 50 walkers x 900 steps (450 burn-in), 2026-09-29, wall clock
2822.1 s. R-hat reached 3.12 and ESS (bulk) as low as 57 on every one of
the eighteen free parameters -- the ensemble under-run about fifty-fold
(:math:`\\tau \\approx 400` steps against 900 sampled) -- and four of the
sixteen qualified parameters' 95 % intervals missed the truth, all on the
disordered ``Tcold0``/``Twarm0`` pair or the abundances that trade off
against them. That under-run, and the coverage run above showing the same
non-convergence persisting after the ordering was fixed, is the lesson: this
model's multimodality has (at least) two independent sources, of which only
one was in this item's scope to fix.
"""

from __future__ import annotations

import argparse
import math
import sys
import time
from collections.abc import Mapping
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st
from scipy.interpolate import interp1d

from ampere.core import (
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    QuasisepGP,
    Spectrum,
)
from ampere.backends.reference import CalibrationScale, Resample
from ampere.inference import EmceeEngine, ZeusEngine

from . import dust_mass, generators

__all__ = [
    "CALIBRATION_PRIOR",
    "DEFAULT_BURN_IN",
    "DEFAULT_GRID",
    "DEFAULT_STEPS",
    "DEFAULT_WALKERS",
    "DERIVED_TRUTH",
    "ENGINES",
    "GP_AMPLITUDE_PRIOR",
    "GP_LENGTH_SCALE_PRIOR",
    "PHYSICAL_TRUTH",
    "QUALIFIED_TRUTH",
    "QUICK_BURN_IN",
    "QUICK_STEPS",
    "KemperTwoShell",
    "build_instrument",
    "build_model",
    "build_problem",
    "derived_temperatures",
    "fit",
    "main",
    "recovers_truth",
    "report",
]

# Legacy's own wavelength grid (``examples/NGC6302.py`` lines 295-297): a
# dense synthetic grid, not the observed data's own sampling -- the model's
# fallback grid before an instrument negotiates the real one.
_WAVE1 = np.linspace(2.3603, 35.0603, 327)
_WAVE2 = np.linspace(1.0 / 196.6261, 1.0 / 35.1, 117)
DEFAULT_GRID = np.concatenate((_WAVE1, 1.0 / _WAVE2[::-1]))

#: ``calUnc``'s v2 counterpart -- see the module docstring.
CALIBRATION_PRIOR = st.lognorm(0.05, scale=1.0)

#: legacy's own ``scalelengthPrior=0.1`` (micron, half-normal).
GP_LENGTH_SCALE_PRIOR = st.halfnorm(scale=0.1)
#: No legacy number to match; weakly informative, scaled to this spectrum's
#: own flux magnitude (see the module docstring).
GP_AMPLITUDE_PRIOR = st.halfnorm(scale=100.0)

ENGINES = ("emcee", "zeus")

#: The two solvers ``_likelihood``/``build_problem`` accept -- ruling 1 makes
#: ``quasisep`` (:class:`~ampere.core.QuasisepGP`, O(N)) the default in place
#: of ``dense`` (:class:`~ampere.core.DenseGP`, O(N^3)); both are kept so the
#: agreement row in :mod:`tests.examples.test_ngc6302` can build the problem
#: both ways.
_SOLVER_CLASSES: dict[str, type] = {"quasisep": QuasisepGP, "dense": DenseGP}

#: The upper edge of each temperature pair's shared box (legacy's own
#: ``lims``, see the module docstring) -- the ``80``/``180`` in
#: :func:`derived_temperatures`.
_TCOLD_UPPER = 80.0
_TWARM_UPPER = 180.0


def derived_temperatures(values: Mapping[str, Any]) -> dict[str, Any]:
    """``Tcold1``/``Twarm1`` from ``Tcold0``/``Tcold_fraction`` (and warm).

    See the module docstring's "Temperature ordering" section for why this
    recovers legacy's own physical temperatures exactly, as derived
    quantities rather than fitted parameters.

    Parameters
    ----------
    values
        Maps (at least) ``Tcold0``, ``Tcold_fraction``, ``Twarm0`` and
        ``Twarm_fraction`` to scalars or equal-shaped arrays (e.g. a fit's
        posterior draws, already flattened or not).

    Returns
    -------
    dict
        ``{"Tcold1": ..., "Twarm1": ...}``, in *values*' own shape.
    """
    tcold0 = np.asarray(values["Tcold0"], dtype=float)
    tcold_fraction = np.asarray(values["Tcold_fraction"], dtype=float)
    twarm0 = np.asarray(values["Twarm0"], dtype=float)
    twarm_fraction = np.asarray(values["Twarm_fraction"], dtype=float)
    return {
        "Tcold1": tcold0 + tcold_fraction * (_TCOLD_UPPER - tcold0),
        "Twarm1": twarm0 + twarm_fraction * (_TWARM_UPPER - twarm0),
    }


#: The fifteen model parameters plus the calibration factor -- the sixteen
#: *declared* parameters :func:`recovers_truth` checks.
QUALIFIED_TRUTH: dict[str, float] = {
    f"model.{name}": value for name, value in generators.TRUTH.items()
}
QUALIFIED_TRUTH["iso.instrument.calibration_scale.scale"] = generators.CALIBRATION_TRUTH

#: The two *derived* physical temperatures' truth (ruling 2) -- scored by
#: :func:`recovers_truth` and printed by :func:`report` alongside the sixteen
#: qualified (declared) parameters, but not part of :data:`QUALIFIED_TRUTH`
#: (they are not a posterior variable name).
DERIVED_TRUTH: dict[str, float] = {
    name: float(value) for name, value in derived_temperatures(generators.TRUTH).items()
}

#: All four physical temperatures' truth, in the "cool -> hot" reading order
#: -- :func:`report`'s combined block.
PHYSICAL_TRUTH: dict[str, float] = {
    "Tcold0": generators.TRUTH["Tcold0"],
    "Tcold1": DERIVED_TRUTH["Tcold1"],
    "Twarm0": generators.TRUTH["Twarm0"],
    "Twarm1": DERIVED_TRUTH["Twarm1"],
}

# Budgets. Legacy's own emcee script (50 walkers, 50 000 steps, 40 000
# burn-in -- ``examples/NGC6302.py`` lines 399, 456) is kept as the shared
# full-budget default for both engines (see the module docstring).
DEFAULT_WALKERS = 50
DEFAULT_STEPS = 50_000
DEFAULT_BURN_IN = 40_000
QUICK_STEPS = 2_000
QUICK_BURN_IN = 1_000

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
    Tcold0, Tcold_fraction, Twarm0, Twarm_fraction
        The cold shell's outer (cooler) temperature and the fraction of the
        remaining ``80 - Tcold0`` K the inner (hotter) one sits at (and
        likewise, over ``180 - Twarm0``, for the warm shell) -- kelvin and
        dimensionless respectively. :meth:`evaluate` derives the physical
        ``Tcold1``/``Twarm1`` (``ckmodbb``'s ``tin``) from these; see the
        module docstring's "Temperature ordering" section for why this
        parameterisation makes ``Tcold0 < Tcold1`` (and ``Twarm0 < Twarm1``)
        exact by construction rather than a rejected or unenforced box.
        :func:`derived_temperatures` recovers the physical pair from a
        point or a set of posterior draws.
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
        Tcold_fraction: Any,
        Twarm0: Any,
        Twarm_fraction: Any,
        channel: str = "sed",
    ) -> None:
        grid = np.asarray(wavelength, dtype=float)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(f"KemperTwoShell needs a 1-D, non-empty grid, got shape {grid.shape}.")
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
            wavelength_buffer = self.register_buffer(
                f"{name}_wavelength", table_wavelength, unit=u.um
            )
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
        self.register_parameter(
            _as_parameter("Tcold_fraction", Tcold_fraction, unit=u.dimensionless_unscaled)
        )
        self.register_parameter(_as_parameter("Twarm0", Twarm0, unit=u.K))
        self.register_parameter(
            _as_parameter("Twarm_fraction", Twarm_fraction, unit=u.dimensionless_unscaled)
        )

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

        # The exact ordered reparameterisation (module docstring, "Temperature
        # ordering", ruling 2): Tcold1/Twarm1 are derived here, not read from
        # ctx, and passed to _ckmodbb exactly as legacy calls it (Tcold1/
        # Twarm1 as tin, Tcold0/Twarm0 as tout).
        tcold0 = float(ctx["Tcold0"])
        tcold1 = tcold0 + float(ctx["Tcold_fraction"]) * (_TCOLD_UPPER - tcold0)
        twarm0 = float(ctx["Twarm0"])
        twarm1 = twarm0 + float(ctx["Twarm_fraction"]) * (_TWARM_UPPER - twarm0)

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
                tin=tcold1,
                tout=tcold0,
                n0=n0,
                grid=grid,
            )
        warm = np.zeros(grid.size)
        for index, n0 in enumerate(awarm):
            warm = warm + _ckmodbb(
                opacity[:, index],
                tin=twarm1,
                tout=twarm0,
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
    """The one :class:`KemperTwoShell`: legacy's own box priors on the eleven
    log-abundances, and the exact ordered-triangle prior (ruling 2, see the
    module docstring's "Temperature ordering") on the two temperature pairs.
    """
    abundance_prior = st.uniform(-6.0, 6.0)
    fraction_prior = st.uniform(0.0, 1.0)
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
        Tcold0=st.triang(c=0, loc=10.0, scale=70.0),
        Tcold_fraction=fraction_prior,
        Twarm0=st.triang(c=0, loc=80.0, scale=100.0),
        Twarm_fraction=fraction_prior,
    )


def build_instrument(observed_wavelength: Any) -> Instrument:
    """The ISO spectrum's instrument: resample onto the observed grid, then calibrate."""
    return Instrument(
        [Resample(observed_wavelength), CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed",
        label="iso",
    )


def _likelihood(*, gp: bool, solver: str = "quasisep") -> Likelihood:
    if not gp:
        return Likelihood(GaussianFamily(), IndependentNoise())
    kernel = Matern32(
        GP_AMPLITUDE_PRIOR,
        GP_LENGTH_SCALE_PRIOR,
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
        axes=("spectral_axis",),
    )
    # Ruling 1: QuasisepGP (O(N)) in place of DenseGP (O(N^3)) -- exact for a
    # Matern-3/2 kernel on sorted, strictly increasing 1-D coordinates, which
    # generators.load_observed_spectrum now guarantees. ``solver="dense"`` is
    # kept for the agreement row (tests/examples/test_ngc6302.py) that checks
    # the two give the same log_prob.
    solver_cls = _SOLVER_CLASSES[solver]
    return Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, solver_cls()))


def build_problem(
    *,
    synthetic: bool = False,
    gp: bool = True,
    seed: int = generators.SEED,
    solver: str = "quasisep",
) -> FittingProblem:
    """The composed problem: one model, one dataset (real, or ``--synthetic``)."""
    model = build_model()
    observed = generators.load_observed_spectrum()
    instrument = build_instrument(observed.spectral_axis.values)
    if synthetic:
        observed = generators.synthetic_data(model, instrument, seed=seed)
    datasets = DatasetCollection(
        {"iso": Dataset(observed, instrument, likelihood=_likelihood(gp=gp, solver=solver))}
    )
    return FittingProblem(model, datasets, seed=seed)


def fit(
    problem: FittingProblem,
    *,
    engine: str = "emcee",
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the named engine, at its own legacy-matched default."""
    if engine == "emcee":
        run_engine = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
    elif engine == "zeus":
        run_engine = ZeusEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
    else:
        raise SystemExit(f"unknown engine {engine!r}; choose {', '.join(ENGINES)}.")
    return run_engine.run(
        DEFAULT_STEPS if steps is None else steps,
        burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
        progress=progress,
    )


#: The four posterior variable names :func:`_derived_draws` needs.
_TEMPERATURE_FAMILY_NAMES = (
    "model.Tcold0",
    "model.Tcold_fraction",
    "model.Twarm0",
    "model.Twarm_fraction",
)


def _derived_draws(posterior: Any) -> dict[str, Any] | None:
    """:func:`derived_temperatures` on *posterior*'s draws, or ``None`` if the
    four declared temperature-family variables are not all present (e.g. a
    run built with a fixed, non-fitted temperature)."""
    if not all(name in posterior.data_vars for name in _TEMPERATURE_FAMILY_NAMES):
        return None
    return derived_temperatures(
        {
            "Tcold0": np.asarray(posterior["model.Tcold0"], dtype=float).ravel(),
            "Tcold_fraction": np.asarray(posterior["model.Tcold_fraction"], dtype=float).ravel(),
            "Twarm0": np.asarray(posterior["model.Twarm0"], dtype=float).ravel(),
            "Twarm_fraction": np.asarray(posterior["model.Twarm_fraction"], dtype=float).ravel(),
        }
    )


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each qualified parameter's, and each derived temperature's,
    central *level* interval covers its truth (ruling 2).
    """
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    derived = _derived_draws(posterior)
    if derived is not None:
        for name, truth in DERIVED_TRUTH.items():
            lower, upper = np.percentile(derived[name], [tail, 100.0 - tail])
            covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any) -> str:
    """A human-readable posterior summary, truth in brackets, 95 % coverage
    flagged; the four physical temperatures (``Tcold0``, the derived
    ``Tcold1``, ``Twarm0``, the derived ``Twarm1``) are printed together in
    one block after the declared parameters (ruling 2).
    """
    attrs = run.attrs
    posterior = run["posterior"].dataset
    has_truth = any(name in posterior.data_vars for name in QUALIFIED_TRUTH)
    covered = recovers_truth(run) if has_truth else {}
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        "  posterior (truth in brackets, 95 % coverage flagged, --synthetic only):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [2.5, 97.5])
        truth = QUALIFIED_TRUTH.get(name)
        flag = "" if truth is None else ("  ok" if covered.get(name) else "  MISS")
        bracket = "" if truth is None else f"  (truth {truth:+.6g})"
        lines.append(
            f"    {name:45s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )

    derived = _derived_draws(posterior)
    if derived is not None:
        tcold0 = np.asarray(posterior["model.Tcold0"], dtype=float).ravel()
        twarm0 = np.asarray(posterior["model.Twarm0"], dtype=float).ravel()
        physical = {
            "Tcold0": (tcold0, "model.Tcold0"),
            "Tcold1": (derived["Tcold1"], "Tcold1"),
            "Twarm0": (twarm0, "model.Twarm0"),
            "Twarm1": (derived["Twarm1"], "Twarm1"),
        }
        lines.append(
            "  the four physical temperatures (Tcold1/Twarm1 derived from "
            'Tcold_fraction/Twarm_fraction -- see "Temperature ordering"):'
        )
        for name, (draws, flag_key) in physical.items():
            lower, upper = np.percentile(draws, [2.5, 97.5])
            truth = PHYSICAL_TRUTH.get(name) if has_truth else None
            flag = "" if truth is None else ("  ok" if covered.get(flag_key) else "  MISS")
            bracket = "" if truth is None else f"  (truth {truth:+.6g})"
            lines.append(
                f"    {name:45s} {draws.mean():+.6g} +- {draws.std():.3g}   "
                f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
            )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="emcee", choices=list(ENGINES))
    parser.add_argument("--synthetic", action="store_true", help="fit the 2002-solution truth")
    parser.add_argument("--no-gp", dest="gp", action="store_false", help="IndependentNoise instead")
    parser.add_argument(
        "--quick", action="store_true", help=f"{QUICK_STEPS}/{QUICK_BURN_IN} budget"
    )
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--burn-in", type=int, default=None)
    parser.add_argument("--dust-mass", action="store_true", help="print the dust-mass table too")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem(synthetic=args.synthetic, gp=args.gp, seed=args.seed)
    print(f"free parameters: {problem.parameters.free_names}")

    steps = args.steps
    burn_in = args.burn_in
    if args.quick:
        steps = QUICK_STEPS if steps is None else steps
        burn_in = QUICK_BURN_IN if burn_in is None else burn_in

    started = time.perf_counter()
    run = fit(problem, engine=args.engine, walkers=args.walkers, steps=steps, burn_in=burn_in)
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")

    if args.dust_mass:
        table = dust_mass.dust_masses(run)
        print(dust_mass.format_table(table))
    return 0
