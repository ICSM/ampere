"""The v2 twin of ``examples/star_disc.py`` -- W6.13 (6).

The legacy script (unchanged) fits HD 105, a Sun-like star with a cold debris
belt, using ``QuickSEDModel`` (``ampere/legacy/models/QuickSED.py``): a
Starfish-emulated PHOENIX photosphere plus a modified blackbody. The data are
the tracked ``HD105_SED.vot`` photometry, a Spitzer IRS spectrum and a
synthetic Gaia-RVS spectrum, fitted with emcee. ``QuickSED`` is the one legacy
model class with no v2 counterpart yet; :class:`StarDisc` is that
counterpart. The star comes from ampere's own PHOENIX emulator
(:mod:`examples.phoenix_star.emulator`), not Starfish.

The model
---------
Parameters, with the legacy labels and ``lims`` as flat priors:

==================  ============  =====================================
twin                legacy        prior
==================  ============  =====================================
``luminosity``      ``lstar``     U(0.7, 1.3) L_sun, **linear** (legacy)
``teff``            ``tstar``     U(5500, 6500) K
``logg``            ``log_g``     U(4.1, 4.9)
``feh``             ``[fe/h]``    U(0.01, 0.1)
``log_area``        ``Adust``     U(-3, 3), log10 au
``t_dust``          ``tdust``     U(30, 100) K
``lambda_0``        ``lam0``      U(50, 500) micron
``beta``            ``beta``      U(0, 4)
==================  ============  =====================================

``distance`` is a buffer at ``1000 / 25.7534`` pc (the Gaia parallax the
legacy fixes). The emulator's [Fe/H] axis has only two nodes, 0.0 and +0.5,
so ``feh`` in 0.01-0.1 is GP interpolation between two planes. That is an
honest limit, and the legacy shared it through Starfish.

The star is :mod:`examples.phoenix_star`'s arithmetic: the emulator's
unit-bolometric shape scaled to ``F_nu`` in Jy for ``luminosity`` at
``distance``, with the correct inverse-square law (see "what changed"). It is
placed on the grid by log-log interpolation and extrapolated beyond 5.5
micron as ``F_nu ~ lambda**-2``, which is the legacy's ``F_lambda ~
lambda**-4``. There is no extinction, because ``QuickSED`` has none. The dust
is transcribed exactly from ``QuickSED.__call__`` (lines 158-164):
``B_nu(T_dust)`` in Jy sr^-1 (astropy's ``BlackBody().evaluate`` with bare
numbers returns cgs, and the legacy's ``1e23`` turns it into Jy), times
``(lambda_0 / lambda)**beta`` for ``lambda >= lambda_0``, times ``pi
((10**log_area au) / (d pc))**2``. The model's ``"sed"`` channel is the
legacy grid ``np.logspace(log10(0.31), 4, 1000)`` micron. A second channel,
``"rvs"``, carries the same star plus dust on the RVS window at R = 11 000,
read from the emulator's own R = 11 000 segment.

The data
--------
* ``"photometry"``: the votable's nineteen points as one
  :class:`~ampere.core.PhotometricPoints`. Seventeen filters come from the
  bundled library through ``SyntheticPhotometry.from_library``. The other two,
  ``ALMA/ALMA.B6`` and ``ATCA/ATCA.9mm``, are not in the library and are
  top-hats (``detector="energy"``): ALMA band 6 as 211-275 GHz (1.090-1.421
  mm), ATCA 9 mm as 30-38 GHz (7.889-9.993 mm). The legacy comment
  approximates the ATCA point this way, and the twin treats both the same.
* ``"rvs"``: the synthetic Gaia-RVS spectrum of ``star_disc.py`` lines 96-125,
  drawn from the twin's own emulator at the Marshall et al. (2018) parameters
  (1.216 L_sun, 6034 K, log g 4.478, [Fe/H] 0.02) with no dust. It spans
  0.847-0.871 micron at R = 11 000 with 5 % noise, ``calUnc=0.10`` as
  ``CalibrationScale(st.lognorm(0.10, scale=1.0))``, and ``scaleLengthPrior=
  0.01`` as the ``Matern32`` length-scale prior ``st.halfnorm(scale=0.01)``
  micron. The amplitude prior is weakly informative, ``st.halfnorm(scale=0.1)``
  Jy, a few per cent of the RVS flux. ``--no-gp`` uses ``IndependentNoise``.
* **IRS.** The legacy reads ``examples/star_disc/cassis_yaaar_spcfw_5295616t.fits``,
  which is **not in the repository**. The twin fits the photometry plus the
  RVS spectrum by default. ``--irs PATH`` adds the SL and LL chunks of a
  CASSIS file the user has fetched, read by :func:`.generators.read_irs`
  (the :mod:`examples.linear_sed` SL/LL split, deduplicated), each with
  ``Resample`` and a ``CalibrationScale`` and independent noise.
* ``HD105_SED.csv`` also lists two upper limits, SPIRE 500 micron and LABOCA
  870 micron. They are not in the votable the legacy reads, so the legacy
  never used them and neither does the twin. A ``Censoring`` declaration
  (``ampere/core/likelihood.py``) is where they would go, which is a follow-up.

``--synthetic`` replaces the votable's fluxes by the model at
:data:`.generators.SYNTHETIC_TRUTH` (the Marshall star with ``log_area =
0.5``, ``t_dust = 60``, ``lambda_0 = 150``, ``beta = 1.0``) plus Gaussian
noise at the votable's own uncertainties. That is the coverage run's truth.

What changed in the translation
-------------------------------
* **The inverse-square law.** The legacy divides the star by ``4 pi d**2``
  (``QuickSED.py`` line 148) although ``fbol_1l1p`` already carries the
  ``4 pi``, so its photosphere is ``4 pi`` too faint. The twin does not
  reproduce that. A luminosity fitted by the legacy is ``4 pi`` times the
  twin's.
* **The bolometric normalisation** integrates the emulated spectrum over
  500 A-5.5 micron. The legacy integrates over an extension to 1 cm, and the
  Rayleigh-Jeans tail beyond 5.5 micron adds about 0.1 % for these stars.
* **Vacuum wavelengths**, and ``Resample`` in place of ``spectres``, as in
  :mod:`examples.phoenix_star`.

::

    python -m examples.star_disc                   # the votable + RVS, emcee
    python -m examples.star_disc --synthetic       # the coverage run
    python -m examples.star_disc --irs PATH        # plus a fetched CASSIS file

Coverage run (Accept criterion)
-------------------------------
``pixi run -e dev python -m examples.star_disc --synthetic``: the votable's
nineteen points and the synthetic RVS spectrum, replaced by the model at the
Marshall star plus the synthetic dust (``generators.SYNTHETIC_TRUTH``) at
``generators.SEED`` (20260930); emcee on the reference backend at the legacy
budget, 40 walkers, 4000 steps with 1000 burn-in (3000 kept draws per
walker); 631 s wall clock, 2026-09-30. All eight truths fall inside the 95 %
intervals, **but the run has not converged**: R-hat is 1.10-1.42, above 1.05
on every parameter, and the minimum bulk ESS is 82, against the criterion's
R-hat < 1.05 and ESS > 400. The coverage verdict therefore does not count at
this budget. It was not re-run; the budget that converges is an open question
(W6.13 (6) report).

==========================  =========  ==============================  =========  =====
parameter                   mean       95 % interval                   truth      R-hat
==========================  =========  ==============================  =========  =====
``beta``                    1.0025     [0.9521, 1.1034]                1.0 ok     1.194
``feh``                     0.0178     [0.0103, 0.0593]                0.02 ok    1.191
``lambda_0``                150.02     [144.87, 161.79]                150 ok     1.187
``log_area``                0.49800    [0.49412, 0.50219]              0.5 ok     1.107
``logg``                    4.4807     [4.4710, 4.4981]                4.478 ok   1.418
``luminosity``              1.21570    [1.21447, 1.21781]              1.216 ok   1.193
``t_dust``                  60.12      [59.78, 60.43]                  60 ok      1.115
``teff``                    6033.9     [5959.1, 6085.7]                6034 ok    1.378
calibration scale           0.9959     [0.9765, 1.0036]                           1.190
GP amplitude                0.083      [0.0044, 0.219]                            1.096
GP length scale             8.4e-3     [3.0e-4, 2.3e-2]                           1.141
==========================  =========  ==============================  =========  =====
"""

from __future__ import annotations

import argparse
import sys
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st
from astropy import constants

from ampere.backends.reference import CalibrationScale, Resample, SyntheticPhotometry
from ampere.core import (
    DTYPE,
    ChannelRequirements,
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
    Spectrum,
)
from examples.phoenix_star.emulator import (
    EMULATOR_FILE,
    NUMPY_OPS,
    ChannelPlan,
    as_parameter,
    channel_plan,
    emulator_arrays,
    load_emulator,
    log_shape,
    plan_arrays,
    star_log_flux,
)

from . import generators

__all__ = [
    "GRID",
    "PARALLAX",
    "PRIORS",
    "QUALIFIED_TRUTH",
    "TOPHATS",
    "StarDisc",
    "build_instruments",
    "build_model",
    "build_problem",
    "dust_flux",
    "fit",
    "main",
    "photometry_step",
    "planck_jy",
    "recovers_truth",
    "report",
]

#: The legacy model grid (``star_disc.py`` line 139), micron.
GRID = np.logspace(np.log10(0.31), np.log10(10000), 1000, endpoint=True)

#: The Gaia DR2 parallax the legacy fixes, mas.
PARALLAX = 25.7534

PRIORS: dict[str, Any] = {
    "luminosity": st.uniform(0.7, 0.6),
    "teff": st.uniform(5500.0, 1000.0),
    "logg": st.uniform(4.1, 0.8),
    "feh": st.uniform(0.01, 0.09),
    "log_area": st.uniform(-3.0, 6.0),
    "t_dust": st.uniform(30.0, 70.0),
    "lambda_0": st.uniform(50.0, 450.0),
    "beta": st.uniform(0.0, 4.0),
}

#: The two votable points with no library filter, as top-hats: name ->
#: (blue edge, red edge) in micron, from 211-275 GHz and 30-38 GHz.
TOPHATS: dict[str, tuple[float, float]] = {
    "ALMA/ALMA.B6": (
        float((275 * u.GHz).to(u.um, u.spectral()).value),
        float((211 * u.GHz).to(u.um, u.spectral()).value),
    ),
    "ATCA/ATCA.9mm": (
        float((38 * u.GHz).to(u.um, u.spectral()).value),
        float((30 * u.GHz).to(u.um, u.spectral()).value),
    ),
}

CALIBRATION_PRIOR = st.lognorm(0.10, scale=1.0)
IRS_CALIBRATION_PRIOR = st.lognorm(0.10, scale=1.0)
GP_LENGTH_SCALE_PRIOR = st.halfnorm(scale=0.01)
GP_AMPLITUDE_PRIOR = st.halfnorm(scale=0.1)

QUALIFIED_TRUTH: dict[str, float] = {f"model.{k}": v for k, v in generators.SYNTHETIC_TRUTH.items()}

# The legacy budget (40 walkers, 4000 samples, 1000 burn-in) and --quick's.
DEFAULT_WALKERS = 40
DEFAULT_STEPS = 4000
DEFAULT_BURN_IN = 1000
QUICK_STEPS = 800
QUICK_BURN_IN = 200

_H = constants.h.cgs.value
_C = constants.c.cgs.value
_K = constants.k_B.cgs.value
#: The legacy's own constants (``QuickSED.py`` lines 57-58), metres; only
#: their ratio enters.
_AU = 1.495978707e11
_PC = 3.0857e16


def planck_jy(wavelength: Any, temperature: Any) -> np.ndarray:
    """``B_nu(T)`` in Jy sr^-1 at *wavelength* (micron): cgs Planck x 1e23."""
    nu = _C / (np.asarray(wavelength, dtype=float) * 1e-4)
    # Deep in the Wien tail expm1 overflows to inf and the flux is exactly 0.
    with np.errstate(over="ignore"):
        return 1e23 * 2.0 * _H * nu**3 / _C**2 / np.expm1(_H * nu / (_K * temperature))


def dust_flux(
    wavelength: Any, log_area: float, t_dust: float, lambda_0: float, beta: float, distance: float
) -> np.ndarray:
    """``QuickSED``'s dust term, Jy (lines 158-164, transcribed)."""
    wavelength = np.asarray(wavelength, dtype=float)
    emission = planck_jy(wavelength, t_dust)
    modified = wavelength >= lambda_0
    emission = np.where(modified, emission * (lambda_0 / wavelength) ** beta, emission)
    return emission * np.pi * ((10**log_area * _AU) / (distance * _PC)) ** 2


class StarDisc(Model):
    """``QuickSEDModel``'s twin: an emulated PHOENIX star plus a modified blackbody.

    Parameters
    ----------
    wavelength
        The ``"sed"`` channel's grid, micron (default :data:`GRID`).
    rvs_wavelength
        The ``"rvs"`` channel's grid (default :func:`.generators.rvs_grid`).
    luminosity, teff, logg, feh, log_area, t_dust, lambda_0, beta
        A prior to fit it, or a number to hold it fixed.
    distance
        Buffer, parsec.
    """

    CHANNELS: ClassVar[tuple[str, ...]] = ("sed", "rvs")
    AXIS: ClassVar[str] = "spectral_axis"

    def __init__(
        self,
        wavelength: Any = GRID,
        *,
        rvs_wavelength: Any = None,
        distance: float = 1000.0 / PARALLAX,
        path: str | Path = EMULATOR_FILE,
        **priors: Any,
    ) -> None:
        unknown = set(priors) - set(PRIORS)
        if unknown:
            raise TypeError(f"StarDisc has no parameter(s) {sorted(unknown)}")
        self.data = load_emulator(path)
        for name, default in PRIORS.items():
            self.register_parameter(as_parameter(name, priors.get(name, default)))
        self.register_buffer("distance", np.asarray(distance, dtype=DTYPE), unit=u.pc)
        self.arrays = emulator_arrays(NUMPY_OPS, self.data)
        self.plans: dict[str, ChannelPlan] = {}
        self.plan_arrays: dict[str, dict[str, Any]] = {}
        self.templates: dict[str, Spectrum] = {}
        rvs = generators.rvs_grid() if rvs_wavelength is None else rvs_wavelength
        self._adopt("sed", np.asarray(wavelength, dtype=DTYPE))
        self._adopt("rvs", np.asarray(rvs, dtype=DTYPE))

    def _adopt(self, channel: str, grid: np.ndarray) -> None:
        plan = channel_plan(self.data, channel, grid, extinction=False)
        self.plans[channel] = plan
        self.plan_arrays[channel] = plan_arrays(NUMPY_OPS, plan)
        self.templates[channel] = Spectrum(grid * u.um, np.zeros(grid.size, dtype=DTYPE), unit=u.Jy)

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each channel, once."""
        for channel in self.CHANNELS:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            self._adopt(channel, np.asarray(asked[self.AXIS].coordinates(), dtype=DTYPE))
        return self

    def star(self, channel: str, ctx: Mapping[str, Any]) -> np.ndarray:
        """The photosphere on *channel*, Jy."""
        shape = log_shape(NUMPY_OPS, self.arrays, ctx["teff"], ctx["logg"], ctx["feh"])
        log_flux = star_log_flux(
            NUMPY_OPS,
            self.plans[channel],
            self.plan_arrays[channel],
            shape,
            np.log10(float(ctx["luminosity"])),
            ctx["distance"],
        )
        return 10.0**log_flux

    def dust(self, channel: str, ctx: Mapping[str, Any]) -> np.ndarray:
        """The dust on *channel*, Jy."""
        return dust_flux(
            self.plans[channel].wavelength,
            float(ctx["log_area"]),
            float(ctx["t_dust"]),
            float(ctx["lambda_0"]),
            float(ctx["beta"]),
            float(ctx["distance"]),
        )

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        emitted = {}
        for channel in self.CHANNELS:
            total = self.star(channel, ctx) + self.dust(channel, ctx)
            emitted[channel] = self.templates[channel].with_values(np.asarray(total, dtype=DTYPE))
        return ModelResult(emitted)


def build_model() -> StarDisc:
    """The star-plus-disc model with the legacy ``lims`` as flat priors."""
    return StarDisc()


def photometry_step(names: list[str]) -> SyntheticPhotometry:
    """One step for all the votable's filters: the library's plus the top-hats."""
    library = [name for name in names if name not in TOPHATS]
    step = SyntheticPhotometry.from_library(library, GRID)
    response = np.asarray(step.buffers["response"].array, dtype=float)
    rows = {name: response[i] for i, name in enumerate(library)}
    pivots = {name: float(p) for name, p in zip(library, step.pivots(), strict=True)}
    detectors = dict(zip(library, step.detectors, strict=True))
    for name, (blue, red) in TOPHATS.items():
        if name not in names:
            continue
        curve = ((GRID >= blue) & (GRID <= red)).astype(float)
        rows[name] = curve
        pivots[name] = float(
            np.sqrt(np.trapezoid(curve * GRID, GRID) / np.trapezoid(curve / GRID, GRID))
        )
        detectors[name] = "energy"
    return SyntheticPhotometry(
        names,
        GRID,
        np.array([rows[name] for name in names]),
        detector=[detectors[name] for name in names],
        pivots=np.array([pivots[name] for name in names]),
    )


def _gp_likelihood(*, gp: bool) -> Likelihood:
    if not gp:
        return Likelihood(GaussianFamily(), IndependentNoise())
    kernel = Matern32(
        GP_AMPLITUDE_PRIOR,
        GP_LENGTH_SCALE_PRIOR,
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
        axes=("spectral_axis",),
    )
    return Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))


def build_instruments(
    names: list[str], irs: list[tuple[np.ndarray, np.ndarray, np.ndarray]] | None = None
) -> dict[str, Instrument]:
    """The photometry (``"sed"``), the RVS spectrograph (``"rvs"``) and any IRS chunks."""
    instruments = {
        "photometry": Instrument([photometry_step(names)], channel="sed", label="photometry"),
        "rvs": Instrument(
            [Resample(generators.rvs_grid()), CalibrationScale(CALIBRATION_PRIOR)],
            channel="rvs",
            label="rvs",
        ),
    }
    for label, (wavelength, _, _) in zip(("sl", "ll"), irs or [], strict=False):
        instruments[label] = Instrument(
            [Resample(wavelength), CalibrationScale(IRS_CALIBRATION_PRIOR)],
            channel="sed",
            label=label,
        )
    return instruments


def build_problem(
    *,
    synthetic: bool = False,
    gp: bool = True,
    irs: str | Path | None = None,
    seed: int = generators.SEED,
) -> FittingProblem:
    """The composed problem: the votable (or its synthetic twin), RVS, optional IRS."""
    rng = np.random.default_rng(seed)
    names, flux, error = generators.read_votable()
    chunks = generators.read_irs(irs) if irs is not None else None
    instruments = build_instruments(names, chunks)
    unity = {"calibration_scale.scale": 1.0}
    calibration = {label: unity for label in instruments if label != "photometry"}

    # The RVS spectrum: the Marshall star, no dust, 5 % noise (both modes).
    star_only = {**generators.MARSHALL, "log_area": -3.0, "t_dust": 30.0}
    star_only |= {"lambda_0": 50.0, "beta": 0.0}
    rvs_truth = generators.observe(
        build_model(), {"rvs": instruments["rvs"]}, star_only, {"rvs": unity}
    )["rvs"]
    rvs_sigma = generators.RVS_FRACTIONAL_NOISE * np.abs(rvs_truth.values)
    observed_rvs = generators.with_values(
        rvs_truth, rvs_truth.values + rng.normal(0.0, rvs_sigma), rvs_sigma
    )

    truth = generators.SYNTHETIC_TRUTH
    shapes = generators.observe(build_model(), instruments, truth, calibration)
    phot = shapes["photometry"]
    values = phot.values + rng.normal(0.0, error) if synthetic else flux
    observed = {
        "photometry": generators.with_values(phot, values, error),
        "rvs": observed_rvs,
    }
    for label, (_, irs_flux, irs_error) in zip(("sl", "ll"), chunks or [], strict=False):
        observed[label] = generators.with_values(shapes[label], irs_flux, irs_error)

    datasets = {
        label: Dataset(
            observed[label],
            instrument,
            likelihood=_gp_likelihood(gp=gp)
            if label == "rvs"
            else Likelihood(GaussianFamily(), IndependentNoise()),
        )
        for label, instrument in instruments.items()
    }
    return FittingProblem(build_model(), DatasetCollection(datasets), seed=seed)


def fit(
    problem: FittingProblem,
    *,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    progress: bool = False,
) -> Any:
    """emcee, the legacy's engine, at the legacy budget by default."""
    from ampere.inference import EmceeEngine

    engine = EmceeEngine(problem, walkers=walkers or DEFAULT_WALKERS)
    return engine.run(
        steps or DEFAULT_STEPS,
        burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
        progress=progress,
    )


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each of the eight parameters' central *level* interval covers the synthetic truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any, *, synthetic: bool = True) -> str:
    """The posterior summary (truth flagged under ``--synthetic``), R-hat and ESS."""
    from examples.phoenix_star.phoenix_star import report as _report

    text = _report(run, QUALIFIED_TRUTH if synthetic else {})
    try:
        import arviz as az

        ess = az.ess(run["posterior"].dataset)
        worst = min(float(ess[name]) for name in ess.data_vars)
        text += f"\n  minimum bulk ESS {worst:.0f}"
    except Exception:
        pass
    return text


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--synthetic", action="store_true", help="the coverage-run data")
    parser.add_argument("--irs", type=Path, default=None, help="a fetched CASSIS file")
    parser.add_argument("--no-gp", dest="gp", action="store_false", help="IndependentNoise")
    parser.add_argument("--quick", action="store_true", help="800 steps, 200 burn-in")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--burn-in", type=int, default=None)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)
    problem = build_problem(synthetic=args.synthetic, gp=args.gp, irs=args.irs, seed=args.seed)
    print(f"free parameters: {problem.parameters.free_names}")
    steps = args.steps or (QUICK_STEPS if args.quick else None)
    burn_in = args.burn_in if args.burn_in is not None else (QUICK_BURN_IN if args.quick else None)
    started = time.perf_counter()
    run = fit(problem, walkers=args.walkers, steps=steps, burn_in=burn_in)
    elapsed = time.perf_counter() - started
    print(report(run, synthetic=args.synthetic))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
