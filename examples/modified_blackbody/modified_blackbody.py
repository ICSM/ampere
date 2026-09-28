"""The v2 twin of ``examples/examples_paper/modifiedblackbody.py`` -- W6.13 (3).

The legacy file stays exactly as it is (the characterisation anchor this
module does **not** touch); this package is the same four-parameter model,
the same ten-band photometry and the same question, written once against
:mod:`ampere.core` and the reference backend, with every engine the legacy
script ran by hand now behind one ``--engine``/``--all`` flag.

The model
---------
:class:`ModifiedBlackBody` transcribes the legacy ``ModifiedBlackBody.__call__``
formula exactly (``examples/examples_paper/modifiedblackbody.py`` lines 39-43):
``astropy.modeling.physical_models.BlackBody`` at ``(temperature)``, scaled by
``1 / (distance in cm)**2``, by ``10**logmass`` grams, by the fixed
``kappa=10`` and by ``(lambda / kappa_wavelength=250 micron)**beta``. Two
things about that formula are kept faithfully rather than "fixed", because
this item transcribes the legacy model rather than correcting it:

* The legacy script's own comment calls ``d_true = 0.11`` "kpc", but the code
  converts it with ``u.pc.to(u.cm)`` -- a **parsec**-to-cm factor -- so the
  number is actually treated as 0.11 parsec, not 0.11 kiloparsec. The
  docstring flags the mismatch; the twin reproduces the code's behaviour
  (parsec), not the comment's claim.
* ``BlackBody().evaluate`` is called with a bare (unitless) frequency array
  rather than a ``Quantity``, so its ``scale=1 * u.Jy / u.sr`` never actually
  attaches a unit to the result -- the legacy script's own synthetic flux is
  therefore a dimensionless number of order ``1e-11`` at the model's own
  ``kappa``/``distance``/``logmass``, not a physical flux density in Jy. This
  twin keeps the same numeric recipe (so the same truth reproduces the same
  numbers) and simply *declares* the result's container unit as Jy, exactly
  as the legacy script's ``Photometry(..., photunits="Jy")`` did -- neither
  version's absolute scale is physically meaningful, only its recovery under
  noise is being demonstrated.

The data
--------
One instrument, ``"catalogue"`` --
:meth:`~ampere.backends.reference.SyntheticPhotometry.from_library` on the
ten AKARI/Herschel bands the legacy script names (line 81), 10 % noise.
Every one of the ten loads from ampere's bundled filter library (checked in
:mod:`tests.examples.test_modified_blackbody`); none needed dropping.

Engines
-------
``--engine emcee|zeus|dynesty|sbi`` runs one engine at the legacy script's
own budget for it (``examples/examples_paper/modifiedblackbody.py`` lines
131-145): emcee 100 walkers / 1000 steps / 900 burn-in, zeus 100 walkers /
150 steps / 100 burn-in, dynesty ``dlogz=1.0``, sbi 50 000 simulations / 10
000 posterior draws. sbi needs the ``sbi`` extra and is refused by name,
naming the install command, when it is absent. ``--all`` runs every engine
available in the current environment (sbi included only where the extra is
installed) and draws :func:`overlay_figure`: one composite figure, saved
under ``--figures DIR`` (nothing this module produces is committed), with
each engine's own :func:`~ampere.results.plot_posterior_predictive` panel as
a labelled subplot -- see that function's docstring for why an embedded
image rather than one shared axis.

Coverage
--------
``recovers_truth`` checks all four parameters, the only ones this model has.

::

    python -m examples.modified_blackbody --engine emcee
    python -m examples.modified_blackbody --all --figures /tmp/mbb-figures

Coverage run (Accept criterion)
--------------------------------
Each engine once at the legacy script's own budget, ``--seed 20260928``, run
2026-09-28 (the dev engines in one process through the machine's gate lock,
the sbi arm in ``-e sbi``). The central 95 % interval covers the truth on
all four parameters under every engine, on the first run (no reseeding
needed)::

    emcee   100 walkers / 1000 steps / 900 burn-in       wall clock 356.8 s
    model.temperature   +28.6731  +/- 1.5     95%[+26.0632,  +31.6419]  (truth +30)    ok
    model.logmass       +0.898954 +/- 0.292   95%[+0.378884, +1.35598]  (truth +1)     ok
    model.beta          -2.10879  +/- 0.163   95%[-2.4326,   -1.8001]   (truth -2)     ok
    model.distance      +0.0969118 +/- 0.0297 95%[+0.0517114, +0.147434] (truth +0.11) ok

    zeus    100 walkers / 150 steps / 100 burn-in        wall clock 424.6 s
    model.temperature   +29.1521  +/- 2.78    95%[+25.6915,  +38.6073]  (truth +30)    ok
    model.logmass       +0.875356 +/- 0.32    95%[+0.255603, +1.33331]  (truth +1)     ok
    model.beta          -2.06893  +/- 0.291   95%[-2.48119,  -1.22664]  (truth -2)     ok
    model.distance      +0.0982483 +/- 0.0298 95%[+0.0523269, +0.148705] (truth +0.11) ok

    dynesty dlogz=1.0 (1476 equal-weight draws)         wall clock 42.7 s
    model.temperature   +28.8348  +/- 1.57    95%[+26.0089,  +32.1987]  (truth +30)    ok
    model.logmass       +0.9435   +/- 0.27    95%[+0.412237, +1.34536]  (truth +1)     ok
    model.beta          -2.09291  +/- 0.171   95%[-2.42789,  -1.75449]  (truth -2)     ok
    model.distance      +0.101998 +/- 0.0284  95%[+0.0540334, +0.147008] (truth +0.11) ok

    sbi     50 000 simulations / 10 000 draws (197 epochs) wall clock 891.4 s
    model.temperature   +38.061   +/- 11.2    95%[+15.6214,  +49.9983]  (truth +30)    ok
    model.logmass       +0.457771 +/- 0.956   95%[-1.73619,  +2.3576]   (truth +1)     ok
    model.beta          -1.56658  +/- 0.821   95%[-2.91674,  -0.15714]  (truth -2)     ok
    model.distance      +0.0925207 +/- 0.0296 95%[+0.0507863, +0.1454]  (truth +0.11)  ok

The sbi posterior is several times wider than the three samplers' on
temperature, log-mass and beta: the model's flux spans many decades across
the ten bands, and the network's z-scoring warns of extreme outliers at
training time (see the smoke test's finding) -- a calibrated but loose
amortised posterior, not a miss.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st
from astropy.modeling.physical_models import BlackBody

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Instrument,
    Likelihood,
    Model,
    ModelResult,
    Parameter,
    Spectrum,
)
from ampere.backends.reference import SyntheticPhotometry
from ampere.inference import DynestyEngine, EmceeEngine, SBIEngine, ZeusEngine
from ampere.results import add_posterior_predictive, plot_posterior_predictive

from . import generators

__all__ = [
    "DEFAULT_DYNESTY_DLOGZ",
    "DEFAULT_EMCEE_BURN_IN",
    "DEFAULT_EMCEE_STEPS",
    "DEFAULT_EMCEE_WALKERS",
    "DEFAULT_SBI_BUDGET",
    "DEFAULT_SBI_DRAWS",
    "DEFAULT_ZEUS_BURN_IN",
    "DEFAULT_ZEUS_STEPS",
    "DEFAULT_ZEUS_WALKERS",
    "ENGINES",
    "FILTERS",
    "GRID",
    "KAPPA",
    "KAPPA_WAVELENGTH",
    "QUALIFIED_TRUTH",
    "ModifiedBlackBody",
    "available_engines",
    "build_instrument",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "overlay_figure",
    "recovers_truth",
    "report",
    "run_all",
]

#: The legacy model's own grid (``examples/examples_paper/modifiedblackbody.py``
#: line 63), micron.
GRID = 10 ** np.linspace(0.0, 2.9, 2000)

#: Fixed, not fitted -- the legacy script's own defaults (line 12).
KAPPA = 10.0
KAPPA_WAVELENGTH = 250.0

#: The legacy script's ten AKARI/Herschel bands (line 81).
FILTERS = (
    "AKARI_FIS_N160",
    "AKARI_FIS_N60",
    "AKARI_FIS_WIDEL",
    "AKARI_FIS_WIDES",
    "HERSCHEL_PACS_100",
    "HERSCHEL_PACS_160",
    "HERSCHEL_PACS_70",
    "HERSCHEL_SPIRE_250",
    "HERSCHEL_SPIRE_350",
    "HERSCHEL_SPIRE_500",
)

#: The response tabulation, micron -- within :data:`GRID`'s own coverage.
PHOTOMETRY_TABULATION = np.geomspace(1.0, 794.0, 500)

#: The legacy script's own flat priors (line 18-21): temperature, log10
#: mass, beta, distance.
PRIOR_LIMITS: dict[str, tuple[float, float]] = {
    "temperature": (10.0, 50.0),
    "logmass": (-3.0, 3.0),
    "beta": (-3.0, 0.0),
    "distance": (0.05, 0.15),
}

ENGINES = ("emcee", "zeus", "dynesty", "sbi")

QUALIFIED_TRUTH: dict[str, float] = {
    f"model.{name}": value for name, value in generators.TRUTH.items()
}

# Budgets, one per engine, all the legacy script's own (lines 131-145).
DEFAULT_EMCEE_WALKERS = 100
DEFAULT_EMCEE_STEPS = 1000
DEFAULT_EMCEE_BURN_IN = 900
DEFAULT_ZEUS_WALKERS = 100
DEFAULT_ZEUS_STEPS = 150
DEFAULT_ZEUS_BURN_IN = 100
DEFAULT_DYNESTY_DLOGZ = 1.0
DEFAULT_SBI_BUDGET = 50_000
DEFAULT_SBI_DRAWS = 10_000


def _as_parameter(name: str, spec: Any, *, unit: u.UnitBase | None = None) -> Parameter:
    """A prior (a frozen scipy distribution) or a bare number, as a Parameter."""
    if isinstance(spec, Parameter):
        return spec
    if hasattr(spec, "logpdf") or hasattr(spec, "rvs"):
        return Parameter(name, spec, unit=unit)
    return Parameter(name, None, value=float(spec), fixed=True, unit=unit)


class ModifiedBlackBody(Model):
    """The paper's four-parameter modified blackbody -- see the module docstring.

    ``ampere.backends.reference.ModifiedBlackBody`` is the shipped, three-
    parameter form (temperature, beta, a dimensionless ``scale``): it folds
    distance and mass into one multiplicative factor, which is the right
    parameterisation for a source whose distance is not independently known.
    This twin keeps the paper's own four parameters (distance and log-mass
    separate) because that is what the legacy script demonstrates and what
    ``recovers_truth`` is accepted against -- not because the shipped form is
    wrong for its own purpose.

    Parameters
    ----------
    wavelength
        The model's own grid, micron (a fallback -- negotiation replaces it
        the moment an instrument is bound).
    temperature, logmass, beta, distance
        A prior to fit each, or a number to hold it fixed. Distance is in
        the same numeric convention the legacy script used -- see the module
        docstring's note on ``u.pc.to(u.cm)``.
    kappa, kappa_wavelength
        Fixed (buffers, never fitted), as in the legacy script.
    channel
        Name of the channel the emitted :class:`~ampere.core.Spectrum`
        appears under.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        temperature: Any,
        logmass: Any,
        beta: Any,
        distance: Any,
        kappa: float = KAPPA,
        kappa_wavelength: float = KAPPA_WAVELENGTH,
        channel: str = "sed",
    ) -> None:
        grid = np.asarray(wavelength, dtype=float)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(
                f"ModifiedBlackBody needs a 1-D, non-empty grid, got shape {grid.shape}."
            )
        self.channel = str(channel)
        self.register_buffer("wavelength", grid, unit=u.um)
        self.register_buffer("kappa", float(kappa))
        self.register_buffer("kappa_wavelength", float(kappa_wavelength), unit=u.um)
        self.register_parameter(_as_parameter("temperature", temperature, unit=u.K))
        self.register_parameter(_as_parameter("logmass", logmass))
        self.register_parameter(_as_parameter("beta", beta))
        self.register_parameter(_as_parameter("distance", distance))
        self._template: Spectrum | None = None

    def compile_for(self, requirements: Any) -> ModifiedBlackBody:
        """Adopt the negotiated grid for :attr:`channel`, once (W2.12's contract)."""
        asked = requirements.get(self.channel)
        if asked is not None and "spectral_axis" in asked:
            grid = np.asarray(asked["spectral_axis"].coordinates(), dtype=float)
            self._template = Spectrum(grid * u.um, np.zeros(grid.size), unit=u.Jy)
        return self

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        if self._template is not None:
            grid = self._template.spectral_axis.values
        else:
            grid = np.asarray(ctx["wavelength"], dtype=float)

        # The legacy formula, transcribed exactly -- see the module docstring
        # for the two quirks (parsec vs. "kpc"; the unitless BlackBody call)
        # kept faithfully rather than corrected. The Wien tail overflows
        # `expm1` harmlessly at short wavelengths / cold temperatures (the
        # correctly-computed answer there is zero flux); silenced rather than
        # left to print on every evaluation.
        freq = (grid * u.um).to(u.Hz, equivalencies=u.spectral()).value
        with np.errstate(over="ignore"):
            flux = BlackBody().evaluate(freq, float(ctx["temperature"]), 1.0 * u.Jy / u.sr)
        flux = flux / (float(ctx["distance"]) * u.pc.to(u.cm)) ** 2
        flux = (
            flux
            * (10.0 ** float(ctx["logmass"]))
            * u.Msun.to(u.g)
            * float(ctx["kappa"])
            * (grid / float(ctx["kappa_wavelength"])) ** float(ctx["beta"])
        )
        flux = np.asarray(flux, dtype=float)

        spectrum = (
            self._template.with_values(flux)
            if self._template is not None
            else Spectrum(grid * u.um, flux, unit=u.Jy)
        )
        return ModelResult({self.channel: spectrum})


def build_model() -> ModifiedBlackBody:
    """The one :class:`ModifiedBlackBody`, with the legacy script's flat priors."""
    return ModifiedBlackBody(
        GRID,
        temperature=st.uniform(*_width(PRIOR_LIMITS["temperature"])),
        logmass=st.uniform(*_width(PRIOR_LIMITS["logmass"])),
        beta=st.uniform(*_width(PRIOR_LIMITS["beta"])),
        distance=st.uniform(*_width(PRIOR_LIMITS["distance"])),
    )


def _width(limits: tuple[float, float]) -> tuple[float, float]:
    """``(low, high)`` as ``scipy.stats.uniform``'s ``(loc, scale)``."""
    low, high = limits
    return low, high - low


def build_instrument() -> Instrument:
    """The photometric catalogue, ten AKARI/Herschel bands, one channel."""
    return Instrument(
        [SyntheticPhotometry.from_library(list(FILTERS), PHOTOMETRY_TABULATION)],
        channel="sed",
        label="catalogue",
    )


def build_problem(*, seed: int = generators.SEED) -> FittingProblem:
    """The composed problem: one model, one dataset."""
    model = build_model()
    instrument = build_instrument()
    observed_photometry = generators.synthetic_data(model, instrument, seed=seed)
    datasets = DatasetCollection(
        {
            "catalogue": Dataset(
                observed_photometry,
                instrument,
                likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
            ),
        }
    )
    return FittingProblem(model, datasets, seed=seed)


def available_engines() -> tuple[str, ...]:
    """Every engine name this environment can actually run."""
    if importlib.util.find_spec("sbi") is None:
        return tuple(name for name in ENGINES if name != "sbi")
    return ENGINES


def _needs_sbi() -> None:
    if importlib.util.find_spec("sbi") is None:
        raise SystemExit(
            "engine='sbi' needs the 'sbi' extra, which is not installed in this environment. "
            "Install it with `pixi install -e sbi --frozen` and run this module again with "
            "`pixi run -e sbi python -m examples.modified_blackbody --engine sbi`."
        )


def fit(
    problem: FittingProblem,
    *,
    engine: str = "emcee",
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    live_points: int | None = None,
    dlogz: float | None = None,
    budget: int | None = None,
    draws: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the named engine, at its own legacy-matched default."""
    if engine == "emcee":
        run_engine = EmceeEngine(
            problem, walkers=DEFAULT_EMCEE_WALKERS if walkers is None else walkers
        )
        return run_engine.run(
            DEFAULT_EMCEE_STEPS if steps is None else steps,
            burn_in=DEFAULT_EMCEE_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "zeus":
        run_engine = ZeusEngine(
            problem, walkers=DEFAULT_ZEUS_WALKERS if walkers is None else walkers
        )
        return run_engine.run(
            DEFAULT_ZEUS_STEPS if steps is None else steps,
            burn_in=DEFAULT_ZEUS_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "dynesty":
        run_engine = DynestyEngine(problem, live_points=live_points)
        return run_engine.run(
            progress=progress, dlogz=DEFAULT_DYNESTY_DLOGZ if dlogz is None else dlogz
        )
    if engine == "sbi":
        _needs_sbi()
        run_engine = SBIEngine(problem, budget=DEFAULT_SBI_BUDGET if budget is None else budget)
        return run_engine.run(DEFAULT_SBI_DRAWS if draws is None else draws, progress=progress)
    raise SystemExit(f"unknown engine {engine!r}; choose {', '.join(ENGINES)}.")


def run_all(problem: FittingProblem, **tiny: Any) -> dict[str, Any]:
    """Fit *problem* on every engine :func:`available_engines` names.

    ``**tiny`` is forwarded to every engine's own call to :func:`fit`
    (walkers/steps/burn_in/live_points/dlogz/budget/draws) -- pass nothing
    for the legacy-matched defaults.
    """
    return {name: fit(problem, engine=name, **tiny) for name in available_engines()}


def overlay_figure(runs: dict[str, Any], problem: FittingProblem, *, band: float = 0.68) -> Any:
    """One composite figure, one labelled subplot per engine's own predictive band.

    :func:`~ampere.results.plot_posterior_predictive` builds its own new
    :class:`~matplotlib.figure.Figure` (``ampere/results/plots.py``: it takes
    no ``axes=`` to draw several runs onto), so there is no single shared
    axis to overlay four coloured bands on. Instead, each engine's own
    figure is rendered once and embedded as an image in one subplot of a new
    composite figure -- still "one figure" with all four engines' predictive
    checks on it, titled by engine, at the cost of a legend's colour key.
    """
    import matplotlib.pyplot as plt

    names = list(runs)
    figure, axes = plt.subplots(1, len(names), figsize=(5.0 * len(names), 4.5), squeeze=False)
    for axis, name in zip(axes[0], names, strict=True):
        tree = add_posterior_predictive(runs[name], problem)
        sub_figure = plot_posterior_predictive(tree, band=band)
        sub_figure.canvas.draw()
        image = np.asarray(sub_figure.canvas.buffer_rgba())
        axis.imshow(image)
        axis.set_title(name)
        axis.axis("off")
        plt.close(sub_figure)
    figure.suptitle("posterior-predictive overlay, by engine")
    return figure


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each of the four parameters' central *level* interval covers its truth."""
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in QUALIFIED_TRUTH.items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any) -> str:
    """A human-readable posterior summary, truth in brackets, 95 % coverage flagged."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    covered = recovers_truth(run)
    lines = [
        (
            f"{attrs['ampere_engine']} on {attrs['ampere_backend']}: "
            f"{posterior.sizes['chain']} chain(s) x {posterior.sizes['draw']} draw(s)"
        ),
        "  posterior (truth in brackets, 95 % coverage flagged):",
    ]
    for name in sorted(posterior.data_vars):
        values = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(values, [2.5, 97.5])
        truth = QUALIFIED_TRUTH.get(name)
        flag = "" if truth is None else ("  ok" if covered[name] else "  MISS")
        bracket = "" if truth is None else f"  (truth {truth:+.6g})"
        lines.append(
            f"    {name:45s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )
    return "\n".join(lines)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", default="emcee", choices=list(ENGINES))
    parser.add_argument("--all", action="store_true", help="every engine this environment has")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--steps", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--burn-in", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--live-points", type=int, default=None, help="dynesty only")
    parser.add_argument("--dlogz", type=float, default=None, help="dynesty only")
    parser.add_argument("--budget", type=int, default=None, help="sbi only")
    parser.add_argument("--draws", type=int, default=None, help="sbi only")
    parser.add_argument(
        "--figures", type=Path, default=None, help="--all: directory for the overlay"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem(seed=args.seed)
    print(f"free parameters: {problem.parameters.free_names}")

    if args.all:
        started = time.perf_counter()
        runs = run_all(
            problem,
            walkers=args.walkers,
            steps=args.steps,
            burn_in=args.burn_in,
            live_points=args.live_points,
            dlogz=args.dlogz,
            budget=args.budget,
            draws=args.draws,
        )
        elapsed = time.perf_counter() - started
        for run in runs.values():
            print(report(run))
        print(f"  {elapsed:.1f} s wall clock, {len(runs)} engine(s)")
        if args.figures is not None:
            args.figures.mkdir(parents=True, exist_ok=True)
            figure = overlay_figure(runs, problem)
            figure.savefig(args.figures / "modified_blackbody_overlay.png")
        return 0

    started = time.perf_counter()
    run = fit(
        problem,
        engine=args.engine,
        walkers=args.walkers,
        steps=args.steps,
        burn_in=args.burn_in,
        live_points=args.live_points,
        dlogz=args.dlogz,
        budget=args.budget,
        draws=args.draws,
    )
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
