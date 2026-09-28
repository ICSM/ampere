"""The v2 twin of ``examples/minimal_working_example.py`` and its four
variants (``_dynesty``, ``_zeus``, ``_sbi``, ``_sbi_embedding``) -- W6.13 (1).

The legacy files stay exactly as they are (the characterisation anchors this
module does **not** touch); this package is the same model, the same data
and the same question, written once against :mod:`ampere.core` and the
reference backend, with the engine chosen by ``--engine`` rather than by
which file you ran.

The model
---------
:class:`LinearModel` is ``ASimpleModel``'s twin: ``F_nu = slope * lambda +
intercept``, in Jy, on the legacy grid ``10 ** np.linspace(0.0, 1.9, 2000)``
micron. Legacy's ``lims=[[-10, 10], [-10, 10]]`` flat prior becomes
``st.uniform(-10, 20)`` on each parameter (``scipy.stats.uniform``'s second
argument is a *width*, not an upper bound).

The data
--------
Three instruments bind the model's one channel, ``"sed"`` (distinct labels,
as :mod:`examples.sed_composition` requires once more than one instrument
shares a channel):

* ``"catalogue"`` -- :meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`
  on ``WISE_RSR_W1`` and ``SPITZER_MIPS_70``, the legacy script's two bands,
  10 % noise.
* ``"sl"`` and ``"ll"`` -- the two Spitzer IRS chunks
  (:mod:`.generators`'s :func:`~.generators.deduplicated_grids`, read from the
  tracked CASSIS file with :mod:`astropy.io.fits` rather than the legacy
  reader), each an :class:`~ampere.core.Instrument` of
  :class:`~ampere.backends.reference.Resample` onto that chunk's grid then
  :class:`~ampere.backends.reference.CalibrationScale`, 10 % noise.

The legacy noise model on each ``Spectrum`` was two numbers:
``calUnc=0.0025`` (the legacy docstring: "the sigma for a LogNormal
distribution" on the calibration scale factor) and ``scaleLengthPrior=0.01``
("the sigma, in micron, for a half-normal distribution" on the correlated
noise's RBF length scale -- the legacy default is a broader 0.1 micron half
normal; the script's own ``0.01`` is ten times tighter, i.e. it expects
almost no correlation beyond adjacent samples). Both map onto v2 pieces
directly: ``calUnc`` becomes each spectrum's own
:class:`~ampere.backends.reference.CalibrationScale`\\ (``st.lognorm(0.0025,
scale=1.0)``\\ ), and ``scaleLengthPrior`` becomes the flexible likelihood's
:class:`~ampere.core.GaussianProcessNoise`\\ (:class:`~ampere.core.Matern32`)
length-scale prior, ``st.halfnorm(scale=0.01)`` -- the same distribution
family and the same number, read directly off the legacy docstring, on the
axis v2 spells ``"spectral_axis"`` where legacy called it ``"wavelength"``.
The kernel's amplitude has no legacy number to match (``covWeightPrior``'s
units differ from a Matern32 amplitude's), so it carries a weakly informative
``st.halfnorm(scale=1.0)`` (Jy) instead. ``--no-gp`` swaps in
:class:`~ampere.core.IndependentNoise` on both spectra (the catalogue has no
correlated-noise counterpart in legacy either way, and stays
``IndependentNoise`` throughout).

Engines
-------
``--engine emcee|dynesty|zeus|sbi`` picks
:class:`~ampere.inference.EmceeEngine`, :class:`~ampere.inference.DynestyEngine`,
:class:`~ampere.inference.ZeusEngine` or :class:`~ampere.inference.SBIEngine` --
the twin of the legacy ``minimal_working_example{,_dynesty,_zeus,_sbi}.py``
choice of search. zeus ships in ``dev``; sbi needs the ``sbi`` extra
(``pixi run -e sbi``), and asking for it without that extra installed is
refused by name, naming the install command, rather than failing on an
``ImportError`` a caller has to trace back. ``--embedding`` adds the legacy
``_sbi_embedding`` variant's fully connected embedding network -- the same
dict, ``{"type": "FC", "num_hiddens": 100, "n_layers": 3, "output_dim":
20}`` -- to the ``sbi`` arm only (``ampere/inference/_sbi.py``'s
``embedding=`` vocabulary is the legacy dict's own spelling, read verbatim).

Coverage
--------
``recovers_truth`` checks the two model parameters and the two calibration
factors -- the four parameters with an injected truth; the GP's own
hyperparameters (when ``--no-gp`` is not given) have a prior but no truth to
recover against, exactly as legacy's ``calUnc``/``scaleLengthPrior`` numbers
were never "recovered" either. The full-budget run this module is accepted
on, and its recovered intervals, are recorded at the foot of this docstring.

::

    python -m examples.linear_sed                          # emcee, reference
    python -m examples.linear_sed --engine sbi --embedding  # pixi run -e sbi

Coverage run (Accept criterion)
--------------------------------
TODO: filled in by the coverage-run commit.
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

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
    Spectrum,
)
from ampere.backends.reference import CalibrationScale, Resample, SyntheticPhotometry
from ampere.inference import DynestyEngine, EmceeEngine, SBIEngine, ZeusEngine

from . import generators

__all__ = [
    "CALIBRATION_PRIOR",
    "DEFAULT_BURN_IN",
    "DEFAULT_STEPS",
    "DEFAULT_WALKERS",
    "EMBEDDING",
    "FILTERS",
    "GRID",
    "PHOTOMETRY_TABULATION",
    "QUALIFIED_TRUTH",
    "SBI_DEFAULT_BUDGET",
    "SBI_DEFAULT_DRAWS",
    "LinearModel",
    "build_instruments",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "recovers_truth",
    "report",
]

#: The legacy model's own grid (``examples/minimal_working_example.py``
#: line 94), micron.
GRID = 10 ** np.linspace(0.0, 1.9, 2000)

#: The legacy script's two well-separated bands (line 108).
FILTERS = ("WISE_RSR_W1", "SPITZER_MIPS_70")

#: The catalogue's response tabulation, micron -- within :data:`GRID`'s own
#: coverage (1 to ~79.4 micron), unlike :mod:`examples.sed_composition`'s,
#: which has a much broader model grid to work with.
PHOTOMETRY_TABULATION = np.geomspace(1.0, 79.0, 500)

#: ``calUnc=0.0025``'s v2 counterpart -- see the module docstring.
CALIBRATION_PRIOR = st.lognorm(0.0025, scale=1.0)

#: ``scaleLengthPrior=0.01``'s v2 counterpart -- see the module docstring.
GP_LENGTH_SCALE_PRIOR = st.halfnorm(scale=0.01)
#: No legacy number to match (see the module docstring); weakly informative.
GP_AMPLITUDE_PRIOR = st.halfnorm(scale=1.0)

#: The legacy ``_sbi_embedding`` variant's dict, verbatim (line ~202 there).
EMBEDDING: dict[str, Any] = {"type": "FC", "num_hiddens": 100, "n_layers": 3, "output_dim": 20}

#: The truth, qualified by the names a built ``FittingProblem`` samples.
QUALIFIED_TRUTH: dict[str, float] = {
    "model.slope": generators.TRUTH["slope"],
    "model.intercept": generators.TRUTH["intercept"],
    "sl.instrument.calibration_scale.scale": generators.CALIBRATION_TRUTH,
    "ll.instrument.calibration_scale.scale": generators.CALIBRATION_TRUTH,
}

# Budgets. The legacy emcee script's own (100 walkers, 150 steps, 100
# burn-in -- ``examples/minimal_working_example.py`` line ~226) is kept as
# the default for every ensemble engine, so a caller comparing engines is
# comparing them at one budget, not four.
DEFAULT_WALKERS = 100
DEFAULT_STEPS = 150
DEFAULT_BURN_IN = 100
DEFAULT_LIVE_POINTS = 200
SBI_DEFAULT_BUDGET = 1000
SBI_DEFAULT_DRAWS = 1000


def _as_parameter(name: str, spec: Any) -> Parameter:
    """A prior (a frozen scipy distribution) or a bare number, as a Parameter."""
    if isinstance(spec, Parameter):
        return spec
    if hasattr(spec, "logpdf") or hasattr(spec, "rvs"):
        return Parameter(name, spec)
    return Parameter(name, None, value=float(spec), fixed=True)


class LinearModel(Model):
    """``F_nu = slope * lambda + intercept`` -- the twin of ``ASimpleModel``.

    Parameters
    ----------
    wavelength
        The model's own grid, micron (a fallback -- negotiation replaces it
        the moment an instrument is bound; see :mod:`.generators`).
    slope, intercept
        A prior to fit it, or a number to hold it fixed.
    channel
        Name of the channel the emitted :class:`~ampere.core.Spectrum`
        appears under. Every instrument in this package binds ``"sed"``.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        slope: Any,
        intercept: Any,
        channel: str = "sed",
    ) -> None:
        grid = np.asarray(wavelength, dtype=float)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(f"LinearModel needs a 1-D, non-empty grid, got shape {grid.shape}.")
        self.channel = str(channel)
        self.register_buffer("wavelength", grid, unit=u.um)
        self.register_parameter(_as_parameter("slope", slope))
        self.register_parameter(_as_parameter("intercept", intercept))
        self._template: Spectrum | None = None

    def compile_for(self, requirements: Any) -> LinearModel:
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
        flux = float(ctx["slope"]) * grid + float(ctx["intercept"])
        spectrum = (
            self._template.with_values(flux)
            if self._template is not None
            else Spectrum(grid * u.um, flux, unit=u.Jy)
        )
        return ModelResult({self.channel: spectrum})


def build_model() -> LinearModel:
    """The one :class:`LinearModel`, with the legacy ``lims`` as flat priors."""
    return LinearModel(GRID, slope=st.uniform(-10.0, 20.0), intercept=st.uniform(-10.0, 20.0))


def build_instruments(*, gp: bool = True) -> tuple[Instrument, Instrument, Instrument]:
    """The catalogue and the two IRS chunks, distinct labels, one channel each."""
    sl_grid, ll_grid = generators.deduplicated_grids()
    catalogue = Instrument(
        [SyntheticPhotometry.from_library(list(FILTERS), PHOTOMETRY_TABULATION)],
        channel="sed",
        label="catalogue",
    )
    sl_instrument = Instrument(
        [Resample(sl_grid), CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed",
        label="sl",
    )
    ll_instrument = Instrument(
        [Resample(ll_grid), CalibrationScale(CALIBRATION_PRIOR)],
        channel="sed",
        label="ll",
    )
    return catalogue, sl_instrument, ll_instrument


def _likelihood(*, gp: bool) -> Likelihood:
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


def build_problem(*, gp: bool = True, seed: int = generators.SEED) -> FittingProblem:
    """The composed problem: one model, three datasets, distinct labels."""
    model = build_model()
    catalogue, sl_instrument, ll_instrument = build_instruments(gp=gp)
    observed_sl, observed_ll, observed_photometry = generators.synthetic_data(
        model, catalogue, sl_instrument, ll_instrument, seed=seed
    )
    datasets = DatasetCollection(
        {
            "catalogue": Dataset(
                observed_photometry,
                catalogue,
                likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
            ),
            "sl": Dataset(observed_sl, sl_instrument, likelihood=_likelihood(gp=gp)),
            "ll": Dataset(observed_ll, ll_instrument, likelihood=_likelihood(gp=gp)),
        }
    )
    return FittingProblem(model, datasets, seed=seed)


def _needs_sbi() -> None:
    if importlib.util.find_spec("sbi") is None:
        raise SystemExit(
            "engine='sbi' needs the 'sbi' extra, which is not installed in this environment. "
            "Install it with `pixi install -e sbi --frozen` and run this module again with "
            "`pixi run -e sbi python -m examples.linear_sed --engine sbi`."
        )


def fit(
    problem: FittingProblem,
    *,
    engine: str = "emcee",
    embedding: bool = False,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    live_points: int | None = None,
    budget: int | None = None,
    draws: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the named engine."""
    if embedding and engine != "sbi":
        raise SystemExit(f"--embedding only applies to --engine sbi, got --engine {engine!r}.")
    if engine == "emcee":
        run_engine = EmceeEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return run_engine.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "zeus":
        run_engine = ZeusEngine(problem, walkers=DEFAULT_WALKERS if walkers is None else walkers)
        return run_engine.run(
            DEFAULT_STEPS if steps is None else steps,
            burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
            progress=progress,
        )
    if engine == "dynesty":
        run_engine = DynestyEngine(
            problem, live_points=DEFAULT_LIVE_POINTS if live_points is None else live_points
        )
        return run_engine.run(progress=progress)
    if engine == "sbi":
        _needs_sbi()
        run_engine = SBIEngine(
            problem,
            budget=SBI_DEFAULT_BUDGET if budget is None else budget,
            embedding=EMBEDDING if embedding else None,
        )
        return run_engine.run(SBI_DEFAULT_DRAWS if draws is None else draws, progress=progress)
    raise SystemExit(f"unknown engine {engine!r}; choose emcee, dynesty, zeus or sbi.")


def recovers_truth(run: Any, *, level: float = 0.95) -> dict[str, bool]:
    """Whether each of the four qualified parameters' central *level* interval covers its truth."""
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
    parser.add_argument("--engine", default="emcee", choices=["emcee", "dynesty", "zeus", "sbi"])
    parser.add_argument("--embedding", action="store_true", help="sbi only")
    parser.add_argument("--no-gp", dest="gp", action="store_false", help="IndependentNoise instead")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--steps", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--burn-in", type=int, default=None, help="emcee/zeus only")
    parser.add_argument("--live-points", type=int, default=None, help="dynesty only")
    parser.add_argument("--budget", type=int, default=None, help="sbi only")
    parser.add_argument("--draws", type=int, default=None, help="sbi only")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem(gp=args.gp, seed=args.seed)
    print(f"free parameters: {problem.parameters.free_names}")

    started = time.perf_counter()
    run = fit(
        problem,
        engine=args.engine,
        embedding=args.embedding,
        walkers=args.walkers,
        steps=args.steps,
        burn_in=args.burn_in,
        live_points=args.live_points,
        budget=args.budget,
        draws=args.draws,
    )
    elapsed = time.perf_counter() - started

    print(report(run))
    print(f"  {elapsed:.1f} s wall clock")
    return 0
