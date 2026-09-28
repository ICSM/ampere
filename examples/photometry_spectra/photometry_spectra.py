"""One :class:`~ampere.backends.reference.ModifiedBlackBody`, fitted to a
photometric catalogue *and* two spectra of different resolving power through
:class:`~ampere.core.Instrument`'s LSF-then-resample-then-calibrate chain --
memo `docs/design/modalities/spectrum_photometry.md` §4's "several
observations" case, made runnable (W6.2). See :doc:`the tutorial page
</photometry_spectra>` for the walk-through; this module is the code it
walks through. :mod:`examples.sed_composition` is this package's sibling and
its own tutorial page is the one-model/two-instrument groundwork this page
does not repeat.

The physics
-----------
The model is exactly :mod:`examples.sed_composition.sed_composition`'s
:class:`~ampere.backends.reference.ModifiedBlackBody` -- same priors, same
``channels="sed"``. What is new is the third observation and what it forces:
a spectrum's absolute flux calibration is never perfectly known, so every
real spectrograph carries an uncertain multiplicative factor relative to a
photometric catalogue (or to another spectrograph) -- :data:`FILTERS` sees
no such factor because a catalogue *is* the calibration reference here.

The three instruments
----------------------
**sl** -- a short-wavelength slit spectrograph:
:class:`~ampere.backends.reference.LSFConvolution` at constant resolving
power :data:`SL_RESOLVING_POWER` (~100), :class:`~ampere.backends.reference.Resample`
onto :data:`SL_OBSERVED_WAVELENGTH` (5-14 micron), and a
:class:`~ampere.backends.reference.CalibrationScale` nuisance parameter --
the IRS "SL" module's shape, synthetic.

**ll** -- a long-wavelength slit spectrograph: the same chain at
:data:`LL_RESOLVING_POWER` (~60) over :data:`LL_OBSERVED_WAVELENGTH`
(14-38 micron) -- the IRS "LL" module's shape. The synthetic data
(:mod:`.generators`) also injects a smooth Gaussian bump into this spectrum
alone -- a feature no modified blackbody can express.

**catalogue** -- :meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`,
six bundled mid/far-infrared filters spanning roughly 20-160 micron (see
:data:`FILTERS`) tabulated on a 500-point grid.

All three bind the model's one channel, ``"sed"``, with distinct labels
(``"sl"``, ``"ll"``, ``"catalogue"``) -- :doc:`the sibling tutorial page
</sed_composition>` walks through what happens (and why) without them.

Per spectrum, or shared
------------------------
:func:`build_problem` takes ``tie: bool``. Untied (the default), each
spectrograph's calibration factor is its own free parameter --
``sl.instrument.calibration_scale.scale`` and
``ll.instrument.calibration_scale.scale``. Tied, one
:class:`~ampere.core.Tie` collapses the two into a single
``"calibration"`` parameter, registered on the
:class:`~ampere.core.dataset.FittingProblem` through its ``ties=`` keyword
(``dataset.py``'s own docstring: "Declared here because this is the only
level at which every site is visible; a tie names sites by their full merged
path"). :data:`examples.photometry_spectra.generators.CALIBRATION_TRUTH`
gives the two spectrographs *different* factors (0.92, 1.08) on purpose, so
the tied fit has a genuine conflict to resolve rather than two numbers that
happen to agree -- the shared factor lands between the two truths, and
:doc:`the tutorial page </photometry_spectra>` shows what that costs each
spectrum's residuals.

The complement: what calibration does not explain
----------------------------------------------------
:func:`build_problem` also takes ``gp: bool``. ``False`` (the default) gives
every dataset :class:`~ampere.core.IndependentNoise` -- the ordinary
chi-square. ``True`` keeps :class:`~ampere.core.IndependentNoise` on the
photometry but gives each spectrum its own
:class:`~ampere.core.GaussianProcessNoise` (a Matérn-3/2 in wavelength, log-uniform
"shrinkage-free" priors on the amplitude and length scale, following
:class:`~ampere.core.GaussianProcessNoise`'s own docstring example -- a prior
with mass concentrated near zero would pull a genuinely present feature back
towards "no GP needed", which is the opposite of what this page's point
needs). The tutorial page's measured contrast: under ``gp=False`` the
injected bump has nowhere to go but the calibration factor and the
temperature, biasing both; under ``gp=True`` the physical parameters and the
calibration factors recover, and the GP's posterior mean localises the bump
on the ``"ll"`` dataset (:func:`write_figures`'s ``gp_localisation.png`).

Run it::

    python -m examples.photometry_spectra                  # untied, independent noise, ~1-2 minutes
    python -m examples.photometry_spectra --tie
    python -m examples.photometry_spectra --gp
    python -m examples.photometry_spectra --tie --gp
    python -m examples.photometry_spectra --figures /tmp/photometry_spectra_figures

Engines
-------
emcee on the reference backend only -- this is a tutorial about composing a
third observation and a tie, not a backend study (:mod:`examples.sed_composition`
already carries the backend-flag story). :class:`~ampere.inference.EmceeEngine`'s
ensemble size is left at its own default (four per free dimension, at least
eight -- ``ampere/inference/engine.py``'s ``_default_walkers``) rather than a
single hardcoded constant, because the four ``tie``/``gp`` combinations have
different dimensionality (4, 5, 8 or 9 free parameters) and a walker count
tuned to one would violate emcee's own ``2 x n_dim`` floor on another.

:mod:`tests.examples.test_photometry_spectra` is this module's own coverage:
a fast, always-on suite at a tiny budget. The full-budget recovery this
item's tutorial page quotes is not part of that suite -- see its docstring
for why -- and is instead verified by running this module directly, once per
``tie``/``gp`` combination, through the gate lock; the branch report has the
numbers from doing exactly that.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
import time
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.backends.reference import (
    CalibrationScale,
    LSFConvolution,
    ModifiedBlackBody,
    Resample,
    SyntheticPhotometry,
)
from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    Likelihood,
    Matern32,
    Tie,
)
from ampere.inference import EmceeEngine

from . import generators

__all__ = [
    "CALIBRATION_PRIOR",
    "DEFAULT_BURN_IN",
    "DEFAULT_STEPS",
    "FILTERS",
    "GP_AMPLITUDE_PRIOR",
    "GP_LENGTH_SCALE_PRIOR",
    "GRID",
    "LL_OBSERVED_WAVELENGTH",
    "LL_RESOLVING_POWER",
    "PHOTOMETRY_TABULATION",
    "SL_OBSERVED_WAVELENGTH",
    "SL_RESOLVING_POWER",
    "build_instruments",
    "build_model",
    "build_problem",
    "fit",
    "main",
    "qualified_truth",
    "recovers_truth",
    "report",
    "write_figures",
]

#: The model's own grid -- a fallback only: negotiation replaces it the
#: moment an instrument is bound (see :mod:`.generators`).
GRID = np.geomspace(1.0, 200.0, 2000)

#: "sl"'s own observed grid, micron -- the Spitzer IRS "Short-Low" module's
#: shape, synthetic.
SL_OBSERVED_WAVELENGTH = np.linspace(5.0, 14.0, 90)

#: "ll"'s own observed grid, micron -- the IRS "Long-Low" module's shape.
LL_OBSERVED_WAVELENGTH = np.linspace(14.0, 38.0, 130)

SL_RESOLVING_POWER = 100.0
LL_RESOLVING_POWER = 60.0

#: The photometry step's response tabulation, micron -- becomes its exact
#: ``points=`` requirement. Wide enough to cover every chosen filter's
#: passband core with margin (checked once against the loaded library --
#: see :data:`FILTERS`).
PHOTOMETRY_TABULATION = np.geomspace(10.0, 200.0, 500)

#: Six bundled mid/far-infrared filters spanning roughly 20-160 micron --
#: the sibling example's ``FILTERS`` (dropping ``WISE_RSR_W3`` at 12 micron,
#: below this page's 20 micron floor) plus two AKARI/Herschel names reaching
#: further into the far infrared. Checked once by loading
#: ``ampere_allfilters.hd5`` directly: pivot wavelengths (micron) are
#: WISE_RSR_W4 22.2, IRAS_60 58.2, SPITZER_MIPS_70 70.9, HERSCHEL_PACS_100
#: 100.6, AKARI_FIS_WIDEL 146.4, HERSCHEL_PACS_160 161.0 -- the bundled
#: library's longest-wavelength filters short of Herschel SPIRE (~250-500
#: micron, well past this source's flux), so 161 micron is as close to the
#: nominal 200 micron upper bound as the library's mid/far-infrared bands get.
FILTERS = (
    "WISE_RSR_W4",
    "IRAS_60",
    "SPITZER_MIPS_70",
    "HERSCHEL_PACS_100",
    "AKARI_FIS_WIDEL",
    "HERSCHEL_PACS_160",
)

#: The natural calibration prior: log-normal about 1, since a calibration
#: factor is positive and multiplicative
#: (:class:`~ampere.backends.reference.CalibrationScale`'s own docstring).
#: Shared by both spectrographs' own steps and, when tying, by the
#: :class:`~ampere.core.Tie` itself.
CALIBRATION_PRIOR = st.lognorm(0.05, scale=1.0)

#: The GP kernel's amplitude prior, Jy -- log-uniform ("shrinkage-free": flat
#: in log-amplitude, so it has no pull towards zero the way a half-normal or
#: half-Cauchy would), following :class:`~ampere.core.GaussianProcessNoise`'s
#: own docstring example. Bracketing the injected bump's peak (about 0.43 Jy
#: -- :data:`generators.BUMP_FRACTION` applied to the "ll" continuum, 2-5.5
#: Jy) without reaching so high that the GP could absorb the continuum's own
#: scale.
GP_AMPLITUDE_PRIOR = st.loguniform(1.0e-2, 3.0)

#: The GP kernel's length-scale prior, micron -- log-uniform for the same
#: reason, bracketing the injected bump's 4 micron width
#: (:data:`generators.BUMP_WIDTH`) while capped well below either
#: spectrograph's own span (9 micron for "sl", 24 for "ll"): an unbounded
#: upper end would let a long-length-scale, moderate-amplitude draw mimic a
#: broadband multiplicative recalibration over the whole spectrum -- exactly
#: the degeneracy with the calibration factor this page's GP arm is meant to
#: avoid, not reproduce (found empirically: an upper bound of 100 micron let
#: the "ll" calibration factor's 95 % interval miss its truth).
GP_LENGTH_SCALE_PRIOR = st.loguniform(5.0e-1, 1.2e1)

# Budgets. Four free-parameter counts are possible (4, 5, 8 or 9, depending
# on tie/gp), so walkers are left to EmceeEngine's own default rather than a
# single constant tuned to one of them -- see the module docstring's
# "Engines" section.
DEFAULT_STEPS = 600
DEFAULT_BURN_IN = 250


def build_model() -> ModifiedBlackBody:
    """The one :class:`~ampere.backends.reference.ModifiedBlackBody`, exactly as
    :func:`examples.sed_composition.sed_composition.build_model` declares it
    on the reference backend."""
    return ModifiedBlackBody(
        GRID,
        temperature=st.uniform(50.0, 400.0),
        beta=st.uniform(0.5, 2.5),
        scale=st.loguniform(1e-16, 1e-14),
        channels="sed",
    )


def build_instruments() -> tuple[Instrument, Instrument, Instrument]:
    """The two spectrographs (``"sl"``, ``"ll"``) and the catalogue camera.

    All three bind channel ``"sed"``; distinct labels are passed explicitly,
    the fix :doc:`the sibling tutorial page </sed_composition>` walks
    through for exactly this situation.
    """
    sl = Instrument(
        [
            LSFConvolution(resolving_power=SL_RESOLVING_POWER),
            Resample(SL_OBSERVED_WAVELENGTH),
            CalibrationScale(CALIBRATION_PRIOR),
        ],
        channel="sed",
        label="sl",
    )
    ll = Instrument(
        [
            LSFConvolution(resolving_power=LL_RESOLVING_POWER),
            Resample(LL_OBSERVED_WAVELENGTH),
            CalibrationScale(CALIBRATION_PRIOR),
        ],
        channel="sed",
        label="ll",
    )
    camera = Instrument(
        [SyntheticPhotometry.from_library(list(FILTERS), PHOTOMETRY_TABULATION)],
        channel="sed",
        label="catalogue",
    )
    return sl, ll, camera


def _spectrum_likelihood(*, gp: bool) -> Likelihood:
    """Independent noise, or the flexible likelihood, for one spectrum dataset.

    Each call builds a fresh :class:`~ampere.core.GaussianProcessNoise` (and
    so fresh ``Parameter`` objects) -- "sl" and "ll" get their own amplitude
    and length scale, not a shared pair, since nothing says the two
    spectrographs' unmodelled structure has the same size or shape.
    """
    if not gp:
        return Likelihood(GaussianFamily(), IndependentNoise())
    kernel = Matern32(
        GP_AMPLITUDE_PRIOR,
        GP_LENGTH_SCALE_PRIOR,
        amplitude_unit=u.Jy,
        length_scale_unit=u.um,
    )
    return Likelihood(GaussianFamily(), GaussianProcessNoise(kernel))


def build_problem(
    *, tie: bool = False, gp: bool = False, seed: int = generators.SEED
) -> FittingProblem:
    """The composed problem: one model, three datasets, distinct labels.

    Parameters
    ----------
    tie
        ``False`` (the default): ``sl.instrument.calibration_scale.scale``
        and ``ll.instrument.calibration_scale.scale`` are independent free
        parameters. ``True``: one :class:`~ampere.core.Tie` named
        ``"calibration"`` collapses them into a single shared factor,
        registered on the :class:`~ampere.core.dataset.FittingProblem`
        through ``ties=`` -- the only level at which both qualified sites are
        visible.
    gp
        ``False`` (the default): every dataset gets
        :class:`~ampere.core.IndependentNoise`. ``True``: the two spectra
        get :class:`~ampere.core.GaussianProcessNoise` instead (the
        photometry keeps :class:`~ampere.core.IndependentNoise` regardless --
        three or four broadband points carry no exploitable correlation
        structure at this resolution, the same reasoning
        ``spectrum_photometry.md`` §4 gives).
    seed
        Seeds both the synthetic data's noise draw and the problem's own
        sampler streams.
    """
    model = build_model()
    sl, ll, camera = build_instruments()
    observed_sl, observed_ll, observed_photometry = generators.synthetic_data(
        model, sl, ll, camera, seed=seed
    )

    datasets = DatasetCollection(
        {
            "sl": Dataset(observed_sl, sl, likelihood=_spectrum_likelihood(gp=gp)),
            "ll": Dataset(observed_ll, ll, likelihood=_spectrum_likelihood(gp=gp)),
            "catalogue": Dataset(
                observed_photometry,
                camera,
                likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
            ),
        }
    )

    ties = (
        (
            Tie(
                "calibration",
                (
                    "sl.instrument.calibration_scale.scale",
                    "ll.instrument.calibration_scale.scale",
                ),
                prior=CALIBRATION_PRIOR,
            ),
        )
        if tie
        else ()
    )
    return FittingProblem(model, datasets, ties=ties, seed=seed)


#: Rough starting guesses for the GP hyperparameters -- not a "truth" (the
#: injected bump is not itself a draw from a Matérn-3/2), just a plausible
#: point inside :data:`GP_AMPLITUDE_PRIOR`/:data:`GP_LENGTH_SCALE_PRIOR`'s
#: support to start the ball from.
_GP_AMPLITUDE_GUESS = 0.3
_GP_LENGTH_SCALE_GUESS = 3.0

#: The initial ball's per-walker relative jitter (log-normal sigma).
_INITIAL_BALL_JITTER = 0.05


def _initial_positions(problem: FittingProblem, walkers: int) -> np.ndarray:
    """A tight ball around a good starting guess, not the raw prior span.

    :meth:`~ampere.inference.engine.Engine.initial_positions` -- what
    ``fit`` would otherwise fall back on -- draws each walker from the joint
    prior and keeps the first draw that scores *finitely*, which accepts a
    technically-scoreable but astronomically improbable point exactly when a
    prior spans much more range than the posterior's own footprint (here,
    ``model.scale`` alone spans two decades). Found empirically on this
    composition: about one walker in ten started that way never accepts a
    single proposal in thousands of subsequent steps -- the affine-invariant
    stretch move cannot climb back from so far off -- which then corrupts
    every flattened statistic :func:`report`/:func:`recovers_truth` compute.

    Starting near a known-good point instead is ordinary MCMC practice, and
    this tutorial is exactly the case where it costs nothing to be honest
    about: the synthetic truth is *known* (:mod:`.generators`), so each
    walker starts within a few percent of it (or, for the two GP
    hyperparameters, of a plausible point in their prior) rather than
    wherever an unconstrained prior draw happened to land.
    """
    names = problem.parameters.free_names
    centre: dict[str, float] = {
        "model.temperature": generators.TRUTH["temperature"],
        "model.beta": generators.TRUTH["beta"],
        "model.scale": generators.TRUTH["scale"],
    }
    if "calibration" in names:
        centre["calibration"] = (
            generators.CALIBRATION_TRUTH["sl"] + generators.CALIBRATION_TRUTH["ll"]
        ) / 2.0
    else:
        centre["sl.instrument.calibration_scale.scale"] = generators.CALIBRATION_TRUTH["sl"]
        centre["ll.instrument.calibration_scale.scale"] = generators.CALIBRATION_TRUTH["ll"]
    for label in ("sl", "ll"):
        amplitude_name = f"{label}.likelihood.amplitude"
        length_scale_name = f"{label}.likelihood.length_scale"
        if amplitude_name in names:
            centre[amplitude_name] = _GP_AMPLITUDE_GUESS
        if length_scale_name in names:
            centre[length_scale_name] = _GP_LENGTH_SCALE_GUESS

    if set(centre) != set(names):
        raise AssertionError(
            f"_initial_positions does not know every free parameter: has {sorted(centre)}, "
            f"problem declares {sorted(names)}."
        )

    # Every centred quantity here is positive (temperatures, scales,
    # calibration factors, GP hyperparameters), so a log-normal jitter keeps
    # every draw inside its parameter's positive support without a rejection
    # loop.
    rng = problem.rng("photometry_spectra.initial_ball")
    positions = np.empty((walkers, len(names)), dtype=float)
    for row in range(walkers):
        jittered = {
            name: centre[name] * float(rng.lognormal(0.0, _INITIAL_BALL_JITTER)) for name in names
        }
        positions[row] = problem.parameters.pack(jittered)
    return positions


def fit(
    problem: FittingProblem,
    *,
    walkers: int | None = None,
    steps: int | None = None,
    burn_in: int | None = None,
    progress: bool = False,
) -> Any:
    """Sample *problem* with :class:`~ampere.inference.EmceeEngine`.

    ``walkers=None`` (the default) leaves the ensemble size to the engine's
    own ``_default_walkers`` -- see the module docstring's "Engines" section
    for why this module does not hardcode one, unlike
    :func:`examples.sed_composition.sed_composition.fit`. Walkers start in
    :func:`_initial_positions`'s tight ball rather than the engine's own
    prior-draw default -- see that function's docstring for why.
    """
    engine = EmceeEngine(problem, walkers=walkers)
    return engine.run(
        DEFAULT_STEPS if steps is None else steps,
        burn_in=DEFAULT_BURN_IN if burn_in is None else burn_in,
        initial=_initial_positions(problem, engine.walkers),
        progress=progress,
    )


def qualified_truth(*, tie: bool) -> dict[str, float]:
    """The truth, qualified by the names a built :class:`~ampere.core.dataset.FittingProblem`
    actually samples -- ``problem.parameters.free_names`` once the three
    datasets and the model are merged.

    Under ``tie=True`` the one shared ``"calibration"`` parameter's coverage
    target is the *mean* of the two spectrographs' truths (0.92, 1.08) --
    not either one alone, since a single shared factor cannot recover both.
    The tutorial page's point is exactly that this is the best a tied fit can
    do, and that "recovers" the mean is a weaker, different claim from
    "recovers" either spectrograph's own miscalibration.
    """
    truth = {
        "model.temperature": generators.TRUTH["temperature"],
        "model.beta": generators.TRUTH["beta"],
        "model.scale": generators.TRUTH["scale"],
    }
    if tie:
        truth["calibration"] = (
            generators.CALIBRATION_TRUTH["sl"] + generators.CALIBRATION_TRUTH["ll"]
        ) / 2.0
    else:
        truth["sl.instrument.calibration_scale.scale"] = generators.CALIBRATION_TRUTH["sl"]
        truth["ll.instrument.calibration_scale.scale"] = generators.CALIBRATION_TRUTH["ll"]
    return truth


def recovers_truth(run: Any, *, tie: bool, level: float = 0.95) -> dict[str, bool]:
    """Whether each parameter's central *level* interval contains its truth.

    Only the model's three parameters and the calibration factor(s) have a
    truth to check -- a fitted GP's amplitude and length scale are not
    checked against a "truth" here, because nothing in :mod:`.generators`
    claims the injected bump's shape *is* a Matérn-3/2 draw; the GP's job is
    to absorb it, not to recover parameters that were never simulated from
    its own prior.
    """
    posterior = run["posterior"].dataset
    tail = (1.0 - level) / 2.0 * 100.0
    covered: dict[str, bool] = {}
    for name, truth in qualified_truth(tie=tie).items():
        draws = np.asarray(posterior[name], dtype=float).ravel()
        lower, upper = np.percentile(draws, [tail, 100.0 - tail])
        covered[name] = bool(lower <= truth <= upper)
    return covered


def report(run: Any, *, tie: bool) -> str:
    """A human-readable posterior summary, truth in brackets, 95 % coverage flagged."""
    attrs = run.attrs
    posterior = run["posterior"].dataset
    truth_map = qualified_truth(tie=tie)
    covered = recovers_truth(run, tie=tie)
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
        truth = truth_map.get(name)
        flag = "" if truth is None else ("  ok" if covered[name] else "  MISS")
        bracket = "" if truth is None else f"  (truth {truth:+.6g})"
        lines.append(
            f"    {name:45s} {values.mean():+.6g} +- {values.std():.3g}   "
            f"95%[{lower:+.6g}, {upper:+.6g}]{bracket}{flag}"
        )
    return "\n".join(lines)


def _pyplot() -> Any:
    """``matplotlib.pyplot``, imported here rather than at module scope --
    see :mod:`examples.m2_misspecification.figures`'s ``_pyplot`` for why."""
    import matplotlib.pyplot as plt

    return plt


def write_figures(
    run: Any, problem: FittingProblem, directory: str | pathlib.Path, *, gp: bool, thin: int = 10
) -> list[pathlib.Path]:
    """Render ``corner.png``, ``posterior_predictive.png`` and (if *gp*)
    ``gp_localisation.png`` into *directory*.

    Nothing here is committed (AGENTS.md ground rule 7) -- the caller decides
    where these go, exactly as :mod:`examples.m2_misspecification.figures`'s
    ``save_result_figures`` does. ``gp_localisation.png`` is only written
    under ``gp=True``: :func:`~ampere.results.gp_localisation` refuses a
    dataset with no :class:`~ampere.core.GaussianProcessNoise` by name, and
    there is nothing to localise when every dataset is independent noise.
    """
    from ampere.results import add_posterior_predictive, plot_corner, plot_posterior_predictive

    plt = _pyplot()
    target = pathlib.Path(directory)
    target.mkdir(parents=True, exist_ok=True)
    written: list[pathlib.Path] = []

    def _save(figure: Any, stem: str) -> None:
        path = target / f"{stem}.png"
        figure.savefig(path, bbox_inches="tight", dpi=120)
        plt.close(figure)
        written.append(path)

    _save(plot_corner(run), "corner")

    run = add_posterior_predictive(run, problem, thin=thin)
    _save(
        plot_posterior_predictive(run, datasets=["sl", "ll", "catalogue"]), "posterior_predictive"
    )

    if gp:
        from ampere.results import gp_localisation, plot_gp_localisation

        run = gp_localisation(run, problem, datasets=["ll"], thin=thin)
        _save(plot_gp_localisation(run, datasets=["ll"]), "gp_localisation")

    return written


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tie", action="store_true", help="one shared calibration factor")
    parser.add_argument("--gp", action="store_true", help="the flexible likelihood on both spectra")
    parser.add_argument("--seed", type=int, default=generators.SEED)
    parser.add_argument("--walkers", type=int, default=None)
    parser.add_argument("--steps", type=int, default=None)
    parser.add_argument("--burn-in", type=int, default=None)
    parser.add_argument(
        "--figures", type=str, default=None, metavar="DIR", help="write the three figures here"
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(sys.argv[1:] if argv is None else argv)

    problem = build_problem(tie=args.tie, gp=args.gp, seed=args.seed)
    req = problem.requirements["model"]["sed"]
    print(f"negotiated channels: {list(problem.requirements['model'])}")
    print(f"sources asking of channel 'sed': {req.sources}")
    for axis, axis_req in req.axes.items():
        print(f"  axis {axis}: {axis_req}")
    print(f"free parameters: {problem.parameters.free_names}")

    started = time.perf_counter()
    run = fit(problem, walkers=args.walkers, steps=args.steps, burn_in=args.burn_in)
    elapsed = time.perf_counter() - started

    print(report(run, tie=args.tie))
    print(f"  {elapsed:.1f} s wall clock")

    if args.figures is not None:
        written = write_figures(run, problem, args.figures, gp=args.gp)
        print(f"figures written to {args.figures}: {[p.name for p in written]}")
    return 0
