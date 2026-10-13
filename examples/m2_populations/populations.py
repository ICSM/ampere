"""The scenario, the three arms and the driver (**W7.2**).

One object class, :data:`N_SPECTRA` spectra of it, and the question of how
their GP amplitudes should be declared. Everything a spectrum *is* comes from
the milestone study: the generator
(:func:`~examples.m2_misspecification.generators.generate`, at distinct noise
seeds), the toy and its truth, the flexible likelihood
(:func:`~examples.m2_misspecification.study.build_likelihood`, the study's
Matérn-3/2 with its half-normal priors), the engine dispatch
(:func:`~examples.m2_misspecification.study.run`: emcee on the reference path,
NUTS on torch and jax) and the reduction to
:class:`~examples.m2_misspecification.study.Summary`. What this module adds is
the *collection*: N datasets, one model instance each, and three declarations
of their amplitudes over the same data.

What varies between the spectra
-------------------------------
The truth does not — one object class, the M2 toy at
:data:`~examples.m2_misspecification.model.TRUTH` — and neither does the grid.
Three things do:

* the **noise realisation**: spectrum *i* is drawn at seed
  ``DATA_SEED + i``;
* the **deviation**: :data:`DEVIATED` (three of eight) carry the study's
  ``strong_sharp`` scenario — a 12 % Gaussian emission line at 0.8630 µm the
  model has no term for — and the rest none;
* the **noise level** of one clean member, :data:`POORLY_CONSTRAINED`, observed
  at :data:`POOR_NOISE_FRACTION` (five times the study's 1 %), so that its own
  data say little about its GP amplitude. That member is where a population
  earns its keep, and :func:`shrinkage` is the statistic for it.

The three arms
--------------
``independent``
    The study's own declaration, once per spectrum: ``amplitude ~
    HalfNormal(0.3 Jy)``, N of them.
``tied``
    A :class:`~ampere.core.Tie` collapsing every ``d*.likelihood.amplitude``
    into one ``amplitude`` with the same prior.
``population``
    The memo's §4, verbatim but for the unit: ``z_i ~ Normal(0, 1)`` internal,
    ``amplitude_i = exp(mu + sigma * z_i)`` derived (in Jy, the unit every
    spectrum's amplitude carries), ``mu ~ Normal(-3, 2)``, ``sigma ~
    HalfNormal(1)``. ``mu = -3`` is ``exp(-3) = 0.05`` Jy, the order of the
    amplitudes M2 fits where there is something to absorb; a prior width of 2
    in the log reaches the near-zero amplitudes of a clean spectrum
    (``exp(-7) = 0.001`` Jy) inside two standard deviations.

Every arm keeps each spectrum's own length scale (the study's
``HalfNormal(0.003 µm)``) and its own four model parameters, so the arms differ
in the amplitude declaration and nothing else.
"""

from __future__ import annotations

import argparse
import dataclasses
import time
from collections.abc import Sequence
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    Derived,
    FittingProblem,
    Parameter,
    Tie,
)
from ampere.results import gp_localisation, gp_localisation_score
from examples.m2_misspecification.generators import (
    NOISE_FRACTION,
    SEED as GENERATOR_SEED,
    SyntheticSpectrum,
    generate,
    scenario_named,
)
from examples.m2_misspecification.model import PARAMETER_NAMES, TRUTH
from examples.m2_misspecification.study import (
    BACKENDS,
    SEED,
    EmceeBudget,
    NutsBudget,
    Summary,
    build_likelihood,
    model_for,
    run,
    summarise,
)

__all__ = [
    "ARMS",
    "DATA_SEED",
    "DEVIATED",
    "DEVIATION",
    "N_SPECTRA",
    "POORLY_CONSTRAINED",
    "POOR_NOISE_FRACTION",
    "POPULATION_EMCEE",
    "POPULATION_FULL_NUTS",
    "POPULATION_LABEL",
    "POPULATION_NUTS",
    "QUICK_EMCEE",
    "Comparison",
    "Fit",
    "Localisation",
    "amplitude_draws",
    "budget_for",
    "build_arm",
    "compare",
    "dataset_label",
    "generate_population",
    "hyperpriors",
    "localise",
    "main",
    "members",
    "model_label",
    "physical_summaries",
    "sample_spread",
]

#: How many spectra the scenario has.
N_SPECTRA = 8

#: The spectra that carry the deviation: three of eight, spread through the
#: list so that no arm's result can depend on the deviated ones being first.
DEVIATED: tuple[int, ...] = (1, 4, 6)

#: The deviation they carry: the study's own sharp line (a 12 % Gaussian
#: emission line at 0.8630 µm, sigma 0.00035 µm), the one scenario of the four
#: with a location, so that GP localisation has an answer to be held to.
DEVIATION = "strong_sharp"

#: The clean member observed at five times the noise: its own data barely
#: constrain its GP amplitude, so its amplitude posterior is where a
#: population's shrinkage shows.
POORLY_CONSTRAINED = 7
POOR_NOISE_FRACTION = 5.0 * NOISE_FRACTION

#: The data seed of spectrum 0; spectrum *i* is drawn at ``DATA_SEED + i``.
#: The study's own generator seed, so spectrum 0 is the study's 200-point
#: control spectrum exactly.
DATA_SEED = GENERATOR_SEED

#: The three declarations, in the order every table uses them.
ARMS: tuple[str, ...] = ("independent", "tied", "population")

#: The population's (and therefore the plate dimension's) name.
POPULATION_LABEL = "gp"

#: The tie's name in the ``tied`` arm.
TIED_NAME = "amplitude"

#: The population's hyperprior on the log-amplitude's location and spread.
MU_LOC, MU_SCALE = -3.0, 2.0
SIGMA_SCALE = 1.0

# ---------------------------------------------------------------------------
# Budgets. The widest arm has 50 free dimensions (8 x four model parameters +
# 8 length scales + 8 z + mu and sigma; the independent arm 48, the tied 41),
# which sets every number below.
# ---------------------------------------------------------------------------

#: The per-PR budget of ``tests/m2/test_populations.py``'s native rows: the
#: study's own :data:`~examples.m2_misspecification.study.TEST_NUTS`. NUTS'
#: cost grows with the dimension far more gently than an ensemble's, which is
#: why the science rows are the gradient path's: about two and a half minutes
#: per arm on jax at 200 points.
POPULATION_NUTS = NutsBudget(draws=400, warmup=400, chains=2)

#: The marked (``m2_full``) and CLI budget on the gradient path.
POPULATION_FULL_NUTS = NutsBudget(draws=1_000, warmup=1_000, chains=2)

#: The reference path's budget, for the ``m2_full`` row and the CLI. emcee's
#: stretch move needs more than twice the dimension in walkers, so 112 covers
#: the widest arm with room. An evaluation of the 50-dimensional problem costs
#: about 6.5 ms on the numpy path (eight GP solves, and fifty scipy prior
#: densities, which are half of it), so this is about twenty minutes per arm:
#: the reason the per-PR science rows are NUTS' and this is behind the marker.
POPULATION_EMCEE = EmceeBudget(walkers=112, steps=1_600, burn_in=800)

#: ``--quick`` on the reference path: a first look, not a measurement.
QUICK_EMCEE = EmceeBudget(walkers=112, steps=300, burn_in=150)


def budget_for(backend: str, *, quick: bool = False) -> EmceeBudget | NutsBudget:
    """The budget for *backend*: :data:`POPULATION_NUTS` or :data:`POPULATION_EMCEE` and kin."""
    if backend == "reference":
        return QUICK_EMCEE if quick else POPULATION_EMCEE
    return POPULATION_NUTS if quick else POPULATION_FULL_NUTS


def dataset_label(index: int) -> str:
    """Spectrum *index*'s dataset label: ``d0``, ``d1``, ..."""
    return f"d{index}"


def model_label(index: int) -> str:
    """Spectrum *index*'s model label: ``obj0``, ``obj1``, ..."""
    return f"obj{index}"


def generate_population(
    *,
    size: int = 200,
    seed: int = DATA_SEED,
    count: int = N_SPECTRA,
    deviated: Sequence[int] = DEVIATED,
    poorly_constrained: int | None = POORLY_CONSTRAINED,
) -> tuple[SyntheticSpectrum, ...]:
    """The scenario's spectra: the study's generator at seeds ``seed + i``.

    Deterministic in its arguments, as the generator is. Spectrum *i* carries
    :data:`DEVIATION` if *i* is in *deviated* and nothing otherwise, at the
    study's 1 % noise — except *poorly_constrained*, at
    :data:`POOR_NOISE_FRACTION`.
    """
    chosen = set(deviated)
    if poorly_constrained is not None and poorly_constrained in chosen:
        raise ValueError(
            f"spectrum {poorly_constrained} cannot be both deviated and the poorly constrained "
            f"clean member: the shrinkage statistic is about a member with nothing to absorb."
        )
    return tuple(
        generate(
            DEVIATION if index in chosen else "none",
            size=size,
            seed=seed + index,
            noise_fraction=POOR_NOISE_FRACTION if index == poorly_constrained else NOISE_FRACTION,
        )
        for index in range(count)
    )


def members() -> list[Parameter]:
    """The memo's §4 members: ``z`` internal, ``amplitude`` derived from it.

    ``z`` is declared by no spectrum's likelihood, so it is internal: sampled
    on the population's own component and routed nowhere (W7.0's rule). The
    amplitude carries **Jy** because each spectrum's ``likelihood.amplitude``
    does — the study's kernel checks it against the observed flux — and a
    population member must agree in unit with the leaf it replaces; the memo's
    §4 spells the declaration without a unit because its sketch's kernel had
    none.
    """
    return [
        Parameter("z", st.norm(0.0, 1.0)),
        Parameter("amplitude", Derived("exp(mu + sigma * z)"), unit=u.Jy),
    ]


def hyperpriors() -> list[Parameter]:
    """``mu ~ Normal(-3, 2)``, ``sigma ~ HalfNormal(1)``: the memo's numbers."""
    return [
        Parameter("mu", st.norm(MU_LOC, MU_SCALE)),
        Parameter("sigma", st.halfnorm(0.0, SIGMA_SCALE)),
    ]


def _datasets(data: Sequence[SyntheticSpectrum], backend: str) -> list[Dataset]:
    return [
        Dataset(
            spectrum.container(),
            likelihood=build_likelihood("flexible", backend=backend),
            model=model_label(index),
            label=dataset_label(index),
        )
        for index, spectrum in enumerate(data)
    ]


def _models(data: Sequence[SyntheticSpectrum], backend: str) -> dict[str, Any]:
    return {
        model_label(index): model_for(backend, spectrum.wavelength)
        for index, spectrum in enumerate(data)
    }


def build_arm(
    data: Sequence[SyntheticSpectrum],
    arm: str,
    *,
    backend: str = "reference",
    seed: int = SEED,
) -> FittingProblem:
    """One of :data:`ARMS` over *data*, composed entirely from *backend*'s parts.

    Examples
    --------
    >>> data = generate_population(size=50)
    >>> [build_arm(data, arm).free_size for arm in ARMS]
    [48, 41, 50]
    """
    if arm not in ARMS:
        raise ValueError(f"unknown arm {arm!r}; the three are {list(ARMS)}.")
    models = _models(data, backend)
    datasets = _datasets(data, backend)
    if arm == "independent":
        return FittingProblem(models, datasets, seed=seed)
    if arm == "tied":
        sites = [f"{dataset.label}.likelihood.amplitude" for dataset in datasets]
        return FittingProblem(models, datasets, ties=[Tie(TIED_NAME, sites)], seed=seed)
    collection = DatasetCollection.plate(
        POPULATION_LABEL,
        datasets,
        within="likelihood",
        members=members(),
        hyperpriors=hyperpriors(),
    )
    return FittingProblem(models, collection, seed=seed)


# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------


def _flat(values: Any) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array.reshape(array.shape[0] * array.shape[1], *array.shape[2:])


def amplitude_draws(stored: Any, arm: str, index: int) -> np.ndarray:
    """Spectrum *index*'s GP amplitude, one value per posterior draw, Jy.

    The three arms store it under three names: ``d{i}.likelihood.amplitude``
    (independent), the tie's ``amplitude`` (tied: every spectrum's is the same
    draws) and element *i* of the derived ``gp.amplitude`` (population).
    """
    posterior = stored["posterior"]
    if arm == "independent":
        return _flat(posterior[f"{dataset_label(index)}.likelihood.amplitude"].values)
    if arm == "tied":
        return _flat(posterior[TIED_NAME].values)
    return _flat(posterior[f"{POPULATION_LABEL}.amplitude"].values)[:, index]


def _summary_of(name: str, draws: np.ndarray, truth: float | None = None) -> Summary:
    low, median, high = (float(v) for v in np.percentile(draws, (16.0, 50.0, 84.0)))
    return Summary(name=name, median=median, low=low, high=high, truth=truth, n_draws=draws.size)


def physical_summaries(stored: Any, index: int) -> dict[str, Summary]:
    """Spectrum *index*'s four model parameters, each with its truth attached.

    :func:`~examples.m2_misspecification.study.summarise` attaches the truth
    under the study's single-model names (``model.A``); here every spectrum has
    its own model, ``obj{i}``, so the truth is attached by local name.
    """
    names = [f"{model_label(index)}.{name}" for name in PARAMETER_NAMES]
    found = summarise(stored, names=names)
    return {
        name: dataclasses.replace(found[qualified], truth=TRUTH[name])
        for name, qualified in zip(PARAMETER_NAMES, names, strict=True)
    }


@dataclasses.dataclass(frozen=True)
class Localisation:
    """Where one spectrum's GP was needed, by :func:`ampere.results.gp_localisation_score`."""

    dataset: str
    peak: float
    largest: float


def localise(fit: Fit, *, thin: int = 20) -> dict[str, Localisation]:
    """GP localisation on every spectrum of *fit*, the derived group written back.

    The study's :func:`~examples.m2_misspecification.study.diagnose` answers
    for a one-dataset problem; this is its loop over the collection.
    """
    labels = list(fit.problem.datasets)
    fit.run = gp_localisation(fit.run, fit.problem, thin=thin)
    found: dict[str, Localisation] = {}
    for label in labels:
        score = gp_localisation_score(fit.run, dataset=label)
        values = np.asarray(score.values, dtype=float)
        coordinates = np.asarray(score.coordinates, dtype=float).reshape(values.size, -1)
        best = int(np.argmax(values))
        found[label] = Localisation(label, float(coordinates[best, 0]), float(values[best]))
    return found


@dataclasses.dataclass
class Fit:
    """One arm's problem and its stored run."""

    arm: str
    problem: FittingProblem
    run: Any
    seconds: float


@dataclasses.dataclass
class Comparison:
    """The three arms over one dataset collection, reduced to what is asserted on.

    Attributes
    ----------
    data
        The spectra, in plate order.
    fits
        ``{arm: Fit}``.
    backend, budget, seed
        What the fits were run on, at and with.
    """

    data: tuple[SyntheticSpectrum, ...]
    fits: dict[str, Fit]
    backend: str
    budget: EmceeBudget | NutsBudget | None
    seed: int

    @property
    def deviated(self) -> tuple[int, ...]:
        """The indices whose spectrum carries the deviation."""
        return tuple(i for i, d in enumerate(self.data) if d.scenario.amplitude != 0.0)

    @property
    def clean(self) -> tuple[int, ...]:
        """The indices whose spectrum carries none."""
        return tuple(i for i, d in enumerate(self.data) if d.scenario.amplitude == 0.0)

    def physical(self, arm: str, index: int) -> dict[str, Summary]:
        """Spectrum *index*'s model parameters under *arm*."""
        return physical_summaries(self.fits[arm].run, index)

    def worst_bias(self, arm: str, index: int) -> float:
        """The largest ``bias_in_widths`` of spectrum *index*'s four model parameters."""
        return max(float(s.bias_in_widths or 0.0) for s in self.physical(arm, index).values())

    def amplitude(self, arm: str, index: int) -> Summary:
        """Spectrum *index*'s GP amplitude under *arm*, Jy."""
        draws = amplitude_draws(self.fits[arm].run, arm, index)
        return _summary_of(f"{arm}:{dataset_label(index)}.amplitude", draws)

    def hyperparameter(self, name: str) -> np.ndarray:
        """``gp.mu`` or ``gp.sigma``'s draws, from the population arm."""
        stored = self.fits["population"].run
        return _flat(stored["posterior"][f"{POPULATION_LABEL}.{name}"].values)

    def shrinkage(self, index: int = POORLY_CONSTRAINED) -> float:
        """How much narrower member *index*'s amplitude is under the population.

        The ratio of the 68 % half-widths, independent over population: above
        one means the population constrained the member better than its own
        data could.
        """
        return (
            self.amplitude("independent", index).width / self.amplitude("population", index).width
        )

    def table(self) -> str:
        """The comparison as plain text."""
        arms = [arm for arm in ARMS if arm in self.fits]
        header = "".join(f"{arm:>14s}" for arm in arms)
        rows = [
            (
                f"M2 over {len(self.data)} spectra: {self.backend}, {self.budget}, "
                f"run seed {self.seed}"
            ),
            "",
            "worst |median - truth| of the four model parameters, posterior widths",
            f"{'spectrum':<26s}{header}",
            "-" * (26 + 14 * len(arms)),
        ]
        for index, spectrum in enumerate(self.data):
            kind = spectrum.scenario.key
            if index == POORLY_CONSTRAINED and kind == "none":
                kind = f"none, {POOR_NOISE_FRACTION:.0%} noise"
            cells = "".join(f"{self.worst_bias(arm, index):14.2f}" for arm in arms)
            rows.append(f"{dataset_label(index) + ' (' + kind + ')':<26s}{cells}")
        rows += ["", "GP amplitude median (68 % half-width), Jy", "-" * (26 + 14 * len(arms))]
        for index in range(len(self.data)):
            cells = "".join(
                f"{self.amplitude(arm, index).median:8.4f}({self.amplitude(arm, index).width:.3f})"
                for arm in arms
            )
            rows.append(f"{dataset_label(index):<26s}{cells}")
        if "population" in self.fits:
            mu, sigma = self.hyperparameter("mu"), self.hyperparameter("sigma")
            low, high = np.quantile(sigma, [0.025, 0.975])
            rows += [
                "",
                f"gp.mu    median {np.median(mu):.3f}",
                f"gp.sigma median {np.median(sigma):.3f}, central 95 % [{low:.3f}, {high:.3f}]",
            ]
            if "independent" in self.fits:
                rows += [
                    f"sample spread of log(independent medians): {sample_spread(self):.3f}",
                    (
                        f"shrinkage of {dataset_label(POORLY_CONSTRAINED)}'s amplitude "
                        f"(independent / population width): {self.shrinkage():.2f}"
                    ),
                ]
        rows.append("")
        rows.append(
            "sampled in "
            + ", ".join(f"{arm} {fit.seconds:.0f} s" for arm, fit in self.fits.items())
        )
        return "\n".join(rows)


def sample_spread(comparison: Comparison) -> float:
    """The spread the population is asked to recover: ``std(log a_i)``.

    The truth has no GP amplitude — the deviation is a 12 % line, not a draw
    from a kernel — so the generated amplitudes are *defined* through the
    status-quo fit: ``a_i`` is spectrum *i*'s amplitude posterior median under
    the independent arm, the mapping from "this spectrum carries this
    deviation" to "this is the amplitude a Matérn-3/2 needs to absorb it".
    The deviated spectra map to about 0.02 Jy, the clean ones at 1 % noise to
    about 0.001 Jy and the clean one at 5 % to about 0.01 Jy; the sample
    standard deviation of their logs (``ddof=0``, the population's own
    ``sigma``) is the number.
    """
    medians = np.array(
        [comparison.amplitude("independent", i).median for i in range(len(comparison.data))]
    )
    return float(np.std(np.log(medians)))


def compare(
    *,
    size: int = 200,
    backend: str = "reference",
    budget: EmceeBudget | NutsBudget | None = None,
    seed: int = SEED,
    data_seed: int = DATA_SEED,
    arms: Sequence[str] = ARMS,
    progress: bool = False,
) -> Comparison:
    """Fit every arm in *arms* to one generated collection and keep the runs.

    The arms share one :func:`generate_population` output — the same numbers,
    not the same recipe — so a difference between their posteriors is a
    difference between their declarations. :func:`~examples.m2_misspecification.study.run`
    dispatches on the backend: emcee on the reference path (whose budget must
    have more than twice the widest arm's fifty dimensions in walkers), NUTS on
    torch and jax.
    """
    if backend not in BACKENDS:
        raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")
    data = generate_population(size=size, seed=data_seed)
    fits: dict[str, Fit] = {}
    for arm in arms:
        problem = build_arm(data, arm, backend=backend, seed=seed)
        started = time.perf_counter()
        stored = run(problem, budget, progress=progress)
        fits[arm] = Fit(arm, problem, stored, time.perf_counter() - started)
    return Comparison(data=data, fits=fits, backend=backend, budget=budget, seed=seed)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the comparison and print its table."""
    parser = argparse.ArgumentParser(prog="python -m examples.m2_populations")
    parser.add_argument("--backend", default="reference", choices=list(BACKENDS))
    parser.add_argument("--size", type=int, default=200, help="points per spectrum")
    parser.add_argument("--quick", action="store_true", help="tests/m2's per-PR budget")
    parser.add_argument("--seed", type=int, default=SEED, help="the run seed")
    parser.add_argument("--data-seed", type=int, default=DATA_SEED, help="spectrum 0's data seed")
    parser.add_argument("--arms", nargs="+", default=list(ARMS), choices=list(ARMS))
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--figures", default=None, help="write the figure into this directory")
    arguments = parser.parse_args(argv)
    budget = budget_for(arguments.backend, quick=arguments.quick)
    comparison = compare(
        size=arguments.size,
        backend=arguments.backend,
        budget=budget,
        seed=arguments.seed,
        data_seed=arguments.data_seed,
        arms=arguments.arms,
        progress=arguments.progress,
    )
    print(comparison.table())
    if "population" in comparison.fits:
        print("\nGP localisation peak (score), population arm:")
        for label, found in localise(comparison.fits["population"]).items():
            print(f"  {label}: {found.peak:.5f} um ({found.largest:.2f})")
        centre = scenario_named(DEVIATION).localised_at
        print(f"  the injected line is at {centre} um")
    if arguments.figures is not None:
        from . import figures

        written = figures.save_figure(comparison.data, arguments.figures, comparison)
        print(f"\nWrote {len(written)} figure to {arguments.figures}")
    return 0
