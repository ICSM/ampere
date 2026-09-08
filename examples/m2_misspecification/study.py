"""The study: two likelihoods, four scenarios, three backends, three sizes.

This module is the experiment. Everything above it (:mod:`.generators`,
:mod:`.model`) supplies data and physics; everything below it (:mod:`.figures`,
``tests/m2``, ``tests/benchmarks/test_m2_*.py``) consumes what this module
returns. There is one code path for a study run, whether it is a three-minute
CI assertion, a ten-minute documentation run or the full ladder.

The two likelihoods
-------------------
**Standard** is :class:`~ampere.core.IndependentNoise` under a
:class:`~ampere.core.GaussianFamily`: the ordinary chi-square likelihood every
spectral fit starts from, and nothing else. This is a deliberate departure from
the paper study, which produced its "standard" case by pinning a GP's
hyperparameter priors to ``1e-10``. Pinning is not the same object: it keeps
two extra dimensions in the sampler, keeps the GP solve in the hot loop, and
its "zero" amplitude is a prior whose mass sits at :math:`10^{-10}` rather than
at nothing at all. The honest v2 statement of "no flexible likelihood" is a
likelihood with no GP in it, and it is also the one a reader would actually
write.

**Flexible** is :class:`~ampere.core.GaussianProcessNoise` over the residuals.
The kernel is **Matérn-3/2**, not the paper study's squared exponential. That
is W2.10's binding decision and it has two reasons, both structural: a
once-differentiable sample path represents real model deficiencies better than
an analytic one, and Matérn-3/2 is *exactly* quasiseparable, which is what
makes :class:`~ampere.core.QuasisepGP` an O(N) solve rather than an
approximation — and therefore what makes the 20 000-point rung of the ladder
possible at all. The squared exponential is kept as
``kernel="squared_exponential"``, dense-only, at 200 points, as the cross-check
against what legacy ampere actually did.

The hyperparameter priors are the paper study's numbers, expressed as
half-normals: amplitude ``halfnorm(scale=0.3)`` in Jy, length scale
``halfnorm(scale=0.003)`` in micron. Half-normal rather than uniform because
the study's claim depends on the GP being *able* to switch itself off: a prior
whose mass piles up at zero is what stops the ``none`` scenario's posterior
being inflated by a GP with nothing to do, and that control is what makes the
comparison in the other three scenarios mean anything.

Engines
-------
``reference`` samples with emcee, because the numpy path has no gradient;
``torch`` and ``jax`` sample with NUTS through
:class:`~ampere.inference.NUTSEngine`, which reaches the backend's realised
density through :func:`ampere.core.realise`. :func:`run` dispatches on the
problem's own ``backend`` flag rather than on an argument, so a caller cannot
ask for a sampler the problem cannot feed.

Reading the result
------------------
:func:`summarise` reduces a stored run to one :class:`Summary` per parameter,
and the number the whole milestone turns on is
:attr:`Summary.bias_in_widths` — ``|median - truth|`` divided by the posterior's
own 68 % half-width. It is the right statistic because it is the one that
distinguishes the two failure modes a misspecified fit can have. A likelihood
that is *wrong* puts the truth many posterior widths from its median: it is
confident and mistaken. A likelihood that has *absorbed* the misspecification
may still have a displaced median, but its posterior is wide enough to contain
the truth, and its error bar is then a statement a reader can act on. The
flexible likelihood does not magically recover the right point estimate — it
recovers honest uncertainty, and this module measures exactly that.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping, MutableMapping, Sequence
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Kernel,
    Likelihood,
)
from ampere.results import (
    add_residuals,
    gp_localisation,
    gp_localisation_score,
    residual_whiteness,
)

from .generators import SCENARIOS, Scenario, SyntheticSpectrum, generate, scenario_named
from .model import PARAMETER_NAMES, TRUTH

__all__ = [
    "AGREEMENT_EMCEE",
    "AGREEMENT_NUTS",
    "BACKENDS",
    "CI_INTERVAL_TOLERANCE",
    "CI_MEDIAN_TOLERANCE",
    "DATASET_LABEL",
    "DOC_EMCEE",
    "DOC_NUTS",
    "FLEXIBLE_MAX_BIAS_WIDTHS",
    "GP_AMPLITUDE_SCALE",
    "GP_LENGTH_SCALE",
    "KERNELS",
    "LIKELIHOODS",
    "LOCALISATION_CONTRAST",
    "LOCALISATION_TOLERANCE_POINTS",
    "MEASURED_SCATTER",
    "MILESTONE_EMCEE",
    "MILESTONE_INTERVAL_TOLERANCE",
    "MILESTONE_MEDIAN_TOLERANCE",
    "MILESTONE_NUTS",
    "PAPER_EMCEE",
    "PHYSICAL_NAMES",
    "SEED",
    "SOLVERS",
    "STANDARD_MIN_BIAS_WIDTHS",
    "TEST_EMCEE",
    "TEST_NUTS",
    "WHITENESS_STRUCTURE_LEVEL",
    "WHITENESS_WHITE_LEVEL",
    "Diagnosis",
    "EmceeBudget",
    "NutsBudget",
    "Summary",
    "agreement",
    "build_kernel",
    "build_likelihood",
    "build_problem",
    "diagnose",
    "model_for",
    "prepare",
    "run",
    "run_study",
    "summarise",
]

#: The label the single dataset carries, and therefore the prefix its
#: likelihood's hyperparameters appear under in the posterior.
DATASET_LABEL = "default"

#: The merged posterior names of the four physical parameters.
PHYSICAL_NAMES: tuple[str, ...] = tuple(f"model.{name}" for name in PARAMETER_NAMES)

#: The run seed. Distinct from the *data* seed (``generators.SEED``, 42) on
#: purpose: one identifies the spectrum, the other the chain, and conflating
#: them would make "the same data, a different chain" impossible to ask for.
SEED = 20260910

#: The two likelihoods under comparison.
LIKELIHOODS: tuple[str, ...] = ("standard", "flexible")
#: The backends the study runs on. ``reference`` is always available; the other
#: two need their extra.
BACKENDS: tuple[str, ...] = ("reference", "torch", "jax")
#: The GP solvers. ``quasisep`` is O(N) and Matérn-only; ``dense`` is O(N^3).
SOLVERS: tuple[str, ...] = ("quasisep", "dense")
#: The kernels. ``matern32`` is the reproduction's; ``squared_exponential`` is
#: legacy's, kept for the 200-point cross-check and refused by ``quasisep``.
KERNELS: tuple[str, ...] = ("matern32", "squared_exponential")

#: The flexible likelihood's hyperparameter prior scales — the paper study's
#: ``covWeightPrior`` and ``scaleLengthPrior``, as half-normal scales.
GP_AMPLITUDE_SCALE = 0.3
GP_LENGTH_SCALE = 0.003


# ---------------------------------------------------------------------------
# Budgets
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class EmceeBudget:
    """Walkers, steps and burn-in for an ensemble run."""

    walkers: int
    steps: int
    burn_in: int

    @property
    def kept(self) -> int:
        """Draws that reach the stored run."""
        return self.walkers * (self.steps - self.burn_in)


@dataclasses.dataclass(frozen=True)
class NutsBudget:
    """Draws, warm-up and chains for a NUTS run."""

    draws: int
    warmup: int
    chains: int

    @property
    def kept(self) -> int:
        """Draws that reach the stored run."""
        return self.draws * self.chains


#: The paper study's own budget: 40 walkers, 5 000 steps, 2 000 discarded.
PAPER_EMCEE = EmceeBudget(walkers=40, steps=5_000, burn_in=2_000)
#: What ``tests/m2`` uses: short, and stated as such. Every assertion made at
#: this budget has a margin large compared with the Monte Carlo error it
#: leaves; :data:`MEASURED_MEDIAN_SCATTER` records the measurement.
TEST_EMCEE = EmceeBudget(walkers=20, steps=900, burn_in=450)
TEST_NUTS = NutsBudget(draws=400, warmup=400, chains=2)
#: What ``tests/m2/test_backend_agreement.py`` uses. Longer than
#: :data:`TEST_EMCEE` / :data:`TEST_NUTS` on purpose: a cross-backend
#: comparison at a budget whose Monte Carlo error is comparable with the
#: tolerance is a test of the samplers' luck rather than of the backends, and
#: the flexible likelihood's ``B`` — partly degenerate with the GP — is the
#: parameter that sets the bar. Roughly three minutes for the four runs.
AGREEMENT_EMCEE = EmceeBudget(walkers=24, steps=2_400, burn_in=1_200)
AGREEMENT_NUTS = NutsBudget(draws=500, warmup=500, chains=2)
#: The **milestone** budget: what the numbers recorded in
#: ``docs/source/m2_misspecification.rst`` and in W2.10's report were measured
#: at, and what ``pytest -m m2_full`` uses. Long enough that the cross-backend
#: comparison measures the two posteriors rather than the two samplers.
MILESTONE_EMCEE = EmceeBudget(walkers=32, steps=4_000, burn_in=2_000)
MILESTONE_NUTS = NutsBudget(draws=1_000, warmup=1_000, chains=4)
#: What the CLI and the documentation page run at.
DOC_EMCEE = MILESTONE_EMCEE
DOC_NUTS = MILESTONE_NUTS


# ---------------------------------------------------------------------------
# The thresholds every assertion in ``tests/m2`` is written against. They live
# here, once, so the tests, the documentation page and W2.10's report cite one
# set of numbers rather than three transcriptions of it.
# ---------------------------------------------------------------------------

#: The flexible likelihood's parameter recovery must stay inside this many
#: posterior widths of the truth, in **every** scenario. Measured worst case at
#: :data:`TEST_EMCEE`, over all four scenarios and all four parameters: 0.98.
FLEXIBLE_MAX_BIAS_WIDTHS = 1.5

#: Under misspecification the standard likelihood must miss the truth by at
#: least this many posterior widths, on at least one parameter. Measured worst
#: case (i.e. the *smallest* such excess, in the mild scenario) at
#: :data:`TEST_EMCEE`: 3.16; the strong scenarios reach 7.2 and 4.8.
STANDARD_MIN_BIAS_WIDTHS = 2.0

#: Residual whiteness (``ampere.results.residual_whiteness``, 199 permutations)
#: on a **standard**-likelihood fit. A p-value at or below the first level says
#: the residuals carry structure the model did not capture; at or above the
#: second, they are consistent with white. Measured: 0.005 (the 199-permutation
#: floor) in all three misspecified scenarios, 0.195 in the control.
WHITENESS_STRUCTURE_LEVEL = 0.01
WHITENESS_WHITE_LEVEL = 0.05

#: GP localisation must peak within this many grid spacings of an injected
#: *localised* deviation. Measured at 200 points: 0.33 spacings.
LOCALISATION_TOLERANCE_POINTS = 3.0

#: ... and its peak score must exceed the control scenario's by this factor,
#: so that "it peaked somewhere" is not mistaken for "it found something".
#: Measured: 36.6x.
LOCALISATION_CONTRAST = 5.0

#: Cross-backend agreement, in units of the posterior's own 68 % half-width.
#: The first pair is what ``tests/m2/test_backend_agreement.py`` asserts at
#: :data:`AGREEMENT_EMCEE` / :data:`AGREEMENT_NUTS`; the second is what
#: ``pytest -m m2_full`` asserts at :data:`MILESTONE_EMCEE` /
#: :data:`MILESTONE_NUTS`, and what the milestone runs recorded in
#: ``docs/source/m2_misspecification.rst`` actually achieved (worst observed:
#: 0.072 on a median, 0.109 on an interval endpoint).
#:
#: Both pairs are set at three to five times the Monte Carlo error of the
#: *difference* of two runs at their budget, derived from
#: :data:`MEASURED_SCATTER` below. That is what makes them assertions about the
#: backends rather than about the samplers' luck: a tolerance comparable with
#: the Monte Carlo error would fail a few per cent of the time with nothing
#: wrong, and a test that does that gets deleted rather than investigated.
CI_MEDIAN_TOLERANCE = 0.5
CI_INTERVAL_TOLERANCE = 0.75
MILESTONE_MEDIAN_TOLERANCE = 0.20
MILESTONE_INTERVAL_TOLERANCE = 0.30

#: The measured run-to-run scatter, in posterior widths: the **range** across
#: three run seeds on identical data, over the four physical parameters, for
#: the flexible likelihood on ``strong_smooth`` at 200 points. That combination
#: is chosen because it is the worst-conditioned one in the study — the GP and
#: the continuum slope ``B`` are partly degenerate, so ``B``'s posterior is
#: broad and mixes slowly, and every tolerance above is really set by it.
#:
#: Read as ``(median, interval endpoint)``. A range across three draws is about
#: 1.7 standard deviations, and the difference of two independent runs has
#: :math:`\\sqrt{2}` times one run's standard deviation, which is the
#: arithmetic behind the tolerances.
MEASURED_SCATTER: dict[str, tuple[float, float]] = {
    "emcee at TEST_EMCEE": (0.31, 0.48),
    "emcee at AGREEMENT_EMCEE": (0.14, 0.21),
    "emcee at MILESTONE_EMCEE": (0.08, 0.11),
    "NUTS at TEST_NUTS": (0.18, 0.39),
    "NUTS at AGREEMENT_NUTS": (0.07, 0.24),
    "NUTS at MILESTONE_NUTS": (0.05, 0.10),
}


# ---------------------------------------------------------------------------
# Composition
# ---------------------------------------------------------------------------


def _backend_module(backend: str) -> Any:
    """This backend's namespace, imported lazily and refused by name."""
    if backend == "reference":
        import ampere.core as module

        return module
    if backend == "torch":
        from ampere.backends import torch as module  # type: ignore[assignment]

        return module
    if backend == "jax":
        from ampere.backends import jax as module  # type: ignore[assignment]

        module.configure_x64()
        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")


def model_for(backend: str, wavelength: Any, **given: Any) -> Any:
    """The toy model on *backend*, with the shared declaration.

    Imports the variant's module lazily, so ``import
    examples.m2_misspecification`` needs neither torch nor jax.
    """
    if backend == "reference":
        from .model import AbsorptionLines

        return AbsorptionLines(wavelength, **given)
    if backend == "torch":
        from .model_torch import AbsorptionLines as TorchAbsorptionLines

        return TorchAbsorptionLines(wavelength, **given)
    if backend == "jax":
        _backend_module("jax")  # configure_x64 before the model refuses
        from .model_jax import AbsorptionLines as JaxAbsorptionLines

        return JaxAbsorptionLines(wavelength, **given)
    raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")


#: The class name each choice resolves to. Every backend's namespace —
#: ``ampere.core`` included, which *is* the reference backend's toolkit — spells
#: these identically (W2.12: "one name per backend, everywhere"), so the study
#: reaches for a name and never for a branch on the backend.
_KERNEL_CLASSES = {"matern32": "Matern32", "squared_exponential": "SquaredExponential"}
_SOLVER_CLASSES = {"quasisep": "QuasisepGP", "dense": "DenseGP"}


def build_kernel(backend: str = "reference", kernel: str = "matern32") -> Kernel:
    """The named kernel, on *backend*, with the study's hyperparameter priors.

    The amplitude carries Jy and the length scale micron because
    :class:`~ampere.core.GaussianProcessNoise` checks them against the observed
    container: a length scale silently in the wrong unit is the sort of error
    that produces a plausible number rather than a crash.
    """
    if kernel not in _KERNEL_CLASSES:
        raise ValueError(f"unknown kernel {kernel!r}; the two are {list(KERNELS)}.")
    cls = getattr(_backend_module(backend), _KERNEL_CLASSES[kernel])
    return cls(
        st.halfnorm(scale=GP_AMPLITUDE_SCALE),
        st.halfnorm(scale=GP_LENGTH_SCALE),
        amplitude_unit=u.Jy,
        length_scale_unit=u.micron,
    )


def build_likelihood(
    kind: str = "flexible",
    *,
    backend: str = "reference",
    kernel: str = "matern32",
    solver: str = "quasisep",
) -> Likelihood:
    """The standard or the flexible likelihood, composed from *backend*'s parts.

    Every part of a native problem must come from one backend — W2.13's
    capability widening made a core :class:`~ampere.core.IndependentNoise`
    inside a torch problem a declared disagreement rather than a silent numpy
    island — so this function is the one place that knows which namespace to
    reach into.
    """
    module = _backend_module(backend)
    family = GaussianFamily()
    if kind == "standard":
        return Likelihood(family, module.IndependentNoise())
    if kind != "flexible":
        raise ValueError(f"unknown likelihood {kind!r}; the two are {list(LIKELIHOODS)}.")
    if solver not in _SOLVER_CLASSES:
        raise ValueError(f"unknown solver {solver!r}; the two are {list(SOLVERS)}.")
    if kernel == "squared_exponential" and solver == "quasisep":
        raise ValueError(
            "the squared exponential is not quasiseparable, so QuasisepGP refuses it. That "
            "asymmetry is the concrete reason Matern-3/2 is v2's default; run the legacy "
            "cross-check with solver='dense'."
        )
    solver_cls = getattr(module, _SOLVER_CLASSES[solver])
    return Likelihood(
        family, module.GaussianProcessNoise(build_kernel(backend, kernel), solver_cls())
    )


def build_problem(
    data: SyntheticSpectrum | Scenario | str,
    *,
    size: int = 200,
    backend: str = "reference",
    likelihood: str = "flexible",
    kernel: str = "matern32",
    solver: str = "quasisep",
    seed: int = SEED,
    strict: bool = False,
) -> FittingProblem:
    """One fitting problem: this scenario, this likelihood, on this backend.

    *data* may be a ready-made :class:`~examples.m2_misspecification.generators.SyntheticSpectrum`
    (so several problems can share one spectrum, which is what makes the
    standard/flexible comparison a comparison) or a scenario to generate one
    from at *size* points.

    Examples
    --------
    >>> problem = build_problem("strong_smooth", size=200, likelihood="standard")
    >>> problem.free_size
    4
    >>> problem.backend
    'reference'
    >>> flexible = build_problem("strong_smooth", size=200, likelihood="flexible")
    >>> flexible.free_size
    6
    """
    spectrum = data if isinstance(data, SyntheticSpectrum) else generate(data, size=size)
    composed = build_likelihood(likelihood, backend=backend, kernel=kernel, solver=solver)
    dataset = Dataset(spectrum.container(), likelihood=composed, label=DATASET_LABEL)
    return FittingProblem(
        model_for(backend, spectrum.wavelength),
        [dataset],
        seed=seed,
        strict=strict,
    )


# ---------------------------------------------------------------------------
# Running
# ---------------------------------------------------------------------------


def run(
    problem: FittingProblem,
    budget: EmceeBudget | NutsBudget | None = None,
    *,
    progress: bool = False,
) -> Any:
    """Sample *problem* with the engine its backend can feed, and return the run.

    Dispatches on ``problem.backend``, never on an argument: the reference path
    has no gradient, so emcee is not a preference there but the only option,
    and the two differentiable backends get NUTS for the same structural
    reason.
    """
    from ampere.inference import EmceeEngine, NUTSEngine

    if problem.backend == "reference":
        chosen = budget if isinstance(budget, EmceeBudget) else TEST_EMCEE
        engine = EmceeEngine(problem, walkers=chosen.walkers)
        return engine.run(chosen.steps, burn_in=chosen.burn_in, progress=progress)
    nuts = budget if isinstance(budget, NutsBudget) else TEST_NUTS
    return NUTSEngine(problem).run(
        nuts.draws, warmup=nuts.warmup, chains=nuts.chains, progress=progress
    )


def run_study(
    *,
    size: int = 200,
    backend: str = "reference",
    scenarios: Sequence[Scenario | str] | None = None,
    likelihoods: Sequence[str] = LIKELIHOODS,
    budget: EmceeBudget | NutsBudget | None = None,
    kernel: str = "matern32",
    solver: str = "quasisep",
    seed: int = SEED,
    progress: bool = False,
) -> dict[tuple[str, str], dict[str, Any]]:
    """Run every ``(scenario, likelihood)`` pair and keep what a figure needs.

    Returns a mapping from ``(scenario_key, likelihood)`` to a dictionary
    holding the ``run`` (the ArviZ ``DataTree``), the ``problem`` it was over
    (needed to derive residuals and GP localisation, which are not stored by
    default) and the ``data`` it was fitted to.

    The two likelihoods in a scenario share one
    :class:`~examples.m2_misspecification.generators.SyntheticSpectrum` — the
    same numbers, not merely the same recipe — because otherwise the difference
    between their posteriors would include a difference in their data.
    """
    chosen = tuple(SCENARIOS if scenarios is None else scenarios)
    results: dict[tuple[str, str], dict[str, Any]] = {}
    for entry in chosen:
        scenario = scenario_named(entry) if isinstance(entry, str) else entry
        data = generate(scenario, size=size)
        for kind in likelihoods:
            problem = build_problem(
                data,
                backend=backend,
                likelihood=kind,
                kernel=kernel,
                solver=solver,
                seed=seed,
            )
            results[(scenario.key, kind)] = {
                "run": run(problem, budget, progress=progress),
                "problem": problem,
                "data": data,
                "scenario": scenario,
                "likelihood": kind,
            }
    return results


# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Summary:
    """One parameter's posterior, reduced to the numbers the study asserts on.

    Attributes
    ----------
    name
        The merged posterior name (``model.d1``,
        ``default.likelihood.amplitude``).
    median, low, high
        The posterior median and the 16th and 84th percentiles.
    truth
        The generating value, where there is one. The GP hyperparameters have
        none — there is no true amplitude for a deficiency the model does not
        contain — so theirs is ``None``.
    n_draws
        How many draws the summary was computed from.
    """

    name: str
    median: float
    low: float
    high: float
    truth: float | None
    n_draws: int

    @property
    def width(self) -> float:
        """The 68 % half-width: ``(q84 - q16) / 2``, the posterior's own scale."""
        return 0.5 * (self.high - self.low)

    @property
    def bias(self) -> float | None:
        """``median - truth``, in the parameter's own units."""
        return None if self.truth is None else self.median - self.truth

    @property
    def bias_in_widths(self) -> float | None:
        """``|median - truth|`` in units of :attr:`width` — the study's statistic.

        A degenerate (zero-width) posterior returns ``inf`` rather than raising:
        a fit that is infinitely confident and wrong is exactly the failure this
        number exists to report, and it should not be a ``ZeroDivisionError``.
        """
        if self.truth is None:
            return None
        if self.width <= 0.0:
            return float("inf")
        return abs(self.median - self.truth) / self.width

    @property
    def covers_truth(self) -> bool | None:
        """Whether the 68 % interval contains the truth."""
        if self.truth is None:
            return None
        return bool(self.low <= self.truth <= self.high)

    def __repr__(self) -> str:
        offset = self.bias_in_widths
        bias = "" if offset is None else f", {offset:.2f} widths from truth"
        return (
            f"<Summary {self.name!r}: {self.median:.5g} "
            f"+{self.high - self.median:.3g} -{self.median - self.low:.3g}{bias}>"
        )


def summarise(run: Any, *, names: Iterable[str] | None = None) -> dict[str, Summary]:
    """Reduce a stored run's posterior to one :class:`Summary` per variable.

    Reads the ``posterior`` group and nothing else, so it works on a run
    restored from netCDF as readily as on one just emitted.
    """
    posterior = run["posterior"]
    available = [str(name) for name in posterior.dataset.data_vars]
    wanted = available if names is None else [str(name) for name in names]
    truths = {f"model.{name}": value for name, value in TRUTH.items()}
    summaries: dict[str, Summary] = {}
    for name in wanted:
        if name not in available:
            raise KeyError(f"the run's posterior has no variable {name!r}; it has {available}.")
        draws = np.asarray(posterior[name].values, dtype=float).ravel()
        low, median, high = (float(v) for v in np.percentile(draws, (16.0, 50.0, 84.0)))
        summaries[name] = Summary(
            name=name,
            median=median,
            low=low,
            high=high,
            truth=truths.get(name),
            n_draws=int(draws.size),
        )
    return summaries


@dataclasses.dataclass(frozen=True)
class Diagnosis:
    """What the two post-fit diagnostics said about one run.

    Attributes
    ----------
    whiteness_statistic, whiteness_p_value
        ``ampere.results.residual_whiteness``'s pooled answer. **A small
        p-value means the residuals are not white**, so a misspecified fit with
        a standard likelihood should produce one. ``diagnostics.md`` §3.1
        scopes this family to standard-likelihood fits, and :func:`diagnose`
        computes it for those only: a GP fit's residuals still contain
        everything the GP explained, so testing *them* would answer a different
        question than the one asked.
    localisation_peak
        The coordinate at which the GP-localisation anomaly score is largest,
        for a flexible fit; ``None`` for a standard one, which has no GP.
    localisation_max
        The score there, in units of the posterior standard deviation of the
        conditioned GP mean ("how many sigma of GP the fit needed here").
    """

    dataset: str
    whiteness_statistic: float | None
    whiteness_p_value: float | None
    localisation_peak: float | None
    localisation_max: float | None


def diagnose(
    entry: MutableMapping[str, Any],
    *,
    thin: int = 1,
    n_permutations: int = 199,
    seed: int = 20260910,
) -> Diagnosis:
    """Run the applicable post-fit diagnostic on one :func:`run_study` entry.

    The two families answer different questions and are scoped accordingly
    (``diagnostics.md`` §3.1 and §4): residual whiteness asks "did a
    *standard*-likelihood fit leave structure behind?", GP localisation asks
    "*where* did a flexible fit need its GP?". Neither is computed where it
    would be circular — a GP fit's residuals still contain everything the GP
    explained, so testing those for whiteness would answer a different question
    from the one asked.

    Both derive the group they read (``residuals``, ``gp_localisation``) from
    the run and the problem, because ``results.md`` §7 keeps ``N_draws x
    N_obs`` groups out of the default emission; *thin* is how a large run stays
    affordable. **The derived group is written back into** ``entry["run"]``, so
    the figures can read it without deriving it a second time.
    """
    problem = entry["problem"]
    label = next(iter(problem.datasets))
    statistic: float | None = None
    p_value: float | None = None
    peak: float | None = None
    largest: float | None = None
    if entry["likelihood"] == "standard":
        entry["run"] = add_residuals(entry["run"], problem, thin=thin)
        test = residual_whiteness(
            entry["run"], dataset=label, n_permutations=n_permutations, seed=seed
        )
        statistic, p_value = float(test.statistic), float(test.p_value)
    else:
        entry["run"] = gp_localisation(entry["run"], problem, thin=thin)
        score = gp_localisation_score(entry["run"], dataset=label)
        values = np.asarray(score.values, dtype=float)
        coordinates = np.asarray(score.coordinates, dtype=float).reshape(values.size, -1)
        index = int(np.argmax(values))
        peak = float(coordinates[index, 0])
        largest = float(values[index])
    return Diagnosis(
        dataset=str(label),
        whiteness_statistic=statistic,
        whiteness_p_value=p_value,
        localisation_peak=peak,
        localisation_max=largest,
    )


# ---------------------------------------------------------------------------
# Cross-backend agreement
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Agreement:
    """How far apart two backends put one parameter, in posterior widths.

    Every quantity is divided by the **mean of the two posteriors' own 68 %
    half-widths**, which is the only scale on which "do these agree?" has an
    answer: two samplers on the same posterior differ by their Monte Carlo
    error, and Monte Carlo error is measured in units of the posterior, never
    in the parameter's own units. A disagreement of 0.1 means the two medians
    sit a tenth of a standard-deviation-ish apart; a disagreement of 3 means
    they are sampling different distributions.
    """

    name: str
    width: float
    median: float
    low: float
    high: float

    @property
    def worst_interval(self) -> float:
        """The larger of the two 68 % endpoint disagreements."""
        return max(self.low, self.high)

    def __repr__(self) -> str:
        return (
            f"<Agreement {self.name!r}: median {self.median:.3f}, "
            f"interval {self.worst_interval:.3f} (posterior widths)>"
        )


def agreement(
    left: Mapping[str, Summary],
    right: Mapping[str, Summary],
    *,
    names: Iterable[str] = PHYSICAL_NAMES,
) -> dict[str, Agreement]:
    """Compare two backends' summaries of the same problem, parameter by parameter.

    *left* and *right* are :func:`summarise` outputs. Only *names* are compared,
    defaulting to the four physical parameters: the GP hyperparameters are
    compared too where a caller asks, but they are not what the milestone's
    claim is about, and one of them (the length scale) is the parameter whose
    posterior is most obviously non-Gaussian, so a quantile comparison on it
    says less than it appears to.
    """
    result: dict[str, Agreement] = {}
    for name in names:
        first, second = left[name], right[name]
        width = 0.5 * (first.width + second.width)
        if width <= 0.0:
            raise ValueError(
                f"parameter {name!r} has zero posterior width in both runs; there is no scale on "
                f"which to express a disagreement."
            )
        result[name] = Agreement(
            name=name,
            width=width,
            median=abs(first.median - second.median) / width,
            low=abs(first.low - second.low) / width,
            high=abs(first.high - second.high) / width,
        )
    return result


def prepare(
    results: MutableMapping[tuple[str, str], MutableMapping[str, Any]],
    *,
    thin: int = 20,
    n_permutations: int = 199,
) -> dict[tuple[str, str], Diagnosis]:
    """Run :func:`diagnose` over a whole :func:`run_study` result, in place.

    Every entry's ``run`` gains the derived group its likelihood's diagnostic
    needs — ``residuals`` for a standard fit, ``gp_localisation`` for a flexible
    one — which is also what :mod:`examples.m2_misspecification.figures` reads.
    Deriving them once, here, is what stops the figures and the numbers being
    computed from two different thinnings of the same chain.

    *thin* defaults to 20 because both groups are ``N_draws x N_obs``: at the
    milestone budget that is 128 000 x 200 floats undecimated, and the median
    across 6 400 draws is not measurably worse than the median across 128 000.
    """
    return {
        key: diagnose(entry, thin=thin, n_permutations=n_permutations)
        for key, entry in results.items()
    }
