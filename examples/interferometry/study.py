"""The study: three sky models, one flexible likelihood, three engines.

This module is the experiment, in the shape
:mod:`examples.m2_misspecification.study` established: everything above it
(:mod:`.generators`, :mod:`.model`) supplies data and physics, everything
below it (:mod:`.figures`, ``tests/interferometry``, ``tests/examples``)
consumes what this module returns.

The three arms on one dataset
------------------------------
One binary-plus-disc truth (:mod:`.generators`), observed once as
visibilities and closure phases, fitted three ways:

``"correct"``
    The disc is in the model, at its true value, and the noise is
    independent. The control: nothing is wrong, so the fit should simply be
    calibrated.
``"incomplete"``
    The disc is entirely absent from the model; the noise is still
    independent. What a naive fit of an under-modelled sky looks like:
    confident, and wrong, by the M2 measure (``|median - truth|`` in
    posterior widths).
``"flexible"``
    The disc is still absent, but the visibilities carry
    :class:`~ampere.core.GaussianProcessNoise` over ``(u, v)`` — W4.2's
    circular complex Gaussian process. The M2 question, asked of this
    modality: does the flexible likelihood keep the binary's two parameters
    calibrated when the sky model is wrong, by absorbing the disc's smooth,
    unmodelled correlation in the residual visibilities? The closure phases
    stay under independent von Mises noise in every arm — a GP on a
    *wrapped* observable is a latent composition the reference (emcee) path
    cannot run (``phase4_placement_memo.md`` §7.2's "blind alley"), and it is
    not what either W4.2 or this study needs: the flagship GP is on the
    visibilities.

Two ways of asking the same question
--------------------------------------
:func:`build_problem` + :func:`run` fit **one** noisy dataset per arm — cheap,
and what :func:`run_study` uses for the six plots and for the full-study CLI
run, exactly as :mod:`examples.m2_misspecification.study` does. That is a
claim about one draw.

:func:`run_calibration` asks the calibration question properly:
:func:`ampere.results.calibration.sbc` draws a fresh binary from the same
prior the fitted models use, simulates a fresh dataset from the **correct**
(disc-included) truth, refits it under one arm's formulation, and repeats —
so "the flexible fit covers" is a statement about coverage over many draws,
not a lucky one, in the pattern ``tests/m2/test_visibility_calibration.py``
established for this modality's flagship GP. ``tests/interferometry`` pins
its assertions on this route.

Engines
-------
``reference`` samples with emcee (no gradient); ``torch`` and ``jax`` sample
with NUTS through :class:`~ampere.inference.NUTSEngine`. :func:`run`
dispatches on the problem's own ``backend`` flag, never on an argument, for
:mod:`examples.m2_misspecification.study`'s reason: a caller cannot ask for a
sampler the problem cannot feed. The pinned science
(``tests/interferometry``) runs on the always-available reference backend,
deliberately, matching ``tests/m2/conftest.py``'s own reasoning verbatim: the
claim is about the *likelihood*, not about the array library. NUTS (one
modern backend) and :class:`~ampere.inference.SBIEngine` are exercised as
"the modality under every engine" — see :mod:`tests.interferometry` and this
module's ``__main__`` runner for where.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    FittingProblem,
    Product,
    VonMisesFamily,
)
from ampere.results.calibration import sbc

from . import generators as gen
from . import model as _model

__all__ = [
    "ARMS",
    "BACKENDS",
    "BURN_IN",
    "CHROMATIC_BUDGET",
    "CHROMATIC_KERNELS",
    "CI_BUDGET",
    "DOC_BUDGET",
    "FLUX_RATIO_PRIOR",
    "GP_AMPLITUDE_SCALE",
    "GP_LENGTH_SCALE_RANGE",
    "RANK_DRAWS",
    "SEPARATION_PRIOR",
    "SIMULATIONS",
    "STEPS",
    "WALKERS",
    "EmceeBudget",
    "Summary",
    "build_problem",
    "chromatic_problem",
    "coverage_at",
    "model_for",
    "run",
    "run_calibration",
    "run_chromatic_arm",
    "run_study",
    "summarise",
]

#: The three arms, in the order the item's text puts them.
ARMS: tuple[str, ...] = ("correct", "incomplete", "flexible")
#: The three backends. ``reference`` always available; the others need their extra.
BACKENDS: tuple[str, ...] = ("reference", "torch", "jax")

#: The free binary parameters' priors, every arm, every backend —
#: :mod:`interferometry_fixtures.fitted_model`'s own ranges, unchanged, so a
#: reader comparing this study with W4.3's proof sees one convention.
SEPARATION_PRIOR = st.uniform(6.0, 14.0)
FLUX_RATIO_PRIOR = st.uniform(0.1, 0.7)

#: The flexible arm's kernel priors, narrowed around the disc's own
#: correlation scale rather than left at
#: :func:`interferometry_fixtures.two_dataset_problem`'s wide, unidentified
#: range: W4.3's carried note is that a GP over an unidentified hyperparameter
#: takes 209 s and diverges under NUTS, and the reason is structural rather
#: than a sampler setting — a length scale free to roam from kilo- to
#: giga-wavelength lets the process trade off against the companion's own
#: position, which is a genuine second mode in the likelihood, not a mixing
#: problem. The disc's own Gaussian visibility profile
#: (:data:`~examples.interferometry.generators.DISC_FWHM`) has an e-folding
#: scale in (u, v) of order ``1 / (pi fwhm_rad)`` ~ 6e7 wavelengths, so the
#: prior is centred there.
GP_AMPLITUDE_SCALE = 0.1
GP_LENGTH_SCALE_RANGE = (2.0e7, 2.0e8)


@dataclasses.dataclass(frozen=True)
class EmceeBudget:
    """Walkers, steps and burn-in for an ensemble run."""

    walkers: int
    steps: int
    burn_in: int


#: The per-PR budget: small, and what ``tests/interferometry`` uses for both
#: the single-run figures and the calibration study.
CI_BUDGET = EmceeBudget(walkers=16, steps=350, burn_in=150)
#: The documentation / ``__main__`` budget: longer chains, a cleaner corner plot.
DOC_BUDGET = EmceeBudget(walkers=24, steps=1500, burn_in=500)
#: Arm (d) is informational (the item's own word), so its per-PR budget is
#: smaller again: three fits, not pinned, reported as a ranking.
CHROMATIC_BUDGET = EmceeBudget(walkers=12, steps=220, burn_in=80)

#: The calibration study's own numbers, in ``tests/m2/test_visibility_calibration.py``'s
#: shape: cheap enough for a per-PR gate, informative enough to show a
#: direction. See that module's docstring for what a small ``count`` does and
#: does not buy.
SIMULATIONS = 12
RANK_DRAWS = 120
WALKERS = 14
STEPS = 260
BURN_IN = 120

#: Arm (d)'s three noise models, by name, in the order the item reports them.
CHROMATIC_KERNELS: tuple[str, ...] = ("spatial", "spectral", "product")


# ---------------------------------------------------------------------------
# Backend / model resolution
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


def _itf_module(backend: str) -> Any:
    """This backend's ``interferometry`` module, imported lazily."""
    if backend == "reference":
        import ampere.backends.reference.interferometry as module

        return module
    if backend == "torch":
        _backend_module("torch")
        import ampere.backends.torch.interferometry as module

        return module
    if backend == "jax":
        _backend_module("jax")
        import ampere.backends.jax.interferometry as module

        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")


def model_for(backend: str, arm: str, x: Any, y: Any) -> Any:
    """The fitted model for *arm*, on *backend*'s own ``Binary``/composite.

    ``"correct"`` gets :class:`~examples.interferometry.model.BinaryWithDisc`
    with the disc fixed at its true value; ``"incomplete"`` and
    ``"flexible"`` get the plain ``Binary`` (the disc entirely omitted) — the
    misspecification the flexible likelihood is asked to survive.
    """
    itf = _itf_module(backend)
    if arm == "correct":
        return _model.binary_with_disc(
            itf,
            x,
            y,
            disc_flux=gen.DISC_FLUX,
            disc_fwhm=gen.DISC_FWHM,
            separation=SEPARATION_PRIOR,
            position_angle=gen.BINARY["position_angle"],
            flux_ratio=FLUX_RATIO_PRIOR,
            flux=gen.BINARY["flux"],
            component_fwhm=gen.COMPONENT_FWHM,
        )
    if arm in ("incomplete", "flexible"):
        return itf.Binary(
            x,
            y,
            separation=SEPARATION_PRIOR,
            position_angle=gen.BINARY["position_angle"],
            flux_ratio=FLUX_RATIO_PRIOR,
            flux=gen.BINARY["flux"],
            component_fwhm=gen.COMPONENT_FWHM,
        )
    raise ValueError(f"unknown arm {arm!r}; the three are {ARMS!r}.")


def _noise_models(backend: str, arm: str) -> tuple[Any, Any]:
    """``(visibility_noise, closure_phase_noise)`` for *arm*, on *backend*."""
    module = _backend_module(backend)
    t3_noise = module.IndependentNoise()
    if arm != "flexible":
        return module.IndependentNoise(), t3_noise
    kernel = module.Matern32(
        st.halfnorm(scale=GP_AMPLITUDE_SCALE),
        st.loguniform(*GP_LENGTH_SCALE_RANGE),
        axes=("u", "v"),
    )
    return module.GaussianProcessNoise(kernel, module.DenseGP()), t3_noise


# ---------------------------------------------------------------------------
# One noisy dataset, fitted once per arm (the single-run route)
# ---------------------------------------------------------------------------


def build_problem(backend: str, arm: str, *, seed: int = gen.SEED) -> FittingProblem:
    """The study's one binary-plus-disc dataset, fitted under *arm*, on *backend*.

    The truth and the noise draw are always generated on the **reference**
    backend (:func:`~examples.interferometry.generators.synthetic`) — an
    observed container is plain arrays, so the fit itself is what picks a
    backend, not the data.
    """
    reference_itf = _itf_module("reference")
    _, observed_vis, observed_t3 = gen.synthetic(reference_itf, seed=seed)
    itf = _itf_module(backend)
    grid = gen.seed_grid()
    fitted = model_for(backend, arm, grid, grid)
    vis_noise, t3_noise = _noise_models(backend, arm)

    import interferometry_fixtures as fixtures

    vis_instrument = fixtures.chain(itf, observed_vis, "vis")
    t3_instrument = fixtures.chain(itf, observed_t3, "t3")
    datasets = DatasetCollection(
        {
            "vis": Dataset(
                observed_vis,
                vis_instrument,
                likelihood=_likelihood(vis_noise, complex_=True),
                label="vis",
            ),
            "t3": Dataset(
                observed_t3,
                t3_instrument,
                likelihood=_likelihood(t3_noise, complex_=False),
                label="t3",
            ),
        }
    )
    return FittingProblem(fitted, datasets, seed=seed)


def _likelihood(noise: Any, *, complex_: bool) -> Any:
    from ampere.core import Likelihood

    family = ComplexGaussianFamily() if complex_ else VonMisesFamily()
    return Likelihood(family, noise)


def run(problem: FittingProblem, budget: EmceeBudget, *, progress: bool = False) -> Any:
    """Sample *problem* with the engine its backend can feed."""
    from ampere.inference import EmceeEngine, NUTSEngine

    if problem.backend == "reference":
        engine = EmceeEngine(problem, walkers=budget.walkers)
        return engine.run(budget.steps, burn_in=budget.burn_in, progress=progress)
    draws = max(budget.steps - budget.burn_in, 100)
    return NUTSEngine(problem).run(draws, warmup=budget.burn_in, chains=2, progress=progress)


def run_study(
    *,
    backend: str = "reference",
    arms: Sequence[str] = ARMS,
    budget: EmceeBudget = CI_BUDGET,
    seed: int = gen.SEED,
) -> dict[str, dict[str, Any]]:
    """One run per arm, keeping what the figures and the report need."""
    results: dict[str, dict[str, Any]] = {}
    for arm in arms:
        problem = build_problem(backend, arm, seed=seed)
        results[arm] = {"run": run(problem, budget), "problem": problem, "arm": arm}
    return results


# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Summary:
    """One parameter's posterior, reduced to the numbers this study asserts on.

    The same reduction as :class:`examples.m2_misspecification.study.Summary`
    (duplicated rather than imported: this package is meant to stand on its
    own as the template for a new modality, and four small properties are a
    poor reason to depend on a sibling example).
    """

    name: str
    median: float
    low: float
    high: float
    truth: float | None

    @property
    def width(self) -> float:
        return 0.5 * (self.high - self.low)

    @property
    def bias_in_widths(self) -> float | None:
        if self.truth is None:
            return None
        if self.width <= 0.0:
            return float("inf")
        return abs(self.median - self.truth) / self.width

    @property
    def covers_truth(self) -> bool | None:
        if self.truth is None:
            return None
        return bool(self.low <= self.truth <= self.high)


def summarise(run: Any, *, names: Sequence[str] | None = None) -> dict[str, Summary]:
    """Reduce a stored run's posterior to one :class:`Summary` per variable."""
    posterior = run["posterior"]
    available = [str(name) for name in posterior.dataset.data_vars]
    wanted = available if names is None else list(names)
    truths = dict(gen.TRUTH)
    summaries: dict[str, Summary] = {}
    for name in wanted:
        draws = np.asarray(posterior[name].values, dtype=float).ravel()
        low, median, high = (float(v) for v in np.percentile(draws, (16.0, 50.0, 84.0)))
        summaries[name] = Summary(
            name=name, median=median, low=low, high=high, truth=truths.get(name)
        )
    return summaries


# ---------------------------------------------------------------------------
# The calibration route: many small refits, per arm
# ---------------------------------------------------------------------------


def _simulating_problem(*, seed: int = gen.SEED) -> FittingProblem:
    """The correct (disc-included) truth, free binary, independent noise.

    Always the reference backend: :func:`~ampere.results.calibration.sbc`
    calls ``problem.simulate()`` many times, and the reference path's emcee
    fit is what every arm is refitted with, so there is nothing a modern
    backend would add here.
    """
    return build_problem("reference", "correct", seed=seed)


def _calibration_factory(arm: str, backend: str) -> Any:
    """A ``sbc`` ``engine_factory``: refit one replica's data under *arm*."""
    from ampere.inference import EmceeEngine

    itf = _itf_module(backend)
    grid = gen.seed_grid()

    def factory(replica: FittingProblem) -> Any:
        observed_vis = replica.datasets["vis"].observed
        observed_t3 = replica.datasets["t3"].observed
        fitted = model_for(backend, arm, grid, grid)
        vis_noise, t3_noise = _noise_models(backend, arm)

        import interferometry_fixtures as fixtures

        vis_instrument = fixtures.chain(itf, observed_vis, "vis")
        t3_instrument = fixtures.chain(itf, observed_t3, "t3")
        datasets = DatasetCollection(
            {
                "vis": Dataset(
                    observed_vis,
                    vis_instrument,
                    likelihood=_likelihood(vis_noise, complex_=True),
                    label="vis",
                ),
                "t3": Dataset(
                    observed_t3,
                    t3_instrument,
                    likelihood=_likelihood(t3_noise, complex_=False),
                    label="t3",
                ),
            }
        )
        problem = FittingProblem(fitted, datasets, seed=replica.seed)
        return EmceeEngine(problem, walkers=WALKERS)

    return factory


def run_calibration(
    arm: str,
    *,
    backend: str = "reference",
    count: int = SIMULATIONS,
    draws: int = RANK_DRAWS,
    seed: int = 20260913,
) -> Any:
    """``sbc`` over *count* replicas: does *arm*'s credible interval hold up?

    The simulating problem is always the correct, disc-included truth
    (:func:`_simulating_problem`); *arm* names the formulation that refits
    each replica. ``tests/interferometry`` asserts on the ``coverage`` group
    this returns; see that suite for the pinned numbers and margins.
    """
    return sbc(
        _simulating_problem(seed=seed),
        _calibration_factory(arm, backend),
        count=count,
        draws=draws,
        run_options={"steps": STEPS, "burn_in": BURN_IN},
        parameters=list(gen.TRUTH),
        seed=seed,
        label=f"{arm} interferometric fit",
    )


def coverage_at(calibration: Any, nominal: float = 0.9) -> np.ndarray:
    """The empirical coverage of the central *nominal* interval, per parameter.

    ``tests/m2/test_visibility_calibration.py``'s own helper, under this
    module's name: :func:`~ampere.results.calibration.sbc` always lays its
    ``level`` coordinate out on the same grid
    (:data:`~ampere.results.calibration.DEFAULT_LEVELS`), so *nominal* is
    exact on it rather than interpolated.
    """
    levels = np.asarray(calibration["level"].values)
    index = int(np.argmin(np.abs(levels - nominal)))
    if abs(float(levels[index]) - nominal) > 1e-9:
        raise ValueError(f"{nominal} is not one of sbc's levels: {levels.tolist()}.")
    return np.asarray(calibration["coverage"].values)[index]


# ---------------------------------------------------------------------------
# Arm (d): the chromatic case
# ---------------------------------------------------------------------------


def _chromatic_kernel(backend: str, kind: str) -> Any:
    """One of arm (d)'s three noise models, over the visibility axes.

    ``"spatial"`` and ``"spectral"`` fix the amplitude that is not their own
    to ``1.0`` (a product's marginal variance is the product of its terms',
    so both amplitudes free over-parameterises it by one degree of freedom —
    ``ampere.core.kernels.Product``'s own docstring). Units are per-leaf: the
    spectral length scale is in the container's spectral unit, micron.
    """
    module = _backend_module(backend)
    spatial = module.Matern32(
        st.halfnorm(scale=GP_AMPLITUDE_SCALE),
        st.loguniform(*GP_LENGTH_SCALE_RANGE),
        axes=("u", "v"),
    )
    spectral = module.Matern32(
        1.0, st.loguniform(1.0e-3, 0.1), axes=("spectral_axis",), length_scale_unit=u.micron
    )
    if kind == "spatial":
        return spatial
    if kind == "spectral":
        return spectral
    if kind == "product":
        return (
            Product(spatial, spectral)
            if backend == "reference"
            else module.Product(spatial, spectral)
        )
    raise ValueError(f"unknown chromatic kernel {kind!r}; the three are {CHROMATIC_KERNELS!r}.")


def chromatic_problem(backend: str, kind: str, *, seed: int = gen.SEED) -> FittingProblem:
    """Arm (d): the achromatic binary, fitted to the dispersed truth under *kind*."""
    reference_itf = _itf_module("reference")
    _, observed = gen.chromatic_synthetic(reference_itf, seed=seed)
    itf = _itf_module(backend)
    grid = gen.seed_grid()
    fitted = itf.Binary(
        grid,
        grid,
        separation=SEPARATION_PRIOR,
        position_angle=gen.BINARY["position_angle"],
        flux_ratio=FLUX_RATIO_PRIOR,
        flux=gen.BINARY["flux"],
        component_fwhm=gen.COMPONENT_FWHM,
    )
    module = _backend_module(backend)
    kernel = _chromatic_kernel(backend, kind)
    noise = module.GaussianProcessNoise(kernel, module.DenseGP())

    import interferometry_fixtures as fixtures

    vis_instrument = fixtures.chain(itf, observed, "vis")
    datasets = DatasetCollection(
        {
            "vis": Dataset(
                observed, vis_instrument, likelihood=_likelihood(noise, complex_=True), label="vis"
            )
        }
    )
    return FittingProblem(fitted, datasets, seed=seed)


def run_chromatic_arm(
    *, backend: str = "reference", budget: EmceeBudget = CHROMATIC_BUDGET, seed: int = gen.SEED
) -> dict[str, dict[str, Any]]:
    """One run per chromatic kernel: the ranking the item asks be measured."""
    results: dict[str, dict[str, Any]] = {}
    for kind in CHROMATIC_KERNELS:
        problem = chromatic_problem(backend, kind, seed=seed)
        results[kind] = {"run": run(problem, budget), "problem": problem}
    return results
