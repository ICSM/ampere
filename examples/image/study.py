"""The study: three arms on one image, and the solver benchmark the phase asks for.

The shape :mod:`examples.m2_misspecification.study` established and
:mod:`examples.interferometry.study` carried to a second modality, applied to a
third: everything above this module (:mod:`.generators`, :mod:`.model`)
supplies data and physics, everything below it (:mod:`.figures`,
``tests/examples/test_image_study.py``) consumes what it returns.

The three arms
--------------
One source-plus-background truth (:mod:`.generators`), observed once as an
``Image``, fitted three ways:

``"correct"``
    The background is in the model, at its true value, and the noise is
    independent. The control.
``"incomplete"``
    The background is absent; the noise is still independent. What a naive fit
    of an under-modelled image looks like.
``"flexible"``
    The background is still absent, and the likelihood carries
    :class:`~ampere.core.GaussianProcessNoise` over a two-axis Matérn-3/2
    kernel on ``(x, y)``. M2's question, asked of a gridded observable: does
    the flexible likelihood keep the source's parameters calibrated when the
    sky model is wrong, by absorbing the background's smooth spatial
    correlation into the residual?

The benchmark, which is the phase's own
---------------------------------------
``DEVELOPMENT_PLAN.md`` §5's rule for Phase 5 is that a solver is "chosen by
measurement". :func:`benchmark_solvers` is this item's measurement: the same
flexible likelihood, at three image sizes, under the exact
:class:`~ampere.core.DenseGP` and under W5.4's approximate
:class:`~ampere.core.HilbertSpaceGP` with a tensor-product basis, reporting
wall clock and peak memory for one ``log_prob``. It is **informational** — the
item's own word — and nothing in ``tests/`` asserts on the numbers, only that
the table has the shape it claims.

Note what both solvers are and are not doing: neither exploits the grid. A
Kronecker or SKI solver would, and that is W5.6's bake-off; these two treat an
image as ``N`` scattered points that happen to lie on a lattice, which is the
honest baseline the next item's candidates have to beat.

Read :mod:`.grid_gp` first
---------------------------
Every flexible arm here goes through :mod:`.grid_gp`, which lifts — out of
tree, from the public API — the two ``Layout.POINTS`` gates that otherwise
refuse a correlated noise model on an ``Image``. That module explains what is
blocked, why the block is a Phase 5 slot rather than mathematics, and what the
library change it stands in for would be.
"""

from __future__ import annotations

import dataclasses
import gc
import time
import tracemalloc
from collections.abc import Sequence
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    negotiate,
)

from . import generators as gen
from .grid_gp import GridDenseGP, GridHilbertSpaceGP, GridLikelihood

__all__ = [
    "ARMS",
    "BACKENDS",
    "BENCHMARK_SIZES",
    "CI_BUDGET",
    "DOC_BUDGET",
    "SMALL_PIXELS",
    "EmceeBudget",
    "SolverCost",
    "Summary",
    "benchmark_solvers",
    "build_problem",
    "coverage_at",
    "run",
    "run_calibration",
    "run_study",
    "summarise",
]

#: The three arms, in the order the report reads them.
ARMS: tuple[str, ...] = ("correct", "incomplete", "flexible")

#: The three backends. ``reference`` is always available; the others need
#: their extra.
BACKENDS: tuple[str, ...] = ("reference", "torch", "jax")

#: The image size every *fit* in this study runs at. Small, and deliberately:
#: an exact ``DenseGP`` over a 24x24 image is a 576-square Cholesky per
#: evaluation, which an ensemble sampler can afford thousands of times. The
#: larger sizes appear only in :func:`benchmark_solvers`, where one ``log_prob``
#: is measured rather than sampled.
SMALL_PIXELS = 24

#: The three sizes the benchmark reports, as the item asks. 128x128 is 16,384
#: pixels, which is where an exact solve stops being a choice.
BENCHMARK_SIZES: tuple[int, ...] = (24, 64, 128)

#: Basis functions **per axis** for the HSGP arm; the total is the square.
#: Eight per axis is sixty-four, which is ample for a background whose
#: correlation length is a good fraction of the field — and the point of the
#: method is that it does not grow with N.
BASIS_PER_AXIS = 8

#: The Hilbert-space boundary factor. Riutort-Mayol et al.'s rule of thumb is
#: ``c >= 1.2 max(length_scale) / S`` for a data half-extent ``S``; the kernel
#: prior below is centred well inside that at ``c = 2``.
BOUNDARY_FACTOR = 2.0

#: The flexible arm's kernel amplitude prior, **in units of the observation's
#: own declared uncertainty**. That scaling is not a convenience: an image of a
#: source a few milliarcseconds across is measured in Jy/sr, so its values are
#: of order 10^15, and a prior written as a bare number would be fifteen orders
#: of magnitude from anything the data could support — a GP pinned at zero
#: amplitude, indistinguishable at a glance from a flexible likelihood that
#: simply did not help. (It is how this study's first run came out, and the
#: fix is this constant.) The uncertainty is a *declared constant of the
#: observation*, not a statistic of the values, so scaling a nuisance prior to
#: it is the ordinary "a few times the noise" statement rather than a peek at
#: the data.
GP_AMPLITUDE_SIGMAS = 1.5

#: The flexible arm's kernel length-scale prior, mas. Centred on the
#: background's own width rather than left wide: a GP whose length scale is
#: free to roam from a pixel to the whole field can trade off against the
#: source itself, which is a second mode in the likelihood rather than a mixing
#: problem — :mod:`examples.interferometry.study` records the same hazard in
#: (u, v), and the fix is the same, an informed prior.
GP_LENGTH_PRIOR = (4.0, 20.0)


@dataclasses.dataclass(frozen=True)
class EmceeBudget:
    """Walkers, steps and burn-in for an ensemble run."""

    walkers: int
    steps: int
    burn_in: int


#: The per-PR budget: what ``tests/examples/test_image_study.py`` uses.
CI_BUDGET = EmceeBudget(walkers=12, steps=200, burn_in=80)
#: The documentation / ``__main__`` budget: longer chains, a cleaner corner.
DOC_BUDGET = EmceeBudget(walkers=20, steps=1200, burn_in=400)

#: The calibration study's numbers, in ``examples/interferometry``'s shape:
#: cheap enough for a per-PR gate, informative enough to show a direction.
#: Twelve simulations is a smoke budget in Talts et al.'s sense — evidence of a
#: gross effect, not a fine one.
SIMULATIONS = 12
RANK_DRAWS = 120
CALIBRATION_BUDGET = EmceeBudget(walkers=12, steps=180, burn_in=70)


# ---------------------------------------------------------------------------
# Backend resolution
# ---------------------------------------------------------------------------


def _backend_module(backend: str) -> Any:
    """The namespace this backend's ``PSFConvolution`` lives in."""
    if backend == "reference":
        from ampere.backends import reference as module

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
    """The namespace this backend's image *models* live in (W4.1's placement)."""
    if backend == "reference":
        from ampere.backends.reference import interferometry as module

        return module
    if backend == "torch":
        from ampere.backends.torch import interferometry as module  # type: ignore[assignment]

        return module
    if backend == "jax":
        from ampere.backends.jax import configure_x64
        from ampere.backends.jax import interferometry as module  # type: ignore[assignment]

        configure_x64()
        return module
    raise ValueError(f"unknown backend {backend!r}; the three are {list(BACKENDS)}.")


def _noise_module(backend: str) -> Any:
    """Where this backend's ``IndependentNoise``/``GaussianProcessNoise`` live."""
    if backend == "reference":
        import ampere.core as module

        return module
    return _backend_module(backend)


# ---------------------------------------------------------------------------
# The model, per arm
# ---------------------------------------------------------------------------


def model_for(backend: str, arm: str, *, fitted: bool = True) -> Any:
    """The sky model this *arm* fits with, on *backend*.

    ``"correct"`` carries the background at its true, fixed value; the other
    two do not have it at all (see :mod:`.model` for why "not at all" rather
    than "at zero flux").
    """
    if arm not in ARMS:
        raise ValueError(f"unknown arm {arm!r}; the three are {list(ARMS)}.")
    grid = gen.seed_grid()
    kwargs: dict[str, Any] = {
        "flux": st.lognorm(0.5, scale=gen.SOURCE["flux"]) if fitted else gen.SOURCE["flux"],
        "fwhm": st.uniform(1.0, 6.0) if fitted else gen.SOURCE["fwhm"],
    }
    if arm == "correct":
        kwargs["background_flux"] = gen.BACKGROUND["flux"]
        kwargs["background_fwhm"] = gen.BACKGROUND["fwhm"]
    if backend == "torch":
        from . import model_torch

        return model_torch.build(grid, grid, **kwargs)
    if backend == "jax":
        from . import model_jax

        return model_jax.build(grid, grid, **kwargs)
    from .model import source_with_background

    return source_with_background(_itf_module(backend), grid, grid, **kwargs)


def _noise_for(backend: str, arm: str, *, sigma: float, solver: Any = None) -> Any:
    """Independent noise, or the flexible arm's two-axis GP over ``(x, y)``.

    *sigma* is the observation's declared per-pixel uncertainty, and the
    kernel's amplitude prior is written in units of it; see
    :data:`GP_AMPLITUDE_SIGMAS` for why that is required rather than tidy.
    """
    module = _noise_module(backend)
    if arm != "flexible":
        return module.IndependentNoise()
    kernel = module.Matern32(
        st.halfnorm(scale=GP_AMPLITUDE_SIGMAS * float(sigma)),
        st.uniform(GP_LENGTH_PRIOR[0], GP_LENGTH_PRIOR[1] - GP_LENGTH_PRIOR[0]),
        axes=("x", "y"),
    )
    return module.GaussianProcessNoise(kernel, solver or GridDenseGP())


def declared_sigma(observed: Any) -> float:
    """The observation's declared per-pixel uncertainty, as one number."""
    return float(np.median(np.asarray(observed.uncertainty, dtype=float)))


def _likelihood(arm: str, noise: Any) -> Likelihood:
    """A ``GridLikelihood`` for the flexible arm, the shipped one otherwise.

    The distinction is :mod:`.grid_gp`'s subject and should not survive the
    library change that module proposes: an independent noise model never asks
    for coordinates, so only the correlated arm needs the override at all.
    """
    if arm == "flexible":
        return GridLikelihood(GaussianFamily(), noise)
    return Likelihood(GaussianFamily(), noise)


# ---------------------------------------------------------------------------
# The problem
# ---------------------------------------------------------------------------


def build_problem(
    backend: str,
    arm: str,
    *,
    pixels: int = SMALL_PIXELS,
    seed: int = gen.SEED,
    solver: Any = None,
) -> FittingProblem:
    """One arm's :class:`~ampere.core.FittingProblem` on a fresh synthetic image."""
    module = _backend_module(backend)
    _, observed = gen.synthetic(module, _itf_module(backend), pixels, seed=seed)
    instrument = gen.camera(module, observed)
    fitted = model_for(backend, arm)
    noise = _noise_for(backend, arm, sigma=declared_sigma(observed), solver=solver)
    datasets = DatasetCollection(
        {
            "image": Dataset(
                observed,
                instrument,
                likelihood=_likelihood(arm, noise),
                label="image",
            )
        }
    )
    return FittingProblem(fitted, datasets, seed=seed)


def run(problem: FittingProblem, budget: EmceeBudget, *, progress: bool = False) -> Any:
    """Sample *problem* with the engine its backend can feed.

    Dispatches on the problem's own ``backend`` flag rather than on an
    argument, for :mod:`examples.m2_misspecification.study`'s reason: a caller
    cannot ask for a sampler the problem cannot feed.
    """
    from ampere.inference import EmceeEngine, NUTSEngine

    if problem.backend == "reference":
        engine = EmceeEngine(problem, walkers=budget.walkers)
        return engine.run(budget.steps, burn_in=budget.burn_in, progress=progress)
    draws = max(budget.steps - budget.burn_in, 80)
    return NUTSEngine(problem).run(draws, warmup=budget.burn_in, chains=2, progress=progress)


def run_study(
    *,
    backend: str = "reference",
    arms: Sequence[str] = ARMS,
    budget: EmceeBudget = CI_BUDGET,
    pixels: int = SMALL_PIXELS,
    seed: int = gen.SEED,
) -> dict[str, dict[str, Any]]:
    """One run per arm, keeping what the figures and the report need."""
    results: dict[str, dict[str, Any]] = {}
    for arm in arms:
        problem = build_problem(backend, arm, pixels=pixels, seed=seed)
        results[arm] = {"run": run(problem, budget), "problem": problem, "arm": arm}
    return results


# ---------------------------------------------------------------------------
# Reading a run
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Summary:
    """One parameter's posterior, reduced to the numbers this study asserts on.

    The same reduction :mod:`examples.interferometry.study` uses, duplicated
    rather than imported for that module's own reason: each study package is
    meant to stand on its own as the template for a new modality.
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


def summarise(run_result: Any, *, names: Sequence[str] | None = None) -> dict[str, Summary]:
    """Reduce a stored run's posterior to one :class:`Summary` per variable."""
    posterior = run_result["posterior"]
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
# The calibration route: many small refits, one arm
# ---------------------------------------------------------------------------


def _simulating_problem(*, pixels: int = SMALL_PIXELS, seed: int = gen.SEED) -> FittingProblem:
    """The correct (background-included) truth, free source, independent noise."""
    return build_problem("reference", "correct", pixels=pixels, seed=seed)


def _calibration_factory(arm: str, backend: str, pixels: int, budget: EmceeBudget) -> Any:
    """An ``sbc`` ``engine_factory``: refit one replica's image under *arm*."""
    from ampere.inference import EmceeEngine

    module = _backend_module(backend)

    def factory(replica: FittingProblem) -> Any:
        observed = replica.datasets["image"].observed
        instrument = gen.camera(module, observed)
        fitted = model_for(backend, arm)
        noise = _noise_for(backend, arm, sigma=declared_sigma(observed))
        datasets = DatasetCollection(
            {
                "image": Dataset(
                    observed,
                    instrument,
                    likelihood=_likelihood(arm, noise),
                    label="image",
                )
            }
        )
        problem = FittingProblem(fitted, datasets, seed=replica.seed)
        return EmceeEngine(problem, walkers=budget.walkers)

    return factory


def run_calibration(
    arm: str,
    *,
    backend: str = "reference",
    pixels: int = SMALL_PIXELS,
    count: int = SIMULATIONS,
    draws: int = RANK_DRAWS,
    budget: EmceeBudget = CALIBRATION_BUDGET,
    seed: int = gen.SEED,
) -> Any:
    """``sbc`` over *count* replicas: does *arm*'s credible interval hold up?

    The simulating problem is always the correct, background-included truth;
    *arm* names the formulation each replica is refitted under.
    ``tests/examples/test_image_study.py`` asserts on the ``coverage`` group
    this returns.
    """
    from ampere.results.calibration import sbc

    return sbc(
        _simulating_problem(pixels=pixels, seed=seed),
        _calibration_factory(arm, backend, pixels, budget),
        count=count,
        draws=draws,
        run_options={"steps": budget.steps, "burn_in": budget.burn_in},
        parameters=list(gen.TRUTH),
        seed=seed,
        label=f"{arm} image fit",
    )


def coverage_at(calibration: Any, nominal: float = 0.9) -> np.ndarray:
    """The empirical coverage of the central *nominal* interval, per parameter."""
    levels = np.asarray(calibration["level"].values)
    index = int(np.argmin(np.abs(levels - nominal)))
    if abs(float(levels[index]) - nominal) > 1e-9:
        raise ValueError(f"{nominal} is not one of sbc's levels: {levels.tolist()}.")
    return np.asarray(calibration["coverage"].values)[index]


# ---------------------------------------------------------------------------
# The benchmark: DenseGP against HilbertSpaceGP, at three N
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class SolverCost:
    """One ``(N, solver)`` cell of the benchmark table."""

    pixels: int
    solver: str
    samples: int
    seconds: float
    peak_mib: float
    log_prob: float

    def row(self) -> str:
        """One fixed-width line, for pasting into a report."""
        return (
            f"{self.pixels:>3}x{self.pixels:<3} {self.samples:>6d}  {self.solver:<12s} "
            f"{self.seconds:>9.3f}  {self.peak_mib:>9.1f}  {self.log_prob:>14.2f}"
        )


def _measure(callable_: Any) -> tuple[float, float, float]:
    """``(seconds, peak MiB, returned value)`` for one call.

    :mod:`tracemalloc` rather than ``resource.getrusage``: the question is how
    much this *solve* allocates, and a high-water mark for the process would be
    dominated by whatever ran before it. It counts Python-level allocations,
    which for numpy means the arrays themselves, which is what a Cholesky of an
    ``N x N`` covariance is made of.
    """
    gc.collect()
    tracemalloc.start()
    start = time.perf_counter()
    value = callable_()
    seconds = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return seconds, peak / (1024.0 * 1024.0), float(value)


def benchmark_solvers(
    *,
    sizes: Sequence[int] = BENCHMARK_SIZES,
    backend: str = "reference",
    basis_per_axis: int = BASIS_PER_AXIS,
    seed: int = gen.SEED,
    include_dense: bool = True,
) -> list[SolverCost]:
    """One ``log_prob`` per ``(N, solver)``: wall clock and peak memory.

    **Informational.** The plan's Phase 5 rule is that a solver is chosen by
    measurement, and this is the measurement — not an assertion. Both solvers
    score the *same* flexible likelihood on the *same* image, at the same
    hyperparameter values, so the difference is the solve and nothing else; the
    ``log_prob`` column is reported so that a reader can see the approximation's
    size beside its cost, which is the trade being measured.

    ``include_dense=False`` skips the exact solver, for the sizes where an
    ``N x N`` Cholesky is not something to run by accident.
    """
    module = _backend_module(backend)
    itf = _itf_module(backend)
    noise_module = _noise_module(backend)
    costs: list[SolverCost] = []
    for pixels in sizes:
        _, observed = gen.synthetic(module, itf, pixels, seed=seed)
        instrument = gen.camera(module, observed)
        fitted = model_for(backend, "flexible", fitted=False)
        compiled = fitted.compile_for(negotiate([instrument]))
        predicted = instrument(compiled.evaluate())
        # The amplitude is in the data's own units, for :data:`GP_AMPLITUDE_SIGMAS`'
        # reason; held fixed here rather than fitted, because the benchmark is
        # about the cost of one solve and not about a posterior.
        kernel = noise_module.Matern32(0.5 * declared_sigma(observed), 8.0, axes=("x", "y"))
        solvers: list[tuple[str, Any]] = []
        if include_dense:
            solvers.append(("DenseGP", GridDenseGP()))
        solvers.append(
            (
                f"HSGP({basis_per_axis}x{basis_per_axis})",
                GridHilbertSpaceGP(
                    basis_size=(basis_per_axis, basis_per_axis),
                    boundary_factor=BOUNDARY_FACTOR,
                ),
            )
        )
        for name, solver in solvers:
            likelihood = GridLikelihood(
                GaussianFamily(), noise_module.GaussianProcessNoise(kernel, solver)
            )
            seconds, peak, value = _measure(
                lambda lik=likelihood, pred=predicted, obs=observed: lik.log_prob(
                    pred, obs, values={}
                )
            )
            costs.append(
                SolverCost(
                    pixels=pixels,
                    solver=name,
                    samples=pixels * pixels,
                    seconds=seconds,
                    peak_mib=peak,
                    log_prob=value,
                )
            )
    return costs


def benchmark_table(costs: Sequence[SolverCost]) -> str:
    """:func:`benchmark_solvers`' result as a fixed-width table."""
    header = (
        f"{'image':>7} {'N':>6}  {'solver':<12s} {'seconds':>9}  {'peak MiB':>9}  {'log_prob':>14}"
    )
    return "\n".join([header, "-" * len(header), *(cost.row() for cost in costs)])
