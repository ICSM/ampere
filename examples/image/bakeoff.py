"""W5.6's solver bake-off: EFGP and Vecchia against HSGP and DenseGP on the image.

``DEVELOPMENT_PLAN.md`` §5's Phase 5 rule is that a 2-3D solver is **chosen by
measurement**, not by argument, and this module is the measurement.
:mod:`examples.image.study`'s :func:`~examples.image.study.benchmark_solvers`
took the first half of it at W5.5 — ``DenseGP`` against ``HilbertSpaceGP`` on
wall clock and memory — and said in as many words that "neither exploits the
grid ... that is W5.6's bake-off". This module is W5.6's: the same image, the
same misspecification, the same kernel, four solvers, and five quantities.

What is measured, and why each one
----------------------------------
``bias``, ``rmse``, ``coverage``, ``localisation``
    Not against ``DenseGP``'s answer but against the **truth**: the image's
    omitted component is a known field (a PSF-convolved broad Gaussian), so
    "did the GP find it?" has an answer that does not depend on which solver
    is asked. ``bias`` is the mean signed error over pixels in units of the
    declared uncertainty, ``rmse`` the same unsigned, ``coverage`` the
    fraction of pixels whose truth lies inside the solver's own 90% interval —
    an approximation that reports honest uncertainty must keep this near 0.9,
    and one that reports the exact prior beside an approximate correction will
    not — and ``localisation`` the correlation between the conditional mean
    and the true field, which is M2's question ("did the GP put its correction
    where the deviation actually is?") asked of a plane.
``d_log_prob``
    The disagreement with ``DenseGP``'s marginal likelihood, in nats. The
    quantity a sampler actually consumes, and the one the conformance class of
    W5.4 holds an approximate solver to.
``seconds``, ``peak MiB``, ``setup``
    Wall clock and peak allocation for one ``log_prob``, with the **first**
    call reported separately: a Vecchia solver's conditioning sets are a
    constant of the data and are built once (``ampere.core.vecchia``), exactly
    as a spectral solver's grid is, and hiding that build inside a
    per-evaluation number would flatter the spectral solvers as surely as
    charging it every evaluation would flatter nothing at all.

The kernel is held fixed across every arm, at an amplitude taken from the
truth field's own r.m.s. and a length scale of :data:`LENGTH_SCALE`. That is
deliberate and it is what makes the comparison about the *solve*: every solver
is scoring the same declared model on the same data, so a difference in any
column is the approximation and nothing else.

What this module is not
-----------------------
It is not a fit. Bias and coverage here are of the **GP's conditional field**,
which one ``condition`` call produces, rather than of a posterior over the
source's parameters, which would need one MCMC run per (solver, N, kernel)
cell and is not what a bake-off between solvers needs to decide. The image
study's own :func:`~examples.image.study.run_study` and
:func:`~examples.image.study.run_calibration` are where parameter bias and
SBC coverage live, on the arm the bake-off's winner will eventually carry.
"""

from __future__ import annotations

import dataclasses
import gc
import time
import tracemalloc
from collections.abc import Callable, Sequence
from typing import Any

import numpy as np

from ampere.core import (
    EquispacedFourierGP,
    GaussianFamily,
    VecchiaResponseGP,
    negotiate,
    sample_coordinates,
)

from . import generators as gen
from . import study
from .grid_gp import GriddedSolver, GridDenseGP, GridHilbertSpaceGP, GridLikelihood

__all__ = [
    "ACCURACY_SIZES",
    "BOUNDARY_FACTOR",
    "COST_SIZES",
    "LENGTH_SCALE",
    "LENGTH_SCALES",
    "RESOLUTIONS",
    "SMOOTHNESS",
    "GridEquispacedFourierGP",
    "GridVecchiaResponseGP",
    "SolverAccuracy",
    "SolverCost",
    "accuracy_table",
    "cost_table",
    "measure_accuracy",
    "measure_cost",
    "measure_resolution",
    "measure_smoothness",
    "smoothness_table",
    "solver_arms",
    "study_case",
]


class GridEquispacedFourierGP(GriddedSolver, EquispacedFourierGP):
    """:class:`~ampere.core.EquispacedFourierGP`, allowed to see an ``Image``.

    :mod:`examples.image.grid_gp`'s mixin, applied to W5.6's prototype for
    exactly the reason that module records: the ``Layout.GRID`` refusal is a
    closed Phase 5 slot rather than mathematics, and lifting it out of tree
    keeps ``ampere/core/likelihood.py`` untouched while the library change is
    with Peter.
    """


class GridVecchiaResponseGP(GriddedSolver, VecchiaResponseGP):
    """:class:`~ampere.core.VecchiaResponseGP`, allowed to see an ``Image``."""


#: The kernel length scale every arm is measured at, mas. The omitted
#: background is 17 mas across, so this is the scale of the thing the GP is
#: being asked to absorb; ``study.benchmark_solvers`` uses the same number.
LENGTH_SCALE = 8.0

#: The spectral solvers' box/spacing factor. The study's own, so that the HSGP
#: column here and the HSGP column in W5.5's table describe the same solver.
BOUNDARY_FACTOR = study.BOUNDARY_FACTOR

#: Sizes the accuracy panel runs at: every one small enough that ``DenseGP``
#: can supply the exact answer beside the three approximations. 64x64 is where
#: W5.5's benchmark stopped and it is where this one stops too — a dense
#: 128x128 covariance is about 10 GiB across the copies a Cholesky makes.
ACCURACY_SIZES: tuple[int, ...] = (24, 48, 64)

#: Sizes the cost panel runs at. The item's "realistic N" is 10^4-10^5: 128x128
#: is 16,384 pixels and 256x256 is 65,536. ``DenseGP`` appears only at the
#: first, for the reason above.
COST_SIZES: tuple[int, ...] = (64, 128, 256)

#: The smoothness ladder, as the item asks: Matérn-1/2 through 5/2.
SMOOTHNESS: tuple[str, ...] = ("Matern12", "Matern32", "Matern52")


@dataclasses.dataclass(frozen=True)
class Case:
    """One image, with the truth field the flexible arm has to find."""

    pixels: int
    observed: Any
    #: The source-only prediction the flexible arm scores its residual
    #: against, as the instrument produced it — an ``Image``, so
    #: ``Likelihood.log_prob`` takes it unchanged.
    predicted: Any
    coordinates: np.ndarray
    residual: np.ndarray
    variance: np.ndarray
    truth: np.ndarray
    sigma: float

    @property
    def samples(self) -> int:
        """``N``, the pixel count."""
        return int(self.residual.shape[0])


def study_case(pixels: int, *, backend: str = "reference", seed: int = gen.SEED) -> Case:
    """The image at *pixels*, its residual, and the component the model omits.

    The residual is ``observed - source-only model``, which is exactly what the
    flexible arm's GP is handed, and ``truth`` is the PSF-convolved background
    the misspecified model does not have — the deviation the GP is supposed to
    absorb. Both are built from :mod:`examples.image.generators`' own truth, so
    nothing here re-derives the physics.
    """
    module = study._backend_module(backend)
    itf = study._itf_module(backend)
    _, observed = gen.synthetic(module, itf, pixels, seed=seed)
    instrument = gen.camera(module, observed)
    requirements = negotiate([instrument])
    complete = instrument(gen.truth_model(itf, pixels).compile_for(requirements).evaluate())
    partial = instrument(
        study.model_for(backend, "incomplete", fitted=False).compile_for(requirements).evaluate()
    )
    values = np.asarray(observed.values, dtype=float).reshape(-1)
    complete_values = np.asarray(complete.values, dtype=float).reshape(-1)
    partial_values = np.asarray(partial.values, dtype=float).reshape(-1)
    sigma = study.declared_sigma(observed)
    return Case(
        pixels=pixels,
        observed=observed,
        predicted=partial,
        coordinates=np.ascontiguousarray(sample_coordinates(observed)),
        residual=values - partial_values,
        variance=np.full(values.shape[0], sigma * sigma),
        truth=complete_values - partial_values,
        sigma=sigma,
    )


def kernel_for(
    case: Case,
    family: str = "Matern32",
    *,
    length_scale: float = LENGTH_SCALE,
    backend: str = "reference",
) -> Any:
    """The one kernel every arm is measured on: fixed, and the same for all.

    The amplitude is the truth field's own r.m.s. rather than a prior draw,
    because a bake-off asks which *solver* computes a declared model best and
    a hyperparameter that differed between arms would answer a different
    question. It is a constant of the generator, not a statistic of the data.
    """
    module = study._noise_module(backend)
    amplitude = float(np.sqrt(np.mean(case.truth**2)))
    return getattr(module, family)(amplitude, length_scale, axes=("x", "y"))


def solver_arms(*, coarse: bool = False, include_dense: bool = True) -> list[tuple[str, Any]]:
    """The arms, at matched resolution: HSGP, EFGP and Vecchia, plus the anchor.

    ``m`` is matched *by frequency reach* rather than by count. HSGP's ``m``
    members per axis run over positive frequencies up to ``pi m / (2 L)``;
    EFGP's grid of ``2M + 1`` runs over ``+-M pi / L``, so ``2M + 1 = m + 1``
    reaches the same place. That is why the counts below are 8 against 9, and
    16 against 17, rather than being equal numbers that would quietly compare
    two different bandwidths.
    """
    hilbert, fourier, width = (8, 9, 10) if coarse else (16, 17, 30)
    arms: list[tuple[str, Any]] = []
    if include_dense:
        arms.append(("DenseGP", GridDenseGP()))
    arms.append(
        (
            f"HSGP m={hilbert**2}",
            GridHilbertSpaceGP(basis_size=(hilbert, hilbert), boundary_factor=BOUNDARY_FACTOR),
        )
    )
    arms.append(
        (
            f"EFGP m={fourier**2}",
            GridEquispacedFourierGP(basis_size=(fourier, fourier), boundary_factor=BOUNDARY_FACTOR),
        )
    )
    arms.append((f"Vecchia k={width}", GridVecchiaResponseGP(neighbours=width, seed=gen.SEED)))
    return arms


# ---------------------------------------------------------------------------
# Accuracy: against the truth, and against DenseGP
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class SolverAccuracy:
    """One ``(N, solver)`` cell of the accuracy table."""

    pixels: int
    samples: int
    solver: str
    d_log_prob: float
    bias: float
    rmse: float
    coverage: float
    localisation: float

    def row(self) -> str:
        """One fixed-width line, for pasting into a report."""
        return (
            f"{self.pixels:>3}x{self.pixels:<3} {self.samples:>6d}  {self.solver:<24s} "
            f"{self.d_log_prob:>11.3e}  {self.bias:>+8.3f}  {self.rmse:>7.3f}  "
            f"{self.coverage:>8.3f}  {self.localisation:>8.4f}"
        )


def _accuracy(case: Case, kernel: Any, solver: Any, exact_log_prob: float) -> tuple[float, ...]:
    values: dict[str, Any] = {}
    log_prob = solver.log_marginal_likelihood(
        kernel, case.coordinates, case.residual, case.variance, values
    )
    conditional = solver.condition(kernel, case.coordinates, case.residual, case.variance, values)
    error = np.asarray(conditional.mean, dtype=float) - case.truth
    deviation = np.sqrt(np.maximum(np.asarray(conditional.variance, dtype=float), 0.0))
    # 1.6449 is the 95th percentile of the standard normal: a central 90%
    # interval, the same nominal level ``study.coverage_at`` uses.
    inside = np.abs(error) <= 1.6448536269514722 * deviation
    correlation = float(np.corrcoef(np.asarray(conditional.mean, dtype=float), case.truth)[0, 1])
    return (
        abs(log_prob - exact_log_prob),
        float(np.mean(error)) / case.sigma,
        float(np.sqrt(np.mean(error**2))) / case.sigma,
        float(np.mean(inside)),
        correlation,
    )


def measure_accuracy(
    *,
    sizes: Sequence[int] = ACCURACY_SIZES,
    backend: str = "reference",
    coarse: bool = False,
    seed: int = gen.SEED,
) -> list[SolverAccuracy]:
    """Bias, coverage, localisation and ``|Δ log p|`` for every ``(N, solver)``.

    **Informational**, as W5.5's benchmark is: the plan's rule is that a solver
    is chosen by measurement, and this is the measurement rather than an
    assertion. ``tests/benchmarks/test_solver_bakeoff.py`` asserts only that
    each cell is finite and that the exact anchor behaves.
    """
    results: list[SolverAccuracy] = []
    for pixels in sizes:
        case = study_case(pixels, backend=backend, seed=seed)
        kernel = kernel_for(case, "Matern32", backend=backend)
        exact = GridDenseGP().log_marginal_likelihood(
            kernel, case.coordinates, case.residual, case.variance, {}
        )
        for name, solver in solver_arms(coarse=coarse, include_dense=True):
            measured = _accuracy(case, kernel, solver, exact)
            results.append(
                SolverAccuracy(
                    pixels=pixels,
                    samples=case.samples,
                    solver=name,
                    d_log_prob=measured[0],
                    bias=measured[1],
                    rmse=measured[2],
                    coverage=measured[3],
                    localisation=measured[4],
                )
            )
    return results


def accuracy_table(rows: Sequence[SolverAccuracy]) -> str:
    """:func:`measure_accuracy`'s result as a fixed-width table."""
    header = (
        f"{'image':>7} {'N':>6}  {'solver':<24s} {'|dlogp|':>11}  {'bias/s':>8}  "
        f"{'rmse/s':>7}  {'cover90':>8}  {'localis.':>8}"
    )
    return "\n".join([header, "-" * len(header), *(row.row() for row in rows)])


# ---------------------------------------------------------------------------
# Smoothness: the same cell, across the Matérn ladder
# ---------------------------------------------------------------------------


#: The two correlation lengths the smoothness panel is measured at, mas, on a
#: 24 mas field. They are not a refinement of one another: at 8 mas the kernel
#: reaches a third of the image — the regime the flexible likelihood is
#: actually used in, because that is the scale of the background it absorbs —
#: and at 2 mas it reaches a twelfth, which is where a screening-based method
#: has a near field to screen with.
LENGTH_SCALES: tuple[float, ...] = (8.0, 2.0)


def measure_smoothness(
    *,
    pixels: int = 64,
    length_scales: Sequence[float] = LENGTH_SCALES,
    backend: str = "reference",
    coarse: bool = False,
    seed: int = gen.SEED,
) -> list[SolverAccuracy]:
    """The accuracy panel at one ``N``, swept over :data:`SMOOTHNESS`.

    The row the bake-off turns on. A reduced-rank spectral method's error is
    the tail of the spectral density it truncates, which for an isotropic
    Matérn-``nu`` in ``d`` axes decays as ``omega ** -(2 nu + d)``; a Vecchia
    approximation's error is the far field its conditioning sets fail to
    screen, and Stein (2002) shows screening works *better* the rougher the
    process. The two therefore fail in opposite directions, and this is where
    that shows.
    """
    case = study_case(pixels, backend=backend, seed=seed)
    results: list[SolverAccuracy] = []
    for length_scale in length_scales:
        for family in SMOOTHNESS:
            kernel = kernel_for(case, family, length_scale=length_scale, backend=backend)
            exact = GridDenseGP().log_marginal_likelihood(
                kernel, case.coordinates, case.residual, case.variance, {}
            )
            for name, solver in solver_arms(coarse=coarse, include_dense=False):
                measured = _accuracy(case, kernel, solver, exact)
                results.append(
                    SolverAccuracy(
                        pixels=pixels,
                        samples=case.samples,
                        solver=f"M{family[6:]} l={length_scale:g} {name}",
                        d_log_prob=measured[0],
                        bias=measured[1],
                        rmse=measured[2],
                        coverage=measured[3],
                        localisation=measured[4],
                    )
                )
    return results


def smoothness_table(rows: Sequence[SolverAccuracy]) -> str:
    """:func:`measure_smoothness`'s result as a fixed-width table."""
    header = (
        f"{'image':>7} {'N':>6}  {'kernel/solver':<24s} {'|dlogp|':>11}  {'bias/s':>8}  "
        f"{'rmse/s':>7}  {'cover90':>8}  {'localis.':>8}"
    )
    return "\n".join([header, "-" * len(header), *(row.row() for row in rows)])


# ---------------------------------------------------------------------------
# Cost: wall clock, peak memory, and the one-off setup each method carries
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class SolverCost:
    """One ``(N, solver)`` cell of the cost table."""

    pixels: int
    samples: int
    solver: str
    setup_seconds: float
    seconds: float
    peak_mib: float
    log_prob: float

    def row(self) -> str:
        """One fixed-width line, for pasting into a report."""
        return (
            f"{self.pixels:>3}x{self.pixels:<3} {self.samples:>6d}  {self.solver:<14s} "
            f"{self.setup_seconds:>9.3f}  {self.seconds:>9.3f}  {self.peak_mib:>9.1f}  "
            f"{self.log_prob:>14.2f}"
        )


def _measure(call: Callable[[], float]) -> tuple[float, float, float]:
    """``(seconds, peak MiB, value)`` for one call.

    :mod:`tracemalloc` rather than ``resource.getrusage``, for
    :func:`examples.image.study.benchmark_solvers`' reason: the question is
    what this *solve* allocates, and a process high-water mark would be
    dominated by whatever ran before it.
    """
    gc.collect()
    tracemalloc.start()
    start = time.perf_counter()
    value = call()
    seconds = time.perf_counter() - start
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return seconds, peak / (1024.0 * 1024.0), float(value)


def measure_cost(
    *,
    sizes: Sequence[int] = COST_SIZES,
    backend: str = "reference",
    coarse: bool = False,
    dense_below: int = 65,
    seed: int = gen.SEED,
) -> list[SolverCost]:
    """Wall clock and peak memory for one ``log_prob``, per ``(N, solver)``.

    Measured through :class:`examples.image.grid_gp.GridLikelihood`, not
    through the solver method, so the numbers are directly comparable with
    :func:`examples.image.study.benchmark_solvers`' table — the library path a
    sampler actually takes, alignment checks included.

    The **first** call is timed separately as ``setup``: a Vecchia solver
    builds its conditioning sets once per coordinate set and caches them, so
    charging that build to every evaluation would misreport the method by the
    ratio of a fit's length. ``dense_below`` keeps the exact solver out of the
    sizes where an ``N x N`` Cholesky is not something to run by accident.
    """
    results: list[SolverCost] = []
    for pixels in sizes:
        case = study_case(pixels, backend=backend, seed=seed)
        kernel = kernel_for(case, "Matern32", backend=backend)
        noise_module = study._noise_module(backend)
        for name, solver in solver_arms(coarse=coarse, include_dense=pixels < dense_below):
            likelihood = GridLikelihood(
                GaussianFamily(), noise_module.GaussianProcessNoise(kernel, solver)
            )

            def call(
                lik: Any = likelihood, pred: Any = case.predicted, obs: Any = case.observed
            ) -> float:
                return lik.log_prob(pred, obs, values={})

            setup_seconds, _, _ = _measure(call)
            seconds, peak, value = _measure(call)
            results.append(
                SolverCost(
                    pixels=pixels,
                    samples=case.samples,
                    solver=name,
                    setup_seconds=setup_seconds,
                    seconds=seconds,
                    peak_mib=peak,
                    log_prob=value,
                )
            )
    return results


def cost_table(rows: Sequence[SolverCost]) -> str:
    """:func:`measure_cost`'s result as a fixed-width table."""
    header = (
        f"{'image':>7} {'N':>6}  {'solver':<14s} {'setup s':>9}  {'seconds':>9}  "
        f"{'peak MiB':>9}  {'log_prob':>14}"
    )
    return "\n".join([header, "-" * len(header), *(row.row() for row in rows)])


#: The resolution ladder the ``m``-scaling panel sweeps, **per axis**. The
#: total ``m`` is the square, so this is 64, 256, 576 and 1024 basis members.
RESOLUTIONS: tuple[int, ...] = (8, 16, 24, 32)


def measure_resolution(
    *,
    pixels: int = 128,
    resolutions: Sequence[int] = RESOLUTIONS,
    backend: str = "reference",
    seed: int = gen.SEED,
) -> list[SolverCost]:
    """Cost against ``m`` at fixed ``N``: the panel that decides the EFGP column.

    The ``(N, solver)`` table alone cannot settle EFGP against HSGP, because
    the two differ in their dependence on **m** rather than on N: HSGP forms
    ``Phi^T D^-1 Phi`` at ``O(N m²)`` and stores an ``(N, m)`` block, EFGP
    assembles the same matrix from a Toeplitz generator at ``O(N 2^d m)`` and
    stores none. At small ``m`` the difference is swamped by EFGP's complex
    exponentials costing several times what HSGP's sines do; the claim is
    about where the quadratic term takes over, and a table that swept only
    ``N`` would report the small-``m`` regime as though it were the method.

    Both solvers are given the same frequency reach at each rung
    (:func:`solver_arms`' 8-against-9 rule), so a row compares two
    approximations of comparable quality and not merely two array sizes.
    """
    case = study_case(pixels, backend=backend, seed=seed)
    kernel = kernel_for(case, "Matern32", backend=backend)
    noise_module = study._noise_module(backend)
    results: list[SolverCost] = []
    for per_axis in resolutions:
        arms: list[tuple[str, Any]] = [
            (
                f"HSGP m={per_axis**2}",
                GridHilbertSpaceGP(
                    basis_size=(per_axis, per_axis), boundary_factor=BOUNDARY_FACTOR
                ),
            ),
            (
                f"EFGP m={(per_axis + 1) ** 2}",
                GridEquispacedFourierGP(
                    basis_size=(per_axis + 1, per_axis + 1), boundary_factor=BOUNDARY_FACTOR
                ),
            ),
        ]
        for name, solver in arms:
            likelihood = GridLikelihood(
                GaussianFamily(), noise_module.GaussianProcessNoise(kernel, solver)
            )

            def call(
                lik: Any = likelihood, pred: Any = case.predicted, obs: Any = case.observed
            ) -> float:
                return lik.log_prob(pred, obs, values={})

            setup_seconds, _, _ = _measure(call)
            seconds, peak, value = _measure(call)
            results.append(
                SolverCost(
                    pixels=pixels,
                    samples=case.samples,
                    solver=name,
                    setup_seconds=setup_seconds,
                    seconds=seconds,
                    peak_mib=peak,
                    log_prob=value,
                )
            )
    return results
