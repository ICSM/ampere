"""W5.6's solver bake-off, attached to ``pixi run bench``.

``DEVELOPMENT_PLAN.md`` §5's Phase 5 rule is that an approximate 2-3D solver is
**chosen by measurement**. W5.4 landed the first one (``HilbertSpaceGP``), W5.5
measured it against ``DenseGP`` on an image and said the grid-exploiting
candidates were "W5.6's bake-off", and this module is where that bake-off's
numbers are produced on every run::

    pixi run bench                     # the rows below, and benchmark.json
    pixi run -e dev pytest tests/benchmarks -m image_full -s
                                       # the full table, at 128x128 and 256x256

The harness is ``tests/benchmarks/test_gp_solvers.py``'s — pytest-benchmark,
for the reasons that module's docstring records — and the same rule applies:
**nothing here asserts a time**. A benchmark that fails on a slow runner is a
benchmark that gets deleted. These rows fail only when the code under them
raises or returns something non-finite, which is a real failure on any machine.
The *science* the timings feed — bias, coverage, localisation, and the
``|Δ log p|`` against ``DenseGP`` — lives in :mod:`examples.image.bakeoff` and
is reported as a table rather than asserted, exactly as W5.5's benchmark is.

The three questions these rows exist to answer
-----------------------------------------------
1. **What does one evaluation cost, per solver, at realistic N?** The
   ``solver`` rows, at a per-PR size; the ``image_full`` rows at 128x128 and
   256x256, which is the item's 10^4 to 10^5.
2. **Where does EFGP overtake HSGP?** In ``m`` far more than in ``N``: HSGP
   forms ``Phi^T D^-1 Phi`` at ``O(N m²)`` and stores an ``(N, m)`` block,
   EFGP assembles the same matrix from a Toeplitz generator at
   ``O(N 2^d m)`` and stores none. The ``resolution`` rows sweep ``m`` at
   fixed ``N``; the per-PR size sits below the crossover and says so, and the
   ``image_full`` table finds it near ``m = 576`` at 16,384 pixels.
3. **What does the FFT actually buy?** The ``route`` rows, and they are the
   bake-off's central finding: the published ``O(N + m log m)`` is the cost of
   a *solve*, and a marginal likelihood also needs ``log|M|``, for which a
   Toeplitz matrix has no FFT. The two rows time the same normal equations
   through a Cholesky (which gives the determinant) and through
   conjugate gradients on FFT matvecs (which does not).
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
import scipy.linalg

from ampere.core import GaussianFamily
from ampere.core.efgp import (
    fourier_grid,
    solve_iterative,
    spectral_weights,
    toeplitz_generator,
    toeplitz_matrix,
)
from examples.image import bakeoff
from examples.image.grid_gp import GridLikelihood

pytest.importorskip("pytest_benchmark")

#: The size the exact anchor is measured at. 32x32 is 1,024 pixels, which is a
#: 1,024-square Cholesky — tens of milliseconds, and affordable in a per-PR
#: job. 64x64 is four seconds and 640 MiB of it, which is not.
DENSE_PIXELS = 32

#: The size the three approximations are measured at per PR. Large enough that
#: the ``N``-dependent half of each method's cost is doing real work, small
#: enough that the whole file is seconds.
APPROX_PIXELS = 48

#: The sizes the item's own table is measured at, behind ``image_full``.
FULL_PIXELS = (128, 256)

#: The resolution rows' two rungs, **per axis**: ``m = 256`` and ``m = 1024``
#: for HSGP, matched in frequency reach by 17 and 33 for EFGP.
RESOLUTION_RUNGS = (16, 32)

#: The ``route`` rows' grid, per axis: 33 is ``m = 1089``, large enough that
#: the cubic term is visible and small enough to factorise five times.
ROUTE_PER_AXIS = 33

_CASES: dict[int, bakeoff.Case] = {}


def case(pixels: int) -> bakeoff.Case:
    """The image at *pixels*, built once per process.

    Building it is not what is being measured — it is a PSF convolution of a
    truth model — so it is cached rather than timed, which is also what keeps
    a parametrised row from rebuilding it per round.
    """
    found = _CASES.get(pixels)
    if found is None:
        found = bakeoff.study_case(pixels)
        _CASES[pixels] = found
    return found


def _log_prob(pixels: int, solver: Any) -> Callable[[], float]:
    """One ``Likelihood.log_prob`` of the flexible arm, ready to time.

    Through the library path rather than the solver method, so a number here
    is comparable with :func:`examples.image.study.benchmark_solvers`' table
    and with what a sampler actually pays.
    """
    from ampere.backends.reference import noise as noise_module

    built = case(pixels)
    kernel = bakeoff.kernel_for(built)
    likelihood = GridLikelihood(GaussianFamily(), noise_module.GaussianProcessNoise(kernel, solver))

    def call() -> float:
        return likelihood.log_prob(built.predicted, built.observed, values={})

    return call


def _run(benchmark: Any, call: Callable[[], float], group: str, name: str) -> None:
    """Time *call*, and refuse a fast wrong answer."""
    benchmark.group = group
    benchmark.name = name
    # One call outside the timer, so a method with a one-off setup — Vecchia's
    # conditioning sets, which are a constant of the data — is measured on the
    # cost it actually repeats rather than on its first evaluation.
    call()
    value = benchmark(call)
    assert math.isfinite(value), f"{name} returned {value!r}"


# ---------------------------------------------------------------------------
# 1. What one evaluation costs, per solver
# ---------------------------------------------------------------------------


def test_dense_anchor(benchmark: Any) -> None:
    """``DenseGP``: exact, ``O(N³)``, and the reason the other three exist."""
    name, solver = bakeoff.solver_arms(coarse=True, include_dense=True)[0]
    _run(
        benchmark,
        _log_prob(DENSE_PIXELS, solver),
        group=f"image-log-prob-n{DENSE_PIXELS**2}",
        name=f"{name} n={DENSE_PIXELS**2}",
    )


@pytest.mark.parametrize("index", [1, 2, 3], ids=["hsgp", "efgp", "vecchia"])
def test_approximate_solvers(benchmark: Any, index: int) -> None:
    """HSGP, EFGP and Vecchia on the same image, at matched resolution."""
    name, solver = bakeoff.solver_arms(coarse=False, include_dense=True)[index]
    _run(
        benchmark,
        _log_prob(APPROX_PIXELS, solver),
        group=f"image-log-prob-n{APPROX_PIXELS**2}",
        name=f"{name} n={APPROX_PIXELS**2}",
    )


@pytest.mark.image_full
@pytest.mark.parametrize("pixels", FULL_PIXELS)
@pytest.mark.parametrize("index", [1, 2, 3], ids=["hsgp", "efgp", "vecchia"])
def test_approximate_solvers_at_scale(benchmark: Any, pixels: int, index: int) -> None:
    """The item's own sizes: 16,384 and 65,536 pixels, approximate solvers only."""
    name, solver = bakeoff.solver_arms(coarse=False, include_dense=False)[index - 1]
    _run(
        benchmark,
        _log_prob(pixels, solver),
        group=f"image-log-prob-n{pixels**2}",
        name=f"{name} n={pixels**2}",
    )


# ---------------------------------------------------------------------------
# 2. Where EFGP overtakes HSGP: the m crossover
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("per_axis", RESOLUTION_RUNGS)
@pytest.mark.parametrize("family", ["hsgp", "efgp"])
def test_resolution_crossover(benchmark: Any, per_axis: int, family: str) -> None:
    """The same image at two resolutions, for the two reduced-rank solvers.

    HSGP's ``O(N m²)`` assembly and EFGP's ``O(N 2^d m)`` one cross somewhere,
    and a table that swept only ``N`` would never find it. **The crossover is
    not at this size**: it depends on ``N`` as well as ``m`` — the term EFGP
    removes is the one proportional to ``N m²`` — and at 2,304 pixels HSGP is
    ahead at both rungs. These rows are the per-PR record of that ratio; the
    crossover itself is in the ``image_full`` table, measured at 16,384 pixels
    where it falls near ``m = 576``, and reported in ``likelihoods.md`` §7.
    """
    if family == "hsgp":
        solver = bakeoff.GridHilbertSpaceGP(
            basis_size=(per_axis, per_axis), boundary_factor=bakeoff.BOUNDARY_FACTOR
        )
        label = f"HSGP m={per_axis**2}"
    else:
        solver = bakeoff.GridEquispacedFourierGP(
            basis_size=(per_axis + 1, per_axis + 1), boundary_factor=bakeoff.BOUNDARY_FACTOR
        )
        label = f"EFGP m={(per_axis + 1) ** 2}"
    _run(
        benchmark,
        _log_prob(APPROX_PIXELS, solver),
        group=f"image-resolution-m{per_axis**2}",
        name=f"{label} n={APPROX_PIXELS**2}",
    )


# ---------------------------------------------------------------------------
# 3. What the FFT buys, and what it does not
# ---------------------------------------------------------------------------


def _route_pieces() -> tuple[Any, np.ndarray, np.ndarray, np.ndarray]:
    built = case(APPROX_PIXELS)
    kernel = bakeoff.kernel_for(built)
    points = np.ascontiguousarray(kernel.select(built.coordinates))
    grid = fourier_grid(points, (ROUTE_PER_AXIS, ROUTE_PER_AXIS), bakeoff.BOUNDARY_FACTOR)
    root = np.sqrt(spectral_weights(kernel, grid, {}))
    generator = toeplitz_generator(grid, points, 1.0 / built.variance)
    return grid, generator, root, points


def test_route_cholesky(benchmark: Any) -> None:
    """The determinant route: form the Toeplitz matrix and factorise it, ``O(m³)``.

    This is what :class:`~ampere.core.EquispacedFourierGP` does, because
    ``log|M|`` is half of a marginal likelihood and there is no FFT for it.
    """
    grid, generator, root, _ = _route_pieces()

    def call() -> float:
        matrix = toeplitz_matrix(grid, generator) * root[:, None] * root[None, :]
        matrix[np.diag_indices_from(matrix)] += 1.0
        factor = scipy.linalg.cho_factor(matrix, lower=True)
        return 2.0 * float(np.sum(np.log(np.abs(np.diag(factor[0])))))

    _run(
        benchmark,
        call,
        group=f"efgp-normal-equations-m{grid.size}",
        name=f"EFGP Cholesky (value + log-determinant) m={grid.size}",
    )


def test_route_iterative(benchmark: Any) -> None:
    """The published route: CG on FFT matvecs, ``O(m log m)`` an iteration.

    It solves the same system and gives the posterior mean; it gives no
    ``log|M|``, which is why this row sits beside the one above rather than
    replacing it.
    """
    grid, generator, root, _ = _route_pieces()
    right = np.asarray(np.random.default_rng(0).normal(size=grid.size), dtype=complex)

    def call() -> float:
        solution, iterations = solve_iterative(grid, generator, root, right)
        assert iterations > 0
        return float(np.real(np.vdot(right, solution)))

    _run(
        benchmark,
        call,
        group=f"efgp-normal-equations-m{grid.size}",
        name=f"EFGP conjugate gradients (value only) m={grid.size}",
    )


# ---------------------------------------------------------------------------
# The tables the item's decision-log row is written from
# ---------------------------------------------------------------------------


def test_accuracy_panel_has_the_shape_it_claims(capsys: Any) -> None:
    """Bias, coverage and localisation at one small size, printed and checked.

    Not a benchmark row: it asserts the *shape* of the measurement (one cell
    per solver, every number finite, the exact anchor exactly zero against
    itself) and prints the table, which ``-s`` shows. Nothing asserts an
    accuracy, for the reason W5.5's benchmark gives: the plan's rule is that
    the solver is chosen by measurement, so the measurement is reported.
    """
    rows = bakeoff.measure_accuracy(sizes=(24,), coarse=True)
    assert len(rows) == 4
    anchor = rows[0]
    assert anchor.solver == "DenseGP"
    assert anchor.d_log_prob == 0.0
    for row in rows:
        assert math.isfinite(row.d_log_prob)
        assert math.isfinite(row.bias)
        assert math.isfinite(row.rmse)
        assert 0.0 <= row.coverage <= 1.0
        assert -1.0 <= row.localisation <= 1.0
    with capsys.disabled():
        print("\n" + bakeoff.accuracy_table(rows))


@pytest.mark.image_full
def test_full_tables(capsys: Any) -> None:
    """The whole bake-off, as ``likelihoods.md`` §7 records it. Minutes, not seconds."""
    accuracy = bakeoff.measure_accuracy()
    smoothness = bakeoff.measure_smoothness()
    cost = bakeoff.measure_cost()
    resolution = bakeoff.measure_resolution()
    with capsys.disabled():
        print("\n" + bakeoff.accuracy_table(accuracy))
        print("\n" + bakeoff.smoothness_table(smoothness))
        print("\n" + bakeoff.cost_table(cost))
        print("\n" + bakeoff.cost_table(resolution))
    assert accuracy and smoothness and cost and resolution
