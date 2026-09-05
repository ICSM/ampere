"""Wall-clock scaling of the GP solvers: W2.3's "demonstrated, not claimed".

``DEVELOPMENT_PLAN.md`` §4.4 makes the O(N) GP the project's distinguishing
feature, and W2.3's acceptance criterion is deliberately empirical:
"10^3-10^5-point scaling demonstrated (wall-clock, not asymptotics claimed)".
This module measures it.

It is **not** part of any default gate: ``pixi run test-all`` enumerates
``tests/core tests/results tests/conformance tests/backends`` and this
directory is not among them, because a dense 10⁴-point solve costs the better
part of a minute and several gigabytes. Run it deliberately::

    pixi run scaling

The assertions are chosen to survive a shared CI runner. Nothing here asserts
an absolute time: the load-bearing claim is the *slope* of wall-clock against
N on log-log axes, which separates a linear solver from a quadratic one by a
factor that no amount of machine noise closes. ``DenseGP`` is capped at 3000
points — the largest size worth paying for in a test — and the correctness
cross-check runs there rather than at 10⁵, where no dense answer exists to
compare against.
"""

from __future__ import annotations

import math
import time
from collections.abc import Callable

import numpy as np
import pytest

from ampere.core import DenseGP, GPSolver, Matern32, QuasisepGP

#: Sizes the O(N) solver is measured at: the decades W2.3 names, plus two
#: intermediate points so the fitted slope has something to fit.
QUASISEP_SIZES: tuple[int, ...] = (1_000, 3_000, 10_000, 30_000, 100_000)

#: Sizes the O(N³) solver is measured at. 10⁴ is feasible but costs ~80 s and
#: ~2.5 GB in ampere's dense path (the ``(n, n, 1)`` separation array
#: dominates), which is more than a test should spend; 10⁵ would need 80 GB
#: for the covariance alone and is infeasible outright. Capping it here, and
#: saying why, is part of the demonstration.
DENSE_SIZES: tuple[int, ...] = (1_000, 3_000)

AMPLITUDE = 0.4
LENGTH_SCALE = 2.0
SIGMA = 0.1
#: Points per length scale, held fixed as N grows, so every size describes the
#: same physical problem sampled longer rather than a different one.
SAMPLING = 10.0


def problem(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """An irregular, ordered 1-D problem of *n* points, deterministic in *n*."""
    rng = np.random.default_rng(20260905 + n)
    span = n * LENGTH_SCALE / SAMPLING
    coordinates = np.sort(rng.uniform(0.0, span, n))[:, None]
    residual = rng.normal(0.0, 0.3, n)
    variance = np.full(n, SIGMA**2)
    return coordinates, residual, variance


def fastest(call: Callable[[], float], repeats: int) -> tuple[float, float]:
    """Best-of-*repeats* wall-clock seconds, and the value computed."""
    best = math.inf
    value = math.nan
    for _ in range(repeats):
        started = time.perf_counter()
        value = call()
        best = min(best, time.perf_counter() - started)
    return best, value


def measure(solver: GPSolver, sizes: tuple[int, ...], repeats: int) -> dict[int, float]:
    """Wall-clock seconds per size, and a finite value at every one."""
    kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
    values = kernel.resolve(None)
    timings: dict[int, float] = {}
    for n in sizes:
        coordinates, residual, variance = problem(n)

        def solve(
            coordinates: np.ndarray = coordinates,
            residual: np.ndarray = residual,
            variance: np.ndarray = variance,
        ) -> float:
            return solver.log_marginal_likelihood(kernel, coordinates, residual, variance, values)

        elapsed, value = fastest(solve, repeats)
        assert math.isfinite(value), f"{solver.NAME} returned {value!r} at n={n}"
        timings[n] = elapsed
    return timings


def slope(timings: dict[int, float]) -> float:
    """The fitted exponent p in ``time ~ N**p``, over the measured sizes."""
    sizes = np.log(np.array(sorted(timings), dtype=float))
    seconds = np.log(np.array([timings[n] for n in sorted(timings)], dtype=float))
    return float(np.polyfit(sizes, seconds, 1)[0])


def table(quasisep: dict[int, float], dense: dict[int, float]) -> str:
    lines = [
        "",
        f"{'N':>8} {'QuasisepGP':>14} {'DenseGP':>14} {'speed-up':>10}",
        f"{'':>8} {'(ms)':>14} {'(ms)':>14} {'':>10}",
    ]
    for n in sorted(set(quasisep) | set(dense)):
        q = quasisep.get(n)
        d = dense.get(n)
        lines.append(
            f"{n:>8} "
            f"{f'{q * 1e3:.2f}' if q else '—':>14} "
            f"{f'{d * 1e3:.2f}' if d else 'not measured':>14} "
            f"{f'{d / q:.0f}x' if q and d else '—':>10}"
        )
    lines.append(
        f"\nfitted exponent p in time ~ N**p: QuasisepGP {slope(quasisep):.2f}"
        + (f", DenseGP {slope(dense):.2f} (two sizes only)" if len(dense) > 1 else "")
    )
    lines.append(
        f"DenseGP above {max(dense) if dense else 0} points is omitted by choice, not by "
        "failure: 10^4 runs in ~80 s and ~2.5 GB, 10^5 would need ~80 GB for the covariance."
    )
    return "\n".join(lines)


def test_the_quasiseparable_solver_scales_linearly(capsys: pytest.CaptureFixture[str]) -> None:
    """10³ → 10⁵ points, measured. A quadratic solver cannot pass this."""
    quasisep = measure(QuasisepGP(), QUASISEP_SIZES, repeats=5)
    dense = measure(DenseGP(), DENSE_SIZES, repeats=1)
    with capsys.disabled():
        print(table(quasisep, dense))

    exponent = slope(quasisep)
    # True O(N) is p == 1. The bound is generous — cache effects and a shared
    # runner both inflate it — but nowhere near the p ≈ 2 an accidentally
    # quadratic implementation would show, which is the discrimination the
    # row exists to make.
    assert exponent < 1.5, f"QuasisepGP scaled as N**{exponent:.2f}, which is not linear"
    # And the point of it all: at the largest size both solvers can reach, the
    # O(N) path is faster by a margin no runner's noise touches. The dense
    # exponent is reported in the table rather than asserted: over the two
    # sizes a test can afford, matrix construction still competes with the
    # cubic factorisation, so the fitted slope is a weak statistic.
    largest = max(DENSE_SIZES)
    assert dense[largest] / quasisep[largest] > 100.0


def test_the_two_solvers_still_agree_at_the_largest_dense_size() -> None:
    """Speed is worthless if it is speed at the wrong answer."""
    kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
    values = kernel.resolve(None)
    for n in DENSE_SIZES:
        coordinates, residual, variance = problem(n)
        expected = DenseGP().log_marginal_likelihood(
            kernel, coordinates, residual, variance, values
        )
        got = QuasisepGP().log_marginal_likelihood(kernel, coordinates, residual, variance, values)
        # Absolute, and scaled by N: both solvers accumulate rounding over the
        # whole problem, and a fixed absolute bound would be a claim about the
        # dense path's own error rather than about their agreement.
        assert got == pytest.approx(expected, abs=1e-9 * n)
