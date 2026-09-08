"""Milestone M2's other benchmark row: what a short sampling run costs.

W2.10 asks for "one short sampling run per backend at 200 points" alongside the
per-evaluation table in :mod:`tests.benchmarks.test_m2_likelihood`. The two say
different things and both are needed. A log-density evaluation is a property of
the *problem*; a sampling run is a property of the problem **and** the sampler,
and the ratio between them — how many evaluations a sampler spends per stored
draw — is exactly where an ensemble sampler and a gradient sampler differ. NUTS
costs far more per draw and buys far more per draw, and only a row like this one
shows both halves of that trade.

The budgets here are deliberately tiny: a hundred-odd evaluations, enough to
time and nowhere near enough to converge. **These rows are not evidence about
the posterior** — that is ``tests/m2``, and the milestone budget behind
``pytest -m m2_full``. They are evidence about cost.

``rounds`` is fixed rather than calibrated (:meth:`pedantic`), because
pytest-benchmark's calibration would otherwise decide for itself how many
seconds-long sampling runs a per-PR job should pay for.
"""

from __future__ import annotations

import math
import pathlib
import sys
from collections.abc import Callable
from typing import Any

import pytest

pytest.importorskip("pytest_benchmark")

_ROOT = pathlib.Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from examples.m2_misspecification import study  # noqa: E402
from examples.m2_misspecification.generators import generate  # noqa: E402

SCENARIO = "strong_smooth"
SIZE = 200

#: Short on purpose. ``EmceeBudget(16, 150, 50)`` is 2 400 evaluations;
#: ``NutsBudget(50, 50, 1)`` is a hundred draws' worth of leapfrog steps.
BENCH_EMCEE = study.EmceeBudget(walkers=16, steps=150, burn_in=50)
BENCH_NUTS = study.NutsBudget(draws=50, warmup=50, chains=1)


def _run(benchmark: Any, call: Callable[[], float], *, group: str, name: str) -> None:
    benchmark.group = group
    benchmark.name = name
    value = benchmark.pedantic(call, rounds=2, iterations=1, warmup_rounds=0)
    assert math.isfinite(value), f"{name} returned {value!r}"


def _sampling_call(backend: str, likelihood: str) -> Callable[[], float]:
    data = generate(SCENARIO, size=SIZE)

    def call() -> float:
        problem = study.build_problem(data, backend=backend, likelihood=likelihood)
        budget = BENCH_EMCEE if backend == "reference" else BENCH_NUTS
        run = study.run(problem, budget)
        # Something cheap, finite and load-bearing: a run that failed to emit a
        # posterior is a failed benchmark, not a fast one.
        return float(run["posterior"]["model.A"].values.mean())

    return call


@pytest.mark.parametrize("likelihood", list(study.LIKELIHOODS))
def test_reference_emcee(benchmark: Any, likelihood: str) -> None:
    """emcee on the numpy path: the only sampler the reference backend has."""
    _run(
        benchmark,
        _sampling_call("reference", likelihood),
        group=f"m2-sampling-{likelihood}",
        name=f"reference emcee {BENCH_EMCEE.walkers}x{BENCH_EMCEE.steps} ({likelihood})",
    )


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("likelihood", list(study.LIKELIHOODS))
def test_backend_nuts(benchmark: Any, backend: str, likelihood: str) -> None:
    """NUTS over the realised density, warm-up included, compilation included."""
    pytest.importorskip(backend)
    _run(
        benchmark,
        _sampling_call(backend, likelihood),
        group=f"m2-sampling-{likelihood}",
        name=f"{backend} NUTS {BENCH_NUTS.draws}+{BENCH_NUTS.warmup} ({likelihood})",
    )
