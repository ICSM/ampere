"""The tracked benchmark suite: what the GP solve costs, on every backend present.

W2.11 is where ``DEVELOPMENT_PLAN.md`` §6's "benchmark harness (Phase 2):
pytest-benchmark vs asv" was settled, and this module is the answer's shape.
Run it with one command::

    pixi run bench            # the reference backend, in `dev`
    pixi run -e torch bench   # and the torch rows as well
    pixi run -e jax bench     # and the jax rows as well

Each run writes ``benchmark.json`` — the machine-readable record CI uploads as
the run's artefact, and the thing that makes two runs comparable rather than
two tables that scrolled past.

Why pytest-benchmark and not asv
---------------------------------
The two are not close, once the requirement is written out. What this project
needs is: results as **CI artefacts on a PR**, comparable across runs, and one
command locally.

* **Environments.** asv owns its environments — it builds them itself from a
  matrix in ``asv.conf.json``, out of a wheel it builds from the checkout. This
  project's environments are pixi's, and deliberately so: three of them
  (``dev``, ``torch``, ``jax``) differ in exactly the way the benchmarks care
  about, one of them pins torch to a CPU-only index, and AGENTS.md's working
  agreement is that every CI command is a pixi task. An asv matrix would be a
  second, divergent definition of the same three environments, and the first
  time it drifted the benchmark numbers would silently stop describing what CI
  runs. pytest-benchmark is a pytest plugin: it inherits the environment it is
  invoked in, and the rows below skip or run according to the same
  ``importorskip`` the rest of the suite uses.
* **What "tracked" means here.** asv's value is its own history database and
  web front end, which want either a committed ``.asv/`` results tree or a
  separate results repository. AGENTS.md ground rule 7 forbids run outputs in
  git, so the history would have to live outside the repository — real
  infrastructure, for a project whose benchmark requirement today is "the
  numbers are attached to the PR that changed them". ``--benchmark-json``
  attached to a workflow run is exactly that, and it is a superset of what a
  reviewer needs: full per-round statistics, the machine and interpreter, and
  the commit.
* **Reuse.** These rows need ampere's fixtures, its skip logic, and (on jax)
  ``configure_x64``. In pytest that is the code already written. In asv it
  would be rewritten in asv's own class-based API, and the O(N) claims would
  then be measured by code the conformance suite has never run.

asv's genuine advantages — regression detection against history, machine
calibration, bisection — are the things to revisit if and when this project
wants a benchmark *server* rather than benchmark artefacts. Nothing here
forecloses it: the measured quantities live in ordinary functions.

What is measured, and what is not
----------------------------------
The GP marginal likelihood, because ``DEVELOPMENT_PLAN.md`` §4.4 makes the
O(N) GP the project's distinguishing feature and §5's Phase 2 line asks for
that feature's cost to be tracked. On the differentiable backends the row is
the **value and gradient**, jitted where the backend jits, because that — not
the value alone — is what one NUTS leapfrog step costs.

Sizes are small on purpose. This suite runs on every PR, on a shared runner,
so it must cost seconds; ``tests/scaling`` is where the 10³ to 10⁵ demonstration
lives and it stays out of every gate. The two suites answer different
questions: ``tests/scaling`` asserts a *slope* (a claim that survives a noisy
runner), while this one records *absolute* times (a number that does not, and
is therefore compared only against other runs of the same job).

Nothing here asserts a time. A benchmark that fails on a slow runner is a
benchmark that gets deleted; these rows fail only when the code under them
raises or returns something non-finite, which is a real failure on any
machine.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

from ampere.core import DenseGP, GPSolver, Kernel, Matern32, QuasisepGP

pytest.importorskip("pytest_benchmark")

AMPLITUDE = 0.4
LENGTH_SCALE = 2.0
SIGMA = 0.1
#: Points per length scale, held fixed as N grows, so every size is the same
#: physical problem sampled longer. Matches ``tests/scaling``.
SAMPLING = 10.0

#: The one size the O(N³) solver is measured at. 1000 points is ~10 ms of
#: dense Cholesky here; 3000 is nearly a second, which is more than a
#: per-PR job should spend on a solver nothing recommends at that size.
DENSE_SIZE = 1_000

#: The sizes the O(N) solver is measured at: one shared with the dense row, so
#: the two are directly comparable, and one decade up, where the difference is
#: the whole point.
QUASISEP_SIZES: tuple[int, ...] = (1_000, 10_000)


def problem(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """An irregular, ordered 1-D problem of *n* points, deterministic in *n*.

    Identical to ``tests/scaling``'s generator, seed included, so a number
    measured here and a number measured there describe the same problem.
    """
    rng = np.random.default_rng(20260905 + n)
    span = n * LENGTH_SCALE / SAMPLING
    coordinates = np.sort(rng.uniform(0.0, span, n))[:, None]
    residual = rng.normal(0.0, 0.3, n)
    variance = np.full(n, SIGMA**2)
    return coordinates, residual, variance


def _run(benchmark: Any, call: Callable[[], float], group: str, name: str) -> None:
    """Time *call*, and refuse a fast wrong answer."""
    benchmark.group = group
    benchmark.name = name
    value = benchmark(call)
    assert math.isfinite(value), f"{name} returned {value!r}"


# ---------------------------------------------------------------------------
# The reference backend (``ampere.core``): always present, never optional.
# ---------------------------------------------------------------------------


def _reference_call(solver: GPSolver, kernel: Kernel, n: int) -> Callable[[], float]:
    coordinates, residual, variance = problem(n)
    values = kernel.resolve(None)

    def call() -> float:
        return solver.log_marginal_likelihood(kernel, coordinates, residual, variance, values)

    return call


@pytest.mark.parametrize("n", QUASISEP_SIZES)
def test_reference_quasisep(benchmark: Any, n: int) -> None:
    """The O(N) solve, in numpy, through celerite2's compiled kernels."""
    kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
    _run(
        benchmark,
        _reference_call(QuasisepGP(), kernel, n),
        group=f"gp-marginal-likelihood-n{n}",
        name=f"reference QuasisepGP n={n}",
    )


def test_reference_dense(benchmark: Any) -> None:
    """The O(N³) solve at the one size worth paying for, as the baseline."""
    kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
    _run(
        benchmark,
        _reference_call(DenseGP(), kernel, DENSE_SIZE),
        group=f"gp-marginal-likelihood-n{DENSE_SIZE}",
        name=f"reference DenseGP n={DENSE_SIZE}",
    )


# ---------------------------------------------------------------------------
# The torch backend: value **and** gradient, which is what NUTS pays for.
# ---------------------------------------------------------------------------


def _torch_call(solver_name: str, n: int) -> Callable[[], float]:
    torch = pytest.importorskip("torch")
    from ampere.backends import torch as backend

    coordinates, residual, variance = problem(n)
    kernel = backend.Matern32(AMPLITUDE, LENGTH_SCALE)
    solver = getattr(backend, solver_name)()

    def call() -> float:
        amplitude = torch.tensor(AMPLITUDE, dtype=torch.float64, requires_grad=True)
        length_scale = torch.tensor(LENGTH_SCALE, dtype=torch.float64, requires_grad=True)
        value = solver.log_marginal_likelihood_tensor(
            kernel,
            coordinates,
            residual,
            variance,
            {"amplitude": amplitude, "length_scale": length_scale},
        )
        value.backward()
        # ``detach`` because the tensor is still attached to the graph at this
        # point and ``float()`` on it warns; the timed work is already done.
        return float(value.detach())

    return call


@pytest.mark.parametrize("n", QUASISEP_SIZES)
def test_torch_quasisep(benchmark: Any, n: int) -> None:
    """The O(N) solve under ``torch.autograd``, forward and reverse."""
    _run(
        benchmark,
        _torch_call("QuasisepGP", n),
        group=f"gp-marginal-likelihood-n{n}",
        name=f"torch QuasisepGP (value+grad) n={n}",
    )


def test_torch_dense(benchmark: Any) -> None:
    """The dense baseline on the same backend, so the ratio is within-library."""
    _run(
        benchmark,
        _torch_call("DenseGP", DENSE_SIZE),
        group=f"gp-marginal-likelihood-n{DENSE_SIZE}",
        name=f"torch DenseGP (value+grad) n={DENSE_SIZE}",
    )


# ---------------------------------------------------------------------------
# The jax backend: value and gradient, jitted. Compilation is excluded — it is
# paid once per shape, and a sampler pays it once per run.
# ---------------------------------------------------------------------------


def _jax_call(solver_name: str, n: int) -> Callable[[], float]:
    jax = pytest.importorskip("jax")
    from ampere.backends import jax as backend

    # ``lowering.md`` §10.2(a): the *application* turns x64 on, and a test
    # module is an application. Every class below refuses without it.
    backend.configure_x64()

    coordinates, residual, variance = problem(n)
    kernel = backend.Matern32(AMPLITUDE, LENGTH_SCALE)
    solver = getattr(backend, solver_name)()
    values = kernel.resolve(None)
    axis = np.ascontiguousarray(coordinates[:, 0])

    def density(amplitude: Any) -> Any:
        return solver.log_marginal_likelihood_jax(
            kernel, axis, residual, variance, {**values, "amplitude": amplitude}
        )

    compiled = jax.jit(jax.value_and_grad(density))
    argument = jax.numpy.asarray(AMPLITUDE)
    jax.block_until_ready(compiled(argument))  # compile here, not in the timed call

    def call() -> float:
        value, _ = compiled(argument)
        return float(jax.block_until_ready(value))

    return call


@pytest.mark.parametrize("n", QUASISEP_SIZES)
def test_jax_quasisep(benchmark: Any, n: int) -> None:
    """The O(N) solve under XLA, value and gradient, compilation excluded."""
    _run(
        benchmark,
        _jax_call("QuasisepGP", n),
        group=f"gp-marginal-likelihood-n{n}",
        name=f"jax QuasisepGP (value+grad, jitted) n={n}",
    )


def test_jax_dense(benchmark: Any) -> None:
    """The dense baseline on the same backend, for the same within-library reason."""
    _run(
        benchmark,
        _jax_call("DenseGP", DENSE_SIZE),
        group=f"gp-marginal-likelihood-n{DENSE_SIZE}",
        name=f"jax DenseGP (value+grad, jitted) n={DENSE_SIZE}",
    )
