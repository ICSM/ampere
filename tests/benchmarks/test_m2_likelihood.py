"""Milestone M2's benchmark table: what one likelihood evaluation costs.

W2.10's acceptance asks for a wall-clock comparison against legacy, "produced
by CI-runnable code, not by hand". This is that code. It runs under the harness
W2.11 chose (pytest-benchmark, ``pixi run -e <env> bench``), writes
``benchmark.json``, and CI uploads that per environment — so the table in
``docs/source/m2_misspecification.rst`` and in W2.10's report is a transcription
of an artefact rather than a set of numbers somebody timed once.

The quantity, and why it is this one
------------------------------------
**One log-density evaluation of the M2 study's own problem at a fixed
parameter vector**, on the 200 / 2 000 / 20 000-point ladder. Not a sampling
run (that is ``test_m2_sampling.py``, and a sampler's cost is this number times
a step count that depends on the geometry), and not the bare GP solve
(``test_gp_solvers.py`` already tracks that, on a synthetic problem). What is
measured here is the whole composed thing a user pays for: forward model,
alignment, noise model, solve.

Two groups, because there are two honest answers
-------------------------------------------------
``m2-log-prob-nN`` is the **contract path** —
:meth:`ampere.core.FittingProblem.log_prob` — which every backend and legacy
ampere can run, and is therefore the only apples-to-apples comparison
available. On torch and jax it is *not* the fast path: the containers convert
to numpy at every step and nothing is jitted, so those rows measure dispatch
overhead as much as arithmetic, and they are here to be compared with legacy
and with the reference backend rather than with each other.

``m2-realised-nN`` is the **realised path** — the differentiable log-density
``ampere.core.realise`` builds and NUTS actually consumes — timed as *value and
gradient*, which is what one leapfrog step costs. It has no reference or legacy
row because neither has a gradient, which is the point of the redesign rather
than a gap in the table.

Legacy
------
``ampere.data.Spectrum.lnlike`` with the paper study's own model class, at 200
and 2 000 points. It is dense and its kernel is a squared exponential, so it
has no 20 000-point row and cannot have one: that is the comparison, not a
missing measurement. Two caveats a reader of the table needs. Legacy's
covariance is ``(I + w M) * sigma_i sigma_j`` — a *relative* correlated
component with a truncation threshold — where v2's is ``K(theta) + diag(sigma^2)``
with an absolute amplitude, so the two are not the same likelihood and their
values are not comparable; only their costs are. And legacy resamples the model
onto the data grid on every call (``setResampler("fast")``), which v2 does not
need here because the model already emits on the observed grid.

Nothing here asserts a time. These rows fail when the code under them raises or
returns something non-finite, which is a real failure on any machine; a
benchmark that fails on a slow runner is a benchmark that gets deleted.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest

pytest.importorskip("pytest_benchmark")

# examples/ is on sys.path via tests/conftest.py (W5.28(k)).
from examples.m2_misspecification import study
from examples.m2_misspecification.generators import generate

#: The scenario every row is timed on. The cost does not depend on which one —
#: the arithmetic is the same — so one is chosen and named rather than looped.
SCENARIO = "strong_smooth"

#: The ladder. The dense solver is only measured on the first two: a 20 000 x
#: 20 000 Cholesky is minutes, and nothing recommends it at that size.
SIZES: tuple[int, ...] = (200, 2_000, 20_000)
DENSE_SIZES: tuple[int, ...] = (200, 2_000)


def _run(
    benchmark: Any,
    call: Callable[[], float],
    *,
    group: str,
    name: str,
    expensive: bool = False,
) -> None:
    """Time *call*, and refuse a fast wrong answer.

    ``expensive=True`` fixes the round count rather than letting
    pytest-benchmark calibrate it, so that one 200 ms dense-solve row cannot
    decide for itself how many seconds a per-PR benchmark job should spend.
    """
    benchmark.group = group
    benchmark.name = name
    if expensive:
        value = benchmark.pedantic(call, rounds=5, iterations=1, warmup_rounds=1)
    else:
        value = benchmark(call)
    assert math.isfinite(value), f"{name} returned {value!r}"


def _contract_call(backend: str, size: int, likelihood: str, solver: str) -> Callable[[], float]:
    """One ``FittingProblem.log_prob`` at the problem's own reference vector."""
    data = generate(SCENARIO, size=size)
    problem = study.build_problem(data, backend=backend, likelihood=likelihood, solver=solver)
    theta = problem.reference_values

    def call() -> float:
        return float(problem.log_prob(theta))

    return call


# ---------------------------------------------------------------------------
# Legacy: the baseline the redesign is measured against.
# ---------------------------------------------------------------------------


def _legacy_call(size: int) -> Callable[[], float]:
    import scipy.stats as st

    from ampere.data import Spectrum as LegacySpectrum
    from ampere.models import Model as LegacyModel
    from ampere.models.results import ModelResults
    from examples.m2_misspecification.model import (
        LINE1_CENTRE,
        LINE1_WIDTH,
        LINE2_CENTRE,
        LINE2_WIDTH,
        PRIOR_LIMITS,
        REFERENCE_WAVELENGTH,
        TRUTH,
    )

    class LegacyToy(LegacyModel):
        """The paper study's ``SimpleSpectralModel``, on the frozen legacy API."""

        def __init__(self, wavelength: Any, **kwargs: Any) -> None:
            self.wavelength = np.asarray(wavelength, dtype=float)
            self.parLabels = ["A", "B", "d1", "d2"]
            self.npars = 4
            self.npars_ptform = 4
            self.priors = [st.uniform(lo, hi - lo) for lo, hi in PRIOR_LIMITS.values()]

        def __call__(self, A: float, B: float, d1: float, d2: float, **kwargs: Any) -> dict:
            grid = self.wavelength
            continuum = A + B * (grid - REFERENCE_WAVELENGTH)
            first = np.exp(-0.5 * ((grid - LINE1_CENTRE) / LINE1_WIDTH) ** 2)
            second = np.exp(-0.5 * ((grid - LINE2_CENTRE) / LINE2_WIDTH) ** 2)
            flux = continuum * (1.0 - d1 * first - d2 * second)
            return {"spectrum": {"wavelength": grid, "flux": flux}}

        def lnprior(self, theta: Any, **kwargs: Any) -> float:
            return float(sum(p.logpdf(t) for p, t in zip(self.priors, theta, strict=True)))

        def prior_transform(self, u: Any, **kwargs: Any) -> np.ndarray:
            return np.array([p.ppf(v) for p, v in zip(self.priors, u, strict=True)])

    data = generate(SCENARIO, size=size)
    spectrum = LegacySpectrum(
        data.wavelength,
        data.observed,
        data.uncertainty,
        "um",
        "Jy",
        calUnc=1e-10,
        scaleLengthPrior=study.GP_LENGTH_SCALE,
        covWeightPrior=study.GP_AMPLITUDE_SCALE,
    )
    spectrum.setResampler(resampleMethod="fast")
    result = ModelResults(**LegacyToy(data.wavelength)(*TRUTH.values()))
    # [calibration scale factor, correlated weight, correlation length]: legacy's
    # three nuisance parameters, at the same length scale the v2 prior centres on.
    nuisance = np.array([1.0, study.GP_AMPLITUDE_SCALE, study.GP_LENGTH_SCALE])

    def call() -> float:
        return float(spectrum.lnlike(nuisance, result))

    return call


@pytest.mark.parametrize("n", DENSE_SIZES)
def test_legacy_dense_rbf(benchmark: Any, n: int) -> None:
    """Legacy ampere's flexible likelihood: dense, squared-exponential, O(N^3)."""
    _run(
        benchmark,
        _legacy_call(n),
        group=f"m2-log-prob-n{n}",
        name=f"legacy dense RBF n={n}",
        expensive=n > 200,
    )


# ---------------------------------------------------------------------------
# The reference backend: always present, never optional.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("n", SIZES)
def test_reference_standard(benchmark: Any, n: int) -> None:
    """The no-GP floor: what the same fit costs with independent noise."""
    _run(
        benchmark,
        _contract_call("reference", n, "standard", "quasisep"),
        group=f"m2-log-prob-n{n}",
        name=f"reference IndependentNoise n={n}",
    )


@pytest.mark.parametrize("n", SIZES)
def test_reference_quasisep(benchmark: Any, n: int) -> None:
    _run(
        benchmark,
        _contract_call("reference", n, "flexible", "quasisep"),
        group=f"m2-log-prob-n{n}",
        name=f"reference QuasisepGP n={n}",
    )


@pytest.mark.parametrize("n", DENSE_SIZES)
def test_reference_dense(benchmark: Any, n: int) -> None:
    _run(
        benchmark,
        _contract_call("reference", n, "flexible", "dense"),
        group=f"m2-log-prob-n{n}",
        name=f"reference DenseGP n={n}",
        expensive=n > 200,
    )


# ---------------------------------------------------------------------------
# The accelerated backends, contract path: comparable with the rows above.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("n", SIZES)
def test_backend_quasisep_contract(benchmark: Any, backend: str, n: int) -> None:
    pytest.importorskip(backend)
    _run(
        benchmark,
        _contract_call(backend, n, "flexible", "quasisep"),
        group=f"m2-log-prob-n{n}",
        name=f"{backend} QuasisepGP (contract path) n={n}",
    )


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("n", DENSE_SIZES)
def test_backend_dense_contract(benchmark: Any, backend: str, n: int) -> None:
    pytest.importorskip(backend)
    _run(
        benchmark,
        _contract_call(backend, n, "flexible", "dense"),
        group=f"m2-log-prob-n{n}",
        name=f"{backend} DenseGP (contract path) n={n}",
        expensive=n > 200,
    )


# ---------------------------------------------------------------------------
# The realised path: value and gradient, which is what a NUTS step costs.
# ---------------------------------------------------------------------------


def _realised_call(backend: str, size: int, solver: str) -> Callable[[], float]:
    data = generate(SCENARIO, size=size)
    problem = study.build_problem(data, backend=backend, likelihood="flexible", solver=solver)
    if backend == "jax":
        import jax

        from ampere.backends.jax import lower_problem

        lowered = lower_problem(problem)
        compiled = jax.jit(jax.value_and_grad(lowered.log_prob_unconstrained))
        point = jax.numpy.zeros(lowered.free_size)
        jax.block_until_ready(compiled(point))  # compile here, not in the timed call

        def call() -> float:
            value, _ = compiled(point)
            return float(jax.block_until_ready(value))

        return call

    import torch

    from ampere.backends.torch import lower_problem as lower_torch

    lowered = lower_torch(problem)

    def call() -> float:
        point = torch.zeros(lowered.free_size, dtype=torch.float64, requires_grad=True)
        value = lowered.log_prob_unconstrained(point)
        value.backward()
        return float(value.detach())

    return call


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("n", SIZES)
def test_backend_quasisep_realised(benchmark: Any, backend: str, n: int) -> None:
    """The O(N) solve with its gradient — the rung of the ladder legacy cannot reach."""
    pytest.importorskip(backend)
    _run(
        benchmark,
        _realised_call(backend, n, "quasisep"),
        group=f"m2-realised-n{n}",
        name=f"{backend} QuasisepGP (value+grad) n={n}",
    )


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("n", DENSE_SIZES)
def test_backend_dense_realised(benchmark: Any, backend: str, n: int) -> None:
    """The dense baseline on the same backend, so the ratio is within-library."""
    pytest.importorskip(backend)
    _run(
        benchmark,
        _realised_call(backend, n, "dense"),
        group=f"m2-realised-n{n}",
        name=f"{backend} DenseGP (value+grad) n={n}",
        expensive=n > 200,
    )
