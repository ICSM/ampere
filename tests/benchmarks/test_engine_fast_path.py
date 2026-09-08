"""What one gradient-free proposal costs, on each route (W2.5 slice 3).

Run it with the rest of the tracked suite::

    pixi run -e jax bench     # and the jax rows as well
    pixi run -e torch bench   # and the torch rows as well

``tests/benchmarks/test_gp_solvers.py``'s module docstring carries the harness
choice and the rules; this module obeys them, with one deliberate exception
stated below.

What is measured
----------------
One call to :meth:`ampere.inference.Engine.log_prob` — a single proposal, as
emcee, dynesty or zeus would score it — on a problem under the flexible
likelihood, twice: through the numpy **contract** path
(``use_realisation=False``), and through the backend's **realisation**
(``use_realisation=True``, the default since W2.5 slice 3). The two rows share
a benchmark group so the comparison is the thing the table shows.

This is the quantity the fast path exists for. W2.10 measured the jax contract
path at about 28 ms flat per ``log_prob`` whatever the problem's size, because
``FittingProblem.evaluate`` is numpy and a jax model evaluated an operation at
a time pays jax's dispatch on each one — which made the *gradient-free* engines
slower on the backend built for speed. A `log_prob` is what a sampler spends
its life doing, so it is the right unit; a whole run would fold in emcee's own
overhead and the walker count and say less.

The one assertion, and why it is allowed here
----------------------------------------------
The sibling module says, correctly, that nothing in this suite asserts a time:
"a benchmark that fails on a slow runner is a benchmark that gets deleted".
:func:`test_the_fast_path_beats_the_contract_path` asserts a **ratio** instead,
and a ratio is a different kind of claim: both halves are measured in the same
process, on the same machine, within milliseconds of each other, so the runner's
speed divides out. The threshold is set an order of magnitude below the measured
margin for the same reason — the claim being defended is "this is not a
pessimisation", not "this is exactly 50 times faster".
"""

from __future__ import annotations

import dataclasses
import importlib
import math
import time
from collections.abc import Callable
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
from ampere.inference import EmceeEngine

pytest.importorskip("pytest_benchmark")

REFERENCE_WAVELENGTH = 5.0
SEED = 20260908

#: Big enough that the GP solve is real work and small enough that a row costs
#: about a second on the contract path. ``tests/scaling`` is where size is the
#: variable; here the *route* is.
SIZE = 400

GRID = np.linspace(1.0, 20.0, SIZE)


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so a row never names a library."""

    name: str
    module: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns x64 on.
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no backend with a registered realisation is installed here",
)


def _observed() -> Spectrum:
    rng = np.random.default_rng(3)
    values = 2.0 * (GRID / REFERENCE_WAVELENGTH) ** -1.2 + rng.normal(0.0, 0.05, SIZE)
    return Spectrum(GRID * u.micron, values * u.Jy, uncertainty=np.full(SIZE, 0.05) * u.Jy)


OBSERVED = _observed()


def _problem(kit: Kit, solver_name: str) -> FittingProblem:
    module = kit.module
    return FittingProblem(
        module.PowerLaw(
            GRID,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [
            Dataset(
                OBSERVED,
                likelihood=Likelihood(
                    GaussianFamily(),
                    module.GaussianProcessNoise(
                        module.Matern32(st.lognorm(0.5, scale=0.1), st.lognorm(0.5, scale=2.0)),
                        getattr(module, solver_name)(),
                    ),
                ),
            )
        ],
        seed=SEED,
    )


def _proposal(kit: Kit, solver_name: str, *, realised: bool) -> Callable[[], float]:
    """One scored proposal, on one route, warmed up so compilation is not timed."""
    problem = _problem(kit, solver_name)
    engine = EmceeEngine(problem, walkers=8, use_realisation=realised)
    assert engine._cache.realised is realised, (
        f"the {kit.name} backend did not take the "
        f"{'realised' if realised else 'contract'} route; the row would measure the other one"
    )
    theta = problem.parameters.pack(problem.reference_values)

    def call() -> float:
        return engine.log_prob(theta)

    call()
    return call


def _run(benchmark: Any, call: Callable[[], float], group: str, name: str) -> None:
    """Time *call*, and refuse a fast wrong answer."""
    benchmark.group = group
    benchmark.name = name
    value = benchmark(call)
    assert math.isfinite(value), f"{name} returned {value!r}"


@pytest.mark.parametrize("kit", KITS, ids=[kit.name for kit in KITS])
@pytest.mark.parametrize("solver", ["QuasisepGP", "DenseGP"])
@pytest.mark.parametrize("realised", [True, False], ids=["realised", "contract"])
def test_one_proposal(benchmark: Any, kit: Kit, solver: str, realised: bool) -> None:
    """``Engine.log_prob`` once, on one backend, one solver and one route."""
    _run(
        benchmark,
        _proposal(kit, solver, realised=realised),
        group=f"engine-log-prob-{solver}-n{SIZE}",
        name=f"{kit.name} {solver} {'realised' if realised else 'contract'}",
    )


def _median_seconds(call: Callable[[], float], rounds: int = 25) -> float:
    """A median rather than a mean: one stalled round must not decide a ratio."""
    timings = []
    for _ in range(rounds):
        start = time.perf_counter()
        call()
        timings.append(time.perf_counter() - start)
    return float(np.median(timings))


#: How much faster the fast path must be before this suite believes it is one.
#: The margin measured while writing it was about 40x on a jax ``QuasisepGP``
#: problem of this size; three is the floor a shared runner cannot plausibly
#: cross by accident, and crossing it would mean the route had stopped working
#: rather than that the machine was busy.
MINIMUM_SPEED_UP = 3.0


@pytest.mark.skipif(
    not any(kit.name == "jax" for kit in KITS),
    reason="the fast path's acceptance measurement is stated on jax (W2.5 slice 3)",
)
def test_the_fast_path_beats_the_contract_path() -> None:
    """W2.5 slice 3's acceptance criterion, measured rather than asserted.

    A ratio of two medians taken in one process, so the runner's speed divides
    out; see this module's docstring for why this suite tolerates one assertion
    of that shape when it tolerates no assertion about an absolute time.
    """
    kit = next(found for found in KITS if found.name == "jax")
    contract = _median_seconds(_proposal(kit, "QuasisepGP", realised=False))
    fast = _median_seconds(_proposal(kit, "QuasisepGP", realised=True))
    ratio = contract / fast
    print(
        f"\nengine.log_prob on a jax QuasisepGP problem of {SIZE} points: "
        f"contract {contract * 1e3:.2f} ms, realised {fast * 1e3:.2f} ms, "
        f"speed-up {ratio:.1f}x"
    )
    assert ratio > MINIMUM_SPEED_UP, (
        f"the realised route was only {ratio:.1f}x the contract path "
        f"({fast * 1e3:.2f} ms against {contract * 1e3:.2f} ms); the fast path exists "
        f"because the contract path costs jax's per-operation dispatch on every "
        f"proposal, so a margin this small means it is not being taken."
    )
