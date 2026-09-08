"""The gradient-free fast path: scoring through the backend's realisation.

W2.5 slice 3, and the reason it exists is a measurement rather than a taste.
W2.10 timed the jax **contract** path at about 28 ms flat per ``log_prob``,
whatever the problem, because ``FittingProblem.evaluate`` is numpy and a jax
model evaluated an operation at a time pays jax's dispatch on every one of
them. That made emcee, dynesty and zeus *slower* on the backend built for
speed, which is not a trade-off anybody chose.

``inference.md`` §10a already had the answer in it: a realisation may offer
``log_likelihood_terms``, and W2.4 slice 2 taught :meth:`Engine.finish` to
consume it for the per-draw record. Slice 3 lets
:class:`ampere.inference.engine._EvaluationCache` *score* through it as well —
the prior from the declaration, the per-dataset terms from the realisation,
the same ``Evaluation`` out.

What this module holds the path to
-----------------------------------
1. **The same numbers.** Every row that matters here compares the fast path
   against ``problem.evaluate`` on the same θ. The numpy path is the oracle
   and it stays the oracle; ``ampere.core.realise`` already checks agreement
   at the reference point before a driver may use a realisation, and this adds
   the check at scattered points, including outside the prior's support.
2. **The same record.** ``Evaluation``'s prior/likelihood split, its
   per-dataset contributions, and its NaN-not-``-inf`` answer for a point of
   zero prior mass are all part of what ``ampere.results`` writes, so they are
   compared field by field rather than through ``log_prob`` alone.
3. **The declared opt-outs.** ``use_realisation=False`` restores the contract
   path, and ``strict=True`` keeps it whether asked or not — a realised density
   cannot raise or record a reason, and ``strict`` is the declaration that says
   reasons matter more than speed.
4. **Both backends.** The rows are parametrised over whichever of jax and
   torch this environment has, as ``test_nuts.py``'s are, and for the same
   reason: the claim is about ``ampere.inference``, which knows about neither.

The timing itself is **not** asserted here — ``tests/benchmarks`` owns that,
where a ratio measured in one process on one machine is the only honest form
of the claim.
"""

from __future__ import annotations

import dataclasses
import importlib
import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)
from ampere.inference import DynestyEngine, EmceeEngine

REFERENCE_WAVELENGTH = 1.0
SEED = 20260908
GRID = np.geomspace(1.0, 20.0, 40)


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so a test body never names a library."""

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
            # ``lowering.md`` §10.2(a): the application turns x64 on, never
            # ampere. A test suite is an application.
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no backend with a registered realisation is installed here",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


def _power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def _observed() -> Spectrum:
    rng = np.random.default_rng(11)
    values = _power_law(GRID, 2.0, -1.2) + rng.normal(0.0, 0.08, GRID.size)
    return Spectrum(GRID * u.micron, values * u.Jy, uncertainty=np.full(GRID.size, 0.08) * u.Jy)


OBSERVED = _observed()


def _problem(kit: Kit, *, gp: bool = False, strict: bool = False) -> FittingProblem:
    """One dataset, every piece of it from *kit*'s backend."""
    module = kit.module
    noise = (
        module.GaussianProcessNoise(
            module.Matern32(st.lognorm(0.5, scale=0.1), st.lognorm(0.5, scale=2.0)),
            module.DenseGP(),
        )
        if gp
        else module.IndependentNoise()
    )
    return FittingProblem(
        module.PowerLaw(
            GRID,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(OBSERVED, likelihood=Likelihood(GaussianFamily(), noise))],
        seed=SEED,
        strict=strict,
    )


def _points(problem: FittingProblem, count: int = 8) -> list[np.ndarray]:
    """Prior draws plus one deliberately impossible point."""
    rng = np.random.default_rng(4)
    drawn = [problem.parameters.pack(problem.sample_prior(rng)) for _ in range(count)]
    outside = np.array(drawn[0], dtype=float)
    outside[0] = -1.0e3  # below a Log-bijected parameter's support
    return [*drawn, outside]


class TestTheFastPathIsTakenAtAll:
    """It is on by default where a realisation exists, and off where it does not."""

    def test_a_backend_problem_scores_through_its_realisation(self, kit: Kit) -> None:
        engine = EmceeEngine(_problem(kit), walkers=8)
        assert engine._cache.realised is True

    def test_it_can_be_declined(self, kit: Kit) -> None:
        engine = EmceeEngine(_problem(kit), walkers=8, use_realisation=False)
        assert engine._cache.realised is False

    def test_a_strict_problem_keeps_the_contract_path(self, kit: Kit) -> None:
        """``strict=True`` says "do not turn a failure into a number", and a
        realised density can do nothing else (§10a sub-decision 2)."""
        engine = EmceeEngine(_problem(kit, strict=True), walkers=8)
        assert engine._cache.realised is False

    def test_a_reference_problem_has_no_realisation_and_says_so(self) -> None:
        """No realisation is registered for the reference backend, and the
        fall back to the contract path is silent because it is not a fault."""
        from ampere.backends.reference import PowerLaw

        problem = FittingProblem(
            PowerLaw(
                GRID,
                norm=st.lognorm(0.4, scale=2.0),
                index=st.norm(-1.2, 0.3),
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(OBSERVED, likelihood=Likelihood(GaussianFamily()))],
            seed=SEED,
        )
        engine = EmceeEngine(problem, walkers=8)
        assert engine._cache.realised is False


class TestTheNumbersAreUnchanged:
    """The numpy path is the oracle; the fast path must not move a number."""

    @pytest.mark.parametrize("gp", [False, True], ids=["iid", "gp"])
    def test_every_field_of_the_evaluation_agrees(self, kit: Kit, gp: bool) -> None:
        problem = _problem(kit, gp=gp)
        engine = EmceeEngine(problem, walkers=8)
        assert engine._cache.realised is True
        compared = 0
        for theta in _points(problem):
            expected = problem.evaluate(theta)
            got = engine._cache.evaluate(theta)
            assert got.log_prior == pytest.approx(expected.log_prior, abs=1e-9)
            if not math.isfinite(expected.log_prior):
                # Zero prior mass: no model runs, and "not evaluated" is NaN
                # on both paths rather than -inf, which would mean impossible.
                assert math.isnan(got.log_likelihood)
                assert got.log_prob == -math.inf
                continue
            assert got.log_likelihood == pytest.approx(expected.log_likelihood, abs=1e-8)
            assert got.log_prob == pytest.approx(expected.log_prob, abs=1e-8)
            assert set(got.contributions) == set(expected.contributions)
            for label, value in expected.contributions.items():
                assert got.contributions[label] == pytest.approx(value, abs=1e-8)
            compared += 1
        assert compared > 0, "every point was outside the support; the row proved nothing"

    def test_the_two_routes_agree_with_each_other(self, kit: Kit) -> None:
        """Belt and braces: the same driver, twice, differing only in the flag."""
        problem = _problem(kit, gp=True)
        fast = EmceeEngine(problem, walkers=8)
        slow = EmceeEngine(problem, walkers=8, use_realisation=False)
        for theta in _points(problem, count=5):
            assert fast.log_prob(theta) == pytest.approx(slow.log_prob(theta), abs=1e-8)


class TestTheRunRecordsWhichPathItTook:
    """``engine_realised_evaluations`` and ``ampere_realised``."""

    def test_a_realised_run_says_so_and_counts_the_proposals(self, kit: Kit) -> None:
        problem = _problem(kit)
        engine = EmceeEngine(problem, walkers=8)
        run = engine.run(steps=8, burn_in=2)
        assert run.attrs["ampere_realised"] == 1
        assert run.attrs["ampere_engine_realised_evaluations"] > 0
        assert (
            run.attrs["ampere_engine_realised_evaluations"]
            <= run.attrs["ampere_engine_evaluations"]
        )

    def test_a_contract_path_run_still_records_zero(self, kit: Kit) -> None:
        problem = _problem(kit)
        engine = EmceeEngine(problem, walkers=8, use_realisation=False)
        run = engine.run(steps=8, burn_in=2)
        assert run.attrs["ampere_realised"] == 0
        assert run.attrs["ampere_engine_realised_evaluations"] == 0

    def test_the_posterior_is_the_same_run_either_way(self, kit: Kit) -> None:
        """Same seed, same declaration, same draws: the fast path is a change
        of arithmetic route, not of the chain the sampler walks."""
        fast = EmceeEngine(_problem(kit), walkers=8).run(steps=8, burn_in=2)
        slow = EmceeEngine(_problem(kit), walkers=8, use_realisation=False).run(steps=8, burn_in=2)
        for name in ("model.norm", "model.index"):
            assert np.asarray(fast["posterior"][name]) == pytest.approx(
                np.asarray(slow["posterior"][name]), abs=1e-8
            )


class TestTheOtherDrivers:
    """The path is the base class's, so every gradient-free driver has it."""

    def test_dynesty_scores_through_the_realisation_too(self, kit: Kit) -> None:
        pytest.importorskip("dynesty")
        problem = _problem(kit)
        engine = DynestyEngine(problem, live_points=40)
        assert engine._cache.realised is True
        run = engine.run(maxiter=120)
        assert run.attrs["ampere_realised"] == 1
        assert run.attrs["ampere_engine_realised_evaluations"] > 0

    def test_zeus_scores_through_the_realisation_too(self, kit: Kit) -> None:
        pytest.importorskip("zeus")
        from ampere.inference import ZeusEngine

        engine = ZeusEngine(_problem(kit), walkers=8)
        assert engine._cache.realised is True


class TestWhatTheFastPathCosts:
    """The one thing it does change, stated and tested rather than implied."""

    def test_a_failure_loses_its_reason_on_this_path(self, kit: Kit) -> None:
        """``inference.md`` §10a sub-decision 2: inside a realised density
        nothing can raise and no ``Failure`` can be built, so a point the model
        cannot score is a bare ``-inf``. The contract path records a reason for
        the same point; this one does not, and the run's failure summary is
        therefore silent. That is the price of the speed and the reason
        ``use_realisation=False`` and ``strict=True`` exist.
        """
        problem = _problem(kit, gp=True)
        # A kernel amplitude whose square overflows float64: the covariance is
        # non-finite, the factorisation fails, and the prior density is still
        # finite there, so the point is reachable and both paths must answer.
        values = dict(problem.reference_values)
        amplitude = next(name for name in values if name.endswith("amplitude"))
        values[amplitude] = 1.0e200
        theta = problem.parameters.pack(values)

        contract = _problem(kit, gp=True).evaluate(theta)
        assert contract.log_prob == -math.inf
        assert contract.failure is not None

        engine = EmceeEngine(problem, walkers=8)
        assert engine._cache.realised is True
        got = engine._cache.evaluate(theta)
        assert got.log_prob == -math.inf
        assert got.failure is None
        assert problem.failure_summary() == ""
