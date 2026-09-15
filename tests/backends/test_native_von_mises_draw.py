"""W5.2 item 6: the native von_mises twin, on both backends.

``ampere.core.VonMisesFamily.sample`` (W4.1) drew ``rng.vonmises(mean, kappa)``
on the numpy contract path, but ``von_mises`` was missing from both backends'
``_TWINNED_FAMILIES`` — a closure-phase dataset's draws fell back to the numpy
path silently rather than being drawn natively, unlike its density (already
transcribed at W4.3). This file is the twin: it holds each backend's native
draw to circular agreement with the numpy family's own ``sample``, "circular"
because a closure phase's draws live on the circle and an ordinary mean would
average -pi and +pi to zero instead of to pi.

Compared *distributionally*, not draw for draw — exactly as the other twinned
families are (``test_torch_lowering.py``'s ``test_the_poisson_twin_draws_counts``
and friends): each backend's own random stream (torch's ``Generator``,
numpyro's ``VonMises`` under a jax key) is not numpy's, so a bit-identical
draw is not the claim. What is checked is that the two populations describe
the same circular distribution, at a sample size where the Monte Carlo error
of a circular mean/variance estimate (``~ 1/sqrt(N * kappa)`` for the
concentrated ``kappa`` a closure phase carries) is far smaller than the
tolerance.

Built on ``interferometry_fixtures.two_dataset_problem``'s ``"t3"`` dataset --
the same ``VonMisesFamily`` + ``IndependentNoise`` closure-phase declaration
``tests/inference/test_interferometry.py`` fits under NUTS -- rather than a
bespoke model, so this exercises the real native lowering path (instrument
chain, retained mask, ``_sigma``) instead of a hand-rolled shortcut.
"""

from __future__ import annotations

import dataclasses
import importlib
from typing import Any

import numpy as np
import pytest

import interferometry_fixtures as kit_data


@dataclasses.dataclass(frozen=True)
class Kit:
    name: str
    backend: Any
    itf: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            backend = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            backend.configure_x64()
        found.append(
            Kit(name, backend, importlib.import_module(f"ampere.backends.{name}.interferometry"))
        )
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS, reason="no modern backend installed; the native von_mises twin needs one"
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


#: Draws per population. Large enough that the circular mean's Monte Carlo
#: error (``~1/sqrt(N * kappa)``, ``kappa = 1/SIGMA_CLOSURE**2 = 400`` here) is
#: two orders of magnitude below :data:`TOLERANCE`.
DRAWS = 4000

#: Loose enough that a correct twin passes essentially always (the expected
#: error at DRAWS draws is ~1.3e-3 rad) and far tighter than the error a wrong
#: kappa, a wrong sign, or an unwrapped mean would produce.
TOLERANCE = 0.02


def _circular_mean_and_variance(angles: np.ndarray) -> tuple[float, float]:
    """The circular statistics a linear mean/variance would get wrong on angles."""
    resultant = np.mean(np.exp(1j * angles))
    return float(np.angle(resultant)), float(1.0 - np.abs(resultant))


def _numpy_population(problem: Any, theta_row: np.ndarray, triangle: int, draws: int) -> np.ndarray:
    """*draws* i.i.d. closure phases at *theta_row*, through the contract path."""
    values = problem.parameters.unpack(theta_row)
    return np.array(
        [
            float(
                np.asarray(
                    problem.simulate(values, observe=True, stream=f"von_mises_oracle_{index}")
                    .observations["t3"]
                    .values
                )[triangle]
            )
            for index in range(draws)
        ]
    )


def _native_population(
    problem: Any, lowered: Any, theta_row: np.ndarray, triangle: int, draws: int
) -> np.ndarray:
    """*draws* i.i.d. closure phases at the same *theta_row*, natively."""
    theta_batch = np.tile(theta_row, (draws, 1))
    predicted = lowered.simulate_batched(theta_batch).predicted
    seeds = list(range(10_000, 10_000 + draws))
    drawn = lowered.sample_observations(theta_batch, predicted, seeds)["t3"]
    return np.asarray(drawn)[:, triangle]


class TestTheNativeVonMisesTwin:
    def test_the_family_is_twinned_and_no_longer_falls_back(self, kit: Kit) -> None:
        """The membership this item adds, checked directly before the statistics are."""
        module = importlib.import_module(f"ampere.backends.{kit.name}.problem")
        assert "von_mises" in module._TWINNED_FAMILIES

        from ampere.core import VonMisesFamily

        assert module._TWINNED_FAMILIES["von_mises"] is VonMisesFamily

    def test_the_native_draw_agrees_with_the_numpy_family_circularly(self, kit: Kit) -> None:
        problem = kit_data.two_dataset_problem(kit.backend, kit.itf)
        lowered = kit.backend.lower_problem(problem)
        theta_row = problem.parameters.pack(problem.reference_values)
        triangle = 0

        native = _native_population(problem, lowered, theta_row, triangle, DRAWS)
        oracle = _numpy_population(problem, theta_row, triangle, DRAWS)

        native_mean, native_variance = _circular_mean_and_variance(native)
        oracle_mean, oracle_variance = _circular_mean_and_variance(oracle)

        # The mean angle is compared through its own wrapped residual: a mean
        # near +-pi could otherwise fail a naive subtraction that is actually
        # a difference of a few hundredths of a radian across the branch cut.
        wrapped_gap = float(np.angle(np.exp(1j * (native_mean - oracle_mean))))
        assert abs(wrapped_gap) < TOLERANCE
        assert abs(native_variance - oracle_variance) < TOLERANCE

    def test_the_draw_is_a_pure_function_of_its_seed(self, kit: Kit) -> None:
        """The property every twinned family's native draw keeps (W3.14), for von_mises."""
        problem = kit_data.two_dataset_problem(kit.backend, kit.itf)
        lowered = kit.backend.lower_problem(problem)
        theta_row = problem.parameters.pack(problem.reference_values)
        theta_batch = np.tile(theta_row, (4, 1))
        predicted = lowered.simulate_batched(theta_batch).predicted
        first = lowered.sample_observations(theta_batch, predicted, [11, 22, 33, 44])["t3"]
        second = lowered.sample_observations(theta_batch, predicted, [11, 22, 33, 44])["t3"]
        shuffled = lowered.sample_observations(theta_batch, predicted, [11, 22, 33, 45])["t3"]
        assert np.array_equal(first, second)
        assert np.array_equal(first[:3], shuffled[:3])
        assert not np.array_equal(first[3], shuffled[3])

    def test_every_draw_stays_in_the_wrapped_range(self, kit: Kit) -> None:
        """The native twin's own convention, matching ``VonMisesFamily.sample``'s ``(-pi, pi]``."""
        problem = kit_data.two_dataset_problem(kit.backend, kit.itf)
        lowered = kit.backend.lower_problem(problem)
        theta_row = problem.parameters.pack(problem.reference_values)
        theta_batch = np.tile(theta_row, (256, 1))
        predicted = lowered.simulate_batched(theta_batch).predicted
        drawn = lowered.sample_observations(theta_batch, predicted, list(range(256)))["t3"]
        assert bool(np.all(drawn > -np.pi - 1e-9)) and bool(np.all(drawn <= np.pi + 1e-9))
