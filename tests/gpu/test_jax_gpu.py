"""GPU smoke tests for the jax backend — skipped wherever there is no accelerator.

Run them with::

    pixi run -e jax gpu

and expect a skip on any ordinary machine and in CI, which is CPU-only by
ruling (``DEVELOPMENT_PLAN.md`` §2's CI row, enforced since W2.11 by pinning
CPU wheels). That is the point: ``architecture.md`` §5 says ampere **never**
auto-detects a device, so a suite that silently found one and ran on it would
be evidence against the rule it is meant to check.

What these rows claim, and what they deliberately do not
--------------------------------------------------------
They claim that the ``device=`` keyword slice 3 added to every piece of this
backend is *plumbing that works*: the arrays a piece owns really move, a
problem composed entirely on the accelerator composes, and the density it
computes is the CPU's answer to the tolerance a change of device is entitled
to. They claim nothing about speed — a smoke test on unknown hardware cannot —
and nothing about float32, which stays the recorded per-run opt-out it is on
the CPU.

``jax.devices()`` is asked once, by name. A machine with a GPU present but only
a CPU jaxlib installed reports the CPU alone, which is exactly the state this
repository's ``jax`` environment is in, and these rows skip there too — the
absence of an *accelerator* is what they test for, not the absence of
hardware.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

pytest.importorskip("jax")
pytest.importorskip("numpyro")

import astropy.units as u
import jax
import scipy.stats as st

from ampere.backends.jax import (
    CalibrationScale,
    DenseGP,
    GaussianProcessNoise,
    Matern32,
    PowerLaw,
    QuasisepGP,
    configure_x64,
    lower_problem,
)
from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)

configure_x64()


def _accelerators() -> list[Any]:
    """Every non-CPU device this process has. Asked once, never guessed at."""
    return [device for device in jax.devices() if device.platform != "cpu"]


ACCELERATORS = _accelerators()

pytestmark = pytest.mark.skipif(
    not ACCELERATORS,
    reason=(
        "no jax accelerator in this process (jax.devices() reports "
        f"{sorted({d.platform for d in jax.devices()})}); ampere never auto-detects one, "
        "so there is nothing for these rows to place anything on."
    ),
)

#: The platform name the rows ask for. One per run: ampere places by platform
#: and leaves the choice among several of a kind to the caller.
PLATFORM = ACCELERATORS[0].platform if ACCELERATORS else "cpu"

GRID = np.linspace(1.0, 12.0, 48)
REFERENCE = 5.0


def _observed() -> Spectrum:
    truth = 2.0 * (GRID / REFERENCE) ** -1.2
    values = truth + np.random.default_rng(20260908).normal(0.0, 0.05, GRID.size)
    return Spectrum(GRID * u.micron, values * u.Jy, uncertainty=np.full(GRID.size, 0.05) * u.Jy)


OBSERVED = _observed()


def _problem(solver_factory: Any, device: str) -> FittingProblem:
    """One flexible-likelihood problem, every piece of it on *device*."""
    return FittingProblem(
        PowerLaw(
            GRID,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE,
            device=device,
        ),
        [
            Dataset(
                OBSERVED,
                likelihood=Likelihood(
                    GaussianFamily(),
                    GaussianProcessNoise(
                        Matern32(
                            st.lognorm(0.5, scale=0.1),
                            st.lognorm(0.5, scale=2.0),
                            device=device,
                        ),
                        solver_factory(device=device),
                        device=device,
                    ),
                ),
            )
        ],
        seed=20260908,
    )


class TestPlacement:
    """The arrays a piece owns are where it says they are."""

    def test_a_model_grid_lives_on_the_accelerator(self) -> None:
        model = PowerLaw(GRID, device=PLATFORM)
        assert model.DEVICE == PLATFORM
        assert {d.platform for d in model.grid("default").devices()} == {PLATFORM}

    def test_a_step_influence_matrix_lives_on_the_accelerator(self) -> None:
        step = CalibrationScale(1.0, device=PLATFORM)
        assert step.DEVICE == PLATFORM

    def test_a_kernel_covariance_lives_on_the_accelerator(self) -> None:
        kernel = Matern32(0.4, 2.0, device=PLATFORM)
        covariance = kernel.matrix(GRID, GRID, {})
        assert {d.platform for d in covariance.devices()} == {PLATFORM}

    def test_a_solver_places_what_it_is_handed(self) -> None:
        solver = DenseGP(device=PLATFORM)
        assert solver.DEVICE == PLATFORM
        assert dict(solver.provenance_config())["device"] == PLATFORM
        assert {d.platform for d in solver.place(np.ones(4)).devices()} == {PLATFORM}


class TestTheDensityAgrees:
    """The same declaration, on two devices, is the same posterior."""

    @pytest.mark.parametrize("factory", [DenseGP, QuasisepGP], ids=["dense", "quasisep"])
    def test_the_realised_density_matches_the_cpu(self, factory: Any) -> None:
        on_cpu = lower_problem(_problem(factory, "cpu"))
        on_device = lower_problem(_problem(factory, PLATFORM))
        theta = np.random.default_rng(4).normal(0.0, 0.5, on_cpu.free_size)
        expected = float(np.asarray(on_cpu.log_prob_unconstrained(theta)))
        got = float(np.asarray(on_device.log_prob_unconstrained(theta)))
        # A change of device is entitled to a different summation order and a
        # different fused kernel; it is not entitled to a different answer.
        assert got == pytest.approx(expected, rel=1e-9, abs=1e-9)

    def test_the_contract_path_matches_the_cpu(self) -> None:
        on_cpu = _problem(DenseGP, "cpu")
        on_device = _problem(DenseGP, PLATFORM)
        theta = on_cpu.parameters.pack(on_cpu.reference_values)
        assert on_device.log_prob(theta) == pytest.approx(on_cpu.log_prob(theta), rel=1e-9)


class TestTheRefusalsStillHold:
    """Having an accelerator does not license using it by accident."""

    def test_the_default_is_still_the_cpu(self) -> None:
        assert PowerLaw(GRID).DEVICE == "cpu"
        assert Matern32(0.4, 2.0).DEVICE == "cpu"
        assert DenseGP().DEVICE == "cpu"

    def test_a_mixed_device_problem_is_still_refused(self) -> None:
        from ampere.core.exceptions import DatasetError

        with pytest.raises(DatasetError, match="different devices"):
            FittingProblem(
                PowerLaw(GRID, norm=st.lognorm(0.4, scale=2.0), device="cpu"),
                [
                    Dataset(
                        OBSERVED,
                        likelihood=Likelihood(
                            GaussianFamily(),
                            GaussianProcessNoise(
                                Matern32(0.4, 2.0, device=PLATFORM),
                                DenseGP(device=PLATFORM),
                                device=PLATFORM,
                            ),
                        ),
                    )
                ],
                seed=1,
            )
