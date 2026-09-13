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
from ampere.backends.jax import IndependentNoise
from ampere.core import (
    ClosurePhases,
    ComplexGaussianFamily,
    Dataset,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    Spectrum,
    VisibilitySet,
    VonMisesFamily,
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


class TestChunkSharding:
    """W3.1 slice 2's multi-device hook, at the API level.

    Peter's addendum of 2026-09-09 asks for the sharding API to be designed,
    documented and smoke-tested here, with the full exercise left to the GPU
    item — so these rows say what a sharder must be true of and run them on
    whatever accelerators this process reports. The single-device degenerate
    case is exercised on CPU by ``tests/backends``; what needs hardware is the
    claim that splitting a chunk across devices does not change the answer,
    which is the *whole* of a sharder's contract.
    """

    def problem(self) -> FittingProblem:
        return _problem(DenseGP, PLATFORM)

    def theta(self, problem: FittingProblem, draws: int = 8) -> np.ndarray:
        rng = np.random.default_rng(20260909)
        return np.stack(
            [problem.prior_transform(row) for row in rng.uniform(0.05, 0.95, (draws, 2))]
        )

    def test_the_sharder_reports_the_devices_it_will_use(self) -> None:
        from ampere.backends.jax import MeshSharder

        sharder = MeshSharder(ACCELERATORS)
        assert len(sharder.devices()) == len(ACCELERATORS)

    def test_sharding_does_not_change_the_prediction(self) -> None:
        """The value contract: where it was computed is not part of the answer."""
        from ampere.backends.jax import MeshSharder, SingleDeviceSharder

        problem = self.problem()
        lowered = lower_problem(problem)
        theta = self.theta(problem)
        alone = lowered.simulate_batched(theta, sharder=SingleDeviceSharder(PLATFORM))
        spread = lowered.simulate_batched(theta, sharder=MeshSharder(ACCELERATORS))
        for label, values in alone.predicted.items():
            assert np.allclose(values, spread.predicted[label], rtol=0.0, atol=1e-9)

    def test_a_ragged_chunk_is_padded_not_refused(self) -> None:
        """A budget has no reason to be a multiple of the device count."""
        from ampere.backends.jax import MeshSharder

        problem = self.problem()
        lowered = lower_problem(problem)
        theta = self.theta(problem, draws=len(ACCELERATORS) * 2 + 1)
        spread = lowered.simulate_batched(theta, sharder=MeshSharder(ACCELERATORS))
        assert len(spread) == theta.shape[0]


# ---------------------------------------------------------------------------
# The Phase 4 modality's placement (W4.3)
# ---------------------------------------------------------------------------

#: A two-baseline, one-triangle array, inline and tiny. Self-contained on
#: purpose: the *arithmetic* of the interferometric twins is held to the
#: reference backend on the CPU by ``tests/backends/test_native_interferometry.py``
#: and by the conformance battery, so what is left for a GPU row is the
#: **placement** — and a placement row should not also carry a synthetic source,
#: a negotiated grid and a noise realisation it cannot check.
_UV = np.array([1.0e7, -1.6e7, 0.6e7])
_VV = np.array([0.4e7, 1.2e7, -1.6e7])
_WAVES = np.full(3, 2.2)
_FIELD_OF_VIEW = 40.0


def _visibilities() -> VisibilitySet:
    """Three baselines with a real error bar apiece."""
    return VisibilitySet(
        _UV,
        _VV,
        _WAVES * u.micron,
        np.array([0.9 + 0.1j, 0.6 - 0.2j, 0.8 + 0.05j]) * u.Jy,
        uncertainty=np.full(3, 0.02) * u.Jy,
    )


def _closure_phases() -> ClosurePhases:
    """One triangle, whose third baseline is the negated sum of the stored two."""
    return ClosurePhases(
        _UV[:1],
        _VV[:1],
        _UV[1:2],
        _VV[1:2],
        _WAVES[:1] * u.micron,
        np.array([-0.4]) * u.rad,
        uncertainty=np.full(1, 0.05) * u.rad,
    )


def interferometric_problem(device: str = PLATFORM) -> FittingProblem:
    """Visibilities and closure phases from one sky, every piece on *device*."""
    from ampere.backends.jax import interferometry as itf

    visibility = _visibilities()
    phase = _closure_phases()
    model = itf.Binary.on_field(
        _FIELD_OF_VIEW * u.mas,
        8,
        channels="sky",
        separation=st.uniform(6.0, 14.0),
        position_angle=0.7,
        flux_ratio=st.uniform(0.1, 0.7),
        flux=1.7,
        component_fwhm=2.0,
        device=device,
    )
    return FittingProblem(
        model,
        [
            Dataset(
                visibility,
                Instrument(
                    [
                        itf.FourierSample.from_observed(
                            visibility, field_of_view=_FIELD_OF_VIEW * u.mas, device=device
                        )
                    ],
                    channel="sky",
                    label="vis",
                ),
                Likelihood(ComplexGaussianFamily(), IndependentNoise(device=device)),
                label="vis",
            ),
            Dataset(
                phase,
                Instrument(
                    [
                        itf.FourierSample.from_observed(
                            phase, field_of_view=_FIELD_OF_VIEW * u.mas, device=device
                        ),
                        itf.ClosurePhase(device=device),
                    ],
                    channel="sky",
                    label="t3",
                ),
                Likelihood(VonMisesFamily(), IndependentNoise(device=device)),
                label="t3",
            ),
        ],
        seed=20260913,
    )


class TestTheInterferometricTwinsComposeOnTheDevice:
    """``declared_capabilities``' device rule, for Phase 4's proof modality.

    The arithmetic is held to the reference backend on the CPU by
    ``tests/backends/test_native_interferometry.py`` and by the conformance
    battery; what is specific to an accelerator is the **placement**, and that a
    density built from a discrete Fourier transform of a model image evaluates
    and differentiates there.
    """

    def test_every_interferometric_part_declares_the_device(self) -> None:
        problem = interferometric_problem()
        assert problem.device == PLATFORM
        parts = (*problem.models.values(), *problem.datasets.capability_parts)
        assert {str(part.DEVICE) for part in parts} == {PLATFORM}

    def test_the_realised_density_and_its_gradient_live_on_the_device(self) -> None:
        problem = interferometric_problem()
        lowered = lower_problem(problem)
        theta = jax.numpy.asarray(problem.unconstrain(problem.reference_values))
        value = lowered.log_prob_unconstrained(theta)
        assert {device.platform for device in value.devices()} == {PLATFORM}
        assert bool(jax.numpy.isfinite(value))
        gradient = jax.grad(lowered.log_prob_unconstrained)(theta)
        assert {device.platform for device in gradient.devices()} == {PLATFORM}
        assert bool(jax.numpy.all(jax.numpy.isfinite(gradient)))
