"""GPU smoke tests for ``ampere.backends.torch``: the API, not the arithmetic.

What this module is for, and what it deliberately is not
--------------------------------------------------------
W2.4 slice 3 item 2, ruled 2026-09-08. These rows prove that a whole ampere
problem can be **composed, realised, differentiated and sampled on a GPU** —
that the ``device=`` keyword reaches every tensor, that
:func:`~ampere.core.declared_capabilities` composes a problem whose parts all
say ``"cuda"``, and that nothing in the realised path materialises a CPU
constant halfway through and dies on a device mismatch. That is an **API**
claim, and it is the one that can only be made by running the code.

It is *not* a numerical claim. Peter's ruling: the CPU/GPU parity of torch's
own linear algebra is torch's to guarantee, not ampere's to re-measure, so
these rows assert that the numbers are finite, that the gradient is non-zero
and that the sampler produces the draws it was asked for — and leave "the same
answer to 1e-9" to the CPU conformance battery, which runs everywhere.

**CI stays CPU-only.** Every row here skips without
:func:`torch.cuda.is_available`, and the environment ampere's gates run in
installs the CPU wheel (W2.11), so the module is *collected and skipped* on the
development machine and on CI alike. That is deliberate: a GPU suite that
cannot even be imported without a GPU rots silently, and one that is collected
every time cannot.

Running them
------------
On a machine with a CUDA build of torch::

    pytest tests/gpu

``tests/gpu`` is not in the ``test-all`` pixi task, for the same reason
``tests/scaling`` is not: it is a suite whose whole point is a resource the
gate does not have.
"""

from __future__ import annotations

import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

torch = pytest.importorskip("torch", reason="the torch backend is an optional extra")

from ampere.backends.torch import (  # noqa: E402  (after the importorskip, deliberately)
    CalibrationScale,
    DenseGP,
    GaussianProcessNoise,
    IndependentNoise,
    Matern32,
    PowerLaw,
    Resample,
    lower_problem,
)
from ampere.core import (  # noqa: E402
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    Spectrum,
    Tie,
    realise,
)

#: The one device string every part of every problem here is built with. A
#: bare ``"cuda"`` rather than ``"cuda:0"``, and the two are *different*
#: declarations (see ``ampere.backends.torch._config``), so this constant is
#: what keeps the parts of one problem agreeing.
DEVICE = "cuda"

pytestmark = pytest.mark.skipif(
    not torch.cuda.is_available(),
    reason="no CUDA device (the gate's torch is the CPU wheel, W2.11); CI stays CPU-only",
)

REFERENCE_WAVELENGTH = 1.0
FINE = np.geomspace(1.0, 10.0, 24)
COARSE = np.geomspace(1.2, 8.0, 9)
TRUTH = {"norm": 2.0, "index": -1.2}
SEED = 20260908


def _power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def _noisy(grid: np.ndarray, sigma: float, seed: int) -> Spectrum:
    """One observed spectrum. Ordinary numpy: the *data* never move to a GPU.

    A container is backend-neutral (``results_schema.md``), and the realisation
    is what puts its values on a device — once, at composition. So the data are
    built here exactly as they are in ``tests/inference/test_nuts.py``.
    """
    rng = np.random.default_rng(seed)
    values = _power_law(grid, TRUTH["norm"], TRUTH["index"])
    return Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, sigma, values.size)) * u.Jy,
        uncertainty=np.full(values.size, sigma) * u.Jy,
    )


BLUE_DATA = _noisy(FINE, 0.08, seed=11)
RED_DATA = _noisy(COARSE, 0.05, seed=13)


def _start(problem: FittingProblem) -> torch.Tensor:
    """The problem's own reference point, unconstrained, as a device tensor.

    The prior median rather than the unconstrained origin. They are not the
    same point for a ``loguniform`` hyperparameter — ``default_bijection_for``
    gives it ``Log(lower=0)``, whose origin is the constrained value 1.0 and
    may sit outside the prior's support (W2.5 slice 2's carried finding) — and
    a row asserting the density is *finite* has to start somewhere it is.
    """
    return torch.as_tensor(
        problem.unconstrain(problem.reference_values), dtype=torch.float64, device=DEVICE
    )


def joint_problem(device: str = DEVICE) -> FittingProblem:
    """``inference.md`` §15's toy joint problem, composed on *device*.

    The same declaration ``tests/inference/test_nuts.py`` builds on the CPU —
    one power law on two channels, two different instrument chains, the
    calibration factors tied — with ``device=`` threaded through every piece.
    Threaded rather than defaulted, because ``architecture.md`` §5 forbids
    ampere detecting a device: if any one of these five keywords were dropped
    the problem would declare two devices and
    :func:`~ampere.core.declared_capabilities` would refuse it, which is the
    behaviour :class:`TestTheDeviceRuleComposes` pins.
    """
    model = PowerLaw(
        FINE,
        norm=st.lognorm(0.4, scale=2.0),
        index=st.norm(-1.2, 0.3),
        reference_wavelength=REFERENCE_WAVELENGTH,
        channels=("blue", "red"),
        device=device,
    )
    return FittingProblem(
        model,
        DatasetCollection(
            {
                "blue": Dataset(
                    BLUE_DATA,
                    Instrument(
                        [CalibrationScale(st.lognorm(0.05), label="calibration", device=device)],
                        channel="blue",
                    ),
                    Likelihood(GaussianFamily(), IndependentNoise(device=device)),
                ),
                "red": Dataset(
                    RED_DATA,
                    Instrument(
                        [
                            Resample(COARSE, device=device),
                            CalibrationScale(st.lognorm(0.05), label="calibration", device=device),
                        ],
                        channel="red",
                    ),
                    Likelihood(GaussianFamily(), IndependentNoise(device=device)),
                ),
            }
        ),
        ties=[
            Tie(
                "calibration",
                ("blue.instrument.calibration.scale", "red.instrument.calibration.scale"),
            )
        ],
        seed=SEED,
    )


def gp_problem(device: str = DEVICE) -> FittingProblem:
    """One dataset under a dense-GP noise model, every piece on *device*.

    The GP is the case worth a separate row: it is the one that reaches a
    Cholesky, it is where ``architecture.md`` §5's float64-on-the-GPU policy
    actually costs something, and it is the composition with the most pieces
    to get onto one device — a kernel, a solver and a noise model, each taking
    the keyword separately.
    """
    return FittingProblem(
        PowerLaw(
            FINE,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
            device=device,
        ),
        [
            Dataset(
                BLUE_DATA,
                likelihood=Likelihood(
                    GaussianFamily(),
                    GaussianProcessNoise(
                        Matern32(st.loguniform(0.01, 1.0), st.loguniform(0.2, 5.0), device=device),
                        DenseGP().configured(device=device),
                        device=device,
                    ),
                ),
            )
        ],
        seed=SEED,
    )


class TestTheDeviceRuleComposes:
    """``declared_capabilities``' device rule, on a machine that has two devices."""

    def test_every_part_declares_the_device_it_was_built_on(self) -> None:
        problem = joint_problem()
        assert problem.device == DEVICE
        parts = (*problem.models.values(), *problem.datasets.capability_parts)
        assert {str(part.DEVICE) for part in parts} == {DEVICE}

    def test_a_gp_problems_kernel_solver_and_noise_agree(self) -> None:
        problem = gp_problem()
        assert problem.device == DEVICE
        noise = problem.datasets["default"].likelihood.noise
        assert noise.DEVICE == DEVICE
        assert noise.kernel.DEVICE == DEVICE
        assert noise.solver.DEVICE == DEVICE
        assert noise.solver.provenance_config()["device"] == DEVICE

    def test_a_cpu_part_in_a_gpu_problem_is_refused_at_composition(self) -> None:
        """The whole reason the flag is per-instance: the mistake is loud.

        A model on the GPU beside a CPU calibration step is a configuration
        error that would otherwise surface as a torch device mismatch inside
        the first evaluation — or, worse, not at all, since torch broadcasts a
        CPU *scalar* against a CUDA tensor without complaint.
        """
        from ampere.core.exceptions import DatasetError

        with pytest.raises(DatasetError, match="different devices"):
            FittingProblem(
                PowerLaw(FINE, norm=st.lognorm(0.4, scale=2.0), index=-1.2, device=DEVICE),
                [
                    Dataset(
                        BLUE_DATA,
                        Instrument([CalibrationScale(st.lognorm(0.05), label="calibration")]),
                        Likelihood(GaussianFamily(), IndependentNoise(device=DEVICE)),
                    )
                ],
                seed=SEED,
            )


class TestTheRealisedDensityRunsOnTheGpu:
    """Compose, realise, evaluate, differentiate — all on the device."""

    def test_the_realisation_reports_the_device_and_lives_there(self) -> None:
        lowered = lower_problem(joint_problem())
        assert str(lowered.device) == DEVICE
        for dataset in lowered._datasets:
            assert dataset.observed_values.device.type == "cuda"
            assert dataset.observed_coordinates.device.type == "cuda"

    @pytest.mark.parametrize("build", [joint_problem, gp_problem], ids=["joint", "gp"])
    def test_the_density_is_finite_and_on_the_device(self, build: Any) -> None:
        problem = build()
        lowered = realise(problem)
        value = lowered.log_prob_unconstrained(_start(problem))
        assert value.device.type == "cuda"
        assert torch.isfinite(value)

    @pytest.mark.parametrize("build", [joint_problem, gp_problem], ids=["joint", "gp"])
    def test_the_gradient_exists_and_is_not_zero(self, build: Any) -> None:
        """The claim ``DIFFERENTIABLE = True`` makes, taken on the GPU.

        Non-zero rather than merely finite: a gradient that is exactly zero in
        every coordinate is what a graph cut at a device boundary produces, and
        it is the failure this row exists to catch.
        """
        problem = build()
        lowered = realise(problem)
        y = _start(problem).requires_grad_(True)
        lowered.log_prob_unconstrained(y).backward()
        assert y.grad is not None
        assert y.grad.device.type == "cuda"
        assert torch.all(torch.isfinite(y.grad))
        assert float(torch.max(torch.abs(y.grad))) > 0.0

    def test_a_stack_is_evaluated_in_one_vmap_call(self) -> None:
        """W2.4 slice 2's batched surface, on the device slice 3 added."""
        problem = joint_problem()
        lowered = realise(problem)
        stack = _start(problem).reshape(1, -1).repeat(4, 1)
        batched = lowered.log_prob_unconstrained_batched(stack)
        assert batched.shape == (4,)
        assert batched.device.type == "cuda"
        one_at_a_time = torch.stack([lowered.log_prob_unconstrained(row) for row in stack])
        assert torch.allclose(batched, one_at_a_time)


class TestNutsDraws:
    """A handful of draws, through the registered realisation, on the GPU.

    Deliberately a *handful*: the ruling is that these rows prove the API
    works, so the budget is the smallest one that exercises warmup, the mass
    matrix and the emission of a run. Recovering a posterior is
    ``tests/inference/test_nuts.py``'s job and it does it on the CPU.
    """

    def test_it_produces_the_draws_it_was_asked_for(self) -> None:
        from ampere.inference import NUTSEngine

        problem = joint_problem()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = NUTSEngine(problem).run(draws=10, warmup=10, chains=1)
        assert run.attrs["ampere_backend"] == "torch"
        assert bool(run.attrs["ampere_realised"]) is True
        posterior = run["posterior"].ds
        assert posterior.sizes["draw"] == 10
        assert posterior.sizes["chain"] == 1
        for name in problem.parameters.free_names:
            assert np.all(np.isfinite(np.asarray(posterior[name].values)))
