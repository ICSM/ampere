"""W2.4 slice 3: the per-instance ``device=``, exercised without a GPU.

``architecture.md`` §5 makes the device an explicit, per-instance choice —
"CPU by default; GPU strictly opt-in; never auto-detected" — and W2.4 slice 3
threads a ``device=`` keyword through every model, instrument step, kernel,
noise model and solver ``ampere.backends.torch`` ships, reports it back as the
instance's ``DEVICE`` capability flag, and moves buffers with ``.to(...)``.

The interesting claims are about a **second** device, and the gate has one CPU.
``torch.device("meta")`` is what makes them testable here anyway: a meta tensor
has a dtype, a shape and a device and *no storage*, so every piece in this
package will happily declare itself on it and none of them can then be
evaluated. That is exactly the shape of the property under test — this module
is about **declaration and composition**, and ``tests/gpu/`` is about
evaluation, which needs real memory and therefore a real accelerator.

So: ``"meta"`` stands in for ``"cuda"`` in every row that only has to reach
composition, and the two rows that must actually compute run on the CPU. The
one thing meta cannot stand in for is CUDA's arithmetic, and by Peter's ruling
(W2.4 slice 3 item 2) ampere does not re-measure torch's CPU/GPU parity anyway.

The complex-family rows are here too, for the ordinary reason a backend suite
holds a family body: :mod:`ampere.backends.torch._families` is a transcription
of ``ampere.core``'s closed form, and the smallest honest test of a
transcription is the two side by side.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

torch = pytest.importorskip("torch")

from ampere.backends.torch import (  # noqa: E402
    BlackBody,
    CalibrationScale,
    DenseGP,
    FractionalModelGPNoise,
    FractionalModelNoise,
    GaussianProcessNoise,
    IndependentNoise,
    LSFConvolution,
    Matern32,
    ModifiedBlackBody,
    PowerLaw,
    QuasisepGP,
    Resample,
    SquaredExponential,
    SyntheticPhotometry,
    lower_problem,
)
from ampere.backends.torch._config import (  # noqa: E402
    complex_dtype,
    device_name,
    resolve_device,
)
from ampere.core import (  # noqa: E402
    ComplexGaussianFamily,
    Dataset,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    NoiseParams,
    Spectrum,
    realise,
)
from ampere.core.exceptions import DatasetError, LikelihoodError, LoweringError  # noqa: E402

#: The stand-in for a second device. See the module docstring.
OTHER = "meta"

GRID = np.array([1.0, 2.0, 4.0, 8.0])
OBSERVED = Spectrum(
    GRID * u.micron,
    [2.1, 0.9, 0.55, 0.24] * u.Jy,
    uncertainty=np.full(GRID.size, 0.1) * u.Jy,
)


def a_model(**kwargs: Any) -> PowerLaw:
    return PowerLaw(GRID, norm=st.loguniform(0.1, 10.0), index=st.norm(-1.0, 0.5), **kwargs)


# ---------------------------------------------------------------------------
# 1. The keyword reaches every shipped piece
# ---------------------------------------------------------------------------


def _synthetic_photometry(device: Any) -> SyntheticPhotometry:
    """The one shipped step that cannot be built on the meta device.

    Not a gap in its ``device=`` keyword: its ``__init__`` *validates the
    response curves* — finite, non-negative, non-degenerate — which is a
    ``bool()`` of a tensor, and a meta tensor has no data to read. That is
    right for a step whose whole content is a tabulated filter set (a bad curve
    should be refused where it is declared, not discovered in a fit), and it
    means this one piece is exercised on the CPU here rather than against the
    meta stand-in. On a real accelerator it takes the keyword like any other.
    """
    return SyntheticPhotometry(
        ["W1", "W2"],
        GRID,
        np.array([[0.0, 1.0, 1.0, 0.0], [1.0, 1.0, 0.0, 0.0]]),
        detector="photon",
        device=device,
    )


def _pieces(device: Any) -> dict[str, Any]:
    """One of every shipped kind, built on *device*.

    Written as a mapping rather than a list so a failure names the class that
    dropped the keyword rather than an index. :class:`SyntheticPhotometry` is
    the one omission and has its own row; see :func:`_synthetic_photometry`.
    """
    return {
        "PowerLaw": a_model(device=device),
        "BlackBody": BlackBody(GRID, temperature=1000.0, device=device),
        "ModifiedBlackBody": ModifiedBlackBody(GRID, temperature=100.0, device=device),
        "CalibrationScale": CalibrationScale(st.lognorm(0.2), label="cal", device=device),
        "Resample": Resample([1.5, 3.0, 6.0], device=device),
        "LSFConvolution": LSFConvolution(resolving_power=100.0, device=device),
        "Matern32": Matern32(0.4, 2.0, device=device),
        "SquaredExponential": SquaredExponential(0.4, 2.0, device=device),
        "IndependentNoise": IndependentNoise(device=device),
        "GaussianProcessNoise": GaussianProcessNoise(
            Matern32(0.4, 2.0, device=device),
            DenseGP().configured(device=device),
            device=device,
        ),
        "FractionalModelNoise": FractionalModelNoise(0.05, device=device),
        "FractionalModelGPNoise": FractionalModelGPNoise(
            Matern32(0.4, 2.0, device=device),
            DenseGP().configured(device=device),
            f=0.05,
            device=device,
        ),
        "DenseGP": DenseGP().configured(device=device),
    }


class TestEveryPieceTakesTheKeyword:
    """Item 1: models, steps, kernels, noise models and solvers, all of them."""

    @pytest.mark.parametrize("name", sorted(_pieces("cpu")))
    def test_it_reports_the_device_it_was_built_on(self, name: str) -> None:
        assert _pieces(OTHER)[name].DEVICE == OTHER

    @pytest.mark.parametrize("name", sorted(_pieces("cpu")))
    def test_the_default_is_the_cpu_and_is_never_detected(self, name: str) -> None:
        """No argument means the CPU, on a machine with an accelerator or without."""
        assert _pieces("cpu")[name].DEVICE == "cpu"

    @pytest.mark.parametrize("name", sorted(_pieces("cpu")))
    def test_the_class_default_is_untouched(self, name: str) -> None:
        """An *instance* attribute shadows the ClassVar; it does not replace it.

        The mechanism matters as much as the result: if these wrote the class
        attribute instead, building one model on a GPU would silently move
        every other model in the process onto it.
        """
        piece = _pieces(OTHER)[name]
        assert type(piece).DEVICE == "cpu"
        assert piece.DEVICE == OTHER

    def test_synthetic_photometry_takes_it_too(self) -> None:
        """The step ``_pieces`` leaves out, on a device it can validate itself on."""
        assert _synthetic_photometry("cpu").DEVICE == "cpu"
        assert type(_synthetic_photometry("cpu")).DEVICE == "cpu"
        assert all(
            tensor.device.type == "cpu" for tensor in _synthetic_photometry("cpu").tensors.buffers()
        )

    def test_the_quasiseparable_solver_still_refuses_to_be_moved(self) -> None:
        """celerite2's kernels are compiled float64 CPU; the refusal is the truth."""
        with pytest.raises(LikelihoodError, match="configure"):
            QuasisepGP().configured(device=OTHER)
        assert QuasisepGP().DEVICE == "cpu"

    def test_a_solvers_provenance_records_the_device_it_was_put_on(self) -> None:
        assert DenseGP().configured(device=OTHER).provenance_config()["device"] == OTHER


class TestTheDeviceStringIsTorchs:
    """:func:`device_name` is ``str(torch.device(...))``, and nothing cleverer."""

    @pytest.mark.parametrize("given", ["cpu", "meta", "cuda", "cuda:1", torch.device("cpu")])
    def test_it_round_trips_through_torch(self, given: Any) -> None:
        assert device_name(given) == str(torch.device(given))

    def test_cuda_and_cuda_zero_stay_distinct(self) -> None:
        """Not a bug: torch's two spellings differ once a second GPU exists."""
        assert device_name("cuda") != device_name("cuda:0")

    def test_none_is_refused_rather_than_taken_as_choose_one(self) -> None:
        with pytest.raises(TypeError, match="never detects a device"):
            resolve_device(None)


# ---------------------------------------------------------------------------
# 2. Buffers move with .to(...)
# ---------------------------------------------------------------------------


class TestBuffersMove:
    """``architecture.md`` §5: "buffers move with parameters under ``.to(...)``"."""

    def test_a_model_moves_its_grid_and_re_declares_itself(self) -> None:
        model = a_model()
        assert model.tensors.get_buffer("wavelength").device.type == "cpu"
        returned = model.to(device=OTHER)
        assert returned is model, "a piece is held by identity; .to() must be in place"
        assert model.DEVICE == OTHER
        assert model.tensors.get_buffer("wavelength").device.type == OTHER

    def test_a_models_negotiated_grids_move_with_it(self) -> None:
        """The per-channel grids ``compile_for`` caches are buffers in all but name."""
        model = a_model()
        model.compile_for({"default": {}})
        model._grids["default"] = torch.zeros(3, dtype=torch.float64)
        model.to(device=OTHER)
        assert model._grids["default"].device.type == OTHER

    def test_a_step_moves_its_buffers(self) -> None:
        step = Resample([1.5, 3.0, 6.0])
        step.to(device=OTHER)
        assert step.DEVICE == OTHER
        assert all(tensor.device.type == OTHER for tensor in step.tensors.buffers())

    def test_a_dtype_move_leaves_the_device_alone(self) -> None:
        model = a_model()
        model.to(dtype=torch.float32)
        assert model.dtype == torch.float32
        assert model.DEVICE == "cpu"
        assert model.tensors.get_buffer("wavelength").dtype == torch.float32


# ---------------------------------------------------------------------------
# 3. The device rule composes a whole problem
# ---------------------------------------------------------------------------


def _problem(model_device: Any, step_device: Any, noise_device: Any) -> FittingProblem:
    return FittingProblem(
        a_model(device=model_device),
        [
            Dataset(
                OBSERVED,
                Instrument([CalibrationScale(st.lognorm(0.2), label="cal", device=step_device)]),
                Likelihood(GaussianFamily(), IndependentNoise(device=noise_device)),
            )
        ],
        seed=1,
    )


class TestTheDeviceRuleComposesAProblem:
    def test_one_device_everywhere_composes_and_is_reported(self) -> None:
        assert _problem("cpu", "cpu", "cpu").device == "cpu"

    @pytest.mark.parametrize(
        "devices",
        [(OTHER, "cpu", "cpu"), ("cpu", OTHER, "cpu"), ("cpu", "cpu", OTHER)],
        ids=["model", "step", "noise"],
    )
    def test_one_part_elsewhere_is_refused_at_composition(self, devices: Any) -> None:
        """The point of making the flag per-instance: the mistake is loud, and early.

        Every one of these three would otherwise surface as a torch device
        mismatch inside the first evaluation — or not at all, since torch
        broadcasts a CPU scalar against a foreign tensor without complaint.
        """
        with pytest.raises(DatasetError, match="different devices"):
            _problem(*devices)

    def test_a_gp_solver_elsewhere_is_refused_before_the_problem_is_built(self) -> None:
        """W2.13 put the solver among the capability parts; slice 3 made it move.

        The refusal arrives from the noise model rather than from
        ``declared_capabilities``, and earlier: a ``GaussianProcessNoise``
        checks its own kernel and solver at construction, which is before there
        is a ``Dataset`` to attach them to. ``declared_capabilities`` would
        catch the solver a step later and the *kernel* never, so the check is
        worth having at both altitudes.
        """
        with pytest.raises(LikelihoodError, match="DenseGP on 'meta'"):
            GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP().configured(device=OTHER))

    def test_a_kernel_elsewhere_is_refused_by_the_noise_model_itself(self) -> None:
        """A kernel is not a capability part, so the noise model checks it.

        ``Likelihood.capability_parts`` carries the noise model and the solver
        onto the problem, deliberately and by ruling — a kernel is a
        declaration. So a kernel left behind would reach torch rather than
        ``declared_capabilities``, and the message would be about strides.
        """
        with pytest.raises(LikelihoodError, match="Matern32 on 'meta'"):
            GaussianProcessNoise(Matern32(0.4, 2.0, device=OTHER), DenseGP(), device="cpu")

    def test_the_realisation_lives_where_the_problem_says_it_does(self) -> None:
        lowered = lower_problem(_problem("cpu", "cpu", "cpu"))
        assert str(lowered.device) == "cpu"
        assert repr(lowered).endswith("device='cpu'>")
        assert lowered._datasets[0].observed_values.device.type == "cpu"


class TestThePriorDrawIsDeviceAware:
    """``lowering.md`` §9's two routes, and the device each of them draws on.

    Both were CPU-only by accident: route (2) asked ``torch.rand`` for a draw
    on the space's device with a **CPU** generator (torch refuses that
    outright), and route (1) forked ``devices=[]``, which saves and restores
    the CPU generator and nothing else. Nothing in the gate could see either,
    because ``FittingProblem.sample_prior`` goes through ``ampere.core``'s numpy
    path and the only caller of this surface is the conformance parameter
    battery, on the CPU.

    Only the *mechanism* can be pinned here. The CPU-generator choice is
    invisible on a CPU space by construction — the two candidate
    implementations agree there — and the meta stand-in cannot help either,
    because building a ``TorchParameterSpace`` on it fails earlier
    (``lower_prior`` takes ``float()`` of a support bound, and a meta tensor
    has no data). ``tests/gpu`` holds the behaviour, on a device that has a
    generator of its own.
    """

    def test_a_cpu_space_forks_no_extra_device(self) -> None:
        from ampere.backends.torch.parameters import _fork_devices

        assert _fork_devices(torch.device("cpu")) == []
        assert _fork_devices(torch.device("meta")) == []

    def test_a_cuda_space_forks_its_own_generator(self) -> None:
        from ampere.backends.torch.parameters import _fork_devices

        assert _fork_devices(torch.device("cuda")) == [torch.device("cuda")]
        assert _fork_devices(torch.device("cuda:1")) == [torch.device("cuda:1")]

    def test_the_draw_is_still_reproducible_on_the_cpu(self) -> None:
        """The property the change must not have moved."""
        from ampere.backends.torch import TorchParameterSpace

        problem = _problem("cpu", "cpu", "cpu")
        space = TorchParameterSpace(problem.parameters)
        first = space.sample(4242)
        second = space.sample(4242)
        assert set(first) == set(problem.parameters.names)
        for name, value in first.items():
            assert np.array_equal(value, second[name])


# ---------------------------------------------------------------------------
# 4. The sigma_tensor hook has no silent consumers left
# ---------------------------------------------------------------------------


class _SilentlyPredictionAware(IndependentNoise):
    """A noise model that overrides ``sigma`` and supplies no ``sigma_tensor``.

    Exactly the shape W2.4 slice 3 item 4 is about: on the contract path it
    doubles the uncertainties, and before slice 3 the realised path would have
    used the *base* quadrature instead and scored a different likelihood
    without saying anything.
    """

    def sigma(self, observed: Any, retain: Any, values: Any, *, predicted: Any = None) -> Any:
        base = super().sigma(observed, retain, values, predicted=predicted)
        return None if base is None else 2.0 * base


class TestThePredictionAwareHook:
    def test_a_noise_model_that_overrides_sigma_without_the_hook_is_refused(self) -> None:
        problem = FittingProblem(
            a_model(),
            [
                Dataset(
                    OBSERVED,
                    likelihood=Likelihood(GaussianFamily(), _SilentlyPredictionAware()),
                )
            ],
            seed=1,
        )
        with pytest.raises(LoweringError, match="sigma_tensor"):
            lower_problem(problem)

    @pytest.mark.parametrize(
        "noise",
        [IndependentNoise(), GaussianProcessNoise(Matern32(0.4, 2.0))],
        ids=["independent", "gp"],
    )
    def test_the_two_transcribed_noise_models_are_not_caught_by_it(self, noise: Any) -> None:
        problem = FittingProblem(
            a_model(),
            [Dataset(OBSERVED, likelihood=Likelihood(GaussianFamily(), noise))],
            seed=1,
        )
        assert lower_problem(problem) is not None

    def test_the_shipped_prediction_aware_pair_supplies_the_hook(self) -> None:
        pair = (
            FractionalModelNoise(0.05),
            FractionalModelGPNoise(Matern32(0.4, 2.0), f=0.05),
        )
        for noise in pair:
            assert hasattr(noise, "sigma_tensor")

    def test_the_hook_takes_the_modulus_of_a_complex_prediction(self) -> None:
        """The other half of item 4: ``|predicted|``, not ``Re(predicted)``.

        Coercing a complex tensor to float64 first drops the imaginary part;
        the amplitude a prediction-aware model inflates by must be the modulus.
        """
        noise = FractionalModelNoise(0.5)
        predicted = torch.tensor([3.0 + 4.0j], dtype=torch.complex128)
        got = noise.sigma_tensor(None, {"f": 0.5}, predicted=predicted)
        assert float(got[0]) == pytest.approx(2.5)  # 0.5 * |3+4i| = 2.5


# ---------------------------------------------------------------------------
# 5. complex_gaussian
# ---------------------------------------------------------------------------


class TestTheComplexGaussianBody:
    """The transcription against ``ampere.core``'s own closed form."""

    def test_it_agrees_with_the_core_family(self) -> None:
        from ampere.backends.torch._families import FamilyInputs, native_log_prob

        observed = np.array([1.4 + 0.3j, -0.8 + 1.1j, 0.2 - 0.9j, 1.0 + 0.0j])
        predicted = np.array([1.5 + 0.2j, -0.7 + 1.0j, 0.1 - 1.0j, 0.9 + 0.1j])
        sigma = np.array([0.1, 0.2, 0.15, 0.3])
        family = ComplexGaussianFamily()
        expected = family.log_prob(predicted, observed, NoiseParams(sigma=sigma, values={}))
        got = native_log_prob(
            FamilyInputs(
                predicted=torch.as_tensor(predicted, dtype=torch.complex128),
                observed=torch.as_tensor(observed, dtype=torch.complex128),
                sigma=torch.as_tensor(sigma, dtype=torch.float64),
                values={},
                family=family,
            )
        )
        assert float(got) == pytest.approx(expected, abs=1e-12)

    def test_the_observed_container_stays_complex_in_the_lowering(self) -> None:
        """The one line slice 2's ``dtype=float`` made impossible."""
        from tests.conformance.backends.torch_backend import TorchBackend
        from tests.conformance.composition import build_problem
        from tests.conformance.test_inference import COMPLEX

        lowered = lower_problem(build_problem(TorchBackend(), COMPLEX))
        dataset = lowered._datasets[0]
        assert dataset.complex_valued is True
        assert dataset.observed_values.dtype == complex_dtype(torch.float64)

    def test_a_complex_gp_is_no_longer_refused(self) -> None:
        """**W4.2**: the closed form exists, so the staging refusal does not fire.

        Until Phase 4 this row asserted the opposite, and the reason was a rule
        rather than a shortage of time: ``ampere.core`` declared
        ``GP_ANALYTIC_IMPLEMENTED = False``, so there was no numpy oracle and a
        native path would have been this backend inventing a likelihood. The
        closed form is written now, and it lives on the **solver** — the
        ``gp_marginal`` branch of :mod:`ampere.backends.torch.problem`, stacking
        the complex residual into two real columns — not in this module, which is
        why the family body below is still the uncorrelated one.
        """
        from ampere.backends.torch._families import refuse_family

        assert (
            refuse_family(ComplexGaussianFamily(), censored=False, latent=False, correlated=True)
            is None
        )

    def test_the_staging_guard_still_refuses_a_family_that_has_not_written_one(self) -> None:
        """The guard keeps its rule and loses its instance (*W4.2*).

        ``refuse_family``'s ``correlated`` argument exists for a family that
        declares ``ANALYTIC_WITH_GP`` without ``GP_ANALYTIC_IMPLEMENTED``. No
        shipped family is in that position any more, so the rule is held against
        one declared here — otherwise lifting the complex refusal would have
        quietly removed the check along with its subject.
        """
        from ampere.backends.torch._families import refuse_family

        class StagedFamily(ComplexGaussianFamily):
            NAME = "staged_complex"
            GP_ANALYTIC_IMPLEMENTED = False

        refusal = refuse_family(StagedFamily(), censored=False, latent=False, correlated=True)
        assert refusal is not None
        assert "staged_complex" in str(refusal)
        assert "GP_ANALYTIC_IMPLEMENTED is False" in str(refusal)

    def test_the_circular_complex_gp_realisation_agrees_with_the_contract_path(self) -> None:
        """The native closed form against its numpy oracle (*W4.2*).

        ``inference.md`` §10a makes the numpy path the thing a realisation is
        checked against, and this is that check for the circular complex GP. Two
        things only this shape exercises on the realised path: the kernel's
        ``axes=("u", "v")`` selector resolved against the lowering's own
        coordinate block (a ``VisibilitySet`` has three axes, and the lowering
        read only the first until W4.2), and the complex residual stacked into
        the solver's two-column right-hand side.
        """
        from tests.conformance.backends.torch_backend import TorchBackend
        from tests.conformance.composition import build_problem
        from tests.conformance.test_inference import COMPLEX_GP

        problem = build_problem(TorchBackend(), COMPLEX_GP)
        lowered = lower_problem(problem)
        theta = problem.parameters.pack(problem.reference_values)
        expected = problem.evaluate(theta).log_likelihood
        terms = lowered.log_likelihood_terms(problem.unconstrain(theta))
        assert float(terms["sed"].detach()) == pytest.approx(expected, abs=1e-9)

    def test_the_circular_complex_gp_gives_the_hyperparameters_a_gradient(self) -> None:
        """A native solver exists so the kernel's knobs can be fitted.

        W2.4 slice 1's finding, asked of the complex path: a lowering that
        coerced anything to numpy between the kernel and the density would leave
        the amplitude and the length scale with a derivative of exactly zero, and
        the fit would sample them against the prior alone.
        """
        from tests.conformance.backends.torch_backend import TorchBackend
        from tests.conformance.composition import build_problem
        from tests.conformance.protocol import CovarianceSpec, KernelFamily
        from tests.conformance.test_inference import COMPLEX_GP
        import dataclasses

        fitted = dataclasses.replace(
            COMPLEX_GP,
            datasets=(
                dataclasses.replace(
                    COMPLEX_GP.datasets[0],
                    covariance=CovarianceSpec(
                        KernelFamily.MATERN32,
                        st.loguniform(1e-2, 1.0),
                        st.loguniform(0.5, 10.0),
                        axes=("u", "v"),
                    ),
                ),
            ),
        )
        problem = build_problem(TorchBackend(), fitted)
        lowered = lower_problem(problem)
        point = torch.tensor(problem.unconstrain(problem.reference_values), requires_grad=True)
        lowered.log_prob_unconstrained(point).backward()
        gradient = point.grad
        assert gradient is not None
        assert torch.all(torch.isfinite(gradient))
        # The kernel's two hyperparameters are the last columns, in declaration
        # order: the noise model registers them before the family's (none here).
        assert bool(torch.all(gradient[-2:] != 0.0))

    def test_the_gradient_reaches_the_phase(self) -> None:
        """A complex residual whose gradient stops would be silent otherwise.

        ``ModelKind.COMPLEX``'s modulus does not depend on ``index``, so a
        derivative of zero here means the imaginary part was dropped rather
        than that the point was stationary.
        """
        from tests.conformance.backends.torch_backend import TorchBackend
        from tests.conformance.composition import build_problem
        from tests.conformance.test_inference import COMPLEX

        problem = build_problem(TorchBackend(), COMPLEX)
        lowered = realise(problem)
        y = torch.tensor(
            problem.unconstrain(problem.reference_values),
            dtype=torch.float64,
            requires_grad=True,
        )
        lowered.log_prob_unconstrained(y).backward()
        assert y.grad is not None
        index = problem.parameters.free_names.index("model.index")
        assert abs(float(y.grad[index])) > 1e-6


# ---------------------------------------------------------------------------
# 6. tests/gpu is real code, exercised on the CPU
# ---------------------------------------------------------------------------


class TestTheGpuModuleIsExercised:
    """``tests/gpu`` skips everywhere the gate runs, so nothing else would run it.

    A suite that only executes on hardware nobody has in CI is a suite that
    rots. Its problem builders take a ``device=`` argument precisely so that
    this row can compose them on the CPU: everything but the device string is
    then covered by the gate, and what CUDA adds is the part the ruling says
    ampere does not re-measure anyway.
    """

    def test_the_toy_joint_problem_builds_and_realises_on_the_cpu(self) -> None:
        module = pytest.importorskip("tests.gpu.test_torch_gpu")
        problem = module.joint_problem(device="cpu")
        assert problem.device == "cpu"
        assert problem.free_size == 3
        lowered = realise(problem)
        y = problem.unconstrain(problem.reference_values)
        assert np.isfinite(float(lowered.log_prob_unconstrained(y)))

    def test_the_gp_problem_builds_and_realises_on_the_cpu(self) -> None:
        module = pytest.importorskip("tests.gpu.test_torch_gpu")
        problem = module.gp_problem(device="cpu")
        assert problem.device == "cpu"
        noise = problem.datasets["default"].likelihood.noise
        assert noise.kernel.DEVICE == noise.solver.DEVICE == noise.DEVICE == "cpu"
        lowered = realise(problem)
        y = problem.unconstrain(problem.reference_values)
        assert np.isfinite(float(lowered.log_prob_unconstrained(y)))
