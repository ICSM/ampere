"""Unit tests for ``ampere.backends.torch``'s models, steps and solver.

The conformance battery already holds the shipped ``PowerLaw``,
``CalibrationScale``, ``Resample``, ``TorchParameterSpace`` and ``DenseGP`` to
the §4 contracts, once per registered backend, and holds them against the
reference backend directly in ``test_cross_backend.py``. Those rows are not
repeated here. What this file covers is what the battery cannot reach:

* the **physics**, against an independent oracle — the Planck function from
  ``astropy.constants``, and the two limits where the arithmetic is delicate;
* the two steps with no conformance kind, :class:`LSFConvolution` and the real
  :class:`SyntheticPhotometry`;
* the **gradients**, which are the whole reason this backend exists and which
  the neutral battery has no vocabulary for: every tensor-valued entry point is
  differentiated here and checked against a closed form or a finite difference;
* the **agreement with the reference backend piece by piece**, at a tolerance
  tight enough to catch a translation error rather than a rounding one.

``ampere.backends.reference`` is imported as an oracle in a handful of rows.
That is deliberate and is not "comparing ampere against ampere" in the sense
``tests/conformance/README.md`` forbids: the *claim* under test is precisely
that two backends implementing one declaration agree, and the reference path is
what ``architecture.md`` §2 appoints to say what the answer is. Where a closed
form exists, it is used instead.
"""

from __future__ import annotations

import math

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

pytest.importorskip("torch", reason="the torch backend needs ampere[torch]")

import torch

from ampere.backends import reference as ref
from ampere.backends.torch import (
    BACKEND,
    BlackBody,
    CalibrationScale,
    DenseGP,
    LSFConvolution,
    ModifiedBlackBody,
    PowerLaw,
    QuasisepGP,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    planck_jy,
)
from ampere.backends.torch import (
    Matern32 as TorchMatern32,
)
from ampere.backends.torch import (
    SquaredExponential as TorchSquaredExponential,
)
from ampere.core import (
    DenseGP as CoreDenseGP,
)
from ampere.core import (
    QuasisepGP as CoreQuasisepGP,
)
from ampere.core import (
    Matern32,
    Model,
    ParameterError,
    Spectrum,
    Transformation,
    TransformationError,
    negotiate,
)

GRID = np.geomspace(1.0, 30.0, 200)


def tensor(values: object) -> torch.Tensor:
    return torch.as_tensor(np.asarray(values, dtype=float), dtype=torch.float64)


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------


class TestPlanckFunction:
    """Against ``astropy.constants``, assembled independently of the module."""

    def test_it_matches_the_closed_form(self) -> None:
        temperature = 1500.0
        wavelength = np.array([1.0, 5.0, 20.0, 100.0])
        frequency = (const.c / (wavelength * u.micron)).to(u.Hz)
        expected = (
            2.0
            * const.h
            * frequency**3
            / const.c**2
            / np.expm1((const.h * frequency / (const.k_B * temperature * u.K)).to_value(""))
        ).to_value(u.Jy)
        got = planck_jy(tensor(wavelength), tensor(temperature)).numpy()
        assert got == pytest.approx(expected, rel=1e-12)

    def test_the_rayleigh_jeans_tail_survives_the_cancellation(self) -> None:
        """``expm1`` rather than ``exp(x) - 1``: at long wavelengths ``x`` is tiny.

        The Rayleigh-Jeans form is an *approximation* — ``B_nu`` approaches it
        as ``x = h nu / kT -> 0``, with a leading relative error of ``x/2`` —
        so the evidence that no cancellation is happening is not that the two
        agree to machine precision (they must not) but that the deviation
        *shrinks like ``x``* as the wavelength grows, decade for decade. A
        naive ``exp(x) - 1`` would lose its significant digits to cancellation
        long before the third decade, and the deviation would stop tracking
        ``x/2`` and start tracking the rounding.
        """
        wavelength = np.array([1.0e5, 1.0e6, 1.0e7])
        temperature = 30.0
        frequency = (const.c / (wavelength * u.micron)).to(u.Hz)
        expected = (2.0 * frequency**2 * const.k_B * temperature * u.K / const.c**2).to_value(u.Jy)
        got = planck_jy(tensor(wavelength), tensor(temperature)).numpy()
        deviation = np.abs(got / expected - 1.0)
        x = (const.h * frequency / (const.k_B * temperature * u.K)).to_value("")
        assert deviation == pytest.approx(x / 2.0, rel=1e-3)

    def test_the_wien_side_underflows_to_zero_rather_than_to_nan(self) -> None:
        """The clamp, and why it differs from the reference backend's guard.

        numpy can form an ``inf`` and replace it afterwards; ``torch.where``
        evaluates both branches and would carry a ``nan`` gradient out of the
        discarded one, so the exponent is clamped before ``expm1`` instead. The
        observable behaviour must be the same: a vanishing radiance, never a
        ``nan``.
        """
        got = planck_jy(tensor([0.001, 0.01]), tensor(3.0)).numpy()
        assert np.all(np.isfinite(got))
        assert np.all(got >= 0.0)
        assert np.all(got < 1e-100)

    def test_it_agrees_with_the_reference_backend(self) -> None:
        wavelength = np.geomspace(0.5, 500.0, 97)
        for temperature in (30.0, 300.0, 5800.0):
            got = planck_jy(tensor(wavelength), tensor(temperature)).numpy()
            expected = ref.planck_jy(wavelength, temperature)
            assert got == pytest.approx(expected, rel=1e-13)

    def test_it_refuses_an_impossible_temperature(self) -> None:
        with pytest.raises(ValueError, match="finite and positive"):
            planck_jy(tensor([1.0]), tensor(-5.0))

    def test_it_refuses_a_non_positive_wavelength(self) -> None:
        with pytest.raises(ValueError, match="strictly positive"):
            planck_jy(tensor([0.0, 1.0]), tensor(100.0))


class TestNativeModels:
    """The trio, against the reference backend and against their own gradients."""

    @pytest.mark.parametrize(
        ("build_torch", "build_reference", "values"),
        [
            (BlackBody, ref.BlackBody, {"temperature": 1200.0, "scale": 3.0}),
            (
                ModifiedBlackBody,
                ref.ModifiedBlackBody,
                {"temperature": 40.0, "beta": 1.8, "scale": 0.5},
            ),
            (PowerLaw, ref.PowerLaw, {"norm": 2.5, "index": -1.3}),
        ],
    )
    def test_each_model_agrees_with_the_reference_backend(
        self, build_torch: type, build_reference: type, values: dict[str, float]
    ) -> None:
        got = build_torch(GRID)(**values)["default"].values
        expected = build_reference(GRID)(**values)["default"].values
        assert got == pytest.approx(expected, rel=1e-12)

    def test_every_model_declares_the_four_capability_flags(self) -> None:
        for build in (BlackBody, ModifiedBlackBody, PowerLaw):
            model = build(GRID)
            assert model.BACKEND == BACKEND
            assert model.DIFFERENTIABLE is True
            # True since W2.4 slice 2: the realised density is vmap-able.
            assert model.BATCHABLE is True
            assert model.DEVICE == "cpu"

    def test_the_grid_is_a_torch_buffer_and_a_declared_one(self) -> None:
        """``lowering.md`` §7: a constant array is a buffer, never a parameter.

        And it is declared *twice*, deliberately: once in ``ampere.core``'s
        ``BufferSet``, which is what provenance and the neutral model identity
        read, and once as a torch buffer, which is what a device move and a
        ``state_dict`` see.
        """
        model = PowerLaw(GRID)
        assert "wavelength" in model.buffers
        assert "wavelength" in dict(model.tensors.named_buffers())
        assert not model.tensors.get_buffer("wavelength").requires_grad
        assert model.tensors.get_buffer("wavelength").dtype is torch.float64
        assert list(model.tensors.parameters()) == []

    def test_the_tensor_entry_point_is_differentiable(self) -> None:
        """The claim ``DIFFERENTIABLE = True`` makes, checked against a closed form.

        ``d/d index [norm * (x / x_ref)**index] = flux * log(x / x_ref)``, so
        the gradient of the summed flux has an exact value rather than a
        finite-difference one.
        """
        reference_wavelength = 2.0
        model = PowerLaw(GRID, reference_wavelength=reference_wavelength)
        norm = torch.tensor(2.5, dtype=torch.float64, requires_grad=True)
        index = torch.tensor(-1.3, dtype=torch.float64, requires_grad=True)
        flux = model.evaluate_tensor(norm=norm, index=index)["default"]
        flux.sum().backward()
        ratio = GRID / reference_wavelength
        assert float(norm.grad) == pytest.approx(float(np.sum(ratio**-1.3)), rel=1e-12)
        assert float(index.grad) == pytest.approx(
            float(np.sum(2.5 * ratio**-1.3 * np.log(ratio))), rel=1e-12
        )

    def test_the_blackbody_temperature_takes_a_gradient(self) -> None:
        """A finite-difference check, because ``dB/dT`` has no short closed form."""
        model = BlackBody(GRID)
        temperature = torch.tensor(1200.0, dtype=torch.float64, requires_grad=True)
        model.evaluate_tensor(temperature=temperature, scale=tensor(1.0))[
            "default"
        ].sum().backward()
        step = 1e-4
        plus = float(np.sum(ref.planck_jy(GRID, 1200.0 + step)))
        minus = float(np.sum(ref.planck_jy(GRID, 1200.0 - step)))
        assert float(temperature.grad) == pytest.approx((plus - minus) / (2.0 * step), rel=1e-6)

    def test_the_container_boundary_is_where_the_graph_stops(self) -> None:
        """The package's principal caveat, stated as a test rather than prose.

        ``evaluate`` must succeed with tensor inputs that carry a graph — a
        model that raised there would be unusable — and the ``Spectrum`` it
        returns necessarily holds detached numpy. The difference between the
        two entry points is the whole of what W2.4's report escalates.
        """
        model = PowerLaw(GRID)
        norm = torch.tensor(2.5, dtype=torch.float64, requires_grad=True)
        emitted = model(norm=norm, index=tensor(-1.3))["default"]
        assert isinstance(emitted.values, np.ndarray)
        assert model.evaluate_tensor(norm=norm, index=tensor(-1.3))["default"].requires_grad

    def test_compile_for_adopts_the_negotiated_grid(self) -> None:
        """``transformations.md`` §14: one container per channel, refilled."""
        target = np.linspace(2.0, 20.0, 31)
        model = PowerLaw(GRID)
        step = Resample(target)
        requirements = negotiate([_instrument(step)])
        compiled = model.compile_for(requirements)
        first = compiled(norm=1.0, index=-1.0)["default"]
        second = compiled(norm=2.0, index=-1.0)["default"]
        assert first.spectral_axis is second.spectral_axis
        assert first.values.size == requirements["default"]["spectral_axis"].coordinates().size

    def test_a_parameter_may_be_fixed_with_a_tensor(self) -> None:
        model = PowerLaw(GRID, norm=tensor(3.0))
        assert model.parameters["norm"].is_fixed
        assert model.parameters["norm"].value == pytest.approx(3.0)

    def test_a_multi_element_tensor_is_refused_as_a_fixed_value(self) -> None:
        with pytest.raises(ParameterError, match="scalar"):
            PowerLaw(GRID, norm=tensor([1.0, 2.0]))


def _instrument(*steps: Transformation):
    from ampere.core import Instrument

    return Instrument(steps, channel="default", input_kind=Spectrum, label="sed")


# ---------------------------------------------------------------------------
# Instrument steps
# ---------------------------------------------------------------------------


def flat(grid: np.ndarray, value: float = 1.0, mask: np.ndarray | None = None) -> Spectrum:
    return Spectrum(grid * u.micron, np.full(grid.size, value), unit=u.Jy, mask=mask)


class TestBinEdges:
    def test_the_midpoints_bracket_the_centres(self) -> None:
        edges = bin_edges(tensor([1.0, 2.0, 4.0])).numpy()
        assert edges == pytest.approx([0.5, 1.5, 3.0, 5.0])

    def test_a_single_sample_gets_a_unit_bin(self) -> None:
        assert bin_edges(tensor([3.0])).numpy() == pytest.approx([2.5, 3.5])

    def test_it_agrees_with_the_reference_backend(self) -> None:
        grid = np.geomspace(1.0, 50.0, 41)
        assert bin_edges(tensor(grid)).numpy() == pytest.approx(ref.bin_edges(grid), rel=1e-15)


class TestCalibrationScale:
    def test_it_multiplies_and_keeps_the_axes(self) -> None:
        step = CalibrationScale(st.lognorm(0.2))
        samples = flat(np.array([1.0, 2.0, 3.0]), 2.0)
        scaled = step(samples, {"scale": 1.5})
        assert scaled.values == pytest.approx(3.0)
        assert scaled.spectral_axis is samples.spectral_axis

    def test_it_publishes_no_requirements(self) -> None:
        assert CalibrationScale(1.0).requirements() == ()

    def test_the_factor_takes_a_gradient(self) -> None:
        step = CalibrationScale(st.lognorm(0.2))
        factor = torch.tensor(1.5, dtype=torch.float64, requires_grad=True)
        step.apply_tensor(tensor([2.0, 4.0]), factor).sum().backward()
        assert float(factor.grad) == pytest.approx(6.0)


class TestResample:
    def test_the_influence_matrix_matches_the_reference_backend(self) -> None:
        source = np.linspace(1.0, 10.0, 60)
        target = np.linspace(2.0, 8.0, 11)
        got = Resample(target).influence(source)
        assert got == pytest.approx(ref.Resample(target).influence(source), abs=1e-15)

    def test_each_row_sums_to_one(self) -> None:
        """A mean, not a sum: the containers carry a flux density."""
        source = np.linspace(1.0, 10.0, 60)
        weights = Resample(np.linspace(2.0, 8.0, 11)).influence(source)
        assert weights.sum(axis=1) == pytest.approx(1.0)

    def test_a_flat_spectrum_resamples_to_itself(self) -> None:
        target = np.linspace(2.0, 8.0, 11)
        out = Resample(target)(flat(np.linspace(1.0, 10.0, 60), 3.0), None)
        assert out.values == pytest.approx(3.0)

    def test_it_propagates_the_mask_by_the_any_rule(self) -> None:
        source = np.linspace(1.0, 10.0, 60)
        mask = np.zeros(60, dtype=bool)
        mask[30] = True
        out = Resample(np.linspace(2.0, 8.0, 11))(flat(source, 1.0, mask), None)
        assert out.mask is not None
        assert bool(np.any(out.mask))

    def test_it_publishes_a_density_requirement_over_the_target_range(self) -> None:
        target = np.linspace(2.0, 8.0, 11)
        (requirement,) = Resample(target).requirements()
        assert requirement.axis == "spectral_axis"
        assert requirement.intervals[0][0] == pytest.approx(2.0)

    def test_it_refuses_an_unsorted_target(self) -> None:
        with pytest.raises(TransformationError, match="strictly increasing"):
            Resample(np.array([3.0, 1.0, 2.0]))

    def test_the_binning_integral_is_differentiable_in_the_flux(self) -> None:
        source = np.linspace(1.0, 10.0, 30)
        step = Resample(np.linspace(2.0, 8.0, 7))
        values = torch.ones(30, dtype=torch.float64, requires_grad=True)
        weights = step.influence_tensor(tensor(source))
        step.apply_tensor(weights, values).sum().backward()
        assert values.grad is not None
        assert float(values.grad.sum()) == pytest.approx(float(weights.sum()))


class TestLSFConvolution:
    def test_the_kernel_matches_the_reference_backend(self) -> None:
        grid = np.geomspace(2.0, 20.0, 120)
        got = LSFConvolution(resolving_power=200.0).influence(grid)
        assert got == pytest.approx(ref.LSFConvolution(resolving_power=200.0).influence(grid))

    def test_a_constant_width_kernel_matches_too(self) -> None:
        grid = np.linspace(2.0, 20.0, 120)
        got = LSFConvolution(fwhm=0.5).influence(grid)
        assert got == pytest.approx(ref.LSFConvolution(fwhm=0.5).influence(grid))

    def test_a_flat_spectrum_is_unchanged_by_convolution(self) -> None:
        grid = np.linspace(2.0, 20.0, 120)
        out = LSFConvolution(fwhm=0.5)(flat(grid, 4.0), None)
        assert out.values == pytest.approx(4.0)

    def test_it_takes_exactly_one_width_description(self) -> None:
        with pytest.raises(TransformationError, match="both"):
            LSFConvolution(resolving_power=100.0, fwhm=0.1)
        with pytest.raises(TransformationError, match="neither"):
            LSFConvolution()

    def test_it_stays_silent_until_configure_from_finds_a_range(self) -> None:
        step = LSFConvolution(resolving_power=100.0)
        assert step.requirements() == ()

    def test_configure_from_pads_the_downstream_range(self) -> None:
        """gap I-3's mechanism: the LSF learns its range from its successors."""
        target = np.linspace(5.0, 10.0, 21)
        step = LSFConvolution(fwhm=0.5)
        step.configure_from([Resample(target)])
        (requirement,) = step.requirements()
        low, high = requirement.intervals[0]
        assert low < 5.0
        assert high > 10.0


class TestSyntheticPhotometry:
    """The real step: response curves, the two detector conventions, ``points=``."""

    @staticmethod
    def top_hat(grid: np.ndarray, low: float, high: float) -> np.ndarray:
        return ((grid >= low) & (grid <= high)).astype(float)

    def build(self, detector: str) -> SyntheticPhotometry:
        grid = np.linspace(1.0, 3.0, 201)
        response = np.vstack([self.top_hat(grid, 1.2, 1.8), self.top_hat(grid, 2.2, 2.8)])
        return SyntheticPhotometry(("blue", "red"), grid, response, detector=detector)

    @pytest.mark.parametrize("detector", ["photon", "energy"])
    def test_it_agrees_with_the_reference_backend(self, detector: str) -> None:
        grid = np.linspace(1.0, 3.0, 201)
        response = np.vstack([self.top_hat(grid, 1.2, 1.8), self.top_hat(grid, 2.2, 2.8)])
        spectrum = Spectrum(grid * u.micron, 2.0 + np.sin(grid), unit=u.Jy)
        got = self.build(detector)(spectrum, None)
        expected = ref.SyntheticPhotometry(("blue", "red"), grid, response, detector=detector)(
            spectrum, None
        )
        assert got.values == pytest.approx(expected.values, rel=1e-13)
        assert got.spectral_axis.values == pytest.approx(expected.spectral_axis.values, rel=1e-13)

    def test_a_flat_spectrum_gives_its_own_level(self) -> None:
        """The weights are normalised, so a constant flux passes through."""
        grid = np.linspace(1.0, 3.0, 201)
        out = self.build("photon")(flat(grid, 7.0), None)
        assert out.values == pytest.approx(7.0)

    def test_the_two_conventions_differ(self) -> None:
        """Different numbers, not different spellings — which is why there is no default."""
        grid = np.linspace(1.0, 3.0, 201)
        spectrum = Spectrum(grid * u.micron, 2.0 + np.sin(5.0 * grid), unit=u.Jy)
        photon = self.build("photon")(spectrum, None).values
        energy = self.build("energy")(spectrum, None).values
        assert not np.allclose(photon, energy)

    def test_it_publishes_points_at_its_own_tabulation(self) -> None:
        """``spectrum_photometry.md`` Gap 1: coverage is not enough."""
        step = self.build("photon")
        (requirement,) = step.requirements()
        assert requirement.points is not None
        assert np.asarray(requirement.points, dtype=float) == pytest.approx(step.tabulation())

    def test_it_reads_a_wider_negotiated_grid_by_lookup(self) -> None:
        """The failure ``Axis.locate`` exists to prevent, exercised.

        A second instrument widens the negotiated grid; reading the response
        matrix positionally against it would multiply each response column by
        the wrong flux. The output must be unchanged.
        """
        step = self.build("photon")
        narrow = flat(step.tabulation(), 7.0)
        wide_grid = np.union1d(step.tabulation(), np.linspace(0.5, 4.0, 57))
        wide = flat(wide_grid, 7.0)
        assert step(narrow, None).values == pytest.approx(step(wide, None).values)

    def test_the_detector_convention_is_recorded_in_describe(self) -> None:
        assert self.build("energy").describe() == {"detector": ["energy", "energy"]}

    def test_it_refuses_an_unknown_convention(self) -> None:
        with pytest.raises(TransformationError, match="unknown detector"):
            self.build("bolometric")

    def test_it_refuses_repeated_filter_names(self) -> None:
        grid = np.linspace(1.0, 3.0, 51)
        response = np.vstack([self.top_hat(grid, 1.2, 1.8)] * 2)
        with pytest.raises(TransformationError, match="repeated filter names"):
            SyntheticPhotometry(("blue", "blue"), grid, response, detector="photon")

    def test_it_refuses_a_response_that_integrates_to_zero(self) -> None:
        grid = np.linspace(1.0, 3.0, 51)
        response = np.zeros((1, grid.size))
        with pytest.raises(TransformationError, match="integrates to zero"):
            SyntheticPhotometry(("dead",), grid, response, detector="photon")

    def test_every_step_declares_the_four_capability_flags(self) -> None:
        for step in (
            CalibrationScale(1.0),
            Resample(np.linspace(1.0, 2.0, 5)),
            LSFConvolution(fwhm=0.1),
            self.build("photon"),
        ):
            assert step.BACKEND == BACKEND
            assert step.DIFFERENTIABLE is True
            # True since W2.4 slice 2, as for the models.
            assert step.BATCHABLE is True
            assert step.DEVICE == "cpu"


# ---------------------------------------------------------------------------
# The GP solver
# ---------------------------------------------------------------------------


def gp_case(size: int = 30) -> tuple[Matern32, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(20260907)
    coordinates = np.sort(rng.uniform(0.0, 12.0, size))
    residual = rng.normal(0.0, 0.3, size)
    variance = np.full(size, 0.04)
    return Matern32(0.4, 2.0), coordinates, residual, variance


class TestDenseGP:
    """Torch's Cholesky against scipy's, and against the multivariate normal."""

    def test_the_marginal_likelihood_is_the_multivariate_normal(self) -> None:
        """A closed form, not the reference solver: a genuinely external oracle."""
        kernel, coordinates, residual, variance = gp_case()
        covariance = kernel.matrix(coordinates, coordinates, {}) + np.diag(variance)
        expected = st.multivariate_normal(cov=covariance).logpdf(residual)
        got = DenseGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {})
        assert got == pytest.approx(float(expected), abs=1e-8)

    def test_it_agrees_with_the_reference_solver(self) -> None:
        kernel, coordinates, residual, variance = gp_case()
        assert DenseGP().log_marginal_likelihood(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreDenseGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {}),
            abs=1e-9,
        )

    def test_the_leave_one_out_terms_agree_with_the_reference_solver(self) -> None:
        kernel, coordinates, residual, variance = gp_case()
        assert DenseGP().conditional_loo(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreDenseGP().conditional_loo(kernel, coordinates, residual, variance, {}), abs=1e-9
        )

    def test_the_conditional_agrees_with_the_reference_solver(self) -> None:
        kernel, coordinates, residual, variance = gp_case()
        at = np.linspace(-1.0, 13.0, 41)
        got = DenseGP().condition(kernel, coordinates, residual, variance, {}, at=at)
        expected = CoreDenseGP().condition(kernel, coordinates, residual, variance, {}, at=at)
        assert got.mean == pytest.approx(expected.mean, abs=1e-9)
        assert got.variance == pytest.approx(expected.variance, abs=1e-9)

    def test_the_whitening_transform_agrees_with_the_reference_solver(self) -> None:
        kernel, coordinates, _, _ = gp_case()
        whitened = np.random.default_rng(7).normal(size=coordinates.size)
        assert DenseGP().latent_transform(kernel, coordinates, whitened, {}) == pytest.approx(
            CoreDenseGP().latent_transform(kernel, coordinates, whitened, {}), abs=1e-9
        )

    def test_the_marginal_likelihood_is_differentiable_in_the_residual(self) -> None:
        """The reason a torch dense solver exists at all.

        ``d/dr log N(r; 0, C) = -C**-1 r``, so the gradient has a closed form
        rather than a finite-difference one.
        """
        kernel, coordinates, residual, variance = gp_case(20)
        residuals = torch.tensor(residual, dtype=torch.float64, requires_grad=True)
        DenseGP().log_marginal_likelihood_tensor(
            kernel, coordinates, residuals, variance, {}
        ).backward()
        covariance = kernel.matrix(coordinates, coordinates, {}) + np.diag(variance)
        expected = -np.linalg.solve(covariance, residual)
        assert residuals.grad.numpy() == pytest.approx(expected, abs=1e-8)

    def test_a_singular_covariance_raises_rather_than_returning_a_nan(self) -> None:
        """W2.3's finding, applied to this solver: never a silent NaN."""
        from ampere.core import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        coordinates = np.array([1.0, 1.0, 2.0])
        with pytest.raises(LikelihoodError, match="positive definite"):
            DenseGP().log_marginal_likelihood(kernel, coordinates, np.zeros(3), np.zeros(3), {})

    def test_a_non_finite_covariance_raises_by_name(self) -> None:
        from ampere.core import LikelihoodError

        kernel, coordinates, residual, _ = gp_case(5)
        with pytest.raises(LikelihoodError, match="non-finite"):
            DenseGP().log_marginal_likelihood(
                kernel, coordinates, residual, np.full(5, math.inf), {}
            )

    def test_it_refuses_a_negative_jitter(self) -> None:
        from ampere.core import LikelihoodError

        with pytest.raises(LikelihoodError, match="jitter"):
            DenseGP(jitter=-1.0)

    def test_it_declares_the_four_capability_flags(self) -> None:
        solver = DenseGP()
        assert solver.BACKEND == BACKEND
        assert solver.DIFFERENTIABLE is True
        assert solver.DEVICE == "cpu"

    def test_its_spec_config_matches_the_reference_solver(self) -> None:
        """Why precision and device are class-level rather than fields.

        ``Likelihood.to_spec`` records a dataclass solver's fields as its
        config, and ``results.md`` §14 requires ``ampere_spec_hash`` to agree
        across backends — so an extra field here would make two backends'
        *declaration* of one problem differ.
        """
        import dataclasses

        assert [field.name for field in dataclasses.fields(DenseGP())] == [
            field.name for field in dataclasses.fields(CoreDenseGP())
        ]


class TestQuasisepGP:
    """The O(N) solver: exact against the dense one, and differentiable.

    W2.4 slice 2. Every row here compares against something that is *not*
    another celerite recursion — the dense Cholesky, a closed form, or central
    differences — because two implementations of one recursion agreeing proves
    only that they are the same recursion.
    """

    def test_the_marginal_likelihood_is_the_multivariate_normal(self) -> None:
        """A closed form, not another solver: a genuinely external oracle."""
        kernel, coordinates, residual, variance = gp_case()
        covariance = kernel.matrix(coordinates, coordinates, {}) + np.diag(variance)
        expected = st.multivariate_normal(cov=covariance).logpdf(residual)
        got = QuasisepGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {})
        assert got == pytest.approx(float(expected), abs=1e-6)

    def test_it_agrees_with_the_reference_quasiseparable_solver_to_the_last_bits(self) -> None:
        """The same factorisation, through a different Python interface.

        ``ampere.core.QuasisepGP`` reaches celerite2 through
        ``celerite2.GaussianProcess``; this one calls ``celerite2.backprop``
        directly. Both feed **the same exact rank-2 representation** to the
        same compiled kernel, so the only thing left to differ is the order the
        two sum ``alpha**2/d`` and ``log d`` in — numpy's reduction against
        torch's. The tolerance is therefore an accumulation tolerance, three
        orders tighter than ``tolerances.linear_algebra`` and eight tighter
        than ``cross_solver``: anything looser passing here would mean the two
        term builders had drifted, which is the failure the shared table
        exists to prevent.
        """
        kernel, coordinates, residual, variance = gp_case()
        assert QuasisepGP().log_marginal_likelihood(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreQuasisepGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {}),
            rel=1e-14,
        )

    def test_it_agrees_with_the_dense_solver_on_unsorted_coordinates(self) -> None:
        """Permutation invariance: this solver sorts internally and undoes it."""
        rng = np.random.default_rng(3)
        coordinates = rng.uniform(0.0, 12.0, 40)  # deliberately unsorted
        residual = rng.normal(0.0, 0.3, 40)
        variance = np.full(40, 0.04)
        kernel = Matern32(0.4, 2.0)
        assert QuasisepGP().log_marginal_likelihood(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreDenseGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {}),
            rel=1e-9,
        )

    def test_the_marginal_likelihood_is_differentiable_in_the_residual(self) -> None:
        """``d/dr log N(r; 0, C) = -C**-1 r`` — a closed form, not a difference."""
        kernel, coordinates, residual, variance = gp_case(20)
        residuals = torch.tensor(residual, dtype=torch.float64, requires_grad=True)
        QuasisepGP().log_marginal_likelihood_tensor(
            kernel, coordinates, residuals, variance, {}
        ).backward()
        covariance = kernel.matrix(coordinates, coordinates, {}) + np.diag(variance)
        expected = -np.linalg.solve(covariance, residual)
        assert residuals.grad.numpy() == pytest.approx(expected, abs=1e-8)

    def test_the_marginal_likelihood_is_differentiable_in_the_hyperparameters(self) -> None:
        """W2.4 slice 1's principal finding, closed on the O(N) path too.

        The gradient in the amplitude and the length scale is what a
        torch-backed NUTS over a flexible likelihood actually needs, and it is
        the thing celerite2's *numpy* interface cannot give at any speed: it
        has no reverse pass. Checked against central differences taken through
        the **dense** solver, so neither the value nor the derivative is
        confirmed by the recursion under test.
        """
        _, coordinates, residual, variance = gp_case(40)
        kernel = TorchMatern32(0.4, 2.0)
        amplitude = torch.tensor(0.4, dtype=torch.float64, requires_grad=True)
        length_scale = torch.tensor(2.0, dtype=torch.float64, requires_grad=True)
        QuasisepGP().log_marginal_likelihood_tensor(
            kernel,
            coordinates,
            residual,
            variance,
            {"amplitude": amplitude, "length_scale": length_scale},
        ).backward()

        def dense(a: float, ell: float) -> float:
            k = Matern32(a, ell)
            return CoreDenseGP().log_marginal_likelihood(
                k, coordinates, residual, variance, k.resolve(None)
            )

        step = 1e-6
        assert float(amplitude.grad) == pytest.approx(
            (dense(0.4 + step, 2.0) - dense(0.4 - step, 2.0)) / (2 * step), rel=1e-6
        )
        assert float(length_scale.grad) == pytest.approx(
            (dense(0.4, 2.0 + step) - dense(0.4, 2.0 - step)) / (2 * step), rel=1e-6
        )

    def test_the_leave_one_out_terms_agree_with_the_dense_closed_form(self) -> None:
        """The O(N) recursion W2.3 deferred, checked against a Cholesky.

        ``DenseGP.conditional_loo`` forms ``(K + diag)**-1`` explicitly and
        reads its diagonal; this one accumulates it backwards through the
        semiseparable inverse in O(N). Nothing but agreement would show that
        the accumulation is right — a wrong ``A_ii`` is still finite, still
        per-sample and still looks like a leave-one-out term.
        """
        kernel, coordinates, residual, variance = gp_case(60)
        assert QuasisepGP().conditional_loo(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreDenseGP().conditional_loo(kernel, coordinates, residual, variance, {}), abs=1e-8
        )

    def test_the_leave_one_out_terms_survive_unsorted_coordinates(self) -> None:
        rng = np.random.default_rng(19)
        coordinates = rng.uniform(0.0, 12.0, 45)
        residual = rng.normal(0.0, 0.3, 45)
        variance = np.full(45, 0.04)
        kernel = Matern32(0.4, 2.0)
        assert QuasisepGP().conditional_loo(
            kernel, coordinates, residual, variance, {}
        ) == pytest.approx(
            CoreDenseGP().conditional_loo(kernel, coordinates, residual, variance, {}), abs=1e-8
        )

    def test_the_conditional_agrees_with_the_dense_solver(self) -> None:
        kernel, coordinates, residual, variance = gp_case()
        at = np.linspace(-1.0, 13.0, 41)
        got = QuasisepGP().condition(kernel, coordinates, residual, variance, {}, at=at)
        expected = CoreDenseGP().condition(kernel, coordinates, residual, variance, {}, at=at)
        assert got.mean == pytest.approx(expected.mean, abs=1e-8)
        assert got.variance == pytest.approx(expected.variance, abs=1e-8)

    def test_the_whitening_transform_agrees_with_the_reference_quasiseparable_one(self) -> None:
        """``L`` is not unique, so the comparison is with the same factorisation.

        A dense Cholesky and a celerite ``L√D`` are different square roots of
        one matrix; they agree on ``L Lᵀ = K``, not on ``L z``. So this row
        compares with ``ampere.core.QuasisepGP`` — same factorisation, other
        library — and the *statistical* claim is the next row's.
        """
        kernel, coordinates, _, _ = gp_case()
        whitened = np.random.default_rng(7).normal(size=coordinates.size)
        assert QuasisepGP().latent_transform(kernel, coordinates, whitened, {}) == pytest.approx(
            CoreQuasisepGP().latent_transform(kernel, coordinates, whitened, {}), abs=1e-12
        )

    def test_the_whitening_transform_reproduces_the_kernel_covariance(self) -> None:
        """``L Lᵀ = K``, which is the property that actually matters."""
        kernel, coordinates, _, _ = gp_case(12)
        identity = np.eye(coordinates.size)
        factor = np.stack(
            [QuasisepGP().latent_transform(kernel, coordinates, column, {}) for column in identity]
        ).T
        assert factor @ factor.T == pytest.approx(
            kernel.matrix(coordinates, coordinates, {}), abs=1e-8
        )

    def test_a_singular_covariance_raises_rather_than_returning_a_nan(self) -> None:
        """W2.3's carried finding: celerite2 returns quiet NaN, so guard first."""
        from ampere.core import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        coordinates = np.array([1.0, 1.0, 2.0])
        with pytest.raises(LikelihoodError, match="positive definite"):
            QuasisepGP().log_marginal_likelihood(kernel, coordinates, np.zeros(3), np.zeros(3), {})

    def test_a_non_finite_diagonal_raises_by_name(self) -> None:
        from ampere.core import LikelihoodError

        kernel, coordinates, residual, _ = gp_case(5)
        with pytest.raises(LikelihoodError, match="non-finite"):
            QuasisepGP().log_marginal_likelihood(
                kernel, coordinates, residual, np.full(5, math.inf), {}
            )

    def test_a_negative_variance_raises_by_name(self) -> None:
        from ampere.core import LikelihoodError

        kernel, coordinates, residual, _ = gp_case(5)
        with pytest.raises(LikelihoodError, match="negative entries"):
            QuasisepGP().log_marginal_likelihood(
                kernel, coordinates, residual, np.full(5, -1.0), {}
            )

    @pytest.mark.parametrize(
        "broken",
        ["non-finite", "negative", "singular"],
    )
    def test_the_native_surface_never_raises_and_says_minus_infinity(self, broken: str) -> None:
        """``inference.md`` §10a: a realised density returns ``-inf``, never raises.

        The three ways this solver can fail, and each reaches ``-inf`` by a
        different route: a non-finite or negative diagonal is **sanitised**
        before the compiled kernel sees it (celerite2 would return quiet NaN),
        while a genuinely non-positive-definite matrix — duplicated
        coordinates with no observational noise to separate them — makes
        celerite2 raise, which :func:`ampere.backends.torch._celerite.factor`
        converts to a NaN ``d``.
        """
        kernel, coordinates, residual, variance = gp_case()
        if broken == "non-finite":
            variance = np.full(coordinates.size, math.inf)
        elif broken == "negative":
            variance = np.full(coordinates.size, -1.0)
        else:
            coordinates = np.repeat(coordinates[: coordinates.size // 2], 2)
            variance = np.zeros(coordinates.size)
        value = QuasisepGP().log_marginal_likelihood_native(
            kernel, coordinates, residual, variance, {}
        )
        assert float(value) == -math.inf

    def test_the_native_surface_matches_the_raising_one_where_both_work(self) -> None:
        kernel, coordinates, residual, variance = gp_case()
        solver = QuasisepGP()
        assert float(
            solver.log_marginal_likelihood_native(kernel, coordinates, residual, variance, {})
        ) == solver.log_marginal_likelihood(kernel, coordinates, residual, variance, {})

    def test_the_native_surface_keeps_the_graph_on_a_refusal(self) -> None:
        """pyro's NUTS refuses a potential with no ``grad_fn``.

        A ``torch.where`` keeps the ``-inf`` attached to the graph; a fresh
        constant would not, and the sampler would fail at the first divergent
        proposal rather than at the first bad one.
        """
        kernel, coordinates, residual, _ = gp_case(10)
        residuals = torch.tensor(residual, dtype=torch.float64, requires_grad=True)
        value = QuasisepGP().log_marginal_likelihood_native(
            kernel, coordinates, residuals, np.full(10, -1.0), {}
        )
        assert value.requires_grad
        assert value.grad_fn is not None

    def test_it_refuses_a_kernel_with_no_exact_representation(self) -> None:
        from ampere.core import LikelihoodError

        kernel = TorchSquaredExponential(0.4, 2.0)
        assert not kernel.QUASISEPARABLE
        with pytest.raises(LikelihoodError, match="quasiseparable representation"):
            QuasisepGP().check_compatible(kernel, Spectrum(GRID * u.um, np.ones(GRID.size) * u.Jy))

    def test_it_refuses_a_negative_jitter(self) -> None:
        from ampere.core import LikelihoodError

        with pytest.raises(LikelihoodError, match="jitter"):
            QuasisepGP(jitter=-1.0)

    def test_it_declares_the_flags_the_compiled_kernels_actually_support(self) -> None:
        """CPU float64, declared rather than discovered.

        celerite2's kernels are double-precision CPU C++, so an honest
        ``BATCHABLE``/``DEVICE`` here is ``False``/``"cpu"``. Declaring
        otherwise would make ``FittingProblem.capabilities`` a promise the
        solver cannot keep.
        """
        solver = QuasisepGP()
        assert solver.BACKEND == BACKEND
        assert solver.DIFFERENTIABLE is True
        assert solver.BATCHABLE is False
        assert solver.DEVICE == "cpu"
        assert solver.EXACT is True
        assert solver.IMPLEMENTED is True

    def test_its_spec_config_matches_the_reference_solver(self) -> None:
        """Cross-backend spec-hash agreement: the same dataclass fields, exactly."""
        import dataclasses

        assert [field.name for field in dataclasses.fields(QuasisepGP())] == [
            field.name for field in dataclasses.fields(CoreQuasisepGP())
        ]

    def test_the_solver_configuration_is_recorded_but_never_hashed(self) -> None:
        config = QuasisepGP().provenance_config()
        assert config["dtype"] == "torch.float64"
        assert config["device"] == "cpu"
        assert config["library"] == "celerite2"

    def test_it_scales_linearly_where_the_dense_solver_cannot_run(self) -> None:
        """The claim that makes this solver worth its existence, measured cheaply.

        Not a benchmark — a shared runner has no business asserting a wall
        clock — but a *feasibility* assertion: 20 000 points is 3.2 GB of
        covariance for the dense solver and 480 kB of generators for this one,
        so a finite answer here is itself the demonstration. Agreement is
        asserted at a size the dense solver can reach.
        """
        rng = np.random.default_rng(2026)
        size = 20_000
        coordinates = np.sort(rng.uniform(0.0, size / 10.0, size))
        residual = rng.normal(0.0, 0.3, size)
        variance = np.full(size, 0.01)
        value = QuasisepGP().log_marginal_likelihood(
            Matern32(0.4, 2.0), coordinates, residual, variance, {}
        )
        assert math.isfinite(value)


# ---------------------------------------------------------------------------
# Import discipline and documentation
# ---------------------------------------------------------------------------


class TestDocumentation:
    """The package's own examples, executed."""

    @pytest.mark.parametrize(
        "module",
        [
            "ampere.backends.torch",
            "ampere.backends.torch._config",
            "ampere.backends.torch.lowering",
            "ampere.backends.torch.parameters",
            "ampere.backends.torch.rng",
            "ampere.backends.torch.models",
            "ampere.backends.torch.instrument",
            "ampere.backends.torch.gp",
            "ampere.backends.torch._celerite",
        ],
    )
    def test_the_module_docstring_examples_run(self, module: str) -> None:
        import doctest
        import importlib

        results = doctest.testmod(
            importlib.import_module(module),
            optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
            verbose=False,
        )
        assert results.failed == 0

    def test_the_model_base_class_is_a_core_model(self) -> None:
        """A backend supplies models and transformations, and nothing else."""
        assert issubclass(PowerLaw, Model)
        assert issubclass(Resample, Transformation)
