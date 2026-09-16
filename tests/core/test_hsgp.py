"""The Hilbert-space basis and the reduced-rank solver (W5.4).

What is asserted here rather than in the conformance battery: the properties of
``ampere.core``'s own basis construction, which is shared by all three backends
and so has no cross-backend column to compare against, and the closed-form
spectral densities, whose oracle is a **numerical Fourier transform** of the
kernel's own closed form. That oracle is the point — a spectral density
transcribed with a wrong constant is still smooth, still positive and still
converges the approximation to *something*, and only an independent transform
catches it.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.integrate
import scipy.stats as st

from ampere.core import (
    SHO,
    DenseGP,
    HilbertSpaceGP,
    LikelihoodError,
    Matern12,
    Matern32,
    Matern52,
    Product,
    RotationTerm,
    SpectralMixture,
    SquaredExponential,
    Sum,
)
from ampere.core.hsgp import (
    DEFAULT_BASIS_SIZE,
    DEFAULT_BOUNDARY_FACTOR,
    basis_matrix,
    hilbert_basis,
    normalise_counts,
    spectral_values,
)

GRID = np.array([1.0, 1.4, 1.9, 2.5, 3.2, 4.0, 5.0, 6.1, 7.3, 8.6, 10.0, 11.5]).reshape(-1, 1)
SIGMA = 0.1
AMPLITUDE = 0.4
LENGTH_SCALE = 2.0


def residuals(seed: int = 1) -> np.ndarray:
    return np.random.default_rng(seed).normal(0.0, 0.3, GRID.shape[0])


def variances() -> np.ndarray:
    return np.full(GRID.shape[0], SIGMA**2)


def numerical_spectral_density(kernel, frequency: float, values=None) -> float:
    """``∫ k(|tau|) e^{-i omega tau} d tau`` in 1-D, by quadrature.

    The oracle. It calls the kernel's own ``value`` — its closed form in one
    separation — and nothing else of ampere's, so it is a transform of the
    covariance rather than a second copy of the spectral density.
    """
    resolved = kernel.resolve({} if values is None else values)

    def integrand(tau: float) -> float:
        return float(np.asarray(kernel.value(np.array([tau]), resolved))[0]) * math.cos(
            frequency * tau
        )

    total, _ = scipy.integrate.quad(integrand, 0.0, 400.0, limit=4000)
    return 2.0 * total


class TestSpectralDensities:
    """Each closed form against a numerical transform of the same kernel."""

    @pytest.mark.parametrize(
        "kernel",
        [
            Matern12(AMPLITUDE, LENGTH_SCALE),
            Matern32(AMPLITUDE, LENGTH_SCALE),
            Matern52(AMPLITUDE, LENGTH_SCALE),
            SquaredExponential(AMPLITUDE, LENGTH_SCALE),
            SHO(0.5, 1.5, 3.0),
            Sum(Matern32(0.3, 3.0), Matern12(0.15, 0.6)),
        ],
        ids=["matern12", "matern32", "matern52", "squared_exponential", "sho", "sum"],
    )
    @pytest.mark.parametrize("frequency", [0.0, 0.3, 1.0, 2.5, 6.0])
    def test_the_closed_form_is_the_fourier_transform_of_the_covariance(
        self, kernel, frequency: float
    ) -> None:
        got = float(
            np.asarray(kernel.spectral_density(np.array([frequency]), kernel.resolve({})))[0]
        )
        expected = numerical_spectral_density(kernel, frequency)
        assert got == pytest.approx(expected, rel=1e-6, abs=1e-12)

    def test_it_is_non_negative_everywhere_which_is_what_makes_it_a_covariance(self) -> None:
        """Bochner's theorem, asserted rather than assumed: a negative ``S`` would
        make ``Phi diag(S) Phi^T`` not a covariance at all, and the square root the
        solver takes of it would be NaN."""
        frequencies = np.linspace(0.0, 60.0, 601)
        for kernel in (
            Matern12(AMPLITUDE, LENGTH_SCALE),
            Matern32(AMPLITUDE, LENGTH_SCALE),
            Matern52(AMPLITUDE, LENGTH_SCALE),
            SquaredExponential(AMPLITUDE, LENGTH_SCALE),
            SHO(0.5, 1.5, 0.51),
            SHO(0.5, 1.5, 40.0),
        ):
            density = np.asarray(kernel.spectral_density(frequencies, kernel.resolve({})))
            assert np.all(density >= 0.0)

    def test_the_dimension_is_part_of_the_answer(self) -> None:
        """The same covariance has a different spectral density in one axis and
        in two, and a solver that forgot the dimension would build a basis for
        the wrong process rather than fail."""
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        one = float(np.asarray(kernel.spectral_density(np.array([1.0]), kernel.resolve({})))[0])
        two = float(
            np.asarray(kernel.spectral_density(np.array([1.0]), kernel.resolve({}), dimensions=2))[
                0
            ]
        )
        assert one != pytest.approx(two, rel=1e-3)

    def test_the_squared_exponential_integrates_back_to_its_own_marginal_variance(self) -> None:
        """``k(0) = (2 pi)^-d integral S``, the inverse transform at zero separation."""
        kernel = SquaredExponential(AMPLITUDE, LENGTH_SCALE)
        resolved = kernel.resolve({})
        total, _ = scipy.integrate.quad(
            lambda w: float(np.asarray(kernel.spectral_density(np.array([w]), resolved))[0]),
            -60.0,
            60.0,
            limit=4000,
        )
        assert total / (2.0 * math.pi) == pytest.approx(AMPLITUDE**2, rel=1e-9)

    @pytest.mark.parametrize(
        ("kernel", "named"),
        [
            (Product(Matern32(0.3, 2.0), Matern12(0.2, 0.5)), "Product"),
            (RotationTerm(0.4, 2.0, 3.0, 0.5, 0.3), "RotationTerm"),
            (SpectralMixture(0.4, 2.0, 1.0), "SpectralMixture"),
        ],
        ids=["product", "rotation", "spectral_mixture"],
    )
    def test_a_family_without_one_refuses_by_name(self, kernel, named: str) -> None:
        with pytest.raises(LikelihoodError) as refusal:
            kernel.spectral_density(np.array([1.0]), kernel.resolve({}))
        assert named in str(refusal.value)
        assert "spectral density" in str(refusal.value)

    def test_an_sho_refuses_more_than_one_axis(self) -> None:
        kernel = SHO(0.5, 1.5, 3.0)
        with pytest.raises(LikelihoodError) as refusal:
            kernel.spectral_density(np.array([1.0]), kernel.resolve({}), dimensions=2)
        assert "one ordered coordinate" in str(refusal.value)

    def test_a_stationary_subclass_without_a_smoothness_declaration_refuses(self) -> None:
        """A user's own ``StationaryKernel`` must not silently inherit the Matérn
        form: the shared constructor gives it the two hyperparameters, not a
        spectral density."""

        class Triangular(SquaredExponential.__mro__[1]):  # StationaryKernel
            FAMILY = "triangular"
            HYPERPARAMETERS = ("amplitude", "length_scale")
            QUASISEPARABLE = False

            def _covariance(self, separation, values):
                amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
                length_scale = self._hyperparameter(values, "length_scale")
                scaled = self.ops.scalar(separation) / length_scale
                return amplitude * amplitude * np.clip(1.0 - scaled, 0.0, None)

        kernel = Triangular(0.4, 2.0)
        with pytest.raises(LikelihoodError) as refusal:
            kernel.spectral_density(np.array([1.0]), kernel.resolve({}))
        assert "triangular" in str(refusal.value)


class TestTheBasis:
    """The box, the frequencies and the eigenfunctions."""

    def test_the_box_contains_the_data_with_the_declared_room_to_spare(self) -> None:
        basis = hilbert_basis(GRID, (16,), 2.0)
        centre = float(basis.centre[0])
        width = float(basis.half_width[0])
        assert centre == pytest.approx(0.5 * (GRID.min() + GRID.max()))
        assert width == pytest.approx(2.0 * 0.5 * (GRID.max() - GRID.min()))
        assert np.all(np.abs(GRID.ravel() - centre) < width)

    def test_the_eigenfunctions_are_orthonormal_on_the_box(self) -> None:
        """``integral phi_i phi_j = delta_ij`` over ``[-L, L]``, which is what makes
        ``Phi diag(S) Phi^T`` the right factorisation rather than merely a
        plausible one."""
        basis = hilbert_basis(GRID, (12,), 2.0)
        centre = float(basis.centre[0])
        width = float(basis.half_width[0])
        fine = np.linspace(centre - width, centre + width, 20001).reshape(-1, 1)
        matrix = np.asarray(basis_matrix(basis, fine))
        gram = np.trapezoid(matrix[:, :, None] * matrix[:, None, :], fine.ravel(), axis=0)
        assert np.max(np.abs(gram - np.eye(12))) < 1e-6

    def test_the_frequencies_are_the_laplacian_eigenvalues_of_the_box(self) -> None:
        basis = hilbert_basis(GRID, (5,), 1.5)
        width = float(basis.half_width[0])
        expected = math.pi * np.arange(1, 6) / (2.0 * width)
        assert np.allclose(basis.norms, expected)
        assert basis.size == 5
        assert basis.dimensions == 1

    def test_the_tensor_product_orders_its_columns_with_its_frequencies(self) -> None:
        """Column ``j`` of ``Phi`` and row ``j`` of ``frequencies`` must be the same
        basis member, or the spectral density is applied to the wrong column."""
        points = np.column_stack([np.linspace(-1.0, 1.0, 7), np.linspace(-2.0, 2.0, 7)[::-1]])
        basis = hilbert_basis(points, (3, 4), 2.0)
        assert basis.size == 12
        matrix = np.asarray(basis_matrix(basis, points))
        assert matrix.shape == (7, 12)
        for index in range(12):
            frequency = basis.frequencies[index]
            direct = np.ones(7)
            for axis in range(2):
                width = float(basis.half_width[axis])
                centre = float(basis.centre[axis])
                direct = (
                    direct
                    * np.sin(frequency[axis] * (points[:, axis] - centre + width))
                    / math.sqrt(width)
                )
            assert np.allclose(matrix[:, index], direct)

    def test_a_degenerate_axis_gets_a_unit_extent_rather_than_a_zero_box(self) -> None:
        flat = np.zeros((5, 1))
        basis = hilbert_basis(flat, (4,), 2.0)
        assert float(basis.half_width[0]) == pytest.approx(2.0)

    def test_a_count_per_axis_is_required(self) -> None:
        with pytest.raises(LikelihoodError, match="axis count"):
            hilbert_basis(GRID, (8, 8), 2.0)

    def test_the_counts_declaration_is_validated(self) -> None:
        assert normalise_counts(32) == (32,)
        assert normalise_counts([4, 5]) == (4, 5)
        with pytest.raises(LikelihoodError, match="at least one axis"):
            normalise_counts([])
        with pytest.raises(LikelihoodError, match="positive"):
            normalise_counts([4, 0])

    def test_the_spectral_values_are_the_density_at_those_frequencies(self) -> None:
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        basis = hilbert_basis(GRID, (6,), 2.0)
        got = np.asarray(spectral_values(kernel, basis, kernel.resolve({})))
        expected = np.asarray(kernel.spectral_density(basis.norms, kernel.resolve({})))
        assert np.allclose(got, expected)


class TestTheSolver:
    """``HilbertSpaceGP``'s own declarations and internal consistency."""

    def test_the_declaration_is_visible_and_normalised(self) -> None:
        solver = HilbertSpaceGP()
        assert solver.basis_size == (DEFAULT_BASIS_SIZE,)
        assert solver.boundary_factor == DEFAULT_BOUNDARY_FACTOR
        assert not solver.EXACT
        assert solver.IMPLEMENTED
        assert HilbertSpaceGP(basis_size=(8, 8)).latent_size(Matern32(0.3, 2.0), 10_000) == 64

    def test_the_approximation_reaches_the_spec_and_so_the_hash(self) -> None:
        """``m`` and ``c`` change the answer, so two runs at different ``m`` are
        not the same declared model and must not hash alike."""
        import dataclasses

        fields = [field.name for field in dataclasses.fields(HilbertSpaceGP())]
        assert fields == ["basis_size", "boundary_factor", "jitter"]

    @pytest.mark.parametrize(
        ("argument", "message"),
        [
            ({"boundary_factor": 0.0}, "boundary_factor"),
            ({"boundary_factor": math.inf}, "boundary_factor"),
            ({"jitter": -1.0}, "jitter"),
        ],
    )
    def test_an_inadmissible_declaration_is_refused(self, argument, message: str) -> None:
        with pytest.raises(LikelihoodError, match=message):
            HilbertSpaceGP(**argument)

    def test_the_marginal_likelihood_is_the_dense_one_of_its_own_approximation(self) -> None:
        """The solver's Woodbury algebra against ``scipy``, on the covariance the
        solver's *own* whitening defines. This is the row that says the
        factorisation is right independently of how good the approximation is."""
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        values = kernel.resolve({})
        solver = HilbertSpaceGP(basis_size=24, boundary_factor=2.0)
        factor = solver.latent_transform(kernel, GRID, np.eye(24), values)
        covariance = factor @ factor.T + np.diag(variances())
        oracle = float(
            st.multivariate_normal.logpdf(residuals(), mean=np.zeros(GRID.shape[0]), cov=covariance)
        )
        got = solver.log_marginal_likelihood(kernel, GRID, residuals(), variances(), values)
        assert got == pytest.approx(oracle, abs=1e-8)

    def test_the_leave_one_out_terms_are_the_dense_identity_on_that_covariance(self) -> None:
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        values = kernel.resolve({})
        solver = HilbertSpaceGP(basis_size=24, boundary_factor=2.0)
        factor = solver.latent_transform(kernel, GRID, np.eye(24), values)
        precision = np.linalg.inv(factor @ factor.T + np.diag(variances()))
        alpha = precision @ residuals()
        expected = (
            0.5 * np.log(np.diag(precision))
            - 0.5 * math.log(2.0 * math.pi)
            - alpha**2 / (2.0 * np.diag(precision))
        )
        got = solver.conditional_loo(kernel, GRID, residuals(), variances(), values)
        assert np.allclose(got, expected, atol=1e-9)

    def test_the_conditioned_variance_is_non_negative_at_every_basis_size(self) -> None:
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        values = kernel.resolve({})
        at = np.linspace(float(GRID.min()) - 2.0, float(GRID.max()) + 2.0, 61).reshape(-1, 1)
        for size in (4, 8, 16, 64):
            solver = HilbertSpaceGP(basis_size=size, boundary_factor=2.0)
            conditioned = solver.condition(kernel, GRID, residuals(), variances(), values, at=at)
            assert np.all(np.asarray(conditioned.variance) >= 0.0)

    def test_the_error_falls_as_the_basis_is_refined(self) -> None:
        """The convergence claim, in its simplest form. The conformance battery
        states it as a tolerance class over every fixture; this is the reference
        path's own sanity check."""
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        values = kernel.resolve({})
        exact = DenseGP().log_marginal_likelihood(kernel, GRID, residuals(), variances(), values)
        errors = [
            abs(
                HilbertSpaceGP(basis_size=size, boundary_factor=2.0).log_marginal_likelihood(
                    kernel, GRID, residuals(), variances(), values
                )
                - exact
            )
            for size in (16, 32, 64, 128)
        ]
        assert errors[0] > errors[1] > errors[2] > errors[3]
        assert errors[-1] < 1e-2

    def test_a_zero_noise_diagonal_is_refused_rather_than_divided_by(self) -> None:
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        with pytest.raises(LikelihoodError, match="strictly positive noise diagonal"):
            HilbertSpaceGP().log_marginal_likelihood(
                kernel, GRID, residuals(), np.zeros(GRID.shape[0]), kernel.resolve({})
            )

    def test_a_whitened_block_of_the_wrong_length_is_refused_by_name(self) -> None:
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        solver = HilbertSpaceGP(basis_size=16)
        with pytest.raises(LikelihoodError, match="basis size m, not the sample count"):
            solver.latent_transform(kernel, GRID, np.zeros(GRID.shape[0]), kernel.resolve({}))

    def test_a_stacked_right_hand_side_is_the_sum_of_its_columns(self) -> None:
        """W4.2's rule, on this solver: ``k`` independent realisations sharing one
        covariance, so the joint marginal is the sum of the columns'."""
        kernel = Matern32(AMPLITUDE, LENGTH_SCALE)
        values = kernel.resolve({})
        solver = HilbertSpaceGP(basis_size=24)
        left, right = residuals(1), residuals(2)
        joint = solver.log_marginal_likelihood(
            kernel, GRID, np.column_stack([left, right]), variances(), values
        )
        separate = sum(
            solver.log_marginal_likelihood(kernel, GRID, column, variances(), values)
            for column in (left, right)
        )
        assert joint == pytest.approx(separate, abs=1e-9)
