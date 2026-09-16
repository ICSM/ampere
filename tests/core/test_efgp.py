"""``ampere.core.efgp``: the equispaced-Fourier prototype W5.6 measures.

The assertions are in **W5.4's comparison class for an ``EXACT = False``
solver** (``tests/conformance/protocol.py``'s ``approximation_envelope`` and
``tests/conformance/README.md`` §3), not at a fixed tolerance: an
approximation's agreement with :class:`~ampere.core.DenseGP` depends on the
kernel, the data and its own two parameters, so what is asserted is that the
error **falls as the grid is refined**, inside an envelope anchored on the
coarsest setting of the same sweep and floored at the periodisation the
spacing fixes — and that the finest setting is actually close.

The envelope is re-stated here rather than imported because ``tests/`` has no
package root and ``tests/conformance`` is not on this module's import path;
:func:`_envelope` is ``approximation_envelope``'s formula and nothing else.
Rows that hold the *exact* structure of the method — the Toeplitz identity,
the circulant matvec, the real feature pair — are ordinary equalities, because
those are not approximations at all.
"""

from __future__ import annotations

import doctest

import astropy.units as u
import numpy as np
import pytest

import ampere.core.efgp as efgp_module
from ampere.core import (
    DenseGP,
    EquispacedFourierGP,
    GaussianFamily,
    GaussianProcessNoise,
    Matern12,
    Matern32,
    Matern52,
    Product,
    Spectrum,
    SquaredExponential,
)
from ampere.core.efgp import (
    fourier_grid,
    real_features,
    solve_iterative,
    spectral_weights,
    toeplitz_generator,
    toeplitz_matrix,
    toeplitz_matvec,
)
from ampere.core.exceptions import LikelihoodError

#: The floor of the envelope: the box's own truncation, as W5.4's class names
#: it. Here it is the rectangle rule's periodisation, which the spacing fixes
#: and no amount of extra grid removes.
APPROXIMATION_FLOOR = 1e-2


def _envelope(coarsest_error: float, coarsest_size: int, size: int, order: float) -> float:
    """``approximation_envelope``'s formula (``tests/conformance/protocol.py``)."""
    if size <= coarsest_size:
        return float(coarsest_error)
    return float(coarsest_error) * (coarsest_size / size) ** order + APPROXIMATION_FLOOR


def _spectrum() -> Spectrum:
    """A minimal one-axis container, for the composition-time refusals."""
    grid = np.linspace(1.0, 4.0, 40)
    return Spectrum(grid * u.um, np.zeros(40) * u.Jy, uncertainty=np.full(40, 0.1) * u.Jy)


@pytest.fixture
def problem() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """An irregular 1-D problem whose residual is a draw from the process."""
    rng = np.random.default_rng(20260916)
    points = np.sort(rng.uniform(0.0, 12.0, 220))[:, None]
    kernel = Matern32(0.4, 1.5)
    covariance = np.asarray(kernel.matrix(points, points, kernel.resolve(None)))
    lower = np.linalg.cholesky(covariance + 1e-10 * np.eye(points.shape[0]))
    residual = lower @ rng.normal(size=points.shape[0]) + rng.normal(0.0, 0.05, points.shape[0])
    return points, residual, np.full(points.shape[0], 0.05**2)


# ---------------------------------------------------------------------------
# The declaration
# ---------------------------------------------------------------------------


def test_flags_and_declaration() -> None:
    solver = EquispacedFourierGP(basis_size=(17, 17), boundary_factor=3.0)
    assert solver.EXACT is False
    assert solver.IMPLEMENTED is True
    assert solver.STACKED_RESIDUALS is True
    assert solver.DIFFERENTIABLE is False, "the prototype is the reference path only"
    assert solver.counts == (17, 17)
    assert solver.latent_size(Matern32(0.3, 2.0), 10_000) == 289
    assert solver.provenance_config() == {"basis_size": [17, 17], "boundary_factor": 3.0}


@pytest.mark.parametrize(
    ("kwargs", "fragment"),
    [
        ({"basis_size": 32}, "must be odd on every axis"),
        ({"boundary_factor": 0.0}, "boundary_factor must be finite and > 0"),
        ({"jitter": -1.0}, "jitter must be finite and >= 0"),
    ],
)
def test_declaration_refusals(kwargs: dict[str, object], fragment: str) -> None:
    with pytest.raises(LikelihoodError, match=fragment):
        EquispacedFourierGP(**kwargs)  # type: ignore[arg-type]


def test_refuses_a_kernel_with_no_spectral_density() -> None:
    data = _spectrum()
    product = Product(Matern32(0.3, 2.0), Matern12(0.2, 0.5))
    noise = GaussianProcessNoise(product, EquispacedFourierGP())
    with pytest.raises(LikelihoodError, match="closed-form spectral density"):
        noise.check_compatible(GaussianFamily(), data)


def test_refuses_a_count_per_axis_mismatch() -> None:
    data = _spectrum()
    noise = GaussianProcessNoise(Matern32(0.3, 2.0), EquispacedFourierGP(basis_size=(9, 9)))
    with pytest.raises(LikelihoodError, match="basis_size"):
        noise.check_compatible(GaussianFamily(), data)


def test_refuses_a_zero_noise_diagonal(problem: tuple[np.ndarray, ...]) -> None:
    points, residual, variance = problem
    kernel = Matern32(0.4, 1.5)
    with pytest.raises(LikelihoodError, match="strictly positive noise diagonal"):
        EquispacedFourierGP(basis_size=17).log_marginal_likelihood(
            kernel, points, residual, np.zeros_like(variance), kernel.resolve(None)
        )


# ---------------------------------------------------------------------------
# The structure: exact identities, asserted exactly
# ---------------------------------------------------------------------------


def test_toeplitz_generator_reproduces_the_normal_matrix(
    problem: tuple[np.ndarray, ...],
) -> None:
    """The claim the whole method rests on: ``Phi^H D^-1 Phi`` is Toeplitz.

    Formed both ways — the ``O(N m²)`` product HSGP would compute, and the
    generator EFGP assembles at ``O(N 2^d m)`` — and required to agree.
    """
    points, _, variance = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    grid = fourier_grid(points, (33,), 4.0)
    weights = spectral_weights(kernel, grid, values)
    root = np.sqrt(weights)
    features = np.exp(1j * ((points - grid.centre[None, :]) @ grid.frequencies.T)) * root[None, :]
    direct = features.conj().T @ (features / variance[:, None])
    structured = toeplitz_matrix(grid, toeplitz_generator(grid, points, 1.0 / variance))
    structured = structured * root[:, None] * root[None, :]
    assert np.max(np.abs(direct - structured)) < 1e-9


def test_circulant_embedding_matches_the_dense_product(
    problem: tuple[np.ndarray, ...],
) -> None:
    points, _, variance = problem
    kernel = Matern32(0.4, 1.5)
    grid = fourier_grid(points, (33,), 4.0)
    generator = toeplitz_generator(grid, points, 1.0 / variance)
    dense = toeplitz_matrix(grid, generator)
    rng = np.random.default_rng(7)
    vector = rng.normal(size=grid.size) + 1j * rng.normal(size=grid.size)
    assert np.max(np.abs(dense @ vector - toeplitz_matvec(grid, generator, vector))) < 1e-8
    del kernel


def test_conjugate_gradients_solve_the_normal_equations(
    problem: tuple[np.ndarray, ...],
) -> None:
    """``solve_iterative`` is the ``O(m log m)`` half of the published claim."""
    points, _, variance = problem
    kernel = Matern32(0.4, 1.5)
    grid = fourier_grid(points, (33,), 4.0)
    weights = spectral_weights(kernel, grid, kernel.resolve(None))
    generator = toeplitz_generator(grid, points, 1.0 / variance)
    root = np.sqrt(weights)
    matrix = toeplitz_matrix(grid, generator) * root[:, None] * root[None, :]
    matrix[np.diag_indices_from(matrix)] += 1.0
    rng = np.random.default_rng(11)
    right = rng.normal(size=grid.size) + 1j * rng.normal(size=grid.size)
    solution, iterations = solve_iterative(grid, generator, root, right)
    assert iterations < grid.size, "CG on I + PSD must not need a full basis of iterations"
    assert np.max(np.abs(matrix @ solution - right)) < 1e-7


def test_real_features_factorise_the_represented_covariance(
    problem: tuple[np.ndarray, ...],
) -> None:
    """``Psi Psi^T`` is the same covariance the complex grid represents.

    The two representations are what ``latent_transform`` and
    ``log_marginal_likelihood`` respectively use, so ``simulate(observe=True)``
    and the latent path agree only if this holds — exactly, not approximately.
    """
    points, _, _ = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    grid = fourier_grid(points, (33,), 4.0)
    weights = spectral_weights(kernel, grid, values)
    features = real_features(grid, weights, points)
    assert features.shape == (points.shape[0], grid.size)
    separation = points - points.T
    represented = np.real(np.exp(1j * separation[:, :, None] * grid.frequencies[None, None, :, 0]))
    represented = represented @ weights
    assert np.max(np.abs(features @ features.T - represented)) < 1e-10


def test_latent_transform_size_and_refusal(problem: tuple[np.ndarray, ...]) -> None:
    points, _, _ = problem
    kernel = Matern32(0.4, 1.5)
    solver = EquispacedFourierGP(basis_size=17, boundary_factor=4.0)
    draws = np.random.default_rng(2).normal(size=17)
    field = solver.latent_transform(kernel, points, draws, kernel.resolve(None))
    assert field.shape == (points.shape[0],)
    with pytest.raises(LikelihoodError, match="Fourier coefficient"):
        solver.latent_transform(kernel, points, np.zeros(5), kernel.resolve(None))


def test_stacked_residuals_are_the_sum_of_their_columns(
    problem: tuple[np.ndarray, ...],
) -> None:
    points, residual, variance = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    solver = EquispacedFourierGP(basis_size=33, boundary_factor=4.0)
    block = np.column_stack([residual, 0.5 * residual])
    joint = solver.log_marginal_likelihood(kernel, points, block, variance, values)
    apart = solver.log_marginal_likelihood(
        kernel, points, residual, variance, values
    ) + solver.log_marginal_likelihood(kernel, points, 0.5 * residual, variance, values)
    assert joint == pytest.approx(apart, rel=1e-12)


# ---------------------------------------------------------------------------
# The approximation: a convergence, not a number
# ---------------------------------------------------------------------------


#: ``(family, refinement order, final tolerance)``. The rate is a property of
#: the *kernel*: an isotropic Matérn-nu's spectral density decays as
#: ``omega**-(2 nu + d)``, so the truncated tail — and with it the reduced-rank
#: error — falls as ``m**-2nu``, and a squared exponential's falls
#: exponentially. Matérn-1/2 is the roughest process a spectral method can be
#: asked for and its own row says so.
SWEEP = [
    # Matérn-1/2's row states its own rate, **looser** than the default 1.0,
    # and the reason is the bake-off's headline rather than a fudge: the
    # ``m ** -2nu`` law is the rate of the *kernel's* truncated tail, and the
    # likelihood error inherits it only asymptotically. The measured rate on
    # this fixture is 0.84, and the finest setting here is still 8.5 nats from
    # the exact answer — a spectral method does not converge usefully on a
    # rough process at any size worth paying for, which is precisely why
    # ``VecchiaResponseGP`` is in this bake-off beside it.
    (Matern12, 0.75, 2.0e1),
    (Matern32, 3.0, 1.0e-2),
    (Matern52, 5.0, 1.0e-4),
    (SquaredExponential, 5.0, 1.0e-9),
]
SIZES = (65, 129, 257, 513)


@pytest.mark.parametrize(("family", "order", "final"), SWEEP, ids=lambda value: str(value))
def test_marginal_likelihood_converges_to_dense(
    problem: tuple[np.ndarray, ...], family: type, order: float, final: float
) -> None:
    points, residual, variance = problem
    kernel = family(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().log_marginal_likelihood(kernel, points, residual, variance, values)
    errors = [
        abs(
            EquispacedFourierGP(basis_size=size, boundary_factor=2.0).log_marginal_likelihood(
                kernel, points, residual, variance, values
            )
            - exact
        )
        for size in SIZES
    ]
    for size, error in zip(SIZES[1:], errors[1:], strict=True):
        assert error <= _envelope(errors[0], SIZES[0], size, order), (
            f"{family.__name__} at m={size}: {error:.3e} outside the envelope anchored on "
            f"{errors[0]:.3e} at m={SIZES[0]}"
        )
    assert errors[-1] < final


def test_conditional_moments_converge_to_dense(problem: tuple[np.ndarray, ...]) -> None:
    points, residual, variance = problem
    kernel = Matern52(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().condition(kernel, points, residual, variance, values)
    previous_mean = previous_variance = np.inf
    for size in SIZES:
        conditional = EquispacedFourierGP(basis_size=size, boundary_factor=2.0).condition(
            kernel, points, residual, variance, values
        )
        mean_error = float(np.max(np.abs(conditional.mean - exact.mean)))
        variance_error = float(np.max(np.abs(conditional.variance - exact.variance)))
        assert np.all(conditional.variance >= 0.0), "a variance is never negative at any m"
        assert mean_error <= previous_mean + APPROXIMATION_FLOOR
        assert variance_error <= previous_variance + APPROXIMATION_FLOOR
        previous_mean, previous_variance = mean_error, variance_error
    assert previous_mean < 1e-6
    assert previous_variance < 1e-8


def test_leave_one_out_terms_converge_to_dense(problem: tuple[np.ndarray, ...]) -> None:
    points, residual, variance = problem
    kernel = Matern52(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().conditional_loo(kernel, points, residual, variance, values)
    coarse = EquispacedFourierGP(basis_size=65, boundary_factor=2.0).conditional_loo(
        kernel, points, residual, variance, values
    )
    fine = EquispacedFourierGP(basis_size=513, boundary_factor=2.0).conditional_loo(
        kernel, points, residual, variance, values
    )
    assert fine.shape == exact.shape
    assert np.max(np.abs(fine - exact)) < np.max(np.abs(coarse - exact))
    assert np.max(np.abs(fine - exact)) < 1e-4


def test_the_two_spectral_solvers_agree_where_both_converge(
    problem: tuple[np.ndarray, ...],
) -> None:
    """EFGP and HSGP represent the same covariance by different quadratures.

    They are not required to agree at small ``m`` — different bases truncate
    differently — but at a grid fine enough for both, they must, because they
    converge to the same exact answer.
    """
    from ampere.core import HilbertSpaceGP

    points, residual, variance = problem
    kernel = Matern52(0.4, 1.5)
    values = kernel.resolve(None)
    fourier = EquispacedFourierGP(basis_size=513, boundary_factor=2.0).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    hilbert = HilbertSpaceGP(basis_size=513, boundary_factor=2.0).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    assert fourier == pytest.approx(hilbert, abs=1e-3)


def test_module_docstring_examples_run() -> None:
    results = doctest.testmod(efgp_module, optionflags=doctest.ELLIPSIS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


def test_circulant_embedding_in_two_axes() -> None:
    """The 2-D embedding has four corners, and a 1-D test cannot see the other two.

    A multilevel circulant wraps each axis independently, so an embedding that
    filled only the ``(forward, forward)`` and ``(backward, backward)`` blocks
    is exactly right in one axis and wrong in two — and wrong in a way the
    marginal likelihood, which factorises the explicit matrix, never notices:
    only ``condition``'s mean and ``solve_iterative`` go through the FFT. That
    is how it got past the 1-D row above once.
    """
    rng = np.random.default_rng(13)
    points = rng.uniform(-1.0, 1.0, (120, 2))
    variance = np.full(120, 0.04)
    kernel = Matern32(0.4, 0.5, axes=None)
    grid = fourier_grid(points, (7, 5), 2.0)
    generator = toeplitz_generator(grid, points, 1.0 / variance)
    dense = toeplitz_matrix(grid, generator)
    vector = rng.normal(size=grid.size) + 1j * rng.normal(size=grid.size)
    assert np.max(np.abs(dense @ vector - toeplitz_matvec(grid, generator, vector))) < 1e-8
    del kernel


def test_conditional_mean_in_two_axes_matches_dense() -> None:
    """``condition`` takes the FFT route for its mean, so it needs its own 2-D row."""
    rng = np.random.default_rng(17)
    points = rng.uniform(-1.0, 1.0, (200, 2))
    kernel = Matern32(0.4, 0.5)
    values = kernel.resolve(None)
    covariance = np.asarray(kernel.matrix(points, points, values))
    lower = np.linalg.cholesky(covariance + 1e-10 * np.eye(200))
    residual = lower @ rng.normal(size=200) + rng.normal(0.0, 0.05, 200)
    variance = np.full(200, 0.05**2)
    from ampere.core import HilbertSpaceGP

    exact = DenseGP().condition(kernel, points, residual, variance, values)
    coarse = EquispacedFourierGP(basis_size=(9, 9), boundary_factor=2.0).condition(
        kernel, points, residual, variance, values
    )
    approximate = EquispacedFourierGP(basis_size=(33, 33), boundary_factor=2.0).condition(
        kernel, points, residual, variance, values
    )
    hilbert = HilbertSpaceGP(basis_size=(33, 33), boundary_factor=2.0).condition(
        kernel, points, residual, variance, values
    )
    assert np.max(np.abs(approximate.mean - exact.mean)) < np.max(np.abs(coarse.mean - exact.mean))
    assert np.max(np.abs(approximate.mean - exact.mean)) < 2e-2
    assert np.max(np.abs(approximate.variance - exact.variance)) < 1e-3
    # The sharper statement, and the one the broken embedding failed by two
    # orders of magnitude: the two reduced-rank solvers represent the same
    # covariance, so at equal resolution their conditional means agree far
    # more closely than either agrees with the exact one.
    assert np.max(np.abs(approximate.mean - hilbert.mean)) < 0.2 * np.max(
        np.abs(approximate.mean - exact.mean)
    )
