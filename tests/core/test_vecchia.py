"""``ampere.core.vecchia``: the sparse prototype W5.6 measures against the spectral two.

Held to the same **convergence** class as ``test_efgp.py``'s rows
(``tests/conformance/protocol.py``'s ``approximation_envelope``,
``tests/conformance/README.md`` §3) — the error must fall as the conditioning
sets grow, and the largest ``k`` must actually be close — with two differences
that are the method's own and are asserted as such:

* the refinement parameter is ``k``, not a basis size, and the rate is not the
  spectral ``m ** -2nu``; and
* the ordering is part of the approximation (Guinness 2018), so it is swept
  rather than fixed.

The rows that matter most to the bake-off are the two at the ends of the
smoothness range: :func:`test_converges_where_a_spectral_method_cannot`, which
is Vecchia on a Matérn-1/2, and
:func:`test_screening_is_weak_for_a_smooth_process`, which is the price it pays
for it.
"""

from __future__ import annotations

import doctest

import astropy.units as u
import numpy as np
import pytest

import ampere.core.vecchia as vecchia_module
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
    VecchiaResponseGP,
)
from ampere.core.exceptions import LikelihoodError
from ampere.core.vecchia import neighbour_structure, separation_covariance

APPROXIMATION_FLOOR = 1e-2


def _envelope(coarsest_error: float, coarsest_size: int, size: int, order: float) -> float:
    """``approximation_envelope``'s formula (``tests/conformance/protocol.py``)."""
    if size <= coarsest_size:
        return float(coarsest_error)
    return float(coarsest_error) * (coarsest_size / size) ** order + APPROXIMATION_FLOOR


WIDTHS = (4, 8, 16, 32)


def _spectrum() -> Spectrum:
    grid = np.linspace(1.0, 4.0, 40)
    return Spectrum(grid * u.um, np.zeros(40) * u.Jy, uncertainty=np.full(40, 0.1) * u.Jy)


@pytest.fixture
def problem() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``test_efgp.py``'s problem, seed included, so the two files compare."""
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
    solver = VecchiaResponseGP(neighbours=20, ordering="coordinate", seed=3)
    assert solver.EXACT is False
    assert solver.IMPLEMENTED is True
    assert solver.STACKED_RESIDUALS is True
    assert solver.DIFFERENTIABLE is False, "the prototype is the reference path only"
    assert solver.provenance_config() == {
        "neighbours": 20,
        "ordering": "coordinate",
        "seed": 3,
    }


def test_latent_block_is_n_not_a_rank() -> None:
    """The one place the two candidates differ in kind, asserted.

    A reduced-rank solver's latent block is ``m``; Vecchia has no rank to
    shrink to, so it is ``N`` — the contract's default, and a real cost under
    NUTS that the bake-off has to weigh against its accuracy.
    """
    assert VecchiaResponseGP(neighbours=20).latent_size(Matern32(0.3, 2.0), 10_000) == 10_000
    assert EquispacedFourierGP(basis_size=33).latent_size(Matern32(0.3, 2.0), 10_000) == 33


@pytest.mark.parametrize(
    ("kwargs", "fragment"),
    [
        ({"neighbours": 0}, "at least 1"),
        ({"ordering": "maximin"}, "ordering must be one of"),
        ({"jitter": -1.0}, "jitter must be finite"),
    ],
)
def test_declaration_refusals(kwargs: dict[str, object], fragment: str) -> None:
    with pytest.raises(LikelihoodError, match=fragment):
        VecchiaResponseGP(**kwargs)  # type: ignore[arg-type]


def test_accepts_a_kernel_no_spectral_method_can_represent() -> None:
    """A ``Product`` has no spectral density, and Vecchia does not need one."""
    data = _spectrum()
    product = Product(Matern32(0.3, 2.0), Matern12(0.2, 0.5))
    GaussianProcessNoise(product, VecchiaResponseGP()).check_compatible(GaussianFamily(), data)
    with pytest.raises(LikelihoodError, match="closed-form spectral density"):
        GaussianProcessNoise(product, EquispacedFourierGP()).check_compatible(
            GaussianFamily(), data
        )


# ---------------------------------------------------------------------------
# The structure
# ---------------------------------------------------------------------------


def test_neighbour_sets_are_the_true_nearest_predecessors() -> None:
    """The blocked KD-tree search must agree with brute force, exactly."""
    rng = np.random.default_rng(4)
    points = rng.uniform(0.0, 1.0, (400, 2))
    width = 7
    structure = neighbour_structure(points, width, ordering="random", seed=2)
    ordered = points[structure.order]
    for position in (1, 5, 40, 399):
        expected_count = min(position, width)
        distances = np.linalg.norm(ordered[:position] - ordered[position], axis=1)
        expected = set(np.argsort(distances)[:expected_count].tolist())
        found = set(structure.neighbours[position][: structure.counts[position]].tolist())
        assert found == expected, f"position {position}"


def test_neighbour_structure_is_cached_by_content() -> None:
    rng = np.random.default_rng(5)
    points = rng.uniform(0.0, 1.0, (300, 2))
    first = neighbour_structure(points, 6, ordering="random", seed=1)
    again = neighbour_structure(points.copy(), 6, ordering="random", seed=1)
    assert again is first, "the structure is a constant of the data and is reused"
    other = neighbour_structure(points, 6, ordering="coordinate", seed=1)
    assert other is not first


def test_separation_covariance_matches_the_public_matrix() -> None:
    """The batched evaluation must be the kernel's own arithmetic, not a copy."""
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    separations = np.array([[0.0, 0.3], [1.2, 5.0]])
    found = separation_covariance(kernel, separations, values, 1)
    expected = np.asarray(
        kernel.matrix(np.zeros((1, 1)), separations.reshape(-1, 1), values)
    ).reshape(2, 2)
    assert np.allclose(found, expected)


def test_stacked_residuals_are_the_sum_of_their_columns(
    problem: tuple[np.ndarray, ...],
) -> None:
    points, residual, variance = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    solver = VecchiaResponseGP(neighbours=16, seed=5)
    block = np.column_stack([residual, 0.5 * residual])
    joint = solver.log_marginal_likelihood(kernel, points, block, variance, values)
    apart = solver.log_marginal_likelihood(
        kernel, points, residual, variance, values
    ) + solver.log_marginal_likelihood(kernel, points, 0.5 * residual, variance, values)
    assert joint == pytest.approx(apart, rel=1e-12)


def test_latent_transform_draws_the_declared_covariance(
    problem: tuple[np.ndarray, ...],
) -> None:
    """``U^-T`` whitens the **noise-free** process, to within Monte Carlo error."""
    points, _, _ = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    solver = VecchiaResponseGP(neighbours=32, seed=5)
    total = points.shape[0]
    draws = np.random.default_rng(9).normal(size=(total, 8000))
    field = solver.latent_transform(kernel, points, draws, values)
    assert field.shape == (total, 8000)
    exact = np.asarray(kernel.matrix(points, points, values))
    empirical = np.cov(field)
    # Five standard errors of the empirical covariance of 8000 draws.
    tolerance = 5.0 * exact[0, 0] * np.sqrt(2.0 / 8000)
    assert np.max(np.abs(empirical - exact)) < tolerance
    with pytest.raises(LikelihoodError, match="one per retained sample"):
        solver.latent_transform(kernel, points, np.zeros(7), values)


# ---------------------------------------------------------------------------
# The approximation: a convergence in k, and a dependence on the ordering
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("family", "order", "final"),
    # The order is an **empirical** rate here, not a spectral law: Vecchia's
    # error has no closed-form dependence on ``k``, so each row states the
    # slowest segment its own sweep measures. Matérn-1/2 falls at better than
    # 1 per doubling throughout; Matérn-3/2's first step buys only 0.62 before
    # the rate rises past 3, so its envelope is anchored on that.
    [(Matern12, 1.0, 1.0e-6), (Matern32, 0.5, 1.0e-1)],
    ids=["matern12", "matern32"],
)
def test_marginal_likelihood_converges_in_k(
    problem: tuple[np.ndarray, ...], family: type, order: float, final: float
) -> None:
    points, residual, variance = problem
    kernel = family(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().log_marginal_likelihood(kernel, points, residual, variance, values)
    errors = [
        abs(
            VecchiaResponseGP(neighbours=width, seed=5).log_marginal_likelihood(
                kernel, points, residual, variance, values
            )
            - exact
        )
        for width in WIDTHS
    ]
    for width, error in zip(WIDTHS[1:], errors[1:], strict=True):
        assert error <= _envelope(errors[0], WIDTHS[0], width, order), (
            f"{family.__name__} at k={width}: {error:.3e} outside the envelope anchored on "
            f"{errors[0]:.3e} at k={WIDTHS[0]}"
        )
    assert errors[-1] < final


def test_converges_where_a_spectral_method_cannot(problem: tuple[np.ndarray, ...]) -> None:
    """The bake-off's headline row, in one assertion.

    A Matérn-1/2 is the roughest process ampere supports. Its spectral density
    decays as ``omega**-2``, so a reduced-rank method's truncated tail falls
    only as ``m**-1`` and even a 513-point grid is nats away; Vecchia needs no
    tail to converge, because the screening effect is *strongest* exactly where
    the process is roughest.
    """
    points, residual, variance = problem
    kernel = Matern12(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().log_marginal_likelihood(kernel, points, residual, variance, values)
    sparse = VecchiaResponseGP(neighbours=32, seed=5).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    spectral = EquispacedFourierGP(basis_size=513, boundary_factor=2.0).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    assert abs(sparse - exact) < 1e-6
    assert abs(spectral - exact) > 1.0
    assert abs(sparse - exact) < abs(spectral - exact) / 1e6


def test_screening_is_weak_for_a_smooth_process(problem: tuple[np.ndarray, ...]) -> None:
    """And the price, measured: the two methods' strengths are complementary.

    Stein (2002) on the screening effect: conditioning on near neighbours
    screens the far field well for a rough process and badly for a smooth one,
    so a Matérn-5/2 — where the spectral methods reach machine precision — is
    where Vecchia's error stops falling reliably with ``k``. The row asserts
    the *comparison* rather than a rate, because a non-monotone error is
    exactly what it is recording.
    """
    points, residual, variance = problem
    kernel = Matern52(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().log_marginal_likelihood(kernel, points, residual, variance, values)
    sparse = abs(
        VecchiaResponseGP(neighbours=32, seed=5).log_marginal_likelihood(
            kernel, points, residual, variance, values
        )
        - exact
    )
    spectral = abs(
        EquispacedFourierGP(basis_size=513, boundary_factor=2.0).log_marginal_likelihood(
            kernel, points, residual, variance, values
        )
        - exact
    )
    assert spectral < 1e-4
    assert sparse > spectral


def test_the_ordering_is_part_of_the_approximation(
    problem: tuple[np.ndarray, ...],
) -> None:
    """Two orderings are two different numbers, which is why both are declared."""
    points, residual, variance = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    random = VecchiaResponseGP(neighbours=8, ordering="random", seed=5).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    coordinate = VecchiaResponseGP(neighbours=8, ordering="coordinate").log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    assert random != coordinate
    seeded = VecchiaResponseGP(neighbours=8, ordering="random", seed=6).log_marginal_likelihood(
        kernel, points, residual, variance, values
    )
    assert seeded != random, "the seed is a declaration too"


def test_conditional_moments_and_loo_converge(problem: tuple[np.ndarray, ...]) -> None:
    points, residual, variance = problem
    kernel = Matern32(0.4, 1.5)
    values = kernel.resolve(None)
    exact = DenseGP().condition(kernel, points, residual, variance, values)
    exact_loo = DenseGP().conditional_loo(kernel, points, residual, variance, values)
    coarse = VecchiaResponseGP(neighbours=4, seed=5)
    fine = VecchiaResponseGP(neighbours=64, seed=5)
    coarse_conditional = coarse.condition(kernel, points, residual, variance, values)
    fine_conditional = fine.condition(kernel, points, residual, variance, values)
    assert np.all(fine_conditional.variance >= 0.0)
    assert np.max(np.abs(fine_conditional.mean - exact.mean)) < np.max(
        np.abs(coarse_conditional.mean - exact.mean)
    )
    assert np.max(np.abs(fine_conditional.mean - exact.mean)) < 1e-4
    assert np.max(np.abs(fine_conditional.variance - exact.variance)) < 1e-8
    fine_loo = fine.conditional_loo(kernel, points, residual, variance, values)
    coarse_loo = coarse.conditional_loo(kernel, points, residual, variance, values)
    assert fine_loo.shape == exact_loo.shape
    assert np.max(np.abs(fine_loo - exact_loo)) < np.max(np.abs(coarse_loo - exact_loo))
    assert np.max(np.abs(fine_loo - exact_loo)) < 1e-3


def test_module_docstring_examples_run() -> None:
    results = doctest.testmod(vecchia_module, optionflags=doctest.ELLIPSIS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0
