"""Unit tests for the W1.6 Likelihood & NoiseModel contract.

Organised around ``WORK_ITEMS.md`` W1.6's acceptance criteria. The headline
ones come first:

* ``TestDenseGPAgainstAnalyticCases`` — the DenseGP correctness anchor, checked
  against (a) the closed-form i.i.d. Gaussian log-likelihood in the
  zero-amplitude limit, (b) ``scipy.stats.multivariate_normal.logpdf`` on the
  full covariance, and (c) mask excision, which must leave a masked sample with
  *exactly* zero influence;
* ``TestLatentPathDeclaration`` — the analytic-vs-latent declaration, exercised
  by the Poisson family, including the capability flagging and the error path
  on an engine that cannot deliver it.

The supporting behaviour each of those rests on follows.
"""

from __future__ import annotations

import math

import astropy.units as u
import numpy as np
import pytest
import scipy.linalg
import scipy.special
import scipy.stats as st
from scipy.stats import multivariate_normal

from ampere.core import (
    AxisSpec,
    ClosurePhases,
    CauchyFamily,
    Censoring,
    ComplexGaussianFamily,
    DenseGP,
    FunctionSamples,
    GaussianFamily,
    GaussianProcessNoise,
    GPConditional,
    Identity,
    Image,
    IndependentNoise,
    InducingPointGP,
    Layout,
    Likelihood,
    LikelihoodFamily,
    LimitKind,
    Log,
    Marginalisation,
    Matern32,
    NoiseModel,
    NoiseParams,
    Order,
    Parameter,
    PhotometricPoints,
    PoissonFamily,
    QuasisepGP,
    RiceFamily,
    Spectrum,
    SquaredExponential,
    StudentTFamily,
    TimeSeries,
    VisibilitySet,
    VonMisesFamily,
    WindowedSparseGP,
    family_named,
    latent_parameter,
    list_families,
    register_family,
)
from ampere.core.exceptions import ContractError, LikelihoodError

# ---------------------------------------------------------------------------
# Fixtures: one small, irregular, ordered spectrum and a model for it
# ---------------------------------------------------------------------------

N_SAMPLES = 14


@pytest.fixture
def coordinates() -> np.ndarray:
    rng = np.random.default_rng(20260901)
    return np.sort(rng.uniform(1.0, 12.0, N_SAMPLES))


@pytest.fixture
def observed(coordinates: np.ndarray) -> Spectrum:
    rng = np.random.default_rng(7)
    truth = np.sin(coordinates) + 0.3 * coordinates
    return Spectrum(
        coordinates * u.um,
        (truth + rng.normal(0.0, 0.15, coordinates.size)) * u.Jy,
        uncertainty=np.full(coordinates.size, 0.15) * u.Jy,
    )


@pytest.fixture
def predicted(observed: Spectrum, coordinates: np.ndarray) -> Spectrum:
    # A deliberately imperfect model: the right trend, the wrong wiggles. This
    # is the situation the flexible likelihood exists for.
    return observed.with_values(0.3 * coordinates)


def matern32_matrix(x: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """The reference Matern-3/2 covariance, written out independently."""
    separation = np.abs(x[:, None] - x[None, :])
    scaled = math.sqrt(3.0) * separation / length_scale
    return amplitude**2 * (1.0 + scaled) * np.exp(-scaled)


def iid_gaussian_log_prob(residual: np.ndarray, sigma: np.ndarray) -> float:
    return float(np.sum(-0.5 * ((residual / sigma) ** 2 + math.log(2.0 * math.pi)) - np.log(sigma)))


# ---------------------------------------------------------------------------
# Headline criterion 1: DenseGP against analytic marginal-likelihood cases
# ---------------------------------------------------------------------------


class TestDenseGPAgainstAnalyticCases:
    """``Accept``: DenseGP validated against analytic marginal-likelihood cases."""

    def test_zero_amplitude_reduces_exactly_to_iid_gaussian(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        """A GP of zero amplitude is no GP at all, to the last bit."""
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.0, 2.0)))
        plain = Likelihood(GaussianFamily(), IndependentNoise())
        assert gp.log_prob(predicted, observed) == pytest.approx(
            plain.log_prob(predicted, observed), rel=0.0, abs=1e-12
        )

    def test_zero_amplitude_matches_the_closed_form_written_out(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.0, 2.0)))
        expected = iid_gaussian_log_prob(
            observed.values - predicted.values, np.asarray(observed.uncertainty)
        )
        assert gp.log_prob(predicted, observed) == pytest.approx(expected, abs=1e-12)

    @pytest.mark.parametrize(
        ("amplitude", "length_scale"),
        [(0.05, 0.5), (0.3, 2.0), (1.5, 7.0), (0.3, 0.05)],
    )
    def test_agrees_with_scipy_multivariate_normal(
        self,
        predicted: Spectrum,
        observed: Spectrum,
        coordinates: np.ndarray,
        amplitude: float,
        length_scale: float,
    ) -> None:
        """The definition of the right answer, from an independent implementation."""
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(amplitude, length_scale)))
        covariance = matern32_matrix(coordinates, amplitude, length_scale)
        covariance = covariance + np.diag(np.asarray(observed.uncertainty) ** 2)
        expected = multivariate_normal.logpdf(
            observed.values - predicted.values,
            mean=np.zeros(coordinates.size),
            cov=covariance,
        )
        assert gp.log_prob(predicted, observed) == pytest.approx(float(expected), abs=1e-9)

    def test_agrees_with_scipy_for_the_squared_exponential_too(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """The solver is kernel-agnostic; only the kernel changes."""
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(SquaredExponential(0.4, 1.5)))
        separation = np.abs(coordinates[:, None] - coordinates[None, :])
        covariance = 0.16 * np.exp(-0.5 * (separation / 1.5) ** 2)
        covariance = covariance + np.diag(np.asarray(observed.uncertainty) ** 2)
        expected = multivariate_normal.logpdf(
            observed.values - predicted.values, mean=np.zeros(coordinates.size), cov=covariance
        )
        assert gp.log_prob(predicted, observed) == pytest.approx(float(expected), abs=1e-9)

    def test_a_masked_sample_is_excised_not_downweighted(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """Masking equals deleting the sample: same number, exactly."""
        mask = np.zeros(coordinates.size, dtype=bool)
        mask[4] = True
        mask[9] = True
        masked_observed = Spectrum(
            coordinates * u.um,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
            mask=mask,
        )
        masked_predicted = masked_observed.with_values(predicted.values)

        keep = ~mask
        sub_observed = Spectrum(
            coordinates[keep] * u.um,
            observed.values[keep] * u.Jy,
            uncertainty=np.asarray(observed.uncertainty)[keep] * u.Jy,
        )
        sub_predicted = sub_observed.with_values(predicted.values[keep])

        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0)))
        assert gp.log_prob(masked_predicted, masked_observed) == pytest.approx(
            gp.log_prob(sub_predicted, sub_observed), rel=0.0, abs=1e-12
        )

    def test_a_masked_sample_has_no_influence_whatever_its_value(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """Zero information means zero: not a small contribution, none."""
        mask = np.zeros(coordinates.size, dtype=bool)
        mask[6] = True
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0)))

        def build(value: float, sigma: float) -> tuple[Spectrum, Spectrum]:
            values = np.array(observed.values, copy=True)
            values[6] = value
            uncertainty = np.array(observed.uncertainty, copy=True)
            uncertainty[6] = sigma
            data = Spectrum(
                coordinates * u.um, values * u.Jy, uncertainty=uncertainty * u.Jy, mask=mask
            )
            return data.with_values(predicted.values), data

        baseline = gp.log_prob(*build(0.0, 0.15))
        assert gp.log_prob(*build(1e6, 0.15)) == pytest.approx(baseline, rel=0.0, abs=0.0)
        assert gp.log_prob(*build(0.0, 1e-6)) == pytest.approx(baseline, rel=0.0, abs=0.0)

    def test_masking_and_excision_agree_for_independent_noise_too(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """The two rules coincide term by term when there is no covariance."""
        mask = np.zeros(coordinates.size, dtype=bool)
        mask[2] = True
        data = Spectrum(
            coordinates * u.um,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
            mask=mask,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        expected = iid_gaussian_log_prob(
            (observed.values - predicted.values)[~mask],
            np.asarray(observed.uncertainty)[~mask],
        )
        assert like.log_prob(data.with_values(predicted.values), data) == pytest.approx(
            expected, abs=1e-12
        )

    def test_the_infinite_variance_limit_really_does_diverge(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """The argument the mask convention rests on, kept honest.

        ``likelihoods.md`` §8 rejects ``masked_uncertainty()`` on the grounds
        that ``sigma_i -> inf`` in a GP marginal likelihood does not converge to
        the excised value but to ``excised - 0.5*log(2*pi*sigma_i**2)``, which
        diverges. If that ever stopped being true the convention would need
        revisiting, so it is asserted rather than merely asserted *about*.
        """
        residual = observed.values - predicted.values
        covariance = matern32_matrix(coordinates, 0.4, 2.0)
        sigma = np.asarray(observed.uncertainty)
        keep = np.ones(coordinates.size, dtype=bool)
        keep[5] = False

        excised = multivariate_normal.logpdf(
            residual[keep],
            mean=np.zeros(int(keep.sum())),
            cov=covariance[np.ix_(keep, keep)] + np.diag(sigma[keep] ** 2),
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0)))
        masked = Spectrum(
            coordinates * u.um,
            observed.values * u.Jy,
            uncertainty=sigma * u.Jy,
            mask=~keep,
        )
        assert gp.log_prob(masked.with_values(predicted.values), masked) == pytest.approx(
            float(excised), abs=1e-9
        )

        # ... whereas inflating the uncertainty converges to a *different*,
        # divergent quantity, one that gets worse the larger sigma is made.
        for inflated in (1e1, 1e2, 1e3):
            widened = sigma.copy()
            widened[5] = inflated
            value = multivariate_normal.logpdf(
                residual,
                mean=np.zeros(coordinates.size),
                cov=covariance + np.diag(widened**2),
            )
            offset = 0.5 * math.log(2.0 * math.pi * inflated**2)
            assert float(value) == pytest.approx(float(excised) - offset, abs=1e-3)
            assert float(value) < float(excised) - 2.0

    def test_the_diagonal_case_diverges_by_exactly_the_same_amount(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """The half of the argument it would be easy to get wrong.

        A zero weight and an infinite uncertainty are equivalent statements
        about a *chi-square term*, which is what ``results_schema.md`` §7
        claims for them — but not about a normalised log-density, because each
        term carries a ``-log sigma`` that diverges too. The uncorrelated case
        is therefore no more forgiving than the GP case, and the spec must not
        claim it is.
        """
        residual = observed.values - predicted.values
        sigma = np.asarray(observed.uncertainty)
        keep = np.ones(coordinates.size, dtype=bool)
        keep[5] = False
        excised = iid_gaussian_log_prob(residual[keep], sigma[keep])

        masked = Spectrum(
            coordinates * u.um, observed.values * u.Jy, uncertainty=sigma * u.Jy, mask=~keep
        )
        plain = Likelihood(GaussianFamily(), IndependentNoise())
        assert plain.log_prob(masked.with_values(predicted.values), masked) == pytest.approx(
            excised, abs=1e-12
        )

        for inflated in (1e1, 1e2, 1e3, 1e6):
            widened = sigma.copy()
            widened[5] = inflated
            value = iid_gaussian_log_prob(residual, widened)
            offset = 0.5 * math.log(2.0 * math.pi * inflated**2)
            # The remaining discrepancy is the vanishing chi-square, O(1/sigma^2).
            tolerance = max(5.0 / inflated**2, 1e-9)
            assert value == pytest.approx(excised - offset, abs=tolerance)
            # The chi-square term does vanish, which is the true half of §7's
            # claim; it is the normalisation that does not.
            assert (residual[5] / inflated) ** 2 <= (residual[5] / 10.0) ** 2

    def test_a_fully_masked_pair_contributes_nothing(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        data = Spectrum(
            coordinates * u.um,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
            mask=np.ones(coordinates.size, dtype=bool),
        )
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0)))
        assert like.log_prob(data.with_values(predicted.values), data) == 0.0

    def test_the_predicted_containers_mask_counts_too(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """The effective mask is the union: a sample needs both sides valid."""
        mask = np.zeros(coordinates.size, dtype=bool)
        mask[3] = True
        model = Spectrum(
            coordinates * u.um,
            predicted.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
            mask=mask,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        both_ways = Spectrum(
            coordinates * u.um,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
            mask=mask,
        )
        assert like.log_prob(model, observed) == pytest.approx(
            like.log_prob(both_ways.with_values(predicted.values), both_ways), abs=1e-12
        )

    def test_everything_is_float64(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        """``DEVELOPMENT_PLAN.md`` §7: float64 for all GP linear algebra."""
        kernel = Matern32(0.4, 2.0)
        matrix = kernel.matrix(
            coordinates[:, None].astype(np.float32),
            coordinates[:, None].astype(np.float32),
            {"amplitude": 0.4, "length_scale": 2.0},
        )
        assert matrix.dtype == np.float64

    def test_the_kernel_matches_the_closed_form(self, coordinates: np.ndarray) -> None:
        kernel = Matern32(0.7, 3.0)
        got = kernel.matrix(
            coordinates[:, None], coordinates[:, None], {"amplitude": 0.7, "length_scale": 3.0}
        )
        assert np.allclose(got, matern32_matrix(coordinates, 0.7, 3.0), atol=1e-14)
        assert np.allclose(np.diag(got), 0.49)
        assert np.allclose(got, got.T)
        assert np.linalg.eigvalsh(got).min() > -1e-10

    def test_a_singular_covariance_fails_loudly_with_a_fix(self) -> None:
        """Not a silent nan: a message naming the cause and the remedy."""
        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um,
            [1.0, 1.0, 1.0] * u.Jy,
            uncertainty=[0.0, 0.0, 0.0] * u.Jy,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(1.0, 1.0)))
        with pytest.raises(LikelihoodError, match="zero or negative uncertainties"):
            gp.log_prob(data.with_values([0.0, 0.0, 0.0]), data)


class TestConditionedGP:
    """The localisation output ``DEVELOPMENT_PLAN.md`` §4.8's family C consumes."""

    def test_the_conditioned_mean_recovers_the_residual_in_the_noiseless_limit(
        self, coordinates: np.ndarray
    ) -> None:
        residual = np.sin(coordinates)
        data = Spectrum(
            coordinates * u.um,
            residual * u.Jy,
            uncertainty=np.full(coordinates.size, 1e-4) * u.Jy,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(2.0, 1.0)))
        conditional = gp.conditional(data.with_values(np.zeros(coordinates.size)), data)
        assert isinstance(conditional, GPConditional)
        assert np.allclose(conditional.mean, residual, atol=1e-3)
        assert np.all(conditional.variance >= -1e-12)
        assert np.all(conditional.standard_deviation < 0.05)

    def test_the_conditioned_mean_is_signed(self, coordinates: np.ndarray) -> None:
        """A localisation diagnostic must show direction, not just magnitude."""
        residual = np.where(coordinates < coordinates.mean(), -1.0, 1.0)
        data = Spectrum(
            coordinates * u.um,
            residual * u.Jy,
            uncertainty=np.full(coordinates.size, 1e-3) * u.Jy,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(3.0, 1.0)))
        conditional = gp.conditional(data.with_values(np.zeros(coordinates.size)), data)
        assert conditional.mean.min() < 0.0 < conditional.mean.max()

    def test_it_evaluates_at_masked_coordinates_too(self, coordinates: np.ndarray) -> None:
        """Masked points stay available as evaluation locations, by design."""
        mask = np.zeros(coordinates.size, dtype=bool)
        mask[5] = True
        data = Spectrum(
            coordinates * u.um,
            np.sin(coordinates) * u.Jy,
            uncertainty=np.full(coordinates.size, 1e-3) * u.Jy,
            mask=mask,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(2.0, 1.0)))
        conditional = gp.conditional(data.with_values(np.zeros(coordinates.size)), data)
        assert conditional.mean.size == coordinates.size

    def test_a_bare_1d_grid_means_n_points_not_one_n_dimensional_point(
        self, coordinates: np.ndarray
    ) -> None:
        """``np.atleast_2d`` reads a 1-D grid the wrong way round, silently.

        ``atleast_2d`` on shape ``(m,)`` yields ``(1, m)`` — one m-dimensional
        point — which then *broadcasts* through the kernel rather than erroring,
        turning a 101-point localisation grid into a one-point answer. Passing a
        bare list of wavelengths is the natural call, so it must mean what a
        reader thinks it means.
        """
        data = Spectrum(
            coordinates * u.um,
            np.sin(coordinates) * u.Jy,
            uncertainty=np.full(coordinates.size, 0.05) * u.Jy,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(1.0, 1.5)))
        model = data.with_values(np.zeros(coordinates.size))
        flat = np.linspace(coordinates.min(), coordinates.max(), 25)

        one_dimensional = gp.conditional(model, data, at=flat)
        assert one_dimensional.mean.shape == (25,)
        column = gp.conditional(model, data, at=flat[:, None])
        assert np.array_equal(one_dimensional.mean, column.mean)
        assert np.array_equal(one_dimensional.variance, column.variance)

        # An unusable shape is refused rather than broadcast.
        with pytest.raises(LikelihoodError, match="coordinate\\(s\\) per point"):
            gp.conditional(model, data, at=np.zeros((5, 3)))
        with pytest.raises(LikelihoodError, match="must be a 1-D array"):
            gp.conditional(model, data, at=np.zeros((2, 2, 2)))

    def test_the_solver_reads_a_1d_grid_the_same_way(self, coordinates: np.ndarray) -> None:
        """The same guarantee one level down, where a backend will call it."""
        kernel = Matern32(1.0, 1.5)
        values = {"amplitude": 1.0, "length_scale": 1.5}
        residual = np.sin(coordinates)
        variance = np.full(coordinates.size, 0.05**2)
        flat = np.linspace(coordinates.min(), coordinates.max(), 9)
        got = DenseGP().condition(kernel, coordinates, residual, variance, values, at=flat)
        assert got.mean.shape == (9,)
        # The kernel itself reads a bare 1-D array as a column too.
        assert kernel.matrix(flat, flat, values).shape == (9, 9)
        assert kernel.diagonal(flat, values).shape == (9,)

    def test_it_evaluates_on_a_finer_grid_when_asked(self, coordinates: np.ndarray) -> None:
        """W1.12's family C wants a visualisation grid, not only the data's own."""
        residual = np.sin(coordinates)
        data = Spectrum(
            coordinates * u.um,
            residual * u.Jy,
            uncertainty=np.full(coordinates.size, 0.05) * u.Jy,
        )
        gp = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(1.0, 1.5)))
        grid = np.linspace(coordinates.min(), coordinates.max(), 101)[:, None]
        conditional = gp.conditional(data.with_values(np.zeros(coordinates.size)), data, at=grid)
        assert conditional.mean.shape == (101,)
        assert conditional.variance.shape == (101,)

        # Against the textbook formulae, computed with an explicit inverse.
        covariance = matern32_matrix(coordinates, 1.0, 1.5)
        total = covariance + np.diag(np.full(coordinates.size, 0.05**2))
        inverse = np.linalg.inv(total)
        separation = np.abs(grid - coordinates[None, :])
        scaled = math.sqrt(3.0) * separation / 1.5
        cross = (1.0 + scaled) * np.exp(-scaled)
        assert np.allclose(conditional.mean, cross @ inverse @ residual, atol=1e-9)
        assert np.allclose(
            conditional.variance,
            1.0 - np.einsum("ij,ji->i", cross, inverse @ cross.T),
            atol=1e-9,
        )
        # Away from any datum the posterior must revert towards the prior; at a
        # datum it must be tighter than the prior. Both, or the band is wrong.
        assert conditional.variance.max() > conditional.variance.min()
        assert conditional.variance.max() <= 1.0 + 1e-9

    def test_it_refuses_an_uncorrelated_noise_model(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="induces no correlations"):
            like.conditional(predicted, observed)


# ---------------------------------------------------------------------------
# Headline criterion 2: the latent-path declaration, via the Poisson family
# ---------------------------------------------------------------------------


class TestLatentPathDeclaration:
    """``Accept``: latent-path declaration exercised by a stub Poisson family."""

    @pytest.fixture
    def counts(self) -> Spectrum:
        return Spectrum([1.0, 2.0, 3.0, 4.0] * u.um, [4.0, 7.0, 2.0, 9.0])

    def test_poisson_with_independent_noise_is_analytic(self, counts: Spectrum) -> None:
        like = Likelihood(PoissonFamily(), IndependentNoise())
        assert like.marginalisation is Marginalisation.ANALYTIC
        rate = counts.with_values([3.5, 6.0, 2.5, 8.0])
        assert like.log_prob(rate, counts) == pytest.approx(
            float(np.sum(st.poisson.logpmf(counts.values, rate.values))), abs=1e-12
        )

    def test_poisson_with_a_gp_is_latent(self, counts: Spectrum) -> None:
        """The structural consequence ``DEVELOPMENT_PLAN.md`` §4.4 spells out."""
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert like.marginalisation is Marginalisation.LATENT

    def test_gaussian_with_a_gp_stays_analytic(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        """Only the Gaussian family folds a GP covariance in closed form."""
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert like.marginalisation is Marginalisation.ANALYTIC

    @pytest.mark.parametrize("family", [StudentTFamily(), CauchyFamily(), PoissonFamily()])
    def test_every_non_gaussian_family_goes_latent_under_a_gp(
        self, family: LikelihoodFamily
    ) -> None:
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        assert family.marginalisation_with(noise) is Marginalisation.LATENT

    def test_the_circular_complex_gp_is_declared_analytic_and_implemented(self) -> None:
        """Ruled 2026-09-03 (§17 Q6), implemented at **W4.2**.

        ``complex_gaussian`` + ``GaussianProcessNoise`` declares ``ANALYTIC``
        with the circular (equal-component, zero-pseudo-covariance) complex GP
        as the fixed meaning. Until W4.2 the pair was *refused* at composition
        with Phase 4 named, because the closed form was declared rather than
        written; this row is the other side of that, and it is kept here rather
        than deleted because the staging flag is still the mechanism and the
        declaration is still what decides the marginalisation.
        """
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        family = ComplexGaussianFamily()
        assert family.marginalisation_with(noise) is Marginalisation.ANALYTIC
        assert ComplexGaussianFamily.GP_ANALYTIC_IMPLEMENTED
        assert Likelihood(family, noise).marginalisation is Marginalisation.ANALYTIC

    def test_a_family_whose_closed_form_is_staged_is_still_refused(self) -> None:
        """The discipline outlives its first instance (*W4.2*).

        ``GP_ANALYTIC_IMPLEMENTED`` exists so that a family may declare its GP
        marginalisation before writing it and have composition refuse rather
        than fall back on a model the declaration does not describe. Every
        shipped family now either implements the closed form or is latent, so
        the rule is held against a family declared here.
        """

        class StagedFamily(ComplexGaussianFamily):
            NAME = "staged_complex"
            GP_ANALYTIC_IMPLEMENTED = False

        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        family = StagedFamily()
        with pytest.raises(LikelihoodError, match="GP_ANALYTIC_IMPLEMENTED is False"):
            Likelihood(family, noise)

    @pytest.mark.parametrize("family", [StudentTFamily(), CauchyFamily()])
    def test_a_latent_family_that_ignores_the_latent_values_cannot_be_composed(
        self, family: LikelihoodFamily
    ) -> None:
        """The defect this flag exists to stop, pinned.

        A family whose ``log_prob`` never reads ``noise.latent`` would return
        the *uncorrelated* likelihood under a GP — so every GP hyperparameter
        and every latent value an engine sampled would leave the
        log-probability untouched, and the fit would run, converge and be
        wrong. Declaring LATENT is not enough; the family must implement it.
        """
        assert not family.CONSUMES_LATENT_GP
        with pytest.raises(LikelihoodError, match="does not implement the latent-conditional"):
            Likelihood(family, GaussianProcessNoise(Matern32(0.3, 1.0)))

    def test_the_family_that_does_implement_it_composes(self) -> None:
        assert PoissonFamily.CONSUMES_LATENT_GP
        Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))

    def test_a_third_party_family_must_opt_in_to_the_latent_path(self) -> None:
        """The safe default: a new family does not silently inherit the defect."""
        assert LikelihoodFamily.CONSUMES_LATENT_GP is False

        @register_family
        class TestLatentAware(LikelihoodFamily):
            NAME = "test_latent_aware"
            CONSUMES_LATENT_GP = True
            REQUIRES_UNCERTAINTY = False

            def log_prob(
                self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
            ) -> float:
                if noise.correlated and noise.latent is None:
                    raise LikelihoodError("needs latent values")
                offset = 0.0 if noise.latent is None else float(np.sum(noise.latent))
                return float(np.sum(observed - predicted)) + offset

        try:
            data = Spectrum([1.0, 2.0] * u.um, [1.0, 2.0] * u.Jy)
            model = data.with_values([0.0, 0.0])
            like = Likelihood(TestLatentAware(), GaussianProcessNoise(Matern32(0.3, 1.0)))
            assert like.marginalisation is Marginalisation.LATENT
            # ``latent=`` is the whitened ``z``; the noise model applies
            # ``f = L z`` before the family sees it (W2.14), so the offset the
            # family adds is the sum of the *correlated* values.
            whitened = np.array([1.0, 1.0])
            realised = DenseGP().latent_transform(
                Matern32(0.3, 1.0),
                np.array([[1.0], [2.0]]),
                whitened,
                {"amplitude": 0.3, "length_scale": 1.0},
            )
            assert like.log_prob(model, data, latent=whitened) == pytest.approx(
                3.0 + float(np.sum(realised))
            )
        finally:
            from ampere.core import likelihood as module

            module._FAMILIES.pop("test_latent_aware", None)

    def test_a_gradient_free_engine_is_refused_loudly(self, counts: Spectrum) -> None:
        """Capability flagging: the plan's "modern-backend capability", enforced."""
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        with pytest.raises(LikelihoodError, match="emcee cannot run this likelihood"):
            like.check_engine(differentiable=False, engine="emcee")
        # A differentiable engine is fine, and so is any engine on an analytic
        # likelihood: the check must not be a blanket refusal.
        like.check_engine(differentiable=True, engine="NUTS")
        Likelihood(PoissonFamily(), IndependentNoise()).check_engine(
            differentiable=False, engine="emcee"
        )

    def test_the_latent_declaration_is_whitened_and_iid(self, counts: Spectrum) -> None:
        """Not a HierarchicalPrior: an i.i.d. standard normal, correlated by L."""
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        declaration = like.latent_declaration(counts.n_samples)
        parameter = declaration.parameter
        assert parameter.shape == (counts.n_samples,)
        assert parameter.is_free and not parameter.is_hierarchical
        assert parameter.references == ()
        assert isinstance(parameter.unconstraining_bijection(), Identity)
        assert declaration.as_parameter_set().free_size == counts.n_samples

    def test_an_analytic_likelihood_declares_no_latent_values(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        with pytest.raises(LikelihoodError, match="marginalises analytically"):
            like.latent_declaration(observed.n_samples)

    def test_the_whitening_transform_reproduces_the_kernel(self, coordinates: np.ndarray) -> None:
        """``f = L z`` with ``L L^T == K``: the correlation the prior does not carry."""
        kernel = Matern32(0.6, 2.0)
        solver = DenseGP()
        values = {"amplitude": 0.6, "length_scale": 2.0}
        points = coordinates[:, None]
        basis = np.eye(coordinates.size)
        columns = [
            solver.latent_transform(kernel, points, basis[:, i], values)
            for i in range(coordinates.size)
        ]
        lower = np.column_stack(columns)
        assert np.allclose(lower @ lower.T, matern32_matrix(coordinates, 0.6, 2.0), atol=1e-8)
        assert np.allclose(np.triu(lower, k=1), 0.0)

    def test_the_latent_poisson_likelihood_is_conditional_on_f(self, counts: Spectrum) -> None:
        """``latent=`` is ``z``; the family scores at ``f = L(theta) z`` (W2.14)."""
        kernel = Matern32(0.3, 1.0)
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(kernel))
        rate = counts.with_values([3.5, 6.0, 2.5, 8.0])
        with pytest.raises(LikelihoodError, match="latent-variable model"):
            like.log_prob(rate, counts)
        whitened = np.array([0.1, -0.2, 0.05, 0.0])
        points = np.asarray(counts.axes[0].values, dtype=float).reshape(-1, 1)
        realised = DenseGP().latent_transform(
            kernel, points, whitened, {"amplitude": 0.3, "length_scale": 1.0}
        )
        assert like.log_prob(rate, counts, latent=whitened) == pytest.approx(
            float(np.sum(st.poisson.logpmf(counts.values, rate.values * np.exp(realised)))),
            abs=1e-12,
        )

    def test_latent_values_on_an_analytic_likelihood_are_refused(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="marginalises analytically"):
            like.log_prob(predicted, observed, latent=np.zeros(observed.n_samples))

    def test_a_wrong_sized_latent_vector_is_caught(self, counts: Spectrum) -> None:
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        rate = counts.with_values([3.5, 6.0, 2.5, 8.0])
        with pytest.raises(LikelihoodError, match="One latent value per retained sample"):
            like.log_prob(rate, counts, latent=np.zeros(3))

    def test_non_integer_counts_are_refused(self) -> None:
        # At composition since 2026-09-02 (X-3, via the check_observed hook):
        # an O(N) property of the data is checked once, not per draw.
        counts = Spectrum([1.0, 2.0] * u.um, [4.5, 7.0])
        like = Likelihood(PoissonFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="non-negative integer counts"):
            like.check_alignment(counts.with_values([4.0, 7.0]), counts)

    def test_latent_parameter_rejects_a_nonsense_size(self) -> None:
        with pytest.raises(LikelihoodError, match="positive number of latent values"):
            latent_parameter("latent", 0)


# ---------------------------------------------------------------------------
# W2.14: the latent path sees its kernel
# ---------------------------------------------------------------------------


def cholesky_factor(x: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """``L`` with ``L Lt = K + jitter``, written out from scipy, not from ampere.

    The stabiliser is ``GPSolver.latent_transform``'s documented default — a
    relative ``1e-10`` on the diagonal, scaled by the kernel's own mean
    variance — so this really is the same matrix and the comparison is exact
    rather than approximate.
    """
    covariance = matern32_matrix(x, amplitude, length_scale)
    stabilised = covariance + np.eye(x.size) * (1e-10 * float(np.mean(np.diag(covariance))))
    return scipy.linalg.cholesky(stabilised, lower=True)


class TestTheLatentPathSeesItsKernel:
    """W2.14: ``noise.latent`` is ``f = L(theta) z``, not the whitened ``z``.

    The defect this class exists for was found independently by both Phase-2
    backend tracks and reproduced on master: ``latent_parameter`` declares
    ``z`` whitened and says the correlation "enters through ``f = L(theta) z``,
    a deterministic transform owned by the ``GPSolver``", and *nothing on the
    scoring path applied it*. ``Dataset.log_likelihood_of`` handed ``z``
    straight to the family, which read it as ``f``, so the kernel
    hyperparameters entered the likelihood nowhere at all and a latent fit
    sampled the amplitude and the length scale against a flat likelihood while
    reporting nothing wrong.

    Every row here fails on the pre-W2.14 code, which is the point: the rows
    that assert *movement* fail because the values were bit-identical, and the
    rows that assert a closed form fail because the value scored was the one
    at ``f = z``.
    """

    @pytest.fixture
    def counts(self) -> Spectrum:
        return Spectrum([1.0, 2.0, 3.0, 4.0] * u.um, [4.0, 7.0, 2.0, 9.0])

    @pytest.fixture
    def rate(self, counts: Spectrum) -> Spectrum:
        return counts.with_values([3.5, 6.0, 2.5, 8.0])

    @pytest.fixture
    def whitened(self) -> np.ndarray:
        return np.array([0.4, -1.1, 0.25, 0.9])

    @staticmethod
    def grid(counts: Spectrum) -> np.ndarray:
        return np.asarray(counts.axes[0].values, dtype=float)

    # -- the transform itself, where the fix lives --------------------------

    def test_noise_params_returns_f_not_z(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """The one line the fix is: ``noise.latent`` is the correlated draw.

        Checked against ``scipy.linalg.cholesky`` of the same stabilised
        kernel matrix, so nothing in the comparison came from the solver under
        test.
        """
        noise = GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP())
        retain = np.ones(counts.n_samples, dtype=bool)
        points = self.grid(counts).reshape(-1, 1)
        params = noise.noise_params(
            counts,
            retain,
            {"amplitude": 0.7, "length_scale": 1.5},
            predicted=np.asarray(rate.values, dtype=float),
            coordinates=points,
            latent=whitened,
        )
        assert params.latent is not None
        expected = cholesky_factor(self.grid(counts), 0.7, 1.5) @ whitened
        assert params.latent == pytest.approx(expected, abs=1e-12)
        # And it is genuinely a different vector from the one handed in, which
        # is what the pre-W2.14 code returned.
        assert not np.allclose(params.latent, whitened)

    def test_the_gaussian_closed_form_at_f_equals_l_z(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """The whitened-``z`` Gaussian statement, written out.

        With ``f = L(theta) z`` the conditional density of ``y`` is
        ``N(y; mu + f, sigma)``. This row scores that closed form off the
        ``f`` the noise model produced, against one built entirely from
        ``scipy.linalg.cholesky``, at three amplitudes — so the two track each
        other rather than agreeing by accident at one point.
        """
        sigma = np.full(counts.n_samples, 0.4)
        observed = np.asarray(counts.values, dtype=float)
        mean = np.asarray(rate.values, dtype=float)
        points = self.grid(counts).reshape(-1, 1)
        retain = np.ones(counts.n_samples, dtype=bool)
        scored = []
        for amplitude in (0.3, 0.9, 2.5):
            noise = GaussianProcessNoise(Matern32(amplitude, 1.5), DenseGP())
            params = noise.noise_params(
                counts,
                retain,
                {"amplitude": amplitude, "length_scale": 1.5},
                predicted=mean,
                coordinates=points,
                latent=whitened,
            )
            assert params.latent is not None
            reference = cholesky_factor(self.grid(counts), amplitude, 1.5) @ whitened
            got = float(np.sum(st.norm.logpdf(observed, loc=mean + params.latent, scale=sigma)))
            expected = float(np.sum(st.norm.logpdf(observed, loc=mean + reference, scale=sigma)))
            assert got == pytest.approx(expected, abs=1e-12)
            scored.append(got)
        assert len(set(scored)) == 3

    # -- the consequence: the likelihood moves ------------------------------

    def test_the_log_likelihood_moves_with_the_amplitude(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """Bit-identical for 0.5, 5 and 50 before W2.14; three numbers now."""
        like = Likelihood(
            PoissonFamily(), GaussianProcessNoise(Matern32(st.halfnorm(0.0, 5.0), 1.5), DenseGP())
        )
        values = [
            like.log_prob(rate, counts, {"amplitude": amplitude}, latent=whitened)
            for amplitude in (0.5, 5.0, 50.0)
        ]
        assert len(set(values)) == 3
        assert all(math.isfinite(value) for value in values)

    def test_the_log_likelihood_moves_with_the_length_scale(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """The same statement for the other free hyperparameter."""
        like = Likelihood(
            PoissonFamily(),
            GaussianProcessNoise(Matern32(0.7, st.loguniform(0.01, 1e3)), DenseGP()),
        )
        values = [
            like.log_prob(rate, counts, {"length_scale": length_scale}, latent=whitened)
            for length_scale in (0.1, 1.0, 100.0)
        ]
        assert len(set(values)) == 3
        assert all(math.isfinite(value) for value in values)

    def test_it_equals_the_poisson_closed_form_at_f(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """``sum log Poisson(y; rate e^(L z))``, from scipy end to end."""
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP()))
        realised = cholesky_factor(self.grid(counts), 0.7, 1.5) @ whitened
        expected = float(np.sum(st.poisson.logpmf(counts.values, rate.values * np.exp(realised))))
        assert like.log_prob(rate, counts, latent=whitened) == pytest.approx(expected, abs=1e-12)

    def test_the_zero_amplitude_limit_is_the_uncorrelated_likelihood(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """``K = 0`` means ``f = 0``: the sanity check at the edge of the family.

        Not *exactly* zero, and the reason is worth recording rather than
        hiding behind a loose tolerance: ``latent_transform`` scales its
        relative jitter by the kernel's own mean variance and falls back to
        ``1.0`` when that is zero, so a zero-amplitude kernel factorises
        ``1e-10 * I`` and ``f`` is ``1e-5 z`` rather than ``0``. The
        uncorrelated limit is therefore approached to about ``1e-5`` in ``f``,
        which is what this asserts.
        """
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.0, 1.5), DenseGP()))
        assert like.log_prob(rate, counts, latent=whitened) == pytest.approx(
            float(np.sum(st.poisson.logpmf(counts.values, rate.values))), abs=1e-4
        )

    def test_the_two_solvers_agree_on_the_realised_likelihood(
        self, counts: Spectrum, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """§4.6's DenseGP-QuasisepGP row, now reachable on the latent path too."""
        dense = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP()))
        quasisep = Likelihood(
            PoissonFamily(), GaussianProcessNoise(Matern32(0.7, 1.5), QuasisepGP())
        )
        assert quasisep.log_prob(rate, counts, latent=whitened) == pytest.approx(
            dense.log_prob(rate, counts, latent=whitened), abs=1e-9
        )

    # -- the failures the transform makes possible ---------------------------

    def test_a_mis_sized_latent_block_is_refused_by_name(self, counts: Spectrum) -> None:
        """Named before a solver produces a matrix-shape error instead."""
        noise = GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP())
        with pytest.raises(LikelihoodError, match="One latent value per retained sample"):
            noise.noise_params(
                counts,
                np.ones(counts.n_samples, dtype=bool),
                {"amplitude": 0.7, "length_scale": 1.5},
                coordinates=self.grid(counts).reshape(-1, 1),
                latent=np.zeros(counts.n_samples + 1),
            )

    def test_latent_values_without_coordinates_are_refused_by_name(self, counts: Spectrum) -> None:
        """A direct caller that forgets them is told what is missing."""
        noise = GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP())
        with pytest.raises(LikelihoodError, match="without the coordinates"):
            noise.noise_params(
                counts,
                np.ones(counts.n_samples, dtype=bool),
                {"amplitude": 0.7, "length_scale": 1.5},
                latent=np.zeros(counts.n_samples),
            )

    def test_a_masked_sample_shrinks_the_transform_with_the_block(
        self, rate: Spectrum, whitened: np.ndarray
    ) -> None:
        """The transform runs on the *retained* coordinates, not on all of them.

        Masking is excision (``likelihoods.md`` §9), so a masked sample leaves
        the covariance entirely — its row and column never enter ``L``. The
        oracle is therefore the Cholesky of the sub-matrix, which is a
        different factor from any three rows of the full one.
        """
        masked = Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [4.0, 7.0, 2.0, 9.0],
            mask=np.array([False, True, False, False]),
        )
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.7, 1.5), DenseGP()))
        kept = [0, 2, 3]
        realised = cholesky_factor(np.array([1.0, 3.0, 4.0]), 0.7, 1.5) @ whitened[kept]
        expected = float(
            np.sum(
                st.poisson.logpmf(
                    np.asarray(masked.values)[kept],
                    np.asarray(rate.values)[kept] * np.exp(realised),
                )
            )
        )
        got = like.log_prob(rate, masked, latent=whitened[kept])
        assert got == pytest.approx(expected, abs=1e-12)


# ---------------------------------------------------------------------------
# The family registry
# ---------------------------------------------------------------------------


class TestFamilyRegistry:
    def test_the_target_family_set_is_on_record(self) -> None:
        """§4.4's list, whether or not each one is implemented yet."""
        assert set(list_families()) >= {
            "gaussian",
            "student_t",
            "cauchy",
            "complex_gaussian",
            "poisson",
            "rice",
            "von_mises",
        }

    def test_lookup_by_name(self) -> None:
        assert family_named("gaussian") is GaussianFamily
        with pytest.raises(LikelihoodError, match="no likelihood family named 'nonesuch'"):
            family_named("nonesuch")

    def test_a_user_can_register_a_family_without_touching_ampere(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        """§4.4: "so users can add their own without touching ampere"."""

        @register_family
        class LaplaceFamily(LikelihoodFamily):
            NAME = "test_laplace"

            def log_prob(
                self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
            ) -> float:
                sigma = noise.sigma
                assert sigma is not None
                return float(
                    np.sum(st.laplace.logpdf((observed - predicted) / sigma) - np.log(sigma))
                )

        try:
            assert family_named("test_laplace") is LaplaceFamily
            like = Likelihood(LaplaceFamily(), IndependentNoise())
            expected = float(
                np.sum(
                    st.laplace.logpdf(
                        (observed.values - predicted.values) / np.asarray(observed.uncertainty)
                    )
                    - np.log(np.asarray(observed.uncertainty))
                )
            )
            assert like.log_prob(predicted, observed) == pytest.approx(expected, abs=1e-12)
        finally:
            from ampere.core import likelihood as module

            module._FAMILIES.pop("test_laplace", None)

    def test_a_duplicate_name_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="already registered"):

            @register_family
            class Clash(LikelihoodFamily):
                NAME = "gaussian"

                def log_prob(
                    self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
                ) -> float:
                    return 0.0

    def test_an_unnamed_family_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="must declare a non-empty NAME"):

            @register_family
            class Anonymous(LikelihoodFamily):
                def log_prob(
                    self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
                ) -> float:
                    return 0.0

    @pytest.mark.parametrize("family", [RiceFamily()])
    def test_declared_but_unimplemented_families_refuse_composition(
        self, family: LikelihoodFamily
    ) -> None:
        """Declared in the registry, refused at composition — never silently wrong.

        ``von_mises`` left this list at W4.1, when the interferometric
        modality brought the data type that needed it. ``rice`` stays: the
        sampling-form principle (``likelihoods.md`` §3) says a family's
        likelihood and its draw arrive together with the data type that needs
        them, and Rice's data type is polarimetry.
        """
        assert family.NAME in list_families()
        assert not family.IMPLEMENTED
        with pytest.raises(LikelihoodError, match="declared but not implemented"):
            Likelihood(family, IndependentNoise())

    def test_a_non_family_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="needs a LikelihoodFamily"):
            Likelihood("gaussian", IndependentNoise())  # type: ignore[arg-type]

    def test_a_non_noise_model_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="needs a NoiseModel"):
            Likelihood(GaussianFamily(), "white")  # type: ignore[arg-type]


class TestFamilyMathematics:
    def test_student_t_matches_scipy(self, predicted: Spectrum, observed: Spectrum) -> None:
        like = Likelihood(StudentTFamily(nu=3.0), IndependentNoise())
        sigma = np.asarray(observed.uncertainty)
        standardised = (observed.values - predicted.values) / sigma
        expected = float(np.sum(st.t.logpdf(standardised, 3.0) - np.log(sigma)))
        assert like.log_prob(predicted, observed) == pytest.approx(expected, abs=1e-12)

    def test_student_t_degrees_of_freedom_are_an_ordinary_parameter(self) -> None:
        family = StudentTFamily(nu=st.loguniform(2.0, 50.0))
        assert family.parameters.free_names == ("nu",)
        assert isinstance(family.parameters["nu"].unconstraining_bijection(), Log)

    def test_cauchy_matches_scipy(self, predicted: Spectrum, observed: Spectrum) -> None:
        like = Likelihood(CauchyFamily(), IndependentNoise())
        sigma = np.asarray(observed.uncertainty)
        standardised = (observed.values - predicted.values) / sigma
        expected = float(np.sum(st.cauchy.logpdf(standardised) - np.log(sigma)))
        assert like.log_prob(predicted, observed) == pytest.approx(expected, abs=1e-12)

    def test_the_complex_gaussian_is_the_circular_model(self) -> None:
        """``results_schema.md`` §16: the container's real sigma means exactly this."""
        vis = VisibilitySet(
            [120.0, -35.0, 88.0],
            [45.0, 190.0, -66.0],
            [2.2, 2.2, 2.2] * u.um,
            [1.0 + 0.2j, 0.6 - 0.3j, 0.4 + 0.0j],
            uncertainty=[0.02, 0.03, 0.05],
        )
        model = vis.with_values([1.02 + 0.18j, 0.58 - 0.31j, 0.41 + 0.01j])
        like = Likelihood(ComplexGaussianFamily(), IndependentNoise())
        like.check_alignment(model, vis)
        sigma = np.asarray(vis.uncertainty)
        residual = vis.values - model.values
        expected = float(
            np.sum(
                st.norm.logpdf(residual.real, scale=sigma)
                + st.norm.logpdf(residual.imag, scale=sigma)
            )
        )
        assert like.log_prob(model, vis) == pytest.approx(expected, abs=1e-12)

    def test_a_real_family_refuses_complex_data(self) -> None:
        vis = VisibilitySet(
            [1.0, 2.0], [3.0, 4.0], [2.2, 2.2] * u.um, [1.0 + 0j, 0.5 + 0j], uncertainty=[0.1, 0.1]
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="holds real values"):
            like.check_alignment(vis.with_values([1.0 + 0j, 0.5 + 0j]), vis)

    def test_scale_and_jitter_are_ordinary_parameters(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        noise = IndependentNoise(scale=st.loguniform(0.5, 5.0), jitter=0.2)
        like = Likelihood(GaussianFamily(), noise)
        assert like.parameters.free_names == ("scale",)
        assert like.parameters["jitter"].is_fixed
        sigma = np.sqrt((2.0 * np.asarray(observed.uncertainty)) ** 2 + 0.04)
        expected = iid_gaussian_log_prob(observed.values - predicted.values, sigma)
        assert like.log_prob(predicted, observed, {"scale": 2.0}) == pytest.approx(
            expected, abs=1e-12
        )

    def test_a_parameter_free_likelihood_needs_no_values(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        """The simple path stays simple."""
        like = Likelihood(GaussianFamily())
        assert like.parameters.free_size == 0
        assert math.isfinite(like.log_prob(predicted, observed))


# ---------------------------------------------------------------------------
# Censoring (issue #11)
# ---------------------------------------------------------------------------


@pytest.fixture
def photometry() -> PhotometricPoints:
    return PhotometricPoints(
        ["W1", "W2", "W3"],
        [3.4, 4.6, 12.0] * u.um,
        [1.0, 2.0, 0.3] * u.Jy,
        uncertainty=[0.1, 0.1, 0.1] * u.Jy,
        extra_coords={"limit_kind": np.array([0, 0, 1])},
    )


class TestVonMisesFamily:
    """The wrapped family W4.1 implemented, and why it is not a Gaussian.

    ``likelihoods.md`` §17 Q4 fixed the interface at the freeze: angles in
    radians, a *wrapped* residual, and ``kappa = 1/sigma**2`` per sample from
    the container's own uncertainties. These rows are the implementation
    meeting that.
    """

    @pytest.fixture
    def phases(self) -> ClosurePhases:
        return ClosurePhases(
            [1.0e7, 2.0e7, 3.0e7],
            [2.0e7, 1.0e7, 0.5e7],
            [-0.5e7, 1.5e7, -2.0e7],
            [1.0e7, -1.0e7, 2.5e7],
            [2.2, 2.2, 2.2] * u.um,
            [0.31, -0.12, 3.1241] * u.rad,
            uncertainty=[0.02, 0.05, 0.0873] * u.rad,
        )

    def test_it_matches_scipy_von_mises(self, phases: ClosurePhases) -> None:
        """The normalised density, against ``scipy.stats.vonmises``."""
        predicted = phases.with_values([0.30, -0.10, -3.1241])
        like = Likelihood(VonMisesFamily(), IndependentNoise())
        like.check_alignment(predicted, phases)
        kappa = 1.0 / np.asarray(phases.uncertainty) ** 2
        residual = np.angle(np.exp(1j * (phases.values - predicted.values)))
        expected = float(np.sum(st.vonmises.logpdf(residual, kappa)))
        assert like.log_prob(predicted, phases) == pytest.approx(expected, abs=1e-12)

    def test_the_residual_is_wrapped_where_a_gaussian_would_be_catastrophic(
        self, phases: ClosurePhases
    ) -> None:
        """A 2-degree error across the branch cut, charged as 2 degrees.

        ``interferometry.md`` §5 measured the alternative: an unwrapped
        residual turns the third sample's 2-degree error into a 358-degree
        one, a several-thousand-nat penalty on a triangle that fits perfectly,
        silently. No sampler recovers from that; it avoids the region of
        parameter space where the model phase is near pi.
        """
        predicted = phases.with_values([0.31, -0.12, -3.1241])
        wrapped = Likelihood(VonMisesFamily(), IndependentNoise()).log_prob(predicted, phases)
        unwrapped = Likelihood(GaussianFamily(), IndependentNoise()).log_prob(predicted, phases)
        assert wrapped > 4.0
        assert unwrapped < -2000.0

    def test_it_tends_to_the_uniform_distribution_at_large_sigma(self) -> None:
        """The normalisation is not decoration: at sigma = pi it is what keeps
        the density a density. An unnormalised wrapped Gaussian assigns *less*
        mass to a perfect match than the uniform distribution on the circle
        does, which is impossible."""
        wide = ClosurePhases(
            [1.0e7],
            [2.0e7],
            [-0.5e7],
            [1.0e7],
            [2.2] * u.um,
            [0.0] * u.rad,
            uncertainty=[100.0] * u.rad,
        )
        like = Likelihood(VonMisesFamily(), IndependentNoise())
        assert like.log_prob(wide.with_values([0.0]), wide) == pytest.approx(
            -np.log(2.0 * np.pi), abs=1e-3
        )

    def test_degrees_are_refused_at_composition_rather_than_scored(self) -> None:
        """Gap I-5's own example: the family's half of the composition-time hook."""
        degrees = ClosurePhases(
            [1.0e7],
            [2.0e7],
            [-0.5e7],
            [1.0e7],
            [2.2] * u.um,
            [18.0] * u.deg,
            uncertainty=[1.5] * u.deg,
        )
        like = Likelihood(VonMisesFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="scores angles in radians"):
            like.check_alignment(degrees.with_values([17.0]), degrees)

    def test_a_bare_array_is_taken_as_radians(self) -> None:
        bare = ClosurePhases([1.0e7], [2.0e7], [-0.5e7], [1.0e7], [2.2] * u.um, [0.2])
        VonMisesFamily().check_observed(bare)

    def test_a_correlated_noise_model_is_refused_at_composition(self) -> None:
        """A GP on a wrapped observable is latent, and this family does not consume one.

        ``Likelihood`` refuses the pair before anything is evaluated, which is
        the general rule for a ``LATENT`` combination whose family has not
        opted in. The family's own guard below is the second line of the same
        defence, for a caller assembling ``NoiseParams`` by hand.
        """
        with pytest.raises(LikelihoodError, match="CONSUMES_LATENT_GP is False"):
            Likelihood(VonMisesFamily(), GaussianProcessNoise(Matern32(0.3, 1.0e7), DenseGP()))

    def test_the_family_itself_refuses_a_correlated_evaluation(self) -> None:
        noise = NoiseParams(
            sigma=np.ones(3), values={}, kernel=Matern32(0.3, 1.0), solver=DenseGP()
        )
        with pytest.raises(LikelihoodError, match="latent-variable model"):
            VonMisesFamily().log_prob(np.zeros(3), np.zeros(3), noise)

    def test_it_draws_from_the_distribution_it_scores(self, phases: ClosurePhases) -> None:
        """The sampling form, added with the likelihood (``likelihoods.md`` §3).

        Checked distributionally rather than by value: the circular mean of
        many draws sits at the predicted angle, and the circular variance
        matches ``1 - I1(kappa)/I0(kappa)`` for the concentration the density
        uses.
        """
        family = VonMisesFamily()
        rng = np.random.default_rng(20260911)
        mean = np.full(4000, 0.7)
        sigma = np.full(4000, 0.4)
        noise = NoiseParams(sigma=sigma, values={})
        drawn = family.sample(mean, noise, rng)
        assert drawn.shape == mean.shape
        assert np.all(np.abs(drawn) <= np.pi + 1e-12)
        resultant = np.mean(np.exp(1j * drawn))
        assert np.angle(resultant) == pytest.approx(0.7, abs=0.02)
        kappa = 1.0 / 0.4**2
        expected = scipy.special.i1e(kappa) / scipy.special.i0e(kappa)
        assert abs(resultant) == pytest.approx(expected, abs=0.02)

    def test_a_draw_is_one_the_density_can_score(self, phases: ClosurePhases) -> None:
        """The rule ``simulate(observe=True)`` depends on: never hand back a
        draw the fitting likelihood refuses."""
        family = VonMisesFamily()
        rng = np.random.default_rng(7)
        noise = NoiseParams(sigma=np.asarray(phases.uncertainty), values={})
        drawn = family.sample(np.asarray(phases.values), noise, rng)
        observed = phases.with_values(drawn)
        like = Likelihood(family, IndependentNoise())
        like.check_alignment(phases, observed)
        assert np.isfinite(like.log_prob(phases, observed))

    def test_it_refuses_to_draw_under_a_correlated_noise_model(self) -> None:
        family = VonMisesFamily()
        noise = NoiseParams(
            sigma=np.ones(3), values={}, kernel=Matern32(0.3, 1.0), solver=DenseGP()
        )
        with pytest.raises(LikelihoodError, match="cannot draw"):
            family.sample(np.zeros(3), noise, np.random.default_rng(0))


class TestCensoring:
    def test_it_reads_the_hook_results_schema_left_open(
        self, photometry: PhotometricPoints
    ) -> None:
        censoring = Censoring.from_extra_coord(photometry, "limit_kind")
        assert censoring.n_samples == 3
        assert censoring.n_censored == 1
        assert censoring.any_censored
        assert repr(censoring) == "Censoring(DETECTION=2, UPPER_LIMIT=1)"

    def test_a_missing_extra_coordinate_says_what_to_do(
        self, photometry: PhotometricPoints
    ) -> None:
        with pytest.raises(LikelihoodError, match="no extra coordinate 'limits'"):
            Censoring.from_extra_coord(photometry, "limits")

    def test_a_boolean_array_must_say_which_kind_of_limit(self) -> None:
        with pytest.raises(LikelihoodError, match="A boolean array is ambiguous"):
            Censoring(np.array([True, False]))
        assert Censoring.upper_limits(np.array([False, True])).n_censored == 1
        assert Censoring.lower_limits(np.array([True, True])).n_censored == 2

    def test_unknown_codes_are_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="unknown limit codes"):
            Censoring(np.array([0, 7]))

    def test_an_upper_limit_uses_the_cdf_not_the_density(
        self, photometry: PhotometricPoints
    ) -> None:
        censoring = Censoring.from_extra_coord(photometry, "limit_kind")
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
        model = photometry.with_values([1.05, 1.95, 0.1])
        like.check_alignment(model, photometry)
        sigma = np.asarray(photometry.uncertainty)
        expected = float(np.sum(st.norm.logpdf(photometry.values[:2], model.values[:2], sigma[:2])))
        expected += float(st.norm.logcdf((photometry.values[2] - model.values[2]) / sigma[2]))
        assert like.log_prob(model, photometry) == pytest.approx(expected, abs=1e-12)

    def test_an_upper_limit_rewards_a_fainter_model(self, photometry: PhotometricPoints) -> None:
        """The point of doing it properly: below the limit is good, above is bad."""
        censoring = Censoring.from_extra_coord(photometry, "limit_kind")
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
        faint = photometry.with_values([1.0, 2.0, 0.0])
        bright = photometry.with_values([1.0, 2.0, 1.0])
        assert like.log_prob(faint, photometry) > like.log_prob(bright, photometry)

    def test_a_lower_limit_is_the_mirror_image(self) -> None:
        data = Spectrum([1.0, 2.0] * u.um, [1.0, 5.0] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy)
        censoring = Censoring(np.array([LimitKind.DETECTION, LimitKind.LOWER_LIMIT]))
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
        brighter = data.with_values([1.0, 9.0])
        fainter = data.with_values([1.0, 2.0])
        assert like.log_prob(brighter, data) > like.log_prob(fainter, data)

    def test_masking_wins_over_censoring(self, photometry: PhotometricPoints) -> None:
        """A masked limit contributes nothing at all — mask is the stronger word."""
        masked = PhotometricPoints(
            ["W1", "W2", "W3"],
            [3.4, 4.6, 12.0] * u.um,
            [1.0, 2.0, 0.3] * u.Jy,
            uncertainty=[0.1, 0.1, 0.1] * u.Jy,
            mask=np.array([False, False, True]),
        )
        censoring = Censoring(np.array([0, 0, 1]))
        censored = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
        plain = Likelihood(GaussianFamily(), IndependentNoise())
        model = masked.with_values([1.05, 1.95, 0.1])
        assert censored.log_prob(model, masked) == pytest.approx(
            plain.log_prob(model, masked), abs=1e-12
        )

    def test_censoring_composes_with_the_student_t_family(self) -> None:
        data = Spectrum([1.0, 2.0] * u.um, [1.0, 0.5] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy)
        censoring = Censoring(np.array([0, 1]))
        like = Likelihood(StudentTFamily(nu=4.0), IndependentNoise(), censoring=censoring)
        model = data.with_values([1.0, 0.2])
        expected = float(st.t.logpdf(0.0, 4.0) - math.log(0.1))
        expected += float(st.t.logcdf((0.5 - 0.2) / 0.1, 4.0))
        assert like.log_prob(model, data) == pytest.approx(expected, abs=1e-12)

    def test_a_family_that_does_not_consume_limits_refuses_them(self) -> None:
        with pytest.raises(LikelihoodError, match="does not consume a censoring declaration"):
            Likelihood(
                ComplexGaussianFamily(), IndependentNoise(), censoring=Censoring(np.array([0, 1]))
            )

    def test_censoring_plus_a_gp_is_declared_latent(self) -> None:
        """The answer to issue #11's "how the matrix algebra changes"."""
        censoring = Censoring(np.array([0, 0, 1]))
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        assert GaussianFamily().marginalisation_with(noise, censoring) is Marginalisation.LATENT

        like = Likelihood(GaussianFamily(), noise, censoring=censoring)
        assert like.marginalisation is Marginalisation.LATENT
        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um, [1.0, 2.0, 0.3] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        )
        # No family implements the truncated-latent form, so composition with
        # real data is refused rather than silently evaluated.
        with pytest.raises(LikelihoodError, match="does not implement the latent-conditional"):
            like.check_alignment(data.with_values([1.0, 2.0, 0.1]), data)
        # The family-level guard stays as a second line of defence for a caller
        # driving log_prob directly rather than through a Likelihood.
        with pytest.raises(LikelihoodError, match="orthant probability"):
            like.log_prob(data.with_values([1.0, 2.0, 0.1]), data)

    def test_a_masked_limit_does_not_force_the_latent_path(self) -> None:
        """§9's "masking beats censoring" has to hold for the declaration too.

        A limit sitting on a sample the mask excludes contributes nothing to
        the arithmetic, so it must not condemn an otherwise analytic problem to
        a gradient-based engine, nor be refused at composition.
        """
        censoring = Censoring(np.array([0, 1, 0, 0]))
        data = Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [1.0, 2.0, 3.0, 4.0] * u.Jy,
            uncertainty=np.full(4, 0.1) * u.Jy,
            mask=np.array([False, True, False, False]),
        )
        like = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)), censoring=censoring
        )
        # Declaration-time: conservative, because it has not seen the mask.
        assert like.marginalisation is Marginalisation.LATENT
        # Data-aware: the only limit is masked, so nothing is censored in fact.
        assert like.marginalisation_for(data) is Marginalisation.ANALYTIC
        like.check_alignment(data.with_values(np.zeros(4)), data)
        like.check_engine(differentiable=False, engine="emcee", observed=data)
        plain = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert like.log_prob(data.with_values(np.zeros(4)), data) == pytest.approx(
            plain.log_prob(data.with_values(np.zeros(4)), data), abs=1e-12
        )
        # Without the mask the same declaration is genuinely latent, and refused.
        exposed = Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [1.0, 2.0, 3.0, 4.0] * u.Jy,
            uncertainty=np.full(4, 0.1) * u.Jy,
        )
        assert like.marginalisation_for(exposed) is Marginalisation.LATENT
        with pytest.raises(LikelihoodError, match="does not implement the latent-conditional"):
            like.check_alignment(exposed.with_values(np.zeros(4)), exposed)

    def test_codes_and_samples_stay_aligned_on_a_gridded_container(self) -> None:
        """Excision ravels; a censoring declaration must ravel the same way.

        Nothing else in the contract pairs a 1-D per-sample array with a 2-D
        container, so this is the one place a silent index shift could hide.
        """
        codes = np.array([[0, 0, 1], [0, 1, 0]])
        image = Image(
            np.linspace(-1.0, 1.0, 2) * u.arcsec,
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.ones((2, 3)) * u.Jy,
            uncertainty=np.full((2, 3), 0.1) * u.Jy,
            extra_coords={"limit_kind": codes},
        )
        censoring = Censoring.from_extra_coord(image, "limit_kind")
        assert censoring.n_samples == 6
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=censoring)
        model = image.with_values(np.zeros((2, 3)))
        like.check_alignment(model, image)

        flat = codes.ravel()
        detected = flat == int(LimitKind.DETECTION)
        expected = float(np.sum(st.norm.logpdf(np.ones(int(detected.sum())), 0.0, 0.1)))
        expected += float(np.sum(st.norm.logcdf(np.full(int((~detected).sum()), 1.0 / 0.1))))
        assert like.log_prob(model, image) == pytest.approx(expected, abs=1e-9)

    def test_a_misaligned_declaration_is_caught(self, photometry: PhotometricPoints) -> None:
        like = Likelihood(
            GaussianFamily(), IndependentNoise(), censoring=Censoring(np.array([0, 0]))
        )
        with pytest.raises(LikelihoodError, match="covers 2 samples"):
            like.check_alignment(photometry.with_values([1.0, 2.0, 0.3]), photometry)


# ---------------------------------------------------------------------------
# Kernels and solver strategies
# ---------------------------------------------------------------------------


class TestKernelSpecification:
    def test_the_spec_is_neutral_and_serialisable(self) -> None:
        """What W1.9's lowering table consumes: a name and ordered hyperparameters."""
        spec = Matern32(1.0, 1.0).spec()
        assert spec.family == "matern32"
        assert spec.hyperparameters == ("amplitude", "length_scale")
        assert spec.quasiseparable is True
        assert spec.to_dict() == {
            "family": "matern32",
            "hyperparameters": ["amplitude", "length_scale"],
            "quasiseparable": True,
        }

    def test_the_squared_exponential_is_not_quasiseparable(self) -> None:
        assert SquaredExponential(1.0, 1.0).spec().quasiseparable is False

    def test_hyperparameters_are_ordinary_parameters_with_log_bijections(self) -> None:
        """``parameters.md`` §13's instruction to this contract, verified."""
        kernel = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(1e-2, 1e2))
        assert kernel.parameters.free_names == ("amplitude", "length_scale")
        assert kernel.parameters.bijections() == (Log(), Log())
        assert kernel.parameters.free_size == 2

    def test_a_number_fixes_a_hyperparameter(self) -> None:
        kernel = Matern32(0.5, st.loguniform(0.1, 10.0))
        assert kernel.parameters["amplitude"].is_fixed
        assert kernel.parameters.free_names == ("length_scale",)

    def test_a_ready_made_parameter_is_accepted_and_name_checked(self) -> None:
        kernel = Matern32(Parameter("amplitude", st.halfnorm(0.0, 1.0)), 1.0)
        assert kernel.parameters["amplitude"].is_free
        with pytest.raises(LikelihoodError, match="was given a Parameter named 'a'"):
            Matern32(Parameter("a", st.halfnorm(0.0, 1.0)), 1.0)

    def test_nonsense_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match=r"must be a frozen scipy\.stats distribution"):
            Matern32("big", 1.0)

    def test_a_non_positive_length_scale_is_refused_at_evaluation(self) -> None:
        kernel = Matern32(1.0, 1.0)
        with pytest.raises(LikelihoodError, match="'length_scale' must be finite and > 0"):
            kernel.matrix([[0.0]], [[1.0]], {"amplitude": 1.0, "length_scale": 0.0})

    def test_the_amplitude_is_a_standard_deviation(self) -> None:
        """``k(0) == amplitude**2`` — the convention celerite2 uses, stated."""
        kernel = Matern32(3.0, 1.0)
        variance = kernel.diagonal(np.zeros((4, 1)), {"amplitude": 3.0, "length_scale": 1.0})
        assert np.allclose(variance, 9.0)

    def test_hyperparameter_units_must_match_the_data(self, coordinates: np.ndarray) -> None:
        data = Spectrum(
            coordinates * u.um,
            np.zeros(coordinates.size) * u.Jy,
            uncertainty=np.full(coordinates.size, 0.1) * u.Jy,
        )
        wrong = GaussianProcessNoise(Matern32(1.0, 1.0, length_scale_unit=u.nm))
        with pytest.raises(LikelihoodError, match="declared in nm"):
            wrong.check_compatible(GaussianFamily(), data)
        wrong_amplitude = GaussianProcessNoise(Matern32(1.0, 1.0, amplitude_unit=u.mJy))
        with pytest.raises(LikelihoodError, match="declared in mJy"):
            wrong_amplitude.check_compatible(GaussianFamily(), data)
        right = GaussianProcessNoise(
            Matern32(1.0, 1.0, amplitude_unit=u.Jy, length_scale_unit=u.um)
        )
        right.check_compatible(GaussianFamily(), data)


class TestSolverStrategies:
    @pytest.mark.parametrize(
        "solver",
        [WindowedSparseGP(), InducingPointGP()],
    )
    def test_the_slots_are_declared_and_refuse_to_pretend(
        self, solver, predicted: Spectrum, observed: Spectrum
    ) -> None:
        assert solver.NAME
        noise = GaussianProcessNoise(Matern32(0.3, 1.0), solver)
        with pytest.raises(LikelihoodError, match="no implementation yet"):
            noise.check_compatible(GaussianFamily(), observed)
        with pytest.raises(LikelihoodError, match="no implementation yet"):
            Likelihood(GaussianFamily(), noise).log_prob(predicted, observed)

    def test_quasisep_reports_the_kernel_problem_before_the_schedule(
        self, observed: Spectrum
    ) -> None:
        """A permanent incompatibility beats a temporary one in the message."""
        noise = GaussianProcessNoise(SquaredExponential(0.3, 1.0), QuasisepGP())
        with pytest.raises(LikelihoodError, match="exact quasiseparable representation"):
            noise.check_compatible(GaussianFamily(), observed)

    def test_the_windowed_sparse_slot_is_declared_approximate(self) -> None:
        """``prior_art.md`` S2: a third, qualitatively different O(N) family."""
        assert WindowedSparseGP.EXACT is False
        assert QuasisepGP.EXACT is True
        assert DenseGP.EXACT is True

    def test_a_gridded_container_points_at_the_phase_5_slots(self) -> None:
        image = Image(
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.ones((3, 3)) * u.Jy,
            uncertainty=np.full((3, 3), 0.1) * u.Jy,
        )
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        with pytest.raises(LikelihoodError, match="SVGP / SKI / Vecchia"):
            noise.check_compatible(GaussianFamily(), image)

    def test_independent_noise_is_perfectly_happy_on_a_grid(self) -> None:
        """Only the *GP* solvers are 1D-ish; the simple path is not."""
        image = Image(
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.ones((3, 3)) * u.Jy,
            uncertainty=np.full((3, 3), 0.1) * u.Jy,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        like.check_alignment(image.with_values(np.zeros((3, 3))), image)
        expected = iid_gaussian_log_prob(np.ones(9), np.full(9, 0.1))
        assert like.log_prob(image.with_values(np.zeros((3, 3))), image) == pytest.approx(
            expected, abs=1e-12
        )

    def test_mixed_coordinate_units_are_refused(self) -> None:
        """Euclidean separation across mixed units is meaningless; say so.

        A ``VisibilitySet`` is now the standing example rather than a
        contrived one: since W4.1 it carries dimensionless ``(u, v)`` *and* a
        spectral axis, so an isotropic kernel over the whole point set is
        refused by construction. That refusal is correct — a Euclidean
        distance across wavelengths and baselines is meaningless — and the
        flexible likelihood on visibilities therefore waits for a kernel that
        can select the axes it applies to (W4.5).
        """
        vis = VisibilitySet(
            [1.0, 2.0],
            [3.0, 4.0],
            [2.2, 2.2] * u.um,
            [1.0 + 0j, 0.5 + 0j],
            uncertainty=[0.1, 0.1],
        )
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        with pytest.raises(LikelihoodError, match="different units"):
            noise.check_compatible(ComplexGaussianFamily(), vis)

    def test_a_two_axis_point_set_is_allowed_with_one_unit(self) -> None:
        """Several axes are fine; *mixed units* are the objection.

        Written against a user-defined kind since W4.1, because no shipped
        kind with more than one axis is single-unit any more.
        """

        class UVPoints(FunctionSamples):
            AXES = (
                AxisSpec("u", physical_types=("dimensionless",)),
                AxisSpec("v", physical_types=("dimensionless",)),
            )
            LAYOUT = Layout.POINTS

        points = UVPoints(
            {"u": [1.0, 2.0, 5.0], "v": [3.0, 4.0, 1.0]},
            [1.0, 0.5, 0.2],
            uncertainty=[0.1, 0.1, 0.1],
        )
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        noise.check_compatible(GaussianFamily(), points)

    def test_jitter_must_be_sane(self) -> None:
        with pytest.raises(LikelihoodError, match="jitter must be finite"):
            DenseGP(jitter=-1.0)

    def test_an_indefinite_covariance_fails_with_an_actionable_message(self) -> None:
        """Not a nan quietly propagated into a posterior: a named cause and a fix."""

        class Indefinite(Matern32):
            def matrix(self, left, right, values):  # type: ignore[no-untyped-def]
                n = np.atleast_2d(left).shape[0]
                return np.full((n, n), 2.0) + np.eye(n) * -1.0

        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um, [0.0, 0.0, 0.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        )
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Indefinite(1.0, 1.0)))
        with pytest.raises(LikelihoodError, match="not positive definite"):
            like.log_prob(data.with_values([0.0, 0.0, 0.0]), data)
        with pytest.raises(LikelihoodError, match=r"DenseGP\(jitter=\.\.\.\)"):
            like.log_prob(data.with_values([0.0, 0.0, 0.0]), data)

    def test_jitter_enters_the_diagonal_in_quadrature(self, coordinates: np.ndarray) -> None:
        data = Spectrum(
            coordinates * u.um,
            np.sin(coordinates) * u.Jy,
            uncertainty=np.full(coordinates.size, 0.1) * u.Jy,
        )
        model = data.with_values(np.zeros(coordinates.size))
        with_jitter = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP(jitter=0.3))
        )
        covariance = matern32_matrix(coordinates, 0.4, 2.0)
        covariance = covariance + np.diag(np.full(coordinates.size, 0.1**2 + 0.3**2))
        expected = multivariate_normal.logpdf(
            data.values, mean=np.zeros(coordinates.size), cov=covariance
        )
        assert with_jitter.log_prob(model, data) == pytest.approx(float(expected), abs=1e-9)


# ---------------------------------------------------------------------------
# W2.3: the quasiseparable solver on the reference path (celerite2)
# ---------------------------------------------------------------------------


class TestQuasisepGP:
    """The O(N) strategy, against the anchor it must reproduce exactly.

    The conformance battery owns the cross-solver agreement rows; these hold
    the pieces specific to *this* implementation — the exact rank-2 celerite
    representation of Matérn-3/2, the internal sort, the deferred
    leave-one-out recursion, and the refusals.
    """

    def test_the_celerite_generators_rebuild_the_kernel_matrix(
        self, coordinates: np.ndarray
    ) -> None:
        """The correspondence itself, not a log-likelihood that hides it.

        ``K[n, m] = sum_j U[n, j] V[m, j] exp(-c[j] (t_n - t_m))`` for
        ``n > m``, and ``a[n]`` on the diagonal: celerite2's own convention,
        written out here and compared with ampere's dense Matérn-3/2. If the
        amplitude/length-scale correspondence were wrong, this fails first and
        by name.
        """
        from ampere.core.likelihood import _celerite_term

        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        term = _celerite_term(kernel, values)
        c, a, lower, right = term.get_celerite_matrices(coordinates, np.zeros(coordinates.size))
        assert c.shape == (2,) and lower.shape == (coordinates.size, 2)

        rebuilt = np.zeros((coordinates.size, coordinates.size))
        for i in range(coordinates.size):
            rebuilt[i, i] = a[i]
            for j in range(i):
                entry = float(
                    np.sum(lower[i] * right[j] * np.exp(-c * (coordinates[i] - coordinates[j])))
                )
                rebuilt[i, j] = rebuilt[j, i] = entry
        np.testing.assert_allclose(rebuilt, matern32_matrix(coordinates, 0.4, 2.0), atol=1e-14)
        # k(0) == amplitude**2: the standard-deviation convention, in the
        # generator the solver actually factorises.
        assert a[0] == pytest.approx(0.4**2, abs=1e-15)

    def test_it_reproduces_the_dense_marginal_likelihood(
        self, predicted: Spectrum, observed: Spectrum
    ) -> None:
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP()))
        quasisep = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), QuasisepGP())
        )
        assert quasisep.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=1e-9
        )

    def test_the_jitter_matches_the_dense_solver_s(self, coordinates: np.ndarray) -> None:
        data = Spectrum(
            coordinates * u.um,
            np.sin(coordinates) * u.Jy,
            uncertainty=np.full(coordinates.size, 0.1) * u.Jy,
        )
        model = data.with_values(np.zeros(coordinates.size))
        kernel = Matern32(0.4, 2.0)
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP(jitter=0.3)))
        quasisep = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), QuasisepGP(jitter=0.3))
        )
        assert quasisep.log_prob(model, data) == pytest.approx(
            dense.log_prob(model, data), abs=1e-9
        )
        with pytest.raises(LikelihoodError, match="jitter must be finite"):
            QuasisepGP(jitter=-1.0)

    def test_unordered_coordinates_are_sorted_internally(self, coordinates: np.ndarray) -> None:
        """A Gaussian density is permutation invariant; celerite2 is not.

        ``PhotometricPoints`` declares ``Order.ANY``, so a one-axis container
        can legitimately arrive unsorted. The solver must permute, solve and
        permute back rather than hand celerite2 something it will reject.
        """
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        rng = np.random.default_rng(20260905)
        residual = rng.normal(0.0, 0.3, coordinates.size)
        variance = np.full(coordinates.size, 0.1**2)
        shuffle = rng.permutation(coordinates.size)

        ordered = QuasisepGP().log_marginal_likelihood(
            kernel, coordinates[:, None], residual, variance, values
        )
        shuffled = QuasisepGP().log_marginal_likelihood(
            kernel,
            coordinates[shuffle][:, None],
            residual[shuffle],
            variance[shuffle],
            values,
        )
        expected = DenseGP().log_marginal_likelihood(
            kernel, coordinates[:, None], residual, variance, values
        )
        assert ordered == pytest.approx(expected, abs=1e-9)
        assert shuffled == pytest.approx(expected, abs=1e-9)

    def test_the_conditioned_gp_matches_the_dense_one(
        self, predicted: Spectrum, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP()))
        quasisep = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), QuasisepGP())
        )
        at = np.linspace(coordinates[0] - 0.5, coordinates[-1] + 0.5, 23)
        for kwargs in ({}, {"at": at}):
            expected = dense.conditional(predicted, observed, **kwargs)
            got = quasisep.conditional(predicted, observed, **kwargs)
            np.testing.assert_allclose(got.mean, expected.mean, atol=1e-10)
            np.testing.assert_allclose(got.variance, expected.variance, atol=1e-10)

    def test_condition_at_accepts_a_multi_axis_container(self) -> None:
        """W5.28(d): ``at=`` must match the container's own axis count, not 1.

        An ``axes=(...)`` kernel already trains fine on a multi-axis
        container -- ``_axis`` reduces it to the one ordered column through
        ``kernel.select``. Before this fix ``condition(at=...)`` hardcoded
        ``dimensions=1`` for the *target*, so an ``at`` shaped like the
        training container itself (the only shape ``DenseGP.condition``
        accepts, and the one a caller matching its own coordinates would
        naturally pass) was refused as a shape mismatch instead of accepted.
        """
        kernel = Matern32(0.4, 2.0, axes=("t",)).for_axes(["t", "other"])
        coordinates = np.column_stack([np.linspace(0.0, 5.0, 9), np.zeros(9)])
        rng = np.random.default_rng(20260922)
        residual = rng.normal(0.0, 0.3, 9)
        variance = np.full(9, 0.05**2)
        values = kernel.resolve(None)
        at = np.column_stack([np.linspace(0.5, 4.5, 6), np.full(6, 3.0)])

        dense = DenseGP().condition(kernel, coordinates, residual, variance, values, at=at)
        quasisep = QuasisepGP().condition(kernel, coordinates, residual, variance, values, at=at)
        np.testing.assert_allclose(quasisep.mean, dense.mean, atol=1e-9)
        np.testing.assert_allclose(quasisep.variance, dense.variance, atol=1e-9)

        # A target that does not match the training container's own axis
        # count is still refused, by name -- W5.28(d) widens the accepted
        # shape to match DenseGP, it does not drop the check.
        with pytest.raises(LikelihoodError, match=r"coordinate\(s\) per point"):
            QuasisepGP().condition(kernel, coordinates, residual, variance, values, at=at[:, :1])

    def test_the_whitening_transform_factorises_the_kernel(self, coordinates: np.ndarray) -> None:
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        basis = np.eye(coordinates.size)
        lower = np.column_stack(
            [
                QuasisepGP().latent_transform(kernel, coordinates[:, None], basis[:, i], values)
                for i in range(coordinates.size)
            ]
        )
        np.testing.assert_allclose(np.triu(lower, k=1), 0.0, atol=1e-15)
        np.testing.assert_allclose(
            lower @ lower.T, matern32_matrix(coordinates, 0.4, 2.0), atol=1e-9
        )

    def test_the_leave_one_out_recursion_is_deferred_and_says_so(self) -> None:
        """W2.3's recorded deferral, refusing rather than costing O(N**2)."""
        with pytest.raises(LikelihoodError, match=r"DEVELOPMENT_PLAN\.md"):
            QuasisepGP().conditional_loo(
                Matern32(0.3, 1.5), np.array([[1.0]]), np.zeros(1), np.ones(1), {}
            )

    def test_a_quasiseparable_kernel_with_no_celerite_term_is_refused(
        self, observed: Spectrum
    ) -> None:
        """The declaration is not enough; ampere must hold the representation."""

        class Unregistered(Matern32):
            FAMILY = "no_such_family"

        noise = GaussianProcessNoise(Unregistered(0.3, 1.0), QuasisepGP())
        with pytest.raises(LikelihoodError, match="no exact celerite representation"):
            noise.check_compatible(GaussianFamily(), observed)

    def test_two_coordinate_axes_are_refused_at_the_solver_too(self) -> None:
        """``check_compatible`` catches it at composition; a direct call here."""
        kernel = Matern32(0.4, 2.0)
        points = np.array([[1.0, 2.0], [2.0, 3.0], [3.0, 4.0]])
        with pytest.raises(LikelihoodError, match="one ordered coordinate per sample"):
            QuasisepGP().log_marginal_likelihood(
                kernel, points, np.zeros(3), np.full(3, 0.01), kernel.resolve(None)
            )

    def test_an_indefinite_system_names_its_own_jitter(self) -> None:
        """Duplicated coordinates: singular in K, and refused rather than NaN."""
        kernel = Matern32(1.0, 1.0)
        duplicated = np.array([[1.0], [1.0], [2.0]])
        with pytest.raises(LikelihoodError, match=r"QuasisepGP\(jitter=\.\.\.\)"):
            QuasisepGP().log_marginal_likelihood(
                kernel, duplicated, np.zeros(3), np.zeros(3), kernel.resolve(None)
            )
        # And the same system is fine once the jitter separates the rows.
        assert math.isfinite(
            QuasisepGP(jitter=0.1).log_marginal_likelihood(
                kernel, duplicated, np.zeros(3), np.zeros(3), kernel.resolve(None)
            )
        )

    def test_a_negative_variance_is_refused_rather_than_returning_nan(self) -> None:
        """celerite2 returns NaN where a Cholesky raises; the guard is ours."""
        kernel = Matern32(1.0, 1.0)
        points = np.array([[1.0], [2.0], [3.0]])
        with pytest.raises(LikelihoodError, match="not a covariance matrix"):
            QuasisepGP().log_marginal_likelihood(
                kernel, points, np.zeros(3), np.full(3, -1.0), kernel.resolve(None)
            )

    def test_the_spec_records_the_strategy_and_its_configuration(self) -> None:
        """``to_spec`` must distinguish the two solvers — ``results.md`` §15 R7."""
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP()))
        quasisep = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), QuasisepGP(jitter=0.1))
        )
        assert dense.to_spec()["solver"]["name"] == "DenseGP"
        assert quasisep.to_spec()["solver"] == {
            "name": "QuasisepGP",
            "class": "QuasisepGP",
            "exact": True,
            "config": {"jitter": 0.1},
        }

    def test_importing_ampere_core_does_not_import_celerite2(self) -> None:
        """The lazy-import discipline: celerite2 arrives on use, not on import."""
        import subprocess
        import sys

        probe = "import sys; import ampere.core; print('celerite2' in sys.modules)"
        result = subprocess.run(
            [sys.executable, "-c", probe], capture_output=True, text=True, check=True
        )
        assert result.stdout.strip().endswith("False")


# ---------------------------------------------------------------------------
# Composition-time checking
# ---------------------------------------------------------------------------


class TestComposition:
    def test_the_parameter_namespace_is_flat(self) -> None:
        """One ParameterSet, no nested merge — ``parameters.md`` §12.4."""
        gp = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(
                Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
                scale=st.loguniform(0.5, 2.0),
                jitter=st.halfnorm(0.0, 1.0),
            ),
        )
        assert gp.parameters.free_names == ("amplitude", "length_scale", "scale", "jitter")
        assert gp.parameters.free_size == 4
        # A family's own parameters join the same flat set, after the noise
        # model's, in declaration order.
        with_family = Likelihood(
            StudentTFamily(nu=st.loguniform(2.0, 50.0)),
            IndependentNoise(scale=st.loguniform(0.5, 2.0)),
        )
        assert with_family.parameters.free_names == ("scale", "nu")

    def test_a_name_collision_is_loud(self) -> None:
        @register_family
        class ScaledFamily(LikelihoodFamily):
            NAME = "test_scaled"

            def __init__(self) -> None:
                self.register_parameter(Parameter("scale", st.norm(0.0, 1.0)))

            def log_prob(
                self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
            ) -> float:
                return 0.0

        try:
            with pytest.raises(LikelihoodError, match="both declare a parameter named 'scale'"):
                Likelihood(ScaledFamily(), IndependentNoise(scale=st.loguniform(0.5, 2.0)))
        finally:
            from ampere.core import likelihood as module

            module._FAMILIES.pop("test_scaled", None)

    def test_a_shape_mismatch_is_caught(self, observed: Spectrum) -> None:
        short = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy)
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="check_alignment"):
            like.log_prob(short, observed)

    def test_a_unit_mismatch_is_caught_at_composition(self, observed: Spectrum) -> None:
        like = Likelihood(GaussianFamily(), IndependentNoise())
        in_mjy = observed.to_unit(u.mJy)
        with pytest.raises(LikelihoodError, match="Convert once at composition time"):
            like.check_alignment(in_mjy, observed)

    def test_a_kind_mismatch_is_caught_at_composition(self, observed: Spectrum) -> None:
        times = TimeSeries(
            np.arange(float(observed.n_samples)) * u.day,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="compares like with like"):
            like.check_alignment(times, observed)

    def test_differing_coordinates_are_caught_at_composition(
        self, observed: Spectrum, coordinates: np.ndarray
    ) -> None:
        shifted = Spectrum(
            (coordinates + 0.5) * u.um,
            observed.values * u.Jy,
            uncertainty=np.asarray(observed.uncertainty) * u.Jy,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="axes differ"):
            like.check_alignment(shifted, observed)

    def test_a_family_needing_uncertainties_says_so(self, coordinates: np.ndarray) -> None:
        bare = Spectrum(coordinates * u.um, np.zeros(coordinates.size) * u.Jy)
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="needs per-sample uncertainties"):
            like.check_alignment(bare.with_values(np.zeros(coordinates.size)), bare)

    def test_poisson_needs_none(self) -> None:
        counts = Spectrum([1.0, 2.0] * u.um, [3.0, 4.0])
        like = Likelihood(PoissonFamily(), IndependentNoise())
        like.check_alignment(counts.with_values([2.5, 4.5]), counts)

    def test_non_finite_data_are_refused_rather_than_propagated(
        self, coordinates: np.ndarray
    ) -> None:
        """No sentinel NaNs — ``results_schema.md`` §7's rule, enforced here."""
        values = np.zeros(coordinates.size)
        values[3] = np.nan
        data = Spectrum(
            coordinates * u.um,
            values * u.Jy,
            uncertainty=np.full(coordinates.size, 0.1) * u.Jy,
        )
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="non-finite entries"):
            like.log_prob(data.with_values(np.zeros(coordinates.size)), data)

    def test_a_user_defined_container_kind_works_unchanged(self) -> None:
        """§4.3's extensibility requirement, from this contract's side.

        ``results_schema.md`` §13 makes a new container kind three class
        attributes and no changes inside ampere. That promise is only worth
        anything if the likelihood contract consumes one without knowing it
        exists, GP and all.
        """

        class PolarisationCurve(FunctionSamples):
            AXES = (
                AxisSpec(
                    "spectral_axis",
                    physical_types=("length",),
                    order=Order.STRICTLY_INCREASING,
                ),
            )
            LAYOUT = Layout.POINTS

        curve = PolarisationCurve(
            {"spectral_axis": [0.4, 0.6, 0.8] * u.um},
            [0.02, 0.03, 0.01],
            uncertainty=[0.002, 0.002, 0.002],
        )
        model = curve.with_values([0.021, 0.028, 0.011])
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.005, 0.2)))
        like.check_alignment(model, curve)

        coordinates = np.array([0.4, 0.6, 0.8])
        covariance = matern32_matrix(coordinates, 0.005, 0.2) + np.diag(np.full(3, 0.002**2))
        expected = multivariate_normal.logpdf(
            curve.values - model.values, mean=np.zeros(3), cov=covariance
        )
        assert like.log_prob(model, curve) == pytest.approx(float(expected), abs=1e-9)
        assert like.conditional(model, curve).mean.size == 3

    def test_a_container_with_a_different_axis_signature_is_refused(self) -> None:
        """A subclass is comparable with its base only if it kept the axes."""

        class TwoAxisSpectrum(Spectrum):
            AXES = (
                AxisSpec("spectral_axis", physical_types=("length",), order=Order.ANY),
                AxisSpec("epoch", physical_types=("time",), order=Order.ANY),
            )

            __init__ = FunctionSamples.__init__  # bypass Spectrum's one-axis signature

        odd = TwoAxisSpectrum(
            {"spectral_axis": [1.0, 2.0] * u.um, "epoch": [0.0, 1.0] * u.day},
            [1.0, 1.0],
        )
        plain = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy, uncertainty=[0.1, 0.1] * u.Jy)
        like = Likelihood(GaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="compares like with like"):
            like.check_alignment(odd, plain)

    def test_non_finite_complex_data_are_refused_too(self) -> None:
        """The complex path must not be the one that leaks a NaN into a posterior."""
        vis = VisibilitySet(
            [1.0, 2.0, 3.0],
            [4.0, 5.0, 6.0],
            [2.2, 2.2, 2.2] * u.um,
            [1.0 + 0.2j, 0.6 - 0.3j, complex(np.nan, 1.0)],
            uncertainty=[0.05, 0.05, 0.05],
        )
        like = Likelihood(ComplexGaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="non-finite entries"):
            like.log_prob(vis.with_values([1.0 + 0j, 0.6 + 0j, 0.4 + 0j]), vis)

    def test_an_overflowing_amplitude_raises_this_contracts_error(self) -> None:
        """Not a bare ValueError: W1.7 is told to catch ``LikelihoodError``."""
        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um, [0.0, 0.0, 0.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        )
        model = data.with_values([0.0, 0.0, 0.0])
        # amplitude**2 overflows float64 while amplitude itself is finite, so
        # the support check passes and the covariance is full of infs.
        huge = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(1e200, 1.0)))
        with pytest.raises(LikelihoodError, match="non-finite entries"):
            huge.log_prob(model, data)
        with pytest.raises(LikelihoodError, match="non-finite entries"):
            huge.conditional(model, data)
        with pytest.raises(LikelihoodError, match="non-finite entries"):
            DenseGP().latent_transform(
                Matern32(1e200, 1.0),
                np.array([[1.0], [2.0], [3.0]]),
                np.zeros(3),
                {"amplitude": 1e200, "length_scale": 1.0},
            )

    def test_a_correlated_noise_model_must_supply_its_own_noise_params(self) -> None:
        """Two notions of "correlated" must not be able to disagree.

        A third-party subclass that sets ``CORRELATED = True`` but inherits the
        base ``noise_params`` would report ``NoiseParams.correlated is False``,
        so every family would take its uncorrelated branch and the correlations
        would silently do nothing.
        """

        class HalfDeclaredNoise(NoiseModel):
            CORRELATED = True

            def sigma(
                self,
                observed: FunctionSamples,
                retain: np.ndarray,
                values: dict[str, object],
                *,
                predicted: np.ndarray | None = None,
            ) -> np.ndarray:
                return np.asarray(observed.uncertainty).ravel()[retain]

        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um, [1.0, 2.0, 3.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        )
        like = Likelihood(GaussianFamily(), HalfDeclaredNoise())
        with pytest.raises(LikelihoodError, match="declares CORRELATED = True but inherits"):
            like.log_prob(data.with_values([0.0, 0.0, 0.0]), data)

    def test_the_repr_states_the_marginalisation(self) -> None:
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert repr(like) == (
            "Likelihood(PoissonFamily, GaussianProcessNoise, marginalisation=latent)"
        )

    def test_every_error_here_is_a_contract_error(self) -> None:
        """Catchable by the family, and by the builtin, per ``exceptions.py``."""
        assert issubclass(LikelihoodError, ContractError)
        assert issubclass(LikelihoodError, ValueError)


class TestCompositionTimeDataChecks:
    """The three 2026-09-02 rulings on check_alignment and NoiseParams.

    W1.11 gaps I-1 (value-dtype comparison), I-5 with X-3 (the family's
    composition-time hook, and Poisson's integrality test moving into it) and
    X-2 (``NoiseParams.retain``), all approved as proposed.
    """

    def test_a_complex_prediction_is_refused_against_real_amplitudes(self) -> None:
        # I-1: VisibilitySet is the only kind legal with both dtypes, so only
        # there could a complex prediction silently be fitted against |V|.
        uv, vv = np.array([10.0, 20.0, 30.0]), np.array([5.0, 15.0, 25.0])
        wave = np.full(3, 2.2) * u.um
        complex_prediction = VisibilitySet(uv, vv, wave, (np.ones(3) + 0.5j) * u.Jy)
        real_amplitudes = VisibilitySet(
            uv, vv, wave, np.abs(np.ones(3) + 0.5j) * u.Jy, uncertainty=0.1 * np.ones(3) * u.Jy
        )
        like = Likelihood(ComplexGaussianFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="take the modulus in the instrument chain"):
            like.check_alignment(complex_prediction, real_amplitudes)

    def test_check_observed_is_called_at_composition(self) -> None:
        # I-5: the family's half of NoiseModel.check_compatible's obligation.
        class Circular(LikelihoodFamily):
            NAME = "test_circular"

            def log_prob(
                self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
            ) -> float:
                return 0.0

            def check_observed(self, observed: FunctionSamples) -> None:
                if np.any(np.abs(np.asarray(observed.values)) > np.pi):
                    raise LikelihoodError("a circular family needs angles in radians.")

        degrees = Spectrum([1.0, 2.0] * u.um, [170.0, -50.0], uncertainty=[5.0, 5.0])
        like = Likelihood(Circular(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="angles in radians"):
            like.check_alignment(degrees.with_values([171.0, -49.0]), degrees)

    def test_poisson_integrality_is_checked_at_composition_not_per_draw(self) -> None:
        # X-3, folded into I-5: an O(N) property of the data belongs in
        # check_alignment, not in the hot loop.
        rates = Spectrum([1.0, 2.0, 3.0] * u.um, [3.5, 6.0, 2.5])
        fractional = Spectrum([1.0, 2.0, 3.0] * u.um, [4.0, 6.51, 2.0])
        like = Likelihood(PoissonFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="not integral"):
            like.check_alignment(rates, fractional)
        # A masked non-integer is exempt: a masked sample carries zero
        # information, so it cannot fail a precondition.
        masked = Spectrum(
            [1.0, 2.0, 3.0] * u.um,
            [4.0, 6.51, 2.0],
            mask=np.array([False, True, False]),
        )
        like.check_alignment(rates, masked)
        # The rate > 0 guard stays per draw: it is a property of the prediction.
        counts = Spectrum([1.0, 2.0, 3.0] * u.um, [4.0, 7.0, 2.0])
        with pytest.raises(LikelihoodError, match="strictly positive expected count"):
            like.log_prob(counts.with_values([3.5, -1.0, 2.5]), counts)

    def test_noise_params_carry_the_retain_indicator(self) -> None:
        # X-2: a family with its own aligned per-sample data excises it with
        # noise.retain, full-container length, and matches deletion.
        captured: dict[str, np.ndarray | None] = {}

        class Background(LikelihoodFamily):
            NAME = "test_background"

            def __init__(self, background: np.ndarray) -> None:
                self.background = np.asarray(background, dtype=float)

            def log_prob(
                self, predicted: np.ndarray, observed: np.ndarray, noise: NoiseParams
            ) -> float:
                captured["retain"] = noise.retain
                assert noise.retain is not None
                aligned = self.background[noise.retain]
                assert aligned.shape == observed.shape
                residual = (observed - predicted - aligned) / noise.sigma
                return float(np.sum(-0.5 * residual**2 - np.log(noise.sigma)))

        data = Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [3.0, 2.5, 2.2, 1.4] * u.Jy,
            uncertainty=[0.1, 0.1, 0.1, 0.1] * u.Jy,
            mask=np.array([False, True, False, False]),
        )
        background = np.array([0.5, 99.0, 0.4, 0.3])
        like = Likelihood(Background(background), IndependentNoise())
        masked_value = like.log_prob(data.with_values([2.6, 0.0, 1.9, 1.2]), data)
        assert captured["retain"] is not None
        assert captured["retain"].tolist() == [True, False, True, True]
        excised = Spectrum(
            [1.0, 3.0, 4.0] * u.um,
            [3.0, 2.2, 1.4] * u.Jy,
            uncertainty=[0.1, 0.1, 0.1] * u.Jy,
        )
        deleted = Likelihood(Background(background[[0, 2, 3]]), IndependentNoise())
        unmasked_value = deleted.log_prob(excised.with_values([2.6, 1.9, 1.2]), excised)
        assert masked_value == pytest.approx(unmasked_value, abs=1e-12)


class TestLikelihoodToSpec:
    """R7's promotion (ruled 2026-09-03): the object that knows itself describes itself."""

    def test_the_spec_distinguishes_what_parameter_specs_cannot(self) -> None:
        # The correctness argument R7 was granted on: same parameters,
        # different kernel family or solver configuration, different spec.
        matern = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0)))
        rbf = Likelihood(GaussianFamily(), GaussianProcessNoise(SquaredExponential(0.3, 2.0)))
        assert matern.parameters.to_spec() == rbf.parameters.to_spec()
        assert matern.to_spec()["kernel"]["family"] == "matern32"
        assert rbf.to_spec()["kernel"]["family"] == "squared_exponential"
        jittered = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), DenseGP(jitter=0.1))
        )
        assert matern.to_spec()["solver"]["config"] == {"jitter": 0.0}
        assert jittered.to_spec()["solver"]["config"] == {"jitter": 0.1}

    def test_the_spec_is_declarative_and_censoring_carries_counts_only(self) -> None:
        codes = np.array([0, 0, int(LimitKind.UPPER_LIMIT)], dtype=np.int8)
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))
        spec = like.to_spec()
        assert spec["family"] == "gaussian"
        assert spec["marginalisation"] == "analytic"
        assert spec["censoring"] == {"n_samples": 3, "n_censored": 1}
        # Per-sample content (the code positions) is provenance's business,
        # deliberately absent from the declarative spec.
        assert "kinds" not in spec["censoring"]
        import json

        json.dumps(spec)  # the whole mapping is JSON-able

    def test_provenance_composes_the_spec_rather_than_reassembling_it(self) -> None:
        from ampere.results.provenance import describe_likelihood

        codes = np.array([0, 0, int(LimitKind.UPPER_LIMIT)], dtype=np.int8)
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))
        described = describe_likelihood(like)
        spec = like.to_spec()
        for key, value in spec.items():
            if key == "censoring":
                continue
            assert described[key] == value
        assert "buffers" in described  # the content-fingerprint layer's addition
        assert "kinds" in described["censoring"]


class TestPointwiseLogProb:
    """The §4.4 addition results.md §15 R2 granted (ruled 2026-09-03).

    Factorised terms for independent noise (exact: they sum to log_prob);
    leave-one-out conditionals for a GP (a different decomposition, which
    deliberately does not sum to the joint value).
    """

    @staticmethod
    def _data(mask: object = None) -> Spectrum:
        return Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [3.0, 2.5, 2.2, 1.4] * u.Jy,
            uncertainty=[0.1, 0.2, 0.1, 0.3] * u.Jy,
            mask=mask,
        )

    def test_factorised_terms_sum_to_log_prob(self) -> None:
        data = self._data(mask=np.array([False, True, False, False]))
        model = data.with_values([3.05, 0.0, 2.25, 1.35])
        like = Likelihood(GaussianFamily(), IndependentNoise())
        terms = like.pointwise_log_prob(model, data)
        assert terms.shape == (3,)  # masked sample excised
        assert float(np.sum(terms)) == pytest.approx(like.log_prob(model, data), abs=1e-12)
        # And each term is the closed form for its own sample.
        sigma = np.array([0.1, 0.1, 0.3])
        residual = np.array([3.0 - 3.05, 2.2 - 2.25, 1.4 - 1.35])
        np.testing.assert_allclose(terms, st.norm(0.0, sigma).logpdf(residual), atol=1e-12)

    def test_a_censored_sample_contributes_its_tobit_term(self) -> None:
        data = self._data()
        model = data.with_values([3.05, 2.40, 2.25, 1.35])
        codes = np.array([0, 0, 0, int(LimitKind.UPPER_LIMIT)], dtype=np.int8)
        like = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))
        terms = like.pointwise_log_prob(model, data)
        assert float(np.sum(terms)) == pytest.approx(like.log_prob(model, data), abs=1e-12)
        assert terms[-1] == pytest.approx(float(st.norm.logcdf((1.4 - 1.35) / 0.3)), abs=1e-12)

    def test_a_parameterised_family_and_noise_resolve_the_same_values(self) -> None:
        data = self._data()
        model = data.with_values([3.05, 2.40, 2.25, 1.35])
        like = Likelihood(
            StudentTFamily(nu=st.loguniform(2.0, 50.0)),
            IndependentNoise(scale=st.loguniform(0.5, 5.0)),
        )
        values = {"nu": 4.0, "scale": 2.0}
        terms = like.pointwise_log_prob(model, data, values)
        assert float(np.sum(terms)) == pytest.approx(like.log_prob(model, data, values), abs=1e-12)

    def test_a_family_with_its_own_aligned_data_excises_per_term(self) -> None:
        # X-2's retain contract holds per single-sample call: the family sees
        # a full-length indicator selecting exactly that sample.
        class Background(LikelihoodFamily):
            NAME = "test_pointwise_background"

            def __init__(self, background: np.ndarray) -> None:
                self.background = background

            def log_prob(self, predicted, observed, noise):
                assert noise.retain is not None
                aligned = self.background[noise.retain]
                residual = (observed - predicted - aligned) / noise.sigma
                return float(np.sum(-0.5 * residual**2 - np.log(noise.sigma)))

        data = self._data(mask=np.array([False, True, False, False]))
        model = data.with_values([2.6, 0.0, 1.9, 1.2])
        like = Likelihood(Background(np.array([0.5, 99.0, 0.4, 0.3])), IndependentNoise())
        terms = like.pointwise_log_prob(model, data)
        assert float(np.sum(terms)) == pytest.approx(like.log_prob(model, data), abs=1e-12)

    def test_gp_terms_are_the_leave_one_out_conditionals(self) -> None:
        data = self._data()
        model = data.with_values([3.05, 2.40, 2.25, 1.35])
        kernel_values = {"amplitude": 0.3, "length_scale": 1.5}
        like = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.5), DenseGP()))
        terms = like.pointwise_log_prob(model, data)
        # Brute force: for each i, condition on the other samples explicitly.
        coordinates = np.array([[1.0], [2.0], [3.0], [4.0]])
        kernel = Matern32(0.3, 1.5)
        covariance = kernel.matrix(coordinates, coordinates, kernel_values) + np.diag(
            np.array([0.1, 0.2, 0.1, 0.3]) ** 2
        )
        residual = np.asarray(model.values - data.values, dtype=float) * -1.0
        expected = np.empty(4)
        for i in range(4):
            others = [j for j in range(4) if j != i]
            solve = np.linalg.solve(covariance[np.ix_(others, others)], residual[others])
            mean = covariance[i, others] @ solve
            var = covariance[i, i] - covariance[i, others] @ np.linalg.solve(
                covariance[np.ix_(others, others)], covariance[others, i]
            )
            expected[i] = st.norm(mean, np.sqrt(var)).logpdf(residual[i])
        np.testing.assert_allclose(terms, expected, atol=1e-9)
        # A different decomposition: the LOO terms do not sum to the joint.
        assert float(np.sum(terms)) != pytest.approx(like.log_prob(model, data), abs=1e-6)

    def test_a_solver_without_the_recursion_refuses(self) -> None:
        with pytest.raises(LikelihoodError, match="conditional_loo"):
            QuasisepGP().conditional_loo(
                Matern32(0.3, 1.5), np.array([[1.0]]), np.zeros(1), np.ones(1), {}
            )

    def test_a_latent_combination_is_refused(self) -> None:
        counts = Spectrum([1.0, 2.0, 3.0] * u.um, [4.0, 7.0, 2.0])
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        with pytest.raises(LikelihoodError, match="latent"):
            like.pointwise_log_prob(counts.with_values([4.0, 7.0, 2.0]), counts)


class TestPredictionAwareNoise:
    """The X-1 ruling (2026-09-03): a noise model sees the prediction.

    ``NoiseModel.sigma`` and ``noise_params`` take the retained predicted
    values as a keyword-only ``predicted`` argument, passed at every call
    site — ``awkward_instrument.md`` §6's detailed design, accepted as
    written. The conformance battery holds the arithmetic
    (``tests/conformance/test_likelihoods.py``); these tests hold the
    plumbing: the argument arrives, excised, at each call site.
    """

    @staticmethod
    def _recording_independent(captured: dict) -> NoiseModel:
        class Recording(IndependentNoise):
            def sigma(self, observed, retain, values, *, predicted=None):
                captured["predicted"] = predicted
                return super().sigma(observed, retain, values, predicted=predicted)

        return Recording()

    def test_log_prob_passes_the_retained_predicted_values(self) -> None:
        captured: dict = {}
        data = Spectrum(
            [1.0, 2.0, 3.0, 4.0] * u.um,
            [3.0, 2.5, 2.2, 1.4] * u.Jy,
            uncertainty=[0.1, 0.1, 0.1, 0.1] * u.Jy,
            mask=np.array([False, True, False, False]),
        )
        like = Likelihood(GaussianFamily(), self._recording_independent(captured))
        like.log_prob(data.with_values([2.6, 0.0, 1.9, 1.2]), data)
        # Excised exactly as the family's own first argument is: the masked
        # second sample is gone.
        assert captured["predicted"] is not None
        np.testing.assert_array_equal(captured["predicted"], [2.6, 1.9, 1.2])

    def test_conditional_passes_the_same_prediction_the_fit_used(self) -> None:
        captured: dict = {}

        class RecordingGP(GaussianProcessNoise):
            def sigma(self, observed, retain, values, *, predicted=None):
                captured["predicted"] = predicted
                return super().sigma(observed, retain, values, predicted=predicted)

        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um,
            [1.0, 2.0, 3.0] * u.Jy,
            uncertainty=[0.1, 0.1, 0.1] * u.Jy,
        )
        like = Likelihood(GaussianFamily(), RecordingGP(Matern32(0.3, 1.0)))
        like.conditional(data.with_values([0.9, 1.8, 2.7]), data)
        assert captured["predicted"] is not None
        np.testing.assert_array_equal(captured["predicted"], [0.9, 1.8, 2.7])

    def test_a_fractional_model_noise_composes_with_the_gaussian_family(self) -> None:
        """The motivating case: sigma_eff**2 = sigma_data**2 + (f * mu)**2."""

        class FractionalModelNoise(NoiseModel):
            def __init__(self, f: float) -> None:
                self._f = float(f)

            def sigma(self, observed, retain, values, *, predicted=None):
                base = np.asarray(observed.uncertainty).ravel()[retain]
                assert predicted is not None
                return np.sqrt(base**2 + (self._f * predicted) ** 2)

        fraction = 0.1
        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um,
            [1.1, 2.1, 2.9] * u.Jy,
            uncertainty=[0.1, 0.2, 0.3] * u.Jy,
        )
        prediction = np.array([1.0, 2.0, 3.0])
        like = Likelihood(GaussianFamily(), FractionalModelNoise(fraction))
        # Diagonal noise: the declaration stays ANALYTIC however sigma was
        # computed (X-1 point 6).
        assert like.marginalisation is Marginalisation.ANALYTIC
        sigma_eff = np.sqrt(np.array([0.1, 0.2, 0.3]) ** 2 + (fraction * prediction) ** 2)
        residual = np.array([1.1, 2.1, 2.9]) - prediction
        expected = float(np.sum(st.norm(0.0, sigma_eff).logpdf(residual)))
        assert like.log_prob(data.with_values(prediction), data) == pytest.approx(
            expected, abs=1e-12
        )
