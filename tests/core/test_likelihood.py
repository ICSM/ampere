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
import scipy.stats as st
from scipy.stats import multivariate_normal

from ampere.core import (
    AxisSpec,
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

    @pytest.mark.parametrize(
        "family", [StudentTFamily(), CauchyFamily(), PoissonFamily(), ComplexGaussianFamily()]
    )
    def test_every_non_gaussian_family_goes_latent_under_a_gp(
        self, family: LikelihoodFamily
    ) -> None:
        like = Likelihood(family, GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert like.marginalisation is Marginalisation.LATENT

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
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        rate = counts.with_values([3.5, 6.0, 2.5, 8.0])
        with pytest.raises(LikelihoodError, match="latent-variable model"):
            like.log_prob(rate, counts)
        latent = np.array([0.1, -0.2, 0.05, 0.0])
        assert like.log_prob(rate, counts, latent=latent) == pytest.approx(
            float(np.sum(st.poisson.logpmf(counts.values, rate.values * np.exp(latent)))),
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
        counts = Spectrum([1.0, 2.0] * u.um, [4.5, 7.0])
        like = Likelihood(PoissonFamily(), IndependentNoise())
        with pytest.raises(LikelihoodError, match="non-negative integer counts"):
            like.log_prob(counts.with_values([4.0, 7.0]), counts)

    def test_latent_parameter_rejects_a_nonsense_size(self) -> None:
        with pytest.raises(LikelihoodError, match="positive number of latent values"):
            latent_parameter("latent", 0)


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

    @pytest.mark.parametrize("family", [RiceFamily(), VonMisesFamily()])
    def test_declared_but_unimplemented_families_refuse_composition(
        self, family: LikelihoodFamily
    ) -> None:
        """Declared in the registry, refused at composition — never silently wrong."""
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
        vis = VisibilitySet([1.0, 2.0], [3.0, 4.0], [1.0 + 0j, 0.5 + 0j], uncertainty=[0.1, 0.1])
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
        like = Likelihood(
            GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)), censoring=censoring
        )
        assert like.marginalisation is Marginalisation.LATENT
        data = Spectrum(
            [1.0, 2.0, 3.0] * u.um, [1.0, 2.0, 0.3] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        )
        with pytest.raises(LikelihoodError, match="orthant probability"):
            like.log_prob(data.with_values([1.0, 2.0, 0.1]), data)

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
        [QuasisepGP(), WindowedSparseGP(), InducingPointGP()],
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
        """Euclidean separation across mixed units is meaningless; say so."""

        class MixedUnits(Spectrum):
            pass

        vis = VisibilitySet(
            [1.0, 2.0] / u.rad, [3.0, 4.0], [1.0 + 0j, 0.5 + 0j], uncertainty=[0.1, 0.1]
        )
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        with pytest.raises(LikelihoodError, match="different units"):
            noise.check_compatible(ComplexGaussianFamily(), vis)

    def test_a_two_axis_point_set_is_allowed_with_one_unit(self) -> None:
        vis = VisibilitySet(
            [1.0, 2.0, 5.0],
            [3.0, 4.0, 1.0],
            [1.0 + 0j, 0.5 + 0j, 0.2 + 0j],
            uncertainty=[0.1, 0.1, 0.1],
        )
        noise = GaussianProcessNoise(Matern32(0.3, 1.0))
        noise.check_compatible(ComplexGaussianFamily(), vis)

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
# Composition-time checking
# ---------------------------------------------------------------------------


class TestComposition:
    def test_the_parameter_namespace_is_flat(self) -> None:
        """One ParameterSet, no nested merge — ``parameters.md`` §12.4."""
        like = Likelihood(
            StudentTFamily(nu=st.loguniform(2.0, 50.0)),
            GaussianProcessNoise(
                Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
                scale=st.loguniform(0.5, 2.0),
            ),
        )
        assert like.parameters.free_names == ("amplitude", "length_scale", "scale", "nu")
        assert like.parameters.free_size == 4

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

    def test_the_repr_states_the_marginalisation(self) -> None:
        like = Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))
        assert repr(like) == (
            "Likelihood(PoissonFamily, GaussianProcessNoise, marginalisation=latent)"
        )

    def test_every_error_here_is_a_contract_error(self) -> None:
        """Catchable by the family, and by the builtin, per ``exceptions.py``."""
        assert issubclass(LikelihoodError, ContractError)
        assert issubclass(LikelihoodError, ValueError)
