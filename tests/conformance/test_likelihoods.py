"""Conformance rows for the likelihood and noise-model contract.

The inventory is ``likelihoods.md`` §16's hand-down, verbatim: "DenseGP against
``scipy.stats.multivariate_normal.logpdf``; the zero-amplitude reduction to the
i.i.d. Gaussian; mask excision equalling deletion of the sample;
DenseGP↔QuasisepGP agreement on Matérn-3/2 once the latter exists (§4.6 names
it); ``L Lᵀ == K`` for the whitening transform; the Tobit censored likelihood
against ``scipy.stats.norm.logcdf``; and the marginalisation declaration for
every family/noise-model pair."

This is the module W1.10 exists for. ``DEVELOPMENT_PLAN.md`` §4.4 makes the
O(N) GP the project's distinguishing feature, and the dense solve is the
correctness anchor every faster strategy must reproduce — so every oracle here
is either ``scipy`` or the kernel written out from its own defining formula.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
import scipy.stats as st
from scipy.stats import multivariate_normal

from ampere.core import (
    Censoring,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Likelihood,
    LikelihoodError,
    LimitKind,
    Marginalisation,
    NoiseModel,
    Spectrum,
    family_named,
    list_families,
)

from .composition import COORDINATE_UNIT, FLUX_UNIT, GP_GRID
from .oracles import analytic_diagonal_gaussian_log_prob, kernel_matrix, matern32_matrix
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    SolverKind,
    Tolerances,
)

AMPLITUDE = 0.4
LENGTH_SCALE = 2.0
SIGMA = 0.1

MATERN32 = CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE)
SQUARED_EXPONENTIAL = CovarianceSpec(KernelFamily.SQUARED_EXPONENTIAL, AMPLITUDE, LENGTH_SCALE)


# ---------------------------------------------------------------------------
# Containers
# ---------------------------------------------------------------------------


def spectra(mask: np.ndarray | None = None, sigma: float = SIGMA) -> tuple[Spectrum, Spectrum]:
    """A fixed ``(predicted, observed)`` pair on :data:`GP_GRID`."""
    grid = np.asarray(GP_GRID, dtype=float)
    rng = np.random.default_rng(20260902)
    truth = 1.0 + 0.3 * np.sin(grid)
    observed = Spectrum(
        grid * COORDINATE_UNIT,
        (truth + rng.normal(0.0, sigma, grid.size)) * FLUX_UNIT,
        uncertainty=np.full(grid.size, sigma) * FLUX_UNIT,
        mask=mask,
    )
    return observed.with_values(truth), observed


def gp_likelihood(
    backend: ConformanceBackend,
    covariance: CovarianceSpec = MATERN32,
    solver: SolverKind = SolverKind.DENSE,
) -> Likelihood:
    """A Gaussian family over *backend*'s kernel and solver."""
    noise = GaussianProcessNoise(backend.kernel(covariance), backend.gp_solver(solver))
    return Likelihood(GaussianFamily(), noise)


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


class TestDenseGPAgainstScipy:
    """The correctness anchor: the dense solve against a general-purpose one."""

    @pytest.mark.parametrize("covariance", [MATERN32, SQUARED_EXPONENTIAL], ids=lambda c: c.family)
    def test_it_agrees_with_multivariate_normal_logpdf(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        predicted, observed = spectra()
        likelihood = gp_likelihood(backend, covariance)

        grid = np.asarray(GP_GRID, dtype=float)
        covariance_matrix = kernel_matrix(
            covariance.family, grid, covariance.amplitude, covariance.length_scale
        )
        total = covariance_matrix + np.diag(np.full(grid.size, SIGMA**2))
        expected = float(
            multivariate_normal.logpdf(
                observed.values - predicted.values, mean=np.zeros(grid.size), cov=total
            )
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )

    def test_the_kernel_matches_its_defining_formula(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        kernel = backend.kernel(MATERN32)
        built = backend.to_numpy(kernel.matrix(grid[:, None], grid[:, None], kernel.resolve(None)))
        expected = matern32_matrix(grid, AMPLITUDE, LENGTH_SCALE)
        assert built == pytest.approx(expected, abs=tolerances.analytic)

    @pytest.mark.parametrize(
        ("covariance", "correlation_at_one_length_scale"),
        [
            # Matérn-3/2: rho(l) = (1 + sqrt(3)) exp(-sqrt(3))
            (MATERN32, (1.0 + math.sqrt(3.0)) * math.exp(-math.sqrt(3.0))),
            # Squared exponential: rho(l) = exp(-1/2)
            (SQUARED_EXPONENTIAL, math.exp(-0.5)),
        ],
        ids=lambda value: getattr(value, "family", "rho"),
    )
    def test_the_kernel_matches_hand_computed_values(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        correlation_at_one_length_scale: float,
        tolerances: Tolerances,
    ) -> None:
        """Two numbers, not a formula — so the whole GP chain has a fixed anchor.

        ``test_the_kernel_matches_its_defining_formula`` compares an
        implementation against a transcription of the same equation, which
        would pass if both carried the same misreading. These two values are
        arithmetic: the marginal variance is ``amplitude**2`` by the
        amplitude-is-a-standard-deviation convention, and the correlation at a
        separation of one length scale is a number that can be checked by
        hand.
        """
        kernel = backend.kernel(covariance)
        points = np.array([[0.0], [covariance.length_scale]])
        built = backend.to_numpy(kernel.matrix(points, points, kernel.resolve(None)))

        assert built[0, 0] == pytest.approx(covariance.amplitude**2, abs=tolerances.analytic)
        assert built[0, 1] / covariance.amplitude**2 == pytest.approx(
            correlation_at_one_length_scale, abs=tolerances.analytic
        )
        assert built[0, 1] == pytest.approx(built[1, 0], abs=tolerances.exact)


class TestZeroAmplitudeReduction:
    """A GP with no amplitude *is* the i.i.d. Gaussian, not merely close to it."""

    def test_it_reduces_exactly_to_the_independent_case(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        predicted, observed = spectra()
        flat = gp_likelihood(backend, CovarianceSpec(KernelFamily.MATERN32, 0.0, LENGTH_SCALE))
        independent = Likelihood(GaussianFamily(), IndependentNoise())
        assert flat.log_prob(predicted, observed) == pytest.approx(
            independent.log_prob(predicted, observed), abs=tolerances.linear_algebra
        )

    def test_it_matches_the_closed_form_written_out(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        predicted, observed = spectra()
        flat = gp_likelihood(backend, CovarianceSpec(KernelFamily.MATERN32, 0.0, LENGTH_SCALE))
        expected = analytic_diagonal_gaussian_log_prob(
            observed.values - predicted.values, np.asarray(observed.uncertainty)
        )
        assert flat.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )


class TestMaskExcision:
    """A masked sample is deleted, not down-weighted."""

    @pytest.mark.parametrize("correlated", [False, True], ids=["independent", "gp"])
    def test_masking_equals_deleting_the_sample(
        self, backend: ConformanceBackend, correlated: bool, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        mask = np.zeros(grid.size, dtype=bool)
        mask[[3, 8]] = True

        masked_predicted, masked_observed = spectra(mask=mask)
        full_predicted, full_observed = spectra()
        keep = ~mask
        sub_observed = Spectrum(
            grid[keep] * COORDINATE_UNIT,
            full_observed.values[keep] * FLUX_UNIT,
            uncertainty=np.asarray(full_observed.uncertainty)[keep] * FLUX_UNIT,
        )
        sub_predicted = sub_observed.with_values(full_predicted.values[keep])

        likelihood = (
            gp_likelihood(backend)
            if correlated
            else Likelihood(GaussianFamily(), IndependentNoise())
        )
        assert likelihood.log_prob(masked_predicted, masked_observed) == pytest.approx(
            likelihood.log_prob(sub_predicted, sub_observed), abs=tolerances.exact
        )

    def test_a_masked_sample_has_no_influence_whatever_its_value(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        mask = np.zeros(grid.size, dtype=bool)
        mask[4] = True
        predicted, observed = spectra(mask=mask)
        likelihood = gp_likelihood(backend)

        baseline = likelihood.log_prob(predicted, observed)
        wrecked = observed.values.copy()
        wrecked[4] = 1e6
        assert likelihood.log_prob(
            predicted, observed.with_values(wrecked * FLUX_UNIT)
        ) == pytest.approx(baseline, abs=tolerances.exact)

    def test_a_fully_masked_dataset_contributes_exactly_nothing(
        self, backend: ConformanceBackend
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        predicted, observed = spectra(mask=np.ones(grid.size, dtype=bool))
        assert gp_likelihood(backend).log_prob(predicted, observed) == 0.0


class TestWhiteningTransform:
    """``f = L z`` with ``L Lᵀ == K`` — the latent path's change of variables."""

    def test_the_recovered_factor_reproduces_the_kernel(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        kernel = backend.kernel(MATERN32)
        solver = backend.gp_solver(SolverKind.DENSE)
        points = grid[:, None]
        values = kernel.resolve(None)

        basis = np.eye(grid.size)
        lower = np.column_stack(
            [
                backend.to_numpy(solver.latent_transform(kernel, points, basis[:, i], values))
                for i in range(grid.size)
            ]
        )
        expected = matern32_matrix(grid, AMPLITUDE, LENGTH_SCALE)
        assert lower @ lower.T == pytest.approx(expected, abs=tolerances.linear_algebra)

    def test_the_factor_is_lower_triangular(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        kernel = backend.kernel(MATERN32)
        solver = backend.gp_solver(SolverKind.DENSE)
        basis = np.eye(grid.size)
        lower = np.column_stack(
            [
                backend.to_numpy(
                    solver.latent_transform(
                        kernel, grid[:, None], basis[:, i], kernel.resolve(None)
                    )
                )
                for i in range(grid.size)
            ]
        )
        assert np.triu(lower, k=1) == pytest.approx(
            np.zeros((grid.size, grid.size)), abs=tolerances.exact
        )


class TestTobitCensoring:
    """A limit is a cumulative probability, and it is ``scipy``'s."""

    @pytest.mark.parametrize(
        ("kind", "tail"),
        [(LimitKind.UPPER_LIMIT, st.norm.logcdf), (LimitKind.LOWER_LIMIT, st.norm.logsf)],
        ids=["upper", "lower"],
    )
    def test_a_limit_is_the_scipy_tail_probability(
        self, kind: LimitKind, tail: object, tolerances: Tolerances
    ) -> None:
        predicted, observed = spectra()
        censored = (1, 7)
        codes = np.zeros(observed.n_samples, dtype=np.int8)
        codes[list(censored)] = int(kind)
        likelihood = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))

        sigma = np.asarray(observed.uncertainty)
        detected = [i for i in range(observed.n_samples) if i not in censored]
        expected = float(
            np.sum(
                st.norm.logpdf(
                    observed.values[detected], predicted.values[detected], sigma[detected]
                )
            )
        )
        for index in censored:
            standardised = (observed.values[index] - predicted.values[index]) / sigma[index]
            expected += float(tail(standardised))  # type: ignore[operator]

        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.analytic
        )

    def test_masking_beats_censoring(self, tolerances: Tolerances) -> None:
        """A masked limit is excised, not scored as a limit (``likelihoods.md`` §9)."""
        grid = np.asarray(GP_GRID, dtype=float)
        mask = np.zeros(grid.size, dtype=bool)
        mask[1] = True
        predicted, observed = spectra(mask=mask)
        codes = np.zeros(grid.size, dtype=np.int8)
        codes[1] = int(LimitKind.UPPER_LIMIT)

        censored = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))
        plain = Likelihood(GaussianFamily(), IndependentNoise())
        assert censored.log_prob(predicted, observed) == pytest.approx(
            plain.log_prob(predicted, observed), abs=tolerances.exact
        )

    def test_a_limit_under_a_correlated_noise_model_is_refused(
        self, backend: ConformanceBackend
    ) -> None:
        """No closed form exists, and the contract says so rather than guessing."""
        predicted, observed = spectra()
        codes = np.zeros(observed.n_samples, dtype=np.int8)
        codes[2] = int(LimitKind.UPPER_LIMIT)
        noise = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        likelihood = Likelihood(GaussianFamily(), noise, censoring=Censoring(codes))
        with pytest.raises(LikelihoodError, match="orthant probability"):
            likelihood.log_prob(predicted, observed)


# ---------------------------------------------------------------------------
# The marginalisation declaration table
# ---------------------------------------------------------------------------

#: The expected marginalisation for every registered family, under an
#: uncorrelated and a correlated noise model. Written out rather than
#: recomputed from ``ANALYTIC_WITH_GP``, so the table is an oracle and not a
#: restatement: an uncorrelated noise model always marginalises analytically,
#: and under a GP only the Gaussian family still does.
EXPECTED_MARGINALISATION: dict[str, tuple[Marginalisation, Marginalisation]] = {
    "gaussian": (Marginalisation.ANALYTIC, Marginalisation.ANALYTIC),
    "student_t": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "cauchy": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "complex_gaussian": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "poisson": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "rice": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "von_mises": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
}


class TestMarginalisationDeclaration:
    """Every (family, noise model) pair declares how it marginalises."""

    def test_the_table_covers_every_registered_family(self) -> None:
        assert set(list_families()) == set(EXPECTED_MARGINALISATION)

    @pytest.mark.parametrize("name", sorted(EXPECTED_MARGINALISATION))
    @pytest.mark.parametrize("correlated", [False, True], ids=["independent", "gp"])
    def test_each_pair_declares_the_expected_marginalisation(
        self, backend: ConformanceBackend, name: str, correlated: bool
    ) -> None:
        family = family_named(name)()
        noise: NoiseModel = (
            GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
            if correlated
            else IndependentNoise()
        )
        expected = EXPECTED_MARGINALISATION[name][int(correlated)]
        assert family.marginalisation_with(noise) is expected

    def test_marginalisation_is_a_property_of_the_combination(
        self, backend: ConformanceBackend
    ) -> None:
        """Neither piece can declare it alone: Poisson is analytic, Poisson+GP is not."""
        poisson = family_named("poisson")()
        assert poisson.marginalisation_with(IndependentNoise()) is Marginalisation.ANALYTIC
        gp = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        assert poisson.marginalisation_with(gp) is Marginalisation.LATENT

    def test_censoring_forces_the_latent_path_under_a_gp(self, backend: ConformanceBackend) -> None:
        codes = np.zeros(len(GP_GRID), dtype=np.int8)
        codes[0] = int(LimitKind.UPPER_LIMIT)
        gp = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        assert GaussianFamily().marginalisation_with(gp, Censoring(codes)) is Marginalisation.LATENT

    def test_the_conservative_declaration_is_relaxed_by_the_data(self) -> None:
        """A masked limit does not force a gradient-based engine on the run."""
        grid = np.asarray(GP_GRID, dtype=float)
        mask = np.zeros(grid.size, dtype=bool)
        mask[0] = True
        _, observed = spectra(mask=mask)
        codes = np.zeros(grid.size, dtype=np.int8)
        codes[0] = int(LimitKind.UPPER_LIMIT)
        likelihood = Likelihood(GaussianFamily(), IndependentNoise(), censoring=Censoring(codes))
        assert likelihood.marginalisation_for(observed) is Marginalisation.ANALYTIC


# ---------------------------------------------------------------------------
# Skeleton rows: the second GP solver
# ---------------------------------------------------------------------------


class TestSolverAgreement:
    """``DenseGP``↔``QuasisepGP``: the row ``DEVELOPMENT_PLAN.md`` §4.6 names.

    A Matérn-3/2 kernel has an exact celerite representation, so the
    quasiseparable state-space recursion computes the *same* marginal
    likelihood the dense Cholesky does, in linear time. That equivalence is
    why ``DenseGP`` exists at all, and it is the row that will catch a
    quasiseparable implementation that is subtly not exact.
    """

    def test_the_quasiseparable_solver_reproduces_the_dense_one(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(
                f"backend {backend.name!r} declares no quasiseparable solver. "
                "Phase 2 must supply one (celerite2 on the numpy and jax sides, "
                "tinygp's QuasisepSolver, GPyTorch or celerite2-torch on the torch "
                "side) and add SolverKind.QUASISEP to its capabilities; this row "
                "then compares it against DenseGP on Matern32 at "
                "tolerances.cross_solver."
            )
        predicted, observed = spectra()
        dense = gp_likelihood(backend, MATERN32, SolverKind.DENSE)
        quasisep = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        assert quasisep.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=tolerances.cross_solver
        )

    def test_an_unimplemented_solver_slot_refuses_rather_than_pretending(
        self, backend: ConformanceBackend
    ) -> None:
        """Until it exists, asking for it must fail loudly and say what to use.

        The other half of the skeleton row: a declared-but-empty strategy slot
        that silently fell back to the dense path would make the agreement row
        above vacuous the day someone forgot to implement it.
        """
        if SolverKind.QUASISEP in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} implements the quasiseparable solver")
        predicted, observed = spectra()
        likelihood = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        with pytest.raises(LikelihoodError, match="no implementation yet"):
            likelihood.log_prob(predicted, observed)

    def test_a_quasiseparable_solver_refuses_a_kernel_that_has_no_such_form(
        self, backend: ConformanceBackend
    ) -> None:
        """The declaration is checked before the implementation is."""
        _, observed = spectra()
        kernel = backend.kernel(SQUARED_EXPONENTIAL)
        assert not kernel.QUASISEPARABLE
        noise = GaussianProcessNoise(kernel, backend.gp_solver(SolverKind.QUASISEP))
        with pytest.raises(LikelihoodError, match="quasiseparable representation"):
            noise.check_compatible(GaussianFamily(), observed)
