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

from ampere.backends.reference import FractionalModelGPNoise, FractionalModelNoise

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
#: and under a GP only the Gaussian family — and, since the 2026-09-03 ruling
#: fixed the circular complex GP as its meaning (``likelihoods.md`` §17 Q6),
#: the complex Gaussian — still does. The complex closed form is Phase 4's to
#: implement; ``TestStagedAnalyticCombination`` holds the refusal meanwhile.
EXPECTED_MARGINALISATION: dict[str, tuple[Marginalisation, Marginalisation]] = {
    "gaussian": (Marginalisation.ANALYTIC, Marginalisation.ANALYTIC),
    "student_t": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "cauchy": (Marginalisation.ANALYTIC, Marginalisation.LATENT),
    "complex_gaussian": (Marginalisation.ANALYTIC, Marginalisation.ANALYTIC),
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


class TestStagedAnalyticCombination:
    """``complex_gaussian`` + GP: declared ``ANALYTIC``, implemented in Phase 4.

    The 2026-09-03 ruling (``likelihoods.md`` §17 Q6) fixed the circular
    (equal-component, zero-pseudo-covariance) complex GP as the combination's
    meaning and landed the declaration at the freeze. Until Phase 4 implements
    the closed form, composition must refuse with the schedule named — the
    same declared-but-staged discipline as the ``QuasisepGP`` slot above.
    Phase 4 replaces the refusal row with agreement rows against the circular
    closed form.
    """

    def test_the_declaration_is_analytic(self, backend: ConformanceBackend) -> None:
        from ampere.core import ComplexGaussianFamily

        gp = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        family = ComplexGaussianFamily()
        assert family.marginalisation_with(gp) is Marginalisation.ANALYTIC

    def test_the_staged_combination_refuses_rather_than_pretending(
        self, backend: ConformanceBackend
    ) -> None:
        from ampere.core import ComplexGaussianFamily

        gp = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        with pytest.raises(LikelihoodError, match="Phase 4"):
            Likelihood(ComplexGaussianFamily(), gp)


# ---------------------------------------------------------------------------
# The second GP solver (live for the reference path since W2.3)
# ---------------------------------------------------------------------------


class TestSolverAgreement:
    """``DenseGP``↔``QuasisepGP``: the row ``DEVELOPMENT_PLAN.md`` §4.6 names.

    A Matérn-3/2 kernel has an exact celerite representation, so the
    quasiseparable state-space recursion computes the *same* marginal
    likelihood the dense Cholesky does, in linear time. That equivalence is
    why ``DenseGP`` exists at all, and it is the row that will catch a
    quasiseparable implementation that is subtly not exact — including one
    that reached for celerite2's own ``Matern32Term``, whose ``eps``
    approximation misses ``tolerances.cross_solver`` by three orders of
    magnitude.
    """

    def test_the_quasiseparable_solver_reproduces_the_dense_one(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(
                f"backend {backend.name!r} declares no quasiseparable solver. "
                "Supply one (celerite2 on the numpy and jax sides, tinygp's "
                "QuasisepSolver, GPyTorch or celerite2-torch on the torch "
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

    def test_the_two_solvers_agree_on_the_conditioned_gp(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The §4.8 localisation diagnostic must not depend on the strategy."""
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        predicted, observed = spectra()
        dense = gp_likelihood(backend, MATERN32, SolverKind.DENSE)
        quasisep = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        at = np.linspace(min(GP_GRID) - 0.5, max(GP_GRID) + 0.5, 29)
        for kwargs in ({}, {"at": at}):
            expected = dense.conditional(predicted, observed, **kwargs)
            got = quasisep.conditional(predicted, observed, **kwargs)
            assert backend.to_numpy(got.mean) == pytest.approx(
                backend.to_numpy(expected.mean), abs=tolerances.cross_solver
            )
            assert backend.to_numpy(got.variance) == pytest.approx(
                backend.to_numpy(expected.variance), abs=tolerances.cross_solver
            )

    def test_the_declared_solver_is_not_the_dense_one_in_disguise(
        self, backend: ConformanceBackend
    ) -> None:
        """The agreement rows above are only worth running if this holds.

        A backend that declared ``QUASISEP`` and then handed back its dense
        strategy would satisfy every agreement row perfectly and prove
        nothing. So the strategy behind the declaration must be a *different*
        one, must say it is implemented and exact, and must compute.

        Before W2.3 this row's other half ran instead: with no quasiseparable
        solver anywhere, it asserted that asking for the empty slot refused.
        That assertion moved to ``tests/core`` when the slot was filled,
        because it is a property of ``ampere.core``'s declared-slot discipline
        rather than of a backend — and because ``gp_solver`` is contractually
        called only for a kind the backend declares (``protocol.py``), so a
        row must not ask a torch or jax fixture for a solver it has said it
        does not have.
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        predicted, observed = spectra()
        solver = backend.gp_solver(SolverKind.QUASISEP)
        assert solver.IMPLEMENTED
        assert solver.EXACT
        assert solver.NAME != backend.gp_solver(SolverKind.DENSE).NAME
        assert type(solver) is not type(backend.gp_solver(SolverKind.DENSE))
        likelihood = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        assert math.isfinite(likelihood.log_prob(predicted, observed))

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

    def test_the_leave_one_out_terms_are_a_dense_only_decomposition(
        self, backend: ConformanceBackend
    ) -> None:
        """W2.3 deferred ``conditional_loo``'s O(N) recursion; it must say so.

        The refusal is the declared-but-staged discipline again: the O(N)
        route to the diagonal of ``(K + diag(sigma**2))**-1`` is recorded as
        deferred in ``DEVELOPMENT_PLAN.md`` §2, and until it exists a
        quasiseparable strategy refuses rather than quietly costing O(N**2).
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        predicted, observed = spectra()
        dense = gp_likelihood(backend, MATERN32, SolverKind.DENSE)
        assert dense.pointwise_log_prob(predicted, observed).size == len(GP_GRID)
        quasisep = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        with pytest.raises(LikelihoodError, match="conditional_loo"):
            quasisep.pointwise_log_prob(predicted, observed)


# ---------------------------------------------------------------------------
# Prediction-aware noise: the X-1 rows (ruled 2026-09-03)
# ---------------------------------------------------------------------------


# W2.1: these rows ran against in-repo doubles until the shipped classes
# landed. They now exercise ``ampere.backends.reference`` itself, which is what
# ``likelihoods.md`` §5 and the README's debt entry always intended -- the
# contract fixes the name and the semantics, and the implementation ships with
# the reference backend. The doubles are gone: a double that shadows a shipped
# class only tests itself.
FRACTION = 0.5


class TestPredictionAwareNoise:
    """The three X-1 conformance rows (``awkward_instrument.md`` §6 point 9).

    ``NoiseModel.sigma``/``noise_params`` receive the retained predicted
    values as a keyword-only ``predicted`` at every call site, so a noise
    whose magnitude depends on the model is a noise model, not a family.
    """

    def test_fractional_sigma_equals_the_manual_quadrature(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Row (a): sigma_eff from a fractional noise, at fixed theta."""
        predicted, observed = spectra()
        likelihood = Likelihood(GaussianFamily(), FractionalModelNoise(FRACTION))
        mean = np.asarray(predicted.values, dtype=float)
        sigma_eff = np.sqrt(SIGMA**2 + (FRACTION * mean) ** 2)
        expected = analytic_diagonal_gaussian_log_prob(
            np.asarray(observed.values, dtype=float) - mean, sigma_eff
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.analytic
        )

    def test_draw_variance_grows_with_the_prediction(self, backend: ConformanceBackend) -> None:
        """Row (b): ``simulate(observe=True)`` draws with variance ``sigma**2 + (f*mu)**2``."""
        from ampere.core import Dataset

        from .composition import DatasetSpec, ProblemSpec, observed_container

        spec = ProblemSpec()
        dataset_spec = DatasetSpec()
        observed = observed_container(spec, dataset_spec)
        dataset = Dataset(
            observed,
            likelihood=Likelihood(GaussianFamily(), FractionalModelNoise(FRACTION)),
            label="d",
        )
        mean = 1.0 + 0.1 * np.arange(observed.n_samples, dtype=float)
        predicted = observed.with_values(mean)
        rng = np.random.default_rng(20260903)
        draws = np.array(
            [
                np.asarray(dataset.draw_observation(predicted, None, rng).values, dtype=float)
                for _ in range(4000)
            ]
        )
        expected = dataset_spec.uncertainty**2 + (FRACTION * mean) ** 2
        # The sample variance of a Gaussian has standard error var * sqrt(2/(n-1)).
        standard_error = expected * math.sqrt(2.0 / (draws.shape[0] - 1))
        limit = backend.capabilities.tolerances.monte_carlo_sigmas * standard_error
        assert np.all(np.abs(np.var(draws, axis=0, ddof=1) - expected) < limit)
        # And the growth itself: the fractional term dominates sigma_data here,
        # so the last sample's draw variance must exceed the first's.
        assert np.var(draws[:, -1], ddof=1) > 2.0 * np.var(draws[:, 0], ddof=1)

    def test_the_gp_composition_agrees_with_a_manual_diagonal(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Row (c): K + diag(sigma_data**2 + (f*mu)**2), against scipy's own solve."""
        predicted, observed = spectra()
        noise = FractionalModelGPNoise(
            backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE), f=FRACTION
        )
        likelihood = Likelihood(GaussianFamily(), noise)
        grid = np.asarray(GP_GRID, dtype=float)
        mean = np.asarray(predicted.values, dtype=float)
        covariance = kernel_matrix(
            MATERN32.family, grid, MATERN32.amplitude, MATERN32.length_scale
        ) + np.diag(SIGMA**2 + (FRACTION * mean) ** 2)
        expected = float(
            multivariate_normal.logpdf(
                np.asarray(observed.values, dtype=float) - mean,
                mean=np.zeros(grid.size),
                cov=covariance,
            )
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )
