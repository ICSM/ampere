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
from collections.abc import Iterator, Sequence
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st
from scipy.stats import multivariate_normal

from ampere.core import (
    AxisSpec,
    Censoring,
    ComplexGaussianFamily,
    FunctionSamples,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Layout,
    Likelihood,
    LikelihoodError,
    LimitKind,
    Marginalisation,
    NoiseModel,
    Order,
    Spectrum,
    VisibilitySet,
    family_named,
    list_families,
    register_quasiseparable_term,
    term_provenance_entries,
)
from ampere.core import Matern12 as CoreMatern12
from ampere.core.kernels import _forget_quasiseparable_term, matern12_representation

from ampere.backends.reference import FractionalModelGPNoise, FractionalModelNoise

from .backends._kernels import USER_FAMILY, user_kernel_type
from .composition import COORDINATE_UNIT, FLUX_UNIT, GP_GRID
from .oracles import (
    analytic_diagonal_gaussian_log_prob,
    covariance_matrix,
    kernel_matrix,
    matern32_matrix,
)
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    SolverKind,
    Tolerances,
    approximation_envelope,
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
    """A combination declared ``ANALYTIC`` whose closed form is not written.

    The 2026-09-03 ruling (``likelihoods.md`` §17 Q6) fixed the circular
    (equal-component, zero-pseudo-covariance) complex GP as
    ``complex_gaussian`` + GP's meaning and landed the *declaration* at the
    freeze, with composition refused and the schedule named.

    **W4.2 wrote it**, so this class keeps the discipline and loses its
    instance: the rows that say the closed form is right are
    :class:`TestTheCircularComplexGP`, and what is held here is that a family
    declaring ``ANALYTIC_WITH_GP`` *without* ``GP_ANALYTIC_IMPLEMENTED`` is
    still refused rather than quietly computed some other way. Held against a
    family declared in the row, because no shipped family is in that position
    any more — which is the point of keeping the row.
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

        class StagedFamily(ComplexGaussianFamily):
            NAME = "staged_complex"
            GP_ANALYTIC_IMPLEMENTED = False

        gp = GaussianProcessNoise(backend.kernel(MATERN32), backend.gp_solver(SolverKind.DENSE))
        with pytest.raises(LikelihoodError, match="GP_ANALYTIC_IMPLEMENTED is False"):
            Likelihood(StagedFamily(), gp)


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
        """The declaration is checked before the implementation is.

        Guarded on the declared solver set like every other row in this class.
        ``protocol.py`` says ``gp_solver`` is "only called for a kind present
        in ``capabilities.solvers``; a backend may raise for anything else",
        and the row above says the same in prose — "a row must not ask a torch
        or jax fixture for a solver it has said it does not have". This row was
        asking anyway, which was invisible while both in-repo fixtures declared
        ``QUASISEP`` and became a spurious failure the moment W2.4's torch
        fixture registered without one. Restoring the guard is a correction to
        this suite, not a relaxation of it: the claim being checked is
        ``GPSolver.check_compatible``'s declarative refusal, which
        ``tests/core`` also exercises on ``ampere.core``'s own solver, so
        nothing goes unchecked for a backend that skips here.
        (Found independently by W2.4 and W2.5: both tracks' fixtures skip
        here until slice 2 supplies a quasiseparable solver.)
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        _, observed = spectra()
        kernel = backend.kernel(SQUARED_EXPONENTIAL)
        assert not kernel.QUASISEPARABLE
        noise = GaussianProcessNoise(kernel, backend.gp_solver(SolverKind.QUASISEP))
        with pytest.raises(LikelihoodError, match="quasiseparable representation"):
            noise.check_compatible(GaussianFamily(), observed)

    def test_the_leave_one_out_terms_either_agree_exactly_or_refuse_by_name(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``conditional_loo`` on the quasiseparable path: agreement, or a refusal.

        This row used to assert only the refusal, because W2.3 **deferred**
        the O(N) recursion for the diagonal of ``(K + diag(sigma**2))**-1``
        (``DEVELOPMENT_PLAN.md`` §2, 2026-09-05): celerite2's public numpy
        interface exposes no route to that diagonal, and shipping an O(N**2)
        fallback under an O(N) name was rejected outright. W2.4 slice 2
        changed the circumstances rather than the ruling — the torch backend
        calls celerite2's compiled kernels directly, so it already knows the
        factorisation convention the recursion needs, and it supplies the
        terms (decision-log row, 2026-09-07).

        So the row states the *whole* contract, which is what the README's
        debt entry always said it would become: a quasiseparable strategy
        either refuses by name, or computes the same decomposition the dense
        solver does. The one thing it may not do is return something else —
        a wrong ``A_ii`` would still be finite, still be per-sample, and
        still look exactly like a leave-one-out term.
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        predicted, observed = spectra()
        dense = gp_likelihood(backend, MATERN32, SolverKind.DENSE)
        expected = dense.pointwise_log_prob(predicted, observed)
        assert expected.size == len(GP_GRID)
        quasisep = gp_likelihood(backend, MATERN32, SolverKind.QUASISEP)
        try:
            got = quasisep.pointwise_log_prob(predicted, observed)
        except LikelihoodError as refusal:
            assert "conditional_loo" in str(refusal)
            return
        assert got.shape == expected.shape
        assert np.allclose(got, expected, rtol=tolerances.cross_solver, atol=0.0)


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


# ---------------------------------------------------------------------------
# W4.5: the kernel algebra, the axis selector and the term registry
# ---------------------------------------------------------------------------


SHO_SPEC = CovarianceSpec(KernelFamily.SHO, AMPLITUDE, period=1.5, quality=3.0)
LEAF_SPECS = [
    CovarianceSpec(KernelFamily.MATERN12, AMPLITUDE, LENGTH_SCALE),
    CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE),
    CovarianceSpec(KernelFamily.MATERN52, AMPLITUDE, LENGTH_SCALE),
    SHO_SPEC,
    CovarianceSpec(KernelFamily.ROTATION, AMPLITUDE, period=1.5, quality=3.0),
]
SUM_SPEC = CovarianceSpec(
    KernelFamily.SUM,
    terms=(
        CovarianceSpec(KernelFamily.MATERN32, 0.3, 2.0),
        CovarianceSpec(KernelFamily.SHO, 0.2, period=1.5, quality=4.0),
    ),
)
MIXTURE_SPEC = CovarianceSpec(
    KernelFamily.SPECTRAL_MIXTURE,
    terms=(
        CovarianceSpec(KernelFamily.SHO, 0.3, period=1.5, quality=3.0),
        CovarianceSpec(KernelFamily.SHO, 0.2, period=0.4, quality=8.0),
    ),
)
USER_SPEC = CovarianceSpec(KernelFamily.USER, AMPLITUDE, LENGTH_SCALE)

#: The chromatic covariance of ``phase4_placement_memo.md`` §3.6: smooth in
#: spatial frequency, sharp in wavelength, on disjoint axis subsets of one
#: three-axis container. The spatial term's amplitude is 1 because a product's
#: marginal variance is the product of its terms' — two free amplitudes would
#: over-parameterise it by one degree of freedom.
CHROMATIC_SPEC = CovarianceSpec(
    KernelFamily.PRODUCT,
    terms=(
        CovarianceSpec(KernelFamily.MATERN32, 1.0, 2.0, axes=("u", "v")),
        CovarianceSpec(
            KernelFamily.MATERN32, 0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um
        ),
    ),
)


class DispersedPoints(FunctionSamples):
    """A three-axis point set: (u, v) dimensionless and a wavelength in micron.

    Local to this module deliberately. ``VisibilitySet`` is being amended to
    three axes by W4.1 *in parallel*, and a conformance row that waited for it
    would be a dependency W4.5 does not have. What it has to be is the smallest
    container on which the axis selector is necessary rather than merely
    available: an isotropic kernel over all three axes is refused by the
    single-unit rule, and a ``Product`` of two selections is not.
    """

    AXES = (
        AxisSpec("u", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("v", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("spectral_axis", physical_types=("length",), order=Order.ANY),
    )
    LAYOUT = Layout.POINTS

    def __init__(self, u_coord: Any, v_coord: Any, spectral_axis: Any, values: Any, **kw: Any):
        super().__init__(
            {"u": u_coord, "v": v_coord, "spectral_axis": spectral_axis},
            values,
            **kw,
        )


def dispersed_pair() -> tuple[DispersedPoints, DispersedPoints]:
    """A fixed ``(predicted, observed)`` pair on a three-axis container."""
    rng = np.random.default_rng(20260911)
    n = 12
    u_points = rng.uniform(-3.0, 3.0, n)
    v_points = rng.uniform(-3.0, 3.0, n)
    wavelength = np.sort(rng.uniform(2.0, 2.4, n))
    truth = 1.0 + 0.1 * u_points
    observed = DispersedPoints(
        u_points * u.dimensionless_unscaled,
        v_points * u.dimensionless_unscaled,
        wavelength * u.um,
        (truth + rng.normal(0.0, SIGMA, n)) * FLUX_UNIT,
        uncertainty=np.full(n, SIGMA) * FLUX_UNIT,
    )
    return observed.with_values(truth), observed


def dispersed_points(container: DispersedPoints) -> np.ndarray:
    """The ``(n, 3)`` coordinate block, in the container's own axis order."""
    return np.column_stack([np.asarray(axis.values, dtype=float) for axis in container.axes])


class TestNewKernelFamilies:
    """W4.5's five new terms, each against its own defining formula.

    ``tolerances.cross_solver`` throughout, which is the item's acceptance
    tolerance and three orders tighter than any approximation reaches — the
    point being that these representations are exact rather than good.
    """

    @pytest.mark.parametrize("covariance", LEAF_SPECS, ids=lambda c: str(c.family))
    def test_the_dense_solve_agrees_with_multivariate_normal_logpdf(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        predicted, observed = spectra()
        likelihood = gp_likelihood(backend, covariance)
        grid = np.asarray(GP_GRID, dtype=float)
        total = covariance_matrix(covariance, grid) + np.diag(np.full(grid.size, SIGMA**2))
        expected = float(
            multivariate_normal.logpdf(
                observed.values - predicted.values, mean=np.zeros(grid.size), cov=total
            )
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )

    @pytest.mark.parametrize("covariance", LEAF_SPECS, ids=lambda c: str(c.family))
    def test_the_kernel_matches_its_defining_formula(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        kernel = backend.kernel(covariance)
        built = backend.to_numpy(kernel.matrix(grid[:, None], grid[:, None], kernel.resolve(None)))
        assert built == pytest.approx(covariance_matrix(covariance, grid), abs=tolerances.analytic)

    @pytest.mark.parametrize("covariance", LEAF_SPECS, ids=lambda c: str(c.family))
    def test_the_amplitude_is_a_marginal_standard_deviation(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        """``k(0) == amplitude**2`` in every family: one convention, no exceptions."""
        kernel = backend.kernel(covariance)
        diagonal = backend.to_numpy(kernel.diagonal(np.zeros(3), kernel.resolve(None)))
        assert diagonal == pytest.approx(covariance.amplitude**2, abs=tolerances.analytic)

    @pytest.mark.parametrize(
        "covariance", [*LEAF_SPECS, SUM_SPEC, MIXTURE_SPEC], ids=lambda c: str(c.family)
    )
    def test_the_quasiseparable_solve_agrees_with_the_dense_one(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        """The item's headline row: every term exact on the O(N) path.

        The two solvers compute the same number by genuinely different
        recursions — a dense Cholesky against a semiseparable factorisation —
        so agreement at ``cross_solver`` is evidence that the representation is
        the kernel, not merely near it.
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        predicted, observed = spectra()
        dense = gp_likelihood(backend, covariance, SolverKind.DENSE)
        quasisep = gp_likelihood(backend, covariance, SolverKind.QUASISEP)
        assert quasisep.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=tolerances.cross_solver
        )


class TestKernelAlgebra:
    """``Sum`` against the sum of matrices, ``Product`` against their product."""

    def test_a_sum_is_the_sum_of_its_terms_matrices(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        grid = np.asarray(GP_GRID, dtype=float)
        composed = backend.kernel(SUM_SPEC)
        built = backend.to_numpy(
            composed.matrix(grid[:, None], grid[:, None], composed.resolve(None))
        )
        blocks = []
        for term in SUM_SPEC.terms:
            kernel = backend.kernel(term)
            blocks.append(
                backend.to_numpy(kernel.matrix(grid[:, None], grid[:, None], kernel.resolve(None)))
            )
        assert built == pytest.approx(blocks[0] + blocks[1], abs=tolerances.exact)
        # And against the oracle, so the row is not two implementations agreeing.
        assert built == pytest.approx(covariance_matrix(SUM_SPEC, grid), abs=tolerances.analytic)

    def test_a_sums_hyperparameters_are_namespaced_by_term(
        self, backend: ConformanceBackend
    ) -> None:
        composed = backend.kernel(SUM_SPEC)
        assert composed.parameters.names == (
            "term0.amplitude",
            "term0.length_scale",
            "term1.amplitude",
            "term1.period",
            "term1.quality",
        )

    def test_the_spec_is_a_tree_in_declaration_order(self, backend: ConformanceBackend) -> None:
        """The declaration a run's hash is taken over, across every backend."""
        spec = backend.kernel(SUM_SPEC).spec()
        assert spec.family == "sum"
        assert [label for label, _ in spec.terms] == ["term0", "term1"]
        assert [child.family for _, child in spec.terms] == ["matern32", "sho"]

    def test_a_product_is_refused_on_the_quasiseparable_path_word_for_word(
        self, backend: ConformanceBackend
    ) -> None:
        """Not "this kernel is not quasiseparable" — the structural reason, by name."""
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        _, observed = spectra()
        product = CovarianceSpec(
            KernelFamily.PRODUCT,
            terms=(
                CovarianceSpec(KernelFamily.MATERN32, 1.0, 2.0),
                CovarianceSpec(KernelFamily.MATERN12, 0.4, 0.5),
            ),
        )
        noise = GaussianProcessNoise(
            backend.kernel(product), backend.gp_solver(SolverKind.QUASISEP)
        )
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), observed)
        assert str(excinfo.value) == (
            "QuasisepGP cannot lower a Product: a product of quasiseparable kernels is not "
            "quasiseparable. Where the factors act on different axes — which is what a Product "
            "is for — the result is not a function of one ordered coordinate at all, and where "
            "they act on the same one the semiseparable rank multiplies and is not recoverable "
            "from the factors' own representations. Use DenseGP, or replace the Product with a "
            "Sum, which is quasiseparable exactly when every term is."
        )


class TestAxisSelector:
    """The selector's unit rule, and the product across disjoint axis subsets."""

    def test_a_bare_kernel_still_refuses_a_mixed_unit_container_word_for_word(
        self, backend: ConformanceBackend
    ) -> None:
        """The pre-W4.5 rule, unchanged: the selector adds a case, it removes none."""
        _, observed = dispersed_pair()
        noise = GaussianProcessNoise(
            backend.kernel(CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE)),
            backend.gp_solver(SolverKind.DENSE),
        )
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), observed)
        # Amended at W4.2: the advice is now the selector rather than "use one
        # axis", because on the modality this refusal actually fires for --
        # a VisibilitySet, whose axes are (u, v, spectral_axis) -- "use one
        # axis" is not the fix and `axes=("u", "v")` is.
        assert str(excinfo.value) == (
            "DenseGP measures separation as a Euclidean distance across a DispersedPoints's "
            "coordinate axes, but they carry different units ['', 'um']. A single isotropic "
            "length-scale is meaningless across mixed units; name the axes this kernel acts on "
            "with axes=(...), as in Matern32(axes=('u',)), and compose kernels on different axes "
            "with Product."
        )

    def test_a_mixed_unit_selection_is_refused_word_for_word(
        self, backend: ConformanceBackend
    ) -> None:
        """The new rule: the single-unit check follows the kernel's own selection."""
        _, observed = dispersed_pair()
        spec = CovarianceSpec(
            KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE, axes=("u", "spectral_axis")
        )
        noise = GaussianProcessNoise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE))
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), observed)
        assert str(excinfo.value) == (
            "Matern32 selects the DispersedPoints's axes ('u', 'spectral_axis'), which carry "
            "different units ['', 'um']. A single isotropic length-scale is meaningless across "
            "mixed units; select a subset whose axes share one unit, and compose the rest with "
            "Product."
        )

    def test_a_selection_reaches_the_spec_and_so_the_hash(
        self, backend: ConformanceBackend
    ) -> None:
        plain = backend.kernel(CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE))
        selected = backend.kernel(
            CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE, axes=("u", "v"))
        )
        assert plain.spec().axes is None
        assert "axes" not in plain.spec().to_dict()
        assert selected.spec().to_dict()["axes"] == ["u", "v"]

    def test_a_product_of_two_selections_is_the_elementwise_product(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The acceptance row: the chromatic case on a three-axis container.

        The oracle is the elementwise product of the two factors' matrices,
        each built from its own defining formula over its own axes — so this
        compares an implementation with an equation, not with itself.
        """
        _, observed = dispersed_pair()
        points = dispersed_points(observed)
        noise = GaussianProcessNoise(
            backend.kernel(CHROMATIC_SPEC), backend.gp_solver(SolverKind.DENSE)
        )
        noise.check_compatible(GaussianFamily(), observed)
        bound = noise.kernel_for(observed)
        built = backend.to_numpy(bound.matrix(points, points, bound.resolve(None)))
        assert built == pytest.approx(
            covariance_matrix(CHROMATIC_SPEC, points), abs=tolerances.analytic
        )

    def test_the_product_reaches_the_log_likelihood(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Composable is not enough: a number comes out, and scipy agrees with it."""
        predicted, observed = dispersed_pair()
        likelihood = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(
                backend.kernel(CHROMATIC_SPEC), backend.gp_solver(SolverKind.DENSE)
            ),
        )
        points = dispersed_points(observed)
        total = covariance_matrix(CHROMATIC_SPEC, points) + np.diag(
            np.full(points.shape[0], SIGMA**2)
        )
        expected = float(
            multivariate_normal.logpdf(
                np.asarray(observed.values, dtype=float)
                - np.asarray(predicted.values, dtype=float),
                mean=np.zeros(points.shape[0]),
                cov=total,
            )
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )

    def test_a_selected_axis_reaches_the_quasiseparable_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """One ordered axis out of three: the O(N) solve on the spectral column."""
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        predicted, observed = dispersed_pair()
        spec = CovarianceSpec(
            KernelFamily.MATERN32,
            0.3,
            0.05,
            axes=("spectral_axis",),
            length_scale_unit=u.um,
        )
        dense = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE)),
        )
        quasisep = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)),
        )
        assert quasisep.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=tolerances.cross_solver
        )


@pytest.fixture
def registered_user_term() -> Iterator[None]:
    """One registration of the battery's user term, torn down after the row.

    The whole point of the row it serves: ``register_quasiseparable_term`` is
    keyed on the kernel *family*, and every backend's user kernel declares the
    same family, so **one** call reaches all three O(N) paths.
    """
    register_quasiseparable_term(
        user_kernel_type(CoreMatern12), matern12_representation, override=True
    )
    try:
        yield
    finally:
        _forget_quasiseparable_term(USER_FAMILY)


class TestUserRegisteredTerm:
    """A kernel declared outside ampere, made fast by one registration."""

    def test_it_is_refused_until_its_representation_is_registered(
        self, backend: ConformanceBackend
    ) -> None:
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        _, observed = spectra()
        noise = GaussianProcessNoise(
            backend.kernel(USER_SPEC), backend.gp_solver(SolverKind.QUASISEP)
        )
        with pytest.raises(LikelihoodError, match=USER_FAMILY):
            noise.check_compatible(GaussianFamily(), observed)

    def test_the_dense_path_never_needed_a_registration(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The half that already worked: the ABC is public and ``matrix`` is all DenseGP needs."""
        predicted, observed = spectra()
        likelihood = gp_likelihood(backend, USER_SPEC)
        grid = np.asarray(GP_GRID, dtype=float)
        total = covariance_matrix(USER_SPEC, grid) + np.diag(np.full(grid.size, SIGMA**2))
        expected = float(
            multivariate_normal.logpdf(
                observed.values - predicted.values, mean=np.zeros(grid.size), cov=total
            )
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )

    def test_one_registration_carries_it_onto_this_backends_o_n_path(
        self,
        backend: ConformanceBackend,
        registered_user_term: None,
        tolerances: Tolerances,
    ) -> None:
        """The acceptance row, once per registered backend column.

        The registration in the fixture names ``ampere.core``'s ``Matern12``
        subclass; this backend's user kernel is a subclass of **its own**
        ``Matern12`` under the same family name, and that is the only thing
        they share. If the registry had been keyed on ``(family, backend)``,
        this row would need three registrations and would be a much weaker
        claim.
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        predicted, observed = spectra()
        dense = gp_likelihood(backend, USER_SPEC, SolverKind.DENSE)
        quasisep = gp_likelihood(backend, USER_SPEC, SolverKind.QUASISEP)
        assert quasisep.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=tolerances.cross_solver
        )

    def test_the_user_row_is_stamped_in_provenance(self, registered_user_term: None) -> None:
        """Ampere's own rows are not news; a user's are (``lowering.md``'s rule)."""
        entries = [entry for entry in term_provenance_entries() if entry["family"] == USER_FAMILY]
        assert len(entries) == 1
        assert entries[0]["builtin"] is False
        assert entries[0]["kind"] == "quasiseparable_term"


# ---------------------------------------------------------------------------
# The circular complex GP (W4.2)
# ---------------------------------------------------------------------------

#: The (u, v) kernel of the visibility rows: an isotropic Matérn-3/2 on the two
#: dimensionless axes, which is the selection ``likelihoods.md`` §4 calls the
#: normal case for this modality — a bare kernel is refused by the per-leaf unit
#: rule, because a ``VisibilitySet``'s third axis is a wavelength in micron.
UV_SPEC = CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE, axes=("u", "v"))

#: The chromatic product of ``phase4_placement_memo.md`` §3.6, on the real kind:
#: smooth in spatial frequency on the scale of the missing patch, sharp in
#: wavelength on the scale of the band. Both amplitudes are fixed numbers, as
#: §6 requires — a product's marginal variance is the product of its terms', so
#: two fitted amplitudes would over-parameterise it by one degree of freedom.
VISIBILITY_CHROMATIC_SPEC = CovarianceSpec(
    KernelFamily.PRODUCT,
    terms=(
        CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE, axes=("u", "v")),
        CovarianceSpec(
            KernelFamily.MATERN32, 1.0, 0.08, axes=("spectral_axis",), length_scale_unit=u.um
        ),
    ),
)

#: How many baselines the visibility rows carry. Small on purpose: the oracle
#: factorises a dense ``2N`` covariance, and the claim being tested is about the
#: arithmetic rather than about scaling.
VISIBILITY_SIZE = 9


def visibility_pair() -> tuple[VisibilitySet, VisibilitySet]:
    """A fixed ``(predicted, observed)`` pair of **complex** visibilities.

    Three axes in mixed units, and the wavelength genuinely varies, which is
    what the chromatic product row needs: a constant spectral column would make
    the spectral factor the constant matrix ``amplitude**2`` and the row would
    prove nothing about the selection.

    The coordinates are order-unity rather than the 10⁷ wavelengths a real
    baseline is measured in, so that one ``length_scale`` serves every row here;
    the modality at its own scale is ``tests/m2/test_visibility_calibration.py``.
    """
    rng = np.random.default_rng(20260913)
    n = VISIBILITY_SIZE
    u_points = rng.uniform(-3.0, 3.0, n)
    v_points = rng.uniform(-3.0, 3.0, n)
    wavelength = np.sort(rng.uniform(2.0, 2.4, n))
    truth = (1.0 + 0.1 * u_points) * np.exp(0.4j * v_points)
    noise = rng.normal(0.0, SIGMA, n) + 1j * rng.normal(0.0, SIGMA, n)
    observed = VisibilitySet(
        u_points * u.dimensionless_unscaled,
        v_points * u.dimensionless_unscaled,
        wavelength * u.um,
        (truth + noise) * FLUX_UNIT,
        uncertainty=np.full(n, SIGMA) * FLUX_UNIT,
    )
    return observed.with_values(truth), observed


def visibility_points(container: VisibilitySet) -> np.ndarray:
    """The ``(n, 3)`` coordinate block, in the container's own axis order."""
    return np.column_stack([np.asarray(axis.values, dtype=float) for axis in container.axes])


def circular_gp_likelihood(
    backend: ConformanceBackend,
    covariance: CovarianceSpec,
    solver: SolverKind = SolverKind.DENSE,
) -> Likelihood:
    """``complex_gaussian`` over *backend*'s kernel and solver."""
    noise = GaussianProcessNoise(backend.kernel(covariance), backend.gp_solver(solver))
    return Likelihood(ComplexGaussianFamily(), noise)


def real_2n_log_prob(covariance: CovarianceSpec, points: np.ndarray, residual: np.ndarray) -> float:
    """The oracle: the dense real ``2N`` Gaussian, materialised and factorised.

    Deliberately the expensive formulation ``ampere.core`` refuses to use. The
    circular model says the real covariance of ``(Re r, Im r)`` is ``diag(S, S)``
    with a zero off-diagonal block; this builds that ``2N`` by ``2N`` matrix and
    hands it to ``scipy``, so what the rows compare is ampere's exploitation of
    the structure against the structure written out. Nothing here calls a
    solver, a family or a noise model.
    """
    block = covariance_matrix(covariance, points) + np.diag(np.full(points.shape[0], SIGMA**2))
    zeros = np.zeros_like(block)
    total = np.block([[block, zeros], [zeros, block]])
    stacked = np.concatenate([np.real(residual), np.imag(residual)])
    return float(multivariate_normal.logpdf(stacked, mean=np.zeros(stacked.size), cov=total))


class TestTheCircularComplexGP:
    """W4.2: ``complex_gaussian`` + ``GaussianProcessNoise``, the closed form.

    The declaration was made at the freeze (``likelihoods.md`` §4, §17 Q6) and
    the combination refused until Phase 4. These are the rows that say the
    implementation computes what the declaration promised, and every oracle is
    the **dense real 2N Gaussian** rather than another ampere path: the whole
    content of the item is that the block structure may be exploited, so the
    thing to compare against is the structure materialised.
    """

    def test_the_pair_is_declared_and_implemented(self, backend: ConformanceBackend) -> None:
        """The declaration, on every column, since a backend could shadow it."""
        noise = GaussianProcessNoise(backend.kernel(UV_SPEC), backend.gp_solver(SolverKind.DENSE))
        family = ComplexGaussianFamily()
        assert family.ANALYTIC_WITH_GP
        assert family.GP_ANALYTIC_IMPLEMENTED
        assert Likelihood(family, noise).marginalisation is Marginalisation.ANALYTIC

    def test_the_dense_solver_declares_that_it_takes_stacked_residuals(
        self, backend: ConformanceBackend
    ) -> None:
        """``STACKED_RESIDUALS`` is what makes the circular solve expressible."""
        assert backend.gp_solver(SolverKind.DENSE).STACKED_RESIDUALS is True
        if SolverKind.QUASISEP in backend.capabilities.solvers:
            assert backend.gp_solver(SolverKind.QUASISEP).STACKED_RESIDUALS is False

    @pytest.mark.parametrize(
        "covariance", [UV_SPEC, VISIBILITY_CHROMATIC_SPEC], ids=["uv", "chromatic_product"]
    )
    def test_the_closed_form_is_the_dense_real_2n_gaussians(
        self,
        backend: ConformanceBackend,
        covariance: CovarianceSpec,
        tolerances: Tolerances,
    ) -> None:
        """The item's headline row, with the (u, v) kernel and with the product.

        ``tolerances.cross_solver`` because the two computations are genuinely
        different recursions: one Cholesky of an ``N`` by ``N`` matrix with two
        right-hand sides against one Cholesky of a ``2N`` by ``2N`` matrix whose
        off-diagonal block is zero.
        """
        predicted, observed = visibility_pair()
        likelihood = circular_gp_likelihood(backend, covariance)
        likelihood.check_alignment(predicted, observed)
        residual = np.asarray(observed.values) - np.asarray(predicted.values)
        expected = real_2n_log_prob(covariance, visibility_points(observed), residual)
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            expected, abs=tolerances.cross_solver
        )

    def test_a_zero_amplitude_reduces_to_the_independent_complex_gaussian(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The reduction that catches a normalisation counted once too often.

        At ``amplitude = 0`` the GP marginal must be the i.i.d. complex
        Gaussian exactly. It is the cheapest detector of the factor-of-two
        errors this density invites: ``log|S|`` counted once rather than twice,
        or ``-log(2 pi)`` rather than ``-log(2 pi) - log(sigma**2)``.
        """
        predicted, observed = visibility_pair()
        flat = CovarianceSpec(KernelFamily.MATERN32, 0.0, LENGTH_SCALE, axes=("u", "v"))
        correlated = circular_gp_likelihood(backend, flat)
        independent = Likelihood(ComplexGaussianFamily(), backend.independent_noise())
        assert correlated.log_prob(predicted, observed) == pytest.approx(
            independent.log_prob(predicted, observed), abs=tolerances.linear_algebra
        )

    def test_the_o_n_path_is_refused_by_name_word_for_word(
        self, backend: ConformanceBackend
    ) -> None:
        """``REQUIRES_ORDERED_1D`` cannot hold on this modality, so say so.

        Not the generic ordered-1D message, whose advice ("select exactly one
        axis to reach the O(N) path") cannot be taken here: a visibility is a
        point of the (u, v) plane at a wavelength, and no ordering of the plane
        makes a stationary kernel a function of one coordinate.
        """
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        _, observed = visibility_pair()
        noise = GaussianProcessNoise(
            backend.kernel(UV_SPEC), backend.gp_solver(SolverKind.QUASISEP)
        )
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(ComplexGaussianFamily(), observed)
        assert str(excinfo.value) == (
            "the complex_gaussian family on a complex VisibilitySet is the circular complex GP: "
            "the real and imaginary parts are two independent real processes sharing one "
            "covariance, so the solve takes a two-column right-hand side, and QuasisepGP "
            "declares STACKED_RESIDUALS = False. QuasisepGP could not take them even in "
            "principle: it needs one ordered one-dimensional coordinate, and a visibility lives "
            "at a point of the (u, v) plane at a wavelength, which no ordering reduces to one "
            "coordinate. Use DenseGP, with the kernel selecting the axes it acts on -- "
            'Matern32(axes=("u", "v")) for an isotropic (u, v) kernel, or '
            'Product(Matern32(axes=("u", "v")), Matern32(axes=("spectral_axis",))) for an error '
            "that is smooth in (u, v) and sharp in wavelength."
        )

    def test_a_bare_kernel_is_still_refused_on_a_visibility_set(
        self, backend: ConformanceBackend
    ) -> None:
        """W4.5's per-leaf unit rule, on the kind it was written for.

        The message names the remedy since W4.2 — the ``axes=`` selector — and
        this row is why that mattered: before the selector existed the advice
        "use one axis" was the fix, and on this container it is not.
        """
        _, observed = visibility_pair()
        noise = GaussianProcessNoise(
            backend.kernel(CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE)),
            backend.gp_solver(SolverKind.DENSE),
        )
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(ComplexGaussianFamily(), observed)
        message = str(excinfo.value)
        assert "carry different units ['', 'um']" in message
        assert "axes=(...)" in message

    def test_a_real_prediction_against_complex_data_is_still_caught(
        self, backend: ConformanceBackend
    ) -> None:
        """Gap I-1's dtype check, under a GP (the item's explicit criterion).

        Kind equality does not imply comparability: a ``VisibilitySet`` is legal
        with real or complex values, so a model that took the modulus somewhere
        would otherwise be fitted against complex data with no warning anywhere.
        The GP path must not have bypassed it — it is the one path that now
        reshapes the residual before anything looks at it.
        """
        predicted, observed = visibility_pair()
        amplitudes = predicted.with_values(np.abs(np.asarray(predicted.values)))
        likelihood = circular_gp_likelihood(backend, UV_SPEC)
        with pytest.raises(LikelihoodError, match="real"):
            likelihood.check_alignment(amplitudes, observed)

    def test_the_pointwise_terms_are_one_per_sample_and_leave_one_out(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``conditional_loo`` for the pointwise group (``results.md`` §6).

        One term per **sample**, not per component: a complex visibility is one
        observation with two parts. The oracle refits each term from scratch —
        the conditional of sample *i* given every other retained sample, each
        component scored against the same conditional variance — so the row
        compares the Cholesky identity against the definition rather than
        against itself.
        """
        predicted, observed = visibility_pair()
        likelihood = circular_gp_likelihood(backend, UV_SPEC)
        terms = likelihood.pointwise_log_prob(predicted, observed)
        assert terms.shape == (VISIBILITY_SIZE,)

        points = visibility_points(observed)
        residual = np.asarray(observed.values) - np.asarray(predicted.values)
        total = covariance_matrix(UV_SPEC, points) + np.diag(np.full(points.shape[0], SIGMA**2))
        for index in range(VISIBILITY_SIZE):
            others = np.array([i for i in range(VISIBILITY_SIZE) if i != index])
            block = total[np.ix_(others, others)]
            cross = total[index, others]
            solved = np.linalg.solve(block, cross)
            variance = total[index, index] - cross @ solved
            expected = 0.0
            for component in (np.real(residual), np.imag(residual)):
                mean = component[others] @ solved
                expected += float(st.norm.logpdf(component[index], mean, math.sqrt(variance)))
            assert terms[index] == pytest.approx(expected, abs=tolerances.linear_algebra)

    def test_the_conditioned_mean_is_complex_and_localises_both_components(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The localisation diagnostic on a complex residual (§4.8, family C).

        The mean comes back **complex**, because a calibration error has a
        direction in the complex plane and the modulus of the conditioned mean
        would hide which way it went. The oracle is the textbook conditional
        mean, per component, on the same matrix.
        """
        predicted, observed = visibility_pair()
        likelihood = circular_gp_likelihood(backend, UV_SPEC)
        conditioned = likelihood.conditional(predicted, observed)
        assert conditioned.mean.dtype.kind == "c"
        assert conditioned.variance.shape == (VISIBILITY_SIZE,)

        points = visibility_points(observed)
        residual = np.asarray(observed.values) - np.asarray(predicted.values)
        kernel = covariance_matrix(UV_SPEC, points)
        total = kernel + np.diag(np.full(points.shape[0], SIGMA**2))
        expected = kernel @ np.linalg.solve(total, residual)
        assert conditioned.mean == pytest.approx(expected, abs=tolerances.linear_algebra)

    def test_a_draw_has_the_covariance_the_density_scores(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``sample`` under a GP, against the moments of the declared model.

        Three claims, and the second is the one a wrong draw would pass the
        others on: each component's covariance is ``K + diag(sigma**2)``, the
        components are **uncorrelated** with each other, and the correlation is
        genuinely present rather than a diagonal. Drawing one realisation and
        using it for both parts would satisfy the first and fail the second.
        """
        predicted, observed = visibility_pair()
        likelihood = circular_gp_likelihood(backend, UV_SPEC)
        likelihood.check_alignment(predicted, observed)
        noise = likelihood.noise.noise_params(
            observed,
            np.ones(VISIBILITY_SIZE, dtype=bool),
            likelihood.context(None),
            coordinates=visibility_points(observed),
        )
        mean = np.zeros(VISIBILITY_SIZE, dtype=complex)
        rng = np.random.default_rng(4242)
        draws = np.stack([likelihood.family.sample(mean, noise, rng) for _ in range(6000)])
        expected = covariance_matrix(UV_SPEC, visibility_points(observed)) + np.diag(
            np.full(VISIBILITY_SIZE, SIGMA**2)
        )
        # The standard error of a covariance entry from L draws of a Gaussian is
        # about sqrt((C_ii C_jj + C_ij**2) / L); the diagonal bounds it.
        scale = float(np.max(np.diag(expected)))
        error = scale * math.sqrt(2.0 / draws.shape[0])
        for component in (np.real(draws), np.imag(draws)):
            empirical = component.T @ component / draws.shape[0]
            assert np.max(np.abs(empirical - expected)) < tolerances.monte_carlo_sigmas * error
        cross = np.real(draws).T @ np.imag(draws) / draws.shape[0]
        assert np.max(np.abs(cross)) < tolerances.monte_carlo_sigmas * error
        # And the draw is correlated: the nearest pair of baselines shares far
        # more than sigma**2 would give them.
        assert np.max(np.abs(expected - np.diag(np.diag(expected)))) > SIGMA**2


# ---------------------------------------------------------------------------
# The first approximate solver: HilbertSpaceGP (W5.4)
# ---------------------------------------------------------------------------

#: The basis sizes the 1-D convergence rows sweep, coarsest first. Chosen for
#: this fixture rather than in general: on ``GP_GRID`` (half-extent 5.25) with
#: ``boundary_factor`` 2 the basis reaches an angular frequency
#: ``pi m / (2 c S)``, so ``m = 16`` barely resolves a length scale of 2 and
#: ``m = 128`` resolves it comfortably. The row asserts the *shape* of that
#: progression, never the numbers.
BASIS_SIZES: tuple[int, ...] = (16, 32, 64, 128)

#: The per-axis counts of the 2-D row, on the ``(u, v)`` extent of
#: :func:`visibility_pair`. The total basis is the square of each — 576
#: members at the finest, for nine baselines, which is the tensor product's
#: whole character: ``m`` grows as the *product* of the per-axis counts, so a
#: spectral method is a low-dimension method and says so.
UV_BASIS_SIZES: tuple[int, ...] = (4, 8, 16, 24)

#: What the 2-D row's finest basis actually reaches. Looser than the 1-D
#: default, and measured rather than chosen: the ``(u, v)`` box spans only
#: about 1.5 length scales, so the *boundary* error dominates again and more
#: basis members buy little — 32 per axis gets to 4e-2 and 24 to 6e-2. That is
#: the regime ``horizon_notes.md`` §2 says a Vecchia approximation is for, and
#: a row that hid it behind the 1-D number would be asserting something untrue
#: about the method in two axes.
UV_FINAL = 0.1

#: The box for each. A wider box needs more basis members to reach the same
#: frequency, so the 2-D row (whose per-axis counts are small) pays for its
#: boundary accuracy with a factor of three rather than two.
BOUNDARY_FACTOR = 2.0
UV_BOUNDARY_FACTOR = 3.0


def hilbert_likelihood(
    backend: ConformanceBackend,
    covariance: CovarianceSpec = MATERN32,
    *,
    basis_size: Any = 32,
    boundary_factor: float = BOUNDARY_FACTOR,
    family: Any = None,
) -> Likelihood:
    """A likelihood over *backend*'s kernel and its reduced-rank solver."""
    solver = backend.gp_solver(
        SolverKind.HILBERT, basis_size=basis_size, boundary_factor=boundary_factor
    )
    noise = GaussianProcessNoise(backend.kernel(covariance), solver)
    return Likelihood(GaussianFamily() if family is None else family, noise)


def converges(
    errors: Sequence[float],
    sizes: Sequence[int],
    tolerances: Tolerances,
    *,
    order: float | None = None,
    final: float | None = None,
) -> None:
    """Assert the ``EXACT = False`` tolerance class on one measured sweep.

    Two claims, both from :func:`approximation_envelope`: every refinement
    sits inside an envelope anchored on the coarsest error and tightening as
    the basis grows, and the finest setting actually reaches
    ``tolerances.approximation_final``. Neither is a fixed number for the
    quantity under test, which is the whole point — an approximate solver
    held to a constant would be held to one fixture's arithmetic.

    ``order`` and ``final`` let a row state what its own **kernel** predicts,
    because the rate is a property of the kernel rather than of the solver: an
    isotropic Matérn-``nu``'s spectral density decays as
    ``omega ** -(2 nu + d)``, so its reduced-rank error falls as ``m ** -2nu``
    — order 1 for Matérn-1/2, 3 for Matérn-3/2, 5 for Matérn-5/2. The defaults
    are the slowest case, so a row that passes neither is asserting the
    weakest honest claim rather than a convenient one.
    """
    coarsest = float(errors[0])
    for size, error in zip(sizes[1:], errors[1:], strict=True):
        allowed = approximation_envelope(coarsest, sizes[0], size, tolerances, order)
        assert error <= allowed, (
            f"the approximation's error at basis size {size} is {error:.3e}, outside the "
            f"envelope {allowed:.3e} that the error {coarsest:.3e} at size {sizes[0]} sets. "
            f"The sweep was {[f'{value:.3e}' for value in errors]}."
        )
    reached = tolerances.approximation_final if final is None else float(final)
    assert float(errors[-1]) <= reached, (
        f"the approximation's error at the finest basis size {sizes[-1]} is "
        f"{errors[-1]:.3e}, which does not reach {reached:.3e}. The sweep was "
        f"{[f'{value:.3e}' for value in errors]}."
    )


class TestApproximateSolverConvergence:
    """``HilbertSpaceGP``↔``DenseGP``: the convergence rows W5.4 adds.

    ``horizon_notes.md`` §2 asked what a conformance row can assert about an
    ``EXACT = False`` solver, given that the battery compares the exact ones
    bit-tightly and an approximation cannot meet that. The answer these rows
    implement is that the assertion is about the *method*, not about a
    number: refine the approximation and the disagreement with ``DenseGP``
    must fall, at a rate the method's own analysis predicts, down to the floor
    its finite box imposes; and the finest setting must actually be close.

    A wrong implementation fails this in a way a fixed tolerance would not
    catch. A basis built on the wrong box, a spectral density with the wrong
    dimension in it, a Woodbury solve missing a factor — each gives an error
    that is *stable* under refinement rather than falling, because more basis
    members converge to the wrong process just as happily as to the right one.
    """

    def _skip_unless_available(self, backend: ConformanceBackend) -> None:
        if SolverKind.HILBERT not in backend.capabilities.solvers:
            pytest.skip(
                f"backend {backend.name!r} declares no reduced-rank spectral solver. Supply one "
                f"and add SolverKind.HILBERT to its capabilities; these rows then compare it "
                f"against DenseGP at a sequence of basis sizes under the approximation "
                f"tolerance class."
            )

    def test_the_solver_declares_itself_approximate(self, backend: ConformanceBackend) -> None:
        """The declaration the rest of this class relies on, on every column."""
        self._skip_unless_available(backend)
        solver = backend.gp_solver(SolverKind.HILBERT, basis_size=32)
        assert solver.NAME == "HilbertSpaceGP"
        assert not solver.EXACT
        assert solver.IMPLEMENTED
        assert solver.STACKED_RESIDUALS
        # The approximation is visible in the run's attrs, per fold-in 10.
        config = solver.provenance_config()
        assert config["basis_size"] == [32]
        assert config["boundary_factor"] == pytest.approx(2.0)

    def test_the_marginal_likelihood_converges_to_the_dense_one(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The headline row: ``log p`` against ``DenseGP``'s, at four basis sizes."""
        self._skip_unless_available(backend)
        predicted, observed = spectra()
        exact = gp_likelihood(backend, MATERN32, SolverKind.DENSE).log_prob(predicted, observed)
        errors = [
            abs(
                hilbert_likelihood(backend, MATERN32, basis_size=size).log_prob(predicted, observed)
                - exact
            )
            for size in BASIS_SIZES
        ]
        converges(errors, BASIS_SIZES, tolerances)

    @pytest.mark.parametrize(
        ("spec", "final"),
        [
            (CovarianceSpec(KernelFamily.MATERN12, AMPLITUDE, LENGTH_SCALE), 0.25),
            (CovarianceSpec(KernelFamily.MATERN52, AMPLITUDE, LENGTH_SCALE), None),
            (SQUARED_EXPONENTIAL, None),
            (CovarianceSpec(KernelFamily.SHO, AMPLITUDE, LENGTH_SCALE), None),
            (
                CovarianceSpec(
                    KernelFamily.SUM,
                    terms=(
                        CovarianceSpec(KernelFamily.MATERN32, 0.3, 3.0),
                        CovarianceSpec(KernelFamily.MATERN52, 0.2, 1.5),
                    ),
                ),
                None,
            ),
        ],
        ids=["matern12", "matern52", "squared_exponential", "sho", "sum"],
    )
    def test_every_family_with_a_spectral_density_converges(
        self,
        backend: ConformanceBackend,
        tolerances: Tolerances,
        spec: CovarianceSpec,
        final: float | None,
    ) -> None:
        """One row per family W5.4 gives a closed-form spectral density.

        The squared exponential is the interesting column: it has no
        quasiseparable form at all, so ``QuasisepGP`` refuses it and this is
        the only strategy that makes it scale. Its spectral density is a
        Gaussian, so it converges fastest of the five.

        Matérn-1/2 is the interesting column in the other direction, and its
        looser ``final`` is a *measurement rather than an excuse*: a
        Matérn-1/2 is the roughest process this method supports, its spectral
        density decays only as ``omega ** -2``, and the reduced-rank error
        therefore falls as ``1/m`` — so 128 basis members over this grid still
        leave a couple of tenths of a nat. That is exactly the regime
        ``horizon_notes.md`` §2 says a Vecchia approximation is for, and a row
        that hid it behind the same number as the smooth families would be
        asserting something untrue about the method.
        """
        self._skip_unless_available(backend)
        predicted, observed = spectra()
        exact = gp_likelihood(backend, spec, SolverKind.DENSE).log_prob(predicted, observed)
        errors = [
            abs(
                hilbert_likelihood(backend, spec, basis_size=size).log_prob(predicted, observed)
                - exact
            )
            for size in BASIS_SIZES
        ]
        converges(errors, BASIS_SIZES, tolerances, final=final)

    def test_the_conditioned_moments_converge(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """§4.8's localisation diagnostic must converge too, on and off the grid."""
        self._skip_unless_available(backend)
        predicted, observed = spectra()
        expected = gp_likelihood(backend, MATERN32, SolverKind.DENSE).conditional(
            predicted, observed
        )
        means: list[float] = []
        variances: list[float] = []
        for size in BASIS_SIZES:
            got = hilbert_likelihood(backend, MATERN32, basis_size=size).conditional(
                predicted, observed
            )
            means.append(float(np.max(np.abs(np.asarray(got.mean) - np.asarray(expected.mean)))))
            variances.append(
                float(np.max(np.abs(np.asarray(got.variance) - np.asarray(expected.variance))))
            )
        converges(means, BASIS_SIZES, tolerances)
        converges(variances, BASIS_SIZES, tolerances)
        # And the posterior variance is a variance at every basis size, which
        # a subtracted quadratic form would not guarantee near the boundary.
        for size in BASIS_SIZES:
            got = hilbert_likelihood(backend, MATERN32, basis_size=size).conditional(
                predicted, observed
            )
            assert np.all(np.asarray(got.variance) >= 0.0)

    def test_the_leave_one_out_terms_converge_rather_than_refusing(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """W5.4 allowed a refusal and did not need one: Woodbury gives ``A_ii``.

        The leave-one-out identity needs the diagonal of ``(K + diag(sigma**2))**-1``,
        which is exactly what ``QuasisepGP`` has no O(N) route to and refuses
        over. A reduced-rank representation has one in closed form, so these
        terms exist and converge with everything else.
        """
        self._skip_unless_available(backend)
        predicted, observed = spectra()
        expected = gp_likelihood(backend, MATERN32, SolverKind.DENSE).pointwise_log_prob(
            predicted, observed
        )
        errors: list[float] = []
        for size in BASIS_SIZES:
            got = hilbert_likelihood(backend, MATERN32, basis_size=size).pointwise_log_prob(
                predicted, observed
            )
            assert got.shape == expected.shape
            errors.append(float(np.max(np.abs(got - expected))))
        converges(errors, BASIS_SIZES, tolerances)

    def test_the_tensor_product_basis_converges_on_a_visibility_set(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The 2-D row, on W4.2's ``(u, v)`` fixture — where ``QuasisepGP`` cannot go.

        ``Matern32(axes=("u", "v"))`` on a three-axis ``VisibilitySet``, scored
        by the circular complex Gaussian, so the right-hand side is the
        ``(n, 2)`` block of ``STACKED_RESIDUALS`` as well. The basis is the
        tensor product over the box spanning the ``(u, v)`` extent, and ``m``
        is the square of the per-axis count — which is why a spectral method
        is a low-dimension method and why this row's counts are small.
        """
        self._skip_unless_available(backend)
        if not backend.capabilities.complex_models:
            pytest.skip(f"backend {backend.name!r} declares no complex models")
        predicted, observed = visibility_pair()
        exact = circular_gp_likelihood(backend, UV_SPEC, SolverKind.DENSE).log_prob(
            predicted, observed
        )
        errors = [
            abs(
                hilbert_likelihood(
                    backend,
                    UV_SPEC,
                    basis_size=(size, size),
                    boundary_factor=UV_BOUNDARY_FACTOR,
                    family=ComplexGaussianFamily(),
                ).log_prob(predicted, observed)
                - exact
            )
            for size in UV_BASIS_SIZES
        ]
        converges(errors, UV_BASIS_SIZES, tolerances, final=UV_FINAL)

    def test_the_latent_block_is_the_basis_size_and_one_whitening_serves_both_paths(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """W5.4's ruling on ``inference.md`` §17.4, asserted rather than asserted about.

        Two halves, and they are the same fact seen twice. The declaration:
        ``Likelihood.latent_declaration(n)`` gives ``m`` whitened variables,
        not ``n``, because the solver says so. The arithmetic: the covariance
        the marginal likelihood **scores** is exactly the covariance the
        whitening **draws** from — ``scipy`` is handed
        ``L Lᵀ + diag(sigma**2)`` built from the solver's own ``L`` and must
        return the solver's own ``log_prob``. That is what ``simulate(observe=True)``
        and the latent-GP path agreeing *means*, written so that a solver
        whose two paths had drifted apart could not pass it.
        """
        self._skip_unless_available(backend)
        predicted, observed = spectra()
        size = 24
        likelihood = hilbert_likelihood(backend, MATERN32, basis_size=size)
        noise = likelihood.noise
        assert isinstance(noise, GaussianProcessNoise)
        assert noise.solver.latent_size(noise.kernel, observed.n_samples) == size
        # The declaration half needs a family that actually takes the latent
        # path: a Gaussian under a GP marginalises analytically and declares no
        # latent values at all, which is §9's table rather than anything to do
        # with this solver. Poisson + GP is the latent shape.
        latent = hilbert_likelihood(
            backend, MATERN32, basis_size=size, family=family_named("poisson")()
        )
        assert latent.marginalisation is Marginalisation.LATENT
        declaration = latent.latent_declaration(observed.n_samples)
        # ``size`` stays the retained-sample count -- it is what the mask
        # invariance is checked against -- while the whitened block the
        # sampler carries is the solver's basis size.
        assert declaration.size == observed.n_samples
        assert declaration.whitened_size == size
        assert declaration.parameter.shape == (size,)

        points = np.asarray(observed.axes[0].values, dtype=float).reshape(-1, 1)
        values = noise.kernel.resolve({})
        # The whitening as a block: L = latent_transform applied to the
        # identity, one column per whitened variable.
        factor = np.asarray(
            noise.solver.latent_transform(noise.kernel, points, np.eye(size), values), dtype=float
        )
        assert factor.shape == (observed.n_samples, size)
        covariance = factor @ factor.T + np.diag(np.full(observed.n_samples, SIGMA**2))
        residual = np.asarray(observed.values, dtype=float) - np.asarray(
            predicted.values, dtype=float
        )
        oracle = float(
            multivariate_normal.logpdf(residual, mean=np.zeros(observed.n_samples), cov=covariance)
        )
        assert likelihood.log_prob(predicted, observed) == pytest.approx(
            oracle, abs=tolerances.linear_algebra
        )

    def test_a_kernel_with_no_spectral_density_is_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        """The declared-slot discipline, one level down from ``QuasisepGP``'s.

        A ``Product``'s leaves each have a closed-form spectral density and
        the product does not — its transform is a convolution — so the refusal
        has to name the *node* rather than the leaf, and this row is what says
        so on every column.
        """
        self._skip_unless_available(backend)
        product = CovarianceSpec(
            KernelFamily.PRODUCT,
            terms=(
                CovarianceSpec(KernelFamily.MATERN32, AMPLITUDE, LENGTH_SCALE),
                CovarianceSpec(KernelFamily.MATERN12, 0.2, 0.5),
            ),
        )
        _, observed = spectra()
        for spec, named in (
            (product, "Product"),
            (CovarianceSpec(KernelFamily.ROTATION, AMPLITUDE, LENGTH_SCALE), "RotationTerm"),
        ):
            likelihood = hilbert_likelihood(backend, spec)
            with pytest.raises(LikelihoodError) as refusal:
                likelihood.noise.check_compatible(likelihood.family, observed)
            message = str(refusal.value)
            assert named in message
            assert "spectral density" in message
            assert "DenseGP" in message

    def test_a_basis_size_of_the_wrong_rank_is_refused_at_composition(
        self, backend: ConformanceBackend
    ) -> None:
        """One count per selected axis, checked against the container, not the arithmetic.

        At **composition**, which is where it has to be: the latent block's
        size comes from this declaration, so a rank that disagreed with the
        container would be discovered as a matrix-shape error inside a solve
        after a sampler had already been built around the wrong dimension.
        """
        self._skip_unless_available(backend)
        _, observed = spectra()
        likelihood = hilbert_likelihood(backend, MATERN32, basis_size=(8, 8))
        with pytest.raises(LikelihoodError) as refusal:
            likelihood.noise.check_compatible(likelihood.family, observed)
        assert "tensor product" in str(refusal.value)
