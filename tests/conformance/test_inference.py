"""Conformance rows for the dataset, fitting-problem and inference contract.

The inventory is ``inference.md`` §18's hand-down, verbatim: "``log_prob ==
log_prior + log_likelihood``; an out-of-support θ returns ``-inf`` without
evaluating the model; ``prior_transform`` at the centre of the cube is the prior
median; ``constrain ∘ unconstrain`` is the identity and
``log_prob_unconstrained`` differs from ``log_prob`` by exactly the summed
Jacobian; a tie costs one dimension and has one binding per site; the joint
log-likelihood equals the sum of ``contributions``; a Gaussian
``simulate(observe=True)`` has empirical covariance ``K + diag(σ²)``; the same
seed gives the same simulation; ``substream`` is stable across processes."
"""

from __future__ import annotations

import os
import subprocess
import sys
from typing import Any

import numpy as np
import pytest

from ampere.core import (
    Dataset,
    FittingProblem,
    Likelihood,
    Tie,
    family_named,
    log_likelihood_terms_of,
    realise,
    registered_realisations,
    substream,
)

from .composition import (
    GP_GRID,
    DatasetSpec,
    NoiseKind,
    ProblemSpec,
    analytic_flux,
    build_instrument,
    build_problem,
    model_context,
    observed_container,
)
from .oracles import kernel_matrix, summed_log_abs_det
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    ModelKind,
    ModelSpec,
    SolverKind,
    Tolerances,
    TransformationKind,
    TransformationSpec,
)

CALIBRATION = TransformationSpec(TransformationKind.SCALE, label="calibration")

SINGLE = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(DatasetSpec(instrument=(CALIBRATION,)),),
)

JOINT = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, channels=("blue", "red"), coordinates=GP_GRID),
    datasets=(
        DatasetSpec(label="blue", channel="blue", instrument=(CALIBRATION,)),
        DatasetSpec(label="red", channel="red", instrument=(CALIBRATION,), data_seed=771),
    ),
    ties=(
        Tie(
            "calibration",
            ("blue.instrument.calibration.scale", "red.instrument.calibration.scale"),
        ),
    ),
)

CORRELATED = ProblemSpec(
    model=ModelSpec(kind=ModelKind.LINEAR, coordinates=GP_GRID),
    datasets=(DatasetSpec(noise=NoiseKind.GP, covariance=CovarianceSpec()),),
)

#: Readable pytest ids for the three shapes above, used where a row is
#: parametrised over all of them.
SPEC_IDS = ("single", "joint", "correlated")


@pytest.fixture
def problem(backend: ConformanceBackend) -> FittingProblem:
    """A one-dataset power-law problem with a calibration nuisance parameter."""
    return build_problem(backend, SINGLE)


@pytest.fixture
def joint(backend: ConformanceBackend) -> FittingProblem:
    """Two datasets on two channels, sharing one calibration through a tie."""
    return build_problem(backend, JOINT)


class TestBackendIdentity:
    """W2.12 item 4: one name per backend, everywhere.

    The string a backend's models and transformations declare as ``BACKEND``
    is the same string as its conformance fixture's ``name``, the same string
    ``lowering.md`` §12.8's registry is keyed on, and the same string a run's
    ``ampere_backend`` carries. A fixture whose parts declared anything else
    would be composing a problem that reports a backend nobody can look up.
    """

    def test_a_composed_problem_reports_this_fixture_name_as_its_backend(
        self, backend: ConformanceBackend
    ) -> None:
        # Every ProblemSpec in this module, so a fixture cannot pass by
        # declaring the flag on its models and forgetting its instrument steps.
        for spec in (SINGLE, JOINT, CORRELATED):
            problem = build_problem(backend, spec)
            assert problem.backend == backend.name
            assert problem.capabilities.backend == backend.name

    def test_every_part_declares_it_rather_than_inheriting_the_default(
        self, backend: ConformanceBackend
    ) -> None:
        # The point of the row above is lost if it passes because everything
        # silently inherited "reference": check the parts themselves, which is
        # what ``declared_capabilities`` reads.
        #
        # CORRELATED is here as well as JOINT since **W2.13**: the parts set
        # widened to the noise model and the GP solver, and a GP problem is the
        # only shape in this module that has a solver at all.
        for spec in (JOINT, CORRELATED):
            problem = build_problem(backend, spec)
            parts = (*problem.models.values(), *problem.datasets.capability_parts)
            assert parts, "a composed problem with no capability parts proves nothing"
            assert {part.BACKEND for part in parts} == {backend.name}

    def test_the_likelihoods_noise_model_and_solver_are_parts(
        self, backend: ConformanceBackend
    ) -> None:
        """W2.13 fold-in 7: the widened parts set, per fixture.

        Before this, a problem could report ``backend="torch"`` and
        ``differentiable=True`` while its noise model and GP solve ran in
        numpy — both Phase 2 tracks recorded it. The flags now cover the whole
        of what one evaluation passes through.
        """
        plain = build_problem(backend, SINGLE)
        for dataset in plain.datasets.values():
            assert dataset.likelihood.noise in dataset.capability_parts

        correlated = build_problem(backend, CORRELATED)
        for dataset in correlated.datasets.values():
            noise = dataset.likelihood.noise
            assert noise in dataset.capability_parts
            assert noise.solver in dataset.capability_parts

    def test_the_family_is_not_a_part(self, backend: ConformanceBackend) -> None:
        """And deliberately so: a family has no backend of its own to declare."""
        problem = build_problem(backend, SINGLE)
        for dataset in problem.datasets.values():
            assert dataset.likelihood.family not in dataset.capability_parts

    def test_the_declared_flags_match_what_the_fixture_claims(
        self, backend: ConformanceBackend
    ) -> None:
        # The other three flags, checked at the same time and for the same
        # reason: BackendCapabilities is the fixture's claim about itself, and
        # the composed problem is what the parts actually declare.
        problem = build_problem(backend, SINGLE)
        claimed = backend.capabilities
        assert problem.differentiable == claimed.differentiable
        assert problem.batchable == claimed.batchable
        assert problem.device == claimed.device


class TestTheDecomposition:
    """``log_prob``, ``log_prior``, ``log_likelihood`` and ``contributions``."""

    def test_log_prob_is_the_sum_of_its_two_halves(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        theta = problem.prior_transform(np.array([0.4, 0.7, 0.55]))
        assert problem.log_prob(theta) == pytest.approx(
            problem.log_prior(theta) + problem.log_likelihood(theta), abs=tolerances.exact
        )

    def test_one_evaluation_reports_the_same_decomposition(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        theta = problem.prior_transform(np.array([0.4, 0.7, 0.55]))
        evaluation = problem.evaluate(theta)
        assert evaluation.failure is None
        assert evaluation.log_prob == pytest.approx(
            evaluation.log_prior + evaluation.log_likelihood, abs=tolerances.exact
        )
        assert evaluation.log_prob == pytest.approx(problem.log_prob(theta), abs=tolerances.exact)

    def test_the_contributions_sum_to_the_joint_log_likelihood(
        self, joint: FittingProblem, tolerances: Tolerances
    ) -> None:
        theta = joint.prior_transform(np.full(joint.free_size, 0.45))
        evaluation = joint.evaluate(theta)
        assert set(evaluation.contributions) == {"blue", "red"}
        assert sum(evaluation.contributions.values()) == pytest.approx(
            evaluation.log_likelihood, abs=tolerances.exact
        )

    def test_each_contribution_is_that_datasets_own_likelihood(
        self, backend: ConformanceBackend, joint: FittingProblem, tolerances: Tolerances
    ) -> None:
        """The per-*dataset* decomposition is well defined even where per-sample is not."""
        theta = joint.prior_transform(np.full(joint.free_size, 0.45))
        contributions = joint.evaluate(theta).contributions
        for label in ("blue", "red"):
            single = build_problem(
                backend,
                ProblemSpec(
                    model=JOINT.model,
                    datasets=tuple(d for d in JOINT.datasets if d.label == label),
                ),
            )
            routed = joint.mapping.distribute(joint.parameters.complete(theta))
            values = {f"model.{name}": value for name, value in model_context(joint, theta).items()}
            values[f"{label}.instrument.calibration.scale"] = routed[label][
                "instrument.calibration.scale"
            ]
            assert single.log_likelihood(values) == pytest.approx(
                contributions[label], abs=tolerances.exact
            )


class TestOutOfSupport:
    """A rejected θ costs nothing, which for an expensive simulator is the point."""

    def test_it_returns_minus_infinity_without_evaluating_the_model(
        self, problem: FittingProblem
    ) -> None:
        model = problem.model
        model.reset_evaluations()
        # `norm` is LogUniform(0.1, 10); 500 is outside its support.
        rejected = dict(problem.reference_values)
        rejected["model.norm"] = 500.0

        assert problem.log_prob(rejected) == -np.inf
        assert model.evaluations == 0

    def test_the_evaluation_records_no_failure_for_a_rejected_draw(
        self, problem: FittingProblem
    ) -> None:
        """Zero prior mass is an answer, not a failure (``results.md`` §14)."""
        rejected = dict(problem.reference_values)
        rejected["model.norm"] = 500.0
        evaluation = problem.evaluate(rejected)

        assert evaluation.log_prior == -np.inf
        assert np.isnan(evaluation.log_likelihood)
        assert evaluation.failure is None
        assert not evaluation.failed
        assert dict(evaluation.contributions) == {}

    def test_an_accepted_draw_does_evaluate_the_model(self, problem: FittingProblem) -> None:
        model = problem.model
        model.reset_evaluations()
        problem.log_prob(problem.reference_values)
        assert model.evaluations == 1


class TestReferencePoint:
    """The cube's centre is the prior median, and it is the default θ."""

    def test_the_centre_of_the_cube_is_the_prior_median(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        centre = np.full(problem.free_size, 0.5)
        expected = [
            parameter.prior.median()
            for parameter in problem.parameters.parameters
            if not parameter.is_fixed
        ]
        assert problem.prior_transform(centre) == pytest.approx(expected, abs=tolerances.analytic)

    def test_the_default_values_are_that_median(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        centre = problem.prior_transform(np.full(problem.free_size, 0.5))
        assert problem.log_prob() == pytest.approx(problem.log_prob(centre), abs=tolerances.exact)


class TestUnconstrainedSpace:
    """The identity, and the exact size of the change-of-variables term."""

    def test_constrain_and_unconstrain_are_inverse(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        y = np.array([0.3, -0.8, 0.45])[: problem.free_size]
        theta = problem.constrain(y)
        assert problem.unconstrain(theta) == pytest.approx(y, abs=tolerances.linear_algebra)

    def test_the_unconstrained_log_prob_differs_by_exactly_the_jacobian(
        self, problem: FittingProblem, tolerances: Tolerances
    ) -> None:
        y = np.array([0.3, -0.8, 0.45])[: problem.free_size]
        theta = problem.constrain(y)
        jacobian = summed_log_abs_det(problem.parameters, y)
        assert problem.log_prob_unconstrained(y) - problem.log_prob(theta) == pytest.approx(
            jacobian, abs=tolerances.analytic
        )


class TestTheRealisation:
    """``inference.md`` §10a: the backend's differentiable native form.

    One row per **registered** realisation, and none at all for a backend that
    has none — the reference backend registers nothing in v1, deliberately, so
    these skip there rather than failing. The claim is the one §10a's "What the
    conformance suite owes" paragraph makes: the realised density agrees with
    the numpy contract path **at many points, including near a support
    boundary**, at ``tolerances.cross_backend``. ``realise``'s own one-point
    check is a guard; this is the proof, and the numpy path is the oracle.
    """

    def realised(self, problem: FittingProblem) -> Any:
        if problem.backend not in registered_realisations():
            pytest.skip(
                f"the {problem.backend!r} backend registers no realisation "
                f"(inference.md §10a: the reference backend has no differentiable path)"
            )
        return realise(problem)

    @staticmethod
    def points(problem: FittingProblem) -> list[np.ndarray]:
        """Unconstrained vectors spanning the interesting part of the space.

        A grid over the unit cube pushed through ``prior_transform`` and back
        into unconstrained space, plus a scatter, **plus two points a long way
        out**. The last group is what "including near a support boundary"
        means for a backend whose bijections are exponentials and logits: a
        constrained parameter approaches its boundary as its unconstrained
        coordinate goes to infinity, so a large ``|y|`` is how a boundary is
        reached at all — and it is where the change-of-variables term dominates
        the density, so an argument-order error in ``log_abs_det_jacobian``
        shows here and nowhere near the median.
        """
        size = problem.free_size
        rng = np.random.default_rng(20260907)
        cube = np.concatenate(
            [
                np.full((1, size), 0.5),
                np.linspace(0.05, 0.95, 9).reshape(-1, 1) * np.ones((1, size)),
                rng.uniform(0.02, 0.98, size=(12, size)),
            ]
        )
        interior = [problem.unconstrain(problem.prior_transform(row)) for row in cube]
        return [*interior, np.full(size, -8.0), np.full(size, 8.0)]

    def test_it_names_the_problems_backend_and_free_size(self, problem: FittingProblem) -> None:
        realised = self.realised(problem)
        assert realised.backend == problem.backend
        assert int(realised.free_size) == problem.free_size

    @pytest.mark.parametrize("spec", [SINGLE, JOINT, CORRELATED], ids=SPEC_IDS)
    def test_the_realised_density_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances, spec: ProblemSpec
    ) -> None:
        """Many points, every declared shape, against the oracle."""
        problem = build_problem(backend, spec)
        realised = self.realised(problem)
        for y in self.points(problem):
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            if not np.isfinite(expected):
                # A bare -inf is all the native path promises (§10a narrows
                # §11): the *reasons* are recovered post hoc, on the numpy
                # path, for stored draws.
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)

    def test_the_optional_decomposition_sums_to_the_joint_likelihood(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """§10a's optional member, where a backend supplies it."""
        joint = build_problem(backend, JOINT)
        realised = self.realised(joint)
        terms = log_likelihood_terms_of(realised)
        if terms is None:
            pytest.skip(
                f"the {joint.backend!r} realisation supplies no log_likelihood_terms; it is "
                f"optional (inference.md §10a) and a driver recomputes the split instead"
            )
        y = joint.unconstrain(joint.reference_values)
        found = {
            label: float(np.asarray(backend.to_numpy(value))) for label, value in terms(y).items()
        }
        assert set(found) == set(joint.datasets)
        evaluation = joint.evaluate(joint.reference_values)
        for label, value in found.items():
            assert value == pytest.approx(
                evaluation.contributions[label], abs=tolerances.cross_backend
            )

    def test_a_numpy_solver_makes_a_native_problem_a_backend_disagreement(
        self, backend: ConformanceBackend
    ) -> None:
        """W2.13 fold-in 7's loud consequence, checked on every fixture.

        A fixture whose own dense solver *is* ``ampere.core``'s skips: there is
        nothing to disagree about, because every part is the reference
        backend's. On a native fixture the core solver is a different backend,
        and the refusal is the point — a nominally differentiable problem whose
        GP solve runs in numpy is one whose GP hyperparameters get no gradient
        at all, which no other check in this suite would notice.
        """
        from ampere.core import DenseGP as CoreDenseGP
        from ampere.core.exceptions import DatasetError

        native = backend.gp_solver(SolverKind.DENSE)
        if type(native) is CoreDenseGP:
            pytest.skip("this fixture's dense solver is ampere.core's, so there is no mismatch")
        assert native.BACKEND == backend.name
        assert CoreDenseGP.BACKEND == "reference"

        declared = CORRELATED.datasets[0]
        observed = observed_container(CORRELATED, declared)
        likelihood = Likelihood(
            family_named(declared.family)(),
            backend.gp_noise(backend.kernel(declared.covariance), CoreDenseGP()),
        )
        with pytest.raises(DatasetError, match="different backends"):
            FittingProblem(
                backend.model(CORRELATED.model),
                [
                    Dataset(
                        observed,
                        build_instrument(backend, declared),
                        likelihood,
                        label=declared.label,
                    )
                ],
                seed=CORRELATED.seed,
            )


class TestTying:
    """One shared quantity, one dimension, one binding per site."""

    def test_a_tie_costs_one_dimension(
        self, backend: ConformanceBackend, joint: FittingProblem
    ) -> None:
        untied = build_problem(backend, ProblemSpec(model=JOINT.model, datasets=JOINT.datasets))
        assert untied.free_size == joint.free_size + 1
        assert joint.tied_names == ("calibration",)

    def test_the_tie_has_one_binding_per_site(self, joint: FittingProblem) -> None:
        sites = joint.sites()["calibration"]
        assert sorted(sites) == [
            "blue.instrument.calibration.scale",
            "red.instrument.calibration.scale",
        ]

    def test_one_value_reaches_both_sites(self, joint: FittingProblem) -> None:
        theta = dict(joint.reference_values)
        theta["calibration"] = 1.3
        routed = joint.mapping.distribute(joint.parameters.complete(theta))
        for label in ("blue", "red"):
            assert routed[label]["instrument.calibration.scale"] == pytest.approx(1.3)


class TestSimulation:
    """The generative path: its covariance, and its reproducibility."""

    def test_a_gaussian_draw_has_the_declared_covariance(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``simulate(observe=True)`` realises ``K + diag(σ²)``, not something near it.

        The assertion is against the Monte-Carlo standard error of the sample
        covariance itself — ``se(Ĉᵢⱼ) = sqrt((Cᵢᵢ Cⱼⱼ + Cᵢⱼ²)/M)`` — so the
        tolerance stays honest as the draw count changes, rather than being a
        number chosen to make the row pass.
        """
        problem = build_problem(backend, CORRELATED)
        grid = np.asarray(GP_GRID, dtype=float)
        theta = dict(problem.reference_values)
        mean = analytic_flux(CORRELATED.model, model_context(problem, theta), grid)

        draws = 2000
        rng = np.random.default_rng(20260902)
        residuals = np.empty((draws, grid.size))
        for index in range(draws):
            simulation = problem.simulate(theta, observe=True, rng=rng)
            assert simulation.observations is not None
            residuals[index] = simulation.observations["sed"].values - mean

        covariance = CORRELATED.datasets[0].covariance
        expected = kernel_matrix(
            covariance.family, grid, covariance.amplitude, covariance.length_scale
        ) + np.diag(np.full(grid.size, CORRELATED.datasets[0].uncertainty ** 2))
        empirical = np.cov(residuals, rowvar=False)

        diagonal = np.diag(expected)
        standard_error = np.sqrt((np.outer(diagonal, diagonal) + expected**2) / draws)
        assert np.all(
            np.abs(empirical - expected) <= tolerances.monte_carlo_sigmas * standard_error
        )

    def test_the_draw_is_centred_on_the_prediction(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_problem(backend, CORRELATED)
        grid = np.asarray(GP_GRID, dtype=float)
        theta = dict(problem.reference_values)
        mean = analytic_flux(CORRELATED.model, model_context(problem, theta), grid)

        draws = 1000
        rng = np.random.default_rng(771)
        realisations = np.array(
            [
                problem.simulate(theta, observe=True, rng=rng).observations["sed"].values
                for _ in range(draws)
            ]
        )
        covariance = CORRELATED.datasets[0].covariance
        variance = (
            np.diag(
                kernel_matrix(
                    covariance.family, grid, covariance.amplitude, covariance.length_scale
                )
            )
            + CORRELATED.datasets[0].uncertainty ** 2
        )
        standard_error = np.sqrt(variance / draws)
        assert np.all(
            np.abs(realisations.mean(axis=0) - mean)
            <= tolerances.monte_carlo_sigmas * standard_error
        )

    def test_the_prediction_itself_is_the_analytic_form(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_problem(backend, CORRELATED)
        grid = np.asarray(GP_GRID, dtype=float)
        simulation = problem.simulate(problem.reference_values)
        expected = analytic_flux(CORRELATED.model, model_context(problem), grid)
        assert backend.to_numpy(simulation.predicted["sed"].values) == pytest.approx(
            expected, abs=tolerances.analytic
        )

    def test_the_same_seed_gives_the_same_simulation(self, backend: ConformanceBackend) -> None:
        first = build_problem(backend, SINGLE).simulate(observe=True)
        second = build_problem(backend, SINGLE).simulate(observe=True)
        assert first.theta.tolist() == second.theta.tolist()
        assert (
            first.observations["sed"].values.tolist() == second.observations["sed"].values.tolist()
        )

    def test_a_different_seed_gives_a_different_simulation(
        self, backend: ConformanceBackend
    ) -> None:
        other = build_problem(backend, ProblemSpec(SINGLE.model, SINGLE.datasets, seed=1))
        assert (
            build_problem(backend, SINGLE).simulate().theta.tolist()
            != other.simulate().theta.tolist()
        )

    def test_the_stream_advances_between_draws(self, backend: ConformanceBackend) -> None:
        problem = build_problem(backend, SINGLE)
        assert problem.simulate().theta.tolist() != problem.simulate().theta.tolist()


class TestSubstreams:
    """A run's derived randomness is a pure function of its seed and a label."""

    def test_distinct_labels_give_distinct_streams(self, problem: FittingProblem) -> None:
        labels = ("prior", "initialisation", "simulate", "posterior_predictive")
        first = {label: problem.rng(label).random() for label in labels}
        assert len(set(first.values())) == len(labels)

    def test_a_stream_is_reproducible_from_the_seed(self, backend: ConformanceBackend) -> None:
        """Two problems built from the same declaration and seed share a stream.

        Within one problem a named stream is a single advancing generator —
        asking for it twice does not rewind it — so reproducibility is a
        statement about two *runs*, which is what it is needed for.
        """
        first = build_problem(backend, SINGLE)
        second = build_problem(backend, SINGLE)
        assert first.rng("prior").random() == second.rng("prior").random()

    def test_a_different_seed_gives_a_different_stream(self, backend: ConformanceBackend) -> None:
        other = build_problem(backend, ProblemSpec(SINGLE.model, SINGLE.datasets, seed=1))
        assert build_problem(backend, SINGLE).rng("prior").random() != other.rng("prior").random()

    def test_substream_is_stable_across_processes(self) -> None:
        """Not merely within one interpreter: across ``PYTHONHASHSEED`` too.

        ``substream`` is BLAKE2b over the seed's bytes and the label's UTF-8,
        precisely so a run reproduces on another machine. A ``hash()``-based
        implementation would pass every in-process check and fail here.
        """
        script = (
            "from ampere.core import substream;"
            "print(substream(20260902, 'simulate'), substream(-7, 'prior.blue'))"
        )
        outputs = set()
        for hash_seed in ("0", "1", "random"):
            environment = dict(os.environ, PYTHONHASHSEED=hash_seed)
            result = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
                text=True,
                check=True,
                env=environment,
            )
            outputs.add(result.stdout.strip())

        assert len(outputs) == 1
        assert outputs == {f"{substream(20260902, 'simulate')} {substream(-7, 'prior.blue')}"}
