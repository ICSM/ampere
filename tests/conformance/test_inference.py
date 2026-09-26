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

import dataclasses
import os
import pickle
import subprocess
import sys
from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FailureReason,
    FittingProblem,
    Likelihood,
    ProcessExecutor,
    SerialExecutor,
    SimulationBatch,
    ThreadExecutor,
    Tie,
    family_named,
    log_likelihood_terms_of,
    realise,
    registered_realisations,
    substream,
)

from .composition import (
    COARSE_GRID,
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
    KernelFamily,
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

#: The correlated shape with a **fixed diagonal floor** on the GP noise model.
#:
#: W3.1 slice 2's parity row. ``GaussianProcessNoise`` has taken ``scale=`` and
#: ``jitter=`` since ``ampere.core`` declared them, and torch's subclass has
#: forwarded them since W2.4 — jax's did not, so the same three-line
#: composition succeeded on one modern backend and raised ``TypeError`` on the
#: other. That is exactly the drift the battery exists to catch, so the check
#: lives here rather than in either backend's own suite.
CORRELATED_JITTER = ProblemSpec(
    model=ModelSpec(kind=ModelKind.LINEAR, coordinates=GP_GRID),
    datasets=(DatasetSpec(noise=NoiseKind.GP, covariance=CovarianceSpec(), gp_jitter=0.05),),
)

#: The correlated shape sized so that a **draw** can be checked against its own
#: covariance without the tolerance having to be generous.
#:
#: W3.1 slice 2's native-sampling rows. ``sigma`` is deliberately large (0.3, so
#: ``sigma² = 0.09``) for the reason ``tests/core/test_dataset.py`` gives about
#: the numpy oracle it copies: with the 0.1 used elsewhere in this file,
#: dropping the diagonal term entirely would move the covariance by 0.01 and no
#: honest tolerance would catch it. The amplitude is raised to 0.5 for the same
#: reason on the other axis.
CORRELATED_DRAWS = ProblemSpec(
    model=ModelSpec(kind=ModelKind.LINEAR, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(
            noise=NoiseKind.GP,
            covariance=CovarianceSpec(amplitude=0.5, length_scale=2.0),
            uncertainty=0.3,
        ),
    ),
)

#: The same problem through the O(N) solver rather than the dense one.
#: Skipped by any backend that does not declare ``SolverKind.QUASISEP``, so it
#: is additive: a track that has not yet made ``DEVELOPMENT_PLAN.md`` §6's
#: measurement is unaffected. It exists because a realisation is a *whole*
#: lowering, and the GP solve is the part of it most likely to differ between
#: the differentiable and the numpy path — the two run genuinely different
#: recursions, which the marginal-agreement rows check in isolation and this
#: one checks inside the composed density.
CORRELATED_QUASISEP = ProblemSpec(
    model=ModelSpec(kind=ModelKind.LINEAR, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(noise=NoiseKind.GP, covariance=CovarianceSpec(), solver=SolverKind.QUASISEP),
    ),
)

#: A chain that **changes kind**: a spectrum integrated through filters, fitted
#: against ``PhotometricPoints``.
#:
#: Added at W2.5 slice 2, and it closed a real hole. Every other shape here
#: ends in a kind-preserving step, so a backend's native photometry path —
#: ``apply_flux``, the one a realised density actually calls — was composed
#: into no problem at all: the jax implementation of it raised
#: ``AttributeError`` on its first evaluation, and had done since it was
#: written, because nothing ever evaluated it. The photometry step the battery
#: declares is the fixture-local one (bare pivots and filter names,
#: ``protocol.py``'s "a chain whose kinds do not compose"), so what this shape
#: proves is that a backend's *chain* survives a kind change with its gradient
#: and its numbers intact; the shipped ``SyntheticPhotometry``'s own response
#: integrals stay ``tests/backends``'.
PHOTOMETRIC = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(
            instrument=(
                CALIBRATION,
                TransformationSpec(
                    TransformationKind.PHOTOMETRY,
                    label="synphot",
                    target=COARSE_GRID,
                    filters=("W1", "W2", "W3", "W4"),
                ),
            ),
        ),
    ),
)

#: A **latent** GP: Poisson counts whose rate is modulated by ``exp(f)``, with
#: ``f`` a GP draw reached through the whitened block the dataset declares.
#:
#: Added at W2.14, and it closed a hole the battery could not see. Every other
#: GP shape here is a Gaussian family, whose covariance is marginalised in
#: closed form by the solver; the latent path is the *other* correlated story
#: — ``DEVELOPMENT_PLAN.md`` §4.4's flexible likelihood for non-Gaussian data,
#: and the one no gradient-free engine can run — and until W2.14 nothing on
#: the scoring path applied ``GPSolver.latent_transform`` at all, so the
#: kernel hyperparameters entered the likelihood nowhere and both backends
#: were left to refuse or to mirror a wrong oracle. ``SolverKind.DENSE``
#: because every backend declares it; the quasiseparable latent path is
#: exercised in each backend's own suite.
LATENT_GP = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(
            family="poisson",
            noise=NoiseKind.GP,
            covariance=CovarianceSpec(amplitude=0.3, length_scale=2.0),
        ),
    ),
)

#: A **complex** dataset: circular complex Gaussian visibilities.
#:
#: Added at W2.4 slice 3. Every other shape here is real, so a backend's
#: realised path could — and on torch did — carry ``dtype=float`` from the
#: observed container to the residual and refuse the ``complex_gaussian``
#: family by name rather than admit it. The family is one of the five
#: ``ampere.core`` implements and the one ``results_schema.md`` §16's
#: :class:`~ampere.core.VisibilitySet` exists for, so a realisation that cannot
#: score it is a realisation the Phase-4 modality cannot be built on.
#:
#: The model is :attr:`ModelKind.COMPLEX`, whose modulus is ``norm`` and whose
#: phase winds with ``index``; the noise is i.i.d., because the correlated case
#: is declared analytic in ``likelihoods.md`` §4 and **not implemented in
#: ``ampere.core``** (``GP_ANALYTIC_IMPLEMENTED`` is ``False``), so there is no
#: oracle for a realisation to be checked against and every backend refuses it
#: by name. Skipped by a backend that declares no complex model, exactly as the
#: quasiseparable shape is skipped by one that declares no ``QUASISEP`` solver.
COMPLEX = ProblemSpec(
    model=ModelSpec(kind=ModelKind.COMPLEX, coordinates=GP_GRID),
    datasets=(DatasetSpec(family="complex_gaussian"),),
)

#: The same shape under the **circular complex GP** (*W4.2*). ``likelihoods.md``
#: §4 declares ``complex_gaussian`` + ``GaussianProcessNoise`` analytic with the
#: circular model as its fixed meaning, and W4.2 implemented the closed form, so
#: this is the shape that was *refused* on every path until Phase 4 and is now a
#: realisation-agreement row like the others.
#:
#: The kernel selects ``("u", "v")``, and it has to: the kind's three axes are in
#: mixed units, so W4.5's per-leaf unit rule refuses a bare one. That makes this
#: shape the battery's only realisation row whose native path must resolve an
#: ``axes=`` selector **and** stack a two-column right-hand side, which is
#: exactly the pair of things a backend can get wrong silently.
COMPLEX_GP = ProblemSpec(
    model=ModelSpec(kind=ModelKind.COMPLEX, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(
            family="complex_gaussian",
            noise=NoiseKind.GP,
            covariance=CovarianceSpec(KernelFamily.MATERN32, 0.3, 2.0, axes=("u", "v")),
        ),
    ),
)

#: The three shapes W3.14's generative rows need beside :data:`LATENT_GP` and
#: :data:`COMPLEX`, which already existed.
#:
#: ``POISSON_COUNTS`` is the analytic half of the Poisson family: counts, no
#: uncertainties (``REQUIRES_UNCERTAINTY`` is ``False``), independent noise.
#: ``STUDENT_T`` carries a deliberately large ``uncertainty`` for the reason
#: :data:`CORRELATED_DRAWS` does — a scale of 0.1 is small enough that a draw
#: at the wrong ``nu`` would still sit inside any honest tolerance. ``CAUCHY``
#: is the still-refusing family, and exists so that §13's refusal row keeps a
#: subject now that ``poisson`` has acquired a ``sample``.
POISSON_COUNTS = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(DatasetSpec(family="poisson"),),
)

STUDENT_T = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(DatasetSpec(family="student_t", uncertainty=0.3),),
)

CAUCHY = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(DatasetSpec(family="cauchy"),),
)

#: Readable pytest ids for the three shapes above, used where a row is
#: parametrised over all of them. The two shapes each slice 2 added carry
#: their own ids at the one row that takes all five.
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

    @pytest.mark.parametrize(
        "spec",
        [SINGLE, JOINT, CORRELATED, CORRELATED_QUASISEP, PHOTOMETRIC],
        ids=(*SPEC_IDS, "correlated-quasisep", "photometric"),
    )
    def test_the_realised_density_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances, spec: ProblemSpec
    ) -> None:
        """Many points, every declared shape, against the oracle."""
        needed = {dataset.solver for dataset in spec.datasets if dataset.noise is NoiseKind.GP}
        if not needed <= backend.capabilities.solvers:
            pytest.skip(
                f"backend {backend.name!r} declares no "
                f"{sorted(kind.name for kind in needed - backend.capabilities.solvers)} solver"
            )
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

    def test_the_complex_realisation_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """W2.4 slice 3: ``complex_gaussian``, realised, against the oracle.

        Its own row rather than a sixth entry in the parametrised list, because
        the *precondition* is different: this shape needs a complex-valued
        model, which is a fixture capability rather than a solver one.

        The check is the usual one — the realised density against the numpy
        contract path at every interior point, at
        ``tolerances.cross_backend`` — and it is worth stating what a failure
        would mean here specifically. ``ModelKind.COMPLEX`` has a constant
        modulus, so ``index`` enters the likelihood **only** through the phase:
        a backend that dropped the imaginary part anywhere between the model's
        ``flux`` and the residual would score a density independent of one of
        its own two parameters, and would disagree with this oracle by whole
        nats rather than in the last digit.
        """
        if not backend.capabilities.complex_models:
            pytest.skip(
                f"backend {backend.name!r} declares no complex model "
                f"(BackendCapabilities.complex_models), so ModelKind.COMPLEX cannot be built"
            )
        problem = build_problem(backend, COMPLEX)
        observed = problem.datasets["sed"].observed
        assert np.asarray(observed.values).dtype.kind == "c"
        realised = self.realised(problem)
        compared = 0
        for y in self.points(problem):
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            if not np.isfinite(expected):
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)
            compared += 1
        assert compared > 0, "every point was outside the support; the row proved nothing"

    def test_the_circular_complex_gp_realisation_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """**W4.2**: the circular complex GP, realised, against the oracle.

        Until Phase 4 this row was a *refusal* — ``GP_ANALYTIC_IMPLEMENTED``
        was ``False``, so the contract path refused the composition and no
        backend could be quietly computing something for it. W4.2 wrote the
        closed form, so the refusal becomes an agreement row, and it is the one
        realisation row that exercises two things no other does: the native path
        resolving an ``axes=("u", "v")`` selector against a three-axis
        container, and the native path stacking a complex residual into the two
        real columns the solver's right-hand side now carries.

        Both are silent failures if wrong. A native path reading only the first
        axis would score a kernel over ``u`` alone; one taking only the real
        column would score half the data. Either disagrees with this oracle by
        whole nats.
        """
        if not backend.capabilities.complex_models:
            pytest.skip(
                f"backend {backend.name!r} declares no complex model "
                f"(BackendCapabilities.complex_models), so ModelKind.COMPLEX cannot be built"
            )
        problem = build_problem(backend, COMPLEX_GP)
        dataset = problem.datasets["sed"]
        assert np.asarray(dataset.observed.values).dtype.kind == "c"
        # Analytic, not latent: no whitened block is declared for this pair.
        assert dataset.latent is None
        assert dataset.likelihood.marginalisation.value == "analytic"

        realised = self.realised(problem)
        compared = 0
        for y in self.points(problem):
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            if not np.isfinite(expected):
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)
            compared += 1
        assert compared > 0, "every point was outside the support; the row proved nothing"

    def test_the_circular_complex_gp_is_refused_on_the_o_n_path(
        self, backend: ConformanceBackend
    ) -> None:
        """The other half of W4.2's declaration: the O(N) solver cannot carry it.

        A circular complex GP hands the solver two columns, and ``QuasisepGP``
        declares ``STACKED_RESIDUALS = False`` — not as a gap but because
        ``REQUIRES_ORDERED_1D`` cannot hold for a point of the (u, v) plane at a
        wavelength. Refused at composition on the contract path, so no backend
        gets the chance to invent a lowering for it.

        Runs wherever a quasiseparable solver exists; it needs a container
        rather than a model, so it does not wait on ``complex_models``.
        """
        from ampere.core import ComplexGaussianFamily
        from ampere.core.exceptions import LikelihoodError

        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} supplies no quasiseparable solver")
        observed = observed_container(COMPLEX_GP, COMPLEX_GP.datasets[0])
        noise = backend.gp_noise(
            backend.kernel(COMPLEX_GP.datasets[0].covariance),
            backend.gp_solver(SolverKind.QUASISEP),
        )
        with pytest.raises(LikelihoodError, match="STACKED_RESIDUALS = False"):
            noise.check_compatible(ComplexGaussianFamily(), observed)

    def test_the_latent_gp_realisation_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """W2.14: the latent path's realised density, against the fixed oracle.

        Its own row rather than a sixth entry in the parametrised list above,
        because the interesting points are different. :meth:`points` includes
        two vectors a long way out, which is right for a bijected
        hyperparameter and wrong for a latent block: the block is ``Identity``
        bijected, so ``y = ±8`` is a whitened draw eight standard deviations
        from its prior mean at *every* one of its elements, and ``rate *
        exp(f)`` there is a number about which the two paths can only agree
        that it is not finite. The interior grid is the whole of what this row
        needs, and it is where a sampler lives.
        """
        problem = build_problem(backend, LATENT_GP)
        assert problem.datasets["sed"].latent is not None
        assert problem.free_size > len(GP_GRID)
        realised = self.realised(problem)
        compared = 0
        for y in self.points(problem)[:-2]:
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            if not np.isfinite(expected):
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)
            compared += 1
        assert compared > 0, "every point was outside the support; the row proved nothing"

    def test_the_latent_numpy_path_moves_with_the_kernel_hyperparameters(
        self, backend: ConformanceBackend
    ) -> None:
        """The claim the row above is only worth making because of.

        A realisation is checked against the numpy path, so an oracle that is
        flat in the GP hyperparameters makes a *perfect* conformance score out
        of a likelihood that cannot fit them. That is exactly what happened
        before W2.14, on every backend at once. This row holds the oracle
        itself: with the whitened block fixed, moving the amplitude and then
        the length scale must move the log-likelihood.

        It runs on **every** backend, the reference one included, because the
        noise model and the solver are the backend's own since W2.13 — so this
        is a statement about each backend's ``latent_transform``, not only
        about ``ampere.core``'s.
        """
        # A fixed, non-zero whitened block. At z = 0 the latent path is flat
        # in the hyperparameters for the honest reason that f = L(theta) 0 is
        # zero whatever L is, so a row evaluated at the prior median would
        # pass against the defect it exists to catch.
        whitened = np.linspace(-1.2, 1.2, len(GP_GRID))
        declared = LATENT_GP.datasets[0]

        def scored(amplitude: float, length_scale: float) -> float:
            spec = dataclasses.replace(
                LATENT_GP,
                datasets=(
                    dataclasses.replace(
                        declared,
                        covariance=CovarianceSpec(
                            family=declared.covariance.family,
                            amplitude=amplitude,
                            length_scale=length_scale,
                        ),
                    ),
                ),
            )
            problem = build_problem(backend, spec)
            latent = problem.datasets["sed"].latent
            assert latent is not None
            theta = dict(problem.reference_values)
            name = next(key for key in theta if key.endswith(f".{latent.parameter.name}"))
            theta[name] = whitened
            return problem.log_likelihood(theta)

        base = declared.covariance
        assert np.isfinite(scored(base.amplitude, base.length_scale))
        amplitudes = [scored(a, base.length_scale) for a in (0.1, 0.3, 1.2)]
        lengths = [scored(base.amplitude, ell) for ell in (0.5, 2.0, 8.0)]
        for name, values in (("amplitude", amplitudes), ("length_scale", lengths)):
            assert len(set(values)) == 3, (
                f"the latent log-likelihood is flat in {name!r}: {values}. The whitening "
                f"transform f = L(theta) z is not reaching the family (W2.14)."
            )

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


class TestTheGPNoiseFloor:
    """W3.1 slice 2: every backend's ``GaussianProcessNoise`` takes ``jitter=``.

    Three claims, and the first is the one that was actually broken: the
    composition **exists** on every backend. The other two are what make the
    keyword mean the same thing everywhere — the parameter is declared and held
    fixed (so no engine dimension appears), and the floor really is in the
    density, in quadrature with the data's own uncertainties, rather than
    accepted and dropped.
    """

    FLOOR = 0.05

    def test_the_floor_is_a_declared_fixed_parameter(self, backend: ConformanceBackend) -> None:
        problem = build_problem(backend, CORRELATED_JITTER)
        noise = problem.datasets["sed"].likelihood.noise
        assert "jitter" in noise.parameters.names
        assert "jitter" not in noise.parameters.free_names
        assert problem.free_size == build_problem(backend, CORRELATED).free_size

    def test_the_floor_reaches_the_density(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Against ``scipy``, not against the no-floor problem: an oracle, not a difference."""
        problem = build_problem(backend, CORRELATED_JITTER)
        theta = problem.unconstrain(problem.reference_values)
        spec = CORRELATED_JITTER.datasets[0]
        grid = np.asarray(GP_GRID, dtype=float)
        mean = analytic_flux(
            CORRELATED_JITTER.model, model_context(problem, dict(problem.reference_values)), grid
        )
        covariance = kernel_matrix(
            spec.covariance.family, grid, spec.covariance.amplitude, spec.covariance.length_scale
        ) + np.diag(np.full(grid.size, spec.uncertainty**2 + self.FLOOR**2))
        observed = np.asarray(problem.datasets["sed"].observed.values, dtype=float)
        expected = float(
            st.multivariate_normal(mean=mean, cov=covariance, allow_singular=False).logpdf(observed)
        )
        assert float(problem.log_likelihood(problem.constrain(theta))) == pytest.approx(
            expected, abs=tolerances.linear_algebra
        )

    def test_it_lowers_where_the_backend_has_a_realisation(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_problem(backend, CORRELATED_JITTER)
        if problem.backend not in registered_realisations():
            pytest.skip(f"the {problem.backend!r} backend registers no realisation")
        realised = realise(problem)
        theta = problem.unconstrain(problem.reference_values)
        assert float(np.asarray(realised.log_prob_unconstrained(theta))) == pytest.approx(
            float(problem.log_prob_unconstrained(theta)), abs=tolerances.cross_backend
        )


class TestNativeObservationSampling:
    """W3.1 slice 2: every backend draws observations, and the numpy path is the oracle.

    Peter's ruling of 2026-09-08. The comparison is **distributional**, not
    draw-for-draw, and that is forced rather than chosen: ``jax.random`` and
    ``torch.Generator`` do not reproduce numpy's stream, and making one do so
    would mean reimplementing a random library inside another. So what is
    asserted is what a sampling distribution *is* — its mean and its covariance
    — against the closed form, with the same two extra assertions the numpy
    oracle in ``tests/core/test_dataset.py`` carries and for the same reason:
    the joint tolerance alone would not catch a dropped ``K`` or a dropped
    ``sigma``, and both of those look entirely plausible in a plot.

    What *is* asserted exactly is the **refusals**. A family ``ampere.core``
    declines to sample is a family no backend may sample, because a backend
    guessing an observation process the contract will not guess is precisely
    how an SBI posterior gets trained on the wrong forward model.
    """

    DRAWS = 4000
    CHUNK = 500

    @staticmethod
    def natively_drawn(problem: FittingProblem, draws: int, chunk: int) -> Any:
        """*draws* observations at the reference θ, all on the native path."""
        theta = np.tile(
            np.asarray(problem.parameters.pack(dict(problem.reference_values)), dtype=float),
            (draws, 1),
        )
        return problem.simulate_many(
            draws, values=theta, observe=True, native=True, chunk_size=chunk
        )

    def problem(self, backend: ConformanceBackend, spec: ProblemSpec) -> FittingProblem:
        problem = build_problem(backend, spec)
        if problem.backend not in registered_realisations() or not problem.batchable:
            pytest.skip(f"the {problem.backend!r} backend has no native batched path")
        return problem

    def test_the_native_draw_has_the_declared_covariance(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``K + diag(σ²)``, from the backend's own random stream.

        The tolerance is the Monte-Carlo standard error of the sample
        covariance itself — ``se(Ĉᵢⱼ) = sqrt((Cᵢᵢ Cⱼⱼ + Cᵢⱼ²)/M)`` — so it stays
        honest as the draw count changes, exactly as the numpy row's does.
        """
        problem = self.problem(backend, CORRELATED_DRAWS)
        batch = self.natively_drawn(problem, self.DRAWS, self.CHUNK)
        assert batch.provenance["sample_backend"] == problem.backend

        spec = CORRELATED_DRAWS.datasets[0]
        grid = np.asarray(GP_GRID, dtype=float)
        covariance = kernel_matrix(
            spec.covariance.family, grid, spec.covariance.amplitude, spec.covariance.length_scale
        )
        expected = covariance + np.diag(np.full(grid.size, spec.uncertainty**2))
        drawn = np.asarray(batch.observations["sed"].values)
        empirical = np.cov(drawn, rowvar=False)

        diagonal = np.diag(expected)
        standard_error = np.sqrt((np.outer(diagonal, diagonal) + expected**2) / self.DRAWS)
        assert np.all(
            np.abs(empirical - expected) <= tolerances.monte_carlo_sigmas * standard_error
        )
        # The two ways this could be wrong and still look plausible, ruled out
        # explicitly rather than left to the joint tolerance: no correlation at
        # all (the off-diagonals collapse to zero) and no diagonal noise at all
        # (the variances are K's alone). The nearest-neighbour off-diagonals
        # are the ones with signal on this grid; the far corners are near zero
        # by construction at this length scale, so asserting on them would be
        # asserting on noise.
        assert np.abs(np.diag(empirical, 1)).min() > 0.05
        assert np.abs(np.diag(empirical) - np.diag(covariance)).min() > 0.04

    def test_the_native_draw_is_centred_on_the_prediction(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = self.problem(backend, CORRELATED_DRAWS)
        batch = self.natively_drawn(problem, self.DRAWS, self.CHUNK)
        drawn = np.asarray(batch.observations["sed"].values)
        mean = np.asarray(batch.predicted["sed"].values[0])
        spec = CORRELATED_DRAWS.datasets[0]
        grid = np.asarray(GP_GRID, dtype=float)
        variance = (
            np.diag(
                kernel_matrix(
                    spec.covariance.family,
                    grid,
                    spec.covariance.amplitude,
                    spec.covariance.length_scale,
                )
            )
            + spec.uncertainty**2
        )
        standard_error = np.sqrt(variance / self.DRAWS)
        assert np.all(
            np.abs(drawn.mean(axis=0) - mean) <= tolerances.monte_carlo_sigmas * standard_error
        )

    def test_an_uncorrelated_native_draw_has_the_declared_variance(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The other branch: ``x = mu + sigma z``, with no covariance to hide in."""
        problem = self.problem(backend, SINGLE)
        batch = self.natively_drawn(problem, self.DRAWS, self.CHUNK)
        drawn = np.asarray(batch.observations["sed"].values)
        sigma = np.asarray(problem.datasets["sed"].observed.uncertainty, dtype=float)
        empirical = drawn.var(axis=0)
        standard_error = sigma**2 * np.sqrt(2.0 / self.DRAWS)
        assert np.all(
            np.abs(empirical - sigma**2) <= tolerances.monte_carlo_sigmas * standard_error
        )
        off_diagonal = np.abs(np.cov(drawn, rowvar=False)[np.triu_indices(len(sigma), k=1)])
        assert off_diagonal.max() < 0.5 * float(sigma.min() ** 2)

    def test_a_refusing_family_refuses_by_name_on_every_backend(
        self, backend: ConformanceBackend
    ) -> None:
        """§13's refusal is the same sentence whichever backend is underneath.

        Asserted against the **text**, not against the exception type, because
        the ruling is that a backend does not get its own paraphrase: what the
        native path may not sample it hands back to the numpy path, and the
        numpy path names the family and the override that would supply the
        observation process.

        The family carrying this row was ``poisson`` until W3.14 moved it (and
        ``student_t`` and ``complex_gaussian``) to the other side of the line;
        it is ``cauchy`` now, which is the shipped family the core still
        declines to guess an observation process for.
        """
        from ampere.core.exceptions import DatasetError

        problem = build_problem(backend, CAUCHY)
        with pytest.raises(DatasetError) as raised:
            problem.simulate_many(2, observe=True)
        message = str(raised.value)
        assert "the cauchy family does not implement sample()" in message
        assert "override" in message and "sample(predicted, noise, rng)" in message

    def test_the_prediction_is_still_native_when_the_draw_is_not(
        self, backend: ConformanceBackend
    ) -> None:
        """Falling back on the noise does not throw away the vectorised forward model."""
        problem = build_problem(backend, CAUCHY)
        native = problem.backend in registered_realisations() and problem.batchable
        batch = problem.simulate_many(3, observe=False)
        assert batch.provenance["simulate_batched"] is native
        assert "sample_backend" not in batch.provenance


class TestTheUnambiguousFamilyDraws:
    """W3.14: ``poisson``, ``student_t`` and ``complex_gaussian`` draw, everywhere.

    One column per registered backend fixture, and **two paths per family**:
    the numpy oracle (``native=False``, which is ``ampere.core``'s own
    ``LikelihoodFamily.sample``) and the backend's native twin, held to the
    *same* assertions and the same tolerances the Gaussian rows above use.
    That is what "the numpy path is the oracle and the twins are compared
    distributionally" means operationally: not that one is checked against the
    other draw for draw — ``jax.random`` and ``torch.Generator`` reproduce
    neither numpy's stream nor each other's — but that both are checked
    against the moments their own ``log_prob`` implies, so a twin that had
    drifted from the density it is the counterpart of fails on its own.

    Each family is checked at whatever statistic actually pins its shape,
    rather than at a uniform mean-and-variance:

    * ``poisson`` — mean *and* variance, both the rate, plus integrality. A
      rounded Gaussian would pass a mean-only row and be a different forward
      model.
    * ``student_t`` — the scale through the **MAD** and the shape through a KS
      test against ``scipy.stats.t``. Not the variance: at ``nu = 4`` (the
      family's default, and what the battery composes) the fourth moment is
      infinite, so the sample variance has no standard error to compare
      against and a variance row would be a coin flip.
    * ``complex_gaussian`` — each component's variance against ``sigma**2``,
      because ``sigma`` is the *per-component* standard deviation, and the
      cross-covariance against zero, because circularity is the whole content
      of the family.
    """

    DRAWS = 4000
    CHUNK = 500

    #: A flat, generous rate: ``index = 0`` makes every sample's expectation
    #: ``norm``, so the Poisson rows have real power at every point rather than
    #: testing the ``k = 0`` term at the long-wavelength end of a falling power
    #: law (``norm``'s prior is ``LogUniform(0.1, 10)``, so 9.0 is inside it).
    POISSON_THETA: ClassVar[dict[str, float]] = {"model.norm": 9.0, "model.index": 0.0}

    def problem(
        self, backend: ConformanceBackend, spec: ProblemSpec, *, native: bool
    ) -> FittingProblem:
        problem = build_problem(backend, spec)
        if native and (problem.backend not in registered_realisations() or not problem.batchable):
            pytest.skip(f"the {problem.backend!r} backend has no native batched path")
        return problem

    def drawn(
        self,
        problem: FittingProblem,
        *,
        native: bool,
        overrides: Mapping[str, Any] | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """``(observations, prediction)`` for :attr:`DRAWS` draws at one θ."""
        values = dict(problem.reference_values)
        values.update(overrides or {})
        theta = np.tile(np.asarray(problem.parameters.pack(values), dtype=float), (self.DRAWS, 1))
        batch = problem.simulate_many(
            self.DRAWS, values=theta, observe=True, native=native, chunk_size=self.CHUNK
        )
        expected = problem.backend if native else "reference"
        assert batch.provenance["sample_backend"] == expected
        assert batch.observations is not None
        return (
            np.asarray(batch.observations["sed"].values),
            np.asarray(batch.predicted["sed"].values[0]),
        )

    @pytest.mark.parametrize("native", [False, True], ids=["numpy", "native"])
    def test_poisson_counts_have_the_rate_as_both_moments(
        self, backend: ConformanceBackend, tolerances: Tolerances, native: bool
    ) -> None:
        problem = self.problem(backend, POISSON_COUNTS, native=native)
        drawn, rate = self.drawn(problem, native=native, overrides=self.POISSON_THETA)
        assert np.all(drawn >= 0.0)
        assert np.all(drawn == np.round(drawn))
        mean_error = np.sqrt(rate / self.DRAWS)
        assert np.all(
            np.abs(drawn.mean(axis=0) - rate) <= tolerances.monte_carlo_sigmas * mean_error
        )
        # A Poisson's fourth central moment is ``lam + 3 lam**2``, so the
        # sample variance's own standard error is ``sqrt((lam + 2 lam**2)/M)``.
        variance_error = np.sqrt((rate + 2.0 * rate**2) / self.DRAWS)
        assert np.all(
            np.abs(drawn.var(axis=0) - rate) <= tolerances.monte_carlo_sigmas * variance_error
        )

    @pytest.mark.parametrize("native", [False, True], ids=["numpy", "native"])
    def test_a_latent_gp_poisson_draw_uses_the_latent_in_theta(
        self, backend: ConformanceBackend, tolerances: Tolerances, native: bool
    ) -> None:
        """The rate is ``predicted * exp(f)``, at the ``f`` the density scores at.

        The property SBC needs: θ carries the whitened ``z``, and if the draw
        invented its own the pair ``(θ, x)`` would come from no model at all.
        Driven at a *constant* ``z`` so the expected rate is a computed number
        rather than a distribution, and the transform is the solver's own —
        ``f = L(θ) z`` is not ``z``, which is exactly why the family may not
        read the whitened block itself.
        """
        problem = self.problem(backend, LATENT_GP, native=native)
        latent = next(name for name in problem.parameters.names if name.startswith("sed.latent."))
        grid = np.asarray(GP_GRID, dtype=float)
        whitened = np.full(grid.size, 0.8)
        drawn, rate = self.drawn(problem, native=native, overrides={latent: whitened})
        spec = LATENT_GP.datasets[0]
        covariance = kernel_matrix(
            spec.covariance.family, grid, spec.covariance.amplitude, spec.covariance.length_scale
        )
        expected = rate * np.exp(np.linalg.cholesky(covariance) @ whitened)
        assert np.all(drawn == np.round(drawn))
        error = np.sqrt(expected / self.DRAWS)
        assert np.all(
            np.abs(drawn.mean(axis=0) - expected) <= tolerances.monte_carlo_sigmas * error
        )

    @pytest.mark.parametrize("native", [False, True], ids=["numpy", "native"])
    def test_student_t_draws_recover_the_scale_and_the_degrees_of_freedom(
        self, backend: ConformanceBackend, tolerances: Tolerances, native: bool
    ) -> None:
        problem = self.problem(backend, STUDENT_T, native=native)
        drawn, location = self.drawn(problem, native=native)
        sigma = float(STUDENT_T.datasets[0].uncertainty)
        nu = float(family_named("student_t")().parameters["nu"].value)
        standardised = (drawn - location) / sigma
        # The MAD, not the variance: at nu = 4 the fourth moment is infinite,
        # so a sample variance has no standard error to be compared against.
        # median|t| is t.ppf(0.75), and the median's own standard error is
        # 1 / (2 f(m) sqrt(M)).
        median = float(st.t(nu).ppf(0.75))
        error = 1.0 / (2.0 * float(st.t(nu).pdf(median)) * np.sqrt(self.DRAWS))
        empirical = np.median(np.abs(standardised), axis=0)
        assert np.all(np.abs(empirical - median) <= tolerances.monte_carlo_sigmas * error)
        # ... and the shape, which no scale statistic can see.
        assert st.kstest(standardised[:, 0], st.t(nu).cdf).pvalue > 1e-4
        matched = st.norm(0.0, np.std(standardised[:, 0]))
        assert st.kstest(standardised[:, 0], matched.cdf).pvalue < 1e-4

    @pytest.mark.parametrize("native", [False, True], ids=["numpy", "native"])
    def test_the_complex_gaussian_draw_is_circular(
        self, backend: ConformanceBackend, tolerances: Tolerances, native: bool
    ) -> None:
        if not backend.capabilities.complex_models:
            pytest.skip(
                f"backend {backend.name!r} declares no complex model "
                f"(BackendCapabilities.complex_models), so ModelKind.COMPLEX cannot be built"
            )
        problem = self.problem(backend, COMPLEX, native=native)
        drawn, mean = self.drawn(problem, native=native)
        assert np.iscomplexobj(drawn)
        sigma = float(COMPLEX.datasets[0].uncertainty)
        variance = sigma**2
        standard_error = variance * np.sqrt(2.0 / self.DRAWS)
        for component in (drawn.real, drawn.imag):
            assert np.all(
                np.abs(component.var(axis=0) - variance)
                <= tolerances.monte_carlo_sigmas * standard_error
            )
        centre_error = sigma / np.sqrt(self.DRAWS)
        assert np.all(
            np.abs(drawn.mean(axis=0) - mean) <= tolerances.monte_carlo_sigmas * centre_error
        )
        # Circularity: the components are independent, so their covariance is
        # zero and its own standard error is ``sigma**2 / sqrt(M)``.
        cross = np.mean(
            (drawn.real - drawn.real.mean(axis=0)) * (drawn.imag - drawn.imag.mean(axis=0)),
            axis=0,
        )
        assert np.all(
            np.abs(cross) <= tolerances.monte_carlo_sigmas * variance / np.sqrt(self.DRAWS)
        )


class TestANativeSamplerFailureIsolatesOneDraw:
    """W4.0 (2): the native sampler's own whole-chunk guard does not abort the chunk.

    A native sampler is vectorised over the whole chunk it is handed — torch's
    Poisson twin (``_poisson_variates``) checks the *stacked* rate array for a
    non-positive entry in one call, so a single bad θ's exception carries no
    row index and would, unfixed, discard or raise for every draw in the
    chunk (W3.14's finding). ``FittingProblem.simulate_many(native=True)``
    retries one row at a time after such a failure, so only the request that
    earned it is flagged — exactly as the numpy loop's own per-draw ``try``
    around ``draw_observation`` isolates a likelihood failure — and the rest
    of the chunk keeps its native draw rather than falling back to the loop.

    ``poisson`` is the first family that can reach this path (*W3.14*): it is
    the only one of the three native-sampled families whose native twin
    refuses a whole batch for one bad row rather than drawing something for
    every θ regardless of sign.
    """

    def problem(self, backend: ConformanceBackend) -> FittingProblem:
        problem = build_problem(backend, POISSON_COUNTS)
        if problem.backend not in registered_realisations() or not problem.batchable:
            pytest.skip(f"the {problem.backend!r} backend has no native batched path")
        return problem

    def theta(self, problem: FittingProblem, overrides: Mapping[str, Any]) -> np.ndarray:
        values = dict(problem.reference_values)
        values.update(overrides)
        return np.asarray(problem.parameters.pack(values), dtype=float)

    def test_one_negative_rate_among_several_flags_only_itself(
        self, backend: ConformanceBackend
    ) -> None:
        problem = self.problem(backend)
        good = TestTheUnambiguousFamilyDraws.POISSON_THETA
        broken = {"model.norm": -9.0, "model.index": 0.0}
        theta = np.stack(
            [
                self.theta(problem, good),
                self.theta(problem, good),
                self.theta(problem, broken),
                self.theta(problem, good),
            ]
        )
        # One chunk (the default), so the batched sampler call really does
        # cover all four draws at once and the retry path is exercised.
        batch = problem.simulate_many(4, values=theta, observe=True, native=True)
        assert batch.provenance["simulate_batched"] is True
        assert batch.provenance["sample_backend"] == problem.backend
        assert list(np.asarray(batch.failed)) == [False, False, True, False]
        offender = batch[2]
        assert offender.failure is not None
        assert offender.failure.reason == FailureReason.LIKELIHOOD_FAILED
        for index in (0, 1, 3):
            simulation = batch[index]
            assert not simulation.failed
            drawn = np.asarray(simulation.observations["sed"].values)
            assert np.all(drawn >= 0.0) and np.all(drawn == np.round(drawn))

    def test_the_offending_draw_still_raises_under_strict(
        self, backend: ConformanceBackend
    ) -> None:
        """``strict=True`` gets the raise at the offending draw, not a retried one.

        The same rule ``_failure_types`` states for the loop (``likelihoods.md``
        §17 Q1's ruling): non-strict isolates and records, strict propagates.
        """
        strict = build_problem(backend, POISSON_COUNTS, strict=True)
        if strict.backend not in registered_realisations() or not strict.batchable:
            pytest.skip(f"the {strict.backend!r} backend has no native batched path")
        good = TestTheUnambiguousFamilyDraws.POISSON_THETA
        broken = {"model.norm": -9.0, "model.index": 0.0}
        theta = np.stack([self.theta(strict, good), self.theta(strict, broken)])
        with pytest.raises(Exception):  # noqa: B017 - the backend's own exception type
            strict.simulate_many(2, values=theta, observe=True, native=True)


class TestBatchedSimulation:
    """W3.1: ``simulate_many`` is the loop, on every backend and every executor.

    ``inference.md`` §13's *batched form* is a statement about equality with a
    sequence of ``simulate`` calls, so that is what these rows assert — order
    and values, not shapes and plausibility. The batch's right-hand side is
    built on a *second* problem with the same seed, because reading
    ``rng("simulate")`` in the test would advance the very stream the batch is
    about to spawn from.

    **Amended at W3.1 slice 2**, when the equality acquired two grades. The
    loop remains the semantics, and ``native=False`` still reproduces it
    **bitwise** — that row is the literal statement slice 1 landed. On a
    backend that can run the chunk through its own ``vmap``, the *prediction*
    agrees to ``tolerances.cross_backend`` rather than bitwise (vectorised
    arithmetic is not scalar arithmetic in the last digits) and the
    *observations* are a draw from the same distribution rather than the same
    draw, because ``jax.random`` and ``torch.Generator`` are not numpy's
    stream and could not be made to be without reimplementing one library
    inside another. Which of the two happened is not left to inference: the
    batch's ``provenance`` records it, and these rows read it.
    """

    COUNT = 6

    @staticmethod
    def loop(backend: ConformanceBackend, count: int, *, observe: bool) -> SimulationBatch:
        """``n`` calls of ``simulate`` on the batch sub-stream's spawned children."""
        problem = build_problem(backend, SINGLE)
        children = problem.rng("simulate").spawn(count)
        return SimulationBatch(
            tuple(problem.simulate(observe=observe, rng=child) for child in children)
        )

    @staticmethod
    def agrees(
        left: SimulationBatch,
        right: SimulationBatch,
        *,
        tolerance: float = 0.0,
        observations: bool = True,
    ) -> bool:
        """Order and values. *tolerance* 0.0 is the bitwise form slice 1 landed.

        ``observations=False`` compares the noise-free half only, which is what
        a natively drawn batch can promise: the same distribution, from a
        different random stream.
        """

        def same(one: Any, other: Any) -> bool:
            first, second = np.asarray(one), np.asarray(other)
            if tolerance == 0.0:
                return bool(np.array_equal(first, second))
            return bool(np.allclose(first, second, rtol=0.0, atol=tolerance))

        if len(left) != len(right) or not same(left.theta, right.theta):
            return False
        for one, other in zip(left, right, strict=True):
            for label in one.predicted:
                if not same(one.predicted[label].values, other.predicted[label].values):
                    return False
            if (one.observations is None) != (other.observations is None):
                return False
            if observations and one.observations is not None and other.observations is not None:
                for label in one.observations:
                    if not same(one.observations[label].values, other.observations[label].values):
                        return False
        return True

    @staticmethod
    def can_run_natively(problem: FittingProblem) -> bool:
        """Whether this backend offers the vectorised forward path for *problem*."""
        return bool(problem.batchable) and problem.backend in registered_realisations()

    @classmethod
    def matches_the_loop(
        cls,
        batch: SimulationBatch,
        reference: SimulationBatch,
        tolerances: Tolerances,
    ) -> bool:
        """Hold *batch* to whatever its own provenance says it is.

        The one place the two grades of the equality are chosen between, so no
        row has to remember which backend does what: a batch that ran the loop
        is held to the bitwise identity, and one that ran a backend's ``vmap``
        to ``cross_backend`` on the prediction — with the observations compared
        only if the loop drew them.
        """
        native = bool(batch.provenance.get("simulate_batched"))
        drew_natively = batch.provenance.get("sample_backend", "reference") != "reference"
        return cls.agrees(
            batch,
            reference,
            tolerance=tolerances.cross_backend if native else 0.0,
            observations=not drew_natively,
        )

    @pytest.mark.parametrize("observe", [False, True], ids=["predicted", "observed"])
    def test_the_batch_is_n_calls_of_simulate(
        self, backend: ConformanceBackend, tolerances: Tolerances, observe: bool
    ) -> None:
        batch = build_problem(backend, SINGLE).simulate_many(self.COUNT, observe=observe)
        loop = self.loop(backend, self.COUNT, observe=observe)
        assert self.matches_the_loop(batch, loop, tolerances)

    def test_the_loop_is_still_reproduced_bitwise(self, backend: ConformanceBackend) -> None:
        """``native=False`` is slice 1's statement, unamended, on every backend.

        Worth a row of its own rather than a note: the native path is a
        throughput optimisation *underneath* an equality, and the way to keep
        that true is to be able to switch it off and get the equality back
        exactly. If this ever fails, something in the native path has leaked
        into the semantics.
        """
        batch = build_problem(backend, SINGLE).simulate_many(self.COUNT, observe=True, native=False)
        assert self.agrees(batch, self.loop(backend, self.COUNT, observe=True))
        assert batch.provenance["simulate_batched"] is False
        assert batch.provenance["sample_backend"] == "reference"

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ThreadExecutor(2), ProcessExecutor(2)],
        ids=["serial", "thread-2", "process-2"],
    )
    def test_every_shipped_executor_gives_the_same_batch(
        self, backend: ConformanceBackend, tolerances: Tolerances, executor: Any
    ) -> None:
        from ampere.core.exceptions import DatasetError

        problem = build_problem(backend, SINGLE)
        if isinstance(executor, ProcessExecutor) and not backend.capabilities.picklable:
            # Not a skip: the contract for a backend whose arrays live on a
            # device is that the pool is refused *by name*, before any worker
            # starts, rather than failing somewhere inside one.
            with pytest.raises(DatasetError, match="cannot be sent to a worker process"):
                problem.simulate_many(self.COUNT, observe=True, executor=executor)
            return
        batch = problem.simulate_many(self.COUNT, observe=True, executor=executor, chunk_size=4)
        loop = self.loop(backend, self.COUNT, observe=True)
        assert self.matches_the_loop(batch, loop, tolerances)

    @pytest.mark.parametrize("chunk_size", [1, 7, None], ids=["chunk-1", "chunk-7", "whole"])
    def test_the_partition_does_not_change_the_answer(
        self, backend: ConformanceBackend, tolerances: Tolerances, chunk_size: int | None
    ) -> None:
        batch = build_problem(backend, SINGLE).simulate_many(
            self.COUNT, observe=True, chunk_size=chunk_size
        )
        loop = self.loop(backend, self.COUNT, observe=True)
        assert self.matches_the_loop(batch, loop, tolerances)

    @pytest.mark.parametrize("chunk_size", [1, 7, None], ids=["chunk-1", "chunk-7", "whole"])
    def test_the_partition_does_not_change_the_native_answer(
        self, backend: ConformanceBackend, tolerances: Tolerances, chunk_size: int | None
    ) -> None:
        """The same claim on the native path, and it is a stronger one there.

        A chunk is what gets ``vmap``ped, so on the native path the partition
        decides the *shape of every array in the trace* — which is exactly the
        thing that could change an answer without changing a line of the model.
        Compared against the whole-budget native batch rather than against the
        loop, so a difference here is the chunking and nothing else.
        Observations included: the per-draw key is derived by index, so native
        sampling is partition-independent too.
        """
        problem = build_problem(backend, SINGLE)
        if not self.can_run_natively(problem):
            pytest.skip(f"the {problem.backend!r} backend has no native batched path")
        whole = problem.simulate_many(self.COUNT, observe=True, native=True)
        chunked = build_problem(backend, SINGLE).simulate_many(
            self.COUNT, observe=True, native=True, chunk_size=chunk_size
        )
        assert self.agrees(chunked, whole, tolerance=tolerances.cross_backend)

    def test_the_provenance_says_which_path_ran(self, backend: ConformanceBackend) -> None:
        """A stored budget must be able to say how it was made, not be guessed at."""
        problem = build_problem(backend, SINGLE)
        batch = problem.simulate_many(self.COUNT, observe=True)
        assert batch.provenance["simulate_batched"] is self.can_run_natively(problem)
        assert batch.provenance["simulation_context"] == "none"
        expected = problem.backend if self.can_run_natively(problem) else "reference"
        assert batch.provenance["sample_backend"] == expected

    #: Draws in the context row: enough standardised residuals (x 12 samples)
    #: that a sigma which failed to reach the realisation -- off by up to a
    #: factor of ten either way under the prior below -- is dozens of Monte
    #: Carlo sigmas out, and a correct one is inside five.
    CONTEXT_COUNT = 200

    def test_a_context_budget_runs_natively_and_matches_the_loop(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """**W5.29**: a batched draw under a context prior is the loop's draw.

        The reference loop is the oracle (the reference backend has no batched
        branch, so the row is skipped there). What is equal **exactly**: θ,
        every draw's context record, and the sigma every observation carries
        -- the contexts are drawn in ``ampere.core`` on the
        ``"<stream>.context"`` sub-stream by index, whichever path then uses
        them. What is equal to ``cross_backend``: the prediction, as in every
        native row. What cannot be equal draw for draw is the noise, for the
        reason :class:`TestBatchedSimulation`'s docstring gives -- a backend's
        random stream is not numpy's -- so the noise is held to its
        distribution instead: each observation standardised by *its own*
        context sigma is ``N(0, 1)``, which is the claim that the per-draw
        sigma reached the realisation rather than the container alone.
        """
        from ampere.core.simulate import ScaledSigma

        problem = build_problem(backend, SINGLE)
        if not self.can_run_natively(problem):
            pytest.skip(f"the {problem.backend!r} backend has no native batched path")
        prior = ScaledSigma(0.1, 10.0)
        count = self.CONTEXT_COUNT
        batch = problem.simulate_many(count, observe=True, native=True, context=prior)
        loop = build_problem(backend, SINGLE).simulate_many(
            count, observe=True, native=False, context=prior
        )
        assert batch.provenance["simulate_batched"] is True
        assert batch.provenance["sample_backend"] == problem.backend
        assert batch.provenance["simulation_context"] == loop.provenance["simulation_context"]
        assert self.matches_the_loop(batch, loop, tolerances)
        residuals = []
        for one, other in zip(batch, loop, strict=True):
            assert one.context is not None and other.context is not None
            assert dict(one.context.record) == dict(other.context.record)
            assert one.observations is not None and other.observations is not None
            for label, drawn in one.observations.items():
                sigma = np.asarray(drawn.uncertainty)
                assert np.array_equal(sigma, np.asarray(other.observations[label].uncertainty))
                assert np.array_equal(sigma, one.context.sigma[label])
                mean = np.asarray(one.predicted[label].values)
                residuals.append((np.asarray(drawn.values) - mean) / sigma)
        standard = np.concatenate(residuals)
        size = standard.size
        assert abs(float(standard.mean())) <= tolerances.monte_carlo_sigmas / np.sqrt(size)
        # The variance of a unit normal's sample variance is 2/n.
        assert abs(float(standard.var()) - 1.0) <= tolerances.monte_carlo_sigmas * np.sqrt(
            2.0 / size
        )

    def test_a_problem_that_cannot_be_run_natively_is_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        """``native=True`` must be a sentence, on every backend, never a slow path.

        The same rule ``log_prob_unconstrained_batched`` follows: a capability
        that is not there is refused with a message naming what withdrew it,
        because a fast path that quietly becomes a slow one is a fast path
        nobody can measure.
        """
        from ampere.core.exceptions import DatasetError, LoweringError

        problem = build_problem(backend, SINGLE)
        if self.can_run_natively(problem):
            if SolverKind.QUASISEP not in backend.capabilities.solvers:
                pytest.skip(f"{problem.backend!r} runs every declared shape natively")
            problem = build_problem(backend, CORRELATED_QUASISEP)
            if problem.batchable:
                pytest.skip(f"{problem.backend!r} declares the quasiseparable solver batchable")
        with pytest.raises((DatasetError, LoweringError)):
            problem.simulate_many(2, native=True)

    def test_a_composed_problem_can_be_sent_to_a_worker(self, backend: ConformanceBackend) -> None:
        """The precondition for the process pool, asserted rather than assumed.

        A backend that declares ``picklable=False`` is held to the other half
        of the contract instead: it must genuinely *fail* to pickle, so that
        the refusal ``simulate_many`` raises is a true statement about the
        problem rather than a stale declaration on a fixture.
        """
        problem = build_problem(backend, SINGLE)
        if not backend.capabilities.picklable:
            with pytest.raises((TypeError, ValueError, AttributeError, pickle.PicklingError)):
                pickle.dumps(problem)
            return
        restored = pickle.loads(pickle.dumps(problem))
        assert restored.free_size == problem.free_size
        assert restored.backend == problem.backend

    def test_the_stacked_views_agree_with_the_per_draw_simulations(
        self, backend: ConformanceBackend
    ) -> None:
        batch = build_problem(backend, SINGLE).simulate_many(4, observe=True)
        assert batch.theta.shape == (4, build_problem(backend, SINGLE).free_size)
        assert batch.observations is not None
        for index, draw in enumerate(batch):
            assert np.array_equal(np.asarray(draw.theta), batch.theta[index])
            for label, stack in batch.predicted.items():
                assert np.array_equal(
                    np.asarray(draw.predicted[label].values), np.asarray(stack.values[index])
                )
                rebuilt = stack[index]
                assert rebuilt is not None
                assert rebuilt.axes == draw.predicted[label].axes

    def test_a_given_theta_table_is_simulated_row_by_row(self, backend: ConformanceBackend) -> None:
        """The SBC idiom: θ chosen outside, in the order they were given."""
        problem = build_problem(backend, SINGLE)
        table = np.stack(
            [problem.prior_transform(np.full(problem.free_size, q)) for q in (0.3, 0.7)]
        )
        batch = problem.simulate_many(2, values=table)
        assert np.allclose(batch.theta, table)
        for index in range(2):
            expected = problem.simulate(table[index])
            for label in expected.predicted:
                assert np.array_equal(
                    np.asarray(batch[index].predicted[label].values),
                    np.asarray(expected.predicted[label].values),
                )


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
