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
from typing import Any

import numpy as np
import pytest

from ampere.core import (
    Dataset,
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

    def test_a_complex_gp_is_refused_by_name_on_every_path(
        self, backend: ConformanceBackend
    ) -> None:
        """The other half of W2.4 slice 3's complex item, and it is a refusal.

        ``ComplexGaussianFamily`` declares ``ANALYTIC_WITH_GP = True`` — the
        circular complex GP does marginalise in closed form — and
        ``GP_ANALYTIC_IMPLEMENTED = False``, because ``ampere.core`` has not
        written it (``likelihoods.md`` §4; the implementation lands with the
        visibility modality in Phase 4). This row holds the invariant that
        makes that pair safe: the *contract* path refuses the composition, so
        no backend can be quietly computing something for it, and the refusal
        names the family.

        It runs on every fixture, complex model or not, because it composes
        nothing but a likelihood.
        """
        from ampere.core import ComplexGaussianFamily
        from ampere.core.exceptions import LikelihoodError

        noise = backend.gp_noise(
            backend.kernel(CovarianceSpec()), backend.gp_solver(SolverKind.DENSE)
        )
        with pytest.raises(LikelihoodError, match="complex_gaussian"):
            Likelihood(ComplexGaussianFamily(), noise)

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


class TestBatchedSimulation:
    """W3.1: ``simulate_many`` is the loop, on every backend and every executor.

    ``inference.md`` §13's *batched form* is a statement about equality with a
    sequence of ``simulate`` calls, so that is what these rows assert — order
    and values, not shapes and plausibility. The batch's right-hand side is
    built on a *second* problem with the same seed, because reading
    ``rng("simulate")`` in the test would advance the very stream the batch is
    about to spawn from.
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
    def identical(left: SimulationBatch, right: SimulationBatch) -> bool:
        if len(left) != len(right) or not np.array_equal(left.theta, right.theta):
            return False
        for one, other in zip(left, right, strict=True):
            for label in one.predicted:
                if not np.array_equal(
                    np.asarray(one.predicted[label].values),
                    np.asarray(other.predicted[label].values),
                ):
                    return False
            if (one.observations is None) != (other.observations is None):
                return False
            if one.observations is not None and other.observations is not None:
                for label in one.observations:
                    if not np.array_equal(
                        np.asarray(one.observations[label].values),
                        np.asarray(other.observations[label].values),
                    ):
                        return False
        return True

    @pytest.mark.parametrize("observe", [False, True], ids=["predicted", "observed"])
    def test_the_batch_is_n_calls_of_simulate(
        self, backend: ConformanceBackend, observe: bool
    ) -> None:
        batch = build_problem(backend, SINGLE).simulate_many(self.COUNT, observe=observe)
        assert self.identical(batch, self.loop(backend, self.COUNT, observe=observe))

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ThreadExecutor(2), ProcessExecutor(2)],
        ids=["serial", "thread-2", "process-2"],
    )
    def test_every_shipped_executor_gives_the_same_batch(
        self, backend: ConformanceBackend, executor: Any
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
        assert self.identical(batch, self.loop(backend, self.COUNT, observe=True))

    @pytest.mark.parametrize("chunk_size", [1, 7, None], ids=["chunk-1", "chunk-7", "whole"])
    def test_the_partition_does_not_change_the_answer(
        self, backend: ConformanceBackend, chunk_size: int | None
    ) -> None:
        batch = build_problem(backend, SINGLE).simulate_many(
            self.COUNT, observe=True, chunk_size=chunk_size
        )
        assert self.identical(batch, self.loop(backend, self.COUNT, observe=True))

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
