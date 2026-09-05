"""Rows that compare two backends rather than one backend against an oracle.

``DEVELOPMENT_PLAN.md`` §3 says torch and jax are implemented "in parallel
against a frozen interface spec by separate agent tracks", and §4.6 names this
suite as "the contract that keeps two lockstep backends aligned". These are the
rows that do that work directly: hand the *same* declaration to two backends and
require the same answer.

Every row here skips when only one backend is installed, and the module's whole
value multiplies as backends register. The two in-repo fixtures — the reference
path and its deliberately differently-routed mirror — make the rows live today
rather than notional, which is what stops the parametrisation rotting before
Phase 2 arrives.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Parameter, ParameterSet

from .composition import GP_GRID, DatasetSpec, NoiseKind, ProblemSpec, build_problem
from .conftest import cross_backend_pairs, pair_ids
from .protocol import ConformanceBackend, ModelKind, ModelSpec

pytestmark = pytest.mark.parametrize(("first", "second"), cross_backend_pairs(), ids=pair_ids())

if not cross_backend_pairs():  # pragma: no cover - only one backend installed
    pytest.skip(
        "cross-backend agreement needs at least two registered backends; "
        "see tests/conformance/backends/__init__.py",
        allow_module_level=True,
    )


SPECS = {
    "linear-iid": ProblemSpec(
        model=ModelSpec(kind=ModelKind.LINEAR, coordinates=GP_GRID),
        datasets=(DatasetSpec(),),
    ),
    "power-law-gp": ProblemSpec(
        model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
        datasets=(DatasetSpec(noise=NoiseKind.GP),),
    ),
    "masked": ProblemSpec(
        model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
        datasets=(DatasetSpec(masked=(1, 6)),),
    ),
}


def shared_tolerance(first: ConformanceBackend, second: ConformanceBackend, name: str) -> float:
    """The looser of two backends' tolerances for the named comparison class.

    A float32 or different-accumulation-order backend widens its own entry;
    taking the maximum means it does not silently tighten the other's, nor
    force everyone else to loosen.
    """
    return max(
        getattr(first.capabilities.tolerances, name),
        getattr(second.capabilities.tolerances, name),
    )


def sample_cube(size: int) -> np.ndarray:
    """A fixed, unremarkable point inside the unit cube."""
    return np.linspace(0.17, 0.83, size)


class TestParameterSpaceAgreement:
    """Two lowerings of one declaration are the same function."""

    @staticmethod
    def declaration() -> ParameterSet:
        return ParameterSet(
            [
                Parameter("bounded", st.uniform(0.0, 4.0)),
                Parameter("unbounded", st.norm(0.0, 1.0)),
                Parameter("positive", st.halfnorm(0.0, 2.0)),
            ]
        )

    def test_the_prior_transform_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        declaration = self.declaration()
        cube = sample_cube(declaration.free_size)
        assert first.parameter_space(declaration).prior_transform(cube) == pytest.approx(
            second.parameter_space(declaration).prior_transform(cube),
            abs=shared_tolerance(first, second, "cross_backend"),
        )

    def test_the_log_prior_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        declaration = self.declaration()
        theta = first.parameter_space(declaration).prior_transform(
            sample_cube(declaration.free_size)
        )
        assert first.parameter_space(declaration).lnprior(theta) == pytest.approx(
            second.parameter_space(declaration).lnprior(theta),
            abs=shared_tolerance(first, second, "cross_backend"),
        )

    def test_the_unconstrained_log_prior_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        declaration = self.declaration()
        y = np.array([0.4, -1.1, 0.6])
        assert first.parameter_space(declaration).lnprior_unconstrained(y) == pytest.approx(
            second.parameter_space(declaration).lnprior_unconstrained(y),
            abs=shared_tolerance(first, second, "cross_backend"),
        )

    def test_the_flat_layout_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        declaration = self.declaration()
        assert (
            first.parameter_space(declaration).free_labels()
            == second.parameter_space(declaration).free_labels()
        )


@pytest.mark.parametrize("name", sorted(SPECS))
class TestLogProbAgreement:
    """The headline row: the same problem scores the same on either backend."""

    def test_log_prob_agrees_at_the_prior_median(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        spec = SPECS[name]
        left, right = build_problem(first, spec), build_problem(second, spec)
        assert left.free_labels() == right.free_labels()
        assert left.log_prob() == pytest.approx(
            right.log_prob(), abs=shared_tolerance(first, second, "cross_backend")
        )

    def test_the_decomposition_agrees_away_from_the_median(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        spec = SPECS[name]
        left, right = build_problem(first, spec), build_problem(second, spec)
        theta = left.prior_transform(sample_cube(left.free_size))
        tolerance = shared_tolerance(first, second, "cross_backend")

        assert left.log_prior(theta) == pytest.approx(right.log_prior(theta), abs=tolerance)
        assert left.log_likelihood(theta) == pytest.approx(
            right.log_likelihood(theta), abs=tolerance
        )
        assert dict(left.evaluate(theta).contributions) == pytest.approx(
            dict(right.evaluate(theta).contributions), abs=tolerance
        )

    def test_the_prediction_itself_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        """Not only the score: the containers the models produce, sample by sample."""
        spec = SPECS[name]
        left, right = build_problem(first, spec), build_problem(second, spec)
        theta = left.prior_transform(sample_cube(left.free_size))
        for label in left.datasets:
            assert first.to_numpy(left.simulate(theta).predicted[label].values) == pytest.approx(
                second.to_numpy(right.simulate(theta).predicted[label].values),
                abs=shared_tolerance(first, second, "cross_backend"),
            )


class TestProvenanceAgreement:
    """The hashes are functions of the declaration and the data, not the arithmetic.

    ``results.md`` §14 states that claim and calls the hash "the cheapest
    possible detector" of a lowering bug: two backends emitting the same
    problem must agree, because a disagreement can only come from one of them
    having lowered the declaration differently.
    """

    @staticmethod
    def attrs(backend: ConformanceBackend, name: str) -> dict[str, object]:
        pytest.importorskip("arviz", reason="the results contract needs ampere[arviz]")
        from ampere.results import provenance_attrs

        return provenance_attrs(build_problem(backend, SPECS[name]))

    @pytest.mark.parametrize("name", sorted(SPECS))
    def test_the_spec_hash_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        left, right = self.attrs(first, name), self.attrs(second, name)
        assert left["ampere_spec_hash"] == right["ampere_spec_hash"]
        assert left["ampere_component_spec_hashes"] == right["ampere_component_spec_hashes"]

    @pytest.mark.parametrize("name", sorted(SPECS))
    def test_the_data_hash_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        left, right = self.attrs(first, name), self.attrs(second, name)
        assert left["ampere_data_hash"] == right["ampere_data_hash"]
        assert left["ampere_data_hashes"] == right["ampere_data_hashes"]

    @pytest.mark.parametrize("name", sorted(SPECS))
    def test_the_problem_hash_tracks_the_model_class_not_only_the_declaration(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        """A recorded contradiction between ``results.md`` §14 and ``provenance.py``.

        §14 promises that ``ampere_spec_hash``, ``ampere_data_hash`` **and**
        ``ampere_problem_hash`` "are functions of the declaration and the data,
        not of the arithmetic", so that any cross-backend disagreement is a
        lowering bug. The first two hold, and are asserted above. The third
        does not: ``provenance.model_fingerprint`` deliberately records the
        model's ``class`` and ``module`` — with a docstring arguing, correctly,
        that "two models of *different classes* can declare the same parameters
        and compute completely different things" — and ``_describe_instrument``
        records each step's class too. Two backends implementing one
        declaration therefore disagree by construction.

        Both behaviours are defensible on their own terms, and W1.10 has no
        standing to change either. So the row asserts the *actual* rule as an
        equivalence rather than xfailing the claim: the problem hash agrees
        exactly when the recorded model identity does. It is green today, it
        stays green for a backend that reuses another's model classes, and it
        fails loudly the moment either the contract or the fingerprint is
        changed — which is precisely when someone should be looking at it.
        """
        spec = SPECS[name]
        identities = tuple(
            (type(model).__name__, type(model).__module__)
            for model in (build_problem(first, spec).model, build_problem(second, spec).model)
        )
        left, right = self.attrs(first, name), self.attrs(second, name)
        agrees = left["ampere_problem_hash"] == right["ampere_problem_hash"]
        assert agrees == (identities[0] == identities[1])

    @pytest.mark.parametrize("name", sorted(SPECS))
    def test_the_declared_shape_of_the_problem_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        """Everything the problem hash *should* be a function of, checked directly."""
        left, right = self.attrs(first, name), self.attrs(second, name)
        for key in (
            "ampere_free_size",
            "ampere_free_names",
            "ampere_plates",
            "ampere_tied_names",
            "ampere_sites",
            "ampere_dataset_labels",
            "ampere_dataset_channels",
            "ampere_model_bindings",
            "ampere_likelihoods",
        ):
            assert left[key] == right[key], key

    @pytest.mark.parametrize("name", sorted(SPECS))
    def test_the_neutral_model_identity_agrees_across_backends(
        self, first: ConformanceBackend, second: ConformanceBackend, name: str
    ) -> None:
        """The "offer" half of ``results.md`` §14 (ruled 2026-09-03, landed W2.1).

        ``ampere_problem_hash`` is deliberately backend-variant — the row above
        asserts that, and it is what stops one backend's trained artefact being
        *served* to another's fit. The derived neutral identity is the other
        half: the fingerprint minus class and module, which two backends
        implementing the same declaration must agree on, so a cross-backend
        emulator can be **offered** with its provenance shown.

        This is a real assertion on these two fixtures rather than a tautology:
        ``reference`` and ``mirror`` compute the same closed forms by
        deliberately different routes, in different classes and modules. They
        disagree on the problem hash (above) and must agree here.
        """
        left, right = self.attrs(first, name), self.attrs(second, name)
        assert left["ampere_model_identity_hashes"] == right["ampere_model_identity_hashes"]


class TestEmissionAgreement:
    """The same draws through the same declaration give the same stored run."""

    def test_the_emitted_posterior_agrees(
        self, first: ConformanceBackend, second: ConformanceBackend
    ) -> None:
        pytest.importorskip("arviz", reason="the results contract needs ampere[arviz]")
        from ampere.results import POSTERIOR_GROUP, SAMPLE_STATS_GROUP, emit

        spec = SPECS["power-law-gp"]
        left, right = build_problem(first, spec), build_problem(second, spec)
        draws = np.array([left.prior_transform(sample_cube(left.free_size))])
        tolerance = shared_tolerance(first, second, "cross_backend")

        trees = [
            emit(problem, draws, [problem.evaluate(draws[0])], backend=backend.name)
            for problem, backend in ((left, first), (right, second))
        ]
        assert set(trees[0].children) == set(trees[1].children)
        for name in trees[0][POSTERIOR_GROUP]:
            assert trees[0][POSTERIOR_GROUP][name].values == pytest.approx(
                trees[1][POSTERIOR_GROUP][name].values, abs=tolerance
            )
        for name in ("lp", "log_prior", "log_likelihood"):
            assert trees[0][SAMPLE_STATS_GROUP][name].values == pytest.approx(
                trees[1][SAMPLE_STATS_GROUP][name].values, abs=tolerance
            )
        assert trees[0].attrs["ampere_backend"] != trees[1].attrs["ampere_backend"]
