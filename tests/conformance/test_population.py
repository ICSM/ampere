"""The population contract, once per registered backend (**W5.12**).

``hierarchical_population.md`` H-1 and H-2, made a lockstep obligation. Three
claims, and the suite owes all three on every fixture:

1. **The joint density decomposes.** A two-level population's joint log
   density is the hyperprior plus the sum over members — written out here from
   ``scipy.stats`` rather than taken from the object under test, so the row is
   an oracle comparison and not a restatement.
2. **The two layouts are one model.** The plate layout (one array-valued site,
   element *i* routed to component *i*) and the flat layout (N scalar sites,
   the tie-based pattern ``parameters.md`` §9 documents) declare the same
   population, and agree on the joint density at matched values and on every
   per-member routed value exactly.
3. **The realisation lowers the plate faithfully.** On a backend with a
   registered realisation — torch and jax; the reference backend has none by
   design — the realised density of the plated problem agrees with the numpy
   contract path at many points, and with the realised density of the
   flattened twin. That is what "the plate lowering is bit-consistent with the
   flattened form" has to mean once it is measurable: jax opens a real
   ``numpyro.plate`` and torch carries an array-valued parameter, and those
   two structures must score one number.

The models are the fixture's own ``LINEAR`` model, one instance per member, so
no test body here names a backend.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    HierarchicalPrior,
    Parameter,
    Population,
    log_density,
    realise,
    registered_realisations,
)
from ampere.core.exceptions import ParameterError

from .composition import (
    ProblemSpec,
    build_instrument,
    build_likelihood,
    observed_container,
)
from .protocol import ConformanceBackend, ModelKind, ModelSpec, Tolerances

#: Members in the population every row here composes. Small, because each one
#: is a dataset and a model: the claims are about structure, not scale.
MEMBERS = 4

#: The population's hyperpriors, and the oracle's copy of them. Declared once
#: so the two cannot drift — the row is worthless if the "independent" oracle
#: is reading the same objects the problem does, so what is shared is the
#: *numbers*, and each side builds its own frozen distribution from them.
MU_LOC, MU_SCALE = 1.0, 0.5
SIGMA_SCALE = 0.4


def member_labels(count: int = MEMBERS) -> tuple[str, ...]:
    """The model component labels the population is declared over."""
    return tuple(f"obj{index}" for index in range(count))


def population(*, layout: str = "plate", count: int = MEMBERS) -> Population:
    """``slope_i ~ Normal(mu, sigma)`` over *count* model components.

    ``slope`` is the fixture's ``LINEAR`` model's own parameter, declared with
    an ordinary prior; the population replaces that prior, which is exactly
    H-1's point — the model is a library model and was not written for this
    fit.
    """
    return Population(
        "objects",
        members=[Parameter("slope", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
        hyperpriors=[
            Parameter("mu", st.norm(MU_LOC, MU_SCALE)),
            Parameter("sigma", st.halfnorm(0.0, SIGMA_SCALE)),
        ],
        over=member_labels(count),
        layout=layout,
    )


def build_population_problem(
    backend: ConformanceBackend,
    *,
    layout: str = "plate",
    count: int = MEMBERS,
    through_collection: bool = False,
) -> FittingProblem:
    """A population of *count* single-dataset members on *backend*.

    One model instance and one dataset per member, the dataset bound to its
    own model by label, and the population declared over the model labels —
    which is what ``DatasetCollection.plate`` does for the caller, and this
    helper exercises both routes so the convenience and the explicit form are
    held to the same claims.
    """
    spec = ProblemSpec(model=ModelSpec(kind=ModelKind.LINEAR))
    labels = member_labels(count)
    models: dict[str, Any] = {}
    datasets = []
    for index, label in enumerate(labels):
        dataset = dataclasses.replace(spec.datasets[0], label=f"d{index}", data_seed=1000 + index)
        observed = observed_container(spec, dataset)
        models[label] = backend.model(spec.model)
        datasets.append(
            Dataset(
                observed,
                build_instrument(backend, dataset),
                build_likelihood(backend, dataset, observed),
                model=label,
                label=dataset.label,
            )
        )
    if through_collection:
        collection = DatasetCollection.plate(
            "objects",
            datasets,
            members=[
                Parameter("slope", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))
            ],
            hyperpriors=[
                Parameter("mu", st.norm(MU_LOC, MU_SCALE)),
                Parameter("sigma", st.halfnorm(0.0, SIGMA_SCALE)),
            ],
            layout=layout,
        )
        return FittingProblem(models, collection, seed=spec.seed)
    return FittingProblem(
        models, datasets, populations=[population(layout=layout, count=count)], seed=spec.seed
    )


def flatten(values: dict[str, Any], count: int = MEMBERS) -> dict[str, Any]:
    """A plate-layout value mapping rewritten onto the flat layout's names."""
    flat: dict[str, Any] = {}
    for name, value in values.items():
        if name == "objects.slope":
            for index, label in enumerate(member_labels(count)):
                flat[f"{label}.slope"] = np.asarray(value)[index]
        else:
            flat[name] = value
    return flat


def points(problem: FittingProblem) -> list[np.ndarray]:
    """Unconstrained vectors spanning the space, the pattern §10a's rows use."""
    size = problem.free_size
    rng = np.random.default_rng(20260918)
    cube = np.concatenate(
        [
            np.full((1, size), 0.5),
            np.linspace(0.1, 0.9, 5).reshape(-1, 1) * np.ones((1, size)),
            rng.uniform(0.05, 0.95, size=(8, size)),
        ]
    )
    interior = [problem.unconstrain(problem.prior_transform(row)) for row in cube]
    return [*interior, np.full(size, -6.0), np.full(size, 6.0)]


class TestTheDeclaration:
    """What the merge produces, on every backend: one plate, N element bindings."""

    def test_the_plate_layout_is_one_array_valued_site(self, backend: ConformanceBackend) -> None:
        problem = build_population_problem(backend)
        plated = problem.parameters["objects.slope"]
        assert plated.shape == (MEMBERS,)
        assert plated.plate == "objects"
        assert plated.references == ("objects.mu", "objects.sigma")
        # The members' own declaration of ``slope`` is gone: the element
        # binding is what reaches them now, which is H-2's whole point.
        assert [name for name in problem.parameters.names if name.endswith("slope")] == [
            "objects.slope"
        ]

    def test_every_member_is_addressed_by_its_own_index(self, backend: ConformanceBackend) -> None:
        problem = build_population_problem(backend)
        elements = {
            binding.component: binding.index
            for binding in problem.mapping.bindings
            if binding.index is not None
        }
        assert elements == dict(zip(member_labels(), range(MEMBERS), strict=True))
        # Addressing is not tying: each element remains its own draw.
        assert "objects.slope" not in problem.tied_names

    def test_the_element_reaches_the_member_it_addresses(self, backend: ConformanceBackend) -> None:
        problem = build_population_problem(backend)
        values = dict(problem.reference_values)
        values["objects.slope"] = np.arange(MEMBERS, dtype=float) + 0.5
        routed = problem.mapping.distribute(problem.parameters.complete(values))
        for index, label in enumerate(member_labels()):
            assert routed[label]["slope"] == pytest.approx(index + 0.5)

    def test_the_collection_factory_declares_the_same_thing(
        self, backend: ConformanceBackend
    ) -> None:
        """``inference.md`` §9's plate of datasets, and the explicit form."""
        through_factory = build_population_problem(backend, through_collection=True)
        explicit = build_population_problem(backend)
        assert through_factory.parameters.names == explicit.parameters.names
        assert through_factory.free_labels() == explicit.free_labels()

    def test_a_flat_population_beyond_the_limit_is_refused_by_name(self) -> None:
        with pytest.raises(ParameterError, match="layout='plate'"):
            Population(
                "objects",
                members=[Parameter("slope", HierarchicalPrior("norm", {"loc": "mu"}))],
                hyperpriors=[Parameter("mu", st.norm(0.0, 1.0))],
                over=[f"obj{index}" for index in range(500)],
                layout="flat",
            )


class TestTheJointDensity:
    """Claim 1: the joint density is the hyperprior plus the sum over members."""

    @pytest.mark.parametrize("layout", ["plate", "flat"])
    def test_the_joint_log_prior_decomposes(
        self, backend: ConformanceBackend, tolerances: Tolerances, layout: str
    ) -> None:
        problem = build_population_problem(backend, layout=layout)
        rng = np.random.default_rng(20260912)
        for _ in range(8):
            values = problem.sample_prior(rng)
            mu = float(values["objects.mu"])
            sigma = float(values["objects.sigma"])
            draws = (
                np.asarray(values["objects.slope"], dtype=float)
                if layout == "plate"
                else np.array([float(values[f"{label}.slope"]) for label in member_labels()])
            )
            hyperprior = float(
                st.norm(MU_LOC, MU_SCALE).logpdf(mu) + st.halfnorm(0.0, SIGMA_SCALE).logpdf(sigma)
            )
            members = float(np.sum(st.norm(mu, sigma).logpdf(draws)))
            # Everything the population does not own: each member model's
            # ``offset``, and whatever the instrument and likelihood declare.
            rest = sum(
                float(np.sum(log_density(problem.parameters[name].prior, values[name])))
                for name in problem.parameters.free_names
                if name not in {"objects.mu", "objects.sigma", "objects.slope"}
                and not name.endswith(".slope")
            )
            assert problem.log_prior(values) == pytest.approx(
                hyperprior + members + rest, abs=tolerances.cross_backend
            )

    def test_a_member_draw_is_scored_against_the_sampled_hyperparameters(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Two levels, not one: moving ``mu`` moves every member's prior.

        The failure this catches is the plausible one the sketch names — a
        population that silently fits *one* θ for the whole survey, or one
        whose members are scored against the hyperprior's own location rather
        than the sampled ``mu``.
        """
        problem = build_population_problem(backend)
        values = dict(problem.parameters.complete(problem.reference_values))
        values["objects.slope"] = np.full(MEMBERS, 1.0)
        shifted = dict(values)
        shifted["objects.mu"] = float(values["objects.mu"]) + 0.3
        sigma = float(values["objects.sigma"])
        expected = float(
            st.norm(MU_LOC, MU_SCALE).logpdf(float(shifted["objects.mu"]))
            - st.norm(MU_LOC, MU_SCALE).logpdf(float(values["objects.mu"]))
            + np.sum(st.norm(float(shifted["objects.mu"]), sigma).logpdf(np.full(MEMBERS, 1.0)))
            - np.sum(st.norm(float(values["objects.mu"]), sigma).logpdf(np.full(MEMBERS, 1.0)))
        )
        got = problem.log_prior(shifted) - problem.log_prior(values)
        assert got == pytest.approx(expected, abs=tolerances.cross_backend)


class TestTheTwoLayoutsAreOneModel:
    """Claim 2: the plate layout and the flattened form declare one density."""

    def test_the_layouts_agree_on_the_joint_density(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        plated = build_population_problem(backend, layout="plate")
        flat = build_population_problem(backend, layout="flat")
        assert plated.free_size == flat.free_size
        rng = np.random.default_rng(20260913)
        for _ in range(6):
            values = plated.parameters.complete(plated.sample_prior(rng))
            twin = flatten(dict(values))
            assert flat.log_prior(twin) == pytest.approx(
                plated.log_prior(values), abs=tolerances.cross_backend
            )
            assert flat.log_prob(twin) == pytest.approx(
                plated.log_prob(values), abs=tolerances.cross_backend
            )

    def test_the_layouts_route_the_same_values(self, backend: ConformanceBackend) -> None:
        """Exactly, not approximately: routing moves numbers, it does not compute."""
        plated = build_population_problem(backend, layout="plate")
        flat = build_population_problem(backend, layout="flat")
        values = plated.parameters.complete(plated.reference_values)
        values = dict(values)
        values["objects.slope"] = np.linspace(0.4, 1.6, MEMBERS)
        twin = flatten(dict(values))
        routed = plated.mapping.distribute(values)
        twin_routed = flat.mapping.distribute(flat.parameters.complete(twin))
        for label in member_labels():
            assert routed[label]["slope"] == twin_routed[label]["slope"]


class TestTheRealisedPlate:
    """Claim 3: the native lowering of a plated population (``inference.md`` §10a)."""

    def realised(self, problem: FittingProblem) -> Any:
        if problem.backend not in registered_realisations():
            pytest.skip(
                f"the {problem.backend!r} backend registers no realisation "
                f"(inference.md §10a: the reference backend has no differentiable path)"
            )
        return realise(problem)

    def test_the_realised_population_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_population_problem(backend)
        realisation = self.realised(problem)
        assert int(realisation.free_size) == problem.free_size
        for y in points(problem):
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(realisation.log_prob_unconstrained(y))))
            if not np.isfinite(expected):
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)

    def test_the_plate_lowering_matches_the_flattened_form(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """One numpyro plate (or one array-valued tensor) against N scalar sites.

        The structures differ deliberately — ``parameters.md`` §9 chose the
        plate *because* it is what a ``numpyro.plate`` is — so this row is
        what keeps the choice from being a change of model.
        """
        plated = build_population_problem(backend, layout="plate")
        flat = build_population_problem(backend, layout="flat")
        plated_realisation = self.realised(plated)
        flat_realisation = self.realised(flat)
        rng = np.random.default_rng(20260914)
        for _ in range(6):
            values = plated.parameters.complete(plated.sample_prior(rng))
            twin = flat.parameters.complete(flatten(dict(values)))
            got = float(
                np.asarray(
                    backend.to_numpy(
                        plated_realisation.log_prob_unconstrained(plated.unconstrain(values))
                    )
                )
            )
            expected = float(
                np.asarray(
                    backend.to_numpy(
                        flat_realisation.log_prob_unconstrained(flat.unconstrain(twin))
                    )
                )
            )
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)
