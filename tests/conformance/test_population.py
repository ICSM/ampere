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

**W7.1** adds the population over a dataset's qualified path
(``TestAPopulationOverADatasetPath``, the ten rows of
``nuisance_populations_and_derived_memo.md`` §9.1): the members reach the
datasets themselves — their GP amplitude, their calibration step — rather
than the models they name.
"""

from __future__ import annotations

import dataclasses
import json
import math
from collections.abc import Sequence
from typing import Any

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    Derived,
    FittingProblem,
    HierarchicalPrior,
    Log,
    Parameter,
    ParameterSet,
    Population,
    Tie,
    log_density,
    realise,
    registered_realisations,
)
from ampere.core.exceptions import ParameterError

from .composition import (
    NoiseKind,
    ProblemSpec,
    build_instrument,
    build_likelihood,
    observed_container,
)
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    ModelKind,
    ModelSpec,
    Tolerances,
    TransformationKind,
    TransformationSpec,
)

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


# ---------------------------------------------------------------------------
# W7.1: a population over a dataset's qualified path
# ---------------------------------------------------------------------------

#: Datasets in the path rows. Each is a flexible-likelihood (GP) dataset with
#: one calibration step, so both of the memo's motivating cases are present.
PATH_MEMBERS = 3

#: The population's spread hyperprior, and the oracle's copy of its scale.
SPREAD_SCALE = 0.5

#: The GP amplitude every dataset declares *free* (the battery holds kernel
#: hyperparameters fixed by default, ``CovarianceSpec``; a fixed site has no
#: draw for a population to replace). The member's prior replaces it.
FREE_AMPLITUDE = st.lognorm(0.5, scale=0.4)

#: One problem seed for the path rows, so a dataset built alone (row 1's
#: oracle) draws exactly the data its twin inside the population does.
PATH_SEED = 20261009


def dataset_labels(count: int = PATH_MEMBERS) -> tuple[str, ...]:
    """The dataset labels the path population is declared over, in plate order."""
    return tuple(f"d{index}" for index in range(count))


def path_datasets(
    backend: ConformanceBackend, count: int = PATH_MEMBERS
) -> tuple[dict[str, Any], list[Dataset]]:
    """*count* GP datasets, each with a ``calibrate`` step and its own model.

    Each dataset's own merged names are ``instrument.calibrate.scale``,
    ``likelihood.amplitude`` and ``likelihood.length_scale``: the amplitude
    one level down, the calibration scale two.
    """
    spec = ProblemSpec(model=ModelSpec(kind=ModelKind.LINEAR))
    models: dict[str, Any] = {}
    datasets: list[Dataset] = []
    for index, label in enumerate(dataset_labels(count)):
        dataset = dataclasses.replace(
            spec.datasets[0],
            label=label,
            data_seed=3000 + index,
            noise=NoiseKind.GP,
            covariance=CovarianceSpec(amplitude=FREE_AMPLITUDE),  # type: ignore[arg-type]
            instrument=(TransformationSpec(kind=TransformationKind.SCALE, label="calibrate"),),
        )
        observed = observed_container(spec, dataset)
        models[f"m{index}"] = backend.model(spec.model)
        datasets.append(
            Dataset(
                observed,
                build_instrument(backend, dataset),
                build_likelihood(backend, dataset, observed),
                model=f"m{index}",
                label=label,
            )
        )
    return models, datasets


def spread_members(member: str = "amplitude") -> tuple[list[Parameter], list[Parameter]]:
    """``member_i ~ LogNormal(s=spread)``, ``spread ~ HalfNormal``: members, hyperpriors.

    The member declares its ``Log`` bijection: ``lognorm`` takes a shape
    argument, so its support cannot be inferred from a hierarchical
    declaration (``parameters.md`` §6, *Amended W5.30*), and the realised rows
    unconstrain it.
    """
    return (
        [Parameter(member, HierarchicalPrior("lognorm", {"s": "spread"}), bijection=Log())],
        [Parameter("spread", st.halfnorm(0.0, SPREAD_SCALE))],
    )


def path_population(
    *,
    within: str = "likelihood",
    member: str = "amplitude",
    layout: str = "plate",
    over: Sequence[str] | None = None,
) -> Population:
    """The memo's §2.1 population, over ``d*.<within>``."""
    members, hyperpriors = spread_members(member)
    return Population(
        "gp",
        members=members,
        hyperpriors=hyperpriors,
        over=[f"{label}.{within}" for label in dataset_labels()] if over is None else over,
        layout=layout,
    )


def build_path_problem(
    backend: ConformanceBackend,
    population: Population | None = None,
    *,
    ties: Sequence[Tie] = (),
    **kwargs: Any,
) -> FittingProblem:
    models, datasets = path_datasets(backend)
    declared = path_population(**kwargs) if population is None else population
    return FittingProblem(models, datasets, ties=list(ties), populations=[declared], seed=PATH_SEED)


def path_values(problem: FittingProblem, **overrides: Any) -> dict[str, Any]:
    values = dict(problem.parameters.complete(problem.reference_values))
    values.update(overrides)
    return dict(problem.parameters.complete(values))


class TestAPopulationOverADatasetPath:
    """W7.1, the memo's §9.1: ``Population.over`` entries of the form ``component[.path]``.

    A per-dataset nuisance — the flexible likelihood's GP amplitude, a
    calibration step's scale — drawn from one shared prior. Routing is
    unchanged: the dataset's retained mapping takes the second hop it already
    takes for every other value, so these rows hold on every fixture with no
    backend code of its own.
    """

    def test_the_leaf_is_stripped_and_the_element_reaches_the_noise_model(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_path_problem(backend)
        assert problem.parameters["gp.amplitude"].shape == (PATH_MEMBERS,)
        assert problem.parameters["gp.amplitude"].plate == "gp"
        assert not [n for n in problem.parameters.names if n.endswith("likelihood.amplitude")]
        amplitudes = np.array([0.25, 0.6, 1.3])
        values = path_values(problem, **{"gp.amplitude": amplitudes})
        contributions = problem.evaluate(values).contributions
        models, datasets = path_datasets(backend)
        for index, label in enumerate(dataset_labels()):
            # The oracle: the same dataset built alone, its own amplitude set
            # to a_i and every other value matched.
            alone = FittingProblem(
                {f"m{index}": models[f"m{index}"]}, [datasets[index]], seed=PATH_SEED
            )
            own = dict(alone.parameters.complete(alone.reference_values))
            own.update({name: values[name] for name in own if name in values})
            own[f"{label}.likelihood.amplitude"] = amplitudes[index]
            expected = alone.evaluate(own).log_likelihood
            assert math.isfinite(expected)
            assert contributions[label] == pytest.approx(expected, abs=tolerances.cross_solver)

    def test_the_joint_log_prior_decomposes_over_the_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_path_problem(backend)
        rng = np.random.default_rng(20261010)
        for _ in range(6):
            values = problem.sample_prior(rng)
            spread = float(values["gp.spread"])
            amplitudes = np.asarray(values["gp.amplitude"], dtype=float)
            hyperprior = float(st.halfnorm(0.0, SPREAD_SCALE).logpdf(spread))
            members = float(np.sum(st.lognorm(spread).logpdf(amplitudes)))
            rest = sum(
                float(np.sum(log_density(problem.parameters[name].prior, values[name])))
                for name in problem.parameters.free_names
                if not name.startswith("gp.")
            )
            assert problem.log_prior(values) == pytest.approx(
                hyperprior + members + rest, abs=tolerances.cross_backend
            )

    def test_a_two_level_path_reaches_an_instrument_step(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_path_problem(backend, within="instrument.calibrate", member="scale")
        assert problem.parameters["gp.scale"].shape == (PATH_MEMBERS,)
        assert not [n for n in problem.parameters.names if n.endswith("calibrate.scale")]
        scales = np.array([0.8, 1.1, 1.4])
        unit = problem.predict(path_values(problem, **{"gp.scale": np.ones(PATH_MEMBERS)}))
        scaled = problem.predict(path_values(problem, **{"gp.scale": scales}))
        for index, label in enumerate(dataset_labels()):
            ratio = np.asarray(backend.to_numpy(scaled[label].values)) / np.asarray(
                backend.to_numpy(unit[label].values)
            )
            np.testing.assert_allclose(ratio, scales[index], rtol=tolerances.cross_backend)

    def test_the_two_layouts_agree_on_the_joint_density(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        plated = build_path_problem(backend, layout="plate")
        flat = build_path_problem(backend, layout="flat")
        assert plated.free_size == flat.free_size
        assert "d0.likelihood.amplitude" in flat.parameters.names
        rng = np.random.default_rng(20261011)
        for _ in range(5):
            values = dict(plated.parameters.complete(plated.sample_prior(rng)))
            twin = {name: value for name, value in values.items() if name != "gp.amplitude"}
            for index, label in enumerate(dataset_labels()):
                twin[f"{label}.likelihood.amplitude"] = np.asarray(values["gp.amplitude"])[index]
            assert flat.log_prior(twin) == pytest.approx(
                plated.log_prior(values), abs=tolerances.cross_backend
            )
            assert flat.log_prob(twin) == pytest.approx(
                plated.log_prob(values), abs=tolerances.cross_backend
            )

    def test_the_realised_population_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_path_problem(backend)
        if problem.backend not in registered_realisations():
            pytest.skip(
                f"the {problem.backend!r} backend registers no realisation "
                f"(inference.md §10a: the reference backend has no differentiable path)"
            )
        realisation = realise(problem)
        assert int(realisation.free_size) == problem.free_size
        for y in points(problem):
            native = realisation.log_prob_unconstrained(y)
            # The element reaches the noise model as the backend's own array:
            # a coercion to numpy on the way would surface here as a numpy
            # scalar rather than a native one.
            assert not isinstance(native, (float, np.floating, np.ndarray))
            expected = problem.log_prob_unconstrained(y)
            got = float(np.asarray(backend.to_numpy(native)))
            if not np.isfinite(expected):
                assert not np.isfinite(got)
                continue
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)

    def test_refusals_by_name(self, backend: ConformanceBackend) -> None:
        # A path that does not resolve names the components at that level.
        with pytest.raises(
            ParameterError, match=r"no component 'noise'.*\['instrument', 'likelihood'\]"
        ):
            build_path_problem(backend, within="noise")
        with pytest.raises(ParameterError, match=r"no component 'calib'.*\['calibrate'\]"):
            build_path_problem(backend, within="instrument.calib", member="scale")
        # A path on a plain-set component (a model).
        with pytest.raises(ParameterError, match=r"'m0' is a ParameterSet, not a composite"):
            build_path_problem(backend, over=[f"m{i}.likelihood" for i in range(PATH_MEMBERS)])
        # A leaf the composite does not declare.
        with pytest.raises(ParameterError, match=r"declares no 'likelihood\.tilt'"):
            build_path_problem(backend, member="tilt")
        # A tied leaf.
        with pytest.raises(ParameterError, match=r"d0\.likelihood\.amplitude.*already shared"):
            build_path_problem(
                backend,
                ties=[
                    Tie("amp", tuple(f"{label}.likelihood.amplitude" for label in dataset_labels()))
                ],
            )
        # A composite without a path: W5.12's refusal, re-worded to the remedy.
        with pytest.raises(
            ParameterError,
            match=r"'d0'.*without a path.*over=\['d0\.likelihood'.*within='likelihood'",
        ):
            build_path_problem(backend, over=list(dataset_labels()))
        # A leaf an inner shared_as collapsed is named by its tie label: each
        # composite here is the real likelihood's declaration beside a twin,
        # their amplitudes collapsed into ``amp`` by the inner merge.
        _, datasets = path_datasets(backend)
        composites = {}
        for dataset in datasets:
            shared = ParameterSet(
                [
                    dataclasses.replace(p, shared_as="amp") if p.name == "amplitude" else p
                    for p in dataset.likelihood.parameters
                ]
            )
            composites[dataset.label] = ParameterSet.merge({"likelihood": shared, "twin": shared})
        assert "amp" in composites["d0"].merged.names
        with pytest.raises(ParameterError, match=r"collapsed it into d0\.amp \(tie label 'amp'\)"):
            ParameterSet.merge(composites, populations=[path_population()])

    def test_the_factory_writes_the_paths(self, backend: ConformanceBackend) -> None:
        models, datasets = path_datasets(backend)
        members, hyperpriors = spread_members()
        collection = DatasetCollection.plate(
            "gp", datasets, members=members, hyperpriors=hyperpriors, within="likelihood"
        )
        explicit = Population(
            "gp",
            members=members,
            hyperpriors=hyperpriors,
            over=[f"{label}.likelihood" for label in dataset_labels()],
        )
        assert collection.populations == (explicit,)
        through_factory = FittingProblem(models, collection, seed=PATH_SEED)
        direct = build_path_problem(backend)
        assert through_factory.parameters.names == direct.parameters.names
        assert through_factory.free_labels() == direct.free_labels()
        assert through_factory.sites() == direct.sites()

    def test_the_population_is_in_provenance(self, backend: ConformanceBackend) -> None:
        from ampere.results import PROVENANCE_SCHEMA_VERSION, provenance_attrs

        problem = build_path_problem(backend)
        attrs = provenance_attrs(problem)
        assert PROVENANCE_SCHEMA_VERSION == 12
        assert attrs["ampere_schema_version"] == 12
        assert json.loads(attrs["ampere_populations"]) == [
            {
                "name": "gp",
                "label": "gp",
                "layout": "plate",
                "over": ["d0.likelihood", "d1.likelihood", "d2.likelihood"],
                "members": ["amplitude"],
                "hyperpriors": ["spread"],
            }
        ]
        # The same problem built again hashes the same; the same population
        # over the datasets in another order does not.
        again = provenance_attrs(build_path_problem(backend))
        assert again["ampere_problem_hash"] == attrs["ampere_problem_hash"]
        reordered = build_path_problem(
            backend, over=[f"{label}.likelihood" for label in reversed(dataset_labels())]
        )
        assert provenance_attrs(reordered)["ampere_problem_hash"] != attrs["ampere_problem_hash"]

    def test_a_non_centred_population_over_a_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The memo's §4: ``z`` internal, ``amplitude = exp(mu + sigma * z)`` derived.

        The centred twin declares ``amplitude ~ LogNormal(s=sigma,
        scale=exp(mu))`` over the same path. At matched θ the likelihoods are
        equal and the priors differ by exactly the Jacobian of the N-fold map
        ``z ↦ exp(mu + sigma z)``, ``Σ log(sigma · a_i)``.
        """
        over = [f"{label}.likelihood" for label in dataset_labels()]
        mu, sigma = Parameter("mu", st.norm(-1.0, 0.5)), Parameter("sigma", st.halfnorm(0.0, 0.5))
        non_centred = build_path_problem(
            backend,
            Population(
                "gp",
                members=[
                    Parameter("z", st.norm(0.0, 1.0)),
                    Parameter("amplitude", Derived("exp(mu + sigma * z)")),
                ],
                hyperpriors=[mu, sigma],
                over=over,
            ),
        )
        centred = build_path_problem(
            backend,
            Population(
                "gp",
                members=[
                    Parameter(
                        "amplitude", HierarchicalPrior("lognorm", {"s": "sigma", "scale": "scale"})
                    )
                ],
                hyperpriors=[mu, sigma, Parameter("scale", Derived("exp(mu)"))],
                over=over,
            ),
        )
        assert "gp.z" in non_centred.parameters.free_names
        assert "gp.amplitude" in non_centred.parameters.derived_names
        elements = {
            b.global_name: (b.component, b.local_name)
            for b in non_centred.mapping.bindings
            if b.index is not None and b.index == 0
        }
        # z is internal and routed nowhere; the derived amplitude is routed by
        # element, under the dataset-level path.
        assert elements == {"gp.amplitude": ("d0", "likelihood.amplitude")}
        rng = np.random.default_rng(20261012)
        for _ in range(5):
            drawn = non_centred.sample_prior(rng)
            values = dict(non_centred.parameters.complete(drawn))
            amplitudes = np.asarray(values["gp.amplitude"], dtype=float)
            np.testing.assert_allclose(
                amplitudes,
                np.exp(float(values["gp.mu"]) + float(values["gp.sigma"]) * values["gp.z"]),
                rtol=1e-12,
            )
            twin = {name: value for name, value in values.items() if name != "gp.z"}
            twin = dict(centred.parameters.complete(twin))
            left, right = non_centred.evaluate(values), centred.evaluate(twin)
            assert left.log_likelihood == pytest.approx(
                right.log_likelihood, abs=tolerances.cross_backend
            )
            jacobian = float(np.sum(np.log(float(values["gp.sigma"]) * amplitudes)))
            assert left.log_prior == pytest.approx(
                right.log_prior + jacobian, abs=tolerances.cross_backend
            )


class TestThePathInResults:
    """Row 8 of the memo's §9.1, dev only: it goes through ``ampere.results``."""

    def test_the_plate_coordinate_is_the_dataset_labels(self) -> None:
        pytest.importorskip("arviz")
        from ampere.results import emit

        from .backends.reference import ReferenceBackend

        problem = build_path_problem(ReferenceBackend())
        rng = np.random.default_rng(20261013)
        rows = [problem.parameters.pack(problem.sample_prior(rng)) for _ in range(4)]
        draws = np.asarray(rows, dtype=float).reshape(2, 2, problem.free_size)
        evaluations = [[problem.evaluate(draw) for draw in chain] for chain in draws]
        tree = emit(problem, draws, evaluations, engine="emcee")
        amplitude = tree["posterior"]["gp.amplitude"]
        assert amplitude.dims[-1] == "gp"
        assert [str(label) for label in np.asarray(amplitude.coords["gp"])] == list(
            dataset_labels()
        )
