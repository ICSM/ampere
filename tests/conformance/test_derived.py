"""The ``Derived`` parameter node (**W7.0**), the memo's §9.2 rows.

``docs/design/nuisance_populations_and_derived_memo.md`` §3: a fourth
parameter state, a deterministic function of others, declared as a closed
grammar over symbols. The rows are named one per claim of §9.2, and each is
an oracle comparison — the expected value is written out here from
``scipy.stats`` and ``numpy``, never read back from the object under test.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    Derived,
    FittingProblem,
    HierarchicalPrior,
    Parameter,
    ParameterSet,
    Population,
    realise,
    registered_realisations,
    shrinkage_horseshoe,
)
from ampere.core.exceptions import ParameterError, TyingError

from .composition import ProblemSpec, build_instrument, build_likelihood, observed_container
from .protocol import ConformanceBackend, ModelKind, ModelSpec, Tolerances

#: The non-centred triple every row composes: ``theta = mu + sigma * z``.
EXPRESSION = "mu + sigma * z"


def non_centred_set(size: int = 3) -> ParameterSet:
    """``mu``, ``sigma``, an array ``z`` and the derived ``theta``, declared out of order.

    ``theta`` is declared *first*, so the rows that check it is formed after
    its inputs are testing the evaluation order and not the declaration order.
    """
    return ParameterSet(
        [
            Parameter("theta", Derived(EXPRESSION), shape=(size,)),
            Parameter("mu", st.norm(0.0, 1.0)),
            Parameter("sigma", st.halfnorm(0.0, 1.0)),
            Parameter("z", st.norm(0.0, 1.0), shape=(size,)),
        ]
    )


def slab_set() -> ParameterSet:
    """A hierarchical prior over a derived scale, declared before what it references."""
    return ParameterSet(
        [
            Parameter("a", HierarchicalPrior("halfnorm", {"scale": "s_eff"})),
            Parameter("s_eff", Derived("sqrt(c**2 * s**2 / (c**2 + s**2))")),
            Parameter("c", st.halfnorm(0.0, 2.0)),
            Parameter("s", st.halfcauchy(0.0, 1.0)),
        ]
    )


def effective_scale(c: float, s: float) -> float:
    """§3.9's slab, written out independently of the expression under test."""
    return math.sqrt(c * c * s * s / (c * c + s * s))


class TestTheState:
    """Rows 1 and 2: not a sampler dimension; computed, idempotently."""

    def test_a_derived_parameter_is_not_a_free_dimension(self) -> None:
        pset = non_centred_set()
        assert pset.free_size == 5
        assert pset.free_names == ("mu", "sigma", "z")
        assert pset.free_labels() == ("mu", "sigma", "z[0]", "z[1]", "z[2]")
        assert pset.derived_names == ("theta",)
        assert "theta" in pset.names
        order = pset.evaluation_order()
        assert order.index("theta") > max(order.index(name) for name in ("mu", "sigma", "z"))
        theta = pset["theta"]
        assert theta.is_derived
        assert not (theta.is_free or theta.is_fixed or theta.is_deferred)
        vector = np.array([0.5, 2.0, -1.0, 0.0, 1.5])
        values = pset.unpack(vector)
        np.testing.assert_array_equal(values["theta"], 0.5 + 2.0 * np.array([-1.0, 0.0, 1.5]))
        np.testing.assert_array_equal(pset.pack(values), vector)
        with pytest.raises(KeyError, match="derived"):
            pset.free_slice("theta")

    def test_complete_is_idempotent_and_computes_it(self) -> None:
        pset = non_centred_set()
        partial = {"mu": 1.0, "sigma": 0.5, "z": np.array([0.0, 2.0, -2.0])}
        once = pset.complete(partial)
        np.testing.assert_array_equal(once["theta"], [1.0, 2.0, 0.0])
        twice = pset.complete(once)
        assert set(twice) == set(once)
        for name in once:
            np.testing.assert_array_equal(twice[name], once[name])
        # pack ignores the derived entry, stale or not.
        stale = dict(once, theta=np.array([9.0, 9.0, 9.0]))
        np.testing.assert_array_equal(pset.pack(stale), pset.pack(once))
        # On the numpy reference path a stale value is refused, by name.
        with pytest.raises(ParameterError, match="derived parameter 'theta'"):
            pset.complete(stale)
        with pytest.raises(ParameterError, match="derived parameter 'theta'"):
            pset.lnprior(stale)
        # Round-off is not staleness.
        jittered = dict(once, theta=once["theta"] * (1.0 + 1e-13))
        np.testing.assert_array_equal(pset.complete(jittered)["theta"], once["theta"])

    @pytest.mark.parametrize(
        ("keywords", "named"),
        [
            ({"fixed": True}, "fixed=True"),
            ({"value": 1.0}, "value"),
            ({"shared_as": "tied"}, "shared_as"),
        ],
    )
    def test_the_derived_state_refuses_what_belongs_to_its_inputs(
        self, keywords: dict[str, object], named: str
    ) -> None:
        with pytest.raises(ParameterError, match=named):
            Parameter("theta", Derived(EXPRESSION), **keywords)  # type: ignore[arg-type]
        derived = Parameter("theta", Derived(EXPRESSION))
        with pytest.raises(ParameterError, match="derived"):
            derived.fix(1.0)
        with pytest.raises(ParameterError, match="derived"):
            derived.release(st.norm())

    def test_a_tie_on_a_derived_parameter_is_refused(self) -> None:
        from ampere.core import Tie

        one = ParameterSet([Parameter("a", st.norm()), Parameter("t", Derived("2 * a"))])
        two = ParameterSet([Parameter("b", st.norm()), Parameter("t", Derived("3 * b"))])
        with pytest.raises(TyingError, match="tie the inputs"):
            ParameterSet.merge({"one": one, "two": two}, ties=[Tie("t", ["one.t", "two.t"])])

    def test_a_result_that_does_not_broadcast_is_refused_by_name(self) -> None:
        with pytest.raises(ParameterError, match=r"'theta'.*does not broadcast"):
            ParameterSet(
                [
                    Parameter("z", st.norm(), value=np.zeros(3)),
                    Parameter("theta", Derived("2 * z")),
                ]
            )
        lazy = ParameterSet(
            [Parameter("z", st.norm(), shape=(3,)), Parameter("theta", Derived("2 * z"))]
        )
        with pytest.raises(ParameterError, match=r"'theta'.*does not broadcast"):
            lazy.complete({"z": np.zeros(3)})


class TestTheGrammar:
    """Row 3: the grammar is closed, and refuses by name."""

    @pytest.mark.parametrize(
        ("expression", "named"),
        [
            ("mu.real", "Attribute"),
            ("z[0]", "Subscript"),
            ("mu < sigma", "Compare"),
            ("lambda: mu", "Lambda"),
            ("max(mu, sigma)", "'max'"),
            ("mu if sigma else z", "IfExp"),
            ("mu % 2", "Mod"),
            ("+mu", "UAdd"),
            ("'text'", "constant"),
            ("exp + 1", "reserved"),
        ],
    )
    def test_the_grammar_is_closed(self, expression: str, named: str) -> None:
        with pytest.raises(ParameterError, match=named):
            Derived(expression)

    def test_the_five_functions_and_the_operators_are_admitted(self) -> None:
        from ampere.core.kernels import NUMPY_OPS

        derived = Derived("-sqrt(abs(a)) + exp(log(b)) * log1p(c) / 2 ** a - 1.5")
        a, b, c = -4.0, 3.0, 0.25
        expected = -math.sqrt(abs(a)) + math.exp(math.log(b)) * math.log1p(c) / 2**a - 1.5
        assert float(derived.evaluate(NUMPY_OPS, {"a": a, "b": b, "c": c})) == pytest.approx(
            expected, rel=1e-15
        )

    def test_a_symbol_spelt_like_a_reserved_name_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="reserved"):
            Derived("2 * x", {"exp": "x"})
        with pytest.raises(ParameterError, match="does not bind"):
            Derived("a + b", {"a": "mu"})
        with pytest.raises(ParameterError, match="does not use"):
            Derived("a", {"a": "mu", "b": "sigma"})
        with pytest.raises(ParameterError, match="non-empty string"):
            Derived(lambda mu: mu)  # type: ignore[arg-type]


class TestTheSpec:
    """Row 4: the spec round-trips after a merge, and hashes stably."""

    def test_the_spec_round_trips_after_a_merge(self) -> None:
        component = ParameterSet(
            [
                Parameter("mu", st.norm()),
                Parameter("sigma", st.halfnorm()),
                Parameter("z", st.norm(), shape=(3,)),
                Parameter("theta", Derived("mu+sigma*z"), shape=(3,)),
            ]
        )
        other = ParameterSet([Parameter("mu", st.norm())])
        merged = ParameterSet.merge({"gp": component, "other": other}).merged
        derived = merged["gp.theta"].prior
        assert isinstance(derived, Derived)
        # The merge renamed the bindings to dotted names and left the expression alone.
        assert derived.expression == EXPRESSION
        assert dict(derived.symbols or {}) == {"mu": "gp.mu", "sigma": "gp.sigma", "z": "gp.z"}
        assert merged["gp.theta"].references == ("gp.mu", "gp.sigma", "gp.z")
        spec = merged.to_spec()
        assert {"name": "gp.theta", "derived": derived.to_dict(), "shape": [3]} in spec[
            "parameters"
        ]
        rebuilt = ParameterSet.from_spec(json.loads(json.dumps(spec)))
        assert rebuilt == merged
        values = merged.unpack(np.arange(merged.free_size, dtype=float))
        np.testing.assert_array_equal(
            rebuilt.unpack(np.arange(merged.free_size, dtype=float))["gp.theta"],
            values["gp.theta"],
        )

    def test_the_normalised_source_hashes_stably(self) -> None:
        assert Derived("a+b") == Derived("a + b")
        assert Derived("(a)+(b)").expression == "a + b"

        def spec_of(expression: str) -> str:
            pset = ParameterSet(
                [
                    Parameter("a", st.norm()),
                    Parameter("b", st.norm()),
                    Parameter("c", Derived(expression)),
                ]
            )
            return json.dumps(pset.to_spec(), sort_keys=True)

        assert spec_of("a+b") == spec_of("a + b") == spec_of("(a) + b")


class TestAHierarchicalPriorOverADerivedValue:
    """Row 5: ``bind`` finds a derived value, in ``lnprior`` and mid-walk."""

    def test_a_hierarchical_prior_may_reference_a_derived_parameter(self) -> None:
        pset = slab_set()
        assert pset.free_names == ("a", "c", "s")
        rng = np.random.default_rng(20261007)
        for _ in range(10):
            cube = rng.uniform(0.05, 0.95, size=3)
            c = float(st.halfnorm(0.0, 2.0).ppf(cube[1]))
            s = float(st.halfcauchy(0.0, 1.0).ppf(cube[2]))
            scale = effective_scale(c, s)
            a = float(st.halfnorm(0.0, scale).ppf(cube[0]))
            np.testing.assert_allclose(pset.prior_transform(cube), [a, c, s], rtol=1e-12)
            expected = (
                st.halfnorm(0.0, 2.0).logpdf(c)
                + st.halfcauchy(0.0, 1.0).logpdf(s)
                + st.halfnorm(0.0, scale).logpdf(a)
            )
            assert pset.lnprior(np.array([a, c, s])) == pytest.approx(expected, rel=1e-12)
            assert pset.unpack(np.array([a, c, s]))["s_eff"] == pytest.approx(scale, rel=1e-14)
        drawn = pset.sample(np.random.default_rng(3))
        assert drawn["s_eff"] == pytest.approx(effective_scale(drawn["c"], drawn["s"]))


#: Members in the population rows 7 and 8 compose.
MEMBERS = 4


def object_labels(count: int = MEMBERS) -> tuple[str, ...]:
    return tuple(f"obj{index}" for index in range(count))


def hyperpriors() -> list[Parameter]:
    return [Parameter("mu", st.norm(1.0, 0.5)), Parameter("sigma", st.halfnorm(0.0, 0.4))]


def centred(*, layout: str = "plate") -> Population:
    """``theta_i ~ Normal(mu, sigma)``, the declaration W5.12 landed."""
    return Population(
        "objects",
        members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
        hyperpriors=hyperpriors(),
        over=object_labels(),
        layout=layout,
    )


def non_centred(*, layout: str = "plate") -> Population:
    """The same density, sampled in ``(mu, sigma, z)``: ``theta = mu + sigma * z``."""
    return Population(
        "objects",
        members=[Parameter("z", st.norm(0.0, 1.0)), Parameter("theta", Derived(EXPRESSION))],
        hyperpriors=hyperpriors(),
        over=object_labels(),
        layout=layout,
    )


def objects() -> dict[str, ParameterSet]:
    """One component per member, each declaring only the routed ``theta``."""
    return {
        label: ParameterSet([Parameter("theta", st.norm(0.0, 1.0))]) for label in object_labels()
    }


class TestTheNonCentredPopulation:
    """Row 7: internal ``z``, derived ``theta``, and the two member rules."""

    def test_the_non_centred_population_declares_the_centred_density(self) -> None:
        nc = ParameterSet.merge(objects(), populations=[non_centred()])
        c = ParameterSet.merge(objects(), populations=[centred()])
        assert nc.merged.free_names == ("objects.mu", "objects.sigma", "objects.z")
        assert nc.merged["objects.z"].shape == (MEMBERS,)
        assert nc.merged["objects.theta"].shape == (MEMBERS,)
        assert nc.merged["objects.theta"].plate == "objects"
        rng = np.random.default_rng(20261007)
        for _ in range(25):
            mu = float(rng.normal(1.0, 0.5))
            sigma = float(rng.uniform(0.05, 1.5))
            z = rng.normal(size=MEMBERS)
            theta = mu + sigma * z
            values = nc.merged.complete({"objects.mu": mu, "objects.sigma": sigma, "objects.z": z})
            np.testing.assert_allclose(values["objects.theta"], theta, rtol=1e-15)
            lnprior_nc = nc.merged.lnprior(values)
            lnprior_c = c.merged.lnprior(
                {"objects.mu": mu, "objects.sigma": sigma, "objects.theta": theta}
            )
            # The Jacobian of the N-fold affine map z -> mu + sigma z.
            assert lnprior_nc == pytest.approx(lnprior_c + MEMBERS * math.log(sigma), abs=1e-10)

    def test_the_internal_member_reaches_no_component(self) -> None:
        nc = ParameterSet.merge(objects(), populations=[non_centred()])
        assert all(
            binding.local_name != "z" for binding in nc.bindings if binding.index is not None
        )
        assert non_centred().plate_bindings(internal=("z",)) == tuple(
            binding for binding in non_centred().plate_bindings() if binding.local_name != "z"
        )
        values = nc.merged.complete(
            {"objects.mu": 0.5, "objects.sigma": 2.0, "objects.z": np.arange(MEMBERS, dtype=float)}
        )
        routed = nc.distribute(values)
        for index, label in enumerate(object_labels()):
            assert routed[label] == {"theta": pytest.approx(0.5 + 2.0 * index)}

    def test_a_plate_may_derive_a_member_from_a_sibling_member(self) -> None:
        """``Plate.expand`` qualifies a derived member's bindings onto its siblings too."""
        plated = non_centred().as_plate().to_parameter_set()
        assert plated["objects.theta"].references == ("objects.mu", "objects.sigma", "objects.z")
        values = plated.complete(
            {"objects.mu": 1.0, "objects.sigma": 2.0, "objects.z": np.arange(MEMBERS, dtype=float)}
        )
        np.testing.assert_array_equal(values["objects.theta"], 1.0 + 2.0 * np.arange(MEMBERS))

    def test_an_undeclared_member_that_nothing_derives_from_is_refused(self) -> None:
        population = Population(
            "objects",
            members=[
                Parameter("z", st.norm(0.0, 1.0)),
                Parameter("w", st.norm(0.0, 1.0)),
                Parameter("theta", Derived(EXPRESSION)),
            ],
            hyperpriors=hyperpriors(),
            over=object_labels(),
        )
        with pytest.raises(ParameterError, match="member 'w' is declared by none"):
            ParameterSet.merge(objects(), populations=[population])

    def test_a_member_declared_by_only_some_components_is_refused(self) -> None:
        partial = objects()
        partial["obj2"] = ParameterSet([Parameter("other", st.norm())])
        with pytest.raises(ParameterError, match=r"\['obj2'\] do not declare it"):
            ParameterSet.merge(partial, populations=[centred()])
        # The flat layout no longer gives an undeclaring component one.
        with pytest.raises(ParameterError, match=r"\['obj2'\] do not declare it"):
            ParameterSet.merge(partial, populations=[centred(layout="flat")])

    def test_derived_and_internal_members_are_refused_in_the_flat_layout(self) -> None:
        with pytest.raises(ParameterError, match=r"derived member\(s\) \['theta'\]"):
            non_centred(layout="flat")
        internal_only = Population(
            "objects",
            members=[Parameter("w", st.norm(0.0, 1.0))],
            hyperpriors=hyperpriors(),
            over=object_labels(),
            layout="flat",
        )
        with pytest.raises(ParameterError, match="layout='flat' has no internal members"):
            ParameterSet.merge(objects(), populations=[internal_only])


# ---------------------------------------------------------------------------
# On every registered backend
# ---------------------------------------------------------------------------


def scalar(value: object) -> float:
    return float(np.asarray(value, dtype=float))


def build_problem(backend: ConformanceBackend, population: Population) -> FittingProblem:
    """One ``LINEAR`` model and one dataset per member of *population*, on *backend*."""
    spec = ProblemSpec(model=ModelSpec(kind=ModelKind.LINEAR))
    models: dict[str, object] = {}
    datasets = []
    for index, label in enumerate(population.over):
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
    return FittingProblem(models, datasets, populations=[population], seed=spec.seed)


def derived_scale_population() -> Population:
    """Row 8's population: a derived hyperprior, a prior over it, a derived element.

    ``sigma`` is derived (``exp(log_sigma)``); ``offset``'s hierarchical prior
    references it; ``slope``, the ``LINEAR`` model's own parameter, is the
    derived ``mu + sigma * z`` over the internal ``z``.
    """
    return Population(
        "objects",
        members=[
            Parameter("z", st.norm(0.0, 1.0)),
            Parameter("slope", Derived("mu + sigma * z")),
            Parameter("offset", HierarchicalPrior("norm", {"scale": "sigma"}, kwds={"loc": 0.0})),
        ],
        hyperpriors=[
            Parameter("mu", st.norm(1.0, 0.5)),
            Parameter("log_sigma", st.norm(-1.0, 0.5)),
            Parameter("sigma", Derived("exp(log_sigma)")),
        ],
        over=object_labels(),
    )


def chained_population() -> Population:
    """Row 10's chain: derived ``slope`` -> hierarchical ``h`` -> derived ``spread``.

    Declared in reverse, so the evaluation order is the only thing that can
    put each value after its inputs (the §3.4 hazard).
    """
    return Population(
        "objects",
        members=[
            Parameter("slope", Derived("a + h ** 2")),
            Parameter("h", HierarchicalPrior("halfnorm", {"scale": "spread"})),
        ],
        hyperpriors=[Parameter("spread", Derived("exp(a)")), Parameter("a", st.norm(-0.5, 0.3))],
        over=object_labels(),
    )


def unit_cubes(size: int, count: int = 8) -> np.ndarray:
    rng = np.random.default_rng(20261008)
    return np.concatenate(
        [np.full((1, size), 0.5), rng.uniform(0.05, 0.95, size=(count - 1, size))]
    )


def realised(problem: FittingProblem) -> object:
    if problem.backend not in registered_realisations():
        pytest.skip(
            f"the {problem.backend!r} backend registers no realisation "
            f"(inference.md §10a: the reference backend has no differentiable path)"
        )
    return realise(problem)


class TestOnEveryBackend:
    """Rows 1, 2 and 5 through each backend's lowered parameter space."""

    def test_a_derived_parameter_is_not_a_free_dimension(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        pset = non_centred_set()
        space = backend.parameter_space(pset)
        assert space.free_size == pset.free_size == 5
        assert space.free_labels() == pset.free_labels()
        vector = np.array([0.5, 2.0, -1.0, 0.0, 1.5])
        unpacked = space.unpack(vector)
        assert set(unpacked) == set(pset.names)
        np.testing.assert_allclose(
            np.asarray(unpacked["theta"], dtype=float),
            pset.unpack(vector)["theta"],
            atol=tolerances.cross_backend,
        )
        np.testing.assert_allclose(space.pack(unpacked), vector, atol=tolerances.cross_backend)

    def test_complete_is_idempotent_and_computes_it(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """A stale derived value: refused on the numpy path, overwritten on a traced one."""
        pset = non_centred_set()
        space = backend.parameter_space(pset)
        vector = np.array([1.0, 0.5, 0.0, 2.0, -2.0])
        fresh = pset.unpack(vector)
        stale = dict(fresh, theta=np.array([9.0, 9.0, 9.0]))
        # Every path: the free vector has no derived slot, so the derived value
        # a lowering hands back is always the recomputation.
        np.testing.assert_allclose(
            np.asarray(space.unpack(space.pack(stale))["theta"], dtype=float),
            [1.0, 2.0, 0.0],
            atol=tolerances.cross_backend,
        )
        if not backend.capabilities.differentiable:
            with pytest.raises(ParameterError, match="derived parameter 'theta'"):
                space.lnprior(stale)
        elif hasattr(space, "log_prior_tensor"):
            # A traced path that takes a mapping cannot branch on a value: it
            # overwrites the supplied derived entry with the recomputation.
            assert scalar(space.lnprior(stale)) == pytest.approx(
                pset.lnprior(fresh), abs=tolerances.cross_backend
            )

    def test_a_hierarchical_prior_may_reference_a_derived_parameter(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        pset = slab_set()
        space = backend.parameter_space(pset)
        for cube in unit_cubes(pset.free_size):
            expected = pset.prior_transform(cube)
            np.testing.assert_allclose(
                space.prior_transform(cube), expected, rtol=1e-9, atol=tolerances.cross_backend
            )
            assert scalar(space.lnprior(expected)) == pytest.approx(
                pset.lnprior(expected), abs=tolerances.cross_backend
            )
            assert scalar(space.unpack(expected)["s_eff"]) == pytest.approx(
                pset.unpack(expected)["s_eff"], abs=tolerances.cross_backend
            )


class TestTheRealisedDerived:
    """Rows 8 and 10: the lowerings agree with the numpy path."""

    def test_the_realised_derived_agrees_with_the_numpy_path(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_problem(backend, derived_scale_population())
        merged = problem.parameters
        assert merged.derived_names == ("objects.sigma", "objects.slope")
        assert merged["objects.offset"].references == ("objects.sigma",)
        # The model consumes the derived element under its own name.
        values = merged.complete(problem.sample_prior(np.random.default_rng(5)))
        routed = problem.mapping.distribute(values)
        for index, label in enumerate(object_labels()):
            assert scalar(routed[label]["slope"]) == pytest.approx(
                float(values["objects.mu"])
                + float(values["objects.sigma"]) * values["objects.z"][index]
            )
        space = backend.parameter_space(merged)
        if hasattr(space, "numpyro_model"):
            import numpyro

            with numpyro.handlers.seed(rng_seed=0):
                trace = numpyro.handlers.trace(space.numpyro_model()).get_trace()
            for name in merged.derived_names:
                assert trace[name]["type"] == "deterministic"
            assert trace["objects.z"]["type"] == "sample"
        realisation = realised(problem)
        assert int(realisation.free_size) == problem.free_size  # type: ignore[attr-defined]
        for cube in unit_cubes(problem.free_size):
            y = problem.unconstrain(problem.prior_transform(cube))
            expected = problem.log_prob_unconstrained(y)
            got = scalar(backend.to_numpy(realisation.log_prob_unconstrained(y)))  # type: ignore[attr-defined]
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)

    def test_the_evaluation_order_is_honoured_under_tracing(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        problem = build_problem(backend, chained_population())
        merged = problem.parameters
        order = merged.evaluation_order()
        assert (
            order.index("objects.a")
            < order.index("objects.spread")
            < order.index("objects.h")
            < order.index("objects.slope")
        )
        space = backend.parameter_space(merged)
        for cube in unit_cubes(merged.free_size):
            vector = merged.prior_transform(cube)
            np.testing.assert_allclose(
                space.prior_transform(cube), vector, rtol=1e-9, atol=tolerances.cross_backend
            )
            numpy_values = merged.unpack(vector)
            lowered_values = space.unpack(vector)
            for name in merged.names:
                np.testing.assert_allclose(
                    np.asarray(lowered_values[name], dtype=float),
                    numpy_values[name],
                    atol=tolerances.cross_backend,
                    rtol=1e-12,
                )
            a = float(numpy_values["objects.a"])
            h = numpy_values["objects.h"]
            np.testing.assert_allclose(numpy_values["objects.slope"], a + h**2, rtol=1e-14)
            assert scalar(space.lnprior(vector)) == pytest.approx(
                merged.lnprior(vector), abs=tolerances.cross_backend
            )
        realisation = realised(problem)
        for cube in unit_cubes(problem.free_size):
            y = problem.unconstrain(problem.prior_transform(cube))
            expected = problem.log_prob_unconstrained(y)
            got = scalar(backend.to_numpy(realisation.log_prob_unconstrained(y)))  # type: ignore[attr-defined]
            assert got == pytest.approx(expected, abs=tolerances.cross_backend)


#: The three amplitudes row 6 shrinks.
AMPLITUDES = ("broad.amplitude", "narrow.amplitude", "mid.amplitude")

#: sha256 of ``json.dumps(ParameterSet(shrinkage_horseshoe(AMPLITUDES,
#: tail=...)).to_spec(), sort_keys=True)``, recorded from the base commit
#: ``5963b2f`` (before W7.0) — the plain and regularised tails are unchanged
#: byte for byte by the slab's arrival.
UNCHANGED_TAILS = {
    "regularised": "5e4921d65b8067850aa9bd291de90fbeee187287318652009b4bff3b9b065ae0",
    "cauchy": "13cc7e4a2a6cdc0ae5061e08e0367729e6ebf7e4c221b14075cb660e730d82e5",
}


class TestTheSlab:
    """Row 6: ``tail="slab"`` is Piironen & Vehtari's, over the helper's own ``s_j``."""

    def test_the_slab_is_piironen_and_vehtari_in_the_helpers_parameterisation(self) -> None:
        c = 1.5
        declaration = shrinkage_horseshoe(AMPLITUDES, tail="slab", slab_scale=c)
        assert [parameter.name for parameter in declaration] == [
            "shrinkage.global_scale",
            "shrinkage.slab_scale",
            "shrinkage.broad",
            "shrinkage.narrow",
            "shrinkage.mid",
            "shrinkage.broad_effective",
            "shrinkage.narrow_effective",
            "shrinkage.mid_effective",
            *AMPLITUDES,
        ]
        pset = ParameterSet(declaration)
        assert pset["shrinkage.slab_scale"].is_fixed
        assert pset.free_names == (
            "shrinkage.global_scale",
            "shrinkage.broad",
            "shrinkage.narrow",
            "shrinkage.mid",
            *AMPLITUDES,
        )
        assert pset["broad.amplitude"].references == ("shrinkage.broad_effective",)
        rng = np.random.default_rng(20261009)
        for _ in range(6):
            values = pset.sample(rng)
            tau = float(values["shrinkage.global_scale"])
            expected = float(st.halfcauchy(0.0, 1.0).logpdf(tau))
            for leaf in ("broad", "narrow", "mid"):
                s = float(values[f"shrinkage.{leaf}"])
                scale = effective_scale(c, s)
                assert values[f"shrinkage.{leaf}_effective"] == pytest.approx(scale, rel=1e-14)
                expected += float(st.halfcauchy(0.0, tau).logpdf(s))
                expected += float(st.halfnorm(0.0, scale).logpdf(values[f"{leaf}.amplitude"]))
            assert pset.lnprior(values) == pytest.approx(expected, rel=1e-12)
        # The plain and regularised tails, byte for byte as before W7.0.
        for tail, digest in UNCHANGED_TAILS.items():
            spec = ParameterSet(shrinkage_horseshoe(AMPLITUDES, tail=tail)).to_spec()
            text = json.dumps(spec, sort_keys=True)
            assert hashlib.sha256(text.encode()).hexdigest() == digest, tail

    def test_the_slab_scale_is_refused_outside_the_slab_tail(self) -> None:
        with pytest.raises(ParameterError, match="needs slab_scale"):
            shrinkage_horseshoe(AMPLITUDES, tail="slab")
        with pytest.raises(ParameterError, match="used only by tail='slab'"):
            shrinkage_horseshoe(AMPLITUDES, slab_scale=1.0)
        with pytest.raises(ParameterError, match="positive, finite"):
            shrinkage_horseshoe(AMPLITUDES, tail="slab", slab_scale=0.0)

    def test_a_slab_scale_with_a_prior_lowers_and_evaluates(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        declaration = shrinkage_horseshoe(
            AMPLITUDES, tail="slab", slab_scale=Parameter("c", st.halfnorm(0.0, 2.0))
        )
        pset = ParameterSet(declaration)
        assert "shrinkage.slab_scale" in pset.free_names
        space = backend.parameter_space(pset)
        for cube in unit_cubes(pset.free_size, count=5):
            vector = pset.prior_transform(cube)
            np.testing.assert_allclose(
                space.prior_transform(cube), vector, rtol=1e-9, atol=tolerances.cross_backend
            )
            assert scalar(space.lnprior(vector)) == pytest.approx(
                pset.lnprior(vector), abs=tolerances.cross_backend
            )
            values = pset.unpack(vector)
            c = float(values["shrinkage.slab_scale"])
            assert values["shrinkage.mid_effective"] == pytest.approx(
                effective_scale(c, float(values["shrinkage.mid"])), rel=1e-14
            )

    def test_the_slab_rides_on_a_kernel_through_with_shrinkage(self) -> None:
        from ampere.core import Matern32, Sum, with_shrinkage

        kernel = Sum(Matern32(0.3, 0.01), Matern32(0.3, 0.001), labels=("broad", "narrow"))
        shrunk = with_shrinkage(
            kernel,
            shrinkage_horseshoe(
                ("broad.amplitude", "narrow.amplitude"), tail="slab", slab_scale=1.0
            ),
        )
        values = shrunk.parameters.sample(np.random.default_rng(1))
        context = shrunk.context(values)
        assert context["shrinkage.broad_effective"] == pytest.approx(
            effective_scale(1.0, float(values["shrinkage.broad"]))
        )


def results_population() -> Population:
    """Row 9's population: every broadcasting case §3.6 names, and a chain.

    ``slope`` (a plate member, ``(N,)``) is derived from scalar ``mu``, the
    scalar derived ``sigma`` (so a derived parameter of a derived one) and the
    internal ``z``; ``tilt`` is a non-plate ``(3,)`` derived from ``w`` beside
    the scalar ``mu``.
    """
    return Population(
        "objects",
        members=[
            Parameter("z", st.norm(0.0, 1.0)),
            Parameter("slope", Derived("mu + sigma * z")),
            Parameter("offset", HierarchicalPrior("norm", {"scale": "sigma"}, kwds={"loc": 0.0})),
        ],
        hyperpriors=[
            Parameter("mu", st.norm(1.0, 0.5)),
            Parameter("log_sigma", st.norm(-1.0, 0.5)),
            Parameter("sigma", Derived("exp(log_sigma)")),
            Parameter("w", st.norm(0.0, 1.0), shape=(3,)),
            Parameter("tilt", Derived("w * mu"), shape=(3,)),
        ],
        over=object_labels(),
    )


class TestThePosterior:
    """Row 9: every derived parameter is a posterior variable (dev, ``ampere.results``)."""

    def test_the_posterior_carries_the_derived_variable(self) -> None:
        pytest.importorskip("emcee")
        from ampere.inference import EmceeEngine
        from ampere.results import PROVENANCE_SCHEMA_VERSION

        from .backends.reference import ReferenceBackend

        problem = build_problem(ReferenceBackend(), results_population())
        assert PROVENANCE_SCHEMA_VERSION == 10
        run = EmceeEngine(problem, walkers=2 * problem.free_size + 2).run(steps=40, burn_in=2)
        assert run.attrs["ampere_schema_version"] == 10
        derived = json.loads(run.attrs["ampere_derived"])
        assert derived == list(problem.parameters.derived_names)
        assert set(derived) == {"objects.sigma", "objects.slope", "objects.tilt"}
        posterior = run["posterior"]
        chains, draws = np.asarray(posterior["objects.mu"]).shape
        assert np.asarray(posterior["objects.slope"]).shape == (chains, draws, MEMBERS)
        assert np.asarray(posterior["objects.tilt"]).shape == (chains, draws, 3)
        assert np.asarray(posterior["objects.sigma"]).shape == (chains, draws)
        merged = problem.parameters
        for chain in range(chains):
            for draw in range(0, draws, 7):
                values = {
                    name: np.asarray(posterior[name])[chain, draw] for name in merged.free_names
                }
                expected = merged.complete(values)
                for name in derived:
                    np.testing.assert_allclose(
                        np.asarray(posterior[name])[chain, draw], expected[name], rtol=1e-14
                    )

    def test_a_problem_without_one_records_an_empty_list(self) -> None:
        from ampere.results import provenance_attrs

        from .backends.reference import ReferenceBackend

        slopes = Population(
            "objects",
            members=[
                Parameter("slope", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))
            ],
            hyperpriors=hyperpriors(),
            over=object_labels(),
        )
        problem = build_problem(ReferenceBackend(), slopes)
        attrs = provenance_attrs(problem, engine="test")
        assert attrs["ampere_derived"] == "[]"
