"""Unit tests for the W1.3 parameter and prior contract.

Organised by the acceptance criteria in ``WORK_ITEMS.md`` W1.3 and the
dispatch brief: round-trip serialisation, tying reducing the free-dimension
count, plate-aware hierarchical grouping, and the parameter/buffer
distinction — followed by the supporting behaviour each of those rests on.
"""

from __future__ import annotations

import json
import math

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Buffer,
    BufferSet,
    HierarchicalPrior,
    Identity,
    Log,
    Logit,
    Parameter,
    Parameterised,
    ParameterMapping,
    ParameterSet,
    Plate,
    PlateBinding,
    PriorSpec,
    Tie,
    default_bijection_for,
    describe_prior,
    log_density,
    prior_from_spec,
)
from ampere.core.exceptions import (
    OptionalDependencyError,
    ParameterError,
    TyingError,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def mixed_set() -> ParameterSet:
    """A set exercising every layout feature: scalar, array, fixed, units."""
    return ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
            Parameter("log_tau", st.norm(0.0, 1.0)),
            Parameter("offset", st.norm(0.0, 0.05), shape=(3,), unit=u.mag),
            Parameter("distance", value=1.5 * u.kpc, fixed=True, unit=u.kpc),
        ]
    )


def _tied_pair() -> tuple[ParameterSet, ParameterSet]:
    """Two sets with five parameters between them, one pair tied."""
    sed = ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
            Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance"),
        ]
    )
    spectrum = ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
            Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance"),
            Parameter("scale", st.norm(1.0, 0.05)),
        ]
    )
    return sed, spectrum


# ===========================================================================
# Acceptance criterion: round-trip prior serialisation
# ===========================================================================


class TestRoundTrip:
    def test_pack_unpack_is_the_identity(self, mixed_set: ParameterSet) -> None:
        theta = mixed_set.prior_transform(np.linspace(0.1, 0.9, mixed_set.free_size))
        values = mixed_set.unpack(theta)
        assert np.array_equal(mixed_set.pack(values), theta)

    def test_unpack_yields_all_parameters_including_fixed(self, mixed_set: ParameterSet) -> None:
        values = mixed_set.unpack(np.zeros(mixed_set.free_size))
        assert set(values) == set(mixed_set.names)
        assert values["distance"] == 1.5  # injected, not sampled
        assert mixed_set.free_size == 5 and len(mixed_set) == 4

    def test_unpack_restores_array_shape(self, mixed_set: ParameterSet) -> None:
        values = mixed_set.unpack(np.arange(5, dtype=float))
        assert isinstance(values["offset"], np.ndarray)
        assert values["offset"].shape == (3,)
        assert np.array_equal(values["offset"], np.array([2.0, 3.0, 4.0]))

    def test_spec_round_trip_preserves_names_values_and_priors(
        self, mixed_set: ParameterSet
    ) -> None:
        payload = json.loads(json.dumps(mixed_set.to_spec()))  # proves it is JSON-able
        rebuilt = ParameterSet.from_spec(payload)

        assert rebuilt.names == mixed_set.names
        assert rebuilt.free_names == mixed_set.free_names
        assert rebuilt.free_size == mixed_set.free_size
        assert rebuilt == mixed_set

        for original in mixed_set:
            copy = rebuilt[original.name]
            assert copy.unit == original.unit
            assert copy.shape == original.shape
            assert copy.fixed == original.fixed
            if original.prior is None:
                assert copy.prior is None
            else:
                assert describe_prior(copy.prior) == describe_prior(original.prior)

        theta = mixed_set.prior_transform(np.full(mixed_set.free_size, 0.3))
        assert np.allclose(rebuilt.prior_transform(np.full(5, 0.3)), theta)
        assert math.isclose(rebuilt.lnprior(theta), mixed_set.lnprior(theta))

    def test_spec_round_trip_survives_bijections_ties_and_hierarchy(self) -> None:
        original = ParameterSet(
            [
                Parameter("mu", st.norm(0.0, 5.0), shared_as="mu", description="population mean"),
                Parameter("sigma", st.halfnorm(0.0, 2.0), bijection=Log(lower=0.0)),
                Parameter(
                    "theta",
                    HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}),
                    shape=(3,),
                    plate="objects",
                ),
            ]
        )
        rebuilt = ParameterSet.from_spec(json.loads(json.dumps(original.to_spec())))
        assert rebuilt == original
        assert rebuilt["mu"].shared_as == "mu"
        assert rebuilt["mu"].description == "population mean"
        assert rebuilt["sigma"].bijection == Log(lower=0.0)
        assert rebuilt["theta"].plate == "objects"
        assert rebuilt["theta"].references == ("mu", "sigma")
        assert rebuilt.plates == {"objects": 3}

    def test_prior_spec_round_trip(self) -> None:
        for prior in (
            st.norm(0.0, 1.0),
            st.uniform(loc=100.0, scale=9900.0),
            st.loguniform(1e-3, 1e3),
            st.poisson(3.0),
        ):
            spec = describe_prior(prior)
            assert PriorSpec.from_dict(json.loads(json.dumps(spec.to_dict()))) == spec
            rebuilt = prior_from_spec(spec)
            assert np.allclose(rebuilt.ppf(0.42), prior.ppf(0.42))
            assert np.allclose(log_density(rebuilt, 1), log_density(prior, 1))

    def test_discrete_priors_are_marked_and_evaluated(self) -> None:
        spec = describe_prior(st.poisson(3.0))
        assert spec.discrete is True and spec.family == "poisson"
        pset = ParameterSet([Parameter("counts", st.poisson(3.0))])
        assert math.isclose(pset.lnprior(np.array([2.0])), float(st.poisson(3.0).logpmf(2)))

    def test_opaque_prior_evaluates_but_refuses_to_serialise(self) -> None:
        class Custom:
            def logpdf(self, x: float) -> float:
                return -0.5 * float(x) ** 2

            def ppf(self, q):
                return st.norm(0.0, 1.0).ppf(q)

            def support(self):
                return (-np.inf, np.inf)

        pset = ParameterSet([Parameter("x", Custom())])
        assert math.isclose(pset.lnprior(np.array([1.0])), -0.5)
        with pytest.raises(ParameterError, match=r"not a frozen scipy\.stats distribution"):
            pset.to_spec()


# ===========================================================================
# Acceptance criterion: tying reduces the free-dimension count
# ===========================================================================


class TestTying:
    def test_five_parameters_one_tied_pair_gives_four_free_dimensions(self) -> None:
        sed, spectrum = _tied_pair()
        assert sed.free_size + spectrum.free_size == 5

        mapping = ParameterSet.merge({"sed": sed, "spectrum": spectrum})

        assert mapping.merged.free_size == 4
        assert len(mapping.merged) == 4
        assert mapping.merged.names == (
            "sed.temperature",
            "distance",
            "spectrum.temperature",
            "spectrum.scale",
        )
        assert mapping.tied_names == ("distance",)

    def test_tied_parameter_has_two_binding_sites(self) -> None:
        sed, spectrum = _tied_pair()
        mapping = ParameterSet.merge({"sed": sed, "spectrum": spectrum})
        sites = mapping.sites_of("distance")
        assert {(b.component, b.local_name) for b in sites} == {
            ("sed", "distance"),
            ("spectrum", "distance"),
        }
        assert len(mapping.sites_of("spectrum.scale")) == 1
        assert mapping.global_name_for("spectrum", "distance") == "distance"

    def test_distribute_routes_one_value_to_every_binding_site(self) -> None:
        sed, spectrum = _tied_pair()
        mapping = ParameterSet.merge({"sed": sed, "spectrum": spectrum})
        theta = mapping.merged.prior_transform(np.full(4, 0.5))
        routed = mapping.distribute(theta)

        assert set(routed) == {"sed", "spectrum"}
        assert set(routed["sed"]) == {"temperature", "distance"}
        assert set(routed["spectrum"]) == {"temperature", "distance", "scale"}
        assert routed["sed"]["distance"] == routed["spectrum"]["distance"]
        # Untied same-named parameters stay independent.
        merged_values = mapping.merged.unpack(theta)
        assert routed["sed"]["temperature"] == merged_values["sed.temperature"]
        assert routed["spectrum"]["temperature"] == merged_values["spectrum.temperature"]

    def test_distribute_accepts_a_mapping_as_well_as_a_vector(self) -> None:
        sed, spectrum = _tied_pair()
        mapping = ParameterSet.merge({"sed": sed, "spectrum": spectrum})
        theta = mapping.merged.prior_transform(np.full(4, 0.5))
        assert mapping.distribute(theta) == mapping.distribute(mapping.merged.unpack(theta))

    def test_composition_time_tie_is_equivalent_to_declaration_time(self) -> None:
        def untied() -> ParameterSet:
            return ParameterSet(
                [
                    Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
                    Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc),
                ]
            )

        mapping = ParameterSet.merge(
            {"a": untied(), "b": untied()},
            ties=[Tie("distance", ("a.distance", "b.distance"))],
        )
        assert mapping.merged.free_size == 3
        assert mapping.tied_names == ("distance",)
        assert len(mapping.sites_of("distance")) == 2

    def test_tie_can_override_the_prior(self) -> None:
        def untied() -> ParameterSet:
            return ParameterSet([Parameter("d", st.norm(1.5, 0.1))])

        mapping = ParameterSet.merge(
            {"a": untied(), "b": untied()},
            ties=[Tie("d", ("a.d", "b.d"), prior=st.uniform(0.0, 10.0))],
        )
        assert describe_prior(mapping.merged["d"].prior).family == "uniform"

    def test_merging_without_ties_only_qualifies_names(self) -> None:
        sed, _spectrum = _tied_pair()
        plain = ParameterSet([p.rename(p.name) for p in sed if p.shared_as is None])
        mapping = ParameterSet.merge({"one": plain, "two": plain})
        assert mapping.merged.names == ("one.temperature", "two.temperature")
        assert mapping.merged.free_size == 2
        assert mapping.tied_names == ()

    def test_fixed_parameters_may_be_tied(self) -> None:
        def fixed_set() -> ParameterSet:
            return ParameterSet([Parameter("d", value=1.5, fixed=True, unit=u.kpc, shared_as="d")])

        mapping = ParameterSet.merge({"a": fixed_set(), "b": fixed_set()})
        assert mapping.merged.free_size == 0
        assert mapping.merged["d"].is_fixed
        assert mapping.distribute(np.array([]))["a"]["d"] == 1.5

    def test_deferred_site_gets_its_prior_from_the_group(self) -> None:
        authority = ParameterSet(
            [Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance")]
        )
        consumer = ParameterSet([Parameter("distance", unit=u.kpc, shared_as="distance")])

        assert authority.is_resolved
        assert not consumer.is_resolved
        assert consumer.deferred_names == ("distance",)
        with pytest.raises(TyingError, match="Merge this set"):
            consumer.lnprior(np.array([1.5]))
        # Layout still works while unresolved: pack/unpack are pure layout.
        assert consumer.unpack(np.array([1.5])) == {"distance": 1.5}

        mapping = ParameterSet.merge({"one": authority, "two": consumer})
        assert mapping.merged.is_resolved
        assert mapping.merged.free_size == 1
        assert describe_prior(mapping.merged["distance"].prior).family == "norm"

    @pytest.mark.parametrize(
        ("left", "right", "match"),
        [
            (
                Parameter("d", st.norm(1.0, 1.0), unit=u.kpc, shared_as="d"),
                Parameter("d", st.norm(1.0, 1.0), unit=u.pc, shared_as="d"),
                "disagree about units",
            ),
            (
                Parameter("d", st.norm(1.0, 1.0), shared_as="d"),
                Parameter("d", st.norm(2.0, 1.0), shared_as="d"),
                "disagreeing priors",
            ),
            (
                Parameter("d", st.norm(1.0, 1.0), shape=(2,), shared_as="d"),
                Parameter("d", st.norm(1.0, 1.0), shape=(3,), shared_as="d"),
                "disagree about shape",
            ),
            (
                Parameter("d", value=1.0, fixed=True, shared_as="d"),
                Parameter("d", value=2.0, fixed=True, shared_as="d"),
                "different values",
            ),
            (
                Parameter("d", value=1.0, fixed=True, shared_as="d"),
                Parameter("d", st.norm(1.0, 1.0), shared_as="d"),
                "mixes fixed and prior-equipped",
            ),
        ],
    )
    def test_incompatible_ties_fail_loudly(
        self, left: Parameter, right: Parameter, match: str
    ) -> None:
        with pytest.raises(TyingError, match=match):
            ParameterSet.merge({"a": ParameterSet([left]), "b": ParameterSet([right])})

    def test_tie_naming_a_missing_site_fails(self) -> None:
        pset = ParameterSet([Parameter("d", st.norm(1.0, 1.0))])
        with pytest.raises(TyingError, match="does not exist"):
            ParameterSet.merge({"a": pset}, ties=[Tie("z", ("a.d", "a.nope"))])

    def test_a_site_may_not_belong_to_two_ties(self) -> None:
        pset = ParameterSet([Parameter("d", st.norm(1.0, 1.0)), Parameter("e", st.norm(1.0, 1.0))])
        with pytest.raises(TyingError, match="at most one tie"):
            ParameterSet.merge(
                {"a": pset},
                ties=[Tie("x", ("a.d", "a.e")), Tie("y", ("a.d", "a.e"))],
            )

    def test_shared_as_and_an_explicit_tie_on_one_site_is_ambiguous(self) -> None:
        pset = ParameterSet(
            [Parameter("d", st.norm(1.0, 1.0), shared_as="d"), Parameter("e", st.norm(1.0, 1.0))]
        )
        with pytest.raises(TyingError, match="ambiguous"):
            ParameterSet.merge({"a": pset}, ties=[Tie("z", ("a.d", "a.e"))])

    def test_a_tie_needs_at_least_two_sites(self) -> None:
        with pytest.raises(TyingError, match="collapses two or more"):
            Tie("d", ("a.d",))

    def test_tie_group_with_no_prior_anywhere_fails(self) -> None:
        deferred = ParameterSet([Parameter("d", shared_as="d")])
        with pytest.raises(TyingError, match="has no prior"):
            ParameterSet.merge({"a": deferred, "b": deferred})

    def test_positional_and_keyword_prior_declarations_tie_cleanly(self) -> None:
        # scipy accepts the same freezing either way; describe_prior
        # canonicalises, so the two declarations are one prior.
        positional = ParameterSet([Parameter("d", st.norm(1.5, 0.1), shared_as="d")])
        keyword = ParameterSet([Parameter("d", st.norm(loc=1.5, scale=0.1), shared_as="d")])
        assert describe_prior(st.norm(1.5, 0.1)) == describe_prior(st.norm(loc=1.5, scale=0.1))
        mapping = ParameterSet.merge({"a": positional, "b": keyword})
        assert mapping.merged.free_size == 1

    def test_tied_hierarchical_priors_must_share_their_hyperparameters(self) -> None:
        # Each component's theta references its *own* mu; collapsing the two
        # thetas would have to pick one component's mu over the other's.
        def component() -> ParameterSet:
            return ParameterSet(
                [
                    Parameter("mu", st.norm(0.0, 5.0)),
                    Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}), shared_as="theta"),
                ]
            )

        with pytest.raises(TyingError, match=r"compared after qualification"):
            ParameterSet.merge({"a": component(), "b": component()})

    def test_tied_hierarchical_priors_with_shared_hyperparameters_collapse(self) -> None:
        def component() -> ParameterSet:
            return ParameterSet(
                [
                    Parameter("mu", st.norm(0.0, 5.0), shared_as="mu"),
                    Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}), shared_as="theta"),
                ]
            )

        mapping = ParameterSet.merge({"a": component(), "b": component()})
        assert mapping.merged.names == ("mu", "theta")
        assert mapping.merged["theta"].references == ("mu",)
        assert mapping.merged.free_size == 2

    def test_tied_sites_may_not_declare_conflicting_bijections(self) -> None:
        left = ParameterSet([Parameter("s", st.uniform(0.0, 1.0), shared_as="s", bijection=Log())])
        right = ParameterSet(
            [Parameter("s", st.uniform(0.0, 1.0), shared_as="s", bijection=Logit(0.0, 1.0))]
        )
        with pytest.raises(TyingError, match="different bijections"):
            ParameterSet.merge({"a": left, "b": right})
        # One site declaring for the group is fine and wins over inference.
        mapping = ParameterSet.merge(
            {"a": left, "b": ParameterSet([Parameter("s", st.uniform(0.0, 1.0), shared_as="s")])}
        )
        assert mapping.merged["s"].bijection == Log()


# ===========================================================================
# Acceptance criterion: plate-aware hierarchical grouping
# ===========================================================================


def _population_plate(size: int = 4) -> Plate:
    return Plate(
        "objects",
        size=size,
        hyperparameters=[
            Parameter("mu", st.norm(0.0, 5.0)),
            Parameter("sigma", st.halfnorm(0.0, 2.0)),
        ],
        members=[
            Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
        ],
    )


class TestPlates:
    def test_plate_declares_shared_hyperparameters_and_per_member_values(self) -> None:
        population = ParameterSet(
            [Parameter("background", st.norm(0.0, 1.0))], plates=[_population_plate()]
        )
        assert population.names == (
            "background",
            "objects.mu",
            "objects.sigma",
            "objects.theta",
        )
        # mu and sigma are shared: one dimension each, whatever N is.
        assert population["objects.mu"].shape == ()
        assert population["objects.sigma"].shape == ()
        # theta_i is per-member: one dimension per object.
        assert population["objects.theta"].shape == (4,)
        assert population["objects.theta"].plate == "objects"
        assert population.plates == {"objects": 4}
        assert population.free_size == 1 + 1 + 1 + 4

    def test_member_prior_references_resolve_to_the_plates_hyperparameters(self) -> None:
        population = _population_plate().to_parameter_set()
        theta = population["objects.theta"]
        assert theta.is_hierarchical
        assert theta.references == ("objects.mu", "objects.sigma")
        assert set(theta.references) <= set(population.names)

    def test_declared_structure_is_numerically_sound(self) -> None:
        population = _population_plate().to_parameter_set()
        cube = np.array([0.5, 0.5, 0.1, 0.3, 0.7, 0.9])
        values = population.unpack(population.prior_transform(cube))
        mu, sigma = values["objects.mu"], values["objects.sigma"]
        assert np.allclose(values["objects.theta"], st.norm(mu, sigma).ppf(cube[2:]))

    def test_lnprior_sums_hyperprior_and_member_contributions(self) -> None:
        population = _population_plate(size=3).to_parameter_set()
        values = {
            "objects.mu": 1.0,
            "objects.sigma": 2.0,
            "objects.theta": np.array([0.5, 1.0, 1.5]),
        }
        expected = (
            float(st.norm(0.0, 5.0).logpdf(1.0))
            + float(st.halfnorm(0.0, 2.0).logpdf(2.0))
            + float(np.sum(st.norm(1.0, 2.0).logpdf(values["objects.theta"])))
        )
        assert math.isclose(population.lnprior(values), expected)

    def test_sampling_respects_dependency_order(self) -> None:
        population = _population_plate(size=5).to_parameter_set()
        drawn = population.sample(np.random.default_rng(0))
        assert drawn["objects.theta"].shape == (5,)
        assert np.isscalar(drawn["objects.mu"])

    def test_labels_name_every_plate_member(self) -> None:
        population = _population_plate(size=3).to_parameter_set()
        assert population.free_labels() == (
            "objects.mu",
            "objects.sigma",
            "objects.theta[0]",
            "objects.theta[1]",
            "objects.theta[2]",
        )

    def test_per_object_pattern_shares_hyperparameters_across_components(self) -> None:
        """The other hierarchical layout: one component per object, tied hypers."""

        def one_object() -> ParameterSet:
            return ParameterSet(
                [
                    Parameter("mu", st.norm(0.0, 5.0), shared_as="mu"),
                    Parameter("sigma", st.halfnorm(0.0, 2.0), shared_as="sigma"),
                    Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})),
                ]
            )

        mapping = ParameterSet.merge(
            {"obj1": one_object(), "obj2": one_object(), "obj3": one_object()}
        )
        # 3 objects x 3 parameters = 9 declared, but mu and sigma are shared.
        assert mapping.merged.free_size == 2 + 3
        assert mapping.merged.names == ("mu", "sigma", "obj1.theta", "obj2.theta", "obj3.theta")
        assert mapping.tied_names == ("mu", "sigma")
        for label in ("obj1", "obj2", "obj3"):
            assert mapping.merged[f"{label}.theta"].references == ("mu", "sigma")
        # ... and it still evaluates.
        theta = mapping.merged.prior_transform(np.full(5, 0.6))
        assert math.isfinite(mapping.merged.lnprior(theta))

    def test_a_plate_needs_members(self) -> None:
        with pytest.raises(ParameterError, match="declares no members"):
            Plate("objects", size=3, hyperparameters=[Parameter("mu", st.norm(0.0, 1.0))])

    def test_unresolvable_hierarchical_reference_fails_at_construction(self) -> None:
        with pytest.raises(ParameterError, match="not in this set"):
            ParameterSet([Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))])

    def test_cyclic_hierarchical_references_fail_at_construction(self) -> None:
        with pytest.raises(ParameterError, match="cycle"):
            ParameterSet(
                [
                    Parameter("a", HierarchicalPrior("norm", {"loc": "b"})),
                    Parameter("b", HierarchicalPrior("norm", {"loc": "a"})),
                ]
            )

    def test_hierarchical_prior_needs_at_least_one_reference(self) -> None:
        with pytest.raises(ParameterError, match="declares no hyperparameters"):
            HierarchicalPrior("norm", {})

    def test_tying_across_plates_is_refused_with_an_explanation(self) -> None:
        plated = _population_plate(size=2).to_parameter_set()
        shared = ParameterSet(
            [
                *[p.rename(p.name) for p in plated][:-1],
                Parameter(
                    "objects.theta",
                    HierarchicalPrior("norm", {"loc": "objects.mu", "scale": "objects.sigma"}),
                    shape=(2,),
                    plate="objects",
                    shared_as="theta",
                ),
            ]
        )
        with pytest.raises(TyingError, match=r"out of scope for v1\.3"):
            ParameterSet.merge({"a": shared, "b": shared})


# ===========================================================================
# Acceptance criterion: explicit buffer declaration
# ===========================================================================


class ToyModel(Parameterised):
    """A model-like object declaring one parameter and two buffers."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.micron)
        self.register_buffer("beta", 1.8, description="dust emissivity index")
        self.register_parameter(Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K))

    def __call__(self, **values: object) -> np.ndarray:
        ctx = self.context(values)
        return ctx["temperature"] * ctx["wavelength"] ** -ctx["beta"]


class TestBuffers:
    def test_parameters_and_buffers_are_separate_declarations(self) -> None:
        model = ToyModel(np.array([1.0, 2.0, 4.0]))
        assert model.parameters.names == ("temperature",)
        assert model.buffers.names == ("wavelength", "beta")
        # Buffers are excluded from the parameter set entirely...
        assert "wavelength" not in model.parameters
        assert "beta" not in model.parameters
        # ... and parameters from the buffer set.
        assert "temperature" not in model.buffers

    def test_buffers_take_no_sampler_dimensions_and_no_prior(self) -> None:
        model = ToyModel(np.array([1.0, 2.0, 4.0]))
        assert model.parameters.free_size == 1
        assert model.parameters.free_labels() == ("temperature",)
        theta = np.array([300.0])
        assert model.parameters.unpack(theta) == {"temperature": 300.0}
        # lnprior sees only the parameter.
        assert math.isclose(
            model.parameters.lnprior(theta), float(st.uniform(100.0, 9900.0).logpdf(300.0))
        )
        assert not hasattr(model.buffers, "lnprior")
        assert not hasattr(model.buffers, "prior_transform")

    def test_context_merges_both_into_one_namespace(self) -> None:
        model = ToyModel(np.array([1.0, 2.0, 4.0]))
        ctx = model.context({"temperature": 100.0})
        assert set(ctx) == {"temperature", "wavelength", "beta"}
        assert np.allclose(model(temperature=100.0), 100.0 * np.array([1.0, 2.0, 4.0]) ** -1.8)
        assert np.allclose(model.context(np.array([100.0]))["temperature"], 100.0)

    def test_promotion_is_configuration_not_a_code_change(self) -> None:
        """architecture.md section 6's stated test of this contract."""
        model = ToyModel(np.array([1.0, 2.0, 4.0]))
        before = model(temperature=100.0)

        promoted = model.promote_buffer("beta", prior=st.norm(1.8, 0.2))

        assert promoted.is_free
        assert model.parameters.names == ("temperature", "beta")
        assert model.buffers.names == ("wavelength",)
        assert model.parameters.free_size == 2
        # __call__ was not touched, and gives the same answer at the same value.
        assert np.allclose(model(temperature=100.0, beta=1.8), before)

    def test_promotion_without_a_prior_yields_a_fixed_parameter(self) -> None:
        model = ToyModel(np.array([1.0, 2.0]))
        promoted = model.promote_buffer("beta")
        assert promoted.is_fixed and promoted.value == 1.8
        assert model.parameters.free_size == 1  # still only temperature
        assert "beta" in model.parameters.fixed_names
        assert np.allclose(model(temperature=100.0), 100.0 * np.array([1.0, 2.0]) ** -1.8)

    def test_demotion_is_the_inverse(self) -> None:
        model = ToyModel(np.array([1.0, 2.0]))
        model.promote_buffer("beta", prior=st.norm(1.8, 0.2))
        model.demote_parameter("beta", value=1.8)
        assert model.buffers.names == ("wavelength", "beta")
        assert model.parameters.names == ("temperature",)
        assert np.allclose(model(temperature=100.0), 100.0 * np.array([1.0, 2.0]) ** -1.8)

    def test_buffer_arrays_are_read_only(self) -> None:
        buffer = Buffer("wavelength", np.array([1.0, 2.0, 3.0]))
        with pytest.raises(ValueError, match="read-only"):
            buffer.value[0] = 9.0

    def test_buffers_keep_their_dtype_and_units(self) -> None:
        indices = Buffer("channel", np.arange(4))
        assert indices.dtype == np.dtype(int)
        wl = Buffer("wavelength", np.array([1.0, 2.0]) * u.micron)
        assert wl.unit == u.micron
        assert np.allclose(wl.value, [1.0, 2.0])
        assert wl.quantity().unit == u.micron

    def test_a_name_cannot_be_both(self) -> None:
        model = ToyModel(np.array([1.0]))
        with pytest.raises(ParameterError, match="already declares a buffer"):
            model.register_parameter(Parameter("wavelength", st.norm(0.0, 1.0)))
        with pytest.raises(ParameterError, match="already declares a parameter"):
            model.register_buffer("temperature", np.array([1.0]))

    def test_a_declaration_cannot_shadow_the_mixins_api(self) -> None:
        model = ToyModel(np.array([1.0]))
        with pytest.raises(ParameterError, match="shadows an attribute"):
            model.register_buffer("context", np.array([1.0]))

    def test_buffer_set_is_immutable_and_ordered(self) -> None:
        original = BufferSet([Buffer("a", np.array([1.0]))])
        extended = original.with_buffer(Buffer("b", np.array([2.0])))
        assert original.names == ("a",)
        assert extended.names == ("a", "b")
        assert extended.without("a").names == ("b",)
        assert list(extended)[1].name == "b"
        assert set(extended.values()) == {"a", "b"}

    def test_buffers_are_not_inferred_from_stray_attributes(self) -> None:
        """Explicit declaration (architecture.md section 6): nothing is adopted."""

        class Sloppy(Parameterised):
            def __init__(self) -> None:
                self.opacity = np.array([1.0, 2.0, 3.0])  # never registered
                self.register_parameter(Parameter("scale", st.norm(1.0, 0.1)))

        model = Sloppy()
        assert model.buffers.names == ()
        assert "opacity" not in model.buffers
        assert "opacity" not in model.context({"scale": 1.0})


# ===========================================================================
# Supporting behaviour
# ===========================================================================


class TestParameterStates:
    def test_fixed_and_prior_together_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="delta-function prior is not the same"):
            Parameter("t", st.norm(0.0, 1.0), value=1.0, fixed=True)

    def test_fixed_without_a_value_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="defined by its value"):
            Parameter("t", fixed=True)

    def test_neither_prior_nor_fixed_nor_shared_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="nothing determines it"):
            Parameter("t")

    def test_fix_and_release_round_trip(self) -> None:
        free = Parameter("t", st.uniform(100.0, 9900.0), value=2500.0, unit=u.K)
        held = free.fix()
        assert held.is_fixed and held.value == 2500.0 and held.prior is None
        released = held.release(st.uniform(100.0, 9900.0))
        assert released.is_free
        assert describe_prior(released.prior) == describe_prior(free.prior)

    @pytest.mark.parametrize("name", ["2x", "a b", "", "a-b", "a..b"])
    def test_unusable_names_are_refused(self, name: str) -> None:
        with pytest.raises(ParameterError):
            Parameter(name, st.norm(0.0, 1.0))

    def test_qualified_names_are_allowed_but_not_declarable_as_ties(self) -> None:
        assert Parameter("objects.mu", st.norm(0.0, 1.0)).name == "objects.mu"
        with pytest.raises(ParameterError, match="must not contain"):
            Parameter("x", st.norm(0.0, 1.0), shared_as="a.b")

    def test_equality_compares_priors_by_description(self) -> None:
        left = Parameter("t", st.norm(0.0, 1.0), unit=u.K)
        right = Parameter("t", st.norm(0.0, 1.0), unit=u.K)
        assert left is not right
        assert left == right
        assert left != Parameter("t", st.norm(0.0, 2.0), unit=u.K)
        assert left != Parameter("t", st.norm(0.0, 1.0), unit=u.Jy)


class TestUnits:
    def test_quantities_are_converted_once_at_declaration(self) -> None:
        p = Parameter("distance", value=1500.0 * u.pc, fixed=True, unit=u.kpc)
        assert p.value == 1.5
        assert not isinstance(p.value, u.Quantity)
        assert p.quantity() == 1.5 * u.kpc

    def test_incompatible_units_fail_loudly(self) -> None:
        with pytest.raises(ParameterError, match="cannot convert"):
            Parameter("distance", value=1500.0 * u.K, fixed=True, unit=u.kpc)

    def test_pack_accepts_quantities_and_stores_plain_floats(self) -> None:
        pset = ParameterSet([Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc)])
        flat = pset.pack({"distance": 2000.0 * u.pc})
        assert flat.dtype == np.float64
        assert np.allclose(flat, [2.0])

    def test_unit_must_be_an_astropy_unit(self) -> None:
        with pytest.raises(ParameterError, match="astropy unit"):
            Parameter("t", st.norm(0.0, 1.0), unit="K")


class TestLayout:
    def test_free_size_and_len_are_different_questions(self, mixed_set: ParameterSet) -> None:
        assert len(mixed_set) == 4
        assert mixed_set.free_size == 5

    def test_free_slices_are_contiguous_and_in_declaration_order(
        self, mixed_set: ParameterSet
    ) -> None:
        assert mixed_set.free_slice("temperature") == slice(0, 1)
        assert mixed_set.free_slice("log_tau") == slice(1, 2)
        assert mixed_set.free_slice("offset") == slice(2, 5)
        with pytest.raises(KeyError, match="fixed"):
            mixed_set.free_slice("distance")

    def test_labels_name_every_flat_dimension(self, mixed_set: ParameterSet) -> None:
        assert mixed_set.free_labels() == (
            "temperature",
            "log_tau",
            "offset[0]",
            "offset[1]",
            "offset[2]",
        )
        assert len(mixed_set.free_labels()) == mixed_set.free_size

    def test_multidimensional_labels_use_c_order(self) -> None:
        pset = ParameterSet([Parameter("m", st.norm(0.0, 1.0), shape=(2, 2))])
        assert pset.free_labels() == ("m[0,0]", "m[0,1]", "m[1,0]", "m[1,1]")

    def test_pack_rejects_unknown_and_missing_names(self, mixed_set: ParameterSet) -> None:
        values = mixed_set.unpack(np.zeros(5))
        with pytest.raises(ParameterError, match="unknown parameter"):
            mixed_set.pack({**values, "temprature": 1.0})
        with pytest.raises(ParameterError, match="missing a value"):
            mixed_set.pack({k: v for k, v in values.items() if k != "log_tau"})

    def test_unpack_rejects_the_wrong_length(self, mixed_set: ParameterSet) -> None:
        with pytest.raises(ParameterError, match="length 5"):
            mixed_set.unpack(np.zeros(4))

    def test_duplicate_names_are_refused(self) -> None:
        with pytest.raises(ParameterError, match="duplicate parameter name"):
            ParameterSet([Parameter("t", st.norm(0.0, 1.0)), Parameter("t", st.norm(0.0, 1.0))])

    def test_sets_are_immutable(self, mixed_set: ParameterSet) -> None:
        extended = mixed_set.with_parameter(Parameter("new", st.norm(0.0, 1.0)))
        assert "new" not in mixed_set
        assert "new" in extended
        assert mixed_set.without("log_tau").names == ("temperature", "offset", "distance")


class TestPriorEvaluation:
    def test_lnprior_matches_summed_scipy_log_densities(self, mixed_set: ParameterSet) -> None:
        theta = np.array([2500.0, 0.5, 0.01, -0.02, 0.03])
        expected = (
            float(st.uniform(100.0, 9900.0).logpdf(2500.0))
            + float(st.norm(0.0, 1.0).logpdf(0.5))
            + float(np.sum(st.norm(0.0, 0.05).logpdf(theta[2:])))
        )
        assert math.isclose(mixed_set.lnprior(theta), expected)

    def test_lnprior_short_circuits_outside_the_support(self, mixed_set: ParameterSet) -> None:
        theta = np.array([50.0, 0.0, 0.0, 0.0, 0.0])
        assert mixed_set.lnprior(theta) == -math.inf

    def test_prior_transform_matches_analytic_quantiles(self) -> None:
        pset = ParameterSet(
            [
                Parameter("a", st.uniform(100.0, 9900.0)),
                Parameter("b", st.loguniform(1e-3, 1e3)),
            ]
        )
        cube = np.array([0.25, 0.75])
        assert np.allclose(
            pset.prior_transform(cube),
            [st.uniform(100.0, 9900.0).ppf(0.25), st.loguniform(1e-3, 1e3).ppf(0.75)],
        )

    def test_prior_transform_rejects_the_wrong_length(self, mixed_set: ParameterSet) -> None:
        with pytest.raises(ParameterError, match="unit-cube vector of length 5"):
            mixed_set.prior_transform(np.zeros(2))

    def test_sampling_is_reproducible_and_complete(self, mixed_set: ParameterSet) -> None:
        first = mixed_set.sample(np.random.default_rng(7))
        second = mixed_set.sample(np.random.default_rng(7))
        assert set(first) == set(mixed_set.names)
        assert first["distance"] == 1.5
        assert math.isclose(first["temperature"], second["temperature"])
        assert np.allclose(first["offset"], second["offset"])

    def test_complete_fills_in_fixed_values(self, mixed_set: ParameterSet) -> None:
        completed = mixed_set.complete({"temperature": 1.0, "log_tau": 0.0, "offset": np.zeros(3)})
        assert completed["distance"] == 1.5
        with pytest.raises(ParameterError, match="no value supplied"):
            mixed_set.complete({"temperature": 1.0})

    def test_a_typo_in_a_value_mapping_is_never_ignored(self, mixed_set: ParameterSet) -> None:
        values = mixed_set.unpack(np.zeros(5))
        with pytest.raises(ParameterError, match="unknown parameter"):
            mixed_set.lnprior({**values, "temprature": 3000.0})
        with pytest.raises(ParameterError, match="unknown parameter"):
            mixed_set.complete({**values, "temprature": 3000.0})


class TestBijections:
    def test_defaults_follow_the_support(self) -> None:
        assert default_bijection_for(st.norm(0.0, 1.0)) == Identity()
        assert default_bijection_for(st.halfnorm(0.0, 1.0)) == Log(lower=0.0)
        assert default_bijection_for(st.uniform(100.0, 9900.0)) == Logit(100.0, 10000.0)

    @pytest.mark.parametrize(
        "bijection", [Identity(), Log(lower=0.0), Log(lower=2.5), Logit(-1.0, 3.0)]
    )
    def test_constrain_and_unconstrain_are_inverse(self, bijection: object) -> None:
        y = np.array([-2.0, -0.3, 0.0, 0.7, 1.9])
        x = bijection.constrain(y)  # type: ignore[attr-defined]
        assert np.allclose(bijection.unconstrain(x), y)  # type: ignore[attr-defined]

    @pytest.mark.parametrize(
        "bijection", [Identity(), Log(lower=0.0), Log(lower=2.5), Logit(-1.0, 3.0)]
    )
    def test_log_det_jacobian_matches_numerical_differentiation(self, bijection: object) -> None:
        y = np.array([-1.3, 0.0, 0.8])
        eps = 1e-6
        numerical = (
            bijection.constrain(y + eps) - bijection.constrain(y - eps)  # type: ignore[attr-defined]
        ) / (2 * eps)
        analytic = np.exp(bijection.log_abs_det_jacobian(y))  # type: ignore[attr-defined]
        assert np.allclose(analytic, np.abs(numerical), rtol=1e-5)

    def test_set_level_round_trip_and_jacobian(self, mixed_set: ParameterSet) -> None:
        theta = mixed_set.prior_transform(np.full(5, 0.3))
        y = mixed_set.unconstrain(theta)
        assert np.allclose(mixed_set.constrain(y), theta)

        expected = mixed_set.lnprior(theta)
        for parameter in mixed_set:
            if parameter.is_fixed:
                continue
            where = mixed_set.free_slice(parameter.name)
            expected += float(
                np.sum(parameter.unconstraining_bijection().log_abs_det_jacobian(y[where]))
            )
        assert math.isclose(mixed_set.lnprior_unconstrained(y), expected)

    def test_declared_bijection_overrides_the_default(self) -> None:
        p = Parameter("t", st.uniform(100.0, 9900.0), bijection=Identity())
        assert p.unconstraining_bijection() == Identity()

    def test_unsupported_support_asks_for_an_explicit_declaration(self) -> None:
        upper_bounded = ParameterSet([Parameter("x", st.norm(0.0, 1.0))])[0].release(
            st.uniform(0.0, 1.0)
        )
        assert upper_bounded.unconstraining_bijection() == Logit(0.0, 1.0)
        # A genuinely (-inf, b] support has no built-in bijection.
        with pytest.raises(ParameterError, match="bounded above but not below"):
            default_bijection_for(_UpperBounded())

    def test_logit_bounds_must_be_ordered(self) -> None:
        with pytest.raises(ParameterError, match="upper > lower"):
            Logit(1.0, 0.0)

    @pytest.mark.parametrize(
        ("prior", "expected"),
        [
            # Location-scale, unbounded: safe whatever the hyperparameters are.
            (HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}), Identity()),
            (HierarchicalPrior("cauchy", {"loc": "mu"}), Identity()),
            # A lower bound at a *constant* loc is unmoved by scale > 0.
            (
                HierarchicalPrior("halfnorm", {"scale": "sigma"}, kwds={"loc": 0.0}),
                Log(lower=0.0),
            ),
            (
                HierarchicalPrior("halfnorm", {"scale": "sigma"}, kwds={"loc": 3.0}),
                Log(lower=3.0),
            ),
        ],
    )
    def test_hierarchical_defaults_are_inferred_only_when_provably_safe(
        self, prior: HierarchicalPrior, expected: object
    ) -> None:
        assert default_bijection_for(prior) == expected

    @pytest.mark.parametrize(
        ("prior", "match"),
        [
            # The lower bound moves with mu, so no static bijection is correct.
            (HierarchicalPrior("halfnorm", {"loc": "mu"}), "bounded support whose position"),
            (
                HierarchicalPrior("uniform", {"loc": "a", "scale": "b"}),
                "bounded support whose position",
            ),
            # Shape arguments mean loc/scale do not determine the support.
            (HierarchicalPrior("gamma", {"scale": "s"}, args=(2.0,)), "takes shape arguments"),
        ],
    )
    def test_hierarchical_defaults_refuse_to_guess(
        self, prior: HierarchicalPrior, match: str
    ) -> None:
        with pytest.raises(ParameterError, match=match):
            default_bijection_for(prior)

    def test_an_explicit_bijection_unblocks_a_constrained_hierarchical_prior(self) -> None:
        pset = ParameterSet(
            [
                Parameter("sigma", st.halfnorm(0.0, 1.0)),
                Parameter(
                    "tau",
                    HierarchicalPrior("halfnorm", {"loc": "mu"}),
                    shape=(2,),
                    bijection=Log(lower=0.0),
                ),
                Parameter("mu", st.norm(0.0, 1.0)),
            ]
        )
        assert pset["tau"].unconstraining_bijection() == Log(lower=0.0)


class _UpperBounded:
    def ppf(self, q):  # pragma: no cover - only support() is exercised
        return q

    def support(self):
        return (-np.inf, 3.0)


class TestDiscreteFamilies:
    """The 2026-09-03 ruling on lowering.md §12.1: refuse the bijection, door ajar.

    The refusal lives only in ``default_bijection_for``; declaration, prior
    sampling, constrained-space log-probabilities and ``prior_transform``
    (scipy's discrete families implement ``ppf``) are untouched, and
    discreteness is queryable from the canonical description.
    """

    def test_the_bijection_is_refused_with_a_typed_capability_error(self) -> None:
        from ampere.core import CapabilityError

        with pytest.raises(CapabilityError, match="discrete") as excinfo:
            default_bijection_for(st.poisson(3.0))
        # A capability refusal, not a malformed declaration: catchable as
        # NotImplementedError, deliberately NOT a ValueError/ContractError.
        assert isinstance(excinfo.value, NotImplementedError)
        assert not isinstance(excinfo.value, ValueError)

    def test_a_discrete_hierarchical_family_is_refused_the_same_way(self) -> None:
        from ampere.core import CapabilityError

        with pytest.raises(CapabilityError, match="discrete"):
            default_bijection_for(HierarchicalPrior("poisson", {"mu": "rate"}))

    def test_everything_but_the_bijection_still_works(self) -> None:
        from ampere.core import CapabilityError

        counts = Parameter("counts", st.poisson(3.0))
        assert describe_prior(counts.prior).discrete is True  # the query
        pset = ParameterSet([counts])
        drawn = pset.sample(np.random.default_rng(7))  # prior sampling
        assert float(drawn["counts"]) == int(drawn["counts"])
        assert np.isfinite(pset.lnprior([2.0]))  # constrained log-prob
        quantile = pset.prior_transform([0.7])  # nested sampling's route
        assert quantile == st.poisson(3.0).ppf(0.7)
        with pytest.raises(CapabilityError):  # and only this refuses
            counts.unconstraining_bijection()

    def test_a_continuous_family_still_infers_normally(self) -> None:
        assert default_bijection_for(st.expon(0.0, 1.0)) == Log(lower=0.0)


class TestExceptionHomes:
    """The W1.13 exception dispositions (ruled 2026-09-03, lowering.md §12.4)."""

    def test_lowering_error_is_an_ampere_error_but_not_a_contract_error(self) -> None:
        # lowering.md §3.4: a family torch does not implement is a capability
        # gap in the backend, not a malformed declaration by the user, so the
        # two must stay distinguishable to a caller catching by type.
        from ampere.core import AmpereError, ContractError, LoweringError

        err = LoweringError("truncnorm", backend="torch", parameter="temperature")
        assert isinstance(err, AmpereError)
        assert not isinstance(err, ContractError)
        assert (err.family, err.parameter, err.backend) == ("truncnorm", "temperature", "torch")
        assert "truncnorm" in str(err) and "torch" in str(err)
        assert "never substitutes an approximation" in str(err)

    def test_optional_dependency_error_stands_ratified_in_place(self) -> None:
        # parameters.md §14 Q5, closed at the freeze: the shape pinned in
        # ampere/core/exceptions.py is the contract.
        err = OptionalDependencyError("paramax", extra="jax", context="lowering a ParameterSet")
        assert isinstance(err, ImportError)
        assert (err.package, err.extra) == ("paramax", "jax")

    def test_results_error_lives_in_core_and_is_reexported(self) -> None:
        # results.md §15 R5, implemented 2026-09-03: one class, two import paths.
        from ampere.core.exceptions import ResultsError as from_core
        from ampere.results.exceptions import ResultsError as from_results

        assert from_core is from_results


class TestLoweringExtensionPoint:
    def test_as_paramax_explains_that_it_belongs_to_the_backend(self) -> None:
        pset = ParameterSet([Parameter("t", st.norm(0.0, 1.0))])
        with pytest.raises(OptionalDependencyError) as excinfo:
            pset.as_paramax()
        assert excinfo.value.package == "paramax"
        assert excinfo.value.extra == "jax"
        assert "ampere.backends.jax" in str(excinfo.value)

    def test_core_imports_no_optional_dependency(self) -> None:
        """ampere.core must never pull in torch, jax or paramax (architecture.md 4).

        Run in a fresh interpreter, so the result cannot be polluted by another
        test in the same session having imported a heavy dependency.
        """
        import subprocess
        import sys

        script = (
            "import sys, ampere.core; "
            "bad = [m for m in ('torch', 'jax', 'paramax', 'numpyro', 'pyro') "
            "if m in sys.modules]; "
            "print(','.join(bad))"
        )
        result = subprocess.run(
            [sys.executable, "-c", script], capture_output=True, text=True, check=True
        )
        assert result.stdout.strip() == "", (
            f"ampere.core imported optional dependencies: {result.stdout.strip()}"
        )


class TestLosslessNesting:
    """merge accepting a ParameterMapping component (ruled 2026-09-02).

    Bindings compose to the leaves; routing stays one level deep.
    """

    @staticmethod
    def _instrument_like() -> ParameterMapping:
        return ParameterSet.merge(
            {
                "a": ParameterSet([Parameter("scale", st.lognorm(0.2), shared_as="gain")]),
                "b": ParameterSet([Parameter("scale", st.lognorm(0.2), shared_as="gain")]),
            }
        )

    def test_an_inner_collapse_surfaces_in_the_outer_bindings(self) -> None:
        inner = self._instrument_like()
        assert inner.tied_names == ("gain",)
        outer = ParameterSet.merge(
            {"instrument": inner, "model": ParameterSet([Parameter("t", st.norm(0, 1))])}
        )
        assert outer.merged.names == ("instrument.gain", "model.t")
        assert outer.tied_names == ("instrument.gain",)
        assert {(b.component, b.local_name) for b in outer.sites_of("instrument.gain")} == {
            ("instrument", "a.scale"),
            ("instrument", "b.scale"),
        }

    def test_routing_stays_one_level_deep(self) -> None:
        inner = self._instrument_like()
        outer = ParameterSet.merge(
            {"instrument": inner, "model": ParameterSet([Parameter("t", st.norm(0, 1))])}
        )
        routed = outer.distribute({"instrument.gain": 2.0, "model.t": 0.1})
        assert routed["instrument"] == {"gain": 2.0}
        assert inner.distribute(routed["instrument"]) == {
            "a": {"scale": 2.0},
            "b": {"scale": 2.0},
        }

    def test_plain_set_components_are_unchanged(self) -> None:
        mapping = ParameterSet.merge(
            {
                "x": ParameterSet([Parameter("p", st.norm(0, 1))]),
                "y": ParameterSet([Parameter("q", st.norm(0, 1))]),
            }
        )
        assert mapping.bindings == mapping.routing
        assert mapping.tied_names == ()

    def test_a_cross_level_tie_composes_to_both_leaf_sets(self) -> None:
        def leaf() -> ParameterSet:
            return ParameterSet([Parameter("scale", st.norm(1.0, 0.1))])

        def dataset_like():
            return ParameterSet.merge({"instr": leaf()})

        outer = ParameterSet.merge(
            {"obj0": dataset_like(), "obj1": dataset_like()},
            ties=[Tie("cal", ("obj0.instr.scale", "obj1.instr.scale"))],
        )
        assert outer.tied_names == ("cal",)
        assert {(b.component, b.local_name) for b in outer.sites_of("cal")} == {
            ("obj0", "instr.scale"),
            ("obj1", "instr.scale"),
        }

    def test_global_name_for_accepts_both_forms(self) -> None:
        inner = self._instrument_like()
        outer = ParameterSet.merge({"instrument": inner})
        assert outer.global_name_for("instrument", "gain") == "instrument.gain"
        assert outer.global_name_for("instrument", "a.scale") == "instrument.gain"


class TestPlateBindings:
    """Binding.index element routing (ruled 2026-09-02, gap H-2)."""

    @staticmethod
    def _survey() -> ParameterSet:
        plate = Plate(
            "objects",
            size=3,
            hyperparameters=[
                Parameter("mu", st.norm(0.0, 5.0)),
                Parameter("sigma", st.halfnorm(0.0, 2.0)),
            ],
            members=[
                Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))
            ],
        )
        return ParameterSet([], plates=[plate])

    def _merged(self):
        return ParameterSet.merge(
            {
                "population": self._survey(),
                "obj0": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
                "obj1": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
            },
            plate_bindings=[
                PlateBinding("population.objects.theta", "obj0", "theta", 0),
                PlateBinding("population.objects.theta", "obj1", "theta", 1),
            ],
        )

    def test_each_component_receives_its_own_element(self) -> None:
        mapping = self._merged()
        values = {
            "population.objects.mu": 0.0,
            "population.objects.sigma": 1.0,
            "population.objects.theta": np.array([10.0, 20.0, 30.0]),
            "obj0.cal": 1.1,
            "obj1.cal": 0.9,
        }
        routed = mapping.distribute(values)
        assert routed["obj0"] == {"cal": 1.1, "theta": 10.0}
        assert routed["obj1"] == {"cal": 0.9, "theta": 20.0}
        # The merged set is unchanged: one array-valued parameter, no extras.
        assert mapping.merged["population.objects.theta"].shape == (3,)

    def test_element_bindings_are_addressing_not_tying(self) -> None:
        mapping = self._merged()
        assert "population.objects.theta" not in mapping.tied_names
        by_index = {
            b.index: b.component
            for b in mapping.sites_of("population.objects.theta")
            if b.index is not None
        }
        assert by_index == {0: "obj0", 1: "obj1"}

    def test_bad_declarations_are_refused_loudly(self) -> None:
        base = {
            "population": self._survey(),
            "obj0": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
        }
        with pytest.raises(ParameterError, match="does not exist"):
            ParameterSet.merge(
                dict(base), plate_bindings=[PlateBinding("population.ghost", "obj0", "t", 0)]
            )
        with pytest.raises(ParameterError, match="scalar"):
            ParameterSet.merge(
                dict(base),
                plate_bindings=[PlateBinding("population.objects.mu", "obj0", "t", 0)],
            )
        with pytest.raises(ParameterError, match="out of range"):
            ParameterSet.merge(
                dict(base),
                plate_bindings=[PlateBinding("population.objects.theta", "obj0", "t", 7)],
            )
        with pytest.raises(ParameterError, match="does not have"):
            ParameterSet.merge(
                dict(base),
                plate_bindings=[PlateBinding("population.objects.theta", "ghost", "t", 0)],
            )
        with pytest.raises(ParameterError, match="already receives"):
            ParameterSet.merge(
                dict(base),
                plate_bindings=[PlateBinding("population.objects.theta", "obj0", "cal", 0)],
            )
        with pytest.raises(ParameterError, match="int or a non-empty tuple"):
            PlateBinding("population.objects.theta", "obj0", "t", 1.5)  # type: ignore[arg-type]


class TestDescribeHook:
    """The opt-in ``describe()`` hook (``results.md`` §13.13, landed W2.1).

    Provenance reaches parameters, buffers and class identity, but not a plain
    Python attribute that changes what a model computes. Hashing an arbitrary
    ``__dict__`` is not a safe general answer, so this is the declared
    extension point instead.
    """

    def test_the_default_declares_nothing(self) -> None:
        class Plain(Parameterised):
            pass

        assert Plain().describe() is None

    def test_a_subclass_may_declare_its_configuration(self) -> None:
        class Redden(Parameterised):
            def __init__(self, law: str) -> None:
                self.law = law

            def describe(self) -> dict[str, str]:
                return {"law": self.law}

        assert Redden("ccm89").describe() == {"law": "ccm89"}
        assert Redden("f99").describe() != Redden("ccm89").describe()

    def test_describe_may_not_be_used_as_a_parameter_name(self) -> None:
        """It is a class attribute now, so the shadowing guard covers it.

        The same already holds for ``parameters``, ``buffers`` and ``context``;
        this row records that ``describe`` joined them, because it is a (small)
        narrowing of what a model may call its parameters.
        """

        class Shadowed(Parameterised):
            pass

        with pytest.raises(ParameterError, match="shadows an attribute"):
            Shadowed().register_parameter(Parameter("describe", st.norm(0.0, 1.0)))

    def test_describe_may_not_be_used_as_a_buffer_name(self) -> None:
        class Shadowed(Parameterised):
            pass

        with pytest.raises(ParameterError, match="shadows an attribute"):
            Shadowed().register_buffer("describe", np.arange(3.0))
