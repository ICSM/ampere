"""Conformance rows for the parameter and prior contract.

The inventory is ``parameters.md`` §13's hand-down, verbatim: "pack/unpack
identity, spec round trip, ``prior_transform`` against analytic quantiles,
``lnprior`` against summed scipy ``logpdf``, ``constrain``∘``unconstrain``
identity, the Jacobian identity in §6, and merge dimension counting."

Every oracle here is ``scipy`` or a closed form transcribed from the contract.
Nothing compares ampere against ampere.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Identity, Parameter, ParameterSet, Tie

from .oracles import (
    analytic_constrain,
    analytic_lnprior,
    analytic_prior_transform,
    free_slices,
    summed_log_abs_det,
)
from .protocol import ConformanceBackend, ParameterSpace, Tolerances

# ---------------------------------------------------------------------------
# Declarations under test
# ---------------------------------------------------------------------------


def mixed_declaration() -> ParameterSet:
    """One declaration exercising all three bijections, plus shape and fixing.

    ``uniform`` is bounded both ways (``Logit``), ``norm`` neither way
    (``Identity``) and ``halfnorm`` below only (``Log``) — which is the whole
    of ``default_bijection_for``'s table. The array-valued and fixed entries
    are there because the flat layout is what ``free_labels`` and ``pack`` are
    a contract about.
    """
    return ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
            Parameter("log_tau", st.norm(0.0, 1.0)),
            Parameter("width", st.halfnorm(0.0, 2.0)),
            Parameter("offset", st.norm(0.0, 0.05), shape=(3,), unit=u.mag),
            Parameter("distance", value=1.5 * u.kpc, fixed=True, unit=u.kpc),
        ]
    )


def tied_component() -> ParameterSet:
    """A component declaring one private and one shared parameter."""
    return ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
            Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc, shared_as="distance"),
        ]
    )


@pytest.fixture
def space(backend: ConformanceBackend) -> ParameterSpace:
    """*backend*'s realisation of :func:`mixed_declaration`."""
    return backend.parameter_space(mixed_declaration())


CUBE = np.array([0.3, 0.62, 0.17, 0.44, 0.81, 0.05])
UNCONSTRAINED = np.array([0.4, -1.2, 0.7, 0.05, -0.3, 1.1])


# ---------------------------------------------------------------------------
# Rows
# ---------------------------------------------------------------------------


class TestFlatLayout:
    """``free_size``, ``free_labels`` and the pack/unpack round trip."""

    def test_the_free_dimension_counts_array_elements_once_each(
        self, space: ParameterSpace
    ) -> None:
        assert space.free_size == 6
        assert len(space.free_labels()) == space.free_size

    def test_free_labels_name_every_flat_entry(self, space: ParameterSpace) -> None:
        assert space.free_labels() == (
            "temperature",
            "log_tau",
            "width",
            "offset[0]",
            "offset[1]",
            "offset[2]",
        )

    def test_pack_and_unpack_are_inverse(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        theta = analytic_prior_transform(space.declaration, CUBE)
        assert space.pack(space.unpack(theta)) == pytest.approx(theta, abs=tolerances.exact)

    def test_unpack_restores_every_parameter_including_the_fixed_one(
        self, space: ParameterSpace
    ) -> None:
        values = space.unpack(analytic_prior_transform(space.declaration, CUBE))
        assert set(values) == {"temperature", "log_tau", "width", "offset", "distance"}
        assert values["distance"] == pytest.approx(1.5)
        assert np.shape(values["offset"]) == (3,)


class TestSpecRoundTrip:
    """The serialisable declaration is lossless."""

    def test_the_declaration_survives_to_spec_and_back(self, space: ParameterSpace) -> None:
        rebuilt = ParameterSet.from_spec(space.declaration.to_spec())
        assert rebuilt.to_spec() == space.declaration.to_spec()
        assert rebuilt.names == space.declaration.names
        assert rebuilt.free_size == space.free_size

    def test_the_round_tripped_declaration_scores_identically(
        self, backend: ConformanceBackend, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        rebuilt = backend.parameter_space(ParameterSet.from_spec(space.declaration.to_spec()))
        theta = analytic_prior_transform(space.declaration, CUBE)
        assert rebuilt.lnprior(theta) == pytest.approx(space.lnprior(theta), abs=tolerances.exact)


class TestPriorsAgainstScipy:
    """``prior_transform`` and ``lnprior`` against the distributions themselves."""

    def test_prior_transform_is_the_analytic_quantile(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        expected = analytic_prior_transform(space.declaration, CUBE)
        assert space.prior_transform(CUBE) == pytest.approx(expected, abs=tolerances.analytic)

    def test_the_centre_of_the_cube_is_the_prior_median(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        centre = np.full(space.free_size, 0.5)
        expected = [
            parameter.prior.median()
            for parameter, where in free_slices(space.declaration)
            for _ in range(where.stop - where.start)
        ]
        assert space.prior_transform(centre) == pytest.approx(expected, abs=tolerances.analytic)

    def test_lnprior_is_the_summed_scipy_logpdf(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        theta = analytic_prior_transform(space.declaration, CUBE)
        expected = analytic_lnprior(space.declaration, theta)
        assert space.lnprior(theta) == pytest.approx(expected, abs=tolerances.analytic)

    def test_a_point_outside_the_support_has_no_prior_mass(self, space: ParameterSpace) -> None:
        theta = analytic_prior_transform(space.declaration, CUBE)
        theta[0] = -1.0  # temperature, whose support starts at 100 K
        assert space.lnprior(theta) == -np.inf


class TestUnconstrainedSpace:
    """``constrain``∘``unconstrain`` and the §6 Jacobian identity."""

    def test_constrain_is_the_analytic_bijection(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        expected = np.empty_like(UNCONSTRAINED)
        for parameter, where in free_slices(space.declaration):
            expected[where] = analytic_constrain(parameter, UNCONSTRAINED[where])
        assert space.constrain(UNCONSTRAINED) == pytest.approx(expected, abs=tolerances.analytic)

    def test_constrain_and_unconstrain_are_inverse(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        theta = space.constrain(UNCONSTRAINED)
        assert space.unconstrain(theta) == pytest.approx(
            UNCONSTRAINED, abs=tolerances.linear_algebra
        )
        assert space.constrain(space.unconstrain(theta)) == pytest.approx(
            theta, abs=tolerances.linear_algebra
        )

    def test_the_jacobian_identity_holds_term_by_term(
        self, space: ParameterSpace, tolerances: Tolerances
    ) -> None:
        """``lnprior_unconstrained(y) == lnprior(constrain(y)) + Σ log|dc/dy|``.

        ``parameters.md`` §6, and the reference semantics ``lowering.md``
        makes the native gradient paths agree with. The Jacobian sum is
        computed here from the priors' supports, never from ampere's own
        bijections.
        """
        jacobian = summed_log_abs_det(space.declaration, UNCONSTRAINED)
        theta = space.constrain(UNCONSTRAINED)
        expected = analytic_lnprior(space.declaration, theta) + jacobian
        assert space.lnprior_unconstrained(UNCONSTRAINED) == pytest.approx(
            expected, abs=tolerances.analytic
        )

    def test_the_default_bijection_makes_the_support_inescapable(
        self, space: ParameterSpace
    ) -> None:
        """No unconstrained vector maps outside the priors' support.

        The corollary a gradient-based engine relies on: with the default
        bijections, ``lnprior(constrain(y))`` is finite for every finite ``y``,
        so a NUTS proposal cannot land on zero prior mass. Checked at a point
        far enough out that a ``Logit`` saturates in float64.
        """
        far = np.full(space.free_size, -40.0)
        theta = space.constrain(far)
        for parameter, where in free_slices(space.declaration):
            low, high = parameter.prior.support()
            assert np.all(theta[where] >= low)
            assert np.all(theta[where] <= high)
        assert np.isfinite(space.lnprior_unconstrained(far))

    def test_a_declared_identity_over_a_bounded_prior_still_reports_no_mass(
        self, backend: ConformanceBackend
    ) -> None:
        """The ``-inf`` short circuit, on the one declaration that can reach it.

        A user may override the default bijection — ``Identity`` over a
        bounded prior is legitimate and is what a fixed-support reparametrisation
        looks like — and then an unconstrained point *can* leave the support.
        ``lnprior_unconstrained`` must answer ``-inf`` rather than adding a
        Jacobian to a ``-inf``.
        """
        escaped = backend.parameter_space(
            ParameterSet([Parameter("bounded", st.uniform(0.0, 1.0), bijection=Identity())])
        )
        assert escaped.lnprior_unconstrained(np.array([-40.0])) == -np.inf
        assert escaped.lnprior_unconstrained(np.array([0.5])) == pytest.approx(0.0)


class TestMergeDimensionCounting:
    """A tie collapses two sites into one dimension, and says which two."""

    def test_a_shared_parameter_costs_one_dimension_not_two(
        self, backend: ConformanceBackend
    ) -> None:
        mapping = ParameterSet.merge({"sed": tied_component(), "spectrum": tied_component()})
        space = backend.parameter_space(mapping.merged)
        assert space.free_size == 3
        assert mapping.tied_names == ("distance",)
        assert space.free_labels() == ("sed.temperature", "distance", "spectrum.temperature")

    def test_the_tie_has_one_binding_per_site(self) -> None:
        mapping = ParameterSet.merge({"sed": tied_component(), "spectrum": tied_component()})
        sites = mapping.sites_of("distance")
        assert len(sites) == 2
        assert {(binding.component, binding.local_name) for binding in sites} == {
            ("sed", "distance"),
            ("spectrum", "distance"),
        }

    def test_a_composition_time_tie_counts_the_same_way(self, backend: ConformanceBackend) -> None:
        untied = ParameterSet(
            [
                Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K),
                Parameter("distance", st.norm(1.5, 0.1), unit=u.kpc),
            ]
        )
        mapping = ParameterSet.merge(
            {"sed": untied, "spectrum": untied},
            ties=[Tie("distance", ("sed.distance", "spectrum.distance"))],
        )
        assert backend.parameter_space(mapping.merged).free_size == 3
        assert mapping.tied_names == ("distance",)

    def test_one_value_reaches_every_tied_site(self) -> None:
        mapping = ParameterSet.merge({"sed": tied_component(), "spectrum": tied_component()})
        routed = mapping.distribute(
            {"sed.temperature": 5000.0, "spectrum.temperature": 6000.0, "distance": 1.7}
        )
        assert routed["sed"]["distance"] == pytest.approx(1.7)
        assert routed["spectrum"]["distance"] == pytest.approx(1.7)
