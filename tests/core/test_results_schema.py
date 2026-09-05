"""Unit tests for the W1.4 ModelResult schema contract.

Organised around the acceptance criteria in ``WORK_ITEMS.md`` W1.4 and the
dispatch brief: the low-res-SED + CO-windows example from
``DEVELOPMENT_PLAN.md`` §4.2 expressed and validated (``TestLowResSedPlusCoWindows``,
the headline criterion), named channels with kind checking and loud mismatches,
default-channel sugar for single-output models, coordinate-indexed containers
with uncertainties and first-class masks, units converted once, fidelity tags,
and complex-valued data — followed by the supporting behaviour each rests on.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    COORDINATE_RTOL,
    DEFAULT_CHANNEL,
    Axis,
    AxisSpec,
    Cube,
    FunctionSamples,
    Image,
    Layout,
    ModelResult,
    Order,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    VisibilitySet,
)
from ampere.core.exceptions import ChannelError, ContractError, ParameterError, SchemaError


@pytest.fixture
def sed() -> Spectrum:
    """A broad, coarsely sampled, logarithmically spaced continuum SED."""
    return Spectrum(
        np.geomspace(1.0, 2000.0, 60) * u.um,
        np.geomspace(50.0, 0.2, 60) * u.Jy,
        fidelity="continuum",
    )


@pytest.fixture
def co_windows() -> Spectrum:
    """Two narrow, finely sampled windows around CO rotational lines."""
    window_32 = np.linspace(866.86, 867.06, 21)
    window_21 = np.linspace(1300.30, 1300.50, 21)
    return Spectrum(
        np.concatenate([window_32, window_21]) * u.um,
        np.concatenate([np.linspace(1.0, 1.4, 21), np.linspace(0.8, 1.1, 21)]) * u.Jy,
        fidelity="co_lte",
    )


# ---------------------------------------------------------------------------
# The headline acceptance criterion
# ---------------------------------------------------------------------------


class TestLowResSedPlusCoWindows:
    """The worked example ``DEVELOPMENT_PLAN.md`` §4.2 uses to justify channels.

    One model, two outputs of the *same kind* at different resolutions over
    different (and non-contiguous) spectral ranges. A kind-keyed schema could
    not express this at all; that is the point.
    """

    def test_both_channels_coexist_under_distinct_names(self, sed, co_windows) -> None:
        result = ModelResult({"sed_lowres": sed, "co_windows": co_windows})
        assert sorted(result) == ["co_windows", "sed_lowres"]
        assert len(result) == 2

    def test_two_channels_of_the_same_kind_stay_distinguishable(self, sed, co_windows) -> None:
        result = ModelResult({"sed_lowres": sed, "co_windows": co_windows})
        assert {name: kind.__name__ for name, kind in result.kinds().items()} == {
            "sed_lowres": "Spectrum",
            "co_windows": "Spectrum",
        }
        assert sorted(result.of_kind(Spectrum)) == ["co_windows", "sed_lowres"]
        # Kind alone cannot tell them apart, so the *name* must carry it.
        assert result.require("sed_lowres", Spectrum) is sed
        assert result.require("co_windows", Spectrum) is co_windows

    def test_the_channels_have_genuinely_different_sampling(self, sed, co_windows) -> None:
        assert sed.n_samples == 60
        assert co_windows.n_samples == 42
        # The SED is log-regular over three decades; the CO channel is neither
        # regular nor log-regular, because of the gap between the two windows.
        assert sed.spectral_axis.log_regular
        assert not co_windows.spectral_axis.regular
        assert not co_windows.spectral_axis.log_regular

    def test_the_co_channel_is_ordered_but_gapped(self, co_windows) -> None:
        """Strictly increasing, so a quasiseparable solver can consume it, while
        being nothing like an evenly spaced grid."""
        diffs = np.diff(co_windows.spectral_axis.values)
        assert np.all(diffs > 0.0)
        assert diffs.max() > 400.0  # the inter-window gap, in microns
        assert diffs.min() < 0.02  # the in-window sampling

    def test_each_channel_carries_its_own_fidelity_tag(self, sed, co_windows) -> None:
        result = ModelResult({"sed_lowres": sed, "co_windows": co_windows})
        assert dict(result.fidelities()) == {
            "sed_lowres": "continuum",
            "co_windows": "co_lte",
        }

    def test_binding_a_missing_channel_fails_loudly_naming_what_exists(
        self, sed, co_windows
    ) -> None:
        result = ModelResult({"sed_lowres": sed, "co_windows": co_windows})
        with pytest.raises(ChannelError) as excinfo:
            result.require("co_lines")
        message = str(excinfo.value)
        assert "no channel named 'co_lines'" in message
        assert "'sed_lowres' (Spectrum)" in message
        assert "'co_windows' (Spectrum)" in message

    def test_the_hot_loop_refills_each_channel_without_revalidating(self, sed, co_windows) -> None:
        result = ModelResult({"sed_lowres": sed, "co_windows": co_windows})
        updated = result.with_channels(
            sed_lowres=sed.with_values(np.geomspace(55.0, 0.18, 60)),
            co_windows=co_windows.with_values(np.ones(42)),
        )
        assert updated["sed_lowres"].axes is sed.axes
        assert updated["co_windows"].axes is co_windows.axes
        assert updated["sed_lowres"].unit == u.Jy
        # The original result is untouched.
        assert float(result["sed_lowres"].values[0]) == pytest.approx(50.0)


# ---------------------------------------------------------------------------
# ModelResult
# ---------------------------------------------------------------------------


class TestNamedChannels:
    def test_bare_container_gets_the_default_channel(self, sed) -> None:
        result = ModelResult(sed)
        assert list(result) == [DEFAULT_CHANNEL]
        assert result.is_single
        assert result.single() is sed

    def test_mapping_protocol(self, sed, co_windows) -> None:
        result = ModelResult({"a": sed, "b": co_windows})
        assert "a" in result and "z" not in result
        assert set(result.keys()) == {"a", "b"}
        assert dict(result.items())["b"] is co_windows
        assert result.get("z") is None
        assert list(result.values()) == [sed, co_windows]

    def test_single_refuses_to_guess_on_a_multi_channel_result(self, sed, co_windows) -> None:
        result = ModelResult({"a": sed, "b": co_windows})
        with pytest.raises(ChannelError, match="only meaningful for a one-channel result"):
            result.single()

    def test_default_only_result_explains_why_a_name_is_missing(self, sed) -> None:
        with pytest.raises(ChannelError, match="returned a bare container"):
            ModelResult(sed).require("sed")

    def test_near_miss_channel_name_is_suggested(self, sed) -> None:
        result = ModelResult({"SED": sed})
        with pytest.raises(ChannelError, match="Did you mean 'SED'"):
            result.require("sed")

    def test_channel_names_must_be_identifiers(self, sed) -> None:
        with pytest.raises(SchemaError, match="must be a valid Python identifier"):
            ModelResult({"sed lowres": sed})

    def test_channel_names_may_be_dot_qualified(self, sed, co_windows) -> None:
        # Ruled 2026-09-02 (results_schema.md §17): a model namespaces grouped
        # output with dot-qualified names — per-object channels or one object's
        # several CO lines — mirroring ParameterSet.merge's qualification.
        result = ModelResult({"co.j3_2": co_windows, "obj1.sed": sed})
        assert result.require("obj1.sed", Spectrum) is sed
        assert "co.j3_2" in result

    def test_a_malformed_dotted_channel_name_is_refused(self, sed) -> None:
        for bad in ("obj1..sed", ".sed", "sed.", "obj 1.sed"):
            with pytest.raises(SchemaError, match=r"'\.'-separated sequence"):
                ModelResult({bad: sed})

    def test_empty_result_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="at least one channel"):
            ModelResult({})

    def test_channel_must_hold_a_container(self) -> None:
        with pytest.raises(SchemaError, match="which is not a container"):
            ModelResult({"sed": np.ones(3)})

    def test_non_mapping_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="built from a single container or a mapping"):
            ModelResult([1, 2, 3])

    def test_immutable_updates(self, sed, co_windows) -> None:
        result = ModelResult({"sed": sed}, meta={"model": "toy"})
        added = result.with_channels(lines=co_windows)
        assert sorted(added) == ["lines", "sed"]
        assert list(result) == ["sed"], "the original must be untouched"
        assert added.meta["model"] == "toy", "metadata survives an update"
        assert list(added.without_channels("lines")) == ["sed"]

    def test_removing_an_absent_channel_fails_loudly(self, sed) -> None:
        with pytest.raises(ChannelError, match="no channel named 'nope'"):
            ModelResult({"sed": sed}).without_channels("nope")

    def test_equality_compares_channels_and_metadata(self, sed) -> None:
        assert ModelResult({"sed": sed}) == ModelResult({"sed": sed})
        assert ModelResult({"sed": sed}) != ModelResult({"other": sed})
        assert ModelResult({"sed": sed}) != ModelResult({"sed": sed}, meta={"x": 1})

    def test_results_are_unhashable_because_containers_are(self, sed) -> None:
        with pytest.raises(TypeError):
            hash(ModelResult(sed))


class TestCarriedParameters:
    """`ModelResult.parameters` — the (θ, result) pairing ruled 2026-09-01."""

    def test_a_result_may_carry_the_theta_that_produced_it(self, sed) -> None:
        tagged = ModelResult(sed, parameters={"temperature": 300.0})
        assert tagged.parameters["temperature"] == 300.0
        assert ModelResult(sed).parameters is None, "absent, not empty"

    def test_the_record_survives_channel_updates(self, sed) -> None:
        tagged = ModelResult({"sed": sed}, parameters={"temperature": 300.0})
        assert tagged.with_channels(line=sed).parameters["temperature"] == 300.0
        assert tagged.with_channels(line=sed).without_channels("line").parameters == {
            "temperature": 300.0
        }

    def test_the_record_is_immutable(self, sed) -> None:
        tagged = ModelResult(sed, parameters={"temperature": 300.0})
        with pytest.raises(TypeError):
            tagged.parameters["temperature"] = 400.0  # type: ignore[index]

    def test_equality_is_array_aware(self, sed) -> None:
        record = {"offset": np.array([1.0, 2.0])}
        assert ModelResult(sed, parameters=record) == ModelResult(sed, parameters=record)
        assert ModelResult(sed, parameters=record) != ModelResult(
            sed, parameters={"offset": np.array([1.0, 3.0])}
        )
        assert ModelResult(sed, parameters=record) != ModelResult(sed)

    def test_a_non_mapping_record_is_refused(self, sed) -> None:
        with pytest.raises(SchemaError, match="parameters record"):
            ModelResult(sed, parameters=[300.0])  # type: ignore[arg-type]


class TestKindChecking:
    def test_matching_kind_returns_the_container(self, sed) -> None:
        assert ModelResult({"sed": sed}).require("sed", Spectrum) is sed

    def test_mismatched_kind_raises_and_names_both_kinds(self, sed) -> None:
        with pytest.raises(ChannelError) as excinfo:
            ModelResult({"sed": sed}).require("sed", VisibilitySet)
        message = str(excinfo.value)
        assert "holds a Spectrum" in message
        assert "a VisibilitySet was required" in message

    def test_mismatch_names_channels_of_the_wanted_kind_when_there_are_some(
        self, sed, co_windows
    ) -> None:
        vis = VisibilitySet([1.0, 2.0], [3.0, 4.0], [1 + 0j, 0 + 1j])
        result = ModelResult({"sed": sed, "uv": vis})
        with pytest.raises(ChannelError, match=r"\['uv'\]"):
            result.require("sed", VisibilitySet)

    def test_subclasses_satisfy_a_base_kind_requirement(self) -> None:
        class NarrowBand(Spectrum):
            pass

        narrow = NarrowBand([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        assert ModelResult({"n": narrow}).require("n", Spectrum) is narrow

    def test_channel_error_is_catchable_as_the_contract_family(self, sed) -> None:
        """A consumer may catch the narrow type, the contract family, or the
        builtin it would naturally reach for."""
        for exception_type in (ChannelError, SchemaError, ContractError, ValueError, KeyError):
            with pytest.raises(exception_type):
                ModelResult({"sed": sed}).require("absent")

    def test_channel_error_keeps_its_message_despite_being_a_key_error(self, sed) -> None:
        """``KeyError`` formats itself as ``repr(args[0])``; this one must not."""
        with pytest.raises(ChannelError) as excinfo:
            ModelResult({"sed": sed}).require("absent")
        assert str(excinfo.value).startswith("no channel named 'absent'")

    def test_mapping_mixins_answer_rather_than_raise(self, sed) -> None:
        """``Mapping.get`` and ``__contains__`` are defined in terms of catching
        ``KeyError``, so ``ChannelError`` must be one for these to work."""
        result = ModelResult({"sed": sed})
        assert result.get("absent") is None
        assert result.get("absent", "fallback") == "fallback"
        assert "absent" not in result
        assert "sed" in result

    def test_schema_errors_are_not_parameter_errors(self, sed) -> None:
        """The two contracts have separate namespaces and separate failure modes."""
        assert not issubclass(SchemaError, ParameterError)
        with pytest.raises(SchemaError):
            ModelResult({"1bad": sed})


class TestChannelAndParameterNamespaces:
    """``parameters.md`` §13's obligation: the two namespaces may collide."""

    def test_a_channel_may_share_a_name_with_a_parameter(self) -> None:
        import scipy.stats as st

        from ampere.core import Parameter, ParameterSet

        parameters = ParameterSet([Parameter("temperature", st.uniform(100.0, 9900.0), unit=u.K)])
        temperature_map = Image(
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.full((3, 3), 2500.0) * u.K,
        )
        result = ModelResult({"temperature": temperature_map})
        # Both exist, neither disambiguates the other, and nothing raises.
        assert "temperature" in parameters.names
        assert "temperature" in result
        assert result["temperature"] is temperature_map


# ---------------------------------------------------------------------------
# Containers as coordinate-indexed function samples
# ---------------------------------------------------------------------------


class TestCoordinatesAndValues:
    def test_coordinates_are_explicit_and_aligned(self) -> None:
        spectrum = Spectrum([1.0, 2.0, 4.0] * u.um, [3.0, 2.5, 1.0] * u.Jy)
        assert spectrum.spectral_axis.values.tolist() == [1.0, 2.0, 4.0]
        assert spectrum.flux.tolist() == [3.0, 2.5, 1.0]
        assert spectrum.shape == (3,)
        assert spectrum.n_samples == 3

    def test_length_mismatch_is_caught_at_construction(self) -> None:
        with pytest.raises(SchemaError, match="aligned index-by-index"):
            Spectrum([1.0, 2.0, 4.0] * u.um, [3.0, 2.5])

    def test_uncertainty_shape_is_checked(self) -> None:
        with pytest.raises(SchemaError, match="uncertainties have shape"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], uncertainty=[0.1, 0.2, 0.3])

    def test_negative_uncertainty_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="negative uncertainties"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], uncertainty=[0.1, -0.2])

    def test_containers_and_their_arrays_are_immutable(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])
        with pytest.raises(ValueError, match="read-only"):
            spectrum.values[0] = 99.0
        with pytest.raises(ValueError, match="read-only"):
            spectrum.spectral_axis.values[0] = 99.0
        with pytest.raises(AttributeError):
            spectrum.unit = u.Jy

    def test_wrong_axis_names_are_refused(self) -> None:
        """The base constructor takes a coordinate mapping keyed by axis name."""
        with pytest.raises(SchemaError) as excinfo:
            FunctionSamples.__init__(
                TimeSeries.__new__(TimeSeries), {"wavelength": [1.0, 2.0] * u.um}, [1.0, 1.0]
            )
        message = str(excinfo.value)
        assert "is defined by the coordinate axes ['time']" in message
        assert "missing ['time']" in message
        assert "unexpected ['wavelength']" in message

    def test_points_layout_requires_equal_length_axes(self) -> None:
        with pytest.raises(SchemaError, match="must have the same length"):
            VisibilitySet([1.0, 2.0, 3.0], [1.0, 2.0], [1 + 0j, 0 + 1j])

    def test_mesh_coordinates_are_refused(self) -> None:
        with pytest.raises(SchemaError, match="must be one-dimensional"):
            Image(
                np.ones((2, 2)) * u.arcsec,
                np.linspace(-1.0, 1.0, 2) * u.arcsec,
                np.ones((2, 2)),
            )

    def test_extra_coords_are_aligned_and_may_be_non_numeric(self) -> None:
        spectrum = Spectrum(
            [1.0, 2.0] * u.um,
            [1.0, 1.0],
            extra_coords={"order": np.array(["blue", "red"])},
        )
        assert spectrum.extra_coords["order"].tolist() == ["blue", "red"]

    def test_extra_coords_shape_is_checked(self) -> None:
        with pytest.raises(SchemaError, match="extra coordinate 'order'"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], extra_coords={"order": np.array([1, 2, 3])})

    def test_extra_coords_refuse_a_quantity_rather_than_stripping_it(self) -> None:
        """v1.4 extra coordinates are unitless labels; dropping a supplied unit
        silently is exactly the quiet loss this contract exists to prevent."""
        with pytest.raises(SchemaError, match="unitless per-sample labels"):
            Spectrum(
                [1.0, 2.0] * u.um,
                [1.0, 1.0],
                extra_coords={"frequency": [230.0, 345.0] * u.GHz},
            )

    def test_extra_coords_may_not_shadow_an_axis(self) -> None:
        with pytest.raises(SchemaError, match="collides with one of its declared axes"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], extra_coords={"spectral_axis": [1.0, 2.0]})

    def test_equality_compares_coordinates_values_and_mask(self) -> None:
        a = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        assert a == Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        assert a != Spectrum([1.0, 3.0] * u.um, [1.0, 1.0] * u.Jy)
        assert a != Spectrum([1.0, 2.0] * u.um, [1.0, 2.0] * u.Jy)
        assert a != Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy, mask=np.array([True, False]))


class TestOrdering:
    """``architecture.md`` §7: validated or documented, never assumed."""

    @pytest.mark.parametrize(
        ("build", "match"),
        [
            (lambda: Spectrum([2.0, 1.0] * u.um, [1.0, 1.0]), "strictly increasing"),
            (lambda: TimeSeries([2.0, 1.0] * u.day, [1.0, 1.0]), "strictly increasing"),
        ],
        ids=["spectrum", "timeseries"],
    )
    def test_unordered_coordinates_are_refused(self, build, match) -> None:
        with pytest.raises(SchemaError, match=match):
            build()

    def test_duplicate_coordinates_are_reported_as_such(self) -> None:
        with pytest.raises(SchemaError, match="1 repeated coordinate"):
            Spectrum([1.0, 1.0, 3.0] * u.um, [1.0, 1.0, 1.0])

    def test_ordering_message_names_the_fix(self) -> None:
        with pytest.raises(SchemaError, match=r"Spectrum\.from_unsorted"):
            Spectrum([2.0, 1.0] * u.um, [1.0, 1.0])

    def test_from_unsorted_sorts_and_carries_everything_along(self) -> None:
        tidy = Spectrum.from_unsorted(
            [3.0, 1.0, 2.0] * u.um,
            [30.0, 10.0, 20.0] * u.Jy,
            uncertainty=[3.0, 1.0, 2.0] * u.Jy,
            mask=np.array([True, False, False]),
        )
        assert tidy.spectral_axis.values.tolist() == [1.0, 2.0, 3.0]
        assert tidy.flux.tolist() == [10.0, 20.0, 30.0]
        assert tidy.uncertainty.tolist() == [1.0, 2.0, 3.0]
        assert tidy.mask.tolist() == [False, False, True]
        # Sorting must not silently strip the units it was handed.
        assert tidy.spectral_axis.unit == u.um
        assert tidy.unit == u.Jy

    def test_from_unsorted_accepts_plain_sequences(self) -> None:
        """It should behave like every other constructor, not fail in numpy."""
        tidy = Spectrum.from_unsorted([3.0, 1.0, 2.0] * u.um, [30.0, 10.0, 20.0])
        assert tidy.flux.tolist() == [10.0, 20.0, 30.0]

    def test_from_unsorted_permutes_extra_coords_too(self) -> None:
        """Per-sample labels must follow their samples through the sort."""
        tidy = Spectrum.from_unsorted(
            [3.0, 1.0, 2.0] * u.um,
            [30.0, 10.0, 20.0],
            extra_coords={"epoch": np.array(["c", "a", "b"])},
        )
        assert tidy.spectral_axis.values.tolist() == [1.0, 2.0, 3.0]
        assert tidy.extra_coords["epoch"].tolist() == ["a", "b", "c"]
        stamped = TimeSeries.from_unsorted(
            [5.0, 1.0] * u.day,
            [0.5, 0.1],
            extra_coords={"visit": np.array([2, 1])},
        )
        assert stamped.extra_coords["visit"].tolist() == [1, 2]

    def test_from_unsorted_converts_a_declared_unit(self) -> None:
        tidy = Spectrum.from_unsorted([3.0, 1.0] * u.um, [2000.0, 1000.0] * u.mJy, unit=u.Jy)
        assert tidy.flux.tolist() == [1.0, 2.0]
        assert tidy.unit == u.Jy

    def test_from_unsorted_still_refuses_duplicates(self) -> None:
        with pytest.raises(SchemaError, match="repeated coordinate"):
            Spectrum.from_unsorted([2.0, 1.0, 2.0] * u.um, [1.0, 1.0, 1.0])

    def test_time_series_from_unsorted(self) -> None:
        curve = TimeSeries.from_unsorted([5.0, 1.0] * u.day, [0.5, 0.1] * u.Jy)
        assert curve.time.values.tolist() == [1.0, 5.0]
        assert curve.values.tolist() == [0.1, 0.5]

    def test_spatial_axes_may_decrease(self) -> None:
        """Sky axes legitimately run either way."""
        sky = Image(
            np.linspace(-2.0, 2.0, 5) * u.arcsec,
            np.linspace(2.0, -2.0, 5) * u.arcsec,
            np.ones((5, 5)),
        )
        assert sky.y.values[0] > sky.y.values[-1]

    def test_non_monotonic_spatial_axis_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="strictly monotonic"):
            Image(
                np.array([0.0, 2.0, 1.0]) * u.arcsec,
                np.linspace(-1.0, 1.0, 2) * u.arcsec,
                np.ones((3, 2)),
            )

    @pytest.mark.parametrize(
        "build",
        [
            lambda: PhotometricPoints(["b", "a"], [3.0, 1.0] * u.um, [1.0, 2.0] * u.Jy),
            lambda: VisibilitySet([5.0, -3.0, 1.0], [2.0, 9.0, -4.0], [1 + 0j, 0 + 1j, 1 + 1j]),
        ],
        ids=["photometry", "visibilities"],
    )
    def test_arbitrary_order_is_tolerated_where_documented(self, build) -> None:
        """These kinds have no natural order; imposing one would be a fiction."""
        assert build().n_samples in (2, 3)


class TestAdvertisedRegularity:
    """Advertised so fast paths can be taken; never required."""

    def test_linear_grid_advertises_regular_and_its_step(self) -> None:
        axis = Spectrum(np.linspace(1.0, 10.0, 10) * u.um, np.ones(10)).spectral_axis
        assert axis.regular
        assert axis.step == pytest.approx(1.0)

    def test_log_grid_advertises_log_regular_only(self) -> None:
        axis = Spectrum(np.geomspace(1.0, 100.0, 20) * u.um, np.ones(20)).spectral_axis
        assert not axis.regular
        assert axis.log_regular
        assert axis.step is None

    def test_irregular_grid_advertises_nothing(self) -> None:
        axis = Spectrum([1.0, 1.5, 9.0] * u.um, np.ones(3)).spectral_axis
        assert not axis.regular
        assert not axis.log_regular

    def test_negative_coordinates_are_never_log_regular(self) -> None:
        axis = Image(
            np.linspace(-2.0, 2.0, 5) * u.arcsec,
            np.linspace(-1.0, 1.0, 2) * u.arcsec,
            np.ones((5, 2)),
        ).x
        assert axis.regular
        assert not axis.log_regular

    def test_container_level_regularity_requires_every_axis(self) -> None:
        regular = Image(
            np.linspace(0.0, 1.0, 3) * u.arcsec,
            np.linspace(0.0, 1.0, 3) * u.arcsec,
            np.ones((3, 3)),
        )
        mixed = Image(
            np.linspace(0.0, 1.0, 3) * u.arcsec,
            np.array([0.0, 0.1, 9.0]) * u.arcsec,
            np.ones((3, 3)),
        )
        assert regular.is_regular
        assert not mixed.is_regular

    def test_short_axes_are_trivially_regular(self) -> None:
        assert Spectrum([1.0] * u.um, [1.0]).spectral_axis.regular


class TestUnits:
    def test_quantity_unit_is_adopted_and_values_become_plain_arrays(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [3.0, 4.0] * u.Jy)
        assert spectrum.unit == u.Jy
        assert type(spectrum.values) is np.ndarray
        assert not isinstance(spectrum.values, u.Quantity)

    def test_explicit_unit_converts_once_at_construction(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1000.0, 2000.0] * u.mJy, unit=u.Jy)
        assert spectrum.values.tolist() == [1.0, 2.0]
        assert spectrum.unit == u.Jy

    def test_incompatible_declared_unit_fails_loudly(self) -> None:
        with pytest.raises(SchemaError, match="not convertible"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 2.0] * u.Jy, unit=u.K)

    def test_axis_units_are_converted_once(self) -> None:
        spectrum = Spectrum([1000.0, 2000.0] * u.nm, [1.0, 1.0])
        assert spectrum.spectral_axis.unit == u.nm
        assert spectrum.spectral_axis.values.tolist() == [1000.0, 2000.0]

    def test_to_unit_converts_values_and_uncertainties(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [3.0, 2.5] * u.Jy, uncertainty=[0.1, 0.2] * u.Jy)
        converted = spectrum.to_unit(u.mJy)
        assert converted.values.tolist() == pytest.approx([3000.0, 2500.0])
        assert converted.uncertainty.tolist() == pytest.approx([100.0, 200.0])
        assert converted.unit == u.mJy
        assert spectrum.unit == u.Jy, "the original is untouched"

    def test_to_unit_without_a_unit_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="has no value unit"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0]).to_unit(u.Jy)

    def test_to_unit_inconvertible_raises_schema_error(self) -> None:
        """The contract's own error type, not a raw astropy exception."""
        with pytest.raises(SchemaError, match="not convertible"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy).to_unit(u.K)

    def test_wrong_physical_type_on_an_axis_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="which this axis does not accept"):
            Spectrum(np.arange(3.0) * u.s, np.ones(3))

    def test_missing_unit_on_a_typed_axis_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="need a unit"):
            Spectrum(np.arange(3.0), np.ones(3))

    @pytest.mark.parametrize(
        "unit", [u.um, u.GHz, u.keV], ids=["wavelength", "frequency", "energy"]
    )
    def test_spectral_axis_accepts_wavelength_frequency_or_energy(self, unit) -> None:
        assert Spectrum([1.0, 2.0] * unit, [1.0, 1.0]).spectral_axis.unit == unit

    def test_value_units_are_deliberately_unconstrained(self) -> None:
        """A model may emit F_lambda, F_nu, a brightness or a bare ratio."""
        for unit in (u.Jy, u.erg / u.s / u.cm**2 / u.AA, u.K, None):
            assert Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], unit=unit).unit == (
                u.Unit(unit) if unit is not None else None
            )

    def test_quantity_helpers_exist_for_reporting(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [3.0, 4.0] * u.Jy)
        assert spectrum.quantity().unit == u.Jy
        assert spectrum.spectral_axis.quantity().unit == u.um


class TestMasks:
    """First-class, and exactly zero information (``prior_art.md`` R1)."""

    @pytest.fixture
    def masked(self) -> Spectrum:
        return Spectrum(
            [1.0, 2.0, 4.0] * u.um,
            [3.0, 2.5, 1.0] * u.Jy,
            uncertainty=[0.1, 0.2, 0.4] * u.Jy,
            mask=np.array([False, True, False]),
        )

    def test_mask_counts_and_complement(self, masked) -> None:
        assert masked.is_masked
        assert masked.n_samples == 3
        assert masked.n_valid == 2
        assert masked.valid.tolist() == [True, False, True]

    def test_zero_weight_representation(self, masked) -> None:
        assert masked.weights().tolist() == [1.0, 0.0, 1.0]

    def test_infinite_uncertainty_representation(self, masked) -> None:
        """R1's other expression of the same statement."""
        assert masked.masked_uncertainty().tolist() == [0.1, np.inf, 0.4]

    def test_the_two_representations_agree(self, masked) -> None:
        weights = masked.weights()
        inflated = masked.masked_uncertainty()
        assert np.array_equal(weights == 0.0, np.isinf(inflated))

    def test_masked_values_stay_in_place_no_sentinel_nans(self, masked) -> None:
        assert masked.values.tolist() == [3.0, 2.5, 1.0]
        assert np.all(np.isfinite(masked.values))

    def test_unmasked_container_is_all_valid(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])
        assert not spectrum.is_masked
        assert spectrum.n_valid == 2
        assert spectrum.weights().tolist() == [1.0, 1.0]

    def test_infinite_uncertainty_form_needs_uncertainties(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], mask=np.array([True, False]))
        with pytest.raises(SchemaError, match="has no uncertainties"):
            spectrum.masked_uncertainty()

    def test_mask_must_be_boolean(self) -> None:
        with pytest.raises(SchemaError, match="must be a boolean array"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], mask=np.array([1, 0]))

    def test_mask_shape_is_checked(self) -> None:
        with pytest.raises(SchemaError, match="mask have shape"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], mask=np.array([True, False, True]))

    def test_masks_work_on_gridded_containers(self) -> None:
        mask = np.zeros((3, 3), dtype=bool)
        mask[1, 1] = True
        sky = Image(
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.ones((3, 3)),
            mask=mask,
        )
        assert sky.n_samples == 9
        assert sky.n_valid == 8
        assert sky.weights().shape == (3, 3)

    def test_a_non_detection_is_not_a_masked_point(self) -> None:
        """§8: censoring carries information, masking does not. The hook is that
        per-sample aligned arrays are supported and ``mask`` is strictly binary."""
        photometry = PhotometricPoints(
            ["WISE_W3", "WISE_W4"],
            [12.1, 22.2] * u.um,
            [0.4, 0.9] * u.Jy,
            uncertainty=[0.05, 0.30] * u.Jy,
            extra_coords={"detected": np.array([True, False])},
        )
        assert not photometry.is_masked
        assert photometry.n_valid == 2
        assert photometry.extra_coords["detected"].tolist() == [True, False]

    def test_mask_is_strictly_binary(self) -> None:
        """No third state and no reserved values, so W1.6 may add censoring
        alongside without a semantic collision."""
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], mask=np.array([True, False]))
        assert spectrum.mask.dtype == np.bool_


class TestHotLoopConstruction:
    def test_with_values_shares_the_validated_axes(self) -> None:
        template = Spectrum(np.geomspace(1.0, 100.0, 50) * u.um, np.zeros(50), unit=u.Jy)
        evaluated = template.with_values(np.ones(50))
        assert evaluated.axes is template.axes
        assert evaluated.unit == u.Jy
        assert evaluated.values.tolist() == [1.0] * 50

    def test_with_values_checks_shape(self) -> None:
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])
        with pytest.raises(SchemaError, match="values have shape"):
            template.with_values(np.ones(3))

    def test_with_values_inherits_mask_and_uncertainty(self) -> None:
        template = Spectrum(
            [1.0, 2.0] * u.um,
            [1.0, 1.0],
            uncertainty=[0.1, 0.2],
            mask=np.array([True, False]),
        )
        evaluated = template.with_values([5.0, 6.0])
        assert evaluated.uncertainty.tolist() == [0.1, 0.2]
        assert evaluated.mask.tolist() == [True, False]

    def test_with_values_can_replace_mask_and_uncertainty(self) -> None:
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], uncertainty=[0.1, 0.2])
        evaluated = template.with_values(
            [5.0, 6.0], uncertainty=[0.5, 0.6], mask=np.array([False, True])
        )
        assert evaluated.uncertainty.tolist() == [0.5, 0.6]
        assert evaluated.mask.tolist() == [False, True]

    def test_with_values_preserves_kind_and_extra_coords(self) -> None:
        template = PhotometricPoints(["a", "b"], [1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        evaluated = template.with_values([3.0, 4.0])
        assert isinstance(evaluated, PhotometricPoints)
        assert evaluated.filters.tolist() == ["a", "b"]

    def test_with_values_converts_a_quantity_rather_than_reinterpreting_it(self) -> None:
        """The units trap in its most damaging form: a silent factor of 1000."""
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        evaluated = template.with_values([2000.0, 3000.0] * u.mJy)
        assert evaluated.values.tolist() == [2.0, 3.0]
        assert evaluated.unit == u.Jy

    def test_with_values_refuses_an_inconvertible_unit(self) -> None:
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0] * u.Jy)
        with pytest.raises(SchemaError, match="not convertible"):
            template.with_values([1.0, 2.0] * u.K)

    def test_with_values_refuses_a_quantity_on_a_unitless_template(self) -> None:
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])
        with pytest.raises(SchemaError, match=r"cannot change .* value unit"):
            template.with_values([1.0, 2.0] * u.Jy)

    def test_with_values_can_retag_fidelity(self) -> None:
        template = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], fidelity="cheap")
        assert template.with_values([2.0, 2.0], fidelity="expensive").fidelity == "expensive"
        assert template.with_values([2.0, 2.0]).fidelity == "cheap"


class TestFidelityTags:
    def test_tag_is_recorded_verbatim(self) -> None:
        assert Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], fidelity="lte").fidelity == "lte"

    def test_absent_by_default(self) -> None:
        assert Spectrum([1.0, 2.0] * u.um, [1.0, 1.0]).fidelity is None

    def test_result_reports_the_fidelity_ladder(self) -> None:
        cheap = Spectrum([1.0, 2.0] * u.um, [1.0, 0.9], fidelity="lte")
        expensive = Spectrum([1.0, 2.0] * u.um, [1.1, 0.8], fidelity="full_nlte")
        result = ModelResult({"a": cheap, "b": expensive})
        assert dict(result.fidelities()) == {"a": "lte", "b": "full_nlte"}

    def test_empty_tag_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="non-empty string"):
            Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], fidelity="")


# ---------------------------------------------------------------------------
# The kinds
# ---------------------------------------------------------------------------


class TestPhotometricPoints:
    def test_filters_are_carried_and_ordered_with_the_samples(self) -> None:
        phot = PhotometricPoints(["2MASS_Ks", "WISE_W1"], [2.16, 3.35] * u.um, [12.0, 8.0] * u.Jy)
        assert phot.filters.tolist() == ["2MASS_Ks", "WISE_W1"]
        assert phot.spectral_axis.values.tolist() == [2.16, 3.35]

    def test_duplicate_filters_are_refused(self) -> None:
        with pytest.raises(SchemaError, match=r"repeated filter names \['W1'\]"):
            PhotometricPoints(["W1", "W1"], [3.4, 3.4] * u.um, [1.0, 1.0] * u.Jy)

    def test_unsorted_pivot_wavelengths_are_fine(self) -> None:
        phot = PhotometricPoints(["b", "a"], [3.0, 1.0] * u.um, [1.0, 2.0] * u.Jy)
        assert phot.spectral_axis.values.tolist() == [3.0, 1.0]

    def test_filter_list_must_be_one_dimensional(self) -> None:
        with pytest.raises(SchemaError, match="one-dimensional"):
            PhotometricPoints(np.array([["a"], ["b"]]), [1.0, 2.0] * u.um, [1.0, 1.0])


class TestGriddedKinds:
    def test_image_shape_follows_its_axes(self) -> None:
        sky = Image(
            np.linspace(-2.0, 2.0, 5) * u.arcsec,
            np.linspace(2.0, -2.0, 4) * u.arcsec,
            np.ones((5, 4)) * u.Jy,
        )
        assert sky.shape == (5, 4)
        assert sky.x.size == 5 and sky.y.size == 4

    def test_image_value_shape_is_checked(self) -> None:
        with pytest.raises(SchemaError, match="coordinates imply"):
            Image(
                np.linspace(-1.0, 1.0, 3) * u.arcsec,
                np.linspace(-1.0, 1.0, 4) * u.arcsec,
                np.ones((4, 3)),
            )

    def test_cube_axis_order_is_x_y_spectral(self) -> None:
        cube = Cube(
            np.linspace(-1.0, 1.0, 3) * u.arcsec,
            np.linspace(-1.0, 1.0, 4) * u.arcsec,
            np.linspace(866.9, 867.0, 6) * u.um,
            np.ones((3, 4, 6)) * u.Jy,
        )
        assert cube.shape == (3, 4, 6)
        assert [axis.name for axis in cube.axes] == ["x", "y", "spectral_axis"]

    def test_cube_spectral_axis_must_increase(self) -> None:
        with pytest.raises(SchemaError, match="strictly increasing"):
            Cube(
                np.linspace(-1.0, 1.0, 2) * u.arcsec,
                np.linspace(-1.0, 1.0, 2) * u.arcsec,
                np.array([867.0, 866.9]) * u.um,
                np.ones((2, 2, 2)),
            )

    def test_grid_axes_need_not_be_evenly_spaced(self) -> None:
        sky = Image(
            np.array([0.0, 0.1, 5.0]) * u.arcsec,
            np.linspace(-1.0, 1.0, 2) * u.arcsec,
            np.ones((3, 2)),
        )
        assert not sky.x.regular
        assert sky.shape == (3, 2)


class TestVisibilitySet:
    """The kind that stresses the schema hardest: complex, scattered, unordered."""

    @pytest.fixture
    def vis(self) -> VisibilitySet:
        return VisibilitySet(
            [120.0, -35.0, 88.0, -210.0],
            [45.0, 190.0, -66.0, 12.0],
            [1.0 + 0.2j, 0.6 - 0.3j, 0.4 + 0.0j, 0.1 - 0.05j],
            uncertainty=[0.02, 0.02, 0.03, 0.05],
            extra_coords={"frequency_ghz": np.array([230.5, 230.5, 345.8, 345.8])},
        )

    def test_values_are_complex(self, vis) -> None:
        assert vis.values.dtype.kind == "c"
        assert vis.n_samples == 4

    def test_scattered_unordered_coordinates_are_accepted(self, vis) -> None:
        assert vis.u.values.tolist() == [120.0, -35.0, 88.0, -210.0]
        assert not vis.u.regular

    def test_amplitude_and_phase_views(self, vis) -> None:
        assert np.allclose(vis.amplitude(), np.abs(vis.values))
        assert np.allclose(vis.phase(), np.angle(vis.values))
        assert vis.amplitude()[0] == pytest.approx(1.019803, rel=1e-5)

    def test_uncertainty_stays_real_on_complex_data(self, vis) -> None:
        """The per-component sigma of a circular complex Gaussian."""
        assert vis.uncertainty.dtype.kind == "f"
        assert vis.uncertainty.tolist() == [0.02, 0.02, 0.03, 0.05]

    def test_complex_uncertainty_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="include complex numbers"):
            VisibilitySet([1.0, 2.0], [3.0, 4.0], [1 + 0j, 0 + 1j], uncertainty=[1 + 1j, 1 + 0j])

    def test_bare_arrays_are_dimensionless_baselines_in_wavelengths(self, vis) -> None:
        assert vis.u.unit == u.dimensionless_unscaled

    def test_spatial_frequency_units_are_accepted(self) -> None:
        vis = VisibilitySet([1.0, 2.0] / u.rad, [3.0, 4.0] / u.rad, [1 + 0j, 0.5 + 0j])
        assert vis.v.unit.is_equivalent(u.rad**-1)

    def test_wrong_axis_unit_is_refused(self) -> None:
        with pytest.raises(SchemaError, match="does not accept"):
            VisibilitySet([1.0, 2.0] * u.m, [3.0, 4.0] * u.m, [1 + 0j, 0 + 1j])

    def test_masks_and_extra_coords_work_on_complex_data(self, vis) -> None:
        masked = vis.with_values(vis.values, mask=np.array([False, False, True, False]))
        assert masked.n_valid == 3
        assert masked.weights().tolist() == [1.0, 1.0, 0.0, 1.0]
        assert masked.extra_coords["frequency_ghz"].tolist() == [230.5, 230.5, 345.8, 345.8]

    def test_complex_values_are_refused_by_real_kinds(self) -> None:
        with pytest.raises(SchemaError, match="holds real values"):
            Spectrum([1.0, 2.0] * u.um, [1.0 + 1j, 2.0 + 0j])

    def test_a_visibility_set_is_an_ordinary_channel(self, vis) -> None:
        result = ModelResult({"uv_230ghz": vis})
        assert result.require("uv_230ghz", VisibilitySet) is vis


class TestUserDefinedKinds:
    """§4.3 makes out-of-tree extensibility a first-class requirement."""

    def test_a_new_kind_is_three_class_attributes(self) -> None:
        class PolarisationCurve(FunctionSamples):
            AXES = (
                AxisSpec(
                    "spectral_axis",
                    physical_types=("length",),
                    order=Order.STRICTLY_INCREASING,
                ),
            )
            LAYOUT = Layout.POINTS

        curve = PolarisationCurve({"spectral_axis": [0.4, 0.6, 0.8] * u.um}, [0.02, 0.03, 0.01])
        assert curve.n_samples == 3
        assert curve.axis("spectral_axis").unit == u.um

    def test_inherited_validation_applies_unchanged(self) -> None:
        class PolarisationCurve(FunctionSamples):
            AXES = (
                AxisSpec(
                    "spectral_axis",
                    physical_types=("length",),
                    order=Order.STRICTLY_INCREASING,
                ),
            )

        with pytest.raises(SchemaError, match="PolarisationCurve requires"):
            PolarisationCurve({"spectral_axis": [0.6, 0.4] * u.um}, [0.02, 0.03])

    def test_a_user_kind_is_an_ordinary_channel(self) -> None:
        class Counts(FunctionSamples):
            AXES = (AxisSpec("channel_number", order=Order.STRICTLY_INCREASING),)

        counts = Counts({"channel_number": [0.0, 1.0, 2.0]}, [10.0, 12.0, 9.0])
        result = ModelResult({"pha": counts})
        assert result.require("pha", Counts) is counts
        assert result.require("pha", FunctionSamples) is counts

    def test_unitless_axes_are_allowed_when_no_physical_type_is_declared(self) -> None:
        class Counts(FunctionSamples):
            AXES = (AxisSpec("channel_number", order=Order.STRICTLY_INCREASING),)

        assert (
            Counts({"channel_number": [0.0, 1.0]}, [1.0, 2.0]).axis("channel_number").unit is None
        )


class TestIntrospection:
    def test_axis_lookup_names_what_exists(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0])
        with pytest.raises(SchemaError, match=r"its axes are \['spectral_axis'\]"):
            spectrum.axis("wavelength")

    def test_reprs_are_informative_and_stable(self) -> None:
        spectrum = Spectrum(
            [1.0, 2.0, 4.0] * u.um,
            [1.0, 1.0, 1.0] * u.Jy,
            mask=np.array([True, False, False]),
            fidelity="cheap",
        )
        text = repr(spectrum)
        assert "Spectrum" in text
        assert "spectral_axis=[1..4] um" in text
        assert "unit=Jy" in text
        assert "masked=1" in text
        assert "fidelity='cheap'" in text

    def test_result_repr_lists_channels_kinds_and_sizes(self, sed, co_windows) -> None:
        text = repr(ModelResult({"sed_lowres": sed, "co_windows": co_windows}))
        assert "sed_lowres: Spectrum@continuum[60]" in text
        assert "co_windows: Spectrum@co_lte[42]" in text

    def test_metadata_is_immutable(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 1.0], meta={"note": "toy"})
        assert spectrum.meta["note"] == "toy"
        with pytest.raises(TypeError):
            spectrum.meta["note"] = "changed"

    def test_integer_inputs_are_promoted_to_float(self) -> None:
        spectrum = Spectrum(np.array([1, 2, 4]) * u.um, np.array([3, 2, 1]))
        assert spectrum.values.dtype.kind == "f"
        assert spectrum.spectral_axis.values.dtype.kind == "f"

    def test_non_numeric_values_are_refused(self) -> None:
        with pytest.raises(SchemaError, match="must be numeric"):
            Spectrum([1.0, 2.0] * u.um, np.array(["a", "b"]))


class TestCoreDependencyFloor:
    """``architecture.md`` §3-4: core is numpy/scipy/astropy/stdlib only."""

    def test_module_imports_no_optional_dependency(self) -> None:
        import sys

        import ampere.core.results_schema  # noqa: F401

        for forbidden in ("torch", "jax", "numpyro", "paramax", "equinox"):
            assert forbidden not in sys.modules, (
                f"importing ampere.core.results_schema pulled in {forbidden}"
            )


class TestAnomalyScore:
    """diagnostics.md §5's shared container, landed at the freeze (R4)."""

    @staticmethod
    def _score(**overrides: object) -> object:
        from ampere.core import AnomalyScore

        settings: dict = {
            "coordinates": np.array([1.0, 2.0, 3.0]),
            "values": np.array([0.1, 2.4, 0.3]),
            "provenance": "gp_localisation_postfit",
            "interpretation_notes": "Amplitude localises deficiency; see the docs.",
        }
        settings.update(overrides)
        return AnomalyScore(**settings)

    def test_construction_and_the_repr(self) -> None:
        score = self._score()
        assert score.n_samples == 3
        assert "gp_localisation_postfit" in repr(score)

    def test_provenance_and_notes_are_required_non_empty(self) -> None:
        # The comparability guard: two differently-computed scores must never
        # travel without saying which family produced them.
        with pytest.raises(SchemaError, match="provenance"):
            self._score(provenance="  ")
        with pytest.raises(SchemaError, match="interpretation_notes"):
            self._score(interpretation_notes="")

    def test_misaligned_shapes_are_refused(self) -> None:
        with pytest.raises(SchemaError, match="indexed by its coordinates"):
            self._score(values=np.array([0.1, 2.4]))
        with pytest.raises(SchemaError, match="mask covers"):
            self._score(mask=np.array([True, False]))

    def test_a_masked_sample_may_carry_a_non_finite_score(self) -> None:
        score = self._score(
            values=np.array([0.1, np.nan, 0.3]), mask=np.array([False, True, False])
        )
        assert score.n_samples == 3
        with pytest.raises(SchemaError, match="finite where retained"):
            self._score(values=np.array([0.1, np.nan, 0.3]))

    def test_it_satisfies_the_results_protocol(self) -> None:
        # ampere.results stays typed against the shape; the class is the shape.
        from ampere.results.plots import AnomalyScoreLike

        assert isinstance(self._score(), AnomalyScoreLike)

    def test_two_dimensional_coordinates_are_legal(self) -> None:
        score = self._score(coordinates=np.array([[1.0, 0.0], [2.0, 1.0], [3.0, 2.0]]))
        assert score.coordinates.shape == (3, 2)


class TestAxisLocate:
    """``Axis.locate``: the lookup ``spectrum_photometry.md`` Gap 1 asks for.

    Ruled by Peter 2026-09-03 and landed with W2.1, whose synthetic-photometry
    step is its first consumer. The rule that matters is the tolerance: the
    negotiated union collapses coordinates coinciding to within
    ``COORDINATE_RTOL`` and keeps one representative, so a step's own published
    coordinate may differ from the survivor by up to that much and must still
    be found.
    """

    def spectral(self, values: object = (1.0, 1.25, 1.65, 2.5, 4.0, 10.0)) -> Axis:
        return Axis.build("spectral_axis", np.asarray(values, dtype=float), u.micron)

    def test_exact_coordinates_are_found_in_the_order_given(self) -> None:
        axis = self.spectral()
        assert axis.locate([2.5, 1.0, 10.0]).tolist() == [3, 0, 5]

    def test_the_result_indexes_the_axis(self) -> None:
        axis = self.spectral()
        wanted = np.array([4.0, 1.25])
        assert axis.values[axis.locate(wanted)].tolist() == wanted.tolist()

    def test_a_coordinate_inside_the_tolerance_still_matches(self) -> None:
        """The whole point: a survivor of the union's dedupe is not bit-identical."""
        axis = self.spectral()
        nudged = 1.65 * (1.0 + COORDINATE_RTOL / 2.0)
        assert axis.locate([nudged]).tolist() == [2]

    def test_a_coordinate_outside_the_tolerance_is_refused(self) -> None:
        axis = self.spectral()
        with pytest.raises(SchemaError, match="COORDINATE_RTOL"):
            axis.locate([1.65 * (1.0 + 1e-6)])

    def test_the_message_names_every_unmatched_value(self) -> None:
        axis = self.spectral()
        with pytest.raises(SchemaError, match=r"99\.0"):
            axis.locate([1.0, 99.0])

    def test_an_unsorted_axis_is_handled(self) -> None:
        """``PhotometricPoints`` declares ``Order.ANY``, and Gap 1 is about photometry."""
        axis = Axis.build("spectral_axis", np.array([3.4, 1.25, 22.0, 4.6]), u.micron)
        assert axis.locate([22.0, 1.25, 3.4]).tolist() == [2, 1, 0]

    def test_locating_nothing_returns_an_empty_index(self) -> None:
        assert self.spectral().locate([]).shape == (0,)

    def test_the_tolerance_is_the_one_negotiation_collapses_with(self) -> None:
        """``locate`` must invert ``_dedupe``, or Gap 1 is only half closed.

        Two coordinates that the union would merge into one must both find that
        survivor; two it would keep apart must not be confused for each other.
        """
        from ampere.core.transform import _dedupe

        published = np.array([2.0, 2.0 * (1.0 + COORDINATE_RTOL / 4.0)])
        survivors = _dedupe(published)
        assert survivors.size == 1
        axis = Axis.build("spectral_axis", survivors, u.micron)
        assert axis.locate(published).tolist() == [0, 0]

    def test_an_empty_axis_refuses_rather_than_returning_nothing(self) -> None:
        axis = Axis.build("spectral_axis", np.array([]), u.micron)
        with pytest.raises(SchemaError, match="empty"):
            axis.locate([1.0])
