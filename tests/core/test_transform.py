"""Unit tests for the W1.5 Transformation & Instrument contract.

Organised around the acceptance criteria in ``WORK_ITEMS.md`` W1.5: the simple
path (fixed grid, no negotiation) is demonstrably trivial
(``TestTheSimplePath``), the out-of-tree extension works
(``TestOutOfTreeExtension``, the headline criterion, loading
``tests/core/thirdparty_polarimeter.py`` the way an installed third-party
package would be loaded), and the pieces those two rest on: kind declaration
and checking, chain composition, parameter composition through
``ParameterSet.merge``, the ratified mask-propagation rule, requirements and
their union, ``negotiate``, and the ``Model`` ABC with its ``compile_for``
hook.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    DEFAULT_CHANNEL,
    AxisRequirement,
    ChannelRequirements,
    Image,
    Instrument,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    Parameterised,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    Transformation,
    negotiate,
    propagate_mask,
)
from ampere.core.exceptions import (
    ChannelError,
    CompositionError,
    ContractError,
    ParameterError,
    SchemaError,
    TransformationError,
)

# ---------------------------------------------------------------------------
# Fixtures and small in-tree transformations used by several tests
# ---------------------------------------------------------------------------


class GreyBody(Model):
    """T * lambda**-beta on a fixed grid: the simple-path model."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.micron)
        self.register_buffer("beta", 1.8)
        self.register_parameter(
            Parameter("temperature", st.uniform(100.0, 900.0), unit=u.K, value=300.0)
        )

    def evaluate(self, **values):
        context = self.context(values)
        flux = context["temperature"] * context["wavelength"] ** -context["beta"]
        return Spectrum(context["wavelength"] * u.micron, flux * u.Jy)


class CalibrationScale(Transformation):
    """A multiplicative flux calibration factor: one free nuisance parameter."""

    ACCEPTS = (Spectrum,)

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.register_parameter(Parameter("scale", st.lognorm(0.1), value=1.0))

    def apply(self, samples, values):
        return samples.with_values(samples.values * self.context(values)["scale"])


class Binner(Transformation):
    """A many-to-one binning resampler; exercises mask propagation."""

    ACCEPTS = (Spectrum,)

    def __init__(self, target, **kwargs) -> None:
        super().__init__(**kwargs)
        self.target = np.asarray(target, dtype=float)

    def weights(self, source: np.ndarray) -> np.ndarray:
        nearest = np.abs(self.target[None, :] - source[:, None]).argmin(axis=1)
        hit = nearest[None, :] == np.arange(self.target.size)[:, None]
        return hit / np.maximum(hit.sum(axis=1, keepdims=True), 1)

    def requirements(self):
        return (
            AxisRequirement(
                "spectral_axis",
                unit=u.um,
                intervals=(float(self.target[0]), float(self.target[-1])),
                max_step=float(np.diff(self.target).min()) / 2.0,
            ),
        )

    def apply(self, samples, values):
        weights = self.weights(samples.spectral_axis.values)
        return Spectrum(
            self.target * u.um,
            weights @ samples.values * samples.unit,
            mask=propagate_mask(samples, weights),
        )


class ToPhotometry(Transformation):
    """A kind-changing step: Spectrum in, PhotometricPoints out."""

    ACCEPTS = (Spectrum,)
    PRODUCES = PhotometricPoints

    def __init__(self, response: np.ndarray, **kwargs) -> None:
        super().__init__(**kwargs)
        self.register_buffer("response", response)

    def apply(self, samples, values):
        weights = self.context(values)["response"]
        weights = weights / weights.sum(axis=1, keepdims=True)
        return PhotometricPoints(
            [f"band{index}" for index in range(weights.shape[0])],
            np.arange(1.0, weights.shape[0] + 1.0) * u.um,
            weights @ samples.values * samples.unit,
            mask=propagate_mask(samples, weights),
        )


@pytest.fixture
def model() -> GreyBody:
    return GreyBody(np.geomspace(1.0, 100.0, 6))


@pytest.fixture
def spectrum() -> Spectrum:
    return Spectrum([1.0, 2.0, 4.0, 8.0] * u.um, [8.0, 4.0, 2.0, 1.0] * u.Jy)


@pytest.fixture
def masked_spectrum() -> Spectrum:
    return Spectrum(
        [1.0, 2.0, 4.0, 8.0] * u.um,
        [8.0, 4.0, 2.0, 1.0] * u.Jy,
        mask=np.array([False, True, False, False]),
    )


# ---------------------------------------------------------------------------
# Acceptance criterion: the simple path is trivial
# ---------------------------------------------------------------------------


class TestTheSimplePath:
    """`DEVELOPMENT_PLAN.md` §4.3: fixed grid, no negotiation, must stay simple."""

    def test_a_whole_instrument_is_one_call(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        predicted = instrument(model(temperature=300.0), {"calibration_scale.scale": 2.0})
        assert isinstance(predicted, Spectrum)
        assert predicted.values[0] == pytest.approx(600.0)

    def test_no_channel_names_are_needed(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        assert instrument.channel == DEFAULT_CHANNEL
        assert model(temperature=300.0).is_single

    def test_nothing_is_negotiated(self) -> None:
        instrument = Instrument([CalibrationScale()])
        assert instrument.requirements() == ()
        assert negotiate([instrument])[DEFAULT_CHANNEL].axes == {}

    def test_a_model_that_ignores_negotiation_still_works(self, model: GreyBody) -> None:
        instrument = Instrument([Binner(np.linspace(2.0, 8.0, 4))])
        requirements = negotiate([instrument])
        assert model.compile_for(requirements) is model
        assert model(temperature=300.0).single().n_samples == 6

    def test_an_instrument_may_have_no_steps(self, model: GreyBody) -> None:
        pass_through = Instrument(channel=DEFAULT_CHANNEL, input_kind=Spectrum)
        result = model(temperature=300.0)
        assert pass_through(result) is result.single()

    def test_declared_values_are_enough(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        assert instrument(model()).values[0] == pytest.approx(300.0)


# ---------------------------------------------------------------------------
# Transformation: declaration, checks, parameters
# ---------------------------------------------------------------------------


class TestTransformationKinds:
    def test_kind_preserving_is_the_default(self) -> None:
        assert CalibrationScale.PRODUCES is None
        assert CalibrationScale.output_kind(Spectrum) is Spectrum

    def test_declared_output_kind_wins(self) -> None:
        assert ToPhotometry.output_kind(Spectrum) is PhotometricPoints

    def test_accepts_kind_honours_subclasses(self) -> None:
        class AnyKind(Transformation):
            def apply(self, samples, values):
                return samples

        assert AnyKind.accepts_kind(Spectrum)
        assert AnyKind.accepts_kind(Image)
        assert not CalibrationScale.accepts_kind(TimeSeries)

    def test_wrong_input_kind_raises_at_the_call(self) -> None:
        series = TimeSeries([0.0, 1.0] * u.day, [1.0, 2.0] * u.Jy)
        with pytest.raises(CompositionError, match="accepts Spectrum, but was given a TimeSeries"):
            CalibrationScale()(series)

    def test_a_non_container_input_raises(self) -> None:
        with pytest.raises(TransformationError, match="not a container"):
            CalibrationScale()(np.arange(3.0))

    def test_wrong_output_kind_raises(self, spectrum: Spectrum) -> None:
        class Miscoded(Transformation):
            ACCEPTS = (Spectrum,)
            PRODUCES = PhotometricPoints

            def apply(self, samples, values):
                return samples

        with pytest.raises(TransformationError, match="declares PRODUCES=PhotometricPoints"):
            Miscoded()(spectrum)

    def test_apply_is_abstract(self) -> None:
        with pytest.raises(TypeError):
            Transformation()  # type: ignore[abstract]

    def test_errors_are_all_contract_errors(self) -> None:
        assert issubclass(TransformationError, ContractError)
        assert issubclass(CompositionError, TransformationError)


class TestTransformationParameters:
    """`parameters.md` §13: declare through `Parameterised`, nothing new."""

    def test_a_transformation_is_parameterised(self) -> None:
        step = CalibrationScale()
        assert isinstance(step, Parameterised)
        assert step.parameters.free_names == ("scale",)

    def test_buffers_work_as_they_do_on_models(self) -> None:
        step = ToPhotometry(np.ones((2, 4)))
        assert step.buffers.names == ("response",)
        assert step.parameters.names == ()

    def test_a_buffer_can_be_promoted_without_touching_apply(self, spectrum: Spectrum) -> None:
        class Convolve(Transformation):
            ACCEPTS = (Spectrum,)

            def __init__(self, width, **kwargs):
                super().__init__(**kwargs)
                self.register_buffer("width", width)

            def apply(self, samples, values):
                return samples.with_values(samples.values * self.context(values)["width"])

        step = Convolve(2.0)
        assert step(spectrum).values[0] == pytest.approx(16.0)
        step.promote_buffer("width", prior=st.uniform(1.0, 3.0))
        assert step.parameters.free_names == ("width",)
        assert step(spectrum, {"width": 3.0}).values[0] == pytest.approx(24.0)

    def test_default_label_is_the_snake_case_class_name(self) -> None:
        assert CalibrationScale().label == "calibration_scale"
        assert ToPhotometry(np.ones((1, 4))).label == "to_photometry"

    def test_an_explicit_label_is_checked(self) -> None:
        assert CalibrationScale(label="detector").label == "detector"
        with pytest.raises(CompositionError, match="not usable"):
            CalibrationScale(label="two words")

    def test_repr_names_the_label_and_the_free_count(self) -> None:
        assert (
            repr(CalibrationScale(label="det")) == "<CalibrationScale 'det', 1 free parameter(s)>"
        )


# ---------------------------------------------------------------------------
# Instrument: composition-time checking
# ---------------------------------------------------------------------------


class TestInstrumentComposition:
    def test_a_valid_chain_reports_its_end_kinds(self) -> None:
        instrument = Instrument([CalibrationScale(), ToPhotometry(np.ones((2, 6)))])
        assert instrument.input_kind is Spectrum
        assert instrument.output_kind is PhotometricPoints

    def test_a_mismatched_chain_raises_when_it_is_built(self) -> None:
        with pytest.raises(CompositionError, match="cannot be composed"):
            Instrument([ToPhotometry(np.ones((2, 6))), CalibrationScale()], label="wrong")

    def test_the_first_step_is_checked_against_the_channel_kind(self) -> None:
        with pytest.raises(CompositionError, match="the channel 'sky' holds a Image"):
            Instrument([CalibrationScale()], channel="sky", input_kind=Image)

    def test_duplicate_step_labels_raise(self) -> None:
        with pytest.raises(CompositionError, match="two steps labelled 'calibration_scale'"):
            Instrument([CalibrationScale(), CalibrationScale()])

    def test_distinct_labels_fix_it(self) -> None:
        instrument = Instrument(
            [CalibrationScale(label="detector"), CalibrationScale(label="aperture")]
        )
        assert instrument.parameters.free_names == ("detector.scale", "aperture.scale")

    def test_input_kind_is_inferred_from_a_single_accepts(self) -> None:
        assert Instrument([CalibrationScale()]).input_kind is Spectrum

    def test_an_ambiguous_first_step_must_be_told(self) -> None:
        class Either(Transformation):
            ACCEPTS = (Spectrum, PhotometricPoints)

            def apply(self, samples, values):
                return samples

        with pytest.raises(CompositionError, match="cannot infer the kind it binds"):
            Instrument([Either()])
        assert Instrument([Either()], input_kind=Spectrum).output_kind is Spectrum

    def test_an_empty_chain_is_the_identity(self, model: GreyBody) -> None:
        instrument = Instrument(input_kind=Spectrum)
        assert instrument.output_kind is Spectrum
        assert instrument(model(temperature=300.0)).n_samples == 6

    def test_a_non_transformation_step_raises(self) -> None:
        with pytest.raises(CompositionError, match="chain of Transformations"):
            Instrument([lambda x: x])  # type: ignore[list-item]

    def test_a_bad_input_kind_raises(self) -> None:
        with pytest.raises(CompositionError, match="FunctionSamples subclass"):
            Instrument(input_kind=float)  # type: ignore[arg-type]

    def test_labels_default_to_the_channel(self) -> None:
        assert Instrument([CalibrationScale()], channel="sed").label == "sed"
        assert Instrument([CalibrationScale()], channel="sed", label="wise").label == "wise"

    def test_repr_shows_the_chain(self) -> None:
        instrument = Instrument([CalibrationScale()], channel="sed")
        assert repr(instrument) == (
            "<Instrument 'sed' on channel 'sed': Spectrum -> calibration_scale -> Spectrum>"
        )


class TestChannelBinding:
    """`results_schema.md` §16: bind with `require`, never `[]`."""

    def test_binding_uses_require_so_a_missing_channel_is_loud(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()], channel="sed_lowres")
        with pytest.raises(ChannelError, match="no channel named 'sed_lowres'"):
            instrument(model(temperature=300.0))

    def test_a_kind_mismatch_on_the_channel_is_loud(self) -> None:
        result = ModelResult(
            {"sky": Image([0.0, 1.0] * u.arcsec, [0.0, 1.0] * u.arcsec, np.ones((2, 2)))}
        )
        instrument = Instrument([CalibrationScale()], channel="sky")
        with pytest.raises(ChannelError, match="holds an Image, but a Spectrum was required"):
            instrument(result)

    def test_channel_error_is_catchable_as_a_schema_error(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()], channel="absent")
        with pytest.raises(SchemaError):
            instrument(model(temperature=300.0))

    def test_binding_a_named_channel(self, model: GreyBody) -> None:
        result = ModelResult({"sed": model(temperature=300.0).single()})
        instrument = Instrument([CalibrationScale()], channel="sed")
        assert instrument(result, {"calibration_scale.scale": 1.0}).n_samples == 6

    def test_bind_refuses_a_bare_container(self, spectrum: Spectrum) -> None:
        with pytest.raises(TransformationError, match="evaluated on a ModelResult"):
            Instrument([CalibrationScale()])(spectrum)  # type: ignore[arg-type]

    def test_the_channel_name_must_be_an_identifier(self) -> None:
        with pytest.raises(CompositionError, match="channel name '2 channels' is not usable"):
            Instrument([CalibrationScale()], channel="2 channels")


class TestParameterComposition:
    """`parameters.md` §13: compose via `merge`, step label as component."""

    def test_names_are_qualified_by_step_label(self) -> None:
        instrument = Instrument([CalibrationScale(), Binner(np.linspace(1.0, 8.0, 4))])
        assert instrument.parameters.names == ("calibration_scale.scale",)

    def test_every_step_is_a_component_even_without_parameters(self) -> None:
        instrument = Instrument([CalibrationScale(), Binner(np.linspace(1.0, 8.0, 4))])
        assert instrument.mapping.components == ("calibration_scale", "binner")

    def test_values_are_routed_to_their_step(self) -> None:
        instrument = Instrument([CalibrationScale(), Binner(np.linspace(1.0, 8.0, 4))])
        routed = instrument.mapping.distribute({"calibration_scale.scale": 1.5})
        assert routed == {"calibration_scale": {"scale": 1.5}, "binner": {}}

    def test_a_flat_vector_is_accepted(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        predicted = instrument(model(temperature=300.0), [2.0])
        assert predicted.values[0] == pytest.approx(600.0)

    def test_a_mapping_is_accepted(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        predicted = instrument(model(temperature=300.0), {"calibration_scale.scale": 2.0})
        assert predicted.values[0] == pytest.approx(600.0)

    def test_the_merge_reflects_later_reconfiguration(self) -> None:
        step = CalibrationScale()
        instrument = Instrument([step])
        assert instrument.parameters.free_size == 1
        step.register_parameter(Parameter("offset", st.norm(0.0, 1.0), value=0.0))
        assert instrument.parameters.free_names == (
            "calibration_scale.scale",
            "calibration_scale.offset",
        )

    def test_the_merged_set_is_an_ordinary_parameter_set(self) -> None:
        instrument = Instrument([CalibrationScale()])
        assert isinstance(instrument.parameters, ParameterSet)
        assert instrument.parameters.lnprior([1.0]) == pytest.approx(
            float(st.lognorm(0.1).logpdf(1.0))
        )

    def test_an_unknown_parameter_name_is_loud(self, model: GreyBody) -> None:
        instrument = Instrument([CalibrationScale()])
        with pytest.raises(ParameterError, match="unknown parameter"):
            instrument(model(temperature=300.0), {"calibration_scale.gain": 2.0})

    def test_nesting_a_chain_inside_a_joint_merge(self) -> None:
        """W1.7's route: merge instruments, then re-distribute per level."""
        instrument = Instrument([CalibrationScale()], channel="sed", label="wise")
        joint = ParameterSet.merge({instrument.label: instrument.parameters})
        assert joint.merged.names == ("wise.calibration_scale.scale",)
        outer = joint.distribute({"wise.calibration_scale.scale": 3.0})
        inner = instrument.mapping.distribute(outer["wise"])
        assert inner["calibration_scale"] == {"scale": 3.0}


# ---------------------------------------------------------------------------
# Masks
# ---------------------------------------------------------------------------


class TestMaskPropagation:
    """`results_schema.md` §16's proposed rule, ratified here."""

    def test_an_unmasked_input_propagates_to_nothing(self, spectrum: Spectrum) -> None:
        assert propagate_mask(spectrum) is None
        assert propagate_mask(None) is None

    def test_one_to_one_passes_the_mask_through(self, masked_spectrum: Spectrum) -> None:
        assert propagate_mask(masked_spectrum).tolist() == [False, True, False, False]

    def test_any_masked_contributor_masks_the_output(self) -> None:
        weights = np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]])
        mask = np.array([False, False, True, False])
        assert propagate_mask(mask, weights).tolist() == [False, True]

    def test_a_zero_weight_is_not_influence(self) -> None:
        weights = np.array([[1.0, 0.0]])
        assert propagate_mask(np.array([False, True]), weights).tolist() == [False]

    def test_the_rule_is_any_not_all(self) -> None:
        """The decision this contract had to make, asserted explicitly."""
        weights = np.ones((1, 3))
        assert propagate_mask(np.array([True, False, False]), weights).tolist() == [True]

    def test_a_mismatched_influence_matrix_is_loud(self) -> None:
        with pytest.raises(TransformationError, match=r"\(n_out, n_in\) influence matrix"):
            propagate_mask(np.array([True, False]), np.ones((2, 3)))

    def test_with_values_propagates_for_free(self, masked_spectrum: Spectrum) -> None:
        out = CalibrationScale()(masked_spectrum, {"scale": 2.0})
        assert out.mask.tolist() == [False, True, False, False]

    def test_dropping_a_mask_raises(self, masked_spectrum: Spectrum) -> None:
        class Forgetful(Transformation):
            ACCEPTS = (Spectrum,)

            def apply(self, samples, values):
                return Spectrum(samples.spectral_axis.quantity(), samples.values * u.Jy)

        with pytest.raises(TransformationError, match="dropped the mask"):
            Forgetful()(masked_spectrum)

    def test_an_unmasked_input_needs_no_mask_on_output(self, spectrum: Spectrum) -> None:
        class Forgetful(Transformation):
            ACCEPTS = (Spectrum,)

            def apply(self, samples, values):
                return Spectrum(samples.spectral_axis.quantity(), samples.values * u.Jy)

        assert Forgetful()(spectrum).mask is None

    def test_an_explicit_all_false_mask_is_the_escape_hatch(
        self, masked_spectrum: Spectrum
    ) -> None:
        class Independent(Transformation):
            ACCEPTS = (Spectrum,)

            def apply(self, samples, values):
                return Spectrum(
                    samples.spectral_axis.quantity(),
                    samples.values * u.Jy,
                    mask=np.zeros(samples.n_samples, dtype=bool),
                )

        assert Independent()(masked_spectrum).is_masked is False

    def test_a_binning_step_masks_the_touched_bin(self, masked_spectrum: Spectrum) -> None:
        binner = Binner(np.array([1.5, 6.0]))
        out = binner(masked_spectrum)
        assert out.mask.tolist() == [True, False]


# ---------------------------------------------------------------------------
# Requirements
# ---------------------------------------------------------------------------


class TestAxisRequirement:
    def test_intervals_adopt_a_quantity_unit(self) -> None:
        requirement = AxisRequirement("spectral_axis", intervals=(1.0, 30.0) * u.um)
        assert requirement.unit == u.um
        assert requirement.intervals == ((1.0, 30.0),)

    def test_several_intervals_are_merged_when_they_overlap(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis", intervals=[(1.0, 5.0), (4.0, 9.0), (20.0, 21.0)] * u.um
        )
        assert requirement.intervals == ((1.0, 9.0), (20.0, 21.0))

    def test_disjoint_intervals_stay_disjoint(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis", intervals=[(866.9, 867.0), (1300.3, 1300.5)] * u.um
        )
        assert len(requirement.intervals) == 2

    def test_a_reversed_interval_is_a_typo_and_raises(self) -> None:
        with pytest.raises(TransformationError, match="wrong way round"):
            AxisRequirement("spectral_axis", intervals=(30.0, 1.0) * u.um)

    def test_a_degenerate_interval_points_at_points(self) -> None:
        with pytest.raises(TransformationError, match="pass points="):
            AxisRequirement("spectral_axis", intervals=(2.0, 2.0) * u.um)

    def test_points_are_sorted_and_deduplicated(self) -> None:
        requirement = AxisRequirement("spectral_axis", points=[3.0, 1.0, 3.0, 2.0] * u.um)
        assert requirement.points.tolist() == [1.0, 2.0, 3.0]

    def test_density_must_be_positive(self) -> None:
        with pytest.raises(TransformationError, match="positive finite"):
            AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, max_step=0.0)
        with pytest.raises(TransformationError, match="positive finite"):
            AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, min_resolving_power=-1.0)

    def test_an_axis_name_must_be_an_identifier(self) -> None:
        with pytest.raises(CompositionError, match="axis name"):
            AxisRequirement("spectral axis")

    def test_conversion_is_exact(self) -> None:
        nanometres = AxisRequirement(
            "spectral_axis", intervals=(2000.0, 30000.0) * u.nm, max_step=100.0 * u.nm
        )
        microns = nanometres.convert_to(u.um)
        assert microns.intervals == ((2.0, 30.0),)
        assert microns.segments()[0][2] == pytest.approx(0.1)

    def test_inconvertible_units_are_refused(self) -> None:
        with pytest.raises(CompositionError, match="cannot be converted"):
            AxisRequirement("spectral_axis", unit=u.um, intervals=(1.0, 2.0) * u.s)

    def test_a_malformed_interval_shape_is_loud(self) -> None:
        with pytest.raises(TransformationError, match=r"shape \(2,\) or \(n, 2\)"):
            AxisRequirement("spectral_axis", intervals=(1.0, 2.0, 3.0) * u.um)
        with pytest.raises(TransformationError, match=r"shape \(2,\) or \(n, 2\)"):
            AxisRequirement("spectral_axis", intervals=5.0)

    def test_per_element_quantities_are_refused_rather_than_stripped(self) -> None:
        with pytest.raises(TransformationError, match="one astropy Quantity covering all of them"):
            AxisRequirement("spectral_axis", intervals=[(1.0 * u.um, 2.0 * u.um)])

    def test_an_empty_requirement_is_legal(self) -> None:
        blank = AxisRequirement("spectral_axis")
        assert blank.intervals == ()
        assert blank.points is None
        assert not blank.constrains_density

    def test_max_step_grid_satisfies_the_step(self) -> None:
        requirement = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, max_step=0.1)
        grid = requirement.coordinates().to_value(u.um)
        assert np.max(np.diff(grid)) <= 0.1 + 1e-12
        assert grid[0] == pytest.approx(1.0)
        assert grid[-1] == pytest.approx(2.0)

    def test_resolving_power_grid_satisfies_the_ratio(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis", intervals=(1.0, 30.0) * u.um, min_resolving_power=40.0
        )
        grid = requirement.coordinates().to_value(u.um)
        assert np.max(grid[1:] / grid[:-1]) <= 1.0 + 1.0 / 40.0 + 1e-12

    def test_both_density_constraints_are_satisfied_together(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis",
            intervals=(1.0, 30.0) * u.um,
            min_resolving_power=40.0,
            max_step=0.2,
        )
        grid = requirement.coordinates().to_value(u.um)
        assert np.max(np.diff(grid)) <= 0.2 + 1e-12
        assert np.max(grid[1:] / grid[:-1]) <= 1.0 + 1.0 / 40.0 + 1e-12

    def test_required_points_appear_in_the_grid(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis", intervals=(1.0, 2.0) * u.um, points=[1.234] * u.um, max_step=0.5
        )
        grid = requirement.coordinates().to_value(u.um)
        assert np.any(np.isclose(grid, 1.234))

    def test_the_grid_is_strictly_increasing_so_a_container_accepts_it(self) -> None:
        requirement = AxisRequirement(
            "spectral_axis", intervals=[(1.0, 2.0), (2.0, 3.0)] * u.um, max_step=0.25
        )
        grid = requirement.coordinates()
        assert Spectrum(grid, np.zeros(grid.size) * u.Jy).n_samples == grid.size

    def test_a_requirement_with_no_geometry_cannot_build_a_grid(self) -> None:
        with pytest.raises(TransformationError, match="no intervals and no points"):
            AxisRequirement("spectral_axis").coordinates()

    def test_a_density_without_coverage_is_refused_at_construction(self) -> None:
        """A density with nowhere to apply would be lost, or over-applied, in a union."""
        with pytest.raises(TransformationError, match="no coverage for it to apply to"):
            AxisRequirement("spectral_axis", max_step=0.1)

    def test_resolving_power_needs_positive_coordinates(self) -> None:
        with pytest.raises(TransformationError, match="strictly positive"):
            AxisRequirement("time", intervals=(0.0, 1.0), min_resolving_power=10.0).coordinates()

    def test_a_unitless_axis_returns_a_bare_array(self) -> None:
        requirement = AxisRequirement("u", intervals=(0.0, 100.0), max_step=50.0)
        grid = requirement.coordinates()
        assert not isinstance(grid, u.Quantity)
        assert grid.tolist() == [0.0, 50.0, 100.0]


class TestRequirementUnion:
    def test_coverage_accumulates(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 30.0) * u.um)
        right = AxisRequirement("spectral_axis", intervals=(20.0, 200.0) * u.um)
        assert left.union(right).intervals == ((1.0, 200.0),)

    def test_gaps_survive(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(866.9, 867.0) * u.um)
        right = AxisRequirement("spectral_axis", intervals=(1300.3, 1300.5) * u.um)
        assert left.union(right).intervals == ((866.9, 867.0), (1300.3, 1300.5))

    def test_the_stricter_density_wins(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, max_step=0.5)
        right = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, max_step=0.1)
        both = left.union(right)
        assert both.segments() == ((1.0, 2.0, 0.5, None), (1.0, 2.0, 0.1, None))
        assert np.diff(np.asarray(both.coordinates())).max() <= 0.1 + 1e-12

    def test_the_higher_resolving_power_wins(self) -> None:
        left = AxisRequirement(
            "spectral_axis", intervals=(1.0, 2.0) * u.um, min_resolving_power=10.0
        )
        right = AxisRequirement(
            "spectral_axis", intervals=(1.0, 2.0) * u.um, min_resolving_power=1000.0
        )
        grid = np.asarray(left.union(right).coordinates())
        assert (grid[:-1] / np.diff(grid)).min() >= 1000.0 - 1e-6

    def test_points_accumulate(self) -> None:
        left = AxisRequirement("spectral_axis", points=[1.0, 2.0] * u.um)
        right = AxisRequirement("spectral_axis", points=[2.0, 3.0] * u.um)
        assert left.union(right).points.tolist() == [1.0, 2.0, 3.0]

    def test_each_input_is_satisfied_over_its_own_coverage(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 10.0) * u.um, max_step=1.0)
        right = AxisRequirement("spectral_axis", intervals=(5.0, 20.0) * u.um, max_step=0.5)
        grid = left.union(right).coordinates().to_value(u.um)
        assert grid[0] <= 1.0
        assert grid[-1] >= 20.0
        coarse = grid[(grid >= 1.0) & (grid <= 10.0)]
        fine = grid[(grid >= 5.0) & (grid <= 20.0)]
        assert np.max(np.diff(coarse)) <= 1.0 + 1e-12
        assert np.max(np.diff(fine)) <= 0.5 + 1e-12

    def test_a_fine_window_does_not_refine_a_coarse_neighbour(self) -> None:
        """The defect the per-interval density exists to prevent."""
        broad = AxisRequirement("spectral_axis", intervals=(1.0, 200.0) * u.um, max_step=1.0)
        window = AxisRequirement("spectral_axis", intervals=(100.0, 100.1) * u.um, max_step=0.001)
        both = broad.union(window)
        assert both.intervals == ((1.0, 200.0), (100.0, 100.1))
        assert both.segments() == (
            (1.0, 200.0, 1.0, None),
            (100.0, 100.1, 0.001, None),
        )
        # ~200 from the broad interval plus ~101 from the window, not 200 000.
        assert both.coordinates().size < 400

    def test_intervals_of_equal_density_still_merge(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 10.0) * u.um, max_step=0.5)
        right = AxisRequirement("spectral_axis", intervals=(5.0, 20.0) * u.um, max_step=0.5)
        assert left.union(right).intervals == ((1.0, 20.0),)

    def test_the_density_scalars_do_not_survive_as_attributes(self) -> None:
        # Ruled 2026-09-02 (spec §15.8): after a union a scalar could only be
        # a strictest-anywhere summary that misreads as global; segments() is
        # the only public statement of density.
        broad = AxisRequirement("spectral_axis", intervals=(1.0, 200.0) * u.um, max_step=1.0)
        window = AxisRequirement("spectral_axis", intervals=(100.0, 100.1) * u.um, max_step=0.001)
        both = broad.union(window)
        assert both.segments() == ((1.0, 200.0, 1.0, None), (100.0, 100.1, 0.001, None))
        with pytest.raises(AttributeError):
            both.max_step  # noqa: B018
        with pytest.raises(AttributeError):
            both.min_resolving_power  # noqa: B018

    def test_units_are_reconciled(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um)
        right = AxisRequirement("spectral_axis", intervals=(1500.0, 4000.0) * u.nm)
        assert left.union(right).intervals == ((1.0, 4.0),)

    def test_a_unitless_requirement_never_mixes_with_a_unit_bearing_one(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um)
        right = AxisRequirement("spectral_axis", intervals=(1.0, 2.0))
        with pytest.raises(CompositionError, match="disagree about units"):
            left.union(right)

    def test_different_axes_cannot_be_unioned(self) -> None:
        with pytest.raises(CompositionError, match="different axes"):
            AxisRequirement("time").union(AxisRequirement("spectral_axis"))

    def test_sources_are_carried_and_deduplicated(self) -> None:
        left = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, source="a")
        right = AxisRequirement("spectral_axis", intervals=(1.0, 2.0) * u.um, source="b")
        assert left.union(right).source == "a, b"


class TestNegotiate:
    def test_every_bound_channel_appears(self) -> None:
        first = Instrument([CalibrationScale()], channel="sed", label="one")
        second = Instrument([CalibrationScale()], channel="lines", label="two")
        assert sorted(negotiate([first, second])) == ["lines", "sed"]

    def test_requirements_are_unioned_across_instruments(self) -> None:
        left = Instrument([Binner(np.linspace(1.0, 5.0, 5))], channel="sed", label="a")
        right = Instrument([Binner(np.linspace(4.0, 20.0, 9))], channel="sed", label="b")
        asked = negotiate([left, right])["sed"]["spectral_axis"]
        assert asked.segments() == ((1.0, 5.0, 0.5, None), (4.0, 20.0, 1.0, None))

    def test_requirements_are_unioned_across_steps_of_one_instrument(self) -> None:
        instrument = Instrument(
            [
                Binner(np.linspace(1.0, 5.0, 5), label="coarse"),
                Binner(np.linspace(4.0, 20.0, 9), label="fine"),
            ],
            channel="sed",
        )
        asked = instrument.requirements()
        assert len(asked) == 2
        assert negotiate([instrument])["sed"]["spectral_axis"].intervals == (
            (1.0, 5.0),
            (4.0, 20.0),
        )

    def test_the_source_defaults_to_instrument_dot_step(self) -> None:
        instrument = Instrument([Binner(np.linspace(1.0, 5.0, 5))], channel="sed", label="irs")
        assert instrument.requirements()[0].source == "irs.binner"

    def test_the_most_specific_kind_wins(self) -> None:
        class AnyKind(Transformation):
            def apply(self, samples, values):
                return samples

        loose = Instrument([AnyKind()], channel="sed", label="loose", input_kind=Spectrum)
        general = Instrument(channel="sed", label="general")
        assert negotiate([general, loose])["sed"].kind is Spectrum

    def test_unrelated_kinds_on_one_channel_raise(self) -> None:
        spectral = Instrument([CalibrationScale()], channel="sed", label="a")
        imaging = Instrument(channel="sed", label="b", input_kind=Image)
        with pytest.raises(CompositionError, match="expect unrelated kinds"):
            negotiate([spectral, imaging])

    def test_a_requirement_on_a_nonexistent_axis_raises(self) -> None:
        class WrongAxis(Transformation):
            ACCEPTS = (Spectrum,)

            def requirements(self):
                return (AxisRequirement("time", intervals=(0.0, 1.0) * u.day),)

            def apply(self, samples, values):
                return samples

        with pytest.raises(CompositionError, match="is indexed by"):
            negotiate([Instrument([WrongAxis()], channel="sed", label="clock")])

    def test_a_non_requirement_publication_raises(self) -> None:
        class Sloppy(Transformation):
            ACCEPTS = (Spectrum,)

            def requirements(self):
                return ("1 to 30 microns",)  # type: ignore[return-value]

            def apply(self, samples, values):
                return samples

        with pytest.raises(CompositionError, match="not an AxisRequirement"):
            Instrument([Sloppy()], channel="sed").requirements()

    def test_negotiate_takes_instruments(self) -> None:
        with pytest.raises(CompositionError, match="takes Instruments"):
            negotiate([CalibrationScale()])  # type: ignore[list-item]

    def test_channel_requirements_names_a_missing_axis(self) -> None:
        asked = negotiate([Instrument([CalibrationScale()], channel="sed")])["sed"]
        assert isinstance(asked, ChannelRequirements)
        assert "spectral_axis" not in asked
        with pytest.raises(CompositionError, match="no requirement was published"):
            asked["spectral_axis"]

    def test_coordinates_are_produced_per_constrained_axis(self) -> None:
        instrument = Instrument([Binner(np.linspace(1.0, 5.0, 5))], channel="sed")
        coordinates = negotiate([instrument])["sed"].coordinates()
        assert set(coordinates) == {"spectral_axis"}
        assert coordinates["spectral_axis"].unit == u.um


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class NegotiatingModel(Model):
    """A model that honours requests: templates built once, refilled thereafter."""

    def __init__(self, channels: dict[str, np.ndarray]) -> None:
        self.register_parameter(
            Parameter("temperature", st.uniform(50.0, 450.0), unit=u.K, value=300.0)
        )
        self.templates = {name: self._template(grid) for name, grid in channels.items()}

    @staticmethod
    def _template(wavelength: np.ndarray) -> Spectrum:
        return Spectrum(np.asarray(wavelength) * u.um, np.zeros(len(wavelength)) * u.Jy)

    def compile_for(self, requirements):
        for channel, asked in requirements.items():
            if "spectral_axis" in asked:
                self.templates[channel] = self._template(
                    asked["spectral_axis"].coordinates().to_value(u.um)
                )
        return self

    def evaluate(self, **values):
        return ModelResult(
            {
                name: template.with_values(
                    values["temperature"] * template.spectral_axis.values**-1.8
                )
                for name, template in self.templates.items()
            }
        )


class TestModelABC:
    def test_a_bare_container_is_filed_under_the_default_channel(self, model: GreyBody) -> None:
        result = model(temperature=300.0)
        assert isinstance(result, ModelResult)
        assert list(result) == [DEFAULT_CHANNEL]

    def test_a_flat_vector_is_accepted(self, model: GreyBody) -> None:
        assert model([400.0]).single().values[0] == pytest.approx(400.0)

    def test_a_mapping_is_accepted(self, model: GreyBody) -> None:
        assert model({"temperature": 400.0}) == model(temperature=400.0)

    def test_declared_values_are_the_fallback(self, model: GreyBody) -> None:
        assert model().single().values[0] == pytest.approx(300.0)

    def test_a_vector_and_keywords_together_are_refused(self, model: GreyBody) -> None:
        with pytest.raises(TransformationError, match="pass one or the other"):
            model([400.0], temperature=500.0)

    def test_an_unknown_keyword_is_loud(self, model: GreyBody) -> None:
        with pytest.raises(ParameterError, match="unknown parameter"):
            model(temperture=300.0)

    def test_evaluate_is_abstract(self) -> None:
        with pytest.raises(TypeError):
            Model()  # type: ignore[abstract]

    def test_a_model_returning_rubbish_is_loud(self) -> None:
        class Rubbish(Model):
            def evaluate(self, **values):
                return np.arange(3.0)

        with pytest.raises(TransformationError, match="A model returns a ModelResult"):
            Rubbish()()

    def test_compile_for_defaults_to_the_identity(self, model: GreyBody) -> None:
        assert model.compile_for({}) is model

    def test_call_attaches_the_resolved_theta(self, model: GreyBody) -> None:
        """results_schema.md §17 question 7, ruled 2026-09-01."""
        record = model(temperature=400.0).parameters
        assert record is not None and record["temperature"] == 400.0
        vector_record = model([250.0]).parameters
        assert vector_record is not None and vector_record["temperature"] == 250.0

    def test_a_record_evaluate_attached_is_respected(self) -> None:
        class SelfTagging(Model):
            def evaluate(self, **values):
                return ModelResult(
                    Spectrum([1.0, 2.0] * u.um, np.ones(2) * u.Jy),
                    parameters={"note": "mine"},
                )

        assert SelfTagging()().parameters == {"note": "mine"}


class TestCompileOnceEvaluateMany:
    """`results_schema.md` §16's compile-once/evaluate-many split."""

    def test_negotiation_reconfigures_the_grids(self) -> None:
        model = NegotiatingModel({"sed": np.geomspace(1.0, 100.0, 6)})
        instrument = Instrument([Binner(np.linspace(2.0, 20.0, 19))], channel="sed")
        before = model(temperature=300.0)["sed"].n_samples
        compiled = model.compile_for(negotiate([instrument]))
        after = compiled(temperature=300.0)["sed"].n_samples
        assert before == 6
        assert after > before

    def test_the_hot_loop_reuses_the_validated_axes(self) -> None:
        model = NegotiatingModel({"sed": np.geomspace(1.0, 100.0, 6)})
        first = model(temperature=300.0)["sed"]
        second = model(temperature=200.0)["sed"]
        assert second.axes is first.axes
        assert second.values[0] / first.values[0] == pytest.approx(2.0 / 3.0)

    def test_the_instrument_consumes_the_negotiated_grid(self) -> None:
        model = NegotiatingModel({"sed": np.geomspace(1.0, 100.0, 6)})
        instrument = Instrument([Binner(np.linspace(2.0, 20.0, 19))], channel="sed")
        compiled = model.compile_for(negotiate([instrument]))
        predicted = instrument(compiled(temperature=300.0))
        assert predicted.n_samples == 19

    def test_disjoint_windows_stay_disjoint_end_to_end(self) -> None:
        model = NegotiatingModel({"co": np.linspace(866.9, 867.0, 5)})
        windows = AxisRequirement(
            "spectral_axis", intervals=[(866.9, 867.0), (1300.3, 1300.5)] * u.um, max_step=0.01
        )

        class Windows(Transformation):
            ACCEPTS = (Spectrum,)

            def requirements(self):
                return (windows,)

            def apply(self, samples, values):
                return samples

        compiled = model.compile_for(negotiate([Instrument([Windows()], channel="co")]))
        grid = compiled(temperature=300.0)["co"].spectral_axis.values
        assert float(np.diff(grid).max()) > 400.0


# ---------------------------------------------------------------------------
# The acceptance criterion: an out-of-tree extension
# ---------------------------------------------------------------------------


def _load_as_third_party(name: str) -> ModuleType:
    """Import a module by path, exactly as an installed package would be imported.

    Deliberately not a package-relative import: the point of the exercise is
    that the module has no relationship to ampere's package structure and needs
    none, so it is loaded from its file and nothing else.
    """
    path = Path(__file__).with_name(f"{name}.py")
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


polarimetry = _load_as_third_party("thirdparty_polarimeter")


class TestOutOfTreeExtension:
    """`DEVELOPMENT_PLAN.md` §4.3's explicit acceptance criterion."""

    @pytest.fixture
    def instrument(self) -> Instrument:
        return polarimetry.build_polarimeter(0.1, [1.0, 4.0, 16.0, 100.0], channel="sed")

    @pytest.fixture
    def result(self, model: GreyBody) -> ModelResult:
        return ModelResult({"sed": model(temperature=300.0).single()})

    def test_the_extension_module_imports_only_the_public_api(self) -> None:
        source = Path(polarimetry.__file__).read_text(encoding="utf-8")
        statements = [line for line in source.splitlines() if line.startswith(("import ", "from "))]
        assert [line for line in statements if "ampere" in line] == ["from ampere.core import ("]

    def test_a_user_container_kind_composes_with_the_chain(
        self, instrument: Instrument, result: ModelResult
    ) -> None:
        predicted = instrument(result)
        assert isinstance(predicted, polarimetry.PolarisationCurve)
        assert predicted.n_samples == 3

    def test_the_chain_declares_its_end_kinds(self, instrument: Instrument) -> None:
        assert instrument.input_kind is Spectrum
        assert instrument.output_kind is polarimetry.PolarisationCurve

    def test_the_users_nuisance_parameters_merge_like_any_other(
        self, instrument: Instrument
    ) -> None:
        assert instrument.parameters.names == (
            "atmospheric_depolarisation.airmass",
            "polarimeter_channels.instrumental",
        )
        assert instrument.parameters.free_size == 0

    def test_a_users_parameter_can_be_freed_at_composition_time(self, result: ModelResult) -> None:
        instrument = polarimetry.build_polarimeter(
            0.1, [1.0, 4.0, 16.0, 100.0], channel="sed", prior=st.norm(0.0, 0.01)
        )
        assert instrument.parameters.free_names == ("polarimeter_channels.instrumental",)
        offset = 0.5
        baseline = instrument(result, {"polarimeter_channels.instrumental": 0.0})
        shifted = instrument(result, {"polarimeter_channels.instrumental": offset})
        assert (shifted.values - baseline.values) == pytest.approx(offset)

    def test_the_users_transformation_publishes_requirements(self, instrument: Instrument) -> None:
        asked = negotiate([instrument])["sed"]["spectral_axis"]
        assert asked.intervals == ((1.0, 100.0),)
        assert asked.segments()[0][2] == pytest.approx(0.75)

    def test_the_users_transformation_propagates_masks(self, model: GreyBody) -> None:
        flux = model(temperature=300.0).single()
        masked = ModelResult(
            {
                "sed": flux.with_values(
                    flux.values, mask=np.array([True, False, False, False, False, False])
                )
            }
        )
        instrument = polarimetry.build_polarimeter(0.1, [1.0, 4.0, 16.0, 100.0], channel="sed")
        predicted = instrument(masked)
        assert predicted.mask.tolist() == [True, False, False]

    def test_a_users_buffer_can_be_promoted_without_touching_their_code(
        self, instrument: Instrument, result: ModelResult
    ) -> None:
        step = instrument.steps[0]
        step.demote_parameter("airmass")
        step.promote_buffer("airmass", prior=st.uniform(1.0, 2.0))
        assert instrument.parameters.free_names == ("atmospheric_depolarisation.airmass",)
        thin = instrument(result, {"atmospheric_depolarisation.airmass": 1.0})
        thick = instrument(result, {"atmospheric_depolarisation.airmass": 2.0})
        assert np.all(thick.values < thin.values)

    def test_nothing_inside_ampere_had_to_change(self) -> None:
        """The interfaces suffice: the extension uses only exported names."""
        import ampere.core

        used = {
            "AxisRequirement",
            "AxisSpec",
            "FunctionSamples",
            "Instrument",
            "Layout",
            "Order",
            "Parameter",
            "Spectrum",
            "Transformation",
            "propagate_mask",
        }
        assert used <= set(ampere.core.__all__)
