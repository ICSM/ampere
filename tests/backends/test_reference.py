"""Unit tests for ``ampere.backends.reference``.

The conformance battery already holds the shipped ``PowerLaw``,
``CalibrationScale``, ``Resample``, ``FractionalModelNoise`` and
``FractionalModelGPNoise`` to the §4 contracts, once per registered backend,
and those rows are not repeated here. What this file covers is what the battery
cannot reach:

* the **physics** — the Planck function against an independent oracle, and the
  limits where the arithmetic is delicate;
* the two steps with no conformance kind — :class:`LSFConvolution` and the real
  :class:`SyntheticPhotometry`, including its pyphot route;
* ``spectrum_photometry.md`` **Gap 1's own scenario**, end to end: a second
  instrument bound to one channel, which is the failure ``Axis.locate`` exists
  to prevent.

Oracles are analytic or ``astropy``/``scipy``, never a second copy of the code
under test.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import (
    BlackBody,
    CalibrationScale,
    FractionalModelGPNoise,
    FractionalModelNoise,
    LSFConvolution,
    ModifiedBlackBody,
    PowerLaw,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    planck_jy,
)
from ampere.core import (
    CompositionError,
    GaussianFamily,
    Instrument,
    Likelihood,
    LikelihoodError,
    Matern32,
    Parameter,
    ParameterError,
    SchemaError,
    Spectrum,
    TransformationError,
    negotiate,
)

GRID = np.geomspace(1.0, 30.0, 400)


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------


class TestPlanckFunction:
    """``planck_jy`` against astropy's own implementation, and at the limits."""

    def test_it_agrees_with_astropys_blackbody(self) -> None:
        """A genuinely independent implementation, in different internal units."""
        from astropy.modeling.physical_models import BlackBody as AstropyBlackBody

        for temperature in (30.0, 300.0, 5000.0):
            oracle = AstropyBlackBody(temperature=temperature * u.K)(GRID * u.micron).to_value(
                u.Jy / u.sr
            )
            assert planck_jy(GRID, temperature) == pytest.approx(oracle, rel=1e-12)

    def test_the_rayleigh_jeans_tail_is_not_lost_to_cancellation(self) -> None:
        """``expm1``, not ``exp(x) - 1``: at ``h nu << kT`` the latter has no digits left.

        Far enough down the tail, ``B_nu -> 2 nu**2 k T / c**2``. The naive
        difference loses precision long before this wavelength.
        """
        far = np.array([1.0e6, 1.0e7])  # micron; h nu / kT is ~1e-4 here
        frequency = (299792458.0 * 1e6) / far
        rayleigh_jeans = 2.0 * frequency**2 * 1.380649e-23 * 100.0 / 299792458.0**2 / 1e-26
        assert planck_jy(far, 100.0) == pytest.approx(rayleigh_jeans, rel=1e-3)

    def test_the_wien_tail_underflows_to_zero_rather_than_nan(self) -> None:
        """``expm1`` overflows to inf there; the ratio must be 0, not a nan."""
        emitted = planck_jy(np.array([1.0e-4, 1.0e-3]), 3.0)
        assert np.all(np.isfinite(emitted))
        assert np.all(emitted == 0.0)

    def test_a_non_positive_temperature_is_refused(self) -> None:
        with pytest.raises(ValueError, match="finite and positive"):
            planck_jy(GRID, -1.0)

    def test_a_non_positive_wavelength_is_refused(self) -> None:
        with pytest.raises(ValueError, match="strictly positive"):
            planck_jy(np.array([0.0, 1.0]), 300.0)


class TestNativeModels:
    """The three models ``DEVELOPMENT_PLAN.md`` §5 has this backend ship."""

    def test_a_blackbody_declares_what_it_was_given(self) -> None:
        model = BlackBody(GRID, temperature=st.norm(300.0, 50.0), scale=1e-6)
        assert model.parameters.free_names == ("temperature",)
        assert model.parameters["scale"].is_fixed
        assert model.buffers["wavelength"].unit == u.micron

    def test_a_blackbody_emits_jy_on_its_own_grid(self) -> None:
        emitted = BlackBody(GRID, temperature=300.0, scale=1.0)().single()
        assert emitted.unit == u.Jy
        assert emitted.spectral_axis.values == pytest.approx(GRID)

    def test_the_scale_is_a_plain_multiplier(self) -> None:
        one = BlackBody(GRID, temperature=300.0, scale=1.0)().single().values
        two = BlackBody(GRID, temperature=300.0, scale=2.0)().single().values
        assert two == pytest.approx(2.0 * one)

    def test_a_modified_blackbody_is_the_planck_function_times_the_emissivity(self) -> None:
        reference, beta = 250.0, 1.8
        emitted = (
            ModifiedBlackBody(
                GRID, temperature=40.0, beta=beta, scale=3.0, reference_wavelength=reference
            )()
            .single()
            .values
        )
        expected = 3.0 * (reference / GRID) ** beta * planck_jy(GRID, 40.0)
        assert emitted == pytest.approx(expected, rel=1e-12)

    def test_beta_zero_recovers_the_blackbody(self) -> None:
        """A greybody with no emissivity slope is a blackbody — a real degeneracy check."""
        grey = ModifiedBlackBody(GRID, temperature=40.0, beta=0.0, scale=1.0)().single().values
        black = BlackBody(GRID, temperature=40.0, scale=1.0)().single().values
        assert grey == pytest.approx(black, rel=1e-12)

    def test_a_power_law_is_normalised_at_its_reference(self) -> None:
        model = PowerLaw(GRID, norm=7.0, index=-2.0, reference_wavelength=2.0)
        emitted = model().single().values
        assert emitted == pytest.approx(7.0 * (GRID / 2.0) ** -2.0, rel=1e-12)
        assert float(np.interp(2.0, GRID, emitted)) == pytest.approx(7.0, rel=1e-3)

    def test_a_bad_reference_wavelength_is_refused(self) -> None:
        with pytest.raises(ValueError, match="finite and positive"):
            PowerLaw(GRID, reference_wavelength=0.0)

    def test_an_empty_grid_is_refused(self) -> None:
        with pytest.raises(ValueError, match="non-empty"):
            BlackBody(np.array([]))

    def test_several_channels_share_one_declaration(self) -> None:
        model = BlackBody(GRID, temperature=300.0, scale=1.0, channels=("blue", "red"))
        result = model()
        assert set(result.keys()) == {"blue", "red"}
        assert result["blue"].values == pytest.approx(result["red"].values)
        assert model.parameters.free_size == 0

    def test_repeated_channel_names_are_refused(self) -> None:
        with pytest.raises(ValueError, match="distinct"):
            BlackBody(GRID, channels=("sed", "sed"))

    def test_a_compiled_model_reuses_its_axes_by_identity(self) -> None:
        """``transformations.md`` §14's hot-loop contract."""
        model = BlackBody(GRID, temperature=300.0, scale=1.0, channels="sed")
        instrument = Instrument([Resample(np.linspace(2.0, 20.0, 9))], channel="sed")
        compiled = model.compile_for(negotiate([instrument]))
        assert compiled().single().axes is compiled().single().axes

    def test_a_quantity_grid_is_converted_once(self) -> None:
        """``DEVELOPMENT_PLAN.md`` §7's units trap: never per evaluation."""
        model = BlackBody(GRID * 1e4 * u.AA, temperature=300.0, scale=1.0)
        assert model.buffers["wavelength"].value == pytest.approx(GRID)


class TestParameterCoercion:
    """``as_parameter``: a prior fits, a number fixes, a Parameter passes through."""

    def test_a_number_is_held_fixed(self) -> None:
        assert BlackBody(GRID, temperature=300.0).parameters["temperature"].is_fixed

    def test_no_bijection_is_forced_onto_a_fitted_parameter(self) -> None:
        """The default machinery reads it off the prior; a class guess would not.

        Forcing one here was a real bug: ``norm ~ LogUniform(0.1, 10)`` is
        bounded *both* ways and lowers to ``Logit``, so an ``Identity`` (or a
        ``Log``) imposed by the constructor would silently contradict the
        declaration the conformance battery fixes.
        """
        from ampere.core import default_bijection_for

        norm = PowerLaw(GRID, norm=st.loguniform(0.1, 10.0)).parameters["norm"]
        assert norm.bijection is None
        assert type(default_bijection_for(norm.prior)).__name__ == "Logit"
        scale = BlackBody(GRID, scale=st.lognorm(0.2)).parameters["scale"]
        assert scale.bijection is None
        assert type(default_bijection_for(scale.prior)).__name__ == "Log"

    def test_a_ready_made_parameter_is_used_as_given(self) -> None:
        declared = Parameter("temperature", st.norm(300.0, 10.0), unit=u.K)
        assert BlackBody(GRID, temperature=declared).parameters["temperature"] is declared

    def test_a_parameter_under_the_wrong_name_is_refused(self) -> None:
        with pytest.raises(ParameterError, match="temperature"):
            BlackBody(GRID, temperature=Parameter("temp", st.norm(300.0, 10.0)))

    def test_nonsense_is_refused_by_name(self) -> None:
        with pytest.raises(ParameterError, match="must be a frozen"):
            BlackBody(GRID, temperature="warm")


# ---------------------------------------------------------------------------
# Instrument steps
# ---------------------------------------------------------------------------


def flat(grid: np.ndarray, value: float = 1.0, mask: np.ndarray | None = None) -> Spectrum:
    return Spectrum(grid * u.micron, np.full(grid.size, value) * u.Jy, mask=mask)


class TestBinEdges:
    def test_edges_bracket_the_centres(self) -> None:
        edges = bin_edges(np.array([1.0, 2.0, 4.0]))
        assert edges.tolist() == [0.5, 1.5, 3.0, 5.0]

    def test_a_single_sample_gets_a_unit_bin(self) -> None:
        assert bin_edges(np.array([3.0])).tolist() == [2.5, 3.5]


class TestCalibrationScale:
    def test_it_multiplies_and_keeps_the_axes(self) -> None:
        source = flat(GRID, 2.0)
        scaled = CalibrationScale(3.0)(source, None)
        assert scaled.values == pytest.approx(6.0)
        assert scaled.axes is source.axes

    def test_it_publishes_no_requirements(self) -> None:
        assert CalibrationScale(1.0).requirements() == ()

    def test_it_refuses_a_kind_it_does_not_accept(self) -> None:
        from ampere.core import PhotometricPoints

        points = PhotometricPoints(("A",), [1.0] * u.micron, [1.0] * u.Jy)
        with pytest.raises(CompositionError, match="accepts"):
            CalibrationScale(1.0)(points, None)


class TestResample:
    def test_a_constant_survives_resampling(self) -> None:
        """A flux *density* mean, so a flat spectrum stays at its own value."""
        target = np.linspace(2.0, 20.0, 12)
        assert Resample(target)(flat(GRID), None).values == pytest.approx(1.0)

    def test_the_weights_are_a_partition_of_unity(self) -> None:
        weights = Resample(np.linspace(2.0, 20.0, 12)).influence(GRID)
        assert weights.sum(axis=1) == pytest.approx(1.0)
        assert np.all(weights >= 0.0)

    def test_it_lands_on_exactly_the_target_grid(self) -> None:
        """``check_alignment`` compares axes for equality, so this must be exact."""
        target = np.linspace(2.0, 20.0, 12)
        assert Resample(target)(flat(GRID), None).spectral_axis.values.tolist() == target.tolist()

    def test_a_linear_ramp_is_reproduced(self) -> None:
        """A bin mean of a straight line is the line at the bin centre."""
        fine = np.linspace(1.0, 21.0, 4001)
        ramp = Spectrum(fine * u.micron, (3.0 + 2.0 * fine) * u.Jy)
        target = np.linspace(4.0, 18.0, 15)
        assert Resample(target)(ramp, None).values == pytest.approx(3.0 + 2.0 * target, rel=1e-6)

    def test_the_mask_follows_the_any_rule(self) -> None:
        target = np.linspace(2.0, 20.0, 12)
        mask = np.zeros(GRID.size, dtype=bool)
        # Mask samples that actually reach an output bin: GRID starts at 1 um,
        # well below the target range, so masking its first few would propagate
        # nowhere and the row would pass vacuously.
        mask[np.argmin(np.abs(GRID - 10.0))] = True
        rebinned = Resample(target)(flat(GRID, 1.0, mask), None)
        assert rebinned.mask is not None
        assert rebinned.is_masked
        assert rebinned.mask[np.argmin(np.abs(target - 10.0))]

    def test_it_publishes_coverage_finer_than_its_finest_bin(self) -> None:
        target = np.array([2.0, 4.0, 8.0])
        (requirement,) = Resample(target).requirements()
        assert requirement.axis == "spectral_axis"
        coordinates = requirement.coordinates().to_value(u.micron)
        assert coordinates.min() == pytest.approx(2.0)
        assert coordinates.max() == pytest.approx(8.0)
        assert np.diff(coordinates).max() <= 1.0 + 1e-12

    def test_an_unsorted_target_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="strictly increasing"):
            Resample(np.array([3.0, 1.0, 2.0]))


class TestLSFConvolution:
    def test_exactly_one_width_specification_is_required(self) -> None:
        with pytest.raises(TransformationError, match="neither"):
            LSFConvolution()
        with pytest.raises(TransformationError, match="both"):
            LSFConvolution(resolving_power=100.0, fwhm=0.1)

    def test_the_kernel_conserves_the_integral(self) -> None:
        grid = np.linspace(5.0, 15.0, 4001)
        spike = np.zeros(grid.size)
        spike[2000] = 1.0
        smoothed = LSFConvolution(fwhm=0.2)(Spectrum(grid * u.micron, spike * u.Jy), None)
        assert float(np.trapezoid(smoothed.values, grid)) == pytest.approx(
            float(np.trapezoid(spike, grid)), rel=1e-6
        )

    def test_the_smoothed_line_has_the_requested_width(self) -> None:
        grid = np.linspace(5.0, 15.0, 8001)
        spike = np.zeros(grid.size)
        spike[4000] = 1.0
        smoothed = LSFConvolution(fwhm=0.4)(Spectrum(grid * u.micron, spike * u.Jy), None)
        above = grid[smoothed.values >= smoothed.values.max() / 2.0]
        assert float(above[-1] - above[0]) == pytest.approx(0.4, rel=0.02)

    def test_a_constant_spectrum_is_unchanged_away_from_the_edges(self) -> None:
        grid = np.linspace(5.0, 15.0, 2001)
        smoothed = LSFConvolution(fwhm=0.1)(flat(grid), None)
        assert smoothed.values[100:-100] == pytest.approx(1.0, rel=1e-9)

    def test_constant_resolving_power_widens_with_wavelength(self) -> None:
        step = LSFConvolution(resolving_power=100.0)
        assert step.sigma(np.array([20.0]))[0] == pytest.approx(
            10.0 * step.sigma(np.array([2.0]))[0]
        )

    def test_it_publishes_nothing_until_it_knows_the_downstream_range(self) -> None:
        """An LSF alone genuinely cannot know what range to ask for."""
        assert LSFConvolution(resolving_power=100.0).requirements() == ()

    def test_configure_from_learns_the_range_and_pads_it(self) -> None:
        """``transformations.md`` §5's chain-internal negotiation, in use.

        The padding is what stops the outermost output samples being convolved
        against an edge instead of against data.
        """
        step = LSFConvolution(resolving_power=100.0)
        Instrument([step, Resample(np.linspace(2.0, 20.0, 12))], channel="sed")
        (requirement,) = step.requirements()
        low, high, _, power = requirement.segments()[0]
        assert low < 2.0
        assert high > 20.0
        assert power is not None and power > 100.0

    def test_the_padding_is_several_kernel_widths(self) -> None:
        step = LSFConvolution(fwhm=0.5)
        Instrument([step, Resample(np.linspace(10.0, 20.0, 5))], channel="sed")
        low, high, max_step, _ = step.requirements()[0].segments()[0]
        sigma = float(step.sigma(np.array([10.0]))[0])
        assert 10.0 - low == pytest.approx(5.0 * sigma)
        assert high - 20.0 == pytest.approx(5.0 * sigma)
        assert max_step is not None and max_step < 0.5


class TestSyntheticPhotometry:
    """The real step, and Gap 1's scenario."""

    def tophats(self) -> SyntheticPhotometry:
        tabulation = np.linspace(1.0, 5.0, 41)
        response = np.zeros((2, tabulation.size))
        response[0, (tabulation >= 1.5) & (tabulation <= 2.5)] = 1.0
        response[1, (tabulation >= 3.5) & (tabulation <= 4.5)] = 1.0
        return SyntheticPhotometry(["A", "B"], tabulation, response)

    def test_a_flat_spectrum_measures_its_own_value(self) -> None:
        step = self.tophats()
        measured = step(flat(step.tabulation(), 4.0), None)
        assert measured.values == pytest.approx(4.0)
        assert tuple(measured.filters) == ("A", "B")

    def test_it_changes_the_container_kind(self) -> None:
        from ampere.core import PhotometricPoints

        step = self.tophats()
        assert isinstance(step(flat(step.tabulation()), None), PhotometricPoints)

    def test_the_pivot_lies_inside_the_response(self) -> None:
        pivots = self.tophats().pivots()
        assert 1.5 < pivots[0] < 2.5
        assert 3.5 < pivots[1] < 4.5

    def test_it_publishes_points_at_its_own_tabulation(self) -> None:
        """A density requirement would only guarantee coverage — §3's lesson."""
        step = self.tophats()
        (requirement,) = step.requirements()
        assert requirement.points is not None
        assert requirement.points == pytest.approx(step.tabulation())

    def test_a_second_instrument_on_the_channel_does_not_break_it(self) -> None:
        """``spectrum_photometry.md`` Gap 1, end to end.

        The union hands the step a longer, reordered grid. Read positionally,
        the response matrix would multiply the wrong fluxes — or, with luck,
        raise a bare ``matmul`` shape error naming neither channel nor
        negotiation. ``Axis.locate`` is what makes it come out right.
        """
        step = self.tophats()
        requirements = negotiate(
            [
                Instrument([step], channel="sed", label="phot"),
                Instrument(
                    [Resample(np.array([1.234, 2.345, 3.456]))], channel="sed", label="spec"
                ),
            ]
        )
        grid = requirements["sed"]["spectral_axis"].coordinates().to_value(u.micron)
        assert grid.size > step.tabulation().size
        assert step(flat(grid, 4.0), None).values == pytest.approx(4.0)

    def test_a_grid_missing_its_tabulation_is_refused_loudly(self) -> None:
        step = self.tophats()
        with pytest.raises(SchemaError, match="COORDINATE_RTOL"):
            step(flat(np.linspace(1.0, 5.0, 37)), None)

    def test_the_mask_propagates_through_the_located_columns(self) -> None:
        step = self.tophats()
        tabulation = step.tabulation()
        mask = np.zeros(tabulation.size, dtype=bool)
        mask[(tabulation >= 1.8) & (tabulation <= 1.9)] = True
        measured = step(flat(tabulation, 1.0, mask), None)
        assert measured.mask is not None
        assert measured.mask.tolist() == [True, False]

    def test_repeated_filter_names_are_refused(self) -> None:
        with pytest.raises(TransformationError, match="unique"):
            SyntheticPhotometry(["A", "A"], np.linspace(1.0, 2.0, 5), np.ones((2, 5)))

    def test_a_mis_shaped_response_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="response array"):
            SyntheticPhotometry(["A"], np.linspace(1.0, 2.0, 5), np.ones((1, 4)))

    def test_a_filter_with_no_response_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="integrates to zero"):
            SyntheticPhotometry(["A"], np.linspace(1.0, 2.0, 5), np.zeros((1, 5)))

    def test_a_negative_response_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="non-negative"):
            SyntheticPhotometry(["A"], np.linspace(1.0, 2.0, 5), -np.ones((1, 5)))


class TestSyntheticPhotometryFromLibrary:
    """The pyphot route. pyphot is a base dependency, so this always runs."""

    def test_it_loads_ampere_s_own_filter_set(self) -> None:
        grid = np.geomspace(0.5, 30.0, 2000)
        step = SyntheticPhotometry.from_library(["2MASS_J", "2MASS_H", "2MASS_Ks"], grid)
        assert step.filters == ("2MASS_J", "2MASS_H", "2MASS_Ks")

    def test_the_pivots_are_the_published_ones(self) -> None:
        """Known 2MASS pivot wavelengths, to the precision they are quoted at."""
        grid = np.geomspace(0.5, 30.0, 4000)
        step = SyntheticPhotometry.from_library(["2MASS_J", "2MASS_H", "2MASS_Ks"], grid)
        assert step.pivots() == pytest.approx([1.239, 1.649, 2.164], abs=2e-3)

    def test_a_flat_spectrum_measures_its_own_value(self) -> None:
        grid = np.geomspace(0.5, 30.0, 2000)
        step = SyntheticPhotometry.from_library(["2MASS_J", "WISE_RSR_W1"], grid)
        assert step(flat(grid, 2.5), None).values == pytest.approx(2.5, rel=1e-6)

    def test_it_composes_onto_a_model(self) -> None:
        grid = np.geomspace(0.5, 30.0, 2000)
        step = SyntheticPhotometry.from_library(["2MASS_J", "WISE_RSR_W1"], grid)
        instrument = Instrument([step], channel="sed", label="phot")
        model = BlackBody(grid, temperature=5000.0, scale=1e-14, channels="sed")
        measured = instrument(model())
        assert measured.values.shape == (2,)
        assert np.all(measured.values > 0.0)


# ---------------------------------------------------------------------------
# Noise
# ---------------------------------------------------------------------------


class TestFractionalModelNoise:
    """Beyond the X-1 conformance rows: the edges and the refusals."""

    def spectra(self) -> tuple[Spectrum, Spectrum]:
        grid = np.linspace(1.0, 5.0, 9)
        mean = 1.0 + 0.3 * np.sin(grid)
        predicted = Spectrum(grid * u.micron, mean * u.Jy)
        observed = Spectrum(
            grid * u.micron, (mean + 0.05) * u.Jy, uncertainty=np.full(grid.size, 0.1) * u.Jy
        )
        return predicted, observed

    def test_f_is_an_ordinary_fitted_parameter(self) -> None:
        """W2.1's requirement, and the difference from the spec's fixed-``f`` sketch."""
        noise = FractionalModelNoise(st.loguniform(0.01, 1.0))
        assert noise.parameters.free_names == ("f",)
        assert Likelihood(GaussianFamily(), noise).parameters.free_names == ("f",)

    def test_a_number_holds_f_fixed(self) -> None:
        assert FractionalModelNoise(0.1).parameters["f"].is_fixed

    def test_the_data_scale_multiplies_only_the_data_term(self) -> None:
        """``sigma_eff**2 = (s*sigma_data)**2 + (f*mu)**2`` — ``s`` is inside the first."""
        predicted, observed = self.spectra()
        retain = np.ones(observed.n_samples, dtype=bool)
        mean = np.asarray(predicted.values, dtype=float)
        noise = FractionalModelNoise(0.4, scale=2.0)
        sigma = noise.sigma(observed, retain, {}, predicted=mean)
        assert sigma == pytest.approx(np.sqrt((2.0 * 0.1) ** 2 + (0.4 * mean) ** 2))

    def test_without_data_uncertainties_the_model_error_is_the_whole_noise(self) -> None:
        """A degenerate but well-defined case: ``sigma_eff = f * |mu|``."""
        grid = np.linspace(1.0, 5.0, 9)
        mean = 1.0 + 0.3 * np.sin(grid)
        bare = Spectrum(grid * u.micron, mean * u.Jy)
        sigma = FractionalModelNoise(0.2).sigma(
            bare, np.ones(grid.size, dtype=bool), {}, predicted=mean
        )
        assert sigma == pytest.approx(0.2 * mean)

    def test_it_refuses_to_guess_when_given_no_prediction(self) -> None:
        _, observed = self.spectra()
        with pytest.raises(LikelihoodError, match="needs the model prediction"):
            FractionalModelNoise(0.1).sigma(observed, np.ones(observed.n_samples, dtype=bool), {})

    def test_a_negative_fraction_is_refused(self) -> None:
        predicted, observed = self.spectra()
        with pytest.raises(LikelihoodError, match="non-negative"):
            FractionalModelNoise(-0.1).sigma(
                observed,
                np.ones(observed.n_samples, dtype=bool),
                {},
                predicted=np.asarray(predicted.values),
            )

    def test_f_zero_reduces_to_independent_noise(self) -> None:
        predicted, observed = self.spectra()
        from ampere.core import IndependentNoise

        plain = Likelihood(GaussianFamily(), IndependentNoise())
        fractional = Likelihood(GaussianFamily(), FractionalModelNoise(0.0))
        assert fractional.log_prob(predicted, observed) == pytest.approx(
            plain.log_prob(predicted, observed)
        )

    def test_the_pointwise_terms_sum_to_the_joint(self) -> None:
        """``pointwise_log_prob`` slices ``predicted`` per sample; the lengths must line up."""
        predicted, observed = self.spectra()
        likelihood = Likelihood(GaussianFamily(), FractionalModelNoise(0.3))
        terms = likelihood.pointwise_log_prob(predicted, observed)
        assert float(terms.sum()) == pytest.approx(likelihood.log_prob(predicted, observed))

    def test_a_masked_sample_is_excised_from_both_arrays(self) -> None:
        predicted, observed = self.spectra()
        mask = np.zeros(observed.n_samples, dtype=bool)
        mask[2] = True
        masked = Spectrum(
            observed.spectral_axis.quantity(),
            observed.values * u.Jy,
            uncertainty=observed.uncertainty * u.Jy,
            mask=mask,
        )
        likelihood = Likelihood(GaussianFamily(), FractionalModelNoise(0.3))
        kept = np.ones(observed.n_samples, dtype=bool)
        kept[2] = False
        trimmed_predicted = Spectrum(
            predicted.spectral_axis.values[kept] * u.micron, predicted.values[kept] * u.Jy
        )
        trimmed_observed = Spectrum(
            observed.spectral_axis.values[kept] * u.micron,
            observed.values[kept] * u.Jy,
            uncertainty=observed.uncertainty[kept] * u.Jy,
        )
        assert likelihood.log_prob(predicted, masked) == pytest.approx(
            likelihood.log_prob(trimmed_predicted, trimmed_observed)
        )


class TestFractionalModelGPNoise:
    def test_it_declares_the_kernel_hyperparameters_and_f(self) -> None:
        noise = FractionalModelGPNoise(Matern32(0.4, 2.0), f=st.loguniform(0.01, 1.0))
        assert "f" in noise.parameters

    def test_it_stays_conditionable(self) -> None:
        """``Likelihood.conditional`` refuses anything that is not a GP noise."""
        grid = np.linspace(1.0, 5.0, 9)
        mean = 1.0 + 0.3 * np.sin(grid)
        predicted = Spectrum(grid * u.micron, mean * u.Jy)
        observed = Spectrum(
            grid * u.micron, (mean + 0.05) * u.Jy, uncertainty=np.full(grid.size, 0.1) * u.Jy
        )
        likelihood = Likelihood(GaussianFamily(), FractionalModelGPNoise(Matern32(0.4, 2.0), f=0.2))
        conditional = likelihood.conditional(predicted, observed)
        assert conditional.mean.shape == (grid.size,)
