"""``ampere.backends.reference.interferometry``: the steps, the models, and a fit.

The conformance battery (``tests/conformance/test_interferometry.py``) holds
the arithmetic to closed forms, once per registered backend. This file is the
reference path's own: the construction-time refusals, the messages a user
actually meets, and the end-to-end fit — an emcee run on the synthetic binary
that has to recover the separation and flux ratio it was made from.

The fit runs through the **analytic** route (a model emitting visibilities
directly) rather than through the image-plus-transform one. That is not a
shortcut around the transform: the conformance rows hold the transform to the
same closed form this model evaluates, to a part in 10^10, so a fit through
either route is a fit to the same prediction — and a direct transform of a
100-by-100 image at every one of forty thousand likelihood evaluations is
twenty minutes of gate time to learn nothing the rows have not already said.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.backends.reference import (
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    ClosurePhase,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    TimeSmearing,
    UniformDisc,
    cell_solid_angle,
)
from ampere.backends.reference.interferometry import MAS_PER_RAD, _quadrature
from ampere.core import (
    ClosurePhases,
    ComplexGaussianFamily,
    CompositionError,
    Dataset,
    DatasetCollection,
    DatasetError,
    FittingProblem,
    Image,
    IndependentNoise,
    Instrument,
    Likelihood,
    Spectrum,
    TransformationError,
    VisibilitySet,
    VonMisesFamily,
    negotiate,
)
from ampere.inference import EmceeEngine

WAVELENGTH = 2.2
FIELD_OF_VIEW = 40.0
STATIONS = np.array([[0.0, 0.0], [-32.4, 39.8], [6.5, -30.1], [-14.7, -48.9]])
BINARY = {
    "separation": 12.0,
    "position_angle": 0.7,
    "flux_ratio": 0.42,
    "flux": 1.7,
    "component_fwhm": 2.0,
}


def coverage() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(u, v, lambda)`` for the array's six baselines."""
    metres = WAVELENGTH * 1e-6
    pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
    table = np.array([(STATIONS[j] - STATIONS[i]) / metres for i, j in pairs])
    return table[:, 0], table[:, 1], np.full(len(pairs), WAVELENGTH)


def triangles() -> tuple[np.ndarray, ...]:
    """``(u1, v1, u2, v2, lambda)`` in the canonical ordering."""
    metres = WAVELENGTH * 1e-6
    triples = [(i, j, k) for i in range(4) for j in range(i + 1, 4) for k in range(j + 1, 4)]
    first = np.array([(STATIONS[j] - STATIONS[i]) / metres for i, j, _ in triples])
    second = np.array([(STATIONS[k] - STATIONS[j]) / metres for _, j, k in triples])
    return first[:, 0], first[:, 1], second[:, 0], second[:, 1], np.full(len(triples), WAVELENGTH)


def visibilities(values: np.ndarray | None = None, **kwargs: Any) -> VisibilitySet:
    u_pts, v_pts, waves = coverage()
    filled = np.zeros(u_pts.size, dtype=np.complex128) if values is None else values
    return VisibilitySet(u_pts, v_pts, waves * u.micron, filled * u.Jy, **kwargs)


def closure_phases(values: np.ndarray | None = None, **kwargs: Any) -> ClosurePhases:
    u1, v1, u2, v2, waves = triangles()
    filled = np.zeros(u1.size) if values is None else values
    return ClosurePhases(u1, v1, u2, v2, waves * u.micron, filled * u.rad, **kwargs)


class TestQuadrature:
    """The node set every smearing step averages over."""

    def test_the_centre_node_is_first_and_exactly_zero(self) -> None:
        """The whole reason the layout works: sub-sample 0 is the sample itself.

        A smearing step recovers its output coordinates by slicing its input
        rather than recomputing them, and that is only bit-identical if the
        centre offset is exactly ``0.0``.
        """
        for count in (1, 3, 5, 7):
            offsets, weights = _quadrature(count)
            assert offsets.size == count
            assert offsets[0] == 0.0
            assert weights.sum() == pytest.approx(1.0, abs=1e-15)
            assert np.all(np.abs(offsets) <= 0.5)

    def test_an_even_node_count_is_refused_by_name(self) -> None:
        with pytest.raises(TransformationError, match="odd number of quadrature nodes"):
            _quadrature(4)


class TestCellSolidAngle:
    """The quadrature weight of the Fourier sum."""

    def test_a_uniform_grid_gives_the_pixel_area(self) -> None:
        grid = np.linspace(-10.0, 10.0, 21)
        cells = cell_solid_angle(grid, grid)
        assert cells.shape == (21, 21)
        assert np.allclose(cells, (1.0 / MAS_PER_RAD) ** 2)

    def test_a_decreasing_axis_still_gives_a_positive_area(self) -> None:
        """An ``Image`` axis is monotonic in *either* direction (sky axes run both ways)."""
        rising = np.linspace(-5.0, 5.0, 11)
        falling = rising[::-1]
        assert np.allclose(cell_solid_angle(rising, falling), cell_solid_angle(rising, rising))

    def test_an_uneven_grid_is_integrated_rather_than_assumed_regular(self) -> None:
        """Negotiation may return any grid satisfying a requirement, not a regular one."""
        grid = np.array([-4.0, -1.0, 0.0, 3.0])
        cells = cell_solid_angle(grid, np.array([0.0, 1.0]))
        assert cells.shape == (4, 2)
        # The total solid angle is the product of the two axes' full spans,
        # reflected outer cells included: 9 mas by 1 mas.
        assert cells.sum() == pytest.approx(9.0 / MAS_PER_RAD**2, rel=1e-12)


class TestFourierSampleConstruction:
    """The refusals a user meets before any arithmetic runs."""

    def test_from_observed_takes_a_visibility_set_unchanged(self) -> None:
        observed = visibilities()
        step = FourierSample.from_observed(observed, field_of_view=FIELD_OF_VIEW * u.mas)
        got_u, got_v, got_w = step.expanded_coverage
        assert np.array_equal(got_u, observed.u.values)
        assert np.array_equal(got_v, observed.v.values)
        assert np.array_equal(got_w, observed.spectral_axis.values)

    def test_from_observed_expands_a_triangle_into_three_baselines(self) -> None:
        """The canonical ordering, laid out ``3t``, ``3t+1``, ``3t+2``."""
        observed = closure_phases()
        step = FourierSample.from_observed(observed, field_of_view=FIELD_OF_VIEW * u.mas)
        got_u, got_v, got_w = step.expanded_coverage
        assert got_u.size == 3 * observed.n_samples
        assert np.array_equal(got_u[0::3], observed.u1.values)
        assert np.array_equal(got_u[1::3], observed.u2.values)
        assert np.allclose(got_u[2::3], -(observed.u1.values + observed.u2.values))
        assert np.allclose(got_v[2::3], -(observed.v1.values + observed.v2.values))
        assert np.array_equal(got_w[0::3], observed.spectral_axis.values)

    def test_another_kind_is_refused_by_name(self) -> None:
        spectrum = Spectrum([1.0, 2.0] * u.micron, [1.0, 1.0] * u.Jy)
        with pytest.raises(TransformationError, match="VisibilitySet or a ClosurePhases"):
            FourierSample.from_observed(spectrum, field_of_view=10.0 * u.mas)

    @pytest.mark.parametrize(
        ("kwargs", "message"),
        [
            ({"field_of_view": -1.0}, "finite and positive"),
            ({"field_of_view": 10.0, "oversampling": 0.5}, "at least 1"),
        ],
    )
    def test_bad_configuration_is_refused(self, kwargs: dict[str, Any], message: str) -> None:
        u_pts, v_pts, waves = coverage()
        with pytest.raises(TransformationError, match=message):
            FourierSample(u_pts, v_pts, waves, **kwargs)

    def test_mismatched_coverage_lengths_are_refused(self) -> None:
        with pytest.raises(TransformationError, match="one u, one v and one wavelength"):
            FourierSample([1.0, 2.0], [1.0], [2.2, 2.2], field_of_view=10.0)

    def test_a_changed_image_unit_is_refused_rather_than_silently_rescaled(self) -> None:
        """The units trap: a template built for one unit is not refilled from another."""
        step = FourierSample.from_observed(visibilities(), field_of_view=FIELD_OF_VIEW * u.mas)
        grid = np.linspace(-5.0, 5.0, 9) * u.mas
        step(Image(grid, grid, np.ones((9, 9)) * u.Jy / u.sr))
        with pytest.raises(TransformationError, match="compiled for an image in"):
            step(Image(grid, grid, np.ones((9, 9)) * u.mJy / u.sr))


class TestImageModels:
    """The trio that emits an ``Image``, and what they do with a negotiated grid."""

    def _requirements(self, **kwargs: Any) -> Any:
        instrument = Instrument(
            [
                FourierSample.from_observed(
                    visibilities(), field_of_view=FIELD_OF_VIEW * u.mas, **kwargs
                )
            ],
            channel="sky",
            label="array",
        )
        return negotiate([instrument])

    def test_the_negotiated_grid_is_adopted_and_shared_by_identity(self) -> None:
        model = GaussianSource.on_field(FIELD_OF_VIEW * u.mas, 4, fwhm=4.0, flux=1.4)
        compiled = model.compile_for(self._requirements())
        first = compiled.evaluate().single()
        second = compiled.evaluate().single()
        assert first.axes is second.axes
        assert first.unit == u.Jy / u.sr
        assert first.shape[0] > 4

    def test_a_gaussian_image_integrates_to_its_total_flux(self) -> None:
        """The brightness convention: Jy/sr, so the solid-angle sum is the flux."""
        model = GaussianSource.on_field(FIELD_OF_VIEW * u.mas, 401, fwhm=4.0, flux=1.4)
        image = model.evaluate().single()
        cells = cell_solid_angle(image.x.values, image.y.values)
        assert float(np.sum(image.values * cells)) == pytest.approx(1.4, rel=1e-6)

    def test_a_uniform_disc_image_integrates_to_its_total_flux(self) -> None:
        model = UniformDisc.on_field(FIELD_OF_VIEW * u.mas, 801, diameter=5.0, flux=1.1)
        image = model.evaluate().single()
        cells = cell_solid_angle(image.x.values, image.y.values)
        assert float(np.sum(image.values * cells)) == pytest.approx(1.1, rel=2e-3)

    def test_a_binary_puts_its_secondary_east_of_north(self) -> None:
        """The sky convention, asserted rather than described.

        At a position angle of zero the secondary sits at ``+y`` (north); at
        ``pi/2`` it sits at ``+x`` (east). Get this wrong and every closure
        phase is a statement about a different source.
        """
        for angle, axis in ((0.0, 1), (0.5 * math.pi, 0)):
            model = Binary.on_field(
                FIELD_OF_VIEW * u.mas, 201, **{**BINARY, "position_angle": angle}
            )
            image = model.evaluate().single()
            cells = cell_solid_angle(image.x.values, image.y.values)
            weight = np.asarray(image.values) * cells
            centroid = [
                float(np.sum(weight * image.x.values[:, None]) / np.sum(weight)),
                float(np.sum(weight * image.y.values[None, :]) / np.sum(weight)),
            ]
            expected = BINARY["separation"] * BINARY["flux_ratio"] / (1.0 + BINARY["flux_ratio"])
            assert centroid[axis] == pytest.approx(expected, rel=1e-3)
            assert centroid[1 - axis] == pytest.approx(0.0, abs=1e-6)

    def test_a_locked_grid_that_is_too_narrow_refuses_by_name(self) -> None:
        narrow = np.linspace(-2.0, 2.0, 401)
        model = GaussianSource(narrow, narrow, fwhm=4.0, flux=1.4, adopt_grid=False)
        with pytest.raises(CompositionError, match="still contributes to every visibility"):
            model.compile_for(self._requirements())

    def test_a_locked_grid_that_is_too_coarse_names_the_aliasing(self) -> None:
        coarse = np.linspace(-20.0, 20.0, 9)
        model = GaussianSource(coarse, coarse, fwhm=4.0, flux=1.4, adopt_grid=False)
        with pytest.raises(CompositionError, match="aliases"):
            model.compile_for(self._requirements())

    def test_a_locked_grid_that_satisfies_the_requirement_is_kept(self) -> None:
        fine = np.linspace(-25.0, 25.0, 1001)
        model = GaussianSource(fine, fine, fwhm=4.0, flux=1.4, adopt_grid=False)
        compiled = model.compile_for(self._requirements())
        assert np.array_equal(compiled.evaluate().single().x.values, fine)

    def test_a_degenerate_grid_is_refused_at_construction(self) -> None:
        with pytest.raises(ValueError, match="strictly increasing"):
            GaussianSource([1.0], [1.0], fwhm=2.0)


class TestClosurePhaseStep:
    """The three-to-one step's own refusals."""

    def test_a_sample_count_that_is_not_a_multiple_of_three_is_refused(self) -> None:
        source = visibilities(np.ones(6, dtype=np.complex128))
        # Six is a multiple of three, so slice to four.
        four = VisibilitySet(
            source.u.values[:4],
            source.v.values[:4],
            source.spectral_axis.values[:4] * u.micron,
            np.ones(4, dtype=np.complex128) * u.Jy,
        )
        with pytest.raises(TransformationError, match="multiple of three"):
            ClosurePhase()(four)

    def test_a_triangle_spanning_two_wavelengths_is_refused(self) -> None:
        u1, v1, u2, v2, _ = triangles()
        u_pts = np.stack([u1, u2, -(u1 + u2)], axis=1).reshape(-1)
        v_pts = np.stack([v1, v2, -(v1 + v2)], axis=1).reshape(-1)
        waves = np.full(u_pts.size, WAVELENGTH)
        waves[1] = 1.6
        source = VisibilitySet(
            u_pts,
            v_pts,
            waves * u.micron,
            np.ones(u_pts.size, dtype=np.complex128) * u.Jy,
        )
        with pytest.raises(TransformationError, match="different wavelengths"):
            ClosurePhase()(source)

    def test_the_output_is_in_radians_and_wrapped(self) -> None:
        observed = closure_phases()
        instrument = Instrument(
            [
                FourierSample.from_observed(
                    observed, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=10.0
                ),
                ClosurePhase(),
            ],
            channel="sky",
            label="t3",
        )
        model = Binary.on_field(FIELD_OF_VIEW * u.mas, 4, channels="sky", **BINARY)
        got = instrument(model.compile_for(negotiate([instrument])).evaluate())
        assert isinstance(got, ClosurePhases)
        assert got.unit == u.rad
        assert np.all(np.abs(got.values) <= math.pi)
        assert np.array_equal(got.u1.values, observed.u1.values)


class TestSmearingConstruction:
    """The uv rates that cannot be derived, and the refusal that says so."""

    def test_time_smearing_refuses_a_container_with_no_rates(self) -> None:
        with pytest.raises(TransformationError, match="how fast the baseline was moving"):
            TimeSmearing.from_observed(visibilities(), integration=60.0 * u.s)

    def test_time_smearing_reads_the_rates_from_extra_coords(self) -> None:
        u_pts, v_pts, _ = coverage()
        observed = visibilities(extra_coords={"du_dt": -v_pts * 7e-5, "dv_dt": u_pts * 7e-5})
        step = TimeSmearing.from_observed(observed, integration=60.0 * u.s)
        assert np.array_equal(step.buffers["du_dt"].value, observed.extra_coords["du_dt"])

    def test_a_rate_count_that_does_not_match_the_coverage_is_refused(self) -> None:
        step = TimeSmearing(integration=60.0 * u.s, du_dt=[1.0, 2.0], dv_dt=[1.0, 2.0])
        fourier = FourierSample.from_observed(visibilities(), field_of_view=FIELD_OF_VIEW * u.mas)
        with pytest.raises(TransformationError, match="uv rate"):
            Instrument([fourier, step], channel="sky", label="array")

    def test_a_non_positive_resolving_power_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="finite and positive"):
            BandwidthSmearing(resolving_power=0.0)

    def test_the_smearing_step_reduces_what_the_fourier_step_expanded(self) -> None:
        fourier = FourierSample.from_observed(
            visibilities(), field_of_view=FIELD_OF_VIEW * u.mas, oversampling=4.0
        )
        smearing = BandwidthSmearing(resolving_power=5.0, nodes=3)
        instrument = Instrument([fourier, smearing], channel="sky", label="array")
        assert fourier.expanded_coverage[0].size == 6 * 3
        model = GaussianSource.on_field(FIELD_OF_VIEW * u.mas, 4, fwhm=4.0, flux=1.4)
        got = instrument(model.compile_for(negotiate([instrument])).evaluate())
        assert got.n_samples == 6

    def test_a_smearing_step_handed_the_wrong_block_size_says_so(self) -> None:
        """A step composed with one chain and evaluated on another's output."""
        smearing = BandwidthSmearing(resolving_power=5.0, nodes=3)
        source = visibilities(np.ones(6, dtype=np.complex128))
        four = VisibilitySet(
            source.u.values[:4],
            source.v.values[:4],
            source.spectral_axis.values[:4] * u.micron,
            np.ones(4, dtype=np.complex128) * u.Jy,
        )
        with pytest.raises(TransformationError, match="not a multiple"):
            smearing(four)


class TestAmplitudeStep:
    def test_the_modulus_preserves_the_coordinates_and_the_mask(self) -> None:
        source = visibilities(
            np.array([1.0 + 1.0j, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.complex128),
            mask=np.array([False, True, False, False, False, False]),
        )
        got = Amplitude()(source)
        assert got.values.dtype.kind == "f"
        assert got.values[0] == pytest.approx(math.sqrt(2.0))
        assert got.mask is not None and got.mask.tolist() == source.mask.tolist()
        assert got.axes is source.axes


class TestAnalyticModels:
    """The route with no Fourier step in it."""

    def test_a_point_like_source_has_unit_visibility_everywhere(self) -> None:
        u_pts, v_pts, waves = coverage()
        model = GaussianSourceVisibilities(u_pts, v_pts, waves * u.micron, fwhm=1.0e-6, flux=2.5)
        got = model.evaluate().single()
        assert np.allclose(got.values, 2.5 + 0j, atol=1e-9)

    def test_from_observed_takes_a_triangle_apart_the_same_way_the_step_does(self) -> None:
        observed = closure_phases()
        model = BinaryVisibilities.from_observed(observed, channels="vis", **BINARY)
        instrument = Instrument([ClosurePhase()], channel="vis", label="t3")
        got = instrument(model.evaluate())
        assert isinstance(got, ClosurePhases)
        assert np.array_equal(got.u1.values, observed.u1.values)

    def test_the_zero_baseline_visibility_is_the_total_flux(self) -> None:
        model = BinaryVisibilities([0.0], [0.0], [WAVELENGTH], **BINARY)
        assert model.evaluate().single().values[0] == pytest.approx(BINARY["flux"] + 0j)


class TestTheSyntheticBinaryIsRecovered:
    """The end-to-end acceptance: an emcee fit that finds the source it was made from."""

    #: Enough coverage for two parameters to be well determined: the array
    #: observed at four hour angles, so the uv plane is filled rather than
    #: sampled at six points.
    HOUR_ANGLES = np.linspace(-0.6, 0.6, 4)

    @classmethod
    def _coverage(cls) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        u_pts, v_pts, _ = coverage()
        rotated_u = np.concatenate(
            [u_pts * math.cos(h) - v_pts * math.sin(h) for h in cls.HOUR_ANGLES]
        )
        rotated_v = np.concatenate(
            [u_pts * math.sin(h) + v_pts * math.cos(h) for h in cls.HOUR_ANGLES]
        )
        return rotated_u, rotated_v, np.full(rotated_u.size, WAVELENGTH)

    @pytest.fixture(scope="class")
    @classmethod
    def run(cls) -> Any:
        u_pts, v_pts, waves = cls._coverage()
        truth = BinaryVisibilities(u_pts, v_pts, waves * u.micron, channels="vis", **BINARY)
        exact = truth.evaluate().single().values
        sigma = 0.03
        rng = np.random.default_rng(20260911)
        noisy = exact + sigma * (
            rng.standard_normal(exact.size) + 1j * rng.standard_normal(exact.size)
        )
        observed = VisibilitySet(
            u_pts,
            v_pts,
            waves * u.micron,
            noisy * u.Jy,
            uncertainty=np.full(exact.size, sigma) * u.Jy,
        )
        model = BinaryVisibilities(
            u_pts,
            v_pts,
            waves * u.micron,
            channels="vis",
            separation=st.uniform(4.0, 16.0),
            position_angle=BINARY["position_angle"],
            flux_ratio=st.uniform(0.05, 0.9),
            flux=BINARY["flux"],
            component_fwhm=BINARY["component_fwhm"],
        )
        problem = FittingProblem(
            model,
            DatasetCollection(
                {
                    "vis": Dataset(
                        observed,
                        Instrument([], channel="vis", input_kind=VisibilitySet, label="array"),
                        likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()),
                        label="vis",
                    )
                }
            ),
            seed=20260911,
        )
        return EmceeEngine(problem, walkers=16).run(steps=700, burn_in=300)

    @pytest.mark.parametrize(
        ("name", "injected"),
        [("model.separation", BINARY["separation"]), ("model.flux_ratio", BINARY["flux_ratio"])],
    )
    def test_the_injected_value_is_inside_the_central_95_per_cent(
        self, run: Any, name: str, injected: float
    ) -> None:
        draws = np.asarray(run["posterior"][name].values).reshape(-1)
        low, high = np.percentile(draws, [2.5, 97.5])
        assert low <= injected <= high
        # And the posterior is actually informative, not merely wide enough to
        # contain the truth by accident.
        assert (high - low) < 0.35 * abs(injected)


class TestTwoDatasetsShareOneSky:
    """The reference path's own copy of the composition, for the error messages."""

    def test_two_instruments_on_one_channel_need_distinct_labels(self) -> None:
        instrument = Instrument(
            [FourierSample.from_observed(visibilities(), field_of_view=FIELD_OF_VIEW * u.mas)],
            channel="sky",
        )
        other = Instrument(
            [
                FourierSample.from_observed(closure_phases(), field_of_view=FIELD_OF_VIEW * u.mas),
                ClosurePhase(),
            ],
            channel="sky",
        )
        model = Binary.on_field(FIELD_OF_VIEW * u.mas, 4, channels="sky", **BINARY)
        compiled = model.compile_for(negotiate([instrument, other]))
        result = compiled.evaluate()
        datasets = {
            "vis": Dataset(
                visibilities(instrument(result).values, uncertainty=np.full(6, 0.02) * u.Jy),
                instrument,
                likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()),
                label="vis",
            ),
            "t3": Dataset(
                closure_phases(other(result).values, uncertainty=np.full(4, 0.05) * u.rad),
                other,
                likelihood=Likelihood(VonMisesFamily(), IndependentNoise()),
                label="t3",
            ),
        }
        with pytest.raises(DatasetError, match="instrument label"):
            FittingProblem(compiled, DatasetCollection(datasets), seed=1)
