"""Conformance rows for the ``ModelResult`` schema.

W1.10's own statement of scope names "schema validation" alongside the prior
round trips and the log-prob agreement, and ``results_schema.md`` §16 hands
down the inventory: "``with_values`` preserves axes identity and unit;
mask/weights round trip; unit conversion is exact; ordering validation fires
per kind; complex values are confined to ``ALLOW_COMPLEX`` kinds; a default
channel is created for a bare container; kind mismatch raises ``ChannelError``."

These look like container unit tests, and ``tests/core/test_results_schema.py``
already holds them as such. They are here as well, and parametrised over the
backend, because ``lowering.md`` makes the containers the boundary a backend
lowers *to*: the schema is the two-conventions bug class's cure
(``DEVELOPMENT_PLAN.md`` §7), and a backend that produced containers violating
it would reintroduce exactly that. So each row runs against a container the
*backend's own model* produced, not only against one the suite built.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    DEFAULT_CHANNEL,
    ChannelError,
    Cube,
    Image,
    ModelResult,
    PhotometricPoints,
    SchemaError,
    Spectrum,
    TimeSeries,
    VisibilitySet,
)

from .composition import GP_GRID
from .protocol import ConformanceBackend, ModelSpec, Tolerances


@pytest.fixture
def produced(backend: ConformanceBackend) -> Spectrum:
    """A container the backend's own model produced, at its prior median."""
    model = backend.model(ModelSpec(coordinates=GP_GRID))
    centre = np.full(model.parameters.free_size, 0.5)
    return model(model.parameters.prior_transform(centre)).single()


class TestWhatABackendProduces:
    """The schema a backend's own output must satisfy."""

    def test_the_container_is_well_formed(self, produced: Spectrum) -> None:
        assert isinstance(produced, Spectrum)
        assert produced.spectral_axis.unit == u.micron
        assert produced.unit == u.Jy
        assert produced.shape == (len(GP_GRID),)
        assert produced.spectral_axis.values == pytest.approx(GP_GRID)

    def test_a_bare_container_gets_the_default_channel(self, produced: Spectrum) -> None:
        result = ModelResult(produced)
        assert set(result) == {DEFAULT_CHANNEL}
        assert result.is_single
        assert result.single() is produced

    def test_the_result_carries_the_parameters_it_was_evaluated_at(
        self, backend: ConformanceBackend
    ) -> None:
        """The W1.4 ruling: ``Model.__call__`` attaches θ to what it returns."""
        model = backend.model(ModelSpec(coordinates=GP_GRID))
        centre = np.full(model.parameters.free_size, 0.5)
        result = model(model.parameters.prior_transform(centre))
        assert result.parameters is not None
        assert set(result.parameters) == set(model.parameters.names)

    def test_requiring_the_wrong_kind_is_a_channel_error(self, produced: Spectrum) -> None:
        result = ModelResult({"sed": produced})
        assert result.require("sed", Spectrum) is produced
        with pytest.raises(ChannelError, match="required"):
            result.require("sed", VisibilitySet)
        with pytest.raises(ChannelError, match="no channel named"):
            result.require("missing")


class TestWithValues:
    """The hot-loop refill: axes by identity, unit unchanged."""

    def test_it_reuses_the_axes_object(self, produced: Spectrum) -> None:
        refilled = produced.with_values(produced.values * 2.0)
        assert refilled.axes is produced.axes
        assert refilled.unit == produced.unit
        assert refilled.values == pytest.approx(produced.values * 2.0)

    def test_it_refuses_to_change_the_value_unit(self, produced: Spectrum) -> None:
        with pytest.raises(SchemaError, match="not convertible"):
            produced.with_values(produced.values * u.K)

    def test_it_converts_a_quantity_in_a_convertible_unit(
        self, produced: Spectrum, tolerances: Tolerances
    ) -> None:
        refilled = produced.with_values(produced.values * 1000.0 * u.mJy)
        assert refilled.unit == u.Jy
        assert refilled.values == pytest.approx(produced.values, abs=tolerances.analytic)

    def test_it_inherits_the_mask_unless_told_otherwise(self, produced: Spectrum) -> None:
        mask = np.zeros(produced.n_samples, dtype=bool)
        mask[2] = True
        masked = produced.with_values(produced.values, mask=mask)
        assert masked.with_values(produced.values).mask.tolist() == mask.tolist()


class TestMaskAndWeights:
    """A mask is a declaration about samples; weights are its arithmetic form."""

    def test_the_mask_and_its_weights_agree(self, produced: Spectrum) -> None:
        mask = np.zeros(produced.n_samples, dtype=bool)
        mask[[1, 4]] = True
        masked = produced.with_values(produced.values, mask=mask)

        assert masked.is_masked
        assert masked.n_samples == produced.n_samples
        assert masked.n_valid == produced.n_samples - 2
        assert masked.valid.tolist() == (~mask).tolist()
        assert masked.weights().tolist() == (~mask).astype(float).tolist()

    def test_an_unmasked_container_declares_no_mask(self, produced: Spectrum) -> None:
        assert produced.mask is None
        assert not produced.is_masked
        assert produced.weights().tolist() == [1.0] * produced.n_samples

    def test_a_non_boolean_mask_is_refused(self, produced: Spectrum) -> None:
        with pytest.raises(SchemaError):
            produced.with_values(produced.values, mask=np.arange(produced.n_samples))


class TestUnitConversion:
    """Conversion happens once, at composition, and it is exact."""

    def test_it_scales_values_and_uncertainties_together(self, tolerances: Tolerances) -> None:
        grid = np.asarray(GP_GRID[:4], dtype=float)
        spectrum = Spectrum(
            grid * u.micron,
            np.array([1.0, 2.0, 3.0, 4.0]) * u.Jy,
            uncertainty=np.full(4, 0.1) * u.Jy,
        )
        converted = spectrum.to_unit(u.mJy)
        assert converted.unit == u.mJy
        assert converted.values == pytest.approx(spectrum.values * 1000.0, abs=tolerances.analytic)
        assert converted.uncertainty == pytest.approx(
            np.asarray(spectrum.uncertainty) * 1000.0, abs=tolerances.analytic
        )

    def test_a_round_trip_returns_the_original(self, tolerances: Tolerances) -> None:
        grid = np.asarray(GP_GRID[:4], dtype=float)
        spectrum = Spectrum(grid * u.micron, np.array([1.0, 2.0, 3.0, 4.0]) * u.Jy)
        assert spectrum.to_unit(u.mJy).to_unit(u.Jy).values == pytest.approx(
            spectrum.values, abs=tolerances.analytic
        )

    def test_an_inconvertible_unit_is_refused(self) -> None:
        grid = np.asarray(GP_GRID[:4], dtype=float)
        spectrum = Spectrum(grid * u.micron, np.array([1.0, 2.0, 3.0, 4.0]) * u.Jy)
        with pytest.raises(SchemaError):
            spectrum.to_unit(u.K)


class TestOrderingAndComplexity:
    """Per-kind validation: what each container will and will not accept."""

    def test_a_spectrum_requires_a_strictly_increasing_axis(self) -> None:
        with pytest.raises(SchemaError):
            Spectrum(np.array([2.0, 1.0, 3.0]) * u.micron, np.ones(3) * u.Jy)

    def test_a_time_series_requires_one_too(self) -> None:
        with pytest.raises(SchemaError):
            TimeSeries(np.array([1.0, 0.0]) * u.day, np.ones(2) * u.mag)

    def test_photometry_accepts_any_order(self) -> None:
        points = PhotometricPoints(
            ("W2", "W1"), np.array([4.6, 3.4]) * u.micron, np.array([2.0, 1.0]) * u.Jy
        )
        assert points.spectral_axis.values.tolist() == [4.6, 3.4]

    def test_an_image_accepts_a_monotonically_decreasing_axis(self) -> None:
        image = Image(
            np.array([1.0, 0.0]) * u.arcsec,
            np.array([0.0, 1.0, 2.0]) * u.arcsec,
            np.arange(6.0).reshape(2, 3),
        )
        assert image.shape == (2, 3)

    @pytest.mark.parametrize(
        "build",
        [
            lambda values: Spectrum(np.array([1.0, 2.0]) * u.micron, values),
            lambda values: TimeSeries(np.array([0.0, 1.0]) * u.day, values),
            lambda values: PhotometricPoints(("A", "B"), np.array([1.0, 2.0]) * u.micron, values),
        ],
        ids=["spectrum", "timeseries", "photometry"],
    )
    def test_complex_values_are_refused_by_the_real_kinds(self, build: object) -> None:
        with pytest.raises(SchemaError, match="complex"):
            build(np.array([1 + 2j, 3 - 1j]))  # type: ignore[operator]

    def test_only_visibilities_allow_complex(self) -> None:
        assert VisibilitySet.ALLOW_COMPLEX
        for kind in (Spectrum, PhotometricPoints, TimeSeries, Image, Cube):
            assert not kind.ALLOW_COMPLEX
        visibilities = VisibilitySet(
            np.array([1.0, 2.0]),
            np.array([3.0, 4.0]),
            np.array([2.2, 2.2]) * u.um,
            np.array([1 + 2j, 3 - 1j]),
        )
        assert visibilities.values.dtype.kind == "c"
