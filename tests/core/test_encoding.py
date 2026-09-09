"""The W3.3 coordinate--value--mask encoding, against ``encoding.md``.

Organised around the acceptance criteria of ``WORK_ITEMS.md`` W3.3, and around
the one property the whole contract rests on: **the layout is a function of the
problem, not of the data it is later shown**. Every statistic is computed once,
from the observed containers, and frozen with the mask into a hashable record;
the observation the posterior is conditioned on and the rows a network trains on
are therefore standardised identically by construction rather than by care.

Four things are checked here and nowhere else:

1. **The round trip**, on one container of every kind the schema ships,
   including a masked one, an irregularly sampled one and a complex one. A
   packing that cannot be inverted is a packing whose columns nobody can check.
2. **Padding and the mask column**, across two datasets of different lengths in
   one tensor -- the case that makes a *collection* encodable at all.
3. **The refusals**: a layout that does not describe a problem, named field by
   field; a row cap the observation exceeds; a differently sampled container.
4. **``unpack``**, and in particular ``per_dataset`` on a grid kind, whose slice
   must reshape back to the container's own axes exactly.

The wrappers that turn an unpacked view into a torch embedding are
``tests/inference/test_sbi.py``'s: they need the ``sbi`` extra, and nothing
here does.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    Cube,
    EncodingError,
    EncodingLayout,
    Image,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    VisibilitySet,
    decode,
    encode,
    encode_observations,
    unpack,
)
from ampere.core.encoding import DEFAULT_FOURIER_BANDS


class _Holder:
    """The smallest thing :func:`encode_observations` reads.

    The encoding takes ``observed`` and ``effective_mask`` off each dataset and
    nothing else, so a container fixture can be exercised without composing a
    whole :class:`~ampere.core.dataset.FittingProblem` around it -- which for a
    complex ``VisibilitySet`` would mean choosing a likelihood family this test
    has no opinion about.
    """

    def __init__(self, observed: Any) -> None:
        self.observed = observed
        self.effective_mask = None if observed.mask is None else np.asarray(observed.mask).ravel()


# ---------------------------------------------------------------------------
# One container of every kind the schema ships
# ---------------------------------------------------------------------------


def a_spectrum() -> Spectrum:
    """Regularly sampled, uncertainties, nothing masked."""
    grid = np.linspace(1.0, 10.0, 12)
    return Spectrum(
        grid * u.um,
        (2.0 * grid) * u.Jy,
        uncertainty=np.full(grid.size, 0.3) * u.Jy,
    )


def a_masked_spectrum() -> Spectrum:
    """Two samples excluded, one of them a NaN the mask exists to remove."""
    grid = np.linspace(1.0, 10.0, 12)
    values = 2.0 * grid
    values[3] = np.nan
    mask = np.zeros(grid.size, dtype=bool)
    mask[[3, 7]] = True
    return Spectrum(
        grid * u.um,
        values * u.Jy,
        uncertainty=np.full(grid.size, 0.3) * u.Jy,
        mask=mask,
    )


def an_irregular_spectrum() -> Spectrum:
    """Two narrow windows far apart: the §4.2 sampling no grid describes."""
    windows = np.concatenate([np.linspace(866.9, 867.1, 7), np.linspace(1300.3, 1300.5, 5)])
    return Spectrum(
        windows * u.um,
        np.linspace(1.0, 1.4, windows.size) * u.Jy,
        uncertainty=np.full(windows.size, 0.05) * u.Jy,
    )


def photometry() -> PhotometricPoints:
    """Arbitrary coordinate order, and a filter name the encoding does not carry."""
    return PhotometricPoints(
        ["J", "H", "K", "W1"],
        [1.25, 1.65, 2.2, 3.4] * u.um,
        [10.0, 8.0, 6.0, 3.0] * u.Jy,
        uncertainty=[0.5, 0.4, 0.3, 0.2] * u.Jy,
    )


def a_light_curve() -> TimeSeries:
    """Uneven cadence, no uncertainties: ``has_sigma`` is False and says so."""
    times = np.array([0.0, 0.5, 0.6, 2.0, 5.0, 5.1])
    return TimeSeries(times * u.day, (1.0 + 0.1 * times) * u.mag)


def an_image() -> Image:
    """A grid kind: ``per_dataset`` must reshape its rows back to (x, y)."""
    x = np.linspace(-1.0, 1.0, 4)
    y = np.linspace(-2.0, 2.0, 3)
    values = np.arange(12.0).reshape(4, 3) + 1.0
    return Image(
        x * u.arcsec,
        y * u.arcsec,
        values * u.Jy,
        uncertainty=np.full((4, 3), 0.1) * u.Jy,
    )


def a_cube() -> Cube:
    """Three axes, so the layout's coordinate width is set by this one."""
    x = np.linspace(-1.0, 1.0, 2)
    y = np.linspace(-1.0, 1.0, 3)
    spectral = np.linspace(1.0, 2.0, 4)
    values = np.arange(24.0).reshape(2, 3, 4) + 1.0
    return Cube(x * u.arcsec, y * u.arcsec, spectral * u.um, values * u.Jy)


def visibilities() -> VisibilitySet:
    """Complex values on a scattered point set: two value columns, and a mask."""
    u_coord = np.array([10.0, -30.0, 55.0, 5.0, -80.0])
    v_coord = np.array([-20.0, 40.0, 5.0, 65.0, 10.0])
    values = np.array([1.0 + 0.5j, 0.8 - 0.3j, 0.2 + 0.9j, -0.4 + 0.1j, 0.6 + 0.6j])
    return VisibilitySet(
        u_coord * u.dimensionless_unscaled,
        v_coord * u.dimensionless_unscaled,
        values * u.Jy,
        uncertainty=np.full(5, 0.05) * u.Jy,
        mask=np.array([False, False, True, False, False]),
    )


CONTAINERS = {
    "spectrum": a_spectrum,
    "masked_spectrum": a_masked_spectrum,
    "irregular_spectrum": an_irregular_spectrum,
    "photometry": photometry,
    "light_curve": a_light_curve,
    "image": an_image,
    "cube": a_cube,
    "visibilities": visibilities,
}


def one(name: str) -> dict[str, _Holder]:
    """A one-dataset collection around the named container fixture."""
    return {name: _Holder(CONTAINERS[name]())}


# ---------------------------------------------------------------------------
# 1. The round trip
# ---------------------------------------------------------------------------


class TestTheRoundTrip:
    """Encode, decode, and get the container back."""

    @pytest.mark.parametrize("name", sorted(CONTAINERS))
    def test_every_container_kind_survives_encode_and_decode(self, name: str) -> None:
        datasets = one(name)
        original = datasets[name].observed
        layout = EncodingLayout.from_datasets(datasets)
        back = decode(encode_observations(datasets, layout=layout))[name]

        assert back.values.shape == original.values.shape
        assert np.asarray(back.values) == pytest.approx(
            np.nan_to_num(np.asarray(original.values)), rel=1e-9, abs=1e-12
        )
        for axis in original.axes:
            assert back.coordinates[axis.name] == pytest.approx(
                np.asarray(axis.values), rel=1e-9, abs=1e-12
            )
        if original.uncertainty is None:
            assert back.uncertainty is None
        else:
            assert back.uncertainty == pytest.approx(np.asarray(original.uncertainty), rel=1e-9)
        if original.mask is None:
            assert back.mask is None
        else:
            assert np.array_equal(back.mask, np.asarray(original.mask))

    def test_a_decoded_dataset_rebuilds_a_container_from_its_template(self) -> None:
        """``as_container`` borrows the axes, the unit and the extra coordinates."""
        datasets = one("photometry")
        original = datasets["photometry"].observed
        layout = EncodingLayout.from_datasets(datasets)
        decoded = decode(encode_observations(datasets, layout=layout))["photometry"]
        rebuilt = decoded.as_container(original)
        assert type(rebuilt) is type(original)
        assert rebuilt.unit == original.unit
        assert list(rebuilt.extra_coords["filters"]) == list(original.extra_coords["filters"])
        assert np.asarray(rebuilt.values) == pytest.approx(np.asarray(original.values))

    def test_a_masked_sample_may_be_a_nan_and_the_encoding_stays_finite(self) -> None:
        """The one amendment W3.3 made to §3, checked rather than asserted."""
        datasets = one("masked_spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        assert np.all(np.isfinite(encoded.values))
        assert math.isnan(float(np.asarray(datasets["masked_spectrum"].observed.values)[3]))

    def test_a_non_finite_retained_sample_is_refused_by_name(self) -> None:
        """Masked, a NaN is fine; retained, it is a mistake the packing names."""
        grid = np.linspace(1.0, 5.0, 5)
        values = np.array([1.0, 2.0, np.nan, 4.0, 5.0])
        datasets = {"bad": _Holder(Spectrum(grid * u.um, values * u.Jy))}
        layout = EncodingLayout.from_datasets(datasets)
        with pytest.raises(EncodingError, match="non-finite"):
            encode_observations(datasets, layout=layout)

    def test_a_non_positive_uncertainty_on_a_retained_sample_is_refused(self) -> None:
        grid = np.linspace(1.0, 5.0, 5)
        datasets = {
            "bad": _Holder(
                Spectrum(
                    grid * u.um,
                    grid * u.Jy,
                    uncertainty=np.array([0.1, 0.1, 0.0, 0.1, 0.1]) * u.Jy,
                )
            )
        }
        layout = EncodingLayout.from_datasets(datasets)
        with pytest.raises(EncodingError, match="strictly positive"):
            encode_observations(datasets, layout=layout)


# ---------------------------------------------------------------------------
# 2. A collection: padding, the mask column, and the dataset column
# ---------------------------------------------------------------------------


def two_datasets() -> dict[str, _Holder]:
    """A twelve-sample masked spectrum beside four photometric points."""
    return {"spec": _Holder(a_masked_spectrum()), "phot": _Holder(photometry())}


class TestACollectionInOneTensor:
    """Accept criterion 2: different lengths, one padded tensor, an exact mask."""

    def test_two_datasets_of_different_lengths_share_one_tensor(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        assert layout.row_cap == 12 + 4
        assert encoded.values.shape == (1, 16, layout.columns_total)
        assert layout.bounds == ((0, 12), (12, 16))

    def test_the_mask_column_is_exact_including_on_padded_rows(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets, row_cap=20)
        encoded = encode_observations(datasets, layout=layout)
        mask = np.asarray(encoded.values)[0, :, layout.group("mask").offset]
        expected = np.ones(20)
        expected[[3, 7]] = 0.0  # the spectrum's two masked samples
        expected[16:] = 0.0  # the four padded rows
        assert mask.tolist() == expected.tolist()
        # A padded row is all zero, not merely masked.
        assert np.all(np.asarray(encoded.values)[0, 16:, :] == 0.0)

    def test_the_dataset_column_carries_the_index_in_declared_order(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        column = np.asarray(encoded.values)[0, :, layout.group("dataset").offset]
        assert column[:12].tolist() == [0.0] * 12
        assert column[12:].tolist() == [1.0] * 4

    def test_the_per_set_features_carry_log_n_has_sigma_and_is_complex(self) -> None:
        datasets = {"phot": _Holder(photometry()), "vis": _Holder(visibilities())}
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        group = layout.group("set_features")
        block = np.asarray(encoded.values)[0, :, group.offset : group.stop]
        assert block[0].tolist() == [math.log(4), 1.0, 0.0]
        # The visibility set has one masked point of five, and is complex.
        assert block[4].tolist() == [math.log(4), 1.0, 1.0]

    def test_a_complex_dataset_takes_two_value_columns_and_a_real_one_takes_one(self) -> None:
        datasets = {"phot": _Holder(photometry()), "vis": _Holder(visibilities())}
        layout = EncodingLayout.from_datasets(datasets)
        assert layout.complex_columns is True
        assert layout.group("value").width == 2
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        group = layout.group("value_asinh")
        # The photometric block's imaginary column is zero throughout; the
        # visibility block's is not.
        assert encoded[0, :4, group.offset + 1].tolist() == [0.0] * 4
        assert np.any(encoded[0, 4:, group.offset + 1] != 0.0)


# ---------------------------------------------------------------------------
# 3. Standardisation, and where the statistics come from
# ---------------------------------------------------------------------------


class TestStandardisation:
    """§4: from the observation, once, and recorded in the layout."""

    def test_the_coordinate_columns_span_minus_one_to_one_per_dataset(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        column = encoded[0, :, layout.group("coordinate").offset]
        assert column[:12].min() == pytest.approx(-1.0)
        assert column[:12].max() == pytest.approx(1.0)
        assert column[12:].min() == pytest.approx(-1.0)
        assert column[12:].max() == pytest.approx(1.0)

    def test_the_value_scale_is_the_median_of_the_valid_magnitudes(self) -> None:
        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        values = np.asarray(datasets["spectrum"].observed.values)
        assert layout.datasets[0].value_scale == pytest.approx(float(np.median(np.abs(values))))

    def test_the_statistics_do_not_move_when_a_different_draw_is_encoded(self) -> None:
        """The whole reason the statistics live in the layout rather than the data."""
        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        original = datasets["spectrum"].observed
        loud = original.with_values(np.asarray(original.values) * 1000.0)
        first = np.asarray(encode_observations(datasets, layout=layout).values)
        second = np.asarray(encode({"spectrum": loud}, layout=layout).values)
        coordinate = layout.group("coordinate")
        assert second[0, :, coordinate.offset].tolist() == first[0, :, coordinate.offset].tolist()
        # The values move, and by the full factor: nothing renormalises them.
        asinh = layout.group("value_asinh").offset
        assert float(second[0, 0, asinh]) > float(first[0, 0, asinh])

    def test_the_fourier_features_are_sin_and_cos_of_the_standardised_coordinate(self) -> None:
        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        coordinate = encoded[0, :, layout.group("coordinate").offset]
        group = layout.group("coordinate_features")
        assert group.width == 1 * 2 * DEFAULT_FOURIER_BANDS
        for band in range(DEFAULT_FOURIER_BANDS):
            angle = math.pi * (2.0**band) * coordinate
            assert encoded[0, :, group.offset + 2 * band] == pytest.approx(np.sin(angle))
            assert encoded[0, :, group.offset + 2 * band + 1] == pytest.approx(np.cos(angle))

    def test_the_group_can_be_switched_off(self) -> None:
        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets, fourier_bands=0)
        assert layout.group("coordinate_features").width == 0
        assert layout.columns_total == 1 + 1 + 0 + 1 + 1 + 1 + 1 + 3 + 0

    def test_a_dataset_with_fewer_axes_than_the_widest_is_padded_with_zeros(self) -> None:
        datasets = {"cube": _Holder(a_cube()), "spec": _Holder(a_spectrum())}
        layout = EncodingLayout.from_datasets(datasets)
        assert layout.coordinates == 3
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        group = layout.group("coordinate")
        spectrum_rows = encoded[0, 24:, group.offset : group.stop]
        assert np.any(spectrum_rows[:, 0] != 0.0)
        assert spectrum_rows[:, 1:].tolist() == [[0.0, 0.0]] * 12
        # And an absent axis contributes zero Fourier features, not cos(0) = 1.
        features = layout.group("coordinate_features")
        absent = encoded[0, 24:, features.offset + 2 * DEFAULT_FOURIER_BANDS : features.stop]
        assert not np.any(absent)


class TestTheSigmaColumn:
    """The horizon notes' one requirement of this item: sigma is not optional."""

    def test_it_is_present_by_default(self) -> None:
        layout = EncodingLayout.from_datasets(one("spectrum"))
        assert layout.has_sigma is True
        assert layout.datasets[0].has_sigma is True
        assert layout.group("log_sigma").width == 1

    def test_its_absence_is_recorded_in_the_layout_rather_than_removing_the_column(
        self,
    ) -> None:
        datasets = one("light_curve")
        layout = EncodingLayout.from_datasets(datasets)
        assert layout.has_sigma is False
        assert layout.datasets[0].has_sigma is False
        assert layout.group("log_sigma").width == 1
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        assert encoded[0, :, layout.group("log_sigma").offset].tolist() == [0.0] * 6
        # With no sigma the asinh value stands in for the whitened one (§3).
        value = layout.group("value").offset
        asinh = layout.group("value_asinh").offset
        assert encoded[0, :, value].tolist() == encoded[0, :, asinh].tolist()

    def test_the_whitened_column_is_y_over_sigma(self) -> None:
        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        encoded = np.asarray(encode_observations(datasets, layout=layout).values)
        observed = datasets["spectrum"].observed
        expected = np.asarray(observed.values) / np.asarray(observed.uncertainty)
        assert encoded[0, :, layout.group("value").offset] == pytest.approx(expected)


# ---------------------------------------------------------------------------
# 4. The layout: a hash, a name, and refusals
# ---------------------------------------------------------------------------


class TestTheLayout:
    """A frozen record whose hash is what says two tensors mean the same thing."""

    def test_the_hash_is_stable_and_covers_the_record(self) -> None:
        first = EncodingLayout.from_datasets(two_datasets())
        second = EncodingLayout.from_datasets(two_datasets())
        assert first.hash == second.hash
        assert len(first.hash) == 32
        wider = EncodingLayout.from_datasets(two_datasets(), fourier_bands=2)
        assert wider.hash != first.hash

    def test_it_round_trips_through_its_own_dict(self) -> None:
        layout = EncodingLayout.from_datasets(two_datasets())
        rebuilt = EncodingLayout.from_dict(layout.to_dict())
        assert rebuilt == layout
        assert rebuilt.hash == layout.hash

    def test_the_flat_layout_is_the_same_object_with_its_own_hash(self) -> None:
        datasets = two_datasets()
        flat = EncodingLayout.from_datasets(datasets, kind="flat")
        assert flat.name == "flat"
        assert flat.row_cap == 1
        # Ten valid samples: twelve less two masked, plus four photometric.
        assert flat.columns_total == 14
        assert flat.hash != EncodingLayout.from_datasets(datasets).hash

    def test_a_flat_layout_refuses_a_complex_container_by_name(self) -> None:
        with pytest.raises(EncodingError, match="complex"):
            EncodingLayout.from_datasets(one("visibilities"), kind="flat")

    def test_a_layout_mismatch_is_refused_and_says_which_field_differs(self) -> None:
        layout = EncodingLayout.from_datasets(two_datasets())
        other = {"spec": _Holder(a_spectrum()), "phot": _Holder(photometry())}
        with pytest.raises(EncodingError) as raised:
            layout.check_against(other)
        message = str(raised.value)
        assert "masked samples" in message
        assert "'spec'" in message
        assert layout.hash in message

    def test_a_different_set_of_labels_is_refused_by_name(self) -> None:
        layout = EncodingLayout.from_datasets(two_datasets())
        with pytest.raises(EncodingError, match="dataset labels"):
            layout.check_against(one("spectrum"))

    def test_a_row_cap_below_the_observation_is_refused_by_name(self) -> None:
        with pytest.raises(EncodingError, match="row cap"):
            EncodingLayout.from_datasets(two_datasets(), row_cap=3)

    def test_a_differently_sampled_container_is_refused_at_encode(self) -> None:
        """Amortisation over sampling is a later item; silence would not be."""
        layout = EncodingLayout.from_datasets(one("spectrum"))
        shorter = Spectrum(
            np.linspace(1.0, 10.0, 6) * u.um,
            np.linspace(2.0, 20.0, 6) * u.Jy,
            uncertainty=np.full(6, 0.3) * u.Jy,
        )
        with pytest.raises(EncodingError, match="row"):
            encode({"spectrum": shorter}, layout=layout)

    def test_a_missing_observation_is_refused_by_name(self) -> None:
        layout = EncodingLayout.from_datasets(one("spectrum"))
        with pytest.raises(EncodingError, match="no observation"):
            encode({}, layout=layout)

    def test_a_fully_masked_problem_is_refused_rather_than_encoded_to_nothing(self) -> None:
        grid = np.linspace(1.0, 5.0, 5)
        datasets = {
            "all_gone": _Holder(Spectrum(grid * u.um, grid * u.Jy, mask=np.ones(5, dtype=bool)))
        }
        with pytest.raises(EncodingError, match="empty"):
            EncodingLayout.from_datasets(datasets)

    def test_a_problem_with_no_datasets_is_refused(self) -> None:
        with pytest.raises(EncodingError, match="no datasets"):
            EncodingLayout.from_datasets({})


# ---------------------------------------------------------------------------
# 5. unpack
# ---------------------------------------------------------------------------


class TestUnpack:
    """The one way a network reads the tensor."""

    def test_the_groups_partition_the_columns(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        view = unpack(encoded.values, layout)
        widths = {
            "dataset": view.dataset,
            "coordinate": view.coordinate,
            "coordinate_features": view.coordinate_features,
            "value": view.value,
            "value_asinh": view.value_asinh,
            "log_sigma": view.log_sigma,
            "mask": view.mask,
            "set_features": view.set_features,
            "context": view.context,
        }
        for name, block in widths.items():
            assert block.shape[-1] == layout.group(name).width
        assert view.features.shape[-1] == layout.columns_total - 1
        assert view.context.shape[-1] == 0

    def test_the_valid_flag_is_the_mask_column(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets, row_cap=18)
        view = unpack(encode_observations(datasets, layout=layout).values, layout)
        assert view.valid.shape == (1, 18)
        assert view.valid[0].tolist() == (
            [True, True, True, False, True, True, True, False]
            + [True] * 4
            + [True] * 4
            + [False] * 2
        )

    def test_a_flattened_tensor_is_reshaped_rather_than_refused(self) -> None:
        """``sbi`` flattens ``x`` on some paths; a wrapper must not have to care."""
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        encoded = encode_observations(datasets, layout=layout)
        view = unpack(encoded.matrix, layout)
        assert view.mask.shape == (1, layout.row_cap, 1)

    def test_a_tensor_of_the_wrong_width_is_refused_by_name(self) -> None:
        layout = EncodingLayout.from_datasets(two_datasets())
        with pytest.raises(EncodingError, match="encoded under a different layout"):
            unpack(np.zeros((1, layout.row_cap, layout.columns_total + 1)), layout)

    def test_unpacking_a_flat_layout_is_refused(self) -> None:
        layout = EncodingLayout.from_datasets(two_datasets(), kind="flat")
        with pytest.raises(EncodingError, match="column groups"):
            unpack(np.zeros((1, 1, layout.columns_total)), layout)

    def test_per_dataset_slices_the_rows_in_declared_order(self) -> None:
        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        view = unpack(encode_observations(datasets, layout=layout).values, layout)
        assert [one.label for one in view.per_dataset] == ["spec", "phot"]
        assert [(one.start, one.stop) for one in view.per_dataset] == [(0, 12), (12, 16)]
        assert view.per_dataset[1].valid[0].tolist() == [True] * 4
        assert view.per_dataset[0].features.shape == (1, 12, layout.columns_total - 1)

    def test_a_grid_kinds_slice_reshapes_to_its_axes_exactly(self) -> None:
        """Accept criterion: ``grid_shape``, and rows in C-order, or nothing works."""
        datasets = {"cube": _Holder(a_cube()), "spec": _Holder(a_spectrum())}
        layout = EncodingLayout.from_datasets(datasets)
        cube = datasets["cube"].observed
        view = unpack(encode_observations(datasets, layout=layout).values, layout)
        block = view.per_dataset[0]
        assert block.grid_shape == (2, 3, 4)
        scale = layout.datasets[0].value_scale
        values = scale * np.sinh(np.asarray(block.value_asinh)[0, :, 0])
        assert values.reshape(block.grid_shape) == pytest.approx(np.asarray(cube.values))
        # And the coordinates come back on the axes they were flattened from.
        coordinate = np.asarray(block.coordinate)[0].reshape(*block.grid_shape, 3)
        assert coordinate[:, 0, 0, 0].tolist() == pytest.approx([-1.0, 1.0])
        assert coordinate[0, :, 0, 1].tolist() == pytest.approx([-1.0, 0.0, 1.0])


# ---------------------------------------------------------------------------
# 6. Batches
# ---------------------------------------------------------------------------


class TestABatch:
    """What a training set is encoded from: many draws, one layout."""

    def test_a_container_batch_encodes_into_one_tensor_per_draw(self) -> None:
        from ampere.core.simulate import ContainerBatch

        datasets = one("spectrum")
        layout = EncodingLayout.from_datasets(datasets)
        observed = datasets["spectrum"].observed
        draws = [observed.with_values(np.asarray(observed.values) * factor) for factor in (1, 2, 3)]
        batch = ContainerBatch.stack(draws)
        encoded = encode({"spectrum": batch}, layout=layout, batched=True)
        assert encoded.values.shape == (3, layout.row_cap, layout.columns_total)
        # Everything but the value columns is identical across draws.
        mask = layout.group("mask").offset
        assert np.array_equal(encoded.values[0, :, mask], encoded.values[2, :, mask])
        first = encoded.values[0, :, layout.group("value").offset]
        third = encoded.values[2, :, layout.group("value").offset]
        assert third == pytest.approx(3.0 * first)

    def test_the_datasets_of_one_batch_must_agree_on_the_draw_count(self) -> None:
        from ampere.core.simulate import ContainerBatch

        datasets = two_datasets()
        layout = EncodingLayout.from_datasets(datasets)
        spectrum = datasets["spec"].observed
        points = datasets["phot"].observed
        with pytest.raises(EncodingError, match="same number of draws"):
            encode(
                {
                    "spec": ContainerBatch.stack([spectrum, spectrum]),
                    "phot": ContainerBatch.stack([points]),
                },
                layout=layout,
                batched=True,
            )
