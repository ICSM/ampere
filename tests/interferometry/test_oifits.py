"""The OIFITS reader against a real file and that file's own truth (W6.12).

The file is the 2008 Imaging Beauty Contest's binary (``tests/data/README.md``
has the provenance, checksum and redistribution basis). Its truth is
published — two uniform discs, 5.0 mas apart at 30° east of north, flux ratio
8.9 — so the rows here hold the reader to **the file's own answer**, written
out in a few lines of numpy below, rather than to ampere's models. A sign error
in the closure-phase ordering, a station-order slip or a metres/wavelengths
confusion each moves the residuals by tens of sigma; the convention row
asserts that they do not, and a mirrored control asserts that the row would
notice if they did.

The file is vendored; if the vendored copy is ever removed (a licence ruling
costs one ``git rm``), the fixture downloads the pinned copy to a cache and
skips by name when offline.
"""

from __future__ import annotations

import hashlib
import itertools
import os
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.special
import scipy.stats as st
from astropy.io import fits

from ampere.core import Dataset, DatasetCollection, FittingProblem, Instrument, Likelihood
from ampere.interferometry import (
    Binary,
    ClosurePhase,
    ClosurePhases,
    FourierSample,
    GaussianFamily,
    OIFITSError,
    SquaredAmplitude,
    VisibilitySet,
    VonMisesFamily,
    read_oifits,
)

NAME = "contest-2008-binary.oifits"
SHA256 = "2476bd412d25ddf9af3ee7002f8998a1e1e6f9fbbfbc60310b7412310c235005"
URL = (
    "https://raw.githubusercontent.com/emmt/OIFITS.jl/"
    "0978576aeb42e25fa56223853997d9ddf79c83ac/test/contest-2008-binary.oifits"
)
VENDORED = Path(__file__).resolve().parents[1] / "data" / NAME

#: The contest's published truth (``CONTENTS.md``).
SEPARATION = 5.0  # mas
POSITION_ANGLE = np.deg2rad(30.0)  # east of north, bright -> faint
BRIGHTNESS_RATIO = 8.9
DIAMETERS = (1.2, 0.75)  # mas, primary and secondary

MAS = np.pi / (180.0 * 3600.0 * 1000.0)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture(scope="module")
def sample() -> Path:
    """The vendored file, or the pinned copy downloaded to a cache; checksum-verified."""
    if VENDORED.is_file():
        assert _sha256(VENDORED) == SHA256, f"{VENDORED} does not match tests/data/README.md"
        return VENDORED
    cache = Path(os.environ.get("AMPERE_TEST_DATA", Path.home() / ".cache" / "ampere-test-data"))
    target = cache / NAME
    if not target.is_file():
        cache.mkdir(parents=True, exist_ok=True)
        try:
            with urllib.request.urlopen(URL, timeout=30) as response:
                target.write_bytes(response.read())
        except (urllib.error.URLError, OSError) as exc:
            pytest.skip(f"{NAME} is not vendored and could not be downloaded from {URL}: {exc}")
    assert _sha256(target) == SHA256, f"the downloaded {target} does not match the pinned SHA-256"
    return target


@pytest.fixture(scope="module")
def data(sample: Path) -> Any:
    return read_oifits(sample)


def truth_visibility(
    u_pts: np.ndarray, v_pts: np.ndarray, position_angle: float = POSITION_ANGLE
) -> np.ndarray:
    """The contest binary's normalised complex visibility, from its published truth.

    Two uniform discs, ``2 J1(pi theta rho) / (pi theta rho)`` each, the
    secondary offset to ``(x, y) = s (sin p, cos p)`` with ``x`` east and ``y``
    north, and the phase factor ``exp(-2 pi i (u x + v y))`` — the sign
    :class:`~ampere.backends.reference.FourierSample` uses (its DFT kernel is
    ``exp(-2j pi outer(u, x))``). Normalised to one at zero spacing.
    """
    rho = np.hypot(u_pts, v_pts)

    def disc(diameter: float) -> np.ndarray:
        z = np.pi * diameter * MAS * rho
        return 2.0 * scipy.special.j1(z) / z

    ratio = 1.0 / BRIGHTNESS_RATIO
    x = SEPARATION * np.sin(position_angle) * MAS
    y = SEPARATION * np.cos(position_angle) * MAS
    phase = np.exp(-2j * np.pi * (u_pts * x + v_pts * y))
    return (disc(DIAMETERS[0]) + ratio * disc(DIAMETERS[1]) * phase) / (1.0 + ratio)


def truth_closure_phase(t3: ClosurePhases, position_angle: float = POSITION_ANGLE) -> np.ndarray:
    """``arg(V_ij V_jk V_ki)`` of the truth on the container's own triangles."""
    u1, v1, u2, v2 = (t3.axis(name).values for name in ("u1", "v1", "u2", "v2"))
    return np.angle(
        truth_visibility(u1, v1, position_angle)
        * truth_visibility(u2, v2, position_angle)
        * truth_visibility(-(u1 + u2), -(v1 + v2), position_angle)
    )


def wrapped(angle: np.ndarray) -> np.ndarray:
    return np.angle(np.exp(1j * angle))


# ---------------------------------------------------------------------------
# (b) The round trip: counts, units, masks, ordering, labels
# ---------------------------------------------------------------------------


class TestTheRoundTrip:
    def test_the_file_and_its_selection(self, data: Any) -> None:
        assert (data.target, data.instrument, data.array, data.revision) == (
            "Gam_Vic",
            "MIRC_H",
            "CHARA",
            1,
        )
        assert data.visibilities is None  # no OI_VIS table: absent, never empty
        assert data.wavelengths.unit == u.m and data.wavelengths.size == 8
        assert data.wavelengths.to_value(u.um) == pytest.approx(np.linspace(1.5, 1.75, 8), 1e-6)
        assert data.bandwidths.unit == u.m

    def test_counts_one_sample_per_row_and_channel(self, data: Any) -> None:
        assert isinstance(data.squared_visibilities, VisibilitySet)
        assert isinstance(data.closure_phases, ClosurePhases)
        assert data.squared_visibilities.n_samples == 75 * 8
        assert data.closure_phases.n_samples == 100 * 8

    def test_units(self, data: Any) -> None:
        v2, t3 = data.squared_visibilities, data.closure_phases
        for axis in (v2.u, v2.v, t3.u1, t3.v1, t3.u2, t3.v2):
            assert axis.unit is None or axis.unit.physical_type == "dimensionless"
        for axis in (v2.spectral_axis, t3.spectral_axis):
            assert axis.unit.physical_type == "length"
        assert v2.unit is None  # VIS2DATA is a pure number
        assert t3.unit == u.rad
        assert np.all(np.abs(t3.values) <= np.pi)
        assert np.all(t3.uncertainty > 0.0) and np.max(t3.uncertainty) < 1.0  # rad, not deg

    def test_the_coordinates_are_metres_over_the_channel_wavelength(
        self, sample: Path, data: Any
    ) -> None:
        with fits.open(sample) as hdus:
            table = hdus["OI_VIS2"].data
            waves = np.asarray(hdus["OI_WAVELENGTH"].data["EFF_WAVE"], dtype=float)
            ucoord = np.asarray(table["UCOORD"], dtype=float)
            vis2 = np.asarray(table["VIS2DATA"], dtype=float)
        v2 = data.squared_visibilities
        assert v2.u.values == pytest.approx((ucoord[:, None] / waves).reshape(-1), rel=1e-12)
        spectral = (v2.spectral_axis.values * v2.spectral_axis.unit).to_value(u.m)
        assert spectral == pytest.approx(np.tile(waves, 75), rel=1e-12)
        # The data stay as measured: squared visibilities, not amplitudes.
        assert v2.values == pytest.approx(vis2.reshape(-1), rel=1e-12)
        assert v2.meta["observable"] == "squared_visibility"

    def test_the_mask_is_the_flag_and_the_meta_counts_it(self, data: Any) -> None:
        for container in (data.squared_visibilities, data.closure_phases):
            assert container.mask is not None and not container.mask.any()
            assert container.meta["n_flagged"] == container.meta["n_nonfinite"] == 0
            assert container.meta["file"] == NAME
            assert container.meta["oi_revn"] == 1
            assert container.meta["date_obs"] == "2007-05-11"
        assert data.closure_phases.meta["tables"] == ["OI_T3 (HDU 5)"]

    def test_the_canonical_triangle_ordering(self, data: Any) -> None:
        t3 = data.closure_phases
        order = {"S1": 0, "S2": 1, "E1": 2, "W1": 3, "W2": 4, "E2": 5}
        triangles = [label.split("-") for label in t3.extra_coords["triangle"]]
        assert all(order[a] < order[b] < order[c] for a, b, c in triangles)
        for (a, b, c), ij, jk, ki in zip(
            triangles,
            t3.extra_coords["baseline_ij"],
            t3.extra_coords["baseline_jk"],
            t3.extra_coords["baseline_ki"],
            strict=True,
        ):
            assert (ij, jk, ki) == (f"{a}-{b}", f"{b}-{c}", f"{c}-{a}")
        u3, v3 = t3.implied_baseline()
        assert u3 == pytest.approx(-(t3.u1.values + t3.u2.values))
        assert v3 == pytest.approx(-(t3.v1.values + t3.v2.values))

    def test_a_triangles_baselines_are_the_visibilities_baselines(self, data: Any) -> None:
        """``b_ij`` and ``b_jk`` of every triangle are the ``OI_VIS2`` rows of the same label."""
        v2, t3 = data.squared_visibilities, data.closure_phases
        coverage = {
            (label, mjd, channel): (u_pt, v_pt)
            for label, mjd, channel, u_pt, v_pt in zip(
                v2.extra_coords["baseline"],
                v2.extra_coords["mjd"],
                v2.extra_coords["channel"],
                v2.u.values,
                v2.v.values,
                strict=True,
            )
        }
        matched = 0
        for stored, (u_name, v_name) in (
            ("baseline_ij", ("u1", "v1")),
            ("baseline_jk", ("u2", "v2")),
        ):
            for label, mjd, channel, u_pt, v_pt in zip(
                t3.extra_coords[stored],
                t3.extra_coords["mjd"],
                t3.extra_coords["channel"],
                t3.axis(u_name).values,
                t3.axis(v_name).values,
                strict=True,
            ):
                key = (label, mjd, channel)
                if key in coverage:
                    assert (u_pt, v_pt) == pytest.approx(coverage[key], rel=1e-9)
                    matched += 1
        assert matched == 2 * t3.n_samples

    def test_the_labels(self, data: Any) -> None:
        v2 = data.squared_visibilities
        assert v2.extra_coords["baseline"][0] == "S1-S2"
        assert sorted(set(v2.extra_coords["channel"].tolist())) == list(range(8))
        assert set(v2.extra_coords) == {"baseline", "channel", "mjd", "eff_band"}
        assert set(data.closure_phases.extra_coords) == {
            "triangle",
            "baseline_ij",
            "baseline_jk",
            "baseline_ki",
            "channel",
            "mjd",
            "eff_band",
        }


# ---------------------------------------------------------------------------
# (c) The convention check: the file's own truth
# ---------------------------------------------------------------------------


class TestTheConventionsAgainstTheFilesTruth:
    """The reader's containers, scored against the contest's published binary.

    Tolerances (measured 2026-10-05): ``|V|**2`` agrees everywhere within
    four sigma (median 0.68). The closure phases agree with median 0.66 sigma
    and a 95th percentile of 2.3 sigma; fourteen of the 800 samples, all in the
    1.50 µm channel on triangles whose triple amplitude is ~0.003 (near a
    visibility null), sit beyond 5 sigma, because the contest simulated the
    data from a pixelised, tapered image rather than from the closed form. So
    the phase row asserts on the bulk — median and 95th percentile — which a
    convention error moves by tens of sigma (the mirrored control: median 47).
    """

    def test_squared_visibilities_match_the_truth(self, data: Any) -> None:
        v2 = data.squared_visibilities
        model = np.abs(truth_visibility(v2.u.values, v2.v.values)) ** 2
        pulls = np.abs(model - v2.values) / v2.uncertainty
        assert np.max(pulls) < 5.0
        assert np.median(pulls) < 1.0

    def test_closure_phases_match_the_truth(self, data: Any) -> None:
        t3 = data.closure_phases
        pulls = np.abs(wrapped(truth_closure_phase(t3) - t3.values)) / t3.uncertainty
        assert np.median(pulls) < 1.0
        assert np.quantile(pulls, 0.95) < 3.0

    def test_the_mirrored_source_fails_the_same_row(self, data: Any) -> None:
        """The control: a sign error is a point reflection, and the row must notice it."""
        t3 = data.closure_phases
        mirrored = truth_closure_phase(t3, POSITION_ANGLE + np.pi)
        pulls = np.abs(wrapped(mirrored - t3.values)) / t3.uncertainty
        assert np.median(pulls) > 10.0


# ---------------------------------------------------------------------------
# Ruling 4 off the identity: permuted triangles and reversed baselines
# ---------------------------------------------------------------------------


def _rewrite(sample: Path, out: Path, edit: Any) -> Path:
    with fits.open(sample) as hdus:
        copies = fits.HDUList([hdu.copy() for hdu in hdus])
    edit(copies)
    copies.writeto(out, overwrite=True)
    return out


class TestTheOrderingOffTheIdentity:
    """Every triangle in the contest file is already sorted; these rows are not.

    Each ``OI_T3`` row is rewritten in another of its six station orders, with
    ``(U1, V1)``, ``(U2, V2)`` the file's ``b_ab``, ``b_bc`` for that order and
    ``T3PHI`` **computed from the truth on those baselines** — not from the
    parity rule the reader implements — so the read-back is held to the truth
    and not to the reader's own arithmetic.
    """

    @pytest.mark.parametrize("permutation", list(itertools.permutations(range(3)))[1:])
    def test_any_station_order_reads_back_canonically(
        self, sample: Path, data: Any, tmp_path: Path, permutation: tuple[int, int, int]
    ) -> None:
        waves = data.wavelengths.to_value(u.m)

        def edit(hdus: fits.HDUList) -> None:
            table = hdus["OI_T3"].data
            for row in range(len(table)):
                ijk = [int(s) for s in table["STA_INDEX"][row]]
                positions = {
                    ijk[0]: np.zeros(2),
                    ijk[1]: np.array([table["U1COORD"][row], table["V1COORD"][row]]),
                }
                positions[ijk[2]] = positions[ijk[1]] + np.array(
                    [table["U2COORD"][row], table["V2COORD"][row]]
                )
                a, b, c = (ijk[p] for p in permutation)
                b_ab, b_bc = positions[b] - positions[a], positions[c] - positions[b]
                b_ca = positions[a] - positions[c]
                table["STA_INDEX"][row] = (a, b, c)
                table["U1COORD"][row], table["V1COORD"][row] = b_ab
                table["U2COORD"][row], table["V2COORD"][row] = b_bc
                bispectrum = (
                    truth_visibility(b_ab[0] / waves, b_ab[1] / waves)
                    * truth_visibility(b_bc[0] / waves, b_bc[1] / waves)
                    * truth_visibility(b_ca[0] / waves, b_ca[1] / waves)
                )
                table["T3PHI"][row] = np.rad2deg(np.angle(bispectrum))

        permuted = read_oifits(_rewrite(sample, tmp_path / "permuted.oifits", edit))
        t3 = permuted.closure_phases
        reference = data.closure_phases
        for name in ("u1", "v1", "u2", "v2"):
            assert t3.axis(name).values == pytest.approx(reference.axis(name).values, rel=1e-9)
        assert list(t3.extra_coords["triangle"]) == list(reference.extra_coords["triangle"])
        assert np.max(np.abs(wrapped(t3.values - truth_closure_phase(t3)))) < 1e-5

    def test_a_reversed_baseline_is_negated(self, sample: Path, data: Any, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            table = hdus["OI_VIS2"].data
            table["STA_INDEX"] = table["STA_INDEX"][:, ::-1]
            table["UCOORD"] = -table["UCOORD"]
            table["VCOORD"] = -table["VCOORD"]

        reversed_ = read_oifits(_rewrite(sample, tmp_path / "reversed.oifits", edit))
        v2, reference = reversed_.squared_visibilities, data.squared_visibilities
        assert v2.u.values == pytest.approx(reference.u.values, rel=1e-12)
        assert v2.v.values == pytest.approx(reference.v.values, rel=1e-12)
        assert list(v2.extra_coords["baseline"]) == list(reference.extra_coords["baseline"])


class TestComplexVisibilities:
    """``OI_VIS``, which the contest file lacks: built from the truth, reversed rows included."""

    def test_an_oi_vis_table_reads_back_as_the_complex_truth(
        self, sample: Path, data: Any, tmp_path: Path
    ) -> None:
        waves = data.wavelengths.to_value(u.m)

        def edit(hdus: fits.HDUList) -> None:
            vis2 = hdus["OI_VIS2"]
            rows = vis2.data
            stations = np.array(rows["STA_INDEX"])
            ucoord, vcoord = np.array(rows["UCOORD"]), np.array(rows["VCOORD"])
            reverse = np.arange(len(rows)) % 2 == 1  # every other row listed (j, i)
            stations[reverse] = stations[reverse][:, ::-1]
            ucoord[reverse], vcoord[reverse] = -ucoord[reverse], -vcoord[reverse]
            truth = truth_visibility(ucoord[:, None] / waves, vcoord[:, None] / waves)
            flags = np.zeros(truth.shape, dtype=bool)
            flags[0, 3] = True
            amplitude = np.abs(truth)
            amplitude[1, 2] = np.nan
            n = truth.shape[1]
            columns = [
                fits.Column("TARGET_ID", "1I", array=rows["TARGET_ID"]),
                fits.Column("TIME", "1D", array=rows["TIME"]),
                fits.Column("MJD", "1D", array=rows["MJD"]),
                fits.Column("INT_TIME", "1D", array=rows["INT_TIME"]),
                fits.Column("VISAMP", f"{n}D", array=amplitude),
                fits.Column("VISAMPERR", f"{n}D", array=np.full(truth.shape, 0.01)),
                fits.Column("VISPHI", f"{n}D", unit="deg", array=np.rad2deg(np.angle(truth))),
                fits.Column("VISPHIERR", f"{n}D", unit="deg", array=np.ones(truth.shape)),
                fits.Column("UCOORD", "1D", unit="m", array=ucoord),
                fits.Column("VCOORD", "1D", unit="m", array=vcoord),
                fits.Column("STA_INDEX", "2I", array=stations),
                fits.Column("FLAG", f"{n}L", array=flags),
            ]
            table = fits.BinTableHDU.from_columns(columns, header=vis2.header.copy(), name="OI_VIS")
            hdus.append(table)

        read = read_oifits(_rewrite(sample, tmp_path / "vis.oifits", edit))
        vis = read.visibilities
        assert isinstance(vis, VisibilitySet) and np.iscomplexobj(vis.values)
        assert vis.n_samples == 600 and vis.meta["observable"] == "complex_visibility"
        assert vis.u.values == pytest.approx(data.squared_visibilities.u.values, rel=1e-12)
        assert vis.mask[3] and vis.mask[8 + 2] and vis.mask.sum() == 2
        assert vis.meta["n_flagged"] == 1 and vis.meta["n_nonfinite"] == 1
        good = ~vis.mask
        expected = truth_visibility(vis.u.values, vis.v.values)
        assert np.max(np.abs(vis.values[good] - expected[good])) < 1e-9


# ---------------------------------------------------------------------------
# (d) The fixtures accept the containers: they compose and evaluate
# ---------------------------------------------------------------------------


class TestTheContainersCompose:
    FIELD_OF_VIEW = 16.0 * u.mas

    def test_a_binary_fit_to_the_file_evaluates(self, data: Any) -> None:
        """Both observables on one ``sky`` channel, the reference ``Binary``, no fit."""
        v2, t3 = data.squared_visibilities, data.closure_phases
        vis2_instrument = Instrument(
            [
                FourierSample.from_observed(v2, field_of_view=self.FIELD_OF_VIEW),
                SquaredAmplitude(normalisation=1.0 * u.Jy),
            ],
            channel="sky",
            label="vis2",
        )
        t3_instrument = Instrument(
            [FourierSample.from_observed(t3, field_of_view=self.FIELD_OF_VIEW), ClosurePhase()],
            channel="sky",
            label="t3",
        )
        model = Binary.on_field(
            self.FIELD_OF_VIEW,
            4,
            channels="sky",
            component_fwhm=0.9,
            separation=st.uniform(2.0, 8.0),
            position_angle=float(POSITION_ANGLE),
            flux_ratio=st.uniform(0.02, 0.5),
            flux=1.0,
        )
        problem = FittingProblem(
            model,
            DatasetCollection(
                {
                    "vis2": Dataset(
                        v2,
                        vis2_instrument,
                        likelihood=Likelihood(GaussianFamily()),
                        label="vis2",
                    ),
                    "t3": Dataset(
                        t3,
                        t3_instrument,
                        likelihood=Likelihood(VonMisesFamily()),
                        label="t3",
                    ),
                }
            ),
            seed=20261005,
        )
        truth = {"model.separation": SEPARATION, "model.flux_ratio": 1.0 / BRIGHTNESS_RATIO}
        at_truth = problem.log_prob(truth)
        assert np.isfinite(at_truth)
        assert at_truth > problem.log_prob({**truth, "model.separation": 8.0})


# ---------------------------------------------------------------------------
# (e) The refusals, each by name
# ---------------------------------------------------------------------------


def _second_instrument(hdus: fits.HDUList) -> None:
    for name in ("OI_WAVELENGTH", "OI_VIS2"):
        copy = hdus[name].copy()
        copy.header["INSNAME"] = "MIRC_K"
        hdus.append(copy)


def _second_target(hdus: fits.HDUList) -> None:
    targets = hdus["OI_TARGET"]
    rows = fits.FITS_rec.from_columns(targets.columns, nrows=2, fill=False)
    rows["TARGET_ID"][1] = 2
    rows["TARGET"][1] = "Other"
    hdus["OI_TARGET"] = fits.BinTableHDU(rows, header=targets.header, name="OI_TARGET")
    hdus["OI_VIS2"].data["TARGET_ID"][::2] = 2


class TestTheRefusals:
    def test_two_instruments_need_insname(self, sample: Path, tmp_path: Path) -> None:
        path = _rewrite(sample, tmp_path / "two.oifits", _second_instrument)
        with pytest.raises(OIFITSError, match=r"2 instruments, \['MIRC_H', 'MIRC_K'\].*insname="):
            read_oifits(path)
        chosen = read_oifits(path, insname="MIRC_K")
        assert chosen.instrument == "MIRC_K" and chosen.closure_phases is None
        assert chosen.squared_visibilities is not None
        with pytest.raises(OIFITSError, match="no instrument 'NOPE'"):
            read_oifits(path, insname="NOPE")

    def test_two_targets_need_target(self, sample: Path, tmp_path: Path) -> None:
        path = _rewrite(sample, tmp_path / "targets.oifits", _second_target)
        with pytest.raises(OIFITSError, match=r"2 targets, \['Gam_Vic', 'Other'\].*target="):
            read_oifits(path)
        other = read_oifits(path, target="Other")
        assert other.squared_visibilities is not None
        assert other.squared_visibilities.n_samples == 38 * 8
        assert other.closure_phases is None  # its OI_T3 rows are all Gam_Vic's

    def test_a_file_with_no_oi_tables(self, tmp_path: Path) -> None:
        path = tmp_path / "empty.fits"
        fits.HDUList(
            [
                fits.PrimaryHDU(),
                fits.BinTableHDU.from_columns([fits.Column("X", "1D", array=[1.0])]),
            ]
        ).writeto(path)
        with pytest.raises(OIFITSError, match="no OI_VIS, OI_VIS2 or OI_T3 table"):
            read_oifits(path)

    def test_a_missing_path(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match=r"no OIFITS file at .*nowhere\.oifits"):
            read_oifits(tmp_path / "nowhere.oifits")

    def test_flags_and_non_finite_values_are_masked_and_counted(
        self, sample: Path, tmp_path: Path
    ) -> None:
        def edit(hdus: fits.HDUList) -> None:
            table = hdus["OI_T3"].data
            table["FLAG"][0, 0] = True
            table["T3PHI"][1, 1] = np.nan
            table["T3PHIERR"][2, 2] = 0.0

        t3 = read_oifits(_rewrite(sample, tmp_path / "flags.oifits", edit)).closure_phases
        assert np.flatnonzero(t3.mask).tolist() == [0, 8 + 1, 16 + 2]
        assert (t3.meta["n_flagged"], t3.meta["n_nonfinite"], t3.meta["n_nonpositive_error"]) == (
            1,
            1,
            1,
        )
        assert np.all(np.isfinite(t3.values)) and np.all(t3.uncertainty > 0.0)
