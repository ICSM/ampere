"""The JWST reader against a real ``x1d`` product and rewritten copies of it (W7.3).

The file is a MIRI low-resolution-spectrometer fixed-slit spectrum of the
flux standard BD+60 1753 from the calibration programme 1536
(``tests/data/README.md`` has the provenance, checksum and redistribution
basis). It is a level-3 ``x1d``: one ``EXTRACT1D`` table of 388 rows, written
from long to short wavelength, whose first and last rows are padding — ``DQ``
``DO_NOT_USE``, a ``NaN`` flux, a zero uncertainty — and whose last row
repeats the wavelength of the row before it. Those three facts are what the
round-trip rows hold the reader to; the selection, the masks' other branches
and the unit conversion are exercised on copies rewritten with
``astropy.io.fits``.

The file is vendored; if the vendored copy is ever removed (one ``git rm``),
the fixture downloads the pinned copy from MAST to a cache and skips by name
when offline.
"""

from __future__ import annotations

import hashlib
import os
import urllib.error
import urllib.request
from collections.abc import Callable
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from ampere.interferometry import OIFITSError, read_oifits
from ampere.spectroscopy.jwst import DO_NOT_USE, JWSTError, read_jwst

NAME = "jw01536-o027_t008_miri_p750l_x1d.fits"
SHA256 = "16407b30bbf4bf86837e894cc39b702bf17207a0b7ca7bd56c26d7f23efba886"
URL = f"https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/{NAME}"
DATA = Path(__file__).resolve().parents[1] / "data"
VENDORED = DATA / NAME
OIFITS = DATA / "contest-2008-binary.oifits"

#: The file's own numbers (``tests/data/README.md``).
N_ROWS = 388


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
            with urllib.request.urlopen(URL, timeout=60) as response:
                target.write_bytes(response.read())
        except (urllib.error.URLError, OSError) as exc:
            pytest.skip(f"{NAME} is not vendored and could not be downloaded from {URL}: {exc}")
    assert _sha256(target) == SHA256, f"the downloaded {target} does not match the pinned SHA-256"
    return target


@pytest.fixture(scope="module")
def table(sample: Path) -> fits.FITS_rec:
    with fits.open(sample, memmap=False) as hdus:
        return hdus["EXTRACT1D"].data.copy()


def _rewrite(sample: Path, out: Path, edit: Callable[[fits.HDUList], None]) -> Path:
    with fits.open(sample, memmap=False) as hdus:
        copies = fits.HDUList([hdu.copy() for hdu in hdus])
    edit(copies)
    copies.writeto(out, overwrite=True)
    return out


# ---------------------------------------------------------------------------
# The round trip
# ---------------------------------------------------------------------------


class TestTheRoundTrip:
    def test_counts_order_and_units(self, sample: Path, table: fits.FITS_rec) -> None:
        spectrum = read_jwst(sample)
        # 388 rows: the last repeats the wavelength of the one before and is
        # padding, so it is the duplicate dropped; 387 samples remain.
        assert spectrum.meta["n_rows"] == N_ROWS
        assert spectrum.meta["n_duplicate"] == 1
        assert spectrum.meta["n_nonfinite_wavelength"] == 0
        assert spectrum.n_samples == N_ROWS - 1
        assert spectrum.meta["reordered"] == "reversed"
        assert spectrum.spectral_axis.unit == u.micron
        assert u.Unit(spectrum.unit) == u.Jy  # the file's TUNIT2, read back
        axis = spectrum.spectral_axis.values
        assert np.all(np.diff(axis) > 0.0)
        assert axis == pytest.approx(np.sort(np.unique(table["WAVELENGTH"])), rel=0, abs=0)
        # The unmasked values are the file's, sample for sample.
        rows = np.argsort(table["WAVELENGTH"][:-1])
        keep = ~spectrum.mask
        assert spectrum.values[keep] == pytest.approx(table["FLUX"][:-1][rows][keep], rel=1e-15)
        assert spectrum.uncertainty[keep] == pytest.approx(
            table["FLUX_ERROR"][:-1][rows][keep], rel=1e-15
        )

    def test_the_mask_is_the_do_not_use_bit(self, sample: Path) -> None:
        spectrum = read_jwst(sample)
        dq = np.asarray(spectrum.extra_coords["dq"])
        assert np.array_equal(spectrum.mask, (dq & DO_NOT_USE) != 0)
        # One padding row survives the duplicate rule: the long-wavelength end.
        assert np.flatnonzero(spectrum.mask).tolist() == [N_ROWS - 2]
        meta = spectrum.meta
        assert (meta["n_do_not_use"], meta["n_nonfinite"], meta["n_nonpositive_error"]) == (1, 1, 1)
        assert meta["n_masked"] == 1
        assert np.all(np.isfinite(spectrum.values)) and np.all(spectrum.uncertainty > 0.0)
        assert "npixels" in spectrum.extra_coords

    def test_the_meta(self, sample: Path) -> None:
        meta = read_jwst(sample).meta
        assert {
            key: meta[key]
            for key in (
                "instrument",
                "detector",
                "filter",
                "grating",
                "exp_type",
                "program",
                "observation",
                "target",
                "date_obs",
                "pipeline_version",
                "product",
                "extension",
                "hdu",
                "slit",
                "srctype",
                "columns",
            )
        } == {
            "instrument": "MIRI",
            "detector": "MIRIMAGE",
            "filter": "P750L",
            "grating": None,
            "exp_type": "MIR_LRS-FIXEDSLIT",
            "program": "01536",
            "observation": "027",
            "target": "BD+60 1753",
            "date_obs": "2022-07-08",
            "pipeline_version": "3.0.0",
            "product": "x1d",
            "extension": 1,
            "hdu": 1,
            "slit": "MIR_LRS-FIXEDSLIT",
            "srctype": "POINT",
            "columns": ["FLUX", "FLUX_ERROR"],
        }

    def test_surface_brightness_reads_its_own_unit(self, sample: Path) -> None:
        spectrum = read_jwst(sample, surface_brightness=True)
        assert u.Unit(spectrum.unit) == u.MJy / u.sr
        assert spectrum.meta["columns"] == ["SURF_BRIGHT", "SB_ERROR"]
        assert spectrum.n_samples == N_ROWS - 1


# ---------------------------------------------------------------------------
# Rewritten copies: selection, masks, units, layouts
# ---------------------------------------------------------------------------


def _second_slit(hdus: fits.HDUList) -> None:
    second = hdus["EXTRACT1D"].copy()
    second.ver = 2
    second.header["SLTNAME"] = "OTHER"
    second.header["SOURCEID"] = 7
    second.data["FLUX"] *= 2.0
    hdus.insert(2, second)


class TestTheSelection:
    def test_two_spectra_are_refused_without_a_keyword(self, sample: Path, tmp_path: Path) -> None:
        path = _rewrite(sample, tmp_path / "two.fits", _second_slit)
        with pytest.raises(
            JWSTError,
            match=r"holds 2 spectra, \(EXTRACT1D EXTVER 1, SLTNAME 'MIR_LRS-FIXEDSLIT'.*"
            r"\(EXTRACT1D EXTVER 2, SLTNAME 'OTHER', SOURCEID '7'.*source=.*extension=.*"
            r"does not pick the first",
        ):
            read_jwst(path)

    def test_selected_by_source_or_extension(self, sample: Path, tmp_path: Path) -> None:
        path = _rewrite(sample, tmp_path / "two.fits", _second_slit)
        first = read_jwst(path, source="MIR_LRS-FIXEDSLIT")
        by_slit = read_jwst(path, source="OTHER")
        by_id = read_jwst(path, source=7)
        by_extension = read_jwst(path, extension=2)
        assert first.meta["extension"] == 1 and first.meta["hdu"] == 1
        for chosen in (by_slit, by_id, by_extension):
            assert chosen.meta["extension"] == 2 and chosen.meta["hdu"] == 2
            assert chosen.meta["slit"] == "OTHER" and chosen.meta["source_id"] == "7"
            keep = ~chosen.mask
            assert chosen.values[keep] == pytest.approx(2.0 * first.values[keep])
        assert read_jwst(path, source="OTHER", extension=2).meta["extension"] == 2

    def test_a_wrong_selection_is_refused_by_name(self, sample: Path, tmp_path: Path) -> None:
        path = _rewrite(sample, tmp_path / "two.fits", _second_slit)
        with pytest.raises(
            JWSTError,
            match=r"no spectrum table with SLTNAME, SOURCEID or SRCNAME 'NOPE'.*it has .*'OTHER'",
        ):
            read_jwst(path, source="NOPE")
        with pytest.raises(JWSTError, match=r"no spectrum table with EXTVER 9; it has .*EXTVER 2"):
            read_jwst(path, extension=9)
        with pytest.raises(JWSTError, match=r"SRCNAME 'OTHER' and EXTVER 1"):
            read_jwst(path, source="OTHER", extension=1)

    def test_one_source_in_two_tables_needs_the_extension(
        self, sample: Path, tmp_path: Path
    ) -> None:
        def edit(hdus: fits.HDUList) -> None:
            _second_slit(hdus)
            hdus["EXTRACT1D", 1].header["SOURCEID"] = 7

        path = _rewrite(sample, tmp_path / "orders.fits", edit)
        with pytest.raises(JWSTError, match=r"holds 2 spectra matching 7.*Pass extension="):
            read_jwst(path, source=7)
        assert read_jwst(path, source=7, extension=1).meta["extension"] == 1


class TestTheMasks:
    def test_only_the_do_not_use_bit_masks(self, sample: Path, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            hdus["EXTRACT1D"].data["DQ"][100] = 4  # another bit: informational
            hdus["EXTRACT1D"].data["DQ"][101] = 1 | 4  # DO_NOT_USE among others

        spectrum = read_jwst(_rewrite(sample, tmp_path / "dq.fits", edit))
        dq = np.asarray(spectrum.extra_coords["dq"])
        assert np.array_equal(spectrum.mask, (dq & DO_NOT_USE) != 0)
        assert sorted(dq[dq != 0].tolist()) == [1, 4, 5]
        assert not spectrum.mask[np.flatnonzero(dq == 4)].any()
        assert spectrum.meta["n_do_not_use"] == 2 and spectrum.meta["n_masked"] == 2

    def test_a_non_finite_flux_and_a_zero_error_are_masked_and_counted(
        self, sample: Path, tmp_path: Path
    ) -> None:
        def edit(hdus: fits.HDUList) -> None:
            hdus["EXTRACT1D"].data["FLUX"][10] = np.inf
            hdus["EXTRACT1D"].data["FLUX_ERROR"][20] = 0.0

        spectrum = read_jwst(_rewrite(sample, tmp_path / "bad.fits", edit))
        meta = spectrum.meta
        assert (meta["n_do_not_use"], meta["n_nonfinite"], meta["n_nonpositive_error"]) == (1, 2, 2)
        assert meta["n_masked"] == 3
        # Rows 10 and 20 of a reversed 387-sample axis.
        assert np.flatnonzero(spectrum.mask).tolist() == [
            N_ROWS - 2 - 20,
            N_ROWS - 2 - 10,
            N_ROWS - 2,
        ]
        assert np.all(np.isfinite(spectrum.values)) and np.all(spectrum.uncertainty > 0.0)


class TestTheUnitsAndLayouts:
    def test_a_wavelength_in_nanometres_is_converted(self, sample: Path, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            hdus["EXTRACT1D"].data["WAVELENGTH"] *= 1000.0
            hdus["EXTRACT1D"].columns["WAVELENGTH"].unit = "nm"

        reference = read_jwst(sample)
        spectrum = read_jwst(_rewrite(sample, tmp_path / "nm.fits", edit))
        assert spectrum.spectral_axis.unit == u.micron
        assert spectrum.spectral_axis.values == pytest.approx(
            reference.spectral_axis.values, rel=1e-12
        )

    def test_a_flux_in_millijansky_stays_in_millijansky(self, sample: Path, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            hdus["EXTRACT1D"].columns["FLUX"].unit = "mJy"

        spectrum = read_jwst(_rewrite(sample, tmp_path / "mjy.fits", edit))
        reference = read_jwst(sample)
        assert u.Unit(spectrum.unit) == u.mJy
        keep = ~spectrum.mask
        # The error, still in Jy, is converted to the value's unit.
        assert spectrum.uncertainty[keep] == pytest.approx(1e3 * reference.uncertainty[keep])

    @pytest.mark.parametrize(
        ("unit", "message"), [("", "declares no unit"), ("not-a-unit", "cannot parse")]
    )
    def test_a_unit_the_reader_cannot_use_is_refused(
        self, sample: Path, tmp_path: Path, unit: str, message: str
    ) -> None:
        def edit(hdus: fits.HDUList) -> None:
            if unit:
                hdus["EXTRACT1D"].columns["FLUX"].unit = unit
            else:
                table = hdus["EXTRACT1D"]
                columns = [
                    fits.Column(c.name, c.format, array=table.data[c.name])
                    if c.name == "FLUX"
                    else c
                    for c in table.columns
                ]
                header = table.header.copy()
                del header["TUNIT2"]
                hdus["EXTRACT1D"] = fits.BinTableHDU.from_columns(
                    columns, header=header, name="EXTRACT1D"
                )

        with pytest.raises(JWSTError, match=f"column FLUX .*{message}"):
            read_jwst(_rewrite(sample, tmp_path / "unit.fits", edit))

    def test_a_c1d_layout(self, sample: Path, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            hdus["EXTRACT1D"].columns.change_name("FLUX_ERROR", "ERROR")
            hdus["EXTRACT1D"].name = "COMBINE1D"

        spectrum = read_jwst(_rewrite(sample, tmp_path / "c1d.fits", edit))
        assert spectrum.meta["product"] == "c1d"
        assert spectrum.meta["columns"] == ["FLUX", "ERROR"]
        assert spectrum.n_samples == N_ROWS - 1

    def test_no_flux_column_names_the_keyword(self, sample: Path, tmp_path: Path) -> None:
        def edit(hdus: fits.HDUList) -> None:
            table = hdus["EXTRACT1D"]
            kept = [c for c in table.columns if c.name not in ("FLUX", "FLUX_ERROR")]
            hdus["EXTRACT1D"] = fits.BinTableHDU.from_columns(
                kept, header=table.header, name="EXTRACT1D"
            )

        path = _rewrite(sample, tmp_path / "sb.fits", edit)
        with pytest.raises(JWSTError, match=r"no FLUX, FLUX_ERROR column.*surface_brightness=True"):
            read_jwst(path)
        assert u.Unit(read_jwst(path, surface_brightness=True).unit) == u.MJy / u.sr

    def test_a_vector_layout_is_refused(self, tmp_path: Path) -> None:
        path = tmp_path / "vector.fits"
        columns = [
            fits.Column(name, "5D", unit=unit, array=np.ones((3, 5)))
            for name, unit in (("WAVELENGTH", "um"), ("FLUX", "Jy"), ("FLUX_ERROR", "Jy"))
        ]
        fits.HDUList(
            [fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="EXTRACT1D")]
        ).writeto(path)
        with pytest.raises(JWSTError, match=r"holds \(5,\) values per row"):
            read_jwst(path)


# ---------------------------------------------------------------------------
# The wrong format, both ways, and a missing path
# ---------------------------------------------------------------------------


class TestTheWrongFormat:
    def test_read_jwst_refuses_an_oifits_file(self) -> None:
        with pytest.raises(
            JWSTError,
            match=r"contest-2008-binary\.oifits has no EXTRACT1D or COMBINE1D table "
            r"\(it has \['OI_ARRAY', 'OI_T3', 'OI_TARGET', 'OI_VIS2', 'OI_WAVELENGTH'\]\)",
        ):
            read_jwst(OIFITS)

    def test_read_oifits_refuses_an_x1d_file(self, sample: Path) -> None:
        with pytest.raises(
            OIFITSError,
            match=rf"{NAME.replace('.', r'\.')} has no OI_VIS, OI_VIS2 or OI_T3 table "
            r"\(it has \['ASDF', 'EXTRACT1D'\]\)",
        ):
            read_oifits(sample)

    def test_a_missing_path(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError, match=r"no JWST file at .*nowhere\.fits"):
            read_jwst(tmp_path / "nowhere.fits")


# ---------------------------------------------------------------------------
# The fit: the photometry_spectra page's front-door block, at a token budget
# ---------------------------------------------------------------------------

#: emcee's floor is twice the four free parameters; a token budget, since the
#: row asserts that the reader's Spectrum composes and fits, not what it finds.
RUN_WALKERS, RUN_STEPS, RUN_BURN_IN = 16, 8, 2


class TestTheFit:
    def test_the_spectrum_fits_end_to_end_with_the_flexible_likelihood(self, sample: Path) -> None:
        # The block docs/source/photometry_spectra.rst shows.
        import scipy.stats as st

        from ampere.backends.reference import PowerLaw
        from ampere.core import Dataset, FittingProblem, GaussianProcessNoise, Instrument
        from ampere.core import Likelihood, Matern32
        from ampere.inference import EmceeEngine
        from ampere.spectroscopy import GaussianFamily, LSFConvolution, Resample, read_jwst

        spectrum = read_jwst(sample)
        model = PowerLaw(
            np.geomspace(3.0, 15.0, 400),
            norm=st.loguniform(1e-3, 1e-1),
            index=st.uniform(-4.0, 4.0),
            reference_wavelength=7.5,
        )
        instrument = Instrument(
            [LSFConvolution(resolving_power=100.0), Resample(spectrum.spectral_axis.values)]
        )
        kernel = Matern32(
            st.loguniform(1e-5, 1e-2),
            st.loguniform(0.1, 10.0),
            amplitude_unit=u.Jy,
            length_scale_unit=u.um,
        )
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel))
        problem = FittingProblem(model, {"miri": Dataset(spectrum, instrument, likelihood)})
        # initial="prior" skips the default optimiser warm start (some forty
        # seconds here): the one departure from the page's block, for budget.
        engine = EmceeEngine(problem, walkers=RUN_WALKERS)
        run = engine.run(RUN_STEPS, burn_in=RUN_BURN_IN, initial="prior")

        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) == {
            "model.norm",
            "model.index",
            "miri.likelihood.amplitude",
            "miri.likelihood.length_scale",
        }
        assert posterior.sizes["chain"] == RUN_WALKERS
        assert posterior.sizes["draw"] == RUN_STEPS - RUN_BURN_IN
        assert all(np.all(np.isfinite(posterior[name].values)) for name in posterior.data_vars)
