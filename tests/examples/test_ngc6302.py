"""``examples/ngc6302`` builds, negotiates and fits end to end (W6.13 (2)).

Fast and always on, like :mod:`tests.examples.test_linear_sed` and
:mod:`tests.examples.test_modified_blackbody`: the opacity buffers load with
the right shape and the two unit conversions; one evaluation equals the
legacy ``SpectrumNGC6302`` at a fixed theta to ``rtol=1e-10``; the real data
load applies the 25-120 micron selection and the five-per-cent uncertainty
rule; the problem builds with the fifteen model parameters plus the
calibration and GP ones; a tiny-budget emcee fit runs; and
``main --synthetic --quick`` (sized down) prints a report. The dust-mass row
is in :mod:`tests.examples.test_ngc6302`'s later commit
(``TestDustMass``, added alongside :mod:`examples.ngc6302.dust_mass`).

The full-budget coverage run this item is accepted on is **not** part of
this suite -- as :mod:`tests.examples.test_linear_sed` and
:mod:`tests.examples.test_modified_blackbody` explain at length, it is one
to two orders of magnitude slower than everything else here (a single
GP-on evaluation of this model costs roughly 100x a bare log_prob call on
the ``linear_sed``/``modified_blackbody`` twins) and is instead verified by
the branch report and :mod:`examples.ngc6302.ngc6302`'s own docstring,
which have the numbers from doing exactly that.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from examples.ngc6302 import generators
from examples.ngc6302.ngc6302 import (
    DEFAULT_GRID,
    QUALIFIED_TRUTH,
    build_instrument,
    build_model,
    build_problem,
    fit,
    main,
    recovers_truth,
    report,
)

TINY = {"walkers": 32, "steps": 10, "burn_in": 3}


class TestOpacityBuffers:
    """The eight opacity tables load as sixteen buffers, converted once."""

    def test_every_species_has_a_wavelength_and_opacity_buffer(self) -> None:
        model = build_model()
        for name in generators.SPECIES:
            assert f"{name}_wavelength" in model.buffers
            assert f"{name}_opacity" in model.buffers
            wavelength = model.buffers[f"{name}_wavelength"].array
            opacity = model.buffers[f"{name}_opacity"].array
            assert wavelength.shape == opacity.shape
            assert wavelength.ndim == 1 and wavelength.size > 0

    def test_enstatite_and_diopside_are_converted_by_1e_4(self) -> None:
        model = build_model()
        filenames = np.loadtxt(generators.OPACITY_FILE_LIST, dtype=str)
        for name, filename in zip(generators.SPECIES, filenames, strict=True):
            if name not in ("enstatite", "diopside"):
                continue
            raw_opacity = np.loadtxt(
                generators.OPACITY_DIRECTORY / str(filename), comments="#", unpack=True
            )[1]
            np.testing.assert_allclose(model.buffers[f"{name}_opacity"].array, raw_opacity * 1e-4)

    def test_calcite_and_dolomite_are_converted_by_the_density_factor(self) -> None:
        model = build_model()
        factors = {"calcite": 2.71 * (4.0 / 3.0) * 1e-4, "dolomite": 2.87 * (4.0 / 3.0) * 1e-4}
        filenames = np.loadtxt(generators.OPACITY_FILE_LIST, dtype=str)
        for name, filename in zip(generators.SPECIES, filenames, strict=True):
            if name not in factors:
                continue
            raw_opacity = np.loadtxt(
                generators.OPACITY_DIRECTORY / str(filename), comments="#", unpack=True
            )[1]
            np.testing.assert_allclose(
                model.buffers[f"{name}_opacity"].array, raw_opacity * factors[name]
            )


class TestModelEqualsLegacy:
    """One evaluation, at a fixed theta, equals ``examples.NGC6302.SpectrumNGC6302``."""

    def test_evaluate_matches_legacy_model(self, monkeypatch: pytest.MonkeyPatch) -> None:
        import examples.NGC6302 as legacy_module

        # The legacy __init__ resolves its opacity directory from os.getcwd(),
        # and __call__ references a bare global `wavelengths` (not
        # self.wavelength) that only exists when the script runs as
        # __main__ -- both reproduced here, inside the test only (ruling 1).
        monkeypatch.chdir(Path(__file__).resolve().parents[2] / "examples")
        legacy_module.wavelengths = DEFAULT_GRID
        legacy_model = legacy_module.SpectrumNGC6302(DEFAULT_GRID)
        legacy_model(
            generators.TRUTH["logacold0"],
            generators.TRUTH["logacold1"],
            generators.TRUTH["logacold2"],
            generators.TRUTH["logacold3"],
            generators.TRUTH["logacold4"],
            generators.TRUTH["logacold6"],
            generators.TRUTH["logacold7"],
            generators.TRUTH["logawarm1"],
            generators.TRUTH["logawarm2"],
            generators.TRUTH["logawarm5"],
            generators.TRUTH["logawarm7"],
            generators.TRUTH["Tcold0"],
            generators.TRUTH["Tcold1"],
            generators.TRUTH["Twarm0"],
            generators.TRUTH["Twarm1"],
        )
        legacy_flux = legacy_model.modelFlux

        model = build_model()
        result = model(**generators.TRUTH).single().values
        np.testing.assert_allclose(result, legacy_flux, rtol=1e-10)


class TestDataLoading:
    """The real spectrum: legacy's own 25-120 micron selection and uncertainty rule."""

    def test_selects_25_to_120_micron(self) -> None:
        spectrum = generators.load_observed_spectrum()
        wavelength = spectrum.spectral_axis.values
        assert wavelength.min() >= 25.0
        assert wavelength.max() <= 120.0
        assert wavelength.size > 0

    def test_uncertainty_is_five_percent_of_the_flux(self) -> None:
        spectrum = generators.load_observed_spectrum()
        np.testing.assert_allclose(spectrum.uncertainty, 0.05 * np.abs(spectrum.values))


class TestTheProblemBuilds:
    """One model, one dataset; the fifteen model parameters plus calibration and GP."""

    def test_default_problem_has_the_calibration_and_gp_parameters(self) -> None:
        problem = build_problem()
        assert problem.backend == "reference"
        names = set(problem.parameters.free_names)
        for name in generators.TRUTH:
            assert f"model.{name}" in names
        assert "iso.instrument.calibration_scale.scale" in names
        assert any("likelihood" in name for name in names)
        assert len(names) == 18

    def test_no_gp_problem_drops_the_gp_parameters(self) -> None:
        problem = build_problem(gp=False)
        names = problem.parameters.free_names
        assert not any("likelihood" in name for name in names)
        assert len(names) == 16

    def test_build_instrument_is_labelled_iso(self) -> None:
        observed = generators.load_observed_spectrum()
        instrument = build_instrument(observed.spectral_axis.values)
        assert instrument.label == "iso"


class TestTheFitRuns:
    """A tiny-budget emcee fit, end to end (ruling 6)."""

    @staticmethod
    @pytest.fixture(scope="class")
    def tiny_run():
        problem = build_problem(synthetic=True, gp=False)
        return fit(problem, engine="emcee", **TINY)

    def test_fit_returns_the_expected_parameters(self, tiny_run) -> None:
        posterior = tiny_run["posterior"].dataset
        assert set(QUALIFIED_TRUTH).issubset(set(posterior.data_vars))
        assert posterior.sizes["chain"] == TINY["walkers"]
        assert posterior.sizes["draw"] == TINY["steps"] - TINY["burn_in"]

    def test_report_names_every_qualified_parameter(self, tiny_run) -> None:
        text = report(tiny_run)
        assert "emcee on reference" in text
        for name in QUALIFIED_TRUTH:
            assert name in text

    def test_recovers_truth_reports_one_bool_per_qualified_parameter(self, tiny_run) -> None:
        covered = recovers_truth(tiny_run)
        assert set(covered) == set(QUALIFIED_TRUTH)
        assert all(isinstance(value, bool) for value in covered.values())


class TestMain:
    """The CLI entry point, synthetic and at a tiny budget."""

    def test_main_synthetic_quick_prints_a_report(self, capsys: pytest.CaptureFixture[str]) -> None:
        assert (
            main(
                [
                    "--synthetic",
                    "--no-gp",
                    "--quick",
                    "--walkers",
                    str(TINY["walkers"]),
                    "--steps",
                    str(TINY["steps"]),
                    "--burn-in",
                    str(TINY["burn_in"]),
                ]
            )
            == 0
        )
        out = capsys.readouterr().out
        assert "free parameters:" in out
        assert "emcee on reference" in out
        assert "wall clock" in out
