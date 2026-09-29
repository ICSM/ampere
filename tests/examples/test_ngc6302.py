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
``TestSolverAgreement`` is ruling 1: the O(N) solver agrees with the O(N^3)
one.

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

from ampere.core import DenseGP, QuasisepGP
from examples.ngc6302 import dust_mass, generators
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

# The legacy dust-mass script's own printed numbers
# (`cd examples && pixi run -e dev python NGC6302-calculate-dust-mass.py`),
# quoted here so the reproduction test does not depend on that script
# continuing to run in CI -- it is the frozen legacy characterisation
# anchor, not touched by this item, and this suite must stay fast.
_LEGACY_DUST_MASSES = {
    "cold": {
        "am. oliv.": 0.04754212255388583,
        "forst.": 0.0005689300735648554,
        "calcite": 1.4866336145406103e-05,
        "ice": 1.4885199428994224e-05,
        "diopside": 3.019907040849098e-05,
        "c-enst.": 3.916360746071945e-06,
        "dolomite": 5.53800058979032e-06,
        "iron": 0.0,
    },
    "warm": {
        "am. oliv.": 6.8280769429978754e-06,
        "forst.": 8.66440471981681e-08,
        "calcite": 0.0,
        "ice": 0.0,
        "diopside": 0.0,
        "c-enst.": 8.767434461047079e-08,
        "dolomite": 0.0,
        "iron": 1.5101493748674502e-05,
    },
}


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

    def test_wavelength_is_sorted_and_strictly_increasing(self) -> None:
        """QuasisepGP's own precondition (ruling 1); the tracked file already
        satisfies it (no duplicate wavelength in the 25-120 micron window)."""
        spectrum = generators.load_observed_spectrum()
        wavelength = spectrum.spectral_axis.values
        assert wavelength.size == 625
        assert np.all(np.diff(wavelength) > 0)


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

    def test_dust_masses_runs_on_the_fit_posterior(self, tiny_run) -> None:
        table = dust_mass.dust_masses(tiny_run)
        assert set(table["cold"]) == set(dust_mass.SPECIES_NAMES)
        assert set(table["warm"]) == set(dust_mass.SPECIES_NAMES)
        for value in table["cold"]["am. oliv."].values():
            assert np.isfinite(value)


class TestDustMass:
    """``dust_masses_at`` reproduces the legacy dust-mass script's printout."""

    def test_matches_the_legacy_printout(self) -> None:
        table = dust_mass.dust_masses_at(generators.TRUTH)
        for component in ("cold", "warm"):
            for name, expected in _LEGACY_DUST_MASSES[component].items():
                if expected == 0.0:
                    assert table[component][name] == 0.0
                else:
                    np.testing.assert_allclose(table[component][name], expected, rtol=1e-6)

    def test_format_table_renders_every_species(self) -> None:
        table = dust_mass.dust_masses_at(generators.TRUTH)
        text = dust_mass.format_table(table)
        for name in dust_mass.SPECIES_NAMES:
            assert name in text
        assert "total cold" in text
        assert "total warm" in text


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


class TestSolverAgreement:
    """Ruling 1: ``QuasisepGP`` (O(N)) is now the default, and agrees with
    ``DenseGP`` (O(N^3)) exactly on this problem's Matern-3/2 likelihood."""

    def test_default_solver_is_quasisep(self) -> None:
        from examples.ngc6302.ngc6302 import _likelihood

        assert isinstance(_likelihood(gp=True).noise.solver, QuasisepGP)
        assert isinstance(_likelihood(gp=True, solver="dense").noise.solver, DenseGP)

    def test_quasisep_and_dense_score_the_same_log_prob_at_truth(self) -> None:
        dense = build_problem(synthetic=True, gp=True, solver="dense")
        quasisep = build_problem(synthetic=True, gp=True, solver="quasisep")
        assert dense.parameters.free_names == quasisep.parameters.free_names

        # generators.TRUTH plus the GP hyperparameters' prior median (no
        # injected truth for those -- module docstring).
        theta = dict(dense.reference_values)
        theta.update(QUALIFIED_TRUTH)

        lp_dense = float(dense.log_prob(theta))
        lp_quasisep = float(quasisep.log_prob(theta))
        assert lp_quasisep == pytest.approx(lp_dense, rel=1e-6)

    def test_quasisep_and_dense_agree_away_from_truth(self) -> None:
        """One point could agree by accident; a few draws from the priors cannot."""
        dense = build_problem(synthetic=True, gp=True, solver="dense", seed=7)
        quasisep = build_problem(synthetic=True, gp=True, solver="quasisep", seed=7)
        rng = np.random.default_rng(11)
        for _ in range(3):
            theta = dense.sample_prior(rng)
            lp_dense = float(dense.log_prob(theta))
            lp_quasisep = float(quasisep.log_prob(theta))
            if np.isfinite(lp_dense) or np.isfinite(lp_quasisep):
                assert lp_quasisep == pytest.approx(lp_dense, rel=1e-6)
