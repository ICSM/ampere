"""``examples/star_disc`` reads HD 105's data, builds and fits (W6.13 (6)).

Fast: the votable's nineteen points (seventeen library filters, two
top-hats), the dust term against a hand-written reference, the photosphere's
continuity where the emulator hands over to its extrapolation, the IRS reader
on the tracked CASSIS file of :mod:`examples.linear_sed`, and a tiny
``--synthetic`` emcee fit. The full-budget coverage run is recorded in
:mod:`examples.star_disc.star_disc`'s docstring, not run here.
"""

from __future__ import annotations

import numpy as np
import pytest
import scipy.stats as st

from examples.linear_sed.generators import IRS_FILE, deduplicated_grids
from examples.star_disc import generators
from examples.star_disc.star_disc import (
    PARALLAX,
    QUALIFIED_TRUTH,
    TOPHATS,
    StarDisc,
    build_problem,
    fit,
    main,
    photometry_step,
    recovers_truth,
)

TINY = {"walkers": 22, "steps": 12, "burn_in": 4}


class TestTheData:
    def test_the_votable_has_nineteen_points_seventeen_in_the_library(self) -> None:
        names, flux, error = generators.read_votable()
        assert len(names) == 19 and flux.shape == error.shape == (19,)
        assert set(TOPHATS) - {"APEX/LABOCA.870"} <= set(names)  # LABOCA is the CSV's
        step = photometry_step(names)
        assert list(step.detectors[-2:]) == ["energy", "energy"]
        assert len([n for n in names if n not in TOPHATS]) == 17

    def test_the_top_hats_have_the_stated_widths(self) -> None:
        (a_blue, a_red), (b_blue, b_red) = TOPHATS["ALMA/ALMA.B6"], TOPHATS["ATCA/ATCA.9mm"]
        assert np.allclose([a_blue, a_red], [1090.2, 1420.8], rtol=1e-3)
        assert np.allclose([b_blue, b_red], [7889.3, 9993.1], rtol=1e-3)
        l_blue, l_red = TOPHATS["APEX/LABOCA.870"]
        assert np.allclose([l_blue, l_red], [799.5, 951.7], rtol=1e-3)

    def test_the_csv_limits_are_the_two_flagged_rows(self) -> None:
        names, wavelength, flux, error = generators.read_limits()
        assert names == ["HERSCHEL_SPIRE_PLW", "APEX/LABOCA.870"]
        assert np.array_equal(wavelength, [500.0, 880.0])
        assert np.array_equal(flux, [2.1e-3, 1.07e-2]) and np.array_equal(error, [7e-3, 5.9e-3])
        vot_names, _, _ = generators.read_votable()
        step = photometry_step([*vot_names, *names])
        assert step.detectors[-1] == "energy"  # LABOCA, the top-hat
        assert "HERSCHEL_SPIRE_PLW" not in TOPHATS  # SPIRE 500 is the library's row

    def test_the_rvs_window_is_the_legacy_one(self) -> None:
        grid = generators.rvs_grid()
        assert grid[0] == 0.847 and grid[-2] < 0.871 <= grid[-1]

    def test_the_irs_reader_splits_and_deduplicates(self) -> None:
        (sl, sl_flux, sl_err), (ll, _, _) = generators.read_irs(IRS_FILE)
        expected_sl, expected_ll = deduplicated_grids()
        assert np.array_equal(sl, expected_sl) and np.array_equal(ll, expected_ll)
        assert sl_flux.shape == sl_err.shape == sl.shape


class TestTheModel:
    def test_the_dust_term_matches_a_hand_written_reference(self) -> None:
        import astropy.units as u
        from astropy.modeling.physical_models import BlackBody

        model = StarDisc()
        theta = {**generators.SYNTHETIC_TRUTH}
        ctx = model.context(theta)
        wavelength = model.plans["sed"].wavelength
        # QuickSED.py lines 158-164, verbatim apart from the names.
        freq = (wavelength * u.micron).to(u.Hz, equivalencies=u.spectral()).value
        emission = BlackBody().evaluate(freq, theta["t_dust"], 1 * u.Jy / u.sr)
        modified = np.where(wavelength >= theta["lambda_0"])
        emission[modified] = emission[modified] * (theta["lambda_0"] / wavelength[modified]) ** 1.0
        emission *= (
            np.pi
            * 1e23
            * ((10 ** theta["log_area"] * 1.495978707e11) / ((1000.0 / PARALLAX) * 3.0857e16)) ** 2
        )
        produced = model.dust("sed", ctx)
        # Relative agreement wherever the dust is not underflowing in the Wien tail.
        live = emission > 1e-200
        assert np.array_equal(produced[~live] <= 1e-200, np.ones((~live).sum(), bool))
        assert np.allclose(produced[live], emission[live], rtol=1e-10, atol=0.0)

    def test_the_photosphere_is_continuous_at_the_emulator_edge(self) -> None:
        edge = 5.5
        grid = np.array([0.5, edge * (1 - 1e-6), edge * (1 + 1e-6), 6.0])
        model = StarDisc(grid)
        star = model.star("sed", model.context(dict(generators.SYNTHETIC_TRUTH)))
        assert abs(star[2] / star[1] - 1.0) < 1e-4
        # Beyond the edge the extrapolation is F_nu ~ lambda**-2.
        assert np.isclose(star[3] / star[2], (6.0 / (edge * (1 + 1e-6))) ** -2, rtol=1e-12)


class TestTheFit:
    def test_the_problem_builds_in_both_modes(self) -> None:
        for synthetic in (False, True):
            problem = build_problem(synthetic=synthetic, gp=False)
            assert set(QUALIFIED_TRUTH) <= set(problem.parameters.free_names)
            censoring = problem.datasets["photometry"].likelihood.censoring
            assert censoring is not None and censoring.n_censored == 2
        bare = build_problem(gp=False, limits=False)
        assert bare.datasets["photometry"].likelihood.censoring is None

    def test_the_limits_enter_the_likelihood_as_tobit_terms(self) -> None:
        # likelihoods.md section 9: an upper limit contributes logcdf(z) with
        # z = (recorded - predicted) / sigma, not a Gaussian density.
        values = {**QUALIFIED_TRUTH, "rvs.instrument.calibration_scale.scale": 1.0}
        with_limits = build_problem(gp=False)
        without = build_problem(gp=False, limits=False)
        dataset = with_limits.datasets["photometry"]
        censoring = dataset.likelihood.censoring
        assert censoring is not None
        assert (censoring.n_samples, censoring.n_censored) == (21, 2)
        assert dataset.likelihood.to_spec()["censoring"] == {"n_samples": 21, "n_censored": 2}
        predicted = with_limits.simulate(values).predicted["photometry"].values[-2:]
        recorded = np.asarray(dataset.observed.values)[-2:]
        sigma = np.asarray(dataset.observed.uncertainty)[-2:]
        tobit = st.norm.logcdf((recorded - predicted) / sigma)
        expected = without.evaluate(values).contributions["photometry"] + tobit.sum()
        got = with_limits.evaluate(values).contributions["photometry"]
        assert got == pytest.approx(expected, rel=1e-9)
        assert got - without.evaluate(values).contributions["photometry"] == pytest.approx(
            tobit.sum(), abs=1e-6
        )
        # Dropping the declaration would score them as detections instead.
        assert (np.log(st.norm.pdf((recorded - predicted) / sigma) / sigma) != tobit).all()

    def test_a_tiny_synthetic_emcee_fit_runs(self) -> None:
        run = fit(build_problem(synthetic=True, gp=False), **TINY)
        assert set(recovers_truth(run)) == set(QUALIFIED_TRUTH)

    def test_main_synthetic_prints_a_report(self, capsys: pytest.CaptureFixture[str]) -> None:
        argv = ["--synthetic", "--no-gp", "--walkers", "22", "--steps", "12", "--burn-in", "4"]
        assert main(argv) == 0
        out = capsys.readouterr().out
        assert "emcee on reference" in out and "wall clock" in out
