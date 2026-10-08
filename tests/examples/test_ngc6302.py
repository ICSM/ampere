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
``TestSolverAgreement`` and ``TestExactOrderedPrior`` are W6.13 (2) tranche
B's two rulings: the O(N) solver agrees with the O(N^3) one, and the
triangle-plus-fraction reparameterisation is the same flat ordered prior
legacy's rejection-sampled box is, exactly.

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

from collections.abc import Iterator
from typing import Any

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import DenseGP, QuasisepGP
from examples.ngc6302 import dust_mass, generators
from examples.ngc6302.ngc6302 import (
    DEFAULT_GRID,
    DERIVED_TRUTH,
    QUALIFIED_TRUTH,
    build_instrument,
    build_model,
    build_problem,
    derived_temperatures,
    fit,
    main,
    recovers_truth,
    report,
)

TINY = {"walkers": 32, "steps": 10, "burn_in": 3}

# The "old" (pre-reparameterisation) synthetic spectrum's own summary
# numbers -- computed once from the branch's first commit
# (``generators.synthetic_data(build_model(), instrument, seed=generators.SEED)``,
# the four independent, unordered temperature boxes) and quoted here so
# ``TestExactOrderedPrior.test_synthetic_spectrum_unchanged_by_the_reparameterisation``
# does not depend on checking out that commit. See ruling 2: the
# reparameterisation changes what is *declared* (the prior), not the
# physical temperatures :data:`generators.TRUTH` encodes, so the synthetic
# spectrum it produces is unchanged.
_OLD_SYNTHETIC_SPECTRUM = {
    "mean": 708.7160959988712,
    "std": 189.98174311949063,
    "first5": [
        320.23823679712285,
        316.58977914448604,
        353.01993692987776,
        318.8355657174629,
        324.29206996207455,
    ],
    "last5": [
        424.19577285548036,
        467.06141415498763,
        445.32801798747,
        430.539285335882,
        426.3037219331683,
    ],
    "sum": 442947.55999929446,
}

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


@pytest.fixture(scope="module", autouse=True)
def _from_the_prior(prior_start: Any) -> Iterator[None]:
    """Every emcee and zeus run in this module starts from the prior (W7.12).

    The default start's one-start optimiser costs about four minutes on the
    twin's sixteen parameters, per fit; these rows test the twin, not the start.
    """
    with prior_start():
        yield


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
        # generators.TRUTH stores Tcold0/Tcold_fraction and Twarm0/Twarm_fraction
        # (ruling 2), not the physical Tcold1/Twarm1 legacy's own constructor
        # takes directly; derived_temperatures recovers them exactly.
        derived = derived_temperatures(generators.TRUTH)
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
            float(derived["Tcold1"]),
            generators.TRUTH["Twarm0"],
            float(derived["Twarm1"]),
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
        # Plus the two derived temperatures, Tcold1/Twarm1 (ruling 2).
        assert set(covered) == set(QUALIFIED_TRUTH) | set(DERIVED_TRUTH)
        assert all(isinstance(value, bool) for value in covered.values())

    def test_dust_masses_runs_on_the_fit_posterior(self, tiny_run) -> None:
        table = dust_mass.dust_masses(tiny_run)
        assert set(table["cold"]) == set(dust_mass.SPECIES_NAMES)
        assert set(table["warm"]) == set(dust_mass.SPECIES_NAMES)
        for value in table["cold"]["am. oliv."].values():
            assert np.isfinite(value)


class TestDustMass:
    """``dust_masses_at`` reproduces the legacy dust-mass script's printout.

    ``dust_masses_at`` keeps taking the physical temperatures directly
    (ruling 2), so *theta* here merges :data:`generators.TRUTH` (which stores
    ``Tcold_fraction``/``Twarm_fraction``, not ``Tcold1``/``Twarm1``) with
    ``derived_temperatures``'s output -- the legacy-printout row itself is
    untouched.
    """

    @staticmethod
    def _theta() -> dict[str, float]:
        return {**generators.TRUTH, **derived_temperatures(generators.TRUTH)}

    def test_matches_the_legacy_printout(self) -> None:
        table = dust_mass.dust_masses_at(self._theta())
        for component in ("cold", "warm"):
            for name, expected in _LEGACY_DUST_MASSES[component].items():
                if expected == 0.0:
                    assert table[component][name] == 0.0
                else:
                    np.testing.assert_allclose(table[component][name], expected, rtol=1e-6)

    def test_format_table_renders_every_species(self) -> None:
        table = dust_mass.dust_masses_at(self._theta())
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


class TestExactOrderedPrior:
    """Ruling 2: the triangle-plus-fraction reparameterisation is legacy's
    own flat ordered prior, exactly -- not an approximation of it."""

    def test_density_is_constant_on_the_cold_triangle(self) -> None:
        rng = np.random.default_rng(20260929)
        t0 = rng.uniform(10.0, 80.0, size=1000)
        f = rng.uniform(0.0, 1.0, size=1000)
        const = (
            st.triang(c=0, loc=10.0, scale=70.0).logpdf(t0)
            + st.uniform(0.0, 1.0).logpdf(f)
            - np.log(80.0 - t0)
        )
        np.testing.assert_allclose(const, const[0], atol=1e-10)

    def test_density_is_constant_on_the_warm_triangle(self) -> None:
        rng = np.random.default_rng(20260929)
        t0 = rng.uniform(80.0, 180.0, size=1000)
        f = rng.uniform(0.0, 1.0, size=1000)
        const = (
            st.triang(c=0, loc=80.0, scale=100.0).logpdf(t0)
            + st.uniform(0.0, 1.0).logpdf(f)
            - np.log(180.0 - t0)
        )
        np.testing.assert_allclose(const, const[0], atol=1e-10)

    def test_derived_temperatures_recovers_the_2002_solution_exactly(self) -> None:
        derived = derived_temperatures(generators.TRUTH)
        np.testing.assert_allclose(derived["Tcold1"], 57.03042, rtol=0, atol=1e-10)
        np.testing.assert_allclose(derived["Twarm1"], 122.69678, rtol=0, atol=1e-10)

    def test_derived_temperatures_handles_arrays(self) -> None:
        draws = {
            "Tcold0": np.full(5, generators.TRUTH["Tcold0"]),
            "Tcold_fraction": np.full(5, generators.TRUTH["Tcold_fraction"]),
            "Twarm0": np.full(5, generators.TRUTH["Twarm0"]),
            "Twarm_fraction": np.full(5, generators.TRUTH["Twarm_fraction"]),
        }
        derived = derived_temperatures(draws)
        assert derived["Tcold1"].shape == (5,)
        assert derived["Twarm1"].shape == (5,)
        np.testing.assert_allclose(derived["Tcold1"], DERIVED_TRUTH["Tcold1"])
        np.testing.assert_allclose(derived["Twarm1"], DERIVED_TRUTH["Twarm1"])

    def test_synthetic_spectrum_unchanged_by_the_reparameterisation(self) -> None:
        """generators.TRUTH's Tcold_fraction/Twarm_fraction were chosen to
        round-trip to the exact 2002-solution physical temperatures, so the
        synthetic spectrum built from the new (declared) parameterisation is
        the same array the old (independent-box) parameterisation produced,
        for the same seed (ruling 2)."""
        model = build_model()
        observed = generators.load_observed_spectrum()
        instrument = build_instrument(observed.spectral_axis.values)
        synthetic = generators.synthetic_data(model, instrument, seed=generators.SEED)
        values = synthetic.values

        assert values.size == 625
        np.testing.assert_allclose(values.mean(), _OLD_SYNTHETIC_SPECTRUM["mean"], rtol=1e-12)
        np.testing.assert_allclose(values.std(), _OLD_SYNTHETIC_SPECTRUM["std"], rtol=1e-12)
        np.testing.assert_allclose(values.sum(), _OLD_SYNTHETIC_SPECTRUM["sum"], rtol=1e-12)
        np.testing.assert_allclose(values[:5], _OLD_SYNTHETIC_SPECTRUM["first5"], rtol=1e-12)
        np.testing.assert_allclose(values[-5:], _OLD_SYNTHETIC_SPECTRUM["last5"], rtol=1e-12)


#: The regression row's pinned seed and its margin below the truth (W7.6).
BOUND_AWARE_SEED = 20261007
BOUND_AWARE_MARGIN = 100.0
#: Coordinates at a bound after the pre-W7.6 scipy route (every coordinate in
#: ``u``) on the same problem and seed, measured at ``429960c``.
OLD_ROUTE_SATURATED = 11


@pytest.mark.heavy
class TestTheBoundAwareScipyRoute:
    """W7.6's regression row: the scipy route on fifteen boxes and a half-line.

    The GP-off synthetic problem (sixteen free parameters: eleven abundance
    boxes ``uniform(-6, 0)``, two fractions ``uniform(0, 1)``, two
    ``triang(c=0)`` temperature boxes, and the ``lognorm`` calibration's
    ``Log`` floor), optimised by ``optimise(method="scipy", starts=1,
    seed=20261007)``. Measured on the same problem and seed with two starts
    (the row runs one: the first start alone reaches the same MAP, and two
    took 610 s — too long for CI's examples job; review, 2026-10-07):

    * the **old route** (every coordinate in ``u``, ``429960c``): log
      posterior ``-4146.0``, eleven coordinates at a bound (seven
      abundances at ``-6``, ``Tcold0 = 10``, ``Tcold_fraction = 1``,
      ``Twarm0 = 80``, ``Twarm_fraction = 0``), 19 691 evaluations;
    * the **new route** (Powell in the bounded coordinates' normalised
      constrained values, reflected into the box): log posterior ``-3151.8``,
      one coordinate at a bound (``logacold2 = -6``), 32 545 evaluations,
      both starts ending at Powell's default ``maxfev`` (so ``converged`` is
      false) — the second start at ``-3203.5``;
    * the truth the data were drawn from: ``-3106.9``.

    On seed 2 the old route gave ``-4947.7`` with twelve coordinates at a
    bound and the new one ``-3173.7`` with three, so the margin of 100 nats
    is pinned over the measured gaps of 45 and 67, against the old route's
    1039 and 1841.
    """

    def test_the_map_is_near_the_truth_and_off_the_corner(self) -> None:
        import warnings

        from ampere.inference import BoundSaturationWarning, optimise, saturated_bounds

        problem = build_problem(synthetic=True, gp=False)
        truth = problem.parameters.pack(QUALIFIED_TRUTH)
        with warnings.catch_warnings():
            # the count is asserted below; whether Twarm0 lands on its bound is not
            warnings.simplefilter("ignore", BoundSaturationWarning)
            optimum = optimise(problem, method="scipy", starts=1, seed=BOUND_AWARE_SEED)
        assert set(optimum.coordinates) == {"constrained"}
        assert optimum.log_prob_constrained >= problem.log_prob(truth) - BOUND_AWARE_MARGIN
        saturated = saturated_bounds(problem, optimum.unconstrained)
        print("saturated:", saturated, "log p:", optimum.log_prob_constrained)
        assert len(saturated) < OLD_ROUTE_SATURATED
