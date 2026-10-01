"""``examples/cstar`` loads its data, declares its model and, with Hyperion, runs (W6.13 (5)).

Two tiers. The rows that need nothing beyond ``dev`` always run: the tracked
votable loads as eleven library filters, the IRS CSV as per-chunk spectra
that are strictly increasing, the model declares the legacy's boxes, the
in-memory size recipe equals the tracked ``.size`` files, and ``--engine
sbi`` is refused by name without the ``sbi`` extra.

The ``needs_hyperion`` rows need the pixi ``hyperion`` environment
(``pixi run -e hyperion pytest tests/examples/test_cstar.py``) and skip
everywhere else, CI included (D12 (b): Hyperion is an example-only
requirement, not an extra). They check miepython's opacities (the
refractive-index sign, the SiC feature, the pure-species end points of the
mixture), build the Hyperion dust, run one ``photons="quick"`` simulation,
pool four through a ``ProcessExecutor`` and fit a 40-simulation NPE end to
end. The coverage run is not here: it is hours of CPU, run once and
recorded in :mod:`examples.cstar.cstar`'s docstring.
"""

from __future__ import annotations

import importlib.util
import time
from pathlib import Path

import numpy as np
import pytest

from examples.cstar import dust, generators
from examples.cstar.cstar import (
    EMBEDDING,
    GRID,
    PRIORS,
    SimulatorFailed,
    build_instruments,
    build_model,
    build_problem,
    fit,
)

HAS_SBI = importlib.util.find_spec("sbi") is not None
HAS_HYPERION = (
    importlib.util.find_spec("hyperion") is not None
    and importlib.util.find_spec("miepython") is not None
)
needs_hyperion = pytest.mark.skipif(
    not HAS_HYPERION,
    reason="needs hyperion and miepython (pixi run -e hyperion ...), an example-only requirement",
)
needs_miepython = pytest.mark.skipif(
    importlib.util.find_spec("miepython") is None, reason="needs miepython (pixi -e hyperion)"
)
needs_sbi = pytest.mark.skipif(not HAS_SBI, reason="needs the 'sbi' extra")

LEGACY_FILTERS = [
    "MCPS_B",
    "MCPS_V",
    "MCPS_I",
    "2MASS_J",
    "2MASS_H",
    "2MASS_Ks",
    "SPITZER_IRAC_36",
    "SPITZER_IRAC_45",
    "SPITZER_IRAC_58",
    "SPITZER_IRAC_80",
    "SPITZER_MIPS_24",
]

#: The legacy ``limits`` (``examples/cstar_model_test_sbi_v2.py`` lines 8-10),
#: minus the redundant second abundance.
LEGACY_BOXES = {
    "envelope_mass": (-10.0, -6.0),
    "envelope_rin": (-2.0, 2.0),
    "envelope_rout": (2.0, 4.0),
    "envelope_r0": (-2.0, 2.0),
    "stellar_mass": (1.0, 3.0),
    "stellar_luminosity": (1.0e3, 1.0e4),
    "sic_fraction": (0.0, 1.0),
}


class TestTheData:
    def test_the_votable_is_eleven_library_filters_in_jy(self) -> None:
        names, flux, error = generators.read_votable()
        assert names == LEGACY_FILTERS
        assert flux.size == error.size == 11
        # The file declares mJy; IRAC 8.0 is 66.98 mJy.
        assert flux[names.index("SPITZER_IRAC_80")] == pytest.approx(0.06698)
        instruments = build_instruments(names, generators.read_irs())
        step = instruments["photometry"].steps[0]
        assert list(step.filters) == LEGACY_FILTERS

    def test_the_irs_csv_is_one_strictly_increasing_spectrum_per_chunk(self) -> None:
        chunks = generators.read_irs()
        assert sorted(chunks) == [1, 2]
        assert sum(w.size for w, _, _ in chunks.values()) == 364
        for wavelength, flux, error in chunks.values():
            assert np.all(np.diff(wavelength) > 0.0)
            assert flux.size == error.size == wavelength.size
            assert np.all(error > 0.0)

    def test_the_photosphere_is_read_as_the_legacy_reads_it(self) -> None:
        nu, _fnu = generators.read_photosphere()
        raw = np.loadtxt(generators.PHOTOSPHERE_FILE, delimiter=",")
        assert nu.size == raw.shape[0] - 1  # skiprows=1 drops a data row, as the legacy does
        assert np.array_equal(nu, raw[1:, 0])

    def test_the_problem_composes_without_hyperion(self) -> None:
        problem = build_problem()
        assert problem.backend == "reference"
        names = problem.parameters.free_names
        assert {f"model.{name}" for name in PRIORS} <= set(names)
        assert {"irs_1.instrument.calibration_scale.scale"} <= set(names)
        assert not any("likelihood" in n for n in build_problem(gp=False).parameters.free_names)


class TestTheModel:
    def test_seven_parameters_with_the_legacy_boxes(self) -> None:
        model = build_model()
        assert tuple(model.parameters.free_names) == tuple(LEGACY_BOXES)
        for name, (low, high) in LEGACY_BOXES.items():
            prior = PRIORS[name]
            assert prior.ppf(0.0) == pytest.approx(low)
            assert prior.ppf(1.0) == pytest.approx(high)

    def test_the_grid_is_the_legacy_one(self) -> None:
        assert np.allclose(GRID, np.logspace(np.log10(0.2), np.log10(200), 1000))

    def test_the_embedding_dict_is_the_legacy_one(self) -> None:
        assert EMBEDDING == {"type": "FC", "num_hiddens": 100, "n_layers": 3, "output_dim": 28}

    def test_the_failure_class_is_declared(self) -> None:
        assert SimulatorFailed in build_problem()._failure_types  # the declared catch set


class TestTheDust:
    @pytest.mark.parametrize("name", ["rouleau91_ac", "SiC_Pegourie1988", "zubko96_ac_acar"])
    def test_the_size_recipe_equals_the_tracked_size_files(self, name: str) -> None:
        """Row (i): the legacy ``.size`` recipe, in memory, to 1e-6."""
        sizes, numbers = dust.size_distribution()
        tracked = np.loadtxt(dust.DATA_DIR / f"{name}.size")
        assert tracked.shape == (101, 2)
        assert np.allclose(sizes, tracked[:, 0], rtol=1e-6, atol=0.0)
        assert np.allclose(numbers, tracked[:, 1], rtol=1e-6, atol=0.0)

    @needs_miepython
    def test_the_sign_convention_keeps_absorbers_absorbing(self) -> None:
        """Row (ii): the SiC feature, and amorphous carbon's falling, absorbing opacity."""
        carbon, sic = dust.species_tables()
        wav = carbon["wavelength"]

        def at(table: dict[str, np.ndarray], lam: float, key: str = "kext") -> float:
            return float(np.interp(np.log(lam), np.log(wav), table[key]))

        assert at(sic, 11.3) > 1.3 * at(sic, 9.5)
        assert at(carbon, 1.0) > at(carbon, 10.0) > at(carbon, 100.0)
        assert at(carbon, 1.0, "ksca") / at(carbon, 1.0) < 0.5

    @needs_miepython
    def test_the_mixture_end_points_are_the_pure_species(self) -> None:
        """Row (iii)."""
        carbon, sic = dust.species_tables()
        for fraction, pure in ((0.0, carbon), (1.0, sic)):
            mixed = dust.mixture(fraction)
            assert np.allclose(mixed["chi"], pure["kext"], rtol=1e-12)
            assert np.allclose(mixed["albedo"], pure["ksca"] / pure["kext"], rtol=1e-12)
            assert np.allclose(mixed["g"], pure["g"], rtol=1e-12)

    @needs_hyperion
    def test_the_hyperion_dust_builds_with_finite_emissivities(self) -> None:
        """Row (iv)."""
        built = dust.build_dust(0.1)
        assert np.all(np.isfinite(built.emissivities.jnu))
        assert np.all(np.diff(built.optical_properties.nu) > 0.0)


class TestTheEngine:
    @pytest.mark.skipif(HAS_SBI, reason="the refusal is what an environment without sbi sees")
    def test_sbi_is_refused_by_name_without_the_extra(self) -> None:
        with pytest.raises(SystemExit, match="needs the 'sbi' extra"):
            fit(build_problem(gp=False), simulations=4, rounds=1, draws=4, cache=None)

    def test_other_engines_are_refused(self) -> None:
        with pytest.raises(SystemExit, match="--engine sbi only"):
            fit(build_problem(gp=False), engine="emcee", cache=None)


@needs_hyperion
class TestTheSimulator:
    def test_one_quick_simulation_is_finite_on_the_legacy_grid(self) -> None:
        model = build_model(photons="quick")
        started = time.perf_counter()
        result = model(**generators.SYNTHETIC_TRUTH)
        elapsed = time.perf_counter() - started
        print(f"\none photons='quick' simulation: {elapsed:.1f} s")
        flux = result.single().values
        assert flux.size == GRID.size == 1000
        assert np.all(np.isfinite(flux))
        assert np.all(flux >= 0.0)

    def test_the_composition_check_passes(self) -> None:
        build_problem(photons="quick").validate()

    def test_a_pooled_budget_runs_without_failures(self) -> None:
        from ampere.core import ProcessExecutor

        problem = build_problem(photons="quick", gp=False)
        batch = problem.simulate_many(4, observe=True, executor=ProcessExecutor(2))
        assert len(batch) == 4
        assert int(batch.failed.sum()) == 0
        assert sum(problem.failure_counts.values()) == 0

    @needs_sbi
    def test_a_forty_simulation_fit_runs_and_its_cache_serves_a_rerun(self, tmp_path: Path) -> None:
        problem = build_problem(photons="quick", gp=False)
        run = fit(
            problem,
            simulations=40,
            rounds=1,
            draws=8,  # each stored draw is scored by one serial Hyperion run
            workers=4,
            cache=tmp_path,
            training={"max_num_epochs": 5},
        )
        posterior = run["posterior"].dataset
        assert {f"model.{name}" for name in PRIORS} <= set(posterior.data_vars)
        assert int(run.attrs["ampere_sbi_cache_hit"]) == 0
        # A rerun at the same settings restores the stored posterior: no simulation.
        again = fit(
            build_problem(photons="quick", gp=False),
            simulations=40,
            rounds=1,
            draws=8,  # each stored draw is scored by one serial Hyperion run
            workers=4,
            cache=tmp_path,
            training={"max_num_epochs": 5},
        )
        assert int(again.attrs["ampere_sbi_cache_hit"]) == 1
