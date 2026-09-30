"""``examples/phoenix_star`` loads its emulator, builds and fits (W6.13 (4)).

Fast and download-free: the committed ``phoenix_emulator.npz`` is checked for
size, provenance and the interpolation property on the node spectra it
stores; the module's CCM89 is checked against ``dust_extinction``; the problem
builds and a tiny emcee fit and ``main`` run. The ``sbi`` arm, the three-backend
agreement and a tiny NUTS run skip where their extra is absent. The
full-budget coverage runs are recorded in :mod:`examples.phoenix_star.phoenix_star`'s
docstring, not run here.
"""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from examples.phoenix_star import emulator, generators
from examples.phoenix_star.phoenix_star import (
    QUALIFIED_TRUTH,
    build_problem,
    fit,
    main,
    recovers_truth,
    report,
    star_class,
)

needs_sbi = pytest.mark.skipif(
    importlib.util.find_spec("sbi") is None, reason="needs the 'sbi' extra (pixi run -e sbi)"
)
needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs ampere[jax]")
needs_dust_extinction = pytest.mark.skipif(
    importlib.util.find_spec("dust_extinction") is None, reason="needs the 'extinction' extra"
)

TINY = {"walkers": 14, "steps": 20, "burn_in": 5}
THETA = {**generators.TRUTH, "teff": 6500.0, "logg": 4.2}


class TestTheCommittedEmulator:
    def test_the_file_is_under_a_megabyte_and_names_its_sources(self) -> None:
        assert emulator.EMULATOR_FILE.stat().st_size < 1_000_000
        data = emulator.load_emulator()
        assert len(data.provenance["files"]) == 156
        assert data.provenance["n_files"] == 157
        assert 1 <= data.components <= 16

    def test_the_segments_have_their_sizes(self) -> None:
        data = emulator.load_emulator()
        assert data.segment_a.size == 582
        assert data.segment_b.size == 387
        assert data.basis.shape == (data.components, 582 + 387)

    def test_the_emulator_reproduces_its_stored_node_spectra(self) -> None:
        # Against the true PHOENIX spectra the file stores for four nodes: the
        # smoothing GP (train_emulator's "The jitter") reproduces them to about
        # a per cent RMS, and never worse than the training run recorded.
        data = emulator.load_emulator()
        worst = data.provenance["node_reproduction"]["max_vs_truth"]
        model = emulator.PhoenixEmulator()
        for index, stored in zip(data.stored_index, data.stored_spectra, strict=True):
            teff, logg, feh = data.grid_parameters[index]
            result = model(teff=teff, logg=logg, feh=feh)
            produced = np.concatenate([result["sed"].values, result["rvs"].values])
            fractional = produced / 10.0**stored - 1.0
            assert np.sqrt(np.mean(fractional**2)) < 0.02
            assert np.max(np.abs(fractional)) <= worst * (1 + 1e-9)

    def test_the_shape_has_unit_bolometric_flux_on_segment_a(self) -> None:
        result = emulator.PhoenixEmulator()(teff=6000.0, logg=4.5, feh=0.0)
        sed = result["sed"]
        # Segment A spans 0.3-5.5 micron, most but not all of the flux (per A).
        fraction = np.trapezoid(sed.values, sed.spectral_axis.values * 1e4)
        assert 0.85 < fraction < 1.0


@needs_dust_extinction
def test_ccm89_matches_dust_extinction() -> None:
    from dust_extinction.parameter_averages import CCM89

    import astropy.units as u

    wavelength = np.geomspace(0.305, 3.3, 10)
    for r_v in (2.9, 3.1, 4.0):
        a, b = emulator.ccm89_ab(wavelength)
        expected = CCM89(Rv=r_v)(wavelength * u.um)
        assert np.allclose(a + b / r_v, expected, rtol=1e-6, atol=0.0)


class TestTheProblem:
    def test_the_problem_has_the_five_parameters_and_the_noise_ones(self) -> None:
        problem = build_problem()
        names = problem.parameters.free_names
        for name in QUALIFIED_TRUTH:
            assert name in names
        assert "rvs.instrument.calibration_scale.scale" in names
        assert sum("likelihood" in name for name in names) == 2
        assert len(names) == 8

    def test_feh_and_distance_are_buffers_until_promoted(self) -> None:
        model = star_class("reference")()
        assert set(model.buffers.names) >= {"feh", "distance"}
        model.promote_buffer("feh", prior=__import__("scipy.stats").stats.uniform(0.0, 0.5))
        assert "feh" in model.parameters.free_names

    def test_the_legacy_luminosity_scaling(self) -> None:
        near = star_class("reference")()(**THETA)
        far = star_class("reference")(distance=200.0)(**THETA)
        assert np.allclose(near["sed"].values / far["sed"].values, 4.0, rtol=1e-12)

    def test_a_tiny_emcee_fit_runs_and_reports(self) -> None:
        run = fit(build_problem(gp=False), engine="emcee", **TINY)
        assert set(QUALIFIED_TRUTH) <= set(run["posterior"].dataset.data_vars)
        assert set(recovers_truth(run)) == set(QUALIFIED_TRUTH)
        assert "emcee on reference" in report(run)

    def test_nuts_refuses_the_reference_backend(self) -> None:
        with pytest.raises(SystemExit, match="needs gradients"):
            fit(build_problem(gp=False), engine="nuts", backend="reference")

    def test_main_prints_a_report(self, capsys: pytest.CaptureFixture[str]) -> None:
        argv = ["--engine", "emcee", "--no-gp", "--walkers", "14", "--steps", "20"]
        assert main([*argv, "--burn-in", "5"]) == 0
        out = capsys.readouterr().out
        assert "emcee on reference" in out and "wall clock" in out


@needs_sbi
def test_a_small_npe_fit_runs() -> None:
    run = fit(build_problem(gp=False), engine="sbi", budget=300, draws=200, cache=False)
    assert run.attrs["ampere_engine"] == "sbi"
    assert set(QUALIFIED_TRUTH) <= set(run["posterior"].dataset.data_vars)


def _agreement(backend: str) -> None:
    reference = star_class("reference")()(**THETA)
    native = star_class(backend)()(**THETA)
    for channel in ("sed", "rvs"):
        assert np.allclose(native[channel].values, reference[channel].values, rtol=1e-8, atol=0)
    if backend == "torch":
        from examples.phoenix_star.emulator_torch import PhoenixEmulator
    else:
        from examples.phoenix_star.emulator_jax import PhoenixEmulator
    bare = emulator.PhoenixEmulator()(teff=6100.0, logg=4.4, feh=0.3)
    twin = PhoenixEmulator()(teff=6100.0, logg=4.4, feh=0.3)
    for channel in ("sed", "rvs"):
        assert np.allclose(twin[channel].values, bare[channel].values, rtol=1e-8, atol=0)


def _tiny_nuts(backend: str) -> None:
    problem = build_problem(backend, gp=False)
    run = fit(problem, engine="nuts", backend=backend, draws=10, warmup=5, chains=1)
    assert run.attrs["ampere_engine"] == "nuts"
    assert set(QUALIFIED_TRUTH) <= set(run["posterior"].dataset.data_vars)


@needs_torch
class TestTorch:
    def test_agrees_with_the_reference(self) -> None:
        import torch

        torch.set_num_threads(4)
        _agreement("torch")

    def test_the_flux_is_differentiable(self) -> None:
        import torch

        model = star_class("torch")()
        teff = torch.tensor(6500.0, dtype=torch.float64, requires_grad=True)
        flux = model.flux("sed", {**THETA, "teff": teff})
        flux.sum().backward()
        assert teff.grad is not None and torch.isfinite(teff.grad)

    def test_a_tiny_nuts_run(self) -> None:
        _tiny_nuts("torch")


@needs_jax
class TestJax:
    def test_agrees_with_the_reference(self) -> None:
        from examples.phoenix_star.phoenix_star import backend_module

        backend_module("jax")
        _agreement("jax")

    def test_a_tiny_nuts_run(self) -> None:
        _tiny_nuts("jax")
