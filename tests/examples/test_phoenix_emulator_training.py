"""The PHOENIX emulator's training logic, on a toy grid -- no download (W6.13 (4)).

``examples/phoenix_star/train_emulator.py`` runs once, by hand, against 157
files from Göttingen; these rows test its pure functions -- the flux-conserving
binning, the line-spread function, the PCA truncation rule and the GP fit --
on a synthetic grid of twelve spectra, so the script's logic is checked in CI
without the download.
"""

from __future__ import annotations

import numpy as np

from examples.phoenix_star import train_emulator as te


def _toy_grid(n: int = 12, pixels: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Twelve ``log10`` spectra spanned by exactly three smooth shapes."""
    rng = np.random.default_rng(3)
    nodes = rng.uniform(-1.0, 1.0, size=(n, 3))
    x = np.linspace(0.0, 1.0, pixels)
    shapes = np.array([np.sin(3 * x), np.cos(5 * x), x**2])
    coefficients = np.column_stack(
        [nodes[:, 0] + 0.3 * nodes[:, 1], np.tanh(nodes[:, 1]), 0.5 * nodes[:, 2] ** 2]
    )
    return nodes, -4.0 + 0.1 * coefficients @ shapes


class TestTheGrids:
    def test_segment_b_is_the_legacy_specwaves(self) -> None:
        grid = te.rvs_grid()
        assert grid.size == 387
        assert grid[0] == 0.842 and grid[-2] < 0.872 <= grid[-1]
        assert np.allclose(grid[1:] / grid[:-1], 1 + 1 / 11000)

    def test_segment_a_and_the_nodes(self) -> None:
        assert te.SEGMENT_A.size == 582
        nodes = te.grid_nodes()
        assert nodes.shape == (156, 3)
        assert len({tuple(row) for row in nodes}) == 156

    def test_urls_follow_the_goettingen_layout(self) -> None:
        assert te.spectrum_url(7200, 4.0, 0.5).endswith(
            "Z+0.5/lte07200-4.00+0.5.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
        )
        assert "Z-0.0/lte05000-4.50-0.0" in te.spectrum_url(5000, 4.5, 0.0)


class TestTheBinning:
    def test_bin_mean_is_exact_for_a_linear_spectrum(self) -> None:
        wave = np.linspace(1.0, 10.0, 5001)
        flux = 2.0 * wave + 1.0
        centres = np.geomspace(2.0, 8.0, 30)
        edges = te.bin_edges(centres)
        means = te.bin_mean(wave, flux, edges)
        assert np.allclose(means, (edges[1:] + edges[:-1]) + 1.0, rtol=1e-12)

    def test_bin_mean_conserves_the_total_integral(self) -> None:
        wave = np.linspace(1.0, 10.0, 20001)
        flux = 1.0 + np.sin(40 * wave) ** 2
        edges = te.bin_edges(np.geomspace(2.0, 8.0, 40))
        total = np.sum(te.bin_mean(wave, flux, edges) * np.diff(edges))
        inside = (wave >= edges[0]) & (wave <= edges[-1])
        exact = np.trapezoid(flux[inside], wave[inside])
        assert abs(total / exact - 1.0) < 1e-4


class TestTheLSF:
    def test_a_constant_stays_constant(self) -> None:
        wave = np.linspace(0.8, 0.9, 200001)
        out = te.lsf_sample(wave, np.full(wave.size, 3.0), te.rvs_grid())
        assert np.allclose(out, 3.0, rtol=1e-12)

    def test_a_narrow_line_takes_the_resolution_width_and_keeps_its_area(self) -> None:
        wave = np.linspace(0.8, 0.9, 400001)
        centre, width = 0.857, 0.857 / 500000
        flux = 1.0 - 0.5 * np.exp(-0.5 * ((wave - centre) / width) ** 2)
        targets = np.linspace(0.855, 0.859, 4001)
        out = te.lsf_sample(wave, flux, targets)
        depth = 1.0 - out
        # Equivalent width is conserved by a normalised kernel.
        ew_in = np.trapezoid(1.0 - flux, wave)
        ew_out = np.trapezoid(depth, targets)
        assert abs(ew_out / ew_in - 1.0) < 0.01
        # And the FWHM is about lambda / 11000 (the intrinsic width is ~2 % of it).
        half = depth >= depth.max() / 2
        fwhm = targets[half][-1] - targets[half][0]
        assert abs(fwhm / (centre / 11000) - 1.0) < 0.03


class TestThePCA:
    def test_the_smallest_sufficient_k_is_chosen(self) -> None:
        _, log_flux = _toy_grid()
        k, mean, basis, weights, errors = te.choose_components(log_flux, 120)
        assert k == 3
        rebuilt = mean + weights @ basis
        assert np.allclose(rebuilt, log_flux, atol=1e-10)
        assert errors["rms_a"] <= te.TARGET_A and errors["rms_b"] <= te.TARGET_B

    def test_the_cap_binds(self) -> None:
        _, log_flux = _toy_grid()
        k, _, basis, weights, _ = te.choose_components(log_flux, 120, cap=2)
        assert k == 2 and basis.shape == (2, 200) and weights.shape == (12, 2)


class TestTheGP:
    def test_the_fit_interpolates_the_nodes(self) -> None:
        nodes, log_flux = _toy_grid()
        _, _, _, weights, _ = te.choose_components(log_flux, 120)
        hyper, alpha = te.fit_all(nodes, weights)
        assert hyper.shape == (3, 5) and alpha.shape == (12, 3)
        assert np.all(hyper[:, 4] >= te.JITTER_FLOOR * (1 - 1e-9))
        predicted = te.gp_mean(nodes, nodes, alpha, hyper)
        # Exact up to the jitter (capped at 1e-6 of the amplitude squared, so a
        # residual of order amplitude * 1e-3): twelve scattered toy nodes are
        # a far sparser design than the real 156-node box.
        bound = 3.0 * hyper[:, 3] * np.sqrt(te.JITTER_CAP)
        assert np.all(np.abs(predicted - weights) <= bound)

    def test_the_predictive_mean_is_the_closed_form(self) -> None:
        rng = np.random.default_rng(0)
        nodes = rng.normal(size=(5, 3))
        alpha = rng.normal(size=(5, 1))
        hyper = np.array([[0.7, 1.3, 2.0, 1.5, 1e-8]])
        x = rng.normal(size=(2, 3))
        d2 = (((x[:, None, :] - nodes[None]) / hyper[0, :3]) ** 2).sum(-1)
        expected = (1.5**2 * np.exp(-0.5 * d2)) @ alpha
        assert np.allclose(te.gp_mean(x, nodes, alpha, hyper), expected, rtol=1e-14)
