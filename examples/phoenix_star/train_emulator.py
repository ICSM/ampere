"""Train the PHOENIX emulator committed beside this file -- a one-off, not CI.

``python -m examples.phoenix_star.train_emulator`` downloads a subset of the
PHOENIX-ACES-AGSS-COND-2011 HiRes grid (Husser et al. 2013) from Göttingen,
reduces every spectrum to two segments, compresses them with a PCA, fits one
Gaussian process per PCA weight over (Teff, log g, [Fe/H]), and writes
``phoenix_emulator.npz`` (well under a megabyte). The example
(:mod:`.phoenix_star`, and :mod:`examples.star_disc`) only ever reads that
file, so the download (157 files, about a gigabyte, into the astropy cache --
outside the repository, idempotent, so an interrupted run resumes) and the
training are an optional, documented step. This replaces the legacy script's
Starfish route (``examples/examples_paper/phoenixstar.py``, lines 84-110):
Starfish is not used anywhere in the twin.

The grid
--------
Teff 5000-7000 K in steps of 100 K plus 7200-8000 K in steps of 200 K (the
PHOENIX grid's own spacing), log g 4.0, 4.5 and 5.0, [Fe/H] 0.0 and +0.5:
26 x 3 x 2 = 156 spectra, plus the shared wavelength file.

The two segments
----------------
* **A, the SED**: :data:`SEGMENT_A`, ``np.geomspace(0.3, 5.5, 582)`` micron
  (R ~ 200 in log wavelength). Each pixel is the mean of ``F_lambda`` over its
  bin (edges at the geometric midpoints between pixel centres), computed from
  the cumulative trapezoid integral, so it conserves the flux of the
  piecewise-linear HiRes spectrum exactly -- flux conservation to far better
  than an emulator needs.
* **B, the Gaia RVS window**: :func:`rvs_grid`, built exactly the legacy way
  (``phoenixstar.py`` lines 222-227: start at 0.842 micron and multiply by
  ``1 + 1/11000`` until past 0.872) -- **387** points. The HiRes spectrum is
  convolved with a Gaussian line-spread function of FWHM ``lambda / 11000``
  (constant R, so a fixed width in ``ln lambda``) and sampled at those points.

Wavelengths are **vacuum** throughout: PHOENIX HiRes is tabulated in vacuum.
The legacy Starfish route converted to air; the twin does not, and because the
twin's synthetic observations come from the same emulator, nothing is
inconsistent.

Every spectrum is normalised to **unit bolometric flux** -- ``F_lambda``
divided by its trapezoid integral over the file's whole range (500 A to
5.5 micron), in the file's surface units (the legacy ``fbol_init``) -- so the
emulated shape has units of per angstrom, and the luminosity scaling happens in
the model. The emulated quantity is ``log10`` of that normalised flux on the
concatenation of segments A and B.

PCA and the regressor
---------------------
A PCA on the mean-subtracted ``log10`` grid (:func:`numpy.linalg.svd`), keeping
the smallest K that reconstructs every grid spectrum to 0.5 % RMS in flux on A
and 1 % on B, capped at 16 (:func:`choose_components`). Each PCA weight is then
regressed on standardised (Teff, log g, [Fe/H]) by a Gaussian process with an
ARD squared-exponential kernel (:func:`fit_gp`: three length scales in
0.05-5 standardised units, an amplitude, and a jitter between 1e-8 and 1e-4 of
the amplitude squared, by maximising the marginal likelihood).

**Why a GP and not the MLP the item prefers** (W6.13 lets the tranche choose):
it trains in seconds with scipy alone in the ``dev`` environment, and its
predictive mean ``k(x, X) alpha`` is a closed form in ``exp`` and ``matmul``
that reads identically in numpy, torch and jax. The emulator is therefore
differentiable on every backend by construction. An MLP needs a training loop
in a framework the ``dev`` environment does not have.

The jitter
----------
The dispatch expected the GP to interpolate every node exactly, up to a jitter
floored at 1e-8. On this grid the marginal likelihood does not agree. The PCA
weights vary with Teff on the scale of the 100 K grid step, and the log g axis
has only three nodes, so a noise-free GP needs length scales near the grid
spacing. It then reproduces the nodes (to 3e-6 in flux) but predicts poorly
between them: leave-one-out RMS 9 % on A and 5 % on B. With the jitter free
up to 1e-4, the likelihood uses it and the GP smooths rather than
interpolates. Leave-one-out improves to 0.6 % RMS on A and 0.3 % on B, but a
node is then reproduced only to the smoothing, several per cent at the worst
pixel (the ultraviolet end of A for the coolest stars). An emulator exists to
predict between nodes, so the committed file uses the smoothing GP. Its
provenance record and the W6.13 (4) report carry the numbers.

Checks
------
(i) leave-one-out on every seventh node (the GPs refitted without it, the PCA
basis kept), against the true PHOENIX spectrum; (ii) the emulator at every node
against its PCA reconstruction and against the true spectrum (the maximum and
RMS fractional error in flux). Both go into the file's provenance record and
are printed. The file also stores four nodes' true spectra, and a smoke row
(:mod:`tests.examples.test_phoenix_star`) checks the emulator against them.
"""

from __future__ import annotations

import argparse
import datetime
import json
import sys
from pathlib import Path

import numpy as np
import scipy.optimize

__all__ = [
    "BASE_URL",
    "FEH",
    "LOGG",
    "OUTPUT",
    "RVS_RESOLVING_POWER",
    "SEGMENT_A",
    "TEFF",
    "bin_edges",
    "bin_mean",
    "choose_components",
    "fit_gp",
    "gp_mean",
    "grid_nodes",
    "lsf_sample",
    "main",
    "rvs_grid",
    "spectrum_url",
]

BASE_URL = "https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/"
WAVE_URL = BASE_URL + "WAVE_PHOENIX-ACES-AGSS-COND-2011.fits"
TEFF = tuple(range(5000, 7001, 100)) + tuple(range(7200, 8001, 200))
LOGG = (4.0, 4.5, 5.0)
FEH = (0.0, 0.5)
RVS_RESOLVING_POWER = 11000.0
#: Segment A, micron.
SEGMENT_A = np.geomspace(0.3, 5.5, 582)
OUTPUT = Path(__file__).resolve().parent / "phoenix_emulator.npz"
#: The truncation targets: RMS fractional flux error per spectrum.
TARGET_A = 0.005
TARGET_B = 0.01
MAX_COMPONENTS = 16
JITTER_FLOOR = 1e-8
#: The jitter's ceiling (relative to amplitude squared). See "The jitter" in
#: the module docstring for why it is not held at the floor.
JITTER_CAP = 1e-4
#: The length scales' range, in standardised input units. The ceiling keeps
#: the kernel matrix conditioned (a length scale far beyond the box makes
#: columns nearly identical).
LENGTH_BOUNDS = (0.05, 5.0)
#: How many nodes' PHOENIX spectra the file stores for the smoke row.
STORED_NODES = 4


def rvs_grid() -> np.ndarray:
    """Segment B: the legacy ``specwaves`` exactly (``phoenixstar.py`` 222-227), micron."""
    waves = [0.842]
    while waves[-1] < 0.872:
        waves.append(waves[-1] * (1 + (1 / RVS_RESOLVING_POWER)))
    return np.array(waves)


def grid_nodes() -> np.ndarray:
    """The 156 (Teff, log g, [Fe/H]) nodes, in download order."""
    return np.array([(t, g, z) for z in FEH for g in LOGG for t in TEFF], dtype=float)


def spectrum_url(teff: float, logg: float, feh: float) -> str:
    """The Göttingen URL of one HiRes spectrum."""
    z = "-0.0" if feh == 0.0 else f"{feh:+.1f}"
    return (
        f"{BASE_URL}PHOENIX-ACES-AGSS-COND-2011/Z{z}/lte{int(teff):05d}-{logg:.2f}{z}"
        ".PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
    )


def bin_edges(centres: np.ndarray) -> np.ndarray:
    """Bin edges at the geometric midpoints of *centres*, outer edges mirrored."""
    centres = np.asarray(centres, dtype=float)
    mid = np.sqrt(centres[1:] * centres[:-1])
    first = centres[0] ** 2 / mid[0]
    last = centres[-1] ** 2 / mid[-1]
    return np.concatenate([[first], mid, [last]])


def bin_mean(wave: np.ndarray, flux: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """Mean of *flux* over each ``[edges[i], edges[i+1]]``, flux-conserving.

    The cumulative trapezoid integral of the piecewise-linear *flux*, read at
    the edges (clipped to the tabulated range), differenced and divided by the
    bin width -- exact for a piecewise-linear input.
    """
    wave = np.asarray(wave, dtype=float)
    flux = np.asarray(flux, dtype=float)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (flux[1:] + flux[:-1]) * np.diff(wave))])
    edges = np.clip(np.asarray(edges, dtype=float), wave[0], wave[-1])
    # The integral up to each edge, exactly: the whole samples before it plus
    # the partial trapezoid of the linear segment the edge falls in.
    i = np.clip(np.searchsorted(wave, edges, side="right") - 1, 0, wave.size - 2)
    width = wave[i + 1] - wave[i]
    dx = edges - wave[i]
    slope = (flux[i + 1] - flux[i]) / width
    integral = cumulative[i] + flux[i] * dx + 0.5 * slope * dx**2
    return np.diff(integral) / np.diff(edges)


def lsf_sample(
    wave: np.ndarray,
    flux: np.ndarray,
    targets: np.ndarray,
    resolving_power: float = RVS_RESOLVING_POWER,
    oversample: float = 20.0,
) -> np.ndarray:
    """Convolve with a Gaussian LSF of FWHM ``lambda / R`` and sample at *targets*.

    Constant R is a constant width in ``ln lambda``, so the spectrum is first
    interpolated onto a uniform ``ln lambda`` grid (*oversample* points per
    FWHM, each the bin mean) spanning *targets* plus five FWHM either side,
    convolved with the
    normalised Gaussian there, and read at *targets*.
    """
    targets = np.asarray(targets, dtype=float)
    fwhm = 1.0 / resolving_power
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    step = fwhm / oversample
    lo = np.log(targets[0]) - 5.0 * fwhm
    hi = np.log(targets[-1]) + 5.0 * fwhm
    log_grid = np.arange(lo, hi + step, step)
    # Bin-average (not point-sample) onto the log grid, so HiRes structure
    # finer than the step is conserved rather than aliased.
    edges = np.exp(np.concatenate([log_grid - 0.5 * step, [log_grid[-1] + 0.5 * step]]))
    fine = bin_mean(np.asarray(wave, dtype=float), np.asarray(flux, dtype=float), edges)
    half = int(np.ceil(5.0 * sigma / step))
    offsets = np.arange(-half, half + 1) * step
    kernel = np.exp(-0.5 * (offsets / sigma) ** 2)
    kernel /= kernel.sum()
    smooth = np.convolve(fine, kernel, mode="same")
    return np.interp(np.log(targets), log_grid, smooth)


def choose_components(
    log_flux: np.ndarray,
    n_a: int,
    *,
    target_a: float = TARGET_A,
    target_b: float = TARGET_B,
    cap: int = MAX_COMPONENTS,
) -> tuple[int, np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    """The smallest K meeting both targets (capped), and the PCA pieces.

    *log_flux* is ``(n_spectra, n_pixels)`` of ``log10`` normalised flux, the
    first *n_a* pixels segment A. Returns ``(K, mean, eigenspectra (K,
    n_pixels), weights (n_spectra, K), errors)`` where *errors* holds the worst
    per-spectrum RMS fractional flux error on each segment at the chosen K.
    """
    mean = log_flux.mean(axis=0)
    centred = log_flux - mean
    _, _, vt = np.linalg.svd(centred, full_matrices=False)
    chosen = min(cap, vt.shape[0])
    errors: dict[str, float] = {}
    for k in range(1, min(cap, vt.shape[0]) + 1):
        basis = vt[:k]
        rebuilt = mean + (centred @ basis.T) @ basis
        frac = 10.0 ** (rebuilt - log_flux) - 1.0
        rms_a = float(np.sqrt(np.mean(frac[:, :n_a] ** 2, axis=1)).max())
        rms_b = float(np.sqrt(np.mean(frac[:, n_a:] ** 2, axis=1)).max())
        errors = {"rms_a": rms_a, "rms_b": rms_b}
        if rms_a <= target_a and rms_b <= target_b:
            chosen = k
            break
    basis = vt[:chosen]
    return chosen, mean, basis, centred @ basis.T, errors


def _kernel(x1: np.ndarray, x2: np.ndarray, lengths: np.ndarray, amplitude: float) -> np.ndarray:
    d = (x1[:, None, :] - x2[None, :, :]) / lengths
    return amplitude**2 * np.exp(-0.5 * np.sum(d**2, axis=-1))


def gp_mean(x: np.ndarray, nodes: np.ndarray, alpha: np.ndarray, hyper: np.ndarray) -> np.ndarray:
    """Predictive means ``k(x, X) alpha`` for every component.

    *x* ``(n, 3)`` and *nodes* ``(m, 3)`` standardised; *alpha* ``(m, K)``;
    *hyper* ``(K, 5)`` -- three length scales, amplitude, jitter. Returns
    ``(n, K)``. The emulator modules implement exactly this, per backend.
    """
    out = np.empty((x.shape[0], alpha.shape[1]))
    for k in range(alpha.shape[1]):
        out[:, k] = _kernel(x, nodes, hyper[k, :3], hyper[k, 3]) @ alpha[:, k]
    return out


def fit_gp(nodes: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Maximise the marginal likelihood of one ARD SE GP; return ``(hyper, alpha)``.

    *hyper* is ``[l1, l2, l3, amplitude, jitter]``; the jitter is a variance
    relative to ``amplitude**2``, floored at :data:`JITTER_FLOOR` and capped at
    :data:`JITTER_CAP`. The fit is on
    log hyperparameters with L-BFGS-B.
    """
    scale = float(np.std(y)) or 1.0
    m = y.shape[0]

    def unpack(theta: np.ndarray) -> tuple[np.ndarray, float, float]:
        return np.exp(theta[:3]), float(np.exp(theta[3])), float(np.exp(theta[4]))

    def nll(theta: np.ndarray) -> float:
        lengths, amp, jitter = unpack(theta)
        cov = _kernel(nodes, nodes, lengths, amp) + amp**2 * jitter * np.eye(m)
        try:
            chol = np.linalg.cholesky(cov)
        except np.linalg.LinAlgError:
            return 1e25
        solved = np.linalg.solve(chol, y)
        return float(0.5 * solved @ solved + np.log(np.diag(chol)).sum())

    start = np.array([0.0, 0.0, 0.0, np.log(scale), np.log(1e-6)])
    bounds = [(np.log(LENGTH_BOUNDS[0]), np.log(LENGTH_BOUNDS[1]))] * 3 + [
        (np.log(scale) - 5.0, np.log(scale) + 5.0),
        (np.log(JITTER_FLOOR), np.log(JITTER_CAP)),
    ]
    best = scipy.optimize.minimize(nll, start, method="L-BFGS-B", bounds=bounds)
    lengths, amp, jitter = unpack(best.x)
    cov = _kernel(nodes, nodes, lengths, amp) + amp**2 * jitter * np.eye(m)
    alpha = np.linalg.solve(cov, y)
    return np.array([*lengths, amp, jitter]), alpha


def fit_all(nodes: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """:func:`fit_gp` per column of *weights*; ``(hyper (K, 5), alpha (m, K))``."""
    hypers, alphas = [], []
    for k in range(weights.shape[1]):
        hyper, alpha = fit_gp(nodes, weights[:, k])
        hypers.append(hyper)
        alphas.append(alpha)
    return np.array(hypers), np.array(alphas).T


def reduce_file(wave_aa: np.ndarray, flux_per_aa: np.ndarray, segment_b: np.ndarray) -> np.ndarray:
    """One spectrum's ``log10`` unit-bolometric shape on segments A and B (per angstrom)."""
    fbol = np.trapezoid(flux_per_aa, wave_aa)
    shape = flux_per_aa / fbol
    wave_um = wave_aa / 1e4
    seg_a = bin_mean(wave_um, shape, bin_edges(SEGMENT_A))
    seg_b = lsf_sample(wave_um, shape, segment_b)
    return np.log10(np.concatenate([seg_a, seg_b]))


def _download_and_reduce(segment_b: np.ndarray) -> tuple[np.ndarray, list[str], str]:
    from astropy.io import fits
    from astropy.utils.data import download_file

    with fits.open(download_file(WAVE_URL, cache=True)) as hdul:
        wave_aa = np.asarray(hdul[0].data, dtype=float)
    rows, names, version = [], [], ""
    for teff, logg, feh in grid_nodes():
        url = spectrum_url(teff, logg, feh)
        with fits.open(download_file(url, cache=True)) as hdul:
            header = hdul[0].header
            unit = str(header.get("BUNIT", "")).replace(" ", "")
            if unit != "erg/s/cm^2/cm":
                raise RuntimeError(f"{url}: unexpected BUNIT {unit!r}")
            version = str(header.get("PHXVER", version))
            flux = np.asarray(hdul[0].data, dtype=float) / 1e8  # per cm -> per angstrom
        rows.append(reduce_file(wave_aa, flux, segment_b))
        names.append(url.rsplit("/", 1)[1])
        print(f"  reduced {names[-1]}", flush=True)
    return np.array(rows), names, version


def _standardise(nodes: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    centre = nodes.mean(axis=0)
    spread = nodes.std(axis=0)
    return (nodes - centre) / spread, centre, spread


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)

    segment_b = rvs_grid()
    assert segment_b.size == 387, segment_b.size
    n_a = SEGMENT_A.size
    print("downloading (astropy cache) and reducing 156 spectra ...", flush=True)
    log_flux, names, version = _download_and_reduce(segment_b)
    nodes = grid_nodes()

    k, mean, basis, weights, pca_errors = choose_components(log_flux, n_a)
    print(
        f"PCA: K = {k}, worst RMS fractional error A {pca_errors['rms_a']:.4g}, "
        f"B {pca_errors['rms_b']:.4g}",
        flush=True,
    )

    x, centre, spread = _standardise(nodes)
    hyper, alpha = fit_all(x, weights)
    print("GP hyperparameters (l_teff, l_logg, l_feh, amplitude, jitter):")
    for row in hyper:
        print("  " + "  ".join(f"{v:.4g}" for v in row))

    # (ii) the emulator at every node, against its PCA reconstruction and the truth.
    rebuilt = mean + weights @ basis
    emulated = mean + gp_mean(x, x, alpha, hyper) @ basis
    to_pca = 10.0 ** (emulated - rebuilt) - 1.0
    to_truth = 10.0 ** (emulated - log_flux) - 1.0
    node_error = {
        "max_vs_pca": float(np.abs(to_pca).max()),
        "max_vs_truth": float(np.abs(to_truth).max()),
        "rms_vs_truth_a": float(np.sqrt(np.mean(to_truth[:, :n_a] ** 2))),
        "rms_vs_truth_b": float(np.sqrt(np.mean(to_truth[:, n_a:] ** 2))),
    }
    print(f"(ii) at the nodes: {node_error}", flush=True)

    # (i) leave-one-out on every seventh node, the basis kept, the GPs refitted.
    frac_a, frac_b = [], []
    for held in range(0, nodes.shape[0], 7):
        keep = np.arange(nodes.shape[0]) != held
        h_hyper, h_alpha = fit_all(x[keep], weights[keep])
        predicted = mean + gp_mean(x[held : held + 1], x[keep], h_alpha, h_hyper)[0] @ basis
        frac = 10.0 ** (predicted - log_flux[held]) - 1.0
        frac_a.append(frac[:n_a])
        frac_b.append(frac[n_a:])
    loo = {
        "held_out": len(frac_a),
        "max_a": float(np.abs(frac_a).max()),
        "rms_a": float(np.sqrt(np.mean(np.square(frac_a)))),
        "max_b": float(np.abs(frac_b).max()),
        "rms_b": float(np.sqrt(np.mean(np.square(frac_b)))),
    }
    print(
        f"(i) leave-one-out on {loo['held_out']} nodes: A max {loo['max_a']:.3g} "
        f"rms {loo['rms_a']:.3g}; B max {loo['max_b']:.3g} rms {loo['rms_b']:.3g}",
        flush=True,
    )

    stored = np.linspace(0, nodes.shape[0] - 1, STORED_NODES).astype(int)
    provenance = {
        "source": "PHOENIX-ACES-AGSS-COND-2011 HiRes (Husser et al. 2013), Göttingen",
        "phoenix_version": version,
        "files": names,
        "n_files": len(names) + 1,
        "wavelength_file": WAVE_URL.rsplit("/", 1)[1],
        "date": datetime.date.today().isoformat(),
        "components": k,
        "pca_reconstruction": pca_errors,
        "node_reproduction": node_error,
        "leave_one_out": loo,
        "vacuum_wavelengths": True,
        "script": "python -m examples.phoenix_star.train_emulator",
    }
    np.savez_compressed(
        args.output,
        segment_a=SEGMENT_A,
        segment_b=segment_b,
        mean=mean,
        basis=basis,
        nodes=x,
        input_centre=centre,
        input_spread=spread,
        alpha=alpha,
        hyper=hyper,
        grid_parameters=nodes,
        stored_index=stored,
        stored_spectra=log_flux[stored],
        provenance=np.array(json.dumps(provenance)),
    )
    print(f"wrote {args.output} ({args.output.stat().st_size} bytes)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
