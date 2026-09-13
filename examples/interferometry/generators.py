"""Synthetic data for the interferometry study: a binary plus a fainter disc.

The truth (:data:`BINARY`, :data:`DISC_FLUX`, :data:`DISC_FWHM`) is a resolved
binary — the same one every interferometry suite in this repository fits —
with a much fainter, much more extended circular Gaussian added underneath
it, standing in for a circumstellar disc or an extended halo a naive model
omits. Three of the study's four arms (``examples.interferometry.study``'s
``"correct"``, ``"incomplete"`` and ``"flexible"``) share this one truth and
one noise realisation, observed as both visibilities and closure phases,
because the M2 question is about what happens when *one* dataset is fitted
with the disc present in the model, absent from it, and absent from it but
covered by a Gaussian process.

Geometry, not re-derived
------------------------
The array (four stations), the hour-angle coverage, the field of view and the
Nyquist oversampling are :mod:`tests.backends.interferometry_fixtures`'s —
imported rather than retyped, exactly as ``tests/inference/test_interferometry.py``
and ``tests/backends/test_native_interferometry.py`` do, by putting
``tests/backends`` on ``sys.path`` (that module's own docstring gives the
reason: ``tests/`` is not a package, so pytest puts a test file's own
directory on ``sys.path`` and not its siblings'). A second, independent
transcription of four station positions and four hour angles would be a
second place for the sign of a closure phase to go quietly wrong.

The chromatic arm's coverage — several wavelengths of the *same* baselines —
is this module's own (:func:`dispersed_visibility_coverage`): the fixtures
module is monochromatic by design (the reference ``Image`` channel is
achromatic, W4.1's carried limitation), and dispersion is exactly what W4.4's
chromatic scenario adds.

Why the disc is a ``GaussianSource``, not a ``UniformDisc``
-------------------------------------------------------------
``phase4_placement_memo.md`` and W4.3's status row both note that
``UniformDisc`` is ``DIFFERENTIABLE = False`` on the modern backends (its
visibility is a Bessel function ampere does not carry a native derivative
of), while ``GaussianSource`` is differentiable everywhere. Using the Gaussian
for the disc means the *same* truth, and the *same* "correct" fitted model,
can run under emcee, NUTS and NPE without a second recipe for the modern
backends — "one truth for all three engines is simpler" (the item's design
guidance, verbatim).
"""

from __future__ import annotations

import pathlib
import sys
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import ClosurePhases, VisibilitySet

# ``tests/backends/interferometry_fixtures.py``'s own docstring gives the
# reason for the explicit sys.path insert: an editable install's path hook
# reaching this far would break the first time this ran against a wheel.
_FIXTURES_DIR = pathlib.Path(__file__).resolve().parents[2] / "tests" / "backends"
if str(_FIXTURES_DIR) not in sys.path:
    sys.path.insert(0, str(_FIXTURES_DIR))

import interferometry_fixtures as fixtures  # noqa: E402

__all__ = [
    "BINARY",
    "COMPONENT_FWHM",
    "DISC_FLUX",
    "DISC_FWHM",
    "FIELD_OF_VIEW",
    "HOUR_ANGLES",
    "LINE_CENTRE",
    "LINE_WIDTH",
    "N_CHANNELS",
    "OVERSAMPLING",
    "PATCH_FLUX",
    "PATCH_OFFSET",
    "SEED",
    "SIGMA_CLOSURE",
    "SIGMA_VISIBILITY",
    "STATIONS",
    "TRUTH",
    "WAVELENGTH",
    "chromatic_synthetic",
    "dispersed_visibility_coverage",
    "synthetic",
]

#: Re-exported for readers of this module and for :mod:`.study`, so neither
#: has to know the constants live one module further down. Deliberately not
#: re-declared: two numbers for one array is exactly the risk the module
#: docstring warns about.
STATIONS = fixtures.STATIONS
HOUR_ANGLES = fixtures.HOUR_ANGLES
FIELD_OF_VIEW = fixtures.FIELD_OF_VIEW
OVERSAMPLING = fixtures.OVERSAMPLING
WAVELENGTH = fixtures.WAVELENGTH
BINARY = fixtures.BINARY
COMPONENT_FWHM = fixtures.COMPONENT_FWHM
SIGMA_VISIBILITY = fixtures.SIGMA_VISIBILITY
SIGMA_CLOSURE = fixtures.SIGMA_CLOSURE
SEED = fixtures.SEED

#: The disc's own truth: a **buffer** on every fitted model in this study,
#: never a free parameter — arm (a)'s "correct" model is correct because it
#: is told the disc's properties exactly, not because it fits them. Total
#: flux a factor of eleven below the binary's, so it is genuinely the fainter
#: component. Its FWHM is sized against the array's own resolution
#: (``1/(2 u_max)``, about 2.5 mas for these baselines at
#: :data:`WAVELENGTH`): three times the binary components' own FWHM, so it is
#: **partially** resolved — a sizeable fraction of its flux survives on the
#: shorter baselines and falls off smoothly on the longer ones, rather than
#: either sitting on top of the binary's own signal (too compact to matter)
#: or resolving out to nothing everywhere the array observes (too extended
#: to matter, which an earlier, ten-mas-wider choice for this constant did).
#: That (u, v)-dependent, smoothly-varying excess is exactly what an
#: independent-noise fit of the binary alone cannot absorb without biasing
#: the binary's own parameters (arm "incomplete"), and what a Gaussian
#: process over (u, v) plausibly can (arm "flexible").
DISC_FLUX = 0.18
DISC_FWHM = 6.0

#: The four physical parameters :mod:`.study` reports on, matching
#: :mod:`interferometry_fixtures`'s own naming.
TRUTH: dict[str, float] = {
    "model.separation": BINARY["separation"],
    "model.flux_ratio": BINARY["flux_ratio"],
}

# ---------------------------------------------------------------------------
# Arms (a)-(c): one truth, one noise draw, with and without the disc.
# ---------------------------------------------------------------------------


def seed_grid(pixels: int = 4) -> np.ndarray:
    """A small, strictly increasing placeholder grid, mas.

    Every image model needs *some* grid at construction, but with
    ``adopt_grid=True`` (the default) it is discarded the moment
    :class:`~ampere.core.dataset.FittingProblem` negotiates and calls
    ``compile_for`` — :mod:`interferometry_fixtures.truth_model` does the same
    four-pixel placeholder for the same reason. Kept to four points here too,
    so a reader comparing the two modules sees one convention rather than two.
    """
    half = 0.5 * FIELD_OF_VIEW
    return np.linspace(-half, half, int(pixels))


def truth_model(itf: Any) -> Any:
    """The binary-plus-disc truth, on *itf* (:mod:`.model` builds the composite)."""
    from . import model as _model

    grid = seed_grid()
    return _model.binary_with_disc(
        itf,
        grid,
        grid,
        disc_flux=DISC_FLUX,
        disc_fwhm=DISC_FWHM,
        component_fwhm=COMPONENT_FWHM,
        **BINARY,
    )


def synthetic(itf: Any, *, seed: int = SEED) -> tuple[Any, VisibilitySet, ClosurePhases]:
    """Negotiate, compile, evaluate and observe the binary-plus-disc truth.

    The reference recipe :mod:`interferometry_fixtures`'s own ``synthetic``
    follows, with the truth swapped for :func:`truth_model`'s composite —
    everything about the noise (circular complex on the visibilities, von
    Mises on the phases) and the coverage is unchanged, so the two studies'
    data are comparable sample for sample. Always built on the **reference**
    backend: the observed containers this returns are plain arrays wrapped in
    ``ampere.core`` containers, and a backend does not enter until
    :mod:`.study` binds them into a problem, so there is exactly one
    synthetic dataset, fitted by every engine, rather than one per backend.
    """
    vis_instrument = fixtures.chain(itf, fixtures.visibilities(), "vis")
    t3_instrument = fixtures.chain(itf, fixtures.closure_phases(), "t3")
    from ampere.core import negotiate

    compiled = truth_model(itf).compile_for(negotiate([vis_instrument, t3_instrument]))
    result = compiled.evaluate()
    clean_visibility = np.asarray(vis_instrument(result).values)
    clean_phase = np.asarray(t3_instrument(result).values)
    rng = np.random.default_rng(seed)
    noisy_visibility = clean_visibility + SIGMA_VISIBILITY * (
        rng.standard_normal(clean_visibility.size) + 1j * rng.standard_normal(clean_visibility.size)
    )
    noisy_phase = np.angle(np.exp(1j * rng.vonmises(clean_phase, 1.0 / SIGMA_CLOSURE**2)))
    return compiled, fixtures.visibilities(noisy_visibility), fixtures.closure_phases(noisy_phase)


# ---------------------------------------------------------------------------
# Arm (d): the chromatic case (memo §3.6). Visibilities only — a compact
# patch does not need a closure phase to make its case, and the fitted model
# stays the achromatic binary throughout.
# ---------------------------------------------------------------------------

#: Number of dispersed channels. Small: the point is several *distinct*
#: wavelengths spanning the line, not a densely sampled spectrum.
N_CHANNELS = 7

#: The band the patch emits in, micron: centred on the array's own
#: monochromatic wavelength so the (u, v) footprint stays close to the rest
#: of the study's, one tenth of it wide across the sampled range.
LINE_CENTRE = WAVELENGTH
LINE_WIDTH = 0.01

#: The patch's peak flux (at line centre) and its offset from the phase
#: centre, mas. Off-centre and unresolved (no width of its own — treated as a
#: point, since only its spectral profile is under test here): its Fourier
#: transform is a pure phase times a constant modulus, which is what makes
#: "evaluate the binary and add the patch's analytic visibility" (the design
#: guidance, verbatim) exact rather than approximate.
PATCH_FLUX = 0.12
PATCH_OFFSET = (6.0, -4.0)


def dispersed_visibility_coverage() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(u, v, lambda)``: the array's baselines, at every hour angle, at every channel.

    The same physical baselines as :func:`interferometry_fixtures.visibility_coverage`
    at :data:`N_CHANNELS` wavelengths across the line — :math:`(u, v) = B/\\lambda`
    per channel, so the *same* station pair traces a different point in the
    (u, v) plane at each wavelength, which is what makes wavelength a genuine
    third coordinate rather than one absorbed into (u, v) (memo §3.6).
    """
    wavelengths = np.linspace(LINE_CENTRE - 3.0 * LINE_WIDTH, LINE_CENTRE + 3.0 * LINE_WIDTH, N_CHANNELS)
    metres = wavelengths[:, None, None] * 1e-6
    pairs = np.array(
        [
            (STATIONS[j] - STATIONS[i])
            for i in range(len(STATIONS))
            for j in range(i + 1, len(STATIONS))
        ]
    )
    # (channel, baseline, 2), then rotated by every hour angle and flattened.
    baselines = pairs[None, :, :] / metres  # (channel, baseline, 2)
    cos_h = np.cos(HOUR_ANGLES)
    sin_h = np.sin(HOUR_ANGLES)
    u_pts = (
        baselines[:, :, 0:1] * cos_h[None, None, :] - baselines[:, :, 1:2] * sin_h[None, None, :]
    )
    v_pts = (
        baselines[:, :, 0:1] * sin_h[None, None, :] + baselines[:, :, 1:2] * cos_h[None, None, :]
    )
    waves = np.broadcast_to(wavelengths[:, None, None], u_pts.shape)
    return u_pts.reshape(-1), v_pts.reshape(-1), waves.reshape(-1)


def band_profile(wavelength: np.ndarray) -> np.ndarray:
    """``S(lambda)``: a Gaussian line, peak :data:`PATCH_FLUX`, sharp in wavelength."""
    return PATCH_FLUX * np.exp(-0.5 * ((wavelength - LINE_CENTRE) / LINE_WIDTH) ** 2)


def patch_visibility(u_pts: np.ndarray, v_pts: np.ndarray, wavelength: np.ndarray) -> np.ndarray:
    """The compact patch's analytic visibility, ``S(lambda) * exp(-2 pi i (u dx + v dy))``.

    A point source at :data:`PATCH_OFFSET`, so its Fourier transform is exact
    on any baseline — no image grid, no Nyquist requirement, which is what
    lets this term be added directly to the achromatic binary's direct-DFT
    visibility (design guidance: "generate the chromatic data by evaluating
    the achromatic binary and adding the patch's analytic visibility").
    """
    from ampere.backends.reference.interferometry import MAS_PER_RAD

    dx, dy = PATCH_OFFSET
    phase = -2j * np.pi * (u_pts * dx + v_pts * dy) / MAS_PER_RAD
    return band_profile(wavelength).astype(np.complex128) * np.exp(phase)


def chromatic_synthetic(itf: Any, *, seed: int = SEED) -> tuple[Any, VisibilitySet]:
    """The achromatic binary plus the chromatic patch, dispersed, noisy.

    The *fitted* model in every one of arm (d)'s three sub-fits is the plain
    achromatic :class:`~ampere.backends.reference.interferometry.Binary` —
    the patch is entirely omitted, exactly as the disc is in arms
    ``"incomplete"`` and ``"flexible"`` — so this truth is generated once,
    on the reference backend, from :data:`BINARY` (no disc component here:
    the chromatic case is its own, simpler misspecification) plus the
    patch's analytic contribution.
    """
    u_pts, v_pts, waves = dispersed_visibility_coverage()
    vis_instrument = fixtures.chain(
        itf,
        VisibilitySet(
            u_pts,
            v_pts,
            waves * u.micron,
            np.zeros(u_pts.size, dtype=np.complex128) * u.Jy,
            uncertainty=np.full(u_pts.size, SIGMA_VISIBILITY) * u.Jy,
        ),
        "vis",
    )
    from ampere.core import negotiate

    grid = seed_grid()
    truth = itf.Binary(grid, grid, channels="sky", component_fwhm=COMPONENT_FWHM, **BINARY)
    compiled = truth.compile_for(negotiate([vis_instrument]))
    result = compiled.evaluate()
    clean = np.asarray(vis_instrument(result).values) + patch_visibility(u_pts, v_pts, waves)
    rng = np.random.default_rng(seed)
    noisy = clean + SIGMA_VISIBILITY * (
        rng.standard_normal(clean.size) + 1j * rng.standard_normal(clean.size)
    )
    observed = VisibilitySet(
        u_pts,
        v_pts,
        waves * u.micron,
        noisy * u.Jy,
        uncertainty=np.full(u_pts.size, SIGMA_VISIBILITY) * u.Jy,
    )
    return compiled, observed
