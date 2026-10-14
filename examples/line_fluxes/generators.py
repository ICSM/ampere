"""The synthetic, JWST-like line list and the line-flux table drawn from it (W7.9).

**The line list.** Ten mid-infrared nebular transitions inside MIRI's
5 to 28 µm range — the lines a planetary-nebula or H II-region spectrum from
MIRI's medium-resolution spectrometer is reduced to (issue #70's use case).
The rest wavelengths are vacuum values rounded to 0.1 nm, from the
fine-structure line compilation of Feuchtgruber et al. (1997, ApJ 487, 962)
and the NIST Atomic Spectra Database (Kramida et al.), and for H I 7-6
(Humphreys alpha) from the Rydberg formula.

**The ladder is a toy.** The upper-level energies are ``E_k = 30000 K µm /
λ_k`` — higher for the shorter-wavelength lines, so the excitation structure
lies *along* rest wavelength — and the strengths are chosen so that the
lines are of comparable brightness at the truth, as the strong nebular lines
in a MIRI spectrum are. Neither is atomic data; the inference pattern is the
point (see :mod:`.model`).

**The truth has two components**: a cool one feeding every line and a hot
one feeding mainly the high-excitation (short-wavelength) lines. The study's
misspecified arms fit one component, so they cannot reproduce the ladder's
curvature, and the residual is correlated across neighbouring transitions in
rest wavelength — what the flexible likelihood's kernel over rest wavelength
is there to absorb.

**Synthetic data, negotiated then compiled then evaluated**, as
:mod:`examples.astrometry.generators` does it: the truth model adopts the
line-integration step's grid before it is evaluated, exactly as
:class:`~ampere.core.FittingProblem` would.

No real line list is used: the JWST file vendored with W7.3 is a MIRI LRS
spectrum of a flux standard, a star with no emission lines, so there is no
table to read from it (the page's "What the template did not say" section).
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import Instrument, LineFluxes, negotiate

__all__ = [
    "ENERGY",
    "LABELS",
    "REST",
    "SEED",
    "STRENGTH",
    "TRUTH",
    "UNCERTAINTY",
    "WIDTH",
    "WINDOW",
    "line_table",
    "synthetic_table",
]

#: Reproducible everywhere this example is run.
SEED = 20261014

#: The transitions, in the order a table might list them (by ion, not by
#: wavelength: the kind tolerates any order).
LABELS: tuple[str, ...] = (
    "[Ar II] 6.99",
    "[Ar III] 8.99",
    "[S IV] 10.51",
    "H I 7-6",
    "[Ne II] 12.81",
    "[Ne V] 14.32",
    "[Ne III] 15.56",
    "[S III] 18.71",
    "[Ne V] 24.32",
    "[O IV] 25.89",
)
#: Their vacuum rest wavelengths, micron (see the module docstring's sources).
REST = np.array(
    [6.9853, 8.9914, 10.5105, 12.3719, 12.8136, 14.3217, 15.5551, 18.7130, 24.3175, 25.8903]
)
#: The toy upper-level energies, K: ``30000 K µm / λ``.
ENERGY = 30000.0 / REST
#: The lines' Gaussian velocity dispersion, km/s.
WIDTH = 50.0
#: The line-integration window's half-width, km/s — six dispersions.
WINDOW = 300.0
#: Fractional uncertainty on each tabulated flux.
UNCERTAINTY = 0.03

#: The two-component truth: a cool component feeding every line and a hot one
#: feeding mainly the short-wavelength, high-excitation lines.
TRUTH: dict[str, float] = {
    "temperature": 600.0,
    "column": 3000.0,
    "hot_temperature": 6000.0,
    "hot_column": 50.0,
}

#: The toy strengths: chosen so the cool component alone gives lines of order
#: unity (in 1e-16 W m⁻²) at the truth's temperature, scattered by a factor of
#: about two, as a real nebular line list's strong lines are.
STRENGTH = (
    np.exp(ENERGY / TRUTH["temperature"])
    * np.array([0.6, 1.4, 0.9, 0.5, 1.8, 0.7, 1.5, 1.0, 0.8, 1.2])
    / np.exp(ENERGY / TRUTH["temperature"]).sum()
)


def line_table(fluxes: Any, sigma: Any, *, order: Any = None) -> LineFluxes:
    """A :class:`~ampere.core.LineFluxes` table over :data:`LABELS`, in W m⁻².

    *order* permutes the rows, as a published table's own ordering would; the
    composition aligns by transition label, so any order fits identically.
    This is also the shape of building one by hand from a real table — a
    label, a rest wavelength, a flux and its uncertainty per row.
    """
    index = np.arange(len(LABELS)) if order is None else np.asarray(order)
    return LineFluxes(
        [LABELS[i] for i in index],
        REST[index] * u.um,
        np.asarray(fluxes, dtype=float)[index] * u.W / u.m**2,
        uncertainty=np.asarray(sigma, dtype=float)[index] * u.W / u.m**2,
    )


def synthetic_table(truth_model: Any, instrument: Instrument, *, seed: int = SEED) -> LineFluxes:
    """The truth's line fluxes through *instrument*, with Gaussian noise of :data:`UNCERTAINTY`.

    Negotiate, compile, evaluate — then draw. The uncertainty is the stated
    fraction of the *true* flux (an observer's error bar does not know the
    draw), and the table is returned in reverse label order to exercise the
    alignment.
    """
    compiled = truth_model.compile_for(negotiate([instrument]))
    clean = np.asarray(instrument(compiled.evaluate()).values, dtype=float)
    sigma = UNCERTAINTY * clean
    rng = np.random.default_rng(seed)
    noisy = clean + sigma * rng.standard_normal(clean.size)
    return line_table(noisy, sigma, order=np.arange(len(LABELS))[::-1])
