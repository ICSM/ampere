"""The truth, and the synthetic image it produces.

Everything about the *data* lives here; :mod:`.model` supplies the physics the
fit is written in, and :mod:`.study` composes the two. The split is
:mod:`examples.m2_misspecification`'s and :mod:`examples.interferometry`'s,
kept so that a reader who has followed either recognises the shape.

The truth, and the thing the model gets wrong
---------------------------------------------
A compact Gaussian source on a smooth, much broader background — an extended
halo, scattered light, an unsubtracted sky, whichever story you prefer; what
matters is that it is **smooth on the scale of the image** and that the
misspecified model omits it entirely. That is the same misspecification shape
M2 uses on an SED and :mod:`examples.interferometry` uses on a visibility set:
a real component, left out of the model, whose residual is correlated rather
than white.

It is chosen to be smooth *in the observed coordinates* deliberately. A
flexible likelihood absorbs an unmodelled component when the kernel can
represent its residual, and a kernel sees a container's axes and nothing else
(``docs/source/interferometry.rst`` §1's chromatic lesson, in two spatial
dimensions rather than one spectral one). An ``Image``'s axes are ``x`` and
``y``, so the component this study omits must be smooth in ``x`` and ``y`` —
which a broad Gaussian is, and which, say, a periodic detector fringe would
not be.

The instrument is one step: a Gaussian PSF of
:data:`PSF_FWHM`, wider than the pixel scale at every size this study runs at,
so that the convolution is doing real work rather than being a near-identity.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np

from ampere.core import Image

__all__ = [
    "BACKGROUND",
    "FIELD_OF_VIEW",
    "PSF_FWHM",
    "SEED",
    "SIGMA_FRACTION",
    "SOURCE",
    "TRUTH",
    "camera",
    "field",
    "observed_template",
    "pixel_scale",
    "psf_step",
    "seed_grid",
    "synthetic",
    "truth_model",
]

#: The observed field, mas, on both axes. Fixed across every image size this
#: study runs at, so that changing ``pixels`` changes the *sampling* and
#: nothing else — which is what makes the three-N benchmark a statement about
#: N rather than about the field.
FIELD_OF_VIEW = 24.0

#: The compact source: the thing every arm fits.
SOURCE: dict[str, float] = {"flux": 2.4, "fwhm": 3.2}

#: The smooth background: present in the data, present in the "correct" arm at
#: exactly this value, and absent from the other two.
BACKGROUND: dict[str, float] = {"flux": 1.1, "fwhm": 17.0}

#: The instrument's point-spread function, mas.
PSF_FWHM = 2.5

#: Per-pixel measurement noise, as a fraction of the brightest observed pixel.
#: Uncorrelated and known, which is the M2 setting: the *only* thing wrong in
#: the misspecified arms is the missing background.
SIGMA_FRACTION = 0.03

#: One seed for the whole study, as every other example here uses.
SEED = 20260916

#: The free parameters and their true values, under the names a posterior
#: records them by.
TRUTH: dict[str, float] = {
    "model.flux": SOURCE["flux"],
    "model.fwhm": SOURCE["fwhm"],
}

#: The surface-brightness unit W4.1's image models emit in, and so the unit an
#: image of one is measured in. ``Likelihood.check_alignment`` compares units,
#: so this is load-bearing rather than cosmetic.
BRIGHTNESS = u.Jy / u.sr


def field(pixels: int) -> np.ndarray:
    """*pixels* centres evenly spaced across :data:`FIELD_OF_VIEW`, mas."""
    half = 0.5 * FIELD_OF_VIEW
    return np.linspace(-half, half, int(pixels))


def pixel_scale(pixels: int) -> float:
    """The pixel scale of :func:`field`, mas."""
    grid = field(pixels)
    return float(grid[1] - grid[0])


def seed_grid(pixels: int = 4) -> np.ndarray:
    """A deliberately useless placeholder grid for a model's constructor.

    Negotiation replaces it on the first ``compile_for`` — the PSF step
    publishes the grid it needs, and every model here adopts it — so a model's
    own grid is a constructor formality. Four coordinates rather than the real
    ones, so that a fit which somehow *failed* to negotiate produces an
    obviously wrong four-pixel image rather than a plausible one.
    (:mod:`examples.interferometry.generators` makes the same argument for the
    same reason.)
    """
    return np.linspace(-0.5 * FIELD_OF_VIEW, 0.5 * FIELD_OF_VIEW, int(pixels))


def observed_template(pixels: int, values: np.ndarray | None = None) -> Image:
    """An ``Image`` on :func:`field`, with the study's uniform uncertainty.

    With *values* ``None`` this is the empty container the PSF step takes its
    target grid from (``from_observed`` — the supported route); with values it
    is the observation itself.
    """
    grid = field(pixels)
    filled = np.zeros((pixels, pixels)) if values is None else np.asarray(values, dtype=float)
    scale = float(np.max(np.abs(filled))) or 1.0
    return Image(
        grid * u.mas,
        grid * u.mas,
        filled * BRIGHTNESS,
        uncertainty=np.full((pixels, pixels), SIGMA_FRACTION * scale) * BRIGHTNESS,
    )


def psf_step(module: Any, observed: Image, *, label: str = "psf") -> Any:
    """This backend's :class:`PSFConvolution`, built from the observed image."""
    return module.PSFConvolution.from_observed(observed, fwhm=PSF_FWHM, label=label)


def camera(module: Any, observed: Image, *, channel: str = "sky") -> Any:
    """The one-step imaging instrument: a PSF convolution and nothing else."""
    from ampere.core import Instrument

    return Instrument([psf_step(module, observed)], channel=channel, label="camera")


def truth_model(itf: Any, pixels: int) -> Any:
    """The generating sky: the source at :data:`SOURCE`, plus the background."""
    from .model import SourceWithBackground

    grid = seed_grid()
    return SourceWithBackground(
        itf,
        grid,
        grid,
        flux=SOURCE["flux"],
        fwhm=SOURCE["fwhm"],
        background_flux=BACKGROUND["flux"],
        background_fwhm=BACKGROUND["fwhm"],
    )


def synthetic(module: Any, itf: Any, pixels: int, *, seed: int = SEED) -> tuple[Any, Image]:
    """One noisy image of the truth, and the compiled model that produced it.

    The noise is drawn here rather than through ``FittingProblem.simulate``
    because the observation has to exist *before* the problem does: the PSF
    step takes its target grid from the observed container, which is
    ``transformations.md`` §10's coordinates-from-the-container rule and the
    reason every modality in this repository builds its data first.
    """
    template = observed_template(pixels)
    instrument = camera(module, template)
    from ampere.core import negotiate

    compiled = truth_model(itf, pixels).compile_for(negotiate([instrument]))
    noiseless = instrument(compiled.evaluate())
    values = np.asarray(noiseless.values, dtype=float)
    sigma = SIGMA_FRACTION * float(np.max(np.abs(values)))
    rng = np.random.default_rng(seed)
    observed = Image(
        template.x.values * u.mas,
        template.y.values * u.mas,
        (values + rng.normal(0.0, sigma, values.shape)) * BRIGHTNESS,
        uncertainty=np.full(values.shape, sigma) * BRIGHTNESS,
    )
    return compiled, observed
