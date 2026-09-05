"""The native spectral models: blackbody, modified blackbody, power law.

``DEVELOPMENT_PLAN.md`` §5's Phase 2 list, on the reference path. All three are
written straight out in numpy against ``astropy.constants``, because
correctness is the reference backend's only goal (``architecture.md`` §2): no
gradients, no batching, no GPU, and no cleverness that would need its own test.

Every model here:

* emits **Jy** on a spectral axis in **micron**, the units the v1 slice works
  in, converting once at construction rather than per evaluation
  (``DEVELOPMENT_PLAN.md`` §7's units trap);
* takes its grid at construction and registers it as a **buffer**, since one
  would never put a prior on it (``parameters.md`` §10);
* accepts each physical quantity as a prior (fitted), a number (held fixed) or
  a ready-made :class:`~ampere.core.Parameter`, so the same class serves a fit
  and a fixed-parameter simulation;
* honours :meth:`~ampere.core.Model.compile_for` by building one container per
  requested channel and refilling it with ``with_values``, so successive
  results share their axes by identity (``transformations.md`` §14).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.constants as const
import astropy.units as u
import numpy as np

from ampere.core import (
    DTYPE,
    ChannelRequirements,
    Model,
    ModelResult,
    Spectrum,
)

from ._declare import as_parameter

__all__ = [
    "COORDINATE_UNIT",
    "FLUX_UNIT",
    "BlackBody",
    "ModifiedBlackBody",
    "PowerLaw",
    "planck_jy",
]

#: The spectral coordinate unit every model here works in.
COORDINATE_UNIT = u.micron
#: The flux unit every model here emits.
FLUX_UNIT = u.Jy

# Converted once, at import, into the bare-float constants the hot loop uses
# (DEVELOPMENT_PLAN.md §7: units are converted at composition time, never per
# evaluation). 2h/c^2 carries J s^3 m^-2, so multiplying it by nu^3 in Hz gives
# J s^-1 m^-2 Hz^-1 sr^-1 -- the SI spectral radiance B_nu.
_TWO_H_OVER_C2 = float((2.0 * const.h / const.c**2).to_value(u.J * u.s**3 / u.m**2))
_H_OVER_K = float((const.h / const.k_B).to_value(u.s * u.K))
_C_MICRON_HZ = float(const.c.to_value(u.micron * u.Hz))
#: 1 Jy in SI (W m^-2 Hz^-1).
_JY = 1.0e-26


def planck_jy(wavelength: np.ndarray, temperature: float) -> np.ndarray:
    """The Planck function ``B_nu(T)`` in Jy/sr, for *wavelength* in micron.

    Written in frequency because the flux unit is Jy — a per-frequency density
    — so evaluating ``B_lambda`` and converting back would only add a
    ``lambda**2`` round trip and its rounding.

    The exponent is evaluated with :func:`numpy.expm1`, which is what keeps the
    Rayleigh-Jeans tail honest: at long wavelengths ``h*nu/kT`` is tiny and
    ``exp(x) - 1`` loses every significant digit to cancellation, while
    ``expm1(x)`` does not.
    """
    grid = np.asarray(wavelength, dtype=DTYPE)
    if not np.all(grid > 0.0):
        raise ValueError("the Planck function needs strictly positive wavelengths (micron).")
    if not np.isfinite(temperature) or temperature <= 0.0:
        raise ValueError(f"blackbody temperature must be finite and positive, got {temperature!r}.")
    frequency = _C_MICRON_HZ / grid
    exponent = _H_OVER_K * frequency / float(temperature)
    # Guard the overflow end explicitly: on the Wien side expm1 overflows to
    # inf and the ratio would be a nan rather than the correct zero.
    with np.errstate(over="ignore"):
        denominator = np.expm1(exponent)
    radiance = _TWO_H_OVER_C2 * frequency**3 / denominator
    return np.where(np.isfinite(radiance), radiance, 0.0) / _JY


class _SpectralModel(Model):
    """Shared plumbing: the grid buffer, the channels, and the template cache.

    One declaration, one set of parameters, and one emitted spectrum per named
    channel. Several channels are not a redundancy: two instruments observing
    the same source may need genuinely different negotiated grids — a coarse
    one for broadband photometry, a fine one for a spectrograph — and binding
    them to separate channels is how each gets what it asked for from the same
    physics (``results_schema.md`` §4).
    """

    #: Name of the spectral-axis requirement these models answer.
    AXIS: ClassVar[str] = "spectral_axis"

    def __init__(self, wavelength: Any, *, channels: str | Sequence[str] = "default") -> None:
        grid = _to_micron(wavelength)
        if grid.ndim != 1 or grid.size == 0:
            raise ValueError(
                f"a spectral model needs a one-dimensional, non-empty wavelength grid, got shape "
                f"{grid.shape}."
            )
        names = (channels,) if isinstance(channels, str) else tuple(str(c) for c in channels)
        if not names or len(set(names)) != len(names):
            raise ValueError(
                f"a spectral model needs distinct, non-empty channel names, got {names!r}."
            )
        self.channels = names
        self.register_buffer("wavelength", grid, unit=COORDINATE_UNIT)
        self.templates: dict[str, Spectrum] = {}

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each of this model's channels, once.

        The containers are built here and refilled with ``with_values`` on
        every evaluation, so the axes are validated once rather than per draw
        and successive results share them by identity — the hot-loop contract
        ``transformations.md`` §14 asks a compiled model for.
        """
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * COORDINATE_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
        return self

    def _grid(self, channel: str, context: Mapping[str, Any]) -> np.ndarray:
        template = self.templates.get(channel)
        if template is None:
            return np.asarray(context["wavelength"], dtype=DTYPE)
        return template.spectral_axis.values

    def _emit(self, channel: str, grid: np.ndarray, flux: np.ndarray) -> Spectrum:
        values = np.asarray(flux, dtype=DTYPE)
        template = self.templates.get(channel)
        if template is None:
            return Spectrum(grid * COORDINATE_UNIT, values, unit=FLUX_UNIT)
        return template.with_values(values)

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        raise NotImplementedError

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        emitted = {
            channel: self._emit(channel, grid, self._flux(grid, context))
            for channel in self.channels
            for grid in (self._grid(channel, context),)
        }
        return ModelResult(emitted)


class BlackBody(_SpectralModel):
    """``F_nu = scale * B_nu(T)`` — an isothermal blackbody in Jy.

    ``scale`` absorbs the solid angle (and any distance dilution): it is the
    dimensionless factor turning the Planck radiance into an observed flux
    density, so a fit of a stellar photosphere varies ``scale`` and
    ``temperature`` together.

    Parameters
    ----------
    wavelength
        The model's own grid, micron (bare, or a :class:`~astropy.units.Quantity`).
    temperature
        Kelvin. A prior to fit it, a number to hold it fixed.
    scale
        Dimensionless multiplier. A prior to fit it, a number to hold it fixed.
    channels
        Name (or names) of the channel the emitted
        :class:`~ampere.core.Spectrum` appears under.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        temperature: Any = 1000.0,
        scale: Any = 1.0,
        channels: str | Sequence[str] = "default",
    ) -> None:
        super().__init__(wavelength, channels=channels)
        self.register_parameter(as_parameter("temperature", temperature, unit=u.K))
        self.register_parameter(as_parameter("scale", scale))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return float(context["scale"]) * planck_jy(grid, float(context["temperature"]))


class ModifiedBlackBody(_SpectralModel):
    """``F_nu = scale * (lambda_0 / lambda)**beta * B_nu(T)`` — optically thin dust.

    The standard greybody: an emissivity rising as ``nu**beta`` multiplying the
    Planck function. ``reference_wavelength`` is where the emissivity is unity,
    so ``scale`` stays interpretable as the flux the source would have if it
    radiated as a pure blackbody at that wavelength — without it, ``scale`` and
    ``beta`` are degenerate in a way that makes the posterior hard to read.

    Parameters
    ----------
    wavelength
        The model's own grid, micron.
    temperature
        Kelvin. A prior to fit it, a number to hold it fixed.
    beta
        Emissivity index. Typically 1 to 2 for interstellar dust.
    scale
        Dimensionless multiplier, at ``reference_wavelength``.
    reference_wavelength
        Micron; a **buffer**, not a parameter — one would not put a prior on
        the definition of one's own normalisation.
    channels
        Name (or names) of the emitted channel.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        temperature: Any = 100.0,
        beta: Any = 1.5,
        scale: Any = 1.0,
        reference_wavelength: float = 250.0,
        channels: str | Sequence[str] = "default",
    ) -> None:
        super().__init__(wavelength, channels=channels)
        reference = float(_to_micron(reference_wavelength))
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError(
                f"reference_wavelength must be finite and positive (micron), got {reference!r}."
            )
        self.register_buffer("reference_wavelength", reference, unit=COORDINATE_UNIT)
        self.register_parameter(as_parameter("temperature", temperature, unit=u.K))
        self.register_parameter(as_parameter("beta", beta))
        self.register_parameter(as_parameter("scale", scale))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        emissivity = (float(context["reference_wavelength"]) / grid) ** float(context["beta"])
        return float(context["scale"]) * emissivity * planck_jy(grid, float(context["temperature"]))


class PowerLaw(_SpectralModel):
    """``F_nu = norm * (lambda / lambda_ref)**index`` — a spectral power law.

    Written in wavelength, which is the axis the v1 slice's containers use, so
    a spectral index quoted in frequency changes sign. ``reference_wavelength``
    is a buffer for the same reason as in :class:`ModifiedBlackBody`: it pins
    what ``norm`` means.

    Parameters
    ----------
    wavelength
        The model's own grid, micron.
    norm
        Flux at ``reference_wavelength``, Jy.
    index
        The power-law index in wavelength.
    reference_wavelength
        Micron. A buffer.
    channels
        Name (or names) of the emitted channel.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        norm: Any = 1.0,
        index: Any = -1.0,
        reference_wavelength: float = 1.0,
        channels: str | Sequence[str] = "default",
    ) -> None:
        super().__init__(wavelength, channels=channels)
        reference = float(_to_micron(reference_wavelength))
        if not np.isfinite(reference) or reference <= 0.0:
            raise ValueError(
                f"reference_wavelength must be finite and positive (micron), got {reference!r}."
            )
        self.register_buffer("reference_wavelength", reference, unit=COORDINATE_UNIT)
        self.register_parameter(as_parameter("norm", norm))
        self.register_parameter(as_parameter("index", index))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        ratio = grid / float(context["reference_wavelength"])
        return float(context["norm"]) * ratio ** float(context["index"])


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=DTYPE)
    return np.asarray(coordinates, dtype=DTYPE)
