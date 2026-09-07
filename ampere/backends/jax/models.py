"""The native spectral models, in ``jax.numpy``: blackbody, modified blackbody, power law.

``DEVELOPMENT_PLAN.md`` §5's Phase 2 list, on the jax path. The *declarations*
are identical to the reference backend's — same class names, same parameter
names, same units, same buffers — because that is what makes the conformance
suite's cross-backend rows a real comparison rather than two spellings of one
implementation. What differs is the arithmetic: every flux here is computed in
``jax.numpy``, so it can be differentiated, jitted and (in slice 2) vmapped.

Two evaluation surfaces, and why there are two
----------------------------------------------
:meth:`~_SpectralModel.flux` is the **native** surface: coordinates and
parameter values in, a ``jax`` array out, nothing else. It is pure, traceable,
and is what :mod:`ampere.backends.jax.problem` composes into the differentiable
log-density that NUTS consumes.

:meth:`~_SpectralModel.evaluate` is the **contract** surface
(``transformations.md``): it wraps ``flux`` in the ``ampere.core`` containers a
``ModelResult`` is made of. Those containers convert their values with
``numpy.asarray`` (``results_schema.md``'s ``_as_array``), so a jax array
becomes a numpy one at that boundary and **a gradient does not survive it**.
That is a fact about the frozen container contract rather than a choice made
here; it is why the native surface exists, and it is recorded in this item's
report as a finding rather than worked around silently.

``DIFFERENTIABLE = True`` refers to the capability the *problem* gains — "whether
``log_prob`` admits a gradient, i.e. whether NUTS/HMC and gradient-based VI or
optimisation can run", which is how ``ampere.core.Capabilities`` defines the
flag — and with this backend's lowered path they can.

Everything else follows the reference backend's rules, for the same reasons:
units are converted **once at construction** and never in the hot loop
(``DEVELOPMENT_PLAN.md`` §7); the wavelength grid is a **buffer**, since one
would never put a prior on it (``parameters.md`` §10); each physical quantity
is accepted as a prior, a number or a ready-made ``Parameter``; and
``compile_for`` builds one container per requested channel and refills it with
``with_values`` so successive results share their axes by identity
(``transformations.md`` §14).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.constants as const
import astropy.units as u
import jax
import jax.numpy as jnp
import numpy as np

from ampere.core import (
    DTYPE,
    ChannelRequirements,
    Model,
    ModelResult,
    Spectrum,
)

from ._config import BACKEND, require_x64
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
# (DEVELOPMENT_PLAN.md §7). Identical values to the reference backend's, read
# from astropy in exactly the same way, so a cross-backend disagreement can
# only come from the arithmetic and never from the constants.
_TWO_H_OVER_C2 = float((2.0 * const.h / const.c**2).to_value(u.J * u.s**3 / u.m**2))
_H_OVER_K = float((const.h / const.k_B).to_value(u.s * u.K))
_C_MICRON_HZ = float(const.c.to_value(u.micron * u.Hz))
#: 1 Jy in SI (W m^-2 Hz^-1).
_JY = 1.0e-26


def planck_jy(wavelength: Any, temperature: Any) -> jax.Array:
    """The Planck function ``B_nu(T)`` in Jy/sr, for *wavelength* in micron.

    Written in frequency because the flux unit is Jy — a per-frequency density
    — so evaluating ``B_lambda`` and converting back would only add a
    ``lambda**2`` round trip and its rounding.

    :func:`jax.numpy.expm1` is what keeps the Rayleigh-Jeans tail honest: at
    long wavelengths ``h*nu/kT`` is tiny and ``exp(x) - 1`` loses every
    significant digit to cancellation.

    **No validation, deliberately.** The reference implementation raises for a
    non-positive wavelength or temperature; this one cannot, because it must
    survive tracing — a Python ``if`` on a traced value is not a branch, it is
    an error. The Wien-side overflow (where ``expm1`` returns ``inf`` and the
    ratio would be ``nan``) is handled with :func:`jax.numpy.where`, which is
    the traced equivalent of the reference's ``errstate`` guard, and a
    non-positive temperature produces a non-finite flux that the likelihood's
    §4.5 failure path turns into ``-inf`` with a recorded reason.
    """
    grid = jnp.asarray(wavelength, dtype=jnp.float64)
    frequency = _C_MICRON_HZ / grid
    exponent = _H_OVER_K * frequency / jnp.asarray(temperature, dtype=jnp.float64)
    denominator = jnp.expm1(exponent)
    radiance = _TWO_H_OVER_C2 * frequency**3 / denominator
    return jnp.where(jnp.isfinite(radiance), radiance, 0.0) / _JY


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

    #: The four capability flags (``DEVELOPMENT_PLAN.md`` §4.5, W2.12),
    #: declared rather than inherited: every piece of this backend says what it
    #: is for itself.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: **True since slice 2** (W2.5). The flag means one thing on this
    #: backend: ``jax.vmap`` over the realised density
    #: (:meth:`~ampere.backends.jax.problem.LoweredProblem.log_prob_unconstrained_batched`)
    #: evaluates a stack of parameter vectors in one call, and it is measured
    #: rather than asserted -- ``tests/backends/test_jax.py`` compares a vmapped
    #: density against the same density in a loop. It is true here because every
    #: operation in this class is whole-array ``jax.numpy``: nothing branches on
    #: a value, nothing indexes by one, so vmap maps it as it maps any pure
    #: function. ``QuasisepGP`` is the one part of this backend that still says
    #: False, and says why.
    BATCHABLE: ClassVar[bool] = True
    #: CPU by default and never auto-detected (``architecture.md`` §5): a
    #: machine with a GPU present must not silently take a different code path
    #: from CI. GPU placement is slice 2's, with an explicit ``device_put``.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, wavelength: Any, *, channels: str | Sequence[str] = "default") -> None:
        require_x64(f"a jax {type(self).__name__}")
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
        #: The evaluation grids as jax arrays, one per channel, converted once.
        #: A buffer is an ordinary array leaf here — never an
        #: ``eqx.field(static=True)``, which would put its *contents* into the
        #: pytree structure and so into every JIT cache key (``lowering.md``
        #: §7). Changing one of these changes leaves, never structure.
        self.grids: dict[str, jax.Array] = {
            channel: jnp.asarray(grid, dtype=jnp.float64) for channel in names
        }

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each of this model's channels, once.

        The containers are built here and refilled with ``with_values`` on
        every evaluation, so the axes are validated once rather than per draw
        and successive results share them by identity — the hot-loop contract
        ``transformations.md`` §14 asks a compiled model for. The jax copy of
        each grid is rebuilt at the same time, so the native path and the
        contract path never work on different coordinates.
        """
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * COORDINATE_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
            self.grids[channel] = jnp.asarray(grid, dtype=jnp.float64)
        return self

    # -- the native surface -------------------------------------------------

    def grid(self, channel: str) -> jax.Array:
        """The jax coordinates *channel* is evaluated on, micron."""
        return self.grids[channel]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's flux on *channel*, in Jy, as a jax array.

        Pure and traceable: the surface :mod:`ampere.backends.jax.problem`
        composes, and the one a gradient actually passes through.
        """
        return self._flux(self.grid(channel), self.context(values))

    def _flux(self, grid: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        raise NotImplementedError

    # -- the contract surface ----------------------------------------------

    def _emit(self, channel: str, flux: jax.Array) -> Spectrum:
        template = self.templates.get(channel)
        if template is None:
            return Spectrum(
                np.asarray(self.grids[channel]) * COORDINATE_UNIT,
                np.asarray(flux),
                unit=FLUX_UNIT,
            )
        return template.with_values(np.asarray(flux))

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        emitted = {
            channel: self._emit(channel, self._flux(self.grids[channel], context))
            for channel in self.channels
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

    def _flux(self, grid: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        return jnp.asarray(context["scale"]) * planck_jy(grid, context["temperature"])


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

    def _flux(self, grid: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        emissivity = (jnp.asarray(context["reference_wavelength"]) / grid) ** jnp.asarray(
            context["beta"]
        )
        return jnp.asarray(context["scale"]) * emissivity * planck_jy(grid, context["temperature"])


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

    def _flux(self, grid: jax.Array, context: Mapping[str, Any]) -> jax.Array:
        ratio = grid / jnp.asarray(context["reference_wavelength"])
        return jnp.asarray(context["norm"]) * ratio ** jnp.asarray(context["index"])


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity.

    numpy rather than jax: this runs at *construction*, where the value is
    validated, stored in a core ``Buffer`` and compared for model identity. The
    jax copy is made once, afterwards, in :meth:`_SpectralModel.__init__`.
    """
    if isinstance(coordinates, u.Quantity):
        return np.asarray(coordinates.to_value(COORDINATE_UNIT), dtype=DTYPE)
    return np.asarray(coordinates, dtype=DTYPE)
