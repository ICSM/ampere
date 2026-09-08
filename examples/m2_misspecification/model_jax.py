"""The toy model in ``jax.numpy``: the same physics, traceably.

The mirror of :mod:`.model_torch`, and imported under the same discipline —
only when the study runs on the jax backend, never from the package's
``__init__``.

:mod:`ampere.backends.jax.problem` composes a problem's differentiable
log-density from its parts' native surfaces and finds a model's by
``hasattr(model, "flux")``, so this hand-written class lowers exactly as
``ampere.backends.jax.PowerLaw`` does. ``flux`` is pure and traceable;
``evaluate`` is the contract surface and converts to numpy, which is where a
gradient stops and why the two exist separately.

float64 is policy rather than a default here (``lowering.md`` §10.2): this
module's classes refuse at construction if ``jax_enable_x64`` is off, by
calling :func:`ampere.backends.jax.require_x64` exactly as the shipped models
do. The *application* — the study driver, or a test module — turns it on with
:func:`ampere.backends.jax.configure_x64`.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np

from ampere.backends.jax import BACKEND, require_x64
from ampere.core import DTYPE, ChannelRequirements, Model, ModelResult, Spectrum

from .model import (
    FLUX_UNIT,
    LINE1_CENTRE,
    LINE1_WIDTH,
    LINE2_CENTRE,
    LINE2_WIDTH,
    REFERENCE_WAVELENGTH,
    WAVELENGTH_UNIT,
    _as_parameter,
    _channels,
    _declared,
    _to_micron,
)

__all__ = ["AbsorptionLines", "flux_jax"]


def flux_jax(wavelength: Any, A: Any, B: Any, d1: Any, d2: Any) -> jax.Array:
    """:func:`examples.m2_misspecification.model.flux_at` in ``jax.numpy``.

    A free function so a test can ``jax.grad`` it without composing a problem.
    Nothing branches on a traced value and nothing indexes by one, so it jits
    and vmaps as any pure function does.
    """
    grid = jnp.asarray(wavelength, dtype=jnp.float64)
    continuum = jnp.asarray(A) + jnp.asarray(B) * (grid - REFERENCE_WAVELENGTH)
    first = jnp.exp(-0.5 * ((grid - LINE1_CENTRE) / LINE1_WIDTH) ** 2)
    second = jnp.exp(-0.5 * ((grid - LINE2_CENTRE) / LINE2_WIDTH) ** 2)
    return continuum * (1.0 - jnp.asarray(d1) * first - jnp.asarray(d2) * second)


class AbsorptionLines(Model):
    """The toy model on the jax path.

    Parameters
    ----------
    wavelength
        The model's own grid, micron.
    A, B, d1, d2
        As :class:`examples.m2_misspecification.model.AbsorptionLines`.
    channels
        Name of the emitted channel.
    """

    AXIS: ClassVar[str] = "spectral_axis"
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        wavelength: Any,
        *,
        A: Any = None,
        B: Any = None,
        d1: Any = None,
        d2: Any = None,
        channels: str | Sequence[str] = "default",
    ) -> None:
        require_x64(f"a jax {type(self).__name__}")
        grid = _to_micron(wavelength)
        self.channels = _channels(channels)
        self.register_buffer("wavelength", grid, unit=WAVELENGTH_UNIT)
        for name, given in _declared(A=A, B=B, d1=d1, d2=d2).items():
            self.register_parameter(_as_parameter(name, given))
        self.templates: dict[str, Spectrum] = {}
        #: An ordinary array leaf, never an ``eqx.field(static=True)``: a
        #: buffer's *contents* in the pytree structure would enter every JIT
        #: cache key (``lowering.md`` §7).
        self.grids: dict[str, jax.Array] = {
            channel: jnp.asarray(grid, dtype=jnp.float64) for channel in self.channels
        }

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each channel, once."""
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * WAVELENGTH_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
            self.grids[channel] = jnp.asarray(grid, dtype=jnp.float64)
        return self

    # -- the native surface --------------------------------------------------

    def grid(self, channel: str) -> jax.Array:
        """The jax coordinates *channel* is evaluated on, micron."""
        return self.grids[channel]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> jax.Array:
        """This model's flux on *channel*, Jy, as a traceable jax array."""
        context = self.context(values)
        return flux_jax(
            self.grid(channel), context["A"], context["B"], context["d1"], context["d2"]
        )

    # -- the contract surface ------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        emitted = {}
        for channel in self.channels:
            flux = np.asarray(self.flux(channel, values), dtype=DTYPE)
            template = self.templates.get(channel)
            emitted[channel] = (
                Spectrum(np.asarray(self.grid(channel)) * WAVELENGTH_UNIT, flux, unit=FLUX_UNIT)
                if template is None
                else template.with_values(flux)
            )
        return ModelResult(emitted)
