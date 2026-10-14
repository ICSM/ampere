"""The forward model: an excitation ladder rendered as Gaussian lines (W7.9).

``EmissionLadder`` emits a spectrum — a flux density in Jy on a micron grid —
because that is what :class:`~ampere.backends.reference.LineIntegration`
integrates, and because a model that emitted the line fluxes directly would
skip the instrument chain the item exists to exercise. Each transition *k*
carries a rest wavelength, an upper-level energy ``E_k`` (in kelvin, ``E/k``)
and a line strength ``S_k`` as **buffers**; the integrated flux of line *k* is

.. math::

    F_k = N\\,S_k\\,\\frac{e^{-E_k/T}}{Z(T)}, \\qquad
    Z(T) = 1 + \\sum_j e^{-E_j/T},

in units of :data:`FLUX_SCALE` (``1e-16 W m⁻²``), for one excitation
temperature ``T`` and one column-like normalisation ``N`` — the two free
parameters. A second component (``hot_temperature``, ``hot_column``) adds a
second Boltzmann term with its own partition function; the truth has one, and
the study's misspecified arms leave it out.

Each line is rendered as a Gaussian in velocity of fixed width (the ``width``
buffer, km/s) normalised so that its integral over frequency *is* ``F_k`` —
so a window several widths wide integrates the line back to its flux, to the
quadrature's accuracy. There is no continuum: the observer's table is
continuum-subtracted, and so the instrument is ``continuum="none"``.

**The physics is a toy.** The energies and strengths in :mod:`.generators`
are an illustrative ladder, not atomic data; a real fit to nebular lines
would carry each ion's level populations and its abundance. The inference
pattern — a ladder model, a line table, a misspecified excitation structure
showing as a correlated residual along rest wavelength — is the point.

Written once, against :class:`ampere.core.PortableModel` (W7.13): the torch
and jax twins in :mod:`.model_torch` and :mod:`.model_jax` are one line each.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ampere.core import Parameter, PortableModel

__all__ = ["FLUX_SCALE", "EmissionLadder"]

#: The unit of ``N`` (and of every ``F_k`` before rendering), W m⁻².
FLUX_SCALE = 1.0e-16
#: The speed of light, km/s.
_C_KMS = 299792.458
#: The speed of light, micron per second (``nu = c / lambda``).
_C_MICRON_S = 2.99792458e14
#: One jansky in W m⁻² Hz⁻¹.
_JANSKY = 1.0e-26


def _parameter(name: str, given: Any) -> Parameter:
    """A frozen ``scipy.stats`` distribution fits *name*; a number fixes it."""
    if hasattr(given, "ppf"):
        return Parameter(name, given)
    return Parameter(name, value=float(given), fixed=True)


class EmissionLadder(PortableModel):
    """One (or two) Boltzmann components over a line list, as a spectrum in Jy.

    Parameters
    ----------
    wavelength
        The model's own placeholder grid, micron — replaced by the
        line-integration step's tabulation the moment negotiation runs.
    rest, energy, strength
        Per transition: rest wavelength (micron), upper-level energy (K) and
        line strength (dimensionless). Buffers.
    width
        The lines' Gaussian velocity dispersion, km/s. A buffer.
    temperature, column
        The component every arm fits: a ``scipy.stats`` distribution fits it,
        a number fixes it.
    hot_temperature, hot_column
        The second component, as for the first; ``None`` (both) builds the
        one-component model. Two components are what the truth has.
    channels
        The channel name.
    """

    def __init__(
        self,
        wavelength: Any,
        *,
        rest: Any,
        energy: Any,
        strength: Any,
        width: float,
        temperature: Any,
        column: Any,
        hot_temperature: Any = None,
        hot_column: Any = None,
        channels: str | Sequence[str] = "lines",
    ) -> None:
        super().__init__(wavelength, channels=channels)
        rest_array = np.asarray(rest, dtype=np.float64)
        self.register_buffer("rest", rest_array)
        self.register_buffer("energy", np.asarray(energy, dtype=np.float64))
        self.register_buffer("strength", np.asarray(strength, dtype=np.float64))
        self.register_buffer("width", np.asarray(float(width), dtype=np.float64))
        self.register_parameter(_parameter("temperature", temperature))
        self.register_parameter(_parameter("column", column))
        if (hot_temperature is None) != (hot_column is None):
            raise ValueError(
                "EmissionLadder's second component needs both hot_temperature and hot_column, "
                "or neither."
            )
        self.components = 1 if hot_temperature is None else 2
        if self.components == 2:
            self.register_parameter(_parameter("hot_temperature", hot_temperature))
            self.register_parameter(_parameter("hot_column", hot_column))

    def line_fluxes(self, context: Mapping[str, Any]) -> Any:
        """``F_k`` in units of :data:`FLUX_SCALE`, from one or two components."""
        ops = self.ops
        energy, strength = context["energy"], context["strength"]

        def component(temperature: Any, column: Any) -> Any:
            boltzmann = ops.exp(-energy / temperature)
            return column * strength * boltzmann / (1.0 + ops.sum(boltzmann))

        fluxes = component(context["temperature"], context["column"])
        if self.components == 2:
            fluxes = fluxes + component(context["hot_temperature"], context["hot_column"])
        return fluxes

    def _flux(self, grid: Any, context: Mapping[str, Any]) -> Any:
        ops = self.ops
        rest = context["rest"]
        width = context["width"]
        # Velocity offset of every grid point from every line, km/s: (n, k).
        velocity = _C_KMS * (grid[:, None] - rest[None, :]) / rest[None, :]
        # A Gaussian in frequency of dispersion nu_k * width / c, normalised to
        # unit integral over nu, times F_k in W m^-2, in Jy.
        sigma_nu = (_C_MICRON_S / rest) * width / _C_KMS
        profile = ops.exp(-0.5 * (velocity / width) ** 2) / (np.sqrt(2.0 * np.pi) * sigma_nu)
        fluxes = self.line_fluxes(context) * (FLUX_SCALE / _JANSKY)
        return ops.sum(profile * fluxes[None, :], axis=1)
