"""Astrometric time series on the numpy path.

Phase 4's second worked modality, following the interferometry page's own
template (``docs/source/interferometry.rst``) — chosen by Peter 2026-09-11
("this will be an essential test of the code") to test the template's claim
that a modality nobody wrote the contracts for follows it with no new
contract surface. The design sketch is
``docs/design/modalities/astrometric_timeseries.md``; this module is that
sketch's composition, shipped, under the placement D1 ruled for
interferometry (``WORK_ITEMS.md``'s Phase 4 header): the kind lives in
``ampere.core`` (:class:`~ampere.core.TimeSeries` already did, needing no
amendment at all — it is one ordered ``time`` axis, exactly what this
modality wants), and the arithmetic lives here, in a per-backend
``astrometry.py``, with no grouping namespace.

The composition this module exists to support, in one picture::

    ReflexOrbit (Model)
      ├── channel "ra"  : TimeSeries    (mas offset, time in days)
      │        └── Instrument "astrom_ra"  [EpochSample]
      │                 → TimeSeries, Likelihood(GaussianFamily(), ...)
      └── channel "dec" : TimeSeries
               └── Instrument "astrom_dec" [EpochSample]
                        → TimeSeries, Likelihood(GaussianFamily(), ...)

``results_schema.md`` §15.2's rule for a vector-valued observable ("one value
array per container", applied to Stokes I/Q/U/V or a model's flux and its
optical depth) is the same shape of decision for a 2D sky position: RA and
Dec are two :class:`~ampere.core.TimeSeries` channels sharing one model
evaluation, not one container with a ``(N, 2)`` value array (which
:class:`~ampere.core.FunctionSamples` does not support in any case). ``period``
and ``phase`` are shared *at the language level* by construction — both
channels read the same ``values["period"]``/``values["phase"]`` inside one
``evaluate()`` call — so no ``Tie`` or ``shared_as`` is needed for the
physical tying between the two coordinates.

**What this module deliberately does not do.** The sketch's own reflex model
is a circular orbit projected onto the sky (a periodic wobble plus a linear
proper-motion drift), not a full eccentric Keplerian orbit with Thiele-Innes
constants, an epoch of periastron or a solved Kepler's equation — the sketch
"closes without needing anything new" with exactly this simpler form, and a
richer parameterisation is not what the design sketch asked the template to
prove. See ``docs/source/astrometry.rst``'s closing section for this as a
question carried to W4.8. **A genuinely joint 2-vector GP over the two
coordinates** (a single correlated noise process over RA and Dec together,
rather than two independent per-coordinate processes) is the sketch's own
named gap, dispositioned at the freeze as Phase 5's ``JointGP`` slot; this
module's likelihood is two independent time-domain flexible processes, one
per channel, which is what the merged contracts can express today.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import astropy.units as u
import numpy as np

from ampere.core import (
    AxisRequirement,
    ChannelRequirements,
    DTYPE,
    Model,
    ModelResult,
    TimeSeries,
    Transformation,
    TransformationError,
)

from ._declare import as_parameter

__all__ = [
    "OFFSET_UNIT",
    "PM_UNIT",
    "TIME_UNIT",
    "EpochSample",
    "ReflexOrbit",
]

#: The unit an epoch is given in throughout this module.
TIME_UNIT = u.day
#: The unit a sky-offset value (an RA or Dec residual from the phase centre)
#: is emitted in.
OFFSET_UNIT = u.mas
#: The unit a proper motion is given in.
PM_UNIT = u.mas / u.yr

#: Julian days per year, for the proper-motion term. Written out rather than
#: taken from astropy's own constant, for the same reason ``oracles.py``
#: transcribes its own formulae: a one-line conversion factor is exactly the
#: kind of thing that should agree by calculation, not by import.
DAYS_PER_YEAR = 365.25


def _to_unit(values: Any, unit: u.UnitBase) -> np.ndarray:
    """*values* as a bare float64 array in *unit*, whether or not they arrive as a Quantity."""
    if isinstance(values, u.Quantity):
        return np.asarray(values.to_value(unit), dtype=DTYPE)
    return np.asarray(values, dtype=DTYPE)


class _Step(Transformation):
    """Shared plumbing: read a declared buffer as a plain array.

    The same shape as ``instrument.py``'s private base and
    ``interferometry.py``'s own copy, and separate from both on purpose — a
    buffer's name may not shadow a class attribute, so this is the readable
    spelling of ``self.buffers[name].value``.
    """

    BACKEND: ClassVar[str] = "reference"

    def _data(self, name: str) -> np.ndarray:
        return np.asarray(self.buffers[name].value, dtype=DTYPE)


class EpochSample(_Step):
    """Pin a model channel onto the observed epochs: no arithmetic at all.

    ``transformations.md`` §10's standard-library table names exactly this
    row: "Epoch sampling | ``TimeSeries`` -> ``TimeSeries`` | none | ``points``
    at the observed epochs." It is the epoch-sampling instrument the design
    sketch's §2 asks for, kind-preserving (``PRODUCES`` is left at its default
    of ``None``) and free of parameters: the observation times are a
    **buffer taken from the observed container**, never recomputed and never
    fitted — the same rule ``FourierSample.from_observed`` keeps for a
    baseline coordinate, applied to a time coordinate instead.

    Because the buffer's job is *only* to publish the requirement, ``apply``
    is the identity: once negotiation has adopted the published ``points=``
    requirement, the container it is handed is already on exactly these
    epochs, so there is nothing left to do (:meth:`~ampere.core.transform.
    Transformation.__call__` still checks the kind and the mask on the way
    through).

    Parameters
    ----------
    epochs
        The observation times, days (or a :class:`~astropy.units.Quantity`).
        Strictly increasing, as a :class:`~ampere.core.TimeSeries`'s own time
        axis always is.
    label
        Component label for this step within a chain.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (TimeSeries,)

    def __init__(self, epochs: Any, *, label: str | None = None) -> None:
        super().__init__(label=label)
        grid = _to_unit(epochs, TIME_UNIT).reshape(-1)
        if grid.size == 0:
            raise TransformationError("EpochSample needs at least one observed epoch.")
        if grid.size > 1 and not bool(np.all(np.diff(grid) > 0.0)):
            raise TransformationError(
                "EpochSample's epochs must be strictly increasing; a TimeSeries's own time axis "
                "is strictly increasing, so build one with EpochSample.from_observed(container) "
                "rather than assembling the array by hand."
            )
        self.register_buffer("epochs", grid, unit=TIME_UNIT)

    @classmethod
    def from_observed(cls, container: TimeSeries, *, label: str | None = None) -> EpochSample:
        """Take the epochs from the observed container. **The supported route.**

        ``transformations.md`` §10 (gap I-2's rule, restated for this
        modality): a step reproducing observed coordinates takes them from
        the observed container rather than recomputing them, so
        ``check_alignment``'s exact comparison is a comparison of two
        independent things rather than a tautology.
        """
        return cls(container.time.values, label=label)

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """``points=`` at the exact epochs — coordinates, not a density.

        Astrometric epochs are irregular by nature (scheduling, weather,
        visibility windows), so ``points=`` is the right shape of requirement
        here, not ``intervals`` plus a cadence — the design sketch's own
        reasoning (§2).
        """
        return (AxisRequirement("time", unit=TIME_UNIT, points=self._data("epochs")),)

    def apply(self, samples: TimeSeries, values: Any) -> TimeSeries:
        return samples


class ReflexOrbit(Model):
    """A linear proper motion plus a periodic reflex wobble, per coordinate.

    The model the sketch's whole composition is built to prove: a source
    whose sky position drifts linearly (proper motion) and wobbles
    periodically (an unseen companion's reflex motion, to leading order a
    projected circular orbit) emits two :class:`~ampere.core.TimeSeries`
    channels, ``"ra"`` and ``"dec"``, sharing one set of orbital parameters —
    ``period`` and ``phase`` are read from the same ``values`` mapping inside
    one :meth:`evaluate` call, so the two channels are tied by construction
    and need no explicit ``Tie``.

    .. math::

        \\Delta\\alpha(t) = \\mu_\\alpha \\frac{t}{365.25}
            + a_\\alpha \\sin\\!\\left(\\frac{2\\pi t}{P} + \\phi\\right)

        \\Delta\\delta(t) = \\mu_\\delta \\frac{t}{365.25}
            + a_\\delta \\cos\\!\\left(\\frac{2\\pi t}{P} + \\phi\\right)

    with ``t`` in days since the reference epoch. This is a **circular**,
    node-aligned reflex orbit (the module docstring says why this is the
    sketch's own choice rather than an omission): ``amp_ra``/``amp_dec`` carry
    the projection of the photocentre orbit's semi-major axis onto the two
    sky coordinates, folding inclination and node into one number per
    coordinate rather than exposing them separately.

    Parameters
    ----------
    time
        This model's own time grid, days — used when :meth:`compile_for` is
        not called (a direct, unnegotiated evaluation) and when a channel's
        requirement carries no ``time`` axis. Strictly increasing.
    pmra, pmdec
        Proper motion, mas/yr.
    period
        Orbital period, days.
    phase
        Phase at the reference epoch, radians.
    amp_ra, amp_dec
        Sky-projected reflex semi-amplitude in right ascension and
        declination, mas.
    """

    BACKEND: ClassVar[str] = "reference"

    #: The two channels this model always emits — ``results_schema.md``
    #: §15.2's "one value array per container" rule applied to a 2D sky
    #: position, so this is not a constructor argument the way an
    #: interferometric source model's ``channels=`` is: RA and Dec are this
    #: model's own physical outputs, not a caller's choice of channel name.
    CHANNELS: ClassVar[tuple[str, str]] = ("ra", "dec")

    def __init__(
        self,
        time: Any,
        *,
        pmra: Any = 0.0,
        pmdec: Any = 0.0,
        period: Any = 400.0,
        phase: Any = 0.5,
        amp_ra: Any = 0.5,
        amp_dec: Any = 0.3,
    ) -> None:
        grid = _to_unit(time, TIME_UNIT).reshape(-1)
        if grid.size == 0 or not bool(np.all(np.diff(grid) > 0.0)):
            raise ValueError(
                f"ReflexOrbit needs a strictly increasing, non-empty time grid (days), got "
                f"{grid.size} point(s)."
            )
        self.register_buffer("time", grid, unit=TIME_UNIT)
        self.register_parameter(as_parameter("pmra", pmra, unit=PM_UNIT))
        self.register_parameter(as_parameter("pmdec", pmdec, unit=PM_UNIT))
        self.register_parameter(as_parameter("period", period, unit=TIME_UNIT))
        self.register_parameter(as_parameter("phase", phase))
        self.register_parameter(as_parameter("amp_ra", amp_ra, unit=OFFSET_UNIT))
        self.register_parameter(as_parameter("amp_dec", amp_dec, unit=OFFSET_UNIT))
        self.templates: dict[str, TimeSeries] = {}

    # -- negotiation -----------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> ReflexOrbit:
        """Adopt each channel's own negotiated epochs, independently.

        RA and Dec need not share one epoch set — a caller free to bind two
        different ``EpochSample`` chains to the two channels gets exactly
        that, computed correctly, because each channel's offsets are a
        function of *its own* ``time`` array (a strict generalisation of the
        design sketch's own worked example, which assumes one shared epoch
        set for both coordinates and reads it off the ``"ra"`` channel
        alone).
        """
        for channel in self.CHANNELS:
            asked = requirements.get(channel)
            if asked is None or "time" not in asked:
                continue
            grid = _to_unit(asked["time"].coordinates(), TIME_UNIT)
            self.templates[channel] = TimeSeries(
                grid * TIME_UNIT, np.zeros(grid.size, dtype=DTYPE) * OFFSET_UNIT
            )
        return self

    # -- evaluation --------------------------------------------------------

    def _data(self, name: str) -> np.ndarray:
        return np.asarray(self.buffers[name].value, dtype=DTYPE)

    def _time(self, channel: str) -> np.ndarray:
        template = self.templates.get(channel)
        if template is not None:
            return np.asarray(template.time.values, dtype=DTYPE)
        return self._data("time")

    def _emit(self, channel: str, t: np.ndarray, values: np.ndarray) -> TimeSeries:
        template = self.templates.get(channel)
        if template is None:
            return TimeSeries(t * TIME_UNIT, values * OFFSET_UNIT)
        return template.with_values(values)

    def _offset(self, channel: str, t: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        cycle = 2.0 * np.pi * t / float(context["period"]) + float(context["phase"])
        if channel == "ra":
            return context["pmra"] * t / DAYS_PER_YEAR + context["amp_ra"] * np.sin(cycle)
        return context["pmdec"] * t / DAYS_PER_YEAR + context["amp_dec"] * np.cos(cycle)

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        emitted = {
            channel: self._emit(channel, t, self._offset(channel, t, context))
            for channel in self.CHANNELS
            for t in (self._time(channel),)
        }
        return ModelResult(emitted)
