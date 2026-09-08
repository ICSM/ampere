"""Transformations, instruments and requirements negotiation (``§4.3``).

A **model** produces a :class:`~ampere.core.results_schema.ModelResult`: named
channels of coordinate-indexed function samples, in whatever space the physics
is natural in. An **instrument** turns one of those channels into a prediction
of what a particular detector would have recorded. This module is the contract
for that second half.

The factoring is deliberate and is the standing decision recorded in
``docs/design/prior_art.md`` §3 and §6 (Tension 3): an ``Instrument`` is an
ordered chain of small, reusable :class:`Transformation`\\ s — line-spread
convolution, resampling, synthetic photometry, a calibration scale factor, a
response matrix — rather than 3ML's single opaque per-instrument plugin that
owns its response *and* its likelihood. That buys reuse across instruments at
the cost of requiring the physics to decompose; W1.11's modality sketches carry
the stress test for the cases where it might not.

Three things make the chain more than a list of callables:

1. **Kinds are checked when the chain is assembled**, not when it is
   evaluated. A transformation declares what it accepts and what it produces,
   an instrument declares which channel it binds, and
   :class:`~ampere.core.exceptions.CompositionError` names the mismatch at
   composition time (``DEVELOPMENT_PLAN.md`` §4.2).
2. **Transformations may carry their own nuisance parameters.** They inherit
   :class:`~ampere.core.parameter.Parameterised`, so a calibration offset is
   declared exactly like a model parameter, and a chain's parameters compose
   through :meth:`~ampere.core.parameter.ParameterSet.merge` with each step's
   label as the component.
3. **Requirements negotiation is available but never required.** An instrument
   may publish what it needs of the model's coordinates; :func:`negotiate`
   takes the union across instruments; a model may honour it in a one-off
   :meth:`Model.compile_for` step and then evaluate onto the resulting grids
   for the rest of the run. A model is free to ignore all of it, and the
   fixed-grid path stays a two-line composition.

Everything here is backend-neutral: numpy, ``astropy.units`` and stdlib only
(``architecture.md`` §3-4).

The simple path, end to end:

>>> import numpy as np, astropy.units as u, scipy.stats as st
>>> from ampere.core import Instrument, Model, ModelResult, Parameter, Spectrum, Transformation
>>> class Calibrate(Transformation):
...     ACCEPTS = (Spectrum,)
...     def __init__(self, **kwargs):
...         super().__init__(**kwargs)
...         self.register_parameter(Parameter("scale", st.lognorm(0.1)))
...     def apply(self, samples, values):
...         return samples.with_values(samples.values * self.context(values)["scale"])
>>> class Flat(Model):
...     def evaluate(self, **values):
...         return Spectrum([1.0, 2.0, 4.0] * u.um, np.ones(3) * u.Jy)
>>> instrument = Instrument([Calibrate()])
>>> instrument(Flat()(), {"calibrate.scale": 2.0}).values.tolist()
[2.0, 2.0, 2.0]
"""

from __future__ import annotations

import abc
import dataclasses
import re
import types
from collections.abc import Iterable, Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np

from .exceptions import CompositionError, TransformationError
from .parameter import ParameterMapping, ParameterSet, Parameterised, Value
from .results_schema import COORDINATE_RTOL, DEFAULT_CHANNEL, FunctionSamples, ModelResult

__all__ = [
    "COORDINATE_RTOL",
    "AxisRequirement",
    "ChannelRequirements",
    "Instrument",
    "Model",
    "Transformation",
    "negotiate",
    "propagate_mask",
]

# COORDINATE_RTOL is defined in results_schema.py and re-exported here, where
# negotiation uses it. It moved down a layer with W2.1 so that Axis.locate --
# the exact inverse of the collapsing _dedupe performs below -- shares the one
# constant rather than a copy of its value (spectrum_photometry.md Gap 1).
# Both import paths keep working.

ArrayLike = Any

_CAMEL_BOUNDARY = re.compile(r"(?<=[a-z0-9])(?=[A-Z])|(?<=[A-Z])(?=[A-Z][a-z])")


def _default_label(cls: type) -> str:
    """``SyntheticPhotometry`` -> ``synthetic_photometry``."""
    return _CAMEL_BOUNDARY.sub("_", cls.__name__).lower()


def _check_label(name: object, kind: str) -> str:
    """Labels become component names in a merge, so they must be identifiers."""
    if not isinstance(name, str) or not name.isidentifier():
        raise CompositionError(
            f"{kind} {name!r} is not usable: it must be a valid Python identifier, because it "
            f"becomes the component label under which this object's parameters are qualified "
            f"(ParameterSet.merge), and qualified names are '.'-separated identifiers."
        )
    return name


def _check_dotted_label(name: object, kind: str) -> str:
    """Channel bindings (and the instrument labels defaulting to them) may be dot-qualified.

    Channel names were relaxed to '.'-separated identifiers on 2026-09-02
    (``results_schema.md`` §17); an instrument label is provenance, never a
    merge component (``inference.md`` §4), so it follows the channel rule.
    Step labels and dataset labels stay bare — they are merge components.
    """
    if not isinstance(name, str) or not name or not all(p.isidentifier() for p in name.split(".")):
        raise CompositionError(
            f"{kind} {name!r} is not usable: it must be a valid Python identifier or a "
            f"'.'-separated sequence of them (a model may namespace grouped channels, e.g. "
            f"'co.j3_2')."
        )
    return name


def _check_unit(unit: object, what: str) -> u.UnitBase | None:
    if unit is None:
        return None
    try:
        return u.Unit(unit)
    except (TypeError, ValueError) as exc:
        raise TransformationError(f"{what} is not a usable astropy unit: {unit!r}") from exc


def _positive(value: object, what: str) -> float | None:
    if value is None:
        return None
    number = float(value)
    if not np.isfinite(number) or number <= 0.0:
        raise TransformationError(f"{what} must be a positive finite number, got {value!r}")
    return number


# ---------------------------------------------------------------------------
# Mask propagation
# ---------------------------------------------------------------------------


def propagate_mask(
    source: FunctionSamples | ArrayLike | None,
    influence: ArrayLike | None = None,
) -> np.ndarray | None:
    """Carry a mask through a transformation, conservatively.

    ``docs/design/contracts/results_schema.md`` §16 makes mask propagation an
    obligation on this contract and proposes a rule for the many-to-one case.
    This function is that rule, ratified: **an output sample is masked if any
    input sample that influences it is masked.** A masked sample carries zero
    information (results_schema §7), so an output bin that integrates over one
    is itself uninformative — and the alternative, being permissive, silently
    feeds invalid model values into a likelihood.

    Parameters
    ----------
    source
        The container whose mask is being propagated, or a boolean mask array,
        or ``None``. An unmasked input propagates to ``None``: there is nothing
        to carry.
    influence
        An ``(n_out, n_in)`` array whose non-zero entries say which input
        samples reach which output sample — a resampling weight matrix, a
        response matrix, or a boolean adjacency. ``None`` means the
        correspondence is one-to-one and the mask passes through unchanged.

    Returns
    -------
    numpy.ndarray or None
        A boolean mask over the output samples, or ``None`` if the input was
        unmasked.

    Examples
    --------
    >>> import numpy as np
    >>> propagate_mask(None) is None
    True
    >>> propagate_mask(np.array([False, True, False])).tolist()
    [False, True, False]

    Two input samples per output bin, and the second bin touches the masked
    sample:

    >>> weights = np.array([[0.5, 0.5, 0.0, 0.0], [0.0, 0.0, 0.5, 0.5]])
    >>> propagate_mask(np.array([False, False, True, False]), weights).tolist()
    [False, True]
    """
    mask = source.mask if isinstance(source, FunctionSamples) else source
    if mask is None:
        return None
    mask = np.asarray(mask, dtype=bool)
    if influence is None:
        return mask
    matrix = np.asarray(influence)
    if matrix.ndim != 2 or matrix.shape[1] != mask.size:
        raise TransformationError(
            f"propagate_mask needs an (n_out, n_in) influence matrix whose second axis matches "
            f"the {mask.size} input sample(s), got shape {matrix.shape}. Pass the same weights "
            f"the transformation applies to the values, or None if the mapping is one-to-one."
        )
    return np.any((matrix != 0) & mask.reshape(1, -1), axis=1)


# ---------------------------------------------------------------------------
# Requirements
# ---------------------------------------------------------------------------


def _dedupe(values: np.ndarray, rtol: float = COORDINATE_RTOL) -> np.ndarray:
    """Sort, and collapse coordinates that coincide to within *rtol*."""
    ordered = np.sort(np.asarray(values, dtype=float).reshape(-1))
    if ordered.size < 2:
        return ordered
    scale = np.maximum(np.abs(ordered[1:]), np.abs(ordered[:-1]))
    keep = np.diff(ordered) > rtol * np.where(scale > 0.0, scale, 1.0)
    return np.concatenate([ordered[:1], ordered[1:][keep]])


Density = tuple[float | None, float | None]


def _merge_runs(intervals: Sequence[tuple[float, float]]) -> list[tuple[float, float]]:
    """Union of closed intervals: sorted, with overlapping or touching pairs joined."""
    ordered = sorted(intervals)
    merged: list[list[float]] = [list(ordered[0])]
    for low, high in ordered[1:]:
        if low <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], high)
        else:
            merged.append([low, high])
    return [(low, high) for low, high in merged]


def _merge_intervals(
    intervals: Sequence[tuple[float, float]], densities: Sequence[Density]
) -> tuple[tuple[tuple[float, float], ...], tuple[Density, ...]]:
    """Union of intervals, **grouped by sampling density**.

    Intervals asking for the same density merge when they overlap; intervals
    asking for different densities never do, because merging them would apply
    the stricter of the two across the whole span. That is the difference
    between "cover 1-200 micron at R=40 *and* one 0.1 micron window at 0.0025
    micron sampling" (a few hundred coordinates) and the same statement with
    0.0025 micron sampling imposed from 1 to 200 micron (eighty thousand) —
    which would defeat the purpose ``DEVELOPMENT_PLAN.md`` §4.3 gives
    negotiation, namely to avoid evaluating expensively everywhere.
    """
    if not intervals:
        return (), ()
    grouped: dict[Density, list[tuple[float, float]]] = {}
    for pair, density in zip(intervals, densities, strict=True):
        grouped.setdefault(density, []).append(pair)
    merged = [(pair, density) for density, items in grouped.items() for pair in _merge_runs(items)]
    merged.sort(key=lambda item: item[0])
    return tuple(pair for pair, _ in merged), tuple(density for _, density in merged)


def _to_unit(raw: object, unit: u.UnitBase | None, what: str) -> np.ndarray:
    """*raw* as a plain float array in *unit*, converting once, here."""
    if isinstance(raw, u.Quantity):
        if unit is None:
            raise CompositionError(
                f"{what} was given in {raw.unit} but the requirement declares no unit. Either "
                f"pass unit={raw.unit!r} or give plain numbers already in the axis's own unit."
            )
        try:
            return np.asarray(raw.to_value(unit), dtype=float)
        except u.UnitConversionError as exc:
            raise CompositionError(
                f"{what} is in {raw.unit}, which cannot be converted to the requirement's "
                f"declared unit {unit}. Declare the requirement in a unit of the same physical "
                f"type; ampere does not apply spectral equivalencies here, because converting "
                f"a wavelength interval to frequency reverses it."
            ) from exc
    try:
        return np.asarray(raw, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TransformationError(
            f"{what} must be plain numbers, or one astropy Quantity covering all of them, got "
            f"{raw!r}. A sequence of per-element Quantities is not accepted: give the whole "
            f"sequence a unit, e.g. [(1.0, 2.0), (10.0, 20.0)] * u.um."
        ) from exc


def _first_quantity_unit(*candidates: object) -> u.UnitBase | None:
    for candidate in candidates:
        if isinstance(candidate, u.Quantity):
            return candidate.unit
    return None


@dataclasses.dataclass(frozen=True)
class AxisRequirement:
    """What one transformation needs of one coordinate axis of its channel.

    A requirement is a statement about *coordinates*, because that is what a
    container is built from (``results_schema.md`` §16). It has two
    independent halves, and either may be omitted:

    *where* — :attr:`intervals` (coverage) and :attr:`points` (coordinates that
    must be present exactly, e.g. the sample positions a response matrix was
    tabulated on);

    *how finely* — ``max_step`` (an absolute spacing) and
    ``min_resolving_power`` (:math:`\\lambda/\\Delta\\lambda`, the natural
    statement for a spectrograph). Giving both is allowed and means both must
    hold. A density applies to the intervals declared **alongside** it, which
    is why a requirement declaring one and no coverage is refused: after a
    :meth:`union` there would be no way to say where it applied.

    Requirements are declared **in the coordinates of the channel the
    instrument binds**, not in the coordinates of whichever step in the chain
    happens to want them; there is no backward propagation through a chain
    (see the spec's deliberate limitations).

    Parameters
    ----------
    axis
        Name of the container axis this constrains (``"spectral_axis"``,
        ``"time"``, ``"u"``, ...). Must match the kind's own axis name.
    unit
        Unit the numbers are in. May be adopted from a
        :class:`~astropy.units.Quantity` argument.
    intervals
        A ``(low, high)`` pair, a sequence of them, or a ``Quantity`` of shape
        ``(2,)`` or ``(n, 2)``. Overlapping intervals asking for the same
        density are merged.
    densities
        Per-interval ``(max_step, min_resolving_power)``, aligned with
        :attr:`intervals`. Not normally passed by hand: it is how
        :meth:`union` keeps one instrument's fine window from imposing its
        sampling on another's broad coverage. Omit it and every interval takes
        the requirement's own ``max_step``/``min_resolving_power``.
    points
        Coordinates that must appear exactly.
    max_step, min_resolving_power
        Sampling density; see above. Constructor-only (ruled 2026-09-02):
        both are folded into :attr:`densities` and do not survive as
        attributes, because after a :meth:`union` a single scalar could only
        be a strictest-anywhere summary, which misreads as applying to the
        whole requirement. Read :meth:`segments` instead.
    source
        Free text naming who asked, quoted back in composition errors.

    Examples
    --------
    >>> import astropy.units as u
    >>> band = AxisRequirement("spectral_axis", intervals=(1.0, 30.0) * u.um,
    ...                        min_resolving_power=50.0, source="photometer")
    >>> band.intervals, band.unit
    (((1.0, 30.0),), Unit("um"))
    >>> band.coordinates().size > 100
    True
    """

    axis: str
    _: dataclasses.KW_ONLY
    unit: u.UnitBase | None = None
    intervals: Any = ()
    densities: Any = None
    points: Any = None
    max_step: dataclasses.InitVar[float | None] = None
    min_resolving_power: dataclasses.InitVar[float | None] = None
    source: str = ""

    def __post_init__(self, max_step: float | None, min_resolving_power: float | None) -> None:
        set_ = object.__setattr__
        axis = _check_label(self.axis, "axis name")
        set_(self, "axis", axis)

        unit = _check_unit(self.unit, f"the unit of the requirement on axis {axis!r}")
        if unit is None:
            unit = _first_quantity_unit(self.intervals, self.points, max_step)
        set_(self, "unit", unit)

        what = f"the requirement on axis {axis!r}"
        if self.intervals is None:
            pairs = np.zeros((0, 2))
        else:
            given = _to_unit(self.intervals, unit, f"the intervals of {what}")
            pairs = np.zeros((0, 2)) if given.size == 0 else np.atleast_2d(given)
        if pairs.ndim != 2 or pairs.shape[1] != 2:
            raise TransformationError(
                f"{what} takes intervals as a (low, high) pair or a sequence of them, so the "
                f"array must have shape (2,) or (n, 2), got {np.shape(self.intervals)}."
            )
        for low, high in pairs:
            if not (np.isfinite(low) and np.isfinite(high)):
                raise TransformationError(f"{what} has a non-finite interval ({low}, {high}).")
            if low > high:
                raise TransformationError(
                    f"{what} has the interval ({low}, {high}) the wrong way round; intervals are "
                    f"(low, high)."
                )
            if low == high:
                raise TransformationError(
                    f"{what} has the degenerate interval ({low}, {high}). A single coordinate is "
                    f"a point, not an interval: pass points=[{low}] instead."
                )
        if self.points is None:
            set_(self, "points", None)
        else:
            points = _dedupe(_to_unit(self.points, unit, f"the points of {what}"))
            if not np.all(np.isfinite(points)):
                raise TransformationError(f"{what} has non-finite points.")
            points.setflags(write=False)
            set_(self, "points", points)

        step = max_step
        if isinstance(step, u.Quantity):
            step = float(_to_unit(step, unit, f"the max_step of {what}"))
        declared_step = _positive(step, f"the max_step of {what}")
        declared_power = _positive(min_resolving_power, f"the min_resolving_power of {what}")

        given_pairs = [(float(low), float(high)) for low, high in pairs]
        if self.densities is None:
            densities: list[Density] = [(declared_step, declared_power)] * len(given_pairs)
        else:
            densities = [
                (
                    _positive(one, f"a per-interval max_step of {what}"),
                    _positive(other, f"a per-interval min_resolving_power of {what}"),
                )
                for one, other in self.densities
            ]
            if len(densities) != len(given_pairs):
                raise TransformationError(
                    f"{what} was given {len(densities)} per-interval densities for "
                    f"{len(given_pairs)} interval(s); they are aligned one to one."
                )
        if (declared_step is not None or declared_power is not None) and not given_pairs:
            raise TransformationError(
                f"{what} declares a sampling density but no coverage for it to apply to. A "
                f"density belongs to the intervals declared alongside it — otherwise, once "
                f"requirements are unioned, there is no way to say where it holds. Add "
                f"intervals=(low, high)."
            )
        merged_pairs, merged_densities = _merge_intervals(given_pairs, densities)
        set_(self, "intervals", merged_pairs)
        set_(self, "densities", merged_densities)

        if not isinstance(self.source, str):
            raise TransformationError(f"the source of {what} must be a string, got {self.source!r}")

    @property
    def constrains_density(self) -> bool:
        """Whether this requirement says anything about how finely to sample."""
        return any(step is not None or power is not None for step, power in self.densities)

    def segments(self) -> tuple[tuple[float, float, float | None, float | None], ...]:
        """``(low, high, max_step, min_resolving_power)`` per interval.

        The requirement's canonical form: what :meth:`coordinates` builds
        from, and the only public statement of density — the constructor's
        scalars do not survive as attributes (see the class docstring).
        """
        return tuple(
            (low, high, step, power)
            for (low, high), (step, power) in zip(self.intervals, self.densities, strict=True)
        )

    def quantity(self, values: np.ndarray) -> Any:
        """*values* as a :class:`~astropy.units.Quantity` when the axis has a unit."""
        return values if self.unit is None else u.Quantity(values, self.unit)

    def convert_to(self, unit: u.UnitBase | None) -> AxisRequirement:
        """This requirement restated in *unit*.

        Refuses to guess: a unitless requirement and a unit-bearing one on the
        same axis are a composition error, not a pair to reconcile silently.
        """
        if unit == self.unit:
            return self
        if unit is None or self.unit is None:
            raise CompositionError(
                f"requirements on axis {self.axis!r} disagree about units: one is in "
                f"{self.unit if self.unit is not None else 'no unit'} and another in "
                f"{unit if unit is not None else 'no unit'}. Declare a unit on both (or neither) "
                f"— ampere will not assume a bare number is in the other's unit."
            )
        factor = float(u.Quantity(1.0, self.unit).to_value(unit))
        return AxisRequirement(
            self.axis,
            unit=unit,
            intervals=np.asarray(self.intervals, dtype=float).reshape(-1, 2) * factor
            if self.intervals
            else (),
            densities=tuple(
                (None if step is None else step * factor, power) for step, power in self.densities
            ),
            points=None if self.points is None else self.points * factor,
            source=self.source,
        )

    def union(self, other: AxisRequirement) -> AxisRequirement:
        """The weakest requirement that satisfies both.

        Coverage and required points accumulate. Density constraints stay
        attached to the intervals that asked for them, so unioning a broad,
        coarse requirement with a narrow, fine one gives a coarse grid with a
        fine window in it — not a fine grid everywhere. Merging them into a
        single strictest-everywhere constraint would defeat the purpose
        ``DEVELOPMENT_PLAN.md`` §4.3 gives negotiation.
        """
        if other.axis != self.axis:
            raise CompositionError(
                f"cannot union requirements on different axes, {self.axis!r} and {other.axis!r}."
            )
        aligned = other.convert_to(self.unit)
        points = [p for p in (self.points, aligned.points) if p is not None]
        sources = [s for s in (self.source, aligned.source) if s]
        return AxisRequirement(
            self.axis,
            unit=self.unit,
            intervals=np.asarray(self.intervals + aligned.intervals, dtype=float).reshape(-1, 2)
            if (self.intervals or aligned.intervals)
            else (),
            densities=self.densities + aligned.densities,
            points=np.concatenate(points) if points else None,
            source=", ".join(dict.fromkeys(sources)),
        )

    def coordinates(self) -> Any:
        """Coordinates satisfying this requirement, sorted and deduplicated.

        A :class:`~astropy.units.Quantity` when the axis has a unit, a bare
        array otherwise — in both cases exactly what a container constructor
        wants. Each interval is sampled at *its own* declared density and the
        results are merged, so every contributing requirement's constraint
        holds over its own coverage and nowhere else.

        This is a *reference* grid builder: a model is free to build its own
        grid, or to ignore the requirement entirely.
        """
        if not self.intervals and self.points is None:
            raise TransformationError(
                f"the requirement on axis {self.axis!r} declares no intervals and no points, so "
                f"there are no coordinates to build from it. Add intervals=(low, high), or ask "
                f"the model for its own grid."
            )
        pieces: list[np.ndarray] = []
        if self.points is not None:
            pieces.append(self.points)
        for low, high, step, power in self.segments():
            pieces.append(self._grid(low, high, step, power))
        return self.quantity(_dedupe(np.concatenate(pieces)))

    def _grid(self, low: float, high: float, step: float | None, power: float | None) -> np.ndarray:
        pieces = [np.array([low, high])]
        if step is not None:
            count = max(int(np.ceil((high - low) / step)) + 1, 2)
            pieces.append(np.linspace(low, high, count))
        if power is not None:
            if low <= 0.0:
                raise TransformationError(
                    f"the requirement on axis {self.axis!r} asks for resolving power {power} over "
                    f"an interval starting at {low}, but resolving power is a ratio and needs "
                    f"strictly positive coordinates. Use max_step for an interval that reaches "
                    f"zero."
                )
            ratio = 1.0 + 1.0 / power
            count = max(int(np.ceil(np.log(high / low) / np.log(ratio))) + 1, 2)
            pieces.append(np.geomspace(low, high, count))
        return np.concatenate(pieces)

    def __repr__(self) -> str:
        bits = [repr(self.axis)]
        if self.unit is not None:
            bits.append(str(self.unit))
        bits.extend(_segment_repr(segment) for segment in self.segments())
        if self.points is not None:
            bits.append(f"{self.points.size} point(s)")
        return f"<AxisRequirement {' '.join(bits)}>"


# The density scalars are constructor-only (ruled 2026-09-02). dataclasses
# leaves an InitVar's default behind as a class attribute, so an instance
# read of .max_step would silently return None instead of failing; remove
# them so the read raises AttributeError.
del AxisRequirement.max_step, AxisRequirement.min_resolving_power


def _segment_repr(segment: tuple[float, float, float | None, float | None]) -> str:
    low, high, step, power = segment
    density = [] if step is None else [f"step<={step:g}"]
    if power is not None:
        density.append(f"R>={power:g}")
    suffix = f"@{','.join(density)}" if density else ""
    return f"[{low:g},{high:g}]{suffix}"


@dataclasses.dataclass(frozen=True)
class ChannelRequirements:
    """Everything the instruments bound to one channel need of it.

    The output of :func:`negotiate`, and the input to :meth:`Model.compile_for`.

    Attributes
    ----------
    channel
        Channel name the requirements apply to.
    kind
        The container kind every instrument binding this channel expects; the
        most specific of them, checked to be compatible with all.
    axes
        Merged :class:`AxisRequirement` per axis name. May be empty: an
        instrument that publishes nothing still declares that it binds the
        channel, which is itself useful information.
    sources
        Labels of the instruments that contributed, for error messages and
        provenance.
    """

    channel: str
    kind: type[FunctionSamples]
    axes: Mapping[str, AxisRequirement]
    sources: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "axes", types.MappingProxyType(dict(self.axes)))

    def __contains__(self, axis: object) -> bool:
        return axis in self.axes

    def __getitem__(self, axis: str) -> AxisRequirement:
        try:
            return self.axes[axis]
        except KeyError:
            raise CompositionError(
                f"no requirement was published for axis {axis!r} of channel {self.channel!r}; "
                f"the axes with requirements are {sorted(self.axes)}. A model may sample an "
                f"unconstrained axis however it likes."
            ) from None

    def coordinates(self) -> dict[str, Any]:
        """Reference coordinates per constrained axis, ready for a container.

        The compile-once half of the compile-once/evaluate-many split: build
        the channel's container from these, then refill it with
        :meth:`~ampere.core.results_schema.FunctionSamples.with_values` on
        every evaluation.
        """
        return {name: requirement.coordinates() for name, requirement in self.axes.items()}

    def __repr__(self) -> str:
        return (
            f"<ChannelRequirements {self.channel!r} {self.kind.__name__} "
            f"axes={sorted(self.axes)} from {list(self.sources)}>"
        )


# ---------------------------------------------------------------------------
# Transformation
# ---------------------------------------------------------------------------


class Transformation(Parameterised, abc.ABC):
    """One step of an instrument: a container in, a container out.

    Subclassing is the whole extension mechanism (``DEVELOPMENT_PLAN.md`` §4.3
    makes that a first-class requirement): declare :attr:`ACCEPTS`, implement
    :meth:`apply`, and optionally register nuisance parameters and publish
    :meth:`requirements`. Nothing inside ampere changes.

    Class attributes
    ----------------
    ACCEPTS
        Container kinds this step can be given. Defaults to *any*.
    PRODUCES
        Container kind it returns, or ``None`` — the default — meaning it
        returns the kind it was given. Kind-preserving is much the commoner
        case (convolution, resampling, calibration); synthetic photometry
        (``Spectrum`` to ``PhotometricPoints``) is the other kind.

    Notes
    -----
    :meth:`__call__` wraps :meth:`apply` with the checks this contract owes its
    callers: the input kind, the output kind, and mask propagation. Call the
    instance, not ``apply``, unless you are deliberately bypassing them.
    """

    #: Container kinds this transformation accepts.
    ACCEPTS: ClassVar[tuple[type[FunctionSamples], ...]] = (FunctionSamples,)
    #: Kind produced, or ``None`` for "the same kind it was given".
    PRODUCES: ClassVar[type[FunctionSamples] | None] = None

    #: Capability flags (``DEVELOPMENT_PLAN.md`` §4.5), promoted into this ABC
    #: at the freeze (ruled 2026-09-03, ``inference.md`` §19.6). The defaults
    #: are the reference path's honest answers, so a step that stays silent
    #: promises nothing; Phase 2's torch/jax subclasses override them.
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Whether this step evaluates a batch of parameter vectors in one call.
    BATCHABLE: ClassVar[bool] = False
    #: Device this step's arrays live on: ``"cpu"``, ``"cuda"``, ``"cuda:0"``, ...
    DEVICE: ClassVar[str] = "cpu"
    #: Backend this step belongs to (the fourth flag, W2.12). ``"reference"``
    #: is the conservative default: a hand-written numpy step runs on the
    #: reference path, which is the base install's whole toolkit.
    BACKEND: ClassVar[str] = "reference"

    _label: str

    def __init__(self, *, label: str | None = None) -> None:
        if label is not None:
            self._label = _check_label(label, "transformation label")

    @property
    def label(self) -> str:
        """Component label for this step's parameters within a chain.

        Defaults to the class name in snake_case, so a chain's merged
        parameters read ``lsf.fwhm``, ``calibrate.scale``. Pass ``label=`` when
        one chain uses two steps of the same class.
        """
        try:
            return self._label
        except AttributeError:
            self._label = _default_label(type(self))
            return self._label

    # -- kinds ---------------------------------------------------------------

    @classmethod
    def accepts_kind(cls, kind: type[FunctionSamples]) -> bool:
        """Whether this step can be handed a container of *kind*."""
        return isinstance(kind, type) and issubclass(kind, cls.ACCEPTS)

    @classmethod
    def output_kind(cls, input_kind: type[FunctionSamples]) -> type[FunctionSamples]:
        """The kind this step returns when handed *input_kind*."""
        return input_kind if cls.PRODUCES is None else cls.PRODUCES

    # -- the one method a subclass must write --------------------------------

    @abc.abstractmethod
    def apply(self, samples: FunctionSamples, values: Mapping[str, Value] | None) -> Any:
        """Transform *samples*, given this step's own parameter *values*.

        Parameters
        ----------
        samples
            A container of a kind in :attr:`ACCEPTS`, already checked.
        values
            Mapping of this step's **local** parameter names to values, or
            ``None`` to use their declared values. Read it through
            ``self.context(values)`` so buffers and parameters arrive in one
            namespace and promoting a buffer stays a configuration change
            (``parameters.md`` §10).

        Returns
        -------
        FunctionSamples
            Of kind ``self.output_kind(type(samples))``. Propagate the input's
            mask — :func:`propagate_mask` is the rule — or
            :meth:`with_values`, which inherits it, will do it for you.
        """

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """What this step needs of the *channel's* coordinates. Empty by default.

        Optional at every level: a step that publishes nothing simply consumes
        whatever the model produces, which is the fixed-grid path.
        """
        return ()

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """Read the declarations of the steps after this one. Default: nothing.

        The chain-internal half of negotiation (ruled 2026-09-03,
        ``transformations.md`` §15 Q2 — gap I-3's mechanism, adopted in place
        of a ``pull_back``). ``requirements()`` are statements about the
        *channel's* coordinates, so a step downstream of a kind-changing step
        cannot publish at all; what it can do is tell the step before it what
        it needs — interferometric bandwidth and time smearing want extra
        (u, v) samples that are the Fourier step's own buffer, not the
        model's. :class:`Instrument` calls this once per step at chain
        construction, with the tuple of that step's successors, before the
        hot loop; a step that wants nothing inherits this no-op.

        The posture is **push-forward-and-raise**: a step that is handed
        something it cannot use fails loudly at evaluation rather than ampere
        inferring requirements backwards through the chain. The first real
        instances (the smearing steps) are Phase 4's.
        """

    # -- evaluation ----------------------------------------------------------

    def __call__(self, samples: FunctionSamples, values: Mapping[str, Value] | None = None) -> Any:
        """Check the kinds, run :meth:`apply`, check the result and its mask."""
        kind = type(samples)
        if not isinstance(samples, FunctionSamples):
            raise TransformationError(
                f"{type(self).__name__} was given {samples!r}, which is not a container. A "
                f"transformation maps one FunctionSamples to another; bind a ModelResult channel "
                f"with result.require(name, kind) first."
            )
        if not self.accepts_kind(kind):
            raise CompositionError(
                f"{type(self).__name__} accepts {self._accepted_names()}, but was given a "
                f"{kind.__name__}. Either bind this step to a channel of an accepted kind, or "
                f"put a transformation that produces one earlier in the chain."
            )
        result = self.apply(samples, values)
        expected = self.output_kind(kind)
        if not isinstance(result, expected):
            got = type(result).__name__
            raise TransformationError(
                f"{type(self).__name__}.apply returned a {got}, but the class declares "
                f"PRODUCES={self.PRODUCES.__name__ if self.PRODUCES else None}, so given a "
                f"{kind.__name__} it must return a {expected.__name__}. Fix apply, or declare "
                f"PRODUCES = {got}."
            )
        if samples.mask is not None and result.mask is None:
            raise TransformationError(
                f"{type(self).__name__}.apply dropped the mask: its input excluded "
                f"{int(np.count_nonzero(samples.mask))} of {samples.n_samples} sample(s) and its "
                f"output carries no mask. A masked sample carries zero information and that must "
                f"survive the chain (results_schema.md §16). Use propagate_mask(samples, weights) "
                f"and pass it as mask=..., or pass an explicit all-False mask if the output "
                f"genuinely does not depend on the excluded samples."
            )
        return result

    @classmethod
    def _accepted_names(cls) -> str:
        return " or ".join(kind.__name__ for kind in cls.ACCEPTS)

    def __repr__(self) -> str:
        free = self.parameters.free_size
        suffix = f", {free} free parameter(s)" if free else ""
        return f"<{type(self).__name__} {self.label!r}{suffix}>"


# ---------------------------------------------------------------------------
# Instrument
# ---------------------------------------------------------------------------


class Instrument:
    """An ordered chain of :class:`Transformation`\\ s bound to one channel.

    The chain is validated when it is built: each step must accept what the
    step before it produces, and the first must accept :attr:`input_kind`. The
    channel binding is checked when the instrument is evaluated, through
    :meth:`~ampere.core.results_schema.ModelResult.require`, so a wrong channel
    name or a wrong kind raises
    :class:`~ampere.core.exceptions.ChannelError` naming what *is* available.

    Parameters
    ----------
    steps
        The transformations, in application order. May be empty: an instrument
        with no steps is a pure channel binding, which is a legitimate thing to
        want when the model already produces the observable.
    channel
        Name of the :class:`~ampere.core.results_schema.ModelResult` channel to
        bind. Defaults to
        :data:`~ampere.core.results_schema.DEFAULT_CHANNEL`, so a single-output
        model needs no channel names at all.
    input_kind
        Kind expected on that channel. Inferred from the first step when it
        accepts exactly one kind; otherwise give it explicitly.
    label
        How a user identifies this instrument — in the requirements
        provenance, in failure reports, and (via the dataset it serves) in
        the joint parameter space. Defaults to the channel name, which is
        right when one instrument reads one channel; when more than one
        instrument reads a channel, distinct labels are **required**, checked
        at problem composition (ruled 2026-09-03, ``transformations.md`` §15
        Q4 — the label and the channel are different concepts, and their
        relationship is not one-to-one).
    meta
        Free-form metadata.

    Notes
    -----
    ``Instrument`` is deliberately **not** :class:`~ampere.core.parameter.
    Parameterised`: its parameters are its steps', merged. An instrument-level
    nuisance parameter is a one-step transformation, which is what a
    calibration scale factor already is.
    """

    def __init__(
        self,
        steps: Iterable[Transformation] = (),
        *,
        channel: str = DEFAULT_CHANNEL,
        input_kind: type[FunctionSamples] | None = None,
        label: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        self.steps: tuple[Transformation, ...] = tuple(steps)
        for step in self.steps:
            if not isinstance(step, Transformation):
                raise CompositionError(
                    f"an Instrument is a chain of Transformations, but one step is {step!r}. "
                    f"Subclass Transformation and implement apply()."
                )
        self.channel = _check_dotted_label(channel, "channel name")
        self.label = _check_dotted_label(channel if label is None else label, "instrument label")
        self.meta: Mapping[str, Any] = types.MappingProxyType(dict(meta) if meta else {})

        seen: dict[str, int] = {}
        for position, step in enumerate(self.steps):
            if step.label in seen:
                raise CompositionError(
                    f"instrument {self.label!r} has two steps labelled {step.label!r} (positions "
                    f"{seen[step.label]} and {position}). Step labels become component names when "
                    f"the chain's parameters are merged, so they must be unique: pass "
                    f"label='...' to one of them. They are not numbered automatically, because "
                    f"inserting a step would then silently rename every parameter after it."
                )
            seen[step.label] = position

        self.input_kind: type[FunctionSamples] = self._resolve_input_kind(input_kind)
        kind = self.input_kind
        for position, step in enumerate(self.steps):
            if not step.accepts_kind(kind):
                previous = (
                    f"step {position - 1} ({self.steps[position - 1].label!r}) produces"
                    if position
                    else f"the channel {self.channel!r} holds"
                )
                raise CompositionError(
                    f"instrument {self.label!r} cannot be composed: {previous} a "
                    f"{kind.__name__}, but step {position} ({step.label!r}, "
                    f"{type(step).__name__}) accepts {step._accepted_names()}. Insert a step that "
                    f"produces a {step._accepted_names()}, or reorder the chain."
                )
            kind = step.output_kind(kind)
        self.output_kind: type[FunctionSamples] = kind

        # Chain-internal negotiation (ruled 2026-09-03, §15 Q2 — gap I-3):
        # each step reads its successors' declarations once, here, before the
        # hot loop. The default configure_from is a no-op.
        #
        # Last step first (ruled 2026-09-07): a step's declarations may
        # themselves depend on what *its* successors asked for — a convolution
        # publishes nothing until it has learned the range downstream of it,
        # and then publishes that range padded. Configuring in forward order
        # handed each step successors that had not been configured yet, so a
        # convolution before another convolution read only the resampler's
        # range and never the second kernel's padding: chained same-axis
        # convolutions were under-padded. Walking the chain from the end
        # guarantees every step reads successors in their final state.
        for position in reversed(range(len(self.steps))):
            self.steps[position].configure_from(self.steps[position + 1 :])

        self._frozen_mapping: ParameterMapping | None = None
        self._frozen_declarations: tuple[tuple[str, tuple[int, ...]], ...] | None = None

    def _resolve_input_kind(self, declared: type[FunctionSamples] | None) -> type[FunctionSamples]:
        if declared is not None:
            if not (isinstance(declared, type) and issubclass(declared, FunctionSamples)):
                raise CompositionError(
                    f"input_kind must be a FunctionSamples subclass, got {declared!r}."
                )
            return declared
        if not self.steps:
            return FunctionSamples
        first = self.steps[0]
        if len(first.ACCEPTS) == 1:
            return first.ACCEPTS[0]
        raise CompositionError(
            f"instrument {self.label!r} cannot infer the kind it binds: its first step "
            f"({first.label!r}, {type(first).__name__}) accepts {first._accepted_names()}, so "
            f"there is no single answer. Pass input_kind=... to say which kind channel "
            f"{self.channel!r} holds."
        )

    # -- parameters ----------------------------------------------------------

    @property
    def mapping(self) -> ParameterMapping:
        """The chain's parameters merged, with each step's label as component.

        Recomputed on access rather than cached: a step is
        :class:`~ampere.core.parameter.Parameterised` and may still be
        reconfigured (``promote_buffer``) after the chain is built, and a stale
        snapshot of that is exactly the silent-drift bug class these contracts
        exist to end. The merge is a few dozen dictionary operations; the hot
        loop is inside :meth:`Transformation.apply`, not here.

        After :meth:`freeze` the snapshot is served instead — with the steps'
        declarations verified first, so a post-freeze reconfiguration is
        **refused**, never silently served stale.
        """
        if self._frozen_mapping is not None:
            if self._declarations() != self._frozen_declarations:
                raise CompositionError(
                    f"instrument {self.label!r} was frozen and one of its steps has been "
                    f"reconfigured since (promote_buffer, register_parameter). A frozen "
                    f"instrument's merged parameters are a snapshot — the joint parameter "
                    f"space built from them would silently disagree with the steps. Configure "
                    f"the steps first and freeze() afterwards, or build a fresh Instrument."
                )
            return self._frozen_mapping
        return ParameterSet.merge({step.label: step.parameters for step in self.steps})

    def _declarations(self) -> tuple[tuple[str, tuple[int, ...]], ...]:
        """A cheap fingerprint of each step's parameter declarations.

        ``Parameter`` is a frozen dataclass, so reconfiguration always
        replaces objects; object identity per step is therefore exactly the
        invariant :meth:`freeze` snapshots.
        """
        return tuple(
            (step.label, tuple(id(parameter) for parameter in step.parameters))
            for step in self.steps
        )

    def freeze(self) -> Instrument:
        """Snapshot the merged parameters; refuse later step reconfiguration.

        Ruled 2026-09-03 (``transformations.md`` §15 Q5, with ``inference.md``
        §19.8 as the consumer evidence): W1.7's ``Dataset`` calls this at
        construction, so the per-``log_prob`` re-merge in :meth:`__call__`
        disappears from the hot loop. Deliberately **not** a silent cache: a
        step reconfigured after the freeze makes the next :attr:`mapping`
        access raise rather than serve the stale snapshot — the silent-drift
        bug class this contract exists to end. Idempotent; returns ``self``.
        """
        self._frozen_mapping = ParameterSet.merge(
            {step.label: step.parameters for step in self.steps}
        )
        self._frozen_declarations = self._declarations()
        return self

    @property
    def parameters(self) -> ParameterSet:
        """The merged :class:`~ampere.core.parameter.ParameterSet` an engine sees.

        Names are qualified with the step label (``"calibrate.scale"``).
        """
        return self.mapping.merged

    # -- requirements --------------------------------------------------------

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """Every step's published requirements, in chain order.

        Subclasses may override to publish something the individual steps
        cannot know on their own.
        """
        published: list[AxisRequirement] = []
        for step in self.steps:
            for requirement in step.requirements():
                if not isinstance(requirement, AxisRequirement):
                    raise CompositionError(
                        f"step {step.label!r} of instrument {self.label!r} published "
                        f"{requirement!r}, which is not an AxisRequirement."
                    )
                published.append(
                    requirement
                    if requirement.source
                    # The InitVars are passed explicitly: replace() would
                    # otherwise fetch their (deleted) class-attribute
                    # defaults. densities already carries the declaration.
                    else dataclasses.replace(
                        requirement,
                        max_step=None,
                        min_resolving_power=None,
                        source=f"{self.label}.{step.label}",
                    )
                )
        return tuple(published)

    # -- evaluation ----------------------------------------------------------

    def bind(self, result: ModelResult) -> FunctionSamples:
        """The channel this instrument reads, kind-checked.

        Uses :meth:`~ampere.core.results_schema.ModelResult.require` rather
        than ``result[name]`` so a missing channel or a wrong kind raises
        :class:`~ampere.core.exceptions.ChannelError` (results_schema.md §16).
        """
        if not isinstance(result, ModelResult):
            raise TransformationError(
                f"an Instrument is evaluated on a ModelResult, got {result!r}. A bare container "
                f"becomes one with ModelResult(container), filed under "
                f"{DEFAULT_CHANNEL!r}."
            )
        return result.require(self.channel, self.input_kind)

    def __call__(
        self, result: ModelResult, values: Mapping[str, Value] | ArrayLike | None = None
    ) -> Any:
        """Predicted data for this instrument, given the model's output.

        Parameters
        ----------
        result
            What the model produced.
        values
            Values for the chain's *merged* parameters — a mapping keyed by
            qualified names (``"calibrate.scale"``), a flat free-parameter
            vector in :attr:`parameters` order, or ``None`` to use each step's
            declared values.
        """
        samples = self.bind(result)
        routed: dict[str, Mapping[str, Value] | None] = {step.label: None for step in self.steps}
        if values is not None:
            mapping = self.mapping
            routed.update(mapping.distribute(mapping.merged.complete(values)))
        for step in self.steps:
            samples = step(samples, routed[step.label])
        return samples

    def __repr__(self) -> str:
        chain = " -> ".join(step.label for step in self.steps) or "(no steps)"
        return (
            f"<Instrument {self.label!r} on channel {self.channel!r}: "
            f"{self.input_kind.__name__} -> {chain} -> {self.output_kind.__name__}>"
        )


# ---------------------------------------------------------------------------
# Negotiation
# ---------------------------------------------------------------------------


def negotiate(instruments: Iterable[Instrument]) -> dict[str, ChannelRequirements]:
    """Collect what every instrument needs, per channel.

    The union half of the negotiation protocol (``DEVELOPMENT_PLAN.md`` §4.3).
    Instruments binding the same channel must agree on its kind — one kind may
    be a subclass of the other, and the most specific wins — and their
    requirements on each axis are unioned: coverage and required points
    accumulate, density constraints take the stricter.

    Returns
    -------
    dict
        Channel name to :class:`ChannelRequirements`, in first-seen order.
        Every bound channel appears, even one with no requirements at all:
        "this channel is consumed" is itself worth telling a model.

    Raises
    ------
    CompositionError
        If two instruments bind one channel with incompatible kinds, or publish
        requirements on one axis in incompatible units.
    """
    channels: dict[str, ChannelRequirements] = {}
    for instrument in instruments:
        if not isinstance(instrument, Instrument):
            raise CompositionError(f"negotiate() takes Instruments, got {instrument!r}.")
        existing = channels.get(instrument.channel)
        kind = _reconcile_kinds(instrument, existing)
        axes: dict[str, AxisRequirement] = dict(existing.axes) if existing else {}
        for requirement in instrument.requirements():
            _check_axis_exists(requirement, kind, instrument)
            previous = axes.get(requirement.axis)
            axes[requirement.axis] = (
                requirement if previous is None else previous.union(requirement)
            )
        sources = (*(existing.sources if existing else ()), instrument.label)
        channels[instrument.channel] = ChannelRequirements(
            channel=instrument.channel, kind=kind, axes=axes, sources=sources
        )
    return channels


def _reconcile_kinds(
    instrument: Instrument, existing: ChannelRequirements | None
) -> type[FunctionSamples]:
    if existing is None:
        return instrument.input_kind
    theirs, ours = existing.kind, instrument.input_kind
    if issubclass(ours, theirs):
        return ours
    if issubclass(theirs, ours):
        return theirs
    raise CompositionError(
        f"instruments {list(existing.sources)} and {instrument.label!r} both bind channel "
        f"{instrument.channel!r} but expect unrelated kinds, {theirs.__name__} and "
        f"{ours.__name__}. One model channel holds one container; give them separate channels."
    )


def _check_axis_exists(
    requirement: AxisRequirement, kind: type[FunctionSamples], instrument: Instrument
) -> None:
    names = tuple(spec.name for spec in kind.AXES)
    if requirement.axis not in names:
        raise CompositionError(
            f"instrument {instrument.label!r} published a requirement on axis "
            f"{requirement.axis!r}, but the {kind.__name__} on channel {instrument.channel!r} is "
            f"indexed by {list(names)}. Requirements name the container's own axes."
        )


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------


class Model(Parameterised, abc.ABC):
    """Something that produces a :class:`~ampere.core.results_schema.ModelResult`.

    This is the smallest object the negotiation protocol needs on the model
    side, and nothing more: a :class:`~ampere.core.parameter.Parameterised`
    that returns named channels, plus the optional
    :meth:`compile_for` hook. The engine-facing surface —
    ``log_prob``, ``prior_transform``, ``simulate``, capability flags — belongs
    to the fitting problem, not to the model, and stays W1.7's (``DEVELOPMENT_
    PLAN.md`` §4.5).

    A subclass implements :meth:`evaluate` and may return a bare container,
    which is filed under
    :data:`~ampere.core.results_schema.DEFAULT_CHANNEL` — the single-output
    case stays trivial, as §4.2 requires.

    Examples
    --------
    >>> import numpy as np, astropy.units as u, scipy.stats as st
    >>> from ampere.core import Parameter, Spectrum
    >>> class Powerlaw(Model):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength, unit=u.micron)
    ...         self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
    ...     def evaluate(self, **values):
    ...         ctx = self.context(values)
    ...         return Spectrum(
    ...             ctx["wavelength"] * u.micron, ctx["wavelength"] ** ctx["index"] * u.Jy
    ...         )
    >>> model = Powerlaw(np.array([1.0, 2.0, 4.0]))
    >>> result = model(index=-2.0)
    >>> list(result), result.single().values.tolist()
    (['default'], [1.0, 0.25, 0.0625])
    >>> model([-1.0]).single().values.tolist()      # a flat free-parameter vector
    [1.0, 0.5, 0.25]
    """

    #: Capability flags (``DEVELOPMENT_PLAN.md`` §4.5), promoted into this ABC
    #: at the freeze (ruled 2026-09-03, ``inference.md`` §19.6). Conservative
    #: defaults — the reference path's honest answers; Phase 2's torch/jax
    #: model subclasses override them, and W1.7's ``declared_capabilities``
    #: reads them directly.
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Whether ``evaluate`` accepts a batch of parameter vectors in one call.
    BATCHABLE: ClassVar[bool] = False
    #: Device this model's arrays live on.
    DEVICE: ClassVar[str] = "cpu"
    #: Backend this model belongs to (the fourth flag, W2.12). ``"reference"``
    #: is the conservative default: a model written by hand in numpy runs on
    #: the reference path.
    BACKEND: ClassVar[str] = "reference"

    @abc.abstractmethod
    def evaluate(self, **values: Value) -> Any:
        """Compute this model's output for a complete set of parameter values.

        Every declared parameter arrives by name, fixed ones included, so model
        code never has to know which parameters this particular fit chose to
        vary. Return a :class:`~ampere.core.results_schema.ModelResult`, or a
        bare container for the single-channel case.
        """

    def __call__(
        self, values: Mapping[str, Value] | ArrayLike | None = None, /, **kwargs: Value
    ) -> ModelResult:
        """Evaluate, normalising *values* and checking what comes back.

        Accepts a mapping, a flat free-parameter vector, keyword values, or
        nothing at all. Anything not supplied falls back to the parameter's
        declared value, exactly as
        :meth:`~ampere.core.parameter.Parameterised.context` does; a free
        parameter with neither a supplied nor a declared value is an error.
        """
        if values is None or isinstance(values, Mapping):
            supplied = dict(values) if values is not None else {}
            supplied.update(kwargs)
            resolved_input = {p.name: p.value for p in self.parameters if p.value is not None}
            resolved_input.update(supplied)
            resolved = self.parameters.complete(resolved_input)
        elif kwargs:
            raise TransformationError(
                f"{type(self).__name__} was called with both a flat parameter vector and the "
                f"keyword value(s) {sorted(kwargs)}. A vector already fixes every free "
                f"parameter; pass one or the other."
            )
        else:
            resolved = self.parameters.complete(values)
        produced = self.evaluate(**resolved)
        # Attach the θ that produced the result (results_schema.md §17, ruled
        # 2026-09-01): a (θ, result) pair is then self-contained for emulator
        # training sets and provenance. A record evaluate() attached itself is
        # respected, never overwritten.
        if isinstance(produced, FunctionSamples):
            return ModelResult(produced, parameters=resolved)
        if not isinstance(produced, ModelResult):
            raise TransformationError(
                f"{type(self).__name__}.evaluate returned {produced!r}. A model returns a "
                f"ModelResult, or a single container which is filed under "
                f"{DEFAULT_CHANNEL!r}."
            )
        return produced if produced.parameters is not None else produced.with_parameters(resolved)

    def evaluate_batch(self, batch: Sequence[Mapping[str, Value]]) -> Sequence[Any]:
        """Compute this model's output for a **table** of parameter sets, in one call.

        The reference backend's reading of ``BATCHABLE`` (W3.1), and the same
        sentence as torch's and jax's: *this part can take a stack of θ*. The
        dialects differ because the batching does. On torch and jax a batchable
        model is one ``vmap`` can push a stack of θ through, and declaring the
        flag is all a backend model does. On the reference backend the case is
        an external code — a Fortran or C simulator, a grid interpolator, a
        vectorised analytic form — that is far cheaper called once with a table
        of *n* parameter sets than *n* times with one, and this is the method
        that offers it the table.

        Implement it **and** set ``BATCHABLE = True``;
        :meth:`~ampere.core.dataset.FittingProblem.simulate_many` uses it, under
        the serial executor, for a chunk at a time. It is not used under a
        process or thread pool, and the reason is not an oversight: a table
        evaluated in one call is by definition not partitioned across workers,
        so the two are alternative ways of spending the same batch. Nothing else
        in ampere calls it, so a model may implement it for
        ``simulate_many`` alone.

        Parameters
        ----------
        batch
            One complete set of parameter values per row, in draw order, every
            declared parameter present — exactly what :meth:`evaluate` receives,
            *n* times over.

        Returns
        -------
        Sequence
            One :class:`~ampere.core.results_schema.ModelResult` (or bare
            container) per row, **in the same order**. Order is the contract:
            row *i* out is row *i* in.

        Raises
        ------
        NotImplementedError
            By default. The base implementation exists so
            :meth:`~ampere.core.dataset.FittingProblem.simulate_many` can tell
            "declared and implemented" from "declared", never as a hook to call
            speculatively.

        Examples
        --------
        >>> import numpy as np, astropy.units as u, scipy.stats as st
        >>> from ampere.core import Parameter, Spectrum
        >>> class BatchedPowerlaw(Model):
        ...     BATCHABLE = True
        ...     def __init__(self, wavelength):
        ...         self.register_buffer("wavelength", wavelength, unit=u.micron)
        ...         self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
        ...     def evaluate(self, **values):
        ...         ctx = self.context(values)
        ...         return Spectrum(
        ...             ctx["wavelength"] * u.micron,
        ...             ctx["wavelength"] ** ctx["index"] * u.Jy,
        ...         )
        ...     def evaluate_batch(self, batch):
        ...         grid = self.context(batch[0])["wavelength"]
        ...         indices = np.asarray([row["index"] for row in batch])[:, None]
        ...         table = grid[None, :] ** indices          # one vectorised call
        ...         return [
        ...             Spectrum(grid * u.micron, row * u.Jy) for row in table
        ...         ]
        >>> model = BatchedPowerlaw(np.array([1.0, 2.0, 4.0]))
        >>> [result.single().values.tolist() for result in model.call_batch(
        ...     [{"index": -1.0}, {"index": -2.0}]
        ... )]
        [[1.0, 0.5, 0.25], [1.0, 0.25, 0.0625]]
        """
        raise NotImplementedError(
            f"{type(self).__name__} does not implement evaluate_batch(). A model that declares "
            f"BATCHABLE on the reference backend implements it to take a table of parameter "
            f"sets in one call; one that does not is simulated by the loop instead."
        )

    def call_batch(
        self, batch: Sequence[Mapping[str, Value] | ArrayLike | None]
    ) -> list[ModelResult]:
        """:meth:`__call__` for a table of θ: normalise in, check what comes back.

        The batched twin of :meth:`__call__`, and it does the same two jobs —
        completing each row against the declared parameters, and turning a bare
        container into a :class:`~ampere.core.results_schema.ModelResult` with
        its θ attached — so an implementer of :meth:`evaluate_batch` writes the
        interesting part and nothing else.
        """
        rows = [self.parameters.complete(_supplied(self, values)) for values in batch]
        produced = self.evaluate_batch(rows)
        results = list(produced)
        if len(results) != len(rows):
            raise TransformationError(
                f"{type(self).__name__}.evaluate_batch was given {len(rows)} parameter set(s) "
                f"and returned {len(results)} result(s). A batch call returns one result per "
                f"row, in row order — a shorter or longer sequence would silently misalign θ "
                f"with the outputs it produced."
            )
        normalised: list[ModelResult] = []
        for resolved, outcome in zip(rows, results, strict=True):
            if isinstance(outcome, FunctionSamples):
                normalised.append(ModelResult(outcome, parameters=resolved))
            elif isinstance(outcome, ModelResult):
                normalised.append(
                    outcome if outcome.parameters is not None else outcome.with_parameters(resolved)
                )
            else:
                raise TransformationError(
                    f"{type(self).__name__}.evaluate_batch returned {outcome!r} for one row. "
                    f"Each row returns a ModelResult, or a single container which is filed "
                    f"under {DEFAULT_CHANNEL!r}."
                )
        return normalised

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """One-off configuration from the instruments' published requirements.

        Called once, before the hot loop, with the output of :func:`negotiate`.
        The default **ignores it and returns self**, which is the contract:
        ``DEVELOPMENT_PLAN.md`` §4.3 says models are free to compute what they
        compute, and a model that does so must stay a working model.

        A model that honours the request builds its per-channel evaluation
        grids from ``requirements[channel].coordinates()``, keeps the resulting
        containers as templates, and refills them with
        :meth:`~ampere.core.results_schema.FunctionSamples.with_values` on
        every evaluation. Return ``self`` (having stored the templates) or a
        new, configured instance.

        **A model that engages with a requirement it cannot honour must
        raise** :class:`~ampere.core.exceptions.CompositionError` (ruled
        2026-09-03, ``transformations.md`` §15 Q3 — gap I-4's loud option).
        The asymmetry with the default is deliberate: ignoring negotiation
        entirely is a *declared* stance, visible in the model's code, and the
        model stays correct on its own grid; accepting a requirement and
        quietly under-sampling it is not — an under-sampled image grid
        aliases, and aliased visibilities look like real source structure
        rather than like an error. Raising is the default posture; the
        explicit opt-out that downgrades the refusal to a warning-and-proceed
        lives on the caller (W1.7's ``FittingProblem(lenient_compile=True)``),
        with the *unconfigured* model used instead.
        """
        return self


def _supplied(
    model: Model, values: Mapping[str, Value] | ArrayLike | None
) -> Mapping[str, Value] | ArrayLike:
    """One row of :meth:`Model.call_batch`'s table, normalised as ``__call__`` does it.

    A mapping (or nothing at all) falls back to the parameters' own declared
    values for anything it does not supply; a flat free-parameter vector already
    fixes every free parameter and is passed through untouched.
    """
    if values is None or isinstance(values, Mapping):
        resolved: dict[str, Value] = {
            parameter.name: parameter.value
            for parameter in model.parameters
            if parameter.value is not None
        }
        if values is not None:
            resolved.update(values)
        return resolved
    return values
