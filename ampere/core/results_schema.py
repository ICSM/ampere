"""The ModelResult schema contract (``DEVELOPMENT_PLAN.md`` §4.2).

This module is the backend-neutral vocabulary in which every ampere model says
*what it produced*. It exists to end the two-conventions bug class the plan's
§7 names first: legacy ampere reports model output both as a
``self.model.modelFlux`` attribute and as a ``result['spectrum']`` dictionary
entry, and code that guesses wrong fails confusingly and late. There is now one
answer — a :class:`ModelResult` — and it is a mapping of *named channels* to
typed containers.

Nothing here knows about torch, jax, numpyro or paramax: it is
numpy/scipy/astropy.units/stdlib only, per ``architecture.md`` §3-4.

The narrative version of this contract, with worked examples, is
``docs/design/contracts/results_schema.md``; that document's examples are
executed as doctests by ``tests/core/test_spec_doctests.py``, so it cannot
drift from this module without the suite going red.

What lives here
---------------
:class:`ModelResult`
    An immutable mapping of channel name to container. One model may produce
    several channels of the *same* kind — a low-resolution SED and a set of
    high-resolution CO windows are both ``Spectrum`` — which is why channels
    are named rather than typed-and-unique. Instruments bind by name and check
    the kind (§4.3); mismatches raise :class:`~ampere.core.exceptions.
    ChannelError`.
:class:`FunctionSamples`
    The base container: explicit coordinates + values (+ uncertainties, +
    mask). Every observable ampere handles is a *function* — flux(lambda),
    V(u,v), position(t) — and a container is a set of samples of it. There is
    no regular-grid assumption anywhere in this module.
:class:`Spectrum`, :class:`PhotometricPoints`, :class:`Image`, :class:`Cube`,
:class:`TimeSeries`, :class:`VisibilitySet`, :class:`ClosurePhases`
    The kinds the plan enumerates, and the one Phase 4 added beside them. Each
    fixes an axis signature, a layout, an ordering rule and whether its values
    may be complex.
:class:`Axis`
    One coordinate axis: name, values, unit, and the *advertised* structure
    (:attr:`Axis.regular`, :attr:`Axis.log_regular`) implementations may use to
    take a fast path. Advertised, never required.

Four conventions worth reading before the API
---------------------------------------------
**Coordinates are explicit and their ordering is checked, not assumed.**
``architecture.md`` §7 is explicit that legacy's implicit "wavelengths are
sorted and deduplicated" assumption is the thing to kill. Every container kind
therefore either *validates* its ordering rule at construction
(:class:`Spectrum`, :class:`TimeSeries`, :class:`Image`, :class:`Cube`) or
*documents* that it tolerates arbitrary order (:class:`PhotometricPoints`,
:class:`VisibilitySet`, :class:`ClosurePhases`). Neither is left to the
reader.

**Regularity is advertised, never required.** :attr:`Axis.regular` and
:attr:`Axis.log_regular` let an implementation choose FFT convolution or a
Toeplitz solve when it legitimately can. Nothing in this module, and nothing
downstream of it, may *require* an evenly spaced grid.

**Validation happens once; the hot loop uses** :meth:`FunctionSamples.
with_values`. Full construction validates coordinates in O(N), which is the
right price to pay once at composition time and the wrong price to pay per
likelihood evaluation. A model that evaluates on a fixed (or negotiated, §4.3)
grid builds its container once and calls :meth:`FunctionSamples.with_values`
thereafter, which reuses the already-validated axes and costs O(1) in them.

**Units are converted once, at construction.** A container stores plain
:class:`numpy.ndarray`\\ s plus an :class:`~astropy.units.UnitBase`; no
:class:`~astropy.units.Quantity` ever reaches the hot loop. This is the units
trap in ``DEVELOPMENT_PLAN.md`` §7, applied here.

Masks, and what they are not
----------------------------
:attr:`FunctionSamples.mask` follows the numpy/astropy convention — ``True``
marks an *excluded* sample — and means **zero information**: the sample is to
be treated exactly as if it had not been observed. ``prior_art.md`` lesson R1
records that RHMF expresses the same idea as zero weight or infinite
uncertainty, and :meth:`FunctionSamples.weights` and
:meth:`FunctionSamples.masked_uncertainty` produce both of those
representations, so a consumer may use whichever suits its algebra. Masked
values stay in place in :attr:`FunctionSamples.values`; ampere never threads
sentinel NaNs through arithmetic.

Masking is **not** censoring. An upper limit carries information and belongs to
the likelihood contract (W1.6, issue #11), not here; see
``docs/design/contracts/results_schema.md`` §8 for the hook this module leaves
open for it.

Examples
--------
A model with one output returns a bare container and is filed under
:data:`DEFAULT_CHANNEL`, so the simple case costs nothing:

>>> import numpy as np
>>> import astropy.units as u
>>> sed = Spectrum([1.0, 2.0, 4.0] * u.um, [3.0, 2.5, 1.0] * u.Jy)
>>> ModelResult(sed).single() is sed
True

A model with several outputs names them, and an instrument binds by name with a
kind check. Two channels of the *same* kind is the case that makes naming
necessary:

>>> co = Spectrum([866.9, 866.95, 867.0] * u.um, [1.0, 1.4, 1.1] * u.Jy)
>>> result = ModelResult({"sed_lowres": sed, "co_windows": co})
>>> result.require("co_windows", Spectrum).n_samples
3
>>> result.require("sed_lowres", VisibilitySet)
Traceback (most recent call last):
    ...
ampere.core.exceptions.ChannelError: channel 'sed_lowres' holds a Spectrum,
but a VisibilitySet was required. ...
"""

from __future__ import annotations

import dataclasses
import enum
import types
from collections.abc import Callable, Iterator, Mapping
from typing import Any, ClassVar

import astropy.units as u
import numpy as np

from .exceptions import ChannelError, SchemaError

__all__ = [
    "COORDINATE_RTOL",
    "DEFAULT_CHANNEL",
    "REGULARITY_RTOL",
    "AnomalyScore",
    "Axis",
    "AxisSpec",
    "ClosurePhases",
    "Cube",
    "FunctionSamples",
    "Image",
    "Layout",
    "ModelResult",
    "Order",
    "PhotometricPoints",
    "Spectrum",
    "TimeSeries",
    "VisibilitySet",
]

#: Channel name given to a model that returns a single bare container. The
#: plan's §4.2 requires the one-channel case to stay trivial; this is how.
#:
#: **The name is deliberately not reserved** (ruled 2026-09-01): a user may
#: name a real channel ``"default"``, and it is then indistinguishable from
#: the automatic one. This is an accepted, documented clash — reserving the
#: obvious word was judged more annoying than the ambiguity it prevents. If
#: an instrument must bind by name, name the channel something distinctive.
DEFAULT_CHANNEL = "default"

#: Relative tolerance for *advertising* an axis as evenly spaced
#: (:attr:`Axis.regular`). Deliberately tight: a false positive would send an
#: implementation down an FFT path the data do not justify.
REGULARITY_RTOL = 1e-9

#: Relative tolerance used when the union of several requirements produces
#: coordinates that coincide (routinely, at the edges of adjacent intervals).
#: Two coordinates closer than this fraction of their magnitude are one
#: coordinate — a container kind that requires strictly increasing coordinates
#: would otherwise reject the negotiated grid, and a covariance built on
#: near-duplicate coordinates is near-singular.
#:
#: It lives here, rather than in ``transform.py`` where negotiation uses it,
#: because :meth:`Axis.locate` is the exact inverse of that collapsing and the
#: two must not be able to drift apart (W2.1; ``spectrum_photometry.md``
#: Gap 1). ``transform.py`` re-exports it, so ``COORDINATE_RTOL`` keeps both of
#: its established import paths.
COORDINATE_RTOL = 1e-12

ArrayLike = Any


# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------


def _check_channel_name(name: object) -> str:
    """Validate a channel name.

    A channel name is a Python identifier or a ``.``-separated sequence of
    them (ruled 2026-09-02, closing `results_schema.md` §17's nested-channel
    question): qualification is how a model namespaces grouped output —
    per-object channels (``obj1.sed``) or grouped data such as one object's
    several CO lines (``co.j3_2``) — mirroring how ``ParameterSet.merge``
    qualifies parameter names. Nothing interprets the dots; ``require`` still
    matches the full name. Channel names are not passed as keyword arguments
    (that is the parameter contract's reason), but they do become group names
    on serialisation and coordinate labels in ``ampere.results``, and
    dot-qualified identifiers are the portable intersection of what those
    tolerate.
    """
    if not isinstance(name, str) or not name:
        raise SchemaError(f"a channel name must be a non-empty string, got {name!r}")
    if not all(piece.isidentifier() for piece in name.split(".")):
        raise SchemaError(
            f"channel name {name!r} is not usable: a channel name must be a valid Python "
            f"identifier or a '.'-separated sequence of them, because channel names become "
            f"group names on serialisation and coordinate labels in ampere.results. Rename "
            f"the channel, e.g. {_suggest_identifier(name)!r}."
        )
    return name


def _suggest_identifier(name: str) -> str:
    """Best-effort repair of a bad channel name, for the error message."""
    cleaned = "".join(char if char.isalnum() else "_" for char in name)
    if not cleaned or not cleaned[0].isalpha():
        cleaned = f"channel_{cleaned.lstrip('_')}"
    return cleaned.strip("_") or "channel"


# ---------------------------------------------------------------------------
# Arrays and units
# ---------------------------------------------------------------------------


def _readonly(array: np.ndarray) -> np.ndarray:
    """Return ``array`` with writing disabled, so containers stay immutable."""
    view = array.view()
    view.flags.writeable = False
    return view


def _as_array(
    raw: ArrayLike,
    what: str,
    *,
    unit: u.UnitBase | None = None,
    allow_complex: bool = False,
) -> tuple[np.ndarray, u.UnitBase | None]:
    """Convert ``raw`` to a read-only array plus a unit, exactly once.

    ``raw`` may be a :class:`~astropy.units.Quantity` (its unit is adopted, or
    converted to ``unit`` if one is given) or anything :func:`numpy.asarray`
    accepts (``unit`` is attached as metadata). The result is a plain array:
    no ``Quantity`` survives into the hot loop (``DEVELOPMENT_PLAN.md`` §7).
    """
    if isinstance(raw, u.Quantity):
        try:
            array = np.asarray(raw.to_value(unit) if unit is not None else raw.value)
        except u.UnitConversionError as exc:
            raise SchemaError(
                f"{what} was given in {raw.unit} but this container declares {unit}, and the "
                f"two are not convertible. Pass the values in a convertible unit, or drop the "
                f"explicit unit= and let the Quantity's own unit be adopted. ({exc})"
            ) from exc
        resolved = unit if unit is not None else raw.unit
    else:
        array = np.asarray(raw)
        resolved = unit

    if array.dtype.kind not in ("f", "c", "i", "u", "b"):
        raise SchemaError(
            f"{what} must be numeric, got dtype {array.dtype}. Coordinates, values and "
            f"uncertainties are numeric arrays; non-numeric labels belong in extra_coords."
        )
    if array.dtype.kind == "c" and not allow_complex:
        raise SchemaError(
            f"{what} include complex numbers, but this container kind holds real values. "
            f"Complex data "
            f"belong in a VisibilitySet (or a container kind that declares "
            f"ALLOW_COMPLEX = True); to keep a real projection of complex data, take its "
            f"amplitude, phase, real or imaginary part explicitly."
        )
    if array.dtype.kind in ("i", "u", "b"):
        array = array.astype(float)
    return _readonly(array), resolved


def _check_unit(unit: object, what: str) -> u.UnitBase | None:
    """Coerce ``unit`` to a unit object, or ``None``."""
    if unit is None:
        return None
    try:
        return u.Unit(unit)
    except (ValueError, TypeError) as exc:
        raise SchemaError(f"{what} is not a usable astropy unit: {unit!r} ({exc})") from exc


# ---------------------------------------------------------------------------
# Sampling geometry
# ---------------------------------------------------------------------------


class Layout(enum.Enum):
    """How a container's axes index its values.

    Neither layout implies an evenly spaced grid; both admit arbitrary
    coordinate values.
    """

    #: Every axis has length N and ``values.shape == (N,)``: an arbitrary point
    #: set in as many dimensions as there are axes. Scattered (u,v) coverage,
    #: a spectrum, a light curve.
    POINTS = "points"

    #: Axes are separable: axis *i* has length ``n_i`` and ``values.shape ==
    #: (n_0, n_1, ...)``. The axis coordinates themselves may be spaced however
    #: they like — separable is not the same claim as regular.
    GRID = "grid"


class Order(enum.Enum):
    """The ordering rule a coordinate axis must satisfy at construction."""

    #: Strictly increasing. Required where downstream algebra depends on it —
    #: quasiseparable GP solvers need ordered 1D coordinates, and duplicate
    #: coordinates make a covariance matrix singular.
    STRICTLY_INCREASING = "strictly_increasing"

    #: Strictly increasing *or* strictly decreasing. Sky axes legitimately run
    #: either way (RA decreases left-to-right in most conventions).
    STRICTLY_MONOTONIC = "strictly_monotonic"

    #: No ordering requirement. Used where arbitrary order is genuinely
    #: meaningful — a (u,v) point set, a bag of photometric points — and
    #: documented as tolerated rather than left unstated.
    ANY = "any"


@dataclasses.dataclass(frozen=True)
class AxisSpec:
    """The declaration of one coordinate axis of a container kind.

    Attributes
    ----------
    name
        Keyword the container's constructor accepts, and the axis's name.
    physical_types
        Accepted :func:`astropy.units.get_physical_type` results. Empty means
        "any physical type", but see ``equivalent_units``.
    equivalent_units
        Accepted alternative units, tested with
        :meth:`~astropy.units.UnitBase.is_equivalent`. An axis passes if it
        matches ``physical_types`` *or* ``equivalent_units``; this exists
        because some legitimate coordinates (spatial frequency in rad⁻¹) have
        no named physical type in astropy.
    order
        The ordering rule, validated at construction (:class:`Order`).
    """

    name: str
    physical_types: tuple[str, ...] = ()
    equivalent_units: tuple[u.UnitBase, ...] = ()
    order: Order = Order.ANY

    @property
    def unit_required(self) -> bool:
        """Whether an axis of this kind must carry a unit."""
        return bool(self.physical_types or self.equivalent_units)


@dataclasses.dataclass(frozen=True, eq=False)
class Axis:
    """One coordinate axis: values, unit, and advertised structure.

    An ``Axis`` is built by a container's constructor, not usually by hand. It
    is immutable and its ``values`` array is read-only.

    Attributes
    ----------
    name, values, unit
        As given. ``values`` is a plain read-only array in ``unit``.
    regular
        ``True`` if the coordinates are evenly spaced to
        :data:`REGULARITY_RTOL`. **Advertised, never required**
        (``architecture.md`` §7): an implementation may use it to choose an FFT
        convolution or a Toeplitz solve, and must work without it.
    log_regular
        ``True`` if the *logarithms* of the coordinates are evenly spaced —
        i.e. constant resolving power lambda/dlambda, the case in which a
        constant-velocity LSF becomes a convolution. Requires strictly positive
        coordinates.
    step
        The common spacing when :attr:`regular`, else ``None``.
    """

    name: str
    values: np.ndarray
    unit: u.UnitBase | None
    regular: bool
    log_regular: bool
    step: float | None

    @classmethod
    def build(cls, name: str, values: np.ndarray, unit: u.UnitBase | None) -> Axis:
        """Compute the advertised structure and return the axis."""
        regular, step = _spacing(values)
        if values.size > 1 and bool(np.all(values > 0.0)):
            log_regular, _ = _spacing(np.log(values))
        else:
            log_regular = False
        return cls(
            name=name, values=values, unit=unit, regular=regular, log_regular=log_regular, step=step
        )

    @property
    def size(self) -> int:
        """Number of coordinate values on this axis."""
        return int(self.values.size)

    def quantity(self) -> u.Quantity:
        """The coordinates as a :class:`~astropy.units.Quantity`.

        A convenience for plotting and reporting — never for the hot loop.
        """
        return u.Quantity(
            self.values, self.unit if self.unit is not None else u.dimensionless_unscaled
        )

    def locate(self, values: ArrayLike) -> np.ndarray:
        """Indices of *values* on this axis, matched within :data:`COORDINATE_RTOL`.

        The lookup a step needs when its buffer is tabulated on coordinates it
        published as ``points=`` (``spectrum_photometry.md`` Gap 1, ruled
        2026-09-03, landed W2.1). As soon as a second instrument binds the same
        channel, negotiation's union hands the step a **larger, possibly
        reordered** grid; a step that reads ``samples.values`` positionally
        against its own tabulation is then silently wrong, and fails as a bare
        ``matmul`` shape error naming neither channel nor negotiation. The
        supported pattern is::

            index = samples.axis("spectral_axis").locate(self.tabulated)
            response @ samples.values[index]

        Matching is within :data:`COORDINATE_RTOL` rather than exact, and that
        is the point rather than a convenience. The union collapses coordinates
        that coincide to within exactly that tolerance and keeps one
        representative, so a step's own published coordinate may legitimately
        differ from the one that survived. The predicate here is that
        collapsing read backwards, sharing the constant so the two cannot
        drift.

        The axis is not assumed sorted: ``PhotometricPoints`` and
        ``VisibilitySet`` declare :attr:`Order.ANY`, and Gap 1's own scenario is
        a photometry step. Where two coordinates are both within tolerance the
        nearest wins.

        Parameters
        ----------
        values
            Coordinates to look up, in this axis's own unit. Any shape; the
            result is flat, in the order given.

        Returns
        -------
        numpy.ndarray
            Integer (``intp``) indices into :attr:`values`, one per requested
            coordinate.

        Raises
        ------
        SchemaError
            If any coordinate has no match within the tolerance. That indicates
            a negotiation defect rather than a usage error — the grid was built
            to satisfy the very requirement the step published — so the message
            names the unmatched values instead of returning a sentinel.
        """
        targets = np.asarray(values, dtype=float).reshape(-1)
        coordinates = np.asarray(self.values, dtype=float).reshape(-1)
        if targets.size == 0:
            return np.empty(0, dtype=np.intp)
        if coordinates.size == 0:
            raise SchemaError(
                f"axis {self.name!r} is empty, so none of the {targets.size} requested "
                f"coordinate(s) can be located on it."
            )
        # Nearest neighbour via a sorted view, so an Order.ANY axis is handled
        # without an O(n*m) scan.
        order = np.argsort(coordinates, kind="stable")
        ordered = coordinates[order]
        upper = np.searchsorted(ordered, targets)
        right = np.clip(upper, 0, ordered.size - 1)
        left = np.clip(upper - 1, 0, ordered.size - 1)
        nearest = np.where(
            np.abs(ordered[left] - targets) <= np.abs(ordered[right] - targets), left, right
        )
        index = order[nearest]
        matched = coordinates[index]
        # The same predicate the union's collapsing uses: two coordinates are
        # one when they differ by no more than COORDINATE_RTOL of the larger
        # magnitude, with the scale floored at 1 so a pair at zero still works.
        scale = np.maximum(np.abs(matched), np.abs(targets))
        tolerance = COORDINATE_RTOL * np.where(scale > 0.0, scale, 1.0)
        missed = ~(np.abs(matched - targets) <= tolerance)
        if bool(np.any(missed)):
            unmatched = targets[missed]
            shown = ", ".join(repr(float(value)) for value in unmatched[:5])
            more = "" if unmatched.size <= 5 else f" (and {unmatched.size - 5} more)"
            raise SchemaError(
                f"axis {self.name!r} has no coordinate within COORDINATE_RTOL "
                f"({COORDINATE_RTOL:g}) of {shown}{more}. The axis spans "
                f"[{coordinates.min():g}, {coordinates.max():g}] over {coordinates.size} "
                f"point(s). A step should only ask for coordinates it published as a "
                f"points= requirement, so this is a negotiation defect rather than a "
                f"usage error."
            )
        return index.astype(np.intp, copy=False)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Axis):
            return NotImplemented
        return (
            self.name == other.name
            and self.unit == other.unit
            and np.array_equal(self.values, other.values)
        )

    __hash__: ClassVar[None] = None

    def __repr__(self) -> str:
        structure = " regular" if self.regular else (" log-regular" if self.log_regular else "")
        return f"<Axis {self.name!r} n={self.size}{_unit_suffix(self.unit)}{structure}>"


def _unit_suffix(unit: u.UnitBase | None) -> str:
    """`` " um"`` for a real unit, ``""`` for none or dimensionless."""
    text = "" if unit is None else str(unit)
    return f" {text}" if text else ""


def format_axis_label(name: str, unit: u.UnitBase | None) -> str:
    """``"baseline length [m]"``, or just *name* where there is no unit.

    **W5.3.** The exact convention :func:`ampere.results._plotting.axis_label`
    uses for a stored axis; a kind's :attr:`FunctionSamples.PLOT_COORDINATE`
    default lives in ``ampere.core`` and cannot import that plotting helper
    (``architecture.md`` §4 rule 1), so the one-line format is duplicated here
    rather than shared, and the two must be kept in step by inspection.
    """
    text = "" if unit is None else str(unit)
    return f"{name} [{text}]" if text else name


def _spacing(values: np.ndarray) -> tuple[bool, float | None]:
    """Return ``(evenly_spaced, step)`` for a 1D coordinate array."""
    if values.size < 2:
        return True, None
    diffs = np.diff(values)
    first = float(diffs.flat[0])
    if first == 0.0 or not np.all(np.isfinite(diffs)):
        return False, None
    if bool(np.allclose(diffs, first, rtol=REGULARITY_RTOL, atol=0.0)):
        return True, first
    return False, None


def _check_order(axis_values: np.ndarray, spec: AxisSpec, kind: str) -> None:
    """Enforce ``spec.order`` on ``axis_values``, loudly."""
    if spec.order is Order.ANY or axis_values.size < 2:
        return
    diffs = np.diff(axis_values)
    if spec.order is Order.STRICTLY_INCREASING:
        if not bool(np.all(diffs > 0.0)):
            raise SchemaError(_ordering_message(axis_values, spec, kind, "strictly increasing"))
    elif not (bool(np.all(diffs > 0.0)) or bool(np.all(diffs < 0.0))):
        raise SchemaError(_ordering_message(axis_values, spec, kind, "strictly monotonic"))


def _ordering_message(values: np.ndarray, spec: AxisSpec, kind: str, requirement: str) -> str:
    """Explain an ordering failure, and say what to do about it."""
    diffs = np.diff(values)
    duplicates = int(np.count_nonzero(diffs == 0.0))
    detail = (
        f"it contains {duplicates} repeated coordinate(s)"
        if duplicates
        else "its coordinates are not sorted"
    )
    fix = (
        f"Use {kind}.from_unsorted(...) to sort them, "
        "or split genuinely overlapping data into separate channels."
    )
    return (
        f"{kind} requires its {spec.name!r} coordinates to be {requirement}, but {detail}. "
        f"This is checked here rather than assumed, because quasiseparable GP solvers and "
        f"resampling both silently misbehave on unordered or duplicated coordinates. {fix}"
    )


# ---------------------------------------------------------------------------
# The container base
# ---------------------------------------------------------------------------

#: The callable form of :attr:`FunctionSamples.PLOT_COORDINATE`: the kind's
#: own axes, by name, to one ordered coordinate and its axis label.
_PlotCoordinateFn = Callable[[Mapping[str, "Axis"]], tuple[np.ndarray, str]]


@dataclasses.dataclass(frozen=True, eq=False, init=False)
class FunctionSamples:
    """Samples of an underlying function: coordinates, values, uncertainty, mask.

    This is the base of every container kind. It is not abstract in the
    ``abc`` sense — it is perfectly usable as a base for a user-defined kind,
    which is the point (``DEVELOPMENT_PLAN.md`` §4.3 makes user extensibility a
    first-class requirement). A subclass declares :attr:`AXES`,
    :attr:`LAYOUT` and :attr:`ALLOW_COMPLEX`, and usually a friendlier
    ``__init__``; everything else is inherited.

    Attributes
    ----------
    axes
        Ordered tuple of :class:`Axis`, one per :attr:`AXES` entry.
    values
        Read-only array. Shape is ``(n_samples,)`` under
        :attr:`Layout.POINTS` and the tuple of axis lengths under
        :attr:`Layout.GRID`.
    unit
        Unit of :attr:`values`, or ``None`` for dimensionless output.
    uncertainty
        Standard deviation, same shape as :attr:`values`, or ``None``. Always
        real and non-negative, including for complex values — see
        :class:`VisibilitySet`.
    mask
        Boolean array, same shape as :attr:`values`, or ``None``. ``True``
        marks an **excluded** sample (numpy/astropy convention); see the module
        docstring.
    extra_coords
        Per-sample auxiliary arrays that do *not* define the sampling geometry
        — filter names, per-visibility frequency, epoch labels. Aligned with
        :attr:`values`; may be non-numeric.
    fidelity
        Optional tag naming which variant of a model produced this channel
        (``"cheap"``, ``"lte"``, ``"full_nlte"``). The multi-fidelity hook the
        plan's design horizon (a) reserves; ampere attaches no meaning to the
        string itself.
    meta
        Free-form immutable metadata.
    """

    #: Axis signature of this container kind, in order.
    AXES: ClassVar[tuple[AxisSpec, ...]] = ()
    #: How the axes index the values.
    LAYOUT: ClassVar[Layout] = Layout.POINTS
    #: Whether values may be complex.
    ALLOW_COMPLEX: ClassVar[bool] = False
    #: The default coordinate a diagnostic plot draws this kind's values
    #: against, for a kind with **several** axes (``results.md`` §4/§13 item
    #: 14, *Amended W5.3*). A single-axis kind needs nothing here — the one
    #: axis is already the plotted coordinate, exactly as before this
    #: attribute existed, which is what keeps every such plot byte-identical.
    #: Either an axis name (one of :attr:`AXES`), or a callable taking a
    #: ``Mapping[str, Axis]`` (the kind's own axes, by name) and returning
    #: ``(coordinates, label)`` — one ordered 1-D array and the string a plot
    #: puts on its axis (:func:`format_axis_label` is the label convention).
    #: ``None`` (the default) means the kind has no default: a caller must
    #: pass ``coordinate=`` to any of the four diagnostics that need one, or
    #: they refuse by name, listing :attr:`AXES`. Resolved once, by
    #: ``ampere.results.derived``, from the live container at the point a
    #: derived group is built — plotting-side infrastructure, not part of the
    #: frozen §4 contract's shape (D1 of Phase 4: a kind is class attributes).
    PLOT_COORDINATE: ClassVar[str | _PlotCoordinateFn | None] = None

    axes: tuple[Axis, ...]
    values: np.ndarray
    unit: u.UnitBase | None
    uncertainty: np.ndarray | None
    mask: np.ndarray | None
    extra_coords: Mapping[str, np.ndarray]
    fidelity: str | None
    meta: Mapping[str, Any]

    def __init__(
        self,
        coordinates: Mapping[str, ArrayLike],
        values: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        kind = type(self).__name__
        set_ = object.__setattr__

        axes = self._build_axes(coordinates, kind)
        set_(self, "axes", axes)

        declared_unit = _check_unit(unit, f"{kind}'s value unit")
        value_array, value_unit = _as_array(
            values, f"{kind} values", unit=declared_unit, allow_complex=self.ALLOW_COMPLEX
        )
        expected = self._expected_shape(axes)
        if value_array.shape != expected:
            raise SchemaError(
                self._shape_message(kind, "values", value_array.shape, expected, axes)
            )
        set_(self, "values", value_array)
        set_(self, "unit", value_unit)

        set_(self, "uncertainty", self._build_uncertainty(uncertainty, value_unit, expected, kind))
        set_(self, "mask", self._build_mask(mask, expected, kind))
        set_(self, "extra_coords", self._build_extra_coords(extra_coords, expected, kind))

        if fidelity is not None and (not isinstance(fidelity, str) or not fidelity):
            raise SchemaError(
                f"{kind}'s fidelity tag must be a non-empty string naming the model variant "
                f"that produced this channel (e.g. 'cheap', 'full_nlte'), got {fidelity!r}."
            )
        set_(self, "fidelity", fidelity)
        set_(self, "meta", types.MappingProxyType(dict(meta) if meta else {}))

    # -- construction helpers ------------------------------------------------

    @classmethod
    def _build_axes(cls, coordinates: Mapping[str, ArrayLike], kind: str) -> tuple[Axis, ...]:
        expected_names = tuple(spec.name for spec in cls.AXES)
        given = tuple(coordinates)
        if set(given) != set(expected_names):
            missing = [name for name in expected_names if name not in coordinates]
            extra = [name for name in given if name not in expected_names]
            problem = []
            if missing:
                problem.append(f"missing {missing}")
            if extra:
                problem.append(f"unexpected {extra}")
            raise SchemaError(
                f"{kind} is defined by the coordinate axes {list(expected_names)}, but the "
                f"coordinates given were {list(given)} ({'; '.join(problem)}). Per-sample data "
                f"that do not define the sampling geometry belong in extra_coords."
            )

        axes: list[Axis] = []
        for spec in cls.AXES:
            raw = coordinates[spec.name]
            array, unit = _as_array(raw, f"{kind}'s {spec.name!r} coordinates")
            if array.ndim != 1:
                raise SchemaError(
                    f"{kind}'s {spec.name!r} coordinates must be one-dimensional, got shape "
                    f"{array.shape}. Each axis carries its own coordinate values; a "
                    f"{cls.LAYOUT.value} layout combines them, so no axis is ever a mesh."
                )
            resolved = cls._check_axis_unit(spec, unit, kind)
            _check_order(array, spec, kind)
            axes.append(Axis.build(spec.name, array, resolved))

        if cls.LAYOUT is Layout.POINTS and len({axis.size for axis in axes}) > 1:
            sizes = {axis.name: axis.size for axis in axes}
            raise SchemaError(
                f"{kind} has a points layout, so all of its coordinate axes must have the same "
                f"length — each index is one sample in {len(axes)} dimensions — but they are "
                f"{sizes}. If the axes are meant to be separable (a grid), use a container kind "
                f"whose LAYOUT is Layout.GRID."
            )
        return tuple(axes)

    @classmethod
    def _check_axis_unit(
        cls, spec: AxisSpec, unit: u.UnitBase | None, kind: str
    ) -> u.UnitBase | None:
        if not spec.unit_required:
            return unit
        if unit is None:
            # A bare array is unambiguous when the axis is dimensionless anyway
            # — (u,v) baselines in wavelengths are plain numbers by convention.
            if "dimensionless" in spec.physical_types:
                return u.dimensionless_unscaled
            wanted = " or ".join(spec.physical_types) or "the right"
            raise SchemaError(
                f"{kind}'s {spec.name!r} coordinates need a unit, but a bare array was given. "
                f"Pass an astropy Quantity of {wanted} physical type, e.g. "
                f"{spec.name}=values * u.{_example_unit(spec)}. Units are converted once here "
                f"and never in the hot loop."
            )
        if any(u.get_physical_type(unit) == wanted for wanted in spec.physical_types):
            return unit
        if any(unit.is_equivalent(other) for other in spec.equivalent_units):
            return unit
        accepted = list(spec.physical_types) + [str(other) for other in spec.equivalent_units]
        raise SchemaError(
            f"{kind}'s {spec.name!r} coordinates are in {unit} (physical type "
            f"'{u.get_physical_type(unit)}'), which this axis does not accept. It accepts "
            f"{accepted}. Check you have not passed the axes in the wrong order."
        )

    @classmethod
    def _expected_shape(cls, axes: tuple[Axis, ...]) -> tuple[int, ...]:
        if cls.LAYOUT is Layout.GRID:
            return tuple(axis.size for axis in axes)
        return (axes[0].size,) if axes else (0,)

    @staticmethod
    def _shape_message(
        kind: str,
        what: str,
        got: tuple[int, ...],
        expected: tuple[int, ...],
        axes: tuple[Axis, ...],
    ) -> str:
        layout = " x ".join(f"{axis.name}={axis.size}" for axis in axes) or "no axes"
        return (
            f"{kind}'s {what} have shape {got}, but its coordinates imply {expected} ({layout}). "
            f"Coordinates and values are aligned index-by-index; a mismatch here is a bug that "
            f"would otherwise surface as an unrelated broadcasting error inside a likelihood."
        )

    def _build_uncertainty(
        self,
        uncertainty: ArrayLike | None,
        value_unit: u.UnitBase | None,
        expected: tuple[int, ...],
        kind: str,
    ) -> np.ndarray | None:
        if uncertainty is None:
            return None
        array, _ = _as_array(uncertainty, f"{kind} uncertainties", unit=value_unit)
        if array.shape != expected:
            raise SchemaError(
                self._shape_message(kind, "uncertainties", array.shape, expected, self.axes)
            )
        if not bool(np.all(array[np.isfinite(array)] >= 0.0)):
            raise SchemaError(
                f"{kind} was given negative uncertainties. An uncertainty is a standard "
                f"deviation; to exclude a sample entirely set its mask, which is the same "
                f"statement as an infinite uncertainty (see prior_art.md lesson R1)."
            )
        return array

    def _build_mask(
        self, mask: ArrayLike | None, expected: tuple[int, ...], kind: str
    ) -> np.ndarray | None:
        if mask is None:
            return None
        array = np.asarray(mask)
        if array.dtype != np.bool_:
            raise SchemaError(
                f"{kind}'s mask must be a boolean array, got dtype {array.dtype}. True marks an "
                f"*excluded* sample, following the numpy and astropy NDData convention; if you "
                f"have an array of good-data flags, pass ~good."
            )
        if array.shape != expected:
            raise SchemaError(self._shape_message(kind, "mask", array.shape, expected, self.axes))
        return _readonly(array)

    def _build_extra_coords(
        self,
        extra_coords: Mapping[str, ArrayLike] | None,
        expected: tuple[int, ...],
        kind: str,
    ) -> Mapping[str, np.ndarray]:
        if not extra_coords:
            return types.MappingProxyType({})
        declared = {spec.name for spec in self.AXES}
        built: dict[str, np.ndarray] = {}
        for name, raw in extra_coords.items():
            if name in declared:
                raise SchemaError(
                    f"{kind}'s extra coordinate {name!r} collides with one of its declared axes "
                    f"{sorted(declared)}. Extra coordinates annotate samples; they never define "
                    f"the sampling geometry."
                )
            if isinstance(raw, u.Quantity):
                # Refused rather than silently stripped: extra coordinates are
                # unitless labels in v1.4, and dropping a unit the caller
                # supplied is exactly the kind of quiet loss this contract is
                # meant to prevent. See results_schema.md §15.4.
                raise SchemaError(
                    f"{kind}'s extra coordinate {name!r} was given as a Quantity in {raw.unit}, "
                    f"but extra coordinates are unitless per-sample labels. Pass "
                    f"{name}=values.to_value({raw.unit}) and record the unit in meta, or make it "
                    f"a declared axis if it defines the sampling geometry."
                )
            array = np.asarray(raw)
            if array.shape != expected:
                raise SchemaError(
                    self._shape_message(
                        kind, f"extra coordinate {name!r}", array.shape, expected, self.axes
                    )
                )
            built[name] = _readonly(array)
        return types.MappingProxyType(built)

    # -- the hot-loop constructor -------------------------------------------

    def with_values(
        self,
        values: ArrayLike,
        *,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        fidelity: str | None = None,
    ) -> Any:
        """Return a copy with new values, reusing the already-validated axes.

        **This is the hot-loop constructor.** Full construction validates
        coordinate ordering, units and physical types in O(N); that is the
        right price once, at composition time, and the wrong price on every
        likelihood evaluation. A model that evaluates on a fixed or negotiated
        grid (``DEVELOPMENT_PLAN.md`` §4.3) should build its container once and
        call this thereafter — the axes, their unit checks and their advertised
        structure are shared, not recomputed.

        The value unit, extra coordinates and metadata are inherited. ``mask``
        and ``uncertainty`` are inherited unless given; pass an explicit
        argument to replace one.
        """
        clone = object.__new__(type(self))
        set_ = object.__setattr__
        kind = type(self).__name__

        # Pass the template's unit through, so a Quantity in a convertible unit
        # is converted rather than silently reinterpreted. Getting this wrong
        # would be the units trap (DEVELOPMENT_PLAN.md §7) in its most damaging
        # form: a factor of 1000 that never announces itself.
        array, resolved = _as_array(
            values, f"{kind} values", unit=self.unit, allow_complex=self.ALLOW_COMPLEX
        )
        if resolved != self.unit:
            raise SchemaError(
                f"with_values() cannot change {kind}'s value unit: the template is in "
                f"{self.unit if self.unit is not None else 'no unit'} but the new values are in "
                f"{resolved}. The unit is a property of the composed container, fixed once at "
                f"construction; build a new container, or pass a plain array in the template's "
                f"unit."
            )
        if array.shape != self.values.shape:
            raise SchemaError(
                self._shape_message(kind, "values", array.shape, self.values.shape, self.axes)
            )
        for name in ("axes", "unit", "extra_coords", "meta"):
            set_(clone, name, getattr(self, name))
        set_(clone, "values", array)
        set_(clone, "fidelity", self.fidelity if fidelity is None else fidelity)
        set_(
            clone,
            "uncertainty",
            self.uncertainty
            if uncertainty is None
            else clone._build_uncertainty(uncertainty, self.unit, array.shape, kind),
        )
        set_(
            clone, "mask", self.mask if mask is None else clone._build_mask(mask, array.shape, kind)
        )
        return clone

    # -- introspection -------------------------------------------------------

    @property
    def shape(self) -> tuple[int, ...]:
        """Shape of :attr:`values`."""
        return self.values.shape

    @property
    def n_samples(self) -> int:
        """Total number of samples, masked ones included."""
        return int(self.values.size)

    @property
    def n_valid(self) -> int:
        """Number of samples not excluded by the mask."""
        return self.n_samples if self.mask is None else int(np.count_nonzero(~self.mask))

    @property
    def is_masked(self) -> bool:
        """Whether any sample is excluded."""
        return self.mask is not None and bool(np.any(self.mask))

    @property
    def valid(self) -> np.ndarray:
        """Boolean array, ``True`` where a sample *is* usable (``~mask``)."""
        if self.mask is None:
            return _readonly(np.ones(self.values.shape, dtype=bool))
        return _readonly(~self.mask)

    @property
    def is_regular(self) -> bool:
        """Whether *every* axis advertises even spacing. Never a requirement."""
        return all(axis.regular for axis in self.axes)

    def axis(self, name: str) -> Axis:
        """The named coordinate axis."""
        for candidate in self.axes:
            if candidate.name == name:
                return candidate
        raise SchemaError(
            f"{type(self).__name__} has no axis {name!r}; its axes are "
            f"{[axis.name for axis in self.axes]}."
        )

    def weights(self) -> np.ndarray:
        """Per-sample weights: 1.0 where usable, 0.0 where masked.

        One of the two equivalent expressions of "excluded means zero
        information" that ``prior_art.md`` lesson R1 records; see
        :meth:`masked_uncertainty` for the other.
        """
        return self.valid.astype(float)

    def masked_uncertainty(self) -> np.ndarray:
        """Uncertainties with ``inf`` at masked samples.

        The infinite-uncertainty expression of the mask (lesson R1), for
        consumers whose algebra divides by sigma² rather than multiplying by a
        weight. Raises if the container has no uncertainties, because there is
        then nothing to inflate.
        """
        if self.uncertainty is None:
            raise SchemaError(
                f"{type(self).__name__} has no uncertainties, so the infinite-uncertainty form "
                f"of its mask is undefined. Use weights() instead, or attach uncertainties."
            )
        inflated = np.array(self.uncertainty, dtype=float, copy=True)
        if self.mask is not None:
            inflated[self.mask] = np.inf
        return _readonly(inflated)

    def quantity(self) -> u.Quantity:
        """The values as a :class:`~astropy.units.Quantity`. Not for the hot loop."""
        return u.Quantity(
            self.values, self.unit if self.unit is not None else u.dimensionless_unscaled
        )

    def to_unit(self, unit: object) -> Any:
        """Return a copy with :attr:`values` converted to ``unit``.

        Converting is a composition-time operation. Doing it per evaluation is
        the units trap in ``DEVELOPMENT_PLAN.md`` §7; do it once, here, and
        hand the converted container to the hot loop.
        """
        target = _check_unit(unit, "the target unit")
        if self.unit is None:
            raise SchemaError(
                f"{type(self).__name__} has no value unit, so it cannot be converted to "
                f"{target}. Declare the unit at construction (unit=...) or pass the values as "
                f"a Quantity."
            )
        try:
            factor = float(self.unit.to(target))
        except u.UnitConversionError as exc:
            raise SchemaError(
                f"{type(self).__name__} values are in {self.unit}, which is not convertible to "
                f"{target}. ({exc})"
            ) from exc
        converted = self.with_values(self.values * factor)
        object.__setattr__(converted, "unit", target)
        if self.uncertainty is not None:
            object.__setattr__(converted, "uncertainty", _readonly(self.uncertainty * factor))
        return converted

    # -- protocol ------------------------------------------------------------

    def __eq__(self, other: object) -> bool:
        if type(other) is not type(self):
            return NotImplemented
        assert isinstance(other, FunctionSamples)
        return (
            self.axes == other.axes
            and self.unit == other.unit
            and self.fidelity == other.fidelity
            and np.array_equal(self.values, other.values)
            and _arrays_equal(self.uncertainty, other.uncertainty)
            and _arrays_equal(self.mask, other.mask)
            and dict(self.meta) == dict(other.meta)
            and set(self.extra_coords) == set(other.extra_coords)
            and all(
                np.array_equal(value, other.extra_coords[name])
                for name, value in self.extra_coords.items()
            )
        )

    __hash__: ClassVar[None] = None

    def __repr__(self) -> str:
        bits = [f"shape={self.shape}"]
        for axis in self.axes:
            span = (
                f"{float(axis.values[0]):g}..{float(axis.values[-1]):g}" if axis.size else "empty"
            )
            bits.append(f"{axis.name}=[{span}]{_unit_suffix(axis.unit)}")
        if self.unit is not None and str(self.unit):
            bits.append(f"unit={self.unit}")
        if self.is_masked:
            bits.append(f"masked={self.n_samples - self.n_valid}")
        if self.fidelity is not None:
            bits.append(f"fidelity={self.fidelity!r}")
        return f"<{type(self).__name__} {' '.join(bits)}>"


def _arrays_equal(left: np.ndarray | None, right: np.ndarray | None) -> bool:
    if left is None or right is None:
        return left is None and right is None
    return bool(np.array_equal(left, right))


def _article(word: str) -> str:
    """``"a"`` or ``"an"``, so error messages read like English."""
    return "an" if word[:1].upper() in "AEIOU" else "a"


def _example_unit(spec: AxisSpec) -> str:
    """A plausible unit to name in a "you forgot the unit" message."""
    examples = {
        "length": "um",
        "frequency": "GHz",
        "energy": "keV",
        "time": "day",
        "angle": "arcsec",
    }
    for physical_type in spec.physical_types:
        if physical_type in examples:
            return examples[physical_type]
    return "dimensionless_unscaled"


def _permute(value: ArrayLike, order: np.ndarray) -> Any:
    """Apply ``order`` to ``value``, keeping a Quantity a Quantity.

    Going through :func:`numpy.asarray` here would strip a
    :class:`~astropy.units.Quantity`'s unit *silently*, so the container would
    be built from bare numbers and lose the unit the caller supplied.
    """
    if isinstance(value, u.Quantity):
        return value[order]
    return np.asarray(value)[order]


def _sorted_copy(
    coordinate: ArrayLike,
    arrays: Mapping[str, ArrayLike | None],
    extra_coords: Mapping[str, ArrayLike] | None = None,
) -> tuple[Any, dict[str, Any], dict[str, Any] | None]:
    """Sort ``coordinate`` ascending and apply the same permutation everywhere.

    Sequences are accepted as well as arrays, so ``from_unsorted`` behaves like
    every other constructor here rather than failing with a bare ``TypeError``
    from numpy's fancy indexing.

    ``extra_coords`` are permuted too: they are per-sample labels, and leaving
    them in declaration order while the samples move would silently misalign
    them with the data they annotate — the exact bug class this module exists
    to end.
    """
    if not isinstance(coordinate, u.Quantity):
        coordinate = np.asarray(coordinate)
    raw = coordinate.value if isinstance(coordinate, u.Quantity) else coordinate
    order = np.argsort(np.asarray(raw), kind="stable")
    permuted = {
        name: (None if value is None else _permute(value, order)) for name, value in arrays.items()
    }
    extras = (
        None
        if extra_coords is None
        else {name: _permute(value, order) for name, value in extra_coords.items()}
    )
    return coordinate[order], permuted, extras


# ---------------------------------------------------------------------------
# The container kinds
# ---------------------------------------------------------------------------


class Spectrum(FunctionSamples):
    """Flux (or any spectral quantity) sampled at explicit spectral coordinates.

    The spectral axis may be a wavelength, a frequency or an energy — whichever
    the model naturally works in — and must be **strictly increasing**, which
    is validated here rather than assumed. Legacy ampere assumed it silently in
    its resampling and covariance construction (``architecture.md`` §7); an
    unsorted or duplicated spectral axis makes a quasiseparable GP solve and a
    dense covariance solve misbehave in ways that look like science problems.

    Spacing is *not* constrained: a channel covering a handful of narrow,
    widely separated line windows is an ordinary ``Spectrum``, and is exactly
    the case §4.2's worked example is built around.
    """

    AXES = (
        AxisSpec(
            "spectral_axis",
            physical_types=("length", "frequency", "energy"),
            order=Order.STRICTLY_INCREASING,
        ),
    )
    LAYOUT = Layout.POINTS

    def __init__(
        self,
        spectral_axis: ArrayLike,
        flux: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"spectral_axis": spectral_axis},
            flux,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @classmethod
    def from_unsorted(
        cls,
        spectral_axis: ArrayLike,
        flux: ArrayLike,
        *,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        **kwargs: Any,
    ) -> Spectrum:
        """Sort the spectral axis ascending, then construct.

        The explicit opt-in for data that arrive in arbitrary order. Duplicated
        coordinates still raise — sorting cannot make them well posed.
        Uncertainties, mask and extra coordinates all follow the permutation.
        """
        axis, reordered, extras = _sorted_copy(
            spectral_axis,
            {"flux": flux, "uncertainty": uncertainty, "mask": mask},
            extra_coords,
        )
        return cls(
            axis,
            reordered["flux"],
            uncertainty=reordered["uncertainty"],
            mask=reordered["mask"],
            extra_coords=extras,
            **kwargs,
        )

    @property
    def spectral_axis(self) -> Axis:
        """The spectral coordinate axis."""
        return self.axis("spectral_axis")

    @property
    def flux(self) -> np.ndarray:
        """Alias for :attr:`~FunctionSamples.values`, in the spectral idiom."""
        return self.values


class PhotometricPoints(FunctionSamples):
    """Broadband photometry: one value per named filter.

    Coordinate order is **explicitly tolerated as arbitrary**. A photometric
    point is identified by its filter, not by its position in an array, and
    real catalogues arrive in whatever order the archive emits. The pivot
    wavelength axis exists so the points can be plotted and so a model can be
    asked for the right spectral range; it is not an index.

    Filter names must be unique — that is the identity of a point — and are
    carried as an extra coordinate, since they are labels rather than
    geometry.
    """

    AXES = (
        AxisSpec(
            "spectral_axis", physical_types=("length", "frequency", "energy"), order=Order.ANY
        ),
    )
    LAYOUT = Layout.POINTS

    def __init__(
        self,
        filters: ArrayLike,
        spectral_axis: ArrayLike,
        flux: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        names = np.asarray(filters, dtype=str)
        if names.ndim != 1:
            raise SchemaError(
                f"PhotometricPoints' filter names must be a one-dimensional sequence, got shape "
                f"{names.shape}."
            )
        unique, counts = np.unique(names, return_counts=True)
        if bool(np.any(counts > 1)):
            repeated = sorted(str(name) for name in unique[counts > 1])
            raise SchemaError(
                f"PhotometricPoints was given repeated filter names {repeated}. A filter name is "
                f"the identity of a photometric point, so duplicates make the channel ambiguous. "
                f"Put repeat observations of one filter in separate channels, or disambiguate the "
                f"names (e.g. 'WISE_W1_epoch1')."
            )
        merged = dict(extra_coords or {})
        merged["filters"] = names
        super().__init__(
            {"spectral_axis": spectral_axis},
            flux,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=merged,
            fidelity=fidelity,
            meta=meta,
        )

    @property
    def filters(self) -> np.ndarray:
        """The filter names, in sample order."""
        return self.extra_coords["filters"]

    @property
    def spectral_axis(self) -> Axis:
        """Pivot wavelengths (or frequencies/energies) of the filters."""
        return self.axis("spectral_axis")

    @property
    def flux(self) -> np.ndarray:
        """Alias for :attr:`~FunctionSamples.values`."""
        return self.values


class TimeSeries(FunctionSamples):
    """A quantity sampled at explicit, strictly increasing times.

    Sampling is arbitrary: real light curves have gaps, uneven cadence and
    bursts of high-rate sampling. Ordering *is* required, for the same
    quasiseparable-solver reason as :class:`Spectrum`.
    """

    AXES = (AxisSpec("time", physical_types=("time",), order=Order.STRICTLY_INCREASING),)
    LAYOUT = Layout.POINTS

    def __init__(
        self,
        time: ArrayLike,
        values: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"time": time},
            values,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @classmethod
    def from_unsorted(
        cls,
        time: ArrayLike,
        values: ArrayLike,
        *,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        **kwargs: Any,
    ) -> TimeSeries:
        """Sort by time ascending, then construct.

        Uncertainties, mask and extra coordinates all follow the permutation.
        """
        axis, reordered, extras = _sorted_copy(
            time,
            {"values": values, "uncertainty": uncertainty, "mask": mask},
            extra_coords,
        )
        return cls(
            axis,
            reordered["values"],
            uncertainty=reordered["uncertainty"],
            mask=reordered["mask"],
            extra_coords=extras,
            **kwargs,
        )

    @property
    def time(self) -> Axis:
        """The time axis."""
        return self.axis("time")


class Image(FunctionSamples):
    """A two-dimensional map on separable, strictly monotonic axes.

    Separable is not the same claim as regular: the axes may be spaced however
    they like, and :attr:`Axis.regular` reports whether a fast path is
    available. Monotonic in *either* direction, because sky axes legitimately
    run both ways.

    ``values`` has shape ``(x, y)`` in the axis order declared here.
    """

    AXES = (
        AxisSpec("x", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
        AxisSpec("y", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
    )
    LAYOUT = Layout.GRID

    def __init__(
        self,
        x: ArrayLike,
        y: ArrayLike,
        values: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"x": x, "y": y},
            values,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @property
    def x(self) -> Axis:
        """The first spatial axis."""
        return self.axis("x")

    @property
    def y(self) -> Axis:
        """The second spatial axis."""
        return self.axis("y")


class Cube(FunctionSamples):
    """A spectral cube: two spatial axes and a spectral axis, separable.

    ``values`` has shape ``(x, y, spectral_axis)``. The spectral axis is
    strictly increasing for the same reason :class:`Spectrum`'s is; the spatial
    axes are strictly monotonic in either direction.
    """

    AXES = (
        AxisSpec("x", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
        AxisSpec("y", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
        AxisSpec(
            "spectral_axis",
            physical_types=("length", "frequency", "energy"),
            order=Order.STRICTLY_INCREASING,
        ),
    )
    LAYOUT = Layout.GRID

    def __init__(
        self,
        x: ArrayLike,
        y: ArrayLike,
        spectral_axis: ArrayLike,
        values: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"x": x, "y": y, "spectral_axis": spectral_axis},
            values,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @property
    def x(self) -> Axis:
        """The first spatial axis."""
        return self.axis("x")

    @property
    def y(self) -> Axis:
        """The second spatial axis."""
        return self.axis("y")

    @property
    def spectral_axis(self) -> Axis:
        """The spectral axis."""
        return self.axis("spectral_axis")


def _visibility_plot_coordinate(axes: Mapping[str, Axis]) -> tuple[np.ndarray, str]:
    """**W5.3.** ``VisibilitySet``'s default plotted coordinate: baseline length.

    ``hypot(u, v)``, labelled with the ``u``/``v`` axes' own unit — they are
    forced equal by construction (both dimensionless-or-``rad**-1``), so
    either's is the right one to show.
    """
    u_axis, v_axis = axes["u"], axes["v"]
    u_values = np.asarray(u_axis.values, dtype=float)
    v_values = np.asarray(v_axis.values, dtype=float)
    return np.hypot(u_values, v_values), format_axis_label("baseline length", u_axis.unit)


class VisibilitySet(FunctionSamples):
    """Complex visibilities at scattered ``(u, v, spectral_axis)`` points.

    The kind that stresses this schema hardest, and the Phase 4 proof modality:
    values are **complex**, the sampling is an irregular point set that no grid
    could describe, and (u,v) coverage is the canonical example of why
    "coordinate + value" beats "array with implicit axes".

    Coordinate order is **explicitly tolerated as arbitrary** — a (u,v) point
    set has no natural order, and imposing one would be a fiction. The same
    tolerance covers the spectral axis, which is why it is :attr:`Order.ANY`
    here and strictly increasing on a :class:`Spectrum`: a visibility table
    arrives sorted by baseline, by time, or by neither, and one wavelength is
    repeated across every baseline of a channel.

    ``uncertainty`` is real and non-negative: it is the per-component standard
    deviation of a circular complex Gaussian, the standard interferometric
    noise model, applying independently to the real and imaginary parts.
    Non-circular noise, and the amplitude/phase (Rice, von Mises) formulations,
    are noise-model concerns and belong to W1.6 — the container deliberately
    does not encode them.

    **The third axis (amended at W4.1, ruled by Peter 2026-09-11.)** ``(u, v)``
    are baselines *in wavelengths*, so the convention already divides by
    lambda; for a **grey** model error that absorbs the wavelength dependence
    entirely, and the two-axis form the freeze shipped was enough. It is not
    enough for a **chromatic** one. Write a missing component as
    ``dI(x, y, lambda) = P(x, y) S(lambda)`` — a patch of sky with a band
    profile — and its effect on a baseline **B** at wavelength lambda is
    ``dV = S(lambda) F(B/lambda)``: sharp in wavelength, smooth in spatial
    frequency. Two samples at the *same* ``(u, v)`` and different wavelengths
    are then near-identical to a kernel that cannot see wavelength and very
    different in truth, so the residual correlation a flexible likelihood is
    meant to model is not expressible without the axis. A kernel sees a
    container's axes and nothing else, which is what makes this a container
    question rather than a kernel one (``phase4_placement_memo.md`` §3.6).

    One wavelength (or frequency, or energy — :class:`Spectrum`'s physical
    types) **per sample**, not a separate grid: the layout is still
    :attr:`Layout.POINTS`, so the third axis has the same length as the other
    two. A monochromatic observation is a constant column, which costs one
    array and keeps every consumer writing one thing rather than two.

    Two consequences worth stating, because neither announces itself:

    * the axes now carry **mixed units** (dimensionless ``u``/``v``, a length
      or frequency ``spectral_axis``), so an isotropic GP over the whole point
      set is refused by :meth:`GPSolver.check_compatible` with its mixed-unit
      message. That refusal is correct — a Euclidean distance across
      wavelengths and baselines is meaningless — and the flexible likelihood
      on visibilities therefore wants a kernel with an ``axes`` selector, so
      that a ``Product`` of a ``(u, v)`` block and a spectral block is
      expressible. That selector is W4.5's.
    * ``visibility`` is the **fourth** positional argument. The axis goes
      before the values so that the positional order matches :attr:`AXES`, as
      it does for every other kind here.
    """

    AXES = (
        AxisSpec("u", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("v", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec(
            "spectral_axis",
            physical_types=("length", "frequency", "energy"),
            order=Order.ANY,
        ),
    )
    LAYOUT = Layout.POINTS
    ALLOW_COMPLEX = True
    PLOT_COORDINATE = staticmethod(_visibility_plot_coordinate)

    def __init__(
        self,
        u_coord: ArrayLike,
        v_coord: ArrayLike,
        spectral_axis: ArrayLike,
        visibility: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"u": u_coord, "v": v_coord, "spectral_axis": spectral_axis},
            visibility,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @property
    def u(self) -> Axis:
        """The u baseline coordinate."""
        return self.axis("u")

    @property
    def v(self) -> Axis:
        """The v baseline coordinate."""
        return self.axis("v")

    @property
    def spectral_axis(self) -> Axis:
        """The per-sample spectral coordinate. One value per visibility."""
        return self.axis("spectral_axis")

    @property
    def visibility(self) -> np.ndarray:
        """Alias for :attr:`~FunctionSamples.values`."""
        return self.values

    def amplitude(self) -> np.ndarray:
        """Visibility amplitudes. A view for plotting, not a likelihood."""
        return _readonly(np.abs(self.values))

    def phase(self) -> np.ndarray:
        """Visibility phases in radians. A view for plotting, not a likelihood."""
        return _readonly(np.angle(self.values))


def _closure_phases_plot_coordinate(axes: Mapping[str, Axis]) -> tuple[np.ndarray, str]:
    """**W5.3.** ``ClosurePhases``'s default plotted coordinate: the longest baseline.

    The triangle's three baseline lengths, ``ij``, ``jk`` and the implied
    ``ki = -(ij + jk)`` (:meth:`ClosurePhases.implied_baseline`'s own
    formula, restated here against an axes mapping rather than a live
    instance), the largest of the three per sample — the one the closure
    phase's own resolution is set by.
    """
    u1, v1 = np.asarray(axes["u1"].values, dtype=float), np.asarray(axes["v1"].values, dtype=float)
    u2, v2 = np.asarray(axes["u2"].values, dtype=float), np.asarray(axes["v2"].values, dtype=float)
    b_ij = np.hypot(u1, v1)
    b_jk = np.hypot(u2, v2)
    b_ki = np.hypot(u1 + u2, v1 + v2)
    longest = np.maximum(np.maximum(b_ij, b_jk), b_ki)
    return longest, format_axis_label("longest baseline", axes["u1"].unit)


class ClosurePhases(FunctionSamples):
    """Closure phase per triangle, in radians, indexed by two of its baselines.

    The sum of the three visibility phases around a closed triangle of
    telescopes. It is the interferometric observable that survives an
    atmosphere: a per-telescope phase error enters two of the three baselines
    with opposite signs and cancels exactly, which is why closure phases are
    fitted where raw phases cannot be.

    A triangle is fixed by **two** of its three baselines — the third is minus
    their sum — so the geometry is four dimensionless axes plus the wavelength
    at which the triangle was measured, one per sample. Wavelength is an axis
    for the reason it is one on :class:`VisibilitySet`, and the same way: see
    that class's "third axis" note.

    **Canonical ordering (required).** A triangle has three representations by
    two of its baselines, and the same closure phase listed as
    ``(u_ij, v_ij, u_jk, v_jk)`` in one file and as
    ``(u_jk, v_jk, u_ki, v_ki)`` in another would look far apart to a kernel
    and would fail the bit-identical coordinate pairing a likelihood needs
    (``transformations.md`` §10). This kind therefore fixes one
    (``phase4_placement_memo.md`` §3.4):

    * telescopes are labelled ``i < j < k`` by the array's own station order;
    * ``(u1, v1)`` is the baseline ``ij`` and ``(u2, v2)`` is the baseline
      ``jk``, both with the sign convention ``b_ij = r_j - r_i`` projected on
      the sky and divided by the wavelength;
    * the third baseline ``ki`` is **implied**, and is
      ``-(u1 + u2, v1 + v2)``; it is never stored, because storing a derived
      quantity is how two files come to disagree;
    * the stored value is
      ``arg(V_ij · V_jk · V_ki) = arg(V(b1) · V(b2) · conj(V(b1 + b2)))``,
      wrapped into ``(-pi, pi]``.

    The sign convention is the second half of the ordering and is just as
    load-bearing: with the opposite convention for ``b_ij`` every closure
    phase changes sign, and a fit to a mirrored source would look just as good
    as a fit to the true one. ``ampere.backends.reference.ClosurePhase``
    produces exactly the convention above, and the conformance battery asserts
    it against a closed form rather than against itself.

    ``uncertainty`` is the per-triangle standard deviation **in radians**, the
    propagated quantity every closure-phase pipeline reports; it becomes
    ``kappa = 1/sigma**2`` under :class:`~ampere.core.VonMisesFamily`.

    The identity labels belong in ``extra_coords``, not in the geometry: a
    ``triangle`` label per sample (``"A0-G1-K0"``) and the per-sample
    ``baseline`` names it was formed from are what a human reads and what
    ties a closure phase back to the visibilities it came from — the pairing
    the plan's design horizon (h) needs kept visible for a joint noise model
    over the two channels.

    Examples
    --------
    >>> import astropy.units as u
    >>> t3 = ClosurePhases(
    ...     [40.0e6, 52.0e6], [10.0e6, -8.0e6],
    ...     [-15.0e6, 11.0e6], [33.0e6, 27.0e6],
    ...     [2.2, 2.2] * u.um,
    ...     [0.31, -0.12] * u.rad,
    ...     uncertainty=[0.02, 0.05] * u.rad,
    ...     extra_coords={"triangle": np.array(["A0-G1-K0", "A0-G1-J3"])},
    ... )
    >>> t3.n_samples, t3.unit
    (2, Unit("rad"))
    >>> t3.implied_baseline()[0].round(1).tolist()
    [-25000000.0, -63000000.0]
    """

    AXES = (
        AxisSpec("u1", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("v1", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("u2", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("v2", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec(
            "spectral_axis",
            physical_types=("length", "frequency", "energy"),
            order=Order.ANY,
        ),
    )
    LAYOUT = Layout.POINTS
    ALLOW_COMPLEX = False
    PLOT_COORDINATE = staticmethod(_closure_phases_plot_coordinate)

    def __init__(
        self,
        u1: ArrayLike,
        v1: ArrayLike,
        u2: ArrayLike,
        v2: ArrayLike,
        spectral_axis: ArrayLike,
        phase: ArrayLike,
        *,
        unit: object = None,
        uncertainty: ArrayLike | None = None,
        mask: ArrayLike | None = None,
        extra_coords: Mapping[str, ArrayLike] | None = None,
        fidelity: str | None = None,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(
            {"u1": u1, "v1": v1, "u2": u2, "v2": v2, "spectral_axis": spectral_axis},
            phase,
            unit=unit,
            uncertainty=uncertainty,
            mask=mask,
            extra_coords=extra_coords,
            fidelity=fidelity,
            meta=meta,
        )

    @property
    def u1(self) -> Axis:
        """The u coordinate of the first baseline, ``ij``."""
        return self.axis("u1")

    @property
    def v1(self) -> Axis:
        """The v coordinate of the first baseline, ``ij``."""
        return self.axis("v1")

    @property
    def u2(self) -> Axis:
        """The u coordinate of the second baseline, ``jk``."""
        return self.axis("u2")

    @property
    def v2(self) -> Axis:
        """The v coordinate of the second baseline, ``jk``."""
        return self.axis("v2")

    @property
    def spectral_axis(self) -> Axis:
        """The per-triangle spectral coordinate."""
        return self.axis("spectral_axis")

    @property
    def phase(self) -> np.ndarray:
        """Alias for :attr:`~FunctionSamples.values`, in radians."""
        return self.values

    def implied_baseline(self) -> tuple[np.ndarray, np.ndarray]:
        """The third baseline ``ki``, ``-(u1 + u2, v1 + v2)``.

        Derived on demand rather than stored, which is the canonical ordering's
        own rule: a stored third baseline is a second statement of the same
        fact, and two statements drift.
        """
        return (
            _readonly(-(self.u1.values + self.u2.values)),
            _readonly(-(self.v1.values + self.v2.values)),
        )


# ---------------------------------------------------------------------------
# ModelResult
# ---------------------------------------------------------------------------


class ModelResult(Mapping[str, FunctionSamples]):
    """What a model produced: an immutable mapping of named channels.

    Channels are *named*, not merely typed, because one model may legitimately
    produce several results of the same kind — a low-resolution SED and a set
    of high-resolution line windows are both :class:`Spectrum`, and an
    instrument must be able to say which one it observes
    (``DEVELOPMENT_PLAN.md`` §4.2).

    The single-channel case stays trivial: pass a bare container and it is
    filed under :data:`DEFAULT_CHANNEL`.

    ``ModelResult`` implements :class:`collections.abc.Mapping`, so ``len``,
    iteration, ``in``, ``keys``/``values``/``items`` and ``.get`` all behave as
    expected. Use :meth:`require` rather than ``[]`` when binding an
    instrument: it checks the kind and explains the mismatch.

    A result may carry the **parameter values that produced it**
    (:attr:`parameters` — ruled 2026-09-01, resolving this contract's open
    question 7): a ``(θ, result)`` pair is then self-contained, which is what
    emulator training sets (``DEVELOPMENT_PLAN.md`` design horizon (c)) and
    provenance want. ``Model.__call__`` (W1.5) attaches the resolved values
    automatically; the coupling is a plain name-to-value mapping, not a
    dependency on the parameter contract's types.
    """

    __slots__ = ("_channels", "_meta", "_parameters")

    def __init__(
        self,
        channels: FunctionSamples | Mapping[str, FunctionSamples],
        *,
        meta: Mapping[str, Any] | None = None,
        parameters: Mapping[str, Any] | None = None,
    ) -> None:
        if isinstance(channels, FunctionSamples):
            channels = {DEFAULT_CHANNEL: channels}
        elif not isinstance(channels, Mapping):
            raise SchemaError(
                f"a ModelResult is built from a single container or a mapping of channel name "
                f"to container, got {type(channels).__name__}. A model that produces one output "
                f"may simply return that container."
            )
        if not channels:
            raise SchemaError(
                "a ModelResult must have at least one channel. A model that produces nothing "
                "cannot be compared with data; if a channel is conditionally absent, say so with "
                "a fully masked container rather than by omitting it."
            )
        built: dict[str, FunctionSamples] = {}
        for name, container in channels.items():
            checked = _check_channel_name(name)
            if not isinstance(container, FunctionSamples):
                raise SchemaError(
                    f"channel {checked!r} holds {type(container).__name__}, which is not a "
                    f"container. Every channel holds a FunctionSamples subclass (Spectrum, "
                    f"PhotometricPoints, Image, Cube, TimeSeries, VisibilitySet, or your own)."
                )
            built[checked] = container
        object.__setattr__(self, "_channels", types.MappingProxyType(built))
        object.__setattr__(self, "_meta", types.MappingProxyType(dict(meta) if meta else {}))
        if parameters is None:
            object.__setattr__(self, "_parameters", None)
        else:
            if not isinstance(parameters, Mapping) or not all(
                isinstance(key, str) for key in parameters
            ):
                raise SchemaError(
                    f"a ModelResult's parameters record is a mapping of parameter name to "
                    f"value — the θ that produced this result — got "
                    f"{type(parameters).__name__}. Model.__call__ attaches it automatically; "
                    f"pass parameters=... only when constructing a result by hand."
                )
            object.__setattr__(self, "_parameters", types.MappingProxyType(dict(parameters)))

    # -- Mapping protocol ----------------------------------------------------

    def __getitem__(self, name: str) -> FunctionSamples:
        try:
            return self._channels[name]
        except KeyError:
            raise ChannelError(self._missing_message(name)) from None

    def __iter__(self) -> Iterator[str]:
        return iter(self._channels)

    def __len__(self) -> int:
        return len(self._channels)

    # -- binding -------------------------------------------------------------

    def require(self, name: str, kind: type[FunctionSamples] | None = None) -> Any:
        """Fetch a channel, checking its kind. The instrument-binding entry point.

        ``DEVELOPMENT_PLAN.md`` §4.2 requires a mismatch here to fail loudly at
        composition time, rather than becoming a confusing shape error inside a
        likelihood later. Both failures raise
        :class:`~ampere.core.exceptions.ChannelError`, so a consumer that wants
        to fall back to an alternative channel can catch one type.
        """
        container = self[name]
        if kind is not None and not isinstance(container, kind):
            got = type(container).__name__
            want = kind.__name__
            raise ChannelError(
                f"channel {name!r} holds {_article(got)} {got}, but {_article(want)} {want} was "
                f"required. Either bind to a channel of the right kind "
                f"({self._names_of_kind(kind)}), or change the model to produce "
                f"{_article(want)} {want} on {name!r}."
            )
        return container

    def _missing_message(self, name: str) -> str:
        available = ", ".join(f"{key!r} ({type(value).__name__})" for key, value in self.items())
        hint = ""
        if isinstance(name, str):
            close = [key for key in self._channels if key.lower() == name.lower()]
            if close:
                hint = f" Did you mean {close[0]!r}?"
            elif set(self._channels) == {DEFAULT_CHANNEL}:
                hint = (
                    f" This result has only the automatic {DEFAULT_CHANNEL!r} channel, which "
                    f"means the model returned a bare container; name its channels explicitly "
                    f"if an instrument needs to bind to one by name."
                )
        return f"no channel named {name!r}. Available channels: {available}.{hint}"

    def _names_of_kind(self, kind: type[FunctionSamples]) -> str:
        names = [key for key, value in self.items() if isinstance(value, kind)]
        return repr(names) if names else "there are none in this result"

    # -- introspection -------------------------------------------------------

    @property
    def meta(self) -> Mapping[str, Any]:
        """Free-form immutable metadata about the evaluation as a whole."""
        return self._meta

    @property
    def parameters(self) -> Mapping[str, Any] | None:
        """The parameter values this result was evaluated at, or ``None``.

        A plain name-to-value mapping — the θ half of the ``(θ, result)``
        pairs emulator training sets and provenance need
        (``DEVELOPMENT_PLAN.md`` design horizon (c); ruled 2026-09-01).
        ``Model.__call__`` attaches it automatically; ``None`` means the
        result was built without one, not that the model has no parameters.
        """
        return self._parameters

    def with_parameters(self, parameters: Mapping[str, Any]) -> ModelResult:
        """Return a copy carrying *parameters* as the values that produced it."""
        return ModelResult(dict(self._channels), meta=self._meta, parameters=parameters)

    @property
    def is_single(self) -> bool:
        """Whether this result carries exactly one channel."""
        return len(self._channels) == 1

    def single(self) -> FunctionSamples:
        """The only container, for the one-channel case.

        Raises if there is more than one, naming them — silently picking the
        first would reintroduce exactly the guess-which-output ambiguity this
        contract exists to remove.
        """
        if not self.is_single:
            raise ChannelError(
                f"single() is only meaningful for a one-channel result, but this one has "
                f"{len(self)} channels: {list(self)}. Ask for the one you mean by name."
            )
        return next(iter(self._channels.values()))

    def kinds(self) -> Mapping[str, type[FunctionSamples]]:
        """Channel name to container type — the composition-time type map."""
        return types.MappingProxyType({key: type(value) for key, value in self.items()})

    def of_kind(self, kind: type[FunctionSamples]) -> Mapping[str, FunctionSamples]:
        """Every channel holding a ``kind`` (or a subclass of it)."""
        return types.MappingProxyType(
            {key: value for key, value in self.items() if isinstance(value, kind)}
        )

    def fidelities(self) -> Mapping[str, str | None]:
        """Channel name to fidelity tag; the design-horizon (a) hook."""
        return types.MappingProxyType({key: value.fidelity for key, value in self.items()})

    def with_channels(self, **channels: FunctionSamples) -> ModelResult:
        """Return a copy with channels added or replaced.

        The immutable-update entry point a transformation chain uses when it
        rewrites one channel and passes the rest through.
        """
        merged = dict(self._channels)
        merged.update(channels)
        return ModelResult(merged, meta=self._meta, parameters=self._parameters)

    def without_channels(self, *names: str) -> ModelResult:
        """Return a copy with the named channels removed."""
        for name in names:
            if name not in self._channels:
                raise ChannelError(self._missing_message(name))
        remaining = {key: value for key, value in self._channels.items() if key not in names}
        return ModelResult(remaining, meta=self._meta, parameters=self._parameters)

    def __repr__(self) -> str:
        entries = ", ".join(
            f"{key}: {type(value).__name__}{'' if value.fidelity is None else f'@{value.fidelity}'}"
            f"[{value.n_samples}]"
            for key, value in self.items()
        )
        return f"<ModelResult {entries}>"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ModelResult):
            return NotImplemented
        return (
            dict(self._channels) == dict(other._channels)
            and dict(self._meta) == dict(other._meta)
            and _parameter_records_equal(self._parameters, other._parameters)
        )

    __hash__: ClassVar[None] = None


def _parameter_records_equal(
    left: Mapping[str, Any] | None, right: Mapping[str, Any] | None
) -> bool:
    """Array-aware comparison of two parameter records (values may be arrays)."""
    if left is None or right is None:
        return left is None and right is None
    if set(left) != set(right):
        return False
    return all(np.array_equal(value, right[name]) for name, value in left.items())


# ---------------------------------------------------------------------------
# The shared diagnostics container
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class AnomalyScore:
    """A coordinate-indexed deficiency map — the shared diagnostics container.

    ``diagnostics.md`` §5's proposal, landed at the freeze (ruled 2026-09-03,
    ``results.md`` §15 R4). It lives in ``ampere.core`` precisely so that
    ``ampere.diagnostics`` (family A, pre-fit RHMF screening) and
    ``ampere.results`` (family C, post-fit GP localisation) can each produce
    one **without either namespace depending on the other** — and plain
    numpy, per ``architecture.md`` §4 rule 1.

    It is deliberately *not* a :class:`FunctionSamples` kind: a score is not
    an observable a model predicts or an instrument transforms — nothing
    binds it to a channel, no likelihood consumes it — it is a statement
    *about* a fit or a collection, indexed by the same coordinates.

    ``provenance`` and ``interpretation_notes`` are **required, non-empty**:
    they are what keeps the shared visual grammar from implying a
    comparability two differently-computed statistics do not have
    (``diagnostics.md`` §5 — two panels are never captioned as
    interchangeable without their provenance shown).

    Parameters
    ----------
    coordinates
        Where each score sits: shape ``(n,)`` for one coordinate axis, or
        ``(n, d)`` for *d* of them. Float64.
    values
        The scores, shape ``(n,)``, float64, finite where retained; a
        documented, comparable range with **higher = more anomalous**.
    mask
        Optional boolean, shape ``(n,)``; ``True`` **excludes** a sample,
        the same convention every container carries.
    provenance
        Which diagnostic family produced it — e.g. ``"rhmf_prefit"``,
        ``"gp_localisation_postfit"``.
    interpretation_notes
        Free text: the caveat that travels with the numbers (family C's
        localisation caveat lives here programmatically).

    Examples
    --------
    >>> score = AnomalyScore(
    ...     coordinates=np.array([1.0, 2.0, 3.0]),
    ...     values=np.array([0.1, 2.4, 0.3]),
    ...     provenance="gp_localisation_postfit",
    ...     interpretation_notes="Amplitude localises deficiency; see docs.",
    ... )
    >>> score.n_samples
    3
    """

    coordinates: np.ndarray
    values: np.ndarray
    provenance: str
    interpretation_notes: str
    mask: np.ndarray | None = None

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        if coordinates.ndim not in (1, 2) or coordinates.shape[0] == 0:
            raise SchemaError(
                f"AnomalyScore coordinates must be (n,) or (n, d) with n >= 1, got shape "
                f"{coordinates.shape}."
            )
        values = np.asarray(self.values, dtype=np.float64).ravel()
        if values.shape[0] != coordinates.shape[0]:
            raise SchemaError(
                f"AnomalyScore holds {coordinates.shape[0]} coordinate(s) but "
                f"{values.shape[0]} value(s); a score is indexed by its coordinates."
            )
        mask = self.mask
        if mask is not None:
            mask = np.asarray(mask, dtype=bool).ravel()
            if mask.shape[0] != values.shape[0]:
                raise SchemaError(
                    f"AnomalyScore mask covers {mask.shape[0]} sample(s) but there are "
                    f"{values.shape[0]}."
                )
        retained = values if mask is None else values[~mask]
        if not np.all(np.isfinite(retained)):
            raise SchemaError(
                "AnomalyScore values must be finite where retained; mask the samples that "
                "have no score."
            )
        for field in ("provenance", "interpretation_notes"):
            text = getattr(self, field)
            if not isinstance(text, str) or not text.strip():
                raise SchemaError(
                    f"AnomalyScore.{field} is required and must be a non-empty string: it is "
                    f"what keeps two differently-computed scores from being read as "
                    f"interchangeable (diagnostics.md §5)."
                )
        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "mask", mask)

    @property
    def n_samples(self) -> int:
        """How many scored samples this map holds."""
        return int(self.values.shape[0])

    def __repr__(self) -> str:
        return (
            f"<AnomalyScore {self.provenance!r}: {self.n_samples} sample(s), "
            f"max {float(np.nanmax(self.values)):.3g}>"
        )
