"""Plain-data forms for containers, ``ModelResult``\\ s and (θ, result) pairs.

``results_schema.md`` §15.7 records that ``ModelResult`` has no ``to_spec`` /
``from_spec``, and routes the decision here: "W1.8 owns results emission and
should decide the format rather than having one imposed here; design horizon
(c)'s *serialisable training sets* for emulators is the requirement to design
against." §17 question 6 repeats it and asks that it be nailed down *before*
Phase 2. This module is that decision, in its minimal form.

Why the helpers live here and not on the containers
---------------------------------------------------
Two reasons, one principled and one procedural. The principled one: a container
is a hot-loop object, and ``results_schema.md`` §10 is emphatic that the way to
make one cheaply is :meth:`~ampere.core.results_schema.FunctionSamples.with_values`
on a template; a serialisation method on it would be the only method on the class
that no evaluation ever calls. The procedural one: ``ampere/core/results_schema.py``
is a merged Phase-1 contract, and adding a method to it is a decision-log matter
(``AGENTS.md`` ground rule 9), whereas a function *over* the public surface is
not. ``docs/design/contracts/results.md`` §12 records the recommendation for
W1.13 either way; nothing below changes if the functions later become methods.

The form
--------
One versioned mapping per container, built only from lists, strings, numbers,
booleans and ``None``:

>>> import numpy as np, astropy.units as u
>>> from ampere.core import Spectrum
>>> spectrum = Spectrum(
...     [1.0, 2.0, 3.0] * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.1] * 3 * u.Jy
... )
>>> encoded = container_to_dict(spectrum)
>>> encoded["kind"], encoded["unit"], encoded["coordinates"]["spectral_axis"]["unit"]
('Spectrum', 'Jy', 'um')
>>> encoded["values"]["data"]
[2.0, 4.0, 6.0]
>>> container_from_dict(encoded) == spectrum
True

Round-tripping is by value, so masks, uncertainties, units, per-sample extra
coordinates and the fidelity tag all survive:

>>> from ampere.core import PhotometricPoints
>>> photometry = PhotometricPoints(
...     ["W1", "W2"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy,
...     uncertainty=[0.1, 0.2] * u.Jy, mask=[False, True], fidelity="cheap",
... )
>>> back = container_from_dict(container_to_dict(photometry))
>>> back == photometry, back.filters.tolist(), back.fidelity
(True, ['W1', 'W2'], 'cheap')

Complex containers are carried as separate real and imaginary parts, because
JSON has no complex number:

>>> from ampere.core import VisibilitySet
>>> visibilities = VisibilitySet(
...     [1.0, 2.0], [3.0, 4.0], [2.2, 2.2] * u.um, [1 + 2j, 3 - 1j]
... )
>>> container_to_dict(visibilities)["values"]["imag"]
[2.0, -1.0]
>>> container_from_dict(container_to_dict(visibilities)) == visibilities
True

A whole :class:`~ampere.core.results_schema.ModelResult` encodes channel by
channel, carrying the θ that produced it (``results_schema.md`` §17 question 7's
ruling) — which is exactly the ``(θ, ModelResult)`` pair design horizon (c) wants:

>>> from ampere.core import ModelResult
>>> result = ModelResult({"blue": spectrum}, parameters={"model.slope": 2.0})
>>> encoded = model_result_to_dict(result)
>>> sorted(encoded["channels"]), encoded["parameters"]
(['blue'], {'model.slope': 2.0})
>>> model_result_from_dict(encoded) == result
True

A kind ampere does not know is refused by name, with the remedy:

>>> container_from_dict({"version": CONTAINER_SCHEMA_VERSION, "kind": "Polarimetry"})
Traceback (most recent call last):
    ...
ampere.core.exceptions.ResultsError: no container kind named 'Polarimetry' is registered...

An out-of-tree kind registers itself once, exactly as it already registers
itself with ``ampere.core``'s other extension points:

>>> from ampere.core import FunctionSamples
>>> @register_kind
... class Polarimetry(FunctionSamples):
...     AXES = Spectrum.AXES
>>> kind_named("Polarimetry") is Polarimetry
True

Notes
-----
Non-finite values (a NaN in a masked sample, an infinite uncertainty) survive
this encoding as Python floats, but :func:`json.dumps` then writes them as the
non-standard ``NaN``/``Infinity`` tokens. That is why
``docs/design/contracts/results.md`` §11 specifies **netCDF**, not JSON, as the
on-disk format for a training set: netCDF carries NaN natively, is already the
project's serialisation (``DEVELOPMENT_PLAN.md`` §4.6), and needs no second
dependency. This module is the in-memory intermediate the writer consumes, and
the format an inline record in a provenance attr uses.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core.results_schema import (
    ClosurePhases,
    Cube,
    FunctionSamples,
    Image,
    ModelResult,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    VisibilitySet,
)

from .exceptions import ResultsError

__all__ = [
    "CONTAINER_SCHEMA_VERSION",
    "container_from_dict",
    "container_to_dict",
    "kind_named",
    "model_result_from_dict",
    "model_result_to_dict",
    "register_kind",
    "registered_kinds",
    "training_pair_from_dict",
    "training_pair_to_dict",
]

#: Bumped whenever the meaning of a key below changes.
#:
#: **2 (W2.8)**: an array-valued θ entry — in a training pair's ``theta`` and in
#: a ``ModelResult``'s ``parameters`` — is now the ``{dtype, data}`` record
#: :func:`_array_to_dict` writes, where it used to be a bare
#: ``ndarray.tolist()``. That is the first of the two losses
#: ``serialisation_review.md`` §4 told the training-set writer to fix: a plate's
#: integer index or a float32 latent block came back as a list of Python floats,
#: and the emulator trained on it never knew. The reader accepts both forms —
#: a bare list still decodes, since a list *is* a plain-data value and refusing
#: one would be pedantry — but the meaning of what ampere **writes** changed,
#: and §5's versioning policy is that this bumps in the same commit.
CONTAINER_SCHEMA_VERSION = 2

_KINDS: dict[str, type[FunctionSamples]] = {}


def register_kind(kind: type[FunctionSamples]) -> type[FunctionSamples]:
    """Make a container kind reconstructible by name. Usable as a decorator.

    Raises
    ------
    ResultsError
        If a *different* class is already registered under that name. Silently
        rebinding a kind would make an archived training set decode into
        something other than what wrote it.
    """
    if not isinstance(kind, type) or not issubclass(kind, FunctionSamples):
        raise ResultsError(
            f"only FunctionSamples subclasses can be registered as container kinds, got {kind!r}."
        )
    existing = _KINDS.get(kind.__name__)
    if existing is not None and existing is not kind:
        raise ResultsError(
            f"a different container kind is already registered as {kind.__name__!r} "
            f"({existing!r}). Rename one of them; a kind name is what an archived record "
            f"decodes through."
        )
    _KINDS[kind.__name__] = kind
    return kind


def kind_named(name: str) -> type[FunctionSamples]:
    """The registered container kind called ``name``."""
    try:
        return _KINDS[name]
    except KeyError:
        raise ResultsError(
            f"no container kind named {name!r} is registered; ampere knows "
            f"{sorted(_KINDS)}. An out-of-tree kind must call "
            f"ampere.results.register_kind(MyKind) before its records can be read back."
        ) from None


def registered_kinds() -> Mapping[str, type[FunctionSamples]]:
    """Every container kind that can currently be reconstructed by name."""
    return dict(_KINDS)


for _kind in (
    FunctionSamples,
    Spectrum,
    PhotometricPoints,
    TimeSeries,
    Image,
    Cube,
    VisibilitySet,
    ClosurePhases,
):
    register_kind(_kind)
del _kind


# ---------------------------------------------------------------------------
# Arrays
# ---------------------------------------------------------------------------


def _array_to_dict(array: np.ndarray) -> dict[str, Any]:
    array = np.asarray(array)
    # ``str(dtype)`` rather than ``dtype.name``: the latter renders a unicode
    # dtype as its bit width (``<U2`` becomes ``"str64"``), which numpy will not
    # parse back. Byte order is normalised to little-endian first, so a record
    # written on one machine decodes on another.
    encoded: dict[str, Any] = {"dtype": str(array.dtype.newbyteorder("<"))}
    if array.dtype.kind == "c":
        encoded["real"] = array.real.tolist()
        encoded["imag"] = array.imag.tolist()
    else:
        encoded["data"] = array.tolist()
    return encoded


def _array_from_dict(encoded: Mapping[str, Any], what: str) -> np.ndarray:
    dtype = np.dtype(encoded["dtype"])
    if dtype.kind == "c":
        try:
            real = np.asarray(encoded["real"], dtype=float)
            imaginary = np.asarray(encoded["imag"], dtype=float)
        except KeyError:
            raise ResultsError(
                f"{what} declares a complex dtype {dtype.name!r} but carries no 'real'/'imag' "
                f"parts."
            ) from None
        return (real + 1j * imaginary).astype(dtype)
    try:
        return np.asarray(encoded["data"], dtype=dtype)
    except KeyError:
        raise ResultsError(f"{what} carries no 'data' entry.") from None


def _values_to_dict(values: Mapping[str, Any]) -> dict[str, Any]:
    """A θ mapping in plain data, **keeping each array's dtype**.

    ``serialisation_review.md`` §4 names this as one of the two losses the
    training-set writer must fix: "θ array values go through ``.tolist()``
    (dtype lost)". A plate's integer index, a boolean indicator or a float32
    latent block all came back as lists of Python floats, and nothing
    downstream could tell. Scalars stay scalars — a Python ``float`` is its own
    plain form and wrapping it in a dtype record would make every record wider
    for nothing.
    """
    encoded: dict[str, Any] = {}
    for name, value in values.items():
        if isinstance(value, np.ndarray):
            encoded[str(name)] = _array_to_dict(value)
        elif isinstance(value, np.generic):
            encoded[str(name)] = _array_to_dict(np.asarray(value))
        else:
            encoded[str(name)] = value
    return encoded


def _values_from_dict(encoded: Mapping[str, Any], what: str) -> dict[str, Any]:
    """The inverse of :func:`_values_to_dict`, tolerant of the pre-W2.8 form.

    A bare list is still decoded — as a list, which is what it was written as —
    because a record that merely stated less than it could is not a corrupt
    record. A ``{dtype, ...}`` mapping is rebuilt as the array it was.
    """
    values: dict[str, Any] = {}
    for name, value in encoded.items():
        if isinstance(value, Mapping) and "dtype" in value:
            values[str(name)] = _array_from_dict(value, f"{what} entry {name!r}")
        else:
            values[str(name)] = value
    return values


def _unit_name(unit: u.UnitBase | None) -> str | None:
    return None if unit is None else str(unit.to_string())


def _quantity(values: np.ndarray, unit: str | None) -> Any:
    return values if unit is None else values * u.Unit(unit)


# ---------------------------------------------------------------------------
# Containers
# ---------------------------------------------------------------------------


def container_to_dict(container: FunctionSamples) -> dict[str, Any]:
    """The plain-data form of one container.

    ``meta`` is carried only when it is already plain data; a container whose
    metadata holds a live object is refused rather than silently stripped, since
    a training set that quietly lost an annotation would compare equal to one
    that never had it.
    """
    if not isinstance(container, FunctionSamples):
        raise ResultsError(
            f"container_to_dict takes a FunctionSamples, got {type(container).__name__}."
        )
    encoded: dict[str, Any] = {
        "version": CONTAINER_SCHEMA_VERSION,
        "kind": type(container).__name__,
        "coordinates": {
            axis.name: {"values": axis.values.tolist(), "unit": _unit_name(axis.unit)}
            for axis in container.axes
        },
        "values": _array_to_dict(container.values),
        "unit": _unit_name(container.unit),
        "uncertainty": (
            None if container.uncertainty is None else _array_to_dict(container.uncertainty)
        ),
        "mask": None if container.mask is None else container.mask.tolist(),
        "extra_coords": {
            name: _array_to_dict(values) for name, values in container.extra_coords.items()
        },
        "fidelity": container.fidelity,
    }
    if container.meta:
        encoded["meta"] = _plain_meta(container.meta, type(container).__name__)
    return encoded


def _plain_meta(meta: Mapping[str, Any], kind: str) -> dict[str, Any]:
    plain: dict[str, Any] = {}
    for key, value in meta.items():
        if value is None or isinstance(value, bool | int | float | str):
            plain[str(key)] = value
        elif isinstance(value, np.ndarray):
            plain[str(key)] = _array_to_dict(value)
        elif isinstance(value, list | tuple):
            plain[str(key)] = list(value)
        else:
            raise ResultsError(
                f"{kind}'s meta entry {key!r} holds a {type(value).__name__}, which has no plain "
                f"form. Serialisable metadata is scalars, strings, lists and arrays; put anything "
                f"else outside the container."
            )
    return plain


def container_from_dict(encoded: Mapping[str, Any]) -> FunctionSamples:
    """Rebuild a container from :func:`container_to_dict`'s output.

    Reconstruction goes through
    :class:`~ampere.core.results_schema.FunctionSamples`'s own ``__init__``, so
    the whole **base** contract is re-checked on the way back in: axis names and
    physical types, coordinate ordering, value and uncertainty shapes, complex
    support, non-negative uncertainties, a strictly boolean mask, and extra
    coordinates aligned with the values. A record that has been tampered with in
    any of those ways fails here rather than inside a likelihood.

    What it does **not** re-check is an invariant a subclass declares in its own
    ``__init__`` rather than through :attr:`~ampere.core.results_schema.FunctionSamples.AXES`
    — today that means exactly one thing, ``PhotometricPoints``' rule that filter
    names are unique. A generic reconstructor cannot call the subclass
    constructors, whose signatures differ per kind by design, so it builds the
    base and inherits the base's checks. That is limitation 12 of
    ``docs/design/contracts/results.md`` §13, and the extension point named
    there is a ``validate()`` classmethod on ``FunctionSamples`` for W1.4 to
    add, which this function would then call. It matters only for a
    hand-edited or corrupted record: anything ampere itself wrote had the
    subclass check applied when it was first built.
    """
    version = encoded.get("version", CONTAINER_SCHEMA_VERSION)
    if version != CONTAINER_SCHEMA_VERSION:
        raise ResultsError(
            f"unsupported container record version {version!r}; this ampere writes and reads "
            f"version {CONTAINER_SCHEMA_VERSION}."
        )
    kind = kind_named(str(encoded["kind"]))
    coordinates = {
        name: _quantity(np.asarray(axis["values"]), axis.get("unit"))
        for name, axis in encoded.get("coordinates", {}).items()
    }
    values = _array_from_dict(encoded["values"], f"{kind.__name__} values")
    uncertainty = encoded.get("uncertainty")
    extra = {
        name: _array_from_dict(payload, f"{kind.__name__} extra coordinate {name!r}")
        for name, payload in encoded.get("extra_coords", {}).items()
    }
    mask = encoded.get("mask")
    container = kind.__new__(kind)
    FunctionSamples.__init__(
        container,
        coordinates,
        values,
        unit=None if encoded.get("unit") is None else u.Unit(encoded["unit"]),
        uncertainty=(
            None
            if uncertainty is None
            else _array_from_dict(uncertainty, f"{kind.__name__} uncertainties")
        ),
        # Deliberately not `dtype=bool`: coercing here would quietly turn a
        # tampered `[0, 2, 0]` into `[False, True, False]`, when the base
        # contract's whole point is that a mask is *strictly* boolean.
        mask=None if mask is None else np.asarray(mask),
        extra_coords=extra or None,
        fidelity=encoded.get("fidelity"),
        meta=encoded.get("meta"),
    )
    return container


# ---------------------------------------------------------------------------
# ModelResult and training pairs
# ---------------------------------------------------------------------------


def model_result_to_dict(result: ModelResult) -> dict[str, Any]:
    """The plain-data form of a whole :class:`~ampere.core.results_schema.ModelResult`."""
    if not isinstance(result, ModelResult):
        raise ResultsError(
            f"model_result_to_dict takes a ModelResult, got {type(result).__name__}."
        )
    encoded: dict[str, Any] = {
        "version": CONTAINER_SCHEMA_VERSION,
        "channels": {name: container_to_dict(result[name]) for name in result},
    }
    if result.parameters is not None:
        encoded["parameters"] = _values_to_dict(result.parameters)
    if result.meta:
        encoded["meta"] = _plain_meta(result.meta, "ModelResult")
    return encoded


def model_result_from_dict(encoded: Mapping[str, Any]) -> ModelResult:
    """Rebuild a :class:`~ampere.core.results_schema.ModelResult`."""
    version = encoded.get("version", CONTAINER_SCHEMA_VERSION)
    if version != CONTAINER_SCHEMA_VERSION:
        raise ResultsError(
            f"unsupported ModelResult record version {version!r}; this ampere writes and reads "
            f"version {CONTAINER_SCHEMA_VERSION}."
        )
    channels = {name: container_from_dict(payload) for name, payload in encoded["channels"].items()}
    parameters = encoded.get("parameters")
    return ModelResult(
        channels,
        meta=encoded.get("meta"),
        parameters=(
            None if parameters is None else _values_from_dict(parameters, "ModelResult parameters")
        ),
    )


def training_pair_to_dict(
    theta: Mapping[str, Any],
    result: ModelResult,
    *,
    failed: bool = False,
    observations: Mapping[str, FunctionSamples] | None = None,
    failure: Any = None,
) -> dict[str, Any]:
    """One ``(θ, ModelResult)`` training pair, as design horizon (c) wants it.

    ``failed`` carries :attr:`~ampere.core.dataset.Simulation.failed` through, so
    an emulator training set records the rejected draws instead of losing them —
    ``inference.md`` §13's "reject-and-record rather than train on garbage" is
    only useful if the record survives to the file.

    ``observations`` and ``failure`` close the second of the two losses
    ``serialisation_review.md`` §4 named: the flag alone said *that* a draw was
    rejected and never *why*, and a budget with a 2 % crash rate is worth
    debugging from its own file rather than by re-running it. ``failure`` is a
    :class:`~ampere.core.dataset.Failure` (anything with ``to_dict``);
    ``observations`` are the noisy draws
    :meth:`~ampere.core.dataset.FittingProblem.simulate` produces with
    ``observe=True``, which an SBI budget needs and an emulator budget does not,
    so they are carried when they exist and absent otherwise.

    >>> import numpy as np, astropy.units as u
    >>> from ampere.core import ModelResult, Spectrum
    >>> pair = training_pair_to_dict(
    ...     {"model.slope": 2.0},
    ...     ModelResult(Spectrum([1.0, 2.0] * u.um, [2.0, 4.0] * u.Jy)),
    ... )
    >>> sorted(pair), pair["failed"]
    (['failed', 'result', 'theta', 'version'], False)

    θ keeps its dtype, which it did not before W2.8:

    >>> pair = training_pair_to_dict(
    ...     {"population.objects.theta": np.array([1, 2], dtype=np.int16)},
    ...     ModelResult(Spectrum([1.0, 2.0] * u.um, [2.0, 4.0] * u.Jy)),
    ... )
    >>> back = training_pair_from_dict(pair)
    >>> back["theta"]["population.objects.theta"].dtype
    dtype('int16')
    """
    encoded: dict[str, Any] = {
        "version": CONTAINER_SCHEMA_VERSION,
        "theta": _values_to_dict(theta),
        "result": model_result_to_dict(result),
        "failed": bool(failed),
    }
    if observations is not None:
        encoded["observations"] = {
            str(label): container_to_dict(container) for label, container in observations.items()
        }
    if failure is not None:
        encoded["failure"] = _failure_to_dict(failure)
    return encoded


def _failure_to_dict(failure: Any) -> dict[str, Any]:
    """A :class:`~ampere.core.dataset.Failure` as plain data, by its own method.

    ``Failure.to_dict`` already exists and is already what a run's provenance
    records (``results.md`` §9), so a training set uses it rather than inventing
    a second encoding of the same five fields. Anything else with a ``to_dict``
    is accepted on the same terms; anything without one is refused, because a
    silently stringified failure record is a failure record nobody can count.
    """
    if not hasattr(failure, "to_dict"):
        raise ResultsError(
            f"a training pair's failure must be an ampere Failure (or anything with to_dict()), "
            f"got {type(failure).__name__}. FittingProblem.simulate puts one on "
            f"Simulation.failure."
        )
    return dict(failure.to_dict())


def training_pair_from_dict(encoded: Mapping[str, Any]) -> dict[str, Any]:
    """Rebuild a training pair from :func:`training_pair_to_dict`'s output.

    The round trip ``serialisation_review.md`` §4 said was owed "when the first
    consumer — the SBI layer or an emulator trainer — lands". Returns a plain
    mapping rather than a class, because that is what the writer consumes and
    what a trainer iterates: ``theta`` (by value, dtypes intact), ``result`` (a
    real :class:`~ampere.core.results_schema.ModelResult`), ``failed``, and —
    where the record carries them — ``observations`` (containers, by value) and
    ``failure`` (the plain-data record; a ``Failure`` is deliberately *not*
    reconstructed, since ``FailureReason`` is a closed vocabulary this module
    would have to import and a record is not a live error).

    >>> import astropy.units as u
    >>> from ampere.core import ModelResult, Spectrum
    >>> spectrum = Spectrum([1.0, 2.0] * u.um, [2.0, 4.0] * u.Jy)
    >>> pair = training_pair_to_dict({"model.slope": 2.0}, ModelResult(spectrum))
    >>> back = training_pair_from_dict(pair)
    >>> back["theta"], back["result"]["default"] == spectrum
    ({'model.slope': 2.0}, True)
    """
    version = encoded.get("version", CONTAINER_SCHEMA_VERSION)
    if version != CONTAINER_SCHEMA_VERSION:
        raise ResultsError(
            f"unsupported training-pair record version {version!r}; this ampere writes and "
            f"reads version {CONTAINER_SCHEMA_VERSION}."
        )
    pair: dict[str, Any] = {
        "theta": _values_from_dict(encoded.get("theta", {}), "training pair theta"),
        "result": model_result_from_dict(encoded["result"]),
        "failed": bool(encoded.get("failed", False)),
    }
    observations = encoded.get("observations")
    if observations is not None:
        pair["observations"] = {
            str(label): container_from_dict(payload) for label, payload in observations.items()
        }
    failure = encoded.get("failure")
    if failure is not None:
        pair["failure"] = dict(failure)
    return pair
