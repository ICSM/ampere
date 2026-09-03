"""What a run has to record about itself to be reproducible.

``DEVELOPMENT_PLAN.md`` §4.6 asks every run to store "run provenance: package
versions, seeds, data hashes, spec hash, in the InferenceData attrs". This
module builds exactly that, as a flat mapping of netCDF-safe scalars, and it
does the hashing that makes those hashes mean something.

Why the spec hash is load-bearing, and not book-keeping
-------------------------------------------------------
``docs/design/lowering.md`` §9.2: numpyro's ``seed`` handler splits its key once
per stochastic site, **in trace order**, rather than folding the site name in.
So the draws a lowered model produces depend on the *order* ampere emits its
parameters in — which means the seed alone does not identify a run. Adding one
parameter to a model changes the draws of every parameter emitted after it. The
merged :meth:`~ampere.core.parameter.ParameterSet.to_spec` is an *ordered* list,
so hashing it captures exactly that order; this is why :func:`provenance_attrs`
hashes the merged set (order preserved) and not a sorted union of the component
sets.

The same hash is the cache key the plan's §7 trap list asks for: "cache keys
must hash the model/prior/data spec so stale artefacts are invalidated
automatically". That is why :func:`problem_fingerprint` covers the likelihood
family, the solver and the kernel as well as the parameters — two runs differing
only in Matérn-3/2 versus squared-exponential, or ``DenseGP`` versus
``QuasisepGP``, have *identical* parameter specs, and a hash that could not tell
them apart would happily serve a stale emulator or SBI posterior.

The recipe, stated once
-----------------------
1. Everything is normalised to a JSON tree by :func:`normalise` — numpy scalars
   become Python scalars, arrays become an ``{dtype, shape, digest}``
   fingerprint, units become their ``to_string()``, enums become their values,
   and a non-finite float becomes one of three sentinel strings so that
   ``allow_nan=False`` can stay on and nothing slips through as the
   non-standard ``NaN`` token.
2. That tree is serialised by :func:`canonical_json`: ``sort_keys=True``,
   no whitespace, UTF-8, ``allow_nan=False``. Mapping order is therefore
   irrelevant and **list order is significant**, which is what the trace-order
   argument above needs.
3. The digest is BLAKE2b with a 16-byte output (32 hex characters) and the
   personalisation string ``b"ampere-prov"`` — :mod:`hashlib`, never
   :func:`hash`, for the same reason ``ampere.core.rng`` gives: Python's
   built-in hash is salted per process and would make a "reproducible" record
   differ between two runs of the same script.
4. Array bytes are normalised to **little-endian, C-contiguous** before
   hashing, and the dtype is recorded by ``dtype.name``, so the digest does not
   depend on the machine's byte order.

Examples
--------
>>> import numpy as np
>>> canonical_json({"b": 1, "a": [2, 3]})
'{"a":[2,3],"b":1}'
>>> canonical_json({"a": [2, 3], "b": 1})     # mapping order is irrelevant
'{"a":[2,3],"b":1}'
>>> canonical_json([1, 2]) == canonical_json([2, 1])   # list order is not
False
>>> digest("hello") == digest("hello")
True
>>> len(digest("hello"))
32
>>> canonical_json(float("inf"))
'"__inf__"'
>>> hash_array(np.arange(3.0)) == hash_array(np.arange(3.0))
True
>>> hash_array(np.arange(3.0)) == hash_array(np.arange(3.0, dtype=np.float32))
False
"""

from __future__ import annotations

import enum
import hashlib
import importlib.metadata as _metadata
import json
import math
import sys
from collections.abc import Mapping, Sequence
from typing import Any

import astropy.units as u
import numpy as np

from ampere.core.dataset import Dataset, FittingProblem
from ampere.core.likelihood import Likelihood
from ampere.core.results_schema import FunctionSamples
from ampere.core.transform import Instrument

from .exceptions import ResultsError

__all__ = [
    "ATTR_PREFIX",
    "DIGEST_BYTES",
    "PROVENANCE_SCHEMA_VERSION",
    "buffer_fingerprint",
    "canonical_json",
    "container_fingerprint",
    "dataset_fingerprint",
    "describe_likelihood",
    "digest",
    "hash_array",
    "hash_container",
    "hash_of",
    "model_fingerprint",
    "normalise",
    "package_versions",
    "problem_fingerprint",
    "provenance_attrs",
    "spec_hashes",
]

#: Bumped whenever the meaning of an ``ampere_*`` attribute changes.
#: 2 (W1.13): ``describe_likelihood`` now fingerprints family- and
#: noise-model-owned buffers too (they were invisible before), so
#: ``ampere_problem_hash`` values differ from schema 1's.
PROVENANCE_SCHEMA_VERSION = 2

#: Every attribute this module writes starts with this, so ampere's provenance
#: never collides with ArviZ's own (``created_at``, ``creation_library``, ...)
#: or with a sampler's.
ATTR_PREFIX = "ampere_"

#: BLAKE2b output length in bytes; 16 gives a 32-character hexadecimal digest.
DIGEST_BYTES = 16

_PERSON = b"ampere-prov"

#: The three non-finite floats, as JSON-legal sentinels. JSON has no ``NaN``;
#: Python's ``json`` will emit the non-standard ``NaN``/``Infinity`` tokens
#: unless ``allow_nan=False``, and those tokens are not portable. Mapping them
#: explicitly keeps ``allow_nan=False`` on, so an unnormalised float is an
#: error rather than a silently non-portable record.
_POSITIVE_INFINITY = "__inf__"
_NEGATIVE_INFINITY = "__-inf__"
_NOT_A_NUMBER = "__nan__"

#: Mapping keys :func:`normalise` reserves for its own encodings. A mapping
#: using one could impersonate the thing it encodes, so one is refused.
_SENTINEL_KEYS = frozenset({"__ndarray__", "__unit__", "__quantity__", "__bytes__", "__str__"})

#: A *string* equal to a float sentinel would impersonate the non-finite
#: float it encodes — normalise(float("nan")) and normalise("__nan__") must
#: not compare equal — so such a string is wrapped rather than passed through.
_FLOAT_SENTINELS = frozenset({_POSITIVE_INFINITY, _NEGATIVE_INFINITY, _NOT_A_NUMBER})

#: Packages whose versions every run records. Anything else the caller adds.
_RECORDED_PACKAGES = (
    "ampere",
    "numpy",
    "scipy",
    "astropy",
    "arviz",
    "arviz-base",
    "xarray",
    "h5netcdf",
    "netCDF4",
    "emcee",
    "dynesty",
    "zeus-mcmc",
)


# ---------------------------------------------------------------------------
# Canonicalisation and hashing
# ---------------------------------------------------------------------------


def hash_array(array: np.ndarray) -> dict[str, Any]:
    """A dtype-, shape- and content-sensitive fingerprint of one array.

    The bytes are taken in **little-endian, C-contiguous** order regardless of
    the machine's own, so the same data hash the same on any platform, and the
    dtype is recorded by :attr:`numpy.dtype.name` for the same reason. Object
    arrays cannot be hashed by their bytes at all — those bytes are pointers —
    so they go through :func:`normalise` on ``tolist()`` instead.

    Examples
    --------
    >>> import numpy as np
    >>> hash_array(np.array([1.0, 2.0]))["shape"]
    [2]
    >>> hash_array(np.array([1.0, 2.0]))["dtype"]
    'float64'
    >>> hash_array(np.array(["W1", "W2"]))["dtype"]
    'str64'
    """
    array = np.asarray(array)
    if array.dtype.kind == "O":
        content = digest(canonical_json(normalise(array.tolist())))
    else:
        ordered = np.ascontiguousarray(array.astype(array.dtype.newbyteorder("<"), copy=False))
        content = digest(ordered.tobytes(order="C"))
    return {"dtype": array.dtype.name, "shape": [int(n) for n in array.shape], "digest": content}


def normalise(obj: object) -> Any:
    """Convert an arbitrary object into a JSON tree with no non-finite floats.

    Raises
    ------
    ResultsError
        If something in the tree has no defined normal form. Refusing is
        deliberate: a provenance record that silently dropped an unrecognised
        field would be worse than no record, because it would compare equal to a
        run that differed in exactly that field.

    Examples
    --------
    >>> import numpy as np, astropy.units as u
    >>> normalise({"n": np.int64(3), "u": u.Jy})
    {'n': 3, 'u': {'__unit__': 'Jy'}}
    >>> normalise(float("-inf"))
    '__-inf__'
    >>> normalise(object())
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.ResultsError: cannot record an object of type 'object'...
    """
    if isinstance(obj, str):
        # A string equal to a float sentinel is wrapped, or it would hash the
        # same as the non-finite float it spells.
        return {"__str__": obj} if obj in _FLOAT_SENTINELS else obj
    if obj is None or isinstance(obj, bool | int):
        return obj
    if isinstance(obj, float):
        return _normalise_float(obj)
    if isinstance(obj, enum.Enum):
        return normalise(obj.value)
    if isinstance(obj, np.ndarray):
        return {"__ndarray__": hash_array(obj)}
    if isinstance(obj, np.generic):
        return normalise(obj.item())
    if isinstance(obj, u.Quantity):
        value = normalise(np.asarray(obj.value))
        return {"__quantity__": {"value": value, "unit": _unit(obj.unit)}}
    if isinstance(obj, u.UnitBase):
        return {"__unit__": _unit(obj)}
    if isinstance(obj, bytes | bytearray):
        return {"__bytes__": bytes(obj).hex()}
    if isinstance(obj, Mapping):
        normalised: dict[str, Any] = {}
        for key, value in obj.items():
            name = str(key)
            if name in _SENTINEL_KEYS:
                # Otherwise a mapping could impersonate an array, a unit or a
                # quantity, and two genuinely different objects would hash the
                # same -- the collision this module's refusal-not-omission rule
                # exists to prevent, arriving by the other door.
                raise ResultsError(
                    f"{name!r} is reserved: it is how a normalised record encodes an array, a "
                    f"unit, a quantity or raw bytes, so a mapping using it as a key could not be "
                    f"told apart from the real thing. Rename the key."
                )
            normalised[name] = normalise(value)
        return normalised
    if isinstance(obj, Sequence):
        return [normalise(item) for item in obj]
    if isinstance(obj, set | frozenset):
        # Sets have no order, so they are recorded as their sorted canonical
        # forms rather than in iteration order, which is not reproducible.
        return sorted(canonical_json(normalise(item)) for item in obj)
    raise ResultsError(
        f"cannot record an object of type {type(obj).__name__!r} in a provenance record: it has "
        f"no defined JSON form. Convert it to a plain mapping, sequence, string or number first, "
        f"or leave it out of the record deliberately."
    )


def _normalise_float(value: float) -> float | str:
    if value == float("inf"):
        return _POSITIVE_INFINITY
    if value == float("-inf"):
        return _NEGATIVE_INFINITY
    if value != value:  # NaN, without importing math for one comparison
        return _NOT_A_NUMBER
    return float(value)


def _unit(unit: u.UnitBase | None) -> str | None:
    return None if unit is None else str(unit.to_string())


def canonical_json(obj: object) -> str:
    """The one serialisation everything in this module hashes.

    ``sort_keys=True`` so mapping order cannot change a digest; no whitespace so
    a digest does not depend on a formatter; ``allow_nan=False`` so a
    non-normalised float is an error rather than a non-portable token.
    """
    return json.dumps(
        normalise(obj),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def digest(payload: str | bytes) -> str:
    """A 32-character BLAKE2b hexadecimal digest of ``payload``."""
    data = payload.encode("utf-8") if isinstance(payload, str) else payload
    return hashlib.blake2b(data, digest_size=DIGEST_BYTES, person=_PERSON).hexdigest()


def hash_of(obj: object) -> str:
    """:func:`digest` of :func:`canonical_json` — the recipe, in one call."""
    return digest(canonical_json(obj))


# ---------------------------------------------------------------------------
# Fingerprints of ampere's own objects
# ---------------------------------------------------------------------------


def container_fingerprint(container: FunctionSamples) -> dict[str, Any]:
    """Everything about an observed container that a run's identity depends on.

    Coordinates, values, uncertainties, mask, per-sample auxiliary arrays, units
    and the fidelity tag — every array by content. ``meta`` is deliberately
    excluded: it is free-form annotation, and a changed comment is not a
    different dataset.
    """
    return {
        "kind": type(container).__name__,
        "axes": [
            {"name": axis.name, "unit": _unit(axis.unit), "values": hash_array(axis.values)}
            for axis in container.axes
        ],
        "values": hash_array(container.values),
        "unit": _unit(container.unit),
        "uncertainty": (
            None if container.uncertainty is None else hash_array(container.uncertainty)
        ),
        "mask": None if container.mask is None else hash_array(container.mask),
        "extra_coords": {
            name: hash_array(values) for name, values in sorted(container.extra_coords.items())
        },
        "fidelity": container.fidelity,
    }


def hash_container(container: FunctionSamples) -> str:
    """The data hash of one observed container."""
    return hash_of(container_fingerprint(container))


def buffer_fingerprint(owner: object) -> list[dict[str, Any]]:
    """Every buffer one :class:`~ampere.core.parameter.Parameterised` declares.

    Buffers are **not** parameters and so appear nowhere in
    :meth:`~ampere.core.parameter.ParameterSet.to_spec` — but they are the
    wavelength grids, opacity tables, filter curves and response matrices a
    model computes with (``architecture.md`` §6), so a run whose buffers moved
    is a different run even though its declaration is identical. Leaving them
    out would make :func:`problem_fingerprint` unusable as the cache key
    ``DEVELOPMENT_PLAN.md`` §7 asks for: an emulator trained against one
    response matrix would be served for a fit against another.

    Hashed by content, in declaration order, which is also the order a backend
    registers them in.
    """
    buffers = getattr(owner, "buffers", None)
    if buffers is None:
        return []
    return [
        {
            "name": buffer.name,
            "unit": _unit(buffer.unit),
            "values": hash_array(buffer.array),
        }
        for buffer in buffers
    ]


# The solver's identity-and-configuration description that used to live here
# (_describe_solver) moved into Likelihood.to_spec with R7's promotion: the
# name alone was never enough — DenseGP's jitter changes the number the same
# theta scores — and the object that knows its own configuration is the one
# that should describe it.


def describe_likelihood(likelihood: Likelihood) -> dict[str, Any]:
    """The family, noise model, solver, kernel and censoring of one likelihood.

    R7 was granted (ruled 2026-09-03, confirmed by W1.13's consolidated
    serialisation review): the declarative assembly this function used to do
    is now :meth:`ampere.core.likelihood.Likelihood.to_spec` — the one
    definition backends and the conformance suite share. What remains here is
    exactly this module's business, per the review's rule that **specs
    describe declarations and provenance fingerprints content**: the buffers'
    content hashes, and the censoring codes' — per-sample data the spec
    deliberately carries only as counts.

    The buffers are fingerprinted on the family and the noise model as well
    as on the ``Likelihood`` itself — a fix the review found: ``Likelihood``
    forwards *parameters* from its pieces, never buffers, so the old
    ``buffer_fingerprint(likelihood)`` alone was blind to a family's
    background template or a noise model's tabulated response, exactly the
    stale-cache trap ``DEVELOPMENT_PLAN.md`` §7 warns about. (The mapping
    change rode a ``PROVENANCE_SCHEMA_VERSION`` bump.)
    """
    described: dict[str, Any] = likelihood.to_spec()
    described["buffers"] = {
        "likelihood": buffer_fingerprint(likelihood),
        "family": buffer_fingerprint(likelihood.family),
        "noise": buffer_fingerprint(likelihood.noise),
    }
    censoring = likelihood.censoring
    if censoring is not None:
        described["censoring"]["kinds"] = hash_array(np.asarray(censoring.kinds))
    return described


def _describe_instrument(instrument: Instrument) -> dict[str, Any]:
    return {
        "label": instrument.label,
        "channel": instrument.channel,
        "input_kind": instrument.input_kind.__name__,
        "steps": [
            {
                "label": step.label,
                "class": type(step).__name__,
                "buffers": buffer_fingerprint(step),
            }
            for step in instrument.steps
        ],
    }


def dataset_fingerprint(dataset: Dataset) -> dict[str, Any]:
    """Everything about one :class:`~ampere.core.dataset.Dataset` a run depends on."""
    return {
        "label": dataset.label,
        "model": dataset.model,
        "channel": dataset.channel,
        "observed": container_fingerprint(dataset.observed),
        "instrument": _describe_instrument(dataset.instrument),
        "likelihood": describe_likelihood(dataset.likelihood),
        "latent_size": None if dataset.latent is None else int(dataset.latent.size),
    }


def model_fingerprint(model: object) -> dict[str, Any]:
    """One model's identity: its class, its declaration and its constant data.

    All three are needed and none is implied by another. Two models of
    *different classes* can declare the same parameters and compute completely
    different things; the same class with the same parameters can be built on a
    different wavelength grid; and the declaration itself is what
    ``lowering.md`` §9.2 makes order-sensitive.
    """
    parameters = getattr(model, "parameters", None)
    return {
        "class": type(model).__name__,
        "module": type(model).__module__,
        "parameters": None if parameters is None else parameters.to_spec(),
        "buffers": buffer_fingerprint(model),
    }


def problem_fingerprint(problem: FittingProblem) -> dict[str, Any]:
    """The whole composition, in the form :func:`hash_of` turns into a cache key.

    Ordered where order matters (the merged parameter spec, the dataset
    sequence) and keyed where it does not.
    """
    return {
        "version": PROVENANCE_SCHEMA_VERSION,
        "parameters": problem.parameters.to_spec(),
        "free_size": int(problem.free_size),
        "models": {label: model_fingerprint(model) for label, model in problem.models.items()},
        "model_bindings": dict(problem.bindings),
        "datasets": [dataset_fingerprint(problem.datasets[label]) for label in problem.datasets],
        "ties": [{"name": tie.name, "sites": list(tie.sites)} for tie in problem.ties],
        "sites": {name: list(paths) for name, paths in problem.sites().items()},
    }


def spec_hashes(problem: FittingProblem) -> dict[str, Any]:
    """The joint spec hash, and one per top-level merge component.

    The ``"spec"`` entry is the load-bearing one — it is the *merged* set, in the
    order ``lowering.md`` §9.2 says a lowered model's seeding depends on. The
    per-component entries are a convenience: if two runs disagree, they say
    where.

    The two live in **separate** namespaces, under ``"spec"`` and
    ``"components"``, rather than in one flat mapping. A component label is a
    user's choice and ``"spec"`` is an ordinary word in this domain, so a
    dataset or model innocently labelled ``spec`` would otherwise overwrite the
    joint entry — and silently, leaving ``ampere_spec_hash`` reporting one
    component's declaration instead of the whole run's, which is exactly the
    property everything downstream relies on it for.
    """
    components: dict[str, str] = {}
    for label, component in problem.datasets.components().items():
        merged = getattr(component, "merged", component)
        components[label] = hash_of(merged.to_spec())
    for label, model in problem.models.items():
        components[label] = hash_of(model.parameters.to_spec())
    return {"spec": hash_of(problem.parameters.to_spec()), "components": components}


# ---------------------------------------------------------------------------
# Versions
# ---------------------------------------------------------------------------


def package_versions(extra: Sequence[str] = ()) -> dict[str, str]:
    """Installed versions of the packages a run's numbers depend on.

    Read from installed distribution metadata rather than by importing, so
    asking what version of jax is present does not import jax.
    """
    versions: dict[str, str] = {"python": ".".join(str(n) for n in sys.version_info[:3])}
    for name in (*_RECORDED_PACKAGES, *extra):
        try:
            versions[name] = _metadata.version(name)
        except _metadata.PackageNotFoundError:
            continue
    return versions


# ---------------------------------------------------------------------------
# The attrs themselves
# ---------------------------------------------------------------------------


def provenance_attrs(
    problem: FittingProblem,
    *,
    engine: str | None = None,
    backend: str = "reference",
    extra: Mapping[str, object] | None = None,
) -> dict[str, Any]:
    """The ``ampere_*`` attributes every emitted run carries.

    Every value is a netCDF-safe scalar: an ``int``, or a ``str`` holding
    canonical JSON for anything structured. Nothing here is large — in
    particular the free-parameter *names* are recorded (one per parameter) and
    the per-element ``free_labels()`` are not, because a 10⁵-element latent
    block's labels are exactly what ``likelihoods.md`` §16 tells this contract
    not to materialise.

    Parameters
    ----------
    problem
        The composed problem the run was over.
    engine
        The sampler or optimiser that produced the draws, e.g. ``"emcee"``.
    backend
        Which rung of the capability ladder executed it.
    extra
        Further entries, JSON-normalised and prefixed like the rest. Use it for
        engine-specific settings (step size, number of live points).

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
    >>> class Line(Model):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength, unit=u.um)
    ...         self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
    ...     def evaluate(self, **values):
    ...         ctx = self.context(values)
    ...         return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)
    >>> grid = np.array([1.0, 2.0, 3.0])
    >>> observed = Spectrum(
    ...     grid * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
    ... )
    >>> problem = FittingProblem(Line(grid), [Dataset(observed)], seed=20260902)
    >>> attrs = provenance_attrs(problem, engine="emcee")
    >>> attrs["ampere_seed"], attrs["ampere_seed_source"]
    (20260902, 'explicit')
    >>> len(attrs["ampere_spec_hash"])
    32
    >>> attrs["ampere_free_names"]
    '["model.slope"]'
    >>> attrs["ampere_log_likelihood_decomposition"]
    'per_dataset'

    The spec hash moves when the declaration moves, which is the whole point —
    a changed prior is a different run even at the same seed:

    >>> class Steeper(Line):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength, unit=u.um)
    ...         self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))
    >>> other = FittingProblem(Steeper(grid), [Dataset(observed)], seed=20260902)
    >>> provenance_attrs(other)["ampere_spec_hash"] == attrs["ampere_spec_hash"]
    False

    and so does the data hash when the data move:

    >>> moved = Spectrum(grid * u.um, [2.0, 4.0, 6.1] * u.Jy, uncertainty=[0.1] * 3 * u.Jy)
    >>> shifted = FittingProblem(Line(grid), [Dataset(moved)], seed=20260902)
    >>> provenance_attrs(shifted)["ampere_data_hash"] == attrs["ampere_data_hash"]
    False
    """
    data_hashes = {
        label: hash_container(problem.datasets[label].observed) for label in problem.datasets
    }
    hashes = spec_hashes(problem)
    attrs: dict[str, Any] = {
        "schema_version": PROVENANCE_SCHEMA_VERSION,
        "backend": backend,
        "engine": "" if engine is None else str(engine),
        "library_versions": canonical_json(package_versions()),
        "spec_hash": hashes["spec"],
        "component_spec_hashes": canonical_json(hashes["components"]),
        "data_hash": hash_of(data_hashes),
        "data_hashes": canonical_json(data_hashes),
        "problem_hash": hash_of(problem_fingerprint(problem)),
        "capabilities": canonical_json(problem.capabilities.to_dict()),
        "free_size": int(problem.free_size),
        "free_names": canonical_json(list(problem.parameters.free_names)),
        "plates": canonical_json(problem.parameters.plates),
        "tied_names": canonical_json(list(problem.tied_names)),
        "sites": canonical_json({name: list(v) for name, v in problem.sites().items()}),
        "dataset_labels": canonical_json(list(problem.datasets)),
        "dataset_channels": canonical_json(
            {label: problem.datasets[label].channel for label in problem.datasets}
        ),
        "model_labels": canonical_json(sorted(problem.models)),
        "model_bindings": canonical_json(dict(problem.bindings)),
        "likelihoods": canonical_json(
            {
                label: describe_likelihood(problem.datasets[label].likelihood)
                for label in problem.datasets
            }
        ),
        "failure_counts": canonical_json(
            {str(reason): int(count) for reason, count in problem.failure_counts.items()}
        ),
        "failures": canonical_json([failure.to_dict() for failure in problem.failures]),
        "log_likelihood_decomposition": "per_dataset",
    }
    if problem.seed is None:
        attrs["seed_source"] = "entropy"
    else:
        attrs["seed"] = int(problem.seed)
        attrs["seed_source"] = "explicit"
    if extra:
        for key, value in extra.items():
            if key in attrs:
                raise ResultsError(
                    f"extra provenance key {key!r} would overwrite the one ampere writes; "
                    f"choose another name."
                )
            attrs[key] = _as_attribute(value)
    return {f"{ATTR_PREFIX}{key}": value for key, value in attrs.items()}


def _as_attribute(value: object) -> int | float | str:
    """Coerce one caller-supplied value into something netCDF can actually store.

    The netCDF types are integers, floats and strings; **a boolean is not one of
    them**, and both engines refuse one — netCDF4 with ``illegal data type for
    attribute``, h5netcdf with ``boolean dtypes are not a supported NetCDF
    feature``. A bare ``isinstance(value, int)`` does not catch this, because
    ``bool`` *is* an ``int`` in Python, so ``adapt=True`` used to pass straight
    through and fail much later, at write time, in a backend traceback. Booleans
    become ``0``/``1``; numpy scalars become their Python equivalents (otherwise
    ``np.int64(64)`` was stringified into ``'64'``, losing its type for no
    reason); a non-finite float is JSON-encoded, since netCDF's handling of one
    in an attribute is not portable; everything else becomes canonical JSON.
    """
    if isinstance(value, bool | np.bool_):
        return int(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, int):
        return value
    if isinstance(value, float):
        return value if math.isfinite(value) else canonical_json(value)
    if isinstance(value, str):
        return value
    return canonical_json(value)
