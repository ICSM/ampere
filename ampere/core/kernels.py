"""Kernels: the covariance algebra, its axis selection, and the term registry.

Split out of :mod:`ampere.core.likelihood` at **W4.5**, where the kernel
surface stopped being two classes and a private one-row translation table and
became an algebra. Everything the flexible likelihood's *covariance* is made
of lives here; :mod:`ampere.core.likelihood` keeps the families, the noise
models and the solver strategies and re-exports these names, so
``from ampere.core import Matern32`` and ``ampere.core.likelihood.Matern32``
both continue to resolve.

Four things are new at W4.5, and each lifts a limitation ``likelihoods.md``
§15 recorded deliberately:

**An algebra** (§15.1). :class:`Sum` and :class:`Product` compose kernels.
A sum of quasiseparable terms *is* quasiseparable — the generators
concatenate and the semiseparable ranks add — so a sum reaches the O(N) path
whenever every term does. A product does not, and is refused there by name.

**An axis selector** (§15.2, ruled 2026-09-11, ``phase4_placement_memo.md``
§3.6 item 3). A kernel may act on a *named subset* of its container's axes:
``Matern32(axes=("u", "v"))`` is an isotropic kernel in the (u, v) plane of a
three-axis ``VisibilitySet``, and
``Product(Matern32(axes=("u", "v")), Matern32(axes=("spectral_axis",)))`` is
the chromatic sky error of the memo's §3.6 — smooth in spatial frequency,
sharp in wavelength. The single-unit rule applies to the *selected* subset,
which is what lets one container carry axes in different units at all. The
selector defaults to ``None`` ("every axis"), so every kernel written before
W4.5 is unchanged, **spec hash included**: :meth:`KernelSpec.to_dict` emits
the ``axes`` key only when a selection was made.

**More terms** (§15.3). :class:`Matern12` (rank 1), :class:`Matern52`
(rank 3), :class:`SHO` (the damped harmonic oscillator, rank 2, ``Q > 1/2``),
:class:`RotationTerm` (a pair of SHOs, rank 4) and :class:`SpectralMixture`
(a sum of SHOs with free frequencies). Each is **exactly** semiseparable —
including Matérn-5/2, which §15.3 recorded as "not exactly quasiseparable"
because it has no exact form in the *celerite basis*
``e^{-c t}(a cos d t + b sin d t)``. It does have one in the *semiseparable*
form the solver underneath actually factorises, by exactly the argument W2.3
used for Matérn-3/2; see :func:`matern52_representation`.

**A public registry** (the horizon note's "users can define kernels; they
cannot make them fast"). :func:`register_quasiseparable_term` is the
lowering/realisation registries' shape — one slot per kernel family, no
silent overwrite, ``override=`` to replace deliberately, ``builtin=`` rows
distinguished, a provenance entry per user row — and it is the *only* route
onto the O(N) path, on every backend.

Nothing here imports torch, jax or celerite2 (``architecture.md`` §4 rule 1).
The backends' kernels subclass these and swap :attr:`Kernel.ops` for their
own array namespace; the built-in term builders are written once, against
that namespace, so the three backends cannot drift.
"""

from __future__ import annotations

import abc
import copy
import dataclasses
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, ClassVar, Protocol

import astropy.units as u
import numpy as np

from .exceptions import LikelihoodError
from .parameter import Identity, Log, Parameter, Parameterised

__all__ = [
    "DTYPE",
    "NUMPY_OPS",
    "SHO",
    "ArrayOps",
    "CeleriteRepresentation",
    "Kernel",
    "KernelSpec",
    "Matern12",
    "Matern32",
    "Matern52",
    "NumpyOps",
    "Product",
    "QuasiseparableTerm",
    "RotationTerm",
    "SpectralMixture",
    "SquaredExponential",
    "StationaryKernel",
    "Sum",
    "TermBuilder",
    "lookup_quasiseparable_term",
    "matern12_representation",
    "matern32_representation",
    "matern52_representation",
    "quasiseparable_families",
    "register_quasiseparable_term",
    "registered_quasiseparable_terms",
    "rotation_representation",
    "sho_representation",
    "sum_representation",
    "term_provenance_entries",
]

#: The dtype every array in the kernel algebra is held in. ``DEVELOPMENT_PLAN.md``
#: §7: GP linear algebra in float32 fails in ways that look like science problems.
#: Defined here rather than in :mod:`ampere.core.likelihood` only because this
#: module is the lower of the two; ``likelihood`` re-exports it unchanged.
DTYPE = np.float64

_SQRT3 = math.sqrt(3.0)
_SQRT5 = math.sqrt(5.0)
_TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------
# Shared array helpers (likelihood.py imports these back)
# ---------------------------------------------------------------------------


def _check_finite(array: np.ndarray, what: str) -> np.ndarray:
    """Refuse NaN or inf. Applies to complex arrays as well as real ones."""
    if not np.all(np.isfinite(array)):
        raise LikelihoodError(
            f"{what} contains non-finite entries. A likelihood cannot be evaluated on NaN or "
            f"inf; mask the affected samples (mask=True excludes them entirely) rather than "
            f"threading sentinels through the arithmetic."
        )
    return array


def _as_float64(array: Any, what: str) -> np.ndarray:
    """Cast to a contiguous float64 array, loudly."""
    return _check_finite(np.ascontiguousarray(np.asarray(array), dtype=DTYPE), what)


def _as_points(array: Any, what: str, *, dimensions: int | None = None) -> np.ndarray:
    """Coerce coordinates to an ``(n, d)`` float64 array of *n points*.

    ``np.atleast_2d`` is the wrong tool here and the reason this helper exists:
    it turns a shape ``(m,)`` array into ``(1, m)`` — **one m-dimensional
    point** — when what a caller passing a bare list of wavelengths means is m
    one-dimensional points. That reading silently broadcasts through the kernel
    and yields a one-element answer instead of an error, so a 1-D input is
    interpreted here as a column, explicitly, and anything ambiguous raises.
    """
    values = np.asarray(array, dtype=DTYPE)
    if values.ndim == 1:
        values = values[:, None]
    elif values.ndim != 2:
        raise LikelihoodError(
            f"{what} must be a 1-D array of coordinates or an (n, d) array of points, but it "
            f"has shape {values.shape}."
        )
    if dimensions is not None and values.shape[1] != dimensions:
        raise LikelihoodError(
            f"{what} has {values.shape[1]} coordinate(s) per point, but the data it is being "
            f"compared against have {dimensions}. Pass an (n, {dimensions}) array"
            + (", or a bare 1-D array of coordinates." if dimensions == 1 else ".")
        )
    return _check_finite(np.ascontiguousarray(values), what)


def _positive(value: Any, name: str, owner: str, *, allow_zero: bool = False) -> float:
    number = float(np.asarray(value, dtype=DTYPE))
    if not math.isfinite(number) or number < 0.0 or (number == 0.0 and not allow_zero):
        bound = ">= 0" if allow_zero else "> 0"
        raise LikelihoodError(
            f"{owner}'s {name!r} must be finite and {bound}, got {number!r}. Declare it with a "
            f"prior supported on the positive half-line and a Log bijection so no sampler can "
            f"propose a value outside it."
        )
    return number


# ---------------------------------------------------------------------------
# The array namespace a kernel builds in
# ---------------------------------------------------------------------------


class ArrayOps(Protocol):
    """The handful of array operations a kernel needs, in one backend's namespace.

    **W4.5's answer to "one builder or three".** Before W4.5 each backend
    carried its own transcription of every kernel's closed form *and* of its
    celerite generators — two copies of the Matérn-3/2 algebra, in
    ``ampere.backends.torch.gp`` and ``ampere.backends.jax.gp``, each with a
    comment saying it "must stay a transcription". Five new families would
    have made that fifteen transcriptions of mathematics that is the same in
    every namespace.

    So the mathematics is written **once**, in terms of this protocol, and a
    backend supplies the namespace. The backend-specific parts are exactly the
    ones that genuinely are backend-specific: how a Python float becomes an
    array of the right dtype on the right device (:meth:`scalar`), and how
    coordinates are read (:meth:`points`). Nothing in a kernel's closed form
    or in its generators is.

    :attr:`Kernel.ops` is an instance attribute rather than a class one,
    because torch places a kernel's arrays per instance (``dtype=``,
    ``device=``).
    """

    def scalar(self, value: Any) -> Any:
        """One number (or tracer, or tensor) in this namespace, ready for arithmetic."""
        ...

    def points(self, coordinates: Any, *, dimensions: int | None = None) -> Any:
        """Coordinates as ``(n, d)``; a bare 1-D input is read as a column of points."""
        ...

    def n_points(self, coordinates: Any) -> int:
        """How many points ``coordinates`` holds."""
        ...

    def separation(self, left: Any, right: Any) -> Any:
        """Euclidean separation between two ``(n, d)`` / ``(m, d)`` point sets."""
        ...

    def zeros(self, n: int) -> Any:
        """``n`` zeros — the separations :meth:`Kernel.diagonal` evaluates at."""
        ...

    def ones_like(self, array: Any) -> Any:
        """An array of ones with ``array``'s shape."""
        ...

    def stack(self, arrays: Sequence[Any], axis: int = -1) -> Any:
        """Stack along a new axis."""
        ...

    def concatenate(self, arrays: Sequence[Any], axis: int = -1) -> Any:
        """Join along an existing axis."""
        ...

    def exp(self, array: Any) -> Any:
        """Elementwise exponential."""
        ...

    def cos(self, array: Any) -> Any:
        """Elementwise cosine."""
        ...

    def sin(self, array: Any) -> Any:
        """Elementwise sine."""
        ...

    def take_columns(self, points: Any, columns: Sequence[int]) -> Any:
        """The named columns of an ``(n, d)`` point set, in the order given."""
        ...


class NumpyOps:
    """:class:`ArrayOps` in numpy. The reference path's namespace, and the default."""

    def scalar(self, value: Any) -> Any:
        return np.asarray(value, dtype=DTYPE)

    def points(self, coordinates: Any, *, dimensions: int | None = None) -> np.ndarray:
        return _as_points(coordinates, "kernel coordinates", dimensions=dimensions)

    def n_points(self, coordinates: Any) -> int:
        return int(self.points(coordinates).shape[0])

    def separation(self, left: Any, right: Any) -> np.ndarray:
        difference = np.asarray(left)[:, None, :] - np.asarray(right)[None, :, :]
        return np.sqrt(np.einsum("ijk,ijk->ij", difference, difference))

    def zeros(self, n: int) -> np.ndarray:
        return np.zeros(n, dtype=DTYPE)

    def ones_like(self, array: Any) -> np.ndarray:
        return np.ones_like(np.asarray(array, dtype=DTYPE))

    def stack(self, arrays: Sequence[Any], axis: int = -1) -> np.ndarray:
        return np.stack([np.asarray(item, dtype=DTYPE) for item in arrays], axis=axis)

    def concatenate(self, arrays: Sequence[Any], axis: int = -1) -> np.ndarray:
        return np.concatenate([np.asarray(item, dtype=DTYPE) for item in arrays], axis=axis)

    def exp(self, array: Any) -> np.ndarray:
        return np.exp(array)

    def cos(self, array: Any) -> np.ndarray:
        return np.cos(array)

    def sin(self, array: Any) -> np.ndarray:
        return np.sin(array)

    def take_columns(self, points: Any, columns: Sequence[int]) -> np.ndarray:
        return np.ascontiguousarray(np.asarray(points)[:, list(columns)])


#: The reference namespace. One instance, because it holds no state.
NUMPY_OPS: NumpyOps = NumpyOps()


# ---------------------------------------------------------------------------
# The neutral declaration
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class KernelSpec:
    """The neutral, serialisable description of a kernel.

    Family name plus ordered hyperparameter names is the minimum W1.9's
    lowering table needs to emit a ``celerite2``/``tinygp``/GPyTorch term, and
    the maximum that translates across all three. It deliberately carries no
    values: the values are :class:`~ampere.core.parameter.Parameter`\\ s, which
    have their own declaration contract.

    **W4.5 adds two optional fields**, and both are omitted from
    :meth:`to_dict` when unused so that every spec hash minted before W4.5 is
    unchanged:

    ``axes``
        The names of the container axes this kernel acts on, or ``None`` for
        "every axis" — the pre-W4.5 meaning and still the default. The
        selection **is** part of the declaration: the same Matérn-3/2 on
        ``("u", "v")`` and on ``("spectral_axis",)`` are different models, so
        they must hash differently.
    ``terms``
        ``(label, spec)`` pairs for a composite (:class:`Sum`,
        :class:`Product`), in **declaration order**. Order-preserving rather
        than canonicalised, even for :class:`Sum`, where addition is
        commutative: the labels are what the children's parameters are
        qualified with (``term0.amplitude``), so reordering a sum renames its
        parameters and *is* a different declaration to everything downstream —
        the sampler's coordinate order, the ArviZ variable names, the stored
        chains. Sorting the tree for the hash would let two runs whose chains
        cannot be compared claim one identity, which is the thing the hash
        exists to prevent.

    Examples
    --------
    >>> KernelSpec("matern32", ("amplitude", "length_scale"), True).to_dict()
    {'family': 'matern32', 'hyperparameters': ['amplitude', 'length_scale'], \
'quasiseparable': True}
    """

    family: str
    hyperparameters: tuple[str, ...]
    quasiseparable: bool
    axes: tuple[str, ...] | None = None
    terms: tuple[tuple[str, KernelSpec], ...] = ()

    def to_dict(self) -> dict[str, Any]:
        """A plain-data form for provenance attrs and spec hashing."""
        described: dict[str, Any] = {
            "family": self.family,
            "hyperparameters": list(self.hyperparameters),
            "quasiseparable": self.quasiseparable,
        }
        if self.axes is not None:
            described["axes"] = list(self.axes)
        if self.terms:
            described["terms"] = [
                {"label": label, "kernel": spec.to_dict()} for label, spec in self.terms
            ]
        return described


def _as_hyperparameter(
    name: str,
    given: Any,
    unit: u.UnitBase | None,
    *,
    positive: bool = True,
) -> Parameter:
    """Coerce a prior, a fixed number or a ready-made Parameter into a Parameter.

    GP hyperparameters are *ordinary parameters*, with ``Log`` bijections —
    ``parameters.md`` §13's instruction to this contract, discharged here in
    one place so no kernel can quietly do it differently.
    """
    if isinstance(given, Parameter):
        if given.name != name:
            raise LikelihoodError(
                f"kernel hyperparameter {name!r} was given a Parameter named {given.name!r}. "
                f"A kernel's hyperparameter names are part of its KernelSpec (W1.9 lowers them "
                f"to term keywords), so they are fixed; rename it with .rename({name!r})."
            )
        return given
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True, unit=unit)
    if hasattr(given, "ppf"):
        bijection = Log() if positive else Identity()
        return Parameter(name, given, unit=unit, bijection=bijection)
    raise LikelihoodError(
        f"kernel hyperparameter {name!r} must be a frozen scipy.stats distribution (a prior), a "
        f"number (held fixed), or an ampere Parameter — got {type(given).__name__}."
    )


# ---------------------------------------------------------------------------
# The kernel ABC
# ---------------------------------------------------------------------------


class Kernel(Parameterised, abc.ABC):
    """A covariance function, declared neutrally.

    A kernel is a *declaration*: a family name (:attr:`FAMILY`) and its
    hyperparameters as ordinary :class:`~ampere.core.parameter.Parameter`\\ s.
    It knows how to build its own dense covariance matrix — that is what
    ``DenseGP`` needs, and what the conformance suite compares every other
    solver against — but a solver is free to ignore ``matrix`` entirely and
    lower :meth:`spec` to a state-space term instead.

    Subclasses declare :attr:`FAMILY`, :attr:`HYPERPARAMETERS`,
    :attr:`QUASISEPARABLE`, and implement :meth:`_covariance`, which receives
    non-negative separations. A non-stationary kernel would override
    :meth:`matrix` instead; the machinery does not assume stationarity anywhere
    outside :meth:`matrix`'s default implementation.

    **The axis selector (W4.5).** ``axes`` names the container axes this
    kernel acts on. ``None`` — the default, and every kernel written before
    W4.5 — means "every axis", so nothing about an existing declaration, its
    spec or its hash changes. A selection restricts both the separation
    :meth:`matrix` measures *and* the single-unit rule a solver enforces, which
    is what lets a three-axis container carry ``u``, ``v`` (dimensionless) and
    ``spectral_axis`` (micron) and still have an isotropic kernel on a
    well-defined subset of them.

    The names are resolved to column indices once, against the container's own
    axis order, by :meth:`for_axes` — a **functional** binding that returns a
    bound copy rather than mutating the declaration, so one kernel object may
    be shared between datasets whose axes differ.

    Examples
    --------
    >>> import scipy.stats as st
    >>> kernel = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))
    >>> kernel.axes is None
    True
    >>> Matern32(0.3, 2.0, axes=("u", "v")).axes
    ('u', 'v')
    """

    #: Neutral family name; the lowering table and the quasiseparable-term
    #: registry are both keyed on it.
    FAMILY: ClassVar[str] = ""
    #: Hyperparameter names, in declaration order.
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ()
    #: Whether this kernel has an exact quasiseparable (celerite-class)
    #: representation, and so admits an exact O(N) solve on ordered 1D data.
    #: A composite sets it per instance from its children (see :class:`Sum`).
    QUASISEPARABLE: ClassVar[bool] = False
    #: Hyperparameters measured in the **coordinate**'s units, checked against
    #: the selected axes' unit at composition. Widened at W4.5 from the
    #: hard-coded ``length_scale`` so :class:`SHO`'s ``period`` is checked too.
    COORDINATE_SCALED: ClassVar[tuple[str, ...]] = ("length_scale",)
    #: Hyperparameters measured in the **data**'s units.
    VALUE_SCALED: ClassVar[tuple[str, ...]] = ("amplitude",)
    #: Whether :meth:`_hyperparameter` validates positivity on the way past.
    #: ``True`` on the reference path, where a bad value must raise; ``False``
    #: on the differentiable backends, where the value may be a tracer and the
    #: declaration (a positive-half-line prior with a ``Log`` bijection) is
    #: what keeps a sampler from proposing one.
    VALIDATES: ClassVar[bool] = True

    # The four capability flags, with the reference path's honest answers.
    # **W3.8** (ruled by Peter 2026-09-08 on W2.4 slice 3's carried finding):
    # a kernel is a capability part now (``Likelihood.capability_parts``),
    # so it declares them like every other composed piece. Until then it
    # declared nothing, and an ``ampere.core`` kernel inside a torch or jax GP
    # noise model was accepted while its amplitude silently got no gradient --
    # the covariance was built in numpy and then converted, which detaches the
    # graph in exactly the hyperparameters a native GP fit exists to fit.

    #: Whether a gradient can be taken through this kernel's covariance.
    DIFFERENTIABLE: ClassVar[bool] = False
    #: Whether it builds a batch of covariances in one call.
    BATCHABLE: ClassVar[bool] = False
    #: Device its arrays live on. Never auto-detected (``architecture.md`` §5).
    DEVICE: ClassVar[str] = "cpu"
    #: Which rung of the capability ladder supplies it (W2.12's fourth flag).
    BACKEND: ClassVar[str] = "reference"

    #: Class-level defaults so a subclass that never calls :meth:`_declare_axes`
    #: (every kernel written before W4.5) behaves exactly as it did.
    _axes: tuple[str, ...] | None = None
    _columns: tuple[int, ...] | None = None

    # -- declaration ---------------------------------------------------------

    def _declare_axes(self, axes: Sequence[str] | None) -> None:
        """Record the axis selection. Called from a concrete kernel's ``__init__``."""
        if axes is None:
            return
        if isinstance(axes, str):
            raise LikelihoodError(
                f"{type(self).__name__}'s axes= takes a sequence of axis names, not the single "
                f"string {axes!r}. Pass a tuple — axes=({axes!r},) — so a one-axis selection "
                f"cannot be confused with a string of one-character names."
            )
        names = tuple(str(name) for name in axes)
        if not names:
            raise LikelihoodError(
                f"{type(self).__name__} was given an empty axes= selection. A kernel acts on at "
                f"least one axis; leave axes unset for 'every axis'."
            )
        if len(set(names)) != len(names):
            raise LikelihoodError(
                f"{type(self).__name__}'s axes= selection {names} repeats an axis name. Each "
                f"selected axis contributes once to the Euclidean separation."
            )
        self._axes = names

    @property
    def axes(self) -> tuple[str, ...] | None:
        """The container axes this kernel acts on, or ``None`` for all of them."""
        return self._axes

    @property
    def ops(self) -> ArrayOps:
        """The array namespace this kernel builds in. numpy unless a backend swaps it."""
        return getattr(self, "_ops", NUMPY_OPS)

    @property
    def terms(self) -> tuple[tuple[str, Kernel], ...]:
        """``(label, kernel)`` children of a composite; empty for a leaf."""
        return ()

    def leaves(self) -> tuple[Kernel, ...]:
        """Every leaf kernel in this tree, in declaration order. ``(self,)`` for a leaf."""
        if not self.terms:
            return (self,)
        found: list[Kernel] = []
        for _, child in self.terms:
            found.extend(child.leaves())
        return tuple(found)

    def spec(self) -> KernelSpec:
        """The neutral description a lowering rule consumes."""
        return KernelSpec(
            family=self.FAMILY,
            hyperparameters=self.HYPERPARAMETERS,
            quasiseparable=self.QUASISEPARABLE,
            axes=self._axes,
            terms=tuple((label, child.spec()) for label, child in self.terms),
        )

    # -- axis binding --------------------------------------------------------

    @property
    def _selective(self) -> bool:
        """Whether anything in this tree names axes, and so needs binding."""
        return self._axes is not None or any(child._selective for _, child in self.terms)

    def for_axes(self, axis_names: Sequence[str]) -> Kernel:
        """This kernel bound to a container's axis order.

        Returns ``self`` unchanged when nothing in the tree names axes, which
        is the pre-W4.5 case and costs nothing; otherwise a bound **copy**
        whose leaves know which columns of the ``(n, d)`` coordinate block they
        act on. Binding is functional rather than in-place so that one kernel
        declaration may be shared between datasets whose containers differ, and
        the result is cached per axis tuple so a per-evaluation call is free
        after the first.
        """
        if not self._selective:
            return self
        names = tuple(axis_names)
        cache: dict[tuple[str, ...], Kernel] = self.__dict__.setdefault("_bound_cache", {})
        bound = cache.get(names)
        if bound is None:
            bound = self._bind(names)
            cache[names] = bound
        return bound

    def _bind(self, axis_names: tuple[str, ...]) -> Kernel:
        """A copy of this kernel with its axis names resolved to column indices."""
        bound = copy.copy(self)
        bound.__dict__["_bound_cache"] = {}
        bound.__dict__["_columns"] = self._resolve_columns(axis_names)
        return bound

    def _resolve_columns(self, axis_names: tuple[str, ...]) -> tuple[int, ...] | None:
        if self._axes is None:
            return None
        missing = [name for name in self._axes if name not in axis_names]
        if missing:
            raise LikelihoodError(
                f"{type(self).__name__} selects axes {self._axes} but the container's axes are "
                f"{axis_names}; it has no {missing[0]!r}. A kernel's axes= names must be the "
                f"container's own axis names."
            )
        return tuple(axis_names.index(name) for name in self._axes)

    def with_ops(self, ops: ArrayOps) -> Kernel:
        """This kernel building in *ops*'s namespace, as a copy.

        A solver, not a kernel, decides where a solve happens: celerite2's
        compiled kernels are float64 on the CPU whatever device a torch kernel
        was placed on, and a reference solver handed a torch kernel must still
        produce numpy. So each solver rebinds the kernel it is given to its own
        namespace before building anything, and a kernel declared in one
        backend used by another's solver computes where that solver can use it,
        rather than one function down where the type error is unintelligible.

        Returns ``self`` when the namespace is already the right one, which is
        the ordinary case (a torch kernel in a torch solver). Otherwise a
        shallow copy: parameters are shared by name, not by identity, and
        nothing here mutates them.
        """
        if self.ops is ops:
            return self
        rebound = copy.copy(self)
        rebound.__dict__["_ops"] = ops
        rebound.__dict__["_bound_cache"] = {}
        if self.terms:
            rebound.__dict__["_terms"] = tuple(
                (label, child.with_ops(ops)) for label, child in self.terms
            )
        return rebound

    def select(self, points: Any) -> Any:
        """The columns of an ``(n, d)`` point set this kernel acts on.

        The identity when no selection was made or when :meth:`for_axes` has
        not been applied — a direct call to :meth:`matrix` with hand-built
        coordinates therefore behaves exactly as it did before W4.5.
        """
        columns = self._columns
        if columns is None:
            return points
        return self.ops.take_columns(points, columns)

    # -- composition-time checks --------------------------------------------

    def selected_axes(self, axis_names: Sequence[str]) -> tuple[str, ...]:
        """The axis names this kernel tree actually uses, in the container's order."""
        names = tuple(axis_names)
        used: set[str] = set()
        for leaf in self.leaves():
            used.update(leaf._axes if leaf._axes is not None else names)
        return tuple(name for name in names if name in used)

    def check_axes(self, observed: Any, *, owner: str) -> None:
        """Composition-time check of the selection and the single-unit rule.

        The rule is unchanged in substance — a Euclidean separation across axes
        in different units is meaningless, so it is refused — and changed in
        *scope*: it applies to each leaf kernel's **own** selected subset. That
        is the whole point of the selector. A container whose axes are ``u``,
        ``v`` (dimensionless) and ``spectral_axis`` (micron) is refused for a
        bare ``Matern32``, exactly as before, and accepted for
        ``Matern32(axes=("u", "v"))``.
        """
        axis_names = tuple(axis.name for axis in observed.axes)
        units = {axis.name: axis.unit for axis in observed.axes}
        kind = type(observed).__name__
        for leaf in self.leaves():
            leaf._resolve_columns(axis_names)
            if leaf._axes is None:
                selected = axis_names
                distinct = {units[name] for name in selected}
                if len(distinct) > 1:
                    named = sorted(str(unit) for unit in distinct)
                    raise LikelihoodError(
                        f"{owner} measures separation as a Euclidean distance across a "
                        f"{kind}'s coordinate axes, but they carry different units {named}. A "
                        f"single isotropic length-scale is meaningless across mixed units; use "
                        f"one axis, or declare a kernel that takes a length-scale per axis."
                    )
                continue
            distinct = {units[name] for name in leaf._axes}
            if len(distinct) > 1:
                named = sorted(str(unit) for unit in distinct)
                raise LikelihoodError(
                    f"{type(leaf).__name__} selects the {kind}'s axes "
                    f"{tuple(leaf._axes)}, which carry different units {named}. A single "
                    f"isotropic length-scale is meaningless across mixed units; select a subset "
                    f"whose axes share one unit, and compose the rest with Product."
                )

    def check_units(self, observed: Any, *, prefix: str = "") -> None:
        """Check hyperparameter units against the axes and values they scale.

        Coordinate-scaled hyperparameters (``length_scale``, and
        :class:`SHO`'s ``period``) are compared against the unit of the axes
        this kernel *selects*, not against the container's first axis: a
        product of an (u, v) kernel and a wavelength kernel has two length
        scales in two units, and both are right.
        """
        axis_unit = self._axis_unit(observed)
        for name in self.COORDINATE_SCALED:
            declared = self._declared_unit(name)
            if declared is not None and declared != axis_unit:
                qualified = f"{prefix}{name}"
                raise LikelihoodError(
                    f"the kernel's {qualified!r} is declared in {declared} but the "
                    f"{type(observed).__name__}'s coordinate axis is in {axis_unit}. Priors are "
                    f"numeric in the declared unit and rescaling a distribution correctly is "
                    f"family-specific, so this contract requires an exact match rather than a "
                    f"conversion (parameters.md §8 makes the same ruling for tying)."
                )
        for name in self.VALUE_SCALED:
            declared = self._declared_unit(name)
            if declared is not None and declared != observed.unit:
                qualified = f"{prefix}{name}"
                raise LikelihoodError(
                    f"the kernel's {qualified!r} is declared in {declared} but the "
                    f"{type(observed).__name__}'s values are in {observed.unit}. The amplitude "
                    f"is a marginal standard deviation in the data's own units; declare it in "
                    f"{observed.unit} or leave its unit unset."
                )
        for label, child in self.terms:
            child.check_units(observed, prefix=f"{prefix}{label}.")

    def _declared_unit(self, name: str) -> u.UnitBase | None:
        if name not in self.parameters:
            return None
        return self.parameters[name].unit

    def _axis_unit(self, observed: Any) -> u.UnitBase | None:
        if self._axes is None:
            return observed.axes[0].unit if observed.axes else None
        units = {axis.name: axis.unit for axis in observed.axes}
        return units.get(self._axes[0])

    # -- evaluation ----------------------------------------------------------

    def resolve(self, values: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """This kernel's own hyperparameter values, out of a wider mapping.

        A ``NoiseModel`` hands round one flat mapping covering itself and its
        kernel; the kernel picks out its own names rather than requiring the
        caller to split them.
        """
        if values is None:
            return self.context(None)
        selected = {name: values[name] for name in self.parameters.names if name in values}
        return self.context(selected)

    def _hyperparameter(self, values: Mapping[str, Any], name: str, **kwargs: Any) -> Any:
        """One resolved hyperparameter, validated where the backend wants it validated.

        The check is run for its exception and its result discarded: what comes
        back is the *given* value in this kernel's namespace, so a kernel whose
        namespace is a differentiable one keeps the graph. ``VALIDATES`` is
        ``False`` on the differentiable backends' kernels, where the value may
        be a tracer and the positive-half-line prior is what rules it out.
        """
        value = values[name]
        if self.VALIDATES:
            _positive(value, name, self.FAMILY, **kwargs)
        return self.ops.scalar(value)

    @abc.abstractmethod
    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        """Covariance at non-negative separations, given resolved values."""

    def value(self, separation: Any, values: Mapping[str, Any]) -> Any:
        """``k(τ)`` at given non-negative separations: the stationary closed form.

        The scalar face of :meth:`matrix`, for a caller that already has
        separations rather than coordinates — celerite2's ``Term.get_value``,
        and any plot of the kernel itself. A composite combines its children's
        (which presumes they act on one coordinate; a product across axes has
        no single ``τ``, and says so).
        """
        return self._covariance(self.ops.scalar(separation), self.resolve(values))

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> Any:
        """Dense covariance between two coordinate sets.

        ``left`` and ``right`` are ``(n, d)`` and ``(m, d)`` arrays; a bare 1-D
        array is read as a column of ``n`` one-dimensional points. The result
        is ``(n, m)``. Separation is Euclidean across the axes this kernel
        **selects**, which is why ``GPSolver.check_compatible`` requires those
        axes to share one unit.
        """
        ops = self.ops
        resolved = self.resolve(values)
        points = ops.points(left)
        other = ops.points(right, dimensions=int(np.shape(points)[1]))
        separation = ops.separation(self.select(points), self.select(other))
        return self._covariance(separation, resolved)

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> Any:
        """The prior variance at each coordinate; ``k(0)`` for a stationary kernel."""
        ops = self.ops
        resolved = self.resolve(values)
        return self._covariance(ops.zeros(ops.n_points(coordinates)), resolved)

    def __repr__(self) -> str:
        declared = ", ".join(repr(self.parameters[name]) for name in self.parameters.names)
        selection = "" if self._axes is None else f", axes={self._axes!r}"
        return f"{type(self).__name__}({declared}{selection})"


# ---------------------------------------------------------------------------
# The Matérn family
# ---------------------------------------------------------------------------


class StationaryKernel(Kernel):
    """A stationary kernel with an ``amplitude`` and a ``length_scale``.

    The shared constructor of the Matérn family and the squared exponential,
    public since **W4.5** because a user writing their own kernel of that
    shape should not have to re-declare the two hyperparameters, their units
    and their ``Log`` bijections: subclass this, set ``FAMILY`` and
    ``QUASISEPARABLE``, implement ``_covariance``, and — if the kernel is
    exactly quasiseparable — register its representation with
    :func:`register_quasiseparable_term`.
    """

    def __init__(
        self,
        amplitude: Any,
        length_scale: Any,
        *,
        amplitude_unit: u.UnitBase | None = None,
        length_scale_unit: u.UnitBase | None = None,
        axes: Sequence[str] | None = None,
    ) -> None:
        self.register_parameters(
            _as_hyperparameter("amplitude", amplitude, amplitude_unit),
            _as_hyperparameter("length_scale", length_scale, length_scale_unit),
        )
        self._declare_axes(axes)


class Matern12(StationaryKernel):
    r"""Matérn-1/2 (the Ornstein-Uhlenbeck / exponential kernel): rank 1, exact.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r}{\ell}\right)

    The roughest of the Matérn family — its sample paths are continuous but
    nowhere differentiable — and the cheapest quasiseparable term there is: it
    *is* a single celerite real term, so its semiseparable representation has
    rank 1 (:func:`matern12_representation`) and needs no centring at all.

    Added at W4.5 as the first of §15.3's "obvious additions". It is the right
    component for a residual with no smoothness to speak of — detector-level
    pixel-to-pixel structure, say — and the wrong one for the structured model
    deficiencies ``DEVELOPMENT_PLAN.md`` §2 made Matérn-3/2 the default for.

    Examples
    --------
    >>> float(Matern12(2.0, 1.0).matrix([[0.0]], [[0.0]], {})[0, 0])
    4.0
    """

    FAMILY: ClassVar[str] = "matern12"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = True

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(values, "length_scale")
        return amplitude * amplitude * self.ops.exp(-self.ops.scalar(separation) / length_scale)


class Matern32(StationaryKernel):
    r"""Matérn-3/2: ampere's canonical flexible-likelihood kernel.

    .. math::
        k(r) = a^2 \left(1 + \frac{\sqrt{3}\,r}{\ell}\right)
               \exp\!\left(-\frac{\sqrt{3}\,r}{\ell}\right)

    ``DEVELOPMENT_PLAN.md`` §2 and §4.4 make this the default throughout,
    replacing legacy's hardcoded RBF, for two independent reasons. It
    represents structured residuals better than a squared exponential (a
    once-differentiable sample path, not an analytic one — real model
    deficiencies are not infinitely smooth); and it is **exactly
    quasiseparable**, which is what makes ``QuasisepGP`` an *exact* O(N) solve
    rather than an approximation. ``prior_art.md`` lesson S1 records that
    Starfish (Czekala et al. 2015) independently arrived at exactly this
    kernel, in velocity separation, for exactly this purpose.

    ``amplitude`` is the marginal **standard deviation** — ``k(0) ==
    amplitude²`` — following celerite2's ``Matern32Term(sigma=..., rho=...)``
    convention, so a prior on it is a prior in the data's own units.

    Parameters
    ----------
    amplitude, length_scale
        A frozen ``scipy.stats`` prior (fitted, with a :class:`Log` bijection),
        a number (held fixed), or a ready-made ``Parameter``.
    amplitude_unit
        Unit of ``amplitude``; must match the observed values' unit.
    length_scale_unit
        Unit of ``length_scale``; must match the coordinate axis's unit.
    axes
        The container axes this kernel acts on (W4.5). ``None`` — the default
        — means every axis, which is what every kernel meant before W4.5.

    Examples
    --------
    >>> import scipy.stats as st
    >>> import astropy.units as u
    >>> kernel = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(1e-2, 1e2))
    >>> kernel.spec()
    KernelSpec(family='matern32', hyperparameters=('amplitude', 'length_scale'),
               quasiseparable=True, axes=None, terms=())
    >>> kernel.parameters.free_names
    ('amplitude', 'length_scale')
    >>> kernel.parameters.bijections()
    (Log(lower=0.0), Log(lower=0.0))
    >>> float(kernel.matrix([[0.0]], [[0.0]], {"amplitude": 2.0, "length_scale": 1.0})[0, 0])
    4.0
    """

    FAMILY: ClassVar[str] = "matern32"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = True

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(values, "length_scale")
        scaled = _SQRT3 * self.ops.scalar(separation) / length_scale
        return amplitude * amplitude * (1.0 + scaled) * self.ops.exp(-scaled)


class Matern52(StationaryKernel):
    r"""Matérn-5/2: rank 3, and **exact** — see :func:`matern52_representation`.

    .. math::
        k(r) = a^2 \left(1 + \frac{\sqrt5\,r}{\ell}
                 + \frac{5 r^2}{3\ell^2}\right)
               \exp\!\left(-\frac{\sqrt5\,r}{\ell}\right)

    Twice differentiable, so smoother than Matérn-3/2 and rougher than the
    squared exponential. ``likelihoods.md`` §15.3 recorded it as *not* exactly
    quasiseparable; W4.5 lifts that, and the correction is worth stating
    plainly because it is the same correction W2.3 made for Matérn-3/2. There
    is no exact Matérn-5/2 in the **celerite basis**
    ``e^{-ct}(a cos dt + b sin dt)``, which has no ``t² e^{-ct}`` member; but
    the *solver* underneath factorises any rank-J **semiseparable** matrix, and
    a degree-2 polynomial in ``t_n - t_m`` is a rank-3 bilinear form in
    ``(1, t, t²)``. §15.3's claim was true of the basis and false of the
    solver.
    """

    FAMILY: ClassVar[str] = "matern52"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = True

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(values, "length_scale")
        scaled = _SQRT5 * self.ops.scalar(separation) / length_scale
        polynomial = 1.0 + scaled + scaled * scaled / 3.0
        return amplitude * amplitude * polynomial * self.ops.exp(-scaled)


class SquaredExponential(StationaryKernel):
    r"""Squared exponential (RBF): legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Provided so the misspecification study in milestone M2 can compare the new
    default against what legacy ampere actually did, and so a user who wants it
    can have it. It is **not quasiseparable** — an analytic sample path has no
    finite-order state-space form — so ``QuasisepGP`` refuses it, and it
    therefore does not scale past ``DenseGP``. That asymmetry is the concrete
    reason ``DEVELOPMENT_PLAN.md`` §2 made Matérn the default.
    """

    FAMILY: ClassVar[str] = "squared_exponential"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "length_scale")
    QUASISEPARABLE: ClassVar[bool] = False

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(values, "length_scale")
        scaled = self.ops.scalar(separation) / length_scale
        return amplitude * amplitude * self.ops.exp(-0.5 * scaled * scaled)


# ---------------------------------------------------------------------------
# The damped oscillators
# ---------------------------------------------------------------------------


class SHO(Kernel):
    r"""A damped simple harmonic oscillator: celerite's ``SHOTerm`` with ``Q > 1/2``.

    .. math::
        k(\tau) = a^2 e^{-\tau \omega_0 / 2Q}
                  \left[\cos(\eta\,\omega_0 \tau)
                        + \frac{1}{\sqrt{4Q^2 - 1}}\sin(\eta\,\omega_0 \tau)\right],
        \qquad \eta = \sqrt{1 - \tfrac{1}{4Q^2}},

    with :math:`\omega_0 = 2\pi/P` and ``k(0) = a²``, so ``amplitude`` is the
    marginal standard deviation in the data's own units as it is for every
    other ampere kernel.

    **This is the component for fringing** — the case the horizon note's §1
    follow-up names — and more generally for any *quasi-periodic* residual: a
    ripple whose period is set by an optical path difference and whose
    coherence decays over a few periods is exactly a damped oscillator, and
    absorbing it into a stationary Matérn costs a length scale that is either
    too short to see the ripple's coherence or too long to localise it.
    ``quality`` is how many periods the ripple stays coherent for, to within a
    factor of :math:`\pi`.

    ``Q > 1/2`` (underdamped) is required, not defaulted: at ``Q = 1/2`` the
    representation is critically damped and its celerite form degenerates
    (:math:`\sqrt{4Q^2-1} \to 0`), and below it the kernel is a sum of two real
    exponentials — a different family, and one a sum of two
    :class:`Matern12`\\ s already expresses.

    Parameters
    ----------
    amplitude
        Marginal standard deviation, in the data's units.
    period
        The **undamped natural** period :math:`P = 2\pi/\omega_0`, in the
        coordinate's units — celerite2's ``rho``. The damped oscillation's
        period is :math:`P/\eta`, longer by a fraction of order
        :math:`1/8Q^2`, and equal to it for all practical purposes once
        ``Q`` is more than a few.
    quality
        The oscillator's quality factor ``Q``, dimensionless, ``> 1/2``.
    axes
        As :class:`Matern32`.

    Examples
    --------
    >>> float(SHO(0.5, 2.0, 4.0).matrix([[0.0]], [[0.0]], {})[0, 0])
    0.25
    """

    FAMILY: ClassVar[str] = "sho"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = ("amplitude", "period", "quality")
    QUASISEPARABLE: ClassVar[bool] = True
    COORDINATE_SCALED: ClassVar[tuple[str, ...]] = ("period",)
    VALUE_SCALED: ClassVar[tuple[str, ...]] = ("amplitude",)

    def __init__(
        self,
        amplitude: Any,
        period: Any,
        quality: Any,
        *,
        amplitude_unit: u.UnitBase | None = None,
        period_unit: u.UnitBase | None = None,
        axes: Sequence[str] | None = None,
    ) -> None:
        self.register_parameters(
            _as_hyperparameter("amplitude", amplitude, amplitude_unit),
            _as_hyperparameter("period", period, period_unit),
            _as_hyperparameter("quality", quality, None),
        )
        self._declare_axes(axes)

    def coefficients(self, values: Mapping[str, Any]) -> tuple[Any, Any, Any, Any]:
        """``(a, b, c, d)`` of the celerite form ``e^{-c t}(a cos d t + b sin d t)``.

        The single place the SHO algebra lives: :meth:`_covariance` evaluates
        the closed form from it and :func:`sho_representation` builds the
        generators from it, so the dense and the O(N) paths cannot drift.
        """
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        period = self._hyperparameter(values, "period")
        quality = self._hyperparameter(values, "quality")
        if self.VALIDATES and float(np.asarray(quality)) <= 0.5:
            raise LikelihoodError(
                f"SHO's 'quality' must be > 1/2 (an underdamped oscillator), got "
                f"{float(np.asarray(quality))!r}. At Q = 1/2 the term is critically damped and "
                f"its celerite representation degenerates; below it the kernel is a sum of two "
                f"real exponentials, which Sum(Matern12(...), Matern12(...)) expresses."
            )
        omega0 = _TWO_PI / period
        # ** 0.5 rather than a namespace sqrt: every backend's arrays carry it,
        # and it keeps ArrayOps to the operations that genuinely differ.
        split = (4.0 * quality * quality - 1.0) ** 0.5
        decay = 0.5 * omega0 / quality
        marginal = amplitude * amplitude
        return marginal, marginal / split, decay, decay * split

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        ops = self.ops
        cosine, sine, decay, frequency = self.coefficients(values)
        tau = ops.scalar(separation)
        oscillation = cosine * ops.cos(frequency * tau) + sine * ops.sin(frequency * tau)
        return ops.exp(-decay * tau) * oscillation


class RotationTerm(Kernel):
    r"""Two SHOs, at a period and at its first harmonic: celerite2's ``RotationTerm``.

    The second form the item asks for beside :class:`SHO`, and celerite2's own
    parameterisation, transcribed: a mixture of two underdamped oscillators
    whose *damped* periods are exactly ``period`` and ``period/2``, with the
    power split between them by ``fraction``. Where a single :class:`SHO`
    describes a sinusoidal ripple, this describes a **non-sinusoidal** periodic
    one — a stellar rotation signal with spots on both hemispheres is the case
    it was invented for, and an interference fringe seen through a
    non-sinusoidal blaze is the same shape.

    ``k(0) = amplitude²``, as everywhere else.

    Parameters
    ----------
    amplitude
        Marginal standard deviation, in the data's units.
    period
        The damped period of the fundamental, in the coordinate's units.
    quality
        ``Q0``: the *excess* quality factor of the harmonic over the critical
        ``1/2``, so the harmonic's ``Q`` is ``1/2 + quality``. Positive.
    delta_quality
        How much more coherent the fundamental is than the harmonic; the
        fundamental's ``Q`` is ``1/2 + quality + delta_quality``. Non-negative.
    fraction
        The harmonic's power as a fraction of the fundamental's. Non-negative;
        zero makes this a single SHO in all but rank.
    axes
        As :class:`Matern32`.
    """

    FAMILY: ClassVar[str] = "rotation"
    HYPERPARAMETERS: ClassVar[tuple[str, ...]] = (
        "amplitude",
        "period",
        "quality",
        "delta_quality",
        "fraction",
    )
    QUASISEPARABLE: ClassVar[bool] = True
    COORDINATE_SCALED: ClassVar[tuple[str, ...]] = ("period",)
    VALUE_SCALED: ClassVar[tuple[str, ...]] = ("amplitude",)

    def __init__(
        self,
        amplitude: Any,
        period: Any,
        quality: Any,
        delta_quality: Any = 0.0,
        fraction: Any = 0.5,
        *,
        amplitude_unit: u.UnitBase | None = None,
        period_unit: u.UnitBase | None = None,
        axes: Sequence[str] | None = None,
    ) -> None:
        self.register_parameters(
            _as_hyperparameter("amplitude", amplitude, amplitude_unit),
            _as_hyperparameter("period", period, period_unit),
            _as_hyperparameter("quality", quality, None),
            _as_hyperparameter("delta_quality", delta_quality, None),
            _as_hyperparameter("fraction", fraction, None),
        )
        self._declare_axes(axes)

    def coefficients(self, values: Mapping[str, Any]) -> tuple[tuple[Any, Any, Any, Any], ...]:
        """The two celerite ``(a, b, c, d)`` quadruples, fundamental first.

        celerite2's ``RotationTerm`` algebra, transcribed. With
        :math:`Q_1 = 1/2 + Q_0 + \\Delta Q` and
        :math:`\\omega_1 = 4\\pi Q_1 / (P\\sqrt{4Q_1^2-1})`, the fundamental's
        damped frequency :math:`\\omega_1\\sqrt{4Q_1^2-1}/2Q_1` is exactly
        :math:`2\\pi/P` — which is why ``period`` is documented as the *damped*
        period here and as the undamped one on :class:`SHO`.
        """
        amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)
        period = self._hyperparameter(values, "period")
        quality = self._hyperparameter(values, "quality")
        delta = self._hyperparameter(values, "delta_quality", allow_zero=True)
        fraction = self._hyperparameter(values, "fraction", allow_zero=True)
        power = amplitude * amplitude / (1.0 + fraction)
        quadruples: list[tuple[Any, Any, Any, Any]] = []
        for factor, harmonic, share in (
            (0.5 + quality + delta, 1.0, 1.0),
            (0.5 + quality, 2.0, fraction),
        ):
            split = (4.0 * factor * factor - 1.0) ** 0.5
            omega = harmonic * 2.0 * _TWO_PI * factor / (period * split)
            decay = 0.5 * omega / factor
            cosine = share * power
            quadruples.append((cosine, cosine / split, decay, decay * split))
        return tuple(quadruples)

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        ops = self.ops
        tau = ops.scalar(separation)
        total = None
        for cosine, sine, decay, frequency in self.coefficients(values):
            term = ops.exp(-decay * tau) * (
                cosine * ops.cos(frequency * tau) + sine * ops.sin(frequency * tau)
            )
            total = term if total is None else total + term
        return total


# ---------------------------------------------------------------------------
# The algebra
# ---------------------------------------------------------------------------


def _default_labels(count: int) -> tuple[str, ...]:
    return tuple(f"term{index}" for index in range(count))


class _Composite(Kernel):
    """Shared machinery for :class:`Sum` and :class:`Product`.

    **Parameter namespacing.** A composite's children are labelled, and their
    hyperparameters are declared under ``label.name`` — the convention
    ``ParameterSet.merge`` already uses for component parameters
    (``parameters.md``; ``"spectrum.temperature"``). ``Sum(Matern32(...),
    Matern32(...))`` therefore declares ``term0.amplitude``,
    ``term0.length_scale``, ``term1.amplitude``, ``term1.length_scale`` and two
    Matérn terms in one sum cannot collide.

    Positional indices — ``terms.0.amplitude`` — were the obvious alternative
    and are **not usable**: ``_check_name`` requires every dot-separated
    segment to be a Python identifier, because a parameter name is handed to a
    model as a keyword argument, and ``0`` is not one. ``term0`` is the same
    information under that rule. ``labels=`` replaces the defaults where a name
    carries meaning (``labels=("broad", "narrow")``), which is worth doing:
    these names reach the sampler's coordinates and the stored chains.

    A composite adopts its children's capability flags rather than declaring
    its own. They go into the instance ``__dict__`` — shadowing the class-level
    declarations — because ``declared_capabilities`` reads them off the
    instance, and a ``Sum`` of two torch kernels is a torch kernel however it
    was built.
    """

    #: What the composite does to its children's matrices, for messages.
    OPERATION: ClassVar[str] = ""

    def __init__(self, *terms: Kernel, labels: Sequence[str] | None = None) -> None:
        if len(terms) < 2:
            raise LikelihoodError(
                f"{type(self).__name__} composes two or more kernels, got {len(terms)}. A "
                f"one-term composite is the term itself; pass it directly."
            )
        for term in terms:
            if not isinstance(term, Kernel):
                raise LikelihoodError(
                    f"{type(self).__name__} composes Kernels, got "
                    f"{type(term).__name__}. Kernels are declared neutrally (family name plus "
                    f"Parameter hyperparameters) so that a lowering rule can translate them."
                )
        chosen = _default_labels(len(terms)) if labels is None else tuple(str(x) for x in labels)
        if len(chosen) != len(terms):
            raise LikelihoodError(
                f"{type(self).__name__} was given {len(terms)} kernel(s) and {len(chosen)} "
                f"label(s); there must be one label per term."
            )
        if len(set(chosen)) != len(chosen):
            raise LikelihoodError(
                f"{type(self).__name__}'s labels {chosen} repeat a name. A term's label "
                f"qualifies its hyperparameters, so two terms cannot share one."
            )
        for label in chosen:
            if not label.isidentifier():
                raise LikelihoodError(
                    f"{type(self).__name__}'s term label {label!r} is not a Python identifier. "
                    f"Labels qualify parameter names (label.amplitude), which are handed to "
                    f"models as keyword arguments."
                )
        self._terms: tuple[tuple[str, Kernel], ...] = tuple(zip(chosen, terms, strict=True))
        for label, term in self._terms:
            for parameter in term.parameters:
                self.register_parameter(parameter.rename(f"{label}.{parameter.name}"))
        self.__dict__["HYPERPARAMETERS"] = tuple(self.parameters.names)
        self._adopt(terms)

    def _adopt(self, terms: Sequence[Kernel]) -> None:
        """Take the children's capability flags, conjunctively, and their namespace."""
        devices = {term.DEVICE for term in terms}
        backends = {term.BACKEND for term in terms}
        if len(devices) > 1 or len(backends) > 1:
            raise LikelihoodError(
                f"{type(self).__name__} was given kernels from different places — devices "
                f"{sorted(devices)}, backends {sorted(backends)}. A composed covariance is one "
                f"array computation; ampere will not move arrays between devices or between "
                f"array libraries on your behalf (architecture.md §5)."
            )
        self.__dict__.update(
            {
                "DIFFERENTIABLE": all(term.DIFFERENTIABLE for term in terms),
                "BATCHABLE": all(term.BATCHABLE for term in terms),
                "DEVICE": devices.pop(),
                "BACKEND": backends.pop(),
                "VALIDATES": all(term.VALIDATES for term in terms),
                "_ops": terms[0].ops,
            }
        )

    @property
    def terms(self) -> tuple[tuple[str, Kernel], ...]:
        return self._terms

    def _bind(self, axis_names: tuple[str, ...]) -> Kernel:
        bound = copy.copy(self)
        bound.__dict__["_bound_cache"] = {}
        bound.__dict__["_columns"] = self._resolve_columns(axis_names)
        bound.__dict__["_terms"] = tuple(
            (label, child.for_axes(axis_names)) for label, child in self._terms
        )
        return bound

    def _child_values(self, label: str, resolved: Mapping[str, Any]) -> dict[str, Any]:
        """The child's own, unqualified hyperparameter mapping."""
        prefix = f"{label}."
        return {
            name[len(prefix) :]: value
            for name, value in resolved.items()
            if name.startswith(prefix)
        }

    def _child_matrices(self, left: Any, right: Any, values: Mapping[str, Any]) -> list[Any]:
        resolved = self.resolve(values)
        return [
            child.matrix(left, right, self._child_values(label, resolved))
            for label, child in self._terms
        ]

    def _child_diagonals(self, coordinates: Any, values: Mapping[str, Any]) -> list[Any]:
        resolved = self.resolve(values)
        return [
            child.diagonal(coordinates, self._child_values(label, resolved))
            for label, child in self._terms
        ]

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        raise LikelihoodError(
            f"{type(self).__name__} has no single closed form in one separation: its terms may "
            f"act on different axes, so each measures its own. Call matrix() or diagonal(), "
            f"which is what every solver does."
        )

    def _child_values_all(self, values: Mapping[str, Any]) -> list[dict[str, Any]]:
        resolved = self.resolve(values)
        return [self._child_values(label, resolved) for label, _ in self._terms]

    def value(self, separation: Any, values: Mapping[str, Any]) -> Any:
        """``k(τ)`` of the composite, from its children's — same-axis terms only.

        A :class:`Product` of kernels on *different* axes has no single ``τ``,
        so this is the one place the algebra presumes the same-axis case. Every
        caller that can be on different axes (``matrix``, ``diagonal``, and so
        every solver) measures each term's own separation instead.
        """
        parts = [
            child.value(separation, child_values)
            for (_, child), child_values in zip(
                self._terms, self._child_values_all(values), strict=True
            )
        ]
        total = parts[0]
        for part in parts[1:]:
            total = total + part if self.OPERATION == "sum" else total * part
        return total

    def __repr__(self) -> str:
        described = ", ".join(f"{label}={child!r}" for label, child in self._terms)
        selection = "" if self._axes is None else f", axes={self._axes!r}"
        return f"{type(self).__name__}({described}{selection})"


class Sum(_Composite):
    """A sum of kernels: ``k(x, x') = Σ_i k_i(x, x')``.

    The composition ``likelihoods.md`` §15.1 named as its own extension point
    ("a long length-scale plus a short one"), and the one that keeps the O(N)
    path: a sum of quasiseparable terms **is** quasiseparable, because the
    semiseparable generators simply concatenate and the ranks add
    (:func:`sum_representation`). :attr:`QUASISEPARABLE` is therefore an
    instance-level conjunction of the children's, not a class-level claim —
    ``Sum(Matern32(...), SquaredExponential(...))`` is honestly not
    quasiseparable, and ``QuasisepGP`` refuses it by naming the term that
    cannot be lowered.

    Order matters for the *declaration* even though addition is commutative:
    the labels qualify the children's parameter names, so ``Sum(a, b)`` and
    ``Sum(b, a)`` produce different parameter names, different sampler
    coordinates and different stored chains. :class:`KernelSpec` records the
    terms in declaration order for that reason.

    Examples
    --------
    >>> kernel = Sum(Matern32(0.3, 2.0), Matern12(0.1, 0.2))
    >>> kernel.parameters.names
    ('term0.amplitude', 'term0.length_scale', 'term1.amplitude', 'term1.length_scale')
    >>> kernel.QUASISEPARABLE
    True
    >>> Sum(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)).QUASISEPARABLE
    False
    >>> named = Sum(Matern32(0.3, 2.0), Matern32(0.1, 0.2), labels=("broad", "narrow"))
    >>> named.parameters.free_names
    ()
    >>> named.parameters.names[0]
    'broad.amplitude'
    """

    FAMILY: ClassVar[str] = "sum"
    OPERATION: ClassVar[str] = "sum"
    #: A class-level *capability*, narrowed per instance in ``__init__``: a sum
    #: can be quasiseparable, and this one is exactly when all its terms are.
    QUASISEPARABLE: ClassVar[bool] = True

    def __init__(self, *terms: Kernel, labels: Sequence[str] | None = None) -> None:
        super().__init__(*terms, labels=labels)
        self.__dict__["QUASISEPARABLE"] = all(term.QUASISEPARABLE for term in terms)

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> Any:
        matrices = self._child_matrices(left, right, values)
        total = matrices[0]
        for block in matrices[1:]:
            total = total + block
        return total

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> Any:
        diagonals = self._child_diagonals(coordinates, values)
        total = diagonals[0]
        for block in diagonals[1:]:
            total = total + block
        return total

    def select(self, points: Any) -> Any:
        """The one coordinate column a quasiseparable solve runs on.

        A sum on the O(N) path must agree with itself about *which* coordinate
        it is a function of — the recursion factorises one ordered axis — so
        the children's selections must coincide. ``check_axes`` refuses a sum
        whose terms select different axes before a solver ever gets here; this
        is the direct-call guard.
        """
        if self._columns is not None:
            return self.ops.take_columns(points, self._columns)
        selections = {child._columns for _, child in self._terms}
        if len(selections) != 1:
            raise LikelihoodError(
                "a Sum whose terms select different axes has no single coordinate to run a "
                "quasiseparable solve on. Give every term the same axes=, or use DenseGP, which "
                "evaluates each term on its own axes."
            )
        columns = selections.pop()
        return points if columns is None else self.ops.take_columns(points, columns)


class Product(_Composite):
    """A product of kernels: ``k(x, x') = Π_i k_i(x, x')``.

    The composition the chromatic interferometric case needs
    (``phase4_placement_memo.md`` §3.6): a residual that is a patch on the sky
    with a spectral profile is ``S(λ) · F(B/λ)``, smooth in spatial frequency
    and sharp in wavelength, so its covariance is
    ``k_uv(u, v) · k_λ(λ)`` — two kernels on **disjoint axis subsets** of one
    container, which is what the ``axes`` selector exists for. An isotropic
    kernel over all three axes cannot express it, and its own single-unit
    refusal says so.

    A product's marginal variance is the product of its terms', so two terms
    each declaring an ``amplitude`` over-parameterise it by one degree of
    freedom; fix all but one amplitude (pass a number rather than a prior),
    exactly as one would for any other multiplicative redundancy.

    **Never quasiseparable.** A product of semiseparable matrices is
    semiseparable only in special cases, and never when the factors act on
    different coordinates — which is the case a product is *for*. The
    declaration is a flat ``False`` rather than a conjunction, and the O(N)
    path refuses it by name.

    Examples
    --------
    >>> spatial = Matern32(0.3, 2.0, axes=("u", "v"))
    >>> spectral = Matern32(1.0, 0.05, axes=("spectral_axis",))
    >>> Product(spatial, spectral).QUASISEPARABLE
    False
    """

    FAMILY: ClassVar[str] = "product"
    OPERATION: ClassVar[str] = "product"
    QUASISEPARABLE: ClassVar[bool] = False

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> Any:
        matrices = self._child_matrices(left, right, values)
        total = matrices[0]
        for block in matrices[1:]:
            total = total * block
        return total

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> Any:
        diagonals = self._child_diagonals(coordinates, values)
        total = diagonals[0]
        for block in diagonals[1:]:
            total = total * block
        return total


class SpectralMixture(Sum):
    """A sum of damped oscillators with free frequencies (Wilson & Adams 2013).

    The horizon note's observation, cashed in: "spectral-mixture kernels are
    sums of damped oscillators — celerite terms already are, so a sum of SHO
    terms *is* a spectral mixture in this codebase". This class is that sum,
    with its own family name so provenance can tell it from a hand-built one,
    and with the component count as the only thing a user has to decide.

    Each component is an :class:`SHO`: a frequency (as a period), a coherence
    (as a quality factor) and a weight (as an amplitude). Free frequencies are
    the point — a mixture with enough components approximates any stationary
    covariance — and they are also the risk, which is why W4.5 stops here: the
    sparsity-inducing prior that would let the data switch components off is
    expressible today with ``HierarchicalPrior`` (``parameters.md`` §9) and is
    a study, not a contract.

    Parameters
    ----------
    amplitudes, periods, qualities
        One entry per component: a prior, a number, or a ``Parameter``.
    component_type
        The :class:`SHO` class to build components from. The default is
        ``ampere.core``'s; a backend passes its own, which is all the
        backend-awareness this class needs (its arithmetic is
        :class:`Sum`'s, which is the children's).
    labels
        Component labels; ``component0``, ``component1``, … by default.

    Examples
    --------
    >>> mixture = SpectralMixture([0.2, 0.1], [1.0, 0.3], [4.0, 8.0])
    >>> mixture.spec().family
    'spectral_mixture'
    >>> mixture.parameters.names[:3]
    ('component0.amplitude', 'component0.period', 'component0.quality')
    >>> mixture.QUASISEPARABLE
    True
    """

    FAMILY: ClassVar[str] = "spectral_mixture"

    def __init__(
        self,
        amplitudes: Sequence[Any],
        periods: Sequence[Any],
        qualities: Sequence[Any],
        *,
        component_type: type[SHO] = SHO,
        labels: Sequence[str] | None = None,
        axes: Sequence[str] | None = None,
        **component_kwargs: Any,
    ) -> None:
        counts = {len(amplitudes), len(periods), len(qualities)}
        if len(counts) != 1:
            raise LikelihoodError(
                f"SpectralMixture needs one amplitude, one period and one quality per "
                f"component, but was given {len(amplitudes)}, {len(periods)} and "
                f"{len(qualities)}."
            )
        components = [
            component_type(amplitude, period, quality, axes=axes, **component_kwargs)
            for amplitude, period, quality in zip(amplitudes, periods, qualities, strict=True)
        ]
        chosen = (
            tuple(f"component{index}" for index in range(len(components)))
            if labels is None
            else tuple(labels)
        )
        # The selection lives on the components, not on the sum: every
        # component is a function of the same coordinate, and a selection
        # declared twice would appear twice in the spec for no extra meaning.
        super().__init__(*components, labels=chosen)


# ---------------------------------------------------------------------------
# The semiseparable representations, and the public registry
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class CeleriteRepresentation:
    r"""The exact semiseparable representation of one kernel, on one coordinate set.

    The solver underneath every O(N) path — celerite2's, on all three backends
    — factorises any rank-J matrix of the form

    .. math::
        K_{nm} = \sum_j U_{nj} V_{mj} e^{-c_j (t_n - t_m)}, \qquad n > m,

    with :math:`K_{nn} = k(0) + \sigma_n^2`. This dataclass is that form, and
    it is deliberately **not** the celerite basis
    :math:`e^{-c\tau}(a\cos d\tau + b\sin d\tau)`: the basis has no
    :math:`\tau e^{-c\tau}` member, so Matérn-3/2 and -5/2 have no exact form
    in it, while both have one here. Declaring the representation at this level
    is what lets ``likelihoods.md`` §6 say "exactly quasiseparable" and mean
    an algebraic identity rather than a limit (W2.3, extended at W4.5).

    Attributes
    ----------
    decay
        ``c``, shape ``(J,)``.
    left
        ``U``, shape ``(n, J)``.
    right
        ``V``, shape ``(n, J)``.
    marginal
        ``k(0)``, a scalar; the solver adds it to the noise diagonal.

    Notes
    -----
    **Centre the coordinates.** Every builder below re-references ``axis`` to
    the midpoint of its own range before building the generators. The products
    only ever involve differences, so this is exact, and it is necessary: the
    Matérn generators grow linearly (and Matérn-5/2's quadratically) in the
    coordinate, so ``U_n · V_m`` is a difference of large numbers when the data
    span many length scales. Centring bounds the cancellation by half that
    span. Measured against the dense Cholesky for Matérn-3/2: ~5e-12 over 10
    length scales, ~3e-10 over 10², ~2e-8 over 10⁴.
    """

    decay: Any
    left: Any
    right: Any
    marginal: Any

    @property
    def rank(self) -> int:
        """``J``: how many exponentials the representation carries."""
        return int(np.shape(self.decay)[0])


#: What a quasiseparable-term builder is: ``(kernel, values, axis)`` to a
#: :class:`CeleriteRepresentation`. ``axis`` is the **sorted** bare coordinate
#: in the backend's own array type, and ``kernel.ops`` is the namespace to
#: build in — which is the whole reason one builder serves three backends.
TermBuilder = Callable[["Kernel", Mapping[str, Any], Any], CeleriteRepresentation]


def _centred(kernel: Kernel, axis: Any) -> Any:
    """``axis`` re-referenced to the midpoint of its own range. Exact; see above."""
    points = kernel.ops.scalar(axis)
    return points - 0.5 * (points[0] + points[-1])


def matern12_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    r"""Rank 1: :math:`k(\Delta) = a^2 e^{-\Delta/\ell}` is one real celerite term.

    :math:`U_n = (a^2)`, :math:`V_m = (1)`, :math:`c = (1/\ell)`. No
    cancellation to bound, so no centring is needed — it is applied anyway, for
    free, so every builder here reads the same way.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    amplitude = kernel._hyperparameter(resolved, "amplitude", allow_zero=True)
    length_scale = kernel._hyperparameter(resolved, "length_scale")
    shifted = _centred(kernel, axis)
    ones = ops.ones_like(shifted)
    marginal = amplitude * amplitude
    return CeleriteRepresentation(
        decay=ops.stack([1.0 / length_scale]),
        left=ops.stack([marginal * ones]),
        right=ops.stack([ones]),
        marginal=marginal,
    )


def matern32_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    r"""Rank 2, exact. W2.3's algebra, now written once for all three backends.

    With :math:`f = \sqrt3/\ell` and :math:`\Delta = t_n - t_m > 0`,

    .. math::
        k(\Delta) = a^2 (1 + f\Delta)e^{-f\Delta}
                  = e^{-f(t_n - t_m)}
                    \big[a^2(1 + f t_n)\cdot 1 + (-a^2 f)\cdot t_m\big],

    so :math:`U_n = (a^2(1 + f t_n),\, -a^2 f)`, :math:`V_m = (1,\, t_m)` and
    :math:`c = (f, f)`. The bracket is :math:`a^2(1 + f t_n - f t_m)`, which is
    an algebraic identity and not celerite2's ``Matern32Term`` ε-limit — the
    latter misses ``tolerances.cross_solver`` by three orders of magnitude at
    its default, which is why ampere carries its own term at all.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    amplitude = kernel._hyperparameter(resolved, "amplitude", allow_zero=True)
    length_scale = kernel._hyperparameter(resolved, "length_scale")
    decay = _SQRT3 / length_scale
    marginal = amplitude * amplitude
    shifted = _centred(kernel, axis)
    ones = ops.ones_like(shifted)
    return CeleriteRepresentation(
        decay=ops.stack([decay, decay]),
        left=ops.stack([marginal * (1.0 + decay * shifted), -marginal * decay * ones]),
        right=ops.stack([ones, shifted]),
        marginal=marginal,
    )


def matern52_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    r"""Rank 3, exact — the correction to ``likelihoods.md`` §15.3.

    With :math:`f = \sqrt5/\ell` and :math:`\Delta = t_n - t_m > 0`,
    :math:`k(\Delta) = a^2(1 + f\Delta + \tfrac13 f^2\Delta^2)e^{-f\Delta}`.
    Expanding the polynomial in :math:`t_n` and :math:`t_m`,

    .. math::
        1 + f(t_n - t_m) + \tfrac13 f^2 (t_n^2 - 2 t_n t_m + t_m^2)
        = \big[1 + f t_n + \tfrac13 f^2 t_n^2\big]\cdot 1
        + \big[-f - \tfrac23 f^2 t_n\big]\cdot t_m
        + \big[\tfrac13 f^2\big]\cdot t_m^2,

    a **rank-3 bilinear form** in :math:`(1, t, t^2)`. So
    :math:`U_n = a^2(1 + f t_n + \tfrac13 f^2 t_n^2,\ -f - \tfrac23 f^2 t_n,\
    \tfrac13 f^2)`, :math:`V_m = (1,\, t_m,\, t_m^2)`, :math:`c = (f, f, f)`.

    §15.3 recorded Matérn-5/2 as "not exactly quasiseparable" because it has no
    exact form in the celerite *basis*, which is true; the solver factorises
    the *semiseparable* form, which is wider, and in that form the identity
    above is algebra. The quadratic generators make the cancellation one order
    worse than Matérn-3/2's — the centring in :func:`_centred` is what keeps it
    inside ``tolerances.cross_solver`` — so this is the term to watch when data
    span very many length scales.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    amplitude = kernel._hyperparameter(resolved, "amplitude", allow_zero=True)
    length_scale = kernel._hyperparameter(resolved, "length_scale")
    decay = _SQRT5 / length_scale
    quadratic = decay * decay / 3.0
    marginal = amplitude * amplitude
    shifted = _centred(kernel, axis)
    ones = ops.ones_like(shifted)
    return CeleriteRepresentation(
        decay=ops.stack([decay, decay, decay]),
        left=ops.stack(
            [
                marginal * (1.0 + decay * shifted + quadratic * shifted * shifted),
                marginal * (-decay * ones - 2.0 * quadratic * shifted),
                marginal * quadratic * ones,
            ]
        ),
        right=ops.stack([ones, shifted, shifted * shifted]),
        marginal=marginal,
    )


def _celerite_pair(
    ops: ArrayOps,
    shifted: Any,
    cosine: Any,
    sine: Any,
    decay: Any,
    frequency: Any,
) -> tuple[Any, Any, Any]:
    r"""``(c, U, V)`` for one celerite term ``e^{-c\tau}(a\cos d\tau + b\sin d\tau)``.

    The standard construction, from the angle-difference identities:
    :math:`U_n = (a\cos d t_n + b\sin d t_n,\ a\sin d t_n - b\cos d t_n)`,
    :math:`V_m = (\cos d t_m,\ \sin d t_m)`, :math:`c = (c, c)`.
    """
    angle = frequency * shifted
    cos_t = ops.cos(angle)
    sin_t = ops.sin(angle)
    return (
        ops.stack([decay, decay]),
        ops.stack([cosine * cos_t + sine * sin_t, cosine * sin_t - sine * cos_t]),
        ops.stack([cos_t, sin_t]),
    )


def sho_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    """Rank 2: one celerite term, from :meth:`SHO.coefficients`.

    The coefficients come from the kernel itself rather than being recomputed
    here, so the dense closed form and the O(N) generators are two readings of
    one expression and the conformance row comparing them is a real test of the
    solver rather than of a copy-paste.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    cosine, sine, decay, frequency = kernel.coefficients(resolved)  # type: ignore[attr-defined]
    shifted = _centred(kernel, axis)
    decays, left, right = _celerite_pair(ops, shifted, cosine, sine, decay, frequency)
    return CeleriteRepresentation(decay=decays, left=left, right=right, marginal=cosine)


def rotation_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    """Rank 4: the two celerite terms of :meth:`RotationTerm.coefficients`, concatenated."""
    ops = kernel.ops
    resolved = kernel.resolve(values)
    quadruples = kernel.coefficients(resolved)  # type: ignore[attr-defined]
    shifted = _centred(kernel, axis)
    blocks = [
        _celerite_pair(ops, shifted, cosine, sine, decay, frequency)
        for cosine, sine, decay, frequency in quadruples
    ]
    marginal = quadruples[0][0]
    for quadruple in quadruples[1:]:
        marginal = marginal + quadruple[0]
    return CeleriteRepresentation(
        decay=ops.concatenate([block[0] for block in blocks], axis=0),
        left=ops.concatenate([block[1] for block in blocks], axis=-1),
        right=ops.concatenate([block[2] for block in blocks], axis=-1),
        marginal=marginal,
    )


def sum_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    """A sum of quasiseparable terms is quasiseparable: concatenate the generators.

    ``K = Σ_i K_i`` with ``K_i`` of rank ``J_i`` is of rank ``Σ_i J_i``, with
    ``c``, ``U`` and ``V`` stacked — because the sum of the bilinear forms *is*
    the bilinear form of the stacked generators. This is the one structural
    claim that makes W4.5's algebra worth having: an arbitrarily rich sum of
    terms still costs O(N), with the constant growing in the total rank.

    Each term is asked for its own representation through the registry, so a
    user-registered term composes with the built-in ones for free, and a term
    with no registered representation is refused **by name** here rather than
    silently approximated.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    blocks: list[CeleriteRepresentation] = []
    for label, child in kernel.terms:
        builder = lookup_quasiseparable_term(child.FAMILY, owner=type(kernel).__name__)
        child_values = kernel._child_values(label, resolved)  # type: ignore[attr-defined]
        blocks.append(builder(child, child_values, axis))
    marginal = blocks[0].marginal
    for block in blocks[1:]:
        marginal = marginal + block.marginal
    return CeleriteRepresentation(
        decay=ops.concatenate([block.decay for block in blocks], axis=0),
        left=ops.concatenate([block.left for block in blocks], axis=-1),
        right=ops.concatenate([block.right for block in blocks], axis=-1),
        marginal=marginal,
    )


@dataclasses.dataclass(frozen=True)
class QuasiseparableTerm:
    """One registered row: a kernel family and the builder of its exact representation.

    The shape of ``lowering.LoweringResolution``, and for the same reasons:
    a registration returns what was stored, a lookup returns what was found,
    and a row knows whether it is built in so ``term_provenance_entries`` can
    stamp the user rows into a run's record without stamping ampere's own.
    """

    family: str
    kernel_class: str
    builder: TermBuilder
    builtin: bool
    builder_module: str
    builder_qualname: str

    def to_provenance_entry(self) -> dict[str, Any]:
        """A JSON-plain, netCDF-safe record of this row for a run's provenance."""
        return {
            "kind": "quasiseparable_term",
            "family": self.family,
            "kernel": self.kernel_class,
            "builtin": self.builtin,
            "builder": f"{self.builder_module}.{self.builder_qualname}",
        }


#: Kernel family -> its exact semiseparable representation. **One slot per
#: family, not per (family, backend)**, which is the W4.5 decision the item
#: asked to be argued: a representation is generators built from the coordinate
#: and the hyperparameters with ``exp``, ``cos``, ``sin`` and ``stack``, and
#: that is the same mathematics in numpy, torch and jax. What differs between
#: backends is the *array namespace*, which :attr:`Kernel.ops` already carries,
#: so the backend is a property of the kernel instance handed to the builder
#: rather than of the slot. Three consequences, all wanted: a user registers
#: once and reaches the O(N) path on all three backends; the backends cannot
#: drift (before W4.5 each carried its own transcription, with a comment saying
#: the tables "must carry exactly the same families" and no way to check it);
#: and a backend that genuinely needs to specialise still can, by registering
#: its own builder for the family under ``override=True`` at import.
_TERMS: dict[str, QuasiseparableTerm] = {}


def register_quasiseparable_term(
    kernel_type: type[Kernel],
    builder: TermBuilder,
    *,
    override: bool = False,
    builtin: bool = False,
) -> QuasiseparableTerm:
    """Declare how a kernel family is represented on the O(N) path.

    The public route the horizon note asked for: before W4.5 a user could
    define a :class:`Kernel` subclass and use it on the dense path — the ABC is
    public and ``matrix``/``diagonal`` are all ``DenseGP`` needs — but could
    not make it fast, because the celerite translation table was private. This
    is that table, opened.

    Parameters
    ----------
    kernel_type
        The kernel class. Its ``FAMILY`` is the registry key, and it must
        declare ``QUASISEPARABLE = True``: a family with a registered
        representation that says it has none would be refused at composition
        anyway, and the contradiction is better caught here.
    builder
        ``(kernel, values, axis) -> CeleriteRepresentation``. ``values`` are
        the kernel's *resolved* hyperparameters, ``axis`` the **sorted** bare
        coordinate in the backend's array type, and ``kernel.ops`` the
        namespace to build in.
    override
        Replace an existing row deliberately. Without it, a second
        registration for one family raises, exactly as the lowering registry
        does — a silently replaced covariance is a changed model.
    builtin
        Mark this as one of ampere's own rows. Built-in rows still require
        ``override=True`` to replace, and are never counted as user-registered
        by :func:`term_provenance_entries`. Ordinary callers leave it ``False``.

    Returns
    -------
    QuasiseparableTerm
        The row just stored.

    Examples
    --------
    >>> class Relabelled(Matern12):
    ...     FAMILY = "relabelled_matern12"
    >>> row = register_quasiseparable_term(Relabelled, matern12_representation)
    >>> row.family, row.builtin
    ('relabelled_matern12', False)
    >>> register_quasiseparable_term(Relabelled, matern12_representation)
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: a user quasiseparable representation is already
    registered for kernel family 'relabelled_matern12' ...
    >>> _forget_quasiseparable_term("relabelled_matern12")
    True
    """
    if not isinstance(kernel_type, type) or not issubclass(kernel_type, Kernel):
        raise LikelihoodError(
            f"register_quasiseparable_term takes a Kernel subclass, got {kernel_type!r}."
        )
    family = kernel_type.FAMILY
    if not family:
        raise LikelihoodError(
            f"{kernel_type.__name__} declares no FAMILY, so there is nothing to key a "
            f"quasiseparable representation on. A kernel's FAMILY is its neutral name in every "
            f"registry and in its KernelSpec."
        )
    if not kernel_type.QUASISEPARABLE:
        raise LikelihoodError(
            f"{kernel_type.__name__} declares QUASISEPARABLE = False, so registering an exact "
            f"representation for it would contradict its own declaration. Set QUASISEPARABLE = "
            f"True if the representation is exact; if it is an approximation, it does not belong "
            f"on a solver whose name says EXACT."
        )
    if not callable(builder):
        raise LikelihoodError(f"a quasiseparable term builder must be callable, got {builder!r}.")
    existing = _TERMS.get(family)
    if existing is not None and not override:
        origin = "a built-in" if existing.builtin else "a user"
        raise LikelihoodError(
            f"{origin} quasiseparable representation is already registered for kernel family "
            f"{family!r} ({existing.builder_module}.{existing.builder_qualname}); pass "
            f"override=True to replace it deliberately. A covariance replaced silently is a "
            f"model changed silently."
        )
    row = QuasiseparableTerm(
        family=family,
        kernel_class=kernel_type.__name__,
        builder=builder,
        builtin=builtin,
        builder_module=getattr(builder, "__module__", "<unknown>"),
        builder_qualname=getattr(builder, "__qualname__", repr(builder)),
    )
    _TERMS[family] = row
    return row


def lookup_quasiseparable_term(family: str, *, owner: str = "QuasisepGP") -> TermBuilder:
    """The builder registered for *family*, or a refusal naming it.

    The **only** route onto the O(N) path, on every backend: a kernel that
    declares ``QUASISEPARABLE`` without a registered representation is refused
    by name rather than silently approximated (W2.3's ruling, widened at W4.5
    from one private table per backend to one public registry).
    """
    row = _TERMS.get(family)
    if row is None:
        known = ", ".join(sorted(_TERMS)) or "(none)"
        raise LikelihoodError(
            f"a kernel of family {family!r} declares QUASISEPARABLE = True, but ampere holds no "
            f"exact celerite representation for it, so {owner} has nothing to lower it to. "
            f"Families with one: {known}. Register one with "
            f"ampere.core.register_quasiseparable_term, or use DenseGP — a wrong representation "
            f"would be an approximation wearing an exact solver's name."
        )
    return row.builder


def _forget_quasiseparable_term(family: str) -> bool:
    """Drop a row. Private: for tests and doctests that register one out of tree.

    Deliberately not public. Registration is a *declaration*, and a public
    un-declaration would invite a "register, use, unregister" idiom in which
    two evaluations of one problem lower differently.
    """
    return _TERMS.pop(family, None) is not None


def quasiseparable_families() -> tuple[str, ...]:
    """Every kernel family that can reach the O(N) path, sorted."""
    return tuple(sorted(_TERMS))


def registered_quasiseparable_terms() -> tuple[QuasiseparableTerm, ...]:
    """Every registered row, built-in and user, in registration order."""
    return tuple(_TERMS.values())


def term_provenance_entries(
    rows: Iterable[QuasiseparableTerm] | None = None,
) -> list[dict[str, Any]]:
    """User-registered rows, as provenance entries. Built-in rows are omitted.

    ``lowering.provenance_entries``'s shape and its rule: ampere's own rows are
    not news, a user's are. Pass to
    ``ampere.results.provenance.provenance_attrs``'s ``extra=`` as
    ``{"registered_quasiseparable_terms": term_provenance_entries()}``.
    """
    chosen = registered_quasiseparable_terms() if rows is None else rows
    return [row.to_provenance_entry() for row in chosen if not row.builtin]


for _kernel_type, _builder in (
    (Matern12, matern12_representation),
    (Matern32, matern32_representation),
    (Matern52, matern52_representation),
    (SHO, sho_representation),
    (RotationTerm, rotation_representation),
    (Sum, sum_representation),
    (SpectralMixture, sum_representation),
):
    register_quasiseparable_term(_kernel_type, _builder, builtin=True)
del _kernel_type, _builder
