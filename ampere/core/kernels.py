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
import itertools
import math
from collections.abc import Callable, Iterable, Mapping, Sequence
from typing import Any, ClassVar, Protocol

import astropy.units as u
import numpy as np
import scipy.stats as st

from .exceptions import LikelihoodError
from .parameter import HierarchicalPrior, Identity, Log, Parameter, Parameterised

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
    "WarpedKernel",
    "find_warped_composite",
    "lookup_quasiseparable_term",
    "matern12_representation",
    "matern32_representation",
    "matern52_representation",
    "quantile_knots",
    "quasiseparable_families",
    "refuse_warped_composite",
    "register_quasiseparable_term",
    "registered_quasiseparable_terms",
    "rotation_representation",
    "sho_representation",
    "sum_representation",
    "term_provenance_entries",
    "warped_representation",
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

    def log1p(self, array: Any) -> Any:
        """Elementwise ``log(1 + x)``, accurate for small ``x``."""
        ...

    def absolute(self, array: Any) -> Any:
        """Elementwise absolute value."""
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

    def log1p(self, array: Any) -> np.ndarray:
        return np.log1p(array)

    def absolute(self, array: Any) -> np.ndarray:
        return np.abs(array)


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

    **W5.7 adds a third**, on the same terms — omitted from :meth:`to_dict`
    when unused, so every spec hash minted before W5.7 is unchanged:

    ``metadata``
        ``(key, value)`` pairs of JSON-plain, **non-parameter** structure a
        family needs in its declaration. :class:`WarpedKernel` is what asked
        for it: a warp's knot *locations* are not hyperparameters (they are
        fixed positions, not things a sampler moves) and they are not
        hyperparameter *names* either, yet two warps with different knots are
        different models and must hash differently. Values are tuples rather
        than lists so a spec stays hashable; :meth:`to_dict` converts.

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
    metadata: tuple[tuple[str, Any], ...] = ()

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
        if self.metadata:
            described["metadata"] = {
                key: list(value) if isinstance(value, tuple) else value
                for key, value in self.metadata
            }
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

    #: Whether :meth:`warped_coordinate` is anything but the identity (W5.7).
    #: A kernel that warps declares it, so a solver and a composite can *ask*
    #: rather than compare coordinate arrays at solve time.
    WARPS: ClassVar[bool] = False

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

    @property
    def warps(self) -> bool:
        """Whether anything in this tree warps the coordinate (W5.7)."""
        return self.WARPS or any(child.warps for _, child in self.terms)

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
        # Plain fancy indexing rather than a namespace method: numpy, torch and
        # jax all read ``x[:, [0, 1]]`` the same way, so the selection is one
        # of the few array operations that needs no ArrayOps entry at all --
        # and, more to the point, it then works whatever namespace the *kernel*
        # carries, which need not be the namespace of the array the solver
        # happens to hand it.
        return points[:, list(columns)]

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
                        f"single isotropic length-scale is meaningless across mixed units; name "
                        f"the axes this kernel acts on with axes=(...), as in "
                        f"{type(leaf).__name__}(axes={axis_names[:1]!r}), and compose kernels on "
                        f"different axes with Product."
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

    def _child_values(self, label: str, resolved: Mapping[str, Any]) -> dict[str, Any]:
        """A labelled child's own, unqualified hyperparameter mapping.

        On :class:`Kernel` rather than on :class:`_Composite` since W5.7, when
        a second kind of wrapping kernel — :class:`WarpedKernel` — needed the
        same unqualification, and :func:`sum_representation` stopped needing a
        ``type: ignore`` to reach it.
        """
        prefix = f"{label}."
        return {
            name[len(prefix) :]: value
            for name, value in resolved.items()
            if name.startswith(prefix)
        }

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

    # -- the warp hooks (W5.7) ----------------------------------------------

    def warped_coordinate(self, axis: Any, values: Mapping[str, Any]) -> Any:
        """The coordinate this kernel is **stationary in**. The identity by default.

        **W5.7's one extension to the term-registry contract**, and the
        smallest that makes :class:`WarpedKernel` reach the O(N) path. A
        quasiseparable solve factorises a matrix whose off-diagonal decay is
        :math:`e^{-c_j (t_n - t_m)}`, and *the solver*, not the generator
        builder, supplies that :math:`t`: celerite2 forms the propagators from
        the coordinate handed to ``compute``/``factor``, so a kernel that is
        stationary in :math:`w(x)` rather than in :math:`x` cannot express
        itself through the generators alone. Folding the difference into
        ``U`` and ``V`` would mean multiplying by :math:`e^{\\pm c t}`, which
        is exactly the overflow celerite's factored form exists to avoid.

        So a solver asks the kernel, once per solve, for the coordinate its
        recursion runs on, and the registered builder still receives the
        **raw** sorted axis its docstring promises. Every warp must be
        **monotone increasing**, so the sorting permutation is unchanged and
        ``QuasisepGP``'s ordering precondition survives; :class:`WarpedKernel`
        is monotone by construction.

        The same map is what ``Likelihood.conditional`` reports in, so the
        whiteness and localisation diagnostics see the coordinate in which the
        residuals are supposed to be stationary.

        Parameters
        ----------
        axis
            The bare ``(n,)`` coordinate, in this kernel's own namespace.
        values
            The hyperparameter mapping, as for :meth:`matrix`.
        """
        return axis

    def warp_provenance(self, values: Mapping[str, Any]) -> dict[str, Any] | None:
        """A JSON-plain record of the warp this kernel applies, or ``None``.

        ``None`` — the default, and every kernel but :class:`WarpedKernel` —
        means "I am stationary in the coordinate you gave me", which is what
        lets ``Likelihood.conditional`` leave its output untouched for every
        declaration written before W5.7. A warping kernel returns enough to
        reconstruct the map: the knot locations, the realised knot values and
        the warped knot images.
        """
        return None

    def spectral_density(
        self, frequency: Any, values: Mapping[str, Any], *, dimensions: int = 1
    ) -> Any:
        r"""The kernel's power spectral density :math:`S(\omega)`, in ``dimensions`` axes.

        **W5.4's addition**, and the one thing a reduced-rank spectral solver
        needs that the dense path never did. The convention is the plain
        Fourier transform of the covariance over :math:`\mathbb{R}^d`,

        .. math::
            S(\boldsymbol\omega) = \int k(\|\mathbf r\|)\,
            e^{-i \boldsymbol\omega \cdot \mathbf r}\, \mathrm d^d \mathbf r,

        so that :math:`k(\mathbf r) = (2\pi)^{-d} \int S e^{i \boldsymbol\omega
        \cdot \mathbf r}\,\mathrm d^d\boldsymbol\omega` and
        :math:`S \ge 0` by Bochner's theorem. For an isotropic kernel it
        depends on :math:`\|\boldsymbol\omega\|` alone, which is what
        ``frequency`` carries, and on the number of axes the kernel acts on,
        which is what ``dimensions`` carries — the same covariance has a
        *different* spectral density in one axis and in two, and a solver that
        forgot the dimension would build a basis for the wrong process.

        Written in :attr:`ops`, so a backend's kernel inherits a
        differentiable one from the same source the closed form comes from:
        :class:`~ampere.core.hsgp.HilbertSpaceGP` is written once and the
        backends contribute only the linear algebra.

        The default **refuses by name**, the same declared-slot discipline an
        unimplemented solver follows: a family with no closed-form spectral
        density (a :class:`Product`, a :class:`RotationTerm`, or a user's own)
        must say so rather than have a reduced-rank solver guess. A
        :class:`SpectralMixture` needs nothing of its own -- it *is* a
        :class:`Sum` of :class:`SHO` terms in this codebase, so it inherits
        the sum's, which is the same structural fact the quasiseparable
        registry already exploits.

        Parameters
        ----------
        frequency
            Non-negative angular frequencies :math:`\|\boldsymbol\omega\|`, in
            radians per unit of the coordinate axis. Any shape.
        values
            The hyperparameter mapping, as for :meth:`matrix`.
        dimensions
            How many axes the kernel acts on — ``len(kernel.selected_axes(...))``.
        """
        raise LikelihoodError(
            f"{type(self).__name__} ({self.FAMILY}) has no closed-form spectral density, so a "
            f"reduced-rank spectral solver (HilbertSpaceGP) cannot represent it. The families "
            f"that do are Matern12, Matern32, Matern52, SquaredExponential, SHO and any Sum of "
            f"them — a SpectralMixture included, being a Sum of SHOs; use DenseGP, which needs "
            f"no spectral density at all."
        )

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

    **Since W5.4** a subclass may also declare :attr:`SPECTRAL_NU`, the Matérn
    smoothness ``nu``, and inherit :meth:`spectral_density` from the one
    closed form below rather than transcribing it — which is what puts
    Matérn-1/2, -3/2 and -5/2 on ``HilbertSpaceGP`` in three lines.
    """

    #: The Matérn smoothness this family is, or ``None`` for a stationary
    #: kernel that is not a Matérn (:class:`SquaredExponential` overrides
    #: :meth:`spectral_density` instead). Read only by
    #: :meth:`spectral_density`.
    SPECTRAL_NU: ClassVar[float | None] = None

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

    def spectral_density(
        self, frequency: Any, values: Mapping[str, Any], *, dimensions: int = 1
    ) -> Any:
        r"""The Matérn spectral density, for a family that declares :attr:`SPECTRAL_NU`.

        With ampere's normalisation (:math:`k(0) = a^2`, length scale
        :math:`\ell`), the isotropic Matérn-:math:`\nu` transform in ``d``
        dimensions is

        .. math::
            S(\omega) = a^2\,
            \frac{2^d \pi^{d/2}\,\Gamma(\nu + d/2)\,(2\nu)^\nu}
                 {\Gamma(\nu)\,\ell^{2\nu}}
            \left(\frac{2\nu}{\ell^2} + \omega^2\right)^{-(\nu + d/2)} ,

        which for :math:`\nu = 1/2`, :math:`d = 1` is the familiar
        :math:`2 a^2 \ell / (1 + \ell^2\omega^2)`. The :math:`\Gamma`\ s and the
        powers of two are ordinary Python floats — they depend on ``nu`` and
        ``d``, never on a fitted value — so only the amplitude and the length
        scale carry a gradient, and they do so in this kernel's own namespace.
        """
        nu = self.SPECTRAL_NU
        if nu is None:
            return super().spectral_density(frequency, values, dimensions=dimensions)
        resolved = self.resolve(values)
        amplitude = self._hyperparameter(resolved, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(resolved, "length_scale")
        omega = self.ops.scalar(frequency)
        half = 0.5 * float(dimensions)
        constant = (
            2.0 ** float(dimensions)
            * math.pi**half
            * math.gamma(nu + half)
            * (2.0 * nu) ** nu
            / math.gamma(nu)
        )
        pole = 2.0 * nu / (length_scale * length_scale) + omega * omega
        return (
            amplitude * amplitude * constant * length_scale ** (-2.0 * nu) * pole ** (-(nu + half))
        )


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
    #: W5.4: the smoothness ``StationaryKernel.spectral_density`` reads.
    SPECTRAL_NU: ClassVar[float | None] = 0.5

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
               quasiseparable=True, axes=None, terms=(), metadata=())
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
    #: W5.4: the smoothness ``StationaryKernel.spectral_density`` reads.
    SPECTRAL_NU: ClassVar[float | None] = 1.5

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
    #: W5.4: the smoothness ``StationaryKernel.spectral_density`` reads.
    SPECTRAL_NU: ClassVar[float | None] = 2.5

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

    def spectral_density(
        self, frequency: Any, values: Mapping[str, Any], *, dimensions: int = 1
    ) -> Any:
        r"""A Gaussian in frequency: :math:`S(\omega) = a^2 (2\pi)^{d/2} \ell^d
        e^{-\ell^2\omega^2/2}`.

        The one family whose reduced-rank approximation converges
        *exponentially* in the basis size, because its spectral density does —
        which is the compensation for its having no quasiseparable form at
        all. ``HilbertSpaceGP`` is therefore the scaling answer for a squared
        exponential in the way ``QuasisepGP`` is for a Matérn.
        """
        resolved = self.resolve(values)
        amplitude = self._hyperparameter(resolved, "amplitude", allow_zero=True)
        length_scale = self._hyperparameter(resolved, "length_scale")
        scaled = length_scale * self.ops.scalar(frequency)
        constant = (2.0 * math.pi) ** (0.5 * float(dimensions))
        return (
            amplitude
            * amplitude
            * constant
            * length_scale ** float(dimensions)
            * self.ops.exp(-0.5 * scaled * scaled)
        )


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

    def spectral_density(
        self, frequency: Any, values: Mapping[str, Any], *, dimensions: int = 1
    ) -> Any:
        r"""The celerite term's own transform, from the same :meth:`coefficients`.

        For the celerite form :math:`k(\tau) = e^{-c|\tau|}(a\cos d\tau +
        b\sin d|\tau|)` the one-sided integrals are elementary and give

        .. math::
            S(\omega) = \frac{a c + b\,(d - \omega)}{c^2 + (d-\omega)^2}
                      + \frac{a c + b\,(d + \omega)}{c^2 + (d+\omega)^2},

        the two Lorentzians a damped oscillator's power sits in, one at each
        signed resonance. It is non-negative for every ``Q > 1/2``, which is
        what makes the term a covariance in the first place, and it is written
        from :meth:`coefficients` rather than from the closed form so the
        dense, quasiseparable and reduced-rank paths cannot drift apart.

        Refused in more than one axis: a damped oscillator is a function of an
        ordered coordinate, and there is no isotropic form of it on a plane.
        """
        if int(dimensions) != 1:
            raise LikelihoodError(
                f"SHO has a spectral density in one axis only, but it was asked for one in "
                f"{int(dimensions)}. A damped oscillator is a process in an ordered coordinate "
                f"(time, wavelength); on a two-axis container select the single axis it runs "
                f"along with axes=(...), or use a Matérn, which is isotropic in any dimension."
            )
        cosine, sine, decay, rate = self.coefficients(self.resolve(values))
        omega = self.ops.scalar(frequency)
        lower = rate - omega
        upper = rate + omega
        return (decay * cosine + sine * lower) / (decay * decay + lower * lower) + (
            decay * cosine + sine * upper
        ) / (decay * decay + upper * upper)


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


def _adopted_capabilities(owner: str, terms: Sequence[Kernel]) -> dict[str, Any]:
    """The capability flags and namespace a wrapping kernel takes from its children.

    Shared by :class:`_Composite` and :class:`WarpedKernel` (W5.7), which face
    the same question for the same reason: neither has arithmetic of its own —
    each builds on what its children compute — so a ``Sum`` of two torch
    kernels, and a ``WarpedKernel`` over one, *are* torch kernels however they
    were constructed. The flags go into the instance ``__dict__``, shadowing
    the class-level declarations, because ``declared_capabilities`` reads them
    off the instance.
    """
    devices = {term.DEVICE for term in terms}
    backends = {term.BACKEND for term in terms}
    if len(devices) > 1 or len(backends) > 1:
        raise LikelihoodError(
            f"{owner} was given kernels from different places — devices "
            f"{sorted(devices)}, backends {sorted(backends)}. A composed covariance is one "
            f"array computation; ampere will not move arrays between devices or between "
            f"array libraries on your behalf (architecture.md §5)."
        )
    return {
        "DIFFERENTIABLE": all(term.DIFFERENTIABLE for term in terms),
        "BATCHABLE": all(term.BATCHABLE for term in terms),
        "DEVICE": devices.pop(),
        "BACKEND": backends.pop(),
        "VALIDATES": all(term.VALIDATES for term in terms),
        "_ops": terms[0].ops,
    }


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
        self.__dict__.update(_adopted_capabilities(type(self).__name__, terms))

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

    def warped_coordinate(self, axis: Any, values: Mapping[str, Any]) -> Any:
        """The identity, or a refusal when a term warps (W5.7).

        A quasiseparable recursion runs on **one** coordinate, and a composite
        whose terms warp differently has none: each term's generators would be
        built at its own ``w_i(x)`` while the propagators came from a single
        axis, giving a matrix that is not the composite's covariance and not
        obviously wrong either. The direct-call guard, in the shape of
        :meth:`Sum.select`'s; ``QuasisepGP.check_compatible`` refuses the same
        declaration at composition, which is where a user should meet it.
        """
        if not self.warps:
            return axis
        raise LikelihoodError(
            f"a {type(self).__name__} with a WarpedKernel among its terms has no single "
            f"coordinate for a quasiseparable recursion to run on: each warped term's "
            f"generators are built at its own w(x), while the recursion's propagators come "
            f"from one axis. Warp the composite instead of its terms — "
            f"WarpedKernel({type(self).__name__}(...), input_warp=...) is one warp of the "
            f"coordinate with several kernels on it, which is what the O(N) path can "
            f"represent — or use DenseGP, which evaluates each term on its own coordinate."
        )

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

    def spectral_density(
        self, frequency: Any, values: Mapping[str, Any], *, dimensions: int = 1
    ) -> Any:
        """The spectral density of a sum is the sum of its terms' (W5.4).

        Linearity of the Fourier transform, and the reason ``Sum`` is the
        composite a reduced-rank spectral solver can take while ``Product``
        is not: a product's transform is a *convolution* of the factors',
        which has no closed form in general and none at all when the factors
        act on different axes. One term with no closed-form density refuses
        the whole sum, by its own name, exactly as the quasiseparable path
        refuses an unregistered term.
        """
        parts = [
            child.spectral_density(frequency, child_values, dimensions=dimensions)
            for (_, child), child_values in zip(
                self._terms, self._child_values_all(values), strict=True
            )
        ]
        total = parts[0]
        for part in parts[1:]:
            total = total + part
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
            return points[:, list(self._columns)]
        selections = {child._columns for _, child in self._terms}
        if len(selections) != 1:
            raise LikelihoodError(
                "a Sum whose terms select different axes has no single coordinate to run a "
                "quasiseparable solve on. Give every term the same axes=, or use DenseGP, which "
                "evaluates each term on its own axes."
            )
        columns = selections.pop()
        return points if columns is None else points[:, list(columns)]


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
# Non-stationarity: the warps (W5.7)
# ---------------------------------------------------------------------------


def _relu(ops: ArrayOps, value: Any) -> Any:
    """``max(value, 0)``, written so that it is exactly ``0.0`` at ``value <= 0``.

    Branch-free and in :class:`ArrayOps`, so one expression serves numpy, torch
    and jax and is differentiable on the latter two. ``0.5 * (v + |v|)`` is the
    standard form; it is exact in IEEE-754 arithmetic (the sum is either
    ``2v`` or ``0``, and halving is exact).
    """
    return 0.5 * (value + ops.absolute(value))


def _softplus(ops: ArrayOps, value: Any) -> Any:
    r""":math:`\log(1 + e^v)`, in the overflow-safe form.

    ``max(v, 0) + log1p(exp(-|v|))``: the naive ``log(1 + exp(v))`` overflows
    for ``v`` above about 709 and returns ``inf`` where the answer is ``v``,
    which a sampler reaches by proposing a wide warp increment. The identity is
    exact rather than approximate.
    """
    magnitude = ops.absolute(value)
    return 0.5 * (value + magnitude) + ops.log1p(ops.exp(-magnitude))


def quantile_knots(coordinates: Any, count: int) -> tuple[float, ...]:
    """*count* knot locations at evenly spaced quantiles of *coordinates*.

    The convenience half of W5.7's knot-placement decision. The knots
    themselves are **fixed numbers in the declaration** — a kernel is a
    declaration, and one whose knots were computed from the data it is about to
    be fitted to would have a different spec hash for every dataset, and would
    silently change model when a sample was masked. So ampere does not place
    knots for you; it gives you this, which you call **once**, on your
    coordinate, and paste (or pass) the result into the declaration, where it
    is recorded.

    Quantiles rather than an even grid because the warp's resolution should
    follow the data's: an even grid over a coordinate with a sparse tail spends
    most of its degrees of freedom where there is nothing to fit.

    Parameters
    ----------
    coordinates
        Any 1-D coordinate set. Not stored; only its quantiles are used.
    count
        How many knots, at least two. The first and last are the coordinate's
        own minimum and maximum, so the warp interpolates rather than
        extrapolating over the data.

    Examples
    --------
    >>> quantile_knots(np.arange(11.0), 3)
    (0.0, 5.0, 10.0)
    """
    if count < 2:
        raise LikelihoodError(
            f"a warp needs at least two knots to have a slope at all, got {count}."
        )
    axis = np.asarray(coordinates, dtype=DTYPE).ravel()
    if axis.size == 0:
        raise LikelihoodError("quantile_knots was given an empty coordinate set.")
    quantiles = np.linspace(0.0, 1.0, count)
    return tuple(float(value) for value in np.quantile(axis, quantiles))


def _warp_knots(given: Any, what: str) -> tuple[float, ...]:
    """Validate and freeze one warp's knot locations. ``None`` means "no warp"."""
    if given is None:
        return ()
    knots = tuple(float(value) for value in np.asarray(given, dtype=DTYPE).ravel())
    if len(knots) < 2:
        raise LikelihoodError(
            f"WarpedKernel's {what}= needs at least two knot locations to have a slope at "
            f"all, got {len(knots)}. Use ampere.core.quantile_knots(coordinate, count) to "
            f"place them once, on your own coordinate, and record the result."
        )
    if not all(math.isfinite(knot) for knot in knots):
        raise LikelihoodError(f"WarpedKernel's {what}= knot locations must all be finite.")
    for lower, upper in itertools.pairwise(knots):
        if upper <= lower:
            raise LikelihoodError(
                f"WarpedKernel's {what}= knot locations {knots} are not strictly increasing "
                f"({upper} follows {lower}). The knots are positions on the coordinate axis and "
                f"the warp interpolates between consecutive ones; an out-of-order or repeated "
                f"knot has no monotone interpolant, so it is refused here rather than producing "
                f"a covariance that is silently not a covariance."
            )
    return knots


def _warp_hyperparameter(name: str, given: Any, *, positive: bool) -> Parameter:
    """:func:`_as_hyperparameter`, widened to accept a :class:`HierarchicalPrior`.

    The shrinkage declaration needs it: ``increment_k ~ Normal(0, s)`` with
    ``s`` an ampere parameter *is* a hierarchical prior, and those carry
    references by name rather than a frozen ``scipy`` object, so they have no
    ``ppf`` for :func:`_as_hyperparameter` to recognise.
    """
    if isinstance(given, HierarchicalPrior):
        return Parameter(name, given, bijection=Log() if positive else Identity())
    return _as_hyperparameter(name, given, None, positive=positive)


def _warp_values(given: Any, count: int, default: Any, what: str) -> list[Any]:
    """One declaration per knot variable: broadcast a single one, or check a sequence."""
    if given is None:
        return [default] * count
    if isinstance(given, (list, tuple, np.ndarray)):
        chosen = list(given)
        if len(chosen) != count:
            raise LikelihoodError(
                f"WarpedKernel's {what}= was given {len(chosen)} declaration(s) for {count} "
                f"knot variable(s). Pass one per variable, or a single prior or number to be "
                f"used for all of them."
            )
        return chosen
    if isinstance(given, Parameter):
        raise LikelihoodError(
            f"WarpedKernel's {what}= cannot broadcast a single Parameter across {count} knot "
            f"variables: a Parameter carries its own name, and each knot variable needs its "
            f"own. Pass a prior (shared by declaration, not by identity), a number, or a "
            f"sequence of {count} Parameters."
        )
    return [given] * count


def _selected_copy(kernel: Kernel, axes: tuple[str, ...] | None) -> Kernel:
    """A copy of *kernel*'s tree with *axes* recorded on every node.

    The warp acts on one coordinate, so every kernel underneath it acts on that
    same coordinate. Recording the selection on the children rather than only
    on the wrapper is what makes :meth:`Kernel.check_units` compare a base
    kernel's ``length_scale`` against the unit of the axis the warp actually
    selects, instead of against the container's first axis.
    """
    clone = copy.copy(kernel)
    clone.__dict__["_axes"] = axes
    clone.__dict__["_bound_cache"] = {}
    if kernel.terms:
        clone.__dict__["_terms"] = tuple(
            (label, _selected_copy(child, axes)) for label, child in kernel.terms
        )
    return clone


class WarpedKernel(Kernel):
    r"""A kernel warped in its input, its amplitude, or both — still O(N).

    The plan's Phase-5 "non-stationary flexible likelihood" bullet, landed at
    **W5.7**. Two wrappers around any base kernel, each of which keeps the
    exact quasiseparable solve:

    **Input warping.** A monotone map :math:`w` of the coordinate, so that
    :math:`k_{\mathrm{warped}}(x, x') = k(w(x), w(x'))`. A Matérn-3/2 whose
    length scale should be short in one band and long in another is *this*,
    not a new family: :math:`w` compresses the coordinate where the process
    varies fast and stretches it where it varies slowly. Because :math:`w` is
    monotone **by construction**, the sorting permutation of the warped axis is
    the sorting permutation of the raw one, so ``QuasisepGP``'s ordering
    precondition survives untouched, and the generators of the base kernel
    evaluated at :math:`w(x)` are the generators of the warped one.

    **Amplitude warping.** :math:`D K D` with :math:`D = \mathrm{diag}(a(x))`:
    a smooth, positive, per-coordinate scaling of the marginal standard
    deviation. A residual whose *size* varies across the band — a noisier
    region, a badly calibrated order — is this. A diagonal congruence of a
    rank-:math:`J` semiseparable matrix is rank-:math:`J` semiseparable, with
    :math:`U \to \mathrm{diag}(a) U` and :math:`V \to \mathrm{diag}(a) V`, so
    this one is free.

    Both together give ``a(x) k(w(x), w(x')) a(x')``, which is the general
    non-stationary form this item buys — and it is **not** a new solver, a new
    approximation or a new representation: it is the registry's own extension
    point (W4.5) used as intended, plus one hook
    (:meth:`Kernel.warped_coordinate`) for the coordinate the recursion runs
    on.

    Parameterisation
    ----------------
    **The input warp**, over ``K`` knot locations :math:`x_0 < \dots <
    x_{K-1}`, is piecewise linear with segment slopes

    .. math::
        m_k = \frac{\zeta(u_k)}{\zeta(0)}, \qquad
        \zeta(u) = \log(1 + e^u),

    anchored at :math:`w(x_0) = x_0` and extended linearly beyond the end
    knots with the end segments' slopes. Every slope is strictly positive
    whatever :math:`u_k` is, so **monotonicity is structural rather than
    checked**; there is no constraint for a sampler to violate and no rejection
    region in the posterior. Normalising by :math:`\zeta(0) = \log 2` puts the
    identity warp at :math:`u_k = 0` — and puts it there *exactly*: the ratio
    is computed as one expression evaluated at the proposal and at zero, so at
    :math:`u_k = 0` it is bit-for-bit ``1.0``, the offsets are bit-for-bit
    ``0.0``, and ``w(x)`` is ``x`` with no rounding at all. The identity warp
    is therefore not merely close to the base kernel, it **is** the base
    kernel, which is what makes the shrinkage prior below a prior on "how far
    from the base model did the data make me go".

    **The amplitude warp**, over ``J`` knot locations, is :math:`\log a`
    piecewise linear through per-knot levels :math:`\ell_j`, likewise extended
    linearly. :math:`\ell_j = 0` gives :math:`a \equiv 1` exactly, the same
    bit-identity at the same identity point.

    The degrees-of-freedom guard
    ----------------------------
    Warps are flexible, and flexibility that is free is flexibility that will
    be used to absorb signal. The guard is in the declaration, not in advice:

    * **Few knots.** ``K`` and ``J`` are yours to choose and should be small —
      three to six. The knots are not sampled, so each one costs exactly one
      dimension.
    * **Shrinkage to the identity.** Every knot variable's default prior is
      hierarchical, with a single scale shared across a warp's knots and its
      own half-normal prior: ``u_k ~ Normal(0, s)``, ``s ~ HalfNormal(0.5)``.
      The identity warp is the point of maximum prior density, so the data must
      *pay* to leave it — the horizon note's answer to §1's
      degrees-of-freedom risk, and the same shape the regularised horseshoe on
      summed noise components takes.
    * **Non-centred by default.** ``non_centred=True`` declares
      ``u_k = s · z_k`` with ``z_k ~ Normal(0, 1)``, which is the same prior
      and a far better posterior geometry: the centred form's funnel is
      exactly the pathology NUTS reports as divergences. ``non_centred=False``
      declares the centred form directly, with
      :class:`~ampere.core.parameter.HierarchicalPrior` — offered because it
      is what the plan names, and because a strongly identified warp samples
      fine either way — its default knot declaration *is* that hierarchical
      prior. In that form the scale reaches the knots only through the
      reference, so overriding the knots with an ordinary prior while the
      scale is still free is refused: a scale nothing depends on is a sampled
      dimension with no posterior.
    * **Warped diagnostics.** ``Likelihood.conditional`` reports in
      :math:`w(x)` and records the warp, so the whiteness (family B) and
      localisation (family C) diagnostics are run in the coordinate the
      residuals are supposed to be stationary in. A warp that has been used to
      hide structure shows up there.

    Two redundancies are worth naming, both the same kind as
    :class:`Product`'s amplitude redundancy and both handled by the shrinkage
    prior rather than by a constraint: a warp whose slopes are all equal is a
    rescaling of the coordinate, degenerate with the base ``length_scale``;
    and a constant ``log a`` is degenerate with the base ``amplitude``.

    Parameters
    ----------
    base
        The kernel to warp. Any :class:`Kernel` — including a :class:`Sum`, and
        including another :class:`WarpedKernel` — provided it declares no
        ``axes`` of its own (declare the selection here instead; the warp acts
        on one coordinate and so does everything under it).
    input_warp
        Knot **locations** for the input warp, in the coordinate's own units,
        strictly increasing; ``None`` for no input warp. See
        :func:`quantile_knots`.
    amplitude_warp
        Knot locations for the amplitude warp; ``None`` for no amplitude warp.
    increments
        The declaration of each of the ``K-1`` input-warp knot variables: a
        prior, a number (held fixed), a
        :class:`~ampere.core.parameter.HierarchicalPrior`, or a sequence of
        ``K-1`` of those. The default is the shrinkage declaration above.
    levels
        The same for the ``J`` amplitude-warp levels.
    input_scale, amplitude_scale
        The shrinkage scale of each warp: a prior, a number, or a
        :class:`~ampere.core.parameter.Parameter`. Defaults to a half-normal.
    non_centred
        Whether the knot variables are declared as standard normals scaled by
        the shrinkage scale (the default) or directly under a
        :class:`~ampere.core.parameter.HierarchicalPrior`.
    axes
        The single container axis this kernel warps, as for every other kernel.

    Examples
    --------
    >>> import scipy.stats as st
    >>> base = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))
    >>> kernel = WarpedKernel(base, input_warp=(0.0, 5.0, 10.0))
    >>> [name for name in kernel.parameters.names if "warp" in name]
    ['input_warp.scale', 'input_warp.increment0', 'input_warp.increment1']
    >>> kernel.QUASISEPARABLE
    True
    >>> WarpedKernel(base, input_warp=(0.0, 5.0, 2.0))
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: WarpedKernel's input_warp= knot locations ...

    The identity warp is the base kernel, exactly:

    >>> import numpy as np
    >>> fixed = Matern32(0.4, 2.0)
    >>> warped = WarpedKernel(
    ...     fixed, input_warp=(0.0, 5.0, 10.0), increments=0.0, input_scale=1.0
    ... )
    >>> grid = np.linspace(0.0, 10.0, 7)[:, None]
    >>> plain = fixed.matrix(grid, grid, fixed.resolve(None))
    >>> bool(np.array_equal(warped.matrix(grid, grid, warped.resolve(None)), plain))
    True
    """

    FAMILY: ClassVar[str] = "warped"
    #: A class-level *capability*, narrowed per instance from the base: warping
    #: preserves an exact semiseparable representation, so a warped kernel is
    #: quasiseparable exactly when the kernel it warps is.
    QUASISEPARABLE: ClassVar[bool] = True
    #: This is the family that warps, so a composite and a solver can ask.
    WARPS: ClassVar[bool] = True
    #: The label the base kernel's hyperparameters are qualified with.
    LABEL: ClassVar[str] = "base"

    #: Default shrinkage scales. Modest rather than vague on purpose: the
    #: identity warp must be the prior's centre of mass, not merely inside it.
    DEFAULT_INPUT_SCALE: ClassVar[float] = 0.5
    DEFAULT_AMPLITUDE_SCALE: ClassVar[float] = 0.3

    def __init__(
        self,
        base: Kernel,
        *,
        input_warp: Any = None,
        amplitude_warp: Any = None,
        increments: Any = None,
        levels: Any = None,
        input_scale: Any = None,
        amplitude_scale: Any = None,
        non_centred: bool = True,
        axes: Sequence[str] | None = None,
    ) -> None:
        if not isinstance(base, Kernel):
            raise LikelihoodError(
                f"WarpedKernel warps a Kernel, got {type(base).__name__}. Kernels are declared "
                f"neutrally (family name plus Parameter hyperparameters) so that a lowering rule "
                f"can translate them."
            )
        if base._selective:
            raise LikelihoodError(
                "WarpedKernel's base kernel declares axes= of its own. A warp is a map of one "
                "ordered coordinate and everything under it acts on that same coordinate, so "
                "the selection belongs on the WarpedKernel: "
                "WarpedKernel(Matern32(...), input_warp=..., axes=('spectral_axis',))."
            )
        self._declare_axes(axes)
        self._input_knots = _warp_knots(input_warp, "input_warp")
        self._amplitude_knots = _warp_knots(amplitude_warp, "amplitude_warp")
        if not self._input_knots and not self._amplitude_knots:
            raise LikelihoodError(
                "WarpedKernel was given neither input_warp= nor amplitude_warp=, so it is its "
                "base kernel wearing an extra parameter namespace. Declare at least one warp, "
                "or use the base kernel directly."
            )
        self._non_centred = bool(non_centred)
        child = _selected_copy(base, self._axes)
        self._terms: tuple[tuple[str, Kernel], ...] = ((self.LABEL, child),)
        for parameter in child.parameters:
            self.register_parameter(parameter.rename(f"{self.LABEL}.{parameter.name}"))
        if self._input_knots:
            self._declare_warp(
                "input_warp",
                "increment",
                len(self._input_knots) - 1,
                increments,
                input_scale,
                self.DEFAULT_INPUT_SCALE,
            )
        if self._amplitude_knots:
            self._declare_warp(
                "amplitude_warp",
                "level",
                len(self._amplitude_knots),
                levels,
                amplitude_scale,
                self.DEFAULT_AMPLITUDE_SCALE,
            )
        self.__dict__["HYPERPARAMETERS"] = tuple(self.parameters.names)
        self.__dict__.update(_adopted_capabilities(type(self).__name__, (base,)))
        self.__dict__["QUASISEPARABLE"] = bool(base.QUASISEPARABLE)

    # -- declaration ---------------------------------------------------------

    def _declare_warp(
        self,
        prefix: str,
        stem: str,
        count: int,
        given: Any,
        scale_given: Any,
        default_scale: float,
    ) -> None:
        """Declare one warp's shrinkage scale and its ``count`` knot variables."""
        scale_name = f"{prefix}.scale"
        self.register_parameter(
            _warp_hyperparameter(
                scale_name,
                st.halfnorm(scale=default_scale) if scale_given is None else scale_given,
                positive=True,
            )
        )
        default = (
            st.norm(0.0, 1.0)
            if self._non_centred
            else HierarchicalPrior("norm", {"scale": scale_name}, kwds={"loc": 0.0})
        )
        hierarchical = False
        for index, declaration in enumerate(
            _warp_values(given, count, default, f"{prefix} {stem}s")
        ):
            parameter = _warp_hyperparameter(f"{prefix}.{stem}{index}", declaration, positive=False)
            hierarchical = hierarchical or isinstance(parameter.prior, HierarchicalPrior)
            self.register_parameter(parameter)
        if self._non_centred or hierarchical or self.parameters[scale_name].fixed:
            return
        raise LikelihoodError(
            f"{scale_name!r} is a free parameter that nothing depends on: with "
            f"non_centred=False the shrinkage scale reaches the knot variables only through a "
            f"HierarchicalPrior referencing it, and none of the {prefix} {stem}s declares one. "
            f"An unidentified sampled dimension has no posterior and costs a sampler real "
            f"work, so it is refused here. Either leave non_centred=True (the default, where "
            f"the kernel forms {stem}_k = {scale_name.split('.')[-1]} * z_k itself), pass "
            f"{stem}s=HierarchicalPrior('norm', {{'scale': {scale_name!r}}}, "
            f"kwds={{'loc': 0.0}}), or hold the scale at a number."
        )

    @property
    def base(self) -> Kernel:
        """The kernel being warped."""
        return self._terms[0][1]

    @property
    def terms(self) -> tuple[tuple[str, Kernel], ...]:
        return self._terms

    @property
    def input_knots(self) -> tuple[float, ...]:
        """The input warp's knot locations; empty when there is no input warp."""
        return self._input_knots

    @property
    def amplitude_knots(self) -> tuple[float, ...]:
        """The amplitude warp's knot locations; empty when there is no amplitude warp."""
        return self._amplitude_knots

    def spec(self) -> KernelSpec:
        """The neutral description, with the knot locations in ``metadata``.

        The knots are not hyperparameters — nothing samples them — but two
        warps with different knots are different models, so they must hash
        differently. ``non_centred`` is recorded for the same reason it is a
        keyword rather than a convention: the two forms declare different
        parameters (``z_k`` against ``u_k``), so a chain stored under one
        cannot be read as the other.
        """
        metadata: list[tuple[str, Any]] = []
        if self._input_knots:
            metadata.append(("input_warp_knots", self._input_knots))
        if self._amplitude_knots:
            metadata.append(("amplitude_warp_knots", self._amplitude_knots))
        metadata.append(("non_centred", self._non_centred))
        return dataclasses.replace(super().spec(), metadata=tuple(metadata))

    # -- composition-time checks --------------------------------------------

    def check_axes(self, observed: Any, *, owner: str) -> None:
        """One coordinate, and the base kernel's own rule on that coordinate.

        The single-unit rule is discharged by the stronger requirement: a warp
        is a monotone map of **one** ordered coordinate — that is what keeps
        the quasiseparable ordering precondition, and what makes "the
        coordinate the residuals are stationary in" a well-posed idea at all —
        so anything but one selected axis is refused by name.
        """
        axis_names = tuple(axis.name for axis in observed.axes)
        self._resolve_columns(axis_names)
        selected = axis_names if self._axes is None else self._axes
        if len(selected) != 1:
            raise LikelihoodError(
                f"{owner} was given a WarpedKernel over a {type(observed).__name__} whose "
                f"selected axes are {selected}. A warp is a monotone map of one ordered "
                f"coordinate; name the axis it acts on with axes=(...), and compose across "
                f"axes with Product."
            )
        self.base.check_axes(observed, owner=owner)

    # -- the warp ------------------------------------------------------------

    def _column(self, points: Any) -> Any:
        """The one coordinate column this kernel warps, as an ``(n,)`` array."""
        selected = self.select(points)
        width = int(np.shape(selected)[1])
        if width != 1:
            raise LikelihoodError(
                f"WarpedKernel warps one ordered coordinate, but was handed points with "
                f"{width} coordinates each. check_axes refuses this at composition time; a "
                f"direct call with hand-built coordinates reaches it here. Pass axes=(...) so "
                f"the kernel knows which column to warp."
            )
        return selected[:, 0]

    def _knot_variables(self, prefix: str, stem: str, count: int, resolved: Any) -> list[Any]:
        """The realised knot variables, non-centring undone where it was declared."""
        ops = self.ops
        raw = [ops.scalar(resolved[f"{prefix}.{stem}{index}"]) for index in range(count)]
        if not self._non_centred:
            return raw
        scale = ops.scalar(resolved[f"{prefix}.scale"])
        return [scale * value for value in raw]

    def _piecewise(self, axis: Any, knots: tuple[float, ...], gradients: list[Any]) -> Any:
        """The piecewise-linear function with these segment gradients, anchored at zero.

        Written in the hinge basis — ``g_0 (x - x_0) + Σ_k (g_k - g_{k-1})
        relu(x - x_k)`` — rather than by bucketing the coordinate, because the
        hinge form needs no ``searchsorted``, is one expression in every array
        namespace, is differentiable in the gradients everywhere and in the
        coordinate away from the knots, and extends linearly past both ends
        with the end segments' own gradients, which is the behaviour a warp
        wants outside its knot range.
        """
        ops = self.ops
        total = gradients[0] * (axis - knots[0])
        previous = gradients[0]
        for index in range(1, len(gradients)):
            total = total + (gradients[index] - previous) * _relu(ops, axis - knots[index])
            previous = gradients[index]
        return total

    def _warp_input(self, axis: Any, resolved: Any) -> Any:
        """``w(x)``: the monotone input warp, the identity when none is declared.

        Written as ``x + δ(x)`` rather than as a cumulative sum of knot images.
        The two are the same function; only the offset form makes the identity
        warp *exact*, because at the identity point every gradient of ``δ`` is
        bit-for-bit zero, so ``δ(x)`` is zero and ``x + 0.0`` is ``x``.
        """
        if not self._input_knots:
            return axis
        ops = self.ops
        normaliser = _softplus(ops, ops.scalar(0.0))
        gradients = [
            _softplus(ops, value) / normaliser - 1.0
            for value in self._knot_variables(
                "input_warp", "increment", len(self._input_knots) - 1, resolved
            )
        ]
        return axis + self._piecewise(axis, self._input_knots, gradients)

    def _warp_amplitude(self, axis: Any, resolved: Any) -> Any | None:
        """``a(x)``, or ``None`` when no amplitude warp is declared."""
        if not self._amplitude_knots:
            return None
        ops = self.ops
        knots = self._amplitude_knots
        levels = self._knot_variables("amplitude_warp", "level", len(knots), resolved)
        gradients = [
            (levels[index + 1] - levels[index]) / (knots[index + 1] - knots[index])
            for index in range(len(knots) - 1)
        ]
        return ops.exp(levels[0] + self._piecewise(axis, knots, gradients))

    def warped_coordinate(self, axis: Any, values: Mapping[str, Any]) -> Any:
        """``w(x)``: the coordinate this kernel is stationary in.

        Composed through the base, so a warp of a warp is a warp: the
        recursion must run where the *innermost* kernel's generators were
        built, which for ``WarpedKernel(WarpedKernel(k, w_i), w_o)`` is
        ``w_i(w_o(x))``. Every base but another warping kernel contributes the
        identity here, so the ordinary case is one map.
        """
        resolved = self.resolve(values)
        moved = self._warp_input(self.ops.scalar(axis), resolved)
        return self.base.warped_coordinate(moved, self._child_values(self.LABEL, resolved))

    def warp_provenance(self, values: Mapping[str, Any]) -> dict[str, Any] | None:
        """The warp as a JSON-plain record: knots in, knots out, amplitudes there.

        Enough to reconstruct and to plot the map, and small enough to travel
        as an attribute on a diagnostics group. Called with concrete values —
        it is a reporting surface, not part of any traced evaluation.
        """
        resolved = self.resolve(values)
        record: dict[str, Any] = {
            "kind": "warped",
            "base": self.base.FAMILY,
            "non_centred": self._non_centred,
        }
        if self._input_knots:
            knots = self.ops.scalar(np.asarray(self._input_knots, dtype=DTYPE))
            record["input_warp_knots"] = list(self._input_knots)
            record["input_warp_images"] = _as_floats(self._warp_input(knots, resolved))
        if self._amplitude_knots:
            knots = self.ops.scalar(np.asarray(self._amplitude_knots, dtype=DTYPE))
            record["amplitude_warp_knots"] = list(self._amplitude_knots)
            record["amplitude_warp_values"] = _as_floats(self._warp_amplitude(knots, resolved))
        return record

    # -- evaluation ----------------------------------------------------------

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> Any:
        raise LikelihoodError(
            "a WarpedKernel is not stationary, so it has no closed form in one separation: "
            "k(x, x') depends on where x and x' are, not only on how far apart they are — "
            "that is what a warp is for. Call matrix() or diagonal(), which is what every "
            "solver does."
        )

    def value(self, separation: Any, values: Mapping[str, Any]) -> Any:
        """Refused: ``k(τ)`` presumes stationarity, which is exactly what a warp drops."""
        return self._covariance(separation, values)

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> Any:
        """``a(x) k(w(x), w(x')) a(x')`` — the dense face of the warp."""
        ops = self.ops
        resolved = self.resolve(values)
        points = ops.points(left)
        other = ops.points(right, dimensions=int(np.shape(points)[1]))
        left_axis = self._column(points)
        right_axis = self._column(other)
        block = self.base.matrix(
            self._warp_input(left_axis, resolved)[:, None],
            self._warp_input(right_axis, resolved)[:, None],
            self._child_values(self.LABEL, resolved),
        )
        left_scale = self._warp_amplitude(left_axis, resolved)
        if left_scale is None:
            return block
        right_scale = self._warp_amplitude(right_axis, resolved)
        return left_scale[:, None] * block * right_scale[None, :]

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> Any:
        """``a(x)² k(0)``: the prior variance, which the amplitude warp moves."""
        ops = self.ops
        resolved = self.resolve(values)
        points = ops.points(coordinates)
        axis = self._column(points)
        prior = self.base.diagonal(
            self._warp_input(axis, resolved)[:, None],
            self._child_values(self.LABEL, resolved),
        )
        scale = self._warp_amplitude(axis, resolved)
        return prior if scale is None else prior * scale * scale

    def __repr__(self) -> str:
        warps = []
        if self._input_knots:
            warps.append(f"input_warp={self._input_knots!r}")
        if self._amplitude_knots:
            warps.append(f"amplitude_warp={self._amplitude_knots!r}")
        if self._axes is not None:
            warps.append(f"axes={self._axes!r}")
        return f"WarpedKernel({self.base!r}, {', '.join(warps)})"


def find_warped_composite(kernel: Kernel, path: tuple[str, ...] = ()) -> tuple[str, ...] | None:
    """The label path to a warp that sits **under a composite**, or ``None``.

    ``_find_nested_product``'s shape, and for the same kind of reason: a
    quasiseparable solve needs one coordinate for its propagators, and a
    :class:`Sum` or :class:`Product` whose terms warp differently has none. The
    supported composition is the other way round —
    ``WarpedKernel(Sum(...), input_warp=...)``, one warp of the coordinate with
    several kernels on it — and a warp of a warp is fine too, because those
    compose into a single monotone map.

    Returned as a *path* so the refusal can name the offending term rather than
    the enclosing kernel, which is what W5.2 established for products.
    """
    composite = isinstance(kernel, _Composite)
    for label, child in kernel.terms:
        # Depth first, so the *deepest* offending term is the one named: in
        # ``Sum(k, Sum(k, WarpedKernel(...)))`` the useful answer is
        # ``term1.term1``, not the enclosing sum that merely contains it.
        found = find_warped_composite(child, (*path, label))
        if found is not None:
            return found
        if composite and child.warps:
            return (*path, label)
    return None


def refuse_warped_composite(kernel: Kernel, owner: str) -> None:
    """Refuse a warp that sits under a composite, naming the term (W5.7).

    Called by every backend's quasiseparable solver from ``check_compatible``,
    so the refusal is one sentence rather than three. ``DenseGP`` evaluates each
    term on its own coordinate and is unaffected, which is what the message
    offers as the way out.
    """
    found = find_warped_composite(kernel)
    if found is None:
        return
    location = ".".join(found)
    raise LikelihoodError(
        f"{owner} cannot lower this {type(kernel).__name__}: its term {location!r} is a "
        f"WarpedKernel, and a composite whose terms warp has no single coordinate for the "
        f"recursion to run on — each warped term's generators are built at its own w(x) while "
        f"the propagators come from one axis. Warp the composite instead of its terms, as in "
        f"WarpedKernel(Sum(...), input_warp=...), which is one warp of the coordinate with "
        f"several kernels on it; or use DenseGP, which evaluates each term on its own "
        f"coordinate."
    )


def _as_floats(array: Any) -> list[float]:
    """A concrete array — numpy, torch or jax — as a list of Python floats."""
    return [float(value) for value in np.asarray(_to_numpy(array), dtype=DTYPE).ravel()]


def _to_numpy(array: Any) -> Any:
    """``array`` as something numpy can read, without importing any backend.

    ``np.asarray`` handles jax arrays and numpy arrays directly; a torch tensor
    needs ``detach``/``cpu`` first, and duck-typing for those two methods is
    how a core module asks for them without importing torch.
    """
    detach = getattr(array, "detach", None)
    if detach is not None:
        array = detach()
        to_cpu = getattr(array, "cpu", None)
        if to_cpu is not None:
            array = to_cpu()
    return array


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
        ``k(0)``, a scalar — or, since **W5.7**, an ``(n,)`` array where the
        kernel's marginal variance varies along the coordinate (an amplitude
        warp's ``a(x)² k(0)``). Every consumer adds it to the noise diagonal,
        so the two shapes broadcast identically and nothing downstream has to
        distinguish them.

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
        blocks.append(builder(child, kernel._child_values(label, resolved), axis))
    marginal = blocks[0].marginal
    for block in blocks[1:]:
        marginal = marginal + block.marginal
    return CeleriteRepresentation(
        decay=ops.concatenate([block.decay for block in blocks], axis=0),
        left=ops.concatenate([block.left for block in blocks], axis=-1),
        right=ops.concatenate([block.right for block in blocks], axis=-1),
        marginal=marginal,
    )


def warped_representation(
    kernel: Kernel, values: Mapping[str, Any], axis: Any
) -> CeleriteRepresentation:
    r"""A warped kernel is the base kernel's generators, moved and scaled (W5.7).

    Two algebraic facts, both exact, and neither costing a rank:

    **Input warping.** If :math:`K_{nm} = \sum_j U_{nj} V_{mj} e^{-c_j (t_n -
    t_m)}` represents :math:`k`, then evaluating the *same* builder at
    :math:`w(t)` represents :math:`k(w(\cdot), w(\cdot))` — provided the
    recursion's own propagators are formed from :math:`w(t)` too, which is
    what :meth:`Kernel.warped_coordinate` tells the solver. Monotonicity is
    what makes this legal: it is the condition under which :math:`w(t)` is
    still sorted, which the recursion requires. The builder here receives the
    **raw** sorted axis (the registry's documented contract, unchanged) and
    applies the warp itself, so the two agree by construction rather than by
    the caller remembering.

    **Amplitude warping.** :math:`D K D` with :math:`D = \mathrm{diag}(a)` has
    entries :math:`a_n K_{nm} a_m = \sum_j (a_n U_{nj})(a_m V_{mj}) e^{-c_j
    (t_n - t_m)}`, so scaling both generator blocks by :math:`a` is the whole
    of it, at unchanged rank :math:`J` and unchanged decays. The marginal
    variance becomes per-point, :math:`a_n^2 k(0)`, which is why
    :class:`CeleriteRepresentation`'s ``marginal`` is allowed to be an
    ``(n,)`` array as well as a scalar: every consumer adds it to the noise
    diagonal, and the two shapes broadcast identically.

    The base term is fetched through the registry like any other, so a warp
    over a user-registered kernel works with no further registration, and a
    warp over a family with no representation is refused **by name** here.
    """
    ops = kernel.ops
    resolved = kernel.resolve(values)
    base = kernel.terms[0][1]
    builder = lookup_quasiseparable_term(base.FAMILY, owner=type(kernel).__name__)
    label = kernel.terms[0][0]
    # **This level's** warp only, not the composed one: the base's builder
    # applies the base's own warp when the base is itself a warping kernel, so
    # handing it the composed coordinate would warp twice. The composition is
    # ``warped_coordinate``'s job, because that is what the *solver* needs.
    warped = kernel._warp_input(ops.scalar(axis), resolved)  # type: ignore[attr-defined]
    block = builder(base, kernel._child_values(label, resolved), warped)
    scale = kernel._warp_amplitude(ops.scalar(axis), resolved)  # type: ignore[attr-defined]
    if scale is None:
        return block
    return CeleriteRepresentation(
        decay=block.decay,
        left=scale[:, None] * block.left,
        right=scale[:, None] * block.right,
        marginal=scale * scale * block.marginal,
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
    (WarpedKernel, warped_representation),
):
    register_quasiseparable_term(_kernel_type, _builder, builtin=True)
del _kernel_type, _builder
