"""A point estimate, stored as what it is rather than as a one-draw posterior.

W6.7 (``docs/design/inference_extensions_memo.md`` §3.2, ``results.md`` §4's
*Amended W6.7* note). A run's ``posterior`` group is the wrong container for a
point: ArviZ would compute an R-hat over it, a corner plot would draw a
histogram of one value, and a reader could not tell a mode from a single
sampled draw. :class:`Optimum` is instead a small, frozen, typed record of what
an optimiser found — the mode in **both** coordinate systems, the density at it
in **both** conventions, the curvature there as a covariance (or a refusal by
name), how the optimiser got there, and the problem's provenance at the time —
with :meth:`Optimum.to_datatree` for the one results format and the existing
netCDF route.

The two coordinate systems and the two conventions
--------------------------------------------------
*unconstrained* is the packed free vector an optimiser actually moved (every
bounded parameter mapped to the real line by its bijection, ``lowering.md``
§2); *constrained* is the same point in the parameters' own units, keyed by
qualified parameter name. The density is recorded twice because the two
numbers answer different questions and a single ``log_prob`` field would be
ambiguous:

``log_prob_constrained``
    ``log p(θ) + log p(D | θ)`` at the mode — the constrained-space posterior
    density, which is **what the optimisers maximise** (``inference.md``
    §10b). It is comparable with a run's ``sample_stats.lp``.
``log_prob_unconstrained``
    The same plus the change-of-variables term ``Σ log|dθ/du|`` — the density
    a NUTS kernel works with. It is *not* maximised at the same point unless
    every bijection is the identity.

The covariance is the inverse Hessian of the (negated) objective **in the
unconstrained coordinates**, because that is where a start ball is drawn
(:meth:`ampere.inference.Engine.initial_positions`'s ``around=``) and where the
Hessian of a bounded parameter is well defined at the boundary's approach. When
the Hessian is not positive definite the covariance is ``None`` and
:attr:`Optimum.covariance_refusal` says why — never a silently regularised
matrix, which would hand a sampler a confident ball around a saddle.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .exceptions import ResultsError
from .provenance import ATTR_PREFIX, canonical_json, hash_of

__all__ = ["OPTIMUM_GROUP", "Optimum", "StartSummary", "aligned"]

#: The one group :meth:`Optimum.to_datatree` writes. Deliberately not
#: ``posterior``: nothing that reads a run's posterior should find a mode there.
OPTIMUM_GROUP = "optimum"

#: The dimension the vectors run along, and the second one the covariance needs.
_DIM = "free_parameter"
_DIM2 = "free_parameter_2"

#: :attr:`Optimum.coordinates`' two entries (W7.6).
_UNCONSTRAINED = "unconstrained"
_COORDINATE_KINDS = frozenset({"constrained", _UNCONSTRAINED})


@dataclass(frozen=True)
class StartSummary:
    """One start of a multi-start optimisation, summarised.

    Attributes
    ----------
    start_hash
        :func:`~ampere.results.provenance.hash_of` of the start vector (packed,
        unconstrained) — enough to tell two starts apart and to check a rerun
        began where the record says, without storing every vector.
    objective
        ``log_prob_constrained`` where that start finished; ``-inf`` if it
        never reached a scoreable point.
    status
        ``"converged"``, ``"not converged: <optimiser message>"`` or
        ``"failed: <reason>"``.
    """

    start_hash: str
    objective: float
    status: str

    def to_dict(self) -> dict[str, Any]:
        return {"start_hash": self.start_hash, "objective": self.objective, "status": self.status}

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> StartSummary:
        return cls(
            start_hash=str(payload["start_hash"]),
            objective=_float(payload["objective"]),
            status=str(payload["status"]),
        )


@dataclass(frozen=True, eq=False)
class Optimum:
    """What an optimiser found, and enough to reproduce and reuse it.

    Built by :func:`ampere.inference.optimise` and
    :func:`ampere.inference.warm_start_gp`; consumed by
    ``Engine.initial_positions(count, around=optimum)`` and every sampling
    engine's ``run(initial=optimum)``. The module docstring explains the two
    coordinate systems and the two density conventions.

    Attributes
    ----------
    route
        ``"scipy"``, ``"map"``, ``"vi"`` or ``"empirical_bayes"``; a
        :meth:`combine`\\ d optimum joins its parts' routes with ``"+"``.
    backend
        The problem's backend at the time.
    free_names
        The qualified names of the free parameters this optimum covers, in
        layout order.
    free_labels
        One label per entry of :attr:`unconstrained` (``"offset[2]"`` for an
        array element; equal to :attr:`free_names` when every parameter is a
        scalar). The coordinate of :meth:`to_datatree`'s variables.
    unconstrained
        The mode as the packed unconstrained vector (read-only).
    constrained
        The mode in the parameters' own units, by qualified name.
    log_prob_constrained, log_prob_unconstrained
        The density at the mode in each convention (module docstring). ``nan``
        for a :meth:`combine`\\ d optimum, whose joint point nobody evaluated.
    covariance
        The inverse Hessian of the negated objective at the mode, in the
        unconstrained coordinates, or ``None``.
    covariance_refusal
        Why :attr:`covariance` is ``None``; ``None`` when it is not.
    converged
        Whether the optimiser reported convergence for the kept start.
    message
        The optimiser's own message for the kept start, plus anything the
        route adds (a fallback taken, a temporary solver built).
    evaluations
        Objective evaluations summed over every start (gradient evaluations
        counted with them on the native routes).
    starts
        One :class:`StartSummary` per start, in the order they were run.
    provenance
        :func:`~ampere.results.provenance.provenance_attrs` of the problem
        when the optimum was found, so the start is reproducible from this
        record alone.
    coordinates
        One entry per :attr:`free_labels` entry, ``"constrained"`` or
        ``"unconstrained"``: the coordinates the minimiser moved that entry
        in (W7.6). The ``"scipy"`` route moves a bounded coordinate in its
        constrained value with the bounds passed, every other route moves
        ``u``. Omitted, it is all ``"unconstrained"`` — what every optimum
        stored before W7.6 was. It does not enter :attr:`identity`: the
        point is the same whichever coordinates reached it.
    """

    route: str
    backend: str
    free_names: tuple[str, ...]
    free_labels: tuple[str, ...]
    unconstrained: np.ndarray
    constrained: Mapping[str, Any]
    log_prob_constrained: float
    log_prob_unconstrained: float
    covariance: np.ndarray | None
    covariance_refusal: str | None
    converged: bool
    message: str
    evaluations: int
    starts: tuple[StartSummary, ...] = ()
    provenance: Mapping[str, Any] = field(default_factory=dict)
    coordinates: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        vector = np.array(self.unconstrained, dtype=float).reshape(-1)
        if vector.size != len(self.free_labels):
            raise ResultsError(
                f"an Optimum's unconstrained vector has {vector.size} entries but "
                f"{len(self.free_labels)} free labels; one label per entry is the layout."
            )
        kinds = tuple(str(k) for k in self.coordinates) or (_UNCONSTRAINED,) * vector.size
        if len(kinds) != vector.size or not set(kinds) <= _COORDINATE_KINDS:
            raise ResultsError(
                f"an Optimum's coordinates must be one of {sorted(_COORDINATE_KINDS)} per free "
                f"label ({vector.size}), got {kinds!r}."
            )
        object.__setattr__(self, "coordinates", kinds)
        vector.setflags(write=False)
        object.__setattr__(self, "unconstrained", vector)
        object.__setattr__(self, "free_names", tuple(self.free_names))
        object.__setattr__(self, "free_labels", tuple(self.free_labels))
        object.__setattr__(self, "starts", tuple(self.starts))
        object.__setattr__(self, "constrained", dict(self.constrained))
        object.__setattr__(self, "provenance", dict(self.provenance))
        if self.covariance is not None:
            matrix = np.array(self.covariance, dtype=float)
            if matrix.shape != (vector.size, vector.size):
                raise ResultsError(
                    f"an Optimum's covariance must be {(vector.size, vector.size)}, got "
                    f"{matrix.shape}."
                )
            matrix.setflags(write=False)
            object.__setattr__(self, "covariance", matrix)
            if self.covariance_refusal is not None:
                raise ResultsError(
                    "an Optimum carries either a covariance or a covariance_refusal, not both."
                )
        elif self.covariance_refusal is None:
            raise ResultsError(
                "an Optimum without a covariance must say why in covariance_refusal: a missing "
                "curvature is a fact about the mode a reader needs, not an empty field."
            )

    # -- identity -------------------------------------------------------------

    @property
    def identity(self) -> str:
        """The optimum's identity hash: the normalised vector plus the route.

        What a run seeded from this optimum records in ``ampere_start``: two
        optima that would start a sampler at the same point by the same route
        share it, and nothing else about them (their message, their evaluation
        count) moves it.
        """
        return hash_of(
            {
                "route": self.route,
                "free_labels": list(self.free_labels),
                "unconstrained": [float(v) for v in self.unconstrained],
            }
        )

    def start_record(self) -> dict[str, Any]:
        """The ``ampere_start`` payload (``results.md`` §9, schema 9)."""
        return {
            "route": self.route,
            "identity": self.identity,
            "log_prob_constrained": float(self.log_prob_constrained),
            "evaluations": int(self.evaluations),
            "converged": bool(self.converged),
        }

    # -- views ----------------------------------------------------------------

    @property
    def constrained_vector(self) -> np.ndarray:
        """:attr:`constrained` flattened in layout order, one entry per label."""
        parts = [
            np.ravel(np.asarray(self.constrained[name], dtype=float)) for name in self.free_names
        ]
        return np.concatenate(parts) if parts else np.empty(0)

    @property
    def standard_deviations(self) -> np.ndarray | None:
        """``sqrt(diag(covariance))`` in the unconstrained coordinates, or ``None``."""
        if self.covariance is None:
            return None
        return np.sqrt(np.diag(self.covariance))

    def summary(self) -> str:
        """A plain-text table of the mode, one row per free-vector entry."""
        sd = self.standard_deviations
        width = max([len("parameter"), *(len(label) for label in self.free_labels)])
        lines = [
            (
                f"Optimum by route {self.route!r} on backend {self.backend!r}: "
                f"{'converged' if self.converged else 'NOT converged'} after "
                f"{self.evaluations} evaluations over {max(len(self.starts), 1)} start(s)"
            ),
            (
                f"log p (constrained) = {self.log_prob_constrained:.6g}; "
                f"log p (unconstrained) = {self.log_prob_unconstrained:.6g}"
            ),
            (
                f"{'parameter':<{width}}  {'constrained':>14}  {'unconstrained':>14}  "
                f"{'sd (unc.)':>12}"
            ),
        ]
        for index, label in enumerate(self.free_labels):
            spread = "-" if sd is None else f"{sd[index]:.6g}"
            lines.append(
                f"{label:<{width}}  {self.constrained_vector[index]:>14.6g}  "
                f"{self.unconstrained[index]:>14.6g}  {spread:>12}"
            )
        if self.covariance_refusal is not None:
            lines.append(f"covariance refused: {self.covariance_refusal}")
        if self.message:
            lines.append(f"message: {self.message}")
        return "\n".join(lines)

    # -- composition ----------------------------------------------------------

    @classmethod
    def combine(cls, *optima: Optimum) -> Optimum:
        """One optimum covering the union of several with disjoint free names.

        The use it exists for is :func:`ampere.inference.warm_start_gp`'s
        hyperparameters beside a model optimum, handed together to
        ``initial_positions(around=)``. The vectors are concatenated in the
        order given (``around=`` reorders by name, so the order is immaterial
        there); the covariance is block diagonal when every part has one and
        ``None`` with a refusal naming the parts otherwise; the densities are
        ``nan``, because nobody evaluated the joint point — a combined optimum
        is a start, not a mode. The provenance is the first part's.
        """
        if not optima:
            raise ResultsError("Optimum.combine needs at least one optimum.")
        if len(optima) == 1:
            return optima[0]
        seen: dict[str, str] = {}
        for optimum in optima:
            for name in optimum.free_names:
                if name in seen:
                    raise ResultsError(
                        f"Optimum.combine: {name!r} is covered by both the {seen[name]!r} and the "
                        f"{optimum.route!r} optimum. Combine merges disjoint free names only; "
                        f"choose which estimate of {name!r} to keep before combining."
                    )
                seen[name] = optimum.route
        backends = {optimum.backend for optimum in optima}
        if len(backends) != 1:
            raise ResultsError(
                f"Optimum.combine: the parts were found on different backends {sorted(backends)}."
            )
        refused = [o.route for o in optima if o.covariance is None]
        if refused:
            covariance = None
            refusal = (
                f"combined: the {', '.join(repr(r) for r in refused)} part(s) had no covariance"
            )
        else:
            size = sum(o.unconstrained.size for o in optima)
            covariance = np.zeros((size, size))
            offset = 0
            for optimum in optima:
                n = optimum.unconstrained.size
                assert optimum.covariance is not None
                covariance[offset : offset + n, offset : offset + n] = optimum.covariance
                offset += n
            refusal = None
        constrained: dict[str, Any] = {}
        for optimum in optima:
            constrained.update(optimum.constrained)
        return cls(
            route="+".join(o.route for o in optima),
            backend=backends.pop(),
            free_names=tuple(n for o in optima for n in o.free_names),
            free_labels=tuple(label for o in optima for label in o.free_labels),
            unconstrained=np.concatenate([o.unconstrained for o in optima]),
            constrained=constrained,
            log_prob_constrained=math.nan,
            log_prob_unconstrained=math.nan,
            covariance=covariance,
            covariance_refusal=refusal,
            converged=all(o.converged for o in optima),
            message="combined from: " + "; ".join(f"{o.route}: {o.message}" for o in optima),
            evaluations=sum(o.evaluations for o in optima),
            starts=tuple(s for o in optima for s in o.starts),
            provenance=optima[0].provenance,
            coordinates=tuple(kind for o in optima for kind in o.coordinates),
        )

    # -- the one results format -----------------------------------------------

    def to_datatree(self) -> Any:
        """A ``DataTree`` with an ``optimum`` group and **no** ``posterior``.

        The group holds ``unconstrained`` and ``constrained`` along
        ``free_parameter`` (the :attr:`free_labels`) and, when there is one,
        ``covariance`` over ``(free_parameter, free_parameter_2)``. The tree's
        attrs are :attr:`provenance` plus this record's own fields under
        ``ampere_optimum_*`` — canonical JSON for anything structured, so the
        tree goes through :func:`ampere.results.to_netcdf` unchanged and
        :meth:`from_datatree` reads it back.
        """
        import xarray as xr

        coords = {_DIM: list(self.free_labels)}
        variables: dict[str, Any] = {
            "unconstrained": ((_DIM,), np.asarray(self.unconstrained)),
            "constrained": ((_DIM,), self.constrained_vector),
        }
        if self.covariance is not None:
            coords[_DIM2] = list(self.free_labels)
            variables["covariance"] = ((_DIM, _DIM2), np.asarray(self.covariance))
        group = xr.Dataset(variables, coords=coords)
        tree = xr.DataTree.from_dict({f"/{OPTIMUM_GROUP}": group})
        attrs: dict[str, Any] = dict(self.provenance)
        prefix = f"{ATTR_PREFIX}optimum_"
        attrs.update(
            {
                f"{prefix}route": self.route,
                f"{prefix}backend": self.backend,
                f"{prefix}identity": self.identity,
                f"{prefix}free_names": canonical_json(list(self.free_names)),
                f"{prefix}shapes": canonical_json(
                    {n: list(np.shape(self.constrained[n])) for n in self.free_names}
                ),
                f"{prefix}log_prob_constrained": canonical_json(self.log_prob_constrained),
                f"{prefix}log_prob_unconstrained": canonical_json(self.log_prob_unconstrained),
                f"{prefix}covariance_refusal": (
                    "" if self.covariance_refusal is None else self.covariance_refusal
                ),
                f"{prefix}converged": int(bool(self.converged)),
                f"{prefix}message": self.message,
                f"{prefix}evaluations": int(self.evaluations),
                f"{prefix}starts": canonical_json([s.to_dict() for s in self.starts]),
                f"{prefix}coordinates": canonical_json(list(self.coordinates)),
            }
        )
        tree.attrs.update(attrs)
        return tree

    @classmethod
    def from_datatree(cls, tree: Any) -> Optimum:
        """Read an :meth:`to_datatree` tree (or its netCDF round trip) back."""
        prefix = f"{ATTR_PREFIX}optimum_"
        attrs = dict(tree.attrs)
        if f"{prefix}route" not in attrs or OPTIMUM_GROUP not in tree.children:
            raise ResultsError(
                "this tree is not an Optimum: it has no 'optimum' group or no "
                f"{prefix}route attribute (Optimum.to_datatree writes both)."
            )
        group = tree[OPTIMUM_GROUP].dataset
        labels = tuple(str(v) for v in np.asarray(group[_DIM].values))
        names = tuple(json.loads(attrs[f"{prefix}free_names"]))
        shapes = json.loads(attrs[f"{prefix}shapes"])
        flat = np.asarray(group["constrained"].values, dtype=float)
        constrained: dict[str, Any] = {}
        offset = 0
        for name in names:
            shape = tuple(int(s) for s in shapes[name])
            size = int(np.prod(shape)) if shape else 1
            chunk = flat[offset : offset + size]
            constrained[name] = float(chunk[0]) if not shape else chunk.reshape(shape).copy()
            offset += size
        refusal = str(attrs[f"{prefix}covariance_refusal"]) or None
        covariance = (
            np.asarray(group["covariance"].values, dtype=float)
            if "covariance" in group.data_vars
            else None
        )
        own = {key for key in attrs if key.startswith(prefix)}
        return cls(
            route=str(attrs[f"{prefix}route"]),
            backend=str(attrs[f"{prefix}backend"]),
            free_names=names,
            free_labels=labels,
            unconstrained=np.asarray(group["unconstrained"].values, dtype=float),
            constrained=constrained,
            log_prob_constrained=_float(json.loads(attrs[f"{prefix}log_prob_constrained"])),
            log_prob_unconstrained=_float(json.loads(attrs[f"{prefix}log_prob_unconstrained"])),
            covariance=covariance,
            covariance_refusal=refusal,
            converged=bool(int(attrs[f"{prefix}converged"])),
            message=str(attrs[f"{prefix}message"]),
            evaluations=int(attrs[f"{prefix}evaluations"]),
            starts=tuple(StartSummary.from_dict(s) for s in json.loads(attrs[f"{prefix}starts"])),
            provenance={k: v for k, v in attrs.items() if k not in own},
            # absent from an optimum stored before W7.6: every one moved in u
            coordinates=tuple(json.loads(attrs.get(f"{prefix}coordinates", "[]"))),
        )

    def __repr__(self) -> str:
        values = ", ".join(
            f"{label}={v:.4g}"
            for label, v in zip(self.free_labels, self.constrained_vector, strict=True)
        )
        return f"<Optimum route={self.route!r} {values}>"


def _float(value: Any) -> float:
    """A JSON-normalised float back to a Python float (sentinels included)."""
    sentinels = {"__inf__": math.inf, "__-inf__": -math.inf, "__nan__": math.nan}
    if isinstance(value, str):
        return sentinels[value]
    return float(value)


def aligned(optimum: Optimum, expected: Sequence[str]) -> tuple[np.ndarray, np.ndarray | None]:
    """*optimum*'s vector and covariance reordered to the label order *expected*.

    Raises :class:`ResultsError` naming the difference when the two label sets
    are not the same — the refusal ``initial_positions(around=)`` relies on.
    """
    have = set(optimum.free_labels)
    want = set(expected)
    if have != want:
        missing = sorted(want - have)
        extra = sorted(have - want)
        raise ResultsError(
            f"this Optimum (route {optimum.route!r}) does not cover this problem's free "
            f"parameters: missing {missing}, not in the problem {extra}. An optimum can seed "
            f"only the problem it was found on; Optimum.combine(...) merges a model optimum with "
            f"warm_start_gp's hyperparameters when each covers part of the problem."
        )
    index = {label: i for i, label in enumerate(optimum.free_labels)}
    order = np.asarray([index[label] for label in expected], dtype=int)
    vector = np.asarray(optimum.unconstrained, dtype=float)[order]
    if optimum.covariance is None:
        return vector, None
    return vector, np.asarray(optimum.covariance, dtype=float)[np.ix_(order, order)]
