"""TMNRE's two ideas, as the pieces :class:`~ampere.inference.SBIEngine` drives.

Private module. Truncated marginal neural ratio estimation (Miller et al. 2021,
*Truncated Marginal Neural Ratio Estimation*, NeurIPS 34) is expressed here
through ``sbi`` 0.27's own ``NRE`` trainers and ``RestrictedPrior`` — W3.4's
ruling of 2026-09-09 — rather than by reviving the archived swyft
implementation, whose last release pins ``pytorch-lightning <= 1.9.5`` and
cannot be installed beside the torch 2.13 the ``sbi`` extra resolves to.

What TMNRE is, in the two halves this module implements
--------------------------------------------------------
**Marginal ratio estimation.** Rather than one classifier for the joint ratio
``r(θ, x) = p(θ|x)/p(θ)``, train one per *marginal* of interest: one for each
1-D ``θ_i`` and, optionally, one for each 2-D pair. Each is the same ``NRE``
trainer fed ``theta[:, idx]`` against the same ``x``, so a marginal estimator
never has to represent the joint's correlations and stays trainable as the
parameter count grows — which is the property that makes the method work at
high dimension, and the corner plot is a set of marginals anyway.

Two consequences worth stating, because both are why this is cheap. A marginal
trainer needs **no prior at all** while it trains: ``NRE``'s loss contrasts a
batch's own ``(θ, x)`` pairs against its own permuted ones, so the "marginal"
side of the classification comes from the *training set's* own θ column and
nothing asks a distribution for a density. And the marginal posterior is
recovered without a closed-form marginal prior — which ampere's unconstrained
joint prior does not have — because

    ``p(θ_i | x) ∝ r(θ_i, x) · q(θ_i)``

where ``q`` is the θ-distribution the rows were **actually drawn from**. That
is a density this module can estimate directly from the same column the
estimator trained on, and :func:`marginal_log_density` does exactly that with a
Gaussian kernel density estimate (:class:`scipy.stats.gaussian_kde`, Scott's
rule).

*Why a KDE and not a histogram* (the choice W3.4 asks to be recorded): the
number this feeds is a **threshold crossing** — the interval where the
estimated marginal exceeds ``ε`` times its own maximum — so what matters is
that the estimate be smooth near its tails, where a histogram is a staircase of
zero-count bins whose crossing point moves by a whole bin width when one draw
moves between bins. A KDE at Scott's bandwidth is smooth by construction, needs
no bin choice to record in the provenance, costs nothing at these sizes (a few
thousand points evaluated on a few hundred grid nodes), and its own
over-smoothing biases the box **outwards**, which is the safe direction for a
truncation: a box slightly too wide costs simulations, a box too narrow throws
away posterior mass that no later round can recover.

**Truncation.** After each round, restrict the *prior* to the hyperrectangle
where the estimated 1-D marginals put their mass — per parameter, the interval
on which the marginal posterior exceeds ``ε`` times its maximum — simulate the
next round inside it, and retrain. The restricted prior is the original prior
renormalised on a subset, **not** a learned proposal, and that is the whole
reason no importance correction appears anywhere here: inside the box the
target is unchanged, the estimate stays amortised *within* the box, and the box
itself is a diagnostic a reader can check (:class:`TruncationBox`, recorded per
round in ``ampere_sbi_truncation``).

What is lost, and it is not small: **amortisation across observations.** The
box is computed from the estimated marginals *at the observed data*, so a
truncated estimator is only meaningful for that observation, exactly as any
multi-round SBI run is. A one-round ``method="nre"`` fit stays amortised; a
TMNRE run does not, and the run says so in ``ampere_sbi_amortised``.

Where the boundary of this module is
-------------------------------------
Everything here is *pure* with respect to the engine: boxes, grids, densities,
the marginal estimator wrapper and the ``marginals`` group. The round loop that
calls them is :meth:`ampere.inference.SBIEngine._run_tmnre`, in ``_sbi.py``
beside the loop it is a variant of, so this module imports nothing from there
and the dependency runs one way.
"""

from __future__ import annotations

import dataclasses
import functools
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

__all__ = [
    "DEFAULT_TRUNCATION_EPSILON",
    "GRID_POINTS_1D",
    "GRID_POINTS_2D",
    "MARGINALS_GROUP",
    "MARGINALS_SCHEMA_VERSION",
    "MARGINAL_ORDERS",
    "MarginalEstimator",
    "MarginalSummary",
    "TMNREArtefact",
    "TruncationBox",
    "attach_marginals",
    "constrained_grid",
    "grid_between",
    "interval_above",
    "marginal_indices",
    "marginal_log_density",
    "pair_labels",
    "pair_mesh",
    "rebuild_restricted_prior",
    "restricted_prior_class",
]

#: The default threshold on the estimated 1-D marginal posterior, relative to
#: that marginal's own maximum: the truncation interval is the span of grid
#: nodes whose density exceeds ``ε · max``. 1e-4 rather than something smaller
#: because the estimate under it is a neural ratio times a kernel density,
#: neither of which is worth believing four decades below the peak; and
#: rather than something larger because a box that clips the posterior's tails
#: is an error no later round can undo.
DEFAULT_TRUNCATION_EPSILON = 1e-4

#: Grid nodes per parameter for the 1-D marginals — the resolution of both the
#: box edges and the stored summary. Odd, so the grid has a centre node.
GRID_POINTS_1D = 129

#: Grid nodes per axis for a 2-D pair. Quadratic in storage and in evaluation,
#: so much coarser than the 1-D grid: a pair's summary is read as a contour,
#: never as a box edge, and a 33-by-33 mesh draws a contour.
GRID_POINTS_2D = 33

#: The ``marginals=`` vocabulary: 1 trains one estimator per parameter, 2 adds
#: one per unordered pair. Nothing above 2 — a corner plot has no third panel,
#: and the count of estimators is already quadratic at 2.
MARGINAL_ORDERS: tuple[int, ...] = (1, 2)

#: The name of the derived group a TMNRE run carries beside its ``posterior``.
MARGINALS_GROUP = "marginals"

#: Dimension names in that group. ``marginal_row``/``marginal_column`` are the
#: two axes of a pair's grid and are deliberately *different* names for the
#: same length, because a square array whose two axes carry different
#: coordinates cannot share one dimension in xarray.
PARAMETER_DIM = "marginal_parameter"
PAIR_DIM = "marginal_pair"
NODE_DIM = "marginal_node"
ROW_DIM = "marginal_row"
COLUMN_DIM = "marginal_column"

#: Bumped when the shape of the ``marginals`` group changes.
MARGINALS_SCHEMA_VERSION = 1


# ---------------------------------------------------------------------------
# The box
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True, slots=True)
class TruncationBox:
    """A hyperrectangle in **unconstrained** parameter space.

    The truncation, and the run's own diagnostic. In unconstrained coordinates
    because that is where the estimator, the prior bridge and every θ this
    driver hands ``sbi`` live (``_sbi.py``'s module docstring); the constrained
    edges are recorded beside them in :meth:`to_dict` so that a reader of an
    archived run sees the box in the coordinates they declared.

    An unbounded box (``-inf``/``+inf`` on every edge) is the state before the
    first round, and :meth:`accepts_everything` says so; a
    :class:`~sbi.utils.RestrictedPrior` is never built from one.
    """

    lower: tuple[float, ...]
    upper: tuple[float, ...]

    def __post_init__(self) -> None:
        if len(self.lower) != len(self.upper):
            raise ValueError(
                f"a truncation box needs one lower edge per upper edge, got "
                f"{len(self.lower)} and {len(self.upper)}."
            )

    @classmethod
    def unbounded(cls, size: int) -> TruncationBox:
        """The whole of ℝⁿ — round 1's proposal, before anything is estimated."""
        return cls(tuple([-math.inf] * int(size)), tuple([math.inf] * int(size)))

    @property
    def size(self) -> int:
        return len(self.lower)

    @property
    def accepts_everything(self) -> bool:
        """True while no edge is finite, which is the untruncated state."""
        return not any(math.isfinite(edge) for edge in (*self.lower, *self.upper))

    @property
    def widths(self) -> np.ndarray:
        return np.asarray(self.upper, dtype=float) - np.asarray(self.lower, dtype=float)

    @property
    def log_volume(self) -> float:
        """``sum(log width)``, not the product.

        The product underflows for a narrow box in more than a few dimensions,
        and the quantity a reader wants is a *ratio* of successive volumes —
        which is a difference of these. ``+inf`` for an unbounded box, so the
        first round's shrinkage is a comparison that holds trivially rather
        than a special case in every reader; ``-inf`` for a box some edge of
        which has collapsed to nothing.
        """
        widths = self.widths
        if not np.all(np.isfinite(widths)):
            return math.inf
        if np.any(widths <= 0.0):
            return -math.inf
        return float(np.sum(np.log(widths)))

    def contains(self, theta: Any) -> bool:
        """Whether one unconstrained θ vector lies inside, edges included."""
        row = np.asarray(theta, dtype=float).reshape(-1)
        if row.size != self.size:
            raise ValueError(
                f"a {self.size}-dimensional truncation box was asked about a "
                f"{row.size}-dimensional θ."
            )
        return bool(
            np.all(row >= np.asarray(self.lower, dtype=float))
            and np.all(row <= np.asarray(self.upper, dtype=float))
        )

    def intersect(self, other: TruncationBox) -> TruncationBox:
        """The tightest box inside both — how a round's box is made nested.

        A round's grid already spans only the previous box, so the raw edges
        are inside it by construction; intersecting anyway makes "the box never
        grows" a property of this class rather than of the caller's arithmetic,
        which is what the acceptance criterion asserts.
        """
        if other.size != self.size:
            raise ValueError(
                f"cannot intersect a {self.size}-dimensional box with a "
                f"{other.size}-dimensional one."
            )
        lower = np.maximum(np.asarray(self.lower, dtype=float), np.asarray(other.lower, float))
        upper = np.minimum(np.asarray(self.upper, dtype=float), np.asarray(other.upper, float))
        return TruncationBox(tuple(lower.tolist()), tuple(upper.tolist()))

    def indicator(self) -> _BoxIndicator:
        """The ``accept_reject_fn`` a ``RestrictedPrior`` takes.

        A module-level callable object rather than a closure, because it is
        reached by two routes that a closure does not survive: a trained
        posterior carries its prior, and W3.5's artefact store **pickles** the
        posterior; and a problem under a process pool pickles everything it
        sends. The same reasoning gave :func:`ampere.inference._sbi._rebuild_prior`
        its ``__reduce__``.
        """
        return _BoxIndicator(np.asarray(self.lower, dtype=float), np.asarray(self.upper, float))

    def to_dict(self, constrain: Any = None) -> dict[str, Any]:
        """JSON-safe, for ``ampere_sbi_truncation``.

        *constrain* is an optional callable mapping one unconstrained vector to
        the constrained coordinates the user declared; given, the box's edges
        are recorded in both parameterisations. The bijections are per
        parameter and therefore diagonal (``tests/inference/test_sbi.py``
        checks that they are), so an edge maps edge-wise.
        """
        recorded: dict[str, Any] = {
            "lower": [float(value) for value in self.lower],
            "upper": [float(value) for value in self.upper],
            "log_volume": self.log_volume,
        }
        if constrain is not None and not self.accepts_everything:
            lower = np.asarray(constrain(np.asarray(self.lower, dtype=float)), dtype=float)
            upper = np.asarray(constrain(np.asarray(self.upper, dtype=float)), dtype=float)
            recorded["lower_constrained"] = [float(value) for value in lower.reshape(-1)]
            recorded["upper_constrained"] = [float(value) for value in upper.reshape(-1)]
        return recorded


@dataclasses.dataclass(frozen=True)
class _BoxIndicator:
    """``theta -> bool`` on a box, as a picklable object with a torch signature.

    ``sbi`` calls this with a ``(batch, dimension)`` tensor and wants a
    ``(batch,)`` boolean tensor back. torch is not imported here — the
    comparison is done on the tensor's own operators, which is both faster than
    a round trip through numpy and the only form that keeps this module free of
    the lazy-import rule ``_sbi.py`` opens with.
    """

    lower: np.ndarray
    upper: np.ndarray

    def __call__(self, theta: Any) -> Any:
        lower = theta.new_tensor(self.lower)
        upper = theta.new_tensor(self.upper)
        return ((theta >= lower) & (theta <= upper)).all(dim=-1)


# ---------------------------------------------------------------------------
# Marginals: which ones, where they are evaluated, and what they say
# ---------------------------------------------------------------------------


def marginal_indices(size: int, order: int) -> tuple[tuple[int, ...], ...]:
    """Every marginal of *order* over *size* parameters, in a fixed order.

    ``order=1`` gives ``(0,), (1,), ...``; ``order=2`` gives the unordered
    pairs in lexicographic order. Fixed, because the index of a marginal is a
    coordinate in the stored group and a run archived today must still be
    readable against a later one.
    """
    if order == 1:
        return tuple((index,) for index in range(int(size)))
    if order == 2:
        return tuple(
            (first, second) for first in range(int(size)) for second in range(first + 1, int(size))
        )
    raise ValueError(
        f"marginals of order {order} are not defined; MARGINAL_ORDERS is {MARGINAL_ORDERS}."
    )


def pair_labels(labels: Sequence[str], pairs: Sequence[tuple[int, ...]]) -> tuple[str, ...]:
    """``"a|b"`` per pair — the ``marginal_pair`` coordinate.

    ``|`` rather than ``,`` or ``-`` because a merged parameter name may
    contain a dot and an array-valued one a bracket, and ``|`` appears in
    neither (``inference.md`` §4.5's name rules).
    """
    return tuple("|".join(labels[index] for index in pair) for pair in pairs)


def marginal_log_density(samples: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """``log q`` on *grid*, from the θ columns the estimator actually trained on.

    *samples* is ``(count, k)`` and *grid* is ``(nodes, k)``; the return is
    ``(nodes,)``. This is the proposal's own marginal, estimated by a Gaussian
    KDE at Scott's bandwidth — the module docstring records why a KDE rather
    than a histogram, and why the estimate's over-smoothing errs in the safe
    direction for a truncation.

    A degenerate column (every draw identical, which a tied or effectively
    fixed parameter can produce) has no kernel density; the return is then a
    flat zero, which makes the marginal posterior the bare ratio and leaves the
    box to the ratio alone rather than failing the run.
    """
    import scipy.stats as stats

    columns = np.asarray(samples, dtype=float)
    if columns.ndim == 1:
        columns = columns.reshape(-1, 1)
    nodes = np.asarray(grid, dtype=float)
    if nodes.ndim == 1:
        nodes = nodes.reshape(-1, 1)
    try:
        kernel = stats.gaussian_kde(columns.T)
        density = np.asarray(kernel(nodes.T), dtype=float).reshape(-1)
    except (np.linalg.LinAlgError, ValueError):
        return np.zeros(nodes.shape[0], dtype=float)
    return np.log(np.clip(density, 1e-300, None))


def interval_above(
    grid: np.ndarray, log_density: np.ndarray, epsilon: float
) -> tuple[float, float]:
    """The span of nodes whose density exceeds ``epsilon`` times the maximum.

    ``min``/``max`` of the crossing nodes rather than the connected component
    around the peak: a multimodal marginal must keep **both** modes, and a
    truncation that quietly dropped one would be a bug whose only symptom is a
    posterior that is confidently wrong. The cost is that a box around two
    well-separated modes also contains the gap between them, which is the right
    trade in a method whose whole safety argument is that the box never cuts
    posterior mass.

    A marginal that is flat, or all non-finite, returns the whole grid.
    """
    nodes = np.asarray(grid, dtype=float).reshape(-1)
    values = np.asarray(log_density, dtype=float).reshape(-1)
    finite = np.isfinite(values)
    if not finite.any():
        return float(nodes.min()), float(nodes.max())
    peak = float(values[finite].max())
    keep = finite & (values - peak > math.log(float(epsilon)))
    if not keep.any():  # pragma: no cover - the peak itself always crosses
        return float(nodes.min()), float(nodes.max())
    return float(nodes[keep].min()), float(nodes[keep].max())


def grid_between(lower: float, upper: float, nodes: int) -> np.ndarray:
    """*nodes* points spanning ``[lower, upper]``, widened if the two coincide.

    A degenerate interval happens when every draw of a parameter landed on one
    value, and a zero-width grid would make every downstream quantity a
    constant with no interval to threshold.
    """
    low, high = float(lower), float(upper)
    if not (math.isfinite(low) and math.isfinite(high)):
        raise ValueError(f"a marginal grid needs finite edges, got [{low}, {high}].")
    if high <= low:
        pad = max(abs(low), 1.0) * 1e-6
        low, high = low - pad, high + pad
    return np.linspace(low, high, int(nodes))


# ---------------------------------------------------------------------------
# The trained marginal estimators
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class MarginalEstimator:
    """One ``sbi`` ratio estimator over a subset of the parameter vector.

    Holds the trainer across rounds — ``append_simulations(..., from_round=r)``
    accumulates, exactly as the joint trainer does — and the estimator the last
    :meth:`train` returned.
    """

    indices: tuple[int, ...]
    trainer: Any
    estimator: Any = None

    def append(self, theta: Any, summary: Any, *, round_index: int) -> None:
        """This round's rows, projected onto :attr:`indices`."""
        columns = theta[:, list(self.indices)]
        self.trainer.append_simulations(columns, summary, from_round=int(round_index))

    def train(self, **options: Any) -> Any:
        self.estimator = self.trainer.train(show_train_summary=False, **options)
        return self.estimator

    def log_ratio(
        self, points: np.ndarray, observation: Any, *, torch: Any, dtype: Any
    ) -> np.ndarray:
        """``log r(θ_idx, x_obs)`` at every row of *points*.

        ``sbi`` 0.27's :meth:`~sbi.neural_nets.ratio_estimators.RatioEstimator.
        unnormalized_log_ratio` refuses to broadcast — it compares the shape
        prefixes of θ and x and raises when they differ — so the observation is
        expanded to one row per grid node here rather than passed once. That is
        cheap (a view, not a copy) and it is the only spelling the estimator
        accepts.
        """
        if self.estimator is None:  # pragma: no cover - the loop always trains first
            raise RuntimeError("this marginal estimator has not been trained yet.")
        nodes = np.asarray(points, dtype=float).reshape(-1, len(self.indices))
        theta = torch.as_tensor(nodes, dtype=dtype, device=observation.device)
        # ``observation`` is one summary: ``(features,)`` under a flat layout
        # and ``(rows, columns)`` under a set one. Reshaping through its own
        # shape rather than to ``(1, -1)`` keeps the set packing intact, which
        # a flatten would silently destroy.
        shape = tuple(observation.shape)
        condition = observation.reshape(1, *shape).expand(nodes.shape[0], *shape)
        with torch.no_grad():
            values = self.estimator.unnormalized_log_ratio(theta, condition)
        return np.asarray(values.detach().cpu().numpy(), dtype=float).reshape(-1)


# ---------------------------------------------------------------------------
# The stored summary
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class MarginalSummary:
    """What the ``marginals`` group holds: the estimators, evaluated on a grid.

    The estimators themselves are torch modules and do not belong in a netCDF
    file; what a reader wants from them is the curve each one draws over the
    final box, which is what a corner plot is made of and what this carries.

    ``log_ratio`` is the network's own output; ``log_density`` is
    ``log_ratio + log q`` normalised to a maximum of zero — the estimated
    marginal posterior up to its own normalisation, which is the quantity the
    truncation thresholds and the one a plot draws. Both are stored because
    only the first is the network's and only the second is interpretable.
    """

    labels: tuple[str, ...]
    order: int
    epsilon: float
    grid: np.ndarray
    grid_constrained: np.ndarray
    log_ratio: np.ndarray
    log_density: np.ndarray
    pairs: tuple[tuple[int, ...], ...] = ()
    pair_grid_row: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros((0, 0)))
    pair_grid_column: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros((0, 0)))
    pair_log_ratio: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros((0, 0, 0)))
    pair_log_density: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros((0, 0, 0)))

    def to_dataset(self, attrs: Mapping[str, Any] | None = None) -> Any:
        """The group as an :class:`xarray.Dataset`, ready to hang off a run.

        Built here rather than through :func:`ampere.results.emit` for the
        reason W3.6's ``calibration`` group records: a group attached after
        emission needs no new hook in a §4 contract this item does not own, and
        ``xarray.DataTree`` assignment is all attaching a child is.
        """
        import xarray

        variables: dict[str, Any] = {
            "grid": ((PARAMETER_DIM, NODE_DIM), self.grid),
            "grid_constrained": ((PARAMETER_DIM, NODE_DIM), self.grid_constrained),
            "log_ratio": ((PARAMETER_DIM, NODE_DIM), self.log_ratio),
            "log_density": ((PARAMETER_DIM, NODE_DIM), self.log_density),
        }
        coords: dict[str, Any] = {PARAMETER_DIM: list(self.labels)}
        if self.pairs:
            variables["pair_grid_row"] = ((PAIR_DIM, ROW_DIM), self.pair_grid_row)
            variables["pair_grid_column"] = ((PAIR_DIM, COLUMN_DIM), self.pair_grid_column)
            variables["pair_log_ratio"] = ((PAIR_DIM, ROW_DIM, COLUMN_DIM), self.pair_log_ratio)
            variables["pair_log_density"] = (
                (PAIR_DIM, ROW_DIM, COLUMN_DIM),
                self.pair_log_density,
            )
            coords[PAIR_DIM] = list(pair_labels(self.labels, self.pairs))
        return xarray.Dataset(variables, coords=coords, attrs=dict(attrs or {}))


@dataclasses.dataclass(frozen=True)
class TMNREArtefact:
    """What W3.5's artefact store holds for a TMNRE run: the whole answer.

    A cache entry has to reproduce the run it replaces. For NPE, NLE and NRE
    that is the trained posterior and nothing else; for TMNRE the emitted run
    also carries the truncation history and the ``marginals`` group, neither of
    which can be recomputed from the posterior alone, so a hit that stored only
    the posterior would silently emit a *different* run under the same key.

    Everything here pickles as it stands: the box is plain floats, its
    indicator a dataclass rather than a closure
    (:meth:`TruncationBox.indicator`), the summary numpy arrays, and the
    posterior is what the store already knew how to hold. The trained marginal
    *networks* are deliberately not here — a run stores their curves, not their
    weights, and reloading weights nothing reads would make every hit pay for
    them.
    """

    posterior: Any
    truncation: Any
    history: list[dict[str, Any]]
    marginals: Any


def attach_marginals(tree: Any, marginals: Any) -> Any:
    """Hang a :meth:`MarginalSummary.to_dataset` off a run as ``marginals``.

    The same two lines :func:`ampere.results.attach_calibration` is, and for
    the same reason: a run holds one set of marginals, the most recent, and
    xarray cannot hold two children under one name anyway.
    """
    import xarray

    tree[MARGINALS_GROUP] = xarray.DataTree(marginals)
    return tree


def constrained_grid(
    problem: Any, grid: np.ndarray, index: int, reference: np.ndarray
) -> np.ndarray:
    """*grid* along parameter *index*, mapped into the user's own coordinates.

    ``FittingProblem.constrain`` maps a whole vector, and the bijections are
    per parameter — diagonal, as ``tests/inference/test_sbi.py``'s Jacobian row
    measures independently — so varying one coordinate of *reference* and
    reading back that same coordinate is the parameter's own bijection.
    """
    nodes = np.asarray(grid, dtype=float).reshape(-1)
    mapped = np.empty_like(nodes)
    row = np.asarray(reference, dtype=float).copy()
    for position, value in enumerate(nodes):
        row[index] = value
        mapped[position] = float(np.asarray(problem.constrain(row), dtype=float)[index])
    return mapped


def pair_mesh(row: np.ndarray, column: np.ndarray) -> np.ndarray:
    """The cartesian product of two grids as ``(n*m, 2)``, row-major.

    A named function rather than an inline ``meshgrid`` because the *order* is
    what lets the flat vector of ratios be reshaped back to ``(n, m)`` and
    stored against ``(marginal_row, marginal_column)`` without ambiguity about
    which axis is which.
    """
    first, second = np.meshgrid(
        np.asarray(row, dtype=float), np.asarray(column, dtype=float), indexing="ij"
    )
    return np.stack([first.reshape(-1), second.reshape(-1)], axis=1)


# ---------------------------------------------------------------------------
# The truncated prior
# ---------------------------------------------------------------------------


@functools.cache
def restricted_prior_class(base: Any) -> Any:
    """``sbi``'s ``RestrictedPrior``, with two defaults changed and stated.

    Built lazily from the class handed in — this module imports no ``sbi``,
    for the reason ``_sbi.py``'s docstring gives — and cached, so repeated runs
    share one subclass and ``isinstance`` checks inside ``sbi`` keep meaning
    what they look like they mean.

    Two changes, both defaults rather than behaviour:

    * ``print_rejected_frac`` defaults to **False**. ``sbi``'s
      ``RestrictedPrior.sample`` prints its rejection rate on every call, and
      a rejection-sampled posterior calls it many times per draw batch. The
      rate is worth having, so it is kept — on :attr:`acceptance_rate`, and in
      the run's ``ampere_sbi_truncation`` — rather than printed. This is the
      same rule the engine applies to ``sbi``'s progress bars.
    * ``log_prob`` defaults to ``norm_restricted_prior=False``. The
      normalising factor is the prior mass of the box: a **constant in θ**,
      which changes neither rejection sampling nor MCMC nor the shape of the
      posterior, and which ``sbi`` estimates by drawing 10 000 *accepted*
      samples through the box every time the factor is not already cached.
      Against ampere's prior — a Python loop over
      :meth:`~ampere.core.dataset.FittingProblem.sample_prior`, because ties
      and hierarchical priors are resolved per draw — that is the single most
      expensive thing a TMNRE run can be made to do, and it buys an additive
      constant. ``norm_restricted_prior=True`` still works if a caller asks
      for it explicitly, so nothing is taken away.

    The consequence is stated where it matters: the log-density this prior
    reports is **unnormalised**, and so, therefore, is the potential a TMNRE
    posterior evaluates — which is what ``ampere_sbi_log_prob_kind`` has always
    said for a ratio-based method anyway.
    """

    class QuietRestrictedPrior(base):  # type: ignore[misc, valid-type]
        """A :class:`~sbi.utils.RestrictedPrior` that is quiet and unnormalised."""

        def sample(self, sample_shape: Any = None, **options: Any) -> Any:
            options.setdefault("print_rejected_frac", False)
            options.setdefault("save_acceptance_rate", True)
            shape = () if sample_shape is None else sample_shape
            return super().sample(shape, **options)

        def log_prob(self, theta: Any, **options: Any) -> Any:
            options.setdefault("norm_restricted_prior", False)
            return super().log_prob(theta, **options)

        def __reduce__(self) -> tuple[Any, tuple[Any, ...]]:
            # The class is built inside this function, so pickle cannot find it
            # by name; a trained TMNRE posterior carries this prior and W3.5's
            # artefact store pickles the posterior. Rebuilt through the
            # module-level factory, exactly as the unconstrained prior is.
            return rebuild_restricted_prior, (self._prior, self._accept_reject_fn, self._device)

    return QuietRestrictedPrior


def rebuild_restricted_prior(prior: Any, accept_reject_fn: Any, device: str) -> Any:
    """Unpickle hook for :func:`restricted_prior_class`'s lazily-built subclass."""
    import sbi.utils  # pyrefly: ignore[missing-import]

    return restricted_prior_class(sbi.utils.RestrictedPrior)(
        prior, accept_reject_fn, sample_with="rejection", device=device
    )
