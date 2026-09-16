r"""The Vecchia approximation: the sparse, non-rank-limited candidate (**W5.6**).

The other half of W5.6's bake-off (:mod:`ampere.core.efgp` is the first), and
the one ``DEVELOPMENT_PLAN.md`` §5 describes as "the one that handles rough
processes without a rank limit". Like its sibling this is a **reference-only
measurement prototype**: numpy and scipy, no torch or jax twin, and distinct
from the declared :class:`~ampere.core.VecchiaGP` strategy slot, which stays
closed until the bake-off's decision-log row is acted on.

The method
----------
Vecchia (1988), in the modern form of Datta et al. (2016) and Guinness (2018).
Any joint density factorises exactly into a chain of conditionals,

.. math::
    p(y_1, \ldots, y_N) = \prod_{i} p\big(y_i \mid y_1, \ldots, y_{i-1}\big),

and the approximation is to condition each term on only the ``k`` **nearest**
of its predecessors:

.. math::
    p(\mathbf y) \;\approx\; \prod_i p\big(y_i \mid y_{g(i)}\big),
    \qquad |g(i)| \le k .

Each factor is a :math:`k`-dimensional Gaussian conditional — one
:math:`k \times k` Cholesky — so the whole density costs :math:`O(N k^3)` with
:math:`O(N k)` storage, at *any* dimension and for *any* kernel. There is no
basis, no rank, no spectral density and no box: nothing here asks the kernel to
be smooth, which is exactly why the method is in this bake-off beside two
spectral ones.

**Response Vecchia**, not latent. The chain is applied to the *observations*,
whose covariance is :math:`A = K + D`, rather than to a latent field with the
noise added afterwards. That is the variant that computes the quantity this
contract asks for — the marginal likelihood of the residuals — in one pass, and
the one Guinness's ordering study measures.

What the approximation actually produces
----------------------------------------
Not a likelihood only: a sparse **inverse Cholesky**. Writing each conditional
as :math:`y_i = \mathbf b_i^{\mathsf T} y_{g(i)} + \sqrt{v_i}\,\varepsilon_i`,
the vectors collect into an :math:`N \times N` matrix :math:`U` with
:math:`k+1` non-zeros per column, and

.. math::
    A^{-1} \approx U U^{\mathsf T},
    \qquad \log|A| \approx \sum_i \log v_i .

Every quantity the :class:`~ampere.core.GPSolver` contract asks for follows
from :math:`U` at :math:`O(N k)`: the quadratic form is
:math:`\|U^{\mathsf T} r\|^2`, and :math:`\operatorname{diag}(A^{-1})` — the
one number the leave-one-out identity needs — is the row-wise sum of
:math:`U^2`. So :meth:`VecchiaResponseGP.conditional_loo` is exact in the
approximation rather than refused, as HSGP's is and for the same reason: the
identity is applied to the covariance this solver actually scores.

The ordering is part of the approximation
------------------------------------------
Guinness (2018) is the reference, and its finding is the one thing about
Vecchia that surprises people: the **ordering** matters more than ``k``. A
coordinate (raster) ordering on an image conditions each pixel on a
half-neighbourhood that is systematically on one side, and is markedly worse
than a random one at equal ``k``; maximum-minimum-distance ("maximin")
ordering is best and random is close behind it. So ``ordering`` and ``seed``
are dataclass fields like ``neighbours`` — they change the value of the
likelihood, so they are a *declaration* that enters ``Likelihood.to_spec``'s
hash, and they are repeated in
:meth:`VecchiaResponseGP.provenance_config`. W5.6's table measures both
orderings so the claim is this repository's own rather than a citation.

What it costs that the spectral methods do not
-----------------------------------------------
The neighbour sets are a function of the coordinates alone — no fitted value
enters — so they are built once and cached (:func:`neighbour_structure`),
exactly as a spectral solver's grid or box is a constant of the data. The
build is the method's one super-linear step and is reported separately in the
bake-off table rather than hidden inside a per-evaluation number.
"""

from __future__ import annotations

import dataclasses
import hashlib
import math
from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import scipy.linalg
import scipy.spatial

from .exceptions import LikelihoodError
from .kernels import Kernel, _as_float64, _as_points
from .likelihood import GPConditional, GPSolver, _components
from .results_schema import FunctionSamples

__all__ = [
    "DEFAULT_NEIGHBOURS",
    "ORDERINGS",
    "NeighbourStructure",
    "VecchiaResponseGP",
    "neighbour_structure",
    "separation_covariance",
]

DTYPE = np.float64

_LOG_2PI = math.log(2.0 * math.pi)

#: Conditioning-set size when nothing says otherwise. Guinness (2018) finds
#: the accuracy gain flattening around 30 for two-dimensional Matérn fields,
#: and the cost is cubic in it, so this is the knee rather than a maximum.
DEFAULT_NEIGHBOURS: int = 30

#: The orderings this prototype implements. ``"maximin"`` — the literature's
#: best — is deliberately absent: a faithful maximin ordering is a
#: farthest-point traversal, and W5.6 measures ``"random"`` against
#: ``"coordinate"`` to establish that the ordering matters at all before
#: paying for the better one.
ORDERINGS: tuple[str, ...] = ("random", "coordinate")

#: How many neighbour structures are kept. Small: each is ``O(N k)`` integers,
#: and a fit reuses one structure for its whole run.
_CACHE_LIMIT = 4
_CACHE: dict[tuple[str, ...], NeighbourStructure] = {}


# ---------------------------------------------------------------------------
# The conditioning sets: a constant of the data
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class NeighbourStructure:
    """The ordering and the conditioning sets, for one coordinate set.

    Attributes
    ----------
    order
        ``(N,)`` permutation: ``order[p]`` is the index, in the caller's own
        numbering, of the sample at position ``p`` of the chain.
    neighbours
        ``(N, k)`` positions in the **chain's** numbering, or ``-1`` where a
        sample has fewer than ``k`` predecessors. Row ``p`` holds the nearest
        predecessors of ``order[p]``, nearest first.
    counts
        ``(N,)`` how many of each row's entries are real.
    ordering
        Which rule produced :attr:`order`, for provenance.
    """

    order: np.ndarray
    neighbours: np.ndarray
    counts: np.ndarray
    ordering: str

    @property
    def size(self) -> int:
        """How many samples the structure covers."""
        return int(self.order.shape[0])

    @property
    def width(self) -> int:
        """``k``, the conditioning-set size it was built for."""
        return int(self.neighbours.shape[1])


def _ordering(points: np.ndarray, ordering: str, seed: int) -> np.ndarray:
    if ordering == "random":
        return np.asarray(np.random.default_rng(seed).permutation(points.shape[0]))
    if ordering == "coordinate":
        return np.asarray(np.lexsort(tuple(points[:, axis] for axis in range(points.shape[1]))))
    raise LikelihoodError(
        f"VecchiaResponseGP's ordering must be one of {list(ORDERINGS)}, got {ordering!r}. The "
        f"ordering changes the value of the approximation (Guinness 2018), which is why it is a "
        f"declared field rather than an internal choice."
    )


def neighbour_structure(
    points: np.ndarray,
    width: int,
    *,
    ordering: str = "random",
    seed: int = 0,
) -> NeighbourStructure:
    """The ordering and nearest-predecessor sets for *points*, cached.

    The search is exact and blocked: for each block of the chain a
    :class:`scipy.spatial.cKDTree` over everything that precedes the block
    supplies candidates, and the block's own predecessors are compared
    directly, so the answer is the true ``k`` nearest and the cost is
    ``O((N/c) N log N)`` rather than ``O(N²)`` distances held at once.

    Cached on the coordinates' bytes together with the three settings, because
    the structure is a **constant of the data** — no hyperparameter enters it —
    and a fit evaluates one likelihood thousands of times over one coordinate
    set. The cache holds :data:`_CACHE_LIMIT` entries and is keyed by content
    rather than by identity, so a caller that rebuilds an equal array still
    hits it.
    """
    block = np.ascontiguousarray(points, dtype=DTYPE)
    key = (
        hashlib.blake2b(block.tobytes(), digest_size=16).hexdigest(),
        str(block.shape),
        str(int(width)),
        str(ordering),
        str(int(seed)),
    )
    cached = _CACHE.get(key)
    if cached is not None:
        return cached
    structure = _build_structure(block, int(width), ordering, int(seed))
    if len(_CACHE) >= _CACHE_LIMIT:
        _CACHE.pop(next(iter(_CACHE)))
    _CACHE[key] = structure
    return structure


def _build_structure(
    points: np.ndarray, width: int, ordering: str, seed: int
) -> NeighbourStructure:
    order = _ordering(points, ordering, seed)
    ordered = np.ascontiguousarray(points[order])
    total = int(ordered.shape[0])
    width = max(1, min(int(width), max(total - 1, 1)))
    neighbours = np.full((total, width), -1, dtype=np.int64)
    counts = np.zeros(total, dtype=np.int64)

    # The first ``width`` samples condition on everything before them.
    head = min(width + 1, total)
    for position in range(1, head):
        neighbours[position, :position] = np.arange(position, dtype=np.int64)
        counts[position] = position

    step = max(256, total // 64)
    for start in range(head, total, step):
        stop = min(start + step, total)
        tree = scipy.spatial.cKDTree(ordered[:start])
        distances, indices = tree.query(ordered[start:stop], k=width)
        distances = np.atleast_2d(distances.T).T.reshape(stop - start, width)
        indices = np.atleast_2d(indices.T).T.reshape(stop - start, width)
        # Candidates inside the block itself, which the tree above cannot see.
        inside = np.linalg.norm(
            ordered[start:stop, None, :] - ordered[start:stop][None, :, :], axis=-1
        )
        rows = np.arange(stop - start)
        inside = np.where(rows[None, :] < rows[:, None], inside, np.inf)
        merged_distances = np.concatenate([distances, inside], axis=1)
        merged_indices = np.concatenate([indices, start + rows[None, :].repeat(stop - start, 0)], 1)
        chosen = np.argpartition(merged_distances, width - 1, axis=1)[:, :width]
        picked = np.take_along_axis(merged_indices, chosen, axis=1)
        order_within = np.argsort(np.take_along_axis(merged_distances, chosen, axis=1), axis=1)
        neighbours[start:stop] = np.take_along_axis(picked, order_within, axis=1)
        counts[start:stop] = width
    return NeighbourStructure(
        order=np.asarray(order, dtype=np.int64),
        neighbours=neighbours,
        counts=counts,
        ordering=ordering,
    )


# ---------------------------------------------------------------------------
# Covariance as a function of separation, through the public kernel API
# ---------------------------------------------------------------------------


def separation_covariance(
    kernel: Kernel,
    separation: np.ndarray,
    values: Mapping[str, Any],
    dimensions: int,
) -> np.ndarray:
    """``k(r)`` for an array of Euclidean separations *r*, of any shape.

    :meth:`~ampere.core.kernels.Kernel.matrix` is written for two coordinate
    *sets* and returns their full cross product, which is the wrong shape for
    a batch of :math:`k \\times k` neighbour blocks — forming it would be the
    dense covariance this solver exists not to form. A stationary kernel is a
    function of separation alone across the axes it selects, so the whole
    batch is evaluated by laying every separation along **one** selected axis
    and asking the kernel once, at ``(1, n)``:

    the answer is the kernel's own arithmetic, from the public method, with no
    private hook and nothing re-derived here.
    """
    flat = np.ascontiguousarray(np.reshape(np.asarray(separation, dtype=DTYPE), -1))
    probe = np.asarray(kernel.select(np.arange(dimensions, dtype=DTYPE)[None, :]))
    axis = round(float(probe[0, 0]))
    left = np.zeros((1, dimensions), dtype=DTYPE)
    right = np.zeros((flat.shape[0], dimensions), dtype=DTYPE)
    right[:, axis] = flat
    covariance = np.asarray(kernel.matrix(left, right, values), dtype=DTYPE)
    return np.asarray(np.reshape(covariance, np.shape(separation)))


# ---------------------------------------------------------------------------
# The solver
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class VecchiaResponseGP(GPSolver):
    """Response-Vecchia GP algebra at ``O(N k³)`` (**W5.6 prototype**).

    An ``EXACT = False`` solver in the comparison class W5.4 established for
    :class:`~ampere.core.HilbertSpaceGP` (``tests/conformance/README.md`` §3)
    — the error falls as ``k`` grows and floors where the conditioning sets
    stop seeing the long-range structure — but the *rate* is the method's own
    rather than the spectral one, and so is the floor: refining ``k`` on a
    rough kernel buys what refining a basis cannot.

    **A measurement prototype, and the reference path only.** W5.6's outcome
    is a decision-log row naming which of this and
    :class:`~ampere.core.EquispacedFourierGP` is promoted to a full
    three-backend solver; until then there is no torch or jax twin, so a
    problem composed with it is a numpy problem and the capability ladder says
    so.

    Parameters
    ----------
    neighbours
        ``k``, the conditioning-set size. Cost is cubic in it.
    ordering
        ``"random"`` (Guinness 2018's near-best, and the default) or
        ``"coordinate"``. It changes the value of the approximation, which is
        why it is declared.
    seed
        The permutation's seed under ``ordering="random"``. Declared for the
        same reason: two seeds are two approximations.
    jitter
        A standard deviation added in quadrature to the noise diagonal, the
        same knob, meaning and default as :class:`~ampere.core.DenseGP`'s.

    Examples
    --------
    >>> from ampere.core import Matern32
    >>> solver = VecchiaResponseGP(neighbours=20)
    >>> solver.latent_size(Matern32(0.3, 2.0), 500)
    500
    >>> VecchiaResponseGP(neighbours=20).provenance_config()
    {'neighbours': 20, 'ordering': 'random', 'seed': 0}
    """

    neighbours: int = DEFAULT_NEIGHBOURS
    ordering: str = "random"
    seed: int = 0
    jitter: float = 0.0

    NAME: ClassVar[str] = "VecchiaResponseGP"
    EXACT: ClassVar[bool] = False
    IMPLEMENTED: ClassVar[bool] = True
    #: The sparse inverse Cholesky applies to as many columns as it is given,
    #: so the circular complex GP's two-column right-hand side costs one
    #: structure build and two triangular products.
    STACKED_RESIDUALS: ClassVar[bool] = True

    def __post_init__(self) -> None:
        if int(self.neighbours) < 1:
            raise LikelihoodError(
                f"VecchiaResponseGP's neighbours must be at least 1, got {self.neighbours!r}. "
                f"A conditioning set of zero is an independent-noise likelihood, which "
                f"IndependentNoise already is."
            )
        if self.ordering not in ORDERINGS:
            raise LikelihoodError(
                f"VecchiaResponseGP's ordering must be one of {list(ORDERINGS)}, got "
                f"{self.ordering!r}."
            )
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"VecchiaResponseGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )

    def provenance_config(self) -> Mapping[str, Any]:
        """The approximation, for the run's attrs — as HSGP's, for its reason."""
        return {
            "neighbours": int(self.neighbours),
            "ordering": str(self.ordering),
            "seed": int(self.seed),
        }

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        """The layout and axis rules, plus this solver's one restriction.

        Every term of the kernel must act on the same axes, for
        :func:`separation_covariance`'s reason: the neighbour blocks are
        evaluated as a function of Euclidean separation across the selected
        axes, which a :class:`~ampere.core.kernels.Product` of terms on
        *different* axes is not. Nothing else is required — no spectral
        density, no quasiseparable form, no smoothness — which is the point of
        the method.
        """
        super().check_compatible(kernel, observed)
        selections = {leaf.axes for leaf in kernel.leaves()}
        if len(selections) > 1:
            raise LikelihoodError(
                f"{self.NAME} evaluates its conditioning blocks as a function of Euclidean "
                f"separation across one axis selection, so every term of the kernel must act on "
                f"the same axes; this one's terms select {sorted(map(str, selections))}. Give "
                f"every term the same axes=, or use DenseGP, which evaluates each term on its "
                f"own axes."
            )

    # -- internals -----------------------------------------------------------

    def _diagonal(self, variance: np.ndarray) -> np.ndarray:
        diagonal = _as_float64(variance, "noise variances") + self.jitter**2
        if not np.all(np.isfinite(diagonal)) or np.any(diagonal < 0.0):
            raise LikelihoodError(
                f"{self.NAME} needs a finite, non-negative noise diagonal, and the conditional "
                f"variances it forms are positive only if K + diag(sigma^2) is. Pass "
                f"{self.NAME}(jitter=...) — a standard deviation in the data's units."
            )
        return diagonal

    def _factorise(
        self,
        kernel: Kernel,
        points: np.ndarray,
        diagonal: np.ndarray,
        values: Mapping[str, Any],
    ) -> tuple[NeighbourStructure, np.ndarray, np.ndarray]:
        r"""``(structure, coefficients, variances)`` — the sparse inverse Cholesky.

        ``coefficients[p]`` is :math:`\mathbf b_p` padded with zeros and
        ``variances[p]`` is :math:`v_p`, both in the chain's own ordering.
        One batched :math:`k \times k` Cholesky does the whole chain; the
        first ``k`` rows, whose conditioning sets are short, are padded with an
        identity block so that they go through the same batch rather than
        through a special case.
        """
        selected = np.ascontiguousarray(kernel.select(points), dtype=DTYPE)
        # The **container's** dimension, not the selection's: it is what
        # ``separation_covariance`` rebuilds a point set in, and what the
        # kernel's own ``select`` will then index into.
        dimensions = int(points.shape[1])
        structure = neighbour_structure(
            selected, int(self.neighbours), ordering=self.ordering, seed=int(self.seed)
        )
        ordered = selected[structure.order]
        noise = diagonal[structure.order]
        total = structure.size
        width = structure.width
        prior = float(np.asarray(separation_covariance(kernel, np.zeros(1), values, dimensions))[0])

        coefficients = np.zeros((total, width), dtype=DTYPE)
        variances = np.empty(total, dtype=DTYPE)
        variances[0] = prior + noise[0]

        # Chunked so the (rows, k, k) block stays inside ~64 MiB.
        budget = max(1, (16 * 1024 * 1024) // (8 * max(width * width, 1)))
        for start in range(1, total, budget):
            stop = min(start + budget, total)
            rows = np.arange(start, stop)
            index = structure.neighbours[start:stop]
            present = index >= 0
            safe = np.where(present, index, 0)
            block = ordered[safe]  # (rows, k, d)
            separation = np.linalg.norm(block[:, :, None, :] - block[:, None, :, :], axis=-1)
            gram = separation_covariance(kernel, separation, values, dimensions)
            gram = gram + np.where(present, noise[safe], 0.0)[:, :, None] * np.eye(width)[None]
            # An absent neighbour becomes an isolated unit row, so the batched
            # Cholesky sees a positive-definite block and contributes nothing.
            absent = ~present
            gram[absent] = 0.0
            gram = np.swapaxes(gram, 1, 2)
            gram[absent] = 0.0
            gram = np.swapaxes(gram, 1, 2)
            diagonal_index = np.arange(width)
            gram[:, diagonal_index, diagonal_index] = np.where(
                present, gram[:, diagonal_index, diagonal_index], 1.0
            )
            cross = separation_covariance(
                kernel,
                np.linalg.norm(block - ordered[rows][:, None, :], axis=-1),
                values,
                dimensions,
            )
            cross = np.where(present, cross, 0.0)
            try:
                lower = np.linalg.cholesky(gram)
            except np.linalg.LinAlgError as error:
                raise LikelihoodError(
                    f"a Vecchia conditioning block is not positive definite ({error}). The "
                    f"blocks are submatrices of K + diag(sigma^2), so this is a numerical "
                    f"rather than a structural failure: raise the jitter, or reduce neighbours."
                ) from error
            weights = np.asarray(
                scipy.linalg.solve_triangular(
                    np.swapaxes(lower, 1, 2),
                    scipy.linalg.solve_triangular(lower, cross[:, :, None], lower=True),
                    lower=False,
                )
            )[:, :, 0]
            coefficients[start:stop] = weights
            variances[start:stop] = prior + noise[rows] - np.sum(cross * weights, axis=1)
        if np.any(variances <= 0.0):
            raise LikelihoodError(
                f"{self.NAME} produced a non-positive conditional variance, so the factorised "
                f"density it forms is not a density. Raise the jitter, or reduce neighbours."
            )
        return structure, coefficients, variances

    def _apply_transpose(
        self,
        structure: NeighbourStructure,
        coefficients: np.ndarray,
        variances: np.ndarray,
        ordered_residual: np.ndarray,
    ) -> np.ndarray:
        r""":math:`U^{\mathsf T} r` in the chain's ordering, ``(N,)`` or ``(N, k)``."""
        index = structure.neighbours
        safe = np.where(index >= 0, index, 0)
        mask = (index >= 0).astype(DTYPE)
        weighted = coefficients * mask
        if ordered_residual.ndim == 1:
            gathered = ordered_residual[safe]
            return (ordered_residual - np.sum(weighted * gathered, axis=1)) / np.sqrt(variances)
        gathered = ordered_residual[safe]  # (N, k, c)
        return (ordered_residual - np.einsum("nk,nkc->nc", weighted, gathered)) / np.sqrt(
            variances
        )[:, None]

    def _apply(
        self,
        structure: NeighbourStructure,
        coefficients: np.ndarray,
        variances: np.ndarray,
        vector: np.ndarray,
    ) -> np.ndarray:
        r""":math:`U v` in the chain's ordering, the transpose of :meth:`_apply_transpose`."""
        index = structure.neighbours
        safe = np.where(index >= 0, index, 0)
        mask = (index >= 0).astype(DTYPE)
        weighted = coefficients * mask
        if vector.ndim == 1:
            scaled = vector / np.sqrt(variances)
            result = scaled.copy()
            np.add.at(result, safe, -weighted * scaled[:, None])
            return result
        scaled = vector / np.sqrt(variances)[:, None]
        result = scaled.copy()
        np.add.at(result, safe, -weighted[:, :, None] * scaled[:, None, :])
        return result

    def _precision_diagonal(
        self,
        structure: NeighbourStructure,
        coefficients: np.ndarray,
        variances: np.ndarray,
    ) -> np.ndarray:
        r""":math:`\operatorname{diag}(U U^{\mathsf T})` in the chain's ordering."""
        index = structure.neighbours
        safe = np.where(index >= 0, index, 0)
        mask = (index >= 0).astype(DTYPE)
        weighted = coefficients * mask
        total = np.zeros(structure.size, dtype=DTYPE)
        total += 1.0 / variances
        np.add.at(total, safe, (weighted * weighted) / variances[:, None])
        return total

    # -- the interface -------------------------------------------------------

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        """The **Vecchia** marginal likelihood, summed over the columns.

        Not the declared kernel's: :attr:`EXACT` is ``False``. It is the exact
        marginal likelihood of the Gaussian whose precision is
        :math:`U U^{\\mathsf T}`, which converges to
        :math:`(K + D)^{-1}` as the conditioning sets grow.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        diagonal = self._diagonal(variance)
        structure, coefficients, variances = self._factorise(kernel, points, diagonal, values)
        ordered = residuals[structure.order]
        whitened = self._apply_transpose(structure, coefficients, variances, ordered)
        columns = _components(residuals)
        return -0.5 * (
            float(np.sum(whitened * whitened))
            + columns * float(np.sum(np.log(variances)))
            + residuals.size * _LOG_2PI
        )

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Sundararajan & Keerthi's identity, exact in the approximation.

        The sparse inverse Cholesky supplies both halves at :math:`O(N k)`:
        :math:`A^{-1} r = U (U^{\\mathsf T} r)`, and
        :math:`(A^{-1})_{ii}` is the row-wise sum of :math:`U^2`. So the
        leave-one-out terms cost no more than the likelihood does — which the
        two spectral solvers, whose own identity needs the ``(N, m)`` block,
        cannot say.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        diagonal = self._diagonal(variance)
        structure, coefficients, variances = self._factorise(kernel, points, diagonal, values)
        ordered = residuals[structure.order]
        alpha = self._apply(
            structure,
            coefficients,
            variances,
            self._apply_transpose(structure, coefficients, variances, ordered),
        )
        precision = self._precision_diagonal(structure, coefficients, variances)
        rows = int(residuals.shape[0])
        columns = _components(residuals)
        quadratic = np.sum(np.reshape(alpha, (rows, columns)) ** 2, axis=1)
        terms = columns * (0.5 * np.log(precision) - 0.5 * _LOG_2PI) - quadratic / (2.0 * precision)
        restored = np.empty_like(terms)
        restored[structure.order] = terms
        return np.asarray(restored, dtype=DTYPE)

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The GP posterior, each target conditioned on its ``k`` nearest data.

        The prediction is the approximation's own rather than the exact one
        applied to an approximate solve, for
        :class:`~ampere.core.HilbertSpaceGP`'s reason: mixing the two gives a
        variance that is the variance of nothing. Conditioning a target on its
        ``k`` nearest *observations* is the Vecchia chain extended by one term,
        it costs :math:`O(m k^3)` for ``m`` targets and — unlike the exact
        expression :math:`K_* (K+D)^{-1} r`, which is :math:`O(N m)` — it never
        touches the full cross-covariance.
        """
        points = _as_points(coordinates, "data coordinates")
        residuals = _as_float64(residual, "residuals")
        diagonal = self._diagonal(variance)
        selected = np.ascontiguousarray(kernel.select(points), dtype=DTYPE)
        dimensions = int(points.shape[1])
        if at is None:
            targets = selected
        else:
            targets = np.ascontiguousarray(
                kernel.select(_as_points(at, "conditioning grid", dimensions=points.shape[1])),
                dtype=DTYPE,
            )
        width = max(1, min(int(self.neighbours), int(selected.shape[0])))
        tree = scipy.spatial.cKDTree(selected)
        _, index = tree.query(targets, k=width)
        index = np.asarray(index, dtype=np.int64).reshape(targets.shape[0], width)

        block = selected[index]
        separation = np.linalg.norm(block[:, :, None, :] - block[:, None, :, :], axis=-1)
        gram = separation_covariance(kernel, separation, values, dimensions)
        gram = gram + diagonal[index][:, :, None] * np.eye(width)[None]
        cross = separation_covariance(
            kernel, np.linalg.norm(block - targets[:, None, :], axis=-1), values, dimensions
        )
        lower = np.linalg.cholesky(gram)
        weights = np.asarray(
            scipy.linalg.solve_triangular(
                np.swapaxes(lower, 1, 2),
                scipy.linalg.solve_triangular(lower, cross[:, :, None], lower=True),
                lower=False,
            )
        )[:, :, 0]
        prior = float(np.asarray(separation_covariance(kernel, np.zeros(1), values, dimensions))[0])
        posterior = prior - np.sum(cross * weights, axis=1)
        if residuals.ndim == 1:
            mean = np.sum(weights * residuals[index], axis=1)
        else:
            mean = np.einsum("nk,nkc->nc", weights, residuals[index])
        return GPConditional(
            mean=np.asarray(mean, dtype=DTYPE),
            variance=np.asarray(np.maximum(posterior, 0.0), dtype=DTYPE),
        )

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        r"""``f = U^{-T} z`` — the sparse whitening, from ``N`` variables.

        The one place the two candidates differ in kind. A reduced-rank solver
        whitens through ``m`` basis coefficients and
        :meth:`~ampere.core.GPSolver.latent_size` shrinks accordingly; Vecchia
        has no rank to shrink to, so its latent block is ``N`` — which is the
        contract's default, and a real cost under NUTS. What it gives back is
        that the factor is *sparse*: :math:`U^{\mathsf T}` has ``k + 1``
        non-zeros per row, so the triangular solve is :math:`O(N k)` rather
        than the dense :math:`O(N^2)`, and the draw is of the **noise-free**
        process, built from the same chain with ``D = 0``.
        """
        del jitter
        points = _as_points(coordinates, "data coordinates")
        draws = _as_float64(whitened, "whitened latent draws")
        total = int(_as_points(coordinates, "data coordinates").shape[0])
        if np.shape(draws)[0] != total:
            raise LikelihoodError(
                f"{self.NAME} whitens {total} value(s) — one per retained sample, this solver "
                f"having no reduced rank — but it was handed {np.shape(draws)[0]}."
            )
        structure, coefficients, variances = self._factorise(
            kernel, points, np.zeros(total, dtype=DTYPE), values
        )
        ordered = draws[structure.order]
        root = np.sqrt(variances)
        index = structure.neighbours
        safe = np.where(index >= 0, index, 0)
        mask = (index >= 0).astype(DTYPE)
        weighted = coefficients * mask
        if ordered.ndim == 1:
            field = np.zeros(total, dtype=DTYPE)
            for position in range(total):
                field[position] = (
                    root[position] * ordered[position] + weighted[position] @ field[safe[position]]
                )
            result = np.empty_like(field)
            result[structure.order] = field
            return result
        field_block = np.zeros((total, ordered.shape[1]), dtype=DTYPE)
        for position in range(total):
            field_block[position] = (
                root[position] * ordered[position]
                + weighted[position] @ field_block[safe[position]]
            )
        result_block = np.empty_like(field_block)
        result_block[structure.order] = field_block
        return result_block
