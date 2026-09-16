"""The Hilbert-space basis a reduced-rank spectral GP solve is built on (W5.4).

``DEVELOPMENT_PLAN.md`` §5's first Phase 5 bullet asks for an approximate
``GPSolver``, and ``docs/design/horizon_notes.md`` §2 names HSGP (Solin &
Särkkä 2020, *Stat. Comput.* 30, 419; Riutort-Mayol et al. 2023, *Stat.
Comput.* 33, 17) as the cheapest candidate and the one most useful under NUTS.
The idea in one paragraph:

On a bounded box :math:`\\Omega = \\prod_k [-L_k, L_k]` the Dirichlet
Laplacian has the separable eigenpairs

.. math::
    \\lambda_{\\mathbf j} = \\sum_k \\Big(\\frac{\\pi j_k}{2 L_k}\\Big)^2,
    \\qquad
    \\phi_{\\mathbf j}(\\mathbf x) = \\prod_k \\frac{1}{\\sqrt{L_k}}
        \\sin\\!\\Big(\\frac{\\pi j_k (x_k + L_k)}{2 L_k}\\Big),

and a stationary kernel, being a function of the Laplacian through its
spectral density :math:`S`, is approximated on :math:`\\Omega` by

.. math::
    k(\\mathbf x, \\mathbf x') \\approx
    \\sum_{\\mathbf j} S\\big(\\sqrt{\\lambda_{\\mathbf j}}\\big)\\,
    \\phi_{\\mathbf j}(\\mathbf x)\\,\\phi_{\\mathbf j}(\\mathbf x')
    \\quad\\Longleftrightarrow\\quad
    K \\approx \\Phi \\operatorname{diag}(S)\\, \\Phi^{\\mathsf T} .

That is a rank-``m`` factorisation of the covariance, with ``m`` the number of
basis members kept, so every GP quantity follows from Woodbury at
:math:`O(N m + m^3)` — and the latent-GP path needs ``m`` whitened variables
rather than ``N``.

**What lives here and what does not.** This module holds the part that is the
same in every namespace: where the box is, which frequencies the basis runs
over, and how :math:`\\Phi` is evaluated. The box is computed from the data in
numpy once, because it is *data* — it is fixed before a fit starts, no fitted
value enters it, and nothing differentiates it; :func:`basis_matrix` then
builds :math:`\\Phi` through :class:`~ampere.core.kernels.ArrayOps`, so a
backend gets its own tensors and (though :math:`\\Phi` itself carries no
gradient) keeps the whole solve in one namespace. The hyperparameter-dependent
half — :math:`S(\\sqrt\\lambda)` — is the kernel's own
:meth:`~ampere.core.kernels.Kernel.spectral_density`, written once per family
in ``ArrayOps`` for exactly the same reason.

The linear algebra proper lives with the solvers: ``ampere.core.likelihood``'s
:class:`~ampere.core.likelihood.HilbertSpaceGP` in scipy, and one class per
differentiable backend.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from .exceptions import LikelihoodError
from .kernels import NUMPY_OPS, ArrayOps, Kernel, StationaryKernel

__all__ = [
    "DEFAULT_BASIS_SIZE",
    "DEFAULT_BOUNDARY_FACTOR",
    "HilbertSpaceBasis",
    "basis_matrix",
    "check_spectral_support",
    "hilbert_basis",
    "spectral_values",
]

DTYPE = np.float64

#: Basis members per axis when nothing says otherwise. Large enough that the
#: 1-D rows of the conformance battery converge on it, and small enough that
#: an ``m³`` factorisation is free; there is no universally right value, which
#: is why it is a declaration rather than a hidden constant (see
#: :class:`~ampere.core.likelihood.HilbertSpaceGP` on choosing it).
DEFAULT_BASIS_SIZE: int = 32

#: How far past the data the box reaches, as a multiple of the data's own
#: half-extent. Riutort-Mayol et al. (2023) §3 show the approximation is poor
#: within roughly one length scale of the boundary, so the box must contain
#: the data with room to spare; ``2.0`` is their "safe" end of the range and
#: the value ampere's own convergence rows use. It is not free: a wider box
#: needs proportionally more basis members to reach the same frequency.
DEFAULT_BOUNDARY_FACTOR: float = 2.0


@dataclasses.dataclass(frozen=True)
class HilbertSpaceBasis:
    """Where the box is and which frequencies the basis runs over.

    Every field is a **constant of the data**, computed once before any
    hyperparameter is seen: a fit that moves the length scale does not move
    the basis, which is what makes the reduced-rank representation a fixed
    linear one and the latent block a fixed size.

    Attributes
    ----------
    centre
        ``(d,)`` midpoint of the data's extent on each selected axis. The
        basis is built on centred coordinates so the box is symmetric.
    half_width
        ``(d,)`` half-width :math:`L_k` of the box, already multiplied by the
        boundary factor.
    counts
        Basis members per axis, :math:`m_k`. The total is their product.
    frequencies
        ``(m, d)`` angular frequency of each basis member on each axis,
        :math:`\\pi j_k / (2 L_k)`, in the tensor-product order
        :func:`basis_matrix` builds.
    norms
        ``(m,)`` :math:`\\sqrt{\\lambda_{\\mathbf j}}`, the isotropic argument
        a kernel's spectral density is evaluated at.
    boundary_factor
        The multiple of the data's half-extent :attr:`half_width` came from,
        carried for provenance.
    """

    centre: np.ndarray
    half_width: np.ndarray
    counts: tuple[int, ...]
    frequencies: np.ndarray
    norms: np.ndarray
    boundary_factor: float

    @property
    def size(self) -> int:
        """How many basis members there are: ``prod(counts)``."""
        return int(self.norms.shape[0])

    @property
    def dimensions(self) -> int:
        """How many axes the basis spans."""
        return len(self.counts)


def basis_size(counts: Sequence[int]) -> int:
    """The total basis size a per-axis count sequence declares.

    Free of any data, because the latent block's size has to be known at
    composition time, before a container is in hand (``inference.md`` §17.4).
    """
    total = 1
    for count in counts:
        total *= int(count)
    return total


def normalise_counts(basis_size_declaration: int | Sequence[int]) -> tuple[int, ...]:
    """Read a ``basis_size`` declaration as a per-axis tuple, validating it.

    An integer is one axis; a sequence is one entry per axis, in the order the
    kernel selects them.
    """
    if isinstance(basis_size_declaration, (int, np.integer)):
        counts: tuple[int, ...] = (int(basis_size_declaration),)
    else:
        counts = tuple(int(count) for count in basis_size_declaration)
    if not counts:
        raise LikelihoodError(
            "HilbertSpaceGP's basis_size must name at least one axis; an empty sequence "
            "declares no basis at all."
        )
    for count in counts:
        if count < 1:
            raise LikelihoodError(
                f"HilbertSpaceGP's basis_size must be positive on every axis, got {counts!r}."
            )
    return counts


def hilbert_basis(
    points: Any,
    counts: Sequence[int],
    boundary_factor: float = DEFAULT_BOUNDARY_FACTOR,
) -> HilbertSpaceBasis:
    """The box and the frequency grid for *points*, in numpy.

    Parameters
    ----------
    points
        ``(n, d)`` coordinates — already reduced to the axes the kernel
        selects. Read in numpy whatever namespace they arrive in: they are
        data, so nothing is lost by it.
    counts
        Basis members per axis, one entry per column of *points*.
    boundary_factor
        The box's half-width as a multiple of the data's own half-extent.

    Notes
    -----
    A degenerate axis — every sample at one coordinate — has no extent to
    scale, and would give ``L = 0`` and a basis of zeros. It is given a unit
    half-extent instead, which is harmless (the kernel contributes nothing
    along an axis it cannot separate) and keeps the failure out of the middle
    of a factorisation.
    """
    block = NUMPY_OPS.points(points)
    dimensions = int(block.shape[1])
    if len(counts) != dimensions:
        raise LikelihoodError(
            f"the Hilbert-space basis was declared with {len(counts)} axis count(s) "
            f"{tuple(counts)!r} but the coordinates it is built on have {dimensions}."
        )
    if not math.isfinite(boundary_factor) or boundary_factor <= 0.0:
        raise LikelihoodError(
            f"HilbertSpaceGP's boundary_factor must be finite and > 0, got {boundary_factor!r}."
        )
    lower = np.min(block, axis=0)
    upper = np.max(block, axis=0)
    centre = 0.5 * (lower + upper)
    extent = 0.5 * (upper - lower)
    extent = np.where(extent > 0.0, extent, 1.0)
    half_width = np.asarray(boundary_factor * extent, dtype=DTYPE)

    per_axis = [
        math.pi * np.arange(1, int(count) + 1, dtype=DTYPE) / (2.0 * float(width))
        for count, width in zip(counts, half_width, strict=True)
    ]
    # The tensor product, in the same C order ``basis_matrix`` builds its
    # columns in, so column j of Phi and row j of ``frequencies`` are the same
    # basis member.
    grids = np.meshgrid(*per_axis, indexing="ij")
    frequencies = np.stack([grid.reshape(-1) for grid in grids], axis=-1)
    norms = np.sqrt(np.sum(frequencies * frequencies, axis=-1))
    return HilbertSpaceBasis(
        centre=np.asarray(centre, dtype=DTYPE),
        half_width=half_width,
        counts=tuple(int(count) for count in counts),
        frequencies=np.asarray(frequencies, dtype=DTYPE),
        norms=np.asarray(norms, dtype=DTYPE),
        boundary_factor=float(boundary_factor),
    )


def basis_matrix(basis: HilbertSpaceBasis, points: Any, ops: ArrayOps = NUMPY_OPS) -> Any:
    """:math:`\\Phi`, the ``(n, m)`` block of eigenfunctions evaluated at *points*.

    Written through *ops* so the whole solve stays in one namespace — numpy
    for the reference path, torch or jax tensors for a backend's. The box and
    the frequencies come in as ordinary floats, because they are data.

    *points* need not be the coordinates the basis was built from: the same
    :math:`\\Phi` construction evaluates the basis anywhere, which is what
    :meth:`~ampere.core.likelihood.GPSolver.condition` needs for its
    evaluation grid. A point outside the box is not refused — the
    eigenfunctions are defined on the whole line — but the approximation it
    carries is meaningless there, and the boundary factor exists so that does
    not happen.
    """
    block = ops.points(points, dimensions=basis.dimensions)
    columns: Any = None
    for axis in range(basis.dimensions):
        width = float(basis.half_width[axis])
        centre = float(basis.centre[axis])
        count = basis.counts[axis]
        # The per-axis frequencies, recovered from the first ``count`` rows of
        # the C-ordered product: axis 0 varies slowest, the last fastest.
        stride = 1
        for later in basis.counts[axis + 1 :]:
            stride *= later
        frequency = np.ascontiguousarray(basis.frequencies[::stride, axis][:count])
        shifted = block[:, axis : axis + 1] - (centre - width)
        wave = ops.sin(ops.scalar(frequency[None, :]) * shifted) * (width**-0.5)
        if columns is None:
            columns = wave
        else:
            rows = ops.n_points(block)
            columns = (columns[:, :, None] * wave[:, None, :]).reshape(rows, -1)
    return columns


def spectral_values(
    kernel: Kernel,
    basis: HilbertSpaceBasis,
    values: Mapping[str, Any],
) -> Any:
    """:math:`S(\\sqrt{\\lambda_{\\mathbf j}})`, the ``(m,)`` diagonal of the factorisation.

    The one hyperparameter-dependent piece, evaluated in the *kernel's* own
    namespace so a fitted amplitude and length scale keep their gradient.
    """
    return kernel.spectral_density(
        kernel.ops.scalar(basis.norms), values, dimensions=basis.dimensions
    )


def _has_spectral_density(leaf: Kernel) -> bool:
    """Whether *leaf*'s class actually supplies a closed-form spectral density.

    Structural rather than a name list, so a user's own kernel that implements
    :meth:`~ampere.core.kernels.Kernel.spectral_density` reaches the
    reduced-rank path with no registration at all — and a user's
    :class:`~ampere.core.kernels.StationaryKernel` subclass that inherits the
    Matérn form without declaring ``SPECTRAL_NU`` is caught here rather than
    inside a factorisation.
    """
    method = type(leaf).spectral_density
    if method is Kernel.spectral_density:
        return False
    if method is StationaryKernel.spectral_density:
        return leaf.SPECTRAL_NU is not None
    return True


def _nodes(kernel: Kernel) -> list[Kernel]:
    """Every node of the tree, root first — composites as well as leaves."""
    found = [kernel]
    for _, child in kernel.terms:
        found.extend(_nodes(child))
    return found


def check_spectral_support(kernel: Kernel, dimensions: int, *, owner: str) -> None:
    """Refuse, by name, a kernel tree with no closed-form spectral density.

    The same composition-time discipline ``QuasisepGP`` applies to the
    quasiseparable registry: the node that cannot be represented is named
    here, loudly, rather than discovered when a solve returns nonsense.

    **Every node, not every leaf**, which is what makes a
    :class:`~ampere.core.kernels.Product` of two Matérns a refusal: its leaves
    each have a spectral density and the product of the two kernels does not
    have theirs. A :class:`~ampere.core.kernels.Sum` passes and its children
    are then checked in turn, because linearity is exactly what a sum's own
    implementation relies on.
    """
    for node in _nodes(kernel):
        if not _has_spectral_density(node):
            raise LikelihoodError(
                f"{owner} needs a stationary kernel with a closed-form spectral density, but "
                f"{type(node).__name__} ({node.FAMILY}) has none. Matern12, Matern32, Matern52, "
                f"SquaredExponential and SHO have one, and so does any Sum of them; a Product's "
                f"is a convolution of its factors' with no closed form — and none at all when "
                f"the factors act on different axes — while RotationTerm and SpectralMixture "
                f"were left out of W5.4 deliberately. Use DenseGP, which needs no spectral "
                f"density at all."
            )
        if node.FAMILY == "sho" and int(dimensions) != 1:
            raise LikelihoodError(
                f"{owner} was given an SHO acting on {int(dimensions)} axes. A damped oscillator "
                f"is a process in one ordered coordinate and has no isotropic form on a plane; "
                f"select the single axis it runs along with axes=(...), or use a Matérn, which "
                f"is isotropic in any dimension."
            )
    selections = {leaf.axes for leaf in kernel.leaves()}
    if len(selections) > 1:
        raise LikelihoodError(
            f"{owner} builds one tensor-product basis over one box, so every term of the kernel "
            f"must act on the same axes; this one's terms select {sorted(map(str, selections))}. "
            f"Give every term the same axes=, or use DenseGP, which evaluates each term on its "
            f"own axes."
        )
