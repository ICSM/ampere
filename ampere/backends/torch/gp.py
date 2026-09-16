"""GP solvers on the torch path: the dense one and the quasiseparable one.

``DEVELOPMENT_PLAN.md`` §4.4 makes the solver a strategy chosen per problem, so
that the scaling story can change without the science code changing, and
``likelihoods.md`` §7 fixes the interface. :class:`DenseGP` here is the same
quantity ``ampere.core.DenseGP`` computes — the Gaussian marginal
log-likelihood of the residuals under ``K(θ) + diag(σ²)`` — through torch's own
Cholesky rather than scipy's.

:class:`Matern32` and :class:`SquaredExponential` are here too, **since
W2.13**. They were not in slice 1, and their absence was that slice's
principal carried finding: the solver built its covariance through
``ampere.core.Kernel.matrix``, which computes in numpy, so the Cholesky was
differentiable and the kernel hyperparameters were not. W2.5 had already
solved it on jax the same way — native kernels subclassing the core
declarations — and this is the mirror of that.

Why a torch dense solver is worth having at all
------------------------------------------------
It is not for speed: a dense ``O(N³)`` Cholesky is a dense ``O(N³)`` Cholesky
in either library, and ``architecture.md`` §2 gives the reference backend no
performance goal to beat. It is worth having for three other reasons.

* It is **differentiable**. ``torch.linalg.cholesky`` and
  ``cholesky_solve`` have derivatives, so the whole marginal likelihood does —
  which is what a torch-backed NUTS or a gradient optimiser over kernel
  hyperparameters needs, and what no amount of speed in the numpy path would
  give.
* It is the **second implementation** the conformance battery's cross-backend
  rows compare against. Two backends that agreed because they called the same
  LAPACK entry point through the same wrapper would prove nothing; these do not
  (scipy's ``cho_factor``/``cho_solve`` against torch's ``linalg.cholesky``/
  ``cholesky_solve``), and the agreement is asserted at
  ``tolerances.cross_backend``.
* It is where the **float64 policy** is not negotiable (``architecture.md``
  §5): a GP solve in float32 fails in ways that read as science bugs, so the
  inputs are converted to float64 on the way in regardless of what the caller
  was working in.

``QuasisepGP``, and the choice ``DEVELOPMENT_PLAN.md`` §6 deferred
-------------------------------------------------------------------
:class:`QuasisepGP` is W2.4 slice 2's, and it is the solver that makes this
backend's scaling story the project's own: exact O(N) for ordered 1D data,
differentiable in the kernel hyperparameters as well as in the residual.

§6 asked for two candidate libraries to be *measured* against this suite
rather than chosen from documentation — "GPyTorch structured solvers vs
celerite2's experimental torch interface" — and the measurement (2026-09-07;
the decision-log row carries the table) found that neither exists in the form
the question supposed. celerite2 0.3.3 ships no torch interface at all, but it
does ship :mod:`celerite2.backprop`, the compiled forward *and reverse* passes
its jax and pymc wrappers are built on; :mod:`ampere.backends.torch._celerite`
is the autograd shim over those, and it is the winner on every axis measured
(bit-identical to the reference solver, 123x faster than the alternative at
10^3 points and 5 657x at 10^4, zero new dependencies). GPyTorch ships no
quasiseparable operator whatsoever — its structured routes are Toeplitz (a
regular grid), Kronecker (a product grid) and SKI (approximate, and measured
at 1e-2 against a 1e-6 tolerance) — so on the irregular coordinates ampere
guarantees, its only *exact* route is a dense Cholesky, which is
:class:`DenseGP` under another name.

One finding from W2.3 shaped both solvers here: **celerite2 returns quiet NaN
where a Cholesky raises.** ``ampere.core.QuasisepGP`` guards its preconditions
(finite diagonal, positive amplitude) before calling celerite2, and this
backend does the same on its raising surfaces and sanitises into ``-inf`` on
its native one — a NaN log-likelihood has to surface through §4.5's
failure-signalling route (``-inf`` when non-strict, a raise under ``strict``),
never as a silent NaN in the chain.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.core import DTYPE, GPConditional, GPSolver, Kernel, Product, Sum, WarpedKernel
from ampere.core import SHO as _CoreSHO
from ampere.core import StationaryKernel as _CoreStationary
from ampere.core import Matern12 as _CoreMatern12
from ampere.core import Matern32 as _CoreMatern32
from ampere.core import Matern52 as _CoreMatern52
from ampere.core import RotationTerm as _CoreRotationTerm
from ampere.core import SpectralMixture as _CoreSpectralMixture
from ampere.core import SquaredExponential as _CoreSquaredExponential
from ampere.core.exceptions import LikelihoodError
from ampere.core.hsgp import (
    DEFAULT_BASIS_SIZE,
    DEFAULT_BOUNDARY_FACTOR,
    HilbertSpaceBasis,
    basis_matrix,
    basis_size,
    check_spectral_support,
    hilbert_basis,
    normalise_counts,
    spectral_values,
)
from ampere.core.kernels import lookup_quasiseparable_term

from . import _celerite
from ._config import (
    BACKEND,
    DEFAULT_DEVICE,
    DEFAULT_DTYPE,
    as_tensor,
    device_name,
    place,
    resolve_device,
    resolve_dtype,
    to_numpy,
)

__all__ = [
    "SHO",
    "DenseGP",
    "HilbertSpaceGP",
    "Matern12",
    "Matern32",
    "Matern52",
    "Product",
    "QuasisepGP",
    "RotationTerm",
    "SpectralMixture",
    "SquaredExponential",
    "Sum",
    "TorchOps",
    "WarpedKernel",
]

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)


def _points(
    coordinates: Any,
    *,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> torch.Tensor:
    """``(n, d)`` coordinates from an ``(n,)`` or ``(n, d)`` array or tensor.

    A bare one-dimensional input is read as a column of ``n`` one-dimensional
    points, matching ``ampere.core.Kernel.matrix``'s own convention.

    *dtype* and *device* are threaded rather than defaulted at the call site
    (W2.4 slice 3): a kernel or solver placed on a GPU handed coordinates that
    were already there must not have them silently copied back to the CPU,
    which is exactly what a hard-wired :data:`DEFAULT_DEVICE` here would do —
    once per evaluation, invisibly.
    """
    tensor = as_tensor(coordinates, dtype=dtype, device=device)
    return tensor[:, None] if tensor.ndim == 1 else tensor


def _separation(
    left: Any,
    right: Any,
    *,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> torch.Tensor:
    """Euclidean separation between two coordinate sets, ``(n, m)``.

    Euclidean in the coordinate space, which is why
    :meth:`~ampere.core.GPSolver.check_compatible` requires every coordinate
    axis to share one unit — a single isotropic length scale is meaningless
    across mixed units.
    """
    a = _points(left, dtype=dtype, device=device)
    b = _points(right, dtype=dtype, device=device)
    difference = a[:, None, :] - b[None, :, :]
    return torch.sqrt(torch.sum(difference * difference, dim=-1))


class TorchOps:
    """:class:`~ampere.core.ArrayOps` in torch, on one dtype and device.

    **W4.5.** Before W4.5 this module transcribed each kernel's closed form and
    each celerite representation into torch by hand, with a comment on every
    one saying it "must stay a transcription". Five more families would have
    made that fifteen copies. The mathematics now lives once in
    :mod:`ampere.core.kernels`, written against this protocol, and what remains
    backend-specific is exactly this class: where a number becomes a tensor,
    and on what device.

    One instance per kernel, not one per module, for the reason W2.4 slice 3
    gave: a kernel placed on a GPU and handed coordinates that are already
    there must not have them silently copied back, which a module-level default
    device would do once per evaluation, invisibly.
    """

    def __init__(self, dtype: torch.dtype, device: torch.device) -> None:
        self.dtype = dtype
        self.device = device

    def scalar(self, value: Any) -> torch.Tensor:
        return as_tensor(value, dtype=self.dtype, device=self.device)

    def points(self, coordinates: Any, *, dimensions: int | None = None) -> torch.Tensor:
        tensor = _points(coordinates, dtype=self.dtype, device=self.device)
        if dimensions is not None and int(tensor.shape[1]) != dimensions:
            raise LikelihoodError(
                f"kernel coordinates have {int(tensor.shape[1])} coordinate(s) per point, but "
                f"the data they are being compared against have {dimensions}."
            )
        return tensor

    def n_points(self, coordinates: Any) -> int:
        return int(self.points(coordinates).shape[0])

    def separation(self, left: Any, right: Any) -> torch.Tensor:
        difference = left[:, None, :] - right[None, :, :]
        return torch.sqrt(torch.sum(difference * difference, dim=-1))

    def zeros(self, n: int) -> torch.Tensor:
        return torch.zeros(n, dtype=self.dtype, device=self.device)

    def ones_like(self, array: Any) -> torch.Tensor:
        return torch.ones_like(self.scalar(array))

    def stack(self, arrays: Sequence[Any], axis: int = -1) -> torch.Tensor:
        return torch.stack([self.scalar(item) for item in arrays], dim=axis)

    def concatenate(self, arrays: Sequence[Any], axis: int = -1) -> torch.Tensor:
        return torch.cat([self.scalar(item) for item in arrays], dim=axis)

    def exp(self, array: Any) -> torch.Tensor:
        return torch.exp(self.scalar(array))

    def cos(self, array: Any) -> torch.Tensor:
        return torch.cos(self.scalar(array))

    def sin(self, array: Any) -> torch.Tensor:
        return torch.sin(self.scalar(array))

    def log1p(self, array: Any) -> torch.Tensor:
        return torch.log1p(self.scalar(array))

    def absolute(self, array: Any) -> torch.Tensor:
        return torch.abs(self.scalar(array))


class _TorchKernel(Kernel):
    """Shared plumbing for the torch kernels: torch separations, torch covariances.

    **W2.13.** W2.4 slice 1 shipped a torch :class:`DenseGP` that built its
    covariance by calling ``ampere.core.Kernel.matrix`` — which computed in
    numpy and coerced every hyperparameter with ``float()``. The Cholesky was
    differentiable and the *hyperparameters were not*: the graph was cut at the
    covariance, so a fitted amplitude or length scale received no gradient at
    all. That was slice 1's principal carried finding, and W2.5 had already
    solved it the same way on jax. These classes are the fix.

    They **subclass the core kernels** rather than redeclaring them, exactly as
    jax's do: same ``FAMILY``, same ``HYPERPARAMETERS``, same
    ``QUASISEPARABLE`` flag, so a problem's declaration — and therefore its
    spec hash — is unchanged by which backend computes it.

    **W4.5 removes the transcriptions.** ``matrix``, ``diagonal`` and every
    family's ``_covariance`` used to be overridden here to rebuild the same
    arithmetic in torch. They are not any more: the core's implementations are
    written against :class:`~ampere.core.ArrayOps`, and this class supplies
    :class:`TorchOps`. The one thing still declared here is that the core's
    positivity checks are *not* applied (``VALIDATES = False``), and its
    absence is not silent: a non-positive length scale gives a non-finite
    covariance, the Cholesky then fails, and :class:`DenseGP`'s two surfaces
    turn that into the failure §4.5 asks for (an exception on the contract
    path, ``-inf`` on the traced one). The declaration is what keeps it from
    arising: a kernel hyperparameter is declared with a prior on the positive
    half-line and a ``Log`` bijection, so no sampler proposes one.
    """

    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    #: Batchable since W2.4 slice 2: a kernel's covariance is elementwise
    #: arithmetic over a fixed separation matrix, which ``vmap`` maps over a
    #: stack of hyperparameters without any special handling.
    BATCHABLE: ClassVar[bool] = True
    #: The class-level default only; :meth:`_place` shadows it per instance
    #: from the ``device=`` keyword (W2.4 slice 3).
    DEVICE: ClassVar[str] = "cpu"
    #: The hyperparameters may be tracers; positivity is the declaration's job.
    VALIDATES: ClassVar[bool] = False

    def _place(self, dtype: Any, device: Any) -> None:
        """Record where this kernel's covariance is built. Called from ``__init__``.

        A kernel holds no buffers -- its hyperparameters are
        :class:`~ampere.core.Parameter` objects, which are neutral
        declarations -- so there is nothing to move and no ``to()``. What it
        does hold is a *place to build in*: every tensor :meth:`matrix` and
        :meth:`diagonal` create is created there, and the separation matrix is
        the largest array in an ordinary dense GP evaluation.
        """
        place(self, dtype, device)
        self._ops = TorchOps(self.dtype, self.device)

    def _scalar(self, value: Any) -> torch.Tensor:
        return as_tensor(value, dtype=self.dtype, device=self.device)


class _TorchStationary(_TorchKernel, _CoreStationary):
    """The two-hyperparameter torch kernels' constructor."""

    def __init__(
        self,
        amplitude: Any,
        length_scale: Any,
        *,
        amplitude_unit: Any = None,
        length_scale_unit: Any = None,
        axes: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            amplitude,
            length_scale,
            amplitude_unit=amplitude_unit,
            length_scale_unit=length_scale_unit,
            axes=axes,
        )
        self._place(dtype, device)


class Matern12(_TorchStationary, _CoreMatern12):
    r"""Matérn-1/2 (Ornstein-Uhlenbeck) in torch. Rank 1 on the O(N) path."""


class Matern32(_TorchStationary, _CoreMatern32):
    r"""Matérn-3/2 in torch: ampere's canonical flexible-likelihood kernel.

    .. math::
        k(r) = a^2 \left(1 + \frac{\sqrt{3}\,r}{\ell}\right)
               \exp\!\left(-\frac{\sqrt{3}\,r}{\ell}\right)

    ``amplitude`` is the marginal **standard deviation** — ``k(0) ==
    amplitude²`` — so a prior on it is a prior in the data's own units. The
    declaration, ``QUASISEPARABLE = True`` included, is inherited from
    ``ampere.core.Matern32``: this kernel does have an exact rank-2
    quasiseparable representation, and saying otherwise here would be a lie
    about the mathematics merely because this backend has not yet shipped a
    solver that exploits it.
    """


class Matern52(_TorchStationary, _CoreMatern52):
    r"""Matérn-5/2 in torch. Rank 3 on the O(N) path, and exactly so (W4.5)."""


class SquaredExponential(_TorchStationary, _CoreSquaredExponential):
    r"""Squared exponential (RBF) in torch: legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Not quasiseparable, and inherits that declaration too, so a quasiseparable
    solver refuses it by name — the asymmetry that made Matérn the default.
    """


class SHO(_TorchKernel, _CoreSHO):
    r"""A damped harmonic oscillator in torch: the fringing component (W4.5)."""

    def __init__(
        self,
        amplitude: Any,
        period: Any,
        quality: Any,
        *,
        amplitude_unit: Any = None,
        period_unit: Any = None,
        axes: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            amplitude,
            period,
            quality,
            amplitude_unit=amplitude_unit,
            period_unit=period_unit,
            axes=axes,
        )
        self._place(dtype, device)


class RotationTerm(_TorchKernel, _CoreRotationTerm):
    r"""celerite2's rotation pair in torch: two SHOs at a period and its harmonic."""

    def __init__(
        self,
        amplitude: Any,
        period: Any,
        quality: Any,
        delta_quality: Any = 0.0,
        fraction: Any = 0.5,
        *,
        amplitude_unit: Any = None,
        period_unit: Any = None,
        axes: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            amplitude,
            period,
            quality,
            delta_quality,
            fraction,
            amplitude_unit=amplitude_unit,
            period_unit=period_unit,
            axes=axes,
        )
        self._place(dtype, device)


class SpectralMixture(_CoreSpectralMixture):
    """A sum of torch :class:`SHO` components (W4.5).

    Not a ``_TorchKernel``: a composite has no arithmetic of its own — it
    combines its children's matrices — so it takes its namespace, its device
    and its capability flags from them, which is what
    ``ampere.core.kernels._Composite`` already does. All this subclass supplies
    is the default component type.
    """

    def __init__(
        self,
        amplitudes: Sequence[Any],
        periods: Sequence[Any],
        qualities: Sequence[Any],
        *,
        component_type: type = SHO,
        **kwargs: Any,
    ) -> None:
        super().__init__(amplitudes, periods, qualities, component_type=component_type, **kwargs)


def _right_hand_side(tensor: torch.Tensor) -> torch.Tensor:
    """A solver's right-hand side, flattened only when it is genuinely 1-D (*W4.2*).

    Every solve in this module used to begin ``.reshape(-1)``, which is right
    for the one realisation a real residual is and catastrophic for the ``(n, k)``
    block ``likelihoods.md`` §7 now permits: flattening two columns would score
    ``2n`` residuals against an ``n x n`` covariance. One helper rather than the
    condition written out six times, because the six places have to agree.
    """
    return tensor.reshape(-1) if tensor.ndim < 2 else tensor


def _columns(tensor: torch.Tensor) -> int:
    """How many realisations a right-hand side carries: ``k``."""
    return 1 if tensor.ndim < 2 else int(tensor.shape[1])


@dataclasses.dataclass(frozen=True)
class DenseGP(GPSolver):
    """Exact dense Cholesky, in torch. Agrees with ``ampere.core.DenseGP`` exactly.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal. Defaults to zero, for the same reason the reference
        solver does: a covariance that will not factorise is a fact about the
        model, and this contract does not hide it behind a silent numerical
        fudge. The error message says when setting it is the right answer.

    Notes
    -----
    The solver interface (``likelihoods.md`` §7) is numpy in and numpy out,
    because it is called from ``ampere.core``'s backend-neutral likelihood.
    :meth:`log_marginal_likelihood_tensor` is the same computation without the
    conversion, for a caller that wants the gradient;
    :meth:`log_marginal_likelihood_native` is that one again without exception
    control flow, which is what a realised density needs.

    **``jitter`` is the only dataclass field, and that is a constraint rather
    than a simplification.** ``Likelihood.to_spec`` records a dataclass
    solver's fields as its ``config``, and ``results.md`` §14 requires
    ``ampere_spec_hash`` to agree across backends — so a field this solver had
    and ``ampere.core.DenseGP`` did not would make two backends' *declaration*
    of the same problem differ, which is exactly the cross-backend
    disagreement §14 calls "the cheapest possible detector" of a lowering bug.
    Precision and device are therefore :class:`~typing.ClassVar` policy
    (:data:`TENSOR_DTYPE`, :data:`TENSOR_DEVICE`), not configuration: they are
    how this backend computes, not part of what the user declared.
    **W2.13 ruled the question slice 1 left open** (``inference.md`` §10a,
    fold-in 10): they stay out of the spec and are reported from
    :meth:`provenance_config`, which a run records under
    ``ampere_solver_config`` and no hash reads. A float32 opt-out — which
    ``architecture.md`` §5 requires to be visible in provenance — has a home
    now without disturbing cross-backend spec-hash agreement.
    """

    jitter: float = 0.0

    #: Where the solve happens. ``architecture.md`` §5's float64 policy for
    #: likelihood and GP linear algebra, on the CPU by default (chosen, never
    #: detected). Deliberately class-level; see the class docstring.
    TENSOR_DTYPE: ClassVar[torch.dtype] = DEFAULT_DTYPE
    TENSOR_DEVICE: ClassVar[torch.device] = DEFAULT_DEVICE

    NAME: ClassVar[str] = "DenseGP"
    EXACT: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    #: W4.2: this solver takes the circular complex GP's two-column right-hand
    #: side, exactly as ``ampere.core.DenseGP`` does.
    STACKED_RESIDUALS: ClassVar[bool] = True

    #: The four capability flags. W2.4 declared them against a parts set that
    #: did not yet include a solver and recorded the gap; **W2.13 closed it**
    #: (``inference.md`` §10a, fold-in 7), so these are now what makes a
    #: composed problem report ``backend="torch"`` — and what makes the same
    #: problem left with ``ampere.core.DenseGP`` a refusal rather than a
    #: silently non-differentiable GP.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: Batchable since W2.4 slice 2. ``torch.linalg.cholesky_ex`` and
    #: ``cholesky_solve`` are batched operations in torch, so ``vmap`` maps a
    #: stack of covariances onto a stack of factorisations natively — which is
    #: the one place batching actually buys something, because a dense solve is
    #: where the arithmetic is.
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __post_init__(self) -> None:
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(f"DenseGP's jitter must be finite and >= 0, got {self.jitter!r}.")

    # -- factorisation ---------------------------------------------------------

    def _tensor(self, value: Any) -> torch.Tensor:
        return as_tensor(value, dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)

    def _covariance(
        self,
        kernel: Kernel,
        coordinates: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """``K(θ) + diag(variance + jitter²)``, with every graph intact.

        Both arguments may carry a gradient now: the covariance because a
        native kernel (:class:`Matern32` here, not ``ampere.core``'s) builds
        it in torch, and the diagonal because a prediction-aware noise model's
        sigma depends on θ. ``as_tensor`` is used for both rather than
        ``np.asarray``, which would have detached them — the numpy coercion
        this line used to do is precisely why W2.4's GP hyperparameters got no
        gradient.
        """
        matrix = self._tensor(kernel.matrix(coordinates, coordinates, values))
        diagonal = self._tensor(variance) + self.jitter**2
        return matrix + torch.diag(diagonal)

    def _cholesky(self, total: torch.Tensor) -> torch.Tensor:
        """The lower Cholesky factor, or the reference solver's own two refusals.

        The messages are deliberately the reference backend's, word for word:
        a user who moves a failing problem from one backend to the other should
        get the same diagnosis, and the causes (a zero uncertainty, coordinates
        closer than float64 can separate, an amplitude far above the data
        scale) are properties of the problem rather than of the library.
        """
        if not bool(torch.all(torch.isfinite(total))):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        factor, info = torch.linalg.cholesky_ex(total)
        if int(info) != 0:
            raise LikelihoodError(
                f"the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                f"Cholesky factorisation failed (leading minor {int(info)} is not positive "
                f"definite). Usual causes: a zero or near-duplicate observational uncertainty, "
                f"coordinates that are closer together than float64 can separate at this "
                f"length-scale, or a kernel amplitude far above the data scale. Pass "
                f"DenseGP(jitter=...) — a standard deviation in the data's units — if the matrix "
                f"is merely ill-conditioned rather than wrong."
            )
        return factor

    @staticmethod
    def _solve(factor: torch.Tensor, right: torch.Tensor) -> torch.Tensor:
        """``(L Lᵀ)⁻¹ right``, with *right* accepted as a vector or a matrix."""
        vector = right.ndim == 1
        columns = right.reshape(-1, 1) if vector else right
        solved = torch.cholesky_solve(columns, factor, upper=False)
        return solved.reshape(-1) if vector else solved

    # -- the GPSolver interface ------------------------------------------------

    def log_marginal_likelihood_tensor(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: Any,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """``log N(residual; 0, K + diag(variance))`` as a differentiable tensor.

        **W4.2**: an ``(n, k)`` residual is the joint marginal of ``k``
        realisations sharing this covariance — one factorisation, ``k`` solves,
        the log-determinant counted ``k`` times. ``(residuals * alpha).sum()``
        rather than ``residuals @ alpha`` because the latter is not defined
        between two blocks and is the sum of the quadratic forms for a vector.
        """
        total = self._covariance(kernel, coordinates, variance, values)
        factor = self._cholesky(total)
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(factor, residuals)
        log_determinant = 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = (residuals * alpha).sum()
        return -0.5 * (
            quadratic + _columns(residuals) * log_determinant + residuals.numel() * _LOG_2PI
        )

    def log_marginal_likelihood_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """The same quantity, **without exception control flow** (W2.13).

        The surface :mod:`ampere.backends.torch.problem` composes, and the
        counterpart of jax's ``log_marginal_likelihood_jax``. A realised
        density must not raise — ``inference.md`` §10a narrows §11 for exactly
        this reason — so a factorisation that fails becomes ``-inf`` through
        :func:`torch.where` rather than a :class:`LikelihoodError`.

        Note that :func:`torch.linalg.cholesky_ex` is what makes this possible
        at all: the plain ``cholesky`` raises, and catching the exception would
        put Python control flow back in the middle of the hot loop. ``info``
        is read as a value.
        """
        total = self._covariance(kernel, coordinates, variance, values)
        factor, info = torch.linalg.cholesky_ex(total)
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(factor, residuals)
        log_determinant = 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = (residuals * alpha).sum()
        value = -0.5 * (
            quadratic + _columns(residuals) * log_determinant + residuals.numel() * _LOG_2PI
        )
        failed = torch.logical_or(info != 0, torch.logical_not(torch.isfinite(value)))
        return torch.where(failed, torch.full_like(value, -math.inf), value)

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        return float(
            self.log_marginal_likelihood_tensor(kernel, coordinates, residual, variance, values)
        )

    def provenance_config(self) -> Mapping[str, Any]:
        """``inference.md`` §10a fold-in 10: how this solver computes, for the attrs.

        The precision and device this backend factorises in. They are
        deliberately **not** dataclass fields — ``Likelihood.to_spec`` records
        those, and ``results.md`` §14 requires the spec hash to agree across
        backends, so a dtype in the spec would make two backends' declaration
        of one problem differ. They are still worth recording, because
        ``architecture.md`` §5 requires a float32 opt-out to be visible in a
        run's provenance, and :meth:`configured` is that opt-out.
        """
        return {"dtype": str(self.TENSOR_DTYPE), "device": str(self.TENSOR_DEVICE)}

    def configured(self, *, dtype: Any = None, device: Any = None) -> DenseGP:
        """A copy of this solver that computes in a different precision or place.

        ``architecture.md`` §5's **per-run opt-out** from the float64 policy,
        and the ``device=`` opt-in, in the one form that does not disturb
        anything else (W2.4 slice 2).

        Why a copy with shadowed class attributes rather than constructor
        arguments: ``Likelihood.to_spec`` records a dataclass solver's
        ``dataclasses.fields`` as its ``config``, and ``results.md`` §14
        requires ``ampere_spec_hash`` to agree across backends — so a ``dtype``
        field here would make the torch *declaration* of a problem differ from
        the reference one, and the cross-backend hash row would fail on a
        difference that is not a difference in the model. Set as instance
        attributes, they shadow the :class:`~typing.ClassVar` policy for this
        solver only, ``dataclasses.fields`` still reports ``jitter`` alone, and
        :meth:`provenance_config` reports what actually happened. Visible,
        never hashed.

        **float32 is a real choice with a real cost**, which is why it is an
        opt-out rather than an option. ``DEVELOPMENT_PLAN.md`` §7 lists it
        among the known traps: "GP Cholesky / quasiseparable solves in float32
        fail in ways that look like science problems". A float32 Cholesky of a
        well-conditioned 10³-point Matérn covariance loses roughly seven digits
        relative to float64, which is larger than every tolerance in the
        conformance table; the reason to take it anyway is GPU throughput,
        where the memory bandwidth saved is the whole point. Nothing here
        stops you; the run says what you did.

        Parameters
        ----------
        dtype
            A ``torch.dtype``, or a name torch understands (``"float32"``).
            ``None`` leaves the policy alone.
        device
            A ``torch.device``, or anything ``torch.device`` accepts
            (``"cpu"``, ``"cuda:0"``). ``None`` leaves it alone. **Never
            auto-detected** (``architecture.md`` §5): a machine with a GPU
            present takes the same path as CI unless a caller says otherwise.

        Notes
        -----
        **The device moves the capability flag with it** (W2.4 slice 3). Until
        slice 3 this method shadowed :data:`TENSOR_DEVICE` alone, so a solver
        configured onto a GPU went on *declaring* ``DEVICE = "cpu"`` and
        composed silently into a problem whose every other part really was on
        the CPU — precisely the mixed-device problem
        :func:`~ampere.core.declared_capabilities` exists to refuse, made
        invisible to it. ``DEVICE`` is now shadowed beside it, so composing a
        GPU solver with CPU models is a composition-time refusal naming both
        devices.

        Examples
        --------
        >>> import torch
        >>> from ampere.backends.torch import DenseGP
        >>> fast = DenseGP(jitter=1e-6).configured(dtype=torch.float32)
        >>> fast.provenance_config()["dtype"]
        'torch.float32'
        >>> fast.jitter
        1e-06

        and the declaration is untouched, which is the point:

        >>> import dataclasses
        >>> [f.name for f in dataclasses.fields(fast)]
        ['jitter']

        while the device it says it is on is the device it was put on:

        >>> DenseGP().configured(device="cpu").DEVICE
        'cpu'
        """
        return _configured(self, dtype, device)

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """The closed form from the same Cholesky (Sundararajan & Keerthi 2001).

        With ``A = (K + diag(variance + jitter²))⁻¹``:
        ``sigma_i^{2,-i} = 1 / A_ii`` and ``mu_i^{-i} = y_i - [A r]_i / A_ii``,
        so ``log p_i = 0.5 log A_ii - [A r]_i² / (2 A_ii) - 0.5 log(2π)``. The
        same identity the reference solver uses, so the two agree exactly on
        the same matrix.
        """
        total = self._covariance(kernel, coordinates, variance, values)
        factor = self._cholesky(total)
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(factor, residuals)
        rows = int(residuals.shape[0])
        columns = _columns(residuals)
        identity = torch.eye(rows, dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)
        precision_diagonal = torch.diagonal(self._solve(factor, identity))
        # W4.2: one term per *sample*, summing the k components, which share
        # A_ii -- the quantity a leave-one-out conditional is built from.
        quadratic = (alpha.reshape(rows, columns) ** 2).sum(dim=1)
        terms = columns * (0.5 * torch.log(precision_diagonal) - 0.5 * _LOG_2PI) - quadratic / (
            2.0 * precision_diagonal
        )
        return to_numpy(terms).astype(DTYPE, copy=False)

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        points = _as_points(coordinates, "data coordinates")
        total = self._covariance(kernel, points, variance, values)
        factor = self._cholesky(total)
        target = points if at is None else _as_points(at, "conditioning grid", points.shape[1])
        cross = self._tensor(kernel.matrix(target, points, values))
        residuals = _right_hand_side(self._tensor(residual))
        # An (n, k) residual conditions to an (m, k) mean beside one (m,)
        # variance: a posterior variance does not depend on the data (W4.2).
        mean = cross @ self._solve(factor, residuals)
        solved = self._solve(factor, cross.transpose(0, 1))
        prior_variance = self._tensor(kernel.diagonal(target, values))
        posterior = prior_variance - (cross * solved.transpose(0, 1)).sum(dim=1)
        return GPConditional(
            mean=to_numpy(mean).astype(DTYPE, copy=False),
            variance=to_numpy(posterior).astype(DTYPE, copy=False),
        )

    def latent_transform_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> torch.Tensor:
        """``f = L(θ) z`` as a **differentiable tensor**, without exceptions (W2.14).

        The counterpart of :meth:`log_marginal_likelihood_native`, and needed
        for the same two reasons. :meth:`latent_transform` returns numpy, so a
        realised density that called it would cut the gradient in exactly the
        hyperparameters the latent path exists to fit — W2.4 slice 1's finding
        arriving a second time, on a second surface. And a realised density
        must not raise (``inference.md`` §10a), so a factorisation that fails
        comes back as NaN through :func:`torch.where` and becomes ``-inf`` at
        the one guard the lowered dataset ends with, rather than as a
        :class:`LikelihoodError` in the middle of a NUTS step.

        ``torch.linalg.cholesky_ex`` is again what makes it possible: ``info``
        is read as a value rather than caught as an exception.
        """
        covariance = self._tensor(kernel.matrix(coordinates, coordinates, values))
        # The reference solver's ``float(...) or 1.0``, kept as a value so a
        # fitted amplitude survives it.
        mean_variance = torch.diagonal(covariance).mean()
        scale = torch.where(mean_variance != 0.0, mean_variance, torch.ones_like(mean_variance))
        stabilised = covariance + torch.eye(
            covariance.shape[0], dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE
        ) * (jitter * scale)
        factor, info = torch.linalg.cholesky_ex(stabilised)
        drawn = factor @ _right_hand_side(self._tensor(whitened))
        failed = torch.logical_or(info != 0, torch.logical_not(torch.isfinite(drawn).all()))
        return torch.where(failed, torch.full_like(drawn, math.nan), drawn)

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """``f = L z`` — the deterministic half of the latent-GP declaration.

        The contract surface, so it raises where
        :meth:`latent_transform_native` returns NaN.
        """
        points = _as_points(coordinates, "data coordinates")
        covariance = self._tensor(kernel.matrix(points, points, values))
        if not bool(torch.all(torch.isfinite(covariance))):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        drawn = self.latent_transform_native(kernel, points, whitened, values, jitter=jitter)
        if not bool(torch.all(torch.isfinite(drawn))):
            raise LikelihoodError(
                "the kernel matrix K is not positive definite, so the whitening transform "
                "f = L z is undefined. Increase the jitter argument, or check the kernel "
                "hyperparameters."
            )
        return to_numpy(drawn).astype(DTYPE, copy=False)


def _configured(solver: Any, dtype: Any, device: Any) -> Any:
    """A copy of *solver* whose dtype/device shadow the class policy.

    Shared by :meth:`DenseGP.configured` so that "a copy, with instance
    attributes, leaving ``dataclasses.fields`` alone" is written once. The
    dataclass is frozen, so the copy is built by ``dataclasses.replace`` and
    the shadows are set through ``object.__setattr__`` — the same route
    ``__post_init__`` would take.
    """
    copy = dataclasses.replace(solver)
    if dtype is not None:
        try:
            resolved = resolve_dtype(dtype)
        except TypeError as error:
            raise LikelihoodError(
                f"{type(solver).__name__}.configured(dtype=...) takes a torch.dtype or the name "
                f"of one, got {dtype!r}."
            ) from error
        if not resolved.is_floating_point:
            raise LikelihoodError(
                f"{type(solver).__name__}.configured(dtype={dtype!r}) was given a "
                f"non-floating-point dtype. A covariance is a real matrix; an integer one "
                f"cannot be factorised."
            )
        object.__setattr__(copy, "TENSOR_DTYPE", resolved)
    if device is not None:
        placed = resolve_device(device)
        object.__setattr__(copy, "TENSOR_DEVICE", placed)
        # And the *capability* flag beside it (W2.4 slice 3). The two were
        # separate until slice 3, which meant a solver configured onto a GPU
        # went on declaring ``DEVICE = "cpu"`` and composed happily into a
        # problem every other part of which was on the CPU -- a device
        # disagreement that ``declared_capabilities`` exists to refuse and
        # could not see. They move together now.
        object.__setattr__(copy, "DEVICE", device_name(placed))
    return copy


def _as_points(coordinates: Any, what: str, dimensions: int | None = None) -> np.ndarray:
    """Coordinates as an ``(n, d)`` array, the shape ``Kernel.matrix`` expects."""
    array = np.asarray(coordinates, dtype=DTYPE)
    if array.ndim == 1:
        array = array.reshape(-1, 1)
    if array.ndim != 2:
        raise LikelihoodError(f"{what} must be one- or two-dimensional, got shape {array.shape}.")
    if dimensions is not None and array.shape[1] != dimensions:
        raise LikelihoodError(
            f"{what} has {array.shape[1]} coordinate dimension(s) but the data have {dimensions}."
        )
    return array


# ---------------------------------------------------------------------------
# The quasiseparable solver
# ---------------------------------------------------------------------------


def _celerite_matrices(
    kernel: Kernel,
    values: Mapping[str, Any],
    axis: torch.Tensor,
    diagonal: torch.Tensor,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
    r"""celerite2's ``(c, a, U, V)`` for any registered kernel, in torch.

    **W4.5 replaced a transcription with a lookup.** Until W4.5 this was
    ``_matern32_matrices``, a hand-written torch copy of the core's rank-2
    Matérn-3/2 algebra whose docstring said in as many words that it "must stay
    a transcription", beside a private ``_QUASISEPARABLE_TERMS`` table whose
    comment said it "must carry exactly the same families" as the core's —
    with nothing able to check either claim. Both obligations are now
    structural rather than editorial: the representation comes from
    :func:`ampere.core.kernels.lookup_quasiseparable_term`, the one registry
    all three backends read, and the generators are built by the *same*
    function the reference path uses, in :class:`TorchOps` instead of numpy.
    A family this backend could lower and the reference could not is no longer
    expressible.

    The conversion from :class:`~ampere.core.CeleriteRepresentation` to
    celerite2's argument order is all that is left: ``a`` is the *noise*
    diagonal plus the kernel's marginal variance.
    """
    # celerite2's compiled kernels are float64 CPU whatever device the kernel
    # was placed on, and a reference (numpy) kernel handed to this solver must
    # still build tensors: the solver, not the kernel, decides the namespace.
    kernel = kernel.with_ops(TorchOps(axis.dtype, axis.device))
    builder = lookup_quasiseparable_term(kernel.FAMILY, owner="QuasisepGP")
    representation = builder(kernel, values, axis)
    return (
        representation.decay,
        diagonal + representation.marginal,
        representation.left,
        representation.right,
    )


@dataclasses.dataclass(frozen=True)
class QuasisepGP(GPSolver):
    """Exact O(N) for ordered 1D data, in torch, differentiable in the hyperparameters.

    The same quantity :class:`~ampere.core.QuasisepGP` computes — the Gaussian
    marginal log-likelihood of the residuals under ``K(θ) + diag(σ²)``, through
    celerite2's semiseparable factorisation over ampere's own **exact** rank-2
    Matérn-3/2 representation — with two differences that are the whole point
    of it existing:

    * it is **differentiable**, in the residual *and* in the kernel
      hyperparameters, because :mod:`ampere.backends.torch._celerite` wraps
      celerite2's compiled forward **and reverse** passes as
      :class:`torch.autograd.Function`\\ s. The numpy solver has no reverse
      pass at all, so a torch problem left with ``ampere.core.QuasisepGP``
      would not merely be slow — it would have no gradient, which is why
      :func:`~ampere.core.declared_capabilities` refuses that composition by
      name;
    * it is the **second implementation** of the recursion the conformance
      suite's ``DenseGP``↔``QuasisepGP`` rows compare, on a different
      autodiff path.

    Which library, and why this one
    -------------------------------
    ``DEVELOPMENT_PLAN.md`` §6 deferred the choice to this track — "GPyTorch
    structured solvers vs celerite2's experimental torch interface" — and W2.4
    slice 2 measured both (the decision-log row of 2026-09-07 has the table).
    Neither existed in the form the plan supposed. celerite2 ships no torch
    interface, but it does ship the compiled kernels its jax and pymc wrappers
    are built on, so ampere writes the shim; GPyTorch ships no quasiseparable
    operator at all, and on irregular 1D coordinates its only *exact* route is
    a dense Cholesky. See :mod:`ampere.backends.torch._celerite`.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal — the same knob, meaning and default as
        :class:`DenseGP`'s and as ``ampere.core.QuasisepGP``'s. Zero by
        default: a covariance that will not factorise is a fact about the
        model, not something to hide.

    Notes
    -----
    **Three surfaces, as** :class:`DenseGP` **has three.**
    :meth:`log_marginal_likelihood` is the numpy-in/numpy-out contract surface
    and *raises*; :meth:`log_marginal_likelihood_tensor` is the same
    computation differentiable and raising;
    :meth:`log_marginal_likelihood_native` is differentiable and **never
    raises**, which is what a realised density needs (``inference.md`` §10a).

    **celerite2 returns quiet NaN where a Cholesky raises** (W2.3's carried
    finding for both Phase-2 tracks). Every precondition — a finite,
    non-negative diagonal, a finite non-negative amplitude, a positive length
    scale — is therefore checked *before* the compiled kernel is called on the
    raising surfaces, and *sanitised into a* ``-inf`` on the native one. A NaN
    that gets through anyway is still converted, so the failure cannot reach a
    sampler as a silent NaN.

    **CPU float64 only, and that is declared rather than discovered.** The
    compiled kernels are double-precision CPU C++, so this solver declares
    ``BATCHABLE = False`` and ``DEVICE = "cpu"``, and
    :meth:`provenance_config` records the dtype and device every run
    factorised in. :class:`DenseGP` is the solver to reach for on another
    device or in reduced precision.
    """

    jitter: float = 0.0

    #: Where the solve happens. Class-level policy, never dataclass fields, for
    #: the reason :class:`DenseGP` states at length: ``Likelihood.to_spec``
    #: records a dataclass solver's fields as its ``config`` and ``results.md``
    #: §14 requires the spec hash to agree across backends, so a dtype in the
    #: spec would make two backends' *declaration* of one problem differ.
    #: Fixed rather than merely defaulted here: celerite2's compiled kernels
    #: are float64 CPU, so these are facts about the library.
    TENSOR_DTYPE: ClassVar[torch.dtype] = torch.float64
    TENSOR_DEVICE: ClassVar[torch.device] = torch.device("cpu")

    NAME: ClassVar[str] = "QuasisepGP"
    EXACT: ClassVar[bool] = True
    REQUIRES_ORDERED_1D: ClassVar[bool] = True
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True

    DIFFERENTIABLE: ClassVar[bool] = True
    #: **Not** batchable, and this is the measured price of the library choice
    #: rather than an omission: the solve happens inside ``celerite2.backprop``,
    #: a compiled extension reached through a :class:`torch.autograd.Function`
    #: that converts to numpy, and ``torch.func.vmap`` cannot see through
    #: either. A problem carrying this solver reports ``batchable=False`` and
    #: :meth:`ampere.backends.torch.LoweredProblem.log_prob_unconstrained_batched`
    #: refuses it by name — which is the honest outcome, since evaluating a
    #: stack one member at a time under a batched name would be a lie about the
    #: cost.
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __post_init__(self) -> None:
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"QuasisepGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )

    def check_compatible(self, kernel: Kernel, observed: Any) -> None:
        super().check_compatible(kernel, observed)
        # Every family in the tree, not just the root: a Sum lowers term by
        # term, so one unregistered term stops it and must be named here.
        for leaf in kernel.leaves():
            lookup_quasiseparable_term(leaf.FAMILY, owner=self.NAME)

    # -- internals -----------------------------------------------------------

    def _tensor(self, value: Any) -> torch.Tensor:
        return as_tensor(value, dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)

    def _axis(
        self, coordinates: Any, kernel: Kernel | None = None
    ) -> tuple[torch.Tensor, torch.Tensor]:
        """The bare coordinate axis and the permutation that sorts it.

        Coordinates need not arrive sorted: a Gaussian density is invariant
        under a simultaneous permutation of residuals, variances and
        coordinates, so this solver sorts internally and undoes the
        permutation on the way out — exactly as ``ampere.core.QuasisepGP``
        does, so the two agree on unsorted input as well as on sorted.

        ``kernel`` is passed since W4.5 so an axis-selecting kernel gets its
        own column out of a multi-axis container's coordinates.
        """
        points = _points(coordinates).to(dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)
        if kernel is not None:
            points = kernel.select(points)
        if points.shape[1] != 1:
            raise LikelihoodError(
                f"{self.NAME} needs one ordered coordinate per sample, but the coordinates have "
                f"{points.shape[1]} per point. check_compatible refuses this at composition "
                f"time; a direct solver call reaches it here. Use DenseGP for 2D+ coordinates."
            )
        axis = points[:, 0]
        return axis, torch.argsort(axis, stable=True)

    def _matrices(
        self,
        kernel: Kernel,
        axis: torch.Tensor,
        diagonal: torch.Tensor,
        values: Mapping[str, Any],
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        return _celerite_matrices(kernel, values, axis, diagonal)

    def _guarded_factor(
        self,
        kernel: Kernel,
        axis: torch.Tensor,
        diagonal: torch.Tensor,
        values: Mapping[str, Any],
        *,
        whitening: bool = False,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """``(c, U, d, W)`` for the sorted problem, **raising** on a bad input.

        The preconditions are checked here rather than left to celerite2
        because celerite2 does not check them: it returns quiet NaN for a
        diagonal that cannot belong to a covariance, where a Cholesky raises
        (W2.3, recorded for both Phase-2 tracks). The messages are
        :class:`DenseGP`'s and the reference solver's, word for word — a user
        who moves a failing problem between solvers or backends should get one
        diagnosis, because the causes are properties of the problem.
        """
        if not bool(torch.all(torch.isfinite(diagonal))):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        if bool(torch.any(diagonal < 0.0)):
            raise LikelihoodError(
                "the diagonal handed to QuasisepGP contains negative entries, so K + "
                "diag(sigma^2) is not a covariance matrix at all. Variances are squares; this "
                "is a caller error rather than an ill-conditioned problem."
            )
        # **W5.7.** ``axis`` is the *raw* sorted coordinate, which is what the
        # registry's builders are documented to receive; the recursion runs on
        # ``kernel.warped_coordinate(axis, values)``, which is the identity for
        # every kernel but a warped one. celerite2 forms its propagators
        # ``exp(-c (t_n - t_m))`` from the coordinate handed to ``factor``, so
        # a kernel stationary in ``w(x)`` cannot express itself through the
        # generators alone.
        solve_axis = kernel.warped_coordinate(axis, values)
        c, a, U, V = self._matrices(kernel, axis, diagonal, values)
        finite = (
            bool(torch.all(torch.isfinite(a)))
            and bool(torch.all(torch.isfinite(U)))
            and bool(torch.all(torch.isfinite(V)))
            and bool(torch.all(torch.isfinite(c)))
        )
        if not finite:
            raise LikelihoodError(
                "the kernel's quasiseparable generators contain non-finite entries, so the "
                "representation cannot be built. The usual cause is a kernel amplitude large "
                "enough that amplitude**2 overflows float64; constrain the amplitude prior to "
                "the data's own scale."
            )
        d, W = _celerite.factor(solve_axis, c, a, U, V)
        if not bool(torch.all(torch.isfinite(d))) or bool(torch.any(d <= 0.0)):
            if whitening:
                raise LikelihoodError(
                    "the kernel matrix K is not positive definite in its quasiseparable "
                    "representation, so the whitening transform f = L z is undefined. Increase "
                    "the jitter argument, or check the kernel hyperparameters."
                )
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                "quasiseparable factorisation failed. Usual causes: a zero or near-duplicate "
                "observational uncertainty, coordinates that are closer together than float64 "
                "can separate at this length-scale, or a kernel amplitude far above the data "
                "scale. Pass QuasisepGP(jitter=...) — a standard deviation in the data's units "
                "— if the matrix is merely ill-conditioned rather than wrong."
            )
        return c, U, d, W

    # -- the marginal likelihood ---------------------------------------------

    @staticmethod
    def _marginal(
        axis: torch.Tensor,
        c: torch.Tensor,
        U: torch.Tensor,
        d: torch.Tensor,
        W: torch.Tensor,
        residual: torch.Tensor,
    ) -> torch.Tensor:
        """``-0.5 (rᵀ A r + log|K + diag| + N log 2π)`` from a finished factorisation.

        ``L⁻¹ r`` is all that is needed for the quadratic form: with
        ``K + diag = L D Lᵀ`` and ``z = L⁻¹ r``, ``rᵀ A r = Σ z_i²/d_i`` and
        ``log|K + diag| = Σ log d_i``, because ``L`` is *unit* triangular.
        celerite2's own ``_do_norm`` is the same two lines.
        """
        z = _celerite.solve_lower(axis, c, U, W, residual.reshape(-1, 1))[:, 0]
        return -0.5 * ((z * z / d).sum() + torch.log(d).sum() + float(residual.numel()) * _LOG_2PI)

    def log_marginal_likelihood_tensor(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """``log N(residual; 0, K + diag(variance))`` as a differentiable tensor."""
        axis, order = self._axis(coordinates, kernel)
        residuals = self._tensor(residual).reshape(-1)
        diagonal = self._tensor(variance).reshape(-1) + self.jitter**2
        sorted_axis = axis[order]
        c, U, d, W = self._guarded_factor(kernel, sorted_axis, diagonal[order], values)
        return self._marginal(
            kernel.warped_coordinate(sorted_axis, values), c, U, d, W, residuals[order]
        )

    def log_marginal_likelihood_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """The same quantity, **without exception control flow** (``inference.md`` §10a).

        The surface :mod:`ampere.backends.torch.problem` composes, and the
        counterpart of :meth:`DenseGP.log_marginal_likelihood_native`. Three
        things could raise on the guarded path and none may here: a
        non-finite or negative diagonal, non-finite generators, and celerite2's
        own ``LinAlgError``. The first two are **sanitised** — replaced by
        values the compiled kernel can consume — and remembered in a flag; the
        third is converted to a NaN ``d`` inside
        :func:`ampere.backends.torch._celerite.factor`. All three, and any NaN
        that survives them, become ``-inf`` through one :func:`torch.where`,
        which keeps the result attached to the graph (a fresh ``-inf`` constant
        would not, and pyro's NUTS refuses a potential with no ``grad_fn``).
        """
        axis, order = self._axis(coordinates, kernel)
        residuals = self._tensor(residual).reshape(-1)[order]
        sorted_axis = axis[order]
        # W5.7: the builders get the raw coordinate, the recursion gets w(x).
        solve_axis = kernel.warped_coordinate(sorted_axis, values)
        diagonal = self._tensor(variance).reshape(-1)[order] + self.jitter**2

        # Sanitise, remembering that we did. A diagonal celerite2 cannot
        # factorise would come back as a quiet NaN rather than an exception, so
        # substituting a value it *can* factorise and forcing -inf afterwards is
        # both cheaper and more honest than hoping the NaN propagates.
        unusable = torch.logical_not(torch.isfinite(diagonal)) | (diagonal < 0.0)
        safe_diagonal = torch.where(unusable, torch.ones_like(diagonal), diagonal)
        refused = unusable.any()

        c, a, U, V = self._matrices(kernel, sorted_axis, safe_diagonal, values)
        finite = (
            torch.isfinite(a).all()
            & torch.isfinite(U).all()
            & torch.isfinite(V).all()
            & torch.isfinite(c).all()
        )
        refused = refused | torch.logical_not(finite)
        one = torch.ones_like(a)
        a = torch.where(finite, a, one)
        U = torch.where(finite, U, torch.zeros_like(U))
        V = torch.where(finite, V, torch.zeros_like(V))
        c = torch.where(finite, c, torch.ones_like(c))

        d, W = _celerite.factor(solve_axis, c, a, U, V)
        refused = refused | torch.logical_not(torch.isfinite(d).all()) | (d <= 0.0).any()
        safe_d = torch.where(torch.isfinite(d) & (d > 0.0), d, torch.ones_like(d))
        safe_W = torch.where(torch.isfinite(W), W, torch.zeros_like(W))

        value = self._marginal(solve_axis, c, U, safe_d, safe_W, residuals)
        failed = refused | torch.logical_not(torch.isfinite(value))
        return torch.where(failed, torch.full_like(value, -math.inf), value)

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        return float(
            self.log_marginal_likelihood_tensor(kernel, coordinates, residual, variance, values)
        )

    def provenance_config(self) -> Mapping[str, Any]:
        """``inference.md`` §10a fold-in 10: how this solver computes, for the attrs.

        Fixed rather than configurable, and recorded for exactly that reason:
        celerite2's compiled kernels are float64 CPU, so a run that used this
        solver factorised in float64 on the CPU whatever the rest of the
        problem did — including a problem that opted into float32 elsewhere.
        A reader of the archived run should not have to know that.
        """
        return {
            "dtype": str(self.TENSOR_DTYPE),
            "device": str(self.TENSOR_DEVICE),
            "library": "celerite2",
        }

    def configured(self, *, dtype: Any = None, device: Any = None) -> QuasisepGP:
        """Refuses, by name: celerite2's compiled kernels are float64 CPU.

        :meth:`DenseGP.configured` exists because a dense Cholesky genuinely
        runs in another precision and on another device. This solver's
        arithmetic happens inside ``celerite2.backprop``, a double-precision
        CPU C++ extension, so a float32 or GPU request here has no
        implementation to reach — and silently ignoring it would make
        ``provenance_config`` record a policy the run did not follow, which is
        worse than refusing. ``DenseGP`` is the solver to reach for; the
        decision-log row of 2026-09-07 records this as the measured price of
        the library choice.
        """
        if dtype is None and device is None:
            return self
        raise LikelihoodError(
            f"{self.NAME} on the torch backend computes inside celerite2's compiled kernels, "
            f"which are float64 on the CPU, so it has no float32 or non-CPU implementation to "
            f"configure (asked for dtype={dtype!r}, device={device!r}). Use DenseGP().configured"
            f"(...) if reduced precision or another device matters more than the O(N) solve, or "
            f"leave this solver as it is — provenance_config() records that the run factorised "
            f"in float64 on the CPU whatever the rest of the problem did."
        )

    # -- the rest of the GPSolver interface -----------------------------------

    def _apply_inverse(
        self,
        axis: torch.Tensor,
        c: torch.Tensor,
        U: torch.Tensor,
        d: torch.Tensor,
        W: torch.Tensor,
        right: torch.Tensor,
    ) -> torch.Tensor:
        """``(K + diag)⁻¹ right`` for an ``(n, nrhs)`` right-hand side, in O(N·nrhs)."""
        z = _celerite.solve_lower(axis, c, U, W, right)
        return _celerite.solve_upper(axis, c, U, W, z / d[:, None])

    @staticmethod
    def _precision_diagonal(
        axis: torch.Tensor,
        c: torch.Tensor,
        U: torch.Tensor,
        d: torch.Tensor,
        W: torch.Tensor,
    ) -> torch.Tensor:
        r"""``diag((K + diag(σ²))⁻¹)`` in O(N J³) — the recursion W2.3 deferred.

        **This is what discharges W2.3's deferral** (``DEVELOPMENT_PLAN.md``
        §2, 2026-09-05): "celerite2's public numpy interface has no O(N) route
        to that diagonal […] An O(N) route does exist (a backward accumulation
        […] over the inverse of celerite's unit-lower semiseparable factor),
        but each reimplements celerite2's internal factorisation convention in
        numpy". This backend has to know that convention anyway — it calls the
        compiled kernels directly rather than through ``GaussianProcess`` — so
        the coupling is already paid for, and it is a *tested* coupling: the
        conformance and backend suites assert these terms against
        :meth:`DenseGP.conditional_loo`, which computes them from a Cholesky.

        The derivation, so a reader need not rebuild it. With ``K + diag =
        L D Lᵀ`` and ``M = L⁻¹`` (unit lower triangular), ``A = Lᵀ⁻¹ D⁻¹ L⁻¹``
        gives

        .. math:: A_{ii} = \sum_{k \ge i} M_{ki}^2 / d_k .

        Column *i* of ``M`` solves ``L z = e_i``, and celerite's forward
        substitution — ``f_k = p_k \odot (f_{k-1} + W_{k-1} z_{k-1})``,
        ``z_k = y_k - U_k \cdot f_k`` with ``p_k = e^{-c(t_k - t_{k-1})}`` —
        becomes, once ``y = e_i`` is substituted in,

        .. math::
            f^{(i)}_{i+1} = p_{i+1} \odot W_i,\qquad
            f^{(i)}_k = G_k f^{(i)}_{k-1},\quad
            G_k = \mathrm{diag}(p_k)\,(I - W_{k-1} U_{k-1}^\top),

        and the transition ``G_k`` **does not depend on i**. That is the whole
        trick: every column of ``M`` is the same linear recursion started from
        a different vector, so the sum over columns collapses into one backward
        accumulation of a ``J x J`` matrix,

        .. math::
            A_{ii} = 1/d_i + w_i^\top R_i w_i,\qquad w_i = p_{i+1} \odot W_i,
            \qquad R_i = \frac{U_{i+1} U_{i+1}^\top}{d_{i+1}}
                       + G_{i+2}^\top R_{i+1} G_{i+2},

        with ``R_{N-1} = 0``. No inverses appear, so it is as stable as the
        factorisation itself.

        It is O(N) with a Python-loop constant — ``J`` is 2 for Matérn-3/2, so
        the arithmetic is trivial and the interpreter dominates. That is
        acceptable where a dense Cholesky is not merely slower but
        *impossible*: ``DenseGP`` needs 80 GB for the covariance alone at 10⁵
        points, and this needs none. It is a post-processing surface
        (``results.md`` §15 R2 does not store the decomposition by default),
        never part of a density.
        """
        n, rank = U.shape
        gaps = axis[1:] - axis[:-1]
        decays = torch.exp(-c[None, :] * gaps[:, None])  # (n-1, J): p_k for k = 1..n-1
        identity = torch.eye(rank, dtype=U.dtype, device=U.device)
        transitions = decays[:, :, None] * (
            identity[None, :, :] - W[:-1, :, None] * U[:-1, None, :]
        )  # (n-1, J, J): G_k for k = 1..n-1, indexed [k-1]
        starts = decays * W[:-1, :]  # (n-1, J): w_i for i = 0..n-2

        precision = torch.empty_like(d)
        precision[n - 1] = 1.0 / d[n - 1]
        accumulated = torch.zeros((rank, rank), dtype=U.dtype, device=U.device)
        for index in range(n - 2, -1, -1):
            if index + 2 <= n - 1:
                step = transitions[index + 1]
                accumulated = step.transpose(0, 1) @ accumulated @ step
            accumulated = accumulated + torch.outer(U[index + 1], U[index + 1]) / d[index + 1]
            precision[index] = 1.0 / d[index] + starts[index] @ accumulated @ starts[index]
        return precision

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """The leave-one-out conditional terms, in O(N) (W2.4 slice 2).

        The same closed form :meth:`DenseGP.conditional_loo` uses
        (Sundararajan & Keerthi 2001) over the same matrix, so the two agree
        exactly rather than approximately: with
        ``A = (K + diag(variance + jitter²))⁻¹``,
        ``log p_i = 0.5 log A_ii - [A r]_i² / (2 A_ii) - 0.5 log 2π``. What
        this solver supplies that the reference one could not is ``A_ii`` in
        linear time — see :meth:`_precision_diagonal`.
        """
        axis, order = self._axis(coordinates, kernel)
        residuals = self._tensor(residual).reshape(-1)
        diagonal = self._tensor(variance).reshape(-1) + self.jitter**2
        sorted_axis = axis[order]
        solve_axis = kernel.warped_coordinate(sorted_axis, values)
        c, U, d, W = self._guarded_factor(kernel, sorted_axis, diagonal[order], values)
        alpha = self._apply_inverse(solve_axis, c, U, d, W, residuals[order].reshape(-1, 1))[:, 0]
        precision = self._precision_diagonal(solve_axis, c, U, d, W)
        terms = 0.5 * torch.log(precision) - alpha**2 / (2.0 * precision) - 0.5 * _LOG_2PI
        restored = torch.empty_like(terms)
        restored[order] = terms
        return to_numpy(restored).astype(DTYPE, copy=False)

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The GP posterior at *at*, in O(N·M).

        O(N·M) rather than O(N), and the docstring says so rather than implying
        otherwise: the cross-covariance block is dense by construction — M
        outputs each need all N inputs — and only the solve against it is
        linear per column. :class:`DenseGP` pays O(N³) for the same answer.
        """
        axis, order = self._axis(coordinates, kernel)
        points = axis.reshape(-1, 1)
        residuals = self._tensor(residual).reshape(-1)
        diagonal = self._tensor(variance).reshape(-1) + self.jitter**2
        sorted_axis = axis[order]
        solve_axis = kernel.warped_coordinate(sorted_axis, values)
        c, U, d, W = self._guarded_factor(kernel, sorted_axis, diagonal[order], values)
        alpha = self._apply_inverse(solve_axis, c, U, d, W, residuals[order].reshape(-1, 1))[:, 0]
        target = points if at is None else _points(at).to(dtype=self.TENSOR_DTYPE)
        if target.shape[1] != 1:
            raise LikelihoodError(
                f"the conditioning grid has {target.shape[1]} coordinate dimension(s) but the "
                f"data have 1."
            )
        cross = self._tensor(kernel.matrix(target, points, values))[:, order]
        solved = self._apply_inverse(solve_axis, c, U, d, W, cross.transpose(0, 1))
        prior_variance = self._tensor(kernel.diagonal(target, values))
        return GPConditional(
            mean=to_numpy(cross @ alpha).astype(DTYPE, copy=False),
            variance=to_numpy(prior_variance - (cross * solved.transpose(0, 1)).sum(dim=1)).astype(
                DTYPE, copy=False
            ),
        )

    def latent_transform_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> torch.Tensor:
        """``f = L(θ) z`` in O(N), as a differentiable tensor (W2.14).

        The quasiseparable counterpart of
        :meth:`DenseGP.latent_transform_native`, and the one that matters at
        the scale the latent path is sized for: forming ``L`` densely is
        O(N³), so a 10⁵-sample latent fit reaches its correlation through this
        recursion or not at all.

        No exception control flow, and the same sanitise-then-force discipline
        :meth:`log_marginal_likelihood_native` uses: celerite's compiled
        kernels are handed values they can factorise whatever the
        hyperparameters were, and a NaN is then forced back in through one
        :func:`torch.where`, which keeps the result attached to the graph.
        """
        axis, order = self._axis(coordinates, kernel)
        points = axis.reshape(-1, 1)
        draws = self._tensor(whitened).reshape(-1)
        # The reference solver's ``float(...) or 1.0``, kept as a value so a
        # fitted amplitude survives it.
        mean_variance = torch.mean(self._tensor(kernel.diagonal(points, values)))
        scale = torch.where(mean_variance != 0.0, mean_variance, torch.ones_like(mean_variance))
        sorted_axis = axis[order]
        solve_axis = kernel.warped_coordinate(sorted_axis, values)
        stabiliser = torch.ones_like(sorted_axis) * (jitter * scale)
        usable = torch.isfinite(stabiliser) & (stabiliser >= 0.0)
        safe_stabiliser = torch.where(usable, stabiliser, torch.zeros_like(stabiliser))
        refused = torch.logical_not(usable.all())

        c, a, U, V = self._matrices(kernel, sorted_axis, safe_stabiliser, values)
        finite = (
            torch.isfinite(a).all()
            & torch.isfinite(U).all()
            & torch.isfinite(V).all()
            & torch.isfinite(c).all()
        )
        refused = refused | torch.logical_not(finite)
        a = torch.where(finite, a, torch.ones_like(a))
        U = torch.where(finite, U, torch.zeros_like(U))
        V = torch.where(finite, V, torch.zeros_like(V))
        c = torch.where(finite, c, torch.ones_like(c))

        d, W = _celerite.factor(solve_axis, c, a, U, V)
        refused = refused | torch.logical_not(torch.isfinite(d).all()) | (d <= 0.0).any()
        safe_d = torch.where(torch.isfinite(d) & (d > 0.0), d, torch.ones_like(d))
        safe_W = torch.where(torch.isfinite(W), W, torch.zeros_like(W))

        scaled = (draws[order] * torch.sqrt(safe_d)).reshape(-1, 1)
        transformed = (scaled + _celerite.matmul_lower(solve_axis, c, U, safe_W, scaled))[:, 0]
        restored = torch.empty_like(transformed)
        restored[order] = transformed
        failed = refused | torch.logical_not(torch.isfinite(restored).all())
        return torch.where(failed, torch.full_like(restored, math.nan), restored)

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """``f = L z`` in O(N) — the deterministic half of the latent-GP declaration.

        ``L`` here is the Cholesky-like factor of ``K`` alone, so the same
        stabilisation :meth:`DenseGP.latent_transform` applies — a jitter
        relative to the kernel's own scale — is applied here, and the two
        solvers therefore factorise the same matrix. celerite's ``L D Lᵀ``
        gives it as ``f = L √D z`` with ``L`` unit lower triangular, which is
        one ``matmul_lower`` rather than a triangular multiply.

        The contract surface: it takes the *guarded* factorisation, so
        celerite's quiet NaN becomes the reference path's loud refusal, which
        :meth:`latent_transform_native` may not raise.
        """
        axis, order = self._axis(coordinates, kernel)
        points = axis.reshape(-1, 1)
        draws = self._tensor(whitened).reshape(-1)
        scale = float(torch.mean(self._tensor(kernel.diagonal(points, values)))) or 1.0
        sorted_axis = axis[order]
        solve_axis = kernel.warped_coordinate(sorted_axis, values)
        stabiliser = torch.full_like(sorted_axis, jitter * scale)
        c, U, d, W = self._guarded_factor(kernel, sorted_axis, stabiliser, values, whitening=True)
        scaled = (draws[order] * torch.sqrt(d)).reshape(-1, 1)
        transformed = (scaled + _celerite.matmul_lower(solve_axis, c, U, W, scaled))[:, 0]
        restored = torch.empty_like(transformed)
        restored[order] = transformed
        return to_numpy(restored).astype(DTYPE, copy=False)


# ---------------------------------------------------------------------------
# The reduced-rank spectral solver (W5.4)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class HilbertSpaceGP(GPSolver):
    r"""Approximate ``O(N m + m³)`` by a Hilbert-space basis, in torch.

    The torch twin of :class:`ampere.core.HilbertSpaceGP`: the same
    declaration, the same basis, the same Woodbury algebra, and the same
    answer to ``tolerances.cross_backend``. What this class adds is the two
    things ``ampere.core``'s cannot have — a solve whose every step carries a
    gradient in the kernel hyperparameters, and a
    :meth:`log_marginal_likelihood_native` surface with no exception control
    flow for :mod:`ampere.backends.torch.problem` to compose into a NUTS
    density.

    **Why this is the solver a reduced-rank latent path wants.** The whitening
    is ``f = Phi z`` with ``z`` of length ``m``, so the sampler's dimension is
    the basis size rather than the sample count; the factorisation is an
    ``m`` by ``m`` Cholesky whose cost does not grow with ``N`` at all; and
    there is nothing to stabilise, because the reduced-rank factor times its
    own transpose is positive semi-definite by construction. A 10^4-point
    spectrum under NUTS is then a few dozen latent dimensions and one small
    dense factorisation per leapfrog step.

    See ``ampere.core.HilbertSpaceGP`` for the mathematics, the parameters and
    the warning about choosing ``basis_size`` and ``boundary_factor``
    together; nothing here changes any of it.
    """

    basis_size: int | tuple[int, ...] = DEFAULT_BASIS_SIZE
    boundary_factor: float = DEFAULT_BOUNDARY_FACTOR
    jitter: float = 0.0

    #: Class-level policy, never dataclass fields — see :class:`DenseGP`.
    TENSOR_DTYPE: ClassVar[torch.dtype] = DEFAULT_DTYPE
    TENSOR_DEVICE: ClassVar[torch.device] = DEFAULT_DEVICE

    NAME: ClassVar[str] = "HilbertSpaceGP"
    EXACT: ClassVar[bool] = False
    IMPLEMENTED: ClassVar[bool] = True
    STACKED_RESIDUALS: ClassVar[bool] = True

    DIFFERENTIABLE: ClassVar[bool] = True
    #: Not batchable. ``torch.linalg.cholesky_ex`` would map over a stack of
    #: capacitance matrices happily enough, but the basis is built from the
    #: *coordinates*, which are shared across a batch of parameter vectors and
    #: would have to be broadcast against them; nothing in the library needs it
    #: yet, and a flag claiming it would be a promise this class does not keep.
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __post_init__(self) -> None:
        object.__setattr__(self, "basis_size", normalise_counts(self.basis_size))
        if not math.isfinite(self.boundary_factor) or self.boundary_factor <= 0.0:
            raise LikelihoodError(
                f"HilbertSpaceGP's boundary_factor must be finite and > 0, got "
                f"{self.boundary_factor!r}."
            )
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"HilbertSpaceGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )

    # -- declarations ----------------------------------------------------------

    @property
    def counts(self) -> tuple[int, ...]:
        """Basis members per axis, normalised to a tuple."""
        return normalise_counts(self.basis_size)

    def latent_size(self, kernel: Kernel, n_samples: int) -> int:
        """``m``, not ``N``. See ``ampere.core.HilbertSpaceGP.latent_size``."""
        del kernel, n_samples
        return basis_size(self.counts)

    def provenance_config(self) -> Mapping[str, Any]:
        """The approximation *and* how this backend computed it.

        ``basis_size`` and ``boundary_factor`` are declarations and so are also
        dataclass fields, which is what keeps this backend's spec hash equal to
        the reference one's for the same problem (``results.md`` §14); the
        dtype and the device are configuration and are not.
        """
        return {
            "basis_size": [int(count) for count in self.counts],
            "boundary_factor": float(self.boundary_factor),
            "dtype": str(self.TENSOR_DTYPE),
            "device": str(self.TENSOR_DEVICE),
        }

    def configured(self, *, dtype: Any = None, device: Any = None) -> HilbertSpaceGP:
        """A copy computing elsewhere or in another precision. See :meth:`DenseGP.configured`."""
        return _configured(self, dtype, device)

    def check_compatible(self, kernel: Kernel, observed: Any) -> None:
        """The core solver's checks, unchanged — they are declarations, not arithmetic."""
        super().check_compatible(kernel, observed)
        selected = kernel.selected_axes([axis.name for axis in observed.axes])
        check_spectral_support(kernel, len(selected), owner=self.NAME)
        if len(self.counts) != len(selected):
            raise LikelihoodError(
                f"{self.NAME} was declared with basis_size={self.basis_size!r} — "
                f"{len(self.counts)} axis count(s) — but the kernel selects {len(selected)} "
                f"axis/axes {selected!r} of this {type(observed).__name__}. The basis is a "
                f"tensor product with one count per axis."
            )

    # -- internals -------------------------------------------------------------

    def _tensor(self, value: Any) -> torch.Tensor:
        return as_tensor(value, dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)

    def _points(self, coordinates: Any) -> torch.Tensor:
        tensor = self._tensor(coordinates)
        return tensor.reshape(-1, 1) if tensor.ndim == 1 else tensor

    def _basis(self, kernel: Kernel, coordinates: Any) -> HilbertSpaceBasis:
        """The box and the frequency grid — numpy, because they are data.

        ``to_numpy`` rather than ``np.asarray``: a coordinate block may be a
        tensor on another device, and it never carries a gradient (a fitted
        coordinate is not a thing ampere has), so nothing is detached that
        anything wanted.
        """
        points = kernel.select(self._points(coordinates))
        return hilbert_basis(to_numpy(points), self.counts, self.boundary_factor)

    def _scaled_basis(
        self,
        kernel: Kernel,
        basis: HilbertSpaceBasis,
        coordinates: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        r"""The basis scaled by the square root of the spectral density, graph intact.

        ``Phi`` is built through :class:`TorchOps`, so it lands on this
        solver's own dtype and device whatever namespace the coordinates
        arrived in; ``S`` comes from the *kernel*'s namespace, which is where
        the gradient in the amplitude and the length scale lives.
        ``clamp(min=0)`` rather than a check: a realised density must not
        raise, and a spectral density is non-negative for every admissible
        hyperparameter, so a negative one is a hyperparameter outside its
        prior and is already on its way to ``-inf``.
        """
        ops = TorchOps(self.TENSOR_DTYPE, self.TENSOR_DEVICE)
        matrix = basis_matrix(basis, kernel.select(self._points(coordinates)), ops)
        density = self._tensor(spectral_values(kernel, basis, values))
        return matrix * torch.sqrt(torch.clamp(density, min=0.0))[None, :]

    def _diagonal(self, variance: Any) -> torch.Tensor:
        return self._tensor(variance) + self.jitter**2

    def _factor(
        self, scaled: torch.Tensor, diagonal: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        r"""``(D^-1 Phi, chol(I + Phi^T D^-1 Phi), info)`` — the one factorisation."""
        weighted = scaled / diagonal[:, None]
        identity = torch.eye(scaled.shape[1], dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)
        capacitance = identity + scaled.transpose(0, 1) @ weighted
        factor, info = torch.linalg.cholesky_ex(capacitance)
        return weighted, factor, info

    @staticmethod
    def _solve(
        scaled: torch.Tensor,
        diagonal: torch.Tensor,
        weighted: torch.Tensor,
        factor: torch.Tensor,
        right: torch.Tensor,
    ) -> torch.Tensor:
        r"""``(K + D)^-1 right`` by Woodbury, for an ``(n,)`` or ``(n, k)`` block."""
        vector = right.ndim == 1
        block = right.reshape(-1, 1) if vector else right
        direct = block / diagonal[:, None]
        inner = torch.cholesky_solve(scaled.transpose(0, 1) @ direct, factor, upper=False)
        solved = direct - weighted @ inner
        return solved.reshape(-1) if vector else solved

    def _guard(self, message: str, ok: Any) -> None:
        if not bool(ok):
            raise LikelihoodError(message)

    # -- the interface ---------------------------------------------------------

    def log_marginal_likelihood_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """The marginal likelihood of the approximation, **without exceptions** (W2.13).

        The surface :mod:`ampere.backends.torch.problem` composes. A
        factorisation that fails becomes ``-inf`` through :func:`torch.where`,
        exactly as :meth:`DenseGP.log_marginal_likelihood_native` does, so a
        NUTS step that wanders into an inadmissible hyperparameter is rejected
        rather than raising in the middle of the hot loop.
        """
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        diagonal = self._diagonal(variance)
        weighted, factor, info = self._factor(scaled, diagonal)
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        log_determinant = torch.log(diagonal).sum() + 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = (residuals * alpha).sum()
        value = -0.5 * (
            quadratic + _columns(residuals) * log_determinant + residuals.numel() * _LOG_2PI
        )
        failed = torch.logical_or(info != 0, torch.logical_not(torch.isfinite(value)))
        failed = torch.logical_or(failed, torch.logical_not(torch.all(diagonal > 0.0)))
        return torch.where(failed, torch.full_like(value, -math.inf), value)

    def log_marginal_likelihood_tensor(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: Any,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """The same quantity as a differentiable tensor, with the loud refusals."""
        diagonal = self._diagonal(variance)
        self._guard(
            f"{self.NAME} needs a strictly positive noise diagonal: the Woodbury identity it "
            f"solves by inverts diag(sigma^2 + jitter^2) directly, so a zero uncertainty is a "
            f"division by zero rather than an ill-conditioned matrix. Pass "
            f"{self.NAME}(jitter=...), or use DenseGP.",
            bool(torch.all(torch.isfinite(diagonal))) and bool(torch.all(diagonal > 0.0)),
        )
        value = self.log_marginal_likelihood_native(kernel, coordinates, residual, variance, values)
        self._guard(
            "the reduced-rank capacitance matrix I + Phi^T D^-1 Phi is not positive definite, "
            "so the Woodbury solve failed. It is positive definite for every admissible "
            "hyperparameter, so this is a numerical rather than a structural failure: reduce "
            "basis_size, or raise the jitter.",
            torch.isfinite(value),
        )
        return value

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> float:
        return float(
            self.log_marginal_likelihood_tensor(kernel, coordinates, residual, variance, values)
        )

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Sundararajan & Keerthi from the Woodbury diagonal. See the core solver."""
        points = _as_points(coordinates, "data coordinates")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, factor, info = self._factor(scaled, diagonal)
        self._guard(
            "the reduced-rank capacitance matrix I + Phi^T D^-1 Phi is not positive definite, "
            "so the Woodbury solve failed.",
            int(info) == 0,
        )
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        triangular = torch.linalg.solve_triangular(factor, weighted.transpose(0, 1), upper=False)
        precision_diagonal = 1.0 / diagonal - (triangular * triangular).sum(dim=0)
        rows = int(residuals.shape[0])
        columns = _columns(residuals)
        quadratic = (alpha.reshape(rows, columns) ** 2).sum(dim=1)
        terms = columns * (0.5 * torch.log(precision_diagonal) - 0.5 * _LOG_2PI) - quadratic / (
            2.0 * precision_diagonal
        )
        return to_numpy(terms).astype(DTYPE, copy=False)

    def condition(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        residual: np.ndarray,
        variance: np.ndarray,
        values: Mapping[str, Any],
        at: np.ndarray | None = None,
    ) -> GPConditional:
        """The posterior of the rank-``m`` process. See the core solver."""
        points = _as_points(coordinates, "data coordinates")
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, factor, info = self._factor(scaled, diagonal)
        self._guard(
            "the reduced-rank capacitance matrix I + Phi^T D^-1 Phi is not positive definite, "
            "so the Woodbury solve failed.",
            int(info) == 0,
        )
        residuals = _right_hand_side(self._tensor(residual))
        alpha = self._solve(scaled, diagonal, weighted, factor, residuals)
        target = points if at is None else _as_points(at, "conditioning grid", points.shape[1])
        target_scaled = self._scaled_basis(kernel, basis, target, values)
        mean = target_scaled @ (scaled.transpose(0, 1) @ alpha)
        triangular = torch.linalg.solve_triangular(
            factor, target_scaled.transpose(0, 1), upper=False
        )
        posterior = (triangular * triangular).sum(dim=0)
        return GPConditional(
            mean=to_numpy(mean).astype(DTYPE, copy=False),
            variance=to_numpy(posterior).astype(DTYPE, copy=False),
        )

    def latent_transform_native(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> torch.Tensor:
        r"""``f = Phi z`` as a differentiable tensor, without exceptions (W2.14).

        ``z`` carries :meth:`latent_size` entries. There is no factorisation
        here at all — the reduced-rank factor *is* the whitening — so the
        ``jitter`` argument the two exact solvers need is accepted and unused,
        and this surface has nothing that can fail to sanitise into NaN.
        """
        del jitter
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        draws = self._tensor(whitened)
        vector = draws.ndim == 1
        drawn = scaled @ (draws.reshape(-1, 1) if vector else draws)
        return drawn.reshape(-1) if vector else drawn

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: np.ndarray,
        whitened: np.ndarray,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """``f = Phi z`` — the contract surface, so it raises where the native one does not."""
        points = _as_points(coordinates, "data coordinates")
        expected = self.latent_size(kernel, points.shape[0])
        if int(np.shape(whitened)[0]) != expected:
            raise LikelihoodError(
                f"{self.NAME} whitens {expected} basis coefficient(s), but it was handed "
                f"{int(np.shape(whitened)[0])} whitened value(s). This solver's latent block is "
                f"the basis size m, not the sample count — see GPSolver.latent_size."
            )
        drawn = self.latent_transform_native(kernel, points, whitened, values, jitter=jitter)
        if not bool(torch.all(torch.isfinite(drawn))):
            raise LikelihoodError(
                "the reduced-rank whitening f = Phi diag(sqrt(S)) z is not finite. The usual "
                "cause is a kernel amplitude large enough that amplitude**2 overflows float64."
            )
        return to_numpy(drawn).astype(DTYPE, copy=False)
