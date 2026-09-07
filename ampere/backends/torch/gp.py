"""GP solvers on the torch path. Slice 1 ships the dense one.

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

Slice 2 owes ``QuasisepGP``
----------------------------
``ampere.core.QuasisepGP`` is celerite2's numpy solver over an **exact** rank-2
Matérn-3/2 representation (W2.3) — never celerite2's approximate
``Matern32Term``. The torch equivalent has two candidate libraries, GPyTorch's
structured solvers and celerite2's own torch interface, and
``DEVELOPMENT_PLAN.md`` §6 asks for them to be *measured* against the
conformance suite rather than chosen from documentation. That measurement is
W2.4 slice 2, and until it is made this backend declares ``DENSE`` only, so the
quasiseparable rows skip with a reason naming what is owed rather than passing
against a solver nobody has checked.

One finding carried into that work, from W2.3: **celerite2 returns quiet NaN
where a Cholesky raises.** ``ampere.core.QuasisepGP`` guards its preconditions
(finite diagonal, positive amplitude) before calling celerite2, and any torch
celerite2 path must do the same — a NaN log-likelihood has to surface through
§4.5's failure-signalling route (``-inf`` with a recorded reason when
non-strict, a raise under ``strict``), never as a silent NaN in the chain.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.core import DTYPE, GPConditional, GPSolver, Kernel
from ampere.core import Matern32 as _CoreMatern32
from ampere.core import SquaredExponential as _CoreSquaredExponential
from ampere.core.exceptions import LikelihoodError

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, to_numpy

__all__ = ["DenseGP", "Matern32", "SquaredExponential"]

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)


def _points(coordinates: Any) -> torch.Tensor:
    """``(n, d)`` coordinates from an ``(n,)`` or ``(n, d)`` array or tensor.

    A bare one-dimensional input is read as a column of ``n`` one-dimensional
    points, matching ``ampere.core.Kernel.matrix``'s own convention.
    """
    tensor = as_tensor(coordinates, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
    return tensor[:, None] if tensor.ndim == 1 else tensor


def _separation(left: Any, right: Any) -> torch.Tensor:
    """Euclidean separation between two coordinate sets, ``(n, m)``.

    Euclidean in the coordinate space, which is why
    :meth:`~ampere.core.GPSolver.check_compatible` requires every coordinate
    axis to share one unit — a single isotropic length scale is meaningless
    across mixed units.
    """
    a, b = _points(left), _points(right)
    difference = a[:, None, :] - b[None, :, :]
    return torch.sqrt(torch.sum(difference * difference, dim=-1))


class _TorchKernel(Kernel):
    """Shared plumbing for the torch kernels: torch separations, torch covariances.

    **W2.13.** W2.4 slice 1 shipped a torch :class:`DenseGP` that built its
    covariance by calling ``ampere.core.Kernel.matrix`` — which computes in
    numpy and coerces every hyperparameter with ``float()``. The Cholesky was
    differentiable and the *hyperparameters were not*: the graph was cut at the
    covariance, so a fitted amplitude or length scale received no gradient at
    all. That was slice 1's principal carried finding, and W2.5 had already
    solved it the same way on jax. These classes are the fix.

    They **subclass the core kernels** rather than redeclaring them, exactly as
    jax's do: same ``FAMILY``, same ``HYPERPARAMETERS``, same
    ``QUASISEPARABLE`` flag, same ``NAME``, so a problem's declaration — and
    therefore its spec hash — is unchanged by which backend computes it.
    ``matrix`` and ``diagonal`` are overridden because the core's build their
    separations in numpy; nothing else is.

    The core's ``_positive`` check on the hyperparameters is not reproduced,
    and its absence is not silent: a non-positive length scale gives a
    non-finite covariance, the Cholesky then fails, and :class:`DenseGP`'s two
    surfaces turn that into the failure §4.5 asks for (an exception on the
    contract path, ``-inf`` on the traced one). The declaration is what keeps
    it from arising: a kernel hyperparameter is declared with a prior on the
    positive half-line and a ``Log`` bijection, so no sampler proposes one.
    """

    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> torch.Tensor:
        raise NotImplementedError

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> torch.Tensor:
        """Dense covariance between two coordinate sets, ``(n, m)``, in torch."""
        return self._covariance(_separation(left, right), self.resolve(values))

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> torch.Tensor:
        """The prior variance at each coordinate; ``k(0)`` for a stationary kernel."""
        n = int(_points(coordinates).shape[0])
        zeros = torch.zeros(n, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        return self._covariance(zeros, self.resolve(values))


class Matern32(_TorchKernel, _CoreMatern32):
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

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> torch.Tensor:
        amplitude = as_tensor(values["amplitude"], dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        length_scale = as_tensor(values["length_scale"], dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        scaled = _SQRT3 * as_tensor(separation) / length_scale
        return amplitude * amplitude * (1.0 + scaled) * torch.exp(-scaled)


class SquaredExponential(_TorchKernel, _CoreSquaredExponential):
    r"""Squared exponential (RBF) in torch: legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Not quasiseparable, and inherits that declaration too, so a quasiseparable
    solver refuses it by name — the asymmetry that made Matérn the default.
    """

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> torch.Tensor:
        amplitude = as_tensor(values["amplitude"], dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        length_scale = as_tensor(values["length_scale"], dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        scaled = as_tensor(separation) / length_scale
        return amplitude * amplitude * torch.exp(-0.5 * scaled * scaled)


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

    #: The four capability flags. W2.4 declared them against a parts set that
    #: did not yet include a solver and recorded the gap; **W2.13 closed it**
    #: (``inference.md`` §10a, fold-in 7), so these are now what makes a
    #: composed problem report ``backend="torch"`` — and what makes the same
    #: problem left with ``ampere.core.DenseGP`` a refusal rather than a
    #: silently non-differentiable GP.
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
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
        """``log N(residual; 0, K + diag(variance))`` as a differentiable tensor."""
        total = self._covariance(kernel, coordinates, variance, values)
        factor = self._cholesky(total)
        residuals = self._tensor(residual).reshape(-1)
        alpha = self._solve(factor, residuals)
        log_determinant = 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = residuals @ alpha
        return -0.5 * (quadratic + log_determinant + residuals.numel() * _LOG_2PI)

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
        residuals = self._tensor(residual).reshape(-1)
        alpha = self._solve(factor, residuals)
        log_determinant = 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = residuals @ alpha
        value = -0.5 * (quadratic + log_determinant + residuals.numel() * _LOG_2PI)
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
        ``architecture.md`` §5 says a float32 opt-out must be visible in a
        run's provenance, and this is where it will be visible when slice 2
        offers one.
        """
        return {"dtype": str(self.TENSOR_DTYPE), "device": str(self.TENSOR_DEVICE)}

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
        residuals = self._tensor(residual).reshape(-1)
        alpha = self._solve(factor, residuals)
        identity = torch.eye(residuals.numel(), dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE)
        precision_diagonal = torch.diagonal(self._solve(factor, identity))
        terms = (
            0.5 * torch.log(precision_diagonal)
            - alpha**2 / (2.0 * precision_diagonal)
            - 0.5 * _LOG_2PI
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
        residuals = self._tensor(residual).reshape(-1)
        mean = cross @ self._solve(factor, residuals)
        solved = self._solve(factor, cross.transpose(0, 1))
        prior_variance = self._tensor(kernel.diagonal(target, values))
        posterior = prior_variance - (cross * solved.transpose(0, 1)).sum(dim=1)
        return GPConditional(
            mean=to_numpy(mean).astype(DTYPE, copy=False),
            variance=to_numpy(posterior).astype(DTYPE, copy=False),
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
        """``f = L z`` — the deterministic half of the latent-GP declaration."""
        points = _as_points(coordinates, "data coordinates")
        covariance = self._tensor(kernel.matrix(points, points, values))
        if not bool(torch.all(torch.isfinite(covariance))):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        scale = float(torch.diagonal(covariance).mean()) or 1.0
        stabilised = covariance + torch.eye(
            covariance.shape[0], dtype=self.TENSOR_DTYPE, device=self.TENSOR_DEVICE
        ) * (jitter * scale)
        factor, info = torch.linalg.cholesky_ex(stabilised)
        if int(info) != 0:
            raise LikelihoodError(
                f"the kernel matrix K is not positive definite, so the whitening transform "
                f"f = L z is undefined (leading minor {int(info)}). Increase the jitter "
                f"argument, or check the kernel hyperparameters."
            )
        drawn = factor @ self._tensor(whitened).reshape(-1)
        return to_numpy(drawn).astype(DTYPE, copy=False)


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
