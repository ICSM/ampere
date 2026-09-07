"""celerite2's semiseparable kernels, under ``torch.autograd``.

What this module is, and why ampere writes it
----------------------------------------------
``DEVELOPMENT_PLAN.md`` §6 deferred one choice to this track: "**Torch GP
solver library**: GPyTorch structured solvers vs celerite2's experimental torch
interface for ``QuasisepGP``. Evaluate both against the conformance suite."
W2.4 slice 2 evaluated both, and the measurement settled it (the decision-log
row of 2026-09-07 carries the table). The short version is two facts:

* **celerite2 has no torch interface.** Version 0.3.3 ships ``celerite2.jax``,
  ``celerite2.pymc``, ``celerite2.pymc3`` and ``celerite2.theano`` — and
  nothing for torch. What it *does* ship, and what all four of those wrappers
  are built on, is :mod:`celerite2.backprop`: the compiled forward **and
  reverse** passes of the semiseparable factorisation and its solves. So the
  "torch interface" is a thin autograd shim over kernels ampere already
  depends on, and this module is it — the exact counterpart of
  ``celerite2/jax/ops.py``, which registers the same C++ entry points as jax
  primitives with the same reverse rules.
* **GPyTorch has no quasiseparable operator at all.** Its structured routes
  are Toeplitz (a *regular* grid), Kronecker (a product grid) and SKI /
  inducing points (approximate). ampere's coordinates are irregular by
  contract (``results_schema.md`` §16 forbids a regular-grid assumption), so
  the only exact route GPyTorch offers on this problem is a dense Cholesky —
  which is :class:`~ampere.backends.torch.DenseGP`, already shipped, under
  another name and a 40 MB dependency.

The convention, stated once
---------------------------
celerite2 factorises ``K + diag(a) = L D Lᵀ`` where ``L`` is **unit** lower
triangular and semiseparable of rank ``J``::

    L_nm = sum_j U_nj W_mj exp(-c_j (t_n - t_m))      n > m
    L_nn = 1

The five inputs are therefore ``t`` (sorted coordinates), ``c`` (the ``J``
decay rates), ``a`` (the diagonal of the full matrix, kernel variance
included), and the generators ``U``, ``V`` (both ``(N, J)``).
:func:`factor` returns ``d`` (the diagonal of ``D``) and ``W`` (the second
generator of ``L``); :func:`solve_lower` and :func:`solve_upper` apply
``L⁻¹`` and ``L⁻ᵀ``; :func:`matmul_lower` applies the strictly lower part of
``L`` itself.

Each is O(N J²) and each has a compiled reverse pass, so the whole marginal
likelihood — and every kernel hyperparameter that reaches it through ``c``,
``a``, ``U`` and ``V`` — is differentiable in linear time. That is the point:
the *residual* was already differentiable through a dense Cholesky, and the
amplitude and length scale were the things W2.4's slice-1 finding said were
not.

Why an autograd ``Function`` and not a differentiable torch recursion
---------------------------------------------------------------------
The factorisation is a sequential scan over N with a J x J state, and torch has
no scan primitive: written in torch it is a Python loop of N iterations
building an N-node autograd graph. Measured, that is ~1 800x slower than this
module at 1 000 points and gets worse with N — the third row of the
decision-log table. Wrapping the compiled kernels is what makes the solver
linear *in practice* as well as in asymptotics.

The price is stated honestly in the capability flags on
:class:`~ampere.backends.torch.QuasisepGP`: these kernels are CPU float64 C++,
so the solver declares ``DEVICE = "cpu"`` and ``BATCHABLE = False``, and a
``device=`` or float32 request is refused by name rather than silently ignored.

Failure signalling
------------------
``celerite2.backprop`` raises :class:`celerite2.backprop.LinAlgError` when the
factorisation meets a non-positive-definite matrix, and returns **quiet NaN**
for a diagonal that could not belong to a covariance at all (W2.3's carried
finding). :func:`factor` therefore converts the exception into a NaN ``d``,
so that a caller on the realised path can turn a NaN into ``-inf`` through
:func:`torch.where` rather than catching an exception inside the hot loop
(``inference.md`` §10a). Callers that *may* raise — the numpy contract
surface — check the preconditions themselves before calling and diagnose the
NaN afterwards; see :meth:`~ampere.backends.torch.QuasisepGP._factorise`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import torch

__all__ = [
    "LinAlgError",
    "factor",
    "general_matmul_lower",
    "general_matmul_upper",
    "matmul_lower",
    "solve_lower",
    "solve_upper",
]


def _backprop() -> Any:
    """``celerite2.backprop``, imported on first use.

    Lazy for the reason ``ampere.core.likelihood`` gives for its own celerite2
    import: it is a compiled extension costing ~20 ms to load, and a torch
    problem with no quasiseparable GP in it should not pay for one. Same idiom,
    one level down.
    """
    # celerite2.backprop is the compiled extension: no stubs, so pyrefly
    # cannot resolve it. The same suppression ``ampere.core.likelihood`` applies
    # to ``celerite2.driver``, for the same reason -- the names are public and
    # documented, and the *shapes* are checked by the definitions.json contract
    # this module transcribes.
    import celerite2.backprop  # type: ignore[missing-import]

    return celerite2.backprop


class LinAlgError(Exception):
    """Raised nowhere here; re-exported so callers can name celerite2's own.

    :func:`factor` swallows ``celerite2.backprop.LinAlgError`` deliberately
    (see the module docstring), so this alias exists for a caller that wants
    to catch the underlying library's exception around a *direct* use of the
    compiled kernels.
    """


def _numpy(tensor: torch.Tensor) -> np.ndarray:
    """A contiguous float64 numpy view of *tensor*, detached.

    Detaching is safe and is not a gradient leak: every function here is the
    forward half of a :class:`torch.autograd.Function`, so the graph is carried
    by the ``Function`` and rebuilt from the compiled reverse pass, exactly as
    ``celerite2/jax/ops.py`` does it with ``custom_vjp``.
    """
    return np.ascontiguousarray(tensor.detach().cpu().numpy(), dtype=np.float64)


def _like(array: np.ndarray, reference: torch.Tensor) -> torch.Tensor:
    return torch.as_tensor(array, dtype=reference.dtype, device=reference.device)


class _Factor(torch.autograd.Function):
    """``factor(t, c, a, U, V) -> (d, W)``: the LDLᵀ of ``K + diag(a)``."""

    @staticmethod
    def forward(  # pyrefly: ignore[bad-override]
        ctx: Any,
        t: torch.Tensor,
        c: torch.Tensor,
        a: torch.Tensor,
        U: torch.Tensor,
        V: torch.Tensor,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        backprop = _backprop()
        points, decay = _numpy(t), _numpy(c)
        diagonal, left, right = _numpy(a), _numpy(U), _numpy(V)
        size, rank = left.shape
        try:
            d, W, S = backprop.factor_fwd(
                points,
                decay,
                diagonal,
                left,
                right,
                np.empty(size, dtype=np.float64),
                np.empty((size, rank), dtype=np.float64),
                np.empty((size, rank, rank), dtype=np.float64),
            )
        except backprop.LinAlgError:
            # Not an error here: a caller on the realised path must not see an
            # exception (inference.md §10a), and a NaN propagates to -inf
            # through the torch.where every caller already applies. A caller
            # that may raise diagnoses the NaN itself, with a message naming
            # the cause rather than repeating celerite2's.
            nan = np.full(size, np.nan, dtype=np.float64)
            return _like(nan, t), _like(np.full((size, rank), np.nan), t)
        ctx.celerite = (points, decay, diagonal, left, right, d, W, S)
        return _like(d, t), _like(W, t)

    @staticmethod
    def backward(  # pyrefly: ignore[bad-override]
        ctx: Any, bd: torch.Tensor, bW: torch.Tensor
    ) -> tuple[torch.Tensor | None, ...]:
        saved = getattr(ctx, "celerite", None)
        if saved is None:  # the factorisation failed; there is no gradient to give
            return (None, None, None, None, None)
        backprop = _backprop()
        points, decay, diagonal, left, right, d, W, S = saved
        size, rank = left.shape
        bt, bc, ba, bU, bV = backprop.factor_rev(
            points,
            decay,
            diagonal,
            left,
            right,
            d,
            W,
            S,
            np.ascontiguousarray(bd.detach().cpu().numpy(), dtype=np.float64),
            np.ascontiguousarray(bW.detach().cpu().numpy(), dtype=np.float64),
            np.empty(size, dtype=np.float64),
            np.empty(rank, dtype=np.float64),
            np.empty(size, dtype=np.float64),
            np.empty((size, rank), dtype=np.float64),
            np.empty((size, rank), dtype=np.float64),
        )
        return (
            _like(bt, bd),
            _like(bc, bd),
            _like(ba, bd),
            _like(bU, bd),
            _like(bV, bd),
        )


def _five_argument_op(name: str) -> type[torch.autograd.Function]:
    """The four ``(t, c, U, V, Y) -> Z`` kernels, which share one shape.

    ``solve_lower``, ``solve_upper``, ``matmul_lower`` and ``matmul_upper``
    have identical signatures, identical reverse signatures and identical
    buffer shapes; only the compiled entry point differs. Writing the shim once
    is not merely economy — four hand-copied buffer lists is exactly the sort
    of place a transposed argument survives review.
    """

    class _Op(torch.autograd.Function):
        @staticmethod
        def forward(  # pyrefly: ignore[bad-override]
            ctx: Any,
            t: torch.Tensor,
            c: torch.Tensor,
            U: torch.Tensor,
            V: torch.Tensor,
            Y: torch.Tensor,
        ) -> torch.Tensor:
            backprop = _backprop()
            points, decay = _numpy(t), _numpy(c)
            left, right, target = _numpy(U), _numpy(V), _numpy(Y)
            size, rank = left.shape
            columns = target.shape[1]
            Z, F = getattr(backprop, f"{name}_fwd")(
                points,
                decay,
                left,
                right,
                target,
                np.empty((size, columns), dtype=np.float64),
                np.empty((size, rank, columns), dtype=np.float64),
            )
            ctx.celerite = (points, decay, left, right, target, Z, F)
            return _like(Z, t)

        @staticmethod
        def backward(  # pyrefly: ignore[bad-override]
            ctx: Any, bZ: torch.Tensor
        ) -> tuple[torch.Tensor | None, ...]:
            backprop = _backprop()
            points, decay, left, right, target, Z, F = ctx.celerite
            size, rank = left.shape
            columns = target.shape[1]
            bt, bc, bU, bV, bY = getattr(backprop, f"{name}_rev")(
                points,
                decay,
                left,
                right,
                target,
                Z,
                F,
                np.ascontiguousarray(bZ.detach().cpu().numpy(), dtype=np.float64),
                np.empty(size, dtype=np.float64),
                np.empty(rank, dtype=np.float64),
                np.empty((size, rank), dtype=np.float64),
                np.empty((size, rank), dtype=np.float64),
                np.empty((size, columns), dtype=np.float64),
            )
            return (_like(bt, bZ), _like(bc, bZ), _like(bU, bZ), _like(bV, bZ), _like(bY, bZ))

    _Op.__name__ = f"_{name.title().replace('_', '')}"
    return _Op


_SolveLower = _five_argument_op("solve_lower")
_SolveUpper = _five_argument_op("solve_upper")
_MatmulLower = _five_argument_op("matmul_lower")
_MatmulUpper = _five_argument_op("matmul_upper")


def factor(
    t: torch.Tensor, c: torch.Tensor, a: torch.Tensor, U: torch.Tensor, V: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """``K + diag(a) = L D Lᵀ``: the diagonal ``d`` of ``D`` and ``L``'s generator ``W``.

    ``d`` comes back all-NaN when the matrix is not positive definite, rather
    than raising; see the module docstring.
    """
    return _Factor.apply(t, c, a, U, V)  # type: ignore[return-value]


def solve_lower(
    t: torch.Tensor, c: torch.Tensor, U: torch.Tensor, W: torch.Tensor, Y: torch.Tensor
) -> torch.Tensor:
    """``L⁻¹ Y`` for the unit lower-triangular factor, in O(N J² · nrhs)."""
    return _SolveLower.apply(t, c, U, W, Y)  # type: ignore[return-value]


def solve_upper(
    t: torch.Tensor, c: torch.Tensor, U: torch.Tensor, W: torch.Tensor, Y: torch.Tensor
) -> torch.Tensor:
    """``L⁻ᵀ Y``, in O(N J² · nrhs)."""
    return _SolveUpper.apply(t, c, U, W, Y)  # type: ignore[return-value]


def matmul_lower(
    t: torch.Tensor, c: torch.Tensor, U: torch.Tensor, W: torch.Tensor, Y: torch.Tensor
) -> torch.Tensor:
    """``tril(L, -1) Y`` — the strictly lower part, so ``L Y = Y + matmul_lower(...)``."""
    return _MatmulLower.apply(t, c, U, W, Y)  # type: ignore[return-value]


def matmul_upper(
    t: torch.Tensor, c: torch.Tensor, U: torch.Tensor, W: torch.Tensor, Y: torch.Tensor
) -> torch.Tensor:
    """``triu(Lᵀ, 1) Y``."""
    return _MatmulUpper.apply(t, c, U, W, Y)  # type: ignore[return-value]


def general_matmul_lower(
    t1: torch.Tensor,
    t2: torch.Tensor,
    c: torch.Tensor,
    U: torch.Tensor,
    V: torch.Tensor,
    Y: torch.Tensor,
) -> torch.Tensor:
    """The cross-covariance product for ``t1`` against ``t2``, lower half.

    No reverse pass exists for this one in ``celerite2.backprop`` (it has
    ``has_rev = False``), so it is exposed as a plain, **non-differentiable**
    helper. That is not a limitation in practice: it is used only by
    :meth:`~ampere.backends.torch.QuasisepGP.condition`, which is a
    post-processing surface returning numpy, never part of a density.
    """
    backprop = _backprop()
    size = int(t1.numel())
    rank = int(c.numel())
    columns = int(Y.shape[1])
    Z, _ = backprop.general_matmul_lower_fwd(
        _numpy(t1),
        _numpy(t2),
        _numpy(c),
        _numpy(U),
        _numpy(V),
        _numpy(Y),
        np.empty((size, columns), dtype=np.float64),
        np.empty((int(t2.numel()), rank, columns), dtype=np.float64),
    )
    return _like(Z, t1)


def general_matmul_upper(
    t1: torch.Tensor,
    t2: torch.Tensor,
    c: torch.Tensor,
    U: torch.Tensor,
    V: torch.Tensor,
    Y: torch.Tensor,
) -> torch.Tensor:
    """The upper half of :func:`general_matmul_lower`'s product. Not differentiable."""
    backprop = _backprop()
    size = int(t1.numel())
    rank = int(c.numel())
    columns = int(Y.shape[1])
    Z, _ = backprop.general_matmul_upper_fwd(
        _numpy(t1),
        _numpy(t2),
        _numpy(c),
        _numpy(U),
        _numpy(V),
        _numpy(Y),
        np.empty((size, columns), dtype=np.float64),
        np.empty((int(t2.numel()), rank, columns), dtype=np.float64),
    )
    return _like(Z, t1)
