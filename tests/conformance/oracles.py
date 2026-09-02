"""Closed forms the battery compares against.

``architecture.md`` §2 makes the reference path the "conformance oracle", but
an oracle that is itself the thing under test proves nothing. So every
comparison in this suite is against one of exactly two things: ``scipy``, or a
formula in this module — each of which is transcribed from the contract's own
statement of it (a kernel's defining equation, ``parameters.md`` §6's change of
variables) rather than lifted from an implementation.

Nothing here imports from ``ampere`` beyond the declaration types it needs to
read; if a function in this module started calling ampere to compute a number,
the row using it would stop being a test.
"""

from __future__ import annotations

import math

import numpy as np

from ampere.core import Parameter, ParameterSet

from .protocol import KernelFamily

__all__ = [
    "analytic_constrain",
    "analytic_diagonal_gaussian_log_prob",
    "analytic_lnprior",
    "analytic_log_abs_det",
    "analytic_prior_transform",
    "free_slices",
    "kernel_matrix",
    "matern32_matrix",
    "squared_exponential_matrix",
    "summed_log_abs_det",
]


# ---------------------------------------------------------------------------
# The flat free-parameter layout
# ---------------------------------------------------------------------------


def free_slices(declaration: ParameterSet) -> list[tuple[Parameter, slice]]:
    """``(parameter, slice)`` for each non-fixed parameter, in declaration order.

    The layout — declaration order, array parameters occupying ``size``
    consecutive entries, fixed parameters occupying none — is part of the
    contract, so the oracle reconstructs it rather than asking ``free_slice``
    for it.
    """
    slices: list[tuple[Parameter, slice]] = []
    offset = 0
    for parameter in declaration.parameters:
        if parameter.is_fixed:
            continue
        slices.append((parameter, slice(offset, offset + parameter.size)))
        offset += parameter.size
    return slices


# ---------------------------------------------------------------------------
# Priors and the change of variables
# ---------------------------------------------------------------------------


def analytic_prior_transform(declaration: ParameterSet, cube: np.ndarray) -> np.ndarray:
    """Every free entry through its own ``scipy`` quantile function."""
    out = np.empty_like(np.asarray(cube, dtype=float))
    for parameter, where in free_slices(declaration):
        out[where] = parameter.prior.ppf(np.asarray(cube, dtype=float)[where])
    return out


def analytic_lnprior(declaration: ParameterSet, theta: np.ndarray) -> float:
    """Summed ``scipy`` ``logpdf`` over every free entry."""
    values = np.asarray(theta, dtype=float)
    return float(
        sum(
            np.sum(parameter.prior.logpdf(values[where]))
            for parameter, where in free_slices(declaration)
        )
    )


def _support(parameter: Parameter) -> tuple[float, float]:
    """The prior's support, refusing the one case ampere itself refuses.

    ``default_bijection_for`` raises for a support bounded above but not
    below, because there is no default map for it. An oracle that silently
    produced ``nan`` there would be worse than useless — it would agree with
    anything — so it raises too.
    """
    low, high = parameter.prior.support()
    if np.isneginf(low) and not np.isposinf(high):
        raise ValueError(
            f"the prior on {parameter.name!r} is bounded above but not below; "
            "default_bijection_for refuses this case and so does the oracle."
        )
    return float(low), float(high)


def analytic_constrain(parameter: Parameter, y: np.ndarray) -> np.ndarray:
    """``parameters.md`` §6's default bijection, chosen from the prior's support.

    Unbounded → the identity; bounded below only → ``lower + exp(y)``; bounded
    both ways → a logistic onto ``[lower, upper]``. That is the whole of
    ``default_bijection_for``'s table, written out.
    """
    low, high = _support(parameter)
    if np.isneginf(low) and np.isposinf(high):
        return y
    if np.isposinf(high):
        return low + np.exp(y)
    return low + (high - low) / (1.0 + np.exp(-y))


def analytic_log_abs_det(parameter: Parameter, y: np.ndarray) -> np.ndarray:
    """``log|dc/dy|`` for :func:`analytic_constrain`, evaluated at *y*.

    Derived, not copied: for the logistic ``x = l + (u - l) s(y)`` with
    ``s(y) = 1/(1 + e^-y)``, ``dx/dy = (u - l) s(y) (1 - s(y))``, and
    ``log s(y) = -logaddexp(0, -y)``, ``log(1 - s(y)) = -logaddexp(0, y)``.
    Always evaluated at the **unconstrained** point, which is what ampere's
    ``Bijection`` contract specifies.
    """
    low, high = _support(parameter)
    if np.isneginf(low) and np.isposinf(high):
        return np.zeros_like(y)
    if np.isposinf(high):
        return y
    return np.log(high - low) - np.logaddexp(0.0, -y) - np.logaddexp(0.0, y)


def summed_log_abs_det(declaration: ParameterSet, unconstrained: np.ndarray) -> float:
    """The Jacobian term ``lnprior_unconstrained`` adds to ``lnprior``."""
    y = np.asarray(unconstrained, dtype=float)
    return float(
        sum(
            np.sum(analytic_log_abs_det(parameter, y[where]))
            for parameter, where in free_slices(declaration)
        )
    )


# ---------------------------------------------------------------------------
# Likelihoods and kernels
# ---------------------------------------------------------------------------


def analytic_diagonal_gaussian_log_prob(residual: np.ndarray, sigma: np.ndarray) -> float:
    """The i.i.d. Gaussian log-likelihood, written out term by term."""
    return float(np.sum(-0.5 * ((residual / sigma) ** 2 + math.log(2.0 * math.pi)) - np.log(sigma)))


def matern32_matrix(coordinates: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """``k(r) = a^2 (1 + sqrt(3) r/L) exp(-sqrt(3) r/L)`` — the defining formula."""
    separation = np.abs(coordinates[:, None] - coordinates[None, :])
    scaled = math.sqrt(3.0) * separation / length_scale
    return amplitude**2 * (1.0 + scaled) * np.exp(-scaled)


def squared_exponential_matrix(
    coordinates: np.ndarray, amplitude: float, length_scale: float
) -> np.ndarray:
    """``k(r) = a^2 exp(-0.5 (r/L)^2)`` — the defining formula."""
    separation = np.abs(coordinates[:, None] - coordinates[None, :])
    return amplitude**2 * np.exp(-0.5 * (separation / length_scale) ** 2)


def kernel_matrix(
    family: KernelFamily, coordinates: np.ndarray, amplitude: float, length_scale: float
) -> np.ndarray:
    """The covariance matrix for *family*, from its defining formula."""
    if family is KernelFamily.MATERN32:
        return matern32_matrix(coordinates, amplitude, length_scale)
    return squared_exponential_matrix(coordinates, amplitude, length_scale)
