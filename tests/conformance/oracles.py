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

from .protocol import CovarianceSpec, KernelFamily

__all__ = [
    "DISPERSED_AXES",
    "analytic_constrain",
    "analytic_diagonal_gaussian_log_prob",
    "analytic_lnprior",
    "analytic_log_abs_det",
    "analytic_prior_transform",
    "covariance_matrix",
    "free_slices",
    "kernel_matrix",
    "matern12_matrix",
    "matern32_matrix",
    "matern52_matrix",
    "rotation_matrix",
    "sho_matrix",
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


def matern12_matrix(coordinates: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """``k(r) = a^2 exp(-r/L)`` — the defining formula (W4.5)."""
    separation = _separations(coordinates)
    return amplitude**2 * np.exp(-separation / length_scale)


def matern32_matrix(coordinates: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """``k(r) = a^2 (1 + sqrt(3) r/L) exp(-sqrt(3) r/L)`` — the defining formula."""
    separation = _separations(coordinates)
    scaled = math.sqrt(3.0) * separation / length_scale
    return amplitude**2 * (1.0 + scaled) * np.exp(-scaled)


def matern52_matrix(coordinates: np.ndarray, amplitude: float, length_scale: float) -> np.ndarray:
    """``k(r) = a^2 (1 + sqrt(5) r/L + 5 r^2/3L^2) exp(-sqrt(5) r/L)`` (W4.5)."""
    separation = _separations(coordinates)
    scaled = math.sqrt(5.0) * separation / length_scale
    return amplitude**2 * (1.0 + scaled + scaled**2 / 3.0) * np.exp(-scaled)


def squared_exponential_matrix(
    coordinates: np.ndarray, amplitude: float, length_scale: float
) -> np.ndarray:
    """``k(r) = a^2 exp(-0.5 (r/L)^2)`` — the defining formula."""
    separation = _separations(coordinates)
    return amplitude**2 * np.exp(-0.5 * (separation / length_scale) ** 2)


def _celerite_pair_matrix(
    separation: np.ndarray, cosine: float, sine: float, decay: float, frequency: float
) -> np.ndarray:
    """``e^{-c r}(a cos d r + b sin d r)`` — one celerite term, written out."""
    return np.exp(-decay * separation) * (
        cosine * np.cos(frequency * separation) + sine * np.sin(frequency * separation)
    )


def sho_matrix(
    coordinates: np.ndarray, amplitude: float, period: float, quality: float
) -> np.ndarray:
    """celerite2's underdamped ``SHOTerm``, from its coefficients (W4.5).

    ``a = sigma^2``, ``b = a / sqrt(4Q^2 - 1)``, ``c = w0 / 2Q``,
    ``d = c sqrt(4Q^2 - 1)`` with ``w0 = 2 pi / period`` — so ``k(0) = a`` and
    the amplitude is a marginal standard deviation, as it is everywhere else in
    ampere.
    """
    separation = _separations(coordinates)
    omega0 = 2.0 * math.pi / period
    split = math.sqrt(4.0 * quality**2 - 1.0)
    decay = 0.5 * omega0 / quality
    return _celerite_pair_matrix(
        separation, amplitude**2, amplitude**2 / split, decay, decay * split
    )


def rotation_matrix(
    coordinates: np.ndarray,
    amplitude: float,
    period: float,
    quality: float,
    delta_quality: float,
    fraction: float,
) -> np.ndarray:
    """celerite2's ``RotationTerm``: two SHOs, at a period and its harmonic (W4.5)."""
    separation = _separations(coordinates)
    power = amplitude**2 / (1.0 + fraction)
    total = np.zeros_like(separation)
    for factor, harmonic, share in (
        (0.5 + quality + delta_quality, 1.0, 1.0),
        (0.5 + quality, 2.0, fraction),
    ):
        split = math.sqrt(4.0 * factor**2 - 1.0)
        omega = harmonic * 4.0 * math.pi * factor / (period * split)
        decay = 0.5 * omega / factor
        cosine = share * power
        total = total + _celerite_pair_matrix(
            separation, cosine, cosine / split, decay, decay * split
        )
    return total


def _separations(coordinates: np.ndarray) -> np.ndarray:
    """Euclidean separations of an ``(n,)`` or ``(n, d)`` coordinate set."""
    points = np.asarray(coordinates, dtype=float)
    if points.ndim == 1:
        points = points[:, None]
    difference = points[:, None, :] - points[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", difference, difference))


def kernel_matrix(
    family: KernelFamily,
    coordinates: np.ndarray,
    amplitude: float,
    length_scale: float,
) -> np.ndarray:
    """The covariance matrix for a leaf *family*, from its defining formula.

    Kept at its W1.10 signature — two hyperparameters, positional — because
    every row written before W4.5 calls it that way. The oscillators and the
    composites go through :func:`covariance_matrix`, which takes the whole
    :class:`~tests.conformance.protocol.CovarianceSpec` and recurses.
    """
    if family is KernelFamily.MATERN12:
        return matern12_matrix(coordinates, amplitude, length_scale)
    if family is KernelFamily.MATERN52:
        return matern52_matrix(coordinates, amplitude, length_scale)
    if family is KernelFamily.MATERN32:
        return matern32_matrix(coordinates, amplitude, length_scale)
    return squared_exponential_matrix(coordinates, amplitude, length_scale)


def covariance_matrix(spec: CovarianceSpec, coordinates: np.ndarray) -> np.ndarray:
    """The covariance matrix a :class:`CovarianceSpec` declares (W4.5).

    The oracle for every W4.5 row: built from the defining formulae here, never
    by calling a kernel, so a row compares an implementation against an
    equation rather than against itself. ``coordinates`` is ``(n,)`` or
    ``(n, d)``; a spec with ``axes`` selects its own columns from the latter by
    name, which is what makes the product-across-axes row an oracle and not a
    restatement.
    """
    points = np.asarray(coordinates, dtype=float)
    if points.ndim == 1:
        points = points[:, None]
    if spec.family in (KernelFamily.SUM, KernelFamily.SPECTRAL_MIXTURE):
        blocks = [covariance_matrix(term, points) for term in spec.terms]
        return np.sum(blocks, axis=0)
    if spec.family is KernelFamily.PRODUCT:
        total = covariance_matrix(spec.terms[0], points)
        for term in spec.terms[1:]:
            total = total * covariance_matrix(term, points)
        return total
    selected = points if spec.axes is None else points[:, list(_columns(spec, points))]
    if spec.family is KernelFamily.SHO:
        return sho_matrix(selected, spec.amplitude, spec.period, spec.quality)
    if spec.family is KernelFamily.ROTATION:
        return rotation_matrix(
            selected,
            spec.amplitude,
            spec.period,
            spec.quality,
            spec.delta_quality,
            spec.fraction,
        )
    if spec.family is KernelFamily.USER:
        # The user term of the battery's registry row is a re-labelled
        # Matérn-1/2, so its oracle is Matérn-1/2's.
        return matern12_matrix(selected, spec.amplitude, spec.length_scale)
    return kernel_matrix(spec.family, selected, spec.amplitude, spec.length_scale)


#: Axis order of the battery's three-axis container, for :func:`covariance_matrix`'s
#: ``axes`` selection. Declared here rather than imported so the oracle module
#: keeps importing nothing from the rows it is an oracle for.
DISPERSED_AXES: tuple[str, ...] = ("u", "v", "spectral_axis")


def _columns(spec: CovarianceSpec, points: np.ndarray) -> tuple[int, ...]:
    names = DISPERSED_AXES[: points.shape[1]]
    assert spec.axes is not None
    return tuple(names.index(name) for name in spec.axes)
