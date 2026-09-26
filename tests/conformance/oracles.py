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
import scipy.special

from ampere.core import Parameter, ParameterSet

from .protocol import CovarianceSpec, KernelFamily

#: Milliarcseconds in one radian. Written out rather than imported from the
#: backend under test: an oracle that asked ampere for its own constants would
#: agree with ampere about them by construction.
MAS_PER_RAD = 180.0 * 3600.0 * 1000.0 / math.pi

#: Gaussian FWHM in units of its standard deviation.
FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))

#: Julian days per year, for the astrometric proper-motion term (W4.9). A
#: second transcription of the same number
#: :mod:`ampere.backends.reference.astrometry` uses, deliberately: an oracle
#: that imported ampere's own constant would agree with ampere about it by
#: construction.
DAYS_PER_YEAR = 365.25

__all__ = [
    "DAYS_PER_YEAR",
    "DISPERSED_AXES",
    "FWHM_PER_SIGMA",
    "MAS_PER_RAD",
    "analytic_constrain",
    "analytic_diagonal_gaussian_log_prob",
    "analytic_lnprior",
    "analytic_log_abs_det",
    "analytic_prior_transform",
    "binary_closure_phase",
    "binary_visibility",
    "covariance_matrix",
    "free_slices",
    "gaussian_visibility",
    "kernel_matrix",
    "matern12_matrix",
    "matern32_matrix",
    "matern52_matrix",
    "reflex_orbit_dec",
    "reflex_orbit_ra",
    "rotation_matrix",
    "sho_matrix",
    "squared_exponential_matrix",
    "summed_log_abs_det",
    "uniform_disc_visibility",
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


# ---------------------------------------------------------------------------
# Interferometry (W4.1)
# ---------------------------------------------------------------------------
#
# The three closed forms every interferometric row compares against, and the
# closure phase derived from the third. They are transcribed from the standard
# statements of them — van Cittert-Zernike applied to a top hat, to a Gaussian
# and to a pair of delta functions — with this suite's own sign convention
# written in explicitly (``exp(-2 pi i (u x + v y))``, ``x`` east, ``y`` north,
# a position angle measured east of north), because the sign is exactly what a
# closure-phase row exists to catch.


def uniform_disc_visibility(
    u_pts: np.ndarray, v_pts: np.ndarray, *, diameter: float, flux: float
) -> np.ndarray:
    """``flux 2 J1(pi theta rho) / (pi theta rho)`` for a disc of *diameter* mas."""
    rho = np.hypot(np.asarray(u_pts, dtype=float), np.asarray(v_pts, dtype=float))
    argument = math.pi * diameter / MAS_PER_RAD * rho
    envelope = np.where(argument == 0.0, 1.0, 2.0 * scipy.special.j1(argument) / argument)
    return (flux * envelope).astype(np.complex128)


def gaussian_visibility(
    u_pts: np.ndarray, v_pts: np.ndarray, *, fwhm: float, flux: float
) -> np.ndarray:
    """``flux exp(-2 pi**2 sigma**2 rho**2)`` for a circular Gaussian of *fwhm* mas."""
    rho = np.hypot(np.asarray(u_pts, dtype=float), np.asarray(v_pts, dtype=float))
    sigma = fwhm / FWHM_PER_SIGMA / MAS_PER_RAD
    return (flux * np.exp(-2.0 * math.pi**2 * sigma**2 * rho**2)).astype(np.complex128)


def binary_visibility(
    u_pts: np.ndarray,
    v_pts: np.ndarray,
    *,
    separation: float,
    position_angle: float,
    flux_ratio: float,
    flux: float,
    component_fwhm: float,
) -> np.ndarray:
    """Two Gaussian components: the primary at the phase centre, the secondary offset.

    ``separation`` and ``component_fwhm`` are mas, ``position_angle`` radians
    east of north, so the secondary sits at
    ``(s sin p, s cos p)``. The ``exp(-2 pi i)`` sign is the suite's.
    """
    offset_x = separation * math.sin(position_angle) / MAS_PER_RAD
    offset_y = separation * math.cos(position_angle) / MAS_PER_RAD
    phase = -2j * math.pi * (np.asarray(u_pts) * offset_x + np.asarray(v_pts) * offset_y)
    envelope = gaussian_visibility(u_pts, v_pts, fwhm=component_fwhm, flux=1.0)
    return envelope * flux / (1.0 + flux_ratio) * (1.0 + flux_ratio * np.exp(phase))


def binary_closure_phase(
    u1: np.ndarray, v1: np.ndarray, u2: np.ndarray, v2: np.ndarray, **source: float
) -> np.ndarray:
    """``arg(V_ij V_jk V_ki)`` with ``ki`` implied as ``-(ij + jk)``.

    The canonical ordering ``ClosurePhases`` fixes, applied to
    :func:`binary_visibility`. The Gaussian envelope is real and positive, so
    it drops out of the argument — which is the reason a closure phase is a
    statement about the *geometry* of a source and not about its size.
    """
    third = (-(np.asarray(u1) + np.asarray(u2)), -(np.asarray(v1) + np.asarray(v2)))
    product = (
        binary_visibility(u1, v1, **source)
        * binary_visibility(u2, v2, **source)
        * binary_visibility(third[0], third[1], **source)
    )
    return np.angle(product)


def reflex_orbit_ra(
    time: np.ndarray, *, pmra: float, period: float, phase: float, amp_ra: float
) -> np.ndarray:
    """The closed-form right-ascension offset of a circular reflex orbit, mas.

    ``time`` in days: a linear proper-motion drift plus a periodic wobble,
    :mod:`ampere.backends.reference.astrometry`'s own formula, transcribed
    independently (W4.9's oracle is an equation, not a call to the model).
    """
    t = np.asarray(time, dtype=float)
    cycle = 2.0 * math.pi * t / period + phase
    return pmra * t / DAYS_PER_YEAR + amp_ra * np.sin(cycle)


def reflex_orbit_dec(
    time: np.ndarray, *, pmdec: float, period: float, phase: float, amp_dec: float
) -> np.ndarray:
    """The closed-form declination offset of a circular reflex orbit, mas.

    The cosine partner of :func:`reflex_orbit_ra`, on the same ``period`` and
    ``phase`` — the two coordinates of one orbit, sharing one cycle.
    """
    t = np.asarray(time, dtype=float)
    cycle = 2.0 * math.pi * t / period + phase
    return pmdec * t / DAYS_PER_YEAR + amp_dec * np.cos(cycle)


# ---------------------------------------------------------------------------
# The latent GP on closure phases (W5.1)
# ---------------------------------------------------------------------------


def von_mises_latent_log_likelihood(
    observed: np.ndarray,
    predicted: np.ndarray,
    sigma: np.ndarray,
    covariance: np.ndarray,
    whitened: np.ndarray,
    *,
    jitter: float = 1e-10,
) -> float:
    """``sum log VonMises(observed | predicted + L z, 1/sigma**2)``, written out.

    The oracle for the von Mises latent composition, which has no closed-form
    marginal to compare against: a GP added to a wrapped observable is only
    defined given the latent draw, so the oracle is the draw itself —
    ``f = L z`` with ``L L^T = K`` (plus the relative stabiliser every
    ``DenseGP.latent_transform`` adds, ``jitter * mean(diag K)``) — and the
    normalised von Mises density around ``predicted + f``. The wrap needs no
    code here: ``cos`` is periodic, and ``log I0(kappa)`` is taken as
    ``log(i0e(kappa)) + kappa`` so that a well-measured triangle's large
    ``kappa`` does not overflow.
    """
    scale = float(np.mean(np.diag(covariance))) or 1.0
    lower = np.linalg.cholesky(covariance + np.eye(covariance.shape[0]) * (jitter * scale))
    latent = lower @ np.asarray(whitened, dtype=float)
    kappa = 1.0 / np.asarray(sigma, dtype=float) ** 2
    delta = np.asarray(observed, dtype=float) - np.asarray(predicted, dtype=float) - latent
    log_i0 = np.log(scipy.special.i0e(kappa)) + kappa
    return float(np.sum(kappa * np.cos(delta) - math.log(2.0 * math.pi) - log_i0))


def rotation_coupling_matrix(
    angle: float, log_variance_0: float, log_variance_1: float
) -> np.ndarray:
    """``B = R(θ) diag(exp(v₀), exp(v₁)) R(θ)ᵀ`` — ``RotationCoupling``'s definition (W5.24).

    Written from the parameterisation's statement (an error ellipse at
    position angle ``θ`` with semi-axis variances ``exp(v)``), not by asking
    the coupling for its matrix.
    """
    cosine, sine = math.cos(angle), math.sin(angle)
    rotation = np.array([[cosine, -sine], [sine, cosine]])
    return rotation @ np.diag([math.exp(log_variance_0), math.exp(log_variance_1)]) @ rotation.T


def coregionalised_covariance(
    coupling: np.ndarray, kernel: np.ndarray, variances: np.ndarray
) -> np.ndarray:
    """``B ⊗ K_x + blockdiag(diag(sigma_t^2))``, channel-major, block by block (W5.24).

    The definition of the joint noise model's covariance under heteroscedastic
    channels, assembled by writing every ``(s, t)`` block out rather than by
    ``np.kron``: block ``(s, t)`` is ``B[s, t] K_x``, and diagonal block ``t``
    adds channel ``t``'s own variances, ``variances[:, t]``.
    """
    channels = coupling.shape[0]
    size = kernel.shape[0]
    covariance = np.zeros((channels * size, channels * size))
    for row in range(channels):
        for column in range(channels):
            block = coupling[row, column] * kernel
            if row == column:
                block = block + np.diag(variances[:, row])
            covariance[row * size : (row + 1) * size, column * size : (column + 1) * size] = block
    return covariance


def coregionalised_log_density(
    coupling: np.ndarray, kernel: np.ndarray, variances: np.ndarray, residuals: np.ndarray
) -> float:
    """``log N(vec(R); 0, B ⊗ K_x + blockdiag(diag(sigma_t^2)))`` by ``scipy`` (W5.24).

    *residuals* is the ``(n, T)`` block; ``vec`` stacks it channel-major, the
    same order :func:`coregionalised_covariance` builds the matrix in. The
    density is ``scipy.stats.multivariate_normal``'s, so neither the
    factorisation nor the assembly is ampere's.
    """
    import scipy.stats

    covariance = coregionalised_covariance(coupling, kernel, variances)
    stacked = np.asarray(residuals, dtype=float).T.reshape(-1)
    return float(
        scipy.stats.multivariate_normal(np.zeros(stacked.size), covariance).logpdf(stacked)
    )
