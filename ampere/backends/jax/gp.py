"""Kernels and the dense GP solve, in jax.

``inference.md`` §18 says a backend supplies "models and transformations […]
and nothing else", and kernels, solvers, families and containers are
backend-neutral and live in ``ampere.core``. That is true of the
*declarations*, which is what that sentence is about — and it is exactly why
the classes here **subclass** the core ones rather than redeclaring them: same
``FAMILY``, same ``HYPERPARAMETERS``, same ``QUASISEPARABLE`` flag, same
``NAME``. What a backend does supply is an implementation of the linear
algebra, and ``architecture.md`` §1 makes the GP solvers part of what
distinguishes one rung of the ladder from another.

What is here, and what is not
-----------------------------
:class:`Matern32` and :class:`SquaredExponential` — the two declared families,
with their covariance evaluated in ``jax.numpy`` so a gradient can be taken
through a fitted amplitude or length scale.

:class:`DenseGP` — the exact O(N³) Cholesky, in ``jax.scipy``. It has no
performance goal: it exists so that "what should the answer be?" has an
implementation on this backend too, and so the cross-backend rows have
something to compare.

``QuasisepGP`` is **not** here. ``DEVELOPMENT_PLAN.md`` §6 leaves the jax GP
library open, and slice 2 of this track chooses between ``tinygp``'s
``QuasisepSolver`` and ``celerite2.jax`` by measuring both against the
conformance suite. Declaring the slot now and filling it with a numpy solve
would satisfy every agreement row and prove nothing — which is precisely the
failure ``tests/conformance``'s "a declared quasiseparable strategy must be a
*different* strategy from the dense one" row exists to catch. The conformance
fixture therefore declares ``DENSE`` only, and the quasiseparable rows skip
with a reason naming what is owed.

What slice 1 measured, so slice 2 does not have to start from nothing
---------------------------------------------------------------------
The plan asks for the maintenance question to be settled at the track's start,
so here is the state on 2026-09-07, checked against PyPI and against the
installed environment rather than recollected:

* **celerite2 0.3.3** (released 2026-07-12) **ships ``celerite2.jax``**, with
  a ``GaussianProcess``, a ``terms`` module and its own ``ops``. celerite2 is
  already a **base** dependency of ampere (W2.3 put it there), so choosing it
  costs no new dependency at all — the strongest argument in its favour, and
  it would also mean the numpy and jax quasiseparable paths shared one
  library's arithmetic, which makes the cross-backend row a comparison of
  ampere's lowering rather than of two third-party solvers.
* **and it flips ``jax_enable_x64`` as an import side effect.** Importing
  ``celerite2.jax`` prints "celerite2.jax only works with dtype float64. We're
  enabling x64 now, but you might run into issues if you've already run some
  jax code" and calls ``config.update`` itself. That is exactly what
  ``lowering.md`` §10.2(a) forbids ampere from doing, done to ampere by a
  dependency: the guard in :mod:`ampere.backends.jax._config` would be
  satisfied by a flag nobody in the user's program set, and the "supported
  configurations are exactly two" promise would quietly become three. Slice 2
  must decide what to do about it — most likely importing it lazily, only
  after ampere's own guard has already refused or passed, so the side effect
  can never be what made the guard pass — and it is a point in ``tinygp``'s
  favour that has nothing to do with either library's numerics.
* **tinygp 0.3.1** (released 2026-03-15, ``requires-python >= 3.11``) is a
  live release but is **not** an existing dependency, so adopting it widens the
  ``jax`` extra. Its ``QuasisepSolver`` is the more idiomatic jax object (a
  pytree, composable with ``jit``/``vmap`` in the ordinary way), which matters
  for slice 2's batching work in a way it does not for correctness.

Neither has been benchmarked here and neither has been run against the
conformance suite; that is slice 2's work, and this note is the starting point
rather than the answer.

Cholesky failure is quiet on jax, and must not be
-------------------------------------------------
W2.3's finding about ``celerite2`` generalises to this backend, in a way worth
stating because the mechanism is different and the symptom is identical:
:func:`jax.scipy.linalg.cho_factor` does **not** raise on a non-positive-definite
matrix. It returns NaN, silently, and a NaN log-likelihood in a chain is not a
failed proposal — it is a number that poisons every diagnostic downstream of it
and attributes the damage to the sampler.

So both surfaces guard, and they guard differently because they must:

* :meth:`DenseGP.log_marginal_likelihood` — the contract surface, called on
  concrete values — checks the factor and raises
  :class:`~ampere.core.exceptions.LikelihoodError`, which is in
  ``FittingProblem``'s default catch set and therefore becomes ``-inf`` with a
  recorded reason (or propagates under ``strict=True``), exactly as the
  reference backend's does;
* :meth:`DenseGP.log_marginal_likelihood_jax` — the native, traced surface —
  cannot raise, because a Python ``if`` on a traced value is an error rather
  than a branch. It returns ``-inf`` computed with :func:`jax.numpy.where`,
  which is the same §4.5 answer arrived at by the only means available inside a
  trace. This is the reason ``inference.md`` §11's non-strict ruling exists.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import jax.scipy.linalg as jsl
import numpy as np

from ampere.core import GPConditional, GPSolver, Kernel
from ampere.core import Matern32 as _CoreMatern32
from ampere.core import SquaredExponential as _CoreSquaredExponential
from ampere.core.exceptions import LikelihoodError

from ._config import BACKEND, require_x64

__all__ = ["DenseGP", "Matern32", "SquaredExponential"]

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)


def _points(coordinates: Any) -> jax.Array:
    """``(n, d)`` coordinates from an ``(n,)`` or ``(n, d)`` array.

    A bare one-dimensional array is read as a column of ``n`` one-dimensional
    points, matching ``ampere.core.Kernel.matrix``'s own convention.
    """
    array = jnp.asarray(coordinates, dtype=jnp.float64)
    return array[:, None] if array.ndim == 1 else array


def _separation(left: Any, right: Any) -> jax.Array:
    """Euclidean separation between two coordinate sets, ``(n, m)``.

    Euclidean in the coordinate space, which is why
    :meth:`~ampere.core.GPSolver.check_compatible` requires every coordinate
    axis to share one unit — a single isotropic length scale is meaningless
    across mixed units.
    """
    a, b = _points(left), _points(right)
    difference = a[:, None, :] - b[None, :, :]
    return jnp.sqrt(jnp.sum(difference * difference, axis=-1))


class _JaxKernel(Kernel):
    """Shared plumbing for the jax kernels: jax separations, jax covariances.

    ``ampere.core.Kernel.matrix`` and ``.diagonal`` build their separations in
    numpy and, through ``_positive``, coerce each hyperparameter with
    ``float()`` — which is right on the reference path and fatal on a traced
    one. Both are overridden here rather than reused.

    The dropped ``_positive`` check is not dropped silently: a non-positive
    length scale gives a non-finite covariance, the Cholesky then yields NaN,
    and the guards in :class:`DenseGP` turn that into the §4.5 failure the
    contract asks for. On the reference path the check can be an exception
    because nothing is traced; here the same information has to travel as a
    value. The declaration is what keeps it from arising at all: a kernel
    hyperparameter is declared with a prior supported on the positive half-line
    and a ``Log`` bijection, so no sampler can propose a value outside it.
    """

    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(*args, **kwargs)

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> jax.Array:
        raise NotImplementedError

    def matrix(self, left: Any, right: Any, values: Mapping[str, Any]) -> jax.Array:
        """Dense covariance between two coordinate sets, ``(n, m)``, in jax."""
        return self._covariance(_separation(left, right), self.resolve(values))

    def diagonal(self, coordinates: Any, values: Mapping[str, Any]) -> jax.Array:
        """The prior variance at each coordinate; ``k(0)`` for a stationary kernel."""
        n = int(_points(coordinates).shape[0])
        return self._covariance(jnp.zeros(n, dtype=jnp.float64), self.resolve(values))


class Matern32(_JaxKernel, _CoreMatern32):
    r"""Matérn-3/2 in ``jax.numpy``: ampere's canonical flexible-likelihood kernel.

    .. math::
        k(r) = a^2 \left(1 + \frac{\sqrt{3}\,r}{\ell}\right)
               \exp\!\left(-\frac{\sqrt{3}\,r}{\ell}\right)

    ``amplitude`` is the marginal **standard deviation** — ``k(0) ==
    amplitude²`` — so a prior on it is a prior in the data's own units. The
    declaration, including ``QUASISEPARABLE = True``, is inherited from
    ``ampere.core.Matern32``: this kernel *does* have an exact rank-2
    quasiseparable representation, and saying otherwise here would be a lie
    about the mathematics merely because this backend has not yet shipped a
    solver that exploits it.
    """

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> jax.Array:
        amplitude = jnp.asarray(values["amplitude"], dtype=jnp.float64)
        length_scale = jnp.asarray(values["length_scale"], dtype=jnp.float64)
        scaled = _SQRT3 * jnp.asarray(separation, dtype=jnp.float64) / length_scale
        return amplitude * amplitude * (1.0 + scaled) * jnp.exp(-scaled)


class SquaredExponential(_JaxKernel, _CoreSquaredExponential):
    r"""Squared exponential (RBF) in ``jax.numpy``: legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Not quasiseparable, and inherits that declaration too, so a quasiseparable
    solver refuses it by name — the asymmetry that made Matérn the default.
    """

    def _covariance(self, separation: Any, values: Mapping[str, Any]) -> jax.Array:
        amplitude = jnp.asarray(values["amplitude"], dtype=jnp.float64)
        length_scale = jnp.asarray(values["length_scale"], dtype=jnp.float64)
        scaled = jnp.asarray(separation, dtype=jnp.float64) / length_scale
        return amplitude * amplitude * jnp.exp(-0.5 * scaled * scaled)


@dataclasses.dataclass(frozen=True)
class DenseGP(GPSolver):
    """Exact O(N³) dense Cholesky, in ``jax.scipy``.

    The same quantity ``ampere.core.DenseGP`` computes, by the same
    factorisation, in a different library — which is what makes the
    cross-backend agreement rows a test of the lowering rather than of the
    linear algebra. It has no performance goal.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal. Defaults to zero: a covariance that will not factorise is
        a fact about the model, and this contract does not hide it behind a
        silent numerical fudge.
    """

    jitter: float = 0.0

    NAME: ClassVar[str] = "DenseGP"
    EXACT: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"

    def __post_init__(self) -> None:
        # A frozen dataclass, and that is load-bearing rather than stylistic:
        # `Likelihood.describe()` records a solver's `config` only when the
        # solver is a dataclass (`dataclasses.fields`), so a backend whose
        # solver were a plain class would report a *differently shaped*
        # declaration for the same problem — and `test_cross_backend`'s
        # `ampere_likelihoods` row compares that shape across backends. The
        # core solver is a frozen dataclass with one `jitter` field; this one
        # matches it exactly.
        require_x64("the jax dense GP solver")
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(f"DenseGP's jitter must be finite and >= 0, got {self.jitter!r}.")

    # -- factorisation ------------------------------------------------------

    def _total(
        self,
        kernel: Kernel,
        coordinates: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        covariance = kernel.matrix(coordinates, coordinates, values)
        diagonal = jnp.asarray(variance, dtype=jnp.float64) + self.jitter**2
        return jnp.asarray(covariance, dtype=jnp.float64) + jnp.diag(diagonal)

    def _factor_jax(self, total: jax.Array) -> jax.Array:
        """The lower Cholesky factor. **NaN** where the matrix is not positive definite.

        Deliberately not guarded here: the two callers need the failure in two
        different forms (an exception on the contract path, ``-inf`` on the
        traced one), and a guard in the shared helper could only serve one of
        them.
        """
        return jnp.linalg.cholesky(total)

    def _checked_factor(
        self,
        kernel: Kernel,
        coordinates: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        """The contract path's factor: a :class:`LikelihoodError` where jax gives NaN."""
        total = self._total(kernel, coordinates, variance, values)
        if not bool(jnp.all(jnp.isfinite(total))):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        lower = self._factor_jax(total)
        if not bool(jnp.all(jnp.isfinite(lower))):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) is not positive definite, so its "
                "Cholesky factorisation failed. jax reports this as NaN rather than by raising, "
                "so ampere checks it here — a NaN log-likelihood in a chain is not a rejected "
                "proposal, it is a number that poisons every diagnostic downstream. Usual "
                "causes: a zero or near-duplicate observational uncertainty, coordinates closer "
                "together than float64 can separate at this length-scale, or a kernel amplitude "
                "far above the data scale. Pass DenseGP(jitter=...) — a standard deviation in "
                "the data's units — if the matrix is merely ill-conditioned rather than wrong."
            )
        return lower

    @staticmethod
    def _solve(lower: jax.Array, right: jax.Array) -> jax.Array:
        return jsl.cho_solve((lower, True), right)

    # -- the native, traced surface ----------------------------------------

    def log_marginal_likelihood_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        """``log N(residual; 0, K + diag(variance))`` as a traceable jax scalar.

        No exception control flow: a factorisation that fails returns ``-inf``
        through :func:`jax.numpy.where`, which is §4.5's failure signal computed
        the only way a traced function can compute it.
        """
        r = jnp.asarray(residual, dtype=jnp.float64)
        lower = self._factor_jax(self._total(kernel, coordinates, variance, values))
        alpha = self._solve(lower, r)
        log_determinant = 2.0 * jnp.sum(jnp.log(jnp.abs(jnp.diag(lower))))
        quadratic = jnp.sum(r * alpha)
        value = -0.5 * (quadratic + log_determinant + r.size * _LOG_2PI)
        return jnp.where(jnp.isfinite(value), value, -jnp.inf)

    # -- the contract surface ----------------------------------------------

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> float:
        """``log N(residual; 0, K(values) + diag(variance))``."""
        r = jnp.asarray(residual, dtype=jnp.float64)
        lower = self._checked_factor(kernel, coordinates, variance, values)
        alpha = self._solve(lower, r)
        log_determinant = 2.0 * float(jnp.sum(jnp.log(jnp.abs(jnp.diag(lower)))))
        quadratic = float(jnp.sum(r * alpha))
        return -0.5 * (quadratic + log_determinant + r.size * _LOG_2PI)

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """The closed form from the same Cholesky (Sundararajan & Keerthi 2001).

        With ``A = (K + diag(variance + jitter²))⁻¹``: ``sigma_i^{2,-i} = 1 / A_ii``
        and ``μ_i^{-i} = y_i - [A r]_i / A_ii``, so
        ``log p_i = ½ log A_ii - [A r]_i² / (2 A_ii) - ½ log 2π``. These are
        what ``arviz.loo``/``waic`` consume; they are a *different*
        decomposition from the factorised pointwise terms of independent noise
        and do not sum to the joint log-likelihood.
        """
        r = jnp.asarray(residual, dtype=jnp.float64)
        lower = self._checked_factor(kernel, coordinates, variance, values)
        alpha = self._solve(lower, r)
        precision_diagonal = jnp.diag(self._solve(lower, jnp.eye(r.size, dtype=jnp.float64)))
        return np.asarray(
            0.5 * jnp.log(precision_diagonal)
            - alpha**2 / (2.0 * precision_diagonal)
            - 0.5 * _LOG_2PI
        )

    def condition(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
        at: Any = None,
    ) -> GPConditional:
        """The GP posterior at *at* (default: the data coordinates)."""
        points = _points(coordinates)
        lower = self._checked_factor(kernel, points, variance, values)
        target = points if at is None else _points(at)
        cross = jnp.asarray(kernel.matrix(target, points, values), dtype=jnp.float64)
        mean = cross @ self._solve(lower, jnp.asarray(residual, dtype=jnp.float64))
        solved = self._solve(lower, cross.T)
        prior_variance = jnp.asarray(kernel.diagonal(target, values), dtype=jnp.float64)
        posterior = prior_variance - jnp.einsum("ij,ji->i", cross, solved)
        return GPConditional(mean=np.asarray(mean), variance=np.asarray(posterior))

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """Map whitened latent draws ``z`` to a GP draw ``f = L(θ) z``.

        The deterministic half of the latent-GP declaration: the *prior* stays
        i.i.d. standard normal and all the correlation lives here, where the
        solver can impose it in whatever representation it uses.
        """
        points = _points(coordinates)
        covariance = jnp.asarray(kernel.matrix(points, points, values), dtype=jnp.float64)
        if not bool(jnp.all(jnp.isfinite(covariance))):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        scale = float(jnp.mean(jnp.diag(covariance))) or 1.0
        stabilised = covariance + jnp.eye(covariance.shape[0], dtype=jnp.float64) * (jitter * scale)
        lower = self._factor_jax(stabilised)
        if not bool(jnp.all(jnp.isfinite(lower))):
            raise LikelihoodError(
                "the kernel matrix K is not positive definite, so the whitening transform "
                "f = L z is undefined (jax reports this as NaN rather than by raising). "
                "Increase the jitter argument, or check the kernel hyperparameters."
            )
        return np.asarray(lower @ jnp.asarray(whitened, dtype=jnp.float64))
