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

:class:`QuasisepGP` — the exact O(N) quasiseparable solve, over
**celerite2.jax**, chosen by measurement in slice 2 (``DEVELOPMENT_PLAN.md``
§2's "jax quasiseparable GP library" row, and §6's deferred-choice bullet).

The library choice, and what it cost
------------------------------------
``DEVELOPMENT_PLAN.md`` §6 left the jax GP library open with an instruction to
verify at Phase 2's start which of ``tinygp``'s ``QuasisepSolver`` and
``celerite2.jax`` was the better bet. Slice 1 surveyed both; slice 2
implemented ``QuasisepGP`` over each, against ``ampere.core``'s exact rank-2
Matérn-3/2 representation, and measured. The full table is in the decision-log
row; the two numbers that decided it are these.

**Scaling.** ``celerite2.jax``'s recursions are XLA custom calls into
celerite2's own compiled kernels, and they are linear in N with a constant of
about 0.25 µs per point: 10⁴ points in 2.3 ms, 10⁶ in 0.26 s, a value *and*
gradient at 10⁶ in 0.75 s. ``tinygp``'s ``QuasisepSolver`` is written in jax
itself, as ``lax.scan`` recursions over 2-by-2 blocks, and on XLA's CPU backend
it is **quadratic in practice** — 3.0 µs per point at N = 500, 6.4 at 10³,
18.1 at 4e3, 98.0 at 1.6e4, doubling with every doubling of N. At 10⁴
points that is 0.47 s for one likelihood evaluation against celerite2's
2.3 ms, a factor of 200, and the gap widens with the data. A solver
``DEVELOPMENT_PLAN.md`` §2 calls "the scaling answer" cannot be the one that
stops scaling at 10³ points.

**Accuracy.** Both agree with the dense Cholesky far inside
``tolerances.cross_solver`` (1e-6), and both agree with
``scipy.stats.multivariate_normal`` to 2e-11 on a log-likelihood of 4e3.
tinygp is the more accurate of the two — 1e-13 to 1e-12 whatever the span,
because its state-space form composes *local* transitions — while celerite2
inherits the cancellation ``ampere.core``'s own term docstring describes: the
rank-2 generators grow linearly in the coordinate, so their products are
differences of large numbers and the error grows with the number of length
scales the data span (1.6e-10 over 10, 1.3e-09 over 10², 2.7e-08 over 10³).
Midpoint-centring the coordinates, which this backend inherits from the core's
representation, is what bounds it; 2.7e-08 on a log-likelihood of order 10³ is
two orders inside the tolerance and far inside anything a sampler can tell
apart.

Three things celerite2.jax costs, all recorded here rather than hidden:

* **it flips ``jax_enable_x64`` as an import side effect** — the very thing
  ``lowering.md`` §10.2(a) forbids ampere from doing. So the import is
  **deferred** to first use, inside :func:`_celerite2_jax`, which can only run
  after :func:`~ampere.backends.jax._config.require_x64` has already passed at
  the solver's construction. The flag is therefore always on before celerite2
  looks at it, its ``config.update`` is a no-op and its warning never fires:
  the side effect can never be what made ampere's guard pass.
  :func:`_celerite2_jax` asserts that invariant rather than trusting it.
* **it has no ``vmap`` batching rule** — ``jax.vmap`` over a density
  containing this solver raises ``NotImplementedError: Batching rule for
  'celerite2_factor' not implemented``. So this class declares
  ``BATCHABLE = False`` where the rest of the backend declares ``True``, and
  ``ampere.core.declared_capabilities`` aggregates that down to a problem
  which honestly says it cannot be batched. tinygp, being ordinary jax, vmaps.
* **it exposes no O(N) route to the diagonal of ``(K + diag(σ²))⁻¹``**, so
  :meth:`QuasisepGP.conditional_loo` stays refused here exactly as it is on
  the reference path (``DEVELOPMENT_PLAN.md`` §2, 2026-09-05). tinygp *does*:
  ``L = solver.factor``, ``L.inv().transpose() @ L.inv()`` is a quasiseparable
  matrix whose ``.diag.d`` is that diagonal, measured here exact to 7e-15
  against a dense inverse. That is a real loss, and it is the price of the
  factor of 200.

**Cost as a dependency.** celerite2 is already a *base* dependency (W2.3), so
choosing it adds nothing at all to the ``jax`` extra; tinygp would have added
one package, needing only jax and equinox, both already there. Both are live
releases — celerite2 0.3.3 (2026-07-12), tinygp 0.3.1 (2026-03-15) — and
neither is unmaintained, so maintenance did not decide this.

Choosing celerite2 also means the numpy and jax quasiseparable paths share one
library's arithmetic over one representation ampere wrote twice, which makes
``TestSolverAgreement``'s cross-solver rows a test of ampere's lowering rather
than of two third-party solvers. That is a smaller argument than the scaling
one, and it points the same way.

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
import functools
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

from ._config import BACKEND, require_x64, x64_enabled

__all__ = ["DenseGP", "Matern32", "QuasisepGP", "SquaredExponential"]

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)

#: The default device for every part of this backend. **Chosen, never
#: detected** (``architecture.md`` §5): a library that silently moved a fit
#: onto whatever accelerator it found would make "why is this run slow?" and
#: "why did this run give a different answer?" both unanswerable.
DEVICE = "cpu"

#: The precisions :func:`_solve_dtype` accepts, mapped to their jax dtypes.
#:
#: ``architecture.md`` §5 makes float64 the policy for every likelihood and GP
#: solve and reserves a **per-run opt-out** for GPU throughput. On jax that
#: opt-out is a deliberate departure from the x64 policy ``lowering.md`` §10.2
#: otherwise enforces process-wide, so it is narrow (this solve only), explicit
#: (a constructor argument, never a default) and recorded
#: (:meth:`DenseGP.provenance_config`). A GP Cholesky in float32 fails in ways
#: that read as science problems — ``DEVELOPMENT_PLAN.md`` §7's trap — which is
#: why it is opt-in and why the run says it happened.
PRECISIONS: dict[str, Any] = {"float64": jnp.float64, "float32": jnp.float32}


def _solve_dtype(precision: str, owner: str) -> Any:
    """The jax dtype for *precision*, or a refusal naming the two that exist."""
    try:
        return PRECISIONS[str(precision)]
    except KeyError:
        known = ", ".join(sorted(PRECISIONS))
        raise LikelihoodError(
            f"{owner} was asked for precision {precision!r}; the choices are {known}. float64 is "
            f"the policy (architecture.md §5) and float32 is the recorded per-run opt-out for "
            f"throughput — a GP factorisation in float32 fails in ways that look like science "
            f"problems, so it is never a default."
        ) from None


def _resolve_device(device: str, owner: str) -> Any:
    """A jax device for *device*, by platform name. Never auto-detected.

    ``architecture.md`` §5: ampere chooses nothing here. The name is looked up
    among the devices jax actually has in this process, and an absent platform
    is a refusal naming what is present rather than a silent fallback to the
    CPU — a fit the user asked to run on a GPU and which quietly did not is
    the failure this rule exists to prevent.

    Returns the first device of that platform. Picking *which* GPU is the
    user's business, and is done by passing a ``jax.Device`` straight through.
    """
    if not isinstance(device, str):
        return device  # an explicit jax.Device, used as given
    available = jax.devices()
    for found in available:
        if found.platform == device:
            return found
    platforms = ", ".join(sorted({d.platform for d in available}))
    raise LikelihoodError(
        f"{owner} was asked for device {device!r}, but this jax process has no such platform. "
        f"Available: {platforms}. ampere never falls back to another device — a fit that was "
        f"asked for a GPU and quietly ran on the CPU is a fit whose timings mean nothing."
    )


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
    #: **True since slice 2** (W2.5). The flag means one thing on this
    #: backend: ``jax.vmap`` over the realised density
    #: (:meth:`~ampere.backends.jax.problem.LoweredProblem.log_prob_unconstrained_batched`)
    #: evaluates a stack of parameter vectors in one call, and it is measured
    #: rather than asserted -- ``tests/backends/test_jax.py`` compares a vmapped
    #: density against the same density in a loop. It is true here because every
    #: operation in this class is whole-array ``jax.numpy``: nothing branches on
    #: a value, nothing indexes by one, so vmap maps it as it maps any pure
    #: function. ``QuasisepGP`` is the one part of this backend that still says
    #: False, and says why.
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE

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
    precision
        ``"float64"`` (the policy, and the default) or ``"float32"`` — the
        per-run opt-out ``architecture.md`` §5 reserves for GPU throughput,
        applied to **this solve only**. See :data:`PRECISIONS` and
        :meth:`provenance_config`; it is an ``InitVar`` rather than a field so
        that it stays out of the spec hash.
    device
        Platform name (``"cpu"``, the default) or an explicit ``jax.Device``.
        **Never auto-detected** (``architecture.md`` §5). It sets this
        instance's ``DEVICE`` capability flag, so a solver placed elsewhere
        than the models beside it is a device disagreement
        ``ampere.core.declared_capabilities`` refuses at composition, loudly —
        which is the intended behaviour and also the reason a GPU run needs
        the same keyword on the models and instrument steps. Adding it there
        is what a later GPU item does; the shape it has to fit is this one.
        Also an ``InitVar``, for the same spec-hash reason as *precision*.

    Batching
    --------
    ``BATCHABLE = True`` since slice 2: every operation here is
    ``jnp``/``jax.scipy`` on whole arrays, so ``jax.vmap`` over a batch of
    parameter vectors maps the factorisation as it maps everything else, and
    :func:`~ampere.backends.jax.problem.lower_problem`'s
    ``log_prob_unconstrained_batched`` is the surface that uses it. The flag
    is measured rather than asserted — ``tests/backends/test_jax.py`` vmaps a
    GP density and compares against the loop.
    """

    jitter: float = 0.0
    precision: dataclasses.InitVar[str] = "float64"
    device: dataclasses.InitVar[Any] = DEVICE

    NAME: ClassVar[str] = "DenseGP"
    EXACT: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE

    def __post_init__(self, precision: str, device: Any) -> None:
        # A frozen dataclass, and that is load-bearing rather than stylistic:
        # `Likelihood.describe()` records a solver's `config` only when the
        # solver is a dataclass (`dataclasses.fields`), so a backend whose
        # solver were a plain class would report a *differently shaped*
        # declaration for the same problem — and `test_cross_backend`'s
        # `ampere_likelihoods` row compares that shape across backends. The
        # core solver is a frozen dataclass with one `jitter` field; this one
        # matches it exactly. `precision` is an InitVar precisely so that it
        # does *not* join that list: `dataclasses.fields()` omits InitVars, so
        # the opt-out is invisible to the spec hash and visible in
        # `provenance_config()`, which is fold-in 10's whole point.
        require_x64("the jax dense GP solver")
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(f"DenseGP's jitter must be finite and >= 0, got {self.jitter!r}.")
        object.__setattr__(self, "_dtype", _solve_dtype(precision, self.NAME))
        object.__setattr__(self, "_precision", str(precision))
        resolved = _resolve_device(device, self.NAME)
        object.__setattr__(self, "_device", resolved)
        # The capability flag is per *instance* here rather than per class, so
        # that a solver placed on an accelerator says so and
        # `declared_capabilities` can catch it beside CPU models.
        object.__setattr__(
            self, "DEVICE", device if isinstance(device, str) else str(resolved.platform)
        )

    @property
    def precision_name(self) -> str:
        """``"float64"`` or ``"float32"``: which precision this solve runs in."""
        return str(getattr(self, "_precision", "float64"))

    @property
    def dtype(self) -> Any:
        """The jax dtype :meth:`_total` factorises in."""
        return getattr(self, "_dtype", jnp.float64)

    def place(self, array: Any) -> jax.Array:
        """*array* on this solver's device, as its solve dtype."""
        return jax.device_put(jnp.asarray(array, dtype=self.dtype), getattr(self, "_device", None))

    def provenance_config(self) -> Mapping[str, Any]:
        """``inference.md`` §10a fold-in 10: how this solver computes, for the attrs.

        Precision and device are configuration, not declaration: they change
        how the same declared model is computed, two backends legitimately
        differ on them, and a dtype in the spec hash would break
        ``results.md`` §14's cross-backend agreement. So they are reported
        here, recorded under ``ampere_solver_config``, and hashed nowhere.

        ``architecture.md`` §5 requires the float32 opt-out to be *visible in
        provenance*, and on jax it is more than a dtype: it is a deliberate,
        per-run departure from the x64 policy ``lowering.md`` §10.2 otherwise
        enforces for the whole process. So a run that took it says so.
        """
        return {
            "library": "jax",
            "dtype": self.precision_name,
            "device": self.DEVICE,
            "x64_policy_opt_out": self.precision_name != "float64",
        }

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
        total = jnp.asarray(covariance, dtype=jnp.float64) + jnp.diag(diagonal)
        # The opt-out applies to the *solve*, and only to the solve: the
        # kernel and the noise quadrature are still built in float64, and the
        # scalar that comes back out is float64 again. Narrowing here is what
        # makes it "the GP linear algebra" rather than "the whole likelihood".
        return self.place(total)

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
        r = self.place(residual)
        lower = self._factor_jax(self._total(kernel, coordinates, variance, values))
        alpha = self._solve(lower, r)
        log_determinant = 2.0 * jnp.sum(jnp.log(jnp.abs(jnp.diag(lower))))
        quadratic = jnp.sum(r * alpha)
        value = jnp.asarray(
            -0.5 * (quadratic + log_determinant + r.size * _LOG_2PI), dtype=jnp.float64
        )
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
        r = self.place(residual)
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
        r = self.place(residual)
        lower = self._checked_factor(kernel, coordinates, variance, values)
        alpha = self._solve(lower, r)
        precision_diagonal = jnp.diag(self._solve(lower, jnp.eye(r.size, dtype=self.dtype)))
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

    def latent_transform_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> jax.Array:
        """:meth:`latent_transform` without the exception control flow.

        The native surface a latent-GP dataset needs *inside* the traced
        density: :mod:`ampere.backends.jax.problem` samples the whitened
        vector ``z`` as an ordinary i.i.d. standard-normal block and this
        turns it into ``f = L(θ) z``, so the correlation is imposed where a
        gradient can still reach the hyperparameters. A factorisation that
        fails gives NaN, which propagates to the log-density and is mapped to
        ``-inf`` there — the same §4.5 answer by the only means a trace has.
        """
        points = _points(coordinates)
        covariance = jnp.asarray(kernel.matrix(points, points, values), dtype=jnp.float64)
        scale = jnp.mean(jnp.diag(covariance))
        scale = jnp.where(scale > 0.0, scale, 1.0)
        stabilised = covariance + jnp.eye(covariance.shape[0], dtype=jnp.float64) * (jitter * scale)
        lower = self._factor_jax(self.place(stabilised))
        return jnp.asarray(lower @ self.place(whitened), dtype=jnp.float64)


# ---------------------------------------------------------------------------
# The quasiseparable solve: celerite2.jax, over ampere's own exact term
# ---------------------------------------------------------------------------


@functools.cache
def _celerite2_jax() -> Any:
    """``celerite2.jax``, imported **only after** x64 is already on.

    The deferral is the whole point (see this module's docstring).
    ``celerite2.jax``'s ``__init__`` reads ``jax_enable_x64`` and, if it is
    off, warns and calls ``config.update`` itself — which is exactly what
    ``lowering.md`` §10.2(a) forbids ampere from doing, done to ampere by a
    dependency. Two consequences would follow from importing it at module
    level: a user who imported ``ampere.backends.jax`` would have the
    process-global flag flipped by an ``import`` statement they did not write,
    and :func:`~ampere.backends.jax._config.require_x64` — the guard that is
    supposed to *refuse* when the flag is off — would be satisfied by a flag
    celerite2 had just set on their behalf.

    So the import happens here, on first use, and every caller has already
    passed ``require_x64`` at construction. The check below states that
    invariant rather than trusting it: if it ever fires, the deferral has been
    routed around and the guard has stopped meaning anything.

    :func:`functools.cache` makes it a singleton, so the module object and the
    term class built from it are stable for ``isinstance`` and for celerite2's
    own caches.
    """
    if not x64_enabled():
        raise LikelihoodError(
            "ampere tried to import celerite2.jax while jax's x64 mode was off. It must never "
            "do that: celerite2.jax turns x64 on as an import side effect, which would make "
            "ampere.backends.jax.require_x64() pass on a flag nobody in your program set "
            "(lowering.md §10.2(a)). Call ampere.backends.jax.configure_x64() first. Reaching "
            "this message means the deferral in ampere.backends.jax.gp has been routed around, "
            "which is a bug in ampere rather than in your program."
        )
    import celerite2.jax

    return celerite2.jax


@functools.cache
def _matern32_term_type() -> Any:
    r"""The ``celerite2.jax`` ``Term`` subclass for an **exact** Matérn-3/2.

    ``ampere.core``'s ``_matern32_term_type`` in ``jax.numpy``: the same rank-2
    semiseparable representation, the same midpoint centring, the same
    algebra — transcribed rather than reused, because the core's builds its
    matrices in numpy and coerces the hyperparameters with ``float()``, and a
    traced amplitude or length scale survives neither.

    Never ``celerite2.jax.terms.Matern32Term``: that one is an *approximation*
    (a celerite pair with a small ``eps``, since the celerite basis has no
    :math:`\tau e^{-c\tau}` member), and at its default ``eps = 0.01`` it costs
    about 5e-3 in the log-likelihood on a realistic spectrum — three orders
    outside ``tolerances.cross_solver``. The solver underneath does not need
    the celerite basis at all: it factorises any rank-J semiseparable matrix

    .. math::
        K_{nm} = \sum_j U_{nj} V_{mj} e^{-c_j (t_n - t_m)}\quad (n > m),

    and for :math:`f = \sqrt3/\ell`, :math:`\Delta = t_n - t_m`,

    .. math::
        k(\Delta) = a^2 (1 + f\Delta)e^{-f\Delta}
                  = e^{-f(t_n - t_m)}
                    \big[a^2(1 + f t_n)\cdot 1 + (-a^2 f)\cdot t_m\big],

    with :math:`c = (f, f)` and diagonal :math:`a_n = a^2 + \mathrm{diag}_n`,
    which is algebra rather than a limit.

    The generators grow linearly in the coordinate, so :math:`U_n \cdot V_m` is
    a difference of two large numbers when :math:`f t \gg 1`. The coordinates
    are therefore re-referenced to the midpoint of their own range — the
    products only ever involve differences, so this is exact — which bounds
    the cancellation by half the number of length scales the data span. This
    module's docstring records what that costs, measured.
    """
    terms = _celerite2_jax().terms

    class _ExactMatern32Term(terms.Term):
        """``k(tau) = amplitude**2 (1 + sqrt(3) tau / ell) exp(-sqrt(3) tau / ell)``."""

        def __init__(self, amplitude: Any, length_scale: Any) -> None:
            # Deliberately *not* coerced with float(): these are the values a
            # gradient flows through, and on the traced path they are tracers.
            self.amplitude = jnp.asarray(amplitude, dtype=jnp.float64)
            self.length_scale = jnp.asarray(length_scale, dtype=jnp.float64)

        def get_value(self, tau: Any) -> jax.Array:
            separation = jnp.abs(jnp.atleast_1d(jnp.asarray(tau, dtype=jnp.float64)))
            scaled = _SQRT3 * separation / self.length_scale
            return self.amplitude**2 * (1.0 + scaled) * jnp.exp(-scaled)

        def get_celerite_matrices(
            self,
            x: Any,
            diag: Any,
            **kwargs: Any,
        ) -> tuple[jax.Array, jax.Array, jax.Array, jax.Array]:
            points = jnp.atleast_1d(jnp.asarray(x, dtype=jnp.float64))
            diagonal = jnp.atleast_1d(jnp.asarray(diag, dtype=jnp.float64))
            decay = jnp.asarray(_SQRT3 / self.length_scale, dtype=jnp.float64)
            marginal = self.amplitude * self.amplitude
            # x arrives sorted (QuasisepGP sorts before calling), so the
            # midpoint of the range is (first + last) / 2.
            shifted = points - 0.5 * (points[0] + points[-1])
            return (
                jnp.stack([decay, decay]),
                diagonal + marginal,
                jnp.stack(
                    [
                        marginal * (1.0 + decay * shifted),
                        -marginal * decay * jnp.ones_like(shifted),
                    ],
                    axis=-1,
                ),
                jnp.stack([jnp.ones_like(shifted), shifted], axis=-1),
            )

    return _ExactMatern32Term


def _matern32_term(kernel: Kernel, values: Mapping[str, Any]) -> Any:
    """Build the exact Matérn-3/2 celerite term from a resolved kernel."""
    resolved = kernel.resolve(values)
    return _matern32_term_type()(resolved["amplitude"], resolved["length_scale"])


#: Kernel family -> the builder for its **exact** celerite representation, as
#: ``ampere.core``'s ``_QUASISEPARABLE_TERMS`` is for the numpy path and for
#: the same reason: :class:`QuasisepGP` routes a kernel to the O(N) path only
#: through this table, so a kernel that declares ``QUASISEPARABLE`` without an
#: entry here is refused **by name** rather than silently approximated.
_QUASISEPARABLE_TERMS: dict[str, Any] = {_CoreMatern32.FAMILY: _matern32_term}


@dataclasses.dataclass(frozen=True)
class QuasisepGP(GPSolver):
    """Exact O(N) for ordered 1D data, in jax, over ``celerite2.jax``.

    The jax counterpart of ``ampere.core.QuasisepGP`` — same declaration, same
    exact rank-2 representation of Matérn-3/2, same answer — with the
    recursions run through ``celerite2.jax``'s XLA custom calls so a gradient
    flows to the kernel hyperparameters. This is what makes the flexible
    likelihood a NUTS target on long spectra: 10⁶ points, value and gradient,
    in under a second, where :class:`DenseGP` would need a 10⁶ by 10⁶ matrix.

    A kernel reaches this path only if it declares ``QUASISEPARABLE`` **and**
    this backend holds an exact representation for its family
    (:data:`_QUASISEPARABLE_TERMS`); anything else is refused by name at
    composition time. Coordinates need not arrive sorted — a Gaussian density
    is invariant under a simultaneous permutation of residuals, variances and
    coordinates — so this solver sorts internally and undoes the permutation
    on the way out. The sort is ``jnp.argsort`` rather than numpy's, because
    on the traced path the coordinates may be a traced array; where they are a
    constant, XLA folds it away.

    What is O(N) and what is not, stated plainly, as the core's is:

    * :meth:`log_marginal_likelihood`, :meth:`log_marginal_likelihood_jax`,
      :meth:`latent_transform` and :meth:`latent_transform_jax` are O(N),
      which is what a sampler calls;
    * :meth:`condition` is O(N·M) for M evaluation points, because the
      cross-covariance block is dense by construction;
    * :meth:`conditional_loo` is **refused**, exactly as on the reference path
      (``DEVELOPMENT_PLAN.md`` §2, 2026-09-05): celerite2 exposes no O(N)
      route to the diagonal of ``(K + diag(σ²))⁻¹`` in either interface. See
      that method, and this module's docstring for what tinygp would have
      given here and why it was not taken.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal — the same knob, meaning and default as
        :class:`DenseGP`'s and as the core solver's. Zero by default: a
        covariance that will not factorise is a fact about the model.
    device
        Platform name (``"cpu"``, the default) or an explicit ``jax.Device``.
        Never auto-detected; see :class:`DenseGP`.

    Notes
    -----
    **No float32 opt-out here.** ``celerite2.jax`` coerces every input to
    float64 in its own ``_as_tensor`` and registers its custom calls for
    float64 only, so :class:`DenseGP`'s ``precision=`` argument has no
    counterpart: there is nothing to opt into, and offering the argument and
    ignoring it would be worse than not offering it.

    **BATCHABLE is False**, alone in this backend, because celerite2's
    primitives register no ``vmap`` batching rule — ``jax.vmap`` over a
    density containing one raises ``NotImplementedError: Batching rule for
    'celerite2_factor' not implemented``. The flag is what
    ``ampere.core.declared_capabilities`` reads, so a problem using this
    solver reports ``batchable=False`` and
    :meth:`~ampere.backends.jax.problem.LoweredProblem.log_prob_unconstrained_batched`
    refuses by name rather than failing inside a trace.

    Examples
    --------
    >>> import numpy as np
    >>> from ampere.backends.jax import Matern32, QuasisepGP, configure_x64
    >>> configure_x64()
    >>> kernel = Matern32(0.7, 1.3)
    >>> rng = np.random.default_rng(11)
    >>> t = np.sort(rng.uniform(0.0, 10.0, 200))
    >>> r = rng.normal(0.0, 1.0, 200)
    >>> v = np.full(200, 0.04)
    >>> exact = QuasisepGP().log_marginal_likelihood(kernel, t, r, v, {})
    >>> from ampere.backends.jax import DenseGP
    >>> dense = DenseGP().log_marginal_likelihood(kernel, t, r, v, {})
    >>> bool(abs(exact - dense) < 1e-6)
    True
    """

    jitter: float = 0.0
    device: dataclasses.InitVar[Any] = DEVICE

    NAME: ClassVar[str] = "QuasisepGP"
    EXACT: ClassVar[bool] = True
    REQUIRES_ORDERED_1D: ClassVar[bool] = True
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    #: See the class docstring: celerite2's primitives have no batching rule.
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = DEVICE

    def __post_init__(self, device: Any) -> None:
        require_x64("the jax quasiseparable GP solver")
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"QuasisepGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )
        resolved = _resolve_device(device, self.NAME)
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(
            self, "DEVICE", device if isinstance(device, str) else str(resolved.platform)
        )

    def check_compatible(self, kernel: Kernel, observed: Any) -> None:
        super().check_compatible(kernel, observed)
        if kernel.FAMILY not in _QUASISEPARABLE_TERMS:
            known = ", ".join(sorted(_QUASISEPARABLE_TERMS)) or "(none)"
            raise LikelihoodError(
                f"{type(kernel).__name__} declares QUASISEPARABLE = True, but the jax backend "
                f"holds no exact celerite representation for the {kernel.FAMILY!r} family, so "
                f"{self.NAME} has nothing to lower it to. Families with one: {known}. Add a "
                f"builder to ampere.backends.jax.gp._QUASISEPARABLE_TERMS, or use DenseGP — a "
                f"wrong representation would be an approximation wearing an exact solver's name."
            )

    def provenance_config(self) -> Mapping[str, Any]:
        """``inference.md`` §10a fold-in 10: how this solver computes.

        Records the *library* as well as the dtype and device, because "which
        third-party recursion produced these numbers?" is exactly the question
        a run archived across the jax/numpy split will be asked, and the two
        quasiseparable paths do not use the same one.
        """
        return {
            "library": "celerite2.jax",
            "dtype": "float64",
            "device": self.DEVICE,
            "x64_policy_opt_out": False,
        }

    # -- internals -----------------------------------------------------------

    @staticmethod
    def _axis(coordinates: Any) -> jax.Array:
        """The bare ``(n,)`` coordinate axis, refusing 2D+ as the core does."""
        points = _points(coordinates)
        if points.shape[1] != 1:
            raise LikelihoodError(
                f"QuasisepGP needs one ordered coordinate per sample, but the coordinates have "
                f"{points.shape[1]} per point. check_compatible refuses this at composition "
                f"time; a direct solver call reaches it here. Use DenseGP for 2D+ coordinates."
            )
        return points[:, 0]

    def place(self, array: Any) -> jax.Array:
        """*array* as float64 on this solver's device."""
        return jax.device_put(jnp.asarray(array, dtype=jnp.float64), getattr(self, "_device", None))

    def _sorted(self, coordinates: Any) -> tuple[jax.Array, jax.Array]:
        """The coordinate axis and the permutation that sorts it."""
        axis = self._axis(coordinates)
        return axis, jnp.argsort(axis, stable=True)

    def _factorise(
        self,
        kernel: Kernel,
        ordered_axis: jax.Array,
        ordered_diagonal: jax.Array,
        values: Mapping[str, Any],
    ) -> Any:
        """A ``celerite2.jax`` ``GaussianProcess`` factorised on sorted coordinates.

        No precondition checking here, and that is deliberate: the two callers
        need the failure in two different forms — an exception on the contract
        path, ``-inf`` on the traced one — and a check in the shared helper
        could serve only one of them. celerite2 itself checks nothing: it
        returns quiet NaN where a Cholesky raises (measured), which is the
        reason both callers guard.
        """
        term = _QUASISEPARABLE_TERMS[kernel.FAMILY](kernel, values)
        gp = _celerite2_jax().GaussianProcess(term, mean=0.0)
        gp.compute(ordered_axis, diag=ordered_diagonal, check_sorted=False)
        return gp

    @staticmethod
    def _unsort(ordered: jax.Array, order: jax.Array) -> jax.Array:
        """Undo the sorting permutation on a ``(n,)`` result."""
        return jnp.zeros_like(ordered).at[order].set(ordered)

    # -- the native, traced surface -----------------------------------------

    def log_marginal_likelihood_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        """``log N(residual; 0, K + diag(variance))`` as a traceable jax scalar.

        No exception control flow. Every precondition the contract path raises
        for is tested here as a *value* and folded into one ``jnp.where``: a
        non-finite or negative diagonal, and a factorisation that failed —
        which celerite2 signals as a quiet NaN rather than by raising, so the
        NaN is what this catches. §4.5's ``-inf``, arrived at by the only
        means a trace allows.
        """
        axis, order = self._sorted(coordinates)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        ordered = jnp.take(diagonal, order)
        gp = self._factorise(kernel, jnp.take(axis, order), ordered, values)
        value = jnp.asarray(gp.log_likelihood(jnp.take(r, order)), dtype=jnp.float64)
        usable = jnp.all(jnp.isfinite(ordered)) & jnp.all(ordered >= 0.0) & jnp.isfinite(value)
        return jnp.where(usable, value, -jnp.inf)

    def latent_transform_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> jax.Array:
        """``f = L(θ) z`` in O(N), traceable. See :meth:`DenseGP.latent_transform_jax`."""
        axis, order = self._sorted(coordinates)
        draws = self.place(whitened)
        scale = jnp.mean(jnp.asarray(kernel.diagonal(axis, values), dtype=jnp.float64))
        scale = jnp.where(scale > 0.0, scale, 1.0)
        gp = self._factorise(
            kernel,
            jnp.take(axis, order),
            jnp.full(draws.shape, jitter, dtype=jnp.float64) * scale,
            values,
        )
        transformed = jnp.asarray(gp.dot_tril(jnp.take(draws, order)), dtype=jnp.float64)
        return self._unsort(transformed, order)

    # -- the contract surface -----------------------------------------------

    @staticmethod
    def _check_diagonal(diagonal: jax.Array) -> None:
        """The contract path's preconditions, raising as the core solver's do.

        celerite2's factorisation does not notice a diagonal that cannot
        belong to a covariance — it returns NaN quietly, where a Cholesky
        raises — so the precondition is checked instead, and the messages are
        the reference path's word for word, so both strategies refuse the same
        inputs the same way.
        """
        if not bool(jnp.all(jnp.isfinite(diagonal))):
            raise LikelihoodError(
                "the covariance matrix K + diag(sigma^2) contains non-finite entries, so it "
                "cannot be factorised. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64; constrain the amplitude prior to the data's "
                "own scale."
            )
        if bool(jnp.any(diagonal < 0.0)):
            raise LikelihoodError(
                "the diagonal handed to QuasisepGP contains negative entries, so K + "
                "diag(sigma^2) is not a covariance matrix at all. Variances are squares; this "
                "is a caller error rather than an ill-conditioned problem."
            )

    @staticmethod
    def _refuse_nan(value: Any, *, whitening: bool = False) -> None:
        """Turn celerite2's quiet NaN into the reference path's loud refusal."""
        if bool(jnp.all(jnp.isfinite(jnp.asarray(value)))):
            return
        if whitening:
            raise LikelihoodError(
                "the kernel matrix K is not positive definite in its quasiseparable "
                "representation, so the whitening transform f = L z is undefined. celerite2.jax "
                "reports this as NaN rather than by raising, so ampere checks it here. Increase "
                "the jitter argument, or check the kernel hyperparameters."
            )
        raise LikelihoodError(
            "the covariance matrix K + diag(sigma^2) is not positive definite, so its "
            "quasiseparable factorisation failed. celerite2.jax reports this as NaN rather "
            "than by raising, so ampere checks it here — a NaN log-likelihood in a chain is "
            "not a rejected proposal, it is a number that poisons every diagnostic downstream. "
            "Usual causes: a zero or near-duplicate observational uncertainty, coordinates "
            "closer together than float64 can separate at this length-scale, or a kernel "
            "amplitude far above the data scale. Pass QuasisepGP(jitter=...) — a standard "
            "deviation in the data's units — if the matrix is merely ill-conditioned rather "
            "than wrong."
        )

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> float:
        """``log N(residual; 0, K(values) + diag(variance))``, in O(N)."""
        axis, order = self._sorted(coordinates)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        self._check_diagonal(diagonal)
        gp = self._factorise(kernel, jnp.take(axis, order), jnp.take(diagonal, order), values)
        value = gp.log_likelihood(jnp.take(r, order))
        self._refuse_nan(value)
        return float(value)

    def condition(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
        at: Any = None,
    ) -> GPConditional:
        """The GP posterior at *at* (default: the data coordinates).

        O(N·M): the cross-covariance is dense whatever the solver — M outputs
        each need all N inputs — and only the solve against it is O(N) per
        column.
        """
        axis, order = self._sorted(coordinates)
        points = _points(coordinates)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        self._check_diagonal(diagonal)
        gp = self._factorise(kernel, jnp.take(axis, order), jnp.take(diagonal, order), values)

        solved_residual = jnp.asarray(gp.apply_inverse(jnp.take(r, order)))
        self._refuse_nan(solved_residual)
        alpha = self._unsort(solved_residual, order)
        target = points if at is None else _points(at)
        cross = jnp.asarray(kernel.matrix(target, points, values), dtype=jnp.float64)
        columns = jnp.asarray(gp.apply_inverse(jnp.take(cross.T, order, axis=0)))
        solved = jnp.zeros_like(columns).at[order, :].set(columns)
        prior_variance = jnp.asarray(kernel.diagonal(target, values), dtype=jnp.float64)
        return GPConditional(
            mean=np.asarray(cross @ alpha),
            variance=np.asarray(prior_variance - jnp.einsum("ij,ji->i", cross, solved)),
        )

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """Map whitened latent draws ``z`` to a GP draw ``f = L(θ) z``, in O(N)."""
        axis, order = self._sorted(coordinates)
        draws = self.place(whitened)
        points = _points(coordinates)
        # The same stabilisation DenseGP applies: a jitter relative to the
        # kernel's own scale, so the two solvers factorise the same matrix.
        scale = float(jnp.mean(jnp.asarray(kernel.diagonal(points, values)))) or 1.0
        gp = self._factorise(
            kernel,
            jnp.take(axis, order),
            jnp.full(int(order.size), jitter * scale, dtype=jnp.float64),
            values,
        )
        transformed = jnp.asarray(gp.dot_tril(jnp.take(draws, order)), dtype=jnp.float64)
        self._refuse_nan(transformed, whitening=True)
        return np.asarray(self._unsort(transformed, order))

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Refused, as on the reference path, and for the same missing primitive.

        Every leave-one-out term needs ``A_ii`` for ``A = (K + diag(σ²))⁻¹``,
        and neither of celerite2's interfaces exposes an O(N) route to that
        diagonal: the public surface is ``log_likelihood``, ``apply_inverse``,
        ``dot_tril``, ``predict``, ``condition`` and ``sample``, and the two
        that could give the diagonal — ``predict(..., return_var=True)`` and
        ``condition`` — form the cross-covariance densely and cost O(N·M).
        Supplying one means reimplementing celerite2's internal factorisation
        convention on top of its private ``_d``/``_W``, which is a coupling
        this backend declined for a decomposition nothing stores by default.

        This is the *same* refusal ``ampere.core.QuasisepGP`` makes (W2.3,
        ``DEVELOPMENT_PLAN.md`` §2, 2026-09-05), so the conformance refusal row
        is one row and not two. The alternative was real and was measured:
        tinygp's quasiseparable factor exposes an exact O(N) inverse
        (``L.inv().transpose() @ L.inv()``, whose ``.diag.d`` matched a dense
        inverse to 7e-15), and it was rejected on the factor of 200 in the
        marginal likelihood itself — see this module's docstring.
        """
        raise LikelihoodError(
            f"{self.NAME} does not implement the leave-one-out conditional terms "
            f"(GPSolver.conditional_loo): neither celerite2's numpy nor its jax interface "
            f"exposes an O(N) route to the diagonal of (K + diag(sigma^2))^-1, so this was "
            f"deferred at W2.3 and stays deferred on the jax path for the same reason "
            f"(DEVELOPMENT_PLAN.md §2). DenseGP computes them exactly from its Cholesky — use "
            f"it for pointwise_log_prob under a GP, or the per-dataset log_likelihood group."
        )
