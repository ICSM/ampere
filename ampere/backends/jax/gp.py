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
**celerite2.jax** on the CPU, chosen by measurement in slice 2
(``DEVELOPMENT_PLAN.md`` §2's "jax quasiseparable GP library" row, and §6's
deferred-choice bullet), and over ampere's own ``lax.scan`` recursions on an
accelerator since W7.10 (:mod:`ampere.backends.jax._quasisep`).

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

**W7.10 found where that cost lives, and wrote the scan anyway.** It is not the
recursion's: XLA's CPU backend makes each iteration of a ``scan`` that stacks
per-step outputs cost time proportional to N, and unrolling the loop by two
removes it (``performance_memo.md`` §9). ampere's scan provider
(:mod:`ampere.backends.jax._quasisep`), the same recursions in pure ``jnp``
over the same rank-J arrays, is linear on the CPU — about 2.8 µs per
point for a value and gradient against celerite2's 0.4, so seven times
slower — and it is what :class:`QuasisepGP` runs on an accelerator, where
celerite2 cannot run at all. On the CPU celerite2 stays the default; the scan
is there by name (``provider="scan"``) for batching and for CI.

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
  'celerite2_factor' not implemented``. So the celerite2 provider declares
  ``BATCHABLE = False`` where the rest of the backend declares ``True``, and
  ``ampere.core.declared_capabilities`` aggregates that down to a problem
  which honestly says it cannot be batched. tinygp, being ordinary jax,
  vmaps — and so, since W7.10, does the scan provider, which lifts this cost:
  ``QuasisepGP(provider="scan")`` declares ``BATCHABLE = True`` on itself.
* **its high-level surface exposes no O(N) route to the diagonal of
  ``(K + diag(σ²))⁻¹``**, which is what every leave-one-out term needs. Slice 2
  read that as a loss and left :meth:`QuasisepGP.conditional_loo` refused, as
  the reference path leaves it (``DEVELOPMENT_PLAN.md`` §2, 2026-09-05), while
  noting that tinygp *does* supply one — ``L = solver.factor``,
  ``L.inv().transpose() @ L.inv()`` is a quasiseparable matrix whose
  ``.diag.d`` matched a dense inverse to 7e-15. **Slice 3 closes that loss
  without the swap**: ``celerite2.jax.ops`` is the public entry point to the
  same compiled kernels the torch backend calls through
  ``celerite2.backprop``, so this backend can hold ``c, U, d, W`` in hand and
  run the O(N) backward recursion W2.3 deferred — see
  :meth:`QuasisepGP._precision_diagonal`, which derives it, and W2.4 slice 2's
  decision-log row, where it was derived. The deferral was always about a
  coupling to celerite2's factorisation convention rather than about the
  mathematics; a backend that calls the kernels directly has paid that coupling
  already.

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
import math
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import jax.scipy.linalg as jsl
import numpy as np

from ampere.core import GPConditional, GPSolver, Kernel, Product, Sum, WarpedKernel
from ampere.core import SHO as _CoreSHO
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
from ampere.core.kernels import lookup_quasiseparable_term, refuse_warped_composite

from ._config import BACKEND, require_x64
from ._device import DEVICE, device_flag, place_on, resolve_device
from ._quasisep import CELERITE2, PROVIDERS, QuasisepProvider

__all__ = [
    "DEVICE",
    "SHO",
    "DenseGP",
    "HilbertSpaceGP",
    "JaxOps",
    "Matern12",
    "Matern32",
    "Matern52",
    "Product",
    "QuasisepGP",
    "RotationTerm",
    "SpectralMixture",
    "SquaredExponential",
    "Sum",
    "WarpedKernel",
]

_LOG_2PI = math.log(2.0 * math.pi)
_SQRT3 = math.sqrt(3.0)


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


class JaxOps:
    """:class:`~ampere.core.ArrayOps` in ``jax.numpy``, on one device.

    **W4.5.** This backend used to transcribe each kernel's closed form and
    each celerite representation into ``jax.numpy`` by hand; the torch backend
    did the same, and both carried comments saying the copies had to stay in
    step with the core's, with nothing able to check it. The mathematics now
    lives once in :mod:`ampere.core.kernels`, written against this protocol,
    and what remains backend-specific is this class: float64 coercion (the x64
    policy of ``lowering.md`` §10.2), and placement on the kernel's device.

    One instance per kernel, because the device is per instance since W2.5
    slice 3.

    **Public since W7.13**, and a model's namespace as well as a kernel's: the
    eleven methods after :meth:`log` are what
    :class:`ampere.backends.jax.PortableModel` writes a forward model against.
    """

    def __init__(self, device: Any = None) -> None:
        self.device = device

    def scalar(self, value: Any) -> jax.Array:
        return place_on(jnp.asarray(value, dtype=jnp.float64), self.device)

    def points(self, coordinates: Any, *, dimensions: int | None = None) -> jax.Array:
        array = _points(coordinates)
        if dimensions is not None and int(array.shape[1]) != dimensions:
            raise LikelihoodError(
                f"kernel coordinates have {int(array.shape[1])} coordinate(s) per point, but "
                f"the data they are being compared against have {dimensions}."
            )
        return array

    def n_points(self, coordinates: Any) -> int:
        return int(_points(coordinates).shape[0])

    def separation(self, left: Any, right: Any) -> jax.Array:
        difference = left[:, None, :] - right[None, :, :]
        return jnp.sqrt(jnp.sum(difference * difference, axis=-1))

    def zeros(self, n: int) -> jax.Array:
        return jnp.zeros(n, dtype=jnp.float64)

    def ones_like(self, array: Any) -> jax.Array:
        return jnp.ones_like(jnp.asarray(array, dtype=jnp.float64))

    def stack(self, arrays: Sequence[Any], axis: int = -1) -> jax.Array:
        return jnp.stack([jnp.asarray(item, dtype=jnp.float64) for item in arrays], axis=axis)

    def concatenate(self, arrays: Sequence[Any], axis: int = -1) -> jax.Array:
        return jnp.concatenate([jnp.asarray(item, dtype=jnp.float64) for item in arrays], axis=axis)

    def exp(self, array: Any) -> jax.Array:
        return jnp.exp(jnp.asarray(array, dtype=jnp.float64))

    def cos(self, array: Any) -> jax.Array:
        return jnp.cos(jnp.asarray(array, dtype=jnp.float64))

    def sin(self, array: Any) -> jax.Array:
        return jnp.sin(jnp.asarray(array, dtype=jnp.float64))

    def log1p(self, array: Any) -> jax.Array:
        return jnp.log1p(jnp.asarray(array, dtype=jnp.float64))

    def absolute(self, array: Any) -> jax.Array:
        return jnp.abs(jnp.asarray(array, dtype=jnp.float64))

    def sqrt(self, array: Any) -> jax.Array:
        return jnp.sqrt(jnp.asarray(array, dtype=jnp.float64))

    def log(self, array: Any) -> jax.Array:
        return jnp.log(jnp.asarray(array, dtype=jnp.float64))

    # -- W7.13: the model half of the protocol ----------------------------

    def asarray(self, value: Any) -> jax.Array:
        return place_on(jnp.asarray(value, dtype=jnp.float64), self.device)

    def asindex(self, value: Any) -> jax.Array:
        return place_on(jnp.asarray(value, dtype=jnp.int64), self.device)

    def to_numpy(self, array: Any) -> np.ndarray:
        return np.asarray(array)

    def where(self, condition: Any, if_true: Any, if_false: Any) -> jax.Array:
        return jnp.where(condition, if_true, if_false)

    def interp(self, x: Any, xp: Any, fp: Any) -> jax.Array:
        return jnp.interp(x, xp, fp)

    def cumsum(self, array: Any, axis: int = -1) -> jax.Array:
        return jnp.cumsum(jnp.asarray(array, dtype=jnp.float64), axis=axis)

    def trapezoid(self, y: Any, x: Any, axis: int = -1) -> jax.Array:
        return jnp.trapezoid(y, x, axis=axis)

    def power(self, base: Any, exponent: Any) -> jax.Array:
        return jnp.power(jnp.asarray(base, dtype=jnp.float64), exponent)

    def clip(self, array: Any, low: Any, high: Any) -> jax.Array:
        return jnp.clip(array, low, high)

    def log10(self, array: Any) -> jax.Array:
        return jnp.log10(jnp.asarray(array, dtype=jnp.float64))

    def sum(self, array: Any, axis: int | None = None) -> jax.Array:
        return jnp.sum(array, axis=axis)


class _JaxKernel(Kernel):
    """Shared plumbing for the jax kernels: jax separations, jax covariances.

    ``ampere.core.Kernel``'s arithmetic is written against
    :class:`~ampere.core.ArrayOps`, and this class supplies :class:`JaxOps`;
    ``matrix``, ``diagonal`` and every family's ``_covariance`` are therefore
    the core's, computing in ``jax.numpy``. Before **W4.5** they were all
    overridden here, in transcriptions that had to be kept in step by hand.

    What is still declared here is that the core's ``_positive`` checks do not
    run (``VALIDATES = False``) — right on the reference path, fatal on a
    traced one, where a hyperparameter may be a tracer. The check is not
    dropped silently: a non-positive length scale gives a non-finite
    covariance, the Cholesky then yields NaN, and the guards in
    :class:`DenseGP` turn that into the §4.5 failure the contract asks for. On
    the reference path the check can be an exception because nothing is traced;
    here the same information has to travel as a value. The declaration is what
    keeps it from arising at all: a kernel hyperparameter is declared with a
    prior supported on the positive half-line and a ``Log`` bijection, so no
    sampler can propose a value outside it.
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
    #: Where this kernel's covariance is built. **Per instance since slice 3**
    #: (W2.5): the class default is the CPU and a ``device=`` keyword shadows
    #: it, so a kernel placed on an accelerator beside CPU models is a device
    #: disagreement ``ampere.core.declared_capabilities`` refuses at
    #: composition rather than a mixed-device failure inside a trace. Never
    #: auto-detected (``architecture.md`` §5); see :mod:`ampere.backends.jax._device`.
    DEVICE: ClassVar[str] = DEVICE
    #: The hyperparameters may be tracers; positivity is the declaration's job.
    VALIDATES: ClassVar[bool] = False

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}")
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        object.__setattr__(self, "_ops", JaxOps(resolved))

    def place(self, array: Any) -> jax.Array:
        """*array* as float64 on this kernel's device."""
        return place_on(jnp.asarray(array, dtype=jnp.float64), getattr(self, "_device", None))


class Matern12(_JaxKernel, _CoreMatern12):
    r"""Matérn-1/2 (Ornstein-Uhlenbeck) in ``jax.numpy``. Rank 1 on the O(N) path."""


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


class Matern52(_JaxKernel, _CoreMatern52):
    r"""Matérn-5/2 in ``jax.numpy``. Rank 3 on the O(N) path, and exactly so (W4.5)."""


class SquaredExponential(_JaxKernel, _CoreSquaredExponential):
    r"""Squared exponential (RBF) in ``jax.numpy``: legacy's kernel, kept for comparison.

    .. math::
        k(r) = a^2 \exp\!\left(-\frac{r^2}{2\ell^2}\right)

    Not quasiseparable, and inherits that declaration too, so a quasiseparable
    solver refuses it by name — the asymmetry that made Matérn the default.
    """


class SHO(_JaxKernel, _CoreSHO):
    r"""A damped harmonic oscillator in ``jax.numpy``: the fringing component (W4.5)."""


class RotationTerm(_JaxKernel, _CoreRotationTerm):
    r"""celerite2's rotation pair in ``jax.numpy``: two SHOs, a period and its harmonic."""


class SpectralMixture(_CoreSpectralMixture):
    """A sum of jax :class:`SHO` components (W4.5).

    Not a ``_JaxKernel``: a composite has no arithmetic of its own — it
    combines its children's matrices — so it takes its namespace, its device
    and its capability flags from them, which ``ampere.core.kernels._Composite``
    already does. All this subclass supplies is the default component type.
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


def _columns(array: Any) -> int:
    """How many realisations a right-hand side carries: ``k`` (*W4.2*).

    ``likelihoods.md`` §7 lets a ``residual`` be an ``(n, k)`` block of ``k``
    realisations sharing one covariance, which is how the circular complex GP
    scores its real and imaginary parts against one factorisation. jax's dense
    solver needed almost nothing for it — ``jnp.sum(r * alpha)`` and ``r.size``
    are already the right expressions for a block — so this helper exists for
    the two places that are not: the log-determinant's multiplicity and the
    size of the identity ``conditional_loo`` inverts against.
    """
    shape = jnp.shape(array)
    return 1 if len(shape) < 2 else int(shape[1])


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

    Notes
    -----
    **Batching.** ``BATCHABLE = True`` since slice 2: every operation here is
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
    #: W4.2: this solver takes the circular complex GP's two-column right-hand
    #: side, exactly as ``ampere.core.DenseGP`` does.
    STACKED_RESIDUALS: ClassVar[bool] = True
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
        resolved = resolve_device(device, self.NAME)
        object.__setattr__(self, "_device", resolved)
        # The capability flag is per *instance* here rather than per class, so
        # that a solver placed on an accelerator says so and
        # `declared_capabilities` can catch it beside CPU models.
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))

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
        return place_on(jnp.asarray(array, dtype=self.dtype), getattr(self, "_device", None))

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
            -0.5 * (quadratic + _columns(r) * log_determinant + r.size * _LOG_2PI),
            dtype=jnp.float64,
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
        return -0.5 * (quadratic + _columns(r) * log_determinant + r.size * _LOG_2PI)

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
        rows = int(jnp.shape(r)[0])
        columns = _columns(r)
        precision_diagonal = jnp.diag(self._solve(lower, jnp.eye(rows, dtype=self.dtype)))
        # W4.2: one term per *sample*, summing the k components, which share
        # A_ii -- the quantity a leave-one-out conditional is built from.
        quadratic = jnp.sum(jnp.reshape(alpha, (rows, columns)) ** 2, axis=1)
        return np.asarray(
            columns * (0.5 * jnp.log(precision_diagonal) - 0.5 * _LOG_2PI)
            - quadratic / (2.0 * precision_diagonal)
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
        # An (n, k) residual conditions to an (m, k) mean beside one (m,)
        # variance: a posterior variance does not depend on the data (W4.2).
        mean = cross @ self._solve(lower, jnp.asarray(residual, dtype=jnp.float64))
        solved = self._solve(lower, cross.T)
        prior_variance = jnp.asarray(kernel.diagonal(target, values), dtype=jnp.float64)
        posterior = prior_variance - jnp.einsum("ij,ji->i", cross, solved)
        return GPConditional(mean=np.asarray(mean), variance=np.asarray(posterior))

    def latent_transform_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> jax.Array:
        """``f = L(θ) z`` as a **traceable, differentiable** jax array (W2.14).

        The counterpart of :meth:`log_marginal_likelihood_jax`, and needed for
        the same reason: :meth:`latent_transform` returns numpy, so a realised
        density that called it would cut the gradient in exactly the
        hyperparameters the latent path exists to fit. Since W2.14 the
        contract path applies this transform on every latent evaluation
        (``GaussianProcessNoise.noise_params``), so the realised path has to
        apply it too, natively, or disagree with its own oracle.

        No exception control flow (``inference.md`` §10a): a kernel matrix
        that will not factorise leaves ``jnp.linalg.cholesky`` as NaN, which
        propagates into the family's term and becomes ``-inf`` at the one
        ``jnp.where`` the lowered dataset ends with.
        """
        points = _points(coordinates)
        covariance = jnp.asarray(kernel.matrix(points, points, values), dtype=jnp.float64)
        # The reference solver's ``float(...) or 1.0``, written as a value so
        # a traced amplitude survives it: a zero-variance kernel would give a
        # zero stabiliser and no factorisation at all.
        mean_variance = jnp.mean(jnp.diag(covariance))
        scale = jnp.where(mean_variance != 0.0, mean_variance, 1.0)
        stabilised = covariance + jnp.eye(covariance.shape[0], dtype=jnp.float64) * (jitter * scale)
        lower = self._factor_jax(stabilised)
        return lower @ jnp.asarray(whitened, dtype=jnp.float64)

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
        solver can impose it in whatever representation it uses. The contract
        surface, so it raises where :meth:`latent_transform_jax` returns NaN.
        """
        points = _points(coordinates)
        covariance = jnp.asarray(kernel.matrix(points, points, values), dtype=jnp.float64)
        if not bool(jnp.all(jnp.isfinite(covariance))):
            raise LikelihoodError(
                "the kernel matrix K contains non-finite entries, so the whitening transform "
                "f = L z is undefined. The usual cause is a kernel amplitude large enough that "
                "amplitude**2 overflows float64."
            )
        realised = self.latent_transform_jax(kernel, points, whitened, values, jitter=jitter)
        if not bool(jnp.all(jnp.isfinite(realised))):
            raise LikelihoodError(
                "the kernel matrix K is not positive definite, so the whitening transform "
                "f = L z is undefined (jax reports this as NaN rather than by raising). "
                "Increase the jitter argument, or check the kernel hyperparameters."
            )
        return np.asarray(realised)


# ---------------------------------------------------------------------------
# The quasiseparable solve: ampere's own exact term, over a recursion provider
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _Factorisation:
    """``K + diag(a) = L D Lᵀ`` on sorted coordinates, as the provider returned it.

    ``t`` is the coordinate the recursion ran over — w(x) for a warped kernel
    (W5.7), the sorted axis otherwise — and ``c``, ``U`` the generators the
    solves need beside ``d`` and ``W``. Everything :class:`QuasisepGP` computes
    is a composition of these five arrays and the provider's four operations.
    """

    t: jax.Array
    c: jax.Array
    U: jax.Array
    d: jax.Array
    W: jax.Array


def _celerite_matrices(
    kernel: Kernel,
    values: Mapping[str, Any],
    ordered_axis: Any,
    ordered_diagonal: Any,
) -> tuple[jax.Array, jax.Array, jax.Array, jax.Array]:
    r"""``(c, a, U, V)``: the exact celerite representation of a kernel, in jax.

    Straight from the public registry — :func:`ampere.core.kernels.
    lookup_quasiseparable_term`, which all three backends read — and built by
    the *same* function the reference path uses, in :class:`JaxOps` instead of
    numpy (**W4.5**: a family this backend could lower and another could not is
    no longer expressible). ``a`` is the full diagonal, ``ordered_diagonal``
    plus the kernel's own marginal variance.

    Never ``celerite2.jax.terms.Matern32Term``, and the reason is unchanged
    since slice 2: that one is an *approximation* (a celerite pair with a small
    ``eps``, since the celerite basis has no :math:`\tau e^{-c\tau}` member),
    and at its default ``eps = 0.01`` it costs about 5e-3 in the log-likelihood
    on a realistic spectrum — three orders outside ``tolerances.cross_solver``.
    The recursions do not need the celerite basis at all: they factorise any
    rank-J semiseparable matrix

    .. math::
        K_{nm} = \sum_j U_{nj} V_{mj} e^{-c_j (t_n - t_m)}\quad (n > m),

    which is what :class:`~ampere.core.CeleriteRepresentation` is, and in which
    Matérn-1/2, -3/2 and -5/2, the SHO and the rotation pair are all exact.

    **W7.10** removed the ``celerite2.jax`` ``Term`` subclass this used to be
    wrapped in: the arrays are what either provider consumes, and building them
    without celerite2 is what lets the scan provider run without importing it.
    The arithmetic is the wrapper's own, operation for operation, so the
    celerite2 provider's numbers did not move (checked bit for bit).

    The builders receive the **raw** sorted coordinate and apply any warp
    themselves (W5.7); they re-reference it to the midpoint of its own range,
    which is exact since only differences enter, and which bounds the
    cancellation this module's docstring measures.
    """
    resolved = kernel.resolve(values)
    builder = lookup_quasiseparable_term(kernel.FAMILY)
    points = jnp.atleast_1d(jnp.asarray(ordered_axis, dtype=jnp.float64))
    diagonal = jnp.atleast_1d(jnp.asarray(ordered_diagonal, dtype=jnp.float64))
    representation = builder(kernel.with_ops(JaxOps(None)), resolved, points)
    return (
        representation.decay,
        diagonal + representation.marginal,
        representation.left,
        representation.right,
    )


def _bare_coordinates(coordinates: Any, kernel: Kernel | None = None) -> Any:
    """An ``(n,)`` numpy view of concrete coordinates, refusing 2D+ as the core does.

    The numpy twin of :meth:`QuasisepGP._axis`, used only where the
    coordinates are already concrete. ``kernel`` selects its own column out of
    a multi-axis container's coordinates (W4.5).
    """
    array = np.asarray(coordinates, dtype=float)
    if kernel is not None and array.ndim == 2 and array.shape[1] > 1:
        array = np.asarray(kernel.select(array), dtype=float)
    if array.ndim == 1:
        return array
    if array.ndim == 2 and array.shape[1] == 1:
        return array[:, 0]
    raise LikelihoodError(
        f"QuasisepGP needs one ordered coordinate per sample, but the coordinates have shape "
        f"{array.shape}. check_compatible refuses this at composition time; a direct solver "
        f"call reaches it here. Use DenseGP for 2D+ coordinates."
    )


@dataclasses.dataclass(frozen=True)
class QuasisepGP(GPSolver):
    """Exact O(N) for ordered 1D data, in jax, over a choice of two recursions.

    The jax counterpart of ``ampere.core.QuasisepGP`` — same declaration, same
    exact rank-2 representation of Matérn-3/2, same answer — with the
    recursions run differentiably, so a gradient flows to the kernel
    hyperparameters. This is what makes the flexible likelihood a NUTS target
    on long spectra: 10⁶ points, value and gradient, in under a second on the
    CPU, where :class:`DenseGP` would need a 10⁶ by 10⁶ matrix.

    **Two providers** (W7.10, :mod:`ampere.backends.jax._quasisep`) compute
    the same factorisation: ``"celerite2"``, ``celerite2.jax``'s XLA custom
    calls into compiled C++, CPU only and the faster there; and ``"scan"``,
    ampere's own ``lax.scan`` recursions in pure ``jnp``, which run on any
    device and batch. ``provider="auto"`` (the default) chooses by the
    *asked-for* device — celerite2 on the CPU, the scan anywhere else — which
    is the rule ``performance_memo.md`` §9's table gives; the choice is
    recorded in :meth:`provenance_config`.

    A kernel reaches this path only if it declares ``QUASISEPARABLE`` **and**
    every family in its tree has a representation in the public registry
    (:func:`~ampere.core.kernels.register_quasiseparable_term`, which all three
    backends read); anything else is refused by name at
    composition time. Coordinates need not arrive sorted — a Gaussian density
    is invariant under a simultaneous permutation of residuals, variances and
    coordinates — so this solver sorts internally and undoes the permutation
    on the way out. The sort is ``jnp.argsort`` rather than numpy's, because
    on the traced path the coordinates may be a traced array; where they are a
    constant, XLA folds it away.

    What is O(N) and what is not, stated plainly, as the core's is:

    * :meth:`log_marginal_likelihood`, :meth:`log_marginal_likelihood_jax`
      and :meth:`latent_transform` are O(N), which is what a sampler calls;
    * :meth:`condition` is O(N·M) for M evaluation points, because the
      cross-covariance block is dense by construction;
    * :meth:`conditional_loo` is O(N J³) since slice 3 — the recursion W2.3
      deferred and W2.4 slice 2 derived, over the provider's own ``factor``,
      so the terms agree with :meth:`DenseGP.conditional_loo` exactly rather
      than approximately. ``ampere.core.QuasisepGP`` keeps its
      refusal: the numpy path's circumstances are unchanged.

    Parameters
    ----------
    jitter
        A standard deviation, in the data's own units, added in quadrature to
        the diagonal — the same knob, meaning and default as
        :class:`DenseGP`'s and as the core solver's. Zero by default: a
        covariance that will not factorise is a fact about the model.
    device
        Platform name (``"cpu"``, the default) or an explicit ``jax.Device``.
        Never auto-detected; see :class:`DenseGP`. Off the CPU the solver
        runs the scan provider: celerite2's primitives have CPU lowerings
        alone, so ``provider="celerite2"`` on an accelerator is refused by
        name (W6.18's refusal, now naming the scan as the way out).
    provider
        ``"auto"`` (the default), ``"celerite2"`` or ``"scan"``; anything else
        is refused naming the three. Like ``device``, an ``InitVar``:
        configuration rather than declaration, so it stays out of the spec
        hash, and a problem declared once fits on either.

    Notes
    -----
    **No float32 opt-out here.** ``celerite2.jax`` coerces every input to
    float64 in its own ``_as_tensor`` and registers its custom calls for
    float64 only, and the scan provider is float64 by the same policy, so
    :class:`DenseGP`'s ``precision=`` argument has no counterpart: offering
    the argument and ignoring it would be worse than not offering it.

    **BATCHABLE is per provider.** The class flag is ``False``, the celerite2
    provider's, because its primitives register no ``vmap`` batching rule —
    ``jax.vmap`` over a density containing one raises
    ``NotImplementedError: Batching rule for 'celerite2_factor' not
    implemented``. A scan instance sets ``BATCHABLE = True`` on itself, as
    ``DEVICE`` is set per instance. The flag is what
    ``ampere.core.declared_capabilities`` reads off the instance, so a
    problem on the celerite2 provider reports ``batchable=False`` and
    :meth:`~ampere.backends.jax.problem.LoweredProblem.log_prob_unconstrained_batched`
    refuses by name rather than failing inside a trace, while one on the scan
    provider batches.

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
    provider: dataclasses.InitVar[str] = "auto"

    NAME: ClassVar[str] = "QuasisepGP"
    EXACT: ClassVar[bool] = True
    REQUIRES_ORDERED_1D: ClassVar[bool] = True
    REQUIRES_QUASISEPARABLE: ClassVar[bool] = True
    IMPLEMENTED: ClassVar[bool] = True
    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    #: The celerite2 provider's flag, which is the CPU default: its primitives
    #: have no batching rule. A scan instance sets ``True`` on itself (W7.10).
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = DEVICE
    #: The recursion provider (:mod:`ampere.backends.jax._quasisep`); an
    #: instance attribute set at construction shadows this default.
    _provider: ClassVar[QuasisepProvider] = CELERITE2

    def __post_init__(self, device: Any, provider: str) -> None:
        require_x64("the jax quasiseparable GP solver")
        if not math.isfinite(self.jitter) or self.jitter < 0.0:
            raise LikelihoodError(
                f"QuasisepGP's jitter must be finite and >= 0, got {self.jitter!r}."
            )
        if provider != "auto" and provider not in PROVIDERS:
            raise LikelihoodError(
                f"{self.NAME} has no recursion provider {provider!r}. The three choices are "
                f"'auto' (celerite2 on the CPU, the scan on an accelerator), 'celerite2' "
                f"(celerite2's compiled kernels, CPU only) and 'scan' (ampere's lax.scan "
                f"recursions, any device, batchable)."
            )
        resolved = resolve_device(device, self.NAME)
        platform = str(resolved.platform)
        if provider == "auto":
            # The device rule (W7.10, ``performance_memo.md`` §9): celerite2 on
            # the CPU, where it is the faster of the two at every size measured;
            # the scan wherever else, since celerite2 cannot run there at all.
            # Chosen from the *asked-for* device, never detected (§5).
            provider = CELERITE2.NAME if platform == "cpu" else "scan"
        chosen = PROVIDERS[provider]
        if chosen.PLATFORMS is not None and platform not in chosen.PLATFORMS:
            # celerite2.jax registers its primitives' MLIR lowerings for the
            # CPU alone (celerite2 0.3.3, ``celerite2/jax/ops.py``), so on an
            # accelerator the first factorisation dies inside the trace with
            # "MLIR translation rule for primitive 'celerite2_factor' not found
            # for platform cuda" — the first GPU run's finding (W6.18,
            # 2026-10-05). Refused here, by name, where the user can read it;
            # ampere never substitutes another provider silently.
            raise LikelihoodError(
                f"{self.NAME}(provider={provider!r}) runs on the CPU only: its O(N) "
                f"recursions are celerite2's compiled kernels, whose jax primitives have CPU "
                f"lowerings alone, and {platform!r} has none. Use provider='scan' (or the "
                f"default, provider='auto', which chooses it off the CPU) for ampere's own "
                f"lax.scan recursions, which run on any device."
            )
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        object.__setattr__(self, "_provider", chosen)
        object.__setattr__(self, "BATCHABLE", chosen.BATCHABLE)

    def check_compatible(self, kernel: Kernel, observed: Any) -> None:
        super().check_compatible(kernel, observed)
        # W5.7: a warp under a composite has no single recursion coordinate.
        refuse_warped_composite(kernel, self.NAME)
        # Every family in the tree, not just the root: a Sum lowers term by
        # term, so one unregistered term stops it and must be named here.
        for leaf in kernel.leaves():
            lookup_quasiseparable_term(leaf.FAMILY, owner=self.NAME)

    def provenance_config(self) -> Mapping[str, Any]:
        """``inference.md`` §10a fold-in 10: how this solver computes.

        Records the *library* as well as the dtype and device, because "which
        recursion produced these numbers?" is exactly the question a run
        archived across the jax/numpy split will be asked. Since W7.10 the
        answer can differ within this class too: ``library`` is
        ``"celerite2.jax"`` or ``"ampere.jax.scan"``, and ``provider`` the
        ``provider=`` spelling that selects it — what ``"auto"`` resolved to,
        never ``"auto"`` itself.
        """
        return {
            "library": self._provider.LIBRARY,
            "provider": self._provider.NAME,
            "dtype": "float64",
            "device": self.DEVICE,
            "x64_policy_opt_out": False,
        }

    # -- internals -----------------------------------------------------------

    @staticmethod
    def _axis(coordinates: Any, kernel: Kernel | None = None) -> jax.Array:
        """The bare ``(n,)`` coordinate axis, refusing 2D+ as the core does.

        ``kernel`` is passed since W4.5 so an axis-selecting kernel gets its own
        column out of a multi-axis container's coordinates.
        """
        points = _points(coordinates)
        if kernel is not None:
            points = kernel.select(points)
        if points.shape[1] != 1:
            raise LikelihoodError(
                f"QuasisepGP needs one ordered coordinate per sample, but the coordinates have "
                f"{points.shape[1]} per point. check_compatible refuses this at composition "
                f"time; a direct solver call reaches it here. Use DenseGP for 2D+ coordinates."
            )
        return points[:, 0]

    def place(self, array: Any) -> jax.Array:
        """*array* as float64 on this solver's device."""
        return place_on(jnp.asarray(array, dtype=jnp.float64), getattr(self, "_device", None))

    def _sorted(
        self, coordinates: Any, kernel: Kernel | None = None
    ) -> tuple[jax.Array, jax.Array]:
        """The coordinate axis and the permutation that sorts it.

        Sorted in **numpy** when the coordinates are concrete, which in
        ampere they always are: a dataset's coordinates come from its observed
        container and cannot depend on θ (``Dataset`` resolves the effective
        mask once and forbids a parameter-dependent one). Leaving the sort to
        ``jnp.argsort`` was correct but not free — XLA constant-folds it at
        *compile* time, which at 3e4 points takes over a second and prints
        its own slow-operation alarm at the user. numpy does the same sort in
        milliseconds, once, at lowering time.

        The traced branch is kept rather than assumed away: a caller may hand
        this solver a traced coordinate array (a fitted grid offset, say), and
        refusing one merely because ampere does not do it today would be a
        limitation invented here. ``stable=True`` in both branches, so the two
        agree on ties.
        """
        if isinstance(coordinates, jax.core.Tracer):
            axis = self._axis(coordinates, kernel)
            return axis, jnp.argsort(axis, stable=True)
        # The *argument* is tested, not the converted axis: jax stages
        # operations on a captured constant into the jaxpr, so by the time
        # `_axis` has run `jnp.asarray` the value looks traced whether or not
        # it ever was. Asking the caller's own object is the only test that
        # distinguishes "this is data" from "this is being traced".
        raw = np.asarray(_bare_coordinates(coordinates, kernel), dtype=float)
        # Placed on this solver's device (W7.10): on a node with an accelerator
        # an unplaced array lands on jax's *default* device, which is not the
        # one a CPU solver was asked for.
        order = place_on(
            jnp.asarray(np.argsort(raw, kind="stable")), getattr(self, "_device", None)
        )
        return self.place(raw), order

    def _factorise(
        self,
        kernel: Kernel,
        ordered_axis: jax.Array,
        ordered_diagonal: jax.Array,
        values: Mapping[str, Any],
    ) -> _Factorisation:
        """``K + diag = L D Lᵀ`` on sorted coordinates, by this solver's provider.

        No precondition checking here, and that is deliberate: the callers
        need the failure in two different forms — an exception on the contract
        path, ``-inf`` on the traced one — and a check in the shared helper
        could serve only one of them. Neither provider raises: a failed pivot
        comes back as a ``d`` that is not positive and finite (see
        :mod:`ampere.backends.jax._quasisep`), which is the reason every caller
        guards.

        **W7.10.** Until W7.10 this built a ``celerite2.jax.GaussianProcess``
        (W5.17's subclass of it, whose ``lax.cond`` compiled once). The raw
        arrays and the provider's operations replace it: there is no
        ``lax.cond`` on this route at all, so W5.17's compile-once property
        holds by construction rather than by keeping branch functions stable.
        """
        c, a, U, V = _celerite_matrices(kernel, values, ordered_axis, ordered_diagonal)
        # **W5.7.** The recursion's propagators come from the coordinate it runs
        # over, and a warped kernel is stationary in w(x) rather than in x.
        # ``warped_coordinate`` is the identity for every other kernel; a warp is
        # monotone by construction, so the sorted axis stays sorted. The
        # generators above came from the raw axis, as the registry documents.
        t = jnp.asarray(kernel.warped_coordinate(ordered_axis, values), dtype=jnp.float64)
        d, W = self._provider.factor(t, c, a, U, V)
        return _Factorisation(t=t, c=c, U=U, d=d, W=W)

    def _log_likelihood(self, factorisation: _Factorisation, residual: jax.Array) -> jax.Array:
        """``log N(residual; 0, L D Lᵀ)`` from a factorisation, sorted residual in.

        celerite2's ``GaussianProcess.log_likelihood`` written out over the
        raw arrays, operation for operation (so the celerite2 provider's
        numbers are bit for bit those of the wrapper W7.10 removed): the
        log-determinant ``Σ log d``, a forward substitution, and the
        quadratic form ``Σ z² / d``. A failed factorisation folds to ``-inf``
        through one ``jnp.where``, where celerite2 used a ``lax.cond``.
        """
        f = factorisation
        log_det = jnp.sum(jnp.log(f.d))
        bad = jnp.any(f.d <= 0) | (~jnp.isfinite(log_det))
        size = f.d.shape[0]
        norm = jnp.where(bad, -jnp.inf, -0.5 * (log_det + size * jnp.log(2 * jnp.pi)))
        alpha = self._provider.solve_lower(f.t, f.c, f.U, f.W, residual[:, None])[:, 0]
        return norm - 0.5 * jnp.sum(alpha**2 / f.d)

    def _solve(self, factorisation: _Factorisation, right: jax.Array) -> jax.Array:
        """``(K + diag)⁻¹ right`` for an ``(n,)`` or ``(n, m)`` right-hand side, in O(N).

        A forward substitution, a division by ``D``, a back substitution —
        celerite2's own ``_do_solve``.
        """
        f = factorisation
        columns = right[:, None] if right.ndim == 1 else right
        z = self._provider.solve_lower(f.t, f.c, f.U, f.W, columns)
        z = self._provider.solve_upper(f.t, f.c, f.U, f.W, z / f.d[:, None])
        return z[:, 0] if right.ndim == 1 else z

    def _dot_tril(self, factorisation: _Factorisation, draws: jax.Array) -> jax.Array:
        """``L D^{1/2} z`` for an ``(n,)`` *draws*: celerite2's ``dot_tril``."""
        f = factorisation
        z = draws[:, None] * jnp.sqrt(f.d)[:, None]
        z = z + self._provider.matmul_lower(f.t, f.c, f.U, f.W, z)
        return z[:, 0]

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
        axis, order = self._sorted(coordinates, kernel)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        ordered = jnp.take(diagonal, order)
        factorisation = self._factorise(kernel, jnp.take(axis, order), ordered, values)
        value = jnp.asarray(self._log_likelihood(factorisation, jnp.take(r, order)))
        usable = jnp.all(jnp.isfinite(ordered)) & jnp.all(ordered >= 0.0) & jnp.isfinite(value)
        return jnp.where(usable, value, -jnp.inf)

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
        axis, order = self._sorted(coordinates, kernel)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        self._check_diagonal(diagonal)
        factorisation = self._factorise(
            kernel, jnp.take(axis, order), jnp.take(diagonal, order), values
        )
        value = self._log_likelihood(factorisation, jnp.take(r, order))
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
        axis, order = self._sorted(coordinates, kernel)
        points = _points(coordinates)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        self._check_diagonal(diagonal)
        factorisation = self._factorise(
            kernel, jnp.take(axis, order), jnp.take(diagonal, order), values
        )

        solved_residual = self._solve(factorisation, jnp.take(r, order))
        self._refuse_nan(solved_residual)
        alpha = self._unsort(solved_residual, order)
        target = points if at is None else _points(at)
        cross = jnp.asarray(kernel.matrix(target, points, values), dtype=jnp.float64)
        columns = self._solve(factorisation, jnp.take(cross.T, order, axis=0))
        solved = jnp.zeros_like(columns).at[order, :].set(columns)
        prior_variance = jnp.asarray(kernel.diagonal(target, values), dtype=jnp.float64)
        return GPConditional(
            mean=np.asarray(cross @ alpha),
            variance=np.asarray(prior_variance - jnp.einsum("ij,ji->i", cross, solved)),
        )

    def latent_transform_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> jax.Array:
        """``f = L(θ) z`` in O(N), traceable and differentiable (W2.14).

        The quasiseparable counterpart of :meth:`DenseGP.latent_transform_jax`,
        and the one that matters at the scale the latent path is sized for:
        forming ``L`` densely is O(N³), so a 10⁵-sample latent fit reaches its
        correlation through this recursion or not at all. celerite2's
        ``dot_tril`` carries its own differentiation rules, so the amplitude
        and the length scale keep their gradients through it.
        """
        axis, order = self._sorted(coordinates, kernel)
        draws = self.place(whitened)
        points = _points(coordinates)
        # The same stabilisation DenseGP applies: a jitter relative to the
        # kernel's own scale, so the two solvers factorise the same matrix.
        # Kept as a *value* rather than a Python float, so a traced amplitude
        # survives it (the reference path's ``float(...) or 1.0``).
        mean_variance = jnp.mean(jnp.asarray(kernel.diagonal(points, values), dtype=jnp.float64))
        scale = jnp.where(mean_variance != 0.0, mean_variance, 1.0)
        factorisation = self._factorise(
            kernel,
            jnp.take(axis, order),
            jnp.full(int(order.size), jitter, dtype=jnp.float64) * scale,
            values,
        )
        transformed = self._dot_tril(factorisation, jnp.take(draws, order))
        return self._unsort(transformed, order)

    def latent_transform(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> np.ndarray:
        """Map whitened latent draws ``z`` to a GP draw ``f = L(θ) z``, in O(N).

        The contract surface: celerite2's quiet NaN becomes the reference
        path's loud refusal, which :meth:`latent_transform_jax` cannot do.
        """
        realised = self.latent_transform_jax(kernel, coordinates, whitened, values, jitter=jitter)
        self._refuse_nan(realised, whitening=True)
        return np.asarray(realised)

    # -- the leave-one-out decomposition, in O(N) (W2.5 slice 3) --------------

    @staticmethod
    def _precision_diagonal(
        axis: jax.Array,
        c: jax.Array,
        U: jax.Array,
        d: jax.Array,
        W: jax.Array,
    ) -> jax.Array:
        r"""``diag((K + diag(σ²))⁻¹)`` in O(N J³) — the recursion W2.3 deferred.

        **This is what lifts the refusal on this backend** (W2.5 slice 3,
        mirroring the torch path's W2.4 slice 2; ``DEVELOPMENT_PLAN.md`` §2).
        W2.3's reason was precise and was about a *coupling*, not about the
        mathematics: the leave-one-out terms need ``A_ii`` for
        ``A = (K + diag(σ²))⁻¹``; celerite2's public **numpy** interface
        exposes no O(N) route to it; and the route that exists "reimplements
        celerite2's internal factorisation convention". This backend pays that
        coupling anyway — :meth:`_factorise` calls the provider's ``factor``
        directly — and the coupling is *tested*, against
        :meth:`DenseGP.conditional_loo`'s Cholesky, both here and in the
        conformance row.

        The derivation, so a reader need not rebuild it. With
        ``K + diag = L D Lᵀ`` and ``M = L⁻¹`` (unit lower triangular),
        ``A = L⁻ᵀ D⁻¹ L⁻¹`` gives

        .. math:: A_{ii} = \sum_{k \ge i} M_{ki}^2 / d_k .

        Column *i* of ``M`` solves ``L z = e_i``, and celerite's forward
        substitution — ``f_k = p_k ⊙ (f_{k-1} + W_{k-1} z_{k-1})``,
        ``z_k = y_k - U_k · f_k``, ``p_k = e^{-c(t_k - t_{k-1})}`` — becomes,
        once ``y = e_i`` is substituted in,

        .. math::
            f^{(i)}_{i+1} = p_{i+1} ⊙ W_i,\qquad
            f^{(i)}_k = G_k f^{(i)}_{k-1},\quad
            G_k = \mathrm{diag}(p_k)\,(I - W_{k-1} U_{k-1}^\top),

        and ``G_k`` **does not depend on i**. That is the whole trick: every
        column of ``M`` is the same linear recursion started from a different
        vector, so the sum over columns collapses into one backward
        accumulation of a ``J x J`` matrix,

        .. math::
            A_{ii} = 1/d_i + w_i^\top R_i w_i,\qquad w_i = p_{i+1} ⊙ W_i,
            \qquad R_i = \frac{U_{i+1} U_{i+1}^\top}{d_{i+1}}
                       + G_{i+2}^\top R_{i+1} G_{i+2},

        with ``R_{N-1} = 0``. No inverses appear, so it is as stable as the
        factorisation itself.

        Written as one :func:`jax.lax.scan` rather than the torch path's
        Python loop, which is the one thing that had to change in the
        transcription: an N-step Python loop over ``jnp`` operations stages N
        copies of the body into the jaxpr, and at 10⁵ points that is a trace
        long enough to be the dominant cost. ``scan`` stages the body once.
        The ``reverse=True`` scan runs ``i = N-2 … 0`` and stores each term at
        its own index, so no ``.at[].set()`` bookkeeping is needed either.

        It is a post-processing surface (``results.md`` §15 R2 stores this
        decomposition only on request), never part of a density — but it is
        pure ``jnp`` and traceable all the same, because writing it any other
        way would have meant a second implementation to keep in step.
        """
        size, rank = U.shape
        if size == 1:
            return 1.0 / d
        gaps = axis[1:] - axis[:-1]
        # decays[j] = p_{j+1}; transitions[j] = G_{j+1}; starts[i] = w_i.
        decays = jnp.exp(-c[None, :] * gaps[:, None])
        identity = jnp.eye(rank, dtype=U.dtype)
        transitions = decays[:, :, None] * (
            identity[None, :, :] - W[:-1, :, None] * U[:-1, None, :]
        )
        starts = decays * W[:-1, :]
        # steps[i] = G_{i+2}. The final entry is the identity: it multiplies
        # R = 0 on the first (i = N-2) step, so it stands in for the G that
        # would be out of range without a branch inside the scan.
        steps = jnp.concatenate([transitions, identity[None, :, :]], axis=0)[1:]

        def step(carry: jax.Array, row: Any) -> tuple[jax.Array, jax.Array]:
            transition, u_next, d_next, start, d_here = row
            accumulated = transition.T @ carry @ transition + jnp.outer(u_next, u_next) / d_next
            return accumulated, 1.0 / d_here + start @ accumulated @ start

        _, terms = jax.lax.scan(
            step,
            jnp.zeros((rank, rank), dtype=U.dtype),
            (steps, U[1:], d[1:], starts, d[:-1]),
            reverse=True,
        )
        return jnp.concatenate([terms, 1.0 / d[-1:]])

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """The leave-one-out conditional terms, in O(N) (W2.5 slice 3).

        The same closed form :meth:`DenseGP.conditional_loo` uses
        (Sundararajan & Keerthi 2001) over the same matrix, so the two agree
        exactly rather than approximately: with
        ``A = (K + diag(variance + jitter²))⁻¹``,
        ``sigma_i^{2,-i} = 1 / A_ii``, ``mu_i^{-i} = y_i - [A r]_i / A_ii`` and

            ``log p_i = ½ log A_ii - [A r]_i² / (2 A_ii) - ½ log 2π``.

        What this solver supplies that the reference one still cannot is
        ``A_ii`` in linear time — see :meth:`_precision_diagonal`, which
        records the recursion and why this backend is entitled to it where
        ``ampere.core.QuasisepGP`` is not.

        A contract surface, so it raises where celerite2 goes quietly NaN.
        """
        axis, order = self._sorted(coordinates, kernel)
        r = self.place(residual)
        diagonal = self.place(variance) + self.jitter**2
        self._check_diagonal(diagonal)
        f = self._factorise(kernel, jnp.take(axis, order), jnp.take(diagonal, order), values)
        # Either provider signals a failed factorisation as a non-positive or
        # non-finite D, quietly, exactly as it does for the log-likelihood.
        if not bool(jnp.all(jnp.isfinite(f.d)) & jnp.all(f.d > 0.0)):
            self._refuse_nan(jnp.asarray(jnp.nan))
        # W5.7: the recursion below walks the factorisation, so it runs on the
        # same coordinate the factorisation did — f.t is w(x), the identity
        # unwarped. It consumes (c, U, d, W) from either provider unchanged.
        alpha = self._solve(f, jnp.take(r, order))
        precision = self._precision_diagonal(f.t, f.c, f.U, f.d, f.W)
        terms = 0.5 * jnp.log(precision) - alpha**2 / (2.0 * precision) - 0.5 * _LOG_2PI
        return np.asarray(self._unsort(terms, order))


# ---------------------------------------------------------------------------
# The reduced-rank spectral solve: a Hilbert-space basis in jax (W5.4)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class HilbertSpaceGP(GPSolver):
    r"""Approximate ``O(N m + m³)`` by a Hilbert-space basis, in ``jax.numpy``.

    The jax twin of :class:`ampere.core.HilbertSpaceGP` — same declaration,
    same basis, same Woodbury algebra — with the traced surfaces
    :mod:`ampere.backends.jax.problem` composes:
    :meth:`log_marginal_likelihood_jax` and :meth:`latent_transform_jax`.

    Two properties make it the reduced-rank latent path's natural home here.
    The whitening ``f = Phi z`` takes ``m`` variables rather than ``N``, so a
    NUTS chain over a long spectrum has a manageable dimension; and every step
    is ``jnp`` on whole arrays with an ``m`` by ``m`` Cholesky at the end, so
    :func:`jax.vmap` maps it as it maps the dense solve — :attr:`BATCHABLE` is
    ``True`` for the same measured reason :class:`DenseGP`'s is.

    The box is computed in **numpy**, from concrete coordinates: it is data,
    fixed before the fit, and never differentiated. Coordinates reach a
    lowered problem as closed-over constants, so this holds under ``jit``;
    ``vmap``\ ping over the *coordinates* themselves is the one thing it would
    not survive, and nothing in ampere does that.

    See ``ampere.core.HilbertSpaceGP`` for the mathematics and for the warning
    about choosing ``basis_size`` and ``boundary_factor`` together.
    """

    basis_size: int | tuple[int, ...] = DEFAULT_BASIS_SIZE
    boundary_factor: float = DEFAULT_BOUNDARY_FACTOR
    jitter: float = 0.0
    precision: dataclasses.InitVar[str] = "float64"
    device: dataclasses.InitVar[Any] = DEVICE

    NAME: ClassVar[str] = "HilbertSpaceGP"
    EXACT: ClassVar[bool] = False
    IMPLEMENTED: ClassVar[bool] = True
    STACKED_RESIDUALS: ClassVar[bool] = True
    BACKEND: ClassVar[str] = BACKEND
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = DEVICE

    def __post_init__(self, precision: str, device: Any) -> None:
        require_x64("the jax reduced-rank GP solver")
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
        object.__setattr__(self, "_dtype", _solve_dtype(precision, self.NAME))
        object.__setattr__(self, "_precision", str(precision))
        resolved = resolve_device(device, self.NAME)
        object.__setattr__(self, "_device", resolved)
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))

    # -- declarations -------------------------------------------------------

    @property
    def counts(self) -> tuple[int, ...]:
        """Basis members per axis, normalised to a tuple."""
        return normalise_counts(self.basis_size)

    @property
    def precision_name(self) -> str:
        """``"float64"`` or ``"float32"``: which precision this solve runs in."""
        return str(getattr(self, "_precision", "float64"))

    @property
    def dtype(self) -> Any:
        """The jax dtype the capacitance matrix is factorised in."""
        return getattr(self, "_dtype", jnp.float64)

    def place(self, array: Any) -> jax.Array:
        """*array* on this solver's device, as its solve dtype."""
        return place_on(jnp.asarray(array, dtype=self.dtype), getattr(self, "_device", None))

    def latent_size(self, kernel: Kernel, n_samples: int) -> int:
        """``m``, not ``N``. See ``ampere.core.HilbertSpaceGP.latent_size``."""
        del kernel, n_samples
        return basis_size(self.counts)

    def provenance_config(self) -> Mapping[str, Any]:
        """The approximation, and how this backend computed it.

        ``basis_size`` and ``boundary_factor`` are dataclass fields as well,
        deliberately: they are declaration, not configuration, and two
        backends computing one problem must agree on them or they are not
        computing the same quantity (``results.md`` §14). Precision and device
        are configuration and stay out of the spec, as they do for
        :class:`DenseGP`.
        """
        return {
            "library": "jax",
            "basis_size": [int(count) for count in self.counts],
            "boundary_factor": float(self.boundary_factor),
            "dtype": self.precision_name,
            "device": self.DEVICE,
            "x64_policy_opt_out": self.precision_name != "float64",
        }

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

    # -- internals ----------------------------------------------------------

    def _basis(self, kernel: Kernel, coordinates: Any) -> HilbertSpaceBasis:
        """The box and the frequency grid, in numpy, from the *raw* coordinates.

        ``np.asarray`` on what arrives rather than on ``_points(...)``, and
        that is load-bearing under ``jit``: a ``jnp`` operation applied to a
        closed-over constant inside a trace returns a **tracer**, which numpy
        then refuses to convert, whereas the concrete array the trace closed
        over converts happily. Coordinates are data — a lowered problem holds
        them as a device array built once, before any tracing — so reading
        them here costs nothing and detaches nothing.
        """
        block = np.asarray(coordinates, dtype=np.float64)
        if block.ndim == 1:
            block = block[:, None]
        return hilbert_basis(kernel.select(block), self.counts, self.boundary_factor)

    def _scaled_basis(
        self,
        kernel: Kernel,
        basis: HilbertSpaceBasis,
        coordinates: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        r"""The basis scaled by the square root of the spectral density.

        ``Phi`` is built through :class:`JaxOps`, so it is a jax array on this
        solver's device; ``S`` comes from the *kernel*'s namespace, which is
        where the gradient in the amplitude and the length scale lives. The
        clamp is the traced counterpart of the reference solver's check: a
        realised density must not raise, so an inadmissible hyperparameter
        travels to ``-inf`` rather than to an exception.
        """
        ops = JaxOps(getattr(self, "_device", None))
        matrix = basis_matrix(basis, kernel.select(_points(coordinates)), ops)
        density = jnp.asarray(spectral_values(kernel, basis, values), dtype=jnp.float64)
        return self.place(matrix * jnp.sqrt(jnp.clip(density, 0.0, None))[None, :])

    def _diagonal(self, variance: Any) -> jax.Array:
        return self.place(variance) + self.jitter**2

    def _factor_jax(self, scaled: jax.Array, diagonal: jax.Array) -> tuple[jax.Array, jax.Array]:
        r"""``(D^-1 Phi, chol(I + Phi^T D^-1 Phi))``. **NaN** where it is not definite."""
        weighted = scaled / diagonal[:, None]
        identity = jnp.eye(int(scaled.shape[1]), dtype=self.dtype)
        return weighted, jnp.linalg.cholesky(identity + scaled.T @ weighted)

    @staticmethod
    def _solve(
        scaled: jax.Array,
        diagonal: jax.Array,
        weighted: jax.Array,
        lower: jax.Array,
        right: jax.Array,
    ) -> jax.Array:
        r"""``(K + D)^-1 right`` by Woodbury, for an ``(n,)`` or ``(n, k)`` block."""
        divisor = diagonal if jnp.ndim(right) == 1 else diagonal[:, None]
        direct = right / divisor
        return direct - weighted @ jsl.cho_solve((lower, True), scaled.T @ direct)

    def _refuse(self, lower: jax.Array, diagonal: jax.Array) -> None:
        """The contract path's guard: a :class:`LikelihoodError` where jax gives NaN."""
        if not bool(jnp.all(jnp.isfinite(diagonal)) & jnp.all(diagonal > 0.0)):
            raise LikelihoodError(
                f"{self.NAME} needs a strictly positive noise diagonal: the Woodbury identity it "
                f"solves by inverts diag(sigma^2 + jitter^2) directly, so a zero uncertainty is "
                f"a division by zero rather than an ill-conditioned matrix. Pass "
                f"{self.NAME}(jitter=...), or use DenseGP."
            )
        if not bool(jnp.all(jnp.isfinite(lower))):
            raise LikelihoodError(
                "the reduced-rank capacitance matrix I + Phi^T D^-1 Phi is not positive "
                "definite, so the Woodbury solve failed (jax reports this as NaN rather than by "
                "raising). It is positive definite for every admissible hyperparameter, so this "
                "is a numerical rather than a structural failure: reduce basis_size, or raise "
                "the jitter."
            )

    # -- the native, traced surface -----------------------------------------

    def log_marginal_likelihood_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> jax.Array:
        """The marginal likelihood of the approximation, as a traceable jax scalar.

        No exception control flow: a factorisation that fails leaves NaN in
        ``lower`` and the final :func:`jax.numpy.where` turns it into ``-inf``,
        which is §4.5's failure signal computed the only way a traced function
        can compute it.
        """
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        diagonal = self._diagonal(variance)
        weighted, lower = self._factor_jax(scaled, diagonal)
        r = self.place(residual)
        alpha = self._solve(scaled, diagonal, weighted, lower, r)
        log_determinant = jnp.sum(jnp.log(diagonal)) + 2.0 * jnp.sum(
            jnp.log(jnp.abs(jnp.diag(lower)))
        )
        quadratic = jnp.sum(r * alpha)
        value = jnp.asarray(
            -0.5 * (quadratic + _columns(r) * log_determinant + r.size * _LOG_2PI),
            dtype=jnp.float64,
        )
        admissible = jnp.all(diagonal > 0.0) & jnp.isfinite(value)
        return jnp.where(admissible, value, -jnp.inf)

    def latent_transform_jax(
        self,
        kernel: Kernel,
        coordinates: Any,
        whitened: Any,
        values: Mapping[str, Any],
        *,
        jitter: float = 1e-10,
    ) -> jax.Array:
        r"""``f = Phi z`` as a traceable, differentiable jax array.

        ``z`` carries :meth:`latent_size` entries. The reduced-rank factor *is*
        the whitening, so there is no factorisation to stabilise and the
        ``jitter`` argument the two exact solvers need is accepted and unused.
        """
        del jitter
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        return scaled @ jnp.asarray(whitened, dtype=jnp.float64)

    # -- the contract surface -----------------------------------------------

    def log_marginal_likelihood(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> float:
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        diagonal = self._diagonal(variance)
        weighted, lower = self._factor_jax(scaled, diagonal)
        self._refuse(lower, diagonal)
        r = self.place(residual)
        alpha = self._solve(scaled, diagonal, weighted, lower, r)
        log_determinant = float(jnp.sum(jnp.log(diagonal))) + 2.0 * float(
            jnp.sum(jnp.log(jnp.abs(jnp.diag(lower))))
        )
        quadratic = float(jnp.sum(r * alpha))
        return -0.5 * (quadratic + _columns(r) * log_determinant + r.size * _LOG_2PI)

    def conditional_loo(
        self,
        kernel: Kernel,
        coordinates: Any,
        residual: Any,
        variance: Any,
        values: Mapping[str, Any],
    ) -> np.ndarray:
        """Sundararajan & Keerthi from the Woodbury diagonal. See the core solver."""
        basis = self._basis(kernel, coordinates)
        scaled = self._scaled_basis(kernel, basis, coordinates, values)
        diagonal = self._diagonal(variance)
        weighted, lower = self._factor_jax(scaled, diagonal)
        self._refuse(lower, diagonal)
        r = self.place(residual)
        alpha = self._solve(scaled, diagonal, weighted, lower, r)
        triangular = jsl.solve_triangular(lower, weighted.T, lower=True)
        precision_diagonal = 1.0 / diagonal - jnp.sum(triangular * triangular, axis=0)
        rows = int(jnp.shape(r)[0])
        columns = _columns(r)
        quadratic = jnp.sum(jnp.reshape(alpha, (rows, columns)) ** 2, axis=1)
        return np.asarray(
            columns * (0.5 * jnp.log(precision_diagonal) - 0.5 * _LOG_2PI)
            - quadratic / (2.0 * precision_diagonal)
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
        """The posterior of the rank-``m`` process. See the core solver."""
        points = _points(coordinates)
        basis = self._basis(kernel, points)
        scaled = self._scaled_basis(kernel, basis, points, values)
        diagonal = self._diagonal(variance)
        weighted, lower = self._factor_jax(scaled, diagonal)
        self._refuse(lower, diagonal)
        r = self.place(residual)
        alpha = self._solve(scaled, diagonal, weighted, lower, r)
        target = points if at is None else _points(at)
        target_scaled = self._scaled_basis(kernel, basis, target, values)
        mean = target_scaled @ (scaled.T @ alpha)
        triangular = jsl.solve_triangular(lower, target_scaled.T, lower=True)
        posterior = jnp.sum(triangular * triangular, axis=0)
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
        """``f = Phi z`` — the contract surface, so it raises where the traced one does not."""
        points = _points(coordinates)
        expected = self.latent_size(kernel, int(jnp.shape(points)[0]))
        if int(np.shape(whitened)[0]) != expected:
            raise LikelihoodError(
                f"{self.NAME} whitens {expected} basis coefficient(s), but it was handed "
                f"{int(np.shape(whitened)[0])} whitened value(s). This solver's latent block is "
                f"the basis size m, not the sample count — see GPSolver.latent_size."
            )
        realised = self.latent_transform_jax(kernel, points, whitened, values, jitter=jitter)
        if not bool(jnp.all(jnp.isfinite(realised))):
            raise LikelihoodError(
                "the reduced-rank whitening f = Phi diag(sqrt(S)) z is not finite. The usual "
                "cause is a kernel amplitude large enough that amplitude**2 overflows float64."
            )
        return np.asarray(realised)
