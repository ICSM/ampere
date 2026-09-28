"""Point estimates and warm starts: ``optimise`` and ``warm_start_gp`` (W6.7).

``inference.md`` §10b. Every engine starts from prior draws that merely score
finite (:meth:`~ampere.inference.engine.Engine.initial_positions`), and a user
with an expensive or badly conditioned problem had no cheap way to ask "where
should I be looking?". This module answers with a **function**, not an
:class:`~ampere.inference.engine.Engine`: an optimiser produces no posterior,
so it has no ``run()`` that emits a ``DataTree`` — it returns an
:class:`~ampere.results.Optimum`, which a sampler then takes as its start
(``Engine.initial_positions(count, around=optimum)``, ``run(initial=optimum)``).

The objective
-------------
**The constrained-space posterior density**, ``log p(θ) + log p(D | θ)`` at
``θ = constrain(u)``, optimised over the packed unconstrained vector ``u`` so
that every iterate stays inside the support. It is *not*
``problem.log_prob_unconstrained(u)``, which adds the change-of-variables term
``Σ log|dθ/du|``: that density is what NUTS samples, and its maximum moves with
the choice of bijection (a ``Log`` on a scale parameter pulls it towards larger
values by exactly ``log θ``), so a "MAP" of it is a property of the
parametrisation rather than of the posterior. Every route here maximises the
same function, and a conformance row holds the ``"scipy"`` and ``"map"``
modes to ``1e-3`` in every free parameter on each fixture both run on.

On the native path the realisation only offers ``log_prob_unconstrained``, so
the Jacobian term is **subtracted on the realised side**, and its own gradient
supplied by central differences on the numpy contract path
(:func:`~ampere.inference.engine.unconstrained_jacobian_correction`). The
alternative — rebuilding each bijection natively inside the backend's graph —
was rejected: ``ampere.inference`` imports no backend, and the bijection a
backend lowers is not always a registry lookup (torch's inferred default goes
through ``biject_to(support)``, ``lowering.md`` §4), so a native rebuild here
would restate each backend's lowering table and could drift from it. The
Jacobian term touches no model and no data, so its finite differences cost
microseconds beside one model evaluation, and they are exact to ~1e-9 because
the term is a smooth, elementwise function of ``u``.

The routes
----------
``"scipy"``
    Any backend, gradient-free: multi-start :func:`scipy.optimize.minimize`
    (Powell by default; ``minimiser="L-BFGS-B"`` selectable) on the numpy
    contract path — the harvested ``ScipyMinOpt`` shape
    (``docs/design/harvest/optim_only/``) re-expressed over
    :class:`~ampere.core.dataset.FittingProblem`. The covariance is the
    inverse of a central-difference Hessian.
``"map"``
    torch and jax: a gradient MAP through :func:`ampere.core.realise` —
    ``torch.optim.LBFGS`` (Adam as the fallback) or
    ``jax.scipy.optimize.minimize(method="BFGS")`` (``optax.adam`` as the
    fallback) — with an autodiff Hessian.
``"vi"``
    torch and jax: :class:`~ampere.inference.VIEngine`'s fitted ``laplace``
    guide, its mean as the point and its covariance as the covariance.
``"auto"``
    ``"map"`` where the problem is native and differentiable, else
    ``"scipy"``.

Multi-start draws come from the same prior-draw helper the engines use
(:func:`~ampere.inference.engine.draw_prior_positions`), and the best start by
objective is kept, with every start summarised in the ``Optimum``.

Not here, by ruling (2026-09-24): Bayesian optimisation. Its linear surrogate
on a thin shell assumes a high-dimensional standard-normal latent and a locally
linear objective, which a posterior of ten to fifty curved, often multimodal
parameters does not supply; an expensive simulator is the SBI layer's case.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.results import Optimum, StartSummary, hash_of, provenance_attrs

from .engine import draw_prior_positions
from .exceptions import EngineError

__all__ = [
    "OPTIMISE_METHODS",
    "constrained_objective",
    "covariance_from_hessian",
    "finite_difference_hessian",
    "optimise",
]

#: The routes :func:`optimise` accepts, ``"auto"`` included.
OPTIMISE_METHODS = ("auto", "scipy", "map", "vi")

#: The scipy route's default local minimiser: derivative-free, robust to a
#: ``-inf`` wall at a support boundary the bijections did not remove (a model
#: that refuses a region), and the harvested ``ScipyMinOpt``'s own choice.
DEFAULT_MINIMISER = "Powell"

#: The options each route understands; anything else is refused by name.
_SCIPY_OPTIONS = frozenset({"minimiser", "tol", "minimiser_options"})

#: Relative floor on the Hessian's eigenvalues: below it the curvature is
#: indistinguishable from zero at finite-difference precision, and a
#: covariance built from it would be a confident ball around a flat direction.
_EIGENVALUE_FLOOR = 1e-8


# ---------------------------------------------------------------------------
# The objective
# ---------------------------------------------------------------------------


def constrained_objective(problem: FittingProblem) -> Callable[[np.ndarray], float]:
    """``u -> log p(constrain(u)) + log p(D | constrain(u))`` on the numpy path.

    The constrained-space posterior density as a function of the packed
    unconstrained vector — what every route here maximises (module
    docstring). ``-inf`` outside the support or where the model cannot score.
    """

    def objective(u: np.ndarray) -> float:
        theta = problem.constrain(np.asarray(u, dtype=float))
        value = float(problem.log_prob(theta))
        return value if math.isfinite(value) else -math.inf

    return objective


def finite_difference_hessian(
    function: Callable[[np.ndarray], float], point: np.ndarray
) -> np.ndarray:
    """The Hessian of a scalar *function* at *point*, by central differences.

    Two passes: the diagonal at a fixed step first, to learn each coordinate's
    curvature scale ``1/sqrt(|H_ii|)``, then the full matrix at a step of 5 %
    of that scale (clipped to ``[1e-6, 1e-1]``) — so a posterior a thousand
    times narrower than unit width in ``u`` is not differenced across its
    whole width, and a broad one is not differenced into rounding noise. No
    new dependency, which is the ruling (``inference.md`` §10b).
    """
    x = np.asarray(point, dtype=float).reshape(-1)
    n = x.size
    f0 = function(x)

    def diagonal(step: np.ndarray) -> np.ndarray:
        out = np.empty(n)
        for i in range(n):
            e = np.zeros(n)
            e[i] = step[i]
            out[i] = (function(x + e) - 2.0 * f0 + function(x - e)) / step[i] ** 2
        return out

    first = diagonal(np.full(n, 1e-3))
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(
            np.isfinite(first) & (np.abs(first) > 0), 1.0 / np.sqrt(np.abs(first)), 1.0
        )
    step = np.clip(0.05 * scale, 1e-6, 1e-1)
    hessian = np.empty((n, n))
    hessian[np.diag_indices(n)] = diagonal(step)
    for i in range(n):
        for j in range(i + 1, n):
            ei = np.zeros(n)
            ej = np.zeros(n)
            ei[i] = step[i]
            ej[j] = step[j]
            value = (
                function(x + ei + ej)
                - function(x + ei - ej)
                - function(x - ei + ej)
                + function(x - ei - ej)
            ) / (4.0 * step[i] * step[j])
            hessian[i, j] = hessian[j, i] = value
    return hessian


def covariance_from_hessian(hessian: Any) -> tuple[np.ndarray | None, str | None]:
    """``(inverse, None)`` for a positive-definite Hessian, else ``(None, reason)``.

    *hessian* is that of the **negated** objective (so a mode has a positive
    definite one). A non-finite entry, a non-positive eigenvalue, or one below
    :data:`_EIGENVALUE_FLOOR` times the largest is refused by name — the
    memo's "refusal by name" — rather than regularised, because a regularised
    inverse around a saddle or a ridge would hand a sampler a confident ball
    in a direction the posterior does not constrain.
    """
    matrix = np.asarray(hessian, dtype=float)
    if not np.all(np.isfinite(matrix)):
        return None, (
            "the Hessian of the negative log-posterior at the mode has non-finite entries (the "
            "mode sits against a region the model cannot score, or the density overflowed)"
        )
    matrix = 0.5 * (matrix + matrix.T)
    eigenvalues = np.linalg.eigvalsh(matrix)
    largest = float(np.max(np.abs(eigenvalues))) if eigenvalues.size else 0.0
    smallest = float(np.min(eigenvalues)) if eigenvalues.size else 0.0
    if largest == 0.0 or smallest <= _EIGENVALUE_FLOOR * largest:
        kind = "negative" if smallest < 0 else "zero (to finite precision)"
        return None, (
            f"the Hessian of the negative log-posterior at the mode is not positive definite: "
            f"its smallest eigenvalue is {smallest:.3g} ({kind}) against a largest of "
            f"{largest:.3g}, so the point is a saddle or lies on a ridge the posterior does not "
            f"constrain, and no covariance is reported rather than a regularised one"
        )
    return np.linalg.inv(matrix), None


# ---------------------------------------------------------------------------
# The entry point
# ---------------------------------------------------------------------------


def optimise(
    problem: FittingProblem,
    *,
    method: str = "auto",
    starts: int = 8,
    seed: int | None = None,
    **options: Any,
) -> Optimum:
    """A point estimate of *problem*'s constrained-space posterior mode.

    Parameters
    ----------
    problem
        The composed problem.
    method
        ``"scipy"``, ``"map"``, ``"vi"`` or ``"auto"`` (module docstring).
    starts
        Independent starts drawn from the joint prior; the best by objective
        is kept and every one summarised in :attr:`Optimum.starts`. The
        ``"vi"`` route fits one guide from the first start.
    seed
        Seed for the start draws. ``None`` takes them from the problem's own
        seed, on the ``"optimise.initialisation"`` stream, so a seeded problem
        optimises reproducibly without one.
    **options
        Route-specific settings, refused by name when the route does not know
        them. ``"scipy"``: ``minimiser`` (default ``"Powell"``; any
        :func:`scipy.optimize.minimize` method, ``"L-BFGS-B"`` the usual
        alternative), ``tol``, ``minimiser_options`` (passed as ``options=``).
        ``"map"``: ``steps`` (the iteration budget, default 500),
        ``learning_rate`` (Adam's, default 0.05). ``"vi"``: ``steps``
        (default 2000), ``learning_rate`` (default 0.01), ``draws`` (default
        200).

    Returns
    -------
    Optimum
        The mode in both coordinate systems, the density at it in both
        conventions, the covariance or a named refusal, and provenance.
    """
    if method not in OPTIMISE_METHODS:
        raise EngineError(
            f"optimise does not know the method {method!r}. Available: "
            f"{', '.join(OPTIMISE_METHODS)}. 'scipy' runs on every backend without gradients; "
            f"'map' and 'vi' need a native differentiable (torch or jax) problem."
        )
    if int(starts) < 1:
        raise EngineError(f"optimise needs at least one start, got {starts}.")
    if problem.free_size == 0:
        raise EngineError(
            f"optimise has nothing to optimise: this problem declares no free parameters "
            f"(every one of {list(problem.parameters.names)} is fixed)."
        )
    route = _resolve(problem, method)
    rng = (
        problem.rng("optimise.initialisation") if seed is None else np.random.default_rng(int(seed))
    )
    positions = draw_prior_positions(
        problem, int(starts), rng, problem.log_prob, who=f"optimise({route!r})"
    )
    if route == "scipy":
        return _scipy_route(problem, positions, dict(options))
    if route == "map":
        return _map_route(problem, positions, dict(options))
    return _vi_route(problem, positions, dict(options))


def _resolve(problem: FittingProblem, method: str) -> str:
    """Which route *method* means for *problem*, refusing a native one by name."""
    from ampere.core.realisation import registered_realisations

    native = problem.differentiable and problem.backend in registered_realisations()
    if method == "auto":
        return "map" if native else "scipy"
    if method in ("map", "vi") and not native:
        why = (
            f"its backend is {problem.backend!r}, which has no differentiable realisation"
            if problem.backend not in registered_realisations()
            else "it is not differentiable (a piece declares DIFFERENTIABLE = False)"
        )
        raise EngineError(
            f"optimise(method={method!r}) needs a gradient through ampere.core.realise, and this "
            f"problem cannot supply one: {why}. Use method='scipy', which optimises the same "
            f"constrained-space density without gradients on every backend, or build the problem "
            f"from ampere.backends.torch or ampere.backends.jax pieces."
        )
    return method


def _check_options(route: str, options: Mapping[str, Any], known: frozenset[str]) -> None:
    unknown = sorted(set(options) - known)
    if unknown:
        raise EngineError(
            f"optimise(method={route!r}) does not take the option(s) {unknown}. It takes: "
            f"{', '.join(sorted(known)) or '(none)'}."
        )


def start_hash(u: np.ndarray) -> str:
    """:func:`hash_of` of a start vector, as :class:`StartSummary` records it."""
    return hash_of([float(v) for v in np.asarray(u, dtype=float).reshape(-1)])


def build_optimum(
    problem: FittingProblem,
    *,
    route: str,
    unconstrained: np.ndarray,
    covariance: np.ndarray | None,
    refusal: str | None,
    converged: bool,
    message: str,
    evaluations: int,
    starts: tuple[StartSummary, ...],
) -> Optimum:
    """Assemble an :class:`Optimum` at *unconstrained*, scoring both conventions.

    Both densities are recomputed here on the numpy contract path rather than
    taken from the optimiser, so every route reports them the same way and a
    native route's float32/float64 rounding cannot make two routes' numbers
    disagree about the same point.
    """
    u = np.asarray(unconstrained, dtype=float).reshape(-1)
    theta = problem.constrain(u)
    values = problem.parameters.unpack(theta)
    free = problem.parameters.free_names
    constrained = {
        name: (
            float(np.asarray(values[name]))
            if np.ndim(values[name]) == 0
            else np.asarray(values[name], dtype=float)
        )
        for name in free
    }
    return Optimum(
        route=route,
        backend=problem.backend,
        free_names=tuple(free),
        free_labels=tuple(problem.parameters.free_labels()),
        unconstrained=u,
        constrained=constrained,
        log_prob_constrained=float(problem.log_prob(theta)),
        log_prob_unconstrained=float(problem.log_prob_unconstrained(u)),
        covariance=covariance,
        covariance_refusal=refusal,
        converged=bool(converged),
        message=message,
        evaluations=int(evaluations),
        starts=starts,
        provenance=provenance_attrs(problem, engine=f"optimise.{route}"),
    )


# ---------------------------------------------------------------------------
# Route 1: scipy on the numpy path
# ---------------------------------------------------------------------------


def _scipy_route(
    problem: FittingProblem, positions: np.ndarray, options: dict[str, Any]
) -> Optimum:
    from scipy.optimize import minimize

    _check_options("scipy", options, _SCIPY_OPTIONS)
    minimiser = str(options.get("minimiser", DEFAULT_MINIMISER))
    tol = options.get("tol")
    minimiser_options = dict(options.get("minimiser_options") or {})
    objective = constrained_objective(problem)
    calls = 0

    def negative(u: np.ndarray) -> float:
        nonlocal calls
        calls += 1
        value = objective(u)
        # scipy's minimisers read +inf as "worse than anything", which is what
        # a point the model cannot score is; NaN would poison a line search.
        return -value if math.isfinite(value) else math.inf

    summaries: list[StartSummary] = []
    best: Any = None
    for theta in positions:
        u0 = problem.unconstrain(theta)
        result = minimize(negative, u0, method=minimiser, tol=tol, options=minimiser_options)
        value = -float(result.fun)
        status = "converged" if bool(result.success) else f"not converged: {result.message}"
        if not math.isfinite(value):
            status = f"failed: the optimiser ended where the density is {value}"
        summaries.append(StartSummary(start_hash(u0), value, status))
        if math.isfinite(value) and (best is None or value > -float(best.fun)):
            best = result
    if best is None:
        raise EngineError(
            f"optimise('scipy'): none of the {len(positions)} starts ended at a point the "
            f"problem can score. problem.failure_summary() says why:\n{problem.failure_summary()}"
        )
    mode = np.asarray(best.x, dtype=float)
    before = calls
    hessian = finite_difference_hessian(lambda u: -objective(u), mode)
    covariance, refusal = covariance_from_hessian(hessian)
    return build_optimum(
        problem,
        route="scipy",
        unconstrained=mode,
        covariance=covariance,
        refusal=refusal,
        converged=bool(best.success),
        message=f"scipy.optimize.minimize({minimiser!r}): {best.message}",
        evaluations=before + _hessian_calls(mode.size),
        starts=tuple(summaries),
    )


def _map_route(problem: FittingProblem, positions: np.ndarray, options: dict[str, Any]) -> Optimum:
    raise EngineError("optimise(method='map') lands in W6.7 unit 5.")


def _vi_route(problem: FittingProblem, positions: np.ndarray, options: dict[str, Any]) -> Optimum:
    raise EngineError("optimise(method='vi') lands in W6.7 unit 6.")


def _hessian_calls(n: int) -> int:
    """Objective evaluations :func:`finite_difference_hessian` makes in *n* dimensions."""
    return 1 + 2 * n + 2 * n + 4 * (n * (n - 1) // 2)
