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

from .engine import draw_prior_positions, unconstrained_jacobian_correction
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

#: The scipy route's default ``tol``. scipy's own Powell default (``xtol`` and
#: ``ftol`` of ``1e-4``, the latter *relative* to a log-density in the
#: hundreds) stops a few hundredths of a nat short, which is a visible
#: fraction of a posterior width in ``u``; this is what holds the scipy and
#: native modes to the ``1e-3`` the conformance row asks, for a few hundred
#: more evaluations per start.
DEFAULT_TOLERANCE = 1e-8

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
    tol = float(options.get("tol", DEFAULT_TOLERANCE))
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


# ---------------------------------------------------------------------------
# The Jacobian term, for the native routes
# ---------------------------------------------------------------------------


def jacobian_term(problem: FittingProblem, u: np.ndarray) -> float:
    """``Σ log|dθ/du|`` at *u* — what the realised density carries and a MAP must not.

    :func:`~ampere.inference.engine.unconstrained_jacobian_correction`, the one
    place that knows the formula, for one vector. ``0`` where the prior itself
    is not finite, so the term cannot turn the realised ``-inf`` into ``nan``.
    """
    value = float(unconstrained_jacobian_correction(problem, u)[0])
    return value if math.isfinite(value) else 0.0


def jacobian_gradient(problem: FittingProblem, u: np.ndarray) -> np.ndarray:
    """The gradient of :func:`jacobian_term` at *u*, by central differences.

    The term is a smooth elementwise function of ``u`` that touches no model
    and no data, so ``2 n`` evaluations cost microseconds and a step of
    ``1e-5 max(1, |u_i|)`` leaves an error near ``1e-10``.
    """
    x = np.asarray(u, dtype=float).reshape(-1)
    out = np.empty(x.size)
    for i in range(x.size):
        h = 1e-5 * max(1.0, abs(x[i]))
        up = x.copy()
        down = x.copy()
        up[i] += h
        down[i] -= h
        out[i] = (jacobian_term(problem, up) - jacobian_term(problem, down)) / (2.0 * h)
    return out


def _gradient_converged(gradient: np.ndarray, value: float) -> bool:
    """A native optimiser's convergence test: a finite point with a vanishing gradient."""
    return bool(
        math.isfinite(value)
        and np.all(np.isfinite(gradient))
        and float(np.max(np.abs(gradient))) < _GRADIENT_TOLERANCE
    )


#: A native optimiser has converged when every component of the objective's
#: gradient in ``u`` is below this — about a thousandth of a posterior
#: standard deviation's worth of log-density slope for a well-scaled problem.
_GRADIENT_TOLERANCE = 1e-3

#: The native routes' options and their defaults.
_MAP_DEFAULTS: dict[str, Any] = {"steps": 500, "learning_rate": 0.05}


# ---------------------------------------------------------------------------
# Route 2: the gradient MAP through realise
# ---------------------------------------------------------------------------


def _map_route(problem: FittingProblem, positions: np.ndarray, options: dict[str, Any]) -> Optimum:
    from ampere.core.realisation import realise

    _check_options("map", options, frozenset(_MAP_DEFAULTS))
    settings = {**_MAP_DEFAULTS, **options}
    steps = int(settings["steps"])
    learning_rate = float(settings["learning_rate"])
    if steps < 1 or learning_rate <= 0:
        raise EngineError(
            f"optimise('map') needs steps >= 1 and a positive learning_rate, got {steps} and "
            f"{learning_rate}."
        )
    realised = realise(problem)
    if problem.backend == "torch":
        fit, hessian = _torch_fit, _torch_hessian
    elif problem.backend == "jax":
        fit, hessian = _jax_fit, _jax_hessian
    else:  # a third registered realisation: the scipy route is the honest answer
        raise EngineError(
            f"optimise('map') knows how to drive the torch and jax realisations; this problem's "
            f"backend is {problem.backend!r}. Use method='scipy'."
        )
    summaries: list[StartSummary] = []
    best: tuple[float, np.ndarray, bool, str] | None = None
    evaluations = 0
    for theta in positions:
        u0 = problem.unconstrain(theta)
        u, value, converged, message, calls = fit(problem, realised, u0, steps, learning_rate)
        evaluations += calls
        status = "converged" if converged else f"not converged: {message}"
        if not math.isfinite(value):
            status = f"failed: the optimiser ended where the density is {value}"
        summaries.append(StartSummary(start_hash(u0), value, status))
        if math.isfinite(value) and (best is None or value > best[0]):
            best = (value, u, converged, message)
    if best is None:
        raise EngineError(
            f"optimise('map'): none of the {len(positions)} starts ended at a point the "
            f"problem can score. problem.failure_summary() says why:\n{problem.failure_summary()}"
        )
    _, mode, converged, message = best
    # The negated objective's Hessian: -(H[log_prob_unconstrained] - H[J]),
    # the first by the backend's autodiff, the second by finite differences
    # of the numpy-path term (module docstring: the term is cheap and smooth).
    curvature = -(
        hessian(realised, mode)
        - finite_difference_hessian(lambda u: jacobian_term(problem, u), mode)
    )
    covariance, refusal = covariance_from_hessian(curvature)
    return build_optimum(
        problem,
        route="map",
        unconstrained=mode,
        covariance=covariance,
        refusal=refusal,
        converged=converged,
        message=message,
        evaluations=evaluations,
        starts=tuple(summaries),
    )


def _torch_fit(
    problem: FittingProblem, realised: Any, u0: np.ndarray, steps: int, learning_rate: float
) -> tuple[np.ndarray, float, bool, str, int]:
    """``torch.optim.LBFGS`` from *u0*, with Adam as the fallback.

    The Jacobian term enters each closure call **linearised at the current
    point** — ``J(u₀) + (u - u₀)·∇J(u₀)`` with ``u₀`` the detached iterate —
    so the closure's value is the objective exactly and its gradient is the
    objective's exactly, which is all L-BFGS's line search reads.

    The fallback is :func:`_quasi_newton_then_adam`'s: Adam from the best
    finite point L-BFGS reached (or the start), then L-BFGS again to polish.
    """
    import torch  # pyrefly: ignore[missing-import]

    calls = 0

    def negative(u: Any) -> Any:
        nonlocal calls
        calls += 1
        here = u.detach().cpu().numpy().astype(float)
        term = torch.as_tensor(jacobian_term(problem, here), dtype=u.dtype)
        slope = torch.as_tensor(jacobian_gradient(problem, here), dtype=u.dtype)
        density = realised.log_prob_unconstrained(u)
        return -(density - term - torch.dot(u - u.detach(), slope))

    def state(u: np.ndarray) -> tuple[float, np.ndarray]:
        v = torch.tensor(u, dtype=torch.float64, requires_grad=True)
        loss = negative(v)
        if not bool(torch.isfinite(loss)):
            return -math.inf, np.full(u.size, math.nan)
        (grad,) = torch.autograd.grad(loss, v)
        return -float(loss.detach()), grad.detach().cpu().numpy().astype(float)

    def quasi_newton(u: np.ndarray) -> tuple[np.ndarray, str]:
        v = torch.tensor(u, dtype=torch.float64, requires_grad=True)
        lbfgs = torch.optim.LBFGS(
            [v],
            lr=1.0,
            max_iter=steps,
            tolerance_grad=_GRADIENT_TOLERANCE * 1e-2,
            tolerance_change=1e-12,
            history_size=20,
            line_search_fn="strong_wolfe",
        )

        def closure() -> Any:
            lbfgs.zero_grad()
            loss = negative(v)
            loss.backward()
            return loss

        try:
            lbfgs.step(closure)
        except RuntimeError as error:  # a non-finite line search, typically
            return u, f"torch.optim.LBFGS raised {error}"
        return v.detach().cpu().numpy().astype(float), "torch.optim.LBFGS"

    def adam(u: np.ndarray) -> np.ndarray:
        v = torch.tensor(u, dtype=torch.float64, requires_grad=True)
        optimiser = torch.optim.Adam([v], lr=learning_rate)
        for _ in range(steps):
            optimiser.zero_grad()
            loss = negative(v)
            if not bool(torch.isfinite(loss)):
                break
            loss.backward()
            optimiser.step()
        return v.detach().cpu().numpy().astype(float)

    u, value, converged, message = _quasi_newton_then_adam(
        np.asarray(u0, dtype=float),
        state,
        quasi_newton,
        adam,
        "torch.optim.Adam",
        learning_rate,
        steps,
    )
    return u, value, converged, message, calls


def _quasi_newton_then_adam(
    u0: np.ndarray,
    state: Callable[[np.ndarray], tuple[float, np.ndarray]],
    quasi_newton: Callable[[np.ndarray], tuple[np.ndarray, str]],
    adam: Callable[[np.ndarray], np.ndarray],
    adam_name: str,
    learning_rate: float,
    steps: int,
) -> tuple[np.ndarray, float, bool, str]:
    """The native routes' one strategy: quasi-Newton, Adam if it fails, quasi-Newton again.

    A quasi-Newton method started from a prior draw takes its first step with
    an identity inverse Hessian, which on a posterior whose curvature is
    ``1e4`` in some direction is a step of thousands of units — straight into a
    saturated bijection where the density is ``-inf``, after which the line
    search gives up. Adam's step is bounded by its learning rate whatever the
    gradient, so it walks the same start into the bulk reliably; but it never
    *converges* by a gradient test there, because it oscillates at the scale
    of its learning rate. So the fallback is Adam to get close and the
    quasi-Newton method again to finish, and the message records every leg.
    """
    u, name = quasi_newton(u0)
    value, grad = state(u)
    if _gradient_converged(grad, value):
        return u, value, True, f"{name}: converged"
    restart = u if math.isfinite(value) and np.all(np.isfinite(u)) else u0
    walked = adam(restart)
    walked_value, _ = state(walked)
    if not math.isfinite(walked_value):
        walked = restart
    polished, polish_name = quasi_newton(walked)
    polished_value, polished_grad = state(polished)
    if not math.isfinite(polished_value):
        polished, (polished_value, polished_grad) = walked, state(walked)
    converged = _gradient_converged(polished_grad, polished_value)
    message = (
        f"{name} did not converge (max |gradient| above {_GRADIENT_TOLERANCE}); fallback: "
        f"{adam_name}, lr={learning_rate}, {steps} steps, then {polish_name} again: "
        f"{'converged' if converged else 'not converged'}"
    )
    return polished, polished_value, converged, message


def _torch_hessian(realised: Any, mode: np.ndarray) -> np.ndarray:
    import torch  # pyrefly: ignore[missing-import]

    point = torch.tensor(np.asarray(mode, dtype=float), dtype=torch.float64)
    hessian = torch.autograd.functional.hessian(realised.log_prob_unconstrained, point)
    return np.asarray(hessian.detach().cpu().numpy(), dtype=float)


def _jax_objective(problem: FittingProblem, realised: Any) -> Callable[[Any], Any]:
    """The negated objective as a jax function, the Jacobian term by callback.

    ``jax.scipy.optimize.minimize`` traces its function, so the numpy-path
    term goes in through :func:`jax.pure_callback` with a custom JVP whose
    tangent is ``∇J · t`` — linear in ``t``, so reverse mode can transpose it
    and ``jax.value_and_grad`` sees the exact gradient.
    """
    import jax  # pyrefly: ignore[missing-import]
    import jax.numpy as jnp  # pyrefly: ignore[missing-import]

    def value(u: Any) -> np.ndarray:
        return np.asarray(jacobian_term(problem, np.asarray(u, dtype=float)), dtype=np.float64)

    def slope(u: Any) -> np.ndarray:
        return np.asarray(jacobian_gradient(problem, np.asarray(u, dtype=float)), dtype=np.float64)

    @jax.custom_jvp
    def term(u: Any) -> Any:
        return jax.pure_callback(value, jax.ShapeDtypeStruct((), jnp.float64), u)

    @term.defjvp
    def _term_jvp(primals: Any, tangents: Any) -> Any:
        (u,) = primals
        (t,) = tangents
        gradient = jax.pure_callback(slope, jax.ShapeDtypeStruct(u.shape, jnp.float64), u)
        return term(u), jnp.dot(gradient, t)

    def negative(u: Any) -> Any:
        return -(realised.log_prob_unconstrained(u) - term(u))

    return negative


def _jax_fit(
    problem: FittingProblem, realised: Any, u0: np.ndarray, steps: int, learning_rate: float
) -> tuple[np.ndarray, float, bool, str, int]:
    """``jax.scipy.optimize.minimize(method="BFGS")`` from *u0*, ``optax.adam`` as the fallback.

    The strategy is :func:`_quasi_newton_then_adam`'s, as on torch.
    """
    import jax  # pyrefly: ignore[missing-import]
    import jax.numpy as jnp  # pyrefly: ignore[missing-import]
    import optax  # pyrefly: ignore[missing-import]
    from jax.scipy.optimize import minimize  # pyrefly: ignore[missing-import]

    negative = _jax_objective(problem, realised)
    value_and_grad = jax.jit(jax.value_and_grad(negative))
    calls = 0

    def state(u: np.ndarray) -> tuple[float, np.ndarray]:
        nonlocal calls
        calls += 1
        loss, grad = value_and_grad(jnp.asarray(u, dtype=jnp.float64))
        value = -float(loss)
        return (value if math.isfinite(value) else -math.inf), np.asarray(grad, dtype=float)

    def quasi_newton(u: np.ndarray) -> tuple[np.ndarray, str]:
        nonlocal calls
        result = minimize(
            negative, jnp.asarray(u, dtype=jnp.float64), method="BFGS", options={"maxiter": steps}
        )
        calls += int(result.nfev) + int(result.njev)
        name = f"jax.scipy.optimize.minimize('BFGS') (status {int(result.status)})"
        return np.asarray(result.x, dtype=float), name

    optimiser = optax.adam(learning_rate)

    @jax.jit
    def step(params: Any, opt_state: Any) -> Any:
        loss, grads = jax.value_and_grad(negative)(params)
        updates, opt_state = optimiser.update(grads, opt_state, params)
        return optax.apply_updates(params, updates), opt_state, loss

    def adam(u: np.ndarray) -> np.ndarray:
        nonlocal calls
        params = jnp.asarray(u, dtype=jnp.float64)
        opt_state = optimiser.init(params)
        for _ in range(steps):
            candidate, candidate_state, loss = step(params, opt_state)
            calls += 1
            if not bool(jnp.isfinite(loss)):
                break
            params, opt_state = candidate, candidate_state
        return np.asarray(params, dtype=float)

    u, value, converged, message = _quasi_newton_then_adam(
        np.asarray(u0, dtype=float), state, quasi_newton, adam, "optax.adam", learning_rate, steps
    )
    return u, value, converged, message, calls


def _jax_hessian(realised: Any, mode: np.ndarray) -> np.ndarray:
    import jax  # pyrefly: ignore[missing-import]
    import jax.numpy as jnp  # pyrefly: ignore[missing-import]

    point = jnp.asarray(np.asarray(mode, dtype=float), dtype=jnp.float64)
    return np.asarray(jax.hessian(realised.log_prob_unconstrained)(point), dtype=float)


def _vi_route(problem: FittingProblem, positions: np.ndarray, options: dict[str, Any]) -> Optimum:
    raise EngineError("optimise(method='vi') lands in W6.7 unit 6.")


def _hessian_calls(n: int) -> int:
    """Objective evaluations :func:`finite_difference_hessian` makes in *n* dimensions."""
    return 1 + 2 * n + 2 * n + 4 * (n * (n - 1) // 2)
