"""Population inference by reweighting archived single-object fits (W5.13).

Design horizon (b) (``DEVELOPMENT_PLAN.md`` §5, ``results.md`` §14 "Phase 5
(population inference)"): given a collection of *already-fitted* runs that
share a parameter declaration, recover population-level hyperparameters
without ever re-fitting a single object, by importance-reweighting the
draws each run already stored. Hogg, Myers & Bovy (2010) is the reference
construction; what follows is that construction restated against exactly
the columns ``results.md`` §4/§6/§9 guarantees every run carries.

The mathematics
----------------
Object *i* has stored draws :math:`\\theta_{ik}`, :math:`k = 1 \\ldots K_i`,
of the *named* parameter, from a run whose free vector is in general
:math:`(\\theta, \\varphi)` — the named parameter plus every other (nuisance)
parameter the object's own model declares. The run carries that run's stored
``log_prior`` and ``log_likelihood`` per draw (every run, unconditionally —
``results.md`` §6) — the **joint** interim log-prior
:math:`\\log \\pi_0^{\\text{joint}}(\\theta_{ik}, \\varphi_{ik})` over the whole
free vector, not the named parameter alone — and, since W5.0,
``sample_stats.proposal_log_density`` when the run came from an approximate
engine (``ampere_approximation != "none"``: a
:class:`~ampere.inference.VIEngine` guide or an
:class:`~ampere.inference.SBIEngine` density estimator).

**Turning stored draws into (weighted) posterior draws.** For an *exact*
sampler (an ensemble, NUTS, nested sampling — anything whose
``ampere_approximation`` is ``"none"``) the stored :math:`(\\theta_{ik},
\\varphi_{ik})` already **are** draws from the single-object posterior under
that run's *joint* interim prior, so each carries uniform weight
:math:`w_{ik} = 1 / K_i`. For an *approximate* engine the stored draws are
instead from the fitted proposal :math:`q`, and W5.0's contract is exactly
what makes them usable here: the *self-normalised* importance weights

.. math::

    w_{ik} \\propto \\exp\\!\\big(\\text{log\\_prior}_{ik} +
    \\text{log\\_likelihood}_{ik} - \\text{proposal\\_log\\_density}_{ik}\\big),
    \\qquad \\sum_k w_{ik} = 1

turn the proposal draws into (weighted) draws from the same single-object
posterior :math:`\\pi_0^{\\text{joint}}(\\theta, \\varphi) \\, L_i(\\theta, \\varphi)`
the exact-sampler case already has — the **joint** ``log_prior`` column is
the correct and only correct term here. A run whose ``ampere_approximation``
is not ``"none"`` and has no ``proposal_log_density`` column is refused by
name (§9's contract is not optional; a run built before W5.0 lacks the
column it needs and there is no honest number to substitute for it).

**The population likelihood, and why it needs the *marginal* interim
prior.** The population model replaces the interim prior on :math:`\\theta`
alone with :math:`p(\\theta \\mid \\alpha)`; it says nothing about
:math:`\\varphi`, whose interim prior :math:`\\pi(\\varphi)` (declared,
unconditional on :math:`\\alpha`) is carried over unchanged. Since every
object's free parameters are declared independently (§4.1), the joint
interim prior factorises, :math:`\\pi_0^{\\text{joint}}(\\theta, \\varphi) =
\\pi_0(\\theta)\\,\\pi(\\varphi)`, and the standard importance-sampling identity
(the interim prior divided back out, one draw at a time) gives

.. math::

    p(d_i \\mid \\alpha) \\;\\propto\\; \\sum_k w_{ik} \\,
    \\frac{p(\\theta_{ik} \\mid \\alpha)\\,\\pi(\\varphi_{ik})}
         {\\pi_0(\\theta_{ik})\\,\\pi(\\varphi_{ik})}
    \\;=\\; \\sum_k w_{ik} \\,
    \\frac{p(\\theta_{ik} \\mid \\alpha)}{\\pi_0(\\theta_{ik})}

— :math:`\\pi(\\varphi_{ik})` cancels exactly, for every draw, leaving only
:math:`\\pi_0(\\theta)`, the named parameter's own **marginal** interim
prior, in the denominator. Dividing by the *joint* ``log_prior`` column
instead (an earlier version of this module did) leaves an uncancelled
:math:`1 / \\pi(\\varphi_{ik})` factor that varies draw to draw whenever the
object's model has more than one free parameter, biasing the population
posterior; a single-parameter model has no :math:`\\varphi` to cancel, which
is why the mistake is invisible on a toy with one parameter. Because a run's
provenance records only each parameter's *name* and a hash of its
declaration (``results.md`` §9 — the full ``PriorSpec`` is deliberately not
written, for the same "do not materialise what §16 says not to" reason
``free_labels()`` is not either), :math:`\\pi_0` cannot be read back off a
run: :func:`fit_population` takes it as an explicit, required
``interim_prior`` argument instead of guessing at it, and refuses (via
Python's own required-argument mechanism) rather than falling back to the
joint column when one is not supplied. The population posterior this module
samples is

.. math::

    \\log p(\\alpha \\mid \\text{data}) = \\log p(\\alpha) +
    \\sum_i \\log \\Big[ \\sum_k w_{ik} \\,
    \\frac{p(\\theta_{ik} \\mid \\alpha)}{\\pi_0(\\theta_{ik})} \\Big]

computed entirely from stored columns plus the one supplied prior:
:func:`fit_population` builds this once per object (:class:`_PreparedObject`,
independent of :math:`\\alpha`) and re-evaluates only the
:math:`p(\\theta_{ik} \\mid \\alpha)` term inside the ``emcee`` sampler's
log-probability.

**Effective sample size, and the refusal.** At a given :math:`\\alpha` (this
module uses the population posterior mean, one evaluation, stated in
:func:`fit_population`'s docstring rather than repeated per object) the
per-object combined weights
:math:`v_{ik} \\propto w_{ik}\\, p(\\theta_{ik}\\mid\\alpha) / \\pi_0(\\theta_{ik})`,
normalised to sum to 1, give the effective sample size
:math:`\\mathrm{ESS}_i = 1 / \\sum_k v_{ik}^2`: how many of object *i*'s
draws are actually doing the reweighting work at that population. An object
whose interim prior or proposal put its mass somewhere the population model
disfavours (or the reverse) collapses this to a handful of draws — the
reweighted answer for that object would then be controlled by one or two
points ``emcee`` happened to visit, which is not something to report
silently. :func:`fit_population` computes every object's ESS at the fitted
posterior mean and refuses by name, naming the floor and the worst object's
index, when the minimum falls below ``ess_floor`` (default 20, an argument).

The population model as a declaration
--------------------------------------
:class:`PopulationModel` is a small protocol: ``log_density(theta, alpha)``
and a sequence of ``hyperparameters`` — ordinary :class:`~ampere.core.
parameter.Parameter` objects, each with its own prior. The hyperprior
:math:`p(\\alpha)` is therefore *not* a separate method a model author has
to reimplement: it is read generically off ``hyperparameters`` (each
``Parameter.prior.logpdf`` summed, in the same "a declared prior is a
declared prior" spirit as ``ParameterSet.lnprior`` for an ordinary fit), the
same way this module asks nothing bespoke of a "hyper-fit" that an ordinary
one does not already ask of a fit. :class:`GaussianPopulationModel` is the
one instance shipped: :math:`\\theta \\sim \\mathcal{N}(\\mu, \\tau^2)`,
independently, over **one** named archived-fit parameter. Several
independent parameters extend this directly — one
:class:`GaussianPopulationModel` per parameter, its log-densities summed —
but that composition is not built here because nothing in Phase 5 needs it
yet; §13's "what is not here" below is explicit about it.

The reader is a protocol
-------------------------
:class:`RunColumns` is what this module needs from one archived run — the
named parameter's flattened draws, the run's **joint** ``log_prior`` and
``log_likelihood``, ``proposal_log_density`` (or ``None``), and the run's
provenance ``attrs`` — expressed as a :class:`typing.Protocol` rather than a
base class, because horizon (b)'s "buildable entirely on stored files" claim
must not secretly mean "on an ``xarray.DataTree`` and nothing else": a
columnar store that never materialises a ``DataTree`` at all can satisfy
this protocol and be reweighted here unchanged. Two implementations ship:
:class:`DataTreeRunColumns` (a run already in memory — whatever produced
it) and :class:`NetCDFRunColumns` plus :func:`runs_from_netcdf_directory`
(a directory of archived ``.nc`` files, read with
:func:`~ampere.results.from_netcdf`).

What is not here
-----------------
No per-object nuisance re-sampling (an object's *other* free parameters are
never touched — only the one named parameter is read back and reweighted).
No multi-level populations (populations of populations; ``parameters.md``
§12.2's nested plates are the core-side prerequisite and remain deferred).
No selection function (a population inferred this way is only ever a
population of the objects that were fit, never corrected for which objects
were selected to be fit). No joint fit across objects sharing the
population prior as a single sampler — that is design horizon (b)'s sibling
route, **W5.12**, ruled to land after this item (D4, 2026-09-15) precisely
because it needs nothing this module builds and this module needs nothing
it builds; the two are cross-checked against each other where both exist,
never merged into one code path.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import numpy as np

from ampere.core.exceptions import OptionalDependencyError, ResultsError
from ampere.core.parameter import Parameter, Prior, PriorSpec, prior_from_spec
from ampere.core.parameter import log_density as prior_log_density

from .emission import from_netcdf
from .provenance import ATTR_PREFIX, canonical_json

__all__ = [
    "DEFAULT_ESS_FLOOR",
    "DataTreeRunColumns",
    "GaussianPopulationModel",
    "NetCDFRunColumns",
    "PopulationModel",
    "RunColumns",
    "effective_sample_sizes",
    "fit_population",
    "runs_from_netcdf_directory",
]

#: The ESS floor :func:`fit_population` refuses below, absent an explicit
#: ``ess_floor=``. Not derived from anything principled — it is the "this
#: object's reweighted contribution is a handful of draws, not a posterior"
#: order of magnitude — which is why it stays an argument rather than a
#: constant baked into the refusal.
DEFAULT_ESS_FLOOR = 20.0


def _require_arviz() -> Any:
    """Import arviz on use, never on import — the same reasoning every sibling module gives."""
    try:
        import arviz
    except ImportError as error:  # pragma: no cover - exercised by the minimal-install job
        raise OptionalDependencyError(
            "arviz",
            context="building a population-hyperparameter run's xarray.DataTree "
            "(ampere.results.population, W5.13). arviz is a base dependency of ampere, so this "
            "environment is incomplete rather than merely missing an extra",
        ) from error
    return arviz


# ---------------------------------------------------------------------------
# The population model
# ---------------------------------------------------------------------------


@runtime_checkable
class PopulationModel(Protocol):
    """A population prior over one archived-fit parameter, declared like any other prior.

    ``hyperparameters`` fixes the order :func:`fit_population`'s ``alpha``
    vector uses (the same convention :attr:`~ampere.core.parameter.
    ParameterSet.free_names` fixes for an ordinary fit's flat vector) and is
    also where the hyperprior lives: each is an ordinary
    :class:`~ampere.core.parameter.Parameter` with its own ``prior``, so
    :math:`\\log p(\\alpha)` is computed generically
    (:func:`_log_hyperprior`) rather than by a method every model author
    would otherwise have to reimplement identically.
    """

    hyperparameters: Sequence[Parameter]

    def log_density(self, theta: np.ndarray, alpha: np.ndarray) -> np.ndarray:
        """:math:`\\log p(\\theta \\mid \\alpha)`, elementwise over *theta*."""
        ...


@dataclasses.dataclass(frozen=True)
class GaussianPopulationModel:
    """:math:`\\theta \\sim \\mathcal{N}(\\mu, \\tau^2)`, over one named archived-fit parameter.

    Parameters
    ----------
    mu, tau
        Free :class:`~ampere.core.parameter.Parameter`\\ s (a prior, not
        fixed or deferred) — the population mean and standard deviation.
        ``tau``'s prior should put no mass at or below zero (a
        half-normal or a log-uniform, say): :meth:`log_density` returns
        ``-inf`` for a non-positive ``tau`` regardless, so an ``emcee``
        walker that proposes one is rejected rather than crashing
        ``scipy.stats.norm``, but a prior that already excludes it gives a
        cleaner hyperprior surface.

    Examples
    --------
    >>> import scipy.stats as st
    >>> model = GaussianPopulationModel(
    ...     Parameter("mu", st.norm(0.0, 5.0)),
    ...     Parameter("tau", st.halfnorm(scale=2.0)),
    ... )
    >>> import numpy as np
    >>> float(model.log_density(np.array([1.0, 1.2]), np.array([1.0, 0.3]))[0])  # doctest: +SKIP
    """

    mu: Parameter
    tau: Parameter

    def __post_init__(self) -> None:
        if not (self.mu.is_free and self.tau.is_free):
            raise ResultsError(
                "GaussianPopulationModel needs mu and tau declared as free Parameters (a prior, "
                "not fixed or deferred) -- that is what makes them hyperparameters with a "
                "hyperprior, the same distinction Parameter draws for an ordinary fit."
            )

    @property
    def hyperparameters(self) -> tuple[Parameter, Parameter]:
        return (self.mu, self.tau)

    def log_density(self, theta: np.ndarray, alpha: np.ndarray) -> np.ndarray:
        mu, tau = float(alpha[0]), float(alpha[1])
        values = np.asarray(theta, dtype=float)
        if not (tau > 0.0):
            return np.full(values.shape, -np.inf)
        variance = tau * tau
        return -0.5 * (((values - mu) ** 2) / variance + math.log(2.0 * math.pi * variance))


def _log_hyperprior(model: PopulationModel, alpha: np.ndarray) -> float:
    """:math:`\\log p(\\alpha)`, generically: the sum of each hyperparameter's own prior."""
    total = 0.0
    for parameter, value in zip(model.hyperparameters, alpha, strict=True):
        density = float(parameter.prior.logpdf(value))
        if not math.isfinite(density):
            return -math.inf
        total += density
    return total


# ---------------------------------------------------------------------------
# The reader: a protocol, and two implementations
# ---------------------------------------------------------------------------


@runtime_checkable
class RunColumns(Protocol):
    """The per-draw columns this module needs from one archived run.

    Horizon (b)'s promise is "buildable entirely on stored files"; expressed
    as a base class this would quietly narrow to "stored as an
    ``xarray.DataTree``". A :class:`typing.Protocol` keeps the promise
    literal: anything with these four columns and these attrs — a columnar
    store that never builds a ``DataTree`` at all — is reweightable here.
    """

    @property
    def attrs(self) -> Mapping[str, Any]:
        """This run's root provenance attrs (``ampere_spec_hash``, ``ampere_approximation``, …)."""
        ...

    def parameter_draws(self, name: str) -> np.ndarray:
        """The named parameter's flattened draws, one value per stored draw."""
        ...

    @property
    def log_prior(self) -> np.ndarray:
        """The stored ``sample_stats.log_prior``, flattened — this run's interim prior, per draw."""
        ...

    @property
    def log_likelihood(self) -> np.ndarray:
        """The stored ``sample_stats.log_likelihood``, flattened."""
        ...

    @property
    def proposal_log_density(self) -> np.ndarray | None:
        """``sample_stats.proposal_log_density``, flattened, or ``None`` if the run has none."""
        ...


def _group(tree: Any, group: str) -> Any:
    try:
        node = tree[group]
    except KeyError as exc:
        raise ResultsError(
            f"this run has no {group!r} group -- was it built by ampere.results.emit (or read "
            f"back with ampere.results.from_netcdf)?"
        ) from exc
    return getattr(node, "dataset", node)


@dataclasses.dataclass(frozen=True)
class DataTreeRunColumns:
    """:class:`RunColumns` over a run already held in memory as an :class:`xarray.DataTree`.

    Accepts anything :func:`~ampere.results.emit` or
    :func:`~ampere.results.from_netcdf` returns; does not care which.
    """

    tree: Any

    @property
    def attrs(self) -> Mapping[str, Any]:
        return dict(self.tree.attrs)

    def parameter_draws(self, name: str) -> np.ndarray:
        posterior = _group(self.tree, "posterior")
        if name not in posterior:
            raise ResultsError(
                f"parameter {name!r} is not in this run's posterior group; it has "
                f"{sorted(str(key) for key in posterior.data_vars)}."
            )
        return np.asarray(posterior[name].values, dtype=float).reshape(-1)

    @property
    def log_prior(self) -> np.ndarray:
        return self._stat("log_prior")

    @property
    def log_likelihood(self) -> np.ndarray:
        return self._stat("log_likelihood")

    @property
    def proposal_log_density(self) -> np.ndarray | None:
        stats = _group(self.tree, "sample_stats")
        if "proposal_log_density" not in stats:
            return None
        return np.asarray(stats["proposal_log_density"].values, dtype=float).reshape(-1)

    def _stat(self, name: str) -> np.ndarray:
        stats = _group(self.tree, "sample_stats")
        if name not in stats:
            raise ResultsError(f"this run's sample_stats has no {name!r} column.")
        return np.asarray(stats[name].values, dtype=float).reshape(-1)


@dataclasses.dataclass(frozen=True)
class NetCDFRunColumns:
    """:class:`RunColumns` over one archived run, read from a netCDF file on disk.

    Read eagerly, once, at construction: an archived single-object run is
    small, so nothing is bought by deferring the read, and a bad path fails
    at construction rather than partway through a population fit.
    """

    path: Path
    _columns: DataTreeRunColumns = dataclasses.field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "_columns", DataTreeRunColumns(from_netcdf(self.path)))

    @property
    def attrs(self) -> Mapping[str, Any]:
        return self._columns.attrs

    def parameter_draws(self, name: str) -> np.ndarray:
        return self._columns.parameter_draws(name)

    @property
    def log_prior(self) -> np.ndarray:
        return self._columns.log_prior

    @property
    def log_likelihood(self) -> np.ndarray:
        return self._columns.log_likelihood

    @property
    def proposal_log_density(self) -> np.ndarray | None:
        return self._columns.proposal_log_density


def runs_from_netcdf_directory(directory: str | Path, *, pattern: str = "*.nc") -> list[RunColumns]:
    """Every archived run in *directory*, one :class:`NetCDFRunColumns` per file.

    Sorted by filename, so a population fit built from a directory is
    reproducible from the directory alone rather than from whatever order
    the filesystem happens to hand back.
    """
    root = Path(directory)
    paths = sorted(root.glob(pattern))
    if not paths:
        raise ResultsError(f"no files matching {pattern!r} in {root}; nothing to reweight.")
    return [NetCDFRunColumns(path) for path in paths]


# ---------------------------------------------------------------------------
# Per-object preparation: everything that does not depend on alpha
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _PreparedObject:
    """One object's columns, reduced to what :func:`_log_posterior` needs, computed once."""

    theta: np.ndarray
    marginal_log_prior: np.ndarray  # log pi_0(theta), the *named parameter's own* interim prior
    log_w: np.ndarray  # self-normalised: sums to 1 in linear space
    problem_hash: str
    spec_hash: str


def _as_prior(interim_prior: Prior | PriorSpec) -> Prior:
    """Accept either a live prior or its neutral :class:`PriorSpec` description."""
    if isinstance(interim_prior, PriorSpec):
        return prior_from_spec(interim_prior)
    return interim_prior


def _self_normalised_log_weights(run: RunColumns, joint_log_prior: np.ndarray) -> np.ndarray:
    """Per-draw log importance weights, normalised to sum to 1 in linear space.

    Exact samplers (``ampere_approximation == "none"``): the stored draws are
    already posterior draws under the *joint* interim prior, so the weights
    are uniform. Approximate engines (W5.0): ``exp(joint_log_prior +
    log_likelihood - proposal_log_density)`` self-normalised is the
    importance weight that turns proposal draws into (weighted) posterior
    draws first, per the module docstring. This is the one place the run's
    **joint** ``log_prior`` column is the correct term -- unlike the
    marginal :math:`\\pi_0(\\theta)` :func:`_log_object_evidence` needs.
    """
    from scipy.special import logsumexp

    approximation = str(run.attrs.get(f"{ATTR_PREFIX}approximation", "none"))
    if approximation == "none":
        raw = np.zeros_like(joint_log_prior)
    else:
        proposal = run.proposal_log_density
        if proposal is None:
            raise ResultsError(
                f"this run's {ATTR_PREFIX}approximation is {approximation!r}, so its stored "
                f"draws are from a proposal rather than the target and need "
                f"sample_stats.proposal_log_density (results.md §9, W5.0) to become importance "
                f"weights; this run has no such column. Refusing rather than treating an "
                f"approximate run's draws as if they were exact posterior draws."
            )
        raw = joint_log_prior + run.log_likelihood - proposal
    return raw - logsumexp(raw)


def _prepare_objects(
    runs: Sequence[RunColumns], parameter: str, interim_prior: Prior
) -> list[_PreparedObject]:
    """Read every run's columns once, and refuse a mixture this module cannot reweight.

    *interim_prior* is the named parameter's own **marginal** interim prior
    :math:`\\pi_0(\\theta)` -- not the run's stored (joint) ``log_prior``
    column, which is what :func:`_self_normalised_log_weights` alone uses.
    See the module docstring's "why it needs the marginal interim prior".
    """
    if not runs:
        raise ResultsError("fit_population needs at least one run.")
    prepared: list[_PreparedObject] = []
    spec_hashes: dict[str, list[int]] = {}
    for index, run in enumerate(runs):
        attrs = run.attrs
        spec_hash = attrs.get(f"{ATTR_PREFIX}spec_hash")
        if spec_hash is None:
            raise ResultsError(
                f"run {index} has no {ATTR_PREFIX}spec_hash attribute; is it an ampere run "
                f"(ampere.results.emit or from_netcdf)?"
            )
        spec_hashes.setdefault(str(spec_hash), []).append(index)
        theta = run.parameter_draws(parameter)
        joint_log_prior = run.log_prior
        if theta.shape != joint_log_prior.shape:
            raise ResultsError(
                f"run {index}: parameter {parameter!r} has {theta.size} draws but log_prior has "
                f"{joint_log_prior.size}; they must come from the same run's columns."
            )
        log_w = _self_normalised_log_weights(run, joint_log_prior)
        marginal_log_prior = prior_log_density(interim_prior, theta)
        prepared.append(
            _PreparedObject(
                theta=theta,
                marginal_log_prior=marginal_log_prior,
                log_w=log_w,
                problem_hash=str(attrs.get(f"{ATTR_PREFIX}problem_hash", "")),
                spec_hash=str(spec_hash),
            )
        )
    if len(spec_hashes) > 1:
        raise ResultsError(
            f"these runs do not share {ATTR_PREFIX}spec_hash ({sorted(spec_hashes)}); reweighting "
            f"together assumes the same parameter declaration (and hence the same interim prior "
            f"convention) for every object. Refusing rather than silently mixing declarations."
        )
    return prepared


def _log_object_evidence(
    obj: _PreparedObject, model: PopulationModel, alpha: np.ndarray
) -> tuple[float, np.ndarray]:
    """``log p(d_i | alpha)`` and the per-draw normalised weights the ESS needs."""
    from scipy.special import logsumexp

    log_c = obj.log_w + model.log_density(obj.theta, alpha) - obj.marginal_log_prior
    log_z = float(logsumexp(log_c))
    if not math.isfinite(log_z):
        return -math.inf, np.full(log_c.shape, 1.0 / log_c.size)
    return log_z, np.exp(log_c - log_z)


def effective_sample_sizes(
    prepared: Sequence[_PreparedObject], model: PopulationModel, alpha: np.ndarray
) -> np.ndarray:
    """Per-object effective sample size of the reweighted draws, at this ``alpha``."""
    sizes = np.empty(len(prepared), dtype=float)
    for index, obj in enumerate(prepared):
        _, normalised = _log_object_evidence(obj, model, alpha)
        sizes[index] = 1.0 / float(np.sum(normalised**2))
    return sizes


def _log_posterior(
    alpha: np.ndarray, prepared: Sequence[_PreparedObject], model: PopulationModel
) -> float:
    hyperprior = _log_hyperprior(model, alpha)
    if not math.isfinite(hyperprior):
        return -math.inf
    total = hyperprior
    for obj in prepared:
        log_z, _ = _log_object_evidence(obj, model, alpha)
        if not math.isfinite(log_z):
            return -math.inf
        total += log_z
    return total


# ---------------------------------------------------------------------------
# emcee conventions, restated for a hyperparameter vector rather than a problem
# ---------------------------------------------------------------------------


def _default_walkers(n_hyper: int) -> int:
    """Four per hyperparameter, at least eight — ``ampere.inference.engine``'s own default."""
    return max(8, 4 * n_hyper)


def _check_ensemble(walkers: int, n_hyper: int) -> int:
    if walkers < 2:
        raise ResultsError(f"fit_population needs at least 2 walkers, got {walkers}.")
    if walkers % 2:
        raise ResultsError(
            f"emcee splits its ensemble in half each step, so the number of walkers must be "
            f"even; got {walkers}."
        )
    if walkers < 2 * n_hyper:
        raise ResultsError(
            f"fit_population needs at least 2 x n_hyper = {2 * n_hyper} walkers to span a "
            f"{n_hyper}-dimensional hyperparameter space, got {walkers}."
        )
    return int(walkers)


def _kept(total: int, burn_in: int, thin: int) -> None:
    if total < 1:
        raise ResultsError(f"fit_population needs at least one step, got {total}.")
    if thin < 1:
        raise ResultsError(f"fit_population's thinning factor must be at least 1, got {thin}.")
    if burn_in < 0:
        raise ResultsError(f"fit_population's burn-in cannot be negative, got {burn_in}.")
    if burn_in >= total:
        raise ResultsError(
            f"fit_population was asked to discard {burn_in} of {total} step(s) as burn-in, which "
            f"leaves nothing to emit."
        )
    if (total - burn_in) // thin < 1:
        raise ResultsError(
            f"fit_population was asked to thin {total - burn_in} kept step(s) by {thin}, which "
            f"leaves nothing to emit."
        )


def _initial_positions(
    model: PopulationModel, walkers: int, rng: np.random.Generator
) -> np.ndarray:
    """*walkers* draws from the hyperprior.

    The same "start where the prior says" convention
    :meth:`~ampere.inference.engine.Engine.initial_positions` uses for an
    ordinary fit.
    """
    positions = np.empty((walkers, len(model.hyperparameters)), dtype=float)
    for column, parameter in enumerate(model.hyperparameters):
        positions[:, column] = parameter.prior.rvs(size=walkers, random_state=rng)
    return positions


# ---------------------------------------------------------------------------
# The entry point
# ---------------------------------------------------------------------------


def fit_population(
    runs: Sequence[RunColumns],
    parameter: str,
    model: PopulationModel,
    interim_prior: Prior | PriorSpec,
    *,
    walkers: int | None = None,
    steps: int = 3000,
    burn_in: int = 1000,
    thin: int = 1,
    ess_floor: float = DEFAULT_ESS_FLOOR,
    seed: int | None = None,
) -> Any:
    """Sample the population posterior over ``model``'s hyperparameters with ``emcee``.

    Parameters
    ----------
    runs
        One :class:`RunColumns` per object, sharing ``ampere_spec_hash``
        (refused by name otherwise) and each holding draws of *parameter*.
    parameter
        The archived-fit parameter's name (e.g. ``"model.theta"``), read
        back from every run's ``posterior`` group.
    model
        The population model — see :class:`GaussianPopulationModel`.
    interim_prior
        The named parameter's own **marginal** interim prior
        :math:`\\pi_0(\\theta)` — the prior it was declared with in every
        input run (guaranteed identical across them by the shared
        ``ampere_spec_hash`` refusal above). A frozen ``scipy.stats``
        distribution (a :class:`~ampere.core.parameter.Prior`) or a
        :class:`~ampere.core.parameter.PriorSpec`. Required, and not read
        off the runs themselves: a run's provenance records only each
        parameter's name and a hash of its declaration, not the declaration
        itself (``results.md`` §9), so there is nothing here to read back
        from — see the module docstring's "why it needs the marginal
        interim prior" for the identity this divides out and why the run's
        stored (joint) ``log_prior`` column is the wrong thing to divide by
        whenever the object's model has more than one free parameter.
    walkers, steps, burn_in, thin
        ``emcee.EnsembleSampler`` settings, in the same sense
        :class:`~ampere.inference.EmceeEngine` uses them; ``walkers``
        defaults to four per hyperparameter (at least eight).
    ess_floor
        Refuse, by name, if any object's effective sample size at the
        population posterior mean falls below this (default
        :data:`DEFAULT_ESS_FLOOR`). Evaluated **once**, at the posterior
        mean of ``alpha`` after sampling — not per MCMC step, which would
        cost one ``effective_sample_sizes`` call per proposal for a number
        this module only needs to report once a fit exists.
    seed
        Seeds both the walkers' start points and ``emcee``'s own
        randomness; ``None`` means neither is reproducible.

    Returns
    -------
    xarray.DataTree
        A ``posterior`` group over ``model.hyperparameters``' names, a
        ``sample_stats`` group holding ``lp``, and root attrs
        ``ampere_population_runs`` (each input run's ``ampere_problem_hash``
        and ``ampere_spec_hash``), ``ampere_population_model`` (the model's
        class and hyperparameter names), ``ampere_population_ess`` (the
        per-object effective sample size at the posterior mean, and its
        min/median across objects) and ``ampere_population_parameter``
        (*parameter*). Not a change to ``PROVENANCE_SCHEMA_VERSION``: this
        is a new derived product, not a change to what a single-object run
        stores.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        A spec-hash mixture, a missing parameter or column, an approximate
        run with no ``proposal_log_density``, or an effective sample size
        collapsed below ``ess_floor``.
    """
    import emcee

    prepared = _prepare_objects(runs, parameter, _as_prior(interim_prior))
    hyperparameters = list(model.hyperparameters)
    n_hyper = len(hyperparameters)
    chosen_walkers = _default_walkers(n_hyper) if walkers is None else int(walkers)
    chosen_walkers = _check_ensemble(chosen_walkers, n_hyper)
    _kept(int(steps), int(burn_in), int(thin))

    rng = np.random.default_rng(seed)
    initial = _initial_positions(model, chosen_walkers, rng)

    sampler = emcee.EnsembleSampler(chosen_walkers, n_hyper, _log_posterior, args=(prepared, model))
    if seed is not None:
        sampler.random_state = np.random.RandomState(int(seed)).get_state()
    sampler.run_mcmc(initial, int(steps), progress=False)

    # emcee stores (step, walker, dim); ArviZ wants (chain, draw, dim).
    chain = np.swapaxes(sampler.get_chain(discard=int(burn_in), thin=int(thin)), 0, 1)
    log_prob = np.swapaxes(sampler.get_log_prob(discard=int(burn_in), thin=int(thin)), 0, 1)

    alpha_mean = chain.reshape(-1, n_hyper).mean(axis=0)
    ess = effective_sample_sizes(prepared, model, alpha_mean)
    if float(np.min(ess)) < float(ess_floor):
        worst = int(np.argmin(ess))
        raise ResultsError(
            f"the effective sample size for object {worst} (of {len(prepared)}) collapsed to "
            f"{ess[worst]:.2f} at the population posterior mean, below the floor {ess_floor}. Its "
            f"stored draws carry almost no weight under this population model, so the reweighted "
            f"answer for that object would be controlled by a handful of draws. Refusing rather "
            f"than reporting a population posterior one under-sampled object secretly controls -- "
            f"re-fit that object with a proposal or interim prior closer to the population, widen "
            f"its run's draw count, or drop it."
        )

    names = [parameter_.name for parameter_ in hyperparameters]
    posterior_variables = {name: chain[..., index] for index, name in enumerate(names)}
    sample_stats = {"lp": log_prob}

    root_attrs = {
        f"{ATTR_PREFIX}population_runs": canonical_json(
            [{"problem_hash": obj.problem_hash, "spec_hash": obj.spec_hash} for obj in prepared]
        ),
        f"{ATTR_PREFIX}population_model": canonical_json(
            {"kind": type(model).__name__, "hyperparameters": names}
        ),
        f"{ATTR_PREFIX}population_ess": canonical_json(
            {
                "per_object": [float(value) for value in ess],
                "min": float(np.min(ess)),
                "median": float(np.median(ess)),
            }
        ),
        f"{ATTR_PREFIX}population_parameter": parameter,
    }

    arviz = _require_arviz()
    tree = arviz.from_dict({"posterior": posterior_variables, "sample_stats": sample_stats})
    tree.attrs.update(root_attrs)
    return tree
