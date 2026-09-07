"""What every engine driver shares, and the rule that keeps them portable.

``inference.md`` §10 states the bet this module cashes: an engine that uses
**only** ``DEVELOPMENT_PLAN.md`` §4.5's surface — ``log_prob``, the
``log_likelihood``/``log_prior`` split, ``prior_transform``, ``simulate``, the
capability flags — "works with the reference backend, with torch, with jax, and
with a legacy black-box model behind a thin adapter, and never knows which".

So the enforcement is a rule about imports, and it is the one thing a reviewer
should check first: **nothing under ``ampere.inference`` may import
``ampere.backends``**, and nothing here branches on which backend built the
problem. The only imports are ``ampere.core`` (the problem surface) and
``ampere.results`` (the one results format). A driver that needed to know more
than that would be evidence against §4.5, not a reason to widen it.

What the base class owns
------------------------
Everything that is the same for emcee, dynesty and zeus, which turns out to be
almost everything except the sampler call itself:

* **the engine check** — :meth:`~ampere.core.dataset.FittingProblem.
  check_engine` at construction, so a likelihood no gradient-free engine can
  marginalise is refused before any sampling happens rather than after;
* **seeds** — every stream derived from the problem's own seed through
  :meth:`~ampere.core.dataset.FittingProblem.rng`, which is ``lowering.md``
  §9.2's policy realised (``inference.md`` §12, ratified 2026-09-03);
* **start points** — drawn from the joint prior, retried until one is inside
  the support;
* **the per-draw record** — one :class:`~ampere.core.dataset.Evaluation` per
  *stored* draw, which is what ``ampere.results`` needs and what
  ``inference.md`` §18 shaped ``Evaluation`` to be;
* **emission** — the ArviZ ``DataTree``, through :func:`ampere.results.emit`,
  every run without exception (``DEVELOPMENT_PLAN.md`` §4.6);
* **failure signalling** — the non-strict path consumed exactly as
  ``inference.md`` §11 declares it, and the aggregate surfaced three ways
  (see :meth:`Engine.run`'s note).

Where the per-draw evaluations come from
----------------------------------------
A sampler calls ``log_prob`` on *proposals* and stores *positions*, and the two
sets are not the same: an ensemble walker that rejects keeps a position it was
scored at several steps ago, and a nested sampler's posterior is a subset of the
points it killed. Recomputing ``problem.evaluate`` over the stored draws
afterwards would be correct and would also double the model evaluations of the
kept half of the run — and for a *stochastic* model it would not even be
correct, because the recomputed log-likelihood is not the one the sampler
accepted on.

So the drivers score through :class:`_EvaluationCache`, which keeps the
``Evaluation`` for each θ it was asked about, keyed by the exact bytes of the
float64 vector, and looks them up again when the run is assembled. Hits are the
normal case; a miss (the cache is bounded, so a long run evicts) falls back to
re-evaluating that draw, which is the honest thing and is counted in the
provenance attrs.

A **gradient-based** driver cannot use that cache at all, and the reason is the
same one that made ``inference.md`` §10a necessary: it scores through the
backend's realisation, which returns a scalar in the backend's own array type
rather than an :class:`~ampere.core.dataset.Evaluation`. So the cache is never
populated by the sampler, and every stored draw used to be decomposed
afterwards on the numpy contract path — a full model evaluation per draw, of a
quantity the backend had just computed. §10a's optional
``log_likelihood_terms`` is the answer, and :meth:`Engine.finish` consumes it
(W2.4 slice 2): a driver that has it hands the decomposition straight in and
``engine_draws_recomputed`` is zero.
"""

from __future__ import annotations

import abc
import math
import warnings
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import Evaluation, FittingProblem
from ampere.results import emit

from .exceptions import EngineError, SamplingFailureWarning

__all__ = ["DEFAULT_CACHE_SIZE", "Engine"]

#: How many scored θ an engine remembers so that the stored draws need not be
#: re-evaluated. An ``Evaluation`` is a handful of floats plus one term per
#: dataset, so 100 000 of them is tens of megabytes — the same order as the
#: emitted run itself, and small beside the model evaluations it saves. Past
#: the bound the oldest go and the affected draws are recomputed; nothing is
#: lost but time.
DEFAULT_CACHE_SIZE = 100_000


class _EvaluationCache:
    """Remember one :class:`~ampere.core.dataset.Evaluation` per scored θ.

    Keyed on ``vector.tobytes()`` of the float64 free-parameter vector, which
    is an *exact* key: two θ that are equal to the last bit share an entry and
    nothing else does. That is what makes it safe — there is no tolerance, so
    there is no way for a lookup to return the evaluation of a nearby point.

    Bounded, and evicted oldest-first (``dict`` preserves insertion order), so
    that a 10⁶-proposal run is a bounded diagnostic rather than a memory leak —
    the same reasoning ``inference.md`` §11 applies to the failure history.
    """

    def __init__(self, problem: FittingProblem, maximum: int) -> None:
        if maximum < 0:
            raise EngineError(f"the evaluation cache cannot hold {maximum} entries.")
        self._problem = problem
        self._maximum = int(maximum)
        self._entries: dict[bytes, Evaluation] = {}
        #: How many times the problem was actually evaluated.
        self.calls = 0
        #: How many stored draws had to be evaluated a second time.
        self.recomputed = 0

    @staticmethod
    def _key(theta: Any) -> tuple[np.ndarray, bytes]:
        vector = np.asarray(theta, dtype=float).reshape(-1)
        return vector, vector.tobytes()

    def evaluate(self, theta: Any) -> Evaluation:
        """Score θ, remembering the result. What the sampler's callback calls."""
        vector, key = self._key(theta)
        evaluation = self._problem.evaluate(vector)
        self.calls += 1
        if self._maximum:
            if len(self._entries) >= self._maximum and key not in self._entries:
                del self._entries[next(iter(self._entries))]
            self._entries[key] = evaluation
        return evaluation

    def lookup(self, theta: Any) -> Evaluation:
        """The evaluation of a *stored* draw, recomputing it if it has been evicted."""
        vector, key = self._key(theta)
        found = self._entries.get(key)
        if found is not None:
            return found
        self.recomputed += 1
        self.calls += 1
        return self._problem.evaluate(vector)


class Engine(abc.ABC):
    """Base class for the gradient-free engine drivers.

    Subclasses supply :attr:`NAME`, whatever settings their sampler takes, and
    a ``run`` method that ends by calling :meth:`finish`. They do not touch
    provenance, emission, seeding or failure reporting — those are here, once,
    which is the same economy ``inference.md`` §10 claims for writing the
    engines against §4.5 only.

    Parameters
    ----------
    problem
        The composed :class:`~ampere.core.dataset.FittingProblem`. Its
        ``check_engine`` is called immediately, so a likelihood this class of
        engine cannot run is refused here rather than thousands of draws later.
    cache_size
        Entries in the evaluation cache; see :data:`DEFAULT_CACHE_SIZE`. Zero
        disables it, at the cost of re-evaluating every stored draw.

    Attributes
    ----------
    sampler
        The underlying engine's own object after a run — ``emcee``'s
        ``EnsembleSampler``, ``dynesty``'s sampler — for anything this driver
        does not expose. ``None`` before the first run.
    last_failure_summary
        :meth:`~ampere.core.dataset.FittingProblem.failure_summary` as of the
        end of the last run; ``""`` when nothing failed.

    Notes
    -----
    **Prior support is the user's responsibility, and it bites.** The default
    catch set is narrow on purpose (``inference.md`` §11, confirmed by Peter
    2026-09-03): a :class:`~ampere.core.exceptions.LikelihoodError` and the
    exception types named in ``simulator_failures=`` become ``-inf`` with a
    recorded reason, and *everything else propagates*, because a composition
    bug turned into ``-inf`` is a fit that runs, converges and is wrong.

    The sharp edge that follows is worth stating before a user meets it at
    proposal 40 000. Models raise :class:`ValueError` for physically
    meaningless parameter values — the shipped
    :class:`~ampere.backends.reference.BlackBody` does so for a non-positive
    temperature — and a ``ValueError`` is *not* in the default catch set. So a
    full-real-support prior on a positivity-constrained parameter
    (``scipy.stats.norm`` on a temperature) will eventually propose a negative
    value and kill the run mid-flight. Two remedies, and the first is usually
    the right one:

    * give the parameter a prior whose support is the parameter's support —
      ``loguniform``, ``lognorm``, ``truncnorm``, ``halfnorm``. The sampler
      then never proposes an impossible value at all, and the fit is stating
      what the user actually believes;
    * or declare the exception: ``FittingProblem(..., simulator_failures=(
      ValueError,))``, which turns those proposals into recorded failures. This
      is the right answer for a wrapped external code whose refusals are not
      known in advance, and the wrong one for a prior that should have been
      bounded, because it hides the mis-declaration behind a failure count.

    **Failures are per-process** (``inference.md`` limitation 17.7). These
    drivers are single-process; running a sampler over a ``Pool`` would leave
    each worker with its own counts and the aggregate incomplete. None of them
    takes a pool for that reason.
    """

    #: What ``check_engine`` and the provenance attrs call this engine.
    NAME: ClassVar[str]

    #: Whether the engine can use a gradient. All three drivers here are
    #: gradient-free, and pass this to ``check_engine`` explicitly rather than
    #: letting it default to what the *problem* can supply: the question being
    #: asked is "can emcee run this?", not "is this problem differentiable?".
    OFFERS_GRADIENTS: ClassVar[bool] = False

    def __init__(
        self,
        problem: FittingProblem,
        *,
        cache_size: int = DEFAULT_CACHE_SIZE,
    ) -> None:
        if problem.free_size == 0:
            raise EngineError(
                f"{self.NAME} has nothing to sample: this problem declares no free parameters "
                f"(every one of {list(problem.parameters.names)} is fixed)."
            )
        # likelihoods.md §16's second obligation, discharged per dataset with
        # `observed=` -- so a censored sample the mask excludes is not counted
        # against a gradient-free engine (inference.md §10).
        problem.check_engine(self.NAME, differentiable=self.OFFERS_GRADIENTS)
        self.problem = problem
        self.sampler: Any = None
        self.last_failure_summary: str = ""
        self._cache = _EvaluationCache(problem, cache_size)

    # -- the §4.5 surface, and nothing else -----------------------------------

    @property
    def backend(self) -> str:
        """Which rung of the capability ladder supplies this problem's pieces.

        Read from the problem, never declared here (W2.12): the backend is
        §4.5's fourth capability flag, aggregated from what the models and
        transformations themselves say, so ``ampere_backend`` in the emitted
        provenance is a fact about the run rather than an assertion by whoever
        constructed the driver. A driver that disagreed with its problem used
        to be able to say so and be believed; now it cannot say anything.
        """
        return self.problem.backend

    def log_prob(self, theta: Any) -> float:
        """log p(θ) + log p(data | θ) — what an MCMC engine's callback wants."""
        return self._cache.evaluate(theta).log_prob

    def log_likelihood(self, theta: Any) -> float:
        """log p(data | θ) — what a nested sampler's callback wants.

        ``-inf`` rather than NaN for a point outside the prior's support.
        ``evaluate`` reports NaN there, deliberately ("not evaluated" is not
        "impossible", ``inference.md`` §10), but a nested sampler is fed θ from
        :meth:`prior_transform` and so never sees one; if it somehow does, NaN
        would poison its bookkeeping where ``-inf`` merely rejects the point.
        """
        value = self._cache.evaluate(theta).log_likelihood
        return value if not math.isnan(value) else -math.inf

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        """The unit hypercube to a free-parameter vector — nested sampling."""
        return self.problem.prior_transform(unit_cube)

    # -- shared run machinery -------------------------------------------------

    def stream(self, concern: str) -> np.random.Generator:
        """The generator for one named concern of this engine.

        ``inference.md`` §12: different concerns must not share a stream, or
        adding a diagnostic silently changes a fit's initialisation. Each
        driver's initialisation, its sampler's own randomness and (for nested
        sampling) its resampling are therefore separate labels under the
        engine's name.
        """
        return self.problem.rng(f"{self.NAME}.{concern}")

    def integer_seed(self, concern: str) -> int:
        """A 32-bit seed for a library that wants an integer, not a generator.

        Drawn from :meth:`stream`, so it is derived from the problem's seed
        along with everything else and is entropy-based exactly when the
        problem is (``seed=None`` means nothing is reproducible, which is the
        honest behaviour for a run that did not ask to be).
        """
        return int(self.stream(concern).integers(0, 2**32))

    def initial_positions(self, count: int, *, attempts: int = 200) -> np.ndarray:
        """*count* start points drawn from the joint prior, each inside the support.

        A draw from the prior is the right start for an ensemble: it is where
        the user says the mass is, it needs no tuning, and it is reproducible
        from the problem's seed. Draws that cannot be scored are re-drawn
        rather than kept — an ensemble move from a ``-inf`` walker cannot go
        anywhere — and a prior that cannot produce a scoreable point at all is
        a composition problem the run should stop for, not sample through.
        """
        if count < 1:
            raise EngineError(f"{self.NAME} needs at least one start point, got {count}.")
        positions = np.empty((count, self.problem.free_size), dtype=float)
        rng = self.stream("initialisation")
        for index in range(count):
            for _ in range(attempts):
                theta = self.problem.parameters.pack(self.problem.sample_prior(rng))
                if math.isfinite(self.log_prob(theta)):
                    positions[index] = theta
                    break
            else:
                raise EngineError(
                    f"{self.NAME} could not find a start point with a finite log-probability in "
                    f"{attempts} draws from the joint prior. Either the prior puts (almost) all "
                    f"its mass where the model or likelihood cannot be evaluated, or the data and "
                    f"the model disagree so completely that every prior draw underflows. "
                    f"problem.failure_summary() says which:\n{self.problem.failure_summary()}"
                )
        return positions

    def start(self) -> None:
        """Begin a run: forget the previous one's failures.

        ``failure_counts`` is what the provenance attrs record and what
        :meth:`finish` reports, so it has to mean "this run" rather than "every
        evaluation this problem has ever done". Clearing here rather than after
        emission means the start points count too, which is deliberate: a prior
        half of whose draws cannot be scored is exactly what the user needs to
        be told.
        """
        self.problem.reset_failures()

    def _evaluations_from_terms(
        self,
        draws: np.ndarray,
        terms: Sequence[Sequence[Mapping[str, float]]],
    ) -> list[list[Evaluation]]:
        """Per-draw evaluations from a decomposition the driver supplied.

        ``inference.md`` §10a's optional ``log_likelihood_terms``, consumed —
        W2.13 shipped it on both realisations and left the consumption to
        whichever slice 2 got there first (W2.4).

        What it buys is not a shortcut. The three gradient-free drivers score
        every proposal through :class:`_EvaluationCache`, so a stored draw is
        usually a lookup; a gradient-based driver scores through the *realised*
        density, which returns a scalar rather than an ``Evaluation``, so
        without this the per-dataset split had to be recomputed for every
        stored draw **on the numpy contract path** — a full model evaluation
        per draw, in the slowest available arithmetic, of a quantity the
        backend had just computed. The realisation's own decomposition is the
        same quantity in the backend's own array library, and it is the one the
        sampler actually used.

        The prior is still taken from the declaration, through
        :meth:`~ampere.core.parameter.ParameterSet.lnprior`. That is deliberate
        and it is cheap: a prior evaluation touches no model and no data, so
        there is nothing to save by asking a realisation for it — and §10a's
        mandatory surface does not include a prior/likelihood split, so a
        driver that demanded one would be requiring more of a backend than the
        contract does.

        Every label is checked against the problem's own datasets. A
        realisation whose keys had drifted would otherwise emit a
        ``log_likelihood`` group that silently omitted a dataset, which is a
        quieter failure than it should be.
        """
        expected = set(self.problem.datasets)
        built: list[list[Evaluation]] = []
        for chain, chain_terms in zip(draws, terms, strict=True):
            row: list[Evaluation] = []
            for theta, contributions in zip(chain, chain_terms, strict=True):
                labels = set(contributions)
                if labels != expected:
                    raise EngineError(
                        f"{self.NAME}'s realisation supplied a per-dataset decomposition keyed "
                        f"{sorted(labels)}, but this problem's datasets are "
                        f"{sorted(expected)}. The two must agree exactly: a missing label would "
                        f"be a dataset silently dropped from the run's log_likelihood group, and "
                        f"an extra one a group with no data behind it."
                    )
                log_prior = float(self.problem.parameters.lnprior(theta))
                log_likelihood = float(sum(contributions.values()))
                row.append(
                    Evaluation(
                        log_prior=log_prior,
                        log_likelihood=log_likelihood,
                        log_prob=log_prior + log_likelihood,
                        contributions=dict(contributions),
                    )
                )
            built.append(row)
        return built

    def finish(
        self,
        draws: Any,
        *,
        extra_attrs: Mapping[str, object] | None = None,
        coords: Mapping[str, Sequence[Any]] | None = None,
        realised: bool = False,
        registered_lowerings: Sequence[Mapping[str, Any]] | None = None,
        log_likelihood_terms: Sequence[Sequence[Mapping[str, float]]] | None = None,
    ) -> Any:
        """Assemble the stored draws into the run's ``DataTree``.

        Every run goes through here, which is how "every engine emits it"
        (``DEVELOPMENT_PLAN.md`` §4.6) is a fact about the code rather than a
        convention. The per-draw ``log_likelihood``/``log_prior`` split, the
        per-dataset decomposition, the observed data and the full provenance
        attrs are all :func:`ampere.results.emit`'s to write; what this method
        adds is the evaluations, the engine's own settings, and the failure
        summary.

        *realised* and *registered_lowerings* are the two a driver that
        sampled through a backend's realisation supplies (``inference.md``
        §10a's "Provenance", W2.13). The three gradient-free drivers leave
        them alone and their runs record ``ampere_realised = 0``, which is a
        fact worth having rather than an absence: it says the draws were
        scored on the numpy contract path.

        *log_likelihood_terms* is §10a's optional third, added at W2.4 slice 2:
        the per-dataset decomposition **the realisation computed**, shaped
        ``(chains, draws)`` and keyed by dataset label. Supplied, it replaces
        the evaluation cache entirely for this run and
        ``engine_draws_recomputed`` is zero, because nothing was recomputed —
        see :meth:`_evaluations_from_terms`. Omitted, the cache behaves exactly
        as it always has.

        The summary is surfaced three ways, because the three have different
        audiences: a :class:`~ampere.inference.exceptions.
        SamplingFailureWarning` for the person watching the run,
        ``ampere_failure_summary`` in the attrs for whoever reads the archived
        file, and :attr:`last_failure_summary` for a script that wants to
        branch on it. The counts and the bounded history are in the attrs
        already, by ``inference.md`` §18(b)'s request to W1.8.
        """
        array = np.asarray(draws, dtype=float)
        if array.ndim == 2:
            array = array[np.newaxis, ...]
        if log_likelihood_terms is None:
            evaluations = [[self._cache.lookup(theta) for theta in chain] for chain in array]
        else:
            evaluations = self._evaluations_from_terms(array, log_likelihood_terms)

        summary = self.problem.failure_summary()
        self.last_failure_summary = summary
        attrs: dict[str, object] = {
            "engine_evaluations": self._cache.calls,
            "engine_draws_recomputed": self._cache.recomputed,
        }
        attrs.update(extra_attrs or {})
        if summary:
            attrs["failure_summary"] = summary
        tree = emit(
            self.problem,
            array,
            evaluations,
            engine=self.NAME,
            coords=coords,
            realised=realised,
            registered_lowerings=registered_lowerings,
            extra_attrs=attrs,
        )
        if summary:
            warnings.warn(
                f"{self.NAME}: {summary}",
                SamplingFailureWarning,
                stacklevel=2,
            )
        return tree

    def __repr__(self) -> str:
        return f"<{type(self).__name__} on {self.problem!r}>"


def _check_ensemble(engine: str, walkers: int, free_size: int) -> int:
    """The conditions both ensemble samplers impose on their walker count.

    Shared because emcee and zeus impose the same ones, for the same reason: a
    stretch (or slice) move builds its proposal for one half of the ensemble
    out of the *other* half, so the ensemble must split into two halves at all
    and each half must span the parameter space. emcee warns and continues
    where the last of these fails; zeus raises. Refusing here means the two
    drivers behave the same way and the message names the fix.
    """
    if walkers < 2:
        raise EngineError(f"{engine} needs at least 2 walkers, got {walkers}.")
    if walkers % 2:
        raise EngineError(
            f"{engine} splits its ensemble in half each step, so the number of walkers must be "
            f"even; got {walkers}."
        )
    if walkers < 2 * free_size:
        raise EngineError(
            f"{engine} needs at least 2 x n_dim = {2 * free_size} walkers to span a "
            f"{free_size}-dimensional problem, got {walkers}. Fewer than that and each half of "
            f"the ensemble is confined to a subspace, which the sampler cannot leave."
        )
    return int(walkers)


def _default_walkers(free_size: int) -> int:
    """Four per dimension, at least eight — and even, since ``4 n`` always is.

    Comfortably above the ``2 x n_dim`` floor without making a small problem
    pay for a large ensemble; the usual advice, and cheap to override.
    """
    return max(8, 4 * free_size)


def _kept(total: int, burn_in: int, thin: int, engine: str) -> None:
    """Refuse a burn-in/thin combination that would keep nothing."""
    if total < 1:
        raise EngineError(f"{engine} needs at least one step, got {total}.")
    if thin < 1:
        raise EngineError(f"{engine}'s thinning factor must be at least 1, got {thin}.")
    if burn_in < 0:
        raise EngineError(f"{engine}'s burn-in cannot be negative, got {burn_in}.")
    if burn_in >= total:
        raise EngineError(
            f"{engine} was asked to discard {burn_in} of {total} step(s) as burn-in, which leaves "
            f"nothing to emit. A run has to keep at least one draw."
        )
    if (total - burn_in) // thin < 1:
        raise EngineError(
            f"{engine} was asked to thin {total - burn_in} kept step(s) by {thin}, which leaves "
            f"nothing to emit."
        )
