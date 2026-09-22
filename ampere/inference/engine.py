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

The gradient-free fast path (W2.5 slice 3)
------------------------------------------
The same optional member answers a second question, and this one is about
*scoring* rather than about bookkeeping. ``problem.evaluate`` is the numpy
contract path, and on a backend whose arithmetic is not numpy it pays that
backend's per-operation dispatch on every proposal: W2.10 measured the jax
contract path at about 28 ms flat per ``log_prob``, whatever the problem, which
turns emcee, dynesty and zeus on a jax problem into a pessimisation — the
gradient-free engines were slower on the backend built for speed.

So :class:`_EvaluationCache` grows a second route. When a realisation is
registered for the problem's backend **and** it supplies
``log_likelihood_terms``, the cache scores through it: the prior from
:meth:`~ampere.core.parameter.ParameterSet.lnprior` (which touches no model and
no data, so there is nothing to save by asking a backend for it, and §10a's
mandatory surface does not offer a prior/likelihood split anyway), the
per-dataset terms from the realisation, and the same
:class:`~ampere.core.dataset.Evaluation` record out. Measured on a 400-point
quasiseparable jax problem: 26.5 ms a proposal on the contract path, 0.5 ms on
this one.

Four things about it are deliberate.

* **Nothing here imports a backend.** ``ampere.core.realise`` dispatches on
  ``problem.backend``, so the rule this module opens with is untouched; a
  problem whose backend has no realisation registered simply keeps the contract
  path.
* **The numpy path stays the oracle and the fallback.** ``realise`` checks the
  realised density against ``problem.log_prob_unconstrained`` at the reference
  point before this cache will use it, the labels are checked once at
  attachment, and any exception from the realisation at *evaluation* time falls
  back to ``problem.evaluate`` for that θ rather than failing the run.
* **The failure reason is lost, and that is §10a's own ruling** (sub-decision
  2): inside a realised density nothing can raise and no ``Failure`` record can
  be built, so a point the model cannot score is a bare ``-inf``.
  ``problem.failure_counts`` therefore stays empty on this path and
  ``ampere_failure_summary`` says nothing. A problem declaring
  ``strict=True`` — the declaration that says "do not hide a failure from me" —
  keeps the contract path for exactly that reason, and ``use_realisation=False``
  restores it for anyone else who would rather have the reasons than the speed.
* **The run says which path it took**: ``engine_realised_evaluations`` counts
  the proposals scored through the realisation, beside ``engine_evaluations``,
  and ``ampere_realised`` is 1 for a run that used it — which is what that
  attribute has always meant ("whether the draws were scored through the
  backend's realisation"), now answered by a gradient-free engine too.
"""

from __future__ import annotations

import abc
import contextlib
import math
import random
import warnings
from collections.abc import Callable, Iterator, Mapping, Sequence
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import Evaluation, FittingProblem
from ampere.core.exceptions import LoweringError
from ampere.core.realisation import (
    log_likelihood_terms_of,
    realise,
    registered_realisations,
)
from ampere.results import emit

from .exceptions import EngineError, SamplingFailureWarning

__all__ = [
    "DEFAULT_CACHE_SIZE",
    "Engine",
    "default_live_points",
    "global_seed",
    "unconstrained_jacobian_correction",
]

#: How many scored θ an engine remembers so that the stored draws need not be
#: re-evaluated. An ``Evaluation`` is a handful of floats plus one term per
#: dataset, so 100 000 of them is tens of megabytes — the same order as the
#: emitted run itself, and small beside the model evaluations it saves. Past
#: the bound the oldest go and the affected draws are recomputed; nothing is
#: lost but time.
DEFAULT_CACHE_SIZE = 100_000


def unconstrained_jacobian_correction(problem: FittingProblem, unconstrained: Any) -> np.ndarray:
    """The change-of-variables term alone, per row of an unconstrained draw.

    ``ParameterSet.lnprior_unconstrained(y) = lnprior(constrain(y)) +
    Σ log|d constrain / dy|`` is the reference-path oracle every realisation
    already agrees with (``ampere.core.realisation``); this is that sum term
    on its own, recovered as the difference of two numpy-contract-path calls
    rather than restated, so there is exactly one place that knows the
    formula.

    What it is for (``results.md`` §9, W5.0): a proposal fitted in
    **unconstrained** coordinates — a VI guide, an SBI density estimator —
    has its own log-density there, and that density is *not* what
    ``sample_stats.proposal_log_density`` must hold, because the stored
    ``log_prior``/``log_likelihood`` are in the **constrained** free-parameter
    space (``Engine._evaluations_from_terms`` scores ``lnprior`` at the
    constrained θ, as every gradient-free driver's cache does too). Subtracting
    this term from the proposal's unconstrained-space density moves it into
    the same coordinates, which is what makes
    ``exp(log_prior + log_likelihood - proposal_log_density)`` a valid
    importance weight computed from the stored groups alone — the two
    Jacobian terms cancel by construction, so a driver that skipped this step
    would silently bias every reweighting by exactly the change of variables
    it forgot.

    Backend-neutral: :meth:`~ampere.core.parameter.ParameterSet.constrain`,
    ``.lnprior`` and ``.lnprior_unconstrained`` are all numpy-contract-path
    calls on ``ampere.core``, so a driver may call this without importing a
    backend, wherever in its own routine it happens to hold the unconstrained
    draws.
    """
    parameters = problem.parameters
    array = np.atleast_2d(np.asarray(unconstrained, dtype=float))
    out = np.empty(array.shape[0], dtype=float)
    for row, y in enumerate(array):
        constrained = parameters.constrain(y)
        unconstrained_prior = float(parameters.lnprior_unconstrained(y))
        constrained_prior = float(parameters.lnprior(constrained))
        out[row] = unconstrained_prior - constrained_prior
    return out


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

    def __init__(
        self,
        problem: FittingProblem,
        maximum: int,
        *,
        use_realisation: bool = True,
    ) -> None:
        if maximum < 0:
            raise EngineError(f"the evaluation cache cannot hold {maximum} entries.")
        self._problem = problem
        self._maximum = int(maximum)
        self._entries: dict[bytes, Evaluation] = {}
        #: How many times the problem was actually scored, by either route.
        self.calls = 0
        #: How many stored draws had to be evaluated a second time.
        self.recomputed = 0
        #: How many of :attr:`calls` went through the backend's realisation
        #: rather than the numpy contract path. ``engine_realised_evaluations``.
        self.realised_calls = 0
        #: The realisation this cache scores through, and its per-dataset
        #: decomposition; both ``None`` on the contract path.
        self.realisation: Any = None
        self._terms: Callable[[Any], Mapping[str, Any]] | None = None
        if use_realisation:
            self._attach(problem)

    # -- the fast path -------------------------------------------------------

    def _attach(self, problem: FittingProblem) -> None:
        """Take the backend's realisation, if this problem has a usable one.

        Every "no" here is a silent fall back to the contract path, because
        every one of them is a legitimate configuration rather than a mistake:
        the reference backend registers no realisation, a jax problem using a
        family this backend does not lower cannot be realised at all, and a
        realisation is free by ruling not to offer ``log_likelihood_terms``.

        ``strict=True`` is the one deliberate refusal. It is the declaration
        that says "do not turn a failure into a number" — and a realised
        density can do nothing else (``inference.md`` §10a sub-decision 2:
        inside a trace nothing raises and no ``Failure`` can be built), so
        honouring the flag means declining the speed.
        """
        if getattr(problem, "strict", False):
            return
        if problem.backend not in registered_realisations():
            return
        try:
            realisation = realise(problem)
        except LoweringError:
            return
        terms = log_likelihood_terms_of(realisation)
        if terms is None:
            return
        # One evaluation at the reference point, which checks the labels and
        # also pays whatever compilation the backend does before the sampler
        # starts rather than on its first proposal.
        try:
            supplied = set(terms(problem.unconstrain(problem.reference_values)))
        except Exception:
            return
        expected = set(problem.datasets.contribution_labels())
        if supplied != expected:
            raise EngineError(
                f"the realisation registered for backend {problem.backend!r} supplied a "
                f"per-dataset decomposition keyed {sorted(supplied)}, but this problem's "
                f"contribution labels are {sorted(expected)}. The two must agree exactly: a "
                f"missing label would be a dataset silently dropped from the run's "
                f"log_likelihood group, and an extra one a group with no data behind it."
            )
        self.realisation = realisation
        self._terms = terms

    @property
    def realised(self) -> bool:
        """Whether this cache scores through a realisation (``ampere_realised``)."""
        return self._terms is not None

    def _score(self, vector: np.ndarray) -> Evaluation:
        """One θ, through whichever route this cache has."""
        self.calls += 1
        if self._terms is None:
            return self._problem.evaluate(vector)
        log_prior = float(self._problem.parameters.lnprior(vector))
        if not math.isfinite(log_prior):
            # Zero prior mass is an answer, not a failure, and no model runs
            # for it on either path. NaN rather than -inf for the likelihood,
            # because "not evaluated" is not "impossible" — the contract path
            # says exactly this, and the two must agree.
            return Evaluation(log_prior=-math.inf, log_likelihood=math.nan, log_prob=-math.inf)
        try:
            supplied = self._terms(self._problem.unconstrain(vector))
            contributions = {str(label): float(value) for label, value in supplied.items()}
        except Exception:
            # The oracle is also the fallback: a realisation that cannot score
            # this point hands it back to numpy rather than failing the run.
            return self._problem.evaluate(vector)
        self.realised_calls += 1
        total = math.fsum(contributions.values())
        if not math.isfinite(total):
            # A bare -inf, with no reason: §10a sub-decision 2 says the reason
            # is unrecoverable on this path, and inventing one would be worse
            # than admitting it. `use_realisation=False` is where the reasons are.
            return Evaluation(
                log_prior=log_prior,
                log_likelihood=-math.inf,
                log_prob=-math.inf,
                contributions=contributions,
            )
        log_prob = log_prior + total
        return Evaluation(
            log_prior=log_prior,
            log_likelihood=total,
            # The contract path's own guard, kept identical: two finite halves
            # can still overflow to -inf, and an Evaluation must never carry a
            # log_prob its two halves do not imply.
            log_prob=log_prob if math.isfinite(log_prob) else -math.inf,
            contributions=contributions,
        )

    @staticmethod
    def _key(theta: Any) -> tuple[np.ndarray, bytes]:
        vector = np.asarray(theta, dtype=float).reshape(-1)
        return vector, vector.tobytes()

    def evaluate(self, theta: Any) -> Evaluation:
        """Score θ, remembering the result. What the sampler's callback calls."""
        vector, key = self._key(theta)
        evaluation = self._score(vector)
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
        return self._score(vector)


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
    use_realisation
        Whether to score proposals through the backend's **realisation** when
        one is registered for this problem's backend and supplies
        ``log_likelihood_terms`` (W2.5 slice 3; this module's docstring has the
        shape and the measurement). ``True`` by default, because on a
        differentiable backend the numpy contract path costs that backend's
        per-operation dispatch on every proposal and there is nothing to be
        gained by paying it. It changes no number: ``ampere.core.realise``
        checks the realised density against the contract path before this
        driver will use it, and the conformance suite compares the two at many
        points. It changes one *record*: a proposal the model cannot score
        becomes ``-inf`` with **no recorded reason**, so
        ``problem.failure_counts`` and ``ampere_failure_summary`` stay empty
        (``inference.md`` §10a sub-decision 2). Pass ``False`` to keep the
        reasons; a problem declaring ``strict=True`` keeps them anyway, since
        that flag means precisely "do not hide a failure from me".

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
        use_realisation: bool = True,
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
        self._cache = _EvaluationCache(problem, cache_size, use_realisation=use_realisation)

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

        Every label is checked against the problem's own
        :meth:`~ampere.core.dataset.DatasetCollection.contribution_labels` — a
        dataset's own label, or, for a dataset a **joint noise group** claims,
        the group's (W5.9; ``inference.md`` §4.9's ``"joint"`` decomposition).
        A realisation whose keys had drifted would otherwise emit a
        ``log_likelihood`` group that silently omitted a dataset, which is a
        quieter failure than it should be.
        """
        expected = set(self.problem.datasets.contribution_labels())
        built: list[list[Evaluation]] = []
        for chain, chain_terms in zip(draws, terms, strict=True):
            row: list[Evaluation] = []
            for theta, contributions in zip(chain, chain_terms, strict=True):
                labels = set(contributions)
                if labels != expected:
                    raise EngineError(
                        f"{self.NAME}'s realisation supplied a per-dataset decomposition keyed "
                        f"{sorted(labels)}, but this problem's contribution labels are "
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
        realised: bool | None = None,
        registered_lowerings: Sequence[Mapping[str, Any]] | None = None,
        log_likelihood_terms: Sequence[Sequence[Mapping[str, float]]] | None = None,
        sample_stats: Mapping[str, Any] | None = None,
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
        §10a's "Provenance", W2.13). Left alone (``None``), *realised* is
        answered by the evaluation cache: it is ``True`` exactly when the
        gradient-free fast path scored this run's proposals through the
        realisation, and ``False`` when they went through the numpy contract
        path. Until W2.5 slice 3 the three gradient-free drivers could only
        record ``0``; now the attribute says which of the two actually
        happened, which is what ``ampere_realised`` has always meant.

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

        *sample_stats* is ``ampere.results.emit``'s hook for a per-draw
        quantity beside ``lp``/``log_prior``/``log_likelihood`` — passed
        straight through, shaped ``(chains, draws)`` (``results.md`` §4/§9,
        W5.0). ``proposal_log_density`` is the one every approximate driver
        supplies here.

        ``ampere_approximation`` defaults to ``"none"`` — this base class's
        five gradient-free and gradient-based samplers are exact — and a
        driver whose draws are not from the target overrides it through
        *extra_attrs* (``"mean_field"``/``"multivariate"`` for
        :class:`~ampere.inference.VIEngine`, ``"density_estimator"`` for
        :class:`~ampere.inference.SBIEngine`). Writing the default here rather
        than on each driver is the point: the key ``results.md`` §9 asks every
        run to carry is carried by construction, not by five drivers
        remembering to say "none".
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
            "engine_realised_evaluations": self._cache.realised_calls,
            "approximation": "none",
        }
        attrs.update(extra_attrs or {})
        if summary:
            attrs["failure_summary"] = summary
        if registered_lowerings is None and self._cache.realisation is not None:
            provenance = getattr(self._cache.realisation, "lowering_provenance", None)
            registered_lowerings = provenance() if callable(provenance) else None
        tree = emit(
            self.problem,
            array,
            evaluations,
            engine=self.NAME,
            coords=coords,
            realised=self._cache.realised if realised is None else realised,
            registered_lowerings=registered_lowerings,
            extra_attrs=attrs,
            sample_stats=sample_stats,
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


def _refuse_foreign_parts(engine: str, problem: FittingProblem) -> None:
    """Refuse *problem* for a gradient-based *engine* if it carries foreign parts.

    **W3.8** (ruled by Peter 2026-09-08). ``allow_foreign_parts=True`` already
    makes the problem non-differentiable, so the check below this one would
    refuse it anyway — but it would refuse it as "one of your pieces declares
    ``DIFFERENTIABLE = False``", which sends a user looking at their model when
    the answer is a numpy kernel three levels down. Naming the pieces is the
    whole difference between a refusal a user can act on and one they cannot.

    Shared by :class:`~ampere.inference.NUTSEngine` and
    :class:`~ampere.inference.VIEngine` so the two say the same thing, and
    written here rather than in ``ampere.core`` because it raises this
    package's :class:`~ampere.inference.exceptions.EngineError`.
    """
    if not problem.foreign_parts:
        return
    named = ", ".join(problem.foreign_part_names)
    raise EngineError(
        f"{engine} needs a gradient, and this problem composes {len(problem.foreign_parts)} "
        f"piece(s) from another backend — {named} — whose arrays would have to be converted to "
        f"{problem.backend!r} at every evaluation, which detaches the graph and leaves their "
        f"parameters with no gradient at all. FittingProblem(..., allow_foreign_parts=True) "
        f"declares that a *gradient-free* run may call through a piece that exists only in "
        f"Python; it cannot make one differentiable. Build those pieces from "
        f"{problem.backend}'s own classes, or sample with emcee, dynesty or zeus."
    )


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


@contextlib.contextmanager
def global_seed(seed: int | None) -> Iterator[None]:
    """Seed the two *process-global* generators, then put both back.

    Some sampling libraries take no generator and expose no ``random_state``:
    they draw from process-global state. The only way to make such a run
    reproducible from ampere's own seed is to seed that state around the run —
    and the only way to do it without a side effect on the caller's streams is
    to save and restore, which is what this does.

    **Both** generators, and the second one is the whole reason this is not one
    line. zeus's sampling loop draws from numpy's legacy global
    (``np.random.uniform``/``exponential``/``shuffle``/``choice`` throughout
    ``zeus/ensemble.py``), *and* its default ``DifferentialMove.get_direction``
    picks its walker pairs with the standard library's ``random.sample``
    (``zeus/moves.py``). Seeding numpy alone leaves the pair selection
    entropy-seeded, and a run reproducible in every draw except which walkers
    proposed for which is not reproducible at all — it just looks like it might
    be until someone checks. ultranest (W5.14) needs the numpy half for the
    same reason: its region sampling, its step samplers and its bootstrap all
    call ``np.random.*`` directly.

    With ``seed=None`` both globals are left completely alone, which is the
    honest behaviour: a problem built without a seed asked not to be
    reproducible, and seeding-then-restoring would make its consecutive runs
    identical instead.

    Recorded as a limitation rather than hidden: this is *global* state, so a
    run under this context manager is not thread-safe against other code
    drawing from ``np.random`` or ``random`` at the same time. emcee, dynesty
    and nautilus have per-sampler streams and need none of it.
    """
    if seed is None:
        yield
        return
    numpy_state = np.random.get_state()
    python_state = random.getstate()
    try:
        np.random.seed(seed)
        random.seed(seed)
        yield
    finally:
        np.random.set_state(numpy_state)
        random.setstate(python_state)


def default_live_points(free_size: int) -> int:
    """25 per dimension plus a floor of 100 — every nested sampler's default here.

    Below roughly ``25 (n_dim + 1)`` the ellipsoidal (or neural, or
    region-based) bound is fitted from too few points to be trustworthy, which
    is where nested sampling starts to *under-cover* rather than merely run
    slowly; the floor keeps a one- or two-dimensional problem from being
    sampled by a handful of points.

    Written here rather than in one driver because W5.14 gave ampere three
    nested samplers (dynesty, nautilus, ultranest) and the reasoning is the
    same for all three. Naming them the same default is also what makes the
    engine battery's cross-engine evidence comparison a comparison of the
    *samplers* rather than of three differently-sized live sets.
    """
    return max(100, 25 * (free_size + 1))


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
