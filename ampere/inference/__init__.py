"""Engines: written once against §4.5's surface, run on every backend.

``DEVELOPMENT_PLAN.md`` §3's architectural bet, cashed. A
:class:`~ampere.core.dataset.FittingProblem` exposes ``log_prob``, the
``log_likelihood``/``log_prior`` split, ``prior_transform``, ``simulate`` and
the capability flags, and ``inference.md`` §10 claims that an engine consuming
only those "works with the reference backend, with torch, with jax, and with a
legacy black-box model behind a thin adapter, and never knows which". The
gradient-free drivers here are the first test of that claim, and they are
written to make it checkable rather than merely asserted: **this namespace does
not import ``ampere.backends``**, in any module, at any depth. Its only ampere
imports are ``ampere.core`` and ``ampere.results``.

Both Phase 2 backend tracks found the claim's one limit, and
:class:`NUTSEngine` below is where it shows: §4.5's surface is not
*traceable*, so a gradient-based engine cannot be written against it alone.
``inference.md`` §10a (W2.13) is the answer — a **realisation**, the backend's
registered differentiable native form of the problem — and the rule above is
unchanged by it: that driver still imports no backend, and reaches the density
through :func:`ampere.core.realise`, which dispatches on the problem's own
``backend`` flag.

What is here
------------
:class:`EmceeEngine`
    Goodman & Weare's affine-invariant ensemble sampler. The general-purpose
    default, and a base dependency.
:class:`DynestyEngine`
    Nested sampling, consuming ``prior_transform`` rather than ``log_prior``.
    Multimodal posteriors, and the only one of the three that yields a marginal
    likelihood. Also a base dependency.
:class:`ZeusEngine`
    Ensemble slice sampling. Needs the ``zeus`` extra -- ``pip install
    ".[zeus]"`` from a checkout (PyPI's ``ampere`` package is unrelated); the
    import is lazy and the refusal names the extra.
:class:`NautilusEngine`, :class:`UltranestEngine`
    Two more **nested** samplers, added at W5.14 from the inference-extensions
    memo's tier 1, sharing one driver
    (:mod:`ampere.inference._nested`) with each other and the evidence
    convention with :class:`DynestyEngine`: the engine-neutral triple
    (``ampere_log_evidence``, ``_err``, ``_method``), the equal-weight
    resampling rule, and the live-point default. nautilus is importance nested
    sampling with a neural boundary -- fewer likelihood calls to a given
    effective sample size, which is the budget that matters for an expensive
    forward model -- and refuses a one-dimensional problem. ultranest is
    MLFriends region sampling with a bootstrapped termination criterion and an
    evidence error whose two halves it keeps apart. **One extra each**
    (``nautilus``, ``ultranest``), imported lazily in the constructor, and the
    refusal names the extra.
:class:`BlackjaxEngine`
    The blackjax route, **jax only**, added at W5.14 from the inference-
    extensions memo's §2.1: one dependency, two methods ampere does not
    otherwise have. ``method="mclmc"`` is microcanonical Langevin Monte Carlo
    (Robnik & Seljak) -- a Markov chain with no accept/reject step, tuned for
    you, and reported at several times NUTS's efficiency per gradient on
    smooth high-dimensional posteriors. ``method="pathfinder"`` is Zhang et
    al.'s quasi-Newton approximation, used twice over: as an approximate
    posterior in its own right (writing ``ampere_approximation =
    "pathfinder"`` and the per-draw ``proposal_log_density`` that makes it
    correctable) and, through ``initial="pathfinder"``, as the start-point
    search for an MCLMC chain. Neither estimates an evidence and neither
    pretends to. **One extra** (``blackjax``), imported lazily in the
    constructor, and the refusal names the extra.
:class:`NUTSEngine`
    The No-U-Turn sampler — numpyro's on a jax problem, pyro's on a torch one.
    The first **gradient-based** engine here, and the only one that cannot be
    written against §4.5's surface alone: ``FittingProblem.log_prob_unconstrained``
    is not traceable (the containers coerce with ``numpy.asarray``, and
    ``lnprior`` short-circuits on ``math.isfinite``), so the density comes
    from the backend's **realisation** (``inference.md`` §10a) — obtained
    through :func:`ampere.core.realise`, which importing the backend
    registered. ``NUTSEngine(problem)`` takes no density argument; the driver
    imports no backend, dispatches on the problem's own ``backend`` flag
    between the two samplers, and imports the one it needs lazily inside
    ``run``. Needs the ``jax`` or ``torch`` extra.
:class:`VIEngine`
    Stochastic variational inference — pyro's SVI on a torch problem,
    numpyro's on a jax one, over the same realisation :class:`NUTSEngine`
    samples and dispatched on the problem's own ``backend`` flag in the same
    way. **Approximate**, and the approximation is the guide family the run
    records in ``vi_guide``: ``"normal"`` assumes the posterior factorises
    over parameters, ``"multivariate"`` captures their correlations, and
    neither captures a non-Gaussian posterior. The honest use is a first look
    at a posterior, or the only tractable route when the space is too large
    for MCMC. A VI run's single "chain" of i.i.d. guide draws is deliberate:
    guide draws have no Markov structure, so R-hat has nothing to say about
    them, and the ELBO trace is what a reader looks at instead. Needs the
    ``torch`` or ``jax`` extra.
:class:`SBIEngine`
    **Simulation-based** inference — neural posterior, likelihood or ratio
    estimation through the ``sbi`` package (``method="npe" | "nle" | "nre" |
    "tmnre"``).
    The one driver here that does not consume ``log_prob`` while it fits:
    §4.5's ``simulate(params) -> data``, batched through
    :meth:`~ampere.core.dataset.FittingProblem.simulate_many`, is the whole of
    its training surface, which is why it is the engine for a model whose
    likelihood cannot be written down — a wrapped external simulator composed
    as a black-box model on the reference backend, the case this field
    actually brings. ``executor=`` and ``chunk_size=`` pass straight through,
    so that budget runs under a process pool without this driver knowing.
    The prior is bridged into **unconstrained** coordinates, so bounded priors
    need no ``RestrictedPrior`` and the density estimator sees ℝⁿ; the drawn
    posterior is mapped back through ``constrain`` before anything is stored.
    Every stored draw is then scored on the numpy contract path, so the run
    carries the *true* per-draw split beside the estimator's own log-density
    (``ampere_sbi_log_prob`` in ``sample_stats``) — which is what makes
    calibration and importance reweighting possible later. One chain of
    i.i.d. draws, as :class:`VIEngine` emits and for the same reason.
    ``method="tmnre"`` is truncated marginal ratio estimation (Miller et al.
    2021) over the same trainer: one ratio estimator per 1-D marginal (and per
    pair at ``marginals=2``), a *prior* truncated between rounds to the box
    where those marginals put their mass, and a joint estimator trained last
    on the truncated prior so the run still emits ordinary joint draws. It
    carries a ``marginals`` group beside them and its truncation history in
    the attrs, and it is **not amortised** — the box is chosen at this
    observation. Needs the ``sbi`` extra; both it and torch are imported
    inside ``run``.

The first three are gradient-free, so all three call
:meth:`~ampere.core.dataset.FittingProblem.check_engine` with
``differentiable=False`` at construction — asking "can this engine run this
likelihood?" rather than "is this problem differentiable?" — and a
marginalisation no gradient-free engine can deliver is refused before any
sampling starts. :class:`NUTSEngine` passes ``differentiable=True`` for the
same reason and in the same spirit: it is stating what the *engine* offers,
and :class:`VIEngine` likewise. :class:`SBIEngine` passes ``False``, which is
the honest answer for an engine whose forward model may be a compiled binary.

Every run emits the run
-----------------------
There is no way to sample through these drivers and not get a stored run:
``run()`` returns the ArviZ ``DataTree``, built by :func:`ampere.results.emit`,
carrying per-draw ``log_prior`` and ``log_likelihood``, the per-dataset
log-likelihood decomposition, the observed data, and the full provenance attrs
— spec and data hashes, the backend-neutral model identity, seeds, package
versions, failure counts (``DEVELOPMENT_PLAN.md`` §4.6). Write it with
:func:`ampere.results.to_netcdf` and read it back with
:func:`ampere.results.from_netcdf`; the hashes survive the round trip by value,
which is what makes an archived run identifiable years later.

arviz and ``h5netcdf`` are base dependencies as of this item, for exactly that
reason (``results.md`` §15 R1, ruled 2026-09-03: promoted "with Phase 2's
engine drivers […] when a user can first emit a run").

Failures, and the prior-support responsibility
----------------------------------------------
The non-strict path is consumed as ``inference.md`` §11 declares it: a
proposal the model or likelihood cannot score becomes ``-inf`` with a recorded
reason, never an exception, and the aggregate is surfaced at the end of the run
three ways — a :class:`~ampere.inference.exceptions.SamplingFailureWarning`,
the ``ampere_failure_summary`` attr, and
:attr:`~ampere.inference.engine.Engine.last_failure_summary`. The intended
workflow is both halves together: run non-strict, read the summary, then re-run
with ``FittingProblem(..., strict=True)`` to get the raise at the offending
draw with a full traceback.

**The catch set is narrow, and the sharp edge is the user's priors.** Only
:class:`~ampere.core.exceptions.LikelihoodError` and the types named in
``simulator_failures=`` are caught; everything else propagates, because a
composition bug turned into ``-inf`` is a fit that runs, converges and is
wrong. A model that refuses a physically meaningless value does so with a
:class:`ValueError` — the shipped reference ``BlackBody`` does, for a
non-positive temperature — and a ``ValueError`` is *not* caught. So a
``scipy.stats.norm`` prior on a positivity-constrained parameter will kill a
run mid-flight the first time the sampler proposes a negative value. The remedy
is nearly always to give the parameter a prior whose support **is** the
parameter's support (``loguniform``, ``lognorm``, ``truncnorm``,
``halfnorm``); ``simulator_failures=(ValueError,)`` is the escape hatch, and is
the right answer for a wrapped external code rather than for a prior that
should have been bounded. :class:`~ampere.inference.engine.Engine`'s docstring
says the same at length, because it is the trap a new user meets first.

Seeds
-----
``inference.md`` §12's policy, used rather than re-derived: every stream comes
from :meth:`~ampere.core.dataset.FittingProblem.rng`, under a label naming the
engine and the concern (``"emcee.initialisation"``, ``"dynesty.resample"``), so
initialisation, the sampler's own randomness and any resampling never make each
other irreproducible. With ``problem.seed`` set, a run repeats exactly; with
``seed=None`` nothing is reproducible, which is the honest behaviour for a run
that did not ask to be.

One caveat, and it is zeus's rather than ampere's: zeus takes no generator and
draws from *two* process-global streams — numpy's legacy global for its slice
sampling, and the standard library's ``random`` for the walker pairs its
default move builds its directions from. :class:`ZeusEngine` seeds and restores
both around the run (:func:`ampere.inference.engine.global_seed`); seeding only
the first, which is what an inspection of ``zeus/ensemble.py`` alone suggests,
leaves the run irreproducible in a way that is easy to miss. Being global
state, a zeus run is not thread-safe against other code drawing from either
stream at the same time.

Not here, deliberately
----------------------
No multiprocessing pool **for the samplers**. ``inference.md`` limitation 17.7
records that the failure history is per-process, so a pooled run would leave
each worker with its own counts and the aggregate silently incomplete;
aggregating them is a real piece of work, not a keyword argument.
:class:`SBIEngine` is the one exception, and it is one because it does not
drive its own pool: its ``executor=`` goes to
:meth:`~ampere.core.dataset.FittingProblem.simulate_many`, which owns both
ends of the pool and replays every draw's failure record on the parent in draw
order (``inference.md`` §13, *Execution*). No optimisers —
``architecture.md`` §3 puts them in this namespace too, and they are a later
phase.

Examples
--------
A complete fit, on a model written by hand rather than taken from a backend —
which is the portability claim, exercised:

>>> import numpy as np, scipy.stats as st, astropy.units as u
>>> from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
>>> class Line(Model):
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.um)
...         self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(
...             ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy
...         )
>>> grid = np.array([1.0, 2.0, 3.0])
>>> observed = Spectrum(
...     grid * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.3, 0.3, 0.3] * u.Jy
... )
>>> problem = FittingProblem(Line(grid), [Dataset(observed)], seed=20260905)
>>> run = EmceeEngine(problem, walkers=8).run(steps=300, burn_in=100)

The run is the ArviZ format, with the per-draw split every run stores:

>>> sorted(run.children)
['constant_data', 'log_likelihood', 'observed_data', 'posterior', 'sample_stats']
>>> run["posterior"]["model.slope"].shape
(8, 200)
>>> sorted(run["sample_stats"].dataset.data_vars)
['failed', 'failure_reason', 'failure_where', 'log_likelihood', 'log_prior', 'lp']

and it knows what it was a run of. ``ampere_backend`` is *derived*, not
declared: since W2.12 the backend is §4.5's fourth capability flag, so this
hand-written model's own ``Model.BACKEND`` default is what put ``'reference'``
there, and no driver takes a ``backend=`` argument to say otherwise.

>>> run.attrs["ampere_engine"], run.attrs["ampere_backend"]
('emcee', 'reference')
>>> len(run.attrs["ampere_spec_hash"])
32

The posterior recovers the truth the data were generated at:

>>> bool(abs(float(run["posterior"]["model.slope"].mean()) - 2.0) < 0.1)
True
"""

from __future__ import annotations

from ._blackjax import METHODS as BLACKJAX_METHODS
from ._blackjax import BlackjaxEngine
from ._dynesty import DynestyEngine
from ._emcee import EmceeEngine
from ._nested import NESTED_ENGINES, NautilusEngine, UltranestEngine
from ._nuts import NUTSEngine
from ._sbi import (
    DEFAULT_TRUNCATION_EPSILON,
    EMBEDDINGS,
    LAYOUTS,
    MARGINAL_ORDERS,
    METHODS,
    SET_EMBEDDINGS,
    SUMMARY_LAYOUT,
    TMNRE_SAMPLERS,
    SBIEngine,
)
from ._vi import VIEngine
from ._zeus import ZeusEngine
from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError, SamplingFailureWarning

__all__ = [
    "BLACKJAX_METHODS",
    "DEFAULT_CACHE_SIZE",
    "DEFAULT_TRUNCATION_EPSILON",
    "EMBEDDINGS",
    "LAYOUTS",
    "MARGINAL_ORDERS",
    "METHODS",
    "NESTED_ENGINES",
    "SET_EMBEDDINGS",
    "SUMMARY_LAYOUT",
    "TMNRE_SAMPLERS",
    "BlackjaxEngine",
    "DynestyEngine",
    "EmceeEngine",
    "Engine",
    "EngineError",
    "NUTSEngine",
    "NautilusEngine",
    "SBIEngine",
    "SamplingFailureWarning",
    "UltranestEngine",
    "VIEngine",
    "ZeusEngine",
]
