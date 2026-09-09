"""The simulation-based driver: NPE, NLE, NRE and TMNRE over ``simulate_many``.

Private module; the class is :class:`ampere.inference.SBIEngine`. Named with a
leading underscore for the reason ``_emcee.py``, ``_nuts.py`` and ``_vi.py``
record for theirs.

Why this engine is different from the other five
------------------------------------------------
Every other driver here consumes ``log_prob``. This one does not consume it at
all while it is *fitting*: ``DEVELOPMENT_PLAN.md`` §4.5's ``simulate(params) ->
data`` is the whole of the surface it trains on, which is precisely why it is
the engine for a model whose likelihood cannot be written down — a compiled
radiative-transfer code behind a Python call, the canonical case this field
brings (``DEVELOPMENT_PLAN.md`` §2, *Batched simulation*, note (2)). A
black-box :class:`~ampere.core.transform.Model` composed on the **reference**
backend is therefore the first-class case here rather than an afterthought, and
``executor=`` and ``chunk_size=`` pass straight through to
:meth:`~ampere.core.dataset.FittingProblem.simulate_many` so that such a model
runs under a process pool without this driver knowing anything about it.

The rule the namespace opens with is untouched. ``sbi`` and torch are imported
**lazily, inside** :meth:`SBIEngine.run`, exactly as ``_nuts.py`` imports pyro,
so ``import ampere.inference`` in the base install pulls in neither; the refusal
is an :class:`~ampere.core.exceptions.OptionalDependencyError` naming the
``sbi`` extra (``architecture.md`` §4 rules 2 and 3).

The prior lives in unconstrained space, and so does the estimator
-----------------------------------------------------------------
``sbi`` wants a torch ``Distribution`` over the flat free vector. The obvious
construction — one distribution per parameter over its *constrained* support —
is the one this driver deliberately does not use. A density estimator trained
on bounded coordinates puts mass outside the bounds, which ``sbi`` then has to
correct for with a ``RestrictedPrior`` and a rejection loop whose acceptance
rate is a property of how badly the network overshot; and a normalising flow
asked to model a distribution on a half-line or an interval is being asked to
approximate a density with a hard edge, which is exactly what flows are worst
at.

So the bridge is :class:`_UnconstrainedPrior` over **ℝⁿ**: θ is drawn in the
constrained space the user declared (``FittingProblem.sample_prior``, which
respects ties and dependency order) and mapped through ``unconstrain``, and the
density is ``ParameterSet.lnprior_unconstrained`` — ``lnprior`` plus the
change-of-variables term, the same quantity ``inference.md`` §10 states on the
reference path so that every backend has an oracle. The two agree by
construction, and the suite checks that they do rather than trusting it: a
``sample`` that did not match its own ``log_prob`` would train the estimator on
one distribution and score it under another, which is the kind of error that
produces a plausible posterior that is wrong.

Three consequences, all of them recorded in the run's attrs
(``ampere_sbi_parameterisation``):

* bounded priors need no ``RestrictedPrior`` — the support is all of ℝⁿ, so
  nothing the estimator proposes can fall outside it;
* an NPE posterior has no leakage, so its ``log_prob`` needs no normalisation
  correction and the number the run stores is a normalised log-density;
* the drawn posterior is in unconstrained coordinates and is mapped back
  through ``constrain`` before anything is stored or scored, so the emitted
  ``posterior`` group is in the user's own coordinates like every other run's.

One warning is expected and is not a problem: ``sbi`` looks for ``mean`` and
``stddev`` on the prior to build the network's input standardisation, and a
distribution defined by ``sample`` and ``log_prob`` has neither in closed
form, so it says so and estimates them from draws. That is exactly what this
driver would have done, done by the library that needs the numbers, and it
affects only the affine z-scoring of the network's inputs — never the density.

What a run stores that a sampler's run does not
------------------------------------------------
Every stored draw is scored on the numpy contract path through
:meth:`~ampere.inference.engine.Engine.finish`, so ``lp``, ``log_prior``,
``log_likelihood`` and the per-dataset decomposition are the **true** ones —
computed by ``problem.evaluate``, not by the network. Beside them,
``sample_stats`` carries ``ampere_sbi_log_prob``: the *estimator's* own
log-density at the same draw. Having both is what makes
simulation-based calibration and importance reweighting possible later (W3.6,
design horizon (b)) — the ratio of the two is the importance weight — and it is
cheap to store now and impossible to recover afterwards.

That it is a per-draw **variable** rather than a provenance attribute is this
item's one small departure from the item text's parenthesis, and the reason is
alignment: an importance weight needs one estimator log-density per stored
draw, in draw order, so the values cannot be thinned the way the training-loss
trace is; a variable in ``sample_stats`` is where a ``(chain, draw)``-shaped
quantity belongs, survives the netCDF round trip as an array, and is what
``results.md`` §4 already shapes that group for. ``ampere_sbi_log_prob_kind``
says whether it is normalised — it is for NPE, where the posterior is a direct
density over ℝⁿ; it is not for NLE and NRE, whose posteriors are known only up
to the evidence.

The tensor the network sees, and the layout that fixes it
---------------------------------------------------------
The network sees **one tensor**, not a
:class:`~ampere.core.dataset.DatasetCollection` — that is the whole of what
``sbi``'s interface allows — so what goes into it is a contract in its own
right, and since W3.3 that contract is ``docs/design/contracts/encoding.md``
and lives in :mod:`ampere.core.encoding`. ``layout=`` picks the packing:

* ``"flat"`` (the default) is W3.2's fixed-size summary — each dataset's
  observed values, flattened, masked samples dropped, concatenated in
  ``datasets`` order. It is enough for a single fitting problem, where the data
  layout is fixed anyway, and for a short real vector a flow is both cheaper and
  better conditioned than one behind a set embedding.
* ``"set"`` is the coordinate-value-mask packing: one **row per sample**,
  carrying where it sits, what was measured, how well, whether it counts and
  which dataset it came from. A set-based embedding over those rows is invariant
  to how many there are and where they sit, which is what makes amortisation
  across differently-sampled datasets possible, and its uncertainty columns are
  what makes amortisation across noise realisations possible once the reserved
  observation context varies them.

Either way the run records **which** packing it used, by name and by hash
(``ampere_sbi_summary_layout``, ``ampere_encoding_layout``,
``ampere_encoding_hash``), and so does any training set it writes, because a
network trained on one layout must never be believed about another. The layout
is built from the observed containers **before** the first simulation is drawn
and every statistic in it is frozen there, so the standardisation the network
trains under and the one the observation is shown under are the same object
rather than two computations that agree.

``context=`` is the observation-context slot the horizon notes reserve
(confirmed by Peter 2026-09-09): a per-draw uncertainty pattern, grid or instrument
setting drawn from a context prior, which the embedding would condition on.
Only ``None`` is accepted today and the run records
``ampere_sbi_context = "none"``, so that the signature exists before the
machinery and a stored run from before it says so.

What a run looks like
---------------------
One "chain" of independent draws from the trained posterior, exactly as
:class:`~ampere.inference.VIEngine` emits, and for exactly the same reason:
they are i.i.d. by construction, so there is no chain structure to split and
R-hat has nothing to say about them. The training-loss trace is what an SBI run
has instead of a convergence diagnostic, and it is recorded (thinned, with the
stride).

The fourth method, and what it adds (W3.4)
-------------------------------------------
``method="tmnre"`` is truncated marginal neural ratio estimation, and it is
**not** a fourth family: it drives the same ``NRE`` trainer ``method="nre"``
does, and differs in what happens around it. :mod:`ampere.inference._tmnre`
holds the pieces and states the mathematics; the loop is
:meth:`SBIEngine._run_tmnre`. Three things about it belong here, because they
are what a reader of a *run* meets.

* **A TMNRE run's posterior is an ordinary posterior.** One estimator over the
  whole parameter vector is trained on the last round's truncated prior beside
  the marginal ones, and the run's draws come from it — one chain, i.i.d.,
  scored on the numpy contract path exactly as an NPE run's are. Nothing about
  the ``posterior``, ``sample_stats`` or ``log_likelihood`` groups is special.
* **The marginals are a group of their own.** ``marginals`` carries each 1-D
  (and, at ``marginals=2``, each 2-D) estimator's log-ratio and estimated
  marginal posterior on a grid over the final box, in both parameterisations.
  That is the corner plot the method exists to produce, and it is a *summary*
  of the estimators rather than the estimators themselves, which are torch
  modules with no place in a netCDF file.
* **The truncation is provenance.** ``ampere_sbi_truncation`` is one record per
  round — the box, its log-volume, the round's counts, the proposal's
  acceptance rate — and the boxes are nested by construction, so a run that
  shows a growing volume is a bug and not a judgement call. The box is
  observation-specific, so ``ampere_sbi_amortised`` is ``0`` for every TMNRE
  run whatever its round count: unlike a one-round NPE or NRE fit, a truncated
  estimator must not be re-conditioned on a different observation.
"""

from __future__ import annotations

import contextlib
import dataclasses
import functools
import importlib
import math
import warnings
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.encoding import (
    FLAT_KIND,
    SET_KIND,
    EncodingError,
    EncodingLayout,
    encode,
    encode_observations,
    unpack,
)
from ampere.core.exceptions import OptionalDependencyError
from ampere.core.simulate import Executor, SimulationBatch
from ampere.results.artefacts import ArtefactKey, ArtefactStore, artefact_key
from ampere.results.calibration import (
    PARAMETER_DIM,
    SBI_ROUTE,
    TARP_LEVEL_DIM,
    attach_calibration,
    calibration_dataset,
)
from ampere.results.provenance import ATTR_PREFIX
from ampere.results.training import append_training_set, write_training_set

from ._tmnre import (
    DEFAULT_TRUNCATION_EPSILON,
    GRID_POINTS_1D,
    GRID_POINTS_2D,
    MARGINAL_ORDERS,
    MARGINALS_SCHEMA_VERSION,
    MarginalEstimator,
    MarginalSummary,
    TMNREArtefact,
    TruncationBox,
    attach_marginals,
    constrained_grid,
    grid_between,
    interval_above,
    marginal_indices,
    marginal_log_density,
    pair_mesh,
    restricted_prior_class,
)
from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError, SamplingFailureWarning

__all__ = [
    "DEFAULT_TRUNCATION_EPSILON",
    "EMBEDDINGS",
    "LAYOUTS",
    "MARGINAL_ORDERS",
    "METHODS",
    "SET_EMBEDDINGS",
    "SUMMARY_LAYOUT",
    "TMNRE_SAMPLERS",
    "SBIEngine",
]

#: The three families this driver offers, and the ``sbi`` trainer each names.
#: Keys are **ampere's** vocabulary and are what a run records in
#: ``ampere_sbi_method``; the values are looked up on ``sbi.inference`` at run
#: time, so this table names no import at module level (the same discipline
#: ``_nuts.py``'s ``SAMPLER_LIBRARIES`` and ``_vi.py``'s
#: ``VARIATIONAL_LIBRARIES`` keep).
#:
#: ``sbi`` 0.27 aliases ``NPE`` to ``NPE_C``, ``NLE`` to ``NLE_A`` and ``NRE``
#: to ``NRE_B``; the run records the class that actually ran, in
#: ``ampere_sbi_trainer``, so an archived run says which variant of the family
#: produced it rather than only which family.
#:
#: **W3.4** adds a fourth key that is not a fourth *family*: ``"tmnre"`` drives
#: the same ``NRE`` trainer as ``"nre"`` and differs in what is done around it —
#: one estimator per marginal, a truncated prior between rounds, a joint
#: estimator trained on the last round's truncated prior to supply the run's
#: draws (:mod:`ampere.inference._tmnre`).
METHODS: Mapping[str, str] = {"npe": "NPE", "nle": "NLE", "nre": "NRE", "tmnre": "NRE"}

#: ampere's name for the truncated-marginal method, spelled once.
TMNRE = "tmnre"

#: The methods whose estimator is a **classifier** rather than a density, which
#: is the one thing ``sbi`` spells differently for them (``classifier=`` rather
#: than ``density_estimator=``, and ``classifier_nn`` as the builder).
_RATIO_METHODS: tuple[str, ...] = ("nre", TMNRE)

#: The default network architecture per method. NPE and NLE learn a *density*
#: and take a normalising flow; NRE learns a *classifier* and takes a residual
#: network, which is why the two are not one default.
_DEFAULT_ARCHITECTURE: Mapping[str, str] = {
    "npe": "maf",
    "nle": "maf",
    "nre": "resnet",
    "tmnre": "resnet",
}

#: Which ``sbi.neural_nets`` builder wraps an embedding net for each method,
#: and the keyword each one spells the embedding with. NRE's classifier sees θ
#: and x separately and so has two slots; the embedding is the *data* one.
_BUILDERS: Mapping[str, tuple[str, str]] = {
    "npe": ("posterior_nn", "embedding_net"),
    "nle": ("likelihood_nn", "embedding_net"),
    "nre": ("classifier_nn", "embedding_net_x"),
    "tmnre": ("classifier_nn", "embedding_net_x"),
}

#: How a TMNRE posterior draws: ``"rejection"`` samples the truncated prior and
#: accepts through the ratio, which gives genuinely i.i.d. draws and is the
#: default; ``"mcmc"`` (``sbi``'s vectorised slice sampler) is the fallback for
#: a box so narrow that rejection's acceptance collapses.
TMNRE_SAMPLERS: tuple[str, ...] = ("rejection", "mcmc")

#: The embedding vocabulary carried over from the frozen legacy
#: ``ampere/infer/sbi.py`` (read, never modified): the two shipped nets by
#: name, a user's own ``torch.nn.Module``, and a dict of hyperparameters whose
#: ``type`` key picks one of the two. The legacy spellings are all here,
#: including ``True``/``"default"``/``"Conv"`` for the CNN, so that a script
#: written against the old class keeps meaning what it meant.
EMBEDDINGS: Mapping[str, str] = {
    "CNN": "CNN",
    "Conv": "CNN",
    "default": "CNN",
    "FC": "FC",
    "FullyConnected": "FC",
    # W3.3's two, which read the coordinate-value-mask packing rather than a
    # flat feature vector and are therefore only available under a ``"set"``
    # layout. Both are ``sbi`` 0.27 nets behind an ampere wrapper that unpacks
    # the tensor by the layout and converts the one mask column into whatever
    # that net wants (``encoding.md`` §7).
    "set": "set",
    "PermutationInvariant": "set",
    "transformer": "transformer",
    "Transformer": "transformer",
}

#: The subset of :data:`EMBEDDINGS` that consumes the ``"set"`` packing. A flat
#: layout has no column groups to unpack, so asking for one of these with
#: ``layout="flat"`` is refused by name rather than silently ignored.
SET_EMBEDDINGS: tuple[str, ...] = ("set", "transformer")

#: The layouts ``layout=`` names. ``"flat"`` is W3.2's fixed-size summary and
#: stays the default: for a single fitting problem the data layout is fixed
#: anyway, and a flow over a short real vector is both cheaper and better
#: conditioned than one behind a set embedding. ``"set"`` is the
#: coordinate-value-mask packing of ``encoding.md``.
LAYOUTS: tuple[str, ...] = (FLAT_KIND, SET_KIND)

#: The default layout's name, recorded in ``ampere_sbi_summary_layout``.
SUMMARY_LAYOUT = FLAT_KIND

#: How many loss values a run records per trace. Provenance attributes are
#: JSON-encoded into the archived file, so a trace must stay small whatever the
#: epoch count; the stride is recorded, and the final value separately and
#: exactly, so nothing about the optimisation is lost that a reader would use.
_TRACE_POINTS = 200

#: Below this acceptance rate the truncated prior is warned about rather than
#: silently made expensive: a box holding a thousandth of the prior's mass
#: needs roughly a thousand prior draws per simulation, and ampere's prior is a
#: Python loop over ``sample_prior`` because ties and hierarchies are resolved
#: per draw. It is a warning and not a refusal — a tiny box is the *intended*
#: outcome of a well-behaved run on sharply informative data.
_TRUNCATION_ACCEPTANCE_FLOOR = 1e-3

#: The MCMC settings :meth:`SBIEngine.calibrate` builds its own posterior with
#: when the run's is a rejection-sampled TMNRE one. Short chains on purpose: a
#: calibration check asks for ``count`` times ``posterior_draws`` draws whose *ranks*
#: are the statistic, not for a publication-grade chain, and ``sbi``'s own
#: defaults (twenty chains, two hundred warm-up steps, automatic thinning)
#: multiply that by two orders of magnitude against a prior whose ``log_prob``
#: is a Python loop. A caller who disagrees passes ``posterior=`` their own.
_CALIBRATION_MCMC: Mapping[str, Any] = {
    "num_chains": 8,
    "warmup_steps": 25,
    "thin": 1,
    "init_strategy": "proposal",
}

#: How a TMNRE MCMC posterior is *initialised*, and it is a cost-shaped choice
#: rather than a tuning one. ``sbi``'s default, ``"resample"``, draws
#: ``num_candidate_samples`` (10 000) from the proposal and picks starting
#: points by importance weight; the proposal here is the prior restricted to
#: the box, so those 10 000 draws are themselves rejected out of ampere's
#: Python-loop prior and a narrow box makes the *initialisation* cost more than
#: the chain. ``"proposal"`` draws one starting point per chain from the same
#: distribution, which is already exactly where a truncated posterior's chain
#: should start. Everything a caller usually wants to set — chains, warm-up,
#: thinning — stays theirs, through ``posterior_options=`` at sample time.
_TMNRE_MCMC: Mapping[str, Any] = {"init_strategy": "proposal"}

#: The dtype ``sbi`` 0.27 trains in. Stated once rather than spelled at each
#: tensor: ampere's own arithmetic is float64 throughout, and the single
#: narrowing happens at the boundary, here.
_DTYPE = "float32"


# ---------------------------------------------------------------------------
# The summary: the one function W3.3 replaces
# ---------------------------------------------------------------------------


def _layout_of(problem: FittingProblem, layout: Any) -> EncodingLayout:
    """Resolve ``layout=`` into the frozen record everything downstream uses.

    A string names a kind and the layout is built from *problem*'s own observed
    containers; an :class:`~ampere.core.encoding.EncodingLayout` is taken as
    given and **checked against the problem**, field by field, because a layout
    is exactly the object a caller reuses across runs and one that does not
    describe this problem would otherwise train a network on columns that mean
    something else.
    """
    if isinstance(layout, EncodingLayout):
        try:
            layout.check_against(problem.datasets)
        except EncodingError as error:
            raise EngineError(str(error)) from error
        return layout
    if isinstance(layout, str):
        if layout not in LAYOUTS:
            known = ", ".join(LAYOUTS)
            raise EngineError(
                f"sbi does not know the layout {layout!r}. Available: {known}, or an "
                f"ampere.core.EncodingLayout of your own -- for a wider row cap, a different "
                f"Fourier band count, or a layout carried over from an earlier run."
            )
        try:
            return EncodingLayout.from_datasets(problem.datasets, kind=layout)
        except EncodingError as error:
            raise EngineError(str(error)) from error
    raise EngineError(
        f"sbi's layout= is one of {list(LAYOUTS)} or an ampere.core.EncodingLayout, got {layout!r}."
    )


def _summary_of(
    observations: Mapping[str, Any] | None,
    datasets: Mapping[str, Any],
    *,
    batched: bool,
    layout: EncodingLayout | None = None,
) -> np.ndarray:
    """The tensor the density estimator conditions on, under *layout*.

    W3.3 turned this from the whole of the layout logic into a thin adapter over
    :mod:`ampere.core.encoding`, which is where the packing now lives and where
    its contract (``encoding.md``) can be read. What survives from W3.2 is the
    guarantee, and it is the one that matters: the observation the posterior is
    conditioned on and the rows the network trained on have the same columns in
    the same order, standardised the same way, with the same samples masked.

    That guarantee is now **structural** rather than careful. The layout freezes
    every statistic and every mask at construction, from the observed data
    alone, so no two calls here can disagree about what a column means —
    including the two W3.2 could in principle have disagreed on, since
    :attr:`~ampere.core.dataset.Dataset.effective_mask` is resolved lazily and
    the observation is encoded before the first simulation is drawn.

    Parameters
    ----------
    observations
        Label to container. With *batched* true these are
        :class:`~ampere.core.simulate.ContainerBatch`\\ es, whose leading axis
        is the sample axis; with it false they are single containers.
    datasets
        The problem's own collection, which supplies the order. Reading the
        order from here rather than from *observations* is deliberate: a dict
        built elsewhere could iterate differently and the network would silently
        see its columns permuted.
    batched
        Whether *observations* carries a leading sample axis.
    layout
        The frozen packing. ``None`` builds the default ``"flat"`` one from
        *datasets*, which is what a caller predating ``layout=`` gets.

    Returns
    -------
    numpy.ndarray
        ``(count, features)`` for a ``"flat"`` layout — ``count`` is 1 for the
        unbatched form, so the observation and the training set have the same
        rank and the same code path — and ``(count, rows, columns)`` for a
        ``"set"`` one.
    """
    if not datasets:
        raise EngineError("sbi has nothing to summarise: this problem declares no datasets.")
    try:
        resolved = (
            EncodingLayout.from_datasets(datasets, kind=FLAT_KIND) if layout is None else layout
        )
        encoded = encode(observations, layout=resolved, batched=batched)
    except EncodingError as error:
        raise EngineError(str(error)) from error
    if resolved.kind == FLAT_KIND:
        return encoded.matrix
    return np.asarray(encoded.values)


# ---------------------------------------------------------------------------
# The prior bridge
# ---------------------------------------------------------------------------


@functools.cache
def _prior_class(torch: Any) -> Any:
    """The torch ``Distribution`` subclass, built once per interpreter.

    Defined inside a function because it subclasses
    ``torch.distributions.Distribution``, and this module must import no torch
    at module level. Cached on the module object so that repeated runs share
    one class — ``isinstance`` checks inside ``sbi`` then mean what they look
    like they mean.
    """

    class _UnconstrainedPrior(torch.distributions.Distribution):
        """``FittingProblem``'s joint prior, pushed forward onto ℝⁿ.

        ``sample`` draws θ in the constrained space the user declared and maps
        it through ``unconstrain``; ``log_prob`` is
        ``ParameterSet.lnprior_unconstrained``, which is ``lnprior`` plus the
        log absolute determinant of the bijection's Jacobian. Those are the two
        halves of one change of variables, so the density is the pushforward of
        the sampler by construction rather than by coincidence — and the suite
        checks the agreement numerically, because "by construction" is a claim
        about code that can stop being true.

        Hierarchical priors and ties come along for free: both are resolved by
        ``sample_prior`` and by ``lnprior`` respectively, in dependency order,
        and this class touches neither.
        """

        arg_constraints: ClassVar[dict[str, Any]] = {}
        has_rsample = False

        def __init__(
            self,
            problem: FittingProblem,
            rng: np.random.Generator,
            *,
            dtype: Any,
            device: str,
        ) -> None:
            self._problem = problem
            self._rng = rng
            self._dtype = dtype
            self._device = device
            self._size = int(problem.free_size)
            super().__init__(
                batch_shape=torch.Size([]),
                event_shape=torch.Size([self._size]),
                validate_args=False,
            )

        def __reduce__(self) -> tuple[Any, tuple[Any, ...]]:
            # The class is built inside a function (torch is imported lazily),
            # so pickle cannot find it by name; a trained posterior carries its
            # prior, and W3.5's artefact store pickles the posterior. Rebuild
            # through the module-level factory instead.
            return _rebuild_prior, (self._problem, self._rng, self._dtype, self._device)

        @property
        def support(self) -> Any:
            """All of ℝⁿ — which is the point of working in these coordinates."""
            return torch.distributions.constraints.independent(
                torch.distributions.constraints.real, 1
            )

        def sample(self, sample_shape: Any = None) -> Any:
            shape = torch.Size([]) if sample_shape is None else torch.Size(sample_shape)
            count = 1
            for extent in shape:
                count *= int(extent)
            problem = self._problem
            drawn = np.empty((count, self._size), dtype=float)
            for index in range(count):
                values = problem.parameters.pack(problem.sample_prior(self._rng))
                drawn[index] = problem.unconstrain(values)
            tensor = torch.as_tensor(drawn, dtype=self._dtype, device=self._device)
            return tensor.reshape(*shape, self._size)

        def log_prob(self, value: Any) -> Any:
            array = np.asarray(value.detach().cpu().numpy(), dtype=float)
            rows = array.reshape(-1, self._size)
            densities = np.array(
                [self._problem.parameters.lnprior_unconstrained(row) for row in rows],
                dtype=float,
            )
            return torch.as_tensor(
                densities.reshape(array.shape[:-1]), dtype=self._dtype, device=self._device
            )

        def __repr__(self) -> str:
            return f"<unconstrained joint prior, {self._size} dimension(s)>"

    return _UnconstrainedPrior


def _rebuild_prior(problem: FittingProblem, rng: Any, dtype: Any, device: str) -> Any:
    """Unpickle hook for the lazily-built prior class (see its ``__reduce__``)."""
    _, torch = _require_sbi()
    return _prior_class(torch)(problem, rng, dtype=dtype, device=device)


# ---------------------------------------------------------------------------
# The embedding vocabulary (legacy parity)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _Embedding:
    """One resolved embedding: the module, what to call it, and its width."""

    module: Any
    name: str
    output_dim: int


def _default_output_width(free_size: int, kind: str) -> int:
    """The embedding's default output width, kind-aware since **W3.11**.

    ``"flat"``-layout nets (``"CNN"``/``"FC"``) keep the legacy default,
    ``2 * free_size``: a summary that short is one more layer feeding a flow
    whose own width is unrelated to it. The ``"set"``/``"transformer"``
    embeddings *pool* a whole set of rows into this vector, so for them it is
    the conditioning vector's entire capacity rather than an intermediate
    width, and both shipped nets end in a ReLU: a narrow, randomly-initialised
    one can and does emit all zeros before any gradient step (W3.3's own tests
    set 16 explicitly, for exactly this reason). ``max(2 * free_size, 32)``
    keeps a one- or two-parameter problem's pooled embedding away from that
    floor, at the cost of a wider first layer; a user's own ``output_dim=``
    always overrides this and is never touched here.

    Ruled by Peter 2026-09-09 (WORK_ITEMS.md W3.11) as a default, not a
    finding: how the right width depends on the data, model and problem
    structure is deferred to a later embedding study.
    """
    width = 2 * int(free_size)
    if kind in SET_EMBEDDINGS:
        return max(width, 32)
    return max(width, 1)


def _embedding_of(
    embedding: Any,
    *,
    torch: Any,
    features: int,
    free_size: int,
    layout: EncodingLayout | None = None,
) -> _Embedding:
    """The ``embedding=`` vocabulary, carried over from ``ampere/infer/sbi.py``.

    The legacy class accepted four things and this accepts the same four, with
    the same spellings and the same defaults, so a script written against it
    keeps meaning what it meant: ``True``/``"default"``/``"Conv"``/``"CNN"``
    for a 1-D convolutional net, ``"FC"``/``"FullyConnected"`` for a fully
    connected one, a user's own ``torch.nn.Module``, or a dict of
    hyperparameters whose ``type`` key picks one of the two. Frozen legacy is
    read, never modified.

    The default output dimension is ``2 * free_size``, which is the legacy
    default and a sensible one: a summary narrower than twice the parameter
    count is unlikely to carry enough about the posterior's location *and*
    width. **W3.11** raises that default to ``max(2 * free_size, 32)`` for the
    ``"set"`` and ``"transformer"`` embeddings only (see
    :func:`_default_output_width`): for a pooled set, the output width *is*
    the whole of the conditioning vector's capacity rather than one more
    layer before a wider flat vector, and both shipped nets end in a ReLU, so
    a narrow, randomly-initialised one can emit all zeros before it has seen
    a single gradient step (W3.3's own tests set 16 explicitly for this
    reason). ``"flat"``'s CNN/FC default is unchanged, for legacy parity, and
    an explicit ``output_dim=`` always wins.

    One deliberate difference from the legacy code, and it is a bug fix rather
    than a change of vocabulary: an unknown ``type`` or an unknown
    hyperparameter is **refused by name** instead of silently building nothing
    (the legacy dict branch fell through with ``embedding_net`` left as the
    dict, which ``posterior_nn`` then took as a module).

    One trap is ``sbi``'s and is left to it, because its own message names the
    remedy: ``CNNEmbedding`` defaults to two convolutions with kernel 5 and a
    pool of 2, which drives a summary shorter than roughly twenty features to
    zero width and asserts. The remedy is fewer layers or a smaller kernel,
    through the hyperparameter dict — or, for the low-dimensional summaries
    the ``"flat"`` layout produces, no embedding at all, which is the default
    for exactly this reason.

    **W3.3 adds two spellings and one rule.** ``"set"`` and ``"transformer"``
    build ``sbi``'s ``PermutationInvariantEmbedding`` and ``TransformerEmbedding``
    behind an ampere wrapper that reads the coordinate-value-mask packing
    through :func:`~ampere.core.encoding.unpack` (see :func:`_set_embedding` and
    :func:`_transformer_embedding`). They need column groups, so they need a
    ``"set"`` layout, and asking for one under ``"flat"`` is refused by name
    rather than quietly given a vector with no coordinates in it. The converse
    holds too: a ``"set"`` layout with no embedding at all is refused, because a
    ``(rows, columns)`` tensor fed straight to a flow is a flow over a padded
    matrix and means nothing.
    """
    is_set = layout is not None and layout.kind == SET_KIND
    if embedding is None or embedding is False:
        if is_set:
            raise EngineError(
                "sbi was given layout='set' and no embedding. The set packing is rows of "
                "(coordinate, value, sigma, mask) that only a set-based network can read: "
                "handing it straight to a density estimator would model a padded matrix "
                "rather than an observation. Pass embedding='set' or embedding='transformer', "
                "or a torch.nn.Module of your own over ampere.core.unpack."
            )
        return _Embedding(module=None, name="none", output_dim=0)

    if embedding is True:
        kind = EMBEDDINGS["default"]
        return _built(
            kind,
            {},
            torch=torch,
            features=features,
            width=_default_output_width(free_size, kind),
            layout=layout,
        )
    if isinstance(embedding, str):
        kind = EMBEDDINGS.get(embedding)
        if kind is None:
            known = ", ".join(sorted(EMBEDDINGS))
            raise EngineError(
                f"sbi does not know the embedding {embedding!r}. Available: {known}, a "
                f"torch.nn.Module of your own, or a dict of hyperparameters with a 'type' key."
            )
        return _built(
            kind,
            {},
            torch=torch,
            features=features,
            width=_default_output_width(free_size, kind),
            layout=layout,
        )
    if isinstance(embedding, Mapping):
        settings = dict(embedding)
        named = str(settings.pop("type", "CNN"))
        kind = EMBEDDINGS.get(named)
        if kind is None:
            known = ", ".join(sorted(EMBEDDINGS))
            raise EngineError(
                f"sbi does not know the embedding type {named!r} named in the embedding= dict. "
                f"Available: {known}."
            )
        width = int(settings.pop("output_dim", _default_output_width(free_size, kind)))
        return _built(kind, settings, torch=torch, features=features, width=width, layout=layout)
    if isinstance(embedding, torch.nn.Module):
        shape = (
            (1, layout.row_cap, layout.columns_total) if is_set and layout else (1, int(features))
        )
        return _Embedding(
            module=embedding,
            name="custom",
            output_dim=_probe_width(embedding, torch=torch, shape=shape),
        )
    raise EngineError(
        f"sbi does not know what to do with embedding={embedding!r}. Pass one of "
        f"{sorted(EMBEDDINGS)}, a torch.nn.Module, a dict of hyperparameters, or None."
    )


#: Legacy hyperparameter names to ``sbi`` 0.27's, per embedding kind. The
#: legacy spellings are the contract here — a user's dict was written against
#: them — and ``kernel_size_per_layer`` is the one ``sbi`` itself renamed
#: (0.27 spells it ``kernel_size``), so both are accepted and mean the same.
_CNN_KEYS: Mapping[str, str] = {
    "n_conv_layers": "num_conv_layers",
    "num_conv_layers": "num_conv_layers",
    "out_channels_per_layer": "out_channels_per_layer",
    "kernel_size_per_layer": "kernel_size",
    "kernel_size": "kernel_size",
    "n_linear_layers": "num_linear_layers",
    "num_linear_layers": "num_linear_layers",
    "num_linear_units": "num_linear_units",
    "pool_kernel_size": "pool_kernel_size",
    "in_channels": "in_channels",
}
_FC_KEYS: Mapping[str, str] = {
    "n_layers": "num_layers",
    "num_layers": "num_layers",
    "num_hiddens": "num_hiddens",
    "enable_layer_norm": "enable_layer_norm",
}


#: W3.3's two, whose hyperparameters are ampere's own rather than a legacy
#: dictionary's, so the names are ``sbi``'s where ``sbi`` has one.
_SET_KEYS: Mapping[str, str] = {
    "trial_net_output_dim": "trial_net_output_dim",
    "num_hiddens": "num_hiddens",
    "num_layers": "num_layers",
    "aggregation_fn": "aggregation_fn",
}
_TRANSFORMER_KEYS: Mapping[str, str] = {
    "feature_space_dim": "feature_space_dim",
    "num_hidden_layers": "num_hidden_layers",
    "num_attention_heads": "num_attention_heads",
    "intermediate_size": "intermediate_size",
    "dropout": "dropout",
}


def _built(
    kind: str,
    settings: Mapping[str, Any],
    *,
    torch: Any,
    features: int,
    width: int,
    layout: EncodingLayout | None = None,
) -> _Embedding:
    """One of the four shipped nets, with its hyperparameters translated."""
    from sbi.neural_nets import embedding_nets  # pyrefly: ignore[missing-import]

    table = {
        "CNN": _CNN_KEYS,
        "FC": _FC_KEYS,
        "set": _SET_KEYS,
        "transformer": _TRANSFORMER_KEYS,
    }[kind]
    unknown = sorted(set(settings) - set(table))
    if unknown:
        raise EngineError(
            f"sbi does not know the {kind} embedding hyperparameter(s) {unknown}. "
            f"Available: {sorted(table)}, plus 'type' and 'output_dim'."
        )
    kwargs = {table[key]: value for key, value in settings.items()}
    if kind in SET_EMBEDDINGS:
        if layout is None or layout.kind != SET_KIND:
            raise EngineError(
                f"the {kind!r} embedding reads the coordinate-value-mask packing -- rows of "
                f"(coordinate, value, sigma, mask) with named column groups -- and this run's "
                f"layout is {FLAT_KIND!r}, which is one row of concatenated values and has no "
                f"column groups to read. Pass layout='set' (or an EncodingLayout) alongside "
                f"embedding={kind!r}."
            )
        builder = _set_embedding if kind == "set" else _transformer_embedding
        return builder(embedding_nets, torch=torch, layout=layout, width=int(width), **kwargs)
    if kind == "CNN":
        module = embedding_nets.CNNEmbedding(
            input_shape=(int(features),), output_dim=int(width), **kwargs
        )
    else:
        module = embedding_nets.FCEmbedding(
            input_dim=int(features), output_dim=int(width), **kwargs
        )
    return _Embedding(module=module, name=kind, output_dim=int(width))


def _set_embedding(
    embedding_nets: Any,
    *,
    torch: Any,
    layout: EncodingLayout,
    width: int,
    trial_net_output_dim: int | None = None,
    num_hiddens: int = 40,
    num_layers: int = 2,
    aggregation_fn: str = "mean",
) -> _Embedding:
    """``PermutationInvariantEmbedding`` over the packing, masked by the wrapper.

    ``sbi``'s net already has mask handling: an **all-NaN row** is treated as
    absent, its per-row embedding zeroed, and the aggregation divided by the
    surviving count. The wrapper's job is therefore to write the mask column
    into the tensor as NaN rows and hand over the rest, which is exactly what
    ``encoding.md`` §7 specifies.

    Two things about ``sbi`` 0.27 were verified before relying on any of it, and
    both are why this is a wrapper rather than a bare net.

    * ``aggregation_fn`` defaults to ``"sum"``, which makes the pooled embedding
      scale with the row count. It is set to ``"mean"`` here (W3.3 trap 5); the
      count is not lost, because ``log N`` per dataset is a column of the
      packing's ``set_features`` group.
    * The net's own valid-row count is computed as
      ``isnan(x).sum(dim=1).reshape(-1)[:num_batch]``, which reads the *first*
      batch element's count for every element whenever ``x`` has more than one
      feature column. Under this contract that is harmless and exactly right,
      because the mask is a property of the **layout** and is therefore
      identical across a batch -- but it is harmless by construction rather
      than by luck, and a future packing with a per-draw mask would have to stop
      using this net's aggregation.
    """
    row_features = layout.columns_total - layout.group("mask").width
    trial_width = int(trial_net_output_dim or max(2 * int(width), 8))
    trial = embedding_nets.FCEmbedding(
        input_dim=int(row_features),
        output_dim=trial_width,
        num_layers=int(num_layers),
        num_hiddens=int(num_hiddens),
    )
    net = embedding_nets.PermutationInvariantEmbedding(
        trial_net=trial,
        trial_net_output_dim=trial_width,
        aggregation_fn=str(aggregation_fn),
        num_hiddens=int(num_hiddens),
        num_layers=int(num_layers),
        output_dim=int(width),
        aggregation_dim=1,
    )
    wrappers = _wrapper_classes(torch)
    return _Embedding(module=wrappers[0](layout, net), name="set", output_dim=int(width))


def _transformer_embedding(
    embedding_nets: Any,
    *,
    torch: Any,
    layout: EncodingLayout,
    width: int,
    feature_space_dim: int = 32,
    num_hidden_layers: int = 2,
    num_attention_heads: int = 4,
    intermediate_size: int = 64,
    dropout: float = 0.1,
) -> _Embedding:
    """``TransformerEmbedding`` over the packing, with W3.3's four corrections.

    ``sbi`` 0.27's defaults are wrong for a set and are all overridden here
    (W3.3 trap 4): ``is_causal=True`` would make a spectrum's rows attend only
    to earlier rows, ``pos_emb="rotary"`` would encode the *index* rather than
    the coordinate, and both dropouts default to 0.5. So this builds it with
    ``is_causal=False``, ``pos_emb="none"``, explicit dropout, and the
    coordinate carried as row features (raw plus Fourier bands) by the packing
    itself.

    Two facts about that net were verified and shape the wrapper, and the first
    is the reason it cannot simply pass a mask and trust it:

    * ``forward`` **discards** ``attention_mask`` unless ``is_causal`` is true
      (``else: attention_mask = None``), so on the non-causal path the mask
      alone excludes nothing. The wrapper therefore zeroes masked and padded
      rows' *tokens* after the projection, which makes the output independent of
      what a masked row holds whatever the net does with the mask.
    * ``forward``'s own aggregation is ``hidden_states[:, -1, :]``: the
      **last** token, after full attention. A padded row is a zero token
      whose attention over the sequence is uniform, so that read is still a
      function of every row -- but *which* row it reads depends on row order
      rather than on the observation, so two encodings of one set that differ
      only in row order would train and condition on different summaries.
      **W3.11** replaces it: the wrapper never calls ``net.forward`` at all,
      and instead runs the net's own body modules directly and pools their
      output with a masked mean over every retained token (permutation
      invariant, and using every row rather than whichever lands last) --
      see :func:`_transformer_masked_mean`, which names exactly which ``sbi``
      attributes this depends on.

    The projection from the packing's feature width to the transformer's model
    dimension is ampere's: ``feature_space_dim`` is the model dimension, so
    without it the model dimension would be whatever the column count happened
    to be.
    """
    row_features = layout.columns_total - layout.group("mask").width
    heads = max(int(num_attention_heads), 1)
    model_dim = max(int(feature_space_dim), int(width), heads)
    model_dim += (-model_dim) % heads
    projection = torch.nn.Linear(int(row_features), model_dim)
    net = embedding_nets.TransformerEmbedding(
        feature_space_dim=model_dim,
        final_emb_dimension=int(width),
        is_causal=False,
        pos_emb="none",
        attention_dropout=float(dropout),
        vit_dropout=float(dropout),
        num_hidden_layers=int(num_hidden_layers),
        num_attention_heads=heads,
        num_key_value_heads=heads,
        intermediate_size=int(intermediate_size),
    )
    wrappers = _wrapper_classes(torch)
    return _Embedding(
        module=wrappers[1](layout, projection, net), name="transformer", output_dim=int(width)
    )


#: The ``sbi`` 0.27 ``TransformerEmbedding`` attributes :func:`_transformer_masked_mean`
#: reads directly, because its ``forward`` offers no hook onto the hidden states
#: it pools. Named once so a rename shows up as one clear message rather than a
#: bare ``AttributeError`` three calls deep.
_TRANSFORMER_BODY_ATTRS: tuple[str, ...] = (
    "preprocess",
    "layers",
    "norm",
    "aggregator",
    "is_causal",
)


def _transformer_masked_mean(net: Any, tokens: Any, keep: Any) -> Any:
    """The masked-mean readout head (**W3.11**), over ``sbi``'s transformer body.

    ``TransformerEmbedding.forward`` (``sbi`` 0.27) ends with
    ``self.aggregator(hidden_states[:, -1, :])`` -- the summary is the
    **last** token, after full (non-causal, unmasked) attention. Under this
    wrapper's construction (``is_causal=False``, ``pos_emb="none"``) that
    last token is a function of every row, since nothing restricts what it
    attends to, but *which* row happens to sit last is a property of row
    order, not of the observation: reversing the rows of one set trains and
    conditions the last-token read on a different summary of the same data.
    This head replaces that read with a mean over every retained token's
    hidden state, which is exactly permutation-invariant (the transformer
    body has no positional embedding under ``pos_emb="none"``, so permuting
    the input rows permutes the hidden states the same way, and a mean over
    a permuted set of vectors is the same mean) and uses every row rather
    than one.

    ``sbi`` 0.27 offers no hook for this -- ``forward`` returns a bare
    tensor, not the hidden states it pooled -- so rather than call it and
    discard its answer, this calls the net's own body modules directly, in
    exactly the order and with exactly the arguments ``forward`` does before
    its own last-token slice:

    * ``net.preprocess`` -- identity on the non-ViT path this wrapper always
      builds (``vit`` is never set), so it is called for the day that
      changes rather than assumed away;
    * ``net.layers`` -- the ``ModuleList`` of transformer blocks; each is
      called as ``block(hidden_states, attention_mask=None,
      position_ids=None)``, returning a tuple whose first element is the
      updated hidden state, which is exactly how ``forward``'s own loop
      calls them (``attention_mask`` is always ``None`` on the non-causal
      path regardless of what is passed in, and ``position_ids`` is never
      set by ``forward`` either);
    * ``net.is_causal`` -- checked and required ``False`` here: this wrapper
      only ever constructs the net that way (``encoding.md`` §7), and a mean
      over tokens that attended under a causal mask would not mean what this
      docstring says it means;
    * ``net.norm`` -- the final normalisation layer, applied once after the
      last block, exactly as ``forward`` does before it indexes out the last
      token;
    * ``net.aggregator`` -- the linear layer from ``feature_space_dim`` to
      ``final_emb_dimension``, applied here to the pooled vector instead of
      to the last token.

    If a future ``sbi`` renames or restructures any of these, this raises
    ``AttributeError`` naming the missing one (or the causal one, if that
    flips) rather than silently reading something that is no longer what it
    was -- :class:`TestTheTransformerWrapper`'s
    ``test_the_readout_depends_on_named_sbi_attributes`` in
    ``tests/inference/test_sbi.py`` fails loudly first if either happens.
    """
    missing = [name for name in _TRANSFORMER_BODY_ATTRS if not hasattr(net, name)]
    if missing:
        raise AttributeError(
            f"sbi's TransformerEmbedding no longer has the attribute(s) {missing!r}, which "
            "ampere's masked-mean readout head (docs/design/contracts/encoding.md §7, "
            "Amended W3.11) reads directly because sbi 0.27's own forward() offers no hook "
            "onto the hidden states it pools. This wrapper (_transformer_masked_mean in "
            "ampere/inference/_sbi.py) needs updating for whatever sbi replaced it with."
        )
    if net.is_causal:
        raise AttributeError(
            "sbi's TransformerEmbedding was built with is_causal=True; ampere's masked-mean "
            "readout head (_transformer_masked_mean) is only valid on the non-causal path "
            "this wrapper constructs, and something changed that."
        )
    hidden = net.preprocess(tokens)
    for block in net.layers:
        hidden = block(hidden, attention_mask=None, position_ids=None)[0]
    hidden = net.norm(hidden)
    weights = keep.unsqueeze(-1).to(hidden.dtype)
    counts = weights.sum(dim=1).clamp(min=1.0)
    pooled = (hidden * weights).sum(dim=1) / counts
    return net.aggregator(pooled)


@functools.cache
def _wrapper_classes(torch: Any) -> tuple[Any, Any]:
    """The two ``nn.Module`` wrappers, defined once per interpreter.

    Inside a function for :func:`_prior_class`'s reason: they subclass
    ``torch.nn.Module``, and this module imports torch on use, never on import.

    Both do the same three things and differ only in what they hand the net.
    They unpack ``x`` by the layout (which also reshapes it, since ``sbi``
    flattens ``x`` on some paths), they narrow float64 to the network's float32
    **once, here, at the boundary**, and they convert the packing's single mask
    column into whatever their net wants -- NaN rows for the set embedding,
    zeroed tokens plus a masked-mean readout (:func:`_transformer_masked_mean`,
    **W3.11**) for the transformer. No embedding ever sees the mask column as
    a feature: :attr:`~ampere.core.encoding.Unpacked.features` excludes it by
    construction.
    """

    class _SetEmbedding(torch.nn.Module):
        """Mask as NaN rows, then ``sbi``'s permutation-invariant net."""

        def __init__(self, layout: EncodingLayout, net: Any) -> None:
            super().__init__()
            self.layout = layout
            self.net = net

        def __reduce__(self) -> tuple[Any, tuple[Any, ...]]:
            return _rebuild_set_embedding, (self.layout, self.net)

        def forward(self, x: Any) -> Any:
            view = unpack(x.to(getattr(torch, _DTYPE)), self.layout)
            features = view.features
            absent = torch.full_like(features, float("nan"))
            return self.net(torch.where(view.valid.unsqueeze(-1), features, absent))

    class _TransformerEmbedding(torch.nn.Module):
        """Zeroed tokens, then ``sbi``'s transformer body with a masked-mean readout.

        **W3.11**: ``self.net`` is never called directly (see
        :func:`_transformer_masked_mean`) -- its own ``forward`` reads the
        *last* token as its summary, which is a function of row order rather
        than of the observation alone. This wrapper's tokens for masked and
        padded rows are already zero (below), so the net's body reaches the
        same information either way; only which token stands for the whole
        set changes.
        """

        def __init__(self, layout: EncodingLayout, projection: Any, net: Any) -> None:
            super().__init__()
            self.layout = layout
            self.projection = projection
            self.net = net

        def __reduce__(self) -> tuple[Any, tuple[Any, ...]]:
            return _rebuild_transformer_embedding, (self.layout, self.projection, self.net)

        def forward(self, x: Any) -> Any:
            view = unpack(x.to(getattr(torch, _DTYPE)), self.layout)
            keep = view.valid
            projected = self.projection(view.features)
            tokens = projected * keep.unsqueeze(-1).to(projected.dtype)
            return _transformer_masked_mean(self.net, tokens, keep)

    return _SetEmbedding, _TransformerEmbedding


def _rebuild_set_embedding(layout: EncodingLayout, net: Any) -> Any:
    """Unpickle hook for the lazily-built set wrapper (see its ``__reduce__``)."""
    _, torch = _require_sbi()
    return _wrapper_classes(torch)[0](layout, net)


def _rebuild_transformer_embedding(layout: EncodingLayout, projection: Any, net: Any) -> Any:
    """Unpickle hook for the lazily-built transformer wrapper (see its ``__reduce__``)."""
    _, torch = _require_sbi()
    return _wrapper_classes(torch)[1](layout, projection, net)


def _probe_width(module: Any, *, torch: Any, shape: tuple[int, ...]) -> int:
    """A user module's output width, measured rather than asked for.

    A ``torch.nn.Module`` does not have to declare what it produces, and the
    number is worth recording — it is what a reader of an archived run needs to
    know what the network actually saw. So it is measured with one forward pass
    on a zero input of the shape this layout produces, and a module that will
    not take it reports 0 rather than failing the run here: the real failure, if
    there is one, belongs at training time with ``sbi``'s own message.
    """
    try:
        with torch.no_grad():
            probe = module(torch.zeros(*shape))
        return int(np.asarray(probe.detach().cpu().numpy()).reshape(1, -1).shape[1])
    except Exception:
        return 0


# ---------------------------------------------------------------------------
# The engine
# ---------------------------------------------------------------------------


class SBIEngine(Engine):
    """Fit a :class:`~ampere.core.dataset.FittingProblem` by simulation.

    Neural posterior, likelihood or ratio estimation through the ``sbi``
    package: a budget of ``(θ, x)`` pairs is simulated from the problem's own
    forward model, a network is trained on them, and the trained posterior is
    evaluated at the **observed** data. The likelihood is never written down,
    which is why this is the engine for a wrapped external simulator.

    Requires the ``sbi`` extra -- ``pip install ".[sbi]"`` from a checkout
    (PyPI's ``ampere`` package is unrelated). Both ``sbi`` and torch are
    imported inside :meth:`run`, so constructing this class in the base install
    is fine and the refusal arrives when the run does.

    Parameters
    ----------
    problem
        The composed problem. Any backend: this driver consumes
        ``simulate_many`` and nothing backend-specific, and a black-box model
        on the reference backend is the case it exists for.
    method
        ``"npe"`` (neural posterior estimation -- the posterior directly, and
        the only one of the three that samples without MCMC), ``"nle"`` (the
        likelihood, sampled with MCMC), ``"nre"`` (the likelihood-to-evidence
        ratio, also sampled with MCMC), or ``"tmnre"`` (**W3.4**: the same
        ratio estimator, plus one estimator per marginal and a prior truncated
        between rounds -- see :mod:`ampere.inference._tmnre` and the
        ``marginals``/``truncation_epsilon``/``sample_with`` arguments below).
        See :data:`METHODS`.
    budget
        Simulations **per round**. The single number that decides what this
        engine costs and how good its answer is.
    embedding
        The data embedding network; see :func:`_embedding_of` and
        :data:`EMBEDDINGS`. ``None`` feeds the summary vector to the density
        estimator directly, which is the right choice for the low-dimensional
        summaries the ``"flat"`` layout produces. ``"set"`` and
        ``"transformer"`` are W3.3's two and need ``layout="set"``: they read
        the coordinate-value-mask packing through
        :func:`~ampere.core.encoding.unpack` and mask themselves from its one
        mask column, so no network ever sees that column as a feature.
    layout
        The packing the network is trained on: ``"flat"`` (the default; W3.2's
        fixed-size summary), ``"set"`` (the coordinate-value-mask encoding), or
        an :class:`~ampere.core.encoding.EncodingLayout` built by hand — for a
        wider row cap, a different Fourier band count, or a layout carried over
        from an earlier run so that the two produce comparable tensors. A layout
        that does not describe this problem is refused by name, saying which
        field differs. A ``"set"`` layout passes ``z_score_x="none"`` to
        ``sbi``: the encoding standardises itself, from statistics frozen in the
        layout.
    density_estimator
        The network architecture: a string ``sbi`` knows (``"maf"``,
        ``"nsf"``, ``"mdn"``, ``"made"`` for NPE/NLE; ``"resnet"``, ``"mlp"``
        for NRE), or a builder callable of ``sbi``'s own shape. Defaults to
        ``"maf"`` for NPE and NLE and ``"resnet"`` for NRE. A callable and an
        ``embedding=`` together are refused: the callable already decides what
        the embedding is.
    rounds
        Rounds of simulate-and-train. Round 1 draws θ from the prior; each
        round after it draws θ from the posterior the previous round trained,
        conditioned on the observed data, which concentrates the budget where
        the posterior is. More than one round makes the result **amortised no
        longer** -- the trained posterior is then only valid at this
        observation -- which is a real cost and the reason the default is 1.
        Under ``method="tmnre"`` the proposal is not the trained posterior but
        the *prior truncated to the current box*, which is what makes the
        estimate exact inside the box without an importance correction.
    marginals
        ``method="tmnre"`` only. ``1`` (the default) trains one ratio estimator
        per parameter — which is what the truncation box is built from; ``2``
        adds one per unordered pair, which is what a corner plot's off-diagonal
        panels are. The pairs cost one more estimator per pair and are trained
        only in the final round, since no box depends on them. Refused, by
        name, for the other three methods, which have no marginals to choose.
    truncation_epsilon
        ``method="tmnre"`` only: the threshold, as a fraction of each 1-D
        marginal's own maximum, that decides the truncation interval. Smaller
        keeps a wider box and spends more simulations outside the posterior;
        larger risks cutting posterior mass that no later round can recover.
        See :data:`~ampere.inference._tmnre.DEFAULT_TRUNCATION_EPSILON`.
    sample_with
        ``method="tmnre"`` only: how the final posterior draws. ``"rejection"``
        (the default) proposes from the truncated prior and accepts through the
        ratio, so the draws are genuinely i.i.d. — the property the rest of
        this class's documentation claims for an SBI run's single "chain".
        ``"mcmc"`` is ``sbi``'s vectorised slice sampler, and it takes the usual
        ``num_chains``/``warmup_steps``/``thin`` through ``posterior_options=``.

        **Rejection's cost rises as the method succeeds, and that is not a
        paradox but the arithmetic.** Its acceptance rate is roughly the
        posterior's volume over the box's, and *each* proposal draw is itself
        rejected out of the untruncated prior at the box's own prior mass —
        which is small exactly when the truncation worked. Measured on
        ``examples/sbi/tmnre_fit.py`` at its defaults: a third-round box holding
        1 % of the prior's mass, about 1 % of proposals accepted through the
        ratio, and **353.6 s against 43.8 s** for the same three rounds under
        ``"mcmc"``. ``sbi`` says so itself, in a warning naming the remedy.
        That remedy is this argument, and on a well-truncated problem it is
        usually the right one — the default is ``"rejection"`` because i.i.d.
        draws are what the rest of this class promises, not because it is the
        cheaper of the two.
    device
        ``"cpu"`` (the default), or a torch device string. CI is CPU-only by
        ruling.
    executor
        Passed straight to
        :meth:`~ampere.core.dataset.FittingProblem.simulate_many`. A
        :class:`~ampere.core.simulate.ProcessExecutor` is how an external
        simulator's budget is run in parallel; dask, ray and
        ``MPIPoolExecutor`` satisfy the protocol as they stand. A problem whose
        backend cannot pickle (jax) refuses the process pool by name, in the
        core, before any worker starts.
    chunk_size
        Also passed straight through: how many simulations exist at once. With
        ``training_set=``, it is also the size of one write.
    context
        The reserved per-draw observation-context slot (the horizon notes,
        confirmed by Peter 2026-09-09). Only ``None`` is accepted today and the
        run records ``ampere_sbi_context = "none"``; anything else is refused
        by name, saying that the machinery is a later item.
    training_set
        Write the simulated pairs to this netCDF path
        (:func:`ampere.results.training.write_training_set`), a chunk at a
        time, so a budget larger than memory reaches the file without being
        held. Failed draws are written too: ``inference.md`` §13's
        reject-and-record needs the record.
    cache_size, use_realisation
        See :class:`~ampere.inference.engine.Engine`. They govern the scoring
        of the *stored draws*, which happens on the numpy contract path after
        the fit; the fit itself scores nothing.

    Attributes
    ----------
    posterior
        ``sbi``'s trained posterior object after a run, for anything this
        driver does not expose -- ``sample`` at a different observation,
        ``map()``, ``sbi.diagnostics``. ``None`` before the first run.
    estimator
        The trained network itself.
    sampler
        The ``sbi`` trainer, as the other drivers hold their sampler.
    batch
        The **last** round's :class:`~ampere.core.simulate.SimulationBatch`,
        failures included, for a caller who wants the pairs as well as the fit.
    truncation
        ``method="tmnre"`` only: the final
        :class:`~ampere.inference._tmnre.TruncationBox`, in unconstrained
        coordinates. ``None`` for every other method.
    truncation_history
        One record per round: that round's box in both parameterisations, its
        log-volume, its simulation counts, and the rate at which the box that
        produced those draws accepted prior samples. This is what
        ``ampere_sbi_truncation`` carries, and it is the run's own evidence
        that the truncation behaved — the volumes must not grow.
    marginal_estimators
        The trained :class:`~ampere.inference._tmnre.MarginalEstimator`\\ s by
        index tuple, for a caller who wants to evaluate one somewhere the
        stored grid does not reach. Empty after a cache hit, which restores a
        run's *stored* content and not the networks behind it.
    marginal_summary
        Those estimators evaluated over the final box — what the run's
        ``marginals`` group holds.

    Examples
    --------
    A hand-written model, fitted without ever writing down its likelihood:

    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
    >>> from ampere.inference import SBIEngine
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
    >>> problem = FittingProblem(Line(grid), [Dataset(observed)], seed=20260909)
    >>> run = SBIEngine(problem, budget=2000).run(draws=500)  # doctest: +SKIP
    >>> run.attrs["ampere_engine"], run.attrs["ampere_sbi_method"]  # doctest: +SKIP
    ('sbi', 'npe')

    Notes
    -----
    **The failure path is the reason this engine exists.** A budget of 10⁴
    draws with a 2 % crash rate produces 9 800 usable pairs and a count, and
    the count is in the attrs (``ampere_sbi_failures``, and the per-reason
    breakdown in ``ampere_failure_counts``). Failed draws are dropped from the
    training set the network sees -- training on garbage is exactly what
    ``inference.md`` §13's reject-and-record exists to prevent -- and written
    to the ``training_set=`` file if one was asked for, because a budget's
    failure rate is a measurement of the prior worth keeping.

    **A run is only as good as its budget, and nothing here checks that.** A
    density estimator always returns a posterior; whether it is the right one
    is what the training-loss trace, and W3.6's simulation-based calibration,
    are for. The attrs carry what a reader needs to judge it.
    """

    NAME: ClassVar[str] = "sbi"
    #: Simulation-based inference asks nothing of a gradient: the whole point
    #: is that the forward model may be a compiled binary. Passed to
    #: ``check_engine`` explicitly, as every other driver does, because the
    #: question is "can sbi run this?" and not "is this problem differentiable?"
    OFFERS_GRADIENTS: ClassVar[bool] = False

    def __init__(
        self,
        problem: FittingProblem,
        *,
        method: str = "npe",
        budget: int = 1000,
        embedding: Any = None,
        layout: Any = SUMMARY_LAYOUT,
        density_estimator: Any = None,
        rounds: int = 1,
        marginals: int = 1,
        truncation_epsilon: float = DEFAULT_TRUNCATION_EPSILON,
        sample_with: str | None = None,
        device: str = "cpu",
        executor: Executor | None = None,
        chunk_size: int | None = None,
        context: Any = None,
        training_set: str | Path | None = None,
        cache: ArtefactStore | None = None,
        cache_size: int = DEFAULT_CACHE_SIZE,
        use_realisation: bool = True,
    ) -> None:
        chosen = str(method).lower()
        if chosen not in METHODS:
            known = ", ".join(sorted(METHODS))
            raise EngineError(
                f"sbi does not know the method {method!r}. Available: {known}. 'npe' estimates "
                f"the posterior directly and is the only one that samples without MCMC; 'nle' "
                f"estimates the likelihood and 'nre' the likelihood-to-evidence ratio, and both "
                f"are then sampled with MCMC; 'tmnre' is 'nre' with one estimator per marginal "
                f"and a prior truncated between rounds, and is no longer amortised."
            )
        if int(budget) < 1:
            raise EngineError(f"sbi needs a simulation budget of at least 1, got {budget}.")
        if int(rounds) < 1:
            raise EngineError(f"sbi needs at least one round, got {rounds}.")
        order = int(marginals)
        epsilon = float(truncation_epsilon)
        sampler = None if sample_with is None else str(sample_with).lower()
        if chosen == TMNRE:
            if order not in MARGINAL_ORDERS:
                known = ", ".join(str(value) for value in MARGINAL_ORDERS)
                raise EngineError(
                    f"tmnre's marginals= is the order of the marginals it estimates: {known}. "
                    f"1 trains one ratio estimator per parameter (which is what the truncation "
                    f"box is built from); 2 adds one per unordered pair, which is what a corner "
                    f"plot's off-diagonal panels are. Got {marginals!r}."
                )
            if not 0.0 < epsilon < 1.0:
                raise EngineError(
                    f"tmnre's truncation_epsilon= is a fraction of each 1-D marginal's own "
                    f"maximum, so it lies strictly between 0 and 1 (the default is "
                    f"{DEFAULT_TRUNCATION_EPSILON}). A smaller value keeps a wider box and costs "
                    f"simulations; a larger one risks cutting posterior mass no later round can "
                    f"recover. Got {truncation_epsilon!r}."
                )
            if sampler is not None and sampler not in TMNRE_SAMPLERS:
                known = ", ".join(TMNRE_SAMPLERS)
                raise EngineError(
                    f"tmnre's sample_with= is one of {known}, got {sample_with!r}. 'rejection' "
                    f"draws from the truncated prior and accepts through the ratio, so the draws "
                    f"are i.i.d.; 'mcmc' is the fallback when a narrow box makes rejection's "
                    f"acceptance rate collapse."
                )
        else:
            if order != 1:
                raise EngineError(
                    f"sbi's marginals= belongs to method='tmnre', which trains one ratio "
                    f"estimator per marginal; method={chosen!r} trains one estimator over the "
                    f"whole parameter vector and has no marginals to choose. Got {marginals!r}."
                )
            if epsilon != DEFAULT_TRUNCATION_EPSILON:
                raise EngineError(
                    f"sbi's truncation_epsilon= belongs to method='tmnre', which truncates the "
                    f"prior between rounds; method={chosen!r} does not truncate anything. Got "
                    f"{truncation_epsilon!r}."
                )
            if sampler is not None:
                raise EngineError(
                    f"sbi's sample_with= belongs to method='tmnre', whose posterior is built "
                    f"against a truncated prior; for method={chosen!r} sbi's own default sampler "
                    f"for the family is used and posterior_options= is where its settings go. "
                    f"Got {sample_with!r}."
                )
        if context is not None:
            raise EngineError(
                f"sbi's context= is the reserved per-draw observation-context slot (a sigma "
                f"pattern, a grid, an instrument setting, drawn from a context prior and seen by "
                f"the embedding). The signature exists so that runs recorded before the "
                f"machinery say so -- every run today records ampere_sbi_context = 'none' -- but "
                f"the machinery is a later item, so only None is accepted. Got {context!r}."
            )
        if density_estimator is not None and not isinstance(density_estimator, str):
            if not callable(density_estimator):
                raise EngineError(
                    f"sbi's density_estimator= is an architecture name or a builder callable of "
                    f"sbi's own shape, got {density_estimator!r}."
                )
            if embedding is not None:
                raise EngineError(
                    "sbi was given both a density_estimator= builder and an embedding=. A "
                    "builder already decides what its embedding is, so passing both says two "
                    "different things about one network. Pass the embedding to your builder, or "
                    "name the architecture as a string and let embedding= wrap it."
                )
        super().__init__(problem, cache_size=cache_size, use_realisation=use_realisation)
        self.method = chosen
        self.budget = int(budget)
        self.embedding = embedding
        self.layout = layout
        #: The resolved :class:`~ampere.core.encoding.EncodingLayout` after a
        #: run: what the network was trained under, and what an observation
        #: shown to it later must match by hash.
        self.encoding: EncodingLayout | None = None
        self.density_estimator_spec = density_estimator
        self.rounds = int(rounds)
        #: TMNRE's order (1 or 2), threshold and posterior sampler. Meaningless
        #: for the other three, which the constructor refuses to let a caller
        #: set at all rather than accepting and ignoring.
        self.marginals = order
        self.truncation_epsilon = epsilon
        self.sample_with = (sampler or TMNRE_SAMPLERS[0]) if chosen == TMNRE else None
        self.device = str(device)
        self.executor = executor
        self.chunk_size = None if chunk_size is None else int(chunk_size)
        self.training_set = None if training_set is None else str(training_set)
        #: W3.5: a trained-artefact store keyed on the problem's own hashes, or
        #: ``None`` to train every time. A hit skips simulation and training.
        self.cache = cache
        #: ``sbi``'s trained posterior, network and trainer after a run.
        self.posterior: Any = None
        self.estimator: Any = None
        #: The last round's batch, failures included.
        self.batch: SimulationBatch | None = None
        #: **W3.4**, after a TMNRE run: the final truncation box, one record
        #: per round (what the run's ``ampere_sbi_truncation`` is built from),
        #: the trained marginal estimators by index tuple, and the summary the
        #: ``marginals`` group carries. All ``None``/empty otherwise.
        self.truncation: TruncationBox | None = None
        self.truncation_history: list[dict[str, Any]] = []
        self.marginal_estimators: dict[tuple[int, ...], MarginalEstimator] = {}
        self.marginal_summary: MarginalSummary | None = None
        self._proposal: Any = None
        self._written = False

    # -- the run --------------------------------------------------------------

    def run(
        self,
        draws: int,
        *,
        training: Mapping[str, Any] | None = None,
        posterior_options: Mapping[str, Any] | None = None,
        progress: bool = False,
    ) -> Any:
        """Simulate, train, sample the trained posterior, and emit the run.

        Parameters
        ----------
        draws
            Independent draws from the trained posterior at the observed data.
            They are the run's posterior, in one chain, as
            :class:`~ampere.inference.VIEngine`'s guide draws are and for the
            same reason. For NPE they are cheap; for NLE and NRE each one costs
            an MCMC step, so a large number is not free.

            One cost is easy to overlook and is the same for all three: every
            stored draw is **scored on the numpy contract path** afterwards, so
            ``draws`` is also a count of full ``problem.evaluate`` calls — one
            forward-model evaluation each. That is what buys the true per-draw
            split beside the estimator's own density, and it is why a draw
            count here is not the free number it is for a variational fit.
            ``engine_draws_recomputed`` in the attrs is that count.
        training
            Forwarded verbatim to the ``sbi`` trainer's ``train`` --
            ``max_num_epochs``, ``training_batch_size``, ``learning_rate``,
            ``stop_after_epochs``, ``validation_fraction``. ampere interposes
            no defaults: they are ``sbi``'s and they are good ones. Under
            ``method="tmnre"`` the same settings train every estimator, the
            marginal ones included.
        posterior_options
            Forwarded verbatim to the trained posterior's ``sample``. This is
            where an MCMC posterior's ``num_chains``, ``warmup_steps``,
            ``thin`` and ``init_strategy`` go, which is what makes an NLE or
            NRE run affordable at a small draw count.

            For a TMNRE posterior under the default ``sample_with="rejection"``
            it is also where ``num_samples_to_find_max`` and
            ``max_sampling_batch_size`` go, and they are worth knowing about:
            ``sbi`` defaults both to 10 000 draws *from the proposal*, and this
            driver's proposal is the prior restricted to the box, whose own
            draws are rejected from a prior that is a Python loop over
            ``sample_prior``. A narrow box therefore makes those defaults the
            dominant cost of a run. No default is interposed here — that is
            W3.2 decision (5) and it stands — but a small run should say so
            explicitly.
        progress
            Show ``sbi``'s own progress bars, for training and for sampling.
            Off by default: a driver that prints by default is unusable inside
            a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(1, draws)``, the true
            per-draw ``log_prior``/``log_likelihood`` split scored on the numpy
            path, and ``ampere_sbi_log_prob`` in ``sample_stats`` beside them.
        """
        if int(draws) < 1:
            raise EngineError(f"sbi needs at least one draw, got {draws}.")
        sbi_package, torch = _require_sbi()
        self.start()
        self._written = False

        problem = self.problem
        dtype = getattr(torch, _DTYPE)
        # The layout is built **first**, from the observed containers alone, and
        # everything after it -- the observation, every round's rows, the
        # embedding's input width, the training set's attrs -- is derived from
        # it. That ordering is the contract (``encoding.md`` §4): a statistic
        # computed from the simulations could not be the same at inference.
        layout = _layout_of(problem, self.layout)
        self.encoding = layout
        try:
            observation = np.asarray(encode_observations(problem.datasets, layout=layout).values)
        except EncodingError as error:
            raise EngineError(str(error)) from error
        if layout.kind == FLAT_KIND:
            observation = observation.reshape(1, -1)
        features = int(np.prod(observation.shape[1:]))

        # W3.15: seeded once here, before the first network is built, so that
        # an eager embedding's initial weights are as reproducible as the
        # rest of the run; seeded again before every round's ``train()``
        # (below, and in ``_run_tmnre``) for training's own randomness (batch
        # order, dropout), and once more (a distinct concern) immediately
        # before every posterior draw, so sampling repeats whether this run
        # trained or hit the cache. ``torch_seed`` is this call's own draw --
        # the first from ``"sbi.torch"`` -- and is what the attrs record.
        with self._seeded(torch, "sbi.torch") as torch_seed:
            resolved = _embedding_of(
                self.embedding,
                torch=torch,
                features=features,
                free_size=problem.free_size,
                layout=layout,
            )
            trainer, architecture = self._trainer(
                sbi_package, torch, embedding=resolved, layout=layout, progress=progress
            )
        self.sampler = trainer

        simulated = 0
        usable = 0
        proposal: Any = None
        observation_tensor = torch.as_tensor(observation[0], dtype=dtype, device=self.device)
        # W3.5: the artefact key is the problem's own hashes plus everything
        # that shaped the estimator -- the encoding *hash* stands in for the
        # layout, so a differently-packed observation is a miss by construction.
        # W3.12: marginals=/truncation_epsilon=/sample_with= are now
        # artefact_key's own keywords -- the _key_architecture stopgap that
        # folded TMNRE's settings into the architecture string is gone.
        # sample_with belongs in the key because sbi bakes it into the built
        # posterior, which is exactly what the store holds (_artefact).
        cache_key: ArtefactKey | None = (
            None
            if self.cache is None
            else artefact_key(
                problem,
                layout=layout.hash,
                method=self.method,
                architecture=architecture,
                budget=self.budget,
                rounds=self.rounds,
                marginals=self.marginals if self.method == TMNRE else None,
                truncation_epsilon=self.truncation_epsilon if self.method == TMNRE else None,
                sample_with=self.sample_with if self.method == TMNRE else None,
            )
        )
        cached = None if cache_key is None or self.cache is None else self.cache.get(cache_key)
        cache_hit = cached is not None
        if cache_hit:
            self._restore(cached, observation_tensor)
        elif self.method == TMNRE:
            simulated, usable = self._run_tmnre(
                sbi_package,
                torch,
                trainer=trainer,
                layout=layout,
                features=features,
                dtype=dtype,
                observation=observation_tensor,
                training=training,
                progress=progress,
            )
        else:
            for round_index in range(self.rounds):
                theta, summary, counts = self._simulate_round(
                    round_index, proposal=proposal, torch=torch
                )
                simulated += counts[0]
                usable += counts[1]
                self._append(
                    trainer,
                    theta,
                    summary,
                    torch=torch,
                    dtype=dtype,
                    proposal=proposal,
                    round_index=round_index,
                )
                with self._seeded(torch, "sbi.torch"):
                    self.estimator = trainer.train(show_train_summary=False, **dict(training or {}))
                self.posterior = trainer.build_posterior(self.estimator)
                self.posterior.set_default_x(observation_tensor)
                proposal = self.posterior
        if not cache_hit and cache_key is not None and self.cache is not None:
            self.cache.put(cache_key, self._artefact())

        # A distinct concern from training's, and reseeded here rather than
        # relying on training's own seed to carry through: a cache hit skips
        # training altogether, so sampling must be pinned on its own for the
        # hit and the miss to each repeat from run to run (W3.15).
        with self._seeded(torch, "sbi.torch.sample"):
            drawn = self.posterior.sample(
                (int(draws),), show_progress_bars=progress, **dict(posterior_options or {})
            )
        unconstrained = np.asarray(drawn.detach().cpu().numpy(), dtype=float).reshape(
            int(draws), problem.free_size
        )
        estimator_log_prob = self._estimator_log_prob(drawn, torch=torch)
        chain = np.stack([[problem.constrain(row) for row in unconstrained]])

        attrs = self._attrs(
            sbi_package,
            torch,
            trainer=trainer,
            architecture=architecture,
            embedding=resolved,
            layout=layout,
            features=features,
            simulated=simulated,
            usable=usable,
            draws=int(draws),
            torch_seed=torch_seed,
        )
        if cache_key is not None:
            attrs["sbi_cache_hit"] = int(cache_hit)
            attrs["sbi_cache_key"] = cache_key.digest()
        tree = self.finish(chain, extra_attrs=attrs)
        tree = _with_estimator_log_prob(tree, estimator_log_prob)
        if self.marginal_summary is not None:
            attach_marginals(
                tree,
                self.marginal_summary.to_dataset(
                    attrs={
                        f"{ATTR_PREFIX}marginals_schema_version": MARGINALS_SCHEMA_VERSION,
                        f"{ATTR_PREFIX}marginals_order": self.marginals,
                        f"{ATTR_PREFIX}marginals_epsilon": self.truncation_epsilon,
                        f"{ATTR_PREFIX}marginals_parameterisation": "unconstrained",
                        f"{ATTR_PREFIX}marginals_method": self.method,
                    }
                ),
            )
        return tree

    # -- family D's fast path (W3.6) ------------------------------------------

    def calibrate(
        self,
        *,
        count: int,
        posterior_draws: int = 1000,
        tarp: bool = True,
        levels: Any = None,
        progress: bool = False,
        posterior: Any = None,
        attach_to: Any = None,
        num_workers: int = 1,
    ) -> Any:
        """Family D on this run's trained posterior: SBC ranks and coverage.

        ``diagnostics.md`` §11's calibration family, on the route where it is
        cheap. An amortised posterior can be re-conditioned on a fresh dataset
        for nothing, so the whole of simulation-based calibration costs one
        extra simulation batch and **no retraining at all** — where
        :func:`ampere.results.calibration.sbc` pays for a full fit per
        simulation. The arithmetic is ``sbi``'s own ``run_sbc``/``check_sbc``
        and ``run_tarp``/``check_tarp``, per §2.2's depend-don't-reimplement
        posture; what this method owns is that they are handed the *right*
        tensors.

        Which is the one thing here that could go silently wrong. The fresh
        batch is encoded with **this run's own**
        :class:`~ampere.core.encoding.EncodingLayout` (:attr:`encoding`), not a
        layout rebuilt from the simulated containers: a network is only
        meaningful about a tensor packed the way it was trained, and a
        calibration check that re-derived its standardisation from the
        calibration batch would be testing a different network from the one the
        run produced — and would report *that* one as calibrated.

        **W3.4 adds two things for a truncated run, both recorded.** A TMNRE
        posterior is defined *on its box*, so the distribution it is calibrated
        against is the truncated prior and not the original one — calibrating
        it against the full prior would report miscalibration for every truth
        the box excludes, which is an artefact of the method rather than a
        property of the estimator (``ampere_calibration_reference``). And a
        run whose posterior samples by rejection is checked through an MCMC
        posterior over the *same* trained estimator
        (``ampere_calibration_sampler``), because a rejection posterior pays a
        fixed maximisation stage per conditioning observation and SBC
        re-conditions at ``count`` of them. Pass *posterior* to override either.

        θ is compared in the **unconstrained** parameterisation, which is where
        the estimator lives (see this module's docstring). That costs nothing
        and changes nothing: :meth:`~ampere.core.dataset.FittingProblem.constrain`
        is monotone per parameter, so a rank is the same number in either space,
        and the labels are the problem's own ``free_labels()`` either way.

        Parameters
        ----------
        count
            How many fresh prior draws to simulate and test at. ``sbi``'s own
            checks warn below 100, and so does this: the uniformity test is a
            goodness-of-fit test on this many points.
        posterior_draws
            ``L``: how many posterior draws each rank is taken against.
        tarp
            Also run the TARP expected-coverage test (Lemos et al. 2023) and
            store its curve. It is the *joint* diagnostic the marginal ranks
            cannot replace — a posterior can be calibrated in every margin and
            wrong about the correlations — and it costs a second batch of
            posterior sampling, which for an amortised NPE is cheap and for an
            MCMC-sampled NLE or NRE is not.
        levels
            Nominal levels for the marginal coverage curve; the module default
            (:data:`~ampere.results.calibration.DEFAULT_LEVELS`) otherwise.
        progress
            Show ``sbi``'s progress bars.
        posterior
            Test *this* posterior rather than :attr:`posterior`. The seam a
            calibration study of a **deliberately** miscalibrated posterior
            needs — a temperature-scaled wrapper, say — and the reason the
            suite can prove this check fails when it should.
        attach_to
            A run's ``DataTree`` to write the result into, as its
            ``calibration`` group (``results.md`` §4). Modified in place; the
            dataset is returned either way.
        num_workers
            Forwarded to ``sbi``, which uses it only on the non-batched
            sampling path (an NLE or NRE posterior).

        Returns
        -------
        xarray.Dataset
            The ``calibration`` group: ``ranks``, ``coverage``, ``ks_pvalue``,
            ``c2st_ranks``, and — with *tarp* — ``tarp_coverage`` on its own
            credibility grid, with the checks' scalar verdicts on the attrs.

        Raises
        ------
        ampere.inference.EngineError
            If this engine has not run, so there is no trained posterior and no
            layout to encode the calibration batch with; or if every simulation
            in the calibration batch failed.
        """
        sbi_package, torch = _require_sbi()
        target = self.posterior if posterior is None else posterior
        # W3.4: a rejection-sampled TMNRE posterior is the wrong object to run
        # SBC against, and unusably so rather than merely slowly. A rejection
        # posterior pays a fixed find-the-maximum stage of 10 000 proposal
        # draws *per conditioning observation*, and every one of those draws is
        # itself rejected against the box out of ampere's Python-loop prior;
        # SBC re-conditions the posterior at ``count`` fresh observations, so
        # the check costs that stage ``count`` times over and does not finish.
        # The same trained estimator sampled by slice MCMC answers the same
        # question at a cost that is linear in the draws, so the check builds
        # one -- recorded in ``ampere_calibration_sampler`` -- and a caller who
        # wants something else passes ``posterior=`` as they always could.
        rebuilt = False
        if posterior is None and self.method == TMNRE and self.sample_with == "rejection":
            target = self.sampler.build_posterior(
                self.estimator,
                prior=self._proposal,
                sample_with="mcmc",
                mcmc_parameters=_CALIBRATION_MCMC,
            )
            rebuilt = True
        if target is None or self.encoding is None:
            raise EngineError(
                "sbi cannot calibrate a posterior it has not trained: call run() first, or pass "
                "posterior= a posterior of your own. The encoding layout the calibration batch "
                "must be packed under is recorded by run() and by nothing else."
            )
        simulations = int(count)
        draws = int(posterior_draws)
        if simulations < 1:
            raise EngineError(f"a calibration check needs at least one simulation, got {count}.")
        if draws < 1:
            raise EngineError(f"a rank needs at least one posterior draw, got {posterior_draws}.")

        problem = self.problem
        dtype = getattr(torch, _DTYPE)
        rng = self.stream("calibrate")
        # W3.4: a TMNRE posterior is defined **on its truncation box**, so the
        # distribution it must be calibrated against is the truncated prior and
        # not the original one. Calibrating a truncated posterior against the
        # full prior would report miscalibration for every truth the box
        # excludes -- which is an artefact of the method, not a property of the
        # estimator -- so the batch and the reference draws both come from the
        # box. The attrs say which of the two happened.
        truncated = self.method == TMNRE and self._proposal is not None
        given = self._truncated_theta(simulations) if truncated else None
        thetas, summaries, simulated = self._calibration_batch(simulations, rng, values=given)
        theta_tensor = torch.as_tensor(thetas, dtype=dtype, device=self.device)
        summary_tensor = torch.as_tensor(summaries, dtype=dtype, device=self.device)
        count_kept = int(thetas.shape[0])
        prior = (
            self._truncated_theta(count_kept, unconstrained=True)
            if truncated
            else np.stack(
                [problem.unconstrain(problem.sample_prior(rng)) for _ in range(count_kept)]
            )
        )
        prior_tensor = torch.as_tensor(prior, dtype=dtype, device=self.device)

        diagnostics = importlib.import_module("sbi.diagnostics")
        ranks, dap = diagnostics.run_sbc(
            theta_tensor,
            summary_tensor,
            target,
            num_posterior_samples=draws,
            num_workers=num_workers,
            show_progress_bar=progress,
        )
        checks = diagnostics.check_sbc(ranks, prior_tensor, dap, num_posterior_samples=draws)
        extras: dict[str, tuple[tuple[str, ...], Any]] = {
            "ks_pvalue": ((PARAMETER_DIM,), _numpy(checks["ks_pvals"])),
            "c2st_ranks": ((PARAMETER_DIM,), _numpy(checks["c2st_ranks"])),
        }
        coords: dict[str, Any] = {}
        attrs: dict[str, Any] = {
            f"{ATTR_PREFIX}calibration_uniformity_check": "sbi.diagnostics.check_sbc",
            f"{ATTR_PREFIX}calibration_requested": simulations,
            f"{ATTR_PREFIX}calibration_failures": simulated - int(thetas.shape[0]),
            f"{ATTR_PREFIX}calibration_engine": type(self).__name__,
            f"{ATTR_PREFIX}calibration_method": self.method,
            f"{ATTR_PREFIX}calibration_encoding_layout": self.encoding.name,
            f"{ATTR_PREFIX}calibration_encoding_hash": self.encoding.hash,
            f"{ATTR_PREFIX}calibration_parameterisation": "unconstrained",
            f"{ATTR_PREFIX}calibration_reference": "truncated_prior" if truncated else "prior",
            f"{ATTR_PREFIX}calibration_sampler": "mcmc" if rebuilt else "as_run",
            f"{ATTR_PREFIX}calibration_c2st_dap": float(np.mean(_numpy(checks["c2st_dap"]))),
            f"{ATTR_PREFIX}calibration_sbi_version": str(
                getattr(sbi_package, "__version__", "unknown")
            ),
        }
        if problem.seed is not None:
            attrs[f"{ATTR_PREFIX}calibration_seed"] = int(problem.seed)
        if tarp:
            ecp, alpha = diagnostics.run_tarp(
                theta_tensor,
                summary_tensor,
                target,
                num_posterior_samples=draws,
                num_workers=num_workers,
                show_progress_bar=progress,
            )
            area, ks_pvalue = diagnostics.check_tarp(ecp, alpha)
            extras["tarp_coverage"] = ((TARP_LEVEL_DIM,), _numpy(ecp))
            coords[TARP_LEVEL_DIM] = _numpy(alpha)
            attrs[f"{ATTR_PREFIX}calibration_tarp_atc"] = float(area)
            attrs[f"{ATTR_PREFIX}calibration_tarp_ks_pvalue"] = float(ks_pvalue)

        result = calibration_dataset(
            _numpy(ranks).astype(np.int64),
            problem.free_labels(),
            posterior_draws=draws,
            route=SBI_ROUTE,
            levels=levels,
            extras=extras,
            coords=coords,
            attrs=attrs,
        )
        if attach_to is not None:
            attach_calibration(attach_to, result)
        return result

    def _truncated_theta(self, count: int, *, unconstrained: bool = False) -> np.ndarray:
        """*count* draws from the final truncated prior, for a calibration batch.

        In constrained coordinates by default, which is what
        :meth:`~ampere.core.dataset.FittingProblem.simulate_many`'s ``values=``
        takes; *unconstrained* gives the coordinates ``sbi``'s checks compare
        ranks in.
        """
        drawn = self._proposal.sample((int(count),), show_progress_bars=False)
        rows = np.asarray(drawn.detach().cpu().numpy(), dtype=float).reshape(
            int(count), self.problem.free_size
        )
        if unconstrained:
            return rows
        return np.stack([self.problem.constrain(row) for row in rows])

    def _calibration_batch(
        self, count: int, rng: Any, *, values: Any = None
    ) -> tuple[np.ndarray, np.ndarray, int]:
        """A fresh batch, unconstrained θ and layout-encoded x.

        Streamed in chunks exactly as :meth:`_simulate_round` streams a
        training round, and on its own named sub-stream (``"sbi.calibrate"``),
        because ``lowering.md`` §9.2's rule is that adding a diagnostic must not
        change what a fit simulated.

        *values* is ``None`` for a prior batch and one constrained θ per draw
        for W3.4's truncated one; the rest of the method does not care which.
        """
        problem = self.problem
        thetas: list[np.ndarray] = []
        summaries: list[np.ndarray] = []
        simulated = 0
        for chunk in problem.simulate_many(
            count,
            values=values,
            observe=True,
            rng=rng,
            executor=self.executor,
            chunk_size=self.chunk_size,
            as_chunks=True,
        ):
            simulated += len(chunk)
            keep = chunk.usable
            if not len(keep):
                continue
            thetas.append(np.stack([problem.unconstrain(draw.theta) for draw in keep]))
            summaries.append(
                _summary_of(keep.observations, problem.datasets, batched=True, layout=self.encoding)
            )
        if not thetas:
            raise EngineError(
                f"sbi has nothing to calibrate on: every one of the {simulated} simulation(s) in "
                f"the calibration batch failed. problem.failure_summary() says why:\n"
                f"{problem.failure_summary()}"
            )
        return np.concatenate(thetas, axis=0), np.concatenate(summaries, axis=0), simulated

    # -- W3.15: torch has a global generator this driver does not own ---------

    @contextlib.contextmanager
    def _seeded(self, torch: Any, concern: str) -> Iterator[int | None]:
        """Seed torch's (and, where it matters, numpy's *legacy* global) generator
        for one block, from this run's own sub-stream, then put both back.

        Neither ``ampere`` nor ``sbi`` seeds torch: network initialisation and
        a trainer's batch order are drawn from torch's *global* generator, the
        one piece of randomness ``problem.rng``'s own sub-streams cannot reach,
        so two runs of the same seeded problem disagreed even though
        ``problem.seed`` fixed every simulation. This is ``lowering.md``
        §9.1's route (1) — seed the global stream from this engine's own
        sub-stream — the same one ``_nuts.py`` and ``_vi.py`` use for pyro.

        **numpy's legacy global generator too**, for the same reason
        ``_zeus.py``'s ``_global_seed`` seeds more than one library's global
        state: ``sbi``'s default MCMC method (``"slice_np_vectorized"``, what
        an NLE/NRE posterior and a TMNRE ``sample_with="mcmc"`` one both use)
        draws its slice proposals through ``np.random`` directly rather than
        through anything ``torch.manual_seed`` reaches, so an MCMC-sampled
        posterior stayed irreproducible even after torch was seeded — found
        by this item's own TMNRE acceptance run refusing to repeat. An NPE
        posterior's flow never touches it, so seeding it unconditionally here
        costs that path nothing.

        **Restored on exit**, this driver's own version of the same rule
        ``_nuts.py``'s ``torch.random.fork_rng`` and ``_zeus.py``'s
        ``_global_seed`` both state: ampere does not leave a library's global
        random state changed behind it, so a caller's own unrelated use of
        ``np.random`` or torch after ``run()`` returns is exactly as
        (ir)reproducible as it would have been had this engine done nothing.
        Restoring the *generator's* state is not undoing anything already
        computed — a trained network's weights, once drawn, stay drawn — it
        only affects whatever draws next.

        Called with the *same* ``concern`` more than once in a run, this
        yields a different seed each time: ``FittingProblem.rng``'s generator
        for a label is created once and then advanced, so repeated draws
        differ while the whole sequence is reproducible for a given seed and
        call order (``dataset.py``'s ``rng`` docstring) — which is what makes
        the training seed "distinct per round" for free, from one call
        written once and made before every round's training rather than from
        any per-round label arithmetic here.

        ``problem.seed is None`` is left alone, on purpose: it means this run
        asked not to be reproducible, and seeding torch (or numpy) anyway
        would make it look reproducible without being asked (``integer_seed``'s
        own contract, ``engine.py``).
        """
        if self.problem.seed is None:
            yield None
            return
        seed = self.integer_seed(concern)
        torch_state = torch.random.get_rng_state()
        cuda_states = torch.cuda.get_rng_state_all() if torch.cuda.is_available() else None
        numpy_state = np.random.get_state()
        try:
            torch.manual_seed(seed)
            if "cuda" in self.device:
                torch.cuda.manual_seed_all(seed)
            np.random.seed(seed)
            yield seed
        finally:
            torch.random.set_rng_state(torch_state)
            if cuda_states is not None:
                torch.cuda.set_rng_state_all(cuda_states)
            np.random.set_state(numpy_state)

    # -- the pieces -----------------------------------------------------------

    def _trainer(
        self,
        sbi_package: Any,
        torch: Any,
        *,
        embedding: _Embedding,
        layout: EncodingLayout,
        progress: bool,
    ) -> tuple[Any, str]:
        """``sbi``'s trainer for this method, with the prior and the network.

        The one thing the layout changes here is ``z_score_x``. ``sbi``
        standardises ``x`` column-wise from the training tensor by default
        (``z_score_x="independent"``), and on the set packing that would z-score
        the mask, the dataset index and the coordinate columns, and let padded
        rows enter every column's statistics. So a set layout passes
        ``z_score_x="none"`` and **the encoding does its own standardisation**,
        with every statistic a field of the layout and computed from the
        observation (``encoding.md`` §4). The flat layout keeps ``sbi``'s
        z-scoring, where it is exactly right: one row of real values.
        """
        prior = _prior_class(torch)(
            self.problem,
            self.stream("prior"),
            dtype=getattr(torch, _DTYPE),
            device=self.device,
        )
        network, architecture = self._network(embedding=embedding, layout=layout)
        # NRE spells its network `classifier=`; NPE and NLE spell it
        # `density_estimator=`. One keyword table rather than one branch per
        # method, so a fourth method is a row.
        keyword = "classifier" if self.method in _RATIO_METHODS else "density_estimator"
        trainer = getattr(sbi_package.inference, METHODS[self.method])(
            prior=prior,
            device=self.device,
            show_progress_bars=progress,
            **{keyword: network},
        )
        return trainer, architecture

    def _network(self, *, embedding: _Embedding, layout: EncodingLayout) -> tuple[Any, str]:
        """The estimator ``sbi`` is handed, and the name a run records for it.

        Split out of :meth:`_trainer` at **W3.4** because a TMNRE run builds
        several: one joint estimator and one per marginal, each needing its own
        network **and its own embedding instance** — two trainers sharing one
        ``torch.nn.Module`` would train one net against two objectives.
        """
        spec = self.density_estimator_spec
        if spec is not None and not isinstance(spec, str):
            return spec, getattr(spec, "__name__", type(spec).__name__)
        architecture = str(spec or _DEFAULT_ARCHITECTURE[self.method])
        if embedding.module is None:
            return architecture, architecture
        builder_name, keyword = _BUILDERS[self.method]
        from sbi import neural_nets  # pyrefly: ignore[missing-import]

        builder = getattr(neural_nets, builder_name)
        options: dict[str, Any] = {keyword: embedding.module}
        if layout.kind == SET_KIND:
            options["z_score_x"] = "none"
        return builder(model=architecture, **options), architecture

    # -- W3.4: the truncated marginal loop -------------------------------------

    def _run_tmnre(
        self,
        sbi_package: Any,
        torch: Any,
        *,
        trainer: Any,
        layout: EncodingLayout,
        features: int,
        dtype: Any,
        observation: Any,
        training: Mapping[str, Any] | None,
        progress: bool,
    ) -> tuple[int, int]:
        """TMNRE's rounds: marginals, a box, a truncated prior, a joint fit.

        The loop the item's design paragraph describes, with three decisions
        made here and worth reading before the code.

        **The 1-D estimators are trained every round; the pairs are not.** Only
        the 1-D marginals build the truncation box, so training a pair
        estimator before the last round would be work no decision consumes. The
        pairs still *accumulate* every round's rows, so the estimator they
        finally train on is the whole budget, not the last round's slice.

        **The joint estimator accumulates every round and trains once, last.**
        That is what gives the run ordinary i.i.d. joint draws to emit, scored
        on the numpy path exactly as an NPE run's are, and it is sound over
        mixed rounds for the same reason truncation needs no importance
        correction: a ratio classifier's θ-dependence is ``p(x|θ)`` whatever
        the θ were drawn from, the proposal entering only as a constant this
        method never needs. So the final estimator sees the union of the
        rounds, and is multiplied by the *truncated* prior to give a posterior.

        **Every round simulates through ``values=``.** The proposal is a
        truncated *prior*, not a trained posterior, so θ is drawn from the
        :class:`~sbi.utils.RestrictedPrior`, mapped back through ``constrain``,
        and handed to
        :meth:`~ampere.core.dataset.FittingProblem.simulate_many` as given
        values (``inference.md`` §13) — which is also what keeps every round on
        its own simulation sub-stream.
        """
        problem = self.problem
        free = int(problem.free_size)
        options = dict(training or {})
        singles = marginal_indices(free, 1)
        pairs = marginal_indices(free, 2) if self.marginals == 2 else ()
        self.marginal_estimators = {
            indices: MarginalEstimator(
                indices,
                self._marginal_trainer(
                    sbi_package, torch, layout=layout, features=features, progress=progress
                ),
            )
            for indices in (*singles, *pairs)
        }
        box = TruncationBox.unbounded(free)
        proposal: Any = None
        accumulated: list[np.ndarray] = []
        simulated = 0
        usable = 0
        self.truncation_history = []
        for round_index in range(self.rounds):
            values, acceptance = self._propose(proposal)
            theta, summary, counts = self._simulate_round(round_index, values=values)
            simulated += counts[0]
            usable += counts[1]
            accumulated.append(theta)
            theta_tensor = torch.as_tensor(theta, dtype=dtype, device=self.device)
            summary_tensor = torch.as_tensor(summary, dtype=dtype, device=self.device)
            trainer.append_simulations(theta_tensor, summary_tensor, from_round=round_index)
            final = round_index == self.rounds - 1
            # W3.15: one reseed per round covers every estimator this round
            # trains (one or more 1-D marginals, the pairs on the last round)
            # under a single, round-distinct draw from "sbi.torch" -- the
            # generator behind it advances on every call, so this is already
            # a different seed from the one the previous round used.
            with self._seeded(torch, "sbi.torch"):
                for indices, estimator in self.marginal_estimators.items():
                    estimator.append(theta_tensor, summary_tensor, round_index=round_index)
                    if len(indices) == 1 or final:
                        estimator.train(**options)
            rows = np.concatenate(accumulated, axis=0)
            box = self._truncate(box, rows, observation=observation, torch=torch, dtype=dtype)
            record: dict[str, Any] = {
                "round": round_index + 1,
                "simulations": counts[0],
                "usable": counts[1],
                # The rate at which the proposal that produced *this* round's
                # draws accepted, which is the previous round's box under the
                # untruncated prior -- 1.0 for round 1, which drew from the
                # prior itself.
                "proposal_acceptance": acceptance,
            }
            record.update(box.to_dict(problem.constrain))
            self.truncation_history.append(record)
            proposal = self._restricted_prior(torch, box, round_index)
        self.truncation = box
        self._proposal = proposal
        # Every round's rows, not only the truncated ones. ``sbi``'s
        # ``discard_prior_samples=True`` was measured on this problem and is
        # not an improvement: it costs a third of the budget and leaves the
        # estimate no closer to an emcee reference, because a ratio's
        # θ-dependence is ``p(x|θ)`` whatever the θ were drawn from, so round
        # 1's wide rows are ordinary training data rather than contamination.
        # A caller who wants it can still pass it through ``training=``.
        with self._seeded(torch, "sbi.torch"):
            self.estimator = trainer.train(show_train_summary=False, **options)
        self.posterior = trainer.build_posterior(
            self.estimator,
            prior=proposal,
            sample_with=self.sample_with,
            **({"mcmc_parameters": dict(_TMNRE_MCMC)} if self.sample_with == "mcmc" else {}),
        )
        self.posterior.set_default_x(observation)
        self.marginal_summary = self._marginal_summary(
            box,
            np.concatenate(accumulated, axis=0),
            observation=observation,
            torch=torch,
            dtype=dtype,
        )
        return simulated, usable

    def _marginal_trainer(
        self,
        sbi_package: Any,
        torch: Any,
        *,
        layout: EncodingLayout,
        features: int,
        progress: bool,
    ) -> Any:
        """One ``NRE`` trainer for one marginal, with **no prior at all**.

        ``prior=None`` is not a shortcut: an ``NRE`` classifier is trained by
        contrasting a batch's own ``(θ, x)`` pairs against its own permuted
        ones, so nothing during training asks a distribution for a density, and
        a marginal of ampere's unconstrained joint prior has no closed form to
        give it anyway. The marginal posterior is recovered afterwards as
        ``exp(log r) · q(θ_i)`` with ``q`` estimated from the very θ column the
        estimator trained on (:func:`~ampere.inference._tmnre.
        marginal_log_density`), which is the quantity the ratio is a ratio
        *to*, so the two halves match by construction rather than by assumption.

        The embedding is rebuilt per trainer rather than shared: two trainers
        holding one ``torch.nn.Module`` would optimise one network against two
        objectives. It is built from the same ``embedding=`` and the same
        *feature count* as the joint one, so every estimator in a run reads the
        same summary in the same way; what differs between them is only which
        columns of θ they see.
        """
        rebuilt = _embedding_of(
            self.embedding,
            torch=torch,
            features=features,
            free_size=self.problem.free_size,
            layout=layout,
        )
        network, _ = self._network(embedding=rebuilt, layout=layout)
        return getattr(sbi_package.inference, METHODS[self.method])(
            prior=None,
            device=self.device,
            show_progress_bars=progress,
            classifier=network,
        )

    def _propose(self, proposal: Any) -> tuple[Any, float]:
        """The next round's θ, in **constrained** coordinates, and the rate.

        ``None`` for round 1: ``simulate_many`` draws from the joint prior
        itself, which is the budget idiom and cheaper than rejecting against a
        box that accepts everything. Afterwards the draws come from the
        truncated prior, and the fraction of untruncated prior draws it
        accepted is recorded — it is the box's own prior mass, the honest
        measure of how much a round's budget was concentrated, and the number
        that tells a reader why a later round took longer.
        """
        if proposal is None:
            return None, 1.0
        drawn = proposal.sample((self.budget,), show_progress_bars=False)
        unconstrained = np.asarray(drawn.detach().cpu().numpy(), dtype=float).reshape(
            self.budget, self.problem.free_size
        )
        rate = getattr(proposal, "acceptance_rate", None)
        acceptance = 1.0 if rate is None else float(np.asarray(_numpy(rate)).reshape(-1)[0])
        if acceptance < _TRUNCATION_ACCEPTANCE_FLOOR:
            warnings.warn(
                f"tmnre's truncated prior accepted {acceptance:.3%} of prior draws, so a round's "
                f"budget costs roughly {1.0 / max(acceptance, 1e-12):.0f} prior samples each. The "
                f"box has become very small relative to the prior; lower truncation_epsilon= for "
                f"a wider box, or accept the cost.",
                SamplingFailureWarning,
                stacklevel=3,
            )
        return np.stack([self.problem.constrain(row) for row in unconstrained]), acceptance

    def _truncate(
        self,
        box: TruncationBox,
        rows: np.ndarray,
        *,
        observation: Any,
        torch: Any,
        dtype: Any,
    ) -> TruncationBox:
        """The next box: per parameter, where the 1-D marginal exceeds ``ε·max``.

        The grid spans the intersection of the incoming box with the range of θ
        actually simulated so far, which is what makes the sequence of boxes
        **nested by construction**: no edge can move outward, because no node
        outside the previous box is ever evaluated. :meth:`TruncationBox.
        intersect` then states that as a property rather than leaving it to
        this arithmetic.
        """
        lower: list[float] = []
        upper: list[float] = []
        for index in range(rows.shape[1]):
            column = rows[:, index]
            low = max(float(column.min()), box.lower[index])
            high = min(float(column.max()), box.upper[index])
            grid = grid_between(low, high, GRID_POINTS_1D)
            nodes = grid.reshape(-1, 1)
            estimator = self.marginal_estimators[(index,)]
            density = estimator.log_ratio(
                nodes, observation, torch=torch, dtype=dtype
            ) + marginal_log_density(column.reshape(-1, 1), nodes)
            edges = interval_above(grid, density, self.truncation_epsilon)
            lower.append(edges[0])
            upper.append(edges[1])
        return TruncationBox(tuple(lower), tuple(upper)).intersect(box)

    def _restricted_prior(self, torch: Any, box: TruncationBox, round_index: int) -> Any:
        """The joint prior renormalised on *box*, as ``sbi``'s own object.

        A fresh base prior per round, on its own named sub-stream, for
        ``inference.md`` §12's reason: two rounds sharing a stream would make
        one round's draws depend on how many the other took, and a rejection
        sampler's consumption depends on the box.
        """
        utils = importlib.import_module("sbi.utils")
        base = _prior_class(torch)(
            self.problem,
            self.stream(f"tmnre.propose.{round_index + 2}"),
            dtype=getattr(torch, _DTYPE),
            device=self.device,
        )
        return restricted_prior_class(utils.RestrictedPrior)(
            base, box.indicator(), sample_with="rejection", device=self.device
        )

    def _marginal_summary(
        self,
        box: TruncationBox,
        rows: np.ndarray,
        *,
        observation: Any,
        torch: Any,
        dtype: Any,
    ) -> MarginalSummary:
        """Every trained marginal, evaluated on a grid over the **final** box.

        This is the ``marginals`` group's content and the answer to the results
        question the item names: the estimators are torch modules and have no
        place in a netCDF file, while the curve each one draws over the final
        box is exactly what a corner plot is and what a reader of an archived
        run can act on.

        Evaluated here rather than reused from the round loop because the loop's
        grids span each round's *incoming* box, and the stored summary should
        span the box the run finished with.
        """
        problem = self.problem
        labels = tuple(problem.free_labels())
        free = len(labels)
        centre = 0.5 * (np.asarray(box.lower, dtype=float) + np.asarray(box.upper, dtype=float))
        grid = np.empty((free, GRID_POINTS_1D), dtype=float)
        grid_constrained = np.empty_like(grid)
        log_ratio = np.empty_like(grid)
        log_density = np.empty_like(grid)
        for index in range(free):
            nodes = grid_between(box.lower[index], box.upper[index], GRID_POINTS_1D)
            column = nodes.reshape(-1, 1)
            ratio = self.marginal_estimators[(index,)].log_ratio(
                column, observation, torch=torch, dtype=dtype
            )
            density = ratio + marginal_log_density(rows[:, [index]], column)
            grid[index] = nodes
            grid_constrained[index] = constrained_grid(problem, nodes, index, centre)
            log_ratio[index] = ratio
            log_density[index] = density - np.nanmax(density)

        pairs = marginal_indices(free, 2) if self.marginals == 2 else ()
        pair_row = np.zeros((len(pairs), GRID_POINTS_2D), dtype=float)
        pair_column = np.zeros_like(pair_row)
        pair_ratio = np.zeros((len(pairs), GRID_POINTS_2D, GRID_POINTS_2D), dtype=float)
        pair_density = np.zeros_like(pair_ratio)
        for position, indices in enumerate(pairs):
            first, second = indices
            axis_row = grid_between(box.lower[first], box.upper[first], GRID_POINTS_2D)
            axis_column = grid_between(box.lower[second], box.upper[second], GRID_POINTS_2D)
            mesh = pair_mesh(axis_row, axis_column)
            ratio = self.marginal_estimators[indices].log_ratio(
                mesh, observation, torch=torch, dtype=dtype
            )
            density = ratio + marginal_log_density(rows[:, list(indices)], mesh)
            shape = (GRID_POINTS_2D, GRID_POINTS_2D)
            pair_row[position] = axis_row
            pair_column[position] = axis_column
            pair_ratio[position] = ratio.reshape(shape)
            pair_density[position] = (density - np.nanmax(density)).reshape(shape)

        return MarginalSummary(
            labels=labels,
            order=self.marginals,
            epsilon=self.truncation_epsilon,
            grid=grid,
            grid_constrained=grid_constrained,
            log_ratio=log_ratio,
            log_density=log_density,
            pairs=pairs,
            pair_grid_row=pair_row,
            pair_grid_column=pair_column,
            pair_log_ratio=pair_ratio,
            pair_log_density=pair_density,
        )

    def _artefact(self) -> Any:
        """What the artefact store is handed: the posterior, or TMNRE's bundle.

        A TMNRE run's answer is not only its posterior — the truncation history
        and the marginal summary are part of what the run emits, and a cache
        hit that produced a run missing its ``marginals`` group would be a
        different run under the same key. So the bundle is stored whole
        (:class:`~ampere.inference._tmnre.TMNREArtefact`), and everything in it
        pickles: the box is plain floats, the indicator a dataclass rather than
        a closure, the summary numpy arrays.
        """
        if self.method != TMNRE:
            return self.posterior
        return TMNREArtefact(
            posterior=self.posterior,
            truncation=self.truncation,
            history=list(self.truncation_history),
            marginals=self.marginal_summary,
        )

    def _restore(self, cached: Any, observation: Any) -> None:
        """A cache hit, put back where a run would have left it.

        The *estimators* are not restored — a stored run holds their summary,
        not their weights — so :attr:`marginal_estimators` stays empty after a
        hit. Everything the emitted run carries is restored.
        """
        if isinstance(cached, TMNREArtefact):
            self.posterior = cached.posterior
            self.truncation = cached.truncation
            self.truncation_history = list(cached.history)
            self.marginal_summary = cached.marginals
        else:
            self.posterior = cached
        self.posterior.set_default_x(observation)

    def _simulate_round(
        self, round_index: int, *, proposal: Any = None, values: Any = None, torch: Any = None
    ) -> tuple[np.ndarray, np.ndarray, tuple[int, int]]:
        """One round's budget: θ in unconstrained space, x as a summary matrix.

        Round 1 draws θ from the joint prior (``values=None`` -- the budget
        idiom). Every round after it draws θ from the previous round's
        posterior at the observed data and simulates *there*, which is what
        makes a multi-round run concentrate its budget and what makes it no
        longer amortised.

        The chunk iterator is consumed here rather than materialised, so a
        budget larger than memory is held only as the two small arrays a
        trainer needs -- θ and the summary -- while the containers stream to
        the training-set file and are dropped.

        **W3.4**: *values* is the other way in — one **constrained** θ per
        draw, already drawn elsewhere. TMNRE's rounds arrive that way, because
        its proposal is a truncated *prior* rather than a trained posterior, so
        the draws come from a ``RestrictedPrior`` and reach ``simulate_many``
        through its ``values=`` argument (``inference.md`` §13). *proposal* and
        *values* are alternatives; passing both is a caller error.

        **W3.15**: *torch* is required exactly when *proposal* is given --
        drawing from it is a real, trained posterior's ``sample()`` and needs
        the same reseed the run's final draw does, for round 2 onwards to
        repeat from the problem's seed too.
        """
        problem = self.problem
        if proposal is not None and values is not None:  # pragma: no cover - internal
            raise EngineError("a simulation round takes a proposal or given values, not both.")
        if proposal is not None:
            # W3.15: *proposal* here is a trained posterior (an NPE flow, or
            # an NLE/NRE's MCMC posterior), so drawing from it spends torch's
            # global generator exactly as the final posterior draw does, and
            # needs the same reseed to make round 2's simulated batch (and
            # therefore everything trained on it) repeat from the problem's
            # seed. TMNRE never reaches here with a proposal (its rounds
            # arrive through *values*, from the numpy-seeded prior), so
            # *torch* is only required when *proposal* is given.
            with self._seeded(torch, "sbi.torch.sample"):
                drawn = proposal.sample((self.budget,), show_progress_bars=False)
            unconstrained = np.asarray(drawn.detach().cpu().numpy(), dtype=float).reshape(
                self.budget, problem.free_size
            )
            values = np.stack([problem.constrain(row) for row in unconstrained])

        thetas: list[np.ndarray] = []
        summaries: list[np.ndarray] = []
        simulated = 0
        last: SimulationBatch | None = None
        for chunk in self._chunks(values=values, round_index=round_index):
            simulated += len(chunk)
            last = chunk
            if self.training_set is not None:
                self._write(chunk)
            keep = chunk.usable
            if not len(keep):
                continue
            thetas.append(np.stack([problem.unconstrain(draw.theta) for draw in keep]))
            summaries.append(
                _summary_of(keep.observations, problem.datasets, batched=True, layout=self.encoding)
            )
        self.batch = last
        if not thetas:
            raise EngineError(
                f"sbi has nothing to train on: every one of the {simulated} simulation(s) in "
                f"round {round_index + 1} failed. problem.failure_summary() says why:\n"
                f"{problem.failure_summary()}"
            )
        theta = np.concatenate(thetas, axis=0)
        summary = np.concatenate(summaries, axis=0)
        return theta, summary, (simulated, int(theta.shape[0]))

    def _chunks(self, *, values: Any, round_index: int) -> Iterator[SimulationBatch]:
        """The budget, as chunks, on a per-round sub-stream.

        A distinct stream per round (``"sbi.simulate.1"``, ``".2"``, ...) for
        ``inference.md`` §12's reason: two rounds sharing a stream would make
        round 2's draws depend on round 1's budget, so changing the budget
        would silently change what round 2 simulated.
        """
        return self.problem.simulate_many(
            self.budget,
            values=values,
            observe=True,
            rng=self.stream(f"simulate.{round_index + 1}"),
            executor=self.executor,
            chunk_size=self.chunk_size,
            as_chunks=True,
        )

    def _write(self, chunk: SimulationBatch) -> None:
        """One chunk into the training-set file: create it, then append.

        ``results.md`` limitation 13.9: the append is read-concatenate-rewrite,
        ``O(existing + new)``, so few large chunks beat many small ones. The
        cost is the user's to choose through ``chunk_size=``; what this method
        guarantees is that the whole budget is never held in memory to pay it.

        The layout's name and hash ride along in the chunk's **provenance**,
        which :func:`~ampere.results.training.write_training_set` writes to the
        file's root attributes. That is what makes a stored budget answer the
        question W3.5 and W3.6 both ask of it -- *which packing were these rows
        encoded under?* -- without either of them having to guess, and it needs
        no new argument on a writer this item does not own.
        """
        path = self.training_set
        if path is None:  # pragma: no cover - the caller checks before calling
            return
        stamped = dataclasses.replace(chunk, provenance=self._encoding_provenance(chunk))
        if self._written:
            append_training_set(path, stamped, self.problem)
        else:
            write_training_set(path, stamped, self.problem)
            self._written = True

    def _encoding_provenance(self, chunk: SimulationBatch) -> dict[str, Any]:
        """The chunk's own provenance plus this run's layout name and hash."""
        found = dict(chunk.provenance)
        if self.encoding is not None:
            found["encoding_layout"] = self.encoding.to_dict()
            found["encoding_hash"] = self.encoding.hash
        return found

    def _append(
        self,
        trainer: Any,
        theta: np.ndarray,
        summary: np.ndarray,
        *,
        torch: Any,
        dtype: Any,
        proposal: Any,
        round_index: int,
    ) -> None:
        """Hand one round's pairs to the trainer, in the spelling it wants.

        NPE's ``append_simulations`` takes ``proposal=`` and needs it to know
        that a later round's θ did not come from the prior (which is what its
        atomic loss corrects for); NLE and NRE take ``from_round`` instead, and
        pass no proposal at all. The difference is ``sbi``'s, and it is the one
        place this driver has to know which family it is driving.

        **W3.4**: TMNRE passes the *actual* round index rather than
        ``0``-or-``1``, because its rounds accumulate — a ratio trained on the
        union of every round's rows is still a ratio for whatever prior it is
        later multiplied by, which is precisely why truncation needs no
        importance correction (:mod:`ampere.inference._tmnre`) — and the index
        is what ``sbi``'s own per-round bookkeeping then reports.
        """
        theta_tensor = torch.as_tensor(theta, dtype=dtype, device=self.device)
        summary_tensor = torch.as_tensor(summary, dtype=dtype, device=self.device)
        if self.method == "npe":
            trainer.append_simulations(theta_tensor, summary_tensor, proposal=proposal)
        elif self.method == TMNRE:
            trainer.append_simulations(theta_tensor, summary_tensor, from_round=int(round_index))
        else:
            trainer.append_simulations(
                theta_tensor, summary_tensor, from_round=0 if proposal is None else 1
            )

    def _estimator_log_prob(self, drawn: Any, *, torch: Any) -> np.ndarray:
        """The network's own log-density at each stored draw.

        ``potential`` rather than ``log_prob``, for all three methods, because
        it is the one surface every ``sbi`` posterior has and because
        ``log_prob`` is deprecated on the two that can only answer up to the
        evidence. For NPE the potential *is* the normalised log-posterior --
        there is no leakage to correct for, since the prior's support is all of
        ℝⁿ -- and ``ampere_sbi_log_prob_kind`` records which of the two a run
        has.

        A posterior that cannot evaluate itself gives NaNs rather than failing
        the run: the fit succeeded, and the estimator's own density is a
        convenience beside the true one, not the run's reason for existing.
        """
        try:
            with torch.no_grad():
                values = self.posterior.potential(drawn, track_gradients=False)
            return np.asarray(values.detach().cpu().numpy(), dtype=float).reshape(-1)
        except Exception:
            return np.full(int(drawn.shape[0]), np.nan, dtype=float)

    def _attrs(
        self,
        sbi_package: Any,
        torch: Any,
        *,
        trainer: Any,
        architecture: str,
        embedding: _Embedding,
        layout: EncodingLayout,
        features: int,
        simulated: int,
        usable: int,
        draws: int,
        torch_seed: int | None,
    ) -> dict[str, object]:
        """Everything a reader of an archived SBI run needs to judge it."""
        summary = getattr(trainer, "summary", {}) or {}
        training_loss = [float(value) for value in np.asarray(summary.get("training_loss", []))]
        validation_loss = [float(value) for value in np.asarray(summary.get("validation_loss", []))]
        thinned, stride = _thinned(training_loss)
        attrs: dict[str, object] = {
            "sbi_method": self.method,
            "sbi_trainer": type(trainer).__name__,
            "sbi_density_estimator": architecture,
            "sbi_budget": self.budget,
            "sbi_rounds": self.rounds,
            "sbi_simulations": simulated,
            "sbi_usable_simulations": usable,
            "sbi_failures": simulated - usable,
            "sbi_draws": draws,
            "sbi_embedding": embedding.name,
            "sbi_embedding_output_dim": embedding.output_dim,
            "sbi_summary_layout": layout.name,
            "sbi_summary_features": features,
            # The layout, in full, and its hash. A network is only meaningful
            # about a tensor packed the way it was trained on, so a run that
            # cannot say which packing that was cannot be reused -- by W3.5's
            # artefact cache, by W3.6's coverage batches, or by a reader.
            "encoding_layout": layout.to_dict(),
            "encoding_hash": layout.hash,
            "sbi_encoding_rows": layout.row_cap,
            "sbi_encoding_columns": layout.columns_total,
            "sbi_parameterisation": "unconstrained",
            # Whether the trained estimator is still valid at *another*
            # observation. A single-round fit is; every multi-round one is not,
            # because its proposal (a trained posterior, or W3.4's truncation
            # box) was chosen at this observation. Recorded rather than left to
            # be inferred from the round count, because TMNRE loses it at one
            # round too -- the box is applied to the posterior either way.
            "sbi_amortised": int(self.rounds == 1 and self.method != TMNRE),
            "sbi_context": "none",
            "sbi_device": self.device,
            "sbi_training_loss": thinned,
            "sbi_training_loss_stride": stride,
            "sbi_validation_loss": _thinned(validation_loss)[0],
            "sbi_epochs_trained": _last_int(summary.get("epochs_trained")),
            "sbi_log_prob_kind": "normalised" if self.method == "npe" else "unnormalised",
            # str() rather than the objects themselves, and it is load-bearing:
            # ``torch.__version__`` is a ``TorchVersion``, a str *subclass*, and
            # h5netcdf's attribute writer takes ``numpy.asarray(value).dtype``
            # for anything without a ``dtype`` — which for a str subclass is
            # ``<U10`` rather than h5py's variable-length string type, and h5py
            # then refuses it with "No conversion path for dtype". The run
            # would build and fail only at ``to_netcdf``, hours later.
            "sbi_version": str(getattr(sbi_package, "__version__", "unknown")),
            "torch_version": str(torch.__version__),
        }
        # W3.15: absent, not `None`, for an unseeded problem -- an unseeded
        # run asked for fresh randomness every time, and a recorded `None`
        # would read as "recorded and empty" rather than "not applicable"
        # (the same convention `calibration_seed` already uses, above).
        if torch_seed is not None:
            attrs["sbi_torch_seed"] = torch_seed
        if self.method == TMNRE:
            attrs["sbi_marginals"] = self.marginals
            attrs["sbi_truncation_epsilon"] = self.truncation_epsilon
            attrs["sbi_truncation_sampler"] = str(self.sample_with)
            # One record per round: the box after it, in both parameterisations,
            # its log-volume, the round's simulation counts and the rate at
            # which the box that produced them accepted prior draws. This is
            # the run's own evidence that the truncation behaved -- the volumes
            # must not grow, and the truth (when a study knows one) must stay
            # inside every box.
            attrs["sbi_truncation"] = list(self.truncation_history)
            attrs["sbi_marginal_estimators"] = len(marginal_indices(self.problem.free_size, 1)) + (
                len(marginal_indices(self.problem.free_size, 2)) if self.marginals == 2 else 0
            )
            if self.truncation is not None:
                attrs["sbi_truncation_log_volume"] = self.truncation.log_volume
        if training_loss:
            attrs["sbi_final_training_loss"] = training_loss[-1]
        if validation_loss:
            attrs["sbi_best_validation_loss"] = min(validation_loss)
        if self.training_set is not None:
            attrs["sbi_training_set"] = self.training_set
        if self.executor is not None:
            attrs["sbi_executor"] = type(self.executor).__name__
        if self.chunk_size is not None:
            attrs["sbi_chunk_size"] = self.chunk_size
        return attrs


# ---------------------------------------------------------------------------
# Plumbing
# ---------------------------------------------------------------------------


def _require_sbi() -> tuple[Any, Any]:
    """``sbi`` and torch, imported on use and refused by name.

    ``architecture.md`` §4 rule 3: on use, never on import. Both are behind the
    one ``sbi`` extra (installing it installs torch), so one refusal names one
    extra whichever of the two is missing.
    """
    try:
        import sbi  # pyrefly: ignore[missing-import]
        import sbi.inference  # pyrefly: ignore[missing-import]
        import torch  # pyrefly: ignore[missing-import]
    except ImportError as error:
        missing = getattr(error, "name", None) or "sbi"
        raise OptionalDependencyError(
            missing,
            extra="sbi",
            context="simulation-based inference with SBIEngine",
        ) from error
    return sbi, torch


def _with_estimator_log_prob(tree: Any, values: np.ndarray) -> Any:
    """Put the estimator's own per-draw log-density beside the true one.

    ``sample_stats`` is where ``results.md`` §4 puts a ``(chain, draw)``-shaped
    per-draw quantity, and this is one. Added to the emitted tree rather than
    passed through :func:`ampere.results.emit`, which has no hook for an extra
    sample statistic; a hook there would be a change to a §4 contract, and this
    driver does not need one to put the number where a reader will look for it.

    Nothing is written if the shapes disagree, which they only can if the
    posterior returned a different number of draws than it was asked for --
    in that case the run is still emitted, without the extra column, rather
    than lost.
    """
    try:
        node = tree["sample_stats"]
        dataset = node.dataset
        shape = tuple(dataset["lp"].shape)
        if int(np.prod(shape)) != int(values.size):
            return tree
        node.dataset = dataset.assign(
            {"ampere_sbi_log_prob": (dataset["lp"].dims, values.reshape(shape))}
        )
    except Exception:
        return tree
    return tree


def _numpy(value: Any) -> np.ndarray:
    """A torch tensor as a plain float array, detached and off the device."""
    detach = getattr(value, "detach", None)
    tensor = value if detach is None else detach().cpu()
    return np.asarray(tensor.numpy() if hasattr(tensor, "numpy") else tensor, dtype=float)


def _thinned(values: Sequence[float]) -> tuple[list[float], int]:
    """*values* reduced to at most :data:`_TRACE_POINTS`, with the stride used.

    A real stride rather than evenly spaced indices, because the stride is what
    the attrs record and a reader must be able to reconstruct which epoch each
    recorded value belongs to. The last value is always kept, so a trace that
    is still descending at its end says so.
    """
    kept = list(values)
    if len(kept) <= _TRACE_POINTS:
        return [float(value) for value in kept], 1
    stride = math.ceil(len(kept) / _TRACE_POINTS)
    thinned = [float(value) for value in kept[::stride]]
    if thinned and thinned[-1] != float(kept[-1]):
        thinned.append(float(kept[-1]))
    return thinned, stride


def _last_int(value: Any) -> int:
    """``sbi``'s summary keeps one entry per training call; the last is this run's."""
    if value is None:
        return 0
    array = np.asarray(value).reshape(-1)
    return int(array[-1]) if array.size else 0
