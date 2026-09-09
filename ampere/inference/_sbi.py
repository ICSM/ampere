"""The simulation-based driver: NPE, NLE and NRE over ``simulate_many``.

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
"""

from __future__ import annotations

import dataclasses
import functools
import importlib
import math
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

from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError

__all__ = [
    "EMBEDDINGS",
    "LAYOUTS",
    "METHODS",
    "SET_EMBEDDINGS",
    "SUMMARY_LAYOUT",
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
METHODS: Mapping[str, str] = {"npe": "NPE", "nle": "NLE", "nre": "NRE"}

#: The default network architecture per method. NPE and NLE learn a *density*
#: and take a normalising flow; NRE learns a *classifier* and takes a residual
#: network, which is why the two are not one default.
_DEFAULT_ARCHITECTURE: Mapping[str, str] = {"npe": "maf", "nle": "maf", "nre": "resnet"}

#: Which ``sbi.neural_nets`` builder wraps an embedding net for each method,
#: and the keyword each one spells the embedding with. NRE's classifier sees θ
#: and x separately and so has two slots; the embedding is the *data* one.
_BUILDERS: Mapping[str, tuple[str, str]] = {
    "npe": ("posterior_nn", "embedding_net"),
    "nle": ("likelihood_nn", "embedding_net"),
    "nre": ("classifier_nn", "embedding_net_x"),
}

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
      what a masked row holds whatever the net does with the mask. The mask is
      passed as well, so that a later ``sbi`` honouring it changes nothing here.
    * The aggregation is ``hidden_states[:, -1, :]``: the **last** token, after
      full attention. A padded row is a zero token whose attention over the
      sequence is uniform, so the summary is still a function of every row --
      but it is a last-token read rather than a pooled one, and a caller wanting
      pooling wants a module of their own over :func:`~ampere.core.encoding.unpack`.

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


@functools.cache
def _wrapper_classes(torch: Any) -> tuple[Any, Any]:
    """The two ``nn.Module`` wrappers, defined once per interpreter.

    Inside a function for :func:`_prior_class`'s reason: they subclass
    ``torch.nn.Module``, and this module imports torch on use, never on import.

    Both do the same three things and differ only in what they hand the net.
    They unpack ``x`` by the layout (which also reshapes it, since ``sbi``
    flattens ``x`` on some paths), they narrow float64 to the network's float32
    **once, here, at the boundary**, and they convert the packing's single mask
    column into whatever their net wants -- NaN rows for the set embedding, an
    attention mask and zeroed tokens for the transformer. No embedding ever sees
    the mask column as a feature: :attr:`~ampere.core.encoding.Unpacked.features`
    excludes it by construction.
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
        """Zeroed tokens plus an attention mask, then ``sbi``'s transformer."""

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
            output = self.net(tokens, attention_mask=keep.to(tokens.dtype))
            return output[0] if isinstance(output, tuple) else output

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
        likelihood, sampled with MCMC) or ``"nre"`` (the likelihood-to-evidence
        ratio, also sampled with MCMC). See :data:`METHODS`.
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
                f"are then sampled with MCMC."
            )
        if int(budget) < 1:
            raise EngineError(f"sbi needs a simulation budget of at least 1, got {budget}.")
        if int(rounds) < 1:
            raise EngineError(f"sbi needs at least one round, got {rounds}.")
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
            no defaults: they are ``sbi``'s and they are good ones.
        posterior_options
            Forwarded verbatim to the trained posterior's ``sample``. This is
            where an MCMC posterior's ``num_chains``, ``warmup_steps``,
            ``thin`` and ``init_strategy`` go, which is what makes an NLE or
            NRE run affordable at a small draw count.
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
            )
        )
        cached = None if cache_key is None or self.cache is None else self.cache.get(cache_key)
        cache_hit = cached is not None
        if cache_hit:
            self.posterior = cached
            self.posterior.set_default_x(observation_tensor)
        else:
            for round_index in range(self.rounds):
                theta, summary, counts = self._simulate_round(round_index, proposal=proposal)
                simulated += counts[0]
                usable += counts[1]
                self._append(trainer, theta, summary, torch=torch, dtype=dtype, proposal=proposal)
                self.estimator = trainer.train(show_train_summary=False, **dict(training or {}))
                self.posterior = trainer.build_posterior(self.estimator)
                self.posterior.set_default_x(observation_tensor)
                proposal = self.posterior
            if cache_key is not None and self.cache is not None:
                self.cache.put(cache_key, self.posterior)

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
        )
        if cache_key is not None:
            attrs["sbi_cache_hit"] = int(cache_hit)
            attrs["sbi_cache_key"] = cache_key.digest()
        tree = self.finish(chain, extra_attrs=attrs)
        return _with_estimator_log_prob(tree, estimator_log_prob)

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
        thetas, summaries, simulated = self._calibration_batch(simulations, rng)
        theta_tensor = torch.as_tensor(thetas, dtype=dtype, device=self.device)
        summary_tensor = torch.as_tensor(summaries, dtype=dtype, device=self.device)
        prior = np.stack(
            [problem.unconstrain(problem.sample_prior(rng)) for _ in range(thetas.shape[0])]
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

    def _calibration_batch(self, count: int, rng: Any) -> tuple[np.ndarray, np.ndarray, int]:
        """A fresh prior batch, unconstrained θ and layout-encoded x.

        Streamed in chunks exactly as :meth:`_simulate_round` streams a
        training round, and on its own named sub-stream (``"sbi.calibrate"``),
        because ``lowering.md`` §9.2's rule is that adding a diagnostic must not
        change what a fit simulated.
        """
        problem = self.problem
        thetas: list[np.ndarray] = []
        summaries: list[np.ndarray] = []
        simulated = 0
        for chunk in problem.simulate_many(
            count,
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
        spec = self.density_estimator_spec
        if spec is not None and not isinstance(spec, str):
            network: Any = spec
            architecture = getattr(spec, "__name__", type(spec).__name__)
        else:
            architecture = str(spec or _DEFAULT_ARCHITECTURE[self.method])
            if embedding.module is None:
                network = architecture
            else:
                builder_name, keyword = _BUILDERS[self.method]
                from sbi import neural_nets  # pyrefly: ignore[missing-import]

                builder = getattr(neural_nets, builder_name)
                options: dict[str, Any] = {keyword: embedding.module}
                if layout.kind == SET_KIND:
                    options["z_score_x"] = "none"
                network = builder(model=architecture, **options)
        # NRE spells its network `classifier=`; NPE and NLE spell it
        # `density_estimator=`. One keyword table rather than one branch per
        # method, so a fourth method is a row.
        keyword = "classifier" if self.method == "nre" else "density_estimator"
        trainer = getattr(sbi_package.inference, METHODS[self.method])(
            prior=prior,
            device=self.device,
            show_progress_bars=progress,
            **{keyword: network},
        )
        return trainer, architecture

    def _simulate_round(
        self, round_index: int, *, proposal: Any
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
        """
        problem = self.problem
        values = None
        if proposal is not None:
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
    ) -> None:
        """Hand one round's pairs to the trainer, in the spelling it wants.

        NPE's ``append_simulations`` takes ``proposal=`` and needs it to know
        that a later round's θ did not come from the prior (which is what its
        atomic loss corrects for); NLE and NRE take ``from_round`` instead, and
        pass no proposal at all. The difference is ``sbi``'s, and it is the one
        place this driver has to know which family it is driving.
        """
        theta_tensor = torch.as_tensor(theta, dtype=dtype, device=self.device)
        summary_tensor = torch.as_tensor(summary, dtype=dtype, device=self.device)
        if self.method == "npe":
            trainer.append_simulations(theta_tensor, summary_tensor, proposal=proposal)
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
