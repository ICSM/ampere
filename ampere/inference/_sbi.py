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

The summary, and why it is behind one function
-----------------------------------------------
The network sees a fixed-length vector, not a
:class:`~ampere.core.dataset.DatasetCollection`. Until W3.3 lands the
coordinate-value-mask encoding, the vector is :func:`_summary_of`: each
dataset's observed values, flattened, masked samples dropped, concatenated in
``datasets`` order. That is enough for a single fitting problem, where the data
layout is fixed anyway, and it is deliberately the *whole* of the layout logic
so that W3.3 replaces one function rather than searching the module. The run
records which layout it used (``ampere_sbi_summary_layout``, ``"flat"`` today),
because a network trained on one layout must never be believed about another.

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
import math
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import OptionalDependencyError
from ampere.core.simulate import Executor, SimulationBatch
from ampere.results.training import append_training_set, write_training_set

from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError

__all__ = [
    "EMBEDDINGS",
    "METHODS",
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
}

#: What :func:`_summary_of` builds, recorded in ``ampere_sbi_summary_layout``.
#: W3.3 adds the coordinate-value-mask layouts beside it and makes this one of
#: several; a run that does not say which layout it used cannot be reused.
SUMMARY_LAYOUT = "flat"

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


def _summary_of(
    observations: Mapping[str, Any] | None,
    datasets: Mapping[str, Any],
    *,
    batched: bool,
) -> np.ndarray:
    """The fixed-layout summary vector: ``layout="flat"``.

    Each dataset contributes its observed values, flattened, with the samples
    its :meth:`~ampere.core.dataset.Dataset.effective_mask` excludes dropped;
    the pieces are concatenated in ``datasets`` order. Masked samples are
    dropped **consistently** — the same columns are absent from every draw and
    from the observation the posterior is conditioned on — which is the whole
    of what "fixed layout" means and the reason the mask is read from the
    dataset rather than from the drawn container.

    Parameters
    ----------
    observations
        Label to container. With *batched* true these are
        :class:`~ampere.core.simulate.ContainerBatch`\\ es, whose leading axis
        is the sample axis; with it false they are single containers.
    datasets
        The problem's own collection, which supplies both the order and the
        masks. Reading the order from here rather than from *observations* is
        deliberate: a dict built elsewhere could iterate differently and the
        network would silently see its columns permuted.
    batched
        Whether *observations* carries a leading sample axis.

    Returns
    -------
    numpy.ndarray
        ``(count, features)``, always two-dimensional — ``count`` is 1 for the
        unbatched form, so the observation and the training set have the same
        shape and the same code path.
    """
    if not datasets:
        raise EngineError("sbi has nothing to summarise: this problem declares no datasets.")
    columns: list[np.ndarray] = []
    for label in datasets:
        holding = None if observations is None else observations.get(label)
        if holding is None:
            raise EngineError(
                f"sbi cannot build a summary vector for dataset {label!r}: the simulation "
                f"produced no observation for it. simulate_many(..., observe=True) is what "
                f"draws them, and a dataset whose likelihood family cannot sample refuses by "
                f"name (inference.md §13)."
            )
        values = np.asarray(holding.values)
        if np.iscomplexobj(values):
            raise EngineError(
                f"dataset {label!r} holds complex values, and the {SUMMARY_LAYOUT!r} summary "
                f"layout has no encoding for one — flattening a complex array into a real "
                f"feature vector is a choice (real and imaginary parts as two columns) that "
                f"belongs in the encoding contract rather than in this driver. W3.3's "
                f"coordinate-value-mask encoding defines it; until then, fit a complex "
                f"dataset with a likelihood-based engine."
            )
        flat = np.asarray(values, dtype=float)
        flat = flat.reshape(flat.shape[0], -1) if batched else flat.reshape(1, -1)
        mask = datasets[label].effective_mask
        if mask is not None:
            keep = ~np.asarray(mask, dtype=bool).reshape(-1)
            if keep.size != flat.shape[1]:
                raise EngineError(
                    f"dataset {label!r}: its effective mask covers {keep.size} sample(s) but the "
                    f"container being summarised has {flat.shape[1]}. The mask is resolved at "
                    f"composition and is evaluation-invariant by declaration, so this is a "
                    f"disagreement between the observed container and the simulated one."
                )
            flat = flat[:, keep]
        columns.append(flat)
    summary = np.concatenate(columns, axis=1)
    if summary.shape[1] == 0:
        raise EngineError(
            "sbi's summary vector is empty: every sample of every dataset is masked out, so "
            "there is nothing for the density estimator to condition on."
        )
    return summary


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


# ---------------------------------------------------------------------------
# The embedding vocabulary (legacy parity)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _Embedding:
    """One resolved embedding: the module, what to call it, and its width."""

    module: Any
    name: str
    output_dim: int


def _embedding_of(embedding: Any, *, torch: Any, features: int, free_size: int) -> _Embedding:
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
    width.

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
    """
    if embedding is None or embedding is False:
        return _Embedding(module=None, name="none", output_dim=0)

    default_width = max(2 * int(free_size), 1)

    if embedding is True:
        return _built(
            EMBEDDINGS["default"], {}, torch=torch, features=features, width=default_width
        )
    if isinstance(embedding, str):
        kind = EMBEDDINGS.get(embedding)
        if kind is None:
            known = ", ".join(sorted(EMBEDDINGS))
            raise EngineError(
                f"sbi does not know the embedding {embedding!r}. Available: {known}, a "
                f"torch.nn.Module of your own, or a dict of hyperparameters with a 'type' key."
            )
        return _built(kind, {}, torch=torch, features=features, width=default_width)
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
        width = int(settings.pop("output_dim", default_width))
        return _built(kind, settings, torch=torch, features=features, width=width)
    if isinstance(embedding, torch.nn.Module):
        return _Embedding(
            module=embedding,
            name="custom",
            output_dim=_probe_width(embedding, torch=torch, features=features),
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


def _built(
    kind: str, settings: Mapping[str, Any], *, torch: Any, features: int, width: int
) -> _Embedding:
    """One of the two shipped nets, with its legacy hyperparameters translated."""
    from sbi.neural_nets import embedding_nets  # pyrefly: ignore[missing-import]

    table = _CNN_KEYS if kind == "CNN" else _FC_KEYS
    unknown = sorted(set(settings) - set(table))
    if unknown:
        raise EngineError(
            f"sbi does not know the {kind} embedding hyperparameter(s) {unknown}. "
            f"Available: {sorted(table)}, plus 'type' and 'output_dim'."
        )
    kwargs = {table[key]: value for key, value in settings.items()}
    if kind == "CNN":
        module = embedding_nets.CNNEmbedding(
            input_shape=(int(features),), output_dim=int(width), **kwargs
        )
    else:
        module = embedding_nets.FCEmbedding(
            input_dim=int(features), output_dim=int(width), **kwargs
        )
    return _Embedding(module=module, name=kind, output_dim=int(width))


def _probe_width(module: Any, *, torch: Any, features: int) -> int:
    """A user module's output width, measured rather than asked for.

    A ``torch.nn.Module`` does not have to declare what it produces, and the
    number is worth recording — it is what a reader of an archived run needs to
    know what the network actually saw. So it is measured with one forward pass
    on a zero input, and a module that will not take this problem's summary
    reports 0 rather than failing the run here: the real failure, if there is
    one, belongs at training time with ``sbi``'s own message.
    """
    try:
        with torch.no_grad():
            probe = module(torch.zeros(1, int(features)))
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
        summaries this layout produces.
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
        density_estimator: Any = None,
        rounds: int = 1,
        device: str = "cpu",
        executor: Executor | None = None,
        chunk_size: int | None = None,
        context: Any = None,
        training_set: str | Path | None = None,
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
        self.density_estimator_spec = density_estimator
        self.rounds = int(rounds)
        self.device = str(device)
        self.executor = executor
        self.chunk_size = None if chunk_size is None else int(chunk_size)
        self.training_set = None if training_set is None else str(training_set)
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
        observation = _summary_of(
            {label: problem.datasets[label].observed for label in problem.datasets},
            problem.datasets,
            batched=False,
        )
        features = int(observation.shape[1])

        resolved = _embedding_of(
            self.embedding, torch=torch, features=features, free_size=problem.free_size
        )
        trainer, architecture = self._trainer(
            sbi_package, torch, embedding=resolved, progress=progress
        )
        self.sampler = trainer

        simulated = 0
        usable = 0
        proposal: Any = None
        observation_tensor = torch.as_tensor(observation[0], dtype=dtype, device=self.device)
        for round_index in range(self.rounds):
            theta, summary, counts = self._simulate_round(round_index, proposal=proposal)
            simulated += counts[0]
            usable += counts[1]
            self._append(trainer, theta, summary, torch=torch, dtype=dtype, proposal=proposal)
            self.estimator = trainer.train(show_train_summary=False, **dict(training or {}))
            self.posterior = trainer.build_posterior(self.estimator)
            self.posterior.set_default_x(observation_tensor)
            proposal = self.posterior

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
            features=features,
            simulated=simulated,
            usable=usable,
            draws=int(draws),
        )
        tree = self.finish(chain, extra_attrs=attrs)
        return _with_estimator_log_prob(tree, estimator_log_prob)

    # -- the pieces -----------------------------------------------------------

    def _trainer(
        self, sbi_package: Any, torch: Any, *, embedding: _Embedding, progress: bool
    ) -> tuple[Any, str]:
        """``sbi``'s trainer for this method, with the prior and the network."""
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
                network = builder(model=architecture, **{keyword: embedding.module})
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
            summaries.append(_summary_of(keep.observations, problem.datasets, batched=True))
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
        """
        assert self.training_set is not None
        if self._written:
            append_training_set(self.training_set, chunk, self.problem)
        else:
            write_training_set(self.training_set, chunk, self.problem)
            self._written = True

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
            "sbi_summary_layout": SUMMARY_LAYOUT,
            "sbi_summary_features": features,
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
