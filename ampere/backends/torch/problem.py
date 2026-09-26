"""The native path: a composed ``FittingProblem`` as one differentiable function.

Why this module has to exist
----------------------------
``DEVELOPMENT_PLAN.md`` §4.5's engine-facing surface —
``FittingProblem.log_prob``, ``log_prob_unconstrained``, ``prior_transform`` —
is the one an engine is written against, and ``ampere.inference``'s three
gradient-free drivers consume nothing else. A NUTS driver cannot, and this was
W2.4 slice 1's principal finding: that surface is **not traceable**, for
reasons that are structural rather than incidental and identical on every
backend.

1. ``ampere.core.results_schema``'s containers convert their values with
   ``numpy.asarray``. A ``ModelResult`` is made of those containers, so a
   model that computes its flux in torch hands the likelihood a *numpy* array
   and the graph is cut there. It does not merely stop: a tensor that requires
   grad cannot be converted at all (``RuntimeError: Can't call numpy() on
   Tensor that requires grad``).
2. ``ParameterSet.lnprior``/``lnprior_unconstrained`` short-circuit on
   ``math.isfinite`` — Python control flow on a value — and coerce with
   ``np.asarray(..., dtype=float)``.

Neither is a defect: the containers are what make a heterogeneous joint fit
checkable, and the short circuit is what stops an expensive model being
evaluated where the prior has no mass. Together they mean a gradient cannot be
taken *through the contract path*, on any backend.

``inference.md`` §10a names the answer — a **realisation**, the backend's
one-way translation of a composed problem into a native object exposing the
log density as a pure function of the unconstrained free vector — and this
module is torch's. It composes the same quantity out of this backend's own
surfaces: ``model.native_flux``, ``step.apply_flux``, ``noise.sigma_tensor`` where a
noise model has one, ``DenseGP.log_marginal_likelihood_native``, and
``TorchParameterSpace.log_prior_unconstrained_tensor``.

**It is deliberately the same walk as** :mod:`ampere.backends.jax.problem`,
line for line where the two libraries allow, because the two are one lowering
in two array libraries and a reviewer should be able to see that by reading
them side by side rather than by establishing it.

The two paths must agree, and the conformance suite is what says so
--------------------------------------------------------------------
Nothing here re-derives the mathematics: the Gaussian log-density is
``GaussianFamily.log_prob``'s closed form transcribed, the noise quadrature is
the noise model's own, and the change of variables is
``lnprior_unconstrained``'s. ``ampere.core.realise`` checks agreement with the
contract path at one point when the realisation is used; ``tests/conformance``
compares them at many, including near a support boundary, which is the proof
(``inference.md`` §10a, "What the conformance suite owes").

No exception control flow on the hot path
------------------------------------------
Every refusal this module makes — a model that is not this backend's, a
censored Student-t, a complex family, another backend's GP solver — is made
**at construction**, once, with a
:class:`~ampere.core.exceptions.LoweringError` naming what is unsupported.
That is ``inference.md`` §10a's reading of ``strict=True`` on the native path
and ``likelihoods.md`` §17 Q1's trace-purity ruling.

Inside the density a failure is ``-inf``, computed with :func:`torch.where`.
torch is eager rather than traced, so a Python ``if`` would *work* here in a
way it cannot in jax — and it is still not used, for two reasons that are not
about tracing. A ``torch.where`` keeps the result attached to the graph, which
a fresh ``-inf`` constant would not, and pyro's NUTS requires its potential to
have a ``grad_fn`` (it raises "element 0 of tensors does not require grad"
otherwise); and the two backends' failure semantics should not differ merely
because one of them could get away with it.

What this lowers, after slice 2
-------------------------------
W2.13 set the floor at ``inference.md`` §10a's minimum — this backend's models
and instrument steps, a ``GaussianFamily``, independent or dense-GP noise,
masks, plates and hierarchical priors. Slice 2 widens it to:

* **every likelihood family ``ampere.core`` implements for real data** —
  ``gaussian``, ``student_t``, ``cauchy`` and ``poisson`` — transcribed into
  torch in :mod:`ampere.backends.torch._families`;
* **censoring**, the Tobit form, for the families with a closed log-CDF
  (``gaussian`` and ``cauchy``). Censored ``student_t`` is refused by name:
  its CDF is a regularised incomplete beta function and ``torch.special`` has
  none;
* **the quasiseparable solver** inside the density, differentiable in the
  kernel hyperparameters — which is what makes a 10⁵-point flexible
  likelihood a thing a NUTS run can actually do;
* **prediction-aware noise**, ``FractionalModelNoise`` and
  ``FractionalModelGPNoise``, through the ``sigma_tensor`` hook W2.13 left
  dormant;
* **the latent-GP path** (W2.14): the whitened block ``z`` goes through the
  solver's ``latent_transform_native`` and the family body scores at
  ``f = L(θ) z``. Slice 2 refused this, and said why — ``ampere.core``'s own
  scoring path applied no whitening transform, so the kernel hyperparameters
  entered a latent likelihood nowhere at all, and lowering it here would have
  meant copying that defect or disagreeing with the oracle. W2.14 fixed the
  oracle; this backend follows it rather than leading it.

What slice 3 added
------------------
* **``complex_gaussian``**, and with it complex tensors end to end. A dataset
  whose observed container holds complex values (a
  :class:`~ampere.core.VisibilitySet`) keeps them complex here, and a model
  whose ``flux`` returns a complex tensor composes; the residual is complex,
  the sigma is real, and the density is real. What is *still* refused is the
  correlated case, because ``ampere.core`` refuses it first —
  :mod:`ampere.backends.torch._families` carries the reasoning and the closed
  form Phase 4 will implement.
* **a device**. Every tensor this module builds is built on the device the
  composed problem declares (``problem.device``, aggregated by
  :func:`~ampere.core.declared_capabilities` from parts that now declare it per
  instance), rather than on :data:`~ampere.backends.torch._config.DEFAULT_DEVICE`.
  A realisation therefore lives entirely where its pieces said they do, and a
  problem assembled half on a GPU is refused at composition rather than
  discovered inside a Cholesky.
* **a refusal for a prediction-aware noise model with no native surface.**
  ``_sigma`` consults ``noise.sigma_tensor`` and fell back to the base
  quadrature when there was none — which is right for
  :class:`~ampere.core.IndependentNoise`, whose sigma *is* the base, and
  silently wrong for anything that overrides ``sigma`` to depend on the
  prediction. Such a noise model is refused by name at construction now; see
  :meth:`_LoweredDataset._check_noise_surface`.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import torch

from ampere.core import (
    BatchedPrediction,
    ComplexGaussianFamily,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    LikelihoodFamily,
    PoissonFamily,
    StudentTFamily,
    VonMisesFamily,
    chunk_bounds,
    foreign_parts_refusal,
    sample_coordinates,
)
from ampere.core.dataset import (
    INSTRUMENT_COMPONENT,
    LATENT_COMPONENT,
    LIKELIHOOD_COMPONENT,
)
from ampere.core.exceptions import LoweringError
from ampere.core.likelihood import (
    DENSE_ROUTE,
    REDUCED_RANK_ROUTE,
    REDUCED_RANK_SOLVERS,
    ROTATED_ROUTE,
)

from ._config import (
    BACKEND,
    DEFAULT_DEVICE,
    DEFAULT_DTYPE,
    as_tensor,
    complex_dtype,
    resolve_device,
    to_numpy,
)
from ._families import FamilyInputs, limit_masks, native_log_prob, refuse_family
from .parameters import TorchParameterSpace

__all__ = ["LoweredProblem", "lower_problem"]

_LOG_2PI = math.log(2.0 * math.pi)


def _refuse(what: str, detail: str) -> LoweringError:
    return LoweringError(what, backend=BACKEND, detail=detail)


def _tensor(value: Any, *, device: torch.device = DEFAULT_DEVICE) -> torch.Tensor:
    return as_tensor(value, dtype=DEFAULT_DTYPE, device=device)


def _seed(seed: int) -> int:
    """A per-draw integer as a torch seed: non-negative and inside int64."""
    return int(seed) & 0x7FFFFFFFFFFFFFFF


def _context_rows(sigma: Any, dataset: Any, count: int) -> np.ndarray:
    """A chunk's per-draw context sigmas on *dataset*'s retained samples (**W5.29**).

    ``(count,) + observed.shape`` in, ``(count, n_retained)`` out: each draw's
    sigma flattened as the lowering flattens the container (``.ravel()``
    before the mask, *W5.5*) and restricted to the retained samples. A shape
    that does not match is refused by name, as
    ``Dataset.contextual_observed`` refuses it on the contract path.
    """
    array = np.asarray(sigma, dtype=float)
    shape = np.shape(dataset.dataset.observed.values)
    if array.shape != (count, *shape):
        raise _refuse(
            "context",
            f"dataset {dataset.label!r}: the per-draw context sigma has shape {array.shape}, "
            f"and a chunk of {count} draw(s) on an observed container of shape {shape} needs "
            f"{(count, *shape)} -- one sigma array of the observation's own shape per draw.",
        )
    return array.reshape(count, -1)[:, dataset.retain]


#: The two spellings of a native model's value-and-coordinates surface, in the
#: order they are looked for (*W4.3*, reordered by *W5.20*).
#: ``native_flux``/``native_grid`` is the **canonical** pair; ``flux``/``grid``
#: is the original spelling, kept as a supported **legacy alias** so that
#: models written before the ruling still compose.
#:
#: W4.3 introduced ``native_*`` because a model could not then have a method
#: whose name one of its own parameters used, and ``flux`` is precisely what an
#: interferometric source model calls its total flux density. W5.20 removed
#: that constraint — the reserved set is core's and no longer grows with a
#: backend's class attributes — but kept ``native_*`` as the canonical name,
#: because a surface a *realisation* composes and a quantity a *user* fits
#: should not compete for one word.
#:
#: Either pair composes; a model that offers **both** is refused as ambiguous,
#: naming both methods, rather than silently taken at one of them; one that
#: offers half of either is refused with the missing half named, because the
#: two go together.
_FLUX_NAMES: tuple[str, ...] = ("native_flux", "flux")
_GRID_NAMES: tuple[str, ...] = ("native_grid", "grid")


def _native_surface(model: Any, names: tuple[str, ...], label: str) -> Any | None:
    """*model*'s surface under the first of *names* it offers, or ``None``.

    Refuses a model that offers both spellings (*W5.20*). Before the ruling
    the lookup took the first match and said nothing, so a model carrying a
    legacy ``flux`` beside a canonical ``native_flux`` — the shape a
    half-finished rename leaves behind — composed at whichever the table
    happened to list first. Two methods, one surface, and no way to tell which
    the author meant: that is a question for the author, not a coin toss.
    """
    offered = [name for name in names if callable(getattr(model, name, None))]
    if len(offered) > 1:
        raise _refuse(
            type(model).__name__,
            f"model {label!r} offers both `{names[0]}` and `{names[1]}`, which are two "
            f"spellings of one surface: `{names[0]}` is canonical and `{names[1]}` is the "
            f"legacy alias. A model that defines both leaves nothing to choose between them "
            f"— keep one, and prefer `{names[0]}`.",
        )
    return getattr(model, offered[0]) if offered else None


def _require_flux(model: Any, label: str) -> Any:
    """*model*'s native value surface, or a refusal naming the model (*W4.3*).

    A problem may hold a model no dataset binds, and every channel of every
    model is written to a training set (``results.md`` §11), so this path reaches
    a model the per-dataset refusal never checked. Naming it is the whole point:
    before the two-spelling lookup this was an ``AttributeError`` on ``flux``,
    which at least said which attribute was missing, and a bare ``None`` call
    would say nothing at all.
    """
    found = _native_surface(model, _FLUX_NAMES, label)
    if found is None:
        raise _refuse(
            type(model).__name__,
            f"model {label!r} has no native surface (`{_FLUX_NAMES[0]}` or "
            f"`{_FLUX_NAMES[1]}`), so its channels cannot be produced natively.",
        )
    return found


def _flat_channel(flux: torch.Tensor) -> torch.Tensor:
    """One draw's channel values as a flat vector (*W4.3*).

    ``BatchedPrediction.channels`` is declared ``{model: {channel: (batch, n)}}``
    — flat per draw, in the container's own C-order, because that is what
    ``ampere.core``'s agreement check compares against
    ``result[channel].values.ravel()`` and what a training-set group is written
    from. Every channel before Phase 4 was one-dimensional already, so this was
    a no-op nobody had to write; an ``Image`` channel is ``(nx, ny)`` and needs
    it. Flattening here rather than in the model keeps the native forward
    surface the natural shape for the step that consumes it —
    ``FourierSample.apply_flux`` contracts over two axes — and is applied under
    ``vmap``, where the leading batch axis is hidden, so ``reshape(-1)`` is the
    whole of it.
    """
    return flux.reshape(-1)


#: Neutral family name -> the ``ampere.core`` class whose ``sample`` this
#: backend has a native twin for (*W3.14*). The **one** place this module names
#: them, so the ceiling ``inference.md`` §13 states — a backend samples exactly
#: what the core samples — is one table rather than a condition repeated per
#: family: a core family that gains a ``sample`` joins by adding a row here and
#: a branch to :meth:`_LoweredDataset.draw_chunk`, and one that has no twin
#: here falls back to the numpy path unchanged.
_TWINNED_FAMILIES: dict[str, type[LikelihoodFamily]] = {
    "gaussian": GaussianFamily,
    "poisson": PoissonFamily,
    "student_t": StudentTFamily,
    "complex_gaussian": ComplexGaussianFamily,
    "von_mises": VonMisesFamily,
}


#: The noise models whose ``sigma`` this module transcribes in
#: :meth:`_LoweredDataset._base_sigma`. A subclass that overrides ``sigma`` is
#: computing something else, and unless it also supplies ``sigma_tensor`` the
#: realised density would quietly use the base quadrature instead — a different
#: likelihood from the contract path's, arrived at without a word. See
#: :meth:`_LoweredDataset._check_noise_surface`.
_TRANSCRIBED_SIGMA = (IndependentNoise.sigma, GaussianProcessNoise.sigma)


class _LoweredDataset:
    """One dataset's forward chain and log-likelihood, as torch."""

    def __init__(
        self,
        problem: FittingProblem,
        label: str,
        *,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        dataset = problem.datasets[label]
        self.device = device
        self.label = label
        self.model_label = problem.bindings[label]
        self.model = problem.models[self.model_label]
        self.channel = dataset.instrument.channel
        self.dataset = dataset
        self.likelihood = dataset.likelihood
        self.noise = dataset.likelihood.noise
        self.steps = tuple(dataset.instrument.steps)

        self.model_flux = _native_surface(self.model, _FLUX_NAMES, self.model_label)
        self.model_grid = _native_surface(self.model, _GRID_NAMES, self.model_label)
        if self.model_flux is None or self.model_grid is None:
            missing = ", ".join(
                f"`{names[0]}` (or `{names[1]}`)"
                for names, found in ((_FLUX_NAMES, self.model_flux), (_GRID_NAMES, self.model_grid))
                if found is None
            )
            raise _refuse(
                type(self.model).__name__,
                f"model {self.model_label!r} has no native torch surface ({missing} missing), so "
                f"it cannot be composed into a differentiable log-density. The two go together: "
                f"the first supplies the values and the second the coordinates the instrument "
                f"chain transforms them on, and `predict` calls both. Build the problem from "
                f"ampere.backends.torch's models, or run it on a gradient-free engine.",
            )
        for step in self.steps:
            if not hasattr(step, "apply_flux"):
                raise _refuse(
                    type(step).__name__,
                    f"instrument step {step.label!r} of dataset {label!r} has no native torch "
                    f"surface (an `apply_flux` method). Build the chain from "
                    f"ampere.backends.torch's steps.",
                )
        censoring = dataset.likelihood.censoring
        refusal = refuse_family(
            self.likelihood.family,
            censored=censoring is not None and bool(censoring.any_censored),
            latent=dataset.latent is not None,
            correlated=bool(getattr(self.noise, "CORRELATED", False)),
        )
        if refusal is not None:
            raise refusal
        # The latent block's local name, resolved once at composition. Until
        # W2.14 this was a blanket refusal: ampere.core scored a latent
        # likelihood at the whitened z rather than at f = L(theta) z, so the
        # kernel hyperparameters entered it nowhere, and lowering the path
        # here would have meant either copying that defect or disagreeing with
        # the oracle a realisation is checked against. The core applies the
        # transform now (GaussianProcessNoise.noise_params), so this backend
        # applies its own native one and the conformance row holds the two
        # together.
        self.latent_name: str | None = (
            None if dataset.latent is None else dataset.latent.parameter.name
        )
        correlated = bool(getattr(self.noise, "CORRELATED", False))
        if correlated and getattr(self.noise.solver, "BACKEND", "reference") != BACKEND:
            raise _refuse(
                type(self.noise.solver).__name__,
                f"dataset {label!r} uses the {type(self.noise.solver).__name__} solver, which is "
                f"not this backend's. A torch problem whose GP solve ran in numpy would not be "
                f"differentiable; pass ampere.backends.torch.DenseGP or "
                f"ampere.backends.torch.QuasisepGP.",
            )
        if correlated and not hasattr(self.noise.solver, "log_marginal_likelihood_native"):
            raise _refuse(
                type(self.noise.solver).__name__,
                f"dataset {label!r} uses the {type(self.noise.solver).__name__} solver, which "
                f"declares this backend but has no `log_marginal_likelihood_native` -- the "
                f"non-raising, differentiable surface a realised density needs (inference.md "
                f"§10a). A user-written solver joins the native path by supplying it.",
            )
        if self.latent_name is not None and not callable(
            getattr(self.noise.solver, "latent_transform_native", None)
        ):
            raise _refuse(
                type(self.noise.solver).__name__,
                f"dataset {label!r} declares a latent GP, but the "
                f"{type(self.noise.solver).__name__} solver has no `latent_transform_native` -- "
                f"the non-raising, differentiable form of `f = L(theta) z`. Since W2.14 the "
                f"contract path applies that transform on every latent evaluation, so a "
                f"realisation without it would disagree with its own oracle and would give the "
                f"kernel hyperparameters no gradient. A user-written solver joins the latent "
                f"path by supplying it.",
            )
        if not isinstance(self.noise, (IndependentNoise, GaussianProcessNoise)):
            raise _refuse(
                type(self.noise).__name__,
                f"dataset {label!r} uses a noise model this backend does not know how to "
                f"compose natively.",
            )
        self._check_noise_surface()
        self.correlated = correlated
        # Which of the two correlated stories this dataset tells, decided by
        # the family's own declaration rather than by its name. A family whose
        # GP marginalises in closed form (the Gaussian one) hands the whole
        # covariance to the solver; one that does not (Poisson) needs the
        # latent block instead, and takes its own body with `noise.latent`
        # supplied. Before W2.14 the branch below tested `correlated` alone,
        # which was right only because a latent dataset never got this far.
        self.gp_marginal = correlated and self.likelihood.family.ANALYTIC_WITH_GP

        observed = dataset.observed
        # The effective mask, resolved once at composition: Dataset takes the
        # union of the observed and predicted masks there and forbids a
        # parameter-dependent one, so which samples are retained is a fact
        # about the problem rather than about theta and can be a constant here.
        excluded = dataset.effective_mask
        if excluded is None:
            excluded = ~np.asarray(observed.valid, dtype=bool).ravel()
        self.retain = ~np.asarray(excluded, dtype=bool).ravel()
        # Complex data stay complex (W2.4 slice 3). ``dtype=float`` here was
        # what made ``complex_gaussian`` unreachable, and it would not have
        # failed loudly if it had been reached: numpy 2.5 *warns*
        # (ComplexWarning) and discards the imaginary part rather than raising,
        # so the density would have been the one for the real projection of the
        # data. What kept that from happening was the family refusal at
        # construction, which is exactly why the refusal came first and the
        # dtype second. The container's own dtype kind decides now, once, here.
        # ``.ravel()`` before the mask, not after: ``retain`` is flat, one entry
        # per *sample*, and on a ``Layout.GRID`` container the values are not
        # (W5.5). Without it a 12x12 image indexed by a 144-long mask raises an
        # IndexError about axis 0 that says nothing about layouts.
        #: How many axes the observed container's values carry: one for a
        #: ``Layout.POINTS`` kind, two for an ``Image``, three for a ``Cube``.
        #: What :meth:`predict` needs in order to flatten the container's own
        #: axes without touching a batch axis in front of them (W5.5).
        self._value_ndim = int(np.asarray(observed.values).ndim)
        values = np.asarray(observed.values).ravel()
        self.complex_valued = values.dtype.kind == "c"
        self.observed_values = as_tensor(
            values[self.retain],
            dtype=complex_dtype(DEFAULT_DTYPE) if self.complex_valued else DEFAULT_DTYPE,
            device=device,
        )
        # Every axis, stacked into the `(n, d)` block ``Likelihood._coordinates``
        # builds -- **not** ``axes[0]`` (*W4.2*), and **not** a bare
        # ``column_stack`` either (*W5.5*): on a ``Layout.GRID`` container the
        # axes are the separable grids the samples are the product of, so a
        # 24x24 image has two axes of 24 coordinates and 576 samples.
        # ``ampere.core.sample_coordinates`` is that stacking made layout-aware,
        # and it is the same function ``Dataset.draw_observation`` uses, so the
        # coordinates a kernel sees here and there cannot drift. A one-axis
        # point container is unchanged by either (``(n, 1)`` and ``(n,)`` are
        # the same point set to ``Kernel.matrix``), and a VisibilitySet has
        # three axes, so taking the first would have handed the kernel the `u`
        # column and called it the coordinates. The same stack is what makes
        # ``axes=("u", "v")`` meaningful here, since a selector resolves to
        # *columns of this block*.
        self.observed_coordinates = _tensor(
            sample_coordinates(observed)[self.retain],
            device=device,
        )
        # The kernel bound to this container's axis order, once, at lowering:
        # W4.5 put the binding on the noise model because that is where the
        # container is, and the native path has to use it or an ``axes=``
        # selector would silently do nothing here while working on the contract
        # path it is checked against.
        self.bound_kernel = (
            self.noise.kernel_for(observed) if getattr(self.noise, "CORRELATED", False) else None
        )
        self.uncertainty = (
            None
            if observed.uncertainty is None
            else _tensor(
                np.asarray(observed.uncertainty, dtype=float).ravel()[self.retain], device=device
            )
        )
        if self.gp_marginal and self.likelihood.family.REQUIRES_UNCERTAINTY:
            self._check_gp_uncertainty(observed, label)
        # Which samples are limits is a fact about the data, so the three
        # groups are resolved once here rather than per evaluation. A limit on
        # a masked sample is not a limit at all -- likelihoods.md §9's "masking
        # beats censoring" -- so the declaration is read *through* the
        # effective mask, exactly as Likelihood._retained_limits does.
        self.detection, self.upper, self.lower = (None, None, None)
        if censoring is not None:
            censoring.check_against(observed)
            self.detection, self.upper, self.lower = limit_masks(
                np.asarray(censoring.kinds)[self.retain], device=device
            )

    def _check_gp_uncertainty(self, observed: Any, label: str) -> None:
        """Refuse, at construction, a GP-marginal dataset the contract path cannot normalise.

        The twin of ``ampere.backends.jax.problem._LoweredDataset``'s method of
        the same name (W3.0 landed it there; the identical gap here is the
        carried finding W3.1 slice 2 closes). Same condition, same wording, for
        the same reason: ``GaussianProcessNoise.sigma`` refuses a retained
        uncertainty that is zero or negative -- "an infinitely precise
        measurement, which no likelihood can normalise" -- for any family that
        ``REQUIRES_UNCERTAINTY``, and on the contract path that surfaces only
        when the density is evaluated. ``ampere.core.realise``'s one-point
        agreement check happens to catch it today, because the contract path
        raises while computing the reference log-probability, but that is
        incidental to which one point gets checked rather than a refusal this
        backend makes by name -- and a ``strict`` NUTS run would sample a
        density the contract refuses. So the refusal belongs here, at
        construction, beside the ``REQUIRES_UNCERTAINTY``/sigma-is-None check
        :meth:`log_likelihood` already makes for the non-GP branch.

        **W5.2**: no longer checks ``observed.uncertainty is None`` -- a
        dataset with no observed uncertainty at all and a family that
        ``REQUIRES_UNCERTAINTY`` is already refused unconditionally at
        composition, by ``NoiseModel.check_compatible``
        (``ampere.core.likelihood``), which every ``Likelihood`` calls before
        a ``FittingProblem`` exists at all -- and a ``_LoweredDataset`` is
        never built except from one. That branch (with its own "no jitter, so
        sigma is undefined" wording, which core's unconditional refusal does
        not honour anyway -- jitter never rescues a missing uncertainty at
        composition, so the branch could not even have been reached the way
        it was written) was dead code; removed rather than exercised, since
        making it reachable would mean relaxing the core refusal, which is a
        §4 contract change this item is not scoped to make.
        """
        retained = np.asarray(observed.uncertainty, dtype=float).ravel()[self.retain]
        bad = retained <= 0.0
        if np.any(bad):
            raise _refuse(
                "uncertainty",
                f"dataset {label!r} was given zero or negative uncertainties on "
                f"{int(np.sum(bad))} retained sample(s). A zero uncertainty is an "
                f"infinitely precise measurement, which no likelihood can normalise; "
                f"mask the sample, or give it a real error bar.",
            )

    def _check_noise_surface(self) -> None:
        """Refuse a prediction-aware noise model with no ``sigma_tensor``.

        W2.4 slice 3, item 4: the hook's last consumer gap. :meth:`_sigma`
        computes the base quadrature and hands it to ``noise.sigma_tensor``
        where there is one; where there is not, it used to return the base —
        correct for the two noise models this module transcribes, and silently
        *wrong* for any subclass that overrides ``sigma`` to depend on the
        prediction, because the realised density would then be a different
        likelihood from the contract path's with nothing said.

        ``ampere.core`` has no "prediction-aware" flag to test, and inventing
        one would be a contract change; what it does have is the two ``sigma``
        implementations this module reproduces, so the test is whether the
        noise model still uses one of them. A subclass that overrides ``sigma``
        and supplies ``sigma_tensor`` — this backend's
        :class:`~ampere.backends.torch.FractionalModelNoise` and
        :class:`~ampere.backends.torch.FractionalModelGPNoise` — passes; one
        that overrides it and does not is refused by name, with the remedy in
        the message.
        """
        if hasattr(self.noise, "sigma_tensor"):
            return
        if type(self.noise).sigma in _TRANSCRIBED_SIGMA:
            return
        raise _refuse(
            type(self.noise).__name__,
            f"dataset {self.label!r} uses a noise model that overrides sigma() but supplies no "
            f"native `sigma_tensor`. A realised density transcribes the base quadrature "
            f"sqrt((scale*sigma_data)**2 + jitter**2) in torch and hands it to sigma_tensor; "
            f"with no such method it would silently score the base instead of whatever sigma() "
            f"computes -- a different likelihood from the contract path's. Add sigma_tensor(base, "
            f"values, *, predicted) returning a tensor (see "
            f"ampere.backends.torch.FractionalModelNoise), or run this problem on a "
            f"gradient-free engine, which uses sigma() directly.",
        )

    # -- routing ------------------------------------------------------------

    def _dataset_values(self, routed: Mapping[str, Mapping[str, Any]]) -> Mapping[str, Any]:
        return self.dataset.route(routed[self.label])

    def _instrument_values(self, routed: Mapping[str, Mapping[str, Any]]) -> Mapping[str, Any]:
        chain = self._dataset_values(routed).get(INSTRUMENT_COMPONENT)
        if not chain:
            return {}
        return self.dataset.instrument.mapping.distribute(chain)

    def _step_values(
        self, routed: Mapping[str, Mapping[str, Any]], label: str
    ) -> Mapping[str, Any]:
        """One step's own local values, through the core's two levels of routing.

        Problem to dataset (``routed[self.label]``), dataset to instrument
        (``Dataset.route``), instrument to step (``Instrument.mapping``). Every
        level is the core's own routing metadata, reused rather than
        reimplemented: the lossless-nesting rule (2026-09-02) is what makes an
        inner ``shared_as`` collapse visible at every level above, and a
        backend that walked the chain itself would lose it.
        """
        return self._instrument_values(routed).get(label, {})

    # -- the forward chain --------------------------------------------------

    def predict(self, routed: Mapping[str, Mapping[str, Any]]) -> torch.Tensor:
        """The instrument-transformed prediction, retained samples only.

        The chain is walked in declaration order, each step handed its own
        local values by the routing the core already computed. Masking is
        applied **last** and by index, because ``Dataset`` resolves the
        effective mask once at construction and forbids a parameter-dependent
        one — so which samples are retained is a fact about the problem, not
        about θ, and can be a constant here.

        ``retain`` is flat, one entry per *sample*, and on a ``Layout.GRID``
        container the prediction is not (**W5.5**): a 24x24 image is a
        ``(24, 24)`` tensor and a 576-long mask. The grid branch flattens the
        container's own trailing axes — and only those, so a leading batch axis
        rides through — before applying it. The point-set branch is left exactly
        as it was rather than folded into the same expression, because a
        one-axis container's indexing is what every other modality's rows are
        written against and this is not the item to change it in.
        """
        full = self.predict_full(routed)
        if self._value_ndim == 1:
            return full[self.retain]
        flat = full.reshape(*full.shape[: full.ndim - self._value_ndim], -1)
        return flat[..., self.retain]

    def predict_full(self, routed: Mapping[str, Mapping[str, Any]]) -> torch.Tensor:
        """The same chain, **before** the mask is applied (W3.1 slice 2).

        The density wants the retained samples; a *simulation* wants the whole
        container, because ``Dataset.predict`` returns one and a masked sample
        keeps the observed container's own value rather than vanishing from it
        (``inference.md`` §13). Splitting the two here rather than re-walking
        the chain keeps one implementation of the forward model, which is the
        only way ``simulate_batched`` can be checked against ``log_likelihood``
        at all.
        """
        flux = self.model_flux(self.channel, routed[self.model_label])
        grid = self.model_grid(self.channel)
        for step in self.steps:
            flux, grid = step.apply_flux(flux, grid, self._step_values(routed, step.label))
        return flux

    # -- the log-likelihood -------------------------------------------------

    def _base_sigma(
        self, values: Mapping[str, Any], uncertainty: torch.Tensor | None = None
    ) -> torch.Tensor | None:
        """``sqrt((scale * sigma_data)**2 + jitter**2)`` for the retained samples.

        Transcribed rather than called, because ``NoiseModel.sigma`` coerces
        with ``float()`` and ``scale`` and ``jitter`` are fitted parameters:
        calling it would cut their gradient exactly as calling
        ``Kernel.matrix`` used to cut the GP amplitude's (W2.4 slice 1's
        finding). ``None`` when there are no uncertainties and no jitter — a
        family that cannot live without them has already refused at
        composition through ``NoiseModel.check_compatible``.

        *uncertainty* (**W5.29**) is one draw's observation-context sigma on
        the retained samples, standing in for the cached ``sigma_data``
        exactly as ``Dataset.contextual_observed`` stands a context's sigma in
        for the container's on the contract path — so ``scale``, ``jitter``
        and a fractional inflation apply to it as they apply to the
        observation's own. ``None`` is the observation's own, unchanged.
        """
        own = {key: value for key, value in values.items() if key in self.noise.parameters}
        resolved = self.noise.context(own)
        sigma = self.uncertainty if uncertainty is None else uncertainty
        if sigma is None:
            if "jitter" not in resolved:
                return None
            return (
                _tensor(resolved["jitter"], device=self.device)
                .expand(int(self.retain.sum()))
                .clone()
            )
        if "scale" in resolved:
            sigma = sigma * _tensor(resolved["scale"], device=self.device)
        if "jitter" in resolved:
            floor = _tensor(resolved["jitter"], device=self.device)
            sigma = torch.sqrt(sigma**2 + floor**2)
        return sigma

    def _sigma(
        self,
        predicted: torch.Tensor,
        values: Mapping[str, Any],
        uncertainty: torch.Tensor | None = None,
    ) -> torch.Tensor | None:
        """The effective ``sigma``, prediction-aware noise models included.

        The base quadrature is computed here (once, from the cached
        uncertainties, or from one draw's context sigma when *uncertainty* is
        given — :meth:`_base_sigma`) and *handed* to a noise model that
        declares ``sigma_tensor`` — this backend's ``FractionalModelNoise``
        and ``FractionalModelGPNoise``, which inflate it by ``f * |predicted|``.
        The hook takes the base rather than the container and the mask so that
        the quadrature is written once and every parameter in it keeps its
        graph; see :mod:`ampere.backends.torch.noise`.
        """
        base = self._base_sigma(values, uncertainty)
        native = getattr(self.noise, "sigma_tensor", None)
        if native is None:
            return base
        return native(base, values, predicted=predicted)

    def _latent(
        self, routed: Mapping[str, Mapping[str, Any]], values: Mapping[str, Any]
    ) -> torch.Tensor | None:
        """The dataset's latent function ``f``, or ``None`` when it declares none.

        The whitened values ``z`` arrive as one array-valued parameter under
        the ``latent`` component (``Dataset.route``); the correlation is
        imposed here by the solver's own **native** whitening transform, so
        what the family body receives is ``f = L(θ) z`` on the retained
        coordinates — which is what ``PoissonFamily`` has always read
        ``noise.latent`` as, and what ``ampere.core`` now supplies.

        ``latent_transform_native`` rather than ``latent_transform`` for the
        usual reason: the contract surface returns numpy and would cut the
        gradient in exactly the hyperparameters this path exists to fit.
        """
        if self.latent_name is None:
            return None
        block = self._dataset_values(routed).get(LATENT_COMPONENT, {})
        whitened = _tensor(block[self.latent_name], device=self.device).reshape(-1)
        kernel = self._kernel()
        return self.noise.solver.latent_transform_native(
            kernel,
            self.observed_coordinates,
            whitened,
            kernel.resolve(values),
        )

    # -- the generative half (W3.1 slice 2) ---------------------------------

    def sampling_refusal(self) -> LoweringError | None:
        """Why this dataset cannot be sampled natively, or ``None``.

        The twin of ``ampere.backends.jax.problem._LoweredDataset``'s method of
        the same name, and the same conditions — see it for the reasoning. In
        one sentence: Peter's ruling of 2026-09-08 is that every backend
        samples natively, and what a backend may sample is fixed by what
        ``ampere.core`` samples, because the numpy path is the oracle; a native
        twin for a family the core declines to guess a sampling distribution
        for would be the "silently train an SBI posterior on the wrong forward
        model" §13's refusal exists to prevent. Everything else is handed back
        to the numpy path, where the core's own refusal text is what the caller
        meets.

        **Amended W3.14**: the set the core samples is
        :data:`_TWINNED_FAMILIES` — ``gaussian``, ``poisson``, ``student_t``
        and ``complex_gaussian`` — and the test is that the family's ``sample``
        is the *core's own* method for that family, so a user override still
        falls back to the numpy function its author wrote. The complex
        condition is gone (complex data are what the ``complex_gaussian`` twin
        exists for) and the latent condition now applies only to the families
        that do not read ``noise.latent`` — every twin but ``poisson`` and,
        since **W5.1**, ``von_mises``, whose latent draw is native only (the
        core's own ``sample`` refuses it by name, so there is no numpy draw to
        fall back to).
        """
        family = self.likelihood.family
        core = _TWINNED_FAMILIES.get(family.NAME)
        if core is None or type(family).sample is not core.sample:
            return _refuse(
                family.NAME or type(family).__name__,
                f"dataset {self.label!r} does not draw its observations through one of "
                f"ampere.core's own family sample() methods, so this backend has no native "
                f"twin for it: the numpy path is the oracle, and a backend that guessed a "
                f"sampling distribution the contract declines to guess would train an SBI "
                f"posterior on the wrong forward model. The draw is made on the numpy path "
                f"instead.",
            )
        censoring = self.likelihood.censoring
        if censoring is not None and bool(np.any(np.asarray(censoring.kinds)[self.retain] != 0)):
            return _refuse(
                "censoring",
                f"dataset {self.label!r} declares limits on retained samples, which blocks "
                f"observation drawing on every backend.",
            )
        if self.latent_name is not None and family.NAME not in ("poisson", "von_mises"):
            return _refuse(
                "latent",
                f"dataset {self.label!r} declares a latent GP, whose family reads "
                f"noise.latent; the {family.NAME!r} twin does not consume one.",
            )
        return None

    def draw_chunk(
        self,
        route: Callable[[torch.Tensor], dict[str, dict[str, Any]]],
        stack: torch.Tensor,
        predicted: torch.Tensor,
        seeds: Sequence[int],
        device: torch.device,
        uncertainty: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """A whole chunk's retained observed values, natively (*W3.14*).

        *uncertainty* (**W5.29**) is ``None`` — every draw at the observation's
        own sigma — or a ``(chunk, n_retained)`` stack of per-draw
        observation-context sigmas, mapped in beside θ and the prediction so
        that draw *i*'s noise level is its own context's
        (:meth:`_base_sigma`).

        Two shapes of draw, and which one a family takes is decided by whether
        its variate can be written as arithmetic on standard normals.

        * **One stage** (``gaussian``, ``complex_gaussian``). The normals are
          drawn outside the transform, per draw, from that draw's own seed, and
          mapped in as data; :meth:`sample_retained` then assembles the value.
          ``torch.func.vmap`` has no per-sample random state — a
          ``torch.Generator`` is not a tensor and cannot be mapped over, and
          ``randomness="different"`` would drive the draw from the *global*
          generator and so make it depend on how the chunk was scheduled — so
          this is what makes a draw a pure function of its seed, which is what
          partition independence needs.
        * **Two stage** (``poisson``, ``student_t``, ``von_mises``). Their
          variates are not arithmetic on normals, so the transform computes
          each draw's distribution *parameters* — the Poisson rate,
          ``exp(f)`` and all; the Student-t location, scale and ``nu``; the
          von Mises mean angle and concentration ``kappa = 1/sigma**2``, read
          exactly as ``ampere.core.VonMisesFamily.sample`` reads it — and the
          variate is drawn after it, per draw, from that draw's own seed. Same
          property, one stage later. jax needs none of this: ``jax.random``
          takes a traced rate and a traced ``nu``/``kappa`` inside
          ``jax.vmap``.

        The variates are drawn on the **CPU** and moved, exactly as the normals
        already were: a ``torch.Generator`` bound to a device would make a
        budget's randomness a function of where it ran.
        """
        name = self.likelihood.family.NAME
        size = int(self.retain.sum())
        if name in ("poisson", "student_t", "von_mises"):

            def parameters(vector: torch.Tensor, row: torch.Tensor, of: Any = self) -> Any:
                return of.sample_parameters(route(vector), row)

            def contextual(
                vector: torch.Tensor, row: torch.Tensor, sigma: torch.Tensor, of: Any = self
            ) -> Any:
                return of.sample_parameters(route(vector), row, sigma)

            resolved = (
                torch.func.vmap(parameters)(stack, predicted)
                if uncertainty is None
                else torch.func.vmap(contextual)(stack, predicted, uncertainty)
            )
            if name == "poisson":
                return self._poisson_variates(resolved, seeds, device)
            if name == "von_mises":
                return self._von_mises_variates(resolved, seeds, size, device)
            return self._student_t_variates(resolved, seeds, size, device)

        # Two standard normals per draw, or **four** for the circular complex GP
        # (*W4.2*): that draw needs a correlated and an independent block in each
        # of the real and the imaginary parts, and the two components must share
        # the factor L and nothing else. ``torch.randn`` fills row-major, so the
        # first two rows of a ``(4, size)`` draw are byte-identical to the
        # ``(2, size)`` draw from the same seed -- every stream this did not
        # change is unchanged.
        rows = 4 if (self.correlated and name == "complex_gaussian") else 2
        normals = torch.stack(
            [
                torch.randn(
                    (rows, size),
                    generator=torch.Generator(device="cpu").manual_seed(_seed(seed)),
                    dtype=DEFAULT_DTYPE,
                ).to(device)
                for seed in seeds
            ]
        )

        def one(
            vector: torch.Tensor, row: torch.Tensor, noise: torch.Tensor, of: Any = self
        ) -> torch.Tensor:
            return of.sample_retained(route(vector), row, noise)

        def one_in_context(
            vector: torch.Tensor,
            row: torch.Tensor,
            noise: torch.Tensor,
            sigma: torch.Tensor,
            of: Any = self,
        ) -> torch.Tensor:
            return of.sample_retained(route(vector), row, noise, sigma)

        if uncertainty is None:
            return torch.func.vmap(one)(stack, predicted, normals)
        return torch.func.vmap(one_in_context)(stack, predicted, normals, uncertainty)

    def sample_retained(
        self,
        routed: Mapping[str, Mapping[str, Any]],
        predicted: torch.Tensor,
        normals: torch.Tensor,
        uncertainty: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """One draw of the retained observed values, from two standard normals.

        The one-stage families (:meth:`draw_chunk`). ``normals`` is ``(2, n)``,
        already drawn from this draw's own seed.

        * ``gaussian`` — ``x = mu + sigma z`` uncorrelated, and
          ``x = mu + L z1 + sigma z2`` under a GP, with ``L`` from the solver's
          **own** ``latent_transform_native``, the same whitening the latent
          path uses, and the solver's numerical jitter folded into the diagonal
          because it is part of the covariance the marginal likelihood scores.
          Omitting it would draw from a narrower distribution than the density
          evaluates.
        * ``complex_gaussian`` (*W3.14*) — ``x = mu + sigma (z1 + i z2)``: the
          two normals become the two components, independent and of equal
          variance, which is the circular symmetry the family's density
          assumes. ``sigma`` is the **per-component** standard deviation, as
          ``-|r|**2 / (2 sigma**2) - log 2pi - log sigma**2`` implies.

        *uncertainty* is this draw's context sigma, or ``None`` (**W5.29**).
        """
        values = dict(self._dataset_values(routed).get(LIKELIHOOD_COMPONENT, {}))
        sigma = self._sigma(predicted, values, uncertainty)
        if self.likelihood.family.NAME == "complex_gaussian":
            assert sigma is not None  # REQUIRES_UNCERTAINTY, checked at composition
            if not self.correlated:
                return predicted + sigma * torch.complex(normals[0], normals[1])
            # W4.2: two real GP realisations sharing L, not one complex one.
            # Sharing L is circularity's equal-component half; sharing nothing
            # else is its zero-pseudo-covariance half. Drawing one realisation
            # for both components would give a draw perfectly correlated between
            # the parts, whose modulus statistics the density does not score.
            whitened = torch.stack([normals[0], normals[1]], dim=-1)
            correlated = self._gp_realisation(values, whitened)
            stabiliser = float(getattr(self.noise.solver, "jitter", 0.0) or 0.0)
            scale = sigma if not stabiliser else torch.sqrt(sigma**2 + stabiliser**2)
            independent = torch.stack([normals[2], normals[3]], dim=-1)
            components = correlated + scale.unsqueeze(-1) * independent
            return predicted + torch.complex(components[..., 0], components[..., 1])
        realisation = predicted
        if self.correlated:
            realisation = realisation + self._gp_realisation(values, normals[0])
            stabiliser = float(getattr(self.noise.solver, "jitter", 0.0) or 0.0)
            if stabiliser:
                floor = torch.full_like(normals[1], stabiliser)
                sigma = floor if sigma is None else torch.sqrt(sigma**2 + floor**2)
        if sigma is None:
            return realisation
        return realisation + sigma * normals[1]

    def sample_parameters(
        self,
        routed: Mapping[str, Mapping[str, Any]],
        predicted: torch.Tensor,
        uncertainty: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        """One draw's distribution parameters, for the two-stage families (*W3.14*).

        Everything about the draw that depends on θ, and nothing that depends
        on randomness — so it can go through ``torch.func.vmap`` while the
        variate cannot. ``poisson`` returns its ``rate`` (the prediction, times
        ``exp(f)`` under a GP, with ``f`` the dataset's own latent block, so θ
        and the drawn counts describe one model); ``student_t`` returns the
        ``location``, the ``scale`` — the noise model's sigma, used as the
        scale exactly as the density standardises by it — and ``nu``;
        ``von_mises`` (**W5.2**) returns the ``mean`` angle (the prediction,
        plus the latent phase error ``f`` under a GP since **W5.1**) and the
        concentration ``kappa = 1/sigma**2``, read from the same
        ``sigma`` its density standardises by — exactly
        ``ampere.core.VonMisesFamily.sample``'s own ``rng.vonmises(mean,
        1.0 / sigma**2)``, one family earlier.
        """
        values = dict(self._dataset_values(routed).get(LIKELIHOOD_COMPONENT, {}))
        family = self.likelihood.family
        if family.NAME == "poisson":
            rate = predicted
            latent = self._latent(routed, values)
            if latent is not None:
                rate = rate * torch.exp(latent)
            return {"rate": rate}
        # W5.29: this draw's context sigma, or None for the observation's own.
        # Poisson above reads no sigma, so a context changes nothing there --
        # the same answer ``PoissonFamily.sample`` gives on the contract path.
        sigma = self._sigma(predicted, values, uncertainty)
        assert sigma is not None  # REQUIRES_UNCERTAINTY, checked at composition
        if family.NAME == "von_mises":
            # W5.1: under a latent GP the draw is around predicted + f, f being
            # this draw's own latent block -- the latent first, then the
            # wrapped variate, and the same f the density scores at.
            latent = self._latent(routed, values)
            mean = predicted if latent is None else predicted + latent
            return {"mean": mean, "kappa": 1.0 / sigma**2}
        own = {key: value for key, value in values.items() if key in family.parameters}
        return {
            "location": predicted,
            "scale": sigma,
            "nu": _tensor(family.context(own)["nu"], device=predicted.device),
        }

    def _poisson_variates(
        self, resolved: Mapping[str, torch.Tensor], seeds: Sequence[int], device: torch.device
    ) -> torch.Tensor:
        rates = resolved["rate"].detach().to("cpu", dtype=DEFAULT_DTYPE)
        if not bool(torch.isfinite(rates).all()) or bool((rates < 0.0).any()):
            raise _refuse(
                "poisson",
                f"dataset {self.label!r}: the model predicted a negative or non-finite expected "
                f"count for at least one draw in this chunk, which no Poisson can be drawn "
                f"from. Constrain the prediction to the positive half-line (a Log bijection on "
                f"the norm, or a positive-support prior); ampere.core's own sample() refuses "
                f"the same condition by name.",
            )
        drawn = torch.stack(
            [
                torch.poisson(
                    rates[index], generator=torch.Generator(device="cpu").manual_seed(_seed(seed))
                )
                for index, seed in enumerate(seeds)
            ]
        )
        return drawn.to(device)

    def _student_t_variates(
        self,
        resolved: Mapping[str, torch.Tensor],
        seeds: Sequence[int],
        size: int,
        device: torch.device,
    ) -> torch.Tensor:
        """``mu + sigma t_nu``, one draw at a time from that draw's own seed.

        ``torch.distributions.StudentT`` takes no generator, so the per-draw
        stream comes from ``torch.random.fork_rng`` around a ``manual_seed`` —
        forked so the caller's global RNG state is left exactly as it was
        found, and CPU-only because the variates are drawn on the CPU and moved
        like the normals. It is the public API for the distribution, which is
        the reason it is preferred to the private generator-aware gamma: an
        exact Student-t out of ``torch`` beats a rejection sampler whose
        accuracy would become ampere's problem.
        """
        degrees = resolved["nu"].detach().to("cpu", dtype=DEFAULT_DTYPE)
        variates = []
        with torch.random.fork_rng(devices=[]):
            for index, seed in enumerate(seeds):
                torch.manual_seed(_seed(seed))
                variates.append(torch.distributions.StudentT(degrees[index]).sample((size,)))
        standard = torch.stack(variates).to(device)
        return resolved["location"] + resolved["scale"] * standard

    def _von_mises_variates(
        self,
        resolved: Mapping[str, torch.Tensor],
        seeds: Sequence[int],
        size: int,
        device: torch.device,
    ) -> torch.Tensor:
        """``rng.vonmises(mean, kappa)`` per retained sample, one draw at a time (**W5.2**).

        The same ``torch.random.fork_rng``/``manual_seed`` pattern
        :meth:`_student_t_variates` uses -- ``torch.distributions.VonMises``
        takes no generator either -- but unlike ``nu``, ``mean`` and ``kappa``
        vary **per retained sample** within one draw, the circular analogue of
        the Gaussian family's per-sample ``sigma``. So this distribution's
        batch shape is already ``(size,)`` and a bare ``.sample()`` draws
        exactly one variate per sample, where :meth:`_student_t_variates`'s
        ``.sample((size,))`` instead draws ``size`` i.i.d. variates from one
        shared scalar ``nu``.
        """
        mean = resolved["mean"].detach().to("cpu", dtype=DEFAULT_DTYPE)
        kappa = resolved["kappa"].detach().to("cpu", dtype=DEFAULT_DTYPE)
        variates = []
        with torch.random.fork_rng(devices=[]):
            for index, seed in enumerate(seeds):
                torch.manual_seed(_seed(seed))
                assert mean[index].shape == (size,)  # (chunk, size) after vmap, per (*W3.14*)
                variates.append(torch.distributions.VonMises(mean[index], kappa[index]).sample())
        return torch.stack(variates).to(device)

    def _gp_realisation(self, values: Mapping[str, Any], whitened: torch.Tensor) -> torch.Tensor:
        """``L(θ) z``, from the solver's own native whitening transform."""
        kernel = self._kernel()
        return self.noise.solver.latent_transform_native(
            kernel,
            self.observed_coordinates,
            whitened,
            kernel.resolve(values),
        )

    def _kernel(self) -> Any:
        """This dataset's kernel, bound to its container's axis order (*W4.2*)."""
        assert self.bound_kernel is not None  # only a correlated dataset asks
        return self.bound_kernel

    def log_likelihood(self, routed: Mapping[str, Mapping[str, Any]]) -> torch.Tensor:
        """``log p(data | θ)`` for this dataset alone, as a differentiable scalar.

        Two branches, and which one applies is settled at composition rather
        than here. A correlated noise model with a family that marginalises a
        GP in closed form (the Gaussian one) takes the solver's marginal
        likelihood; everything else takes its family's own body from
        :mod:`ampere.backends.torch._families`, censoring included. The
        combinations that fall in neither — a Student-t with a GP, say — are
        not refused *here* because ``ampere.core.Likelihood`` refuses them
        first, at composition, as an unconsumed latent path.
        """
        predicted = self.predict(routed)
        values = dict(self._dataset_values(routed).get(LIKELIHOOD_COMPONENT, {}))
        # The noise model's own parameters arrive under the likelihood
        # component, flatly -- the same mapping ampere.core hands NoiseModel.
        sigma = self._sigma(predicted, values)
        if self.gp_marginal:
            residual = self.observed_values - predicted
            # W4.2: a complex residual becomes the circular GP's two real
            # columns, and the per-component sigma is the real diagonal of the
            # one covariance both columns are scored against.
            if self.complex_valued:
                residual = torch.stack([residual.real, residual.imag], dim=-1)
            variance = (
                torch.zeros(residual.shape[0], dtype=DEFAULT_DTYPE, device=self.device)
                if sigma is None
                else sigma**2
            )
            kernel = self._kernel()
            value = self.noise.solver.log_marginal_likelihood_native(
                kernel,
                self.observed_coordinates,
                residual,
                variance,
                kernel.resolve(values),
            )
        else:
            if sigma is None and self.likelihood.family.REQUIRES_UNCERTAINTY:
                raise _refuse(
                    "uncertainty",
                    f"dataset {self.label!r} has no observed uncertainties and its noise model "
                    f"declares no jitter, so sigma is undefined.",
                )
            value = native_log_prob(
                FamilyInputs(
                    predicted=predicted,
                    observed=self.observed_values,
                    sigma=sigma,
                    values=values,
                    detection=self.detection,
                    upper=self.upper,
                    lower=self.lower,
                    latent=self._latent(routed, values),
                    family=self.likelihood.family,
                    dtype=DEFAULT_DTYPE,
                    device=self.device,
                )
            )
        return torch.where(torch.isfinite(value), value, torch.full_like(value, -math.inf))


# ---------------------------------------------------------------------------
# W5.9 -- joint noise over a tuple of channels
# ---------------------------------------------------------------------------


class _LoweredJointGroup:
    """One :class:`~ampere.core.JointGaussianProcessNoise` group's term, as torch.

    The twin of ``DatasetCollection.group_log_likelihood`` on the contract
    path, and it exists for one reason: without it, ``B``'s parameters would be
    sampled by NUTS against a density that never saw them. The channels'
    residuals come from the member datasets' own lowered forward chains, so the
    physical model, the instrument steps and the mask rules are shared with
    every other dataset on this backend; what is added here is the rotation and
    the ``T`` rescaled scalar solves.

    The coupling's eigendecomposition is evaluated through
    :meth:`~ampere.core.ChannelCoupling.eigen` with ``xp=torch`` rather than
    transcribed: the parameterisations use only ``cos``, ``sin``, ``exp``,
    ``stack`` and ``linalg.eigh``, which numpy and torch spell identically, so
    there is one implementation for both paths and no transcription to drift.
    The resolved values are coerced to tensors first, because a coupling with
    one parameter fixed and another free hands back a Python float beside a
    tensor and ``torch.stack`` refuses the mixture outright.

    **W5.24: three routes, chosen once at lowering.** The channels' retained
    uncertainties decide, exactly as :meth:`~ampere.core.
    JointGaussianProcessNoise.route` decides on the contract path: identical
    → W5.9's rotated scalar solves, unchanged; unequal under ``DenseGP`` → one
    ``torch.linalg.cholesky_ex`` of the ``TN x TN`` matrix built in torch;
    unequal under ``HilbertSpaceGP`` → Woodbury with the solver's own
    ``(N, m)`` feature factor, read through ``latent_transform_native`` on the
    identity so it carries the gradient in the kernel, and the ``Tm x Tm``
    capacitance. ``B ⊗ K_x`` is never formed on that route. The group's
    ``scale`` and ``jitter`` are group-wide, so they cannot move a group from
    one route to another between lowering and a draw.
    """

    def __init__(
        self,
        problem: FittingProblem,
        label: str,
        lowered: Mapping[str, _LoweredDataset],
        *,
        device: torch.device,
    ) -> None:
        self.label = label
        self.device = device
        self.noise = problem.datasets.joint[label]
        if getattr(self.noise, "BACKEND", None) != BACKEND:
            raise _refuse(
                type(self.noise).__name__,
                f"joint noise group {label!r} is declared by a {type(self.noise).__name__} "
                f"whose BACKEND is {getattr(self.noise, 'BACKEND', None)!r}, not {BACKEND!r}. A "
                f"torch problem whose joint GP solve ran in numpy would not be differentiable "
                f"in B at all; pass ampere.backends.torch.JointGaussianProcessNoise.",
            )
        if not hasattr(self.noise.solver, "log_marginal_likelihood_native"):
            raise _refuse(
                type(self.noise.solver).__name__,
                f"joint noise group {label!r} uses the {type(self.noise.solver).__name__} "
                f"solver, which has no `log_marginal_likelihood_native`. Pass "
                f"ampere.backends.torch.DenseGP, ampere.backends.torch.QuasisepGP or "
                f"ampere.backends.torch.HilbertSpaceGP.",
            )
        self.members = tuple(lowered[member] for member in self.noise.datasets)
        first = self.members[0]
        observed = first.dataset.observed
        self.bound_kernel = self.noise.kernel_for(observed)
        self.coordinates = first.observed_coordinates
        self.size = int(np.count_nonzero(first.retain))
        self.uncertainty = first.uncertainty
        self.uncertainties = tuple(member.uncertainty for member in self.members)
        for member in self.members:
            if member.uncertainty is None:
                continue
            retained = np.asarray(member.dataset.observed.uncertainty, dtype=float).ravel()[
                member.retain
            ]
            if np.any(retained <= 0.0):
                raise _refuse(
                    "uncertainty",
                    f"joint noise group {label!r} was given zero or negative uncertainties on "
                    f"{int(np.sum(retained <= 0.0))} retained sample(s) of channel "
                    f"{member.label!r}. An infinitely precise measurement is one no likelihood "
                    f"can normalise.",
                )
        self.route = self._route()

    def _route(self) -> str:
        """``"rotated"``, ``"dense"`` or ``"reduced_rank"`` — the contract path's rule."""
        present = [uncertainty is not None for uncertainty in self.uncertainties]
        if not any(present):
            return ROTATED_ROUTE
        reference = self.uncertainties[0]
        if all(present) and all(
            bool(torch.equal(uncertainty, reference)) for uncertainty in self.uncertainties
        ):
            return ROTATED_ROUTE
        solver = self.noise.solver
        if solver.NAME == "DenseGP":
            return DENSE_ROUTE
        if solver.NAME in REDUCED_RANK_SOLVERS and hasattr(solver, "latent_transform_native"):
            return REDUCED_RANK_ROUTE
        raise _refuse(
            type(solver).__name__,
            f"joint noise group {self.label!r}'s channels carry different per-sample "
            f"uncertainties, and its {type(solver).__name__} solver has no route for that. "
            f"Pass ampere.backends.torch.DenseGP (the exact dense route) or "
            f"ampere.backends.torch.HilbertSpaceGP (the reduced-rank one).",
        )

    def _sigma_of(
        self, uncertainty: torch.Tensor | None, values: Mapping[str, Any]
    ) -> torch.Tensor | None:
        """One channel's per-sample sigma, with the group's ``scale`` and ``jitter``."""
        own = {key: value for key, value in values.items() if key in self.noise.parameters}
        resolved = self.noise.context(own)
        sigma = uncertainty
        if sigma is None:
            if "jitter" not in resolved:
                return None
            floor = _tensor(resolved["jitter"], device=self.device)
            return floor * torch.ones(self.size, dtype=DEFAULT_DTYPE, device=self.device)
        if "scale" in resolved:
            sigma = sigma * _tensor(resolved["scale"], device=self.device)
        if "jitter" in resolved:
            sigma = torch.sqrt(sigma**2 + _tensor(resolved["jitter"], device=self.device) ** 2)
        return sigma

    def _sigma(self, values: Mapping[str, Any]) -> torch.Tensor | None:
        """The group's shared per-sample sigma, in torch (the rotated route)."""
        return self._sigma_of(self.uncertainty, values)

    def _variances(self, values: Mapping[str, Any]) -> torch.Tensor:
        """Every channel's variances as an ``(n, T)`` block (the two unequal routes)."""
        columns = []
        for uncertainty in self.uncertainties:
            sigma = self._sigma_of(uncertainty, values)
            columns.append(
                torch.zeros(self.size, dtype=DEFAULT_DTYPE, device=self.device)
                if sigma is None
                else sigma**2
            )
        return torch.stack(columns, dim=-1)

    def log_likelihood(self, routed: Mapping[str, Mapping[str, Any]]) -> torch.Tensor:
        """``log N(vec(R); 0, B (x) K_x + blockdiag(diag(sigma_t^2)))``, differentiable."""
        values = dict(routed.get(self.label, {}))
        residuals = torch.stack(
            [member.observed_values - member.predict(routed) for member in self.members],
            dim=-1,
        )
        coupling = self.noise.coupling
        resolved = {
            name: _tensor(value, device=self.device)
            for name, value in coupling.resolved(values).items()
        }
        eigenvalues, rotation = coupling.eigen(resolved, xp=torch)
        kernel = self.bound_kernel
        hyperparameters = kernel.resolve(values)
        if self.route == DENSE_ROUTE:
            total = self._dense(residuals, self._variances(values), eigenvalues, rotation, values)
        elif self.route == REDUCED_RANK_ROUTE:
            total = self._reduced_rank(
                residuals, self._variances(values), eigenvalues, rotation, hyperparameters
            )
        else:
            sigma = self._sigma(values)
            variance = (
                torch.zeros(self.size, dtype=DEFAULT_DTYPE, device=self.device)
                if sigma is None
                else sigma**2
            )
            rotated = residuals @ rotation
            total = _tensor(0.0, device=self.device)
            for index in range(len(self.members)):
                scale = eigenvalues[index]
                total = total + self.noise.solver.log_marginal_likelihood_native(
                    kernel,
                    self.coordinates,
                    rotated[:, index] / torch.sqrt(scale),
                    variance / scale,
                    hyperparameters,
                )
                total = total - 0.5 * self.size * torch.log(scale)
        return torch.where(torch.isfinite(total), total, torch.full_like(total, -math.inf))

    def _dense(
        self,
        residuals: torch.Tensor,
        variances: torch.Tensor,
        eigenvalues: torch.Tensor,
        rotation: torch.Tensor,
        values: Mapping[str, Any],
    ) -> torch.Tensor:
        """The dense route: ``B (x) (K_x + jitter^2 I) + blockdiag(diag(sigma_t^2))``."""
        kernel = self.bound_kernel
        matrix = _tensor(
            kernel.matrix(self.coordinates, self.coordinates, kernel.resolve(values)),
            device=self.device,
        )
        jitter = float(getattr(self.noise.solver, "jitter", 0.0) or 0.0)
        if jitter:
            matrix = matrix + torch.eye(self.size, dtype=DEFAULT_DTYPE, device=self.device) * (
                jitter**2
            )
        coupling = (rotation * eigenvalues[None, :]) @ rotation.transpose(0, 1)
        covariance = torch.kron(coupling, matrix) + torch.diag(
            variances.transpose(0, 1).reshape(-1)
        )
        residual = residuals.transpose(0, 1).reshape(-1)
        factor, info = torch.linalg.cholesky_ex(covariance)
        alpha = torch.cholesky_solve(residual.reshape(-1, 1), factor, upper=False).reshape(-1)
        log_determinant = 2.0 * torch.log(torch.diagonal(factor)).sum()
        value = -0.5 * ((residual * alpha).sum() + log_determinant + residual.numel() * _LOG_2PI)
        return torch.where(info != 0, torch.full_like(value, -math.inf), value)

    def _reduced_rank(
        self,
        residuals: torch.Tensor,
        variances: torch.Tensor,
        eigenvalues: torch.Tensor,
        rotation: torch.Tensor,
        hyperparameters: Mapping[str, Any],
    ) -> torch.Tensor:
        """The reduced-rank route: Woodbury with ``G = (I_T (x) Phi)(F (x) I_m)``, ``B = F F^T``.

        See :meth:`ampere.core.JointGaussianProcessNoise._reduced_rank_log_prob`
        for the algebra; this is the same arithmetic in torch.
        """
        solver = self.noise.solver
        kernel = self.bound_kernel
        count = int(solver.latent_size(kernel, self.size))
        identity = torch.eye(count, dtype=DEFAULT_DTYPE, device=self.device)
        phi = solver.latent_transform_native(kernel, self.coordinates, identity, hyperparameters)
        channels = len(self.members)
        width = channels * count
        diagonal = variances + float(getattr(solver, "jitter", 0.0) or 0.0) ** 2
        inverse = 1.0 / diagonal
        root = rotation * torch.sqrt(eigenvalues)[None, :]
        gram = torch.einsum("nk,nt,nl->tkl", phi, inverse, phi)
        capacitance = torch.einsum("ts,tu,tkl->skul", root, root, gram).reshape(
            width, width
        ) + torch.eye(width, dtype=DEFAULT_DTYPE, device=self.device)
        weighted = residuals * inverse
        projected = ((phi.transpose(0, 1) @ weighted) @ root).transpose(0, 1).reshape(-1, 1)
        factor, info = torch.linalg.cholesky_ex(capacitance)
        solved = torch.cholesky_solve(projected, factor, upper=False)
        log_determinant = torch.log(diagonal).sum() + 2.0 * torch.log(torch.diagonal(factor)).sum()
        quadratic = (residuals * weighted).sum() - (projected * solved).sum()
        value = -0.5 * (quadratic + log_determinant + residuals.numel() * _LOG_2PI)
        failed = torch.logical_or(info != 0, torch.logical_not(torch.all(diagonal > 0.0)))
        return torch.where(failed, torch.full_like(value, -math.inf), value)


class LoweredProblem:
    """A :class:`~ampere.core.FittingProblem`, lowered onto torch.

    Built by :func:`lower_problem`, and registered with ``ampere.core`` as
    this backend's realisation (``inference.md`` §10a) when
    ``ampere.backends.torch`` is imported. Exposes the same three quantities
    §4.5 does — a log prior, a log likelihood and their sum — as
    **differentiable functions of tensors**, plus the unconstrained form a NUTS
    kernel needs.

    Attributes
    ----------
    problem
        The problem this lowers. Kept so a driver can read the seed, the
        parameter names and the capability flags off it rather than
        re-deriving them.
    parameters
        The :class:`~ampere.backends.torch.parameters.TorchParameterSpace` —
        the ``torch.distributions``/``biject_to`` view of the merged
        declaration.
    """

    def __init__(self, problem: FittingProblem) -> None:
        if problem.backend != BACKEND:
            raise _refuse(
                problem.backend,
                f"this problem reports backend {problem.backend!r}, not {BACKEND!r}. The backend "
                f"is derived from what the models, transformations, noise models and GP solvers "
                f"declare (W2.12, widened at W2.13), so a problem built from another backend's "
                f"pieces cannot be lowered here — build it from ampere.backends.torch's.",
            )
        # W3.8: a problem that reports this backend but composes a piece from
        # another one. That composition is refused at construction unless it
        # asked for allow_foreign_parts=True, and the flag buys a gradient-free
        # run and nothing more -- a lowered problem *is* the differentiable
        # form, so it refuses by name here regardless.
        foreign = foreign_parts_refusal(problem, what="a differentiable torch problem")
        if foreign is not None:
            raise foreign
        self.problem = problem
        # Where this realisation lives (W2.4 slice 3). ``problem.device`` is
        # ``declared_capabilities``' aggregate over every part, which refuses a
        # problem whose pieces name two devices -- so by the time we are here
        # there is exactly one answer and every tensor below is built on it.
        self.device = resolve_device(problem.device)
        self.parameters = TorchParameterSpace(
            problem.parameters, strict=problem.strict, device=self.device
        )
        self._mapping = problem.mapping
        self._datasets = tuple(
            _LoweredDataset(problem, label, device=self.device) for label in problem.datasets
        )
        # W5.9: a dataset a joint noise group claims contributes its residual to
        # the group's one term, not a term of its own -- `inference.md` section 4's
        # "joint" decomposition, and the reason `_likelihood_terms` is not
        # simply one entry per dataset any more.
        lowered = {dataset.label: dataset for dataset in self._datasets}
        self._joint = tuple(
            _LoweredJointGroup(problem, label, lowered, device=self.device)
            for label in problem.datasets.joint
        )
        self._grouped = frozenset(
            member for group in self._joint for member in group.noise.datasets
        )

    @property
    def backend(self) -> str:
        """This backend's one name — what ``ampere.core.realise`` checks against the problem's."""
        return BACKEND

    @property
    def free_size(self) -> int:
        return self.parameters.free_size

    def lowering_provenance(self) -> list[dict[str, Any]]:
        """The user-registered lowering rows this problem consulted.

        ``lowering.md`` §12.8's "every registered row is stamped
        user-registered in provenance", ready for ``provenance_attrs``'s
        ``extra=``. Empty for a problem whose priors all took ampere's own
        built-in rows, which is the point: the signal is "did this run depend
        on something outside the conformance suite's guarantees?".
        """
        return self.parameters.lowering_provenance()

    def _route(self, theta: torch.Tensor) -> dict[str, dict[str, Any]]:
        """The constrained free vector, routed to each component's local names.

        ``ParameterMapping.distribute`` is pure routing metadata — it looks a
        value up by name and puts it where the declaration says it goes — so it
        carries tensors through untouched, graph and all, and reusing it is
        what keeps this module from acquiring its own (divergent) idea of how a
        tie is bound.
        """
        return self._mapping.distribute(self.parameters.unpack_tensor(theta))

    def log_prior(self, theta: Any) -> torch.Tensor:
        """``log p(θ)`` in the constrained parameterisation."""
        return self.parameters.log_prior_tensor(theta)

    def _likelihood_terms(self, theta: Any) -> dict[str, torch.Tensor]:
        """Each dataset's own ``log p(data | θ)``, in the **constrained** space."""
        routed = self._route(self.parameters._tensor(theta))
        terms = {
            dataset.label: dataset.log_likelihood(routed)
            for dataset in self._datasets
            if dataset.label not in self._grouped
        }
        for group in self._joint:
            terms[group.label] = group.log_likelihood(routed)
        return terms

    def log_likelihood(self, theta: Any) -> torch.Tensor:
        """``log p(data | θ)``, summed over datasets."""
        total = _tensor(0.0, device=self.device)
        for value in self._likelihood_terms(theta).values():
            total = total + value
        return torch.where(torch.isfinite(total), total, torch.full_like(total, -math.inf))

    def log_likelihood_terms(self, unconstrained: Any) -> dict[str, torch.Tensor]:
        """``inference.md`` §10a's optional member: the per-dataset decomposition.

        Keyed by dataset label, in the **unconstrained** parameterisation — the
        same argument :meth:`log_prob_unconstrained` takes, because that is the
        vector a NUTS driver has in hand for every draw. The terms sum to
        :meth:`log_likelihood` at the corresponding constrained point, and that
        identity is what ``ampere.results.emit``'s per-dataset
        ``log_likelihood`` group records.

        Supplying it is optional by ruling; supplying it costs nothing here,
        because the joint likelihood is already a sum over exactly these terms.
        """
        y = self.parameters._tensor(unconstrained)
        return self._likelihood_terms(self.parameters.constrain_tensor(y))

    def log_prob(self, theta: Any) -> torch.Tensor:
        """``log p(θ) + log p(data | θ)`` — the number a sampler maximises."""
        total = self.log_prior(theta) + self.log_likelihood(theta)
        return torch.where(torch.isfinite(total), total, torch.full_like(total, -math.inf))

    def log_prob_unconstrained(self, unconstrained: Any) -> torch.Tensor:
        """:meth:`log_prob` in unconstrained space, change of variables included.

        The density a NUTS kernel works with, and the torch counterpart of
        ``FittingProblem.log_prob_unconstrained``. The Jacobian term comes from
        :meth:`~ampere.backends.torch.parameters.TorchParameterSpace.log_prior_unconstrained_tensor`,
        which is the one place in this backend that calls
        ``log_abs_det_jacobian``, so the argument order of ``lowering.md``
        §2(b) is got right once.
        """
        y = self.parameters._tensor(unconstrained)
        prior = self.parameters.log_prior_unconstrained_tensor(y)
        constrained = self.parameters.constrain_tensor(y)
        total = prior + self.log_likelihood(constrained)
        return torch.where(torch.isfinite(total), total, torch.full_like(total, -math.inf))

    def log_prob_unconstrained_batched(self, unconstrained: Any) -> torch.Tensor:
        """:meth:`log_prob_unconstrained` over a **stack** of free vectors.

        ``(batch, free_size)`` in, ``(batch,)`` out, in one call, through
        :func:`torch.func.vmap` — W2.4 slice 2's ``BATCHABLE`` claim, and the
        surface that makes it a claim about code rather than about intent.

        It is deliberately **not** a loop wearing a batched name. ``vmap``
        rewrites the whole density as batched operations: one Cholesky over a
        stack of covariances, one matrix multiply through the instrument chain,
        one set of ``log_prob`` calls over stacked distribution parameters. A
        loop would produce the same numbers and none of the benefit, and it
        would make the capability flag a falsehood — which is why this method
        refuses rather than falling back when some part of the problem is not
        batchable.

        Two conditions must hold and both are checked here rather than
        discovered inside torch:

        * **every part must declare** ``BATCHABLE``. That is
          ``problem.batchable``, aggregated conjunctively by
          :func:`~ampere.core.declared_capabilities`, and the refusal names the
          pieces that said no. The commonest one is
          :class:`~ampere.backends.torch.QuasisepGP`, whose solve is a compiled
          extension ``vmap`` cannot see through;
        * the argument must be two-dimensional. A single vector is
          :meth:`log_prob_unconstrained`'s job, and silently accepting one here
          would make the two methods' shapes depend on the input rather than on
          which was called.

        What makes it work at all is that this density contains **no
        data-dependent control flow**: every refusal is at construction, every
        failure is a :func:`torch.where`, and W2.4 slice 2 removed the last
        Python ``if`` on a tensor value — the prior's ``math.isfinite`` short
        circuit in
        :meth:`~ampere.backends.torch.parameters.TorchParameterSpace.log_prior_tensor`.
        That was already the rule this module was written to
        (``inference.md`` §10a, ``likelihoods.md`` §17 Q1); batching is the
        second thing it buys, after tracing.

        Examples
        --------
        >>> import numpy as np, scipy.stats as st, astropy.units as u, torch
        >>> from ampere.backends.torch import IndependentNoise, PowerLaw, lower_problem
        >>> from ampere.core import (
        ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
        ... )
        >>> grid = np.array([1.0, 2.0, 3.0])
        >>> observed = Spectrum(
        ...     grid * u.um, [2.0, 1.0, 0.7] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
        ... )
        >>> problem = FittingProblem(
        ...     PowerLaw(grid, norm=st.lognorm(0.3, scale=2.0), index=-1.0),
        ...     [Dataset(observed, likelihood=Likelihood(GaussianFamily(), IndependentNoise()))],
        ... )
        >>> lowered = lower_problem(problem)
        >>> stack = torch.zeros((4, lowered.free_size), dtype=torch.float64)
        >>> lowered.log_prob_unconstrained_batched(stack).shape
        torch.Size([4])

        and it agrees, point for point, with evaluating them one at a time:

        >>> stack = torch.linspace(-1.0, 1.0, 4).reshape(4, 1).to(torch.float64)
        >>> one_at_a_time = torch.stack([lowered.log_prob_unconstrained(y) for y in stack])
        >>> bool(torch.allclose(lowered.log_prob_unconstrained_batched(stack), one_at_a_time))
        True
        """
        if not self.problem.batchable:
            parts = ", ".join(
                sorted(
                    {
                        type(part).__name__
                        for part in (
                            *self.problem.models.values(),
                            *self.problem.datasets.capability_parts,
                        )
                        if not bool(getattr(part, "BATCHABLE", False))
                    }
                )
            )
            raise _refuse(
                "batching",
                f"this problem declares batchable=False, so its density cannot be evaluated "
                f"over a stack of parameter vectors in one call. The piece(s) that cannot: "
                f"{parts or '(none declared)'}. torch.func.vmap needs every operation in the "
                f"density to be one it can rewrite, and a compiled extension reached through a "
                f"torch.autograd.Function is not — QuasisepGP is the usual answer here, and "
                f"DenseGP is the batchable alternative. Evaluate the vectors one at a time with "
                f"log_prob_unconstrained.",
            )
        stack = as_tensor(unconstrained, dtype=DEFAULT_DTYPE, device=self.device)
        if stack.ndim != 2 or int(stack.shape[1]) != self.free_size:
            raise _refuse(
                "batching",
                f"log_prob_unconstrained_batched takes a (batch, {self.free_size}) stack, got "
                f"shape {tuple(stack.shape)}. A single vector belongs to "
                f"log_prob_unconstrained.",
            )
        return torch.func.vmap(self.log_prob_unconstrained)(stack)

    # -- batched simulation (W3.1 slice 2) ----------------------------------

    def _forward(self, theta: torch.Tensor) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
        """One θ through the whole noise-free forward model, natively.

        Returns what :class:`~ampere.core.BatchedPrediction` holds, as nested
        dictionaries of tensors so ``torch.func.vmap`` maps them in one call:
        every model's every channel, and every dataset's instrument-transformed
        prediction on the *whole* observed grid.

        Every channel, not only the bound ones, because ``results.md`` §11
        writes one training-set group per ``<model>.<channel>`` and a fast path
        that returned fewer would silently write a smaller file than the loop.
        """
        routed = self._route(self.parameters._tensor(theta))
        channels = {
            label: {
                channel: _flat_channel(_require_flux(model, label)(channel, routed.get(label, {})))
                for channel in getattr(model, "channels", ())
            }
            for label, model in self.problem.models.items()
        }
        predicted = {dataset.label: dataset.predict_full(routed) for dataset in self._datasets}
        return channels, predicted

    def simulate_batched(
        self,
        theta: Any,
        *,
        chunk_size: int | None = None,
        sharder: Any = None,
    ) -> BatchedPrediction:
        """``inference.md`` §13's *batched form*, natively: a stack of θ, one ``vmap``.

        ``(batch, free_size)`` **constrained** free vectors in — the coordinates
        ``SimulationBatch.theta`` holds, not the unconstrained ones
        :meth:`log_prob_unconstrained` takes, because a simulation budget is a
        set of parameter values rather than a set of sampler positions.

        **Per chunk, always.** ``chunk_size`` bounds how many simulations are
        vectorised at once and the chunks are *looped*; a ``vmap`` over a whole
        budget is exactly the single-device memory trap Peter's ruling of
        2026-09-08 names, and the reason ``simulate_many`` is built on an
        executor rather than on this. ``None`` means one chunk, which is what
        ``simulate_many`` passes because it has already chunked.

        **What it does not do.** It draws no noise
        (:meth:`sample_observations` does), it flags no failures (there is no
        control flow inside a ``vmap`` to flag with, so a simulator that fails
        produces NaNs, classified by the caller exactly as the loop classifies
        them), and it refuses rather than falling back — ``simulate_many`` owns
        the fallback, and owning it in two places would let a fast path quietly
        become a slow one.

        Refused, by name, when the problem declares ``batchable = False``, for
        the same reason :meth:`log_prob_unconstrained_batched` refuses: a
        compiled extension reached through a ``torch.autograd.Function`` is not
        something ``vmap`` can rewrite.
        """
        if not self.problem.batchable:
            raise _refuse(
                "batched simulation",
                f"this problem declares batchable=False, so its forward model cannot be "
                f"evaluated over a stack of parameter vectors in one call "
                f"({self.problem.capabilities}). torch.func.vmap needs every operation to be "
                f"one it can rewrite, and a compiled extension reached through a "
                f"torch.autograd.Function is not — QuasisepGP is the usual answer here. "
                f"simulate_many falls back to the loop, which is the semantics anyway.",
            )
        stack = as_tensor(theta, dtype=DEFAULT_DTYPE, device=self.device)
        if stack.ndim != 2 or int(stack.shape[1]) != self.free_size:
            raise _refuse(
                "batched simulation",
                f"simulate_batched takes a (batch, {self.free_size}) stack of constrained free "
                f"vectors, got shape {tuple(stack.shape)}.",
            )
        vectorised = torch.func.vmap(self._forward)
        channels: dict[str, dict[str, list[np.ndarray]]] = {}
        predicted: dict[str, list[np.ndarray]] = {}
        for start, stop in chunk_bounds(int(stack.shape[0]), chunk_size):
            piece = stack[start:stop]
            produced, prediction = (
                vectorised(piece) if sharder is None else sharder.shard(vectorised, piece)
            )
            for label, holding in produced.items():
                for channel, values in holding.items():
                    channels.setdefault(label, {}).setdefault(channel, []).append(to_numpy(values))
            for label, values in prediction.items():
                predicted.setdefault(label, []).append(to_numpy(values))
        return BatchedPrediction(
            channels={
                label: {channel: np.concatenate(parts) for channel, parts in holding.items()}
                for label, holding in channels.items()
            },
            predicted={label: np.concatenate(parts) for label, parts in predicted.items()},
        )

    def sample_observations(
        self,
        theta: Any,
        predicted: Mapping[str, Any],
        seeds: Sequence[int],
        *,
        sigma: Mapping[str, Any] | None = None,
    ) -> dict[str, np.ndarray]:
        """Draw the retained observed values for a chunk, natively.

        Peter's ruling of 2026-09-08. The distributions are ``ampere.core``'s
        own families' ``sample`` methods, transcribed in
        :meth:`_LoweredDataset.draw_chunk`; the **stream** is
        ``torch.Generator``'s, which is why the numpy path stays the oracle and
        the two are compared *distributionally* rather than draw for draw.

        *sigma* (**W5.29**) carries a per-draw observation context: dataset
        label to a ``(batch,) + observed.shape`` stack of sigma arrays in the
        observed container's value unit, one per draw — what
        ``Dataset.contextual_observed`` substitutes on the contract path. A
        label it omits draws at the observation's own sigma, as does
        ``sigma=None``. The context's sigma replaces ``sigma_data`` *before*
        the noise model's ``scale``, ``jitter`` and fractional inflation are
        applied, exactly as the contract path's ``NoiseParams`` are built from
        the substituted container.

        *seeds* is one integer per draw, taken from the per-draw child generator
        ``simulate_many`` spawns **by index** — so partition independence
        carries over unchanged: draw *i* gets the same generator whichever chunk
        it ran in, and a budget split 1/7/whole gives the same observations.

        Every variate is drawn **outside** the ``vmap``, per draw, from that
        draw's own seed. ``torch.func.vmap`` has no per-sample random state to
        give a transformed function, and faking one with a global generator
        would make a draw depend on how the chunk was scheduled — which is the
        property this whole path is built to keep. For the Gaussian and complex
        families that means two standard normals mapped in as data; for the
        Poisson and Student-t twins W3.14 added it means the transform computes
        the draw's distribution parameters and the variate follows (see
        :meth:`_LoweredDataset.draw_chunk`).

        Refuses by name, before drawing anything, for any dataset this backend
        may not sample (:meth:`_LoweredDataset.sampling_refusal`); the caller
        then runs the numpy path, where ``ampere.core``'s own refusal text is
        what a user meets.
        """
        if self._joint:
            # W5.9. A joint noise group's channels are correlated, and this
            # method draws dataset by dataset; drawing each channel from its
            # own marginal would write a training set whose cross-covariance is
            # zero -- data from a different model than the one being fitted.
            # The numpy path draws the group correlated
            # (``DatasetCollection.draw_group``) and is what the caller falls
            # back to.
            raise _refuse(
                "joint",
                f"this problem declares joint noise group(s) "
                f"{sorted(group.label for group in self._joint)}, whose channels are correlated "
                f"with one another. The native draw is per dataset, so it would produce "
                f"observations with no cross-covariance at all; the contract path draws the "
                f"group in one correlated call.",
            )
        for dataset in self._datasets:
            refusal = dataset.sampling_refusal()
            if refusal is not None:
                raise refusal
        stack = as_tensor(theta, dtype=DEFAULT_DTYPE, device=self.device)
        drawn: dict[str, np.ndarray] = {}
        for dataset in self._datasets:
            # The prediction's own dtype decides, as it does at lowering time
            # (W2.4 slice 3): a complex_gaussian dataset predicts complex
            # visibilities, and forcing them to float64 here would not fail
            # loudly -- it would draw around the real projection of the model.
            raw = np.asarray(predicted[dataset.label])
            rows = as_tensor(
                raw,
                dtype=complex_dtype(DEFAULT_DTYPE) if raw.dtype.kind == "c" else DEFAULT_DTYPE,
                device=self.device,
            )
            context = None if sigma is None else sigma.get(dataset.label)
            uncertainty = (
                None
                if context is None
                else _tensor(_context_rows(context, dataset, len(seeds)), device=self.device)
            )
            drawn[dataset.label] = to_numpy(
                dataset.draw_chunk(
                    self._route,
                    stack,
                    rows[:, dataset.retain],
                    seeds,
                    self.device,
                    uncertainty=uncertainty,
                )
            )
        return drawn

    def potential(self) -> Callable[[torch.Tensor], torch.Tensor]:
        """The **negative** unconstrained log-density, which is what pyro wants.

        pyro's kernels minimise a potential energy rather than maximise a
        log-density, and getting that sign wrong produces a sampler that
        explores the prior's tails with perfect efficiency. Spelled once, here.
        """

        def potential_energy(unconstrained: torch.Tensor) -> torch.Tensor:
            return -self.log_prob_unconstrained(unconstrained)

        return potential_energy

    def __repr__(self) -> str:
        return (
            f"<LoweredProblem {self.free_size} free dimension(s), "
            f"{len(self._datasets)} dataset(s), backend={BACKEND!r}, "
            f"device={str(self.device)!r}>"
        )


def lower_problem(problem: FittingProblem) -> LoweredProblem:
    """Lower *problem* onto torch: one differentiable log-density, built once.

    ``lowering.md`` §0: "It happens once, when a ``FittingProblem`` is realised
    on a backend — never per evaluation." This is the factory registered as
    this backend's realisation, so ``ampere.core.realise(problem)`` on a torch
    problem arrives here.

    Raises
    ------
    LoweringError
        Naming, by class, whatever part of the problem this backend cannot
        compose natively — a model without a torch surface, a censored
        dataset, a non-Gaussian family, another backend's GP solver. Never a
        silent fallback to the reference path: a nominally torch-backed NUTS
        run that quietly lost its gradients would be a run that cannot work at
        all.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.torch import IndependentNoise, PowerLaw, lower_problem
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
    >>> grid = np.array([1.0, 2.0, 3.0])
    >>> observed = Spectrum(
    ...     grid * u.um, [2.0, 1.0, 0.7] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
    ... )
    >>> # This backend's noise model, not ampere.core's: since W2.13 a noise
    >>> # model declares a backend, and Dataset's default declares "reference".
    >>> likelihood = Likelihood(GaussianFamily(), IndependentNoise())
    >>> problem = FittingProblem(
    ...     PowerLaw(grid, norm=st.lognorm(0.3, scale=2.0), index=-1.0),
    ...     [Dataset(observed, likelihood=likelihood)],
    ... )
    >>> lowered = lower_problem(problem)
    >>> y = np.zeros(lowered.free_size)
    >>> bool(np.isfinite(float(lowered.log_prob_unconstrained(y))))
    True

    and it is differentiable, which is the whole point:

    >>> import torch
    >>> theta = torch.zeros(lowered.free_size, dtype=torch.float64, requires_grad=True)
    >>> lowered.log_prob_unconstrained(theta).backward()
    >>> theta.grad.shape == (lowered.free_size,)
    True

    and it supplies §10a's optional per-dataset decomposition, which sums to
    the joint log-likelihood:

    >>> terms = lowered.log_likelihood_terms(y)
    >>> sorted(terms)
    ['default']
    >>> float(sum(terms.values())) == float(
    ...     lowered.log_likelihood(lowered.parameters.constrain_tensor(y))
    ... )
    True
    """
    return LoweredProblem(problem)
