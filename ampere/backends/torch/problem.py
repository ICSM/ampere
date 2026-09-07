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
surfaces: ``model.flux``, ``step.apply_flux``, ``noise.sigma_tensor`` where a
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
censored dataset, a latent GP, a non-Gaussian family, another backend's GP
solver — is made **at construction**, once, with a
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
  dormant.

Two things stay refused, both by name and both for a stated reason rather
than for want of time. ``complex_gaussian`` needs complex tensors end to end
and belongs with the visibility modality (plan §5, Phase 4). **The latent-GP
path is refused because the numpy path it must agree with does not apply the
whitening transform** — see the refusal's own message, which states the
defect and its measurable consequence; lowering it would mean either copying
the defect or disagreeing with the oracle, and neither is a backend's call.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import numpy as np
import torch

from ampere.core import (
    FittingProblem,
    GaussianProcessNoise,
    IndependentNoise,
)
from ampere.core.dataset import INSTRUMENT_COMPONENT, LIKELIHOOD_COMPONENT
from ampere.core.exceptions import LoweringError

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor
from ._families import FamilyInputs, limit_masks, native_log_prob, refuse_family
from .parameters import TorchParameterSpace

__all__ = ["LoweredProblem", "lower_problem"]

_LOG_2PI = math.log(2.0 * math.pi)


def _refuse(what: str, detail: str) -> LoweringError:
    return LoweringError(what, backend=BACKEND, detail=detail)


def _tensor(value: Any) -> torch.Tensor:
    return as_tensor(value, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)


class _LoweredDataset:
    """One dataset's forward chain and log-likelihood, as torch."""

    def __init__(self, problem: FittingProblem, label: str) -> None:
        dataset = problem.datasets[label]
        self.label = label
        self.model_label = problem.bindings[label]
        self.model = problem.models[self.model_label]
        self.channel = dataset.instrument.channel
        self.dataset = dataset
        self.likelihood = dataset.likelihood
        self.noise = dataset.likelihood.noise
        self.steps = tuple(dataset.instrument.steps)

        if not hasattr(self.model, "flux"):
            raise _refuse(
                type(self.model).__name__,
                f"model {self.model_label!r} has no native torch surface (a `flux` method), so it "
                f"cannot be composed into a differentiable log-density. Build the problem from "
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
        )
        if refusal is not None:
            raise refusal
        if dataset.latent is not None:
            raise _refuse(
                "latent",
                f"dataset {label!r} declares a latent GP, and this backend refuses to lower it "
                f"**because the numpy path it would have to agree with does not apply the "
                f"whitening transform**. ampere.core.likelihood.latent_parameter declares that "
                f"'the covariance enters through f = L(theta) z, a deterministic transform owned "
                f"by the GPSolver', but no scoring path calls GPSolver.latent_transform: "
                f"Dataset.log_likelihood_of hands the whitened z straight to "
                f"Likelihood.log_prob, which hands it to the family as noise.latent, and "
                f"PoissonFamily.log_prob uses it as f. The measurable consequence is that the "
                f"kernel hyperparameters do not enter the likelihood at all -- the value is "
                f"identical for amplitude 0.5, 5 and 50 -- so a latent fit samples them against "
                f"a flat likelihood. Lowering this here would mean either reproducing the defect "
                f"or disagreeing with the oracle ampere.core.realise checks against, and "
                f"neither is a backend's decision to take. Recorded as a W2.4 slice 2 finding; "
                f"the fix belongs in ampere.core.",
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
        if not isinstance(self.noise, (IndependentNoise, GaussianProcessNoise)):
            raise _refuse(
                type(self.noise).__name__,
                f"dataset {label!r} uses a noise model this backend does not know how to "
                f"compose natively.",
            )
        self.correlated = correlated

        observed = dataset.observed
        # The effective mask, resolved once at composition: Dataset takes the
        # union of the observed and predicted masks there and forbids a
        # parameter-dependent one, so which samples are retained is a fact
        # about the problem rather than about theta and can be a constant here.
        excluded = dataset.effective_mask
        if excluded is None:
            excluded = ~np.asarray(observed.valid, dtype=bool).ravel()
        self.retain = ~np.asarray(excluded, dtype=bool).ravel()
        self.observed_values = _tensor(np.asarray(observed.values, dtype=float)[self.retain])
        self.observed_coordinates = _tensor(
            np.asarray(observed.axes[0].values, dtype=float)[self.retain]
        )
        self.uncertainty = (
            None
            if observed.uncertainty is None
            else _tensor(np.asarray(observed.uncertainty, dtype=float)[self.retain])
        )
        # Which samples are limits is a fact about the data, so the three
        # groups are resolved once here rather than per evaluation. A limit on
        # a masked sample is not a limit at all -- likelihoods.md §9's "masking
        # beats censoring" -- so the declaration is read *through* the
        # effective mask, exactly as Likelihood._retained_limits does.
        self.detection, self.upper, self.lower = (None, None, None)
        if censoring is not None:
            censoring.check_against(observed)
            self.detection, self.upper, self.lower = limit_masks(
                np.asarray(censoring.kinds)[self.retain]
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
        """
        flux = self.model.flux(self.channel, routed[self.model_label])
        grid = self.model.grid(self.channel)
        for step in self.steps:
            flux, grid = step.apply_flux(flux, grid, self._step_values(routed, step.label))
        return flux[self.retain]

    # -- the log-likelihood -------------------------------------------------

    def _base_sigma(self, values: Mapping[str, Any]) -> torch.Tensor | None:
        """``sqrt((scale * sigma_data)**2 + jitter**2)`` for the retained samples.

        Transcribed rather than called, because ``NoiseModel.sigma`` coerces
        with ``float()`` and ``scale`` and ``jitter`` are fitted parameters:
        calling it would cut their gradient exactly as calling
        ``Kernel.matrix`` used to cut the GP amplitude's (W2.4 slice 1's
        finding). ``None`` when there are no uncertainties and no jitter — a
        family that cannot live without them has already refused at
        composition through ``NoiseModel.check_compatible``.
        """
        own = {key: value for key, value in values.items() if key in self.noise.parameters}
        resolved = self.noise.context(own)
        sigma = self.uncertainty
        if sigma is None:
            if "jitter" not in resolved:
                return None
            return _tensor(resolved["jitter"]).expand(int(self.retain.sum())).clone()
        if "scale" in resolved:
            sigma = sigma * _tensor(resolved["scale"])
        if "jitter" in resolved:
            floor = _tensor(resolved["jitter"])
            sigma = torch.sqrt(sigma**2 + floor**2)
        return sigma

    def _sigma(self, predicted: torch.Tensor, values: Mapping[str, Any]) -> torch.Tensor | None:
        """The effective ``sigma``, prediction-aware noise models included.

        The base quadrature is computed here (once, from the cached
        uncertainties) and *handed* to a noise model that declares
        ``sigma_tensor`` — this backend's ``FractionalModelNoise`` and
        ``FractionalModelGPNoise``, which inflate it by ``f * |predicted|``.
        The hook takes the base rather than the container and the mask so that
        the quadrature is written once and every parameter in it keeps its
        graph; see :mod:`ampere.backends.torch.noise`.
        """
        base = self._base_sigma(values)
        native = getattr(self.noise, "sigma_tensor", None)
        if native is None:
            return base
        return native(base, values, predicted=predicted)

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
        if self.correlated:
            residual = self.observed_values - predicted
            variance = torch.zeros_like(residual) if sigma is None else sigma**2
            value = self.noise.solver.log_marginal_likelihood_native(
                self.noise.kernel,
                self.observed_coordinates,
                residual,
                variance,
                self.noise.kernel.resolve(values),
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
                    latent=None,
                    family=self.likelihood.family,
                )
            )
        return torch.where(torch.isfinite(value), value, torch.full_like(value, -math.inf))


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
        self.problem = problem
        self.parameters = TorchParameterSpace(problem.parameters, strict=problem.strict)
        self._mapping = problem.mapping
        self._datasets = tuple(_LoweredDataset(problem, label) for label in problem.datasets)

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
        return {dataset.label: dataset.log_likelihood(routed) for dataset in self._datasets}

    def log_likelihood(self, theta: Any) -> torch.Tensor:
        """``log p(data | θ)``, summed over datasets."""
        total = _tensor(0.0)
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
        stack = as_tensor(unconstrained, dtype=DEFAULT_DTYPE, device=DEFAULT_DEVICE)
        if stack.ndim != 2 or int(stack.shape[1]) != self.free_size:
            raise _refuse(
                "batching",
                f"log_prob_unconstrained_batched takes a (batch, {self.free_size}) stack, got "
                f"shape {tuple(stack.shape)}. A single vector belongs to "
                f"log_prob_unconstrained.",
            )
        return torch.func.vmap(self.log_prob_unconstrained)(stack)

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
            f"{len(self._datasets)} dataset(s), backend={BACKEND!r}>"
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
