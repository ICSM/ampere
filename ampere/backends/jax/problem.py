"""The native path: a composed ``FittingProblem`` as one differentiable function.

Why this module has to exist
----------------------------
``DEVELOPMENT_PLAN.md`` §4.5's engine-facing surface —
``FittingProblem.log_prob``, ``log_prob_unconstrained``, ``prior_transform`` —
is the one an engine is written against, and ``ampere.inference``'s three
gradient-free drivers consume nothing else. A NUTS driver cannot: that surface
is **not traceable**, and the reason is structural rather than incidental.

Two links in the chain break a trace, and both are frozen §4 contract:

1. ``ampere.core.results_schema``'s containers convert their values with
   ``numpy.asarray`` (``_as_array``, and ``with_values`` through it). A
   ``ModelResult`` is made of those containers, so a model that computes its
   flux in jax hands the likelihood a *numpy* array and the gradient stops at
   that boundary. Under ``jax.grad`` it does not merely stop, it raises:
   ``np.asarray`` of a tracer is a ``TracerArrayConversionError``.
2. ``ParameterSet.lnprior``/``lnprior_unconstrained`` short-circuit on
   ``math.isfinite`` — Python control flow on a value — and coerce with
   ``np.asarray(..., dtype=float)``.

Neither is a defect: the containers are what make a heterogeneous joint fit
checkable, and the short circuit is what stops an expensive model being
evaluated at a point with no prior mass. But together they mean a gradient
cannot be taken *through the contract path*, on any backend, and a jax backend
that pretended otherwise would be claiming a capability it does not have.

So this module composes the same quantity out of the backend's **native**
surfaces — ``model.flux``, ``step.apply_flux``, ``noise.sigma_jax``,
``DenseGP.log_marginal_likelihood_jax``, ``LoweredParameterSet.log_prior`` —
and hands back a pure jax function of the unconstrained vector. It is the
lowering of a *problem*, in exactly the sense ``lowering.md`` §0 defines: a
one-way translation of declarations, done once when a problem is realised on a
backend, never per evaluation.

**The two paths must agree, and the conformance suite is what says so.**
Nothing here re-derives the mathematics: the Gaussian log-density is
``GaussianFamily.log_prob``'s closed form transcribed, the noise quadrature is
the noise model's own, and the change of variables is
``lnprior_unconstrained``'s. A cross-check against ``problem.log_prob`` at
concrete points is the acceptance test, and it is in ``tests/backends``.

Trace purity
------------
No Python exception is raised on the traced path (``likelihoods.md`` §17 Q1's
ruling, and the reason it exists). Every refusal this module makes — a model
that is not this backend's, a family neither ``ampere.core`` nor this backend
implements, a solver from another rung of the ladder — is made **at lowering
time**, once, with a :class:`~ampere.core.exceptions.LoweringError` naming
what is unsupported. Inside the traced function a failure is ``-inf``,
computed with ``jnp.where``, which is §4.5's answer arrived at by the only
means a trace allows.

What this lowers
----------------
Everything ``ampere.core`` implements, bar the pieces named below:

* this backend's models and instrument steps, on any number of jointly fitted
  datasets, with ties, plates and hierarchical priors;
* every **implemented** likelihood family — ``gaussian``, ``student_t``,
  ``cauchy``, ``complex_gaussian``, ``poisson`` — through
  :mod:`ampere.backends.jax.families`, which transcribes ``ampere.core``'s
  closed forms and refuses by name anything it does not hold;
* **censoring**: the Tobit decomposition, for every family that declares
  ``SUPPORTS_CENSORING``, with the CDFs written out where ``jax.scipy.stats``
  has none;
* **latent GPs**: a dataset declaring one hands its whitened block to the
  family that consumes it (``CONSUMES_LATENT_GP``), which today is Poisson —
  the combination ``DEVELOPMENT_PLAN.md`` §4.4 singles out, and the one no
  gradient-free engine can run at all;
* both GP solves, **dense and quasiseparable**, differentiable in the kernel
  hyperparameters. ``QuasisepGP`` is what makes the flexible likelihood
  tractable at 10⁵ to 10⁶ samples, and it lowers here exactly as ``DenseGP``
  does because the solver interface is the strategy the noise model holds.

Refused, by name, at lowering time: a model or step from another backend, a
GP solver from another backend, a family ``ampere.core`` itself does not
implement (``rice``, ``von_mises``), and the analytic-GP combinations
``ampere.core`` declares but has not implemented (``complex_gaussian`` with a
GP, Phase 4's).

Batching
--------
:meth:`LoweredProblem.log_prob_unconstrained_batched` is ``jax.vmap`` of the
same function, and it is offered **only when the problem declares
``batchable``** — which it does when every part does. That is not caution for
its own sake: ``celerite2.jax``'s primitives register no batching rule, so a
``vmap`` over a quasiseparable density fails *inside* the trace with a message
about a primitive the user never named. Asking the declaration first turns
that into a refusal by name, which is this architecture's rule about
capabilities everywhere else too.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np

from ampere.core import (
    FittingProblem,
    GaussianProcessNoise,
    IndependentNoise,
)
from ampere.core.dataset import (
    INSTRUMENT_COMPONENT,
    LATENT_COMPONENT,
    LIKELIHOOD_COMPONENT,
)
from ampere.core.exceptions import LoweringError

from ._config import BACKEND, require_x64
from .families import lower_family
from .parameters import LoweredParameterSet

__all__ = ["LoweredProblem", "lower_problem"]

_LOG_2PI = math.log(2.0 * math.pi)


def _refuse(what: str, detail: str) -> LoweringError:
    return LoweringError(what, backend=BACKEND, detail=detail)


def _is_native_solver(solver: Any) -> bool:
    """Whether *solver* is one of this backend's, by declaration.

    Two conditions, and both matter. ``BACKEND`` is W2.12's identity flag, so
    a solver from another rung of the ladder is caught even if it happens to
    offer the right method name; ``log_marginal_likelihood_jax`` is the native
    surface this module actually calls, so a jax solver that had not supplied
    it would be caught before the trace rather than inside it. Checked by
    declaration rather than by ``isinstance`` so that a user's own jax solver
    — the strategy interface exists to be extended — composes here without
    subclassing one of ours.
    """
    return getattr(solver, "BACKEND", None) == BACKEND and callable(
        getattr(solver, "log_marginal_likelihood_jax", None)
    )


class _LoweredDataset:
    """One dataset's forward chain and log-likelihood, as jax."""

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
                f"model {self.model_label!r} has no native jax surface (a `flux` method), so it "
                f"cannot be composed into a differentiable log-density. Build the problem from "
                f"ampere.backends.jax's models, or run it on a gradient-free engine.",
            )
        for step in self.steps:
            if not hasattr(step, "apply_flux"):
                raise _refuse(
                    type(step).__name__,
                    f"instrument step {step.label!r} of dataset {label!r} has no native jax "
                    f"surface (an `apply_flux` method). Build the chain from "
                    f"ampere.backends.jax's steps.",
                )
        family = self.likelihood.family
        #: The family's own closed form, in jax (:mod:`ampere.backends.jax.families`).
        #: Refuses by name, here, for anything ``ampere.core`` does not implement
        #: or this backend has not transcribed.
        self.family_log_prob = lower_family(family, label)
        correlated = bool(getattr(self.noise, "CORRELATED", False))
        if correlated and not _is_native_solver(self.noise.solver):
            raise _refuse(
                type(self.noise.solver).__name__,
                f"dataset {label!r} uses the {type(self.noise.solver).__name__} solver, which is "
                f"not this backend's. A jax problem whose GP solve ran in numpy would not be "
                f"differentiable; pass ampere.backends.jax.DenseGP or "
                f"ampere.backends.jax.QuasisepGP.",
            )
        if not isinstance(self.noise, (IndependentNoise, GaussianProcessNoise)):
            raise _refuse(
                type(self.noise).__name__,
                f"dataset {label!r} uses a noise model this backend does not know how to "
                f"compose natively.",
            )
        self.correlated = correlated
        # Which of the two correlated stories this dataset tells. A family
        # whose GP marginalises in closed form (Gaussian) hands the whole
        # covariance to the solver; one that does not (Poisson) needs the
        # latent block instead. The *declaration* decides, not a name check,
        # so a family added to ampere.core lands on the right branch here.
        self.gp_marginal = correlated and family.ANALYTIC_WITH_GP
        if self.gp_marginal and not family.GP_ANALYTIC_IMPLEMENTED:
            raise _refuse(
                family.NAME,
                f"dataset {label!r} composes the {family.NAME!r} family with a correlated noise "
                f"model. That marginalisation is declared analytic but is unimplemented on the "
                f"reference path too (GP_ANALYTIC_IMPLEMENTED = False), so there is nothing for "
                f"this backend to agree with. ampere.core refuses the same composition.",
            )
        censoring = self.likelihood.censoring
        if censoring is not None and not family.SUPPORTS_CENSORING:
            raise _refuse(
                "censoring",
                f"dataset {label!r} declares censored samples under the {family.NAME!r} family, "
                f"which cannot consume a censoring declaration. ampere.core refuses the same "
                f"composition at construction.",
            )
        self.latent_name: str | None = None
        if dataset.latent is not None:
            if not family.CONSUMES_LATENT_GP:
                raise _refuse(
                    "latent",
                    f"dataset {label!r} declares a latent GP that the {family.NAME!r} family "
                    f"does not consume. ampere.core refuses this composition at construction "
                    f"(Likelihood._refuse_unconsumed_latent); reaching it here is a bug.",
                )
            self.latent_name = dataset.latent.parameter.name

        observed = dataset.observed
        # The effective mask, resolved once at composition: Dataset takes the
        # union of the observed and predicted masks there and forbids a
        # parameter-dependent one, so which samples are retained is a fact
        # about the problem rather than about theta and can be a constant here.
        # Public since W2.13 (fold-in 9) precisely because this line read it.
        excluded = dataset.effective_mask
        if excluded is None:
            excluded = ~np.asarray(observed.valid, dtype=bool).ravel()
        self.retain = ~np.asarray(excluded, dtype=bool).ravel()
        # Complex for a complex family, real otherwise, taken from the
        # container rather than assumed: `check_alignment` has already refused
        # a complex observation under a real family, so the container's own
        # dtype is the declaration.
        raw = np.asarray(observed.values)
        self.observed_values = jnp.asarray(
            raw[self.retain],
            dtype=jnp.complex128 if raw.dtype.kind == "c" else jnp.float64,
        )
        self.observed_coordinates = jnp.asarray(
            np.asarray(observed.axes[0].values, dtype=float)[self.retain], dtype=jnp.float64
        )
        self.uncertainty = (
            None
            if observed.uncertainty is None
            else jnp.asarray(
                np.asarray(observed.uncertainty, dtype=float)[self.retain], dtype=jnp.float64
            )
        )
        #: The Tobit codes for the **retained** samples, or ``None``.
        #: ``Likelihood._retained_limits`` computes exactly this on the numpy
        #: path, once per evaluation; a censoring declaration cannot depend on
        #: theta, so here it is a constant resolved at lowering time.
        self.limits = (
            None
            if censoring is None
            else jnp.asarray(np.asarray(censoring.kinds)[self.retain], dtype=jnp.int32)
        )

    # -- routing ------------------------------------------------------------

    def _dataset_values(self, routed: Mapping[str, Mapping[str, Any]]) -> Mapping[str, Any]:
        return self.dataset.route(routed[self.label])

    def _instrument_values(self, routed: Mapping[str, Mapping[str, Any]]) -> Mapping[str, Any]:
        chain = self._dataset_values(routed).get(INSTRUMENT_COMPONENT)
        if not chain:
            return {}
        return self.dataset.instrument.mapping.distribute(chain)

    # -- the forward chain --------------------------------------------------

    def predict(self, routed: Mapping[str, Mapping[str, Any]]) -> jax.Array:
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

    # -- the log-likelihood -------------------------------------------------

    def _sigma(self, predicted: jax.Array, values: Mapping[str, Any]) -> jax.Array | None:
        """``sigma`` for the retained samples, in jax.

        A prediction-aware noise model supplies its own traced surface
        (``sigma_jax``); the plain ones are ``scale``/``jitter`` applied to the
        observed uncertainties, transcribed here rather than called, because
        ``NoiseModel.sigma`` coerces with ``float()``.
        """
        native = getattr(self.noise, "sigma_jax", None)
        if native is not None:
            return native(self.dataset.observed, self.retain, values, predicted=predicted)
        own = {key: value for key, value in values.items() if key in self.noise.parameters}
        resolved = self.noise.context(own)
        sigma = self.uncertainty
        if sigma is None:
            if "jitter" not in resolved:
                return None
            return jnp.full(
                int(self.retain.sum()), jnp.asarray(resolved["jitter"], dtype=jnp.float64)
            )
        if "scale" in resolved:
            sigma = sigma * jnp.asarray(resolved["scale"], dtype=jnp.float64)
        if "jitter" in resolved:
            floor = jnp.asarray(resolved["jitter"], dtype=jnp.float64)
            sigma = jnp.sqrt(sigma**2 + floor**2)
        return sigma

    def _latent(self, routed: Mapping[str, Mapping[str, Any]]) -> jax.Array | None:
        """The dataset's latent block, or ``None`` when it declares none.

        The whitened values ``z`` arrive as one array-valued parameter under
        the ``latent`` component (``Dataset.route``), exactly as they do on
        the numpy path, and are handed to the family unchanged — which is what
        ``Dataset.log_likelihood_of`` does. The correlation would enter through
        ``GPSolver.latent_transform``; the reference path does not apply it and
        neither does this one, because the two must agree and the numpy path is
        the oracle. See :meth:`LoweredProblem.log_likelihood_terms`' note.
        """
        if self.latent_name is None:
            return None
        block = self._dataset_values(routed).get(LATENT_COMPONENT, {})
        return jnp.asarray(block[self.latent_name], dtype=jnp.float64).reshape(-1)

    def log_likelihood(self, routed: Mapping[str, Mapping[str, Any]]) -> jax.Array:
        """``log p(data | θ)`` for this dataset, as a traceable jax scalar."""
        predicted = self.predict(routed)
        values = dict(self._dataset_values(routed).get(LIKELIHOOD_COMPONENT, {}))
        # The noise model's own parameters arrive under the likelihood
        # component, flatly -- the same mapping ampere.core hands NoiseModel.
        sigma = self._sigma(predicted, values)
        if self.gp_marginal:
            # The whole covariance goes to the solver: this is the flexible
            # likelihood, and the branch GaussianFamily.log_prob takes when
            # `noise.correlated`.
            residual = self.observed_values - predicted
            variance = jnp.zeros_like(residual) if sigma is None else sigma**2
            value = self.noise.solver.log_marginal_likelihood_jax(
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
            value = self.family_log_prob(
                predicted,
                self.observed_values,
                sigma,
                self.likelihood.family,
                values,
                self.limits,
                self._latent(routed),
            )
        return jnp.where(jnp.isfinite(value), value, -jnp.inf)


class LoweredProblem:
    """A :class:`~ampere.core.FittingProblem`, lowered onto jax.

    Built by :func:`lower_problem`. Exposes the same three quantities §4.5
    does — a log prior, a log likelihood and their sum — as **pure jax
    functions**, plus the unconstrained form a NUTS kernel needs.

    Attributes
    ----------
    problem
        The problem this lowers. Kept so a driver can read the seed, the
        parameter names and the capability flags off it rather than
        re-deriving them.
    parameters
        The :class:`~ampere.backends.jax.parameters.LoweredParameterSet` — the
        numpyro/biject_to view of the merged declaration.
    """

    def __init__(self, problem: FittingProblem) -> None:
        require_x64("a lowered jax problem")
        if problem.backend != BACKEND:
            raise _refuse(
                problem.backend,
                f"this problem reports backend {problem.backend!r}, not {BACKEND!r}. The backend "
                f"is derived from what the models and transformations declare (W2.12), so a "
                f"problem built from another backend's pieces cannot be lowered here — build it "
                f"from ampere.backends.jax's models and steps.",
            )
        self.problem = problem
        self.parameters = LoweredParameterSet(problem.parameters, strict=problem.strict)
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

        See :meth:`~ampere.backends.jax.parameters.LoweredParameterSet.lowering_provenance`
        — ``lowering.md`` §12.8's "every registered row is stamped
        user-registered in provenance", ready for ``provenance_attrs``'s
        ``extra=``.
        """
        return self.parameters.lowering_provenance()

    def _route(self, theta: jax.Array) -> dict[str, dict[str, Any]]:
        """The constrained free vector, routed to each component's local names.

        ``ParameterMapping.distribute`` is pure routing metadata — it looks a
        value up by name and puts it where the declaration says it goes — so it
        carries tracers through untouched, and reusing it is what keeps this
        module from acquiring its own (divergent) idea of how a tie is bound.
        """
        return self._mapping.distribute(self.parameters._resolved(theta))

    def log_prior(self, theta: Any) -> jax.Array:
        """``log p(θ)`` in the constrained parameterisation."""
        return self.parameters.log_prior(theta)

    def log_likelihood(self, theta: Any) -> jax.Array:
        """``log p(data | θ)``, summed over datasets."""
        total = jnp.asarray(0.0, dtype=jnp.float64)
        for value in self._likelihood_terms(theta).values():
            total = total + value
        return jnp.where(jnp.isfinite(total), total, -jnp.inf)

    def _likelihood_terms(self, theta: Any) -> dict[str, jax.Array]:
        """Each dataset's own ``log p(data | θ)``, in the **constrained** space."""
        vector = jnp.asarray(theta, dtype=jnp.float64).reshape(-1)
        routed = self._route(vector)
        return {dataset.label: dataset.log_likelihood(routed) for dataset in self._datasets}

    def log_likelihood_terms(self, unconstrained: Any) -> dict[str, jax.Array]:
        """``inference.md`` §10a's optional member: the per-dataset decomposition.

        Keyed by dataset label, in the **unconstrained** parameterisation —
        the same argument :meth:`log_prob_unconstrained` takes, because that
        is the vector a NUTS driver has in hand for every draw. The terms sum
        to :meth:`log_likelihood` at the corresponding constrained point, and
        that identity is what ``ampere.results.emit``'s per-dataset
        ``log_likelihood`` group records.

        Supplying it is optional by ruling (sub-decision 1). Supplying it here
        is worth it because it costs nothing: the joint likelihood is already
        a sum over exactly these terms, so returning them instead of the sum
        is a change of return type rather than a second evaluation — and it
        spares a driver recomputing every stored draw on the numpy path.
        """
        y = jnp.asarray(unconstrained, dtype=jnp.float64).reshape(-1)
        return self._likelihood_terms(self.parameters.constrain_jax(y))

    def log_prob(self, theta: Any) -> jax.Array:
        """``log p(θ) + log p(data | θ)`` — the number a sampler maximises."""
        prior = self.log_prior(theta)
        total = prior + self.log_likelihood(theta)
        return jnp.where(jnp.isfinite(total), total, -jnp.inf)

    def log_prob_unconstrained(self, unconstrained: Any) -> jax.Array:
        """:meth:`log_prob` in unconstrained space, change of variables included.

        The density a NUTS kernel works with, and the jax counterpart of
        ``FittingProblem.log_prob_unconstrained``. The Jacobian term comes from
        :meth:`~ampere.backends.jax.parameters.LoweredParameterSet.log_prior_unconstrained`,
        which is the one place in this backend that calls
        ``log_abs_det_jacobian``, so the argument order of ``lowering.md``
        §2(b) is got right once.
        """
        y = jnp.asarray(unconstrained, dtype=jnp.float64).reshape(-1)
        prior = self.parameters.log_prior_unconstrained(y)
        constrained = self.parameters.constrain_jax(y)
        total = prior + self.log_likelihood(constrained)
        return jnp.where(jnp.isfinite(total), total, -jnp.inf)

    @property
    def batchable(self) -> bool:
        """Whether :meth:`log_prob_unconstrained_batched` is offered.

        Read off the *problem*, which aggregates it from what the models,
        instrument steps, noise models and GP solvers declare
        (``ampere.core.declared_capabilities``, conjunctively). One part that
        cannot be batched withdraws the claim for the whole density, which is
        the right answer: a ``vmap`` is over the composed function.
        """
        return bool(self.problem.batchable)

    def log_prob_unconstrained_batched(self, unconstrained: Any) -> jax.Array:
        """:meth:`log_prob_unconstrained` over a ``(batch, n_dim)`` stack.

        ``architecture.md`` §1's ``BATCHABLE`` rung, cashed: one traced
        function evaluated for many parameter vectors at once, which is what
        an ensemble sampler, a population-based optimiser or a vectorised
        importance-sampling step wants, and what a GPU would want even for a
        single chain.

        It is ``jax.vmap`` of the same function and nothing else — there is no
        second implementation to drift, which is the whole reason a jax
        backend can offer batching cheaply where the reference path cannot.

        Refused, by name, when the problem declares ``batchable = False``.
        That refusal is load-bearing rather than defensive: ``celerite2.jax``
        registers no batching rule for its primitives, so a ``vmap`` over a
        quasiseparable density raises ``NotImplementedError: Batching rule for
        'celerite2_factor' not implemented`` from inside the trace — a message
        about a primitive the user never named, at a point that says nothing
        about which part of their problem is at fault. Asking
        :attr:`batchable` first turns it into a sentence naming the solver.
        """
        if not self.batchable:
            raise _refuse(
                "batching",
                f"this problem declares batchable = False, so its density cannot be vmapped. "
                f"The flag is aggregated from what every part declares "
                f"({self.problem.capabilities}), conjunctively, so one part is enough to "
                f"withdraw it — ampere.backends.jax.QuasisepGP is the one that does, because "
                f"celerite2's primitives register no jax batching rule. Use DenseGP if the "
                f"problem is small enough, or evaluate the density in a loop.",
            )
        stacked = jnp.asarray(unconstrained, dtype=jnp.float64)
        if stacked.ndim != 2:
            raise _refuse(
                "batching",
                f"a batched density takes a (batch, n_dim) stack of unconstrained vectors, got "
                f"shape {tuple(stacked.shape)}.",
            )
        return jax.vmap(self.log_prob_unconstrained)(stacked)

    def potential(self) -> Callable[[jax.Array], jax.Array]:
        """The **negative** unconstrained log-density, which is what numpyro wants.

        numpyro's kernels minimise a potential energy rather than maximise a
        log-density, and getting that sign wrong produces a sampler that
        explores the prior's tails with perfect efficiency. Spelled once, here.
        """

        def potential_energy(unconstrained: jax.Array) -> jax.Array:
            return -self.log_prob_unconstrained(unconstrained)

        return potential_energy

    def __repr__(self) -> str:
        return (
            f"<LoweredProblem {self.free_size} free dimension(s), "
            f"{len(self._datasets)} dataset(s), backend={BACKEND!r}>"
        )


def lower_problem(problem: FittingProblem) -> LoweredProblem:
    """Lower *problem* onto jax: one differentiable log-density, built once.

    ``lowering.md`` §0: "It happens once, when a ``FittingProblem`` is realised
    on a backend — never per evaluation."

    Raises
    ------
    LoweringError
        Naming, by class, whatever part of the problem this backend cannot
        compose natively — a model without a jax surface, a censored dataset, a
        non-Gaussian family, another backend's GP solver. Never a silent
        fallback to the reference path: a nominally jax-backed NUTS run that
        quietly lost its gradients would be a run that cannot work at all.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.jax import (
    ...     IndependentNoise, PowerLaw, configure_x64, lower_problem
    ... )
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
    >>> configure_x64()
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
    >>> bool(np.isfinite(lowered.log_prob_unconstrained(y)))
    True

    and it is differentiable, which is the whole point:

    >>> import jax
    >>> gradient = jax.grad(lowered.log_prob_unconstrained)(y)
    >>> gradient.shape == (lowered.free_size,)
    True

    and it supplies §10a's optional per-dataset decomposition, which sums to
    the joint log-likelihood:

    >>> terms = lowered.log_likelihood_terms(y)
    >>> sorted(terms)
    ['default']
    >>> float(sum(terms.values())) == float(
    ...     lowered.log_likelihood(lowered.parameters.constrain_jax(y))
    ... )
    True
    """
    return LoweredProblem(problem)
