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
* **latent GPs**: a dataset declaring one pushes its whitened block ``z``
  through the solver's native whitening transform and hands ``f = L(θ) z`` to
  the family that consumes it (``CONSUMES_LATENT_GP``), which today is Poisson
  — the combination ``DEVELOPMENT_PLAN.md`` §4.4 singles out, and the one no
  gradient-free engine can run at all. The transform is applied here since
  W2.14, when ``ampere.core`` began applying it too;
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

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np

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
    chunk_bounds,
    foreign_parts_refusal,
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


def _refuse(what: str, detail: str) -> LoweringError:
    return LoweringError(what, backend=BACKEND, detail=detail)


#: The two spellings of a native model's value-and-coordinates surface, in the
#: order they are looked for (*W4.3*). ``flux``/``grid`` is the original pair;
#: ``native_flux``/``native_grid`` exists because a model may not have a method
#: whose name one of its own parameters already uses --
#: ``Parameterised._check_free_name`` refuses a parameter that shadows a class
#: attribute, and ``flux`` is precisely what an interferometric source model
#: calls its total flux density. Either pair composes; a model that offers both
#: is taken at the first, and one that offers half of either is refused with the
#: missing half named, because the two go together.
_FLUX_NAMES: tuple[str, ...] = ("flux", "native_flux")
_GRID_NAMES: tuple[str, ...] = ("grid", "native_grid")


def _native_surface(model: Any, names: tuple[str, ...]) -> Any | None:
    """The first callable *model* offers under *names*, or ``None``."""
    for name in names:
        found = getattr(model, name, None)
        if callable(found):
            return found
    return None


#: Neutral family name -> the ``ampere.core`` class whose ``sample`` this
#: backend has a native twin for (*W3.14*). The **one** place this module
#: names them, so the ceiling ``inference.md`` §13 states — a backend samples
#: exactly what the core samples — is one table rather than a condition
#: repeated per family: a core family that gains a ``sample`` joins by adding
#: a row here and a branch to :meth:`_LoweredDataset.sample_retained`, and one
#: that does not have a twin here falls back to the numpy path unchanged.
_TWINNED_FAMILIES: dict[str, type[LikelihoodFamily]] = {
    "gaussian": GaussianFamily,
    "poisson": PoissonFamily,
    "student_t": StudentTFamily,
    "complex_gaussian": ComplexGaussianFamily,
}


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

        self.model_flux = _native_surface(self.model, _FLUX_NAMES)
        self.model_grid = _native_surface(self.model, _GRID_NAMES)
        missing = ", ".join(
            f"`{names[0]}` (or `{names[1]}`)"
            for names, found in ((_FLUX_NAMES, self.model_flux), (_GRID_NAMES, self.model_grid))
            if found is None
        )
        if missing:
            raise _refuse(
                type(self.model).__name__,
                f"model {self.model_label!r} has no native jax surface ({missing} missing), so "
                f"it cannot be composed into a differentiable log-density. The two go together: "
                f"the first supplies the values and the second the coordinates the instrument "
                f"chain transforms them on, and `predict` calls both. Build the problem from "
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
        self.family_log_prob = lower_family(
            family, label, censored=self.likelihood.censoring is not None
        )
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
            if not callable(getattr(self.noise.solver, "latent_transform_jax", None)):
                raise _refuse(
                    type(self.noise.solver).__name__,
                    f"dataset {label!r} declares a latent GP, but the "
                    f"{type(self.noise.solver).__name__} solver has no `latent_transform_jax` — "
                    f"the traceable, differentiable form of `f = L(theta) z`. Since W2.14 the "
                    f"contract path applies that transform on every latent evaluation, so a "
                    f"realisation without it would disagree with its own oracle and give the "
                    f"kernel hyperparameters no gradient. A user-written solver joins the "
                    f"latent path by supplying it.",
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
        # Every axis, stacked into the `(n, d)` block ``Likelihood._coordinates``
        # builds -- **not** ``axes[0]`` (*W4.2*). A one-axis container is
        # unchanged by this (`(n, 1)` and `(n,)` are the same point set to
        # ``Kernel.matrix``), and a VisibilitySet has three axes, so taking the
        # first would have handed the kernel the `u` column and called it the
        # coordinates. The same stack is what makes ``axes=("u", "v")``
        # meaningful here, since a selector resolves to *columns of this block*.
        self.observed_coordinates = jnp.asarray(
            np.column_stack([np.asarray(axis.values, dtype=float) for axis in observed.axes])[
                self.retain
            ],
            dtype=jnp.float64,
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
            else jnp.asarray(
                np.asarray(observed.uncertainty, dtype=float)[self.retain], dtype=jnp.float64
            )
        )
        if self.gp_marginal and family.REQUIRES_UNCERTAINTY:
            self._check_gp_uncertainty(observed, label)
        #: The Tobit codes for the **retained** samples, or ``None``.
        #: ``Likelihood._retained_limits`` computes exactly this on the numpy
        #: path, once per evaluation; a censoring declaration cannot depend on
        #: theta, so here it is a constant resolved at lowering time.
        self.limits = (
            None
            if censoring is None
            else jnp.asarray(np.asarray(censoring.kinds)[self.retain], dtype=jnp.int32)
        )

    def _check_gp_uncertainty(self, observed: Any, label: str) -> None:
        """Refuse, at construction, a GP-marginal dataset the contract path cannot normalise.

        ``GaussianProcessNoise.sigma`` (``ampere.core.likelihood``) refuses a
        dataset with no observed uncertainty at all, or with a retained
        uncertainty that is zero or negative -- "an infinitely precise
        measurement, which no likelihood can normalise" -- for any family
        that ``REQUIRES_UNCERTAINTY``. On the contract path that surfaces
        only when the density is actually evaluated; ``ampere.core.realise``'s
        one-point agreement check happens to catch it today, because the
        contract path itself raises while computing the reference
        log-probability, but that is incidental to which one point gets
        checked, not a refusal this backend makes by name. A NUTS run
        evaluates every other sampled point too, so the refusal belongs here,
        at construction -- mirroring the ``REQUIRES_UNCERTAINTY``/sigma-is-None
        check :meth:`log_likelihood` already makes for the non-GP branch
        (this module's docstring, "No exception control flow on the hot
        path").
        """
        if observed.uncertainty is None:
            if "jitter" not in self.noise.parameters:
                raise _refuse(
                    "uncertainty",
                    f"dataset {label!r} has no observed uncertainties and its noise model "
                    f"declares no jitter, so sigma is undefined.",
                )
            return
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
        return self.predict_full(routed)[self.retain]

    def predict_full(self, routed: Mapping[str, Mapping[str, Any]]) -> jax.Array:
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

    def _latent(
        self, routed: Mapping[str, Mapping[str, Any]], values: Mapping[str, Any]
    ) -> jax.Array | None:
        """The dataset's latent function ``f``, or ``None`` when it declares none.

        The whitened values ``z`` arrive as one array-valued parameter under
        the ``latent`` component (``Dataset.route``); the correlation is then
        imposed here by the solver's own **native** whitening transform, so
        what the family receives is ``f = L(θ) z`` on the retained
        coordinates.

        Until W2.14 this method handed the block over unchanged, and said so
        at length, because the numpy path did the same: nothing on the scoring
        path called ``GPSolver.latent_transform``, and a realisation that
        "fixed" it from a backend would have disagreed with the oracle it is
        checked against (``inference.md`` §10a sub-decision 4). The oracle is
        fixed now — ``GaussianProcessNoise.noise_params`` applies the
        transform — so this applies it too, and the conformance row holds the
        two together. ``latent_transform_jax`` rather than ``latent_transform``
        for the usual reason: the contract surface returns numpy and would cut
        the gradient in exactly the hyperparameters this path exists to fit.
        """
        if self.latent_name is None:
            return None
        block = self._dataset_values(routed).get(LATENT_COMPONENT, {})
        whitened = jnp.asarray(block[self.latent_name], dtype=jnp.float64).reshape(-1)
        kernel = self._kernel()
        return self.noise.solver.latent_transform_jax(
            kernel,
            self.observed_coordinates,
            whitened,
            kernel.resolve(values),
        )

    # -- the generative half (W3.1 slice 2) ---------------------------------

    def sampling_refusal(self) -> LoweringError | None:
        """Why this dataset cannot be sampled natively, or ``None``.

        Peter's ruling of 2026-09-08 is that every backend samples natively;
        what a backend may sample is fixed by what ``ampere.core`` samples,
        because **the numpy path is the oracle**. A native twin for a family
        the core will not sample would be a backend inventing an observation
        process the contract declines to guess — precisely the "silently train
        an SBI posterior on the wrong forward model" the refusal exists to
        prevent — so this backend samples exactly what the core samples, and
        refuses the rest **by handing the draw back to the numpy path**, where
        the core's own refusal text is what the caller meets.

        **Amended W3.14**, when the core acquired three more ``sample``
        implementations. The set is now ``gaussian``, ``poisson``,
        ``student_t`` and ``complex_gaussian`` — :data:`_TWINNED_FAMILIES`,
        which is the one place this backend names them — and the test is still
        the same test: the family's ``sample`` must be the *core's own* method
        for that family. That keeps two properties at once. A family the core
        refuses (``cauchy``, ``rice``, ``von_mises``, a user family with no
        ``sample``) is refused here, so the caller meets the core's refusal
        text. A family whose ``sample`` a *user* has overridden is refused too,
        and deliberately: it is a numpy function, it is the observation process
        its author wrote, and running something else instead would be worse
        than running it slowly (§13's "a user with an exotic observation
        process supplies it by subclassing").

        Two conditions survive the widening:

        * a censoring declaration that survives the mask blocks a draw on every
          backend, because applying the censoring operator to a draw is not
          implemented anywhere;
        * a latent declaration means the family reads ``noise.latent``, which
          only ``poisson`` does. The other three would silently ignore it.

        The complex condition is gone: complex data are exactly what the
        ``complex_gaussian`` twin draws, and a complex container under any
        other family is refused by ``Likelihood.check_alignment`` long before
        a draw is asked for.
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
        if self.latent_name is not None and family.NAME != "poisson":
            return _refuse(
                "latent",
                f"dataset {self.label!r} declares a latent GP, whose family reads "
                f"noise.latent; the {family.NAME!r} twin does not consume one.",
            )
        return None

    def sample_retained(
        self,
        routed: Mapping[str, Mapping[str, Any]],
        predicted: jax.Array,
        key: jax.Array,
    ) -> jax.Array:
        """One draw of the retained observed values, natively.

        The core family's own ``sample`` transcribed into jax, dispatched on
        the neutral family name exactly as :meth:`log_likelihood` dispatches
        the density. Every branch is traceable, so ``jax.vmap`` maps it over a
        chunk and the key is per draw; ``jax.random`` takes a *traced* rate and
        a traced ``nu``, which is why this backend needs none of the
        two-stage machinery torch's twins do.

        * ``gaussian`` — ``x = mu + sigma z`` uncorrelated, and
          ``x = mu + L z1 + sigma z2`` under a GP, with ``L`` from the solver's
          **own** ``latent_transform_jax``, the same whitening the latent path
          uses, and the solver's numerical jitter folded into the diagonal
          because it is part of the covariance the marginal likelihood scores.
          Omitting it would draw from a narrower distribution than the density
          evaluates, and the error is not small at the jitter values the
          library's own error message tells a user to raise.
        * ``poisson`` (*W3.14*) — ``counts ~ Poisson(rate)``, the rate being
          the prediction, or ``prediction * exp(f)`` under a GP with ``f`` the
          dataset's own latent block through :meth:`_latent` — the same ``f``
          the density scores at, so θ and the drawn data describe one model.
          Returned as float64, because that is what the containers hold and
          what ``ampere.core.PoissonFamily.sample`` returns.
        * ``student_t`` (*W3.14*) — ``x = mu + sigma t_nu``, with ``sigma``
          used as the *scale*, exactly as the density standardises by it.
        * ``complex_gaussian`` (*W3.14*) — independent ``Normal(0, sigma**2)``
          on each component, ``sigma`` being the per-component standard
          deviation the density's ``-|r|**2 / (2 sigma**2)`` implies.
        """
        values = dict(self._dataset_values(routed).get(LIKELIHOOD_COMPONENT, {}))
        name = self.likelihood.family.NAME
        if name == "poisson":
            return self._sample_poisson(routed, predicted, values, key)
        if name == "student_t":
            return self._sample_student_t(predicted, values, key)
        if name == "complex_gaussian":
            return self._sample_complex_gaussian(predicted, values, key)
        return self._sample_gaussian(predicted, values, key)

    def _sample_gaussian(
        self, predicted: jax.Array, values: Mapping[str, Any], key: jax.Array
    ) -> jax.Array:
        sigma = self._sigma(predicted, values)
        realisation = predicted
        size = int(self.retain.sum())
        if self.correlated:
            gp_key, noise_key = jax.random.split(key)
            realisation = realisation + self._gp_realisation(values, gp_key, size)
            stabiliser = float(getattr(self.noise.solver, "jitter", 0.0) or 0.0)
            if stabiliser:
                floor = jnp.full((size,), stabiliser, dtype=jnp.float64)
                sigma = floor if sigma is None else jnp.sqrt(sigma**2 + floor**2)
        else:
            noise_key = key
        if sigma is None:
            return realisation
        return realisation + sigma * jax.random.normal(noise_key, (size,), dtype=jnp.float64)

    def _sample_poisson(
        self,
        routed: Mapping[str, Mapping[str, Any]],
        predicted: jax.Array,
        values: Mapping[str, Any],
        key: jax.Array,
    ) -> jax.Array:
        size = int(self.retain.sum())
        rate = predicted
        if self.correlated:
            latent = self._latent(routed, values)
            if latent is None:
                # Unreachable through ``Dataset``, which declares a latent block
                # for every LATENT combination, but the core's own sample() has
                # the same fallback and the two must not diverge.
                gp_key, key = jax.random.split(key)
                latent = self._gp_realisation(values, gp_key, size)
            rate = rate * jnp.exp(latent)
        return jax.random.poisson(key, rate, (size,)).astype(jnp.float64)

    def _sample_student_t(
        self, predicted: jax.Array, values: Mapping[str, Any], key: jax.Array
    ) -> jax.Array:
        sigma = self._sigma(predicted, values)
        assert sigma is not None  # REQUIRES_UNCERTAINTY, checked at composition
        family = self.likelihood.family
        own = {key_: value for key_, value in values.items() if key_ in family.parameters}
        nu = jnp.asarray(family.context(own)["nu"], dtype=jnp.float64)
        size = int(self.retain.sum())
        return predicted + sigma * jax.random.t(key, nu, (size,), dtype=jnp.float64)

    def _sample_complex_gaussian(
        self, predicted: jax.Array, values: Mapping[str, Any], key: jax.Array
    ) -> jax.Array:
        sigma = self._sigma(predicted, values)
        assert sigma is not None  # REQUIRES_UNCERTAINTY, checked at composition
        size = int(self.retain.sum())
        if not self.correlated:
            real_key, imaginary_key = jax.random.split(key)
            real = jax.random.normal(real_key, (size,), dtype=jnp.float64)
            imaginary = jax.random.normal(imaginary_key, (size,), dtype=jnp.float64)
            return predicted + sigma * (real + 1j * imaginary)
        # W4.2: two real GP realisations sharing L, not one complex one. Sharing
        # L is circularity's equal-component half; sharing nothing else is its
        # zero-pseudo-covariance half. Drawing one realisation for both
        # components would give a draw perfectly correlated between the parts,
        # whose modulus statistics the density does not score.
        gp_key, noise_key = jax.random.split(key)
        correlated = self._gp_realisation(values, gp_key, size, components=2)
        stabiliser = float(getattr(self.noise.solver, "jitter", 0.0) or 0.0)
        scale = sigma if not stabiliser else jnp.sqrt(sigma**2 + stabiliser**2)
        independent = jax.random.normal(noise_key, (size, 2), dtype=jnp.float64)
        components = correlated + scale[:, None] * independent
        return predicted + components[:, 0] + 1j * components[:, 1]

    def _gp_realisation(
        self, values: Mapping[str, Any], key: jax.Array, size: int, components: int = 1
    ) -> jax.Array:
        """``L(θ) z`` for a fresh whitened draw, from the solver's own transform.

        *components* is ``k``: one for a real process, two for the circular
        complex GP, whose real and imaginary parts share the factor ``L`` and
        nothing else (*W4.2*).
        """
        shape = (size,) if components == 1 else (size, components)
        whitened = jax.random.normal(key, shape, dtype=jnp.float64)
        kernel = self._kernel()
        return self.noise.solver.latent_transform_jax(
            kernel,
            self.observed_coordinates,
            whitened,
            kernel.resolve(values),
        )

    def _kernel(self) -> Any:
        """This dataset's kernel, bound to its container's axis order (*W4.2*)."""
        assert self.bound_kernel is not None  # only a correlated dataset asks
        return self.bound_kernel

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
            # W4.2: a complex residual becomes the circular GP's two real
            # columns, and the per-component sigma is the real diagonal of the
            # one covariance both columns are scored against.
            if jnp.iscomplexobj(residual):
                residual = jnp.stack([residual.real, residual.imag], axis=-1)
            variance = (
                jnp.zeros(jnp.shape(residual)[0], dtype=jnp.float64) if sigma is None else sigma**2
            )
            kernel = self._kernel()
            value = self.noise.solver.log_marginal_likelihood_jax(
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
            value = self.family_log_prob(
                predicted,
                self.observed_values,
                sigma,
                self.likelihood.family,
                values,
                self.limits,
                self._latent(routed, values),
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
        # W3.8: a problem that reports this backend but composes a piece from
        # another one. That composition is refused at construction unless it
        # asked for allow_foreign_parts=True, and the flag buys a gradient-free
        # run and nothing more -- a lowered problem *is* the differentiable
        # form, so it refuses by name here regardless.
        foreign = foreign_parts_refusal(problem, what="a differentiable jax problem")
        if foreign is not None:
            raise foreign
        self.problem = problem
        self.parameters = LoweredParameterSet(problem.parameters, strict=problem.strict)
        self._mapping = problem.mapping
        self._datasets = tuple(_LoweredDataset(problem, label) for label in problem.datasets)
        #: ``jax.jit`` of :meth:`_terms_unconstrained`, built on first use. See
        #: :meth:`log_likelihood_terms` for why this one member is compiled
        #: here rather than left to whatever transformation the caller applies.
        self._terms_compiled: Callable[[Any], dict[str, jax.Array]] | None = None

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

    def _terms_unconstrained(self, unconstrained: Any) -> dict[str, jax.Array]:
        """:meth:`log_likelihood_terms`, uncompiled. The pure function jit wraps."""
        y = jnp.asarray(unconstrained, dtype=jnp.float64).reshape(-1)
        return self._likelihood_terms(self.parameters.constrain_jax(y))

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

        **Compiled, since W2.5 slice 3.** This member used to be the only part
        of the realised surface a caller might invoke many times *without* a
        surrounding transformation: NUTS reaches the density through
        :meth:`potential`, which numpyro jits, but a driver asking for the
        decomposition — or, since slice 3, ``ampere.inference``'s gradient-free
        fast path asking for it once per proposal — got jax's eager,
        op-at-a-time dispatch. Measured on a 400-point quasiseparable problem
        that is 28.8 ms a call against 0.14 ms compiled, a factor of two
        hundred, and it is the whole of what made the jax contract path
        (26.5 ms) look like a pessimisation for emcee, dynesty and zeus.

        :func:`jax.jit` is applied to the *uncompiled* pure function and cached
        on the instance, so the compilation happens once per lowered problem
        and shapes are fixed by ``free_size``. It composes as everything else
        in jax does: called inside another trace it is an ordinary nested jit,
        and :func:`jax.grad` through it still differentiates.
        """
        if self._terms_compiled is None:
            self._terms_compiled = jax.jit(self._terms_unconstrained)
        return self._terms_compiled(unconstrained)

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

    # -- batched simulation (W3.1 slice 2) ----------------------------------

    def _forward(self, theta: jax.Array) -> tuple[dict[str, dict[str, Any]], dict[str, Any]]:
        """One θ through the whole noise-free forward model, natively.

        Returns what :class:`~ampere.core.BatchedPrediction` holds, as a pytree
        of jax arrays so that ``jax.vmap`` maps it in one traced call: every
        model's every channel, and every dataset's instrument-transformed
        prediction on the *whole* observed grid.

        Every channel, not only the bound ones, because ``results.md`` §11
        writes one training-set group per ``<model>.<channel>`` and a fast path
        that returned fewer would silently write a smaller file than the loop.
        """
        routed = self._route(jnp.asarray(theta, dtype=jnp.float64).reshape(-1))
        channels = {
            label: {
                channel: _native_surface(model, _FLUX_NAMES)(channel, routed.get(label, {}))
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

        ``(batch, free_size)`` **constrained** free vectors in — the
        coordinates ``SimulationBatch.theta`` holds, not the unconstrained ones
        :meth:`log_prob_unconstrained` takes, because a simulation budget is a
        set of parameter values rather than a set of sampler positions.

        **Per chunk, always.** ``chunk_size`` bounds how many simulations are
        vectorised at once and the chunks are *looped*; ``vmap`` over a whole
        budget is exactly the single-device memory trap Peter's ruling of
        2026-09-08 names, and the reason ``simulate_many`` is built on an
        executor rather than on this. ``None`` means one chunk, which is what
        ``simulate_many`` passes because it has already chunked.

        **What it does not do.** It draws no noise (:meth:`sample_observations`
        does), it flags no failures (there is no control flow inside a trace to
        flag with, so a simulator that fails produces NaNs, classified by the
        caller exactly as the loop classifies them), and it refuses rather than
        falling back — ``simulate_many`` owns the fallback, and owning it in
        two places would let a fast path quietly become a slow one.

        Refused, by name, when the problem declares ``batchable = False``. The
        refusal is load-bearing for the same reason
        :meth:`log_prob_unconstrained_batched`'s is: ``celerite2.jax``
        registers no batching rule for its primitives.
        """
        if not self.batchable:
            raise _refuse(
                "batched simulation",
                f"this problem declares batchable = False, so its forward model cannot be "
                f"vmapped over a stack of parameter vectors. The flag is aggregated from what "
                f"every part declares ({self.problem.capabilities}), conjunctively — "
                f"ampere.backends.jax.QuasisepGP is the one that withdraws it, because "
                f"celerite2's primitives register no jax batching rule. simulate_many falls "
                f"back to the loop, which is the semantics anyway.",
            )
        stack = jnp.asarray(theta, dtype=jnp.float64)
        if stack.ndim != 2 or int(stack.shape[1]) != self.free_size:
            raise _refuse(
                "batched simulation",
                f"simulate_batched takes a (batch, {self.free_size}) stack of constrained free "
                f"vectors, got shape {tuple(stack.shape)}.",
            )
        vectorised = jax.vmap(self._forward)
        channels: dict[str, dict[str, list[np.ndarray]]] = {}
        predicted: dict[str, list[np.ndarray]] = {}
        for start, stop in chunk_bounds(int(stack.shape[0]), chunk_size):
            piece = stack[start:stop]
            produced, prediction = (
                vectorised(piece) if sharder is None else sharder.shard(vectorised, piece)
            )
            for label, holding in produced.items():
                for channel, values in holding.items():
                    channels.setdefault(label, {}).setdefault(channel, []).append(
                        np.asarray(values)
                    )
            for label, values in prediction.items():
                predicted.setdefault(label, []).append(np.asarray(values))
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
    ) -> dict[str, np.ndarray]:
        """Draw the retained observed values for a chunk, natively.

        Peter's ruling of 2026-09-08. The distribution is
        :meth:`ampere.core.GaussianFamily.sample`'s, transcribed in
        :meth:`_LoweredDataset.sample_retained`; the **stream** is
        ``jax.random``'s, which is why the numpy path stays the oracle and the
        two are compared *distributionally* rather than draw for draw.

        *seeds* is one integer per draw, taken from the per-draw child
        generator ``simulate_many`` spawns **by index** — so partition
        independence carries over unchanged: draw *i* gets the same key
        whichever chunk it ran in, and a budget split 1/7/whole gives the same
        observations.

        Refuses by name, before drawing anything, for any dataset this backend
        may not sample (:meth:`_LoweredDataset.sampling_refusal`); the caller
        then runs the numpy path, where ``ampere.core``'s own refusal text is
        what a user meets.
        """
        for dataset in self._datasets:
            refusal = dataset.sampling_refusal()
            if refusal is not None:
                raise refusal
        stack = jnp.asarray(theta, dtype=jnp.float64)
        keys = jnp.stack([jax.random.PRNGKey(int(seed) & 0xFFFFFFFF) for seed in seeds])
        drawn: dict[str, np.ndarray] = {}
        for dataset in self._datasets:
            # The prediction's own dtype decides, as it does at lowering time
            # (W2.4 slice 3): a complex_gaussian dataset predicts complex
            # visibilities, and forcing them to float64 here would not fail
            # loudly -- it would draw around the real projection of the model.
            raw = np.asarray(predicted[dataset.label])
            rows = jnp.asarray(raw, dtype=jnp.complex128 if raw.dtype.kind == "c" else jnp.float64)
            retained = rows[:, dataset.retain]

            if dataset.likelihood.family.NAME == "poisson":
                # Eager, before jax.vmap ever runs (W4.0 (2)): unlike the rest
                # of this method, the *base* rate here is already a concrete
                # array -- simulate_batched computed it earlier in the chunk
                # -- so nothing stops checking it the way log_prob does
                # (families.py's `_poisson` refuses rate <= 0, scored to
                # -inf) before drawing from it, exactly as torch's
                # `_poisson_variates` does on its own concrete `rate`. A
                # correlated (GP) draw multiplies this by a strictly positive
                # `exp(latent)` inside the vmap, so checking the pre-latent
                # rate here is sufficient either way, and it is what makes a
                # bad rate one *draw's* failure rather than silently invented
                # counts that only fail much later, at evaluation.
                concrete_rate = np.asarray(retained, dtype=float)
                if not np.all(np.isfinite(concrete_rate)) or np.any(concrete_rate < 0.0):
                    raise _refuse(
                        "poisson",
                        f"dataset {dataset.label!r}: the model predicted a negative or "
                        f"non-finite expected count for at least one draw in this chunk, "
                        f"which no Poisson can be drawn from. Constrain the prediction to the "
                        f"positive half-line (a Log bijection on the norm, or a "
                        f"positive-support prior); ampere.core's own sample() refuses the same "
                        f"condition by name.",
                    )

            def one(vector: jax.Array, row: jax.Array, key: jax.Array, of: Any = dataset) -> Any:
                return of.sample_retained(self._route(vector.reshape(-1)), row, key)

            drawn[dataset.label] = np.asarray(jax.vmap(one)(stack, retained, keys))
        return drawn

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
