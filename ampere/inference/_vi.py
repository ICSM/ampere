"""The variational driver: stochastic variational inference over a realisation.

Private module; the class is :class:`ampere.inference.VIEngine`. Named with a
leading underscore for the reason ``_emcee.py`` and ``_nuts.py`` record for
theirs.

Why VI is here at all
---------------------
``DEVELOPMENT_PLAN.md`` §1 lists "gradient-based inference (NUTS, VI) where
models are differentiable" as one of the three things the redesign is for, and
§5's Phase 2 gives each backend track "NUTS + VI" as its inference deliverable.
The two answer different questions and neither replaces the other: NUTS is
asymptotically exact and expensive, VI is approximate and fast, and the honest
use of VI is as a *first look* at a posterior — or as the only tractable route
when the parameter space is large enough that no MCMC finishes.

The approximation is real and this driver does not disguise it. A run's attrs
record the guide family (``ampere_vi_guide``), and the guide family **is** the approximation: a
mean-field ``AutoNormal`` assumes the posterior factorises over parameters,
which is exactly wrong for the correlated posteriors ampere's fits produce
(a temperature and a scale are anticorrelated in every SED fit), and it
underestimates marginal variances when it is wrong. ``AutoMultivariateNormal``
captures the correlations and costs O(d²) parameters. Neither captures a
non-Gaussian posterior at all. That is why ``ampere_engine`` names the driver
and ``vi_guide`` names the family: a reader of an archived run must be able to
tell what was assumed. Since W5.0 the same fact is also machine-readable, in
the engine-neutral ``ampere_approximation`` root attribute
(``"mean_field"``/``"multivariate"``) every non-exact engine writes
(``results.md`` §9), and the guide's own log-density at each draw is stored
per draw as ``sample_stats.proposal_log_density`` — the fitted guide's density
moved into the same constrained coordinates the stored ``log_prior``/
``log_likelihood`` are, so
``exp(log_prior + log_likelihood - proposal_log_density)`` is an importance
weight an importance-corrected posterior can be built from using nothing but
the stored groups.

How it reaches a backend
------------------------
Exactly as :class:`~ampere.inference.NUTSEngine` does, and for the same
reason: ``ampere.inference`` may import ``ampere.core`` and ``ampere.results``
and nothing else, so the gradient comes from the backend's **realisation**
(``inference.md`` §10a) obtained through :func:`ampere.core.realise`, and the
variational library is imported lazily inside :meth:`VIEngine.run`.

Two libraries, one driver
-------------------------
The *density* is backend-neutral by the time this driver has it, but the
optimiser is not: pyro's ``SVI`` wants a torch model and numpyro's wants a jax
one, and neither will consume the other's array. So :meth:`VIEngine.run`
dispatches on ``problem.backend`` between :meth:`VIEngine._fit_pyro` and
:meth:`VIEngine._fit_numpyro`, exactly as ``_nuts.py`` dispatches between the
two NUTS kernels, and imports the one it needs lazily.
:data:`VARIATIONAL_LIBRARIES` is that table, and
:func:`supported_backends` is the registry intersected with it.

The two routes emit the **same** provenance keys — ``vi_guide``,
``vi_guide_class``, ``vi_steps``, ``vi_draws``, ``vi_optimiser``,
``vi_learning_rate``, ``vi_library``, ``vi_final_elbo``, ``vi_elbo_trace`` —
plus their own library versions. That is not tidiness: a reader of an archived
run should not have to know which library fitted it in order to know what was
fitted, and a key that existed on one backend only would make every
cross-backend comparison of VI runs a special case.

Turning a density into a model SVI can guide
---------------------------------------------
This is the one piece of real work, and it is worth stating plainly because it
looks like a trick. ``SVI`` takes a *model* — a probabilistic program with
sample sites — while a realisation is a bare log density of one flat
unconstrained vector, with no site structure at all (``inference.md`` §10a
fixes it that way, so that ampere's own ``ParameterMapping`` stays the only
thing that routes names). MCMC has ``potential_fn=`` for exactly this; SVI has
no equivalent, in either library.

So the model declares one site over the whole unconstrained vector. On the
pyro route that site is a standard normal, whose support is all of ``R**d``,
and the model then adds a :func:`pyro.factor` of
``density(theta) - base.log_prob(theta)``. The log
joint is therefore ``log N(theta; 0, I) + density(theta) - log N(theta; 0, I)``
= ``density(theta)``, exactly, and the standard normal is a *carrier* for the
site rather than a prior: it cancels identically, term by term, at every point.
Any distribution with full real support would do; the standard normal is
chosen because it is cheap and because its scale matches the unconstrained
coordinates a well-behaved bijection produces, which keeps the subtraction from
being a difference of two large numbers.

numpyro can do better, and does: it ships an ``ImproperUniform``, whose
contribution to the log-density is *identically* zero, so the jax route
declares its site under that and adds ``numpyro.factor("ampere_density",
density(theta))`` with nothing to subtract. Same target, one fewer
cancellation. The two routes therefore differ in the carrier and agree in the
density, which is the right place for a library difference to live.

The autoguides then see a single unconstrained real site of the right shape,
which is what they are best at. ``AutoNormal`` gives a diagonal Gaussian in
that space, ``AutoMultivariateNormal`` a full-covariance one; both are
*unconstrained-space* approximations, so the constrained posterior they imply
is already correctly warped by ampere's own bijections rather than by a second
set of pyro's.

Four guide families, and what each one costs (W5.14)
-----------------------------------------------------
W5.14 adds two families to W2.5's pair, from ``inference_extensions_memo.md``
§6's tier 1. The four now span the choice a user actually has to make, and the
choice **is** the approximation:

* ``normal`` — a diagonal Gaussian. O(d) guide parameters, and wrong in
  exactly the way a correlated posterior is correlated.
* ``multivariate`` — a full-covariance Gaussian, O(d²) guide parameters,
  fitted by maximising the ELBO.
* ``laplace`` — a full-covariance Gaussian too, but fitted differently and
  costing differently. SVI optimises a ``Delta`` guide to the MAP point (so
  what the ELBO trace records along the way is the *log joint*: a ``Delta``
  guide has no entropy, and its ELBO is therefore ``log p(x, z)`` at the
  current iterate), and the covariance is then the inverse Hessian of the
  negative log joint *at that point*, computed once at the end. It is a
  **local** answer — the curvature at one point, not a fit to the mass --
  which makes it cheap where the full-covariance ELBO fit is dear, exact when
  the posterior really is Gaussian, and arbitrarily wrong when it is not. It
  also needs a *second* derivative of the realised density, which is a real
  demand on a backend's lowering rather than a free one.
* ``flow`` — an inverse-autoregressive flow (Kingma et al. 2016) over a
  standard-normal base: the only family here that can represent a skewed,
  heavy-tailed or otherwise non-Gaussian posterior, at the cost of a small
  neural network per transform. Two consequences matter enough to be stated
  where a reader will meet them. It is autoregressive over the coordinates, so
  it needs **at least two** of them and a one-dimensional problem is refused
  by name (:meth:`VIEngine.run`); and **its fit is not a location and a
  scale**, so the three fitted-parameter attributes this driver keeps for a
  Gaussian guide (:attr:`VIEngine.guide_loc`, :attr:`VIEngine.guide_scale`,
  :attr:`VIEngine.guide_scale_tril`) all stay ``None`` — which is why a run's
  attrs carry ``vi_guide_parameters`` saying which of them were kept, rather
  than leaving a reader to infer it from the guide's name.

All four begin at the same place, and the flow gets there differently
-----------------------------------------------------------------------
Every driver in this package starts from a seeded draw from the joint prior,
and the three Gaussian guides are put there by ``init_loc_fn``. A flow has no
location parameter for ``init_loc_fn`` to set — pyro's ``AutoIAFNormal``
documents that it ignores the argument and warns if one is passed, and
numpyro's uses it for nothing the transform reads — because a flow's starting
point *is its base distribution*: a standard normal at the origin, which the
transforms must then learn to carry to wherever the posterior is.

That is not a small matter in ampere's coordinates. The unconstrained
coordinate of a parameter declared with an unbounded prior is the parameter
itself, so a posterior at ``norm = 2.0 ± 0.1`` is twenty base standard
deviations from the origin *and* a tenth of its width; measured while writing
this driver, a flow left at the origin had reached ``norm ≈ 1.1`` after four
thousand steps, while the mean-field guide, started at the prior draw, was
converged in six hundred. A guide family that only works on a problem whose
priors happen to be standardised would be a trap rather than a feature.

So this driver moves the flow's **base** to the start point instead, which is
the library-sanctioned hook for exactly this: ``get_base_dist`` is a method
both ``AutoIAFNormal`` implementations define and both ``get_posterior``
implementations call, and overriding it to return ``Normal(start, 1)`` in
place of ``Normal(0, 1)`` changes where the flow starts without touching what
it can represent or how its density is computed — ``TransformedDistribution``
scores the base at the matching point either way, so ``proposal_log_density``
needs no correction. The two routes make the identical one-line override, for
the identical reason, which is why it is written twice rather than abstracted:
each lives beside the lazy import of the library it overrides.

What a run looks like
---------------------
One "chain" of independent draws from the fitted guide. That is not a
limitation dressed up: guide draws are i.i.d. by construction, so there is no
chain structure to split and R-hat has nothing to say — ArviZ's convergence
diagnostics answer a question about Markov chains, and these are not one. The
ELBO trace is what a VI run has instead, and it is recorded.
"""

from __future__ import annotations

import dataclasses
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import LoweringError
from ampere.core.realisation import (
    Realisation,
    log_likelihood_terms_of,
    realise,
    registered_realisations,
)

from .engine import (
    DEFAULT_CACHE_SIZE,
    Engine,
    _refuse_foreign_parts,
    unconstrained_jacobian_correction,
)
from .exceptions import EngineError

__all__ = [
    "GAUSSIAN_GUIDES",
    "GUIDE_FAMILIES",
    "VARIATIONAL_LIBRARIES",
    "VIEngine",
    "supported_backends",
]

#: Which variational library drives which backend's realisation. Strings, not
#: imports, for the reason ``_nuts.py``'s :data:`SAMPLER_LIBRARIES` gives: this
#: module imports no backend and no library at module level.
#:
#: Both, since W2.5 slice 2 added the numpyro route the torch track's comment
#: reserved. A row here without a route in :meth:`VIEngine.run` would make
#: ``supported_backends()`` claim a capability that does not exist, so the two
#: land together — the same discipline ``_nuts.py``'s
#: :data:`~ampere.inference._nuts.SAMPLER_LIBRARIES` keeps.
VARIATIONAL_LIBRARIES: dict[str, str] = {"torch": "pyro", "jax": "numpyro"}

#: The guide families a user may name, and what each assumes. The **names** are
#: ampere's, not pyro's, because the same two families exist in numpyro under
#: the same names and a run's ``vi_guide`` attribute has to mean one thing
#: across backends — the same discipline W2.12 applied to backend names. The
#: values happen to be the class names both libraries use, which is why one
#: table serves both routes; that is a convenience, not the contract. The
#: contract is the *key*.
GUIDE_FAMILIES: dict[str, str] = {
    "normal": "AutoNormal",
    "multivariate": "AutoMultivariateNormal",
    "laplace": "AutoLaplaceApproximation",
    "flow": "AutoIAFNormal",
}

#: The guides whose fit this driver can write down as plain arrays: a location
#: and one of a diagonal scale or a covariance Cholesky factor. Anything not
#: named here fits something that is not a Gaussian — today only ``"flow"``,
#: whose parameters are a neural network's weights — and leaves all three of
#: :attr:`VIEngine.guide_loc`, :attr:`VIEngine.guide_scale` and
#: :attr:`VIEngine.guide_scale_tril` at ``None``. **W5.14.**
GAUSSIAN_GUIDES: frozenset[str] = frozenset({"normal", "multivariate", "laplace"})

#: What each family leaves on the engine, spelled for the run's attrs
#: (``vi_guide_parameters``) so that "no fitted parameters were kept" is a
#: recorded fact rather than an absence a reader has to interpret. **W5.14.**
_GUIDE_PARAMETERS: dict[str, str] = {
    "normal": "loc, scale",
    "multivariate": "loc, scale_tril",
    "laplace": "loc, scale_tril",
    "flow": "none",
}

#: The smallest free dimension an autoregressive flow can be built over: with
#: one coordinate there is nothing to be autoregressive *about*, and both
#: libraries say so — numpyro raises ``ValueError("latent dim = 1. Consider
#: using AutoDiagonalNormal instead")`` and pyro's autoregressive network warns
#: and degenerates. :meth:`VIEngine.run` refuses first, by name, so that the
#: user reads ampere's diagnosis rather than a library's. **W5.14.**
_FLOW_MIN_DIM = 2

#: ``ampere_approximation`` (``results.md`` §9, W5.0): the family a plot or a
#: summary checks before it reports an R-hat that means nothing for a run
#: that was never a Markov chain. Named for what the guide *assumes* rather
#: than for its class, which is why this is a second table from
#: :data:`GUIDE_FAMILIES` rather than a re-use of it — ``vi_guide_class``
#: already carries the library's own name.
_APPROXIMATION_FAMILIES: dict[str, str] = {
    "normal": "mean_field",
    "multivariate": "multivariate",
    "laplace": "laplace",
    "flow": "normalising_flow",
}

#: The single site the whole unconstrained vector travels under. See the module
#: docstring; it is the same reasoning, and the same shape, as ``_nuts.py``'s
#: ``_PYRO_SITE``.
_SITE = "theta"

#: How many ELBO values a run records. The trace is a provenance *attribute*,
#: JSON-encoded into the archived file, so it must stay small whatever the step
#: count; 200 evenly spaced values describe the shape of any optimisation
#: anyone will read, and the final value is recorded separately and exactly.
_TRACE_POINTS = 200


@dataclasses.dataclass(frozen=True)
class _Settings:
    """One run's optimiser settings, so a route cannot be handed a partial set."""

    draws: int
    steps: int
    guide: str
    learning_rate: float
    progress: bool


def supported_backends() -> frozenset[str]:
    """The backends this driver can fit **right now**, in this interpreter.

    A backend qualifies when it has registered a realisation *and* this driver
    knows a variational library that can consume it. A function rather than a
    constant, for the reason :func:`ampere.inference._nuts.supported_backends`
    gives: a realisation is registered when the user imports the backend, which
    is normally after ``ampere.inference``.
    """
    return frozenset(registered_realisations()) & frozenset(VARIATIONAL_LIBRARIES)


class VIEngine(Engine):
    """Fit a :class:`~ampere.core.dataset.FittingProblem` by variational inference.

    Stochastic variational inference (Hoffman et al. 2013) over the
    *unconstrained* density: a parametrised guide is optimised to maximise the
    evidence lower bound, and the posterior is then drawn from the fitted
    guide. Fast where NUTS is slow, and approximate where NUTS is exact — see
    the module docstring on what the guide family assumes.

    Parameters
    ----------
    problem
        The composed problem. Must report a backend in
        :func:`supported_backends` and must be ``differentiable``; both are
        checked at construction, by name.
    density
        Optional, and normally omitted — as for
        :class:`~ampere.inference.NUTSEngine`. Omitted, the density is the
        backend's registered realisation (``ampere.core.realise``).
    cache_size
        See :class:`~ampere.inference.engine.Engine`. This driver's cache is
        only ever populated by its start-point search: the optimiser scores
        through the realised density, and the stored draws are decomposed from
        the realisation's own ``log_likelihood_terms`` when it offers one.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.torch import IndependentNoise, PowerLaw
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
    >>> from ampere.inference import VIEngine
    >>> grid = np.geomspace(1.0, 10.0, 24)
    >>> truth = 2.0 * grid ** -1.0
    >>> rng = np.random.default_rng(3)
    >>> observed = Spectrum(
    ...     grid * u.um,
    ...     (truth + rng.normal(0.0, 0.03, grid.size)) * u.Jy,
    ...     uncertainty=np.full(grid.size, 0.03) * u.Jy,
    ... )
    >>> problem = FittingProblem(
    ...     PowerLaw(grid, norm=st.norm(2.0, 0.5), index=-1.0),
    ...     [Dataset(observed, likelihood=Likelihood(GaussianFamily(), IndependentNoise()))],
    ...     seed=20260907,
    ... )
    >>> run = VIEngine(problem).run(draws=400, steps=1500)
    >>> run["posterior"]["model.norm"].shape
    (1, 400)
    >>> run.attrs["ampere_engine"], run.attrs["ampere_vi_guide"]
    ('vi', 'normal')
    >>> bool(abs(float(run["posterior"]["model.norm"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "vi"
    #: A variational fit is nothing but gradient descent on the ELBO.
    OFFERS_GRADIENTS: ClassVar[bool] = True

    def __init__(
        self,
        problem: FittingProblem,
        density: Any = None,
        *,
        cache_size: int = DEFAULT_CACHE_SIZE,
    ) -> None:
        # First, ahead of every other refusal: a foreign part is the most
        # specific diagnosis available and the only one that names the piece
        # (W3.8). Ordering it here also means it is the answer on a backend
        # that has no realisation registered at all, where "this driver cannot
        # run on your backend" would send the user to install something that
        # would not have helped.
        _refuse_foreign_parts(self.NAME, problem)
        available = supported_backends()
        if problem.backend not in available:
            known = ", ".join(sorted(available)) or "(none: no backend has been imported here)"
            raise EngineError(
                f"{self.NAME} cannot fit a problem on the {problem.backend!r} backend. This "
                f"driver needs two things: a registered realisation, so there is a "
                f"differentiable density at all (inference.md §10a), and a variational library "
                f"that can consume it — {VARIATIONAL_LIBRARIES}. Backends satisfying both in "
                f"this interpreter: {known}. Import the backend that supplies one and build the "
                f"problem from its models, instrument steps, noise models and GP solver; or use "
                f"a gradient-free engine (emcee, dynesty, zeus), which run every backend on the "
                f"contract path."
            )
        if not problem.differentiable:
            raise EngineError(
                f"{self.NAME} needs a differentiable problem, and this one declares "
                f"differentiable=False. That flag is aggregated from what the models, instrument "
                f"steps, noise models and GP solvers themselves declare, so one "
                f"non-differentiable piece is enough — problem.capabilities says which pieces "
                f"were consulted."
            )
        self.realisation: Realisation | None = None
        #: The fitted guide after a run — pyro's or numpyro's, whichever route
        #: ran — for anything this driver does not expose: the guide's own
        #: ``quantiles``, say. :attr:`sampler` holds the ``SVI`` object beside
        #: it, as it holds the ``MCMC`` object for the other drivers.
        self.guide: Any = None
        #: The fitted guide's own parameters, as plain arrays rather than
        #: library objects, so they outlive pyro's scoped param store and a
        #: caller (or a test) can rebuild the guide's density independently
        #: of which library fitted it. ``guide_scale`` is set for
        #: ``guide="normal"``, ``guide_scale_tril`` for
        #: ``guide="multivariate"`` and for ``guide="laplace"`` (whose fit is
        #: a full-covariance Gaussian too) — never both. **W5.0**, extended
        #: at **W5.14**: for ``guide="flow"`` all three stay ``None``,
        #: because a normalising flow's fit is a neural network's weights and
        #: not a location and a scale at all. A caller reading these must
        #: therefore check for ``None`` rather than assume a Gaussian; the
        #: run's ``ampere_vi_guide_parameters`` attribute records which of
        #: them a given run left behind, so an archived run says it too.
        self.guide_loc: np.ndarray | None = None
        self.guide_scale: np.ndarray | None = None
        self.guide_scale_tril: np.ndarray | None = None
        if density is None:
            try:
                self.realisation = realise(problem)
            except LoweringError as error:
                raise EngineError(
                    f"{self.NAME} could not obtain a differentiable density for this problem: "
                    f"{error}"
                ) from error
            # `use_realisation=False`: this driver already holds the
            # realisation and scores every proposal through it. The base
            # class's gradient-free fast path (W2.5 slice 3) would lower
            # the same problem a second time for the start-point search
            # alone, which is a lowering and a compilation for nothing.
            super().__init__(problem, cache_size=cache_size, use_realisation=False)
            self.density = self.realisation.log_prob_unconstrained
            return
        if not callable(density):
            raise EngineError(
                f"{self.NAME} takes the backend's lowering of the problem as a callable of the "
                f"unconstrained vector, got {density!r}. Omit it to use the backend's registered "
                f"realisation (ampere.core.realise)."
            )
        super().__init__(problem, cache_size=cache_size, use_realisation=False)
        self.density = density

    # -- the run --------------------------------------------------------------

    def run(
        self,
        draws: int,
        *,
        steps: int = 2000,
        guide: str = "normal",
        learning_rate: float = 0.01,
        initial: Any = None,
        progress: bool = False,
    ) -> Any:
        """Optimise the guide, draw from it, and emit the run.

        Parameters
        ----------
        draws
            Independent draws taken from the **fitted** guide. They are the
            run's posterior. Cheap — a guide draw is one reparametrised normal
            sample — so a large number costs almost nothing beside the fit.
        steps
            Optimiser steps. The ELBO trace in the attrs is how a reader tells
            whether this was enough; a trace still climbing at the last step
            means it was not.
        guide
            One of :data:`GUIDE_FAMILIES`: ``"normal"`` (mean-field,
            diagonal), ``"multivariate"`` (full-covariance, ELBO-fitted),
            ``"laplace"`` (full-covariance from the curvature at the MAP
            point) or ``"flow"`` (an inverse-autoregressive flow, the only
            non-Gaussian family, needing at least
            ``_FLOW_MIN_DIM`` free parameters). See the module docstring on
            what each costs: this argument *is* the approximation being made.
        learning_rate
            Adam's step size.
        initial
            A ``(1, n_dim)`` or ``(n_dim,)`` start point **in the constrained
            space**, as for the other drivers. The default draws one from the
            joint prior on this engine's own initialisation stream, so a run
            repeats exactly from the problem's seed.
        progress
            Report progress while optimising — the ELBO every tenth of the
            run on the pyro route, numpyro's own progress bar on the jax one,
            because each library's idiom is what a user of that library
            expects to see. Off by default either way: a driver that prints by
            default is unusable inside a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(1, draws)``.
        """
        if int(draws) < 1:
            raise EngineError(f"{self.NAME} needs at least one draw, got {draws}.")
        if int(steps) < 1:
            raise EngineError(f"{self.NAME} needs at least one optimiser step, got {steps}.")
        if guide not in GUIDE_FAMILIES:
            known = ", ".join(sorted(GUIDE_FAMILIES))
            raise EngineError(
                f"{self.NAME} does not know the guide family {guide!r}. Available: {known}. "
                f"'normal' assumes the posterior factorises over parameters; 'multivariate' "
                f"captures their correlations at O(d**2) guide parameters; 'laplace' takes the "
                f"same full covariance from the curvature at the MAP point rather than from an "
                f"ELBO fit; 'flow' is the only one of the four that is not a Gaussian."
            )
        if guide == "flow" and int(self.problem.free_size) < _FLOW_MIN_DIM:
            raise EngineError(
                f"{self.NAME}'s 'flow' guide is an autoregressive flow over the free parameters, "
                f"and this problem has {self.problem.free_size} of them: there is nothing for the "
                f"flow to be autoregressive over below {_FLOW_MIN_DIM}. Use guide='normal' "
                f"(identical to 'multivariate' in one dimension) or guide='laplace', both of "
                f"which fit a one-dimensional posterior exactly when it is Gaussian."
            )
        if float(learning_rate) <= 0.0:
            raise EngineError(
                f"{self.NAME}'s learning rate must be positive, got {learning_rate!r}."
            )
        self.start()

        start = self.initial_positions(1)[0] if initial is None else self._checked_initial(initial)
        unconstrained = self.problem.unconstrain(start)

        settings = _Settings(
            draws=int(draws),
            steps=int(steps),
            guide=str(guide),
            learning_rate=float(learning_rate),
            progress=bool(progress),
        )
        library = VARIATIONAL_LIBRARIES[self.problem.backend]
        if library == "pyro":
            drawn, log_q, attrs = self._fit_pyro(unconstrained, settings)
        else:
            drawn, log_q, attrs = self._fit_numpyro(unconstrained, settings)

        chain = np.stack([[self.problem.constrain(y) for y in drawn]])
        # The guide's own density is fitted in *unconstrained* coordinates
        # (module docstring); the stored log_prior/log_likelihood are in the
        # *constrained* ones (`_evaluations_from_terms` scores `lnprior` at
        # `chain`, not at `drawn`). Moving log_q into the same coordinates is
        # what makes sample_stats.proposal_log_density usable directly in
        # results.md §9's importance-weight formula (engine.py, W5.0).
        proposal_log_density = log_q - unconstrained_jacobian_correction(self.problem, drawn)
        attrs.update(
            {
                "vi_draws": settings.draws,
                "vi_steps": settings.steps,
                "vi_guide": settings.guide,
                "vi_guide_class": GUIDE_FAMILIES[settings.guide],
                # W5.14: which fitted parameters this run left on the engine.
                # A flow leaves none, and an archived run has to say so rather
                # than leave a reader to infer it from the guide's name.
                "vi_guide_parameters": _GUIDE_PARAMETERS[settings.guide],
                "vi_optimiser": "adam",
                "vi_learning_rate": settings.learning_rate,
                "vi_library": library,
                "approximation": _APPROXIMATION_FAMILIES[settings.guide],
            }
        )
        decomposition = self._decomposition(np.asarray([drawn]))
        return self.finish(
            chain,
            extra_attrs=attrs,
            realised=self.realisation is not None,
            registered_lowerings=self._realisation_provenance(),
            log_likelihood_terms=decomposition,
            sample_stats={"proposal_log_density": proposal_log_density},
        )

    # -- the pyro route -------------------------------------------------------

    def _fit_pyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
        """pyro's SVI over a torch realisation.

        Lazy imports, exactly as ``_nuts.py`` imports pyro: ``import
        ampere.inference`` must not require torch (``architecture.md`` §4 rule
        2). The suppressions are the price of that rule — ``dev``, the
        environment CI typechecks in, deliberately has no torch.

        Three things are deliberate and are the reasons this method is longer
        than a call to ``SVI``.

        * **The parameter store is scoped.** pyro keeps guide parameters in a
          process-global ``ParamStoreDict``; a fit that wrote into it would
          leave state behind for the next fit in the same interpreter to pick
          up under the same site names, which is a silent wrong answer rather
          than a crash. ``scope()`` gives this run its own, exactly as
          ``fork_rng`` gives it its own RNG.
        * **The seed is forked, not set.** pyro's guide draws its
          reparametrisation noise from torch's *global* generator and offers no
          ``generator=``, so this is ``lowering.md`` §9.1's route (1) — seed the
          global stream from this engine's own sub-stream — done inside
          ``torch.random.fork_rng`` so the host process's RNG state is restored
          afterwards. Same rule as the NUTS route.
        * **The guide is initialised at the start point.** ``AutoNormal`` and
          ``AutoMultivariateNormal`` take an ``init_loc_fn``; without one they
          start at zero in unconstrained space, which for a
          ``loguniform``-declared parameter is the geometric middle of a very
          wide prior and can be far enough from the mass that the ELBO's
          gradient is numerically flat. Starting where the other drivers start
          — a seeded draw from the joint prior — makes the run reproducible from
          the problem's seed *and* comparable with an emcee or NUTS run of the
          same problem.

        A fourth thing is new at W5.0: **the guide's own log-density is read
        off the same trace the draw comes from**, inside the scope, because
        the fitted parameters :meth:`~pyro.params.param_store.ParamStoreDict.
        scope` isolated are gone the moment this ``with`` block exits — asking
        for ``guide.log_prob`` afterwards would look up parameters that no
        longer exist in the (now-restored) global store.

        **It is summed over every sample site in the trace, not read off
        ``_SITE`` alone — this is not a stylistic choice, it is the fix for a
        real bug a review caught.** An autoguide does not sample ``_SITE``
        directly: it draws from an *auxiliary* site (``AutoNormal``'s own
        unconstrained latent, ``AutoMultivariateNormal``'s
        ``..._latent``) under the real ``Normal``/``MultivariateNormal`` the
        guide fitted, then records ``_SITE`` itself as a :class:`pyro.
        distributions.Delta` at the transformed value — here the identity
        transform, since the carrier's support is already ``R**d``
        (module docstring). ``trace.nodes[_SITE]["fn"].log_prob(...)`` is
        therefore that ``Delta``'s log-density, which is the change-of-
        variables term between the auxiliary site and ``_SITE`` — **zero**
        for an identity transform, not the guide's density. Reading it alone
        silently stored a constant zero for every draw, whatever the guide
        actually fitted. The auxiliary site's own ``fn.log_prob`` *is* the
        real density, but its distribution and value are only equal to
        ``_SITE``'s **up to** exactly the ``Delta`` term the two together
        cancel — so summing ``fn.log_prob(value)`` over **every** sample-type
        node in the trace gives the guide's total log-density at ``_SITE``'s
        drawn value by construction, without this method needing to know
        which auxiliary site name a given autoguide subclass happens to use.
        """
        import pyro  # pyrefly: ignore[missing-import]
        import pyro.poutine as poutine  # pyrefly: ignore[missing-import]
        import torch  # pyrefly: ignore[missing-import]
        from pyro.distributions import Normal  # pyrefly: ignore[missing-import]
        from pyro.infer import SVI, Trace_ELBO  # pyrefly: ignore[missing-import]
        from pyro.infer import autoguide  # pyrefly: ignore[missing-import]
        from pyro.optim import Adam  # pyrefly: ignore[missing-import]

        size = int(self.problem.free_size)
        start = torch.as_tensor(unconstrained, dtype=torch.float64)
        zeros = torch.zeros(size, dtype=torch.float64)
        ones = torch.ones(size, dtype=torch.float64)

        def model() -> None:
            # The carrier site and the correction that cancels it; see the
            # module docstring. `.to_event(1)` makes the d coordinates one
            # multivariate site rather than d univariate ones, which is what
            # lets AutoMultivariateNormal correlate them.
            base = Normal(zeros, ones).to_event(1)
            theta = pyro.sample(_SITE, base)
            pyro.factor("ampere_density", self.density(theta) - base.log_prob(theta))

        builder = getattr(autoguide, GUIDE_FAMILIES[settings.guide])
        elbo: list[float] = []
        with pyro.get_param_store().scope(), torch.random.fork_rng(devices=[]):
            torch.manual_seed(self.integer_seed("optimiser"))
            if settings.guide == "flow":
                # W5.14. `AutoIAFNormal` ignores `init_loc_fn` and warns if it
                # is given one; a flow is started by moving its *base*, which
                # is what `get_base_dist` exists for. See the module docstring
                # on why starting it at the origin is not an option in
                # ampere's unconstrained coordinates. The subclass is built
                # here, beside the lazy import, because naming pyro's class at
                # module level is exactly what `architecture.md` §4 rule 2
                # forbids.
                started = Normal(start, ones).to_event(1)

                class _StartedFlow(builder):  # type: ignore[misc, valid-type]
                    """``AutoIAFNormal``, started at the start point and in float64."""

                    def get_base_dist(self) -> Any:
                        return started

                    def get_posterior(self, *args: Any, **kwargs: Any) -> Any:
                        # pyro builds the flow's autoregressive network the
                        # first time a posterior is asked for, at torch's
                        # *default* dtype -- float32, while everything else in
                        # this driver is float64, and the two meet as
                        # "mat1 and mat2 must have the same dtype" the moment
                        # the float64 base reaches the network. Building it
                        # here and casting before the assignment registers its
                        # parameters is the float64 that `lowering.md` §10.1
                        # asks to be threaded explicitly rather than obtained
                        # from `torch.set_default_dtype`, which this project
                        # never calls. pyro's own `get_posterior` then finds
                        # the transform already built and leaves it alone.
                        if self.transform is None:
                            self.transform = self._init_transform_fn(self.latent_dim).double()
                        return super().get_posterior(*args, **kwargs)

                guide = _StartedFlow(model)
            else:
                guide = builder(model, init_loc_fn=lambda site: start)
            svi = SVI(model, guide, Adam({"lr": settings.learning_rate}), loss=Trace_ELBO())
            report = max(1, settings.steps // 10)
            for step in range(settings.steps):
                # SVI.step returns the *loss*, which is the negated ELBO.
                elbo.append(-float(svi.step()))
                if settings.progress and step % report == 0:
                    print(f"{self.NAME}: step {step}/{settings.steps}  ELBO {elbo[-1]:.6g}")
            # W5.14, and the one place the laplace family is not just another
            # autoguide: what SVI fitted is a `Delta` at the MAP point, and
            # drawing from *that* would return the same point `draws` times
            # with a meaningless density. `laplace_approximation()` is pyro's
            # own documented second step — it takes the Hessian of the
            # negative log joint at the fitted `loc` and hands back an
            # `AutoMultivariateNormal` carrying `loc`, `scale` and
            # `scale_tril` as buffers. From here on it *is* an
            # `AutoMultivariateNormal`, which is why the draw loop and the
            # parameter extraction below need no further special case. It
            # must run outside `no_grad` (it differentiates twice) and inside
            # the param-store scope (it reads the fitted `loc`).
            fitted = guide.laplace_approximation() if settings.guide == "laplace" else guide
            self.guide = fitted
            self.sampler = svi
            with torch.no_grad():
                values: list[np.ndarray] = []
                log_q = np.empty(settings.draws, dtype=float)
                for i in range(settings.draws):
                    trace = poutine.trace(fitted).get_trace()
                    node = trace.nodes[_SITE]
                    values.append(np.asarray(node["value"].detach().cpu().numpy(), dtype=float))
                    # Every sample-type node, Delta included: see the
                    # docstring's account of why the sum -- not `_SITE`
                    # alone -- is the guide's density at `_SITE`'s value.
                    log_q[i] = float(
                        sum(
                            site["fn"].log_prob(site["value"]).sum()
                            for site in trace.nodes.values()
                            if site["type"] == "sample"
                        )
                    )
                drawn = np.stack(values)
                # Kept as plain arrays (not the guide object, which is
                # useless once the param-store scope below exits) so a
                # caller -- or a test checking this method against an
                # independent scipy density -- can read the fitted guide
                # back without touching pyro's param store at all.
                #
                # The two autoguide classes expose their fit differently, and
                # neither name is `.loc`/`.scale` on both: `AutoNormal` is
                # not an `AutoContinuous` subclass and keeps one `Parameter`
                # per *site name* under `.locs`/`.scales` (`guide.locs.theta`
                # here); `AutoMultivariateNormal` is, and exposes `.loc`
                # directly but its **covariance's** Cholesky factor is not
                # `.scale_tril` alone -- that is a unit-diagonal correlation
                # matrix, row-scaled by the separate `.scale` vector
                # (pyro's own `get_posterior`: `scale[..., None] *
                # scale_tril`). Combining the two here, once, is what let a
                # review's own probe catch this method reconstructing the
                # wrong covariance from `scale_tril` alone.
                # `laplace` joins `multivariate` here rather than needing a
                # branch of its own: `laplace_approximation()` above returned
                # an `AutoMultivariateNormal`, and pyro builds it in exactly
                # the same row-scaled-correlation split (`register_buffer`
                # of `loc`, `scale` and a unit-diagonal `scale_tril`), so the
                # same product is the same covariance's Cholesky factor.
                # `flow` reaches neither branch: there is no location and no
                # scale to keep, and all three attributes stay `None`
                # (W5.14; `GAUSSIAN_GUIDES` is the predicate).
                if settings.guide in {"multivariate", "laplace"}:
                    self.guide_loc = np.asarray(fitted.loc.detach().cpu().numpy(), dtype=float)
                    scale = fitted.scale.detach().cpu().numpy()
                    correlation = fitted.scale_tril.detach().cpu().numpy()
                    self.guide_scale_tril = np.asarray(scale[..., None] * correlation, dtype=float)
                elif settings.guide == "normal":
                    self.guide_loc = np.asarray(
                        getattr(fitted.locs, _SITE).detach().cpu().numpy(), dtype=float
                    )
                    self.guide_scale = np.asarray(
                        getattr(fitted.scales, _SITE).detach().cpu().numpy(), dtype=float
                    )

        attrs: dict[str, object] = {
            "vi_final_elbo": elbo[-1],
            "vi_elbo_trace": _thinned(elbo),
            "pyro_version": str(pyro.__version__),
            "torch_version": str(torch.__version__),
        }
        return drawn, log_q, attrs

    # -- the numpyro route ----------------------------------------------------

    def _fit_numpyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
        """numpyro's SVI over a jax realisation.

        The same three steps as the pyro route, in numpyro's spelling, and
        emitting the same provenance keys — ``vi_final_elbo``,
        ``vi_elbo_trace`` and a pair of library versions — because a reader of
        an archived run should not have to know which library fitted it to
        know what was fitted.

        Lazy imports, exactly as ``_nuts.py`` imports numpyro: ``import
        ampere.inference`` must not require jax (``architecture.md`` §4
        rule 2). The suppressions are the price of that rule — ``dev``, the
        environment CI typechecks in, deliberately has no jax — and the real
        check is ``pixi run -e jax typecheck``.

        Two differences from the pyro route, both forced by the library rather
        than chosen.

        * **The carrier is an improper uniform, not a standard normal.** pyro
          needs a proper site and subtracts its ``log_prob`` back off; numpyro
          offers ``ImproperUniform``, whose contribution to the log-density is
          *identically* zero, so the correction has nothing to cancel and is
          simply absent. That is the better of the two — there is no
          difference of two large numbers to lose precision to — and it is
          available here only because numpyro ships the distribution. The
          support is ``real_vector``, so the autoguide's own ``biject_to`` is
          the identity and the guide is fitted in exactly the coordinates the
          realisation takes, which is the property the pyro route obtains from
          the standard normal.
        * **Nothing is scoped or forked.** numpyro has no process-global
          parameter store and no global RNG: the fitted parameters come back
          in the ``SVIRunResult`` and every stream is an explicit key. So the
          two guards the pyro route needs — ``get_param_store().scope()`` and
          ``fork_rng`` — have no counterpart, and their absence is a property
          of the library rather than an omission. The seed still comes from
          this engine's own sub-stream, so a run repeats from the problem's
          seed exactly as every other driver's does.

        The guide's own log-density (W5.0) is read off a **traced** draw, not
        off ``get_posterior`` — a review's own probe of this environment
        found ``AutoNormal`` has no such method at all
        (``'AutoNormal' object has no attribute 'get_posterior'``): unlike
        ``AutoMultivariateNormal``, numpyro's ``AutoNormal`` is not an
        ``AutoContinuous`` subclass, so the two guide classes this driver
        supports do not share one convenience accessor. Tracing does not need
        one: :func:`numpyro.handlers.substitute` fixes the fitted parameters,
        :func:`numpyro.handlers.seed` gives the draw its own key, and summing
        ``fn.log_prob(value)`` over every ``type == "sample"`` node in the
        resulting trace is the guide's total log-density at ``_SITE``'s drawn
        value regardless of which internal shape a given autoguide happens
        to use — the same reasoning, and the same fix for the same class of
        bug, as the pyro route's identical sum. (``AutoNormal`` samples
        ``_SITE`` directly under a real ``Normal``, so the sum there has one
        term; ``AutoMultivariateNormal`` samples an auxiliary
        ``..._latent`` under the real ``MultivariateNormal`` and reports
        ``_SITE`` as a zero-density ``Delta`` at the identity-transformed
        value, so the sum there has two, one of them zero — either way the
        total is correct without this method needing to know which.)
        """
        import jax  # pyrefly: ignore[missing-import]
        import numpyro  # pyrefly: ignore[missing-import]
        import numpyro.distributions as dist  # pyrefly: ignore[missing-import]
        from numpyro import optim  # pyrefly: ignore[missing-import]
        from numpyro.infer import SVI, Trace_ELBO, autoguide  # pyrefly: ignore[missing-import]
        from numpyro.infer.initialization import (  # pyrefly: ignore[missing-import]
            init_to_value,
        )

        size = int(self.problem.free_size)
        density = self.density
        start = jax.numpy.asarray(unconstrained, dtype=jax.numpy.float64)

        def model() -> None:
            theta = numpyro.sample(
                _SITE, dist.ImproperUniform(dist.constraints.real_vector, (), (size,))
            )
            numpyro.factor("ampere_density", density(theta))

        builder = getattr(autoguide, GUIDE_FAMILIES[settings.guide])
        # The same start point the pyro route uses, and for the same reason:
        # a guide left at zero in unconstrained space can begin far enough
        # from the mass that the ELBO's gradient is numerically flat. The flow
        # is started the same way the pyro route starts it and for the same
        # reason -- by its base rather than by a location parameter it does
        # not have (W5.14; the module docstring has the measurement).
        if settings.guide == "flow":
            started = dist.Normal(start, jax.numpy.ones(size)).to_event(1)

            class _StartedFlow(builder):  # type: ignore[misc, valid-type]
                """``AutoIAFNormal`` whose base sits at the start point."""

                def get_base_dist(self) -> Any:
                    return started

            guide = _StartedFlow(model)
        else:
            guide = builder(model, init_loc_fn=init_to_value(values={_SITE: start}))
        svi = SVI(model, guide, optim.Adam(settings.learning_rate), loss=Trace_ELBO())
        keys = jax.random.split(jax.random.key(self.integer_seed("optimiser")), 2)
        result = svi.run(keys[0], settings.steps, progress_bar=settings.progress)
        self.guide = guide
        self.sampler = svi

        if settings.guide == "laplace":
            drawn, log_q = self._draw_numpyro_laplace(guide, result.params, settings, keys[1])
        else:
            drawn, log_q = self._draw_numpyro_traced(
                numpyro, guide, result.params, settings, keys[1], size
            )

        # Kept as plain arrays, exactly as the pyro route keeps
        # `guide_loc`/`guide_scale`/`guide_scale_tril`, so a caller can
        # rebuild the guide's density independently of which library fitted
        # it. numpyro's own naming for the two classes differs (`AutoNormal`
        # keys one param per site name, `f"{site}_{prefix}_loc"`;
        # `AutoMultivariateNormal` keys the flattened fit as
        # `f"{prefix}_loc"`/`f"{prefix}_scale_tril"` — the Cholesky factor of
        # the covariance directly, unlike pyro's row-scaled-correlation
        # split above).
        # `laplace` is read off the fitted transform rather than off a
        # parameter, because the covariance is not a parameter: SVI fitted
        # only the MAP location, and the Cholesky factor is the inverse
        # Hessian computed there (W5.14). `flow` keeps nothing, for the
        # reason `GAUSSIAN_GUIDES` records.
        if settings.guide == "laplace":
            transform = guide.get_transform(result.params)
            self.guide_loc = np.asarray(transform.loc, dtype=float)
            self.guide_scale_tril = np.asarray(transform.scale_tril, dtype=float)
        elif settings.guide == "multivariate":
            self.guide_loc = np.asarray(result.params[f"{guide.prefix}_loc"], dtype=float)
            self.guide_scale_tril = np.asarray(
                result.params[f"{guide.prefix}_scale_tril"], dtype=float
            )
        elif settings.guide == "normal":
            self.guide_loc = np.asarray(result.params[f"{_SITE}_{guide.prefix}_loc"], dtype=float)
            self.guide_scale = np.asarray(
                result.params[f"{_SITE}_{guide.prefix}_scale"], dtype=float
            )

        # numpyro's SVI *minimises* the negative ELBO, so its losses are -ELBO.
        # Recorded as the ELBO itself, which is what the pyro route records and
        # what a reader expects to see climbing towards a plateau.
        elbo = [-float(value) for value in np.asarray(result.losses, dtype=float)]
        attrs: dict[str, object] = {
            "vi_final_elbo": elbo[-1],
            "vi_elbo_trace": _thinned(elbo),
            "numpyro_version": str(numpyro.__version__),
            "jax_version": str(jax.__version__),
        }
        return drawn, log_q, attrs

    def _draw_numpyro_traced(
        self,
        numpyro: Any,
        guide: Any,
        params: Any,
        settings: _Settings,
        key: Any,
        size: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Draw from a fitted numpyro guide by **tracing** it, and read its density off the trace.

        The route every family but ``laplace`` takes, unchanged since W5.0
        and for W5.0's reason: :func:`numpyro.handlers.substitute` fixes the
        fitted parameters, :func:`numpyro.handlers.seed` gives the draw its
        own key, and summing ``fn.log_prob(value)`` over every ``type ==
        "sample"`` node is the guide's total log-density at ``_SITE``'s drawn
        value whatever internal shape the autoguide happens to use. The flow
        joins the two Gaussian families here without a word of new code,
        which is the property that makes the sum worth having: its auxiliary
        latent is sampled under a
        :class:`~numpyro.distributions.TransformedDistribution` whose
        ``log_prob`` already carries the flow's log-Jacobian, and ``_SITE``
        is the same zero-density ``Delta`` it is for
        ``AutoMultivariateNormal``.
        """
        import jax  # pyrefly: ignore[missing-import]

        draw_keys = jax.random.split(key, settings.draws)
        values: list[np.ndarray] = []
        log_q = np.empty(settings.draws, dtype=float)
        for i in range(settings.draws):
            traced = numpyro.handlers.trace(
                numpyro.handlers.seed(
                    numpyro.handlers.substitute(guide, data=params),
                    rng_seed=draw_keys[i],
                )
            ).get_trace()
            site = traced[_SITE]
            values.append(np.asarray(site["value"], dtype=float))
            # Every sample-type node, Delta included: see the docstring's
            # account of why the sum is the guide's density at `_SITE`'s
            # value regardless of which autoguide class produced the trace.
            log_q[i] = float(
                sum(
                    np.asarray(s["fn"].log_prob(s["value"])).sum()
                    for s in traced.values()
                    if s["type"] == "sample"
                )
            )
        return np.asarray(values, dtype=float).reshape(settings.draws, size), log_q

    def _draw_numpyro_laplace(
        self, guide: Any, params: Any, settings: _Settings, key: Any
    ) -> tuple[np.ndarray, np.ndarray]:
        """Draw from the **Laplace** Gaussian, which is not what tracing the guide would give.

        **W5.14.** numpyro's ``AutoLaplaceApproximation`` is a ``Delta`` guide
        during SVI — it fits the MAP location and nothing else — so the traced
        route above would return the same point ``draws`` times and read a
        ``Delta``'s log-density for it. The Gaussian this family is named for
        exists only in :meth:`~numpyro.infer.autoguide.
        AutoLaplaceApproximation.get_posterior`, which takes the Hessian of
        the negative log joint at the fitted location and returns the
        multivariate normal whose covariance is its inverse. Drawing from
        *that* distribution, and asking *it* for the density, is therefore not
        a shortcut around the trace: it is the only place the approximation
        is. (pyro reaches the same object by a different road —
        ``laplace_approximation()`` hands back a whole
        ``AutoMultivariateNormal``, which the pyro route then traces like any
        other — and the two roads agree because the distribution at the end of
        them is the same one.)

        The distribution is over the *latent* vector, and ``_SITE``'s support
        is ``real_vector``, so the bijection between the two is the identity
        (module docstring) and the drawn latent is the unconstrained vector
        the rest of this driver expects, with no change-of-variables term to
        add.
        """
        posterior = guide.get_posterior(params)
        sample = posterior.sample(key, (settings.draws,))
        drawn = np.asarray(sample, dtype=float)
        log_q = np.asarray(posterior.log_prob(sample), dtype=float)
        return drawn, log_q

    # -- shared plumbing ------------------------------------------------------

    def _realisation_provenance(self) -> list[dict[str, Any]] | None:
        """The user-registered lowering rows the realisation consulted, if any."""
        if self.realisation is None:
            return None
        rows = getattr(self.realisation, "lowering_provenance", None)
        return list(rows()) if callable(rows) else None

    def _decomposition(self, unconstrained: np.ndarray) -> list[list[dict[str, float]]] | None:
        """§10a's optional per-dataset split, for the drawn points.

        The same fetch-don't-require rule as the NUTS driver's: ``None`` when
        the realisation offers no decomposition, in which case
        :meth:`Engine.finish` falls back to the numpy contract path and counts
        what it recomputed.
        """
        if self.realisation is None:
            return None
        terms = log_likelihood_terms_of(self.realisation)
        if terms is None:
            return None
        return [
            [{label: float(np.asarray(value)) for label, value in terms(y).items()} for y in chain]
            for chain in unconstrained
        ]

    def _checked_initial(self, initial: Any) -> np.ndarray:
        positions = np.asarray(initial, dtype=float).reshape(-1)
        if positions.size != self.problem.free_size:
            raise EngineError(
                f"{self.NAME}'s start point must have {self.problem.free_size} element(s), got "
                f"{positions.size}."
            )
        return positions


def _thinned(values: list[float]) -> list[float]:
    """*values* reduced to at most :data:`_TRACE_POINTS`, endpoints kept.

    The trace is stored as a provenance attribute, so it is JSON-encoded into
    the archived file and has to stay small whatever the step count. Evenly
    spaced indices rather than a decimation of the tail: the *shape* of the
    optimisation is what a reader is after, and the last value is recorded
    exactly in ``vi_final_elbo`` regardless.
    """
    if len(values) <= _TRACE_POINTS:
        return [float(value) for value in values]
    indices = np.unique(np.linspace(0, len(values) - 1, _TRACE_POINTS).astype(int))
    return [float(values[index]) for index in indices]
