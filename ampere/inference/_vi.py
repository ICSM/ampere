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
tell what was assumed.

How it reaches a backend
------------------------
Exactly as :class:`~ampere.inference.NUTSEngine` does, and for the same
reason: ``ampere.inference`` may import ``ampere.core`` and ``ampere.results``
and nothing else, so the gradient comes from the backend's **realisation**
(``inference.md`` §10a) obtained through :func:`ampere.core.realise`, and the
variational library is imported lazily inside :meth:`VIEngine.run`.

Turning a density into a model pyro can guide
----------------------------------------------
This is the one piece of real work, and it is worth stating plainly because it
looks like a trick. ``pyro.infer.SVI`` takes a *model* — a probabilistic
program with sample sites — while a realisation is a bare log density of one
flat unconstrained vector, with no site structure at all (``inference.md``
§10a fixes it that way, so that ampere's own ``ParameterMapping`` stays the
only thing that routes names). MCMC has ``potential_fn=`` for exactly this;
SVI has no equivalent.

So the model declares one site over the whole unconstrained vector — a
standard normal, whose support is all of ``R**d`` — and then adds a
:func:`pyro.factor` of ``density(theta) - base.log_prob(theta)``. The log
joint is therefore ``log N(theta; 0, I) + density(theta) - log N(theta; 0, I)``
= ``density(theta)``, exactly, and the standard normal is a *carrier* for the
site rather than a prior: it cancels identically, term by term, at every point.
Any distribution with full real support would do; the standard normal is
chosen because it is cheap and because its scale matches the unconstrained
coordinates a well-behaved bijection produces, which keeps the subtraction from
being a difference of two large numbers.

The autoguides then see a single unconstrained real site of the right shape,
which is what they are best at. ``AutoNormal`` gives a diagonal Gaussian in
that space, ``AutoMultivariateNormal`` a full-covariance one; both are
*unconstrained-space* approximations, so the constrained posterior they imply
is already correctly warped by ampere's own bijections rather than by a second
set of pyro's.

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

from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError

__all__ = ["GUIDE_FAMILIES", "VARIATIONAL_LIBRARIES", "VIEngine", "supported_backends"]

#: Which variational library drives which backend's realisation. Strings, not
#: imports, for the reason ``_nuts.py``'s :data:`SAMPLER_LIBRARIES` gives: this
#: module imports no backend and no library at module level.
#:
#: Only torch today. numpyro's ``SVI`` is the jax counterpart and joins by
#: adding a row here *and* a route in :meth:`VIEngine._fit`; the jax track owns
#: that, and a row here with no route would make ``supported_backends()`` claim
#: a capability that does not exist.
VARIATIONAL_LIBRARIES: dict[str, str] = {"torch": "pyro"}

#: The guide families a user may name, and what each assumes. The **names** are
#: ampere's, not pyro's, because the same two families exist in numpyro under
#: the same names and a run's ``vi_guide`` attribute has to mean one thing
#: across backends — the same discipline W2.12 applied to backend names.
GUIDE_FAMILIES: dict[str, str] = {
    "normal": "AutoNormal",
    "multivariate": "AutoMultivariateNormal",
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
        #: pyro's ``SVI`` object after a run, and the fitted guide, for anything
        #: this driver does not expose — the guide's own ``quantiles``, say.
        self.guide: Any = None
        if density is None:
            try:
                self.realisation = realise(problem)
            except LoweringError as error:
                raise EngineError(
                    f"{self.NAME} could not obtain a differentiable density for this problem: "
                    f"{error}"
                ) from error
            super().__init__(problem, cache_size=cache_size)
            self.density = self.realisation.log_prob_unconstrained
            return
        if not callable(density):
            raise EngineError(
                f"{self.NAME} takes the backend's lowering of the problem as a callable of the "
                f"unconstrained vector, got {density!r}. Omit it to use the backend's registered "
                f"realisation (ampere.core.realise)."
            )
        super().__init__(problem, cache_size=cache_size)
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
            ``"normal"`` (mean-field, diagonal) or ``"multivariate"``
            (full-covariance). See :data:`GUIDE_FAMILIES` and the module
            docstring: this argument *is* the approximation being made.
        learning_rate
            Adam's step size.
        initial
            A ``(1, n_dim)`` or ``(n_dim,)`` start point **in the constrained
            space**, as for the other drivers. The default draws one from the
            joint prior on this engine's own initialisation stream, so a run
            repeats exactly from the problem's seed.
        progress
            Print the ELBO every tenth of the run. Off by default: a driver
            that prints by default is unusable inside a loop or a test suite.

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
                f"captures their correlations at O(d**2) guide parameters."
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
        drawn, attrs = self._fit_pyro(unconstrained, settings)

        chain = np.stack([[self.problem.constrain(y) for y in drawn]])
        attrs.update(
            {
                "vi_draws": settings.draws,
                "vi_steps": settings.steps,
                "vi_guide": settings.guide,
                "vi_guide_class": GUIDE_FAMILIES[settings.guide],
                "vi_optimiser": "adam",
                "vi_learning_rate": settings.learning_rate,
                "vi_library": VARIATIONAL_LIBRARIES[self.problem.backend],
            }
        )
        decomposition = self._decomposition(np.asarray([drawn]))
        return self.finish(
            chain,
            extra_attrs=attrs,
            realised=self.realisation is not None,
            registered_lowerings=self._realisation_provenance(),
            log_likelihood_terms=decomposition,
        )

    # -- the pyro route -------------------------------------------------------

    def _fit_pyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, dict[str, object]]:
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
        """
        import pyro  # pyrefly: ignore[missing-import]
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
            guide = builder(model, init_loc_fn=lambda site: start)
            svi = SVI(model, guide, Adam({"lr": settings.learning_rate}), loss=Trace_ELBO())
            report = max(1, settings.steps // 10)
            for step in range(settings.steps):
                # SVI.step returns the *loss*, which is the negated ELBO.
                elbo.append(-float(svi.step()))
                if settings.progress and step % report == 0:
                    print(f"{self.NAME}: step {step}/{settings.steps}  ELBO {elbo[-1]:.6g}")
            self.guide = guide
            self.sampler = svi
            with torch.no_grad():
                drawn = np.stack(
                    [
                        np.asarray(guide()[_SITE].detach().cpu().numpy(), dtype=float)
                        for _ in range(settings.draws)
                    ]
                )

        attrs: dict[str, object] = {
            "vi_final_elbo": elbo[-1],
            "vi_elbo_trace": _thinned(elbo),
            "pyro_version": pyro.__version__,
            "torch_version": torch.__version__,
        }
        return drawn, attrs

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
