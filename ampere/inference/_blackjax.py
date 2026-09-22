"""The blackjax route on jax: MCLMC as a sampler, Pathfinder twice over.

Private module; the class is :class:`ampere.inference.BlackjaxEngine`.
**W5.14**, from ``docs/design/inference_extensions_memo.md`` §2.1's finding:
"``blackjax`` wants ``logdensity_fn(θ)`` and nothing else — the realisation is
that function", and one dependency then buys six methods. Two of them are
here, because two of them are the ones ampere does not already have.

Why these two and not the other four
------------------------------------
blackjax ships NUTS and HMC as well, and neither is here: numpyro's NUTS is
already :class:`~ampere.inference.NUTSEngine`'s jax route, and a second
implementation of the same algorithm would be a maintenance cost with no new
answer. What blackjax has that nothing else in this package has is:

* **MCLMC** — microcanonical Langevin Monte Carlo (Robnik & Seljak 2023). A
  deterministic isokinetic trajectory with periodic momentum decoherence
  rather than a Metropolis accept/reject: no rejections, no U-turn criterion,
  and a cost per effective sample reported at several times better than NUTS's
  on smooth high-dimensional posteriors. Those are exactly ampere's hard
  cases: a latent GP over a spectrum, a hierarchical plate over a population.
* **Pathfinder** (Zhang et al. 2022) — quasi-Newton on the log density, with
  a Gaussian approximation built from the L-BFGS inverse-Hessian factors at
  every iterate along the path and the ELBO-best one kept. It is used here
  **twice**, which is the memo's own reading of it: as a cheap approximate
  posterior in its own right, and as the best available *initialiser* for a
  sampler that has to start somewhere better than a prior draw.

One engine, two methods, and why they are not two classes
----------------------------------------------------------
:class:`BlackjaxEngine` takes ``method="mclmc"`` or ``method="pathfinder"``,
as ``inference_extensions_memo.md`` §2.1 proposed ("one
``BlackjaxEngine(problem, method=...)`` over the jax realisation"). The two
produce genuinely different runs — a Markov chain and a set of i.i.d. draws
from an approximation — but they share the thing that is worth sharing: one
extra, one lazy import, one library-version record, one refusal for a
non-jax problem, and one account of how a jax array becomes an ampere run.
Splitting them into two classes would duplicate all of it to separate two
``run`` bodies.

What each one writes into a run, and the two attributes that differ
--------------------------------------------------------------------
``ampere_approximation`` (``results.md`` §9) is the difference, and it is
not cosmetic:

* **Pathfinder writes ``"pathfinder"``**, and with it
  ``sample_stats.proposal_log_density`` — the approximation's own log-density
  at each draw, moved into the constrained coordinates the stored
  ``log_prior``/``log_likelihood`` are in, exactly as
  :class:`~ampere.inference.VIEngine` does and for exactly the same reason:
  ``exp(log_prior + log_likelihood - proposal_log_density)`` is then an
  importance weight that can be built from the stored groups alone. blackjax
  hands the density back beside the draws (``blackjax.vi.pathfinder.sample``
  returns ``(positions, logq)``), so recording it costs nothing and *not*
  recording it would throw away the only thing that makes an approximate run
  correctable.
* **MCLMC writes ``"none"``**, and the reason wants stating plainly rather
  than assuming. Unadjusted MCLMC is not a Metropolis chain: there is no
  accept/reject step, so the stationary distribution is the target only up to
  a discretisation bias that the tuner controls by holding the energy
  variance per dimension near ``desired_energy_var`` (blackjax's default,
  5e-4, is recorded in the attrs as ``ampere_blackjax_desired_energy_var``,
  and ``ampere_blackjax_adjusted`` records that no Metropolis correction was
  applied). ``ampere_approximation`` is not a claim about that bias: it is
  the key ``plot_trace`` and ``ampere.results.summary`` read to decide
  whether an R-hat, an ESS or a trace shape means anything, and an MCLMC run
  **is** a Markov chain, so all three do. A run whose draws are not a chain
  at all — VI's, SBI's, Pathfinder's — is what that key exists to flag. The
  bias is a separate fact, it is recorded under the engine's own names, and
  the engine battery SBC-ranks MCLMC precisely so that it is measured rather
  than asserted.

**No evidence.** Neither method estimates a marginal likelihood, so neither
writes ``ampere_log_evidence`` and a reader finds the triple absent rather
than filled with a placeholder. MCLMC is an unnormalised-density sampler like
every MCMC here; Pathfinder's ELBO is a *lower bound* on the log evidence
from a single L-BFGS path, and recording a bound of unknown tightness under
the name a nested sampler's estimate uses would be worse than recording
nothing. The ELBO is kept under Pathfinder's own name
(``ampere_blackjax_pathfinder_elbo``) where nothing will mistake it for one.

Cost accounting (design horizon (g))
-------------------------------------
:meth:`~ampere.inference.engine.Engine.finish` writes
``ampere_engine_evaluations`` for every run from the evaluation cache, so
this driver needs add nothing to satisfy the horizon. It records the
*library's* own count beside it —
``ampere_blackjax_integrator_steps`` for MCLMC (tuning included, because
tuning is gradient evaluations that were really spent) and
``ampere_blackjax_lbfgs_iterations`` for Pathfinder — for the reason
``_nested.py`` records ``ampere_nautilus_likelihood_calls``: the two numbers
are not the same when a cache serves a point twice, and the difference is
itself what a cost model wants.

jax only, and lazily
--------------------
blackjax is a jax library; there is no torch route and the memo's §2.2 says
why (no maintained torch MCLMC or Pathfinder). So this driver refuses a
problem on any other backend by name, in the constructor, and imports
blackjax there too rather than in ``run`` — :mod:`ampere.inference._zeus`'s
reasoning, unchanged: the failure belongs at the moment the user asks for the
engine, not after they have chosen a draw count and pressed go.
"""

from __future__ import annotations

import dataclasses
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import LoweringError, OptionalDependencyError
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

__all__ = ["BLACKJAX_BACKEND", "METHODS", "BlackjaxEngine"]

#: The one backend blackjax can be driven over here. A constant rather than a
#: table like ``_nuts.py``'s :data:`~ampere.inference._nuts.SAMPLER_LIBRARIES`
#: because there is nothing to tabulate: blackjax is a jax library and the
#: memo's §2.2 records that torch has no counterpart to borrow.
BLACKJAX_BACKEND = "jax"

#: What ``method=`` accepts, and what each one is. The **keys** are ampere's
#: and are what a run's ``ampere_blackjax_method`` records; the values are the
#: one-line description a refusal quotes back, so a user who guessed a name
#: learns what the real ones do rather than only that they exist.
METHODS: dict[str, str] = {
    "mclmc": (
        "microcanonical Langevin Monte Carlo: a Markov chain, tuned for you, "
        "strongest on smooth high-dimensional posteriors"
    ),
    "pathfinder": (
        "quasi-Newton variational approximation: i.i.d. draws from a Gaussian "
        "fitted along an L-BFGS path, with an importance-weightable density"
    ),
}

#: The smallest free dimension MCLMC can be run over. blackjax refuses fewer
#: ("The target distribution must have more than 1 dimension for MCLMC"), and
#: the reason is the method rather than the implementation: the isokinetic
#: dynamics decorrelate a *direction* in the momentum, and in one dimension
#: there is no direction to decorrelate. :meth:`BlackjaxEngine.run` refuses
#: first, by name, exactly as :class:`~ampere.inference.NautilusEngine` does
#: for nautilus's own two-parameter floor. Pathfinder has no such floor.
#: **W5.14.**
_MCLMC_MIN_DIM = 2

#: The smallest free dimension MCLMC can be run over. blackjax refuses fewer
#: ("The target distribution must have more than 1 dimension for MCLMC"), and
#: the reason is the method rather than the implementation: the isokinetic
#: dynamics decorrelate a *direction* in the momentum, and in one dimension
#: there is no direction to decorrelate. :meth:`BlackjaxEngine.run` refuses
#: first, by name, exactly as :class:`~ampere.inference.NautilusEngine` does
#: for nautilus's own two-parameter floor. Pathfinder has no such floor.
#: **W5.14.**
_MCLMC_MIN_DIM = 2

#: ``ampere_approximation`` (``results.md`` §9) per method. See the module
#: docstring on why MCLMC's is ``"none"`` and Pathfinder's is not.
_APPROXIMATIONS: dict[str, str] = {"mclmc": "none", "pathfinder": "pathfinder"}


@dataclasses.dataclass(frozen=True)
class _Settings:
    """One run's settings, so a method cannot be handed a partial set."""

    draws: int
    warmup: int
    chains: int
    progress: bool


def supported_backends() -> frozenset[str]:
    """The backends this driver can fit **right now**, in this interpreter.

    A function rather than a constant, for the reason
    :func:`ampere.inference._nuts.supported_backends` gives: a realisation is
    registered when the user imports the backend, which is normally after
    ``ampere.inference``.
    """
    return frozenset(registered_realisations()) & {BLACKJAX_BACKEND}


class BlackjaxEngine(Engine):
    """Fit a :class:`~ampere.core.dataset.FittingProblem` with blackjax, on jax.

    Parameters
    ----------
    problem
        The composed problem. Must be on the ``jax`` backend with a registered
        realisation, and must be ``differentiable``; both are checked at
        construction, by name.
    density
        Optional, and normally omitted — as for
        :class:`~ampere.inference.NUTSEngine`. Omitted, the density is the
        backend's registered realisation (``ampere.core.realise``).
    method
        ``"mclmc"`` (the sampler) or ``"pathfinder"`` (the approximation). See
        :data:`METHODS` and the module docstring: the choice changes what a
        run *is*, and the run says which was made in
        ``ampere_blackjax_method``.
    cache_size
        See :class:`~ampere.inference.engine.Engine`.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> import ampere.backends.jax as backend
    >>> backend.configure_x64()
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
    >>> from ampere.inference import BlackjaxEngine
    >>> grid = np.geomspace(1.0, 10.0, 24)
    >>> rng = np.random.default_rng(3)
    >>> observed = Spectrum(
    ...     grid * u.um,
    ...     (2.0 * grid ** -1.0 + rng.normal(0.0, 0.03, grid.size)) * u.Jy,
    ...     uncertainty=np.full(grid.size, 0.03) * u.Jy,
    ... )
    >>> problem = FittingProblem(
    ...     backend.PowerLaw(grid, norm=st.norm(2.0, 0.5), index=-1.0),
    ...     [Dataset(
    ...         observed,
    ...         likelihood=Likelihood(GaussianFamily(), backend.IndependentNoise()),
    ...     )],
    ...     seed=20260910,
    ... )
    >>> run = BlackjaxEngine(problem, method="pathfinder").run(draws=500)
    >>> run.attrs["ampere_engine"], run.attrs["ampere_approximation"]
    ('blackjax', 'pathfinder')
    >>> bool(abs(float(run["posterior"]["model.norm"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "blackjax"
    #: Everything here is gradient-based; that is what blackjax is for.
    OFFERS_GRADIENTS: ClassVar[bool] = True
    #: The importable name of the library, and the extra that supplies it.
    MODULE: ClassVar[str] = "blackjax"
    EXTRA: ClassVar[str] = "blackjax"

    def __init__(
        self,
        problem: FittingProblem,
        density: Any = None,
        *,
        method: str = "mclmc",
        cache_size: int = DEFAULT_CACHE_SIZE,
    ) -> None:
        # First, ahead of every other refusal: a foreign part is the most
        # specific diagnosis available and the only one that names the piece
        # (W3.8).
        _refuse_foreign_parts(self.NAME, problem)
        if method not in METHODS:
            known = "; ".join(f"{name!r} -- {what}" for name, what in sorted(METHODS.items()))
            raise EngineError(
                f"{self.NAME} does not know the method {method!r}. Available: {known}."
            )
        available = supported_backends()
        if problem.backend not in available:
            known = ", ".join(sorted(available)) or "(none: jax has not been imported here)"
            raise EngineError(
                f"{self.NAME} cannot fit a problem on the {problem.backend!r} backend. blackjax "
                f"is a jax library and there is no torch counterpart to borrow "
                f"(inference_extensions_memo.md §2.2), so this driver needs a problem on the "
                f"{BLACKJAX_BACKEND!r} backend with a registered realisation (inference.md "
                f"§10a). Satisfying that in this interpreter: {known}. Import "
                f"ampere.backends.jax and build the problem from its models, instrument steps, "
                f"noise models and GP solver; or use NUTSEngine or VIEngine, which run the torch "
                f"backend too, or a gradient-free engine, which runs every backend."
            )
        if not problem.differentiable:
            raise EngineError(
                f"{self.NAME} needs a differentiable problem, and this one declares "
                f"differentiable=False. That flag is aggregated from what the models, instrument "
                f"steps, noise models and GP solvers themselves declare, so one "
                f"non-differentiable piece is enough — problem.capabilities says which pieces "
                f"were consulted."
            )
        #: The library, imported at construction so that a missing extra is
        #: reported when the engine is asked for rather than when it is run.
        self._library = self._require()
        self.method = str(method)
        #: The backend's realisation, when this driver obtained one.
        self.realisation: Realisation | None = None
        #: The library's own final object after a run — the tuned MCLMC state
        #: and parameters, or Pathfinder's fitted state — for anything this
        #: driver does not expose.
        self.approximation: Any = None
        if density is None:
            try:
                self.realisation = realise(problem)
            except LoweringError as error:
                raise EngineError(
                    f"{self.NAME} could not obtain a differentiable density for this problem: "
                    f"{error}"
                ) from error
            # `use_realisation=False`: this driver already holds the
            # realisation and scores every proposal through it, so the base
            # class's gradient-free fast path would lower the same problem a
            # second time for the start-point search alone.
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

    @classmethod
    def _require(cls) -> Any:
        """Import the library on use, never on import (``architecture.md`` §4)."""
        import importlib

        try:
            return importlib.import_module(cls.MODULE)
        except ImportError as error:  # pragma: no cover - the minimal-install job exercises it
            raise OptionalDependencyError(
                cls.MODULE,
                extra=cls.EXTRA,
                context=f"running the blackjax engine (ampere.inference.{cls.__name__})",
            ) from error

    # -- the run --------------------------------------------------------------

    def run(
        self,
        draws: int,
        *,
        warmup: int | None = None,
        chains: int = 1,
        initial: Any = None,
        progress: bool = False,
    ) -> Any:
        """Sample or approximate, and emit the run.

        Parameters
        ----------
        draws
            Kept draws **per chain** for ``method="mclmc"`` (which needs at
            least two free parameters — see ``_MCLMC_MIN_DIM``); the number of
            i.i.d. draws taken from the fitted approximation for
            ``method="pathfinder"``, where they are cheap and a large number
            costs almost nothing beside the fit.
        warmup
            MCLMC's tuning budget, in integrator steps. Defaults to *draws*,
            which is the convention both samplers in ``_nuts.py`` use and is
            generous rather than clever. Tuning is what chooses the momentum
            decoherence length ``L``, the step size and the diagonal
            preconditioner, so too small a budget shows up as a step size that
            has not settled. Ignored by ``method="pathfinder"``, whose L-BFGS
            path has its own termination.
        chains
            Independent MCLMC chains. Several is not a luxury: R-hat and the
            split-chain diagnostics need more than one. Pathfinder produces
            one set of i.i.d. draws and refuses more than one chain rather
            than quietly ignoring the argument.
        initial
            ``(chains, n_dim)`` start positions **in the constrained space**,
            as for the other drivers; or the string ``"pathfinder"``, which
            runs Pathfinder first and starts each MCLMC chain at one of its
            draws. That is the memo's second use for it and the reason it is
            here: a quasi-Newton approximation of the posterior is a far
            better start than a prior draw on a posterior the prior barely
            covers. The default draws from the joint prior on this engine's
            own initialisation stream, so a run repeats exactly from the
            problem's seed.
        progress
            Accepted for symmetry with the other drivers and currently unused:
            both routes run inside one ``jax.lax.scan``, which has no place to
            print from. Recorded nowhere, and never on by default.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(chains, draws)`` for MCLMC
            and ``(1, draws)`` for Pathfinder.
        """
        if int(draws) < 1:
            raise EngineError(f"{self.NAME} needs at least one draw, got {draws}.")
        if int(chains) < 1:
            raise EngineError(f"{self.NAME} needs at least one chain, got {chains}.")
        if self.method == "pathfinder" and int(chains) != 1:
            raise EngineError(
                f"{self.NAME}'s 'pathfinder' method draws i.i.d. points from one fitted "
                f"approximation, so there is no such thing as a second chain of them; got "
                f"chains={chains}. Use method='mclmc' for a chain, or initial='pathfinder' to "
                f"start several MCLMC chains from this approximation."
            )
        if self.method == "mclmc" and int(self.problem.free_size) < _MCLMC_MIN_DIM:
            raise EngineError(
                f"{self.NAME}'s 'mclmc' method needs at least {_MCLMC_MIN_DIM} free parameters "
                f"and this problem has {self.problem.free_size}: the microcanonical dynamics "
                f"decorrelate a direction in the momentum, and in one dimension there is no "
                f"direction to decorrelate (blackjax refuses it too, one layer down). Use "
                f"method='pathfinder', which has no such floor, or NUTSEngine, EmceeEngine or "
                f"DynestyEngine, which all sample a one-dimensional posterior."
            )
        if self.method == "mclmc" and int(self.problem.free_size) < _MCLMC_MIN_DIM:
            raise EngineError(
                f"{self.NAME}'s 'mclmc' method needs at least {_MCLMC_MIN_DIM} free parameters "
                f"and this problem has {self.problem.free_size}: the microcanonical dynamics "
                f"decorrelate a direction in the momentum, and in one dimension there is no "
                f"direction to decorrelate (blackjax refuses it too, one layer down). Use "
                f"method='pathfinder', which has no such floor, or NUTSEngine, EmceeEngine or "
                f"DynestyEngine, which all sample a one-dimensional posterior."
            )
        adaptation = int(draws) if warmup is None else int(warmup)
        if adaptation < 1:
            raise EngineError(
                f"{self.NAME}'s tuning budget must be at least one step, got {warmup}. MCLMC has "
                f"no un-tuned default step size to fall back on."
            )
        self.start()

        settings = _Settings(
            draws=int(draws),
            warmup=adaptation,
            chains=int(chains),
            progress=bool(progress),
        )
        unconstrained = self._start_points(initial, settings)

        if self.method == "pathfinder":
            drawn, log_q, attrs = self._run_pathfinder(unconstrained[0], settings)
            sample_stats = {
                # The same coordinate move VIEngine makes, for the same
                # reason: the approximation's density is fitted in
                # unconstrained coordinates and the stored log_prior /
                # log_likelihood are in constrained ones, so §9's importance
                # weight is only buildable from the stored groups if the two
                # are brought together here (engine.py, W5.0).
                "proposal_log_density": log_q
                - unconstrained_jacobian_correction(self.problem, drawn[0])
            }
        else:
            drawn, attrs = self._run_mclmc(unconstrained, settings)
            sample_stats = None

        chain = np.stack(
            [np.stack([self.problem.constrain(y) for y in walker]) for walker in drawn]
        )
        attrs.update(
            {
                "blackjax_method": self.method,
                "blackjax_draws": settings.draws,
                "blackjax_chains": settings.chains,
                "blackjax_version": str(getattr(self._library, "__version__", "unknown")),
                "approximation": _APPROXIMATIONS[self.method],
            }
        )
        return self.finish(
            chain,
            extra_attrs=attrs,
            realised=self.realisation is not None,
            registered_lowerings=self._realisation_provenance(),
            log_likelihood_terms=self._decomposition(drawn),
            sample_stats=sample_stats,
        )

    # -- the two methods ------------------------------------------------------

    def _run_mclmc(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, dict[str, object]]:
        """blackjax's MCLMC over the jax realisation, tuned then run.

        Lazy imports, exactly as ``_nuts.py`` imports numpyro: ``import
        ampere.inference`` must not require jax (``architecture.md`` §4 rule
        2). The suppressions are the price of that rule — ``dev``, the
        environment CI typechecks in, deliberately has no jax — and the real
        check is ``pixi run -e jax typecheck``.

        Three things are worth knowing about what this method does.

        * **Tuning is not optional and is not free.** MCLMC has no default
          step size: the three-phase tuner
          (:func:`blackjax.mclmc_find_L_and_step_size`) estimates the energy
          variance, the step size and the momentum decoherence length ``L``,
          and optionally a diagonal preconditioner, over *warmup* integrator
          steps. Those steps are gradient evaluations that were really spent,
          so they are counted into ``ampere_blackjax_integrator_steps`` rather
          than quietly dropped.
        * **Each chain is tuned separately.** A chain that started somewhere
          else has a different local geometry to adapt to, and sharing one
          tuned step size across chains would make the second chain's
          behaviour depend on where the first one started. The tuned
          parameters are therefore recorded as the mean over chains, with the
          per-chain values kept on the engine.
        * **The seed is a key, not a global.** jax has no global RNG, so every
          stream here is an explicit key split from this engine's own
          ``sampler`` sub-stream, and a run repeats exactly from the problem's
          seed.
        """
        import jax  # pyrefly: ignore[missing-import]
        import blackjax  # pyrefly: ignore[missing-import]
        from blackjax.mcmc.integrators import (  # pyrefly: ignore[missing-import]
            isokinetic_mclachlan,
        )

        density = self.density
        kernel = blackjax.mcmc.mclmc.build_kernel(integrator=isokinetic_mclachlan)
        keys = jax.random.split(jax.random.key(self.integer_seed("sampler")), settings.chains)

        chains: list[np.ndarray] = []
        tuned: list[Any] = []
        tuning_steps = 0
        for index in range(settings.chains):
            init_key, tune_key, run_key = jax.random.split(keys[index], 3)
            state = blackjax.mcmc.mclmc.init(
                position=jax.numpy.asarray(unconstrained[index]),
                logdensity_fn=density,
                rng_key=init_key,
            )
            state, parameters, steps = blackjax.mclmc_find_L_and_step_size(
                mclmc_kernel=kernel,
                num_steps=settings.warmup,
                state=state,
                rng_key=tune_key,
                logdensity_fn=density,
            )
            tuning_steps += int(steps)
            algorithm = blackjax.mclmc(
                density,
                L=parameters.L,
                step_size=parameters.step_size,
                inverse_mass_matrix=parameters.inverse_mass_matrix,
            )
            _, positions = blackjax.util.run_inference_algorithm(
                rng_key=run_key,
                inference_algorithm=algorithm,
                num_steps=settings.draws,
                initial_state=state,
                transform=lambda state, info: state.position,
            )
            chains.append(np.asarray(positions, dtype=float))
            tuned.append(parameters)

        self.sampler = tuned
        attrs: dict[str, object] = {
            "blackjax_warmup": settings.warmup,
            "blackjax_L": float(np.mean([float(p.L) for p in tuned])),
            "blackjax_step_size": float(np.mean([float(p.step_size) for p in tuned])),
            "blackjax_tuning_steps": int(tuning_steps),
            "blackjax_integrator_steps": int(tuning_steps + settings.draws * settings.chains),
            "blackjax_adjusted": False,
            "blackjax_desired_energy_var": _DESIRED_ENERGY_VAR,
            "blackjax_integrator": "isokinetic_mclachlan",
            "jax_version": str(jax.__version__),
        }
        return np.stack(chains), attrs

    def _run_pathfinder(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
        """blackjax's Pathfinder: L-BFGS, then draws from the ELBO-best iterate.

        Returns the draws shaped ``(1, draws, n_dim)`` — one "chain" of i.i.d.
        points, which is the shape :class:`~ampere.inference.VIEngine` emits
        and for the same reason: there is no chain structure to split, so
        R-hat has nothing to say about them.

        The density comes back beside the draws.
        :func:`blackjax.vi.pathfinder.sample` returns ``(positions, logq)``
        despite documenting only the first, and ``logq`` is the approximation's
        own log-density at each drawn point — the quantity
        ``sample_stats.proposal_log_density`` is for. Recomputing it from the
        state's ``alpha``/``beta``/``gamma`` factors would be reimplementing
        the library's own low-rank Gaussian, so it is read rather than rebuilt,
        and the engine battery checks it against the target it was fitted to.
        """
        import jax  # pyrefly: ignore[missing-import]
        import blackjax  # pyrefly: ignore[missing-import]

        approximate_key, sample_key = jax.random.split(
            jax.random.key(self.integer_seed("optimiser")), 2
        )
        state, info = blackjax.vi.pathfinder.approximate(
            rng_key=approximate_key,
            logdensity_fn=self.density,
            initial_position=jax.numpy.asarray(unconstrained),
        )
        positions, log_q = blackjax.vi.pathfinder.sample(sample_key, state, settings.draws)
        self.approximation = state
        drawn = np.asarray(positions, dtype=float).reshape(1, settings.draws, -1)
        attrs: dict[str, object] = {
            "blackjax_pathfinder_elbo": float(np.asarray(state.elbo)),
            "blackjax_lbfgs_iterations": int(np.asarray(info.path.elbo).size),
            "jax_version": str(jax.__version__),
        }
        return drawn, np.asarray(log_q, dtype=float).reshape(settings.draws), attrs

    # -- shared plumbing ------------------------------------------------------

    def _start_points(self, initial: Any, settings: _Settings) -> np.ndarray:
        """``(chains, n_dim)`` unconstrained start points, by whichever route.

        Three routes, and the third is the memo's second use for Pathfinder:
        ``initial="pathfinder"`` fits the approximation and draws one start
        point per chain from it. The Pathfinder fit is *seeded from a prior
        draw* like everything else here, so a run started this way still
        repeats exactly from the problem's seed — the extra quality comes from
        the quasi-Newton path, not from extra randomness.
        """
        if isinstance(initial, str):
            if initial != "pathfinder":
                raise EngineError(
                    f"{self.NAME}'s only named initialiser is 'pathfinder', got {initial!r}. "
                    f"Pass an array of start points in the constrained space, or omit the "
                    f"argument to draw them from the joint prior."
                )
            if self.method == "pathfinder":
                raise EngineError(
                    f"{self.NAME} cannot initialise the 'pathfinder' method with Pathfinder: it "
                    f"would fit the same approximation twice and draw from the second. Omit "
                    f"initial=, or use method='mclmc' to sample from a Pathfinder start."
                )
            return self._pathfinder_starts(settings)
        if initial is None:
            positions = self.initial_positions(settings.chains)
        else:
            positions = self._checked_initial(initial, settings.chains)
        return np.stack([self.problem.unconstrain(theta) for theta in positions])

    def _pathfinder_starts(self, settings: _Settings) -> np.ndarray:
        """*chains* unconstrained start points drawn from a Pathfinder fit."""
        import jax  # pyrefly: ignore[missing-import]
        import blackjax  # pyrefly: ignore[missing-import]

        prior = self.problem.unconstrain(self.initial_positions(1)[0])
        approximate_key, sample_key = jax.random.split(
            jax.random.key(self.integer_seed("initialiser")), 2
        )
        state, _ = blackjax.vi.pathfinder.approximate(
            rng_key=approximate_key,
            logdensity_fn=self.density,
            initial_position=jax.numpy.asarray(prior),
        )
        positions, _ = blackjax.vi.pathfinder.sample(sample_key, state, settings.chains)
        self.approximation = state
        return np.asarray(positions, dtype=float).reshape(settings.chains, -1)

    def _realisation_provenance(self) -> list[dict[str, Any]] | None:
        """The user-registered lowering rows the realisation consulted, if any."""
        if self.realisation is None:
            return None
        rows = getattr(self.realisation, "lowering_provenance", None)
        return list(rows()) if callable(rows) else None

    def _decomposition(self, unconstrained: np.ndarray) -> list[list[dict[str, float]]] | None:
        """§10a's optional per-dataset split, for the drawn points.

        The same fetch-don't-require rule as the NUTS and VI drivers': ``None``
        when the realisation offers no decomposition, in which case
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

    def _checked_initial(self, initial: Any, chains: int) -> np.ndarray:
        positions = np.asarray(initial, dtype=float).reshape(chains, -1)
        if positions.shape != (chains, self.problem.free_size):
            raise EngineError(
                f"{self.NAME}'s start points must have shape ({chains}, "
                f"{self.problem.free_size}), got {positions.shape}."
            )
        return positions


#: blackjax's own default for the tuner's energy-variance target, recorded in
#: the attrs so that a reader of an archived MCLMC run can see the number the
#: discretisation bias was traded against. Repeated here rather than read off
#: the library because it is a *default argument* of
#: :func:`blackjax.mclmc_find_L_and_step_size`, and this driver does not pass
#: one; if that default ever changes, the conformance of this constant with it
#: is exactly what the battery's reproducibility row would notice.
_DESIRED_ENERGY_VAR = 5e-4
