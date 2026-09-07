"""The NUTS driver: numpyro's No-U-Turn sampler over a backend's lowered density.

Private module; the class is :class:`ampere.inference.NUTSEngine`. Named with a
leading underscore for the reason ``_emcee.py`` records for its own: an
``ampere/inference/nuts.py`` would be one directory from the library it
imports, and unambiguous names make that safe in practice as well as in
principle.

Why this driver takes an argument the other three do not
--------------------------------------------------------
``inference.md`` §10's claim is that an engine written against
``DEVELOPMENT_PLAN.md`` §4.5's surface "works with the reference backend, with
torch, with jax […] and never knows which", and ``EmceeEngine``,
``DynestyEngine`` and ``ZeusEngine`` are the proof: they call ``log_prob``,
``prior_transform`` and nothing else.

A gradient-based engine cannot be written that way today, and the reason is
structural rather than a gap in this driver. ``FittingProblem.log_prob`` and
``log_prob_unconstrained`` are **not traceable**: a ``ModelResult`` is built
from ``ampere.core``'s containers, which coerce their values with
``numpy.asarray``, and ``ParameterSet.lnprior`` short-circuits on
``math.isfinite`` — Python control flow on a value. Both are frozen §4
contract, both are there for good reasons, and together they mean a gradient
cannot be taken through the contract path *on any backend*.

So the gradient has to come from the backend's own lowering of the problem —
its **realisation** (``ampere.core.realisation``, the W2.13 prototype). A
backend registers a factory with ``ampere.core`` at import; this driver calls
``ampere.core.realise(problem)``, which dispatches on ``problem.backend``
(W2.12's derived flag) and checks the result against the contract path at one
point::

    import ampere.backends.jax        # registers the jax realisation
    from ampere.inference import NUTSEngine

    run = NUTSEngine(problem).run(draws=1000, warmup=1000)

That keeps the rule this namespace is built on intact — **nothing under
``ampere.inference`` imports ``ampere.backends``**, and this module does not
either; it imports jax and numpyro lazily, inside ``run``, exactly as
``_emcee.py`` imports emcee. The density can still be passed explicitly
(``NUTSEngine(problem, density)``) for a user's own lowering or a test, in which
case the same one-point agreement check runs here.

What the run records
--------------------
Everything a gradient-free run records, through the same :meth:`Engine.finish`:
the per-draw ``log_prior``/``log_likelihood`` split, the per-dataset
decomposition, the observed data and the full provenance attrs, with
``ampere_backend`` read off the problem. Plus NUTS's own diagnostics —
divergences, tree depth, step size and the accept probability — in the attrs.

The per-draw evaluations are **recomputed** rather than cached, and the run
says so in ``engine_draws_recomputed``. The other drivers score every proposal
through :class:`~ampere.inference.engine.Engine`'s evaluation cache, so a
stored draw is usually a lookup; this one scores through the *lowered* density,
which returns a scalar rather than an
:class:`~ampere.core.dataset.Evaluation`, so the decomposition is computed once
per stored draw on the contract path afterwards. That is honest and it is
counted; making the lowered path return a full decomposition is slice 2's, with
the rest of the native likelihood work.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.realisation import realise

from .engine import DEFAULT_CACHE_SIZE, Engine, _kept
from ampere.core.exceptions import LoweringError

from .exceptions import EngineError

__all__ = ["NUTSEngine"]

#: The backends whose lowering this driver knows how to drive. A string, not an
#: import: the driver reads ``problem.backend`` (W2.12's derived capability
#: flag) and compares, so adding torch here is a one-line change once that
#: backend ships a traceable density — and getting it wrong is a refusal by
#: name rather than a crash inside a trace.
SUPPORTED_BACKENDS: frozenset[str] = frozenset({"jax"})


class NUTSEngine(Engine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with numpyro's NUTS.

    Hoffman & Gelman's No-U-Turn sampler over the *unconstrained* density, with
    numpyro's window adaptation for the step size and mass matrix. The engine
    of choice once a problem is differentiable: it scales to far more
    dimensions than an affine-invariant ensemble, and its divergences are a
    diagnostic the gradient-free engines cannot offer.

    Parameters
    ----------
    problem
        The composed problem. Must report a backend in
        :data:`SUPPORTED_BACKENDS` and must be ``differentiable``; both are
        checked at construction, by name, rather than discovered at the first
        gradient.
    density
        Optional. The backend's lowering of *problem*: a **pure jax function**
        of the unconstrained free vector returning the log density, change of
        variables included. **Omitted, it is obtained through
        ``ampere.core.realise(problem)``** — the backend's registered
        realisation (W2.13 prototype), which importing ``ampere.backends.jax``
        registers; this driver still imports no backend. Passing one
        explicitly remains possible (a user's own lowering, or a test), in
        which case it is checked against the contract path at one point here.
    cache_size
        See :class:`~ampere.inference.engine.Engine`. Note that this driver's
        cache is only ever populated by its start-point search, since the
        sampler itself scores through *density*; the stored draws are
        decomposed afterwards on the contract path.

    Notes
    -----
    **The density must agree with the problem.** Nothing here can check that in
    general — they are two implementations of one quantity — so the constructor
    checks the cheapest necessary condition: that *density* and
    ``problem.log_prob_unconstrained`` agree at one point. A backend whose
    lowering had drifted from its contract path would otherwise sample a
    different posterior from the one the run records, and the run would look
    perfectly healthy. ``tests/backends/test_jax.py`` makes the full comparison.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.jax import PowerLaw, configure_x64
    >>> from ampere.core import Dataset, FittingProblem, Spectrum
    >>> from ampere.inference import NUTSEngine
    >>> configure_x64()
    >>> grid = np.geomspace(1.0, 10.0, 12)
    >>> truth = 2.0 * grid ** -1.0
    >>> rng = np.random.default_rng(3)
    >>> observed = Spectrum(
    ...     grid * u.um,
    ...     (truth + rng.normal(0.0, 0.05, grid.size)) * u.Jy,
    ...     uncertainty=np.full(grid.size, 0.05) * u.Jy,
    ... )
    >>> problem = FittingProblem(
    ...     PowerLaw(grid, norm=st.norm(2.0, 0.5), index=-1.0), [Dataset(observed)],
    ...     seed=20260907,
    ... )
    >>> run = NUTSEngine(problem).run(draws=300, warmup=300, chains=2)
    >>> run["posterior"]["model.norm"].shape
    (2, 300)
    >>> run.attrs["ampere_engine"], run.attrs["ampere_backend"]
    ('nuts', 'jax')
    >>> bool(abs(float(run["posterior"]["model.norm"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "nuts"
    #: The whole point of this driver, and what ``check_engine`` is told.
    OFFERS_GRADIENTS: ClassVar[bool] = True

    def __init__(
        self,
        problem: FittingProblem,
        density: Callable[[Any], Any] | None = None,
        *,
        cache_size: int = DEFAULT_CACHE_SIZE,
    ) -> None:
        if problem.backend not in SUPPORTED_BACKENDS:
            raise EngineError(
                f"{self.NAME} cannot sample a problem on the {problem.backend!r} backend: this "
                f"driver runs numpyro's sampler, which needs a jax-traceable log-density, and "
                f"only {sorted(SUPPORTED_BACKENDS)} supply one. Build the model and instrument "
                f"steps from ampere.backends.jax, or use a gradient-free engine (emcee, dynesty, "
                f"zeus), which run every backend."
            )
        if not problem.differentiable:
            raise EngineError(
                f"{self.NAME} needs a differentiable problem, and this one declares "
                f"differentiable=False. That flag is aggregated from what the models and "
                f"instrument steps themselves declare (W2.12), so one non-differentiable step is "
                f"enough to make the whole chain non-differentiable — "
                f"problem.capabilities says which pieces were consulted."
            )
        if density is None:
            # W2.13 prototype: the backend's own realisation, reached through
            # ampere.core rather than by importing the backend -- ``realise``
            # dispatches on problem.backend and has already checked agreement
            # with the contract path at the reference point.
            try:
                realised = realise(problem)
            except LoweringError as error:
                raise EngineError(
                    f"{self.NAME} could not obtain a differentiable density for this problem: "
                    f"{error}"
                ) from error
            density = realised.log_prob_unconstrained
            super().__init__(problem, cache_size=cache_size)
            self.density = density
            return
        if not callable(density):
            raise EngineError(
                f"{self.NAME} takes the backend's lowering of the problem as a callable of the "
                f"unconstrained vector, got {density!r}. Omit it to use the backend's registered "
                f"realisation (ampere.core.realise), or on the jax backend pass "
                f"ampere.backends.jax.lower_problem(problem).log_prob_unconstrained."
            )
        super().__init__(problem, cache_size=cache_size)
        self.density = density
        self._check_density_agrees()

    def _check_density_agrees(self) -> None:
        """The cheapest necessary condition that *density* is this problem's.

        Not a proof — they are two implementations of one quantity, and only
        the conformance and backend suites compare them properly. But a driver
        handed the *wrong* problem's density would sample a healthy-looking
        posterior of something else entirely, and one point costs nothing.
        """
        reference = self.problem.unconstrain(self.problem.reference_values)
        expected = self.problem.log_prob_unconstrained(reference)
        try:
            got = float(np.asarray(self.density(reference)))
        except Exception as error:  # the lowering refusing is not this driver's to interpret
            raise EngineError(
                f"{self.NAME} could not evaluate the supplied density at the problem's reference "
                f"values: {error}. It must be a callable of the unconstrained free vector "
                f"(length {self.problem.free_size}) returning a scalar."
            ) from error
        if not math.isclose(got, expected, rel_tol=1e-6, abs_tol=1e-6):
            raise EngineError(
                f"{self.NAME} was given a density that disagrees with the problem it was given: "
                f"at the problem's reference values the lowered density is {got!r} and "
                f"FittingProblem.log_prob_unconstrained is {expected!r}. The two are meant to be "
                f"the same quantity computed twice; sampling the first and recording the second "
                f"would produce a run that looks healthy and describes a different posterior."
            )

    def run(
        self,
        draws: int,
        *,
        warmup: int | None = None,
        chains: int = 1,
        max_tree_depth: int = 10,
        target_accept_prob: float = 0.8,
        dense_mass: bool = False,
        initial: Any = None,
        progress: bool = False,
    ) -> Any:
        """Adapt, sample, and emit the run.

        Parameters
        ----------
        draws
            Kept draws **per chain**, after warm-up. NUTS keeps everything it
            draws past adaptation, so there is no burn-in argument: *warmup* is
            adaptation and is discarded by construction.
        warmup
            Adaptation steps per chain. Defaults to *draws*, which is
            numpyro's own convention and is generous rather than clever.
        chains
            Independent chains. Several is not a luxury here — R-hat and the
            split-chain diagnostics ArviZ computes need more than one, and NUTS
            chains are cheap enough to make that the default advice.
        max_tree_depth
            NUTS's doubling limit. Hitting it is a diagnostic (the trajectory
            wanted to be longer than 2**depth), recorded in the attrs.
        target_accept_prob
            Window adaptation's target. Raise it towards 0.95 for a posterior
            with strong curvature, at the cost of smaller steps.
        dense_mass
            Adapt a full mass matrix rather than a diagonal one. Worth it for
            a strongly correlated posterior of modest dimension.
        initial
            ``(chains, n_dim)`` start positions **in the constrained space**,
            as for the other drivers. The default draws them from the joint
            prior on this engine's own initialisation stream, so a run repeats
            exactly from the problem's seed.
        progress
            numpyro's progress bar. Off by default: a driver that prints by
            default is unusable inside a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(chains, draws)``.
        """
        # Lazy, exactly as _emcee.py imports emcee: `import ampere.inference`
        # must not require jax (architecture.md §4 rule 2). The suppressions are
        # the price of that rule: `dev` -- the environment CI typechecks in --
        # deliberately has no jax, so pyrefly cannot resolve these and would
        # report three unresolvable imports rather than three type errors. The
        # real check happens in `pixi run -e jax typecheck`, which resolves them.
        import jax  # pyrefly: ignore[missing-import]
        import numpyro  # pyrefly: ignore[missing-import]
        from numpyro.infer import MCMC, NUTS  # pyrefly: ignore[missing-import]

        if int(chains) < 1:
            raise EngineError(f"{self.NAME} needs at least one chain, got {chains}.")
        adaptation = int(draws) if warmup is None else int(warmup)
        if adaptation < 0:
            raise EngineError(f"{self.NAME}'s warm-up cannot be negative, got {warmup}.")
        # `_kept` is the shared refusal of a burn-in/thin combination that keeps
        # nothing; NUTS thins nothing and discards only its adaptation, so the
        # only condition left is that it keeps at least one draw.
        _kept(int(draws), 0, 1, self.NAME)
        self.start()

        positions = (
            self.initial_positions(int(chains))
            if initial is None
            else self._checked_initial(initial, int(chains))
        )
        unconstrained = np.stack([self.problem.unconstrain(theta) for theta in positions])

        def potential(y: Any) -> Any:
            """numpyro minimises an energy; the density is maximised.

            One line, spelled here rather than inferred, because getting the
            sign wrong gives a sampler that explores the prior's tails with
            perfect efficiency and reports no divergences at all.
            """
            return -self.density(y)

        kernel = NUTS(
            potential_fn=potential,
            max_tree_depth=int(max_tree_depth),
            target_accept_prob=float(target_accept_prob),
            dense_mass=bool(dense_mass),
        )
        mcmc = MCMC(
            kernel,
            num_warmup=adaptation,
            num_samples=int(draws),
            num_chains=int(chains),
            chain_method="sequential",
            progress_bar=bool(progress),
        )
        rng = jax.random.key(self.integer_seed("sampler"))
        mcmc.run(
            rng,
            init_params=jax.numpy.asarray(unconstrained),
            extra_fields=("diverging", "num_steps", "adapt_state.step_size", "accept_prob"),
        )
        self.sampler = mcmc

        drawn = np.asarray(mcmc.get_samples(group_by_chain=True), dtype=float)
        if drawn.ndim == 2:  # a single chain comes back flat
            drawn = drawn[np.newaxis, ...]
        # Back to the constrained space the results schema and every other
        # engine record: `constrain` is the problem's own, so the posterior a
        # NUTS run stores is directly comparable with an emcee one.
        chain = np.stack(
            [np.stack([self.problem.constrain(y) for y in walker]) for walker in drawn]
        )

        extra = mcmc.get_extra_fields(group_by_chain=True)
        divergences = int(np.sum(np.asarray(extra["diverging"]))) if "diverging" in extra else 0
        attrs: dict[str, object] = {
            "nuts_draws": int(draws),
            "nuts_warmup": adaptation,
            "nuts_chains": int(chains),
            "nuts_max_tree_depth": int(max_tree_depth),
            "nuts_target_accept_prob": float(target_accept_prob),
            "nuts_dense_mass": bool(dense_mass),
            "nuts_divergences": divergences,
            "numpyro_version": numpyro.__version__,
            "jax_version": jax.__version__,
        }
        if "num_steps" in extra:
            attrs["nuts_mean_tree_size"] = float(np.mean(np.asarray(extra["num_steps"])))
        if "accept_prob" in extra:
            attrs["nuts_mean_accept_prob"] = float(np.mean(np.asarray(extra["accept_prob"])))
        if "adapt_state.step_size" in extra:
            attrs["nuts_step_size"] = float(
                np.mean(np.asarray(extra["adapt_state.step_size"])[..., -1])
            )
        return self.finish(chain, extra_attrs=attrs)

    def _checked_initial(self, initial: Any, chains: int) -> np.ndarray:
        positions = np.asarray(initial, dtype=float)
        expected = (chains, self.problem.free_size)
        if positions.shape != expected:
            raise EngineError(
                f"{self.NAME} start positions must be shaped {expected} (chains, n_dim), got "
                f"{positions.shape}."
            )
        return positions
