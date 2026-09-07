"""The NUTS driver: a No-U-Turn sampler over a backend's realisation.

Private module; the class is :class:`ampere.inference.NUTSEngine`. Named with a
leading underscore for the reason ``_emcee.py`` records for its own: an
``ampere/inference/nuts.py`` would be one directory from the library it
imports, and unambiguous names make that safe in practice as well as in
principle.

Why a gradient-based driver needs something the other three do not
------------------------------------------------------------------
``inference.md`` §10's claim is that an engine written against
``DEVELOPMENT_PLAN.md`` §4.5's surface "works with the reference backend, with
torch, with jax […] and never knows which", and ``EmceeEngine``,
``DynestyEngine`` and ``ZeusEngine`` are the proof: they call ``log_prob``,
``prior_transform`` and nothing else.

A gradient-based engine cannot be written that way, and the reason is
structural rather than a gap in this driver. ``FittingProblem.log_prob`` and
``log_prob_unconstrained`` are **not traceable**: a ``ModelResult`` is built
from ``ampere.core``'s containers, which coerce their values with
``numpy.asarray``, and ``ParameterSet.lnprior`` short-circuits on
``math.isfinite`` — Python control flow on a value. Both are frozen §4
contract, both are there for good reasons, and together they mean a gradient
cannot be taken through the contract path *on any backend*. Both Phase 2
backend tracks found it independently.

So the gradient comes from the backend's own lowering of the problem — its
**realisation** (``inference.md`` §10a, ``ampere.core.realisation``). A
backend registers a factory with ``ampere.core`` at import; this driver calls
``ampere.core.realise(problem)``, which dispatches on ``problem.backend``
(W2.12's derived flag) and checks the result against the contract path at one
point::

    import ampere.backends.jax        # registers the jax realisation
    from ampere.inference import NUTSEngine

    run = NUTSEngine(problem).run(draws=1000, warmup=1000)

and the same three lines with ``ampere.backends.torch`` sample the torch
realisation instead. That keeps the rule this namespace is built on intact —
**nothing under ``ampere.inference`` imports ``ampere.backends``**, and this
module does not either.

Two samplers, one driver
------------------------
The *density* is backend-neutral by the time this driver has it, but the
**sampler** is not: numpyro's NUTS wants a jax function and pyro's wants a
torch one, and neither will consume the other's array. So :meth:`NUTSEngine.run`
dispatches on ``problem.backend`` between the two routes, and imports the
library it needs **lazily, inside ``run``**, exactly as ``_emcee.py`` imports
emcee. Both routes are the same three steps — build the kernel over a
*potential function* (the negated log density), start it from this engine's own
seeded draws from the prior, run chains sequentially — because the potential
form is the one that takes an arbitrary callable rather than a probabilistic
program, and a realisation is a callable.

:data:`SAMPLER_LIBRARIES` is that table. What a backend must do to be sampled
here is register a realisation *and* appear in it;
:func:`supported_backends` reports the intersection, and it is what the
constructor checks — so "which backends can NUTS sample?" is answered from the
registry at the moment it is asked rather than from a constant written when
only jax existed.

What the run records
--------------------
Everything a gradient-free run records, through the same :meth:`Engine.finish`:
the per-draw ``log_prior``/``log_likelihood`` split, the per-dataset
decomposition, the observed data and the full provenance attrs, with
``ampere_backend`` read off the problem. Plus NUTS's own diagnostics —
divergences, tree depth, step size and the accept probability — and, since
W2.13, ``ampere_realised = 1`` and the user-registered lowering rows the
realisation consulted (``inference.md`` §10a, "Provenance").

The per-draw evaluations come from the realisation itself, and the run records
``engine_draws_recomputed = 0``. W2.13 shipped this driver recomputing them:
the other drivers score every proposal through
:class:`~ampere.inference.engine.Engine`'s evaluation cache, so a stored draw
is usually a lookup, while this one scores through the *realised* density,
which returns a scalar rather than an
:class:`~ampere.core.dataset.Evaluation` — so every stored draw was
re-decomposed on the numpy contract path afterwards. That was honest and it
was counted, and it cost one full model evaluation per draw on top of the ones
the sampler had already paid for. §10a's optional ``log_likelihood_terms`` is
what removes it — both shipped realisations supply it — and consuming it meant
teaching ``Engine.finish`` to accept a decomposition it did not compute, which
is a change to the emission path rather than to this driver. Both halves
landed in the slice-2 pass: :meth:`NUTSEngine._decomposition` here, and
``Engine.evaluations_from_terms`` there.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Callable
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

from .engine import DEFAULT_CACHE_SIZE, Engine, _kept
from .exceptions import EngineError

__all__ = ["SAMPLER_LIBRARIES", "NUTSEngine", "supported_backends"]

#: Which sampler library drives which backend's realisation. Strings, not
#: imports: this module must import no backend and no sampler at module level
#: (``architecture.md`` §4 rule 2), so the mapping is names and ``run``
#: imports the one it needs.
#:
#: The two entries are not interchangeable and the table is not indirection
#: for its own sake — numpyro's kernel wants a jax function and pyro's wants a
#: torch one, and a realisation's density is in its own backend's array type
#: by construction (``inference.md`` §10a). A third differentiable backend
#: joins by adding a row here *and* a route in :meth:`NUTSEngine.run`.
SAMPLER_LIBRARIES: dict[str, str] = {"jax": "numpyro", "torch": "pyro"}

#: The one site name pyro's potential-function interface is given.
#:
#: pyro's ``NUTS(potential_fn=...)`` takes a callable of a *dict* of named
#: tensors — the shape a probabilistic program would have — while a
#: realisation is a function of one flat unconstrained vector by construction
#: (``inference.md`` §10a). So the whole free vector travels under one name.
#: Splitting it into per-parameter sites would mean re-deriving ampere's own
#: routing inside the sampler, which is exactly what ``ParameterMapping``
#: exists to prevent two copies of.
_PYRO_SITE = "theta"


@dataclasses.dataclass(frozen=True)
class _Settings:
    """One run's sampler settings, so the two routes take the same arguments.

    Not part of the public surface: it exists so that
    :meth:`NUTSEngine._run_numpyro` and :meth:`NUTSEngine._run_pyro` cannot
    drift apart in what they are given, which is the failure mode a pair of
    long parameter lists invites.
    """

    draws: int
    warmup: int
    chains: int
    max_tree_depth: int
    target_accept_prob: float
    dense_mass: bool
    progress: bool


def supported_backends() -> frozenset[str]:
    """The backends this driver can sample **right now**, in this interpreter.

    A backend qualifies when it has registered a realisation — so
    ``ampere.core.realise`` can build a differentiable density for it — *and*
    this driver knows which sampler drives it (:data:`SAMPLER_LIBRARIES`).

    A function rather than the module constant it replaces (W2.13), because
    the first half of that condition is **not** knowable when this module is
    imported: a realisation is registered when the user imports the backend,
    which is normally after ``ampere.inference``. A frozen set would either
    have to name backends that are not importable here or go stale the moment
    a third one ships. Asking the registry answers the question at the moment
    it is asked, and lets the refusal list remedies that are actually true of
    the interpreter it is raised in.
    """
    return frozenset(registered_realisations()) & frozenset(SAMPLER_LIBRARIES)


class NUTSEngine(Engine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with NUTS.

    Hoffman & Gelman's No-U-Turn sampler over the *unconstrained* density, with
    window adaptation for the step size and mass matrix — numpyro's on a jax
    problem, pyro's on a torch one. The engine of choice once a problem is
    differentiable: it scales to far more dimensions than an affine-invariant
    ensemble, and its divergences are a diagnostic the gradient-free engines
    cannot offer.

    Parameters
    ----------
    problem
        The composed problem. Must report a backend in
        :func:`supported_backends` and must be ``differentiable``; both are
        checked at construction, by name, rather than discovered at the first
        gradient.
    density
        Optional, and normally omitted. The backend's lowering of *problem*: a
        callable of the unconstrained free vector returning the log density in
        that backend's own array type, change of variables included.
        **Omitted, it comes from ``ampere.core.realise(problem)``** — the
        backend's registered realisation (``inference.md`` §10a), which
        importing ``ampere.backends.jax`` or ``ampere.backends.torch``
        registers; this driver still imports no backend. Passing one
        explicitly remains possible (a user's own lowering, or a test), in
        which case it is checked against the contract path at one point here
        and the run records ``ampere_realised = 0``, because it was not the
        registered realisation that produced the draws.
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
    perfectly healthy. ``tests/conformance`` makes the full comparison, once
    per registered realisation.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.jax import IndependentNoise, PowerLaw, configure_x64
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
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
    ...     PowerLaw(grid, norm=st.norm(2.0, 0.5), index=-1.0),
    ...     [Dataset(observed, likelihood=Likelihood(GaussianFamily(), IndependentNoise()))],
    ...     seed=20260907,
    ... )
    >>> run = NUTSEngine(problem).run(draws=300, warmup=300, chains=2)
    >>> run["posterior"]["model.norm"].shape
    (2, 300)
    >>> run.attrs["ampere_engine"], run.attrs["ampere_backend"]
    ('nuts', 'jax')
    >>> run.attrs["ampere_realised"]
    1
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
        available = supported_backends()
        if problem.backend not in available:
            known = ", ".join(sorted(available)) or "(none: no backend has been imported here)"
            raise EngineError(
                f"{self.NAME} cannot sample a problem on the {problem.backend!r} backend. This "
                f"driver needs two things: a registered realisation, so there is a "
                f"differentiable density at all (inference.md §10a), and a sampler that can "
                f"consume it — {SAMPLER_LIBRARIES}. Backends satisfying both in this "
                f"interpreter: {known}. Import the backend that supplies one (importing "
                f"ampere.backends.jax or ampere.backends.torch registers its realisation) and "
                f"build the problem from its models, instrument steps, noise models and GP "
                f"solver; or use a gradient-free engine (emcee, dynesty, zeus), which run every "
                f"backend on the contract path."
            )
        if not problem.differentiable:
            raise EngineError(
                f"{self.NAME} needs a differentiable problem, and this one declares "
                f"differentiable=False. That flag is aggregated from what the models, instrument "
                f"steps, noise models and GP solvers themselves declare (W2.12, widened at "
                f"W2.13), so one non-differentiable piece is enough to make the whole chain "
                f"non-differentiable — problem.capabilities says which pieces were consulted."
            )
        #: The backend's realisation, when this driver obtained one. ``None``
        #: when the caller supplied a density, which is what makes
        #: ``ampere_realised`` a fact about the run rather than about the
        #: backend.
        self.realisation: Realisation | None = None
        if density is None:
            # The backend's own realisation, reached through ampere.core
            # rather than by importing the backend: ``realise`` dispatches on
            # problem.backend and has already checked agreement with the
            # contract path at the reference point.
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
                f"realisation (ampere.core.realise), or pass this backend's own — e.g. "
                f"ampere.backends.jax.lower_problem(problem).log_prob_unconstrained."
            )
        super().__init__(problem, cache_size=cache_size)
        self.density = density
        self._check_density_agrees()

    def _realisation_provenance(self) -> list[dict[str, Any]] | None:
        """The user-registered lowering rows the realisation consulted, if any.

        ``lowering.md`` §12.8's stamping, made first-class at W2.13. Optional
        on the realisation surface — §10a's mandatory three do not include it
        — so it is fetched rather than required, and a realisation that does
        not offer it simply records no rows.
        """
        if self.realisation is None:
            return None
        rows = getattr(self.realisation, "lowering_provenance", None)
        return list(rows()) if callable(rows) else None

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
            Adaptation steps per chain. Defaults to *draws*, which is both
            samplers' own convention and is generous rather than clever.
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
            The sampler's progress bar. Off by default: a driver that prints
            by default is unusable inside a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(chains, draws)``.
        """
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

        settings = _Settings(
            draws=int(draws),
            warmup=adaptation,
            chains=int(chains),
            max_tree_depth=int(max_tree_depth),
            target_accept_prob=float(target_accept_prob),
            dense_mass=bool(dense_mass),
            progress=bool(progress),
        )
        library = SAMPLER_LIBRARIES[self.problem.backend]
        if library == "numpyro":
            drawn, attrs = self._run_numpyro(unconstrained, settings)
        else:
            drawn, attrs = self._run_pyro(unconstrained, settings)

        # Back to the constrained space the results schema and every other
        # engine record: `constrain` is the problem's own, so the posterior a
        # NUTS run stores is directly comparable with an emcee one.
        chain = np.stack(
            [np.stack([self.problem.constrain(y) for y in walker]) for walker in drawn]
        )
        attrs.update(
            {
                "nuts_draws": settings.draws,
                "nuts_warmup": settings.warmup,
                "nuts_chains": settings.chains,
                "nuts_max_tree_depth": settings.max_tree_depth,
                "nuts_target_accept_prob": settings.target_accept_prob,
                "nuts_dense_mass": settings.dense_mass,
                "nuts_sampler": library,
            }
        )
        return self.finish(
            chain,
            extra_attrs=attrs,
            realised=self.realisation is not None,
            registered_lowerings=self._realisation_provenance(),
            log_likelihood_terms=self._decomposition(drawn),
        )

    def _decomposition(self, unconstrained: np.ndarray) -> dict[str, np.ndarray] | None:
        """The per-dataset log-likelihood of every stored draw, natively.

        ``inference.md`` §10a's optional ``log_likelihood_terms``, used
        "when it is there" (sub-decision 1). Until this method existed, every
        stored draw of a NUTS run was re-decomposed on the numpy contract path
        by :meth:`~ampere.inference.engine.Engine.finish` — one full model
        evaluation per draw, on top of the ones the sampler had already paid
        for, and for a stochastic model not even the same number the sampler
        accepted on. W2.13 recorded that and left it, because consuming the
        member is a change to the emission path rather than to this driver;
        the emission path learned it in the slice-2 pass and this is the
        driver half.

        The loop is a Python loop over draws rather than a ``vmap``, and
        deliberately: ``ampere.inference`` may import no backend
        (``architecture.md`` §4 rule 2), so it cannot reach for jax's or
        torch's batching, and the realisations' own batched surfaces are not
        part of §10a's contract. Each call is one traced density evaluation of
        an already-compiled function — cheap beside the model evaluation it
        replaces. A realisation that supplies no decomposition returns
        ``None`` here and the cache-and-recompute path stands.

        Parameters
        ----------
        unconstrained
            ``(chains, draws, n_dim)`` in the **unconstrained** space, which
            is the argument ``log_likelihood_terms`` takes.
        """
        if self.realisation is None:
            return None
        terms = log_likelihood_terms_of(self.realisation)
        if terms is None:
            return None
        chains, draws = unconstrained.shape[0], unconstrained.shape[1]
        collected: dict[str, np.ndarray] = {}
        for chain in range(chains):
            for index in range(draws):
                found = terms(unconstrained[chain, index])
                for label, value in found.items():
                    column = collected.get(label)
                    if column is None:
                        column = np.empty((chains, draws), dtype=float)
                        collected[label] = column
                    column[chain, index] = float(np.asarray(value))
        return collected

    # -- the two sampler routes ----------------------------------------------

    def _run_numpyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, dict[str, object]]:
        """numpyro's NUTS over a jax realisation.

        Lazy imports, exactly as ``_emcee.py`` imports emcee: ``import
        ampere.inference`` must not require jax (``architecture.md`` §4 rule
        2). The suppressions are the price of that rule — ``dev``, the
        environment CI typechecks in, deliberately has no jax, so pyrefly
        cannot resolve these and would report unresolvable imports rather than
        type errors. The real check is ``pixi run -e jax typecheck``.
        """
        import jax  # pyrefly: ignore[missing-import]
        import numpyro  # pyrefly: ignore[missing-import]
        from numpyro.infer import MCMC, NUTS  # pyrefly: ignore[missing-import]

        kernel = NUTS(
            potential_fn=self._potential,
            max_tree_depth=settings.max_tree_depth,
            target_accept_prob=settings.target_accept_prob,
            dense_mass=settings.dense_mass,
        )
        mcmc = MCMC(
            kernel,
            num_warmup=settings.warmup,
            num_samples=settings.draws,
            num_chains=settings.chains,
            chain_method="sequential",
            progress_bar=settings.progress,
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

        extra = mcmc.get_extra_fields(group_by_chain=True)
        attrs: dict[str, object] = {
            "nuts_divergences": (
                int(np.sum(np.asarray(extra["diverging"]))) if "diverging" in extra else 0
            ),
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
        return drawn, attrs

    def _run_pyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, dict[str, object]]:
        """pyro's NUTS over a torch realisation.

        The same three steps as the numpyro route, in pyro's spelling.

        * ``pyro.infer.NUTS(potential_fn=...)`` takes a callable of a **dict**
          of named tensors rather than of a bare vector, and ``MCMC`` takes
          ``initial_params`` in the same shape. Ampere's realisation is a
          function of one flat unconstrained vector, so the whole problem
          travels under a single site name, :data:`_PYRO_SITE`, and the
          adaptation sees exactly the vector every other driver reports. That
          is not a workaround: a realisation deliberately has no site
          structure, because ampere's declaration already routes names and
          the sampler is not the place to re-derive that.
        * **chains are run sequentially, here, rather than by pyro.** pyro's
          own ``num_chains > 1`` forks worker processes through
          ``torch.multiprocessing``, which needs the potential to be
          picklable — it is a closure over a lowered problem, so it is not —
          and would give each worker its own RNG state to reconcile.
          numpyro's route asks for ``chain_method="sequential"`` for the
          comparable reason, so the two agree in behaviour as well as in
          shape, and a chain here starts from this engine's own seeded draw
          from the prior.
        * **the seed is forked, not simply set.** pyro's kernel draws its
          momenta from torch's *global* generator and offers no ``generator=``,
          so ``lowering.md`` §9.1's route (1) — seed the global stream from
          this engine's own ``sampler`` sub-stream — is the only one available,
          and it is what makes a run repeat from the problem's seed. It is done
          inside ``torch.random.fork_rng``, so the host process's RNG state is
          restored afterwards: this backend's rule is that ampere does not
          leave global torch state changed behind it (``lowering.md`` §10.1's
          objection to ``set_default_dtype``, applied to the other piece of
          global state a library can disturb).
        """
        import pyro  # pyrefly: ignore[missing-import]
        import torch  # pyrefly: ignore[missing-import]
        from pyro.infer import MCMC, NUTS  # pyrefly: ignore[missing-import]

        def potential(params: dict[str, Any]) -> Any:
            return -self.density(params[_PYRO_SITE])

        chains: list[np.ndarray] = []
        divergences = 0
        accepted: list[float] = []
        step_sizes: list[float] = []
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(self.integer_seed("sampler"))
            for index in range(settings.chains):
                start = torch.as_tensor(unconstrained[index], dtype=torch.float64)
                kernel = NUTS(
                    potential_fn=potential,
                    max_tree_depth=settings.max_tree_depth,
                    target_accept_prob=settings.target_accept_prob,
                    full_mass=settings.dense_mass,
                )
                mcmc = MCMC(
                    kernel,
                    num_samples=settings.draws,
                    warmup_steps=settings.warmup,
                    initial_params={_PYRO_SITE: start},
                    num_chains=1,
                    disable_progbar=not settings.progress,
                )
                mcmc.run()
                self.sampler = mcmc
                chains.append(
                    np.asarray(mcmc.get_samples()[_PYRO_SITE].detach().cpu().numpy(), dtype=float)
                )
                diagnostics = mcmc.diagnostics()
                divergences += sum(
                    len(found) for found in dict(diagnostics.get("divergences", {})).values()
                )
                accepted.extend(
                    float(value) for value in dict(diagnostics.get("acceptance rate", {})).values()
                )
                step_sizes.append(float(kernel.step_size))

        attrs: dict[str, object] = {
            "nuts_divergences": int(divergences),
            "pyro_version": pyro.__version__,
            "torch_version": torch.__version__,
        }
        if accepted:
            attrs["nuts_mean_accept_prob"] = float(np.mean(accepted))
        if step_sizes:
            attrs["nuts_step_size"] = float(np.mean(step_sizes))
        return np.stack(chains), attrs

    def _potential(self, y: Any) -> Any:
        """The **negated** density: a sampler minimises an energy.

        One line, spelled here rather than inferred at each call site, because
        getting the sign wrong gives a sampler that explores the prior's tails
        with perfect efficiency and reports no divergences at all.
        """
        return -self.density(y)

    def _checked_initial(self, initial: Any, chains: int) -> np.ndarray:
        positions = np.asarray(initial, dtype=float)
        expected = (chains, self.problem.free_size)
        if positions.shape != expected:
            raise EngineError(
                f"{self.NAME} start positions must be shaped {expected} (chains, n_dim), got "
                f"{positions.shape}."
            )
        return positions
