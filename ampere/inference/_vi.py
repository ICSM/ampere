"""The VI driver: stochastic variational inference over a backend's realisation.

Private module; the class is :class:`ampere.inference.VIEngine`. Named with a
leading underscore for the reason ``_emcee.py`` records for its own.

Why a VI driver belongs here at all
-----------------------------------
``DEVELOPMENT_PLAN.md`` §3 lists gradient-based VI beside NUTS as what
differentiability *buys*: NUTS gives the right answer slowly, VI gives an
approximate one fast, and the choice between them is a scientific judgement
about how much a fit is worth. A 10⁵-point spectrum under the flexible
likelihood is exactly the case where the judgement is live.

It is written against the **same** surface :class:`~ampere.inference._nuts.
NUTSEngine` uses — the backend's registered *realisation* (``inference.md``
§10a), reached through :func:`ampere.core.realise` — for the same structural
reason: ``FittingProblem.log_prob_unconstrained`` is not traceable, so a
gradient cannot come from §4.5's surface on any backend. And it keeps the
same rule: **nothing here imports ``ampere.backends``**, numpyro is imported
lazily inside :meth:`VIEngine.run`, and the driver dispatches on the problem's
own ``backend`` flag.

Guides, and what "the posterior" means afterwards
-------------------------------------------------
The guide is an **autoguide** over the whole unconstrained vector:
``AutoNormal`` (mean-field: one independent Normal per dimension) or
``AutoMultivariateNormal`` (a full-covariance Normal, so linear correlations
between parameters survive). Both are numpyro's, both are fitted by maximising
the ELBO, and neither is the posterior — they are a *family* of
approximations, and the fit finds the closest member of it.

That distinction is why this driver records what it does. A VI run's
``posterior`` group holds draws from the fitted guide, not from the posterior,
and the only honest way to say so in an archived file is to record the guide
family, the number of optimisation steps and the ELBO trace, so a reader can
see whether the optimisation converged and what family it converged within.
All three go into the run's provenance attrs (``vi_guide``, ``vi_steps``,
``vi_elbo_trace`` and the rest), and ``ampere_engine`` is ``"vi"``, so no run
can be mistaken for an MCMC one.

**Mean-field VI understates variance**, systematically and by construction,
and ampere does not correct for it. Neither does anything else; the remedy is
to know it, which is why :class:`VIEngine`'s docstring says so and why
``AutoMultivariateNormal`` is one keyword away.

The potential, and the improper prior
--------------------------------------
numpyro's SVI wants a *model* — a probabilistic program — where NUTS wanted a
potential function, and a realisation is neither: it is one callable of one
flat unconstrained vector (``inference.md`` §10a, deliberately, so that
ampere's own routing is not re-derived inside a sampler). The bridge is the
smallest numpyro model that can hold it: one site under an
``ImproperUniform(real_vector)`` prior — which contributes exactly zero to the
log-density — and one ``numpyro.factor`` carrying the realisation's value.
The model's log-density is then the realisation's, exactly, with no term
numpyro added and none ampere has to subtract; and because the site's support
is the whole real vector space, the autoguide's own ``biject_to`` is the
identity and the guide is fitted in precisely the coordinates the realisation
takes.
"""

from __future__ import annotations

import dataclasses
import math
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

__all__ = ["GUIDES", "VI_LIBRARIES", "VIEngine", "supported_backends"]

#: Which variational library drives which backend's realisation. Strings, not
#: imports, for the reason :data:`ampere.inference._nuts.SAMPLER_LIBRARIES`
#: gives: this module imports no backend and no sampler at module level.
#:
#: **jax only, today, and that is a fact about the code rather than about
#: torch.** pyro's ``SVI`` and its ``AutoNormal``/``AutoMultivariateNormal``
#: are the direct counterparts and would slot in beside numpyro's exactly as
#: pyro's NUTS did — one row here and one route in :meth:`VIEngine.run`. It is
#: not written because it cannot be *tested* from this branch: the torch
#: track's slice 2 is a sibling agent's, and a driver route nothing exercises
#: is a claim, not a feature. The refusal names torch by name and says this.
VI_LIBRARIES: dict[str, str] = {"jax": "numpyro"}

#: Guide family -> the numpyro autoguide class name it maps to.
#:
#: Two, and the choice between them is the one that matters: ``normal`` is
#: mean-field (independent per dimension, cheap, understates variance and
#: destroys correlations), ``multivariate_normal`` fits a full covariance
#: (O(d²) parameters, keeps linear correlations, still Gaussian). Names rather
#: than classes so this module imports numpyro only inside ``run``.
GUIDES: dict[str, str] = {
    "normal": "AutoNormal",
    "multivariate_normal": "AutoMultivariateNormal",
}

#: The one site name numpyro's model is given. See this module's docstring:
#: a realisation has no site structure by construction, so the whole
#: unconstrained vector travels under one name.
_SITE = "theta"

#: How many ELBO values a run records at most. A 50 000-step fit's full trace
#: is 50 000 floats of JSON in a netCDF attribute, which is a large thing to
#: put somewhere nothing can index; thinning to this many keeps the *shape* of
#: the optimisation — which is what the trace is read for — at a bounded cost,
#: and the stride is recorded beside it so nobody misreads the x axis.
MAX_ELBO_TRACE = 1000


@dataclasses.dataclass(frozen=True)
class _Settings:
    """One run's optimisation settings, kept together so a route cannot drift."""

    steps: int
    draws: int
    guide: str
    learning_rate: float
    progress: bool


def supported_backends() -> frozenset[str]:
    """The backends this driver can fit **right now**, in this interpreter.

    The registry ∩ :data:`VI_LIBRARIES`, computed when asked rather than when
    imported, for the reason :func:`ampere.inference._nuts.supported_backends`
    gives: a realisation is registered when the user imports the backend,
    which is normally after ``ampere.inference``.
    """
    return frozenset(registered_realisations()) & frozenset(VI_LIBRARIES)


class VIEngine(Engine):
    """Fit a :class:`~ampere.core.dataset.FittingProblem` by variational inference.

    Stochastic variational inference (Hoffman et al. 2013) over the
    *unconstrained* density: an autoguide is fitted by maximising the evidence
    lower bound with Adam, and the run's ``posterior`` group holds draws from
    the fitted guide.

    **This is an approximation, and the run says so.** The draws are from a
    Gaussian family — diagonal by default — not from the posterior, and the
    two differ in ways no diagnostic in the emitted file can fully expose. The
    known bias is one-directional: a mean-field guide **understates the
    posterior variance**, sometimes severely, because it cannot represent
    correlations and pays for a misfit in width. Credible intervals from a
    ``guide="normal"`` run are lower bounds on the true ones. Use it to
    explore, to initialise, or where a fit would otherwise be unaffordable;
    use :class:`~ampere.inference.NUTSEngine` for the number that goes in a
    paper, and compare the two when it matters.

    Parameters
    ----------
    problem
        The composed problem. Must report a backend in
        :func:`supported_backends` and must be ``differentiable``; both are
        checked at construction, by name.
    density
        Optional, and normally omitted — the backend's lowering of *problem*,
        exactly as :class:`~ampere.inference.NUTSEngine` takes it. Omitted, it
        comes from ``ampere.core.realise(problem)``. Passing one explicitly
        makes the run record ``ampere_realised = 0``, because it was not the
        registered realisation that produced the draws.
    cache_size
        See :class:`~ampere.inference.engine.Engine`. As for NUTS, this
        driver's cache is only ever populated by its start-point search.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.backends.jax import IndependentNoise, PowerLaw, configure_x64
    >>> from ampere.core import (
    ...     Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
    ... )
    >>> from ampere.inference import VIEngine
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
    >>> run = VIEngine(problem).run(draws=500, steps=2000)
    >>> run["posterior"]["model.norm"].shape
    (1, 500)
    >>> run.attrs["ampere_engine"], run.attrs["ampere_backend"]
    ('vi', 'jax')
    >>> run.attrs["ampere_vi_guide"]
    'normal'
    >>> bool(abs(float(run["posterior"]["model.norm"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "vi"
    #: Gradient-based, and told to ``check_engine`` as a statement about the
    #: *engine* rather than about the problem — the same reading NUTS takes.
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
            torch_note = (
                " The torch route is not written: pyro's SVI and autoguides are the direct "
                "counterparts and would slot in beside numpyro's, but nothing on this branch "
                "could exercise them, and an untested route is a claim rather than a feature."
                if problem.backend == "torch"
                else ""
            )
            raise EngineError(
                f"{self.NAME} cannot fit a problem on the {problem.backend!r} backend. This "
                f"driver needs two things: a registered realisation, so there is a "
                f"differentiable density at all (inference.md §10a), and a variational library "
                f"that can consume it — {VI_LIBRARIES}. Backends satisfying both in this "
                f"interpreter: {known}.{torch_note} Import the backend that supplies one and "
                f"build the problem from its models, instrument steps, noise models and GP "
                f"solver; or use a gradient-free engine (emcee, dynesty, zeus), which run every "
                f"backend on the contract path."
            )
        if not problem.differentiable:
            raise EngineError(
                f"{self.NAME} needs a differentiable problem, and this one declares "
                f"differentiable=False. That flag is aggregated from what the models, instrument "
                f"steps, noise models and GP solvers themselves declare, so one "
                f"non-differentiable piece is enough to make the whole chain "
                f"non-differentiable — problem.capabilities says which pieces were consulted."
            )
        #: The backend's realisation, when this driver obtained one; ``None``
        #: when the caller supplied a density, which is what makes
        #: ``ampere_realised`` a fact about the run.
        self.realisation: Realisation | None = None
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
        self._check_density_agrees()

    def _check_density_agrees(self) -> None:
        """The cheapest necessary condition that *density* is this problem's.

        The same one-point check :class:`~ampere.inference.NUTSEngine` makes,
        and for the same reason: a driver handed the wrong problem's density
        would fit a healthy-looking guide to something else entirely.
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
                f"the same quantity computed twice."
            )

    def _realisation_provenance(self) -> list[dict[str, Any]] | None:
        """The user-registered lowering rows the realisation consulted, if any."""
        if self.realisation is None:
            return None
        rows = getattr(self.realisation, "lowering_provenance", None)
        return list(rows()) if callable(rows) else None

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
        """Fit the guide, draw from it, and emit the run.

        Parameters
        ----------
        draws
            How many samples to draw from the **fitted guide**. They are
            independent by construction — a guide is a distribution, not a
            chain — so there is no burn-in, no thinning and no
            autocorrelation, and the run reports a single chain. R-hat and ESS
            are meaningless for such draws and ArviZ will still compute them;
            read ``vi_elbo_trace`` instead.
        steps
            ELBO optimisation steps. The default is generous rather than
            clever; watch the recorded trace and raise it if it has not
            flattened.
        guide
            ``"normal"`` (mean-field) or ``"multivariate_normal"``
            (full covariance). See :data:`GUIDES`.
        learning_rate
            Adam's step size.
        initial
            ``(n_dim,)`` start position **in the constrained space**, used to
            initialise the guide's location. The default draws one from the
            joint prior on this engine's own initialisation stream, so a run
            repeats exactly from the problem's seed.
        progress
            The optimiser's progress bar. Off by default: a driver that prints
            by default is unusable inside a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(1, draws)``.
        """
        if guide not in GUIDES:
            known = ", ".join(sorted(GUIDES))
            raise EngineError(
                f"{self.NAME} does not know the guide family {guide!r}; the choices are {known}. "
                f"'normal' is mean-field and understates the posterior variance; "
                f"'multivariate_normal' fits a full covariance and keeps linear correlations."
            )
        if int(steps) < 1:
            raise EngineError(f"{self.NAME} needs at least one optimisation step, got {steps}.")
        if not (float(learning_rate) > 0.0 and math.isfinite(float(learning_rate))):
            raise EngineError(
                f"{self.NAME}'s learning rate must be finite and positive, got {learning_rate!r}."
            )
        # `_kept` is the shared refusal of a combination that keeps nothing;
        # VI thins nothing and discards nothing, so the only condition left is
        # that it draws at least one sample.
        _kept(int(draws), 0, 1, self.NAME)
        self.start()

        start = self.initial_positions(1)[0] if initial is None else self._checked_initial(initial)
        unconstrained = self.problem.unconstrain(start)

        settings = _Settings(
            steps=int(steps),
            draws=int(draws),
            guide=str(guide),
            learning_rate=float(learning_rate),
            progress=bool(progress),
        )
        drawn, attrs = self._run_numpyro(unconstrained, settings)

        chain = np.stack([np.stack([self.problem.constrain(y) for y in drawn])])
        attrs.update(
            {
                "vi_draws": settings.draws,
                "vi_steps": settings.steps,
                "vi_guide": settings.guide,
                "vi_learning_rate": settings.learning_rate,
                "vi_optimiser": "adam",
                "vi_library": VI_LIBRARIES[self.problem.backend],
            }
        )
        return self.finish(
            chain,
            extra_attrs=attrs,
            realised=self.realisation is not None,
            registered_lowerings=self._realisation_provenance(),
            log_likelihood_terms=self._decomposition(drawn[np.newaxis, ...]),
        )

    # -- the one library route -----------------------------------------------

    def _run_numpyro(
        self, unconstrained: np.ndarray, settings: _Settings
    ) -> tuple[np.ndarray, dict[str, object]]:
        """numpyro's SVI over a jax realisation.

        Lazy imports, exactly as ``_nuts.py`` imports numpyro: ``import
        ampere.inference`` must not require jax (``architecture.md`` §4 rule
        2). The suppressions are the price of that rule — ``dev``, the
        environment CI typechecks in, deliberately has no jax — and the real
        check is ``pixi run -e jax typecheck``.

        The model is the smallest one that can carry a bare density: one site
        under ``ImproperUniform(real_vector)``, which adds nothing to the
        log-density, plus a ``factor`` carrying the realisation's value. The
        guide is fitted in exactly the coordinates the realisation takes,
        because that site's support is the whole space and the autoguide's
        ``biject_to`` is therefore the identity.
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

        def model() -> None:
            theta = numpyro.sample(
                _SITE,
                dist.ImproperUniform(dist.constraints.real_vector, (), (size,)),
            )
            numpyro.factor("ampere_log_density", density(theta))

        builder = getattr(autoguide, GUIDES[settings.guide])
        # `init_loc_fn` puts the guide's location at this engine's own seeded
        # draw from the prior, so a VI run repeats from the problem's seed in
        # the same way every other driver's does.
        fitted_guide = builder(
            model,
            init_loc_fn=init_to_value(values={_SITE: jax.numpy.asarray(unconstrained)}),
        )
        svi = SVI(
            model,
            fitted_guide,
            optim.Adam(settings.learning_rate),
            loss=Trace_ELBO(),
        )
        keys = jax.random.split(jax.random.key(self.integer_seed("optimiser")), 2)
        result = svi.run(keys[0], settings.steps, progress_bar=settings.progress)
        self.sampler = svi

        posterior = fitted_guide.sample_posterior(
            keys[1], result.params, sample_shape=(settings.draws,)
        )
        drawn = np.asarray(posterior[_SITE], dtype=float).reshape(settings.draws, size)

        # numpyro's SVI *minimises* the negative ELBO, so its losses are
        # -ELBO. Recorded as the ELBO itself, which is what a reader expects
        # to see rising towards a plateau.
        losses = np.asarray(result.losses, dtype=float)
        elbo = -losses
        stride = max(1, math.ceil(elbo.size / MAX_ELBO_TRACE))
        # Thinned from the **end**, so the last recorded value is always the
        # last step's: a trace whose tail an unlucky stride dropped would be a
        # trace that cannot answer the question it is read for ("had it
        # converged?").
        thinned = elbo[::-1][::stride][::-1]
        attrs: dict[str, object] = {
            "vi_elbo_final": float(elbo[-1]) if elbo.size else float("nan"),
            "vi_elbo_initial": float(elbo[0]) if elbo.size else float("nan"),
            "vi_elbo_trace": [float(value) for value in thinned],
            "vi_elbo_trace_stride": stride,
            "numpyro_version": numpyro.__version__,
            "jax_version": jax.__version__,
        }
        return drawn, attrs

    # -- shared with the NUTS driver, deliberately duplicated ---------------

    def _decomposition(self, unconstrained: np.ndarray) -> dict[str, np.ndarray] | None:
        """The per-dataset log-likelihood of every stored draw, natively.

        ``inference.md`` §10a's optional ``log_likelihood_terms``, used when
        it is there — the same consumption
        :meth:`ampere.inference.NUTSEngine._decomposition` makes, so a VI run
        records ``engine_draws_recomputed = 0`` as a NUTS run does.
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

    def _checked_initial(self, initial: Any) -> np.ndarray:
        position = np.asarray(initial, dtype=float).reshape(-1)
        if position.size != self.problem.free_size:
            raise EngineError(
                f"{self.NAME}'s start position must have {self.problem.free_size} element(s), "
                f"got {position.size}."
            )
        return position
