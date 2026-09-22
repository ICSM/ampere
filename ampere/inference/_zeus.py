"""The zeus driver: ensemble slice sampling over §4.5's ``log_prob``.

Private module; the class is :class:`ampere.inference.ZeusEngine`.

zeus is the one engine here that is **not** a base dependency: it stays behind
the ``zeus`` extra (``architecture.md`` §3's table), so this module imports it
lazily, inside the constructor, and raises
:class:`~ampere.core.exceptions.OptionalDependencyError` naming the extra
(§4 rules 2 and 3). Importing it in ``__init__`` rather than in ``run`` is
deliberate: the failure belongs at the moment the user asks for the engine, not
after they have chosen their walker count and pressed go.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import OptionalDependencyError

from .engine import (
    DEFAULT_CACHE_SIZE,
    Engine,
    _check_ensemble,
    _default_walkers,
    _kept,
    global_seed,
)
from .exceptions import EngineError

__all__ = ["ZeusEngine"]


class ZeusEngine(Engine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with zeus.

    Ensemble *slice* sampling: no acceptance/rejection step and no hand-tuned
    proposal scale, which typically buys a shorter autocorrelation time than
    emcee on a correlated posterior, at the cost of several ``log_prob``
    evaluations per walker per step rather than one.

    Requires the ``zeus`` extra -- ``pip install ".[zeus]"`` from a checkout
    (PyPI's ``ampere`` package is unrelated).

    Parameters
    ----------
    problem
        The composed problem.
    walkers
        Ensemble size; same rule and default as
        :class:`~ampere.inference.EmceeEngine` (even, at least ``2 x n_dim``).
    moves
        Passed to ``zeus.EnsembleSampler`` unchanged.
    cache_size
        See :class:`~ampere.inference.engine.Engine`. There is no ``backend=``:
        W2.12 made the backend §4.5's fourth capability flag, read off the
        problem's own pieces.
    use_realisation
        See :class:`~ampere.inference.engine.Engine`. ``True`` by default: on a
        backend with a registered realisation this driver scores every proposal
        through it rather than through the numpy contract path, which on a jax
        quasiseparable problem is the difference between 27 ms and 0.6 ms a
        proposal. It changes no number and one record — a failure loses its
        reason there (``inference.md`` §10a) — so pass ``False``, or declare
        ``strict=True`` on the problem, to keep the reasons.
    **sampler_settings
        Anything else ``zeus.EnsembleSampler`` takes — ``tune``, ``tolerance``,
        ``maxsteps``, ``maxiter``, ``mu`` — forwarded untouched. ampere has no
        opinion about zeus's tuning and does not interpose defaults for it.

    Notes
    -----
    zeus is more sensitive to its start points than emcee is. Slice sampling
    expands an interval until it brackets the slice, so walkers drawn from a
    prior much wider than the posterior can exhaust the expansion budget before
    they reach the mode — which is a normal situation, not a pathological one,
    whenever the data are far more informative than the prior. The default
    start here is a draw from the joint prior (the engine-neutral choice, and
    the one that needs no tuning); when zeus refuses it, :meth:`run` says so
    and names the remedies, of which passing ``initial=`` with a tight ball
    around a reasonable guess is zeus's own recommendation.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
    >>> class Line(Model):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength, unit=u.um)
    ...         self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))
    ...     def evaluate(self, **values):
    ...         ctx = self.context(values)
    ...         return Spectrum(
    ...             ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy
    ...         )
    >>> grid = np.array([1.0, 2.0, 3.0])
    >>> observed = Spectrum(
    ...     grid * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.3, 0.3, 0.3] * u.Jy
    ... )
    >>> problem = FittingProblem(Line(grid), [Dataset(observed)], seed=20260905)
    >>> run = ZeusEngine(problem, walkers=8).run(steps=80, burn_in=30)
    >>> run["posterior"]["model.slope"].shape
    (8, 50)
    >>> bool(abs(float(run["posterior"]["model.slope"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "zeus"

    def __init__(
        self,
        problem: FittingProblem,
        *,
        walkers: int | None = None,
        moves: Any = None,
        cache_size: int = DEFAULT_CACHE_SIZE,
        use_realisation: bool = True,
        **sampler_settings: Any,
    ) -> None:
        if "backend" in sampler_settings:
            raise EngineError(
                "zeus takes no backend= setting, and neither does this driver any more: since "
                "W2.12 the backend is a capability flag the problem's own models and "
                "transformations declare, so ampere_backend is derived rather than asserted. "
                "Drop the argument; compose the problem from the backend you meant instead."
            )
        self._zeus = _require_zeus()
        super().__init__(problem, cache_size=cache_size, use_realisation=use_realisation)
        chosen = _default_walkers(problem.free_size) if walkers is None else int(walkers)
        self.walkers = _check_ensemble(self.NAME, chosen, problem.free_size)
        self.moves = moves
        self.sampler_settings = sampler_settings

    def run(
        self,
        steps: int,
        *,
        burn_in: int = 0,
        thin: int = 1,
        initial: Any = None,
        progress: bool = False,
    ) -> Any:
        """Run the ensemble for *steps* steps and emit the run.

        Parameters and return value are :meth:`~ampere.inference.EmceeEngine.
        run`'s, deliberately: the two ensemble engines differ in how they move,
        not in what a user has to say to them, and a driver layer that
        re-spelled the same four settings twice would be one more thing to get
        wrong when switching engines to see whether the answer changes.
        """
        zeus = self._zeus
        _kept(int(steps), int(burn_in), int(thin), self.NAME)
        self.start()
        positions = (
            self.initial_positions(self.walkers)
            if initial is None
            else self._checked_initial(initial)
        )

        sampler = zeus.EnsembleSampler(
            self.walkers,
            self.problem.free_size,
            self.log_prob,
            moves=self.moves,
            verbose=False,
            **self.sampler_settings,
        )
        self.sampler = sampler
        with global_seed(None if self.problem.seed is None else self.integer_seed("sampler")):
            try:
                sampler.run_mcmc(positions, int(steps), progress=progress)
            except RuntimeError as error:
                raise self._explain(error) from error

        # zeus stores (step, walker, dim), as emcee does; ArviZ wants
        # (chain, draw, dim).
        chain = np.swapaxes(sampler.get_chain(discard=int(burn_in), thin=int(thin)), 0, 1)
        return self.finish(
            chain,
            extra_attrs={
                "zeus_walkers": self.walkers,
                "zeus_steps": int(steps),
                "zeus_burn_in": int(burn_in),
                "zeus_thin": int(thin),
                "zeus_seeded_numpy_global_rng": self.problem.seed is not None,
                "zeus_version": zeus.__version__,
            },
        )

    def _checked_initial(self, initial: Any) -> np.ndarray:
        positions = np.asarray(initial, dtype=float)
        expected = (self.walkers, self.problem.free_size)
        if positions.shape != expected:
            raise EngineError(
                f"zeus start positions must be shaped {expected} (walkers, n_dim), got "
                f"{positions.shape}."
            )
        return positions

    def _explain(self, error: RuntimeError) -> EngineError:
        """Turn zeus's slice-expansion refusal into advice, keeping the cause attached.

        zeus's own message ("Number of expansions exceeded maximum limit") names
        neither of the two things a user can do about it, and the commoner cause
        is not a malformed pdf at all: an ensemble slice sampler expands its
        interval until it brackets the slice, so walkers started deep in the
        prior's tail of a *sharply* informative posterior — which is the normal
        situation when the data are much more informative than the prior — can
        exhaust the expansion budget before they ever reach the mode. emcee
        tolerates the same start; zeus does not, and that difference is worth
        stating rather than leaving a user to conclude their model is broken.
        """
        if "expansions" not in str(error) and "maxiter" not in str(error):
            return EngineError(f"zeus stopped: {error}")
        return EngineError(
            f"zeus stopped: {error}\n"
            f"This is usually the start points rather than the model. zeus expands its slice "
            f"until it brackets, so walkers drawn from a prior much wider than the posterior "
            f"can exhaust the budget before reaching the mode -- ampere's default start is a "
            f"draw from the joint prior, which emcee tolerates and zeus often does not. Three "
            f"remedies, in order of preference: pass initial= with walkers in a small ball "
            f"around a reasonable guess (zeus's own recommendation); raise the budget with "
            f"ZeusEngine(..., maxsteps=..., maxiter=...), which is passed to "
            f"zeus.EnsembleSampler untouched; or check with problem.failure_summary() and "
            f"FittingProblem(..., strict=True) that the log-probability really is finite and "
            f"well defined over the prior's support."
        )


def _require_zeus() -> Any:
    """Import zeus on use, never on import (``architecture.md`` §4, rules 2 and 3)."""
    try:
        import zeus
    except ImportError as error:  # pragma: no cover - exercised by the minimal-install job
        raise OptionalDependencyError(
            "zeus",
            extra="zeus",
            context="running the zeus ensemble slice sampler (ampere.inference.ZeusEngine)",
        ) from error
    return zeus
