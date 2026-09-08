"""The emcee driver: an affine-invariant ensemble over §4.5's ``log_prob``.

Private module; the class is :class:`ampere.inference.EmceeEngine`. Named with
a leading underscore so that ``ampere/inference/emcee.py`` never sits on the
import path shadowing the package this module imports — absolute imports make
that safe in principle, and unambiguous names make it safe in practice.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem

from .engine import DEFAULT_CACHE_SIZE, Engine, _check_ensemble, _default_walkers, _kept
from .exceptions import EngineError

__all__ = ["EmceeEngine"]


class EmceeEngine(Engine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with emcee.

    Goodman & Weare's affine-invariant ensemble sampler: the default for a
    problem of a few to a few tens of dimensions whose ``log_prob`` is cheap
    enough to call a few hundred thousand times, and the engine legacy ampere
    used most.

    Parameters
    ----------
    problem
        The composed problem.
    walkers
        Ensemble size. The default is four per free dimension (at least
        eight), comfortably above emcee's ``2 x n_dim`` floor. Must be even
        and at least ``2 x n_dim``: each step proposes for one half of the
        ensemble using the other half, so a half that does not span the space
        cannot leave the subspace it starts in.
    moves
        Passed to ``emcee.EnsembleSampler`` unchanged — ampere has no opinion
        about the move mixture and does not interpose one.
    use_realisation
        See :class:`~ampere.inference.engine.Engine`. ``True`` by default: on a
        backend with a registered realisation this driver scores every proposal
        through it rather than through the numpy contract path, which on a jax
        quasiseparable problem is the difference between 27 ms and 0.6 ms a
        proposal. It changes no number and one record — a failure loses its
        reason there (``inference.md`` §10a) — so pass ``False``, or declare
        ``strict=True`` on the problem, to keep the reasons.
    cache_size
        See :class:`~ampere.inference.engine.Engine`. Note that this driver
        takes no ``backend=`` at all — neither ampere's (W2.12 made that
        §4.5's fourth capability flag, read off the problem) nor emcee's own
        (its HDF5 store), which this driver does not use: the run is emitted
        as an ArviZ ``DataTree`` through ``ampere.results``, which is the
        single results format.

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
    >>> run = EmceeEngine(problem, walkers=8).run(steps=300, burn_in=100)
    >>> run["posterior"]["model.slope"].shape
    (8, 200)
    >>> bool(abs(float(run["posterior"]["model.slope"].mean()) - 2.0) < 0.1)
    True
    """

    NAME: ClassVar[str] = "emcee"

    def __init__(
        self,
        problem: FittingProblem,
        *,
        walkers: int | None = None,
        moves: Any = None,
        cache_size: int = DEFAULT_CACHE_SIZE,
        use_realisation: bool = True,
    ) -> None:
        super().__init__(problem, cache_size=cache_size, use_realisation=use_realisation)
        chosen = _default_walkers(problem.free_size) if walkers is None else int(walkers)
        self.walkers = _check_ensemble(self.NAME, chosen, problem.free_size)
        self.moves = moves

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

        Parameters
        ----------
        steps
            Steps per walker, burn-in included.
        burn_in
            Leading steps to discard. Discarded before emission, so they never
            reach the stored run — but their failures are counted, because a
            burn-in that could not be scored is still something the user needs
            told.
        thin
            Keep every *thin*-th retained step.
        initial
            ``(walkers, n_dim)`` start positions. The default draws them from
            the joint prior on this engine's own initialisation stream.
        progress
            emcee's progress bar. Off by default: a driver that prints by
            default is unusable inside a loop or a test suite.

        Returns
        -------
        xarray.DataTree
            The run, with ``(chain, draw)`` = ``(walkers, kept steps)``.
        """
        import emcee

        _kept(int(steps), int(burn_in), int(thin), self.NAME)
        self.start()
        positions = (
            self.initial_positions(self.walkers)
            if initial is None
            else self._checked_initial(initial)
        )

        sampler = emcee.EnsembleSampler(
            self.walkers,
            self.problem.free_size,
            self.log_prob,
            moves=self.moves,
        )
        self._seed(sampler)
        sampler.run_mcmc(positions, int(steps), progress=progress)
        self.sampler = sampler

        # emcee stores (step, walker, dim); ArviZ wants (chain, draw, dim).
        chain = np.swapaxes(sampler.get_chain(discard=int(burn_in), thin=int(thin)), 0, 1)
        return self.finish(
            chain,
            extra_attrs={
                "emcee_walkers": self.walkers,
                "emcee_steps": int(steps),
                "emcee_burn_in": int(burn_in),
                "emcee_thin": int(thin),
                "emcee_mean_acceptance_fraction": float(np.mean(sampler.acceptance_fraction)),
                "emcee_version": emcee.__version__,
            },
        )

    def _checked_initial(self, initial: Any) -> np.ndarray:
        positions = np.asarray(initial, dtype=float)
        expected = (self.walkers, self.problem.free_size)
        if positions.shape != expected:
            raise EngineError(
                f"emcee start positions must be shaped {expected} (walkers, n_dim), got "
                f"{positions.shape}."
            )
        return positions

    def _seed(self, sampler: Any) -> None:
        """Put emcee's own randomness on this run's seed, and check that it took.

        emcee 3 has no seed argument; the supported route is its
        ``random_state`` property, which wraps a legacy
        :class:`numpy.random.RandomState`. That setter **fails silently** by
        design (its own docstring says so), so this verifies the state
        afterwards rather than assuming: a reproducibility guarantee that
        quietly does not hold is worse than none.
        """
        seed = self.integer_seed("sampler")
        intended = np.random.RandomState(seed).get_state()
        sampler.random_state = intended
        applied = sampler.random_state
        if applied is None or not np.array_equal(applied[1], intended[1]):
            raise EngineError(
                "emcee refused the random state derived from this problem's seed, so the run "
                "would not be reproducible. This means emcee's `random_state` setter (which "
                "fails silently by design) rejected a legacy RandomState tuple — check the "
                f"installed emcee version ({getattr(sampler, '__module__', 'unknown')})."
            )
