"""The dynesty driver: nested sampling over §4.5's ``prior_transform``.

Private module; the class is :class:`ampere.inference.DynestyEngine`.

Nested sampling is the one engine here that consumes ``prior_transform``
rather than ``log_prior``, which is exactly why ``inference.md`` §10 lists both
on the surface: the same problem object serves an ensemble sampler that wants a
density and a nested sampler that wants an inverse CDF, and neither knows what
the other needs.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem

from .engine import DEFAULT_CACHE_SIZE, Engine
from .exceptions import EngineError

__all__ = ["DynestyEngine"]


def _default_live_points(free_size: int) -> int:
    """25 per dimension plus a floor of 100.

    Below roughly ``25 (n_dim + 1)`` the ellipsoidal bound is fitted from too
    few points to be trustworthy, which is where nested sampling starts to
    under-cover rather than merely run slowly; the floor keeps a one- or
    two-dimensional problem from being sampled by a handful of points.
    """
    return max(100, 25 * (free_size + 1))


class DynestyEngine(Engine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with dynesty.

    Nested sampling: slower per posterior draw than an ensemble sampler, and
    worth it for a multimodal posterior, for a problem where a chain's
    convergence is hard to judge, or when the marginal likelihood is wanted —
    which the ensemble engines cannot give at all. The log-evidence and its
    uncertainty are recorded in the run's provenance attrs.

    Parameters
    ----------
    problem
        The composed problem.
    live_points
        The live set. Defaults to ``max(100, 25 (n_dim + 1))``.
    dynamic
        Use ``dynesty.DynamicNestedSampler``, which spends its later
        iterations where they most improve the posterior rather than the
        evidence. Better posteriors for the same budget; a less clean evidence
        estimate.
    bound, sample
        Passed to dynesty unchanged.
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

    Notes
    -----
    **The emitted posterior is resampled to equal weight.** Nested sampling
    produces *weighted* dead points, and ArviZ's ``posterior`` group has no
    weight axis — every consumer of it, from ``arviz.summary`` to a corner
    plot, assumes the draws are equally weighted, so storing the raw dead
    points there would silently misreport every posterior summary. The stored
    draws are therefore ``dynesty.utils.resample_equal`` of the dead points,
    drawn on this engine's own ``resample`` stream, and the attrs record the
    original count (``ampere_dynesty_dead_points``) beside the resampled one so
    the reduction is visible. The full weighted output is on :attr:`sampler`
    (``engine.sampler.results``) for anything that wants it.

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
    >>> run = DynestyEngine(problem, live_points=60).run()
    >>> run["posterior"]["model.slope"].sizes["chain"]
    1
    >>> bool(abs(float(run["posterior"]["model.slope"].mean()) - 2.0) < 0.1)
    True
    >>> isinstance(run.attrs["ampere_dynesty_logz"], float)
    True
    """

    NAME: ClassVar[str] = "dynesty"

    def __init__(
        self,
        problem: FittingProblem,
        *,
        live_points: int | None = None,
        dynamic: bool = False,
        bound: str = "multi",
        sample: str = "auto",
        cache_size: int = DEFAULT_CACHE_SIZE,
        use_realisation: bool = True,
    ) -> None:
        super().__init__(problem, cache_size=cache_size, use_realisation=use_realisation)
        chosen = (
            _default_live_points(problem.free_size) if live_points is None else int(live_points)
        )
        if chosen < 2:
            raise EngineError(f"dynesty needs at least 2 live points, got {chosen}.")
        self.live_points = chosen
        self.dynamic = bool(dynamic)
        self.bound = bound
        self.sample = sample

    def run(
        self,
        *,
        maxiter: int | None = None,
        maxcall: int | None = None,
        progress: bool = False,
        **run_nested: Any,
    ) -> Any:
        """Run to termination and emit the run.

        Parameters
        ----------
        maxiter, maxcall
            Hard budgets, passed to dynesty. ``None`` means "run to dynesty's
            own stopping criterion", which is the right default: a nested
            sampler stopped early has not merely fewer draws but a biased
            evidence.
        progress
            dynesty's ``print_progress``. Off by default.
        **run_nested
            Anything else ``run_nested`` takes — ``dlogz`` for the static
            sampler, ``dlogz_init``/``nlive_batch`` for the dynamic one.
            Passed through untouched rather than re-spelled here, because the
            two samplers' stopping criteria are genuinely different and
            inventing one name for both would hide that.

        Returns
        -------
        xarray.DataTree
            The run, with one chain of equal-weight posterior draws.
        """
        import dynesty
        from dynesty.utils import resample_equal

        self.start()
        # Typed as Any because dynesty's two entry points are *factory
        # functions* that build different sampler classes from one keyword
        # surface; a checker resolving the union lands on the shared base
        # class, whose __init__ takes none of the arguments the factories do.
        factory: Any = dynesty.DynamicNestedSampler if self.dynamic else dynesty.NestedSampler
        sampler = factory(
            self.log_likelihood,
            self.prior_transform,
            self.problem.free_size,
            nlive=self.live_points,
            bound=self.bound,
            sample=self.sample,
            rstate=self.stream("sampler"),
        )
        sampler.run_nested(
            maxiter=maxiter,
            maxcall=maxcall,
            print_progress=progress,
            **run_nested,
        )
        self.sampler = sampler

        results = sampler.results
        dead = np.asarray(results.samples, dtype=float)
        weights = np.asarray(results.importance_weights(), dtype=float)
        equal = resample_equal(dead, weights, rstate=self.stream("resample"))
        return self.finish(
            equal[np.newaxis, ...],
            extra_attrs={
                "dynesty_live_points": self.live_points,
                "dynesty_dynamic": self.dynamic,
                "dynesty_bound": self.bound,
                "dynesty_sample": self.sample,
                "dynesty_dead_points": int(dead.shape[0]),
                "dynesty_equal_weight_draws": int(equal.shape[0]),
                "dynesty_logz": float(results.logz[-1]),
                "dynesty_logzerr": float(results.logzerr[-1]),
                "dynesty_niter": int(results.niter),
                "dynesty_ncall": int(np.sum(results.ncall)),
                "dynesty_version": dynesty.__version__,
            },
        )
