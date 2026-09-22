"""Two more nested samplers — nautilus and ultranest — behind one driver.

Private module; the classes are :class:`ampere.inference.NautilusEngine` and
:class:`ampere.inference.UltranestEngine`. **W5.14**, from
``docs/design/inference_extensions_memo.md`` §6's tier 1.

Why a second and a third nested sampler
---------------------------------------
:class:`~ampere.inference.DynestyEngine` has been here since W2.2 and stays the
base install's nested sampler. These two are not replacements for it; they are
the two the memo's §6 table puts at the top of tier 1, and each buys something
dynesty does not:

* **nautilus** (Lange 2023) is *importance* nested sampling with a neural
  boundary: it fits a small ensemble of neural regressors to the likelihood
  surface and samples the shells through them, which typically reaches a given
  effective sample size in far fewer likelihood calls than an ellipsoidal
  decomposition does. That matters most for exactly the problems ampere
  exists for — an expensive radiative-transfer forward model where the
  *number of evaluations* is the budget.
* **ultranest** (Buchner 2021) is MLFriends region sampling with a calibrated,
  bootstrapped termination criterion and an evidence error that separates its
  bootstrap and tail contributions. It is the conservative one: it is designed
  not to miss a mode, its error bars are meant to be believed, and it has an
  insertion-order test that says so when they should not be.

Both are **evidence-producing** and both explore the whole prior volume, so
multimodality is theirs by construction rather than by a tuning choice — the
property the memo's tier table puts beside them.

One driver, two libraries
-------------------------
:class:`_NestedEngine` owns everything nested sampling shares, which after
W5.0 is nearly all of it: the live-point default (:func:`~ampere.inference.
engine.default_live_points`, shared with dynesty since this item), the
equal-weight resampling rule, the engine-neutral evidence triple, the cost
record, and the emission. A subclass supplies three class attributes and one
method, :meth:`_NestedEngine._sample`, which runs its own library and returns
a :class:`_NestedOutput`. Nothing else differs, and that is the point: a run
archived from one of these is comparable with a run archived from the other,
and with a dynesty run, without a reader knowing which produced it.

The weighted-draw rule, applied the same way three times
--------------------------------------------------------
``results.md`` §9 (W5.0, from the memo's §5.2a) promoted dynesty's convention
to the contract: **weighted output is resampled to equal weight for the
``posterior`` group, the original count is recorded, and the raw weighted
output stays on the engine object.** Both drivers here obey it with the same
function dynesty's driver uses — ``dynesty.utils.resample_equal``, drawn on
this engine's own ``resample`` stream. dynesty is a base dependency, so
borrowing its resampler costs nothing and buys the thing that matters: the
rule is one *implementation*, not three, so three nested samplers cannot
drift into three slightly different posteriors from the same dead points.

ultranest offers an equal-weight ``result["samples"]`` of its own and nautilus
a ``posterior(equal_weight=True)``; neither is used, because both draw from
their library's own randomness rather than from the problem's seed, and a run
that is reproducible except for its final resampling is not reproducible.

Evidence, and the one place the two libraries genuinely differ
--------------------------------------------------------------
ultranest reports ``logz`` and ``logzerr`` directly, and its error is the one
to record: it already combines the bootstrap and tail contributions.

nautilus reports ``log_z`` and **no uncertainty at all**. Rather than record a
``nan`` for a quantity a reader will compare across engines, this driver
computes the standard importance-sampling estimate from nautilus's own
effective sample size. For weights ``w_i`` with ``Z = mean(w)``, the relative
standard error of ``Z`` is ``sqrt(Var(w)/n) / mean(w)``, which is
``1 / sqrt(ESS)`` for the Kish effective size ``ESS = (sum w)**2 / sum(w**2)``
— and ``Sampler.n_eff`` is exactly that Kish size, computed shell by shell
(``nautilus/sampler.py``). Since ``d log Z = dZ / Z``, the recorded
``ampere_log_evidence_err`` is ``1 / sqrt(n_eff)``. It is an estimate of the
*sampling* error of an importance sum and says nothing about a mode nautilus
never found, which is true of every nested sampler's error bar and is why the
engine battery checks all three against a closed-form evidence rather than
against each other. The value is recorded under the engine's own name too
(``ampere_nautilus_log_z_err_source`` names the estimator), so a reader is
never told that nautilus reported an uncertainty it did not.

Cost accounting
---------------
The memo's design horizon **(g)** asks for "one cost record every engine
writes", and names ``engine_evaluations`` as the one that already exists for
slot A. :meth:`~ampere.inference.engine.Engine.finish` writes it for every run
from the evaluation cache, so these two engines need add nothing to satisfy
the horizon. They do record the *sampler's own* count beside it —
``ampere_nautilus_likelihood_calls``, ``ampere_ultranest_likelihood_calls``,
exactly as dynesty records ``ampere_dynesty_ncall`` — because the two are not
the same number when a sampler evaluates points the cache then serves again,
and the difference is itself the thing a cost model wants. No contract change:
the memo's §10 ruling records that (e) to (g) need none today.

Extras, and the refusal
-----------------------
One extra per sampler (``inference_extensions_memo.md`` §10, ruled
2026-09-10), imported lazily **in the constructor** rather than in ``run`` —
:mod:`ampere.inference._zeus`'s reasoning, unchanged: the failure belongs at
the moment the user asks for the engine, not after they have chosen a live-
point count and pressed go.
"""

from __future__ import annotations

import contextlib
import dataclasses
import logging
import math
from collections.abc import Iterator, Mapping
from typing import Any, ClassVar

import numpy as np

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import OptionalDependencyError

from .engine import DEFAULT_CACHE_SIZE, Engine, default_live_points, global_seed
from .exceptions import EngineError

__all__ = ["NESTED_ENGINES", "NautilusEngine", "UltranestEngine"]

#: The method family both drivers write into ``ampere_evidence_method``.
#: The same string dynesty's driver writes, deliberately: W5.0's carried
#: question was whether that attribute should name the *family* or the engine,
#: and this item is the second evidence engine it asked to be reviewed
#: against. The family is kept — ``ampere_engine`` already names the engine,
#: so spelling it twice would buy nothing, while the family is what tells a
#: reader whether two archived evidences were estimated the same way.
EVIDENCE_METHOD = "nested_sampling"


@dataclasses.dataclass(frozen=True)
class _NestedOutput:
    """What a nested sampler hands back, in the one shape this driver stores.

    ``points`` and ``weights`` are the *weighted* dead points — the library's
    own output, before the equal-weight resampling ``results.md`` §9 requires
    — and ``weights`` are normalised to sum to one, which is what
    ``dynesty.utils.resample_equal`` expects.
    """

    points: np.ndarray
    weights: np.ndarray
    log_evidence: float
    log_evidence_err: float
    likelihood_calls: int
    attrs: dict[str, object]


def _from_log_weights(log_weights: Any) -> np.ndarray:
    """Importance weights from *log* weights, summing to one and never NaN."""
    values = np.asarray(log_weights, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        raise EngineError(
            "the nested sampler returned no finite importance weight, so there is no posterior "
            "to resample. That normally means every dead point scored -inf, which is a problem "
            "with the likelihood rather than with the sampler."
        )
    shifted = np.exp(values - float(np.max(finite)))
    shifted[~np.isfinite(shifted)] = 0.0
    return _renormalised(shifted)


def _renormalised(weights: Any) -> np.ndarray:
    """Linear importance weights, renormalised to sum to one.

    ``dynesty.utils.resample_equal`` checks that the weights sum to one to
    within a tolerance, and a library's own normalisation is only ever *close*
    to one; renormalising here rather than relying on that tolerance is what
    keeps a long run from failing on the last digit.
    """
    values = np.asarray(weights, dtype=float)
    values = np.where(np.isfinite(values) & (values > 0.0), values, 0.0)
    total = float(np.sum(values))
    if total <= 0.0:
        raise EngineError(
            "the nested sampler's importance weights summed to zero, so there is no posterior "
            "to resample. That normally means every dead point scored -inf, which is a problem "
            "with the likelihood rather than with the sampler."
        )
    return values / total


class _NestedEngine(Engine):
    """What nautilus and ultranest share. Not exported; see the module docstring.

    A subclass declares :attr:`NAME`, :attr:`MODULE` and :attr:`EXTRA`, and
    implements :meth:`_sample`.
    """

    #: The importable name of the library this driver drives.
    MODULE: ClassVar[str]
    #: The ``pip install "ampere[...]"`` extra that supplies it.
    EXTRA: ClassVar[str]

    def __init__(
        self,
        problem: FittingProblem,
        *,
        live_points: int | None = None,
        options: Mapping[str, Any] | None = None,
        cache_size: int = DEFAULT_CACHE_SIZE,
        use_realisation: bool = True,
    ) -> None:
        self._library = self._require()
        super().__init__(problem, cache_size=cache_size, use_realisation=use_realisation)
        chosen = default_live_points(problem.free_size) if live_points is None else int(live_points)
        if chosen < 2:
            raise EngineError(f"{self.NAME} needs at least 2 live points, got {chosen}.")
        self.live_points = chosen
        self.options: dict[str, Any] = dict(options or {})

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
                context=(
                    f"running the {cls.NAME} nested sampler (ampere.inference.{cls.__name__})"
                ),
            ) from error

    # -- the run --------------------------------------------------------------

    def run(self, *, progress: bool = False, **run_options: Any) -> Any:
        """Run to the sampler's own termination and emit the run.

        Parameters
        ----------
        progress
            The library's own progress reporting. Off by default, for the
            reason every driver here has it off by default: a driver that
            prints by default is unusable inside a loop or a test suite.
        **run_options
            Passed to the library's own ``run`` unchanged, and deliberately
            not re-spelled: nautilus stops on ``f_live``/``n_eff`` and
            ultranest on ``dlogz``/``dKL``/``frac_remain``/``min_ess``, and
            inventing one name for two genuinely different stopping criteria
            would hide the difference rather than abstract it. ``DynestyEngine
            .run`` takes ``**run_nested`` for the same reason.

        Returns
        -------
        xarray.DataTree
            The run, with one chain of equal-weight posterior draws.
        """
        from dynesty.utils import resample_equal

        self.start()
        output = self._sample(progress=bool(progress), **run_options)
        equal = np.asarray(
            resample_equal(output.points, output.weights, rstate=self.stream("resample")),
            dtype=float,
        )
        attrs: dict[str, object] = {
            f"{self.NAME}_live_points": self.live_points,
            f"{self.NAME}_dead_points": int(output.points.shape[0]),
            f"{self.NAME}_equal_weight_draws": int(equal.shape[0]),
            f"{self.NAME}_likelihood_calls": int(output.likelihood_calls),
            f"{self.NAME}_version": str(getattr(self._library, "__version__", "unknown")),
        }
        attrs.update(output.attrs)
        # The engine-neutral triple (results.md §9, W5.0). Written last so a
        # subclass cannot shadow it by accident: whatever a library calls its
        # own evidence, a reader finds it here under one name.
        attrs["log_evidence"] = float(output.log_evidence)
        attrs["log_evidence_err"] = float(output.log_evidence_err)
        attrs["evidence_method"] = EVIDENCE_METHOD
        return self.finish(equal[np.newaxis, ...], extra_attrs=attrs)

    def _sample(self, *, progress: bool, **run_options: Any) -> _NestedOutput:
        """Drive this subclass's library. See :class:`_NestedOutput`."""
        raise NotImplementedError


class NautilusEngine(_NestedEngine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with nautilus.

    Importance nested sampling with a neural boundary (Lange 2023): the shells
    are sampled through a small ensemble of neural regressors fitted to the
    likelihood surface, which typically reaches a target effective sample size
    in fewer likelihood calls than an ellipsoidal decomposition. For an
    expensive forward model that is the budget that matters.

    Parameters
    ----------
    problem
        The composed problem. **At least two free parameters**: nautilus
        refuses a one-dimensional problem (``nautilus/sampler.py``: "Cannot run
        Nautilus with less than 2 parameters"), and this driver refuses it
        first, by name, rather than letting a library ``ValueError`` out of the
        constructor.
    live_points
        The live set, nautilus's ``n_live``. Defaults to
        ``max(100, 25 (n_dim + 1))`` — :func:`~ampere.inference.engine.
        default_live_points`, shared with the other two nested samplers so
        that a comparison between them is a comparison of samplers.
        (nautilus's own default is 2000, which is generous for the small
        problems a first look uses and is a long wait for one.)
    options
        Passed to ``nautilus.Sampler`` unchanged — ``n_networks`` (0 disables
        the neural boundary entirely and falls back to plain ellipsoidal
        sampling, which is the right setting for a cheap low-dimensional
        problem), ``enlarge_per_dim``, ``split_threshold``, ``periodic``.
        ``filepath`` is not among them: this driver never writes a checkpoint,
        because a run output on disk is not something an engine should produce
        behind the caller's back (AGENTS.md ground rule 7 makes the same point
        about the repository).
    cache_size, use_realisation
        See :class:`~ampere.inference.engine.Engine`.

    Notes
    -----
    **The evidence uncertainty is ampere's, not nautilus's.** nautilus reports
    ``log_z`` and no error; the recorded ``ampere_log_evidence_err`` is
    ``1/sqrt(n_eff)``, the standard importance-sampling estimate over
    nautilus's own Kish effective sample size. The module docstring derives it,
    and ``ampere_nautilus_log_z_err_source`` records that it was estimated here
    rather than reported by the library.

    Examples
    --------
    >>> import numpy as np, scipy.stats as st, astropy.units as u
    >>> from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
    >>> class Line(Model):
    ...     def __init__(self, wavelength):
    ...         self.register_buffer("wavelength", wavelength, unit=u.um)
    ...         self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))
    ...         self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
    ...     def evaluate(self, **values):
    ...         ctx = self.context(values)
    ...         return Spectrum(
    ...             ctx["wavelength"] * u.um,
    ...             (ctx["slope"] * ctx["wavelength"] + ctx["offset"]) * u.Jy,
    ...         )
    >>> grid = np.array([1.0, 2.0, 3.0, 4.0])
    >>> observed = Spectrum(
    ...     grid * u.um, [2.0, 4.0, 6.0, 8.0] * u.Jy, uncertainty=[0.3] * 4 * u.Jy
    ... )
    >>> problem = FittingProblem(Line(grid), [Dataset(observed)], seed=20260922)
    >>> engine = NautilusEngine(problem, live_points=100, options={"n_networks": 0})
    >>> run = engine.run(n_eff=500)
    >>> run["posterior"]["model.slope"].sizes["chain"]
    1
    >>> bool(abs(float(run["posterior"]["model.slope"].mean()) - 2.0) < 0.4)
    True
    >>> run.attrs["ampere_evidence_method"]
    'nested_sampling'
    """

    NAME: ClassVar[str] = "nautilus"
    MODULE: ClassVar[str] = "nautilus"
    EXTRA: ClassVar[str] = "nautilus"

    def __init__(self, problem: FittingProblem, **kwargs: Any) -> None:
        super().__init__(problem, **kwargs)
        if self.problem.free_size < 2:
            raise EngineError(
                f"nautilus cannot sample a {self.problem.free_size}-dimensional problem: the "
                f"library refuses fewer than 2 free parameters, because its neural boundary is "
                f"fitted in the space it is asked to bound. This problem's only free parameter "
                f"is {next(iter(self.problem.free_labels()))!r}. Use DynestyEngine or "
                f"UltranestEngine, both of which sample one dimension happily, or free a second "
                f"parameter."
            )

    def _sample(self, *, progress: bool, **run_options: Any) -> _NestedOutput:
        # `filepath=None` and no `resume`: no checkpoint file appears beside
        # the caller's working directory. `seed` is an integer drawn from this
        # engine's own `sampler` stream, so the run repeats from the problem's
        # seed (nautilus builds its own Generator from it internally).
        sampler = self._library.Sampler(
            self.prior_transform,
            self.log_likelihood,
            n_dim=self.problem.free_size,
            n_live=self.live_points,
            seed=self.integer_seed("sampler"),
            filepath=None,
            **self.options,
        )
        sampler.run(verbose=bool(progress), **run_options)
        self.sampler = sampler

        points, log_weights, _ = sampler.posterior()
        n_eff = float(sampler.n_eff)
        if not (n_eff > 0.0):  # pragma: no cover - a run that produced nothing
            raise EngineError(
                "nautilus finished with an effective sample size of zero, so it found no "
                "posterior mass at all. Check the prior's support against the likelihood."
            )
        return _NestedOutput(
            points=np.asarray(points, dtype=float),
            weights=_from_log_weights(log_weights),
            log_evidence=float(sampler.log_z),
            # See the module docstring: 1/sqrt(Kish ESS) is the relative
            # standard error of an importance sum, and d(log Z) = dZ / Z.
            log_evidence_err=1.0 / math.sqrt(n_eff),
            likelihood_calls=int(sampler.n_like),
            attrs={
                "nautilus_log_z": float(sampler.log_z),
                "nautilus_n_eff": n_eff,
                "nautilus_log_z_err_source": "1/sqrt(n_eff); nautilus reports no uncertainty",
                "nautilus_n_networks": int(self.options.get("n_networks", 4)),
            },
        )


class UltranestEngine(_NestedEngine):
    """Sample a :class:`~ampere.core.dataset.FittingProblem` with ultranest.

    MLFriends region sampling (Buchner 2021) with a bootstrapped termination
    criterion. The conservative nested sampler of the three: it is built not to
    miss a mode, its evidence error separates the bootstrap and tail
    contributions, and its insertion-order test reports when the sampling was
    not uniform enough for either to be believed.

    Parameters
    ----------
    problem
        The composed problem. One free parameter is fine here.
    live_points
        ultranest's ``min_num_live_points``. Defaults to
        ``max(100, 25 (n_dim + 1))``, shared with the other two nested
        samplers; ultranest's own default is 400.
    options
        Passed to ``ultranest.ReactiveNestedSampler`` unchanged —
        ``wrapped_params`` for a circular parameter, ``num_bootstraps``,
        ``ndraw_min``/``ndraw_max``. ``log_dir`` is not among them: this
        driver never writes a run directory, for the reason
        :class:`NautilusEngine` never writes a checkpoint.
    cache_size, use_realisation
        See :class:`~ampere.inference.engine.Engine`.

    Notes
    -----
    **ultranest draws from numpy's process-global generator**, exactly as zeus
    does (``np.random.uniform``/``randint``/``choice`` throughout its region
    sampling, its step samplers and its bootstrap). So the run is wrapped in
    :func:`~ampere.inference.engine.global_seed`, which seeds and restores that
    state around it; with ``problem.seed=None`` the globals are left alone,
    which is the honest behaviour for a run that did not ask to be
    reproducible. Being global state, an ultranest run is not thread-safe
    against other code drawing from ``np.random`` at the same time —
    ``inference.md`` §12's limitation, now shared by two drivers.

    **The library logs at INFO** through a ``logging`` logger of its own,
    independently of ``show_status``. With ``progress=False`` that logger is
    quietened for the duration of the run and restored afterwards, so a driver
    asked not to print does not print.
    """

    NAME: ClassVar[str] = "ultranest"
    MODULE: ClassVar[str] = "ultranest"
    EXTRA: ClassVar[str] = "ultranest"

    def _sample(self, *, progress: bool, **run_options: Any) -> _NestedOutput:
        seed = None if self.problem.seed is None else self.integer_seed("sampler")
        with global_seed(seed):
            # Constructed *outside* the quietening block on purpose: the
            # constructor is what installs ultranest's own stream handler
            # (`ultranest.utils.create_logger`), so a block entered before it
            # would have no handler to quieten and the run would print anyway.
            sampler = self._library.ReactiveNestedSampler(
                list(self.problem.free_labels()),
                self.log_likelihood,
                transform=self.prior_transform,
                vectorized=False,
                log_dir=None,
                **self.options,
            )
            with _quiet_logger("ultranest", quiet=not progress):
                result = sampler.run(
                    min_num_live_points=self.live_points,
                    show_status=bool(progress),
                    viz_callback=False,
                    **run_options,
                )
        self.sampler = sampler

        weighted = result["weighted_samples"]
        insertion = result.get("insertion_order_MWW_test") or {}
        # `weights`, not `logw`. ultranest's `logw` is the log *prior volume*
        # element of each dead point, not its importance weight: resampling on
        # it returns the prior, which is what an early draft of this driver
        # stored (the battery's posterior row caught it -- the recovered
        # standard deviations were the prior's to three digits). `weights` is
        # the normalised importance weight, and is what `samples` is built
        # from inside ultranest.
        return _NestedOutput(
            points=np.asarray(weighted["points"], dtype=float),
            weights=_renormalised(weighted["weights"]),
            log_evidence=float(result["logz"]),
            log_evidence_err=float(result["logzerr"]),
            likelihood_calls=int(result["ncall"]),
            attrs={
                "ultranest_logz": float(result["logz"]),
                "ultranest_logzerr": float(result["logzerr"]),
                # The two halves ultranest keeps apart, kept apart here too: a
                # bootstrap error much smaller than the tail error means the
                # run stopped on its remaining prior volume rather than on its
                # sampling noise, which changes what a reader should conclude.
                "ultranest_logzerr_bs": float(result["logzerr_bs"]),
                "ultranest_logzerr_tail": float(result["logzerr_tail"]),
                "ultranest_information": float(result["H"]),
                "ultranest_ess": float(result["ess"]),
                "ultranest_niter": int(result["niter"]),
                # ultranest's own answer to "were my error bars trustworthy?".
                # Recorded rather than asserted on: it is a diagnostic for the
                # reader of an archived run, not a reason for this driver to
                # refuse one.
                "ultranest_insertion_order_converged": bool(insertion.get("converged", False)),
            },
        )


@contextlib.contextmanager
def _quiet_logger(name: str, *, quiet: bool) -> Iterator[None]:
    """Silence a library's own ``logging`` logger for the duration, then restore.

    ultranest reports through ``logging.getLogger("ultranest")`` with a stream
    handler it installs itself, which ``show_status=False`` does not touch. The
    handler levels are raised as well as the logger's, because the handler
    ultranest installs carries its own level and a logger level alone would not
    stop it.
    """
    if not quiet:
        yield
        return
    logger = logging.getLogger(name)
    level = logger.level
    handler_levels = [(handler, handler.level) for handler in logger.handlers]
    try:
        logger.setLevel(logging.WARNING)
        for handler, _ in handler_levels:
            handler.setLevel(logging.WARNING)
        yield
    finally:
        logger.setLevel(level)
        for handler, original in handler_levels:
            handler.setLevel(original)


#: The nested samplers this module adds, by ``NAME``. A plain mapping rather
#: than a registry: the memo's design horizon (g) proposes an engine registry
#: shaped like the realisation registry, and until that lands (it needs a
#: contract decision, not a dictionary) a table the battery can enumerate is
#: what the battery needs and all it needs.
NESTED_ENGINES: dict[str, type[_NestedEngine]] = {
    NautilusEngine.NAME: NautilusEngine,
    UltranestEngine.NAME: UltranestEngine,
}
