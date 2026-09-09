"""Family D of ``diagnostics.md``: simulation-based calibration and coverage.

``diagnostics.md`` §11 records posterior calibration as a fourth diagnostic
family and is explicit about what makes it different from families A—C: those
diagnose the **model** — where a smooth model, or the fit's noise budget, is
inadequate for these data — while this one diagnoses the **inference
machinery**. Does the posterior an engine produces actually have the coverage
it claims, whether or not the model is right? An amortised neural posterior
can be silently overconfident in a way no residual test detects: the fit to the
*observed* data can look excellent while the posterior's claimed uncertainties
are fiction.

Where it lives, and why here
----------------------------
§11 left the placement to be "decided when the code lands, per the §8 table's
logic", and that logic settles it. This family consumes ``InferenceData``\\ s
plus ``simulate()``, produces numbers and pictures, and carries **no new
dependency** — the ``sbi`` route needs ``sbi``, which the ``sbi`` extra already
brings for the engine that produced the posterior in the first place, and the
generic route needs nothing at all. That is precisely the argument §8 makes for
families B and C living in :mod:`ampere.results`, so family D lives here too,
with the ``sbi``-specific fast path inside the engine
(:meth:`ampere.inference.SBIEngine.calibrate`) where the trained posterior and
its :class:`~ampere.core.encoding.EncodingLayout` already are.

Its outputs — a rank histogram and a coverage curve — are not
coordinate-indexed deficiency maps, so, like family B and for the same reason,
this family does **not** adopt the :class:`~ampere.core.AnomalyScore`
convention. Forcing a rank statistic into a container shaped for "how bad is
the model *here*" would be exactly the conflation ``diagnostics.md`` Tension 5
warns about.

The two routes
--------------
Both answer the same question and both return the **same**
:class:`xarray.Dataset` (:func:`calibration_dataset`), so one pair of plots and
one stored group serve either.

1. **The ``sbi`` fast path** — :meth:`ampere.inference.SBIEngine.calibrate`.
   An amortised posterior can be re-conditioned on a fresh simulated dataset
   for free, so the whole of SBC costs one extra simulation batch and no
   retraining at all. ``sbi.diagnostics``'s own ``run_sbc``/``check_sbc`` and
   ``run_tarp``/``check_tarp`` do the arithmetic — the
   depend-don't-reimplement posture ``diagnostics.md`` §2.2 takes for RHMF —
   and the batch is encoded **with the run's own layout**, so what is tested is
   the network as trained rather than a network shown columns in another order.

2. **The generic path** — :func:`sbc`. The Talts et al. (2018,
   arXiv:1804.06788) loop: draw θ from the prior, simulate a dataset there,
   **fit it from scratch**, and rank the true θ among the posterior draws. One
   full fit per simulation, so it is expensive by design and budget-controlled
   by the caller. It is also the only route that can validate a *likelihood* —
   the flexible-GP likelihood of milestone M2, or the profiled statistic
   ``examples/wstat_comparison.py`` builds — because there the thing under test
   is the density, not a network.

What a rank is, and what a coverage curve is
--------------------------------------------
For one simulation and one scalar parameter, the **rank** is the number of the
``L`` posterior draws that fall strictly below the true value:

.. code-block:: text

    r = #{ theta_l : theta_l < theta_true },      r in {0, ..., L}

Under a calibrated posterior ``r`` is uniform on ``{0, ..., L}``, and that is
the whole of Talts et al.'s test. Ranks below the mode and above it are the two
signatures worth naming: a **U-shaped** histogram means the posteriors are too
narrow (the truth keeps landing in the tails), a **inverted-U-shaped** one means they
are too wide, and a slope means they are biased.

The **coverage curve** is the same information read the other way, and it is
derived from the ranks rather than computed separately so that the two halves
of one figure can never disagree. The truth lies inside the central credible
interval of nominal level ``alpha`` exactly when its rank fraction ``u = r / L``
lies in ``[(1 - alpha)/2, (1 + alpha)/2]``, so

.. code-block:: text

    coverage(alpha) = mean over simulations of  1[ (1-a)/2 <= r/L <= (1+a)/2 ]

and a calibrated posterior puts that curve on the diagonal. Route 1 stores
TARP's *joint* expected-coverage curve (Lemos et al. 2023, arXiv:2302.03026)
beside it, which the marginal curve cannot replace: a posterior can be
perfectly calibrated in every one-dimensional margin and wrong about the
correlations between them.

The uniformity p-value is a Kolmogorov—Smirnov test of the ranks against
``Uniform(0, L)``, per parameter — ``sbi``'s ``check_sbc`` on route 1, the same
test spelled with :func:`scipy.stats.kstest` on route 2, so the number means
the same thing either way. It is a *frequentist* p-value on a small sample and
should be read as one: a low value is evidence of miscalibration, a high value
at 20 simulations is evidence of very little.

Storage
-------
:func:`attach_calibration` writes the dataset into a run's ``DataTree`` as the
``calibration`` group — ``results.md`` §4's schema table, extended by this item
— so a calibrated run archives its own calibration and a reader needs no second
file to know whether to believe it.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core.dataset import FittingProblem
from ampere.core.exceptions import OptionalDependencyError, ResultsError
from ampere.core.results_schema import FunctionSamples
from ampere.core.rng import generator, substream

from .emission import CHAIN_DIM, DRAW_DIM
from .provenance import ATTR_PREFIX

__all__ = [
    "CALIBRATION_GROUP",
    "CALIBRATION_SCHEMA_VERSION",
    "CALIBRATION_STREAM",
    "DEFAULT_LEVELS",
    "LEVEL_DIM",
    "PARAMETER_DIM",
    "REFIT_ROUTE",
    "SBI_ROUTE",
    "SIMULATION_DIM",
    "TARP_LEVEL_DIM",
    "attach_calibration",
    "calibration_dataset",
    "coverage_from_ranks",
    "replace_observations",
    "sbc",
    "uniformity_pvalues",
]

#: The group a calibration result is stored in (``results.md`` §4).
CALIBRATION_GROUP = "calibration"

#: Bumped when the group's variables or their meaning change. Distinct from
#: ``PROVENANCE_SCHEMA_VERSION`` because a calibration group can be written
#: into a run emitted under an older provenance schema and must still say what
#: shape it is.
CALIBRATION_SCHEMA_VERSION = 1

#: Dimension names. ``simulation`` indexes the fresh datasets the check was run
#: over, ``parameter`` the scalar columns ranked, ``level`` the nominal
#: credible levels of the marginal coverage curve, and ``tarp_level`` TARP's
#: own credibility grid — a separate dimension because ``run_tarp`` chooses its
#: bin count from the sample size and it is not this module's to align.
SIMULATION_DIM = "simulation"
PARAMETER_DIM = "parameter"
LEVEL_DIM = "level"
TARP_LEVEL_DIM = "tarp_level"

#: The nominal credible levels the coverage curve is reported at. Twenty-one
#: points from 0 to 1 inclusive: the endpoints are where a curve that has gone
#: wrong is most visible, and both are exactly attainable from ranks.
DEFAULT_LEVELS = np.linspace(0.0, 1.0, 21)

#: The sub-stream name every random draw in this module derives from
#: (``lowering.md`` §9.2). Distinct from ``"simulate"`` on purpose: running a
#: calibration check must not change what an SBI budget simulated.
CALIBRATION_STREAM = "calibration"

#: The two routes, as recorded in ``ampere_calibration_route``.
SBI_ROUTE = "sbi"
REFIT_ROUTE = "refit"

_SBC_WARNING_FLOOR = 100


def _require_xarray() -> Any:
    """Import xarray on use, never on import.

    xarray arrives with arviz, which W2.2 promoted into the base install, so
    ``extra=None``: an environment without it is incomplete rather than missing
    an extra. Lazy for the reason :mod:`ampere.results`'s docstring gives for
    arviz's own import.
    """
    try:
        import xarray
    except ImportError as error:  # pragma: no cover - exercised by a minimal install
        raise OptionalDependencyError(
            "xarray",
            context="assembling a calibration result (ampere.results.calibration returns an "
            "xarray.Dataset). xarray arrives with arviz, which is a base dependency of ampere, "
            "so this environment is incomplete rather than merely missing an extra",
        ) from error
    return xarray


# ---------------------------------------------------------------------------
# The statistics
# ---------------------------------------------------------------------------


def coverage_from_ranks(
    ranks: Any, posterior_draws: int, levels: Any = None
) -> tuple[np.ndarray, np.ndarray]:
    """The marginal central-interval coverage curve, from SBC ranks alone.

    Derived rather than measured, which is the point: the rank fraction
    ``u = r / L`` is the posterior's own CDF evaluated at the truth, so the
    truth lies inside the central credible interval of nominal level ``alpha``
    exactly when ``(1 - alpha)/2 ≤ u ≤ (1 + alpha)/2``. Computing the curve from the
    same ranks the histogram draws means the two panels of one figure cannot
    disagree about what was found.

    Parameters
    ----------
    ranks
        ``(simulations, parameters)`` integer ranks in ``{0, ..., L}``.
    posterior_draws
        ``L``: how many posterior draws each rank was taken against.
    levels
        The nominal levels to report at; :data:`DEFAULT_LEVELS` by default.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray]
        The levels used, and the ``(levels, parameters)`` empirical coverage.
    """
    grid = np.asarray(DEFAULT_LEVELS if levels is None else levels, dtype=float).reshape(-1)
    if grid.size == 0:
        raise ResultsError("a coverage curve needs at least one nominal level.")
    if np.any(grid < 0.0) or np.any(grid > 1.0):
        raise ResultsError(
            f"nominal credible levels are probabilities in [0, 1]; got "
            f"{np.min(grid):.4g} to {np.max(grid):.4g}."
        )
    draws = int(posterior_draws)
    if draws < 1:
        raise ResultsError(f"a rank is taken against at least one posterior draw, got {draws}.")
    table = np.asarray(ranks, dtype=float)
    if table.ndim != 2:
        raise ResultsError(f"ranks are (simulations, parameters); got shape {table.shape}.")
    fraction = table / float(draws)
    lower = ((1.0 - grid) / 2.0)[:, None, None]
    upper = ((1.0 + grid) / 2.0)[:, None, None]
    inside = (fraction[None, :, :] >= lower) & (fraction[None, :, :] <= upper)
    return grid, np.asarray(inside.mean(axis=1), dtype=float)


def uniformity_pvalues(ranks: Any, posterior_draws: int) -> np.ndarray:
    """One Kolmogorov—Smirnov p-value per parameter, ranks against ``Uniform(0, L)``.

    The same test ``sbi``'s ``check_sbc`` runs, spelled here so that route 2's
    number and route 1's mean the same thing and can sit in one variable. It is
    a frequentist p-value on a small sample: low is evidence of miscalibration,
    high at twenty simulations is evidence of very little.
    """
    table = np.asarray(ranks, dtype=float)
    if table.ndim != 2:
        raise ResultsError(f"ranks are (simulations, parameters); got shape {table.shape}.")
    reference = st.uniform(loc=0.0, scale=float(int(posterior_draws)))
    return np.asarray(
        [float(st.kstest(column, reference.cdf).pvalue) for column in table.T], dtype=float
    )


# ---------------------------------------------------------------------------
# The shared result
# ---------------------------------------------------------------------------


def calibration_dataset(
    ranks: Any,
    names: Sequence[str],
    *,
    posterior_draws: int,
    route: str,
    levels: Any = None,
    extras: Mapping[str, tuple[tuple[str, ...], Any]] | None = None,
    coords: Mapping[str, Any] | None = None,
    attrs: Mapping[str, Any] | None = None,
) -> Any:
    """Assemble the ``calibration`` group both routes return.

    One function so that the two routes cannot drift into two schemas. It
    computes the coverage curve from the ranks (:func:`coverage_from_ranks`)
    and, unless the caller supplies one in *extras*, the uniformity p-value
    too, then hangs whatever else the route measured beside them.

    Parameters
    ----------
    ranks
        ``(simulations, parameters)`` integer ranks in ``{0, ..., L}``.
    names
        One label per ranked column, in the ranks' own order. These are
        **merged parameter names** (``results.md`` §4), with an array-valued
        block expanded as ``name[element]``.
    posterior_draws
        ``L``.
    route
        :data:`SBI_ROUTE` or :data:`REFIT_ROUTE`.
    levels
        Nominal levels for the coverage curve; :data:`DEFAULT_LEVELS` by
        default.
    extras
        Extra variables, as ``{name: (dims, values)}`` — how route 1 carries
        ``c2st_ranks`` and TARP's curve without route 2 having to invent them.
    coords
        Extra coordinates the *extras* need, e.g. TARP's credibility grid.
    attrs
        Extra ``ampere_*`` attributes. The schema version, the route, the
        counts and the variable names are added here.

    Returns
    -------
    xarray.Dataset
        The ``calibration`` group. Every value is netCDF-safe.
    """
    xarray = _require_xarray()
    table = np.asarray(ranks)
    if table.ndim != 2:
        raise ResultsError(f"ranks are (simulations, parameters); got shape {table.shape}.")
    labels = [str(name) for name in names]
    if len(labels) != table.shape[1]:
        raise ResultsError(
            f"a calibration result needs one name per ranked column: {len(labels)} name(s) for "
            f"{table.shape[1]} column(s)."
        )
    if route not in (SBI_ROUTE, REFIT_ROUTE):
        raise ResultsError(
            f"a calibration route is {SBI_ROUTE!r} (an amortised posterior re-conditioned on "
            f"fresh simulations) or {REFIT_ROUTE!r} (a full fit per simulation); got {route!r}."
        )
    draws = int(posterior_draws)
    grid, coverage = coverage_from_ranks(table, draws, levels)

    supplied = dict(extras or {})
    variables: dict[str, Any] = {
        "ranks": ((SIMULATION_DIM, PARAMETER_DIM), table.astype(np.int64)),
        "coverage": ((LEVEL_DIM, PARAMETER_DIM), coverage),
    }
    if "ks_pvalue" not in supplied:
        variables["ks_pvalue"] = ((PARAMETER_DIM,), uniformity_pvalues(table, draws))
    for name, (dims, values) in supplied.items():
        variables[name] = (tuple(str(dim) for dim in dims), np.asarray(values))

    resolved: dict[str, Any] = {
        SIMULATION_DIM: np.arange(table.shape[0], dtype=np.int64),
        PARAMETER_DIM: np.asarray(labels),
        LEVEL_DIM: grid,
    }
    resolved.update({str(key): np.asarray(value) for key, value in (coords or {}).items()})

    recorded: dict[str, Any] = {
        f"{ATTR_PREFIX}calibration_schema_version": int(CALIBRATION_SCHEMA_VERSION),
        f"{ATTR_PREFIX}calibration_route": str(route),
        f"{ATTR_PREFIX}calibration_simulations": int(table.shape[0]),
        f"{ATTR_PREFIX}calibration_posterior_draws": draws,
    }
    recorded.update({str(key): _netcdf_safe(value) for key, value in (attrs or {}).items()})
    return xarray.Dataset(variables, coords=resolved, attrs=recorded)


def _netcdf_safe(value: Any) -> Any:
    """Attributes are netCDF scalars: an int, a float, or a string.

    ``bool`` is written as an ``int`` deliberately — netCDF has no boolean
    attribute type and h5netcdf would round-trip a ``True`` into a ``numpy``
    scalar whose ``repr`` is not ``True``, which is the kind of difference that
    makes a stored run and a live one disagree about a flag.
    """
    if isinstance(value, (bool, np.bool_)):
        return int(value)
    if isinstance(value, (int, np.integer)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value)
    return str(value)


def attach_calibration(tree: Any, calibration: Any) -> Any:
    """Hang a calibration result off a run as its ``calibration`` group.

    ``results.md`` §4's schema table gains this group with W3.6: a calibrated
    run archives the evidence that it is calibrated, so a reader deciding
    whether to believe an amortised posterior needs no second file.

    The tree is modified in place and returned, as attaching a child to a
    :class:`xarray.DataTree` is. An existing ``calibration`` group is replaced:
    a run holds one calibration, the most recent, and silently keeping two
    under one name is not a thing xarray can do anyway.
    """
    xarray = _require_xarray()
    tree[CALIBRATION_GROUP] = xarray.DataTree(calibration)
    return tree


# ---------------------------------------------------------------------------
# Route 2: the Talts et al. loop over any engine
# ---------------------------------------------------------------------------


def replace_observations(
    problem: FittingProblem,
    observations: Mapping[str, FunctionSamples],
    *,
    seed: int | None = None,
) -> FittingProblem:
    """*problem* with each dataset's observed container swapped for a simulated one.

    The replica an SBC iteration is fitted on. Everything else is carried over
    unchanged — the models, the frozen instrument chains, the likelihoods, the
    ties, the labels and the model bindings — so the only difference between
    the replica and the original is the data, which is exactly the difference
    SBC is about. The containers ``simulate(observe=True)`` produces come from
    :meth:`~ampere.core.dataset.Dataset.draw_observation`, which builds them
    with :meth:`~ampere.core.results_schema.FunctionSamples.with_values`, so
    they already carry the original's coordinates, uncertainties and mask.

    *seed* is the replica's own run seed. :func:`sbc` derives one per iteration
    off its calibration sub-stream, so two iterations do not initialise their
    walkers identically; ``None`` leaves the replica unseeded.
    """
    from ampere.core.dataset import Dataset, DatasetCollection

    missing = [label for label in problem.datasets if label not in observations]
    if missing:
        raise ResultsError(
            f"a replica needs one observed container per dataset; {sorted(missing)} "
            f"{'is' if len(missing) == 1 else 'are'} missing. A simulation that failed has no "
            f"observations at all — check Simulation.failed before building a replica."
        )
    rebuilt: dict[str, Any] = {}
    for label, dataset in problem.datasets.items():
        latent = dataset.latent
        extra = {} if latent is None else {"latent_name": latent.parameter.name}
        rebuilt[label] = Dataset(
            observations[label],
            dataset.instrument,
            dataset.likelihood,
            label=dataset.label,
            model=dataset.model,
            meta=dataset.meta,
            **extra,
        )
    return FittingProblem(
        problem.models,
        DatasetCollection(rebuilt),
        ties=problem.ties,
        seed=seed,
        strict=problem.strict,
    )


def sbc(
    problem: FittingProblem,
    engine_factory: Callable[[FittingProblem], Any],
    *,
    count: int,
    draws: int,
    run_options: Mapping[str, Any] | None = None,
    parameters: Mapping[str, str] | Sequence[str] | None = None,
    levels: Any = None,
    seed: int | None = None,
    label: str | None = None,
) -> Any:
    """Simulation-based calibration by refitting: Talts et al. (2018), any engine.

    The honest, expensive route, and the only one that can validate a
    *likelihood* rather than a network. For each of *count* iterations: draw θ
    from *problem*'s joint prior, simulate a dataset there, build a replica
    problem carrying that data, hand it to *engine_factory*, run the fit, and
    rank the true θ among *draws* of the resulting posterior. The cost is one
    full fit per iteration, which is why every budget here is the caller's to
    choose and none of them has a default.

    Parameters
    ----------
    problem
        The **simulating** problem: what θ is drawn from and what generates the
        data. It need not be the problem that is fitted — pointing the two at
        different formulations of the same experiment is how a *statistic* is
        put on trial, which is what ``examples/wstat_comparison.py`` does with
        the profiled Cash-with-background statistic.
    engine_factory
        Called once per iteration with the replica
        :class:`~ampere.core.dataset.FittingProblem` (the same problem carrying
        the simulated data) and returning either an engine — anything with a
        ``run`` method, so that :mod:`ampere.results` still imports no engine —
        or an already-emitted run. A factory that wants to fit something *else*
        reads the simulated containers off the replica's ``datasets`` and
        builds its own problem, which is the WStat case.
    count
        How many simulations. Talts et al.'s uniformity test is a
        goodness-of-fit test on this many points, so a hundred or more is where
        it starts to have power; a handful is a smoke test and is reported as
        one (a warning below :data:`_SBC_WARNING_FLOOR`, matching ``sbi``'s own).
    draws
        ``L``: how many posterior draws each rank is taken against. Thinned
        evenly out of each chain the fit produced, which is what breaks the
        autocorrelation an MCMC run's neighbouring draws carry; a fit that
        produced fewer than this is refused rather than silently ranked against
        a shorter list, because ``L`` is what the rank's uniform null is
        defined against.
    run_options
        Forwarded verbatim to the engine's ``run``. This is where an emcee
        fit's ``steps`` and ``burn_in`` go. Ignored when the factory returns a
        run rather than an engine.
    parameters
        Which columns to rank. ``None`` ranks every variable of the fit's
        ``posterior`` group, matching truth by the same merged name. A sequence
        names a subset of them. A mapping is ``{posterior name: truth name}``,
        for the case where the fitted problem spells a parameter differently
        from the simulating one.
    levels
        Nominal levels for the coverage curve; :data:`DEFAULT_LEVELS` by
        default.
    seed
        Base seed for the simulations and for each replica's fit. *problem*'s
        own seed by default, so a seeded problem gives a reproducible study;
        an unseeded one gets entropy, and the seed actually used is recorded on
        the result.
    label
        A short name for what was calibrated, recorded in
        ``ampere_calibration_label``. Free text for the figure's title.

    Returns
    -------
    xarray.Dataset
        The ``calibration`` group: ``ranks``, ``coverage``, ``ks_pvalue``, and
        the provenance of the study on its attributes.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If no iteration produced a usable fit, or if a fit produced fewer than
        *draws* draws, or if a requested parameter is not in the fit's
        posterior.
    """
    simulations = int(count)
    wanted = int(draws)
    if simulations < 1:
        raise ResultsError(
            f"simulation-based calibration needs at least one simulation, got {count}."
        )
    if wanted < 1:
        raise ResultsError(f"a rank needs at least one posterior draw, got {draws}.")
    if simulations < _SBC_WARNING_FLOOR:
        warnings.warn(
            f"simulation-based calibration over {simulations} simulation(s): the uniformity "
            f"test is a goodness-of-fit test on that many points, so at this budget a passing "
            f"result is weak evidence of calibration and only a gross failure will show. "
            f"{_SBC_WARNING_FLOOR} or more is where it starts to have power.",
            UserWarning,
            stacklevel=2,
        )

    base = _base_seed(problem, seed)
    rng = generator(base, f"{CALIBRATION_STREAM}.simulate")
    options = dict(run_options or {})

    rank_rows: list[np.ndarray] = []
    names: list[str] | None = None
    failures = 0
    engine_names: set[str] = set()
    for index in range(simulations):
        simulation = problem.simulate(observe=True, rng=rng)
        if simulation.failed or simulation.observations is None:
            failures += 1
            continue
        replica = replace_observations(
            problem,
            simulation.observations,
            seed=int(substream(base, f"{CALIBRATION_STREAM}.fit.{index}")),
        )
        candidate = engine_factory(replica)
        engine_names.add(type(candidate).__name__)
        run = candidate if _is_run(candidate) else candidate.run(**options)
        columns = _posterior_columns(run, parameters)
        if names is None:
            names = [column.label for column in columns]
        elif [column.label for column in columns] != names:
            raise ResultsError(
                f"iteration {index} ranked {[column.label for column in columns]} but the first "
                f"ranked {names}. Every iteration must produce the same posterior variables, or "
                f"the ranks are not one statistic."
            )
        truth = _truth_vector(simulation.parameters, columns, index)
        rank_rows.append(_ranks_of(columns, truth, wanted, index))

    if not rank_rows or names is None:
        raise ResultsError(
            f"simulation-based calibration produced no usable iteration: all {simulations} "
            f"simulation(s) failed. problem.failure_summary() says why:\n"
            f"{problem.failure_summary()}"
        )
    return calibration_dataset(
        np.stack(rank_rows),
        names,
        posterior_draws=wanted,
        route=REFIT_ROUTE,
        levels=levels,
        attrs={
            f"{ATTR_PREFIX}calibration_seed": base,
            f"{ATTR_PREFIX}calibration_requested": simulations,
            f"{ATTR_PREFIX}calibration_failures": failures,
            f"{ATTR_PREFIX}calibration_engine": ", ".join(sorted(engine_names)) or "unknown",
            f"{ATTR_PREFIX}calibration_uniformity_check": "scipy.stats.kstest",
            **({} if label is None else {f"{ATTR_PREFIX}calibration_label": str(label)}),
        },
    )


def _is_run(candidate: Any) -> bool:
    """Whether *candidate* is already an emitted run rather than an engine.

    Asked of the object's *shape*, not its type, because :mod:`ampere.results`
    must not import :class:`~ampere.inference.engine.Engine` -- ``ampere.
    inference`` imports this namespace, so the dependency only runs one way. A
    run is a tree with a ``posterior`` child; anything else is asked to
    ``run()``.
    """
    children = getattr(candidate, "children", None)
    return children is not None and "posterior" in children


def _base_seed(problem: FittingProblem, seed: int | None) -> int:
    """The seed every stream here derives from, recorded on the result.

    An explicit argument wins; otherwise the problem's own, so a seeded problem
    gives a reproducible study. An unseeded problem gets entropy — the honest
    behaviour for one that did not ask to be reproducible — and the number is
    stored so the study can be repeated after the fact anyway. The same rule,
    and the same reasoning, as
    :func:`ampere.results.diagnostics.residual_whiteness`'s.
    """
    if seed is not None:
        return int(seed)
    if problem.seed is not None:
        return int(problem.seed)
    return int((np.random.SeedSequence().entropy or 0) % (2**31))


class _Column:
    """One scalar posterior column: its label, its draws, and where truth lives.

    ``results.md`` §4 stores an array-valued parameter as **one** variable with
    a named dimension, never ``N`` scalar names, so a rank per element has to
    expand the block here — and has to remember which element it was, because
    the truth is one entry of the simulating problem's own array.
    """

    __slots__ = ("draws", "element", "label", "source")

    def __init__(self, label: str, source: str, element: int | None, draws: np.ndarray) -> None:
        self.label = label
        self.source = source
        self.element = element
        self.draws = draws


def _posterior_columns(run: Any, parameters: Any) -> list[_Column]:
    """The fit's posterior, as scalar columns in the caller's chosen order."""
    children = getattr(run, "children", None)
    if children is None or "posterior" not in children:
        raise ResultsError(
            "an engine in a calibration loop must return an emitted run with a 'posterior' "
            "group (ampere.results.emit's output). Got an object without one — a factory that "
            "returns a run rather than an engine must return the tree, not the draws."
        )
    group = run["posterior"].dataset
    available = [str(name) for name in group.data_vars]
    if parameters is None:
        wanted = {name: name for name in available}
    elif isinstance(parameters, Mapping):
        wanted = {str(key): str(value) for key, value in parameters.items()}
    else:
        wanted = {str(name): str(name) for name in parameters}
    unknown = [name for name in wanted if name not in available]
    if unknown:
        raise ResultsError(
            f"the fit's posterior has no variable named {sorted(unknown)}; it has "
            f"{sorted(available)}. Variables are keyed by merged parameter name "
            f"(results.md §4), so 'model.index' and 'index' are different names."
        )
    columns: list[_Column] = []
    for name, source in wanted.items():
        data = group[name]
        extra = [str(dim) for dim in data.dims if dim not in (CHAIN_DIM, DRAW_DIM)]
        values = np.asarray(data.values, dtype=float)
        if not extra:
            columns.append(_Column(name, source, None, values))
            continue
        block = values.reshape(values.shape[0], values.shape[1], -1)
        for element, tag in enumerate(_element_labels(data, extra, block.shape[-1])):
            columns.append(_Column(f"{name}[{tag}]", source, element, block[:, :, element]))
    if not columns:
        raise ResultsError("no posterior variable was selected, so there is nothing to rank.")
    return columns


def _element_labels(data: Any, extra: Sequence[str], size: int) -> list[str]:
    """One label per element of an array-valued block, coordinates preferred.

    The same rule :func:`ampere.results._plotting.parameter_columns` uses, so a
    rank histogram's element labels match a corner plot's.
    """
    if len(extra) == 1 and extra[0] in data.coords:
        return [str(value) for value in np.asarray(data.coords[extra[0]].values).ravel()]
    shape = tuple(int(data.sizes[dim]) for dim in extra)
    if len(shape) == 1:
        return [str(index) for index in range(size)]
    return [", ".join(str(part) for part in index) for index in np.ndindex(*shape)]


def _truth_vector(values: Mapping[str, Any], columns: Sequence[_Column], index: int) -> np.ndarray:
    """The simulated θ, one entry per ranked column."""
    truth = np.empty(len(columns), dtype=float)
    for position, column in enumerate(columns):
        if column.source not in values:
            raise ResultsError(
                f"iteration {index}: the simulating problem has no parameter named "
                f"{column.source!r}, so column {column.label!r} has no truth to rank against. "
                f"It has {sorted(values)}. Where the fitted problem spells a parameter "
                f"differently from the simulating one, pass parameters={{posterior name: truth "
                f"name}}."
            )
        found = np.asarray(values[column.source], dtype=float).reshape(-1)
        if column.element is None:
            if found.size != 1:
                raise ResultsError(
                    f"iteration {index}: column {column.label!r} is scalar in the posterior but "
                    f"the simulating problem's {column.source!r} has {found.size} elements."
                )
            truth[position] = float(found[0])
            continue
        if column.element >= found.size:
            raise ResultsError(
                f"iteration {index}: column {column.label!r} indexes element "
                f"{column.element} of the simulating problem's {column.source!r}, which has "
                f"{found.size}."
            )
        truth[position] = float(found[column.element])
    return truth


def _ranks_of(columns: Sequence[_Column], truth: np.ndarray, draws: int, index: int) -> np.ndarray:
    """One rank per column: how many of *draws* thinned posterior draws fall below truth."""
    ranks = np.empty(len(columns), dtype=np.int64)
    for position, column in enumerate(columns):
        thinned = _thin(column.draws, draws, column.label, index)
        ranks[position] = int(np.count_nonzero(thinned < truth[position]))
    return ranks


def _thin(values: np.ndarray, draws: int, label: str, index: int) -> np.ndarray:
    """*draws* draws, thinned evenly **within** each chain and then concatenated.

    Evenly rather than at random, and per chain rather than across the pooled
    array, because the autocorrelation a rank must not inherit is a
    within-chain property: taking every ``m``-th draw of each chain is the
    standard remedy Talts et al. prescribe, while pooling first and thinning
    afterwards would take consecutive draws whenever the chain count does not
    divide the stride.
    """
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        array = array.reshape(1, -1)
    chains, length = array.shape
    per_chain = -(-draws // chains)
    if length < per_chain:
        raise ResultsError(
            f"iteration {index}: ranking {label!r} against {draws} draw(s) needs {per_chain} "
            f"from each of the fit's {chains} chain(s), and each holds {length}. Ask the fit "
            f"for more draws, or lower draws=; a rank's uniform null is defined against a "
            f"fixed number of posterior draws, so it is not shortened silently."
        )
    picked = np.linspace(0, length - 1, per_chain).round().astype(int)
    # Draw-major, so that truncating ``chains * per_chain`` back to ``draws``
    # drops the last thinned draw of a few chains rather than the whole tail of
    # the last chain -- an imbalance that would weight one chain's mixing more
    # heavily than another's for no reason.
    return np.asarray(array[:, picked].ravel(order="F")[:draws], dtype=float)
