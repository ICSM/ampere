"""Groups a run does **not** store by default, and the rule for computing them.

``diagnostics.md`` §7 puts two questions to W1.8, and this module is where the
answers are made concrete rather than left as prose.

**(1) Are posterior-predictive replicates stored by default?** *No.* The cost is
``N_draws * N_obs`` per dataset per run — for a 10⁵-point spectrum and 10⁴ draws
that is 8 GB of float64 in a file whose whole point is being cheap enough to
archive, and family B needs them only when a user actually asks for a check.
They are computed on demand by :func:`add_posterior_predictive`, into the
reserved group :data:`POSTERIOR_PREDICTIVE_GROUP`, from the stored draws and the
problem — which is exactly the re-evaluation ``inference.md`` §13's ``simulate``
already performs. The group name is reserved *now* so that a run which does
store them (a small dataset, a user who asked) is readable by the same code.

**(2) Where do signed per-point residuals live?** In the reserved group
:data:`RESIDUALS_GROUP`, also on demand. ``diagnostics.md`` §3.3 is precise about
why this is not free: §4.6's log-likelihood requirement makes the
posterior-predictive half cheap, but a Ljung-Box whiteness test needs the
**sign**, and for a Gaussian noise model ``log_likelihood_i`` recovers
``|residual_i|`` and not its sign. The derivation is stated here so that no
implementation has to invent it: the standardised residual of dataset *d* at
draw *k* is

.. code-block:: text

    r_dk = (observed_d - predicted_d(theta_k)) / sigma_d

with ``predicted_d(theta_k)`` the instrument chain's output — ``Dataset.predict``
on the model result at ``theta_k`` — and ``sigma_d`` the observed container's own
``masked_uncertainty()``. Masked samples are carried as NaN rather than dropped,
so the coordinate axis stays the container's own and family B's statistic can
decide for itself what to do with a gap.

**Family C's input is not a group at all.** ``likelihoods.md`` §16 fixes it:
:meth:`ampere.core.likelihood.Likelihood.conditional` returns a signed mean and
a variance *on the full coordinate axis, masked samples included* — because
"what would the GP have said here?" is exactly the question asked about an
excluded region. :func:`gp_localisation` evaluates it across posterior draws;
the caveat that must travel with the answer is in :mod:`ampere.results.plots`.

W1.8 declared all three surfaces without implementing them. **W2.7 lands two of
them** — :func:`add_residuals` and :func:`gp_localisation`, the inputs
``diagnostics.md``'s families B and C actually consume — leaving
:func:`add_posterior_predictive` for the item that lands
``plot_posterior_predictive`` beside it. Nothing about the shapes, the group
names or the cost policy has moved: they are what W1.8 fixed, which is the
point of having fixed them before two backend tracks and three diagnostic
families wrote against them.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ampere.core.dataset import LIKELIHOOD_COMPONENT, Dataset, FittingProblem
from ampere.core.exceptions import ResultsError
from ampere.core.likelihood import GaussianProcessNoise
from ampere.core.results_schema import FunctionSamples

from .emission import (
    CHAIN_DIM,
    DRAW_DIM,
    POSTERIOR_GROUP,
    _container_dims,
    _require_arviz,
)
from .provenance import ATTR_PREFIX, hash_of, problem_fingerprint

__all__ = [
    "GP_LOCALISATION_GROUP",
    "POINTWISE_LOG_LIKELIHOOD_GROUP",
    "POSTERIOR_PREDICTIVE_GROUP",
    "RESIDUALS_GROUP",
    "add_posterior_predictive",
    "add_residuals",
    "gp_localisation",
]

#: Replicate observations ``y_rep ~ p(y | θ_k)``, one variable per dataset,
#: dims ``(chain, draw, <dataset dims>)``. Not written by default.
POSTERIOR_PREDICTIVE_GROUP = "posterior_predictive"

#: Signed standardised residuals, one variable per dataset, dims
#: ``(chain, draw, <dataset dims>)``. Not written by default.
RESIDUALS_GROUP = "residuals"

#: The conditioned GP mean and variance on the full coordinate axis, one pair of
#: variables per dataset (``<label>_mean``, ``<label>_variance``). Not written by
#: default; family C's input.
GP_LOCALISATION_GROUP = "gp_localisation"

#: ArviZ's per-**observation** ``log_likelihood`` convention, reserved but not
#: written: computing it needs a method on ``Likelihood`` that does not exist
#: (``results.md`` §6 and its ruling request R2). It is a distinct group from
#: ``log_likelihood`` precisely so the two decompositions cannot be confused —
#: a variable here carries an ``ampere_decomposition`` of ``"factorised"``
#: (independent noise; exact) or ``"conditional_loo"`` (a GP's leave-one-out
#: conditional terms). Reserved now so a run emitted today is
#: forward-compatible with one emitted after it lands.
POINTWISE_LOG_LIKELIHOOD_GROUP = "pointwise_log_likelihood"

_PHASE_2 = (
    "Its shape, group name and cost policy are fixed by W1.8 "
    "(docs/design/contracts/results.md §7); the computation lands in Phase 2, "
    "which is when a backend exists to make it cheap."
)


def add_posterior_predictive(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    seed_stream: str = "posterior_predictive",
) -> Any:
    """Draw ``y_rep`` for the stored posterior and attach them to ``tree``.

    Uses :meth:`~ampere.core.dataset.FittingProblem.simulate` with
    ``observe=True`` at each retained draw, so the replicates come from the same
    ``LikelihoodFamily.sample`` the likelihood scores with — never from a
    Gaussian assumption bolted on here (``inference.md`` §13).

    Parameters
    ----------
    tree
        The run, as :func:`ampere.results.emit` produced it.
    problem
        The problem the run was over. Its ``problem_hash`` must match the
        stored one, or the replicates would come from a different model.
    datasets
        Which datasets to replicate; all of them by default.
    thin
        Take every ``thin``-th draw. The reason this parameter exists is the
        memory cost that keeps the group out of the default emission.
    seed_stream
        The named RNG sub-stream (``lowering.md`` §9.2) the draws come from.
        Separate from ``"simulate"`` on purpose: adding a predictive check must
        not change an SBI budget's draws.
    """
    raise NotImplementedError(
        f"posterior-predictive replicates are not implemented yet. {_PHASE_2}"
    )


def add_residuals(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    standardised: bool = True,
) -> Any:
    """Compute signed per-point residuals for the stored posterior.

    ``standardised=True`` divides by the observed container's own uncertainty,
    which is what family B's whiteness statistic wants; ``False`` keeps them in
    the data's units, which is what a residual panel wants. Masked samples are
    NaN, never dropped.

    The derivation is ``results.md`` §7's, not this function's invention:

    .. code-block:: text

        r_dk = (observed_d - predicted_d(theta_k)) / sigma_d

    with ``predicted_d(theta_k)`` the instrument chain's output at draw ``k``
    and ``sigma_d`` the observed container's own ``masked_uncertainty()``. The
    prediction comes from :meth:`~ampere.core.dataset.FittingProblem.simulate`
    with ``observe=False``, which is the public path to "run the forward model
    once at this theta" and which **flags** rather than raises when a draw
    cannot be evaluated — so a run holding a few prior-rejected or
    model-failing draws yields NaN rows for those and residuals for the rest,
    rather than no residuals at all.

    Parameters
    ----------
    tree
        The run, as :func:`ampere.results.emit` produced it. Returned with the
        ``residuals`` group attached; the tree is modified in place, as
        attaching a child to a :class:`xarray.DataTree` is.
    problem
        The problem the run was over. Its ``problem_hash`` must match the
        stored one, or the residuals would be against a different model.
    datasets
        Which datasets to residualise; all of them by default.
    thin
        Take every ``thin``-th draw. The reason this parameter exists is the
        ``N_draws x N_obs`` memory cost that keeps the group out of the default
        emission; the retained draw indices become the group's own ``draw``
        coordinate, so a thinned group still says which draws it came from.
    standardised
        Divide by ``sigma``. Recorded on the group, because a residual panel
        and a whiteness statistic want different answers and neither should
        have to guess which it was handed.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the problem is not the one the run was over, if a requested dataset
        is not in it, or if ``standardised=True`` and a dataset has no
        uncertainties to standardise by.
    """
    _require_same_problem(tree, problem)
    labels = _requested(problem, datasets)
    if standardised:
        for label in labels:
            observed = problem.datasets[label].observed
            if observed.uncertainty is None:
                raise ResultsError(
                    f"dataset {label!r} has no uncertainties, so a standardised residual is "
                    f"undefined. Pass standardised=False for residuals in the data's own "
                    f"units, or attach uncertainties to the observed container."
                )
    thetas, kept = _stored_thetas(tree, problem, thin)
    chains, draws = thetas.shape[0], thetas.shape[1]
    variables = {
        label: np.full(
            (chains, draws, problem.datasets[label].observed.n_samples), np.nan, dtype=float
        )
        for label in labels
    }
    for chain in range(chains):
        for draw in range(draws):
            simulation = problem.simulate(thetas[chain, draw])
            if simulation.failed:
                continue
            for label in labels:
                dataset = problem.datasets[label]
                predicted = simulation.predicted.get(label)
                if predicted is None:  # pragma: no cover - a failure short-circuits above
                    continue
                variables[label][chain, draw] = _residual_of(dataset, predicted, standardised)
    dims, coords = _observed_axes(problem, labels)
    return _attach(
        tree,
        RESIDUALS_GROUP,
        variables,
        dims,
        coords,
        kept,
        {
            f"{ATTR_PREFIX}standardised": int(bool(standardised)),
            f"{ATTR_PREFIX}derivation": (
                "(observed - predicted(theta)) / sigma, per retained draw; masked samples are "
                "NaN, never dropped (results.md §7)."
                if standardised
                else "observed - predicted(theta), per retained draw, in the data's own units; "
                "masked samples are NaN, never dropped (results.md §7)."
            ),
        },
    )


def _residual_of(dataset: Dataset, predicted: FunctionSamples, standardised: bool) -> np.ndarray:
    """One draw's signed residual for one dataset, masked samples NaN.

    ``masked_uncertainty()`` inflates a masked sample's sigma to ``inf``, which
    would send its residual to a perfectly white **zero** — the one value a
    whiteness statistic must not be handed for a sample nobody measured. The
    mask is therefore applied as NaN afterwards, which is what ``results.md``
    §7 asks for and what lets the statistic "decide for itself what to do with
    a gap".
    """
    observed = dataset.observed
    residual = (
        np.asarray(observed.values, dtype=float).ravel()
        - np.asarray(predicted.values, dtype=float).ravel()
    )
    if standardised:
        residual = residual / np.asarray(observed.masked_uncertainty(), dtype=float).ravel()
    excluded = dataset.effective_mask
    if excluded is not None:
        residual = np.where(np.asarray(excluded).ravel(), np.nan, residual)
    return residual


def gp_localisation(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    at: Any = None,
) -> Any:
    """Evaluate the conditioned GP mean and variance across posterior draws.

    Delegates per draw to :meth:`ampere.core.likelihood.Likelihood.conditional`,
    whose defaults are already the right ones for this diagnostic: the **full**
    coordinate axis, masked samples included, and a **signed** mean, so the
    output shows the direction of the local deficiency and not only its size.
    ``at`` overrides the evaluation grid, for a smoother curve than the data's
    own sampling.

    Only datasets whose likelihood carries a
    :class:`~ampere.core.likelihood.GaussianProcessNoise` can be localised this
    way; ``conditional`` refuses the rest by name, and this function propagates
    that refusal rather than silently skipping them. Concretely: a *named*
    dataset without a GP raises, because asking to localise a dataset that has
    no GP is a mistake worth reporting; ``datasets=None`` means "every dataset
    that has one", and refuses only when the problem has none at all, because
    the alternative would make the default unusable on the mixed
    GP/non-GP joint fit this whole namespace exists to serve.

    Parameters
    ----------
    tree
        The run. Returned with the ``gp_localisation`` group attached, holding
        ``<label>_mean`` and ``<label>_variance`` per dataset.
    problem
        The problem the run was over; its ``problem_hash`` must match.
    datasets
        Which datasets to localise. See above for what ``None`` means.
    thin
        Take every ``thin``-th draw; the retained indices become the group's
        ``draw`` coordinate.
    at
        A 1-D array of coordinates to evaluate on instead of the data's own
        axis. Handed to the solver untouched, so one place decides what a bare
        array of coordinates means.
    """
    _require_same_problem(tree, problem)
    labels = _gp_requested(problem, datasets)
    thetas, kept = _stored_thetas(tree, problem, thin)
    chains, draws = thetas.shape[0], thetas.shape[1]
    grid = None if at is None else np.asarray(at, dtype=float)
    size = {
        label: (problem.datasets[label].observed.n_samples if grid is None else int(grid.shape[0]))
        for label in labels
    }
    variables: dict[str, np.ndarray] = {}
    for label in labels:
        for role in ("mean", "variance"):
            variables[f"{label}_{role}"] = np.full(
                (chains, draws, size[label]), np.nan, dtype=float
            )
    for chain in range(chains):
        for draw in range(draws):
            simulation = problem.simulate(thetas[chain, draw])
            if simulation.failed:
                continue
            routed = problem.mapping.distribute(simulation.parameters)
            for label in labels:
                dataset = problem.datasets[label]
                predicted = simulation.predicted.get(label)
                if predicted is None:  # pragma: no cover - a failure short-circuits above
                    continue
                split = dataset.route(routed.get(label, {}))
                conditional = dataset.likelihood.conditional(
                    predicted,
                    dataset.observed,
                    split.get(LIKELIHOOD_COMPONENT),
                    at=grid,
                )
                variables[f"{label}_mean"][chain, draw] = np.asarray(conditional.mean).ravel()
                variables[f"{label}_variance"][chain, draw] = np.asarray(
                    conditional.variance
                ).ravel()
    if grid is None:
        dims, coords = _observed_axes(problem, labels, roles=("mean", "variance"))
    else:
        dims = {}
        coords = {}
        for label in labels:
            name = f"{label}_gp_axis"
            coords[name] = grid
            dims[f"{label}_mean"] = [name]
            dims[f"{label}_variance"] = [name]
    return _attach(
        tree,
        GP_LOCALISATION_GROUP,
        variables,
        dims,
        coords,
        kept,
        {
            f"{ATTR_PREFIX}evaluation_grid": "data" if grid is None else "explicit",
            f"{ATTR_PREFIX}interpretation_note": (
                "The conditioned GP mean is signed and localises where the model is deficient; "
                "it does not say why. See ampere.results.gp_localisation_caveat()."
            ),
        },
    )


# ---------------------------------------------------------------------------
# Shared machinery: reading a stored run back into thetas
# ---------------------------------------------------------------------------


def _require_same_problem(tree: Any, problem: FittingProblem) -> None:
    """Refuse a problem that is not the one the run was over.

    Recomputing a residual against a different model would produce a group that
    is wrong in a way nothing downstream could detect — the coordinates would
    line up and the numbers would be meaningless. ``ampere_problem_hash``
    covers the models, the data, the instruments and the likelihoods, which is
    exactly the set a derived group depends on.
    """
    stored = getattr(tree, "attrs", {}).get(f"{ATTR_PREFIX}problem_hash")
    if stored is None:
        raise ResultsError(
            "this tree carries no ampere provenance, so there is no way to check that it is a "
            "run over the problem given. Derived groups are computed against the problem, and "
            "computing them against the wrong one is undetectable downstream."
        )
    current = hash_of(problem_fingerprint(problem))
    if stored != current:
        raise ResultsError(
            f"this run was over a different problem: it records problem_hash {stored!r} and the "
            f"problem given hashes to {current!r}. The residual and localisation groups are "
            f"derived by re-running the forward model at the stored draws, so the problem has "
            f"to be the one the draws came from."
        )


def _requested(problem: FittingProblem, datasets: Sequence[str] | None) -> tuple[str, ...]:
    """Which dataset labels to compute for, validated against the problem."""
    if datasets is None:
        return tuple(problem.datasets)
    unknown = [label for label in datasets if label not in problem.datasets]
    if unknown:
        raise ResultsError(
            f"this problem has no dataset(s) {sorted(unknown)}; it has {sorted(problem.datasets)}."
        )
    if not tuple(datasets):
        raise ResultsError("no datasets were requested, so there is nothing to compute.")
    return tuple(datasets)


def _gp_requested(problem: FittingProblem, datasets: Sequence[str] | None) -> tuple[str, ...]:
    """Family C's applicable datasets — see :func:`gp_localisation`'s docstring."""
    if datasets is None:
        chosen = tuple(
            label
            for label in problem.datasets
            if isinstance(problem.datasets[label].likelihood.noise, GaussianProcessNoise)
        )
        if not chosen:
            raise ResultsError(
                "no dataset in this problem was fitted with a GaussianProcessNoise model, so "
                "there is no conditioned GP mean to localise with. Residual whiteness "
                "(ampere.results.add_residuals, then residual_whiteness) is the diagnostic for "
                "a standard-likelihood fit — diagnostics.md family B, which is precisely the "
                "test for whether the flexible likelihood is worth switching on."
            )
        return chosen
    return _requested(problem, datasets)


def _stored_thetas(tree: Any, problem: FittingProblem, thin: int) -> tuple[np.ndarray, np.ndarray]:
    """``(thetas, draw_indices)`` — the posterior group packed back to free vectors.

    The posterior is stored by merged parameter name, which is the
    representation every consumer wants and the one thing a re-evaluation does
    not: :meth:`~ampere.core.parameter.ParameterSet.pack` puts it back into the
    flat order ``simulate`` takes. An array-valued parameter is one variable
    with a named dimension (``results.md`` §4), so it packs straight back
    without ever materialising ``free_labels()``.
    """
    if thin < 1:
        raise ResultsError(f"thin must be at least 1, got {thin}.")
    posterior = getattr(tree, "children", {})
    if POSTERIOR_GROUP not in posterior:
        raise ResultsError(
            "this tree has no 'posterior' group, so there are no draws to re-evaluate at."
        )
    group = tree[POSTERIOR_GROUP].dataset
    names = [str(name) for name in group.data_vars]
    chains = int(group.sizes[CHAIN_DIM])
    count = int(group.sizes[DRAW_DIM])
    kept = np.arange(0, count, thin)
    arrays = {name: np.asarray(group[name].values) for name in names}
    thetas = np.empty((chains, kept.size, problem.free_size), dtype=float)
    for chain in range(chains):
        for position, draw in enumerate(kept):
            values = {name: arrays[name][chain, draw] for name in names}
            thetas[chain, position] = problem.parameters.pack(values)
    return thetas, kept


def _observed_axes(
    problem: FittingProblem,
    labels: Sequence[str],
    roles: Sequence[str] = (),
) -> tuple[dict[str, list[str]], dict[str, Any]]:
    """Dimension names and coordinates matching the run's own ``observed_data``.

    Reuses :func:`ampere.results.emission._container_dims` rather than
    re-deriving the naming, so a derived group's x-axis is the *same*
    coordinate as the observed data it is a residual of — which is what lets a
    plot draw them on one axes without an alignment step.
    """
    dims: dict[str, list[str]] = {}
    coords: dict[str, Any] = {}
    for label in labels:
        observed = problem.datasets[label].observed
        names = list(_container_dims(label, observed))
        for axis, name in zip(observed.axes, names, strict=False):
            coords[name] = np.asarray(axis.values, dtype=float).ravel()
        for role in roles or ("",):
            dims[f"{label}_{role}" if role else label] = names
    return dims, coords


def _attach(
    tree: Any,
    group: str,
    variables: Mapping[str, np.ndarray],
    dims: Mapping[str, list[str]],
    coords: Mapping[str, Any],
    kept: np.ndarray,
    attrs: Mapping[str, Any],
) -> Any:
    """Build one derived group and hang it off the run, coordinates and all.

    The retained draw indices become the group's own ``draw`` coordinate, so a
    thinned group states which draws it holds rather than renumbering them from
    zero and quietly losing the correspondence to the posterior.
    """
    arviz = _require_arviz()
    resolved = dict(coords)
    resolved[DRAW_DIM] = kept
    built = arviz.from_dict({group: dict(variables)}, dims=dict(dims), coords=resolved)
    for name, child in built.children.items():
        tree[name] = child
    tree[group].attrs.update(dict(attrs))
    return tree
