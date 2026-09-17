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

W1.8 declared all three surfaces without implementing them. **W2.7 landed two of
them** — :func:`add_residuals` and :func:`gp_localisation`, the inputs
``diagnostics.md``'s families B and C actually consume — leaving
:func:`add_posterior_predictive` for the item that lands
``plot_posterior_predictive`` beside it. **W2.8 lands that one**, and with it
the fourth on-demand group: :func:`add_pointwise_log_likelihood`, ArviZ's
per-*observation* convention over
:meth:`ampere.core.likelihood.Likelihood.pointwise_log_prob`. Nothing about the
shapes, the group names or the cost policy has moved: they are what W1.8 fixed,
which is the point of having fixed them before two backend tracks and three
diagnostic families wrote against them.

The per-observation group, and why it is still an explicit call
---------------------------------------------------------------
``results.md`` §6 reserved :data:`POINTWISE_LOG_LIKELIHOOD_GROUP` and its
``ampere_decomposition`` attribute — ``"factorised"`` for independent noise,
``"conditional_loo"`` for a GP — and ruled (2026-09-03, §15 R2) that the terms
are "granted at the freeze but **not stored by default**". The three reasons
are unchanged by having a computation available: the group is
``N_draws x N_obs`` per dataset where the per-dataset one is
``N_draws x N_datasets``; a run that does not need LOO should not carry it;
and two decompositions must not share a group name, or ``arviz.loo`` would
mean three different things depending on the fit. So it is written by
:func:`add_pointwise_log_likelihood` and by nothing else.

Where the solver cannot supply the decomposition, this module **refuses by
name and never falls back to a dense solve**:
:meth:`ampere.core.likelihood.QuasisepGP.conditional_loo` is deferred (W2.3),
as is the jax one, and silently substituting :class:`DenseGP` would turn an
``O(N)`` fit's diagnostic into an ``O(N^3)`` one behind the user's back —
which for the 20 000-point spectra the quasiseparable solve exists to serve is
not a slower answer but a different program.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from ampere.core.dataset import (
    LIKELIHOOD_COMPONENT,
    Dataset,
    FittingProblem,
    sample_coordinates,
)
from ampere.core.exceptions import DatasetError, LikelihoodError, ResultsError
from ampere.core.likelihood import GaussianProcessNoise, Marginalisation
from ampere.core.results_schema import FunctionSamples, format_axis_label

from . import _plotting as _p
from .emission import (
    CHAIN_DIM,
    DRAW_DIM,
    LOG_LIKELIHOOD_GROUP,
    POSTERIOR_GROUP,
    _container_dims,
    _require_arviz,
)
from .provenance import ATTR_PREFIX, hash_of, problem_fingerprint

__all__ = [
    "COMPONENTS",
    "CONDITIONAL_LOO_DECOMPOSITION",
    "FACTORISED_DECOMPOSITION",
    "GP_LOCALISATION_GROUP",
    "JOINT_DECOMPOSITION",
    "POINTWISE_LOG_LIKELIHOOD_GROUP",
    "POSTERIOR_PREDICTIVE_GROUP",
    "RESIDUALS_GROUP",
    "add_pointwise_log_likelihood",
    "add_posterior_predictive",
    "add_residuals",
    "base_label",
    "component_variable",
    "gp_localisation",
    "pointwise_as_log_likelihood",
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

#: ArviZ's per-**observation** ``log_likelihood`` convention. Reserved by
#: ``results.md`` §6 and written since W2.8, by :func:`add_pointwise_log_likelihood`
#: and by nothing else. It is a distinct group from ``log_likelihood`` precisely
#: so the two decompositions cannot be confused — a variable here carries an
#: ``ampere_decomposition`` of ``"factorised"`` (independent noise; exact) or
#: ``"conditional_loo"`` (a GP's leave-one-out conditional terms).
POINTWISE_LOG_LIKELIHOOD_GROUP = "pointwise_log_likelihood"

#: The independent-noise decomposition: the family's own ``log_prob`` evaluated
#: pointwise. Exact, and the terms sum to the dataset's joint log-likelihood.
FACTORISED_DECOMPOSITION = "factorised"

#: The GP decomposition: ``log N(y_i | mu_i^{-i}, sigma_i^{2,-i})``, the
#: leave-one-out conditional terms. These are what ``arviz.loo`` consumes, and
#: they deliberately do **not** sum to the joint value — a GP likelihood has no
#: per-observation factorisation (``likelihoods.md`` §16).
CONDITIONAL_LOO_DECOMPOSITION = "conditional_loo"

#: The joint decomposition (**W5.9**): the leave-one-out conditional terms of a
#: :class:`~ampere.core.JointGaussianProcessNoise` group's **rotated outputs**.
#: A group's ``T`` channels are not independent, so a leave-one-out conditional
#: of one channel alone would condition on its sibling's value at the same
#: sample without saying so; the rotated outputs are the things that *are*
#: independent, and they are what this decomposition is indexed by — output
#: ``s`` stored under the group's ``s``-th member's label, on the grid every
#: channel shares. Like ``conditional_loo`` the terms do not sum to the joint
#: value. The rotation itself is parameter-dependent: recover it per draw with
#: ``JointGaussianProcessNoise.eigen(values)``.
JOINT_DECOMPOSITION = "joint"


#: The four views a complex dataset's derived variable may be stored as
#: (**W5.3**). The order matches ``results.md`` §8's list, and every refusal
#: message that names them uses this tuple, so the two cannot drift apart.
COMPONENTS: tuple[str, ...] = ("real", "imag", "abs", "phase")


def _component_of(values: np.ndarray, component: str) -> np.ndarray:
    """One of :data:`COMPONENTS` of complex-valued *values* — the one helper.

    **W5.3.** :func:`add_posterior_predictive`, :func:`add_residuals` and
    :func:`gp_localisation` all reach this, so the four views mean the same
    thing wherever a caller asks for one.
    """
    if component == "real":
        return values.real
    if component == "imag":
        return values.imag
    if component == "abs":
        return np.abs(values)
    return np.angle(values)  # "phase"


def _resolve_complex(
    label: str, complex_valued: bool, component: str | None, group: str, alternative: str
) -> None:
    """Refuse a complex-valued dataset unless *component* names how to view it.

    **W5.2** made the three derived groups refuse a complex-valued dataset the
    same way, deferring the "which component" question. **W5.3** answers it:
    ``component="real" | "imag" | "abs" | "phase"`` names the view, and the
    derived variable is stored as ``<label>_<component>`` — ``results.md`` §4
    splits a complex *observed* container into ``<label>_real`` and
    ``<label>_imag`` at emission time, and a derived group still holds one
    real variable per dataset, so the split is the caller's *choice* of
    which view rather than "both, always". ``component=None`` on a
    complex-valued dataset keeps the plain W5.2 refusal, naming the four
    choices; ``component=`` on a *real*-valued dataset is refused too — there
    is nothing to pick a component of.
    """
    if component is not None and component not in COMPONENTS:
        raise ResultsError(f"component must be one of {COMPONENTS}, got {component!r}.")
    if not complex_valued:
        if component is not None:
            raise ResultsError(
                f"dataset {label!r} is real-valued; component= only applies to a "
                f"complex-valued dataset, and there is nothing to pick a component of here."
            )
        return
    if component is None:
        raise ResultsError(
            f"dataset {label!r} is complex-valued, and a {group} group holds one real "
            f"variable per dataset. results.md §4 splits a complex observed container into "
            f"<label>_real and <label>_imag; the derived groups store one caller-named view "
            f'instead — pass component="real", "imag", "abs" or "phase" (W5.3), '
            f"stored as <label>_<component>. {alternative}"
        )


def component_variable(label: str, component: str | None) -> str:
    """The stored variable name for *label*, given an optional *component*.

    Public (no leading underscore) because :mod:`ampere.results.plots`'s
    three plot functions need the same naming to find what an ``add_*`` call
    with the same ``component=`` stored — the inverse, :func:`base_label`, is
    the other half.
    """
    return f"{label}_{component}" if component else label


def base_label(name: str, component: str | None) -> str:
    """Undo :func:`component_variable`: the dataset label a stored *name* is for.

    ``component=None`` is a no-op (a real dataset's variable is never
    suffixed), so a caller who never asked for a component sees the group's
    variable names unchanged — the byte-identical path W5.3 promises.
    """
    if component and name.endswith(f"_{component}") and len(name) > len(component) + 1:
        return name[: -(len(component) + 1)]
    return name


def _rename_for_component(
    mapping: Mapping[str, Any],
    labels: Sequence[str],
    complex_valued: Mapping[str, bool],
    component: str | None,
) -> dict[str, Any]:
    """*mapping* (keyed by base label), with a complex label's key renamed.

    Applied identically to the ``variables`` payload and to ``dims``, so a
    variable and its own dimension list are found under the same key.
    """
    renamed = dict(mapping)
    for label in labels:
        if complex_valued.get(label, False) and label in renamed:
            renamed[component_variable(label, component)] = renamed.pop(label)
    return renamed


def add_posterior_predictive(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    seed_stream: str = "posterior_predictive",
    component: str | None = None,
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
    component
        For a complex-valued dataset, which view to store —
        ``"real"``, ``"imag"``, ``"abs"`` or ``"phase"`` (**W5.3**); the
        derived variable is named ``<label>_<component>``. Required when any
        requested dataset is complex-valued; refused (there is nothing to
        pick a component of) when none is.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the problem is not the one the run was over, if a requested dataset
        is not in it, if a requested dataset's observed container is
        complex-valued and ``component`` is not one of the four named above,
        or if ``component`` is given and no requested dataset is complex.

    Notes
    -----
    Masked samples are **NaN**, not the observed value.
    :meth:`~ampere.core.dataset.Dataset.draw_observation` deliberately leaves a
    masked sample's observed value in place — it draws no noise for a sample
    nobody measured — and copying that into a group named "replicate" would
    invent a replicate that is exactly the datum, giving a discrepancy
    statistic a free perfect fit at every gap. The mask convention is
    ``add_residuals``', for the same reason it is there.

    A draw the forward model could not complete (a prior-rejected θ, a
    simulator crash) leaves its whole row NaN, exactly as the residual group
    does: ``simulate`` flags rather than raises, so one bad draw costs one row
    rather than the group.
    """
    _require_same_problem(tree, problem)
    labels = _requested(problem, datasets)
    complex_valued: dict[str, bool] = {}
    for label in labels:
        observed = problem.datasets[label].observed
        is_complex = np.asarray(observed.values).dtype.kind == "c"
        _resolve_complex(
            label,
            is_complex,
            component,
            "posterior-predictive",
            "Replicate the real datasets by name, or draw the replicates yourself with "
            "FittingProblem.simulate(observe=True).",
        )
        complex_valued[label] = is_complex
    thetas, kept = _stored_thetas(tree, problem, thin)
    chains, draws = thetas.shape[0], thetas.shape[1]
    variables = {
        label: np.full(
            (chains, draws, problem.datasets[label].observed.n_samples),
            np.nan,
            dtype=complex if complex_valued[label] else float,
        )
        for label in labels
    }
    for chain in range(chains):
        for draw in range(draws):
            simulation = problem.simulate(
                thetas[chain, draw], observe=True, stream=str(seed_stream)
            )
            if simulation.failed or simulation.observations is None:
                continue
            for label in labels:
                drawn = simulation.observations.get(label)
                if drawn is None:  # pragma: no cover - a failure short-circuits above
                    continue
                replicate = np.asarray(
                    drawn.values, dtype=complex if complex_valued[label] else float
                ).ravel()
                excluded = problem.datasets[label].effective_mask
                if excluded is not None:
                    replicate = np.where(np.asarray(excluded).ravel(), np.nan, replicate)
                variables[label][chain, draw] = replicate
    for label in labels:
        if complex_valued[label]:
            variables[label] = _component_of(variables[label], component)
    variables = _rename_for_component(variables, labels, complex_valued, component)
    dims, coords = _observed_axes(problem, labels)
    dims = _rename_for_component(dims, labels, complex_valued, component)
    extra_variables, variable_attrs, extra_attrs = _multi_axis_extras(problem, labels)
    return _attach(
        tree,
        POSTERIOR_PREDICTIVE_GROUP,
        variables,
        dims,
        coords,
        kept,
        {
            f"{ATTR_PREFIX}seed_stream": str(seed_stream),
            f"{ATTR_PREFIX}derivation": (
                "y_rep ~ p(y | theta_k) through FittingProblem.simulate(observe=True), so the "
                "replicates come from the same LikelihoodFamily.sample the likelihood scores "
                "with (inference.md §13); masked samples are NaN, never dropped "
                "(results.md §7)."
            ),
            **extra_attrs,
        },
        variable_attrs,
        extra_variables,
    )


def add_residuals(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    standardised: bool = True,
    component: str | None = None,
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
    component
        For a complex-valued dataset, which view to store —
        ``"real"``, ``"imag"``, ``"abs"`` or ``"phase"`` (**W5.3**); the
        derived variable is named ``<label>_<component>``. Required when any
        requested dataset is complex-valued; refused when none is.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the problem is not the one the run was over, if a requested dataset
        is not in it, if a requested dataset's observed container is
        complex-valued and ``component`` is not one of the four named above,
        if ``component`` is given and no requested dataset is complex, or if
        ``standardised=True`` and a dataset has no uncertainties to
        standardise by.
    """
    _require_same_problem(tree, problem)
    labels = _requested(problem, datasets)
    complex_valued: dict[str, bool] = {}
    for label in labels:
        observed = problem.datasets[label].observed
        is_complex = np.asarray(observed.values).dtype.kind == "c"
        _resolve_complex(
            label,
            is_complex,
            component,
            "residuals",
            "Residualise the real datasets by name, or compute observed - predicted "
            "yourself with FittingProblem.simulate() and Likelihood.conditional's inputs.",
        )
        complex_valued[label] = is_complex
        if standardised and observed.uncertainty is None:
            raise ResultsError(
                f"dataset {label!r} has no uncertainties, so a standardised residual is "
                f"undefined. Pass standardised=False for residuals in the data's own "
                f"units, or attach uncertainties to the observed container."
            )
    thetas, kept = _stored_thetas(tree, problem, thin)
    chains, draws = thetas.shape[0], thetas.shape[1]
    variables = {
        label: np.full(
            (chains, draws, problem.datasets[label].observed.n_samples),
            np.nan,
            dtype=complex if complex_valued[label] else float,
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
                variables[label][chain, draw] = _residual_of(
                    dataset, predicted, standardised, complex_valued=complex_valued[label]
                )
    for label in labels:
        if complex_valued[label]:
            variables[label] = _component_of(variables[label], component)
    variables = _rename_for_component(variables, labels, complex_valued, component)
    dims, coords = _observed_axes(problem, labels)
    dims = _rename_for_component(dims, labels, complex_valued, component)
    extra_variables, variable_attrs, extra_attrs = _multi_axis_extras(problem, labels)
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
            **extra_attrs,
        },
        variable_attrs,
        extra_variables,
    )


def _residual_of(
    dataset: Dataset,
    predicted: FunctionSamples,
    standardised: bool,
    *,
    complex_valued: bool = False,
) -> np.ndarray:
    """One draw's signed residual for one dataset, masked samples NaN.

    ``masked_uncertainty()`` inflates a masked sample's sigma to ``inf``, which
    would send its residual to a perfectly white **zero** — the one value a
    whiteness statistic must not be handed for a sample nobody measured. The
    mask is therefore applied as NaN afterwards, which is what ``results.md``
    §7 asks for and what lets the statistic "decide for itself what to do with
    a gap".

    ``complex_valued`` (**W5.3**) keeps the subtraction complex rather than
    truncating it, for a caller who is about to reduce it with ``component=``.
    """
    dtype = complex if complex_valued else float
    observed = dataset.observed
    residual = (
        np.asarray(observed.values, dtype=dtype).ravel()
        - np.asarray(predicted.values, dtype=dtype).ravel()
    )
    if standardised:
        residual = residual / np.asarray(observed.masked_uncertainty(), dtype=float).ravel()
    excluded = dataset.effective_mask
    if excluded is not None:
        residual = np.where(np.asarray(excluded).ravel(), np.nan, residual)
    return residual


def _rename_roles_for_component(
    mapping: Mapping[str, Any],
    labels: Sequence[str],
    complex_valued: Mapping[str, bool],
    component: str | None,
    roles: Sequence[str] = ("mean", "variance"),
) -> dict[str, Any]:
    """:func:`_rename_for_component`, for a mapping keyed by ``<label>_<role>``.

    ``gp_localisation``'s own naming (``<label>_mean``, ``<label>_variance``)
    is role-suffixed rather than bare, so the component goes **between** the
    label and the role — ``<label>_<component>_mean`` — keeping both roles of
    one label discoverable under the one renamed base
    (:func:`component_variable`), which is what lets a reader recover "the
    label" from either half of the pair the same way.
    """
    renamed = dict(mapping)
    for label in labels:
        if not complex_valued.get(label, False):
            continue
        stored = component_variable(label, component)
        for role in roles:
            old_key = f"{label}_{role}"
            if old_key in renamed:
                renamed[f"{stored}_{role}"] = renamed.pop(old_key)
    return renamed


def gp_localisation(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
    at: Any = None,
    component: str | None = None,
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
        array of coordinates means. Not combined with a multi-axis kind's own
        default coordinate (an explicit grid already says what to plot
        against).
    component
        For a complex-valued dataset's conditioned mean, which view to store
        — ``"real"``, ``"imag"``, ``"abs"`` or ``"phase"`` (**W5.3**); stored
        as ``<label>_<component>_mean`` beside ``<label>_<component>_variance``
        (the variance itself is always real and untouched — see
        :class:`~ampere.core.VisibilitySet`). Required when any requested
        dataset is complex-valued; refused when none is.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the problem is not the one the run was over, if no dataset has a
        ``GaussianProcessNoise`` model, if a requested dataset's observed
        container is complex-valued and ``component`` is not one of the four
        named above, or if ``component`` is given and no requested dataset is
        complex.
    """
    _require_same_problem(tree, problem)
    labels = _gp_requested(problem, datasets)
    complex_valued: dict[str, bool] = {}
    for label in labels:
        observed = problem.datasets[label].observed
        is_complex = np.asarray(observed.values).dtype.kind == "c"
        _resolve_complex(
            label,
            is_complex,
            component,
            "gp_localisation",
            "Localise the real datasets by name, or call Likelihood.conditional yourself "
            "and choose a component of the complex conditioned mean.",
        )
        complex_valued[label] = is_complex
    thetas, kept = _stored_thetas(tree, problem, thin)
    chains, draws = thetas.shape[0], thetas.shape[1]
    grid = None if at is None else np.asarray(at, dtype=float)
    size = {
        label: (problem.datasets[label].observed.n_samples if grid is None else int(grid.shape[0]))
        for label in labels
    }
    variables: dict[str, np.ndarray] = {}
    for label in labels:
        variables[f"{label}_mean"] = np.full(
            (chains, draws, size[label]), np.nan, dtype=complex if complex_valued[label] else float
        )
        variables[f"{label}_variance"] = np.full((chains, draws, size[label]), np.nan, dtype=float)
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
                variables[f"{label}_mean"][chain, draw] = np.asarray(
                    conditional.mean, dtype=complex if complex_valued[label] else float
                ).ravel()
                variables[f"{label}_variance"][chain, draw] = np.asarray(
                    conditional.variance, dtype=float
                ).ravel()
    for label in labels:
        if complex_valued[label]:
            variables[f"{label}_mean"] = _component_of(variables[f"{label}_mean"], component)
    variables = _rename_roles_for_component(variables, labels, complex_valued, component)
    extra_attrs: dict[str, str] = {}
    variable_attrs: dict[str, dict[str, str]] = {}
    extra_variables: dict[str, tuple[list[str], np.ndarray]] = {}
    if grid is None:
        dims, coords = _observed_axes(problem, labels, roles=("mean", "variance"))
        dims = _rename_roles_for_component(dims, labels, complex_valued, component)
        extra_variables, variable_attrs, extra_attrs = _multi_axis_extras(problem, labels)
    else:
        dims = {}
        coords = {}
        for label in labels:
            name = f"{label}_gp_axis"
            coords[name] = grid
            dims[f"{label}_mean"] = [name]
            dims[f"{label}_variance"] = [name]
        dims = _rename_roles_for_component(dims, labels, complex_valued, component)
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
            **extra_attrs,
        },
        variable_attrs,
        extra_variables,
    )


# ---------------------------------------------------------------------------
# ArviZ's per-observation convention, on an explicit call
# ---------------------------------------------------------------------------


def add_pointwise_log_likelihood(
    tree: Any,
    problem: FittingProblem,
    *,
    datasets: Sequence[str] | None = None,
    thin: int = 1,
) -> Any:
    """Compute the per-**observation** log-likelihood terms into the reserved group.

    ``results.md`` §6 names two decompositions and forbids their sharing one
    group, so the terms land in :data:`POINTWISE_LOG_LIKELIHOOD_GROUP` rather
    than in ``log_likelihood`` (which stays the per-*dataset* split), and each
    variable declares which decomposition it is:

    * :data:`FACTORISED_DECOMPOSITION` — independent noise. The family's own
      ``log_prob`` evaluated pointwise; exact, and the terms sum to the
      dataset's stored joint log-likelihood.
    * :data:`CONDITIONAL_LOO_DECOMPOSITION` — a GP. The leave-one-out
      conditional terms ``log N(y_i | mu_i^{-i}, sigma_i^{2,-i})``, from the
      same factorisation the marginal likelihood forms. They are what
      ``arviz.loo`` consumes and they deliberately do **not** sum to the joint
      value.

    Not written by default, and this is the only thing that writes it (§6's
    ruling of 2026-09-03, unchanged by the computation becoming available): the
    group is ``N_draws x N_obs`` per dataset, and a run that will never be
    asked for LOO should not carry it.

    Where the dataset's solver does not supply ``conditional_loo`` — today
    :class:`~ampere.core.likelihood.QuasisepGP` on the reference path and its
    jax twin — this **refuses by name and does not fall back to a dense
    solve**. Substituting :class:`~ampere.core.likelihood.DenseGP` would turn
    an ``O(N)`` fit's diagnostic into an ``O(N^3)`` one silently, which at the
    sizes the quasiseparable solve exists for is a different program rather
    than a slower answer.

    Parameters
    ----------
    tree
        The run, as :func:`ampere.results.emit` produced it. Returned with the
        group attached; the tree is modified in place.
    problem
        The problem the run was over; its ``problem_hash`` must match.
    datasets
        Which datasets to decompose; all of them by default.
    thin
        Take every ``thin``-th draw; the retained indices become the group's
        own ``draw`` coordinate.

    Returns
    -------
    xarray.DataTree
        *tree*, with :data:`POINTWISE_LOG_LIKELIHOOD_GROUP` attached. Masked
        samples are NaN on the container's own coordinate axis, and a draw the
        forward model could not complete is a NaN row —
        :func:`pointwise_as_log_likelihood` is how the group is handed to
        ``arviz.loo``, which has no use for either.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the problem is not the one the run was over, if a requested dataset
        is not in it, if its likelihood marginalises over latent values (whose
        per-observation terms are conditional on values inference owns), or if
        its solver does not implement the leave-one-out conditionals.
    """
    _require_same_problem(tree, problem)
    labels = _requested(problem, datasets)
    decompositions = {label: _decomposition_of(problem, label) for label in labels}
    groups = _requested_groups(problem, labels)
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
            routed = problem.mapping.distribute(simulation.parameters)
            for label in labels:
                dataset = problem.datasets[label]
                predicted = simulation.predicted.get(label)
                if predicted is None:  # pragma: no cover - a failure short-circuits above
                    continue
                # W5.9: a channel of a joint noise group has no per-observation
                # decomposition of its own; the group's rotated outputs are the
                # independent things, and they are taken once per group below.
                if problem.datasets.group_of(label) is not None:
                    continue
                split = dataset.route(routed.get(label, {}))
                variables[label][chain, draw] = _pointwise_of(
                    dataset, predicted, split.get(LIKELIHOOD_COMPONENT)
                )
            for group in groups:
                rotated = _joint_pointwise_of(problem, group, simulation.predicted, routed)
                for member, column in rotated.items():
                    if member in variables:
                        variables[member][chain, draw] = column
    dims, coords = _observed_axes(problem, labels)
    distinct = sorted(set(decompositions.values()))
    attached = _attach(
        tree,
        POINTWISE_LOG_LIKELIHOOD_GROUP,
        variables,
        dims,
        coords,
        kept,
        {
            # One value where the run is unanimous, which is the case §6
            # describes. A joint fit that mixes an independent-noise dataset
            # with a GP one has two decompositions and no single honest answer,
            # so the group says "mixed" and the per-variable attribute — always
            # written — is where a consumer reads which is which. Reporting one
            # of the two would be reporting a falsehood about the other.
            f"{ATTR_PREFIX}decomposition": distinct[0] if len(distinct) == 1 else "mixed",
            f"{ATTR_PREFIX}decomposition_note": (
                "one term per retained observation, on the observed container's own coordinate "
                "axis with masked samples as NaN. 'factorised' terms sum to the dataset's joint "
                "log-likelihood; 'conditional_loo' terms deliberately do not, because a GP "
                "likelihood has no per-observation factorisation (results.md §6). 'joint' terms "
                "are a JointGaussianProcessNoise group's rotated outputs -- output s under the "
                "group's s-th member's label, on the grid every channel shares -- and likewise "
                "do not sum to the joint value; the rotation is parameter-dependent and is "
                "recovered per draw from JointGaussianProcessNoise.eigen(values)."
            ),
        },
    )
    for label in labels:
        attached[POINTWISE_LOG_LIKELIHOOD_GROUP][label].attrs[f"{ATTR_PREFIX}decomposition"] = (
            decompositions[label]
        )
    return attached


def _decomposition_of(problem: FittingProblem, label: str) -> str:
    """Which of §6's two decompositions this dataset's likelihood supplies.

    Refuses the latent case up front rather than after a full pass of forward
    models: ``pointwise_log_prob`` has no unconditional answer there, and
    finding that out on the last draw of a long run is a worse way to be told.
    """
    likelihood = problem.datasets[label].likelihood
    if likelihood.marginalisation is Marginalisation.LATENT:
        raise ResultsError(
            f"dataset {label!r} is fitted with a latent-variable likelihood, whose "
            f"per-observation terms are conditional on latent values inference owns rather "
            f"than this contract (likelihoods.md §7). Use the per-dataset log_likelihood group "
            f"every run already carries."
        )
    if problem.datasets.group_of(label) is not None:
        return JOINT_DECOMPOSITION
    if isinstance(likelihood.noise, GaussianProcessNoise):
        return CONDITIONAL_LOO_DECOMPOSITION
    return FACTORISED_DECOMPOSITION


def _requested_groups(problem: FittingProblem, labels: Sequence[str]) -> tuple[str, ...]:
    """The joint noise groups any requested dataset belongs to (**W5.9**).

    A group's rotated outputs are computed together or not at all, so asking
    for one channel of a group asks for the group; the other channels' columns
    are simply not stored when they were not requested.
    """
    found: list[str] = []
    for label in labels:
        group = problem.datasets.group_of(label)
        if group is not None and group not in found:
            found.append(group)
    return tuple(found)


def _joint_pointwise_of(
    problem: FittingProblem,
    group: str,
    predicted: Mapping[str, FunctionSamples],
    routed: Mapping[str, Any],
) -> dict[str, np.ndarray]:
    """One draw's rotated-output terms for one joint noise group (**W5.9**).

    Returns one full-length row per member label, masked samples NaN, with the
    ``s``-th member carrying the ``s``-th rotated output. That pairing is a
    storage convention rather than a claim that output ``s`` "is" channel
    ``s``: the outputs live in ``B``'s eigenbasis, which rotates with θ, and
    the decomposition attribute says so. What makes it the right convention is
    that there are exactly ``T`` of each and they share one grid, so the group
    needs no coordinate axis of its own.
    """
    noise = problem.datasets.joint[group]
    members = noise.datasets
    first = problem.datasets[members[0]]
    values = dict(routed.get(group, {}))
    rows = {
        label: np.full(problem.datasets[label].observed.n_samples, np.nan, dtype=float)
        for label in members
    }
    missing = [label for label in members if predicted.get(label) is None]
    if missing:  # pragma: no cover - a failure short-circuits the caller
        return rows
    try:
        retain = problem.datasets._group_retained(group, predicted)
        residuals = [
            problem.datasets._group_residual(label, predicted[label], retain) for label in members
        ]
        if not np.any(retain):
            return rows
        sigma = noise.sigma(first.observed, retain, values)
        variance = (
            np.zeros(int(np.count_nonzero(retain)))
            if sigma is None
            else np.asarray(sigma, dtype=float) ** 2
        )
        terms = noise.pointwise_log_prob(
            residuals,
            variance,
            sample_coordinates(first.observed)[retain],
            values,
            kernel=noise.kernel_for(first.observed),
        )
    except (DatasetError, LikelihoodError) as error:
        raise ResultsError(
            f"joint noise group {group!r} cannot supply the per-observation decomposition: {error}"
        ) from error
    for index, label in enumerate(members):
        rows[label][retain] = np.asarray(terms, dtype=float)[:, index]
    return rows


def _pointwise_of(
    dataset: Dataset, predicted: FunctionSamples, values: Mapping[str, Any] | None
) -> np.ndarray:
    """One draw's per-observation terms for one dataset, masked samples NaN.

    :meth:`~ampere.core.likelihood.Likelihood.pointwise_log_prob` returns one
    term per *retained* sample, which is the shortened axis a likelihood works
    on. The group's axis is the container's own, so the terms are scattered
    back into it — the same rule the residual group follows, and the reason
    a plot of either lines up with ``observed_data`` without an alignment step.

    A refusal from the likelihood is re-raised as this contract's error, with
    the solver's own message kept: "QuasisepGP does not implement the
    leave-one-out conditional terms" is the sentence a user needs, and burying
    it under "could not compute the pointwise group" would lose it.
    """
    observed = dataset.observed
    weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
    retain = weights > 0.0
    row = np.full(retain.size, np.nan, dtype=float)
    if not np.any(retain):
        return row
    try:
        terms = np.asarray(
            dataset.likelihood.pointwise_log_prob(predicted, observed, values), dtype=float
        ).ravel()
    except LikelihoodError as error:
        raise ResultsError(
            f"dataset {dataset.label!r} cannot supply the per-observation decomposition: {error}"
        ) from error
    if terms.size != int(np.count_nonzero(retain)):  # pragma: no cover - a contract violation
        raise ResultsError(
            f"dataset {dataset.label!r}'s likelihood returned {terms.size} per-observation "
            f"term(s) for {int(np.count_nonzero(retain))} retained sample(s); "
            f"Likelihood.pointwise_log_prob returns one term per retained sample."
        )
    row[retain] = terms
    return row


def pointwise_as_log_likelihood(tree: Any) -> Any:
    """A view of *tree* whose ``log_likelihood`` group is the per-observation one.

    ``arviz.loo`` and ``arviz.waic`` read the group **named**
    ``log_likelihood``, and ampere deliberately keeps its per-*dataset*
    decomposition there (``results.md`` §15 R6): calling ``arviz.loo`` on a run
    straight out of :func:`~ampere.results.emit` computes
    leave-one-*dataset*-out, which is meaningful for a joint fit and useless
    for a single-dataset fit of 10⁵ points. This function is the one-line
    bridge, so that getting per-observation LOO is an explicit swap rather than
    a hand-edited tree — and so that nobody has to discover by experiment which
    of the two ``arviz.loo`` just gave them.

    The returned tree is a shallow copy: the original run keeps its own
    ``log_likelihood`` group, and both trees share the underlying arrays.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the run has no ``pointwise_log_likelihood`` group, naming
        :func:`add_pointwise_log_likelihood` as the precondition.
    """
    children = getattr(tree, "children", None)
    if children is None or POINTWISE_LOG_LIKELIHOOD_GROUP not in children:
        raise ResultsError(
            f"this run has no {POINTWISE_LOG_LIKELIHOOD_GROUP!r} group, so there is nothing to "
            f"hand to arviz.loo as a per-observation decomposition. It is not stored by default "
            f"— the cost is N_draws x N_obs per dataset (results.md §6) — so compute it first: "
            f"ampere.results.add_pointwise_log_likelihood(tree, problem)."
        )
    view = tree.copy()
    view[LOG_LIKELIHOOD_GROUP] = tree[POINTWISE_LOG_LIKELIHOOD_GROUP].dataset
    view[LOG_LIKELIHOOD_GROUP].attrs.update(dict(tree[POINTWISE_LOG_LIKELIHOOD_GROUP].attrs))
    return view


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
        if len(names) == len(observed.axes):
            # One axis: it doubles as the dimension's own coordinate, exactly
            # as observed_data's does (emission._data_groups). A kind with
            # several axes takes the *other* branch below (_container_dims
            # gives it one joint ``<label>_index`` dimension, and its axes are
            # ordinary variables rather than that dimension's index) — zipping
            # the two mismatched-length sequences here would silently pair
            # only the first axis with the index dimension, which is wrong
            # rather than merely incomplete: it would put ``u``'s values on a
            # dimension that means "position in the sample list, arbitrary
            # order" (results.md §4).
            for axis, name in zip(observed.axes, names, strict=True):
                coords[name] = np.asarray(axis.values, dtype=float).ravel()
        for role in roles or ("",):
            dims[f"{label}_{role}" if role else label] = names
    return dims, coords


#: One extra (non-draw) variable for :func:`_attach` to assign onto a group
#: after arviz has built it: ``(dims, values)``, exactly xarray's own
#: ``Dataset.assign`` shorthand.
_ExtraVariable = tuple[list[str], np.ndarray]


def _multi_axis_extras(
    problem: FittingProblem, labels: Sequence[str]
) -> tuple[dict[str, _ExtraVariable], dict[str, dict[str, str]], dict[str, str]]:
    """A multi-axis kind's raw axes and precomputed default coordinate.

    **W5.3.** Each of the kind's own axes becomes an ordinary variable on the
    shared ``<label>_index`` dimension, mirroring what
    ``emission._data_groups`` already does for ``observed_data``/
    ``constant_data`` (results.md §4), so a derived group carries the same
    information its ``observed_data`` sibling does. Where the kind declares
    :attr:`~ampere.core.results_schema.FunctionSamples.PLOT_COORDINATE`, it is
    resolved **here**, against the *live* container (the only place its real
    axis units are available), and stored as a further variable —
    ``ampere.results._plotting.coordinate_of`` reads it back at plot time
    without needing to know which kind produced it.

    These are returned separately from the per-draw ``variables``/``dims``
    passed to ``arviz.from_dict`` rather than merged into them: arviz assumes
    *every* variable named there carries the group's own leading
    ``(chain, draw)`` dimensions, and a kind's raw axes do not — they are one
    value per sample, not per draw. :func:`_attach` assigns them onto the
    built group afterwards, as plain xarray variables.

    Single-axis kinds return nothing (``_observed_axes`` already gives them
    their one coordinate), which is what keeps every existing plot
    byte-identical.

    Returns ``(extra_variables, variable_attrs, group_attrs)``:
    *extra_variables* for :func:`_attach`'s post-hoc assignment; the per-axis
    unit and the :data:`~ampere.results._plotting.PLOT_AXIS_ATTR` marker for
    each; the default coordinate's axis label, on the group itself.
    """
    extra_variables: dict[str, _ExtraVariable] = {}
    variable_attrs: dict[str, dict[str, str]] = {}
    group_attrs: dict[str, str] = {}
    for label in labels:
        observed = problem.datasets[label].observed
        names = list(_container_dims(label, observed))
        if len(names) == len(observed.axes):
            continue
        index_dim = names[0]
        axes_map = {axis.name: axis for axis in observed.axes}
        for axis in observed.axes:
            axis_variable = f"{label}_{axis.name}"
            extra_variables[axis_variable] = (
                [index_dim],
                np.asarray(axis.values, dtype=float).ravel(),
            )
            attrs = {_p.PLOT_AXIS_ATTR: axis.name}
            if axis.unit is not None:
                attrs["units"] = str(axis.unit.to_string())
            variable_attrs[axis_variable] = attrs
        default = type(observed).PLOT_COORDINATE
        if default is None:
            continue
        if isinstance(default, str):
            axis = axes_map[default]
            coordinate = np.asarray(axis.values, dtype=float).ravel()
            coordinate_label = format_axis_label(axis.name, axis.unit)
        else:
            raw_coordinate, coordinate_label = default(axes_map)
            coordinate = np.asarray(raw_coordinate, dtype=float).ravel()
        coordinate_variable = f"{label}{_p.PLOT_COORDINATE_VARIABLE_SUFFIX}"
        extra_variables[coordinate_variable] = ([index_dim], coordinate)
        group_attrs[f"{_p.PLOT_COORDINATE_LABEL_ATTR_PREFIX}{label}"] = str(coordinate_label)
    return extra_variables, variable_attrs, group_attrs


def _attach(
    tree: Any,
    group: str,
    variables: Mapping[str, np.ndarray],
    dims: Mapping[str, list[str]],
    coords: Mapping[str, Any],
    kept: np.ndarray,
    attrs: Mapping[str, Any],
    variable_attrs: Mapping[str, Mapping[str, str]] | None = None,
    extra_variables: Mapping[str, _ExtraVariable] | None = None,
) -> Any:
    """Build one derived group and hang it off the run, coordinates and all.

    The retained draw indices become the group's own ``draw`` coordinate, so a
    thinned group states which draws it holds rather than renumbering them from
    zero and quietly losing the correspondence to the posterior.

    *variable_attrs* (**W5.3**) sets per-variable attributes afterwards — a
    multi-axis kind's raw axes carry their unit and the
    :data:`~ampere.results._plotting.PLOT_AXIS_ATTR` marker this way, which
    ``arviz.from_dict``'s own ``dims``/``coords`` arguments have no slot for.

    *extra_variables* (**W5.3**) assigns further variables onto the built
    group *after* arviz has built it, as plain xarray variables with their own
    dims — bypassing arviz's assumption that every variable named in the
    ``group`` dict above carries the group's own leading ``(chain, draw)``.
    """
    arviz = _require_arviz()
    resolved = dict(coords)
    resolved[DRAW_DIM] = kept
    built = arviz.from_dict({group: dict(variables)}, dims=dict(dims), coords=resolved)
    for name, child in built.children.items():
        tree[name] = child
    tree[group].attrs.update(dict(attrs))
    if extra_variables:
        dataset = tree[group].dataset.assign(
            {name: tuple(spec) for name, spec in extra_variables.items()}
        )
        tree[group] = dataset
        tree[group].attrs.update(dict(attrs))
    for name, extra in (variable_attrs or {}).items():
        tree[group][name].attrs.update(dict(extra))
    return tree
