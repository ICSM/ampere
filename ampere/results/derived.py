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

Everything below is a declared surface with no implementation: Phase 2 lands
these behind the backends that make them computable, and the point of W1.8 is
that their shape, their group names and their cost policy are settled first.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from ampere.core.dataset import FittingProblem

__all__ = [
    "GP_LOCALISATION_GROUP",
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
    """
    raise NotImplementedError(f"signed residuals are not implemented yet. {_PHASE_2}")


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
    that refusal rather than silently skipping them.
    """
    raise NotImplementedError(f"GP localisation is not implemented yet. {_PHASE_2}")
