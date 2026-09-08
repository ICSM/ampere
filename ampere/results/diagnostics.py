"""Families B and C of ``diagnostics.md``, computed rather than drawn.

``diagnostics.md`` §6 places the two post-fit diagnostic families in
``ampere.results``: family B (residual whiteness and posterior-predictive
checks) because it consumes a run and produces a number, family C
(GP localisation) because it consumes a run and produces a picture of where
the model is deficient. The renderers are in :mod:`ampere.results.plots`; the
arithmetic they draw is here, so that a user who wants the number without the
figure — a pipeline deciding whether to switch the flexible likelihood on —
never has to import matplotlib to get it.

Family B: whiteness, on irregular coordinates
---------------------------------------------
``diagnostics.md`` §3.2 named two honest options for testing whiteness when
§4.2 forbids assuming a regular grid, and deliberately picked neither:
resample onto a nominal grid (simple, lossy), or "a separation-binned
autocorrelation / structure-function statistic native to irregular coordinate
spacing (more faithful to §4.2, more design work, no off-the-shelf
implementation to point to)". ``likelihoods.md`` §12 and this item settle it:
**separation-binned, permutation-calibrated**, hand-rolled in numpy. The
statistic is defined in :func:`separation_binned_autocorrelation` and
:func:`residual_whiteness`, and it is worth being exact about what it is,
because "Ljung-Box" is a name for a thing this deliberately is not.

Let ``r_i`` be the signed standardised residual at coordinate ``x_i``
(:func:`~ampere.results.derived.add_residuals`; signed because
autocorrelation is sign-sensitive and per-observation log-likelihood is not —
``diagnostics.md`` §3.3's precise gap). For every unordered pair ``i < j``
whose separation ``s_ij = |x_i - x_j|`` falls inside the tested window, bin
the pair by ``s_ij`` and average the product:

.. code-block:: text

    rho_b = (1 / n_b) * sum_{(i,j) in bin b} r_i r_j / v,    v = mean_i r_i^2

``rho_b`` is the autocorrelation of the residual field at separation ``b``, and
``v`` is the zero-separation normalisation, so ``rho_b`` sits in ``[-1, 1]``
with the same reading as an ordinary autocorrelation coefficient. The window is
``bins`` bins wide, each one median-nearest-neighbour-spacing wide by default,
which is exactly what makes the reduction honest: **on an evenly spaced axis
bin b contains precisely the lag-b pairs, so rho_b is the lag-b sample
autocorrelation and nothing has been given up by generalising**. On an
irregular axis the same bins hold whatever pairs actually have that
separation, which is the generalisation §4.2 asks for and which resampling
would have thrown away.

The aggregate is Ljung-Box-*shaped* — a pair-count-weighted sum of squared
autocorrelations —

.. code-block:: text

    Q = sum_b n_b * rho_b^2

with ``n_b`` the number of pairs in bin ``b``. Weighting by pair count is the
irregular-coordinate analogue of Ljung-Box's ``n (n + 2) / (n - k)``: a bin
with more pairs pins its ``rho_b`` down better and should count for more.

**But Q has no tabulated null distribution**, and that is the whole reason the
classical test cannot simply be reused. Ljung-Box's asymptotic chi-square rests
on equal lag spacing and on ``n_k`` being a deterministic function of ``n`` and
``k``; here the pair counts depend on the sampling geometry, neighbouring bins
share points, and the residuals of a fitted model are not independent of one
another anyway. So the calibration is **permutation**: shuffle the residual
values across the coordinates — which destroys any relationship between value
and position while preserving the residual marginal *and* the sampling
geometry exactly — recompute ``Q``, and report

.. code-block:: text

    p = (1 + #{Q_perm >= Q_obs}) / (1 + n_permutations)

the add-one estimator, which is never zero and is exact under the
exchangeability null rather than asymptotic. Randomness comes from
:func:`ampere.core.rng.substream` off the run's own recorded seed, so the same
stored run gives the same p-value tomorrow — the reason ``lowering.md`` §9.2
insists a diagnostic draws from its own named stream and not from the fit's.

Family B: the posterior-predictive half
---------------------------------------
``diagnostics.md`` §3.3 observes that ``DEVELOPMENT_PLAN.md`` §4.6's stored
per-draw log-likelihood makes the posterior-predictive discrepancy check cheap
*already*, while giving the whiteness test nothing. :func:`chi_square_pvalue`
is that cheap half, and it is deliberately narrow: for a Gaussian family with
independent noise and no noise nuisance parameters, the stored per-dataset
log-likelihood is an exact affine function of the chi-square discrepancy
``T(y, theta) = sum_i r_i(theta)^2``, whose reference distribution under the
model is ``chi^2`` on the retained count — so the Bayesian p-value needs no
replicate draws at all. Where that does not hold (a GP, a fitted error scale,
censoring, a non-Gaussian family) the function refuses and names
:func:`~ampere.results.derived.add_posterior_predictive`, because sampling
replicates is then the honest route and pretending otherwise would be exactly
the "Gaussian assumption bolted on at the diagnostic layer" ``inference.md``
§13 forbids.

Family C: the score
-------------------
:func:`gp_localisation_score` converts the conditioned GP mean and variance
(:func:`~ampere.results.derived.gp_localisation`) into the shared
:class:`ampere.core.AnomalyScore`, with ``provenance`` and
``interpretation_notes`` populated — the latter carrying
:data:`~ampere.results.plots.GP_LOCALISATION_CAVEAT`, which is how
``diagnostics.md`` §4.3's caveat reaches a user who extracts the numbers and
never looks at a figure.
"""

from __future__ import annotations

import dataclasses
from typing import Any

import numpy as np
import scipy.stats as st

from ampere.core.exceptions import ResultsError
from ampere.core.results_schema import AnomalyScore
from ampere.core.rng import generator

from ._plotting import (
    coordinate_of,
    gp_datasets,
    likelihood_specs,
    require_group,
    run_seed,
    select_datasets,
)
from .derived import GP_LOCALISATION_GROUP, RESIDUALS_GROUP
from .emission import LOG_LIKELIHOOD_GROUP

__all__ = [
    "GP_LOCALISATION_PROVENANCE",
    "WHITENESS_STREAM",
    "ChiSquareCheck",
    "WhitenessTest",
    "chi_square_pvalue",
    "gp_localisation_datasets",
    "gp_localisation_score",
    "residual_whiteness",
    "separation_binned_autocorrelation",
]

#: ``AnomalyScore.provenance`` for family C, as ``diagnostics.md`` §5 names it.
GP_LOCALISATION_PROVENANCE = "gp_localisation_postfit"

#: The named RNG sub-stream the permutation calibration draws from
#: (``lowering.md`` §9.2). Distinct from ``"simulate"`` and from
#: ``"posterior_predictive"`` for the reason that policy exists: running a
#: diagnostic must not change any other consumer's draws.
WHITENESS_STREAM = "residual_whiteness"

#: Cap on the pair budget, so a 10^5-point spectrum does not silently ask for
#: a gigabyte of permutation products. Pairs above it are subsampled uniformly
#: (seeded), and the record says so.
DEFAULT_MAX_PAIRS = 200_000

#: How many posterior draws the test uses by default. Every retained draw gets
#: its own statistic and its own permutation calibration, so the cost is linear
#: in this; the pooled answer is the median across them.
DEFAULT_MAX_DRAWS = 32


# ---------------------------------------------------------------------------
# The statistic
# ---------------------------------------------------------------------------


def _pairs_within(
    coordinates: np.ndarray, max_separation: float, *, max_pairs: int, rng: np.random.Generator
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    """``(i, j, separation, subsampled)`` for every pair closer than *max_separation*.

    Found by sorting and :func:`numpy.searchsorted` rather than by forming the
    full ``n x n`` separation matrix: the window is a few median spacings wide,
    so the number of pairs is linear in ``n`` and the matrix would be the only
    quadratic thing in the computation.
    """
    order = np.argsort(coordinates, kind="stable")
    sorted_coordinates = coordinates[order]
    upper = np.searchsorted(sorted_coordinates, sorted_coordinates + max_separation, side="right")
    counts = upper - np.arange(sorted_coordinates.size) - 1
    counts = np.maximum(counts, 0)
    total = int(counts.sum())
    if total == 0:
        empty = np.empty(0, dtype=np.intp)
        return empty, empty, np.empty(0, dtype=float), False
    left = np.repeat(np.arange(sorted_coordinates.size), counts)
    offsets = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
    right = left + 1 + offsets
    separations = sorted_coordinates[right] - sorted_coordinates[left]
    subsampled = False
    if total > max_pairs:
        keep = rng.choice(total, size=max_pairs, replace=False)
        keep.sort()
        left, right, separations = left[keep], right[keep], separations[keep]
        subsampled = True
    return order[left], order[right], separations, subsampled


def _bin_pairs(
    separations: np.ndarray, bins: int, max_separation: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(order, segment_starts, pair_counts)`` for the non-empty separation bins.

    Pairs are sorted by bin so that the per-bin sums are one
    :func:`numpy.add.reduceat` over the last axis, which vectorises across
    draws and across permutations at once. Empty bins are dropped rather than
    reported as zero: a bin no pair falls in is a statement about the sampling,
    not about the residuals, and ``reduceat`` has no sane answer for one.
    """
    edges = np.linspace(0.0, max_separation, bins + 1)
    # Right-closed intervals, ``(b w, (b+1) w]``: a separation of exactly one
    # median spacing is lag 1, not a boundary case shared with lag 2. That is
    # what makes bin ``b`` *be* lag ``b + 1`` on an evenly spaced axis.
    index = np.clip(np.digitize(separations, edges[1:-1], right=True), 0, bins - 1)
    order = np.argsort(index, kind="stable")
    counts = np.bincount(index, minlength=bins)
    occupied = np.flatnonzero(counts > 0)
    starts = np.concatenate([[0], np.cumsum(counts)])[occupied]
    return order, starts.astype(np.intp), counts[occupied]


def separation_binned_autocorrelation(
    coordinates: Any,
    residuals: Any,
    *,
    bins: int = 8,
    max_separation: float | None = None,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(separation, rho, pair_count)`` — the module docstring's ``rho_b``.

    The statistic native to irregular coordinates that ``diagnostics.md`` §3.2
    left open and this item settles. On an evenly spaced axis with the default
    window it is the ordinary lag-1..lag-``bins`` sample autocorrelation; on an
    irregular one it is the same quantity computed from whatever pairs actually
    have each separation.

    Parameters
    ----------
    coordinates
        The data's own 1-D coordinate axis. Never resampled.
    residuals
        Signed standardised residuals, same length. NaN marks a sample that
        was masked or could not be scored; such samples take no part.
    bins
        How many separation bins, i.e. how many lags on a regular axis.
    max_separation
        Width of the tested window. ``None`` (the default) uses ``bins`` times
        the median nearest-neighbour spacing, which is what makes the regular
        case reduce exactly.
    max_pairs
        Pair budget; above it, pairs are subsampled uniformly at *seed*.
    seed
        Only used for that subsampling.

    Returns
    -------
    tuple
        Bin centres, ``rho_b``, and the number of pairs in each bin — for the
        **non-empty** bins only.

    Examples
    --------
    An anticorrelated saw-tooth on a regular axis: lag 1 negative, lag 2
    positive, exactly as the ordinary autocorrelation would say.

    >>> x = np.arange(40.0)
    >>> r = np.where(np.arange(40) % 2 == 0, 1.0, -1.0)
    >>> sep, rho, counts = separation_binned_autocorrelation(x, r, bins=2)
    >>> np.round(rho, 6).tolist()
    [-1.0, 1.0]
    >>> counts.tolist()
    [39, 38]
    """
    x = np.asarray(coordinates, dtype=float).ravel()
    r = np.asarray(residuals, dtype=float).ravel()
    if x.shape != r.shape:
        raise ResultsError(
            f"the whiteness statistic needs one residual per coordinate; got {x.size} "
            f"coordinate(s) and {r.size} residual(s)."
        )
    finite = np.isfinite(x) & np.isfinite(r)
    x, r = x[finite], r[finite]
    if x.size < 3:
        raise ResultsError(
            f"the whiteness statistic needs at least three unmasked, scored samples; this "
            f"dataset has {x.size}."
        )
    if bins < 1:
        raise ResultsError(f"a whiteness test needs at least one separation bin, got {bins}.")
    window = _default_window(x, bins) if max_separation is None else float(max_separation)
    if not np.isfinite(window) or window <= 0.0:
        raise ResultsError(f"max_separation must be finite and positive, got {max_separation!r}.")
    rng = np.random.default_rng(seed)
    left, right, separations, _ = _pairs_within(x, window, max_pairs=max_pairs, rng=rng)
    if left.size == 0:
        raise ResultsError(
            f"no pair of samples lies within {window:g} of another, so there is nothing to test "
            f"for correlation. Widen max_separation."
        )
    order, starts, counts = _bin_pairs(separations, bins, window)
    products = (r[left] * r[right])[order]
    variance = float(np.mean(r**2))
    if variance <= 0.0:
        raise ResultsError("every residual is exactly zero; there is no correlation to measure.")
    rho = np.add.reduceat(products, starts) / counts / variance
    centres = _bin_centres(separations[order], starts, counts)
    return centres, rho, counts


def _default_window(coordinates: np.ndarray, bins: int) -> float:
    """``bins`` times the median nearest-neighbour spacing.

    Chosen so that on an evenly spaced axis the bins *are* lags 1..``bins`` and
    the statistic reduces to the classical one, and so that on an irregular
    axis the window still adapts to how finely the data are actually sampled
    rather than to the arbitrary total span.
    """
    ordered = np.sort(coordinates)
    spacing = np.diff(ordered)
    spacing = spacing[spacing > 0.0]
    if spacing.size == 0:  # pragma: no cover - guarded by the caller's size check
        raise ResultsError("every coordinate is identical; there are no separations to bin.")
    return float(bins * np.median(spacing))


def _bin_centres(
    sorted_separations: np.ndarray, starts: np.ndarray, counts: np.ndarray
) -> np.ndarray:
    """Mean separation actually realised in each bin, not the nominal centre.

    On an irregular axis the pairs in a bin need not sit symmetrically inside
    it, and a plot's x-axis should say where the measurement came from.
    """
    return np.add.reduceat(sorted_separations, starts) / counts


def _statistic(rho: np.ndarray, counts: np.ndarray) -> float:
    """``Q = sum_b n_b rho_b^2`` — the pair-count-weighted aggregate."""
    return float(np.sum(counts * rho**2))


# ---------------------------------------------------------------------------
# Family B: the whiteness test
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class WhitenessTest:
    """The ``(statistic, p_value)`` record ``diagnostics.md`` §3.3 asks for.

    Deliberately **not** an :class:`~ampere.core.AnomalyScore`: §5 is explicit
    that family B declines the shared convention, because a statistic and a
    p-value are not a coordinate-indexed deficiency map and forcing them into
    that container shape would be exactly the conflation Tension 5 warns about.

    Attributes
    ----------
    dataset
        Which dataset was tested.
    statistic, p_value
        The pooled answer: the median across the tested posterior draws of the
        per-draw ``Q`` and of the per-draw permutation p-value. The median
        rather than the mean because a single badly-scored draw should not
        drag the summary, and because the per-draw arrays are kept anyway.
    separations, autocorrelation, pair_counts
        The binned autocorrelation itself — bin centres, the across-draw median
        ``rho_b``, and how many pairs each bin holds.
    statistics, p_values
        Per tested draw, so the spread is visible rather than summarised away.
    n_pairs, n_draws, n_permutations
        What the answer cost, and therefore how finely ``p_value`` is resolved:
        a permutation p-value cannot be smaller than ``1 / (1 + n_permutations)``.
    seed
        The sub-stream seed the permutations were drawn from, so the number is
        reproducible from the stored run alone.
    subsampled
        Whether the pair budget was hit and pairs were sampled rather than
        exhausted.
    """

    dataset: str
    statistic: float
    p_value: float
    separations: np.ndarray
    autocorrelation: np.ndarray
    pair_counts: np.ndarray
    statistics: np.ndarray
    p_values: np.ndarray
    n_pairs: int
    n_draws: int
    n_permutations: int
    seed: int
    subsampled: bool = False

    @property
    def resolution(self) -> float:
        """The smallest p-value this many permutations can report."""
        return 1.0 / (1.0 + self.n_permutations)

    def __repr__(self) -> str:
        return (
            f"<WhitenessTest {self.dataset!r}: Q={self.statistic:.4g}, p={self.p_value:.3g} "
            f"({self.n_draws} draw(s), {self.n_permutations} permutations)>"
        )


def residual_whiteness(
    tree: Any,
    *,
    dataset: str | None = None,
    bins: int = 8,
    max_separation: float | None = None,
    n_permutations: int = 199,
    max_draws: int = DEFAULT_MAX_DRAWS,
    max_pairs: int = DEFAULT_MAX_PAIRS,
    seed: int | None = None,
) -> WhitenessTest:
    """Test one dataset's signed residuals for leftover structure.

    The separation-binned, permutation-calibrated statistic this module's
    docstring defines, run on the ``residuals`` group
    (:func:`~ampere.results.derived.add_residuals`) — never on a resampled
    copy, and never on ``|residual|`` recovered from the per-observation
    log-likelihood, which has no sign (``diagnostics.md`` §3.3).

    A **small p-value means the residuals are not white**: there is structure
    at some separation the smooth model did not capture, which is the concrete
    trigger ``DEVELOPMENT_PLAN.md`` §4.8 wants in place of switching the
    flexible likelihood on as an act of faith.

    Parameters
    ----------
    tree
        A run carrying a ``residuals`` group.
    dataset
        Which dataset to test; required when the group holds more than one.
    bins, max_separation, max_pairs
        Passed to :func:`separation_binned_autocorrelation`.
    n_permutations
        Size of the permutation null. The p-value cannot resolve below
        ``1 / (1 + n_permutations)``.
    max_draws
        At most this many posterior draws are tested, taken evenly across the
        stored ones. Each costs a full permutation calibration.
    seed
        Overrides the run's own recorded seed. With neither, the permutations
        are entropy-seeded and the p-value will move between calls — the
        honest behaviour for a run that did not ask to be reproducible.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the ``residuals`` group is absent (naming
        :func:`~ampere.results.derived.add_residuals`), if the dataset is
        ambiguous, or if the coordinate axis is not 1-D.
    """
    group = require_group(
        tree,
        RESIDUALS_GROUP,
        remedy="call ampere.results.add_residuals(tree, problem) first, which derives the "
        "signed standardised residuals from the stored draws and the problem "
        "(results.md §7).",
    )
    available = [str(name) for name in group.data_vars]
    label = _one_dataset(available, dataset, RESIDUALS_GROUP)
    axis, coordinates = coordinate_of(group, label)
    values = np.asarray(group[label].values, dtype=float)
    flat = values.reshape(-1, values.shape[-1])
    scored = flat[np.any(np.isfinite(flat), axis=1)]
    if scored.size == 0:
        raise ResultsError(
            f"every stored draw of dataset {label!r} has unscored residuals, so there is nothing "
            f"to test. A run whose draws were all rejected by the prior looks like this."
        )
    if max_draws >= 1 and scored.shape[0] > max_draws:
        scored = scored[np.linspace(0, scored.shape[0] - 1, max_draws).astype(int)]
    resolved_seed = _resolve_seed(tree, seed)
    # One named sub-stream per dataset axis, so testing two datasets of one run
    # does not reuse a permutation sequence (``lowering.md`` §9.2's whole point).
    rng = generator(resolved_seed, f"{WHITENESS_STREAM}.{axis}")

    finite = np.all(np.isfinite(scored), axis=0)
    x = coordinates[finite]
    if x.size < 3:
        raise ResultsError(
            f"dataset {label!r} has {x.size} sample(s) scored in every tested draw; the "
            f"whiteness statistic needs at least three."
        )
    window = _default_window(x, bins) if max_separation is None else float(max_separation)
    left, right, separations, subsampled = _pairs_within(x, window, max_pairs=max_pairs, rng=rng)
    if left.size == 0:
        raise ResultsError(
            f"no pair of dataset {label!r}'s samples lies within {window:g} of another, so there "
            f"is nothing to test for correlation. Widen max_separation."
        )
    order, starts, counts = _bin_pairs(separations, bins, window)
    left, right = left[order], right[order]
    centres = _bin_centres(separations[order], starts, counts)

    rows = scored[:, finite]
    per_draw_rho = np.empty((rows.shape[0], counts.size), dtype=float)
    per_draw_q = np.empty(rows.shape[0], dtype=float)
    per_draw_p = np.empty(rows.shape[0], dtype=float)
    for index, residual in enumerate(rows):
        rho = _rho_of(residual[np.newaxis, :], left, right, starts, counts)[0]
        observed = _statistic(rho, counts)
        null = _permutation_statistics(residual, left, right, starts, counts, n_permutations, rng)
        per_draw_rho[index] = rho
        per_draw_q[index] = observed
        per_draw_p[index] = (1.0 + float(np.count_nonzero(null >= observed))) / (
            1.0 + n_permutations
        )
    return WhitenessTest(
        dataset=label,
        statistic=float(np.median(per_draw_q)),
        p_value=float(np.median(per_draw_p)),
        separations=centres,
        autocorrelation=np.median(per_draw_rho, axis=0),
        pair_counts=counts,
        statistics=per_draw_q,
        p_values=per_draw_p,
        n_pairs=int(left.size),
        n_draws=int(rows.shape[0]),
        n_permutations=int(n_permutations),
        seed=resolved_seed,
        subsampled=subsampled,
    )


def _rho_of(
    residuals: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    starts: np.ndarray,
    counts: np.ndarray,
) -> np.ndarray:
    """``rho_b`` for a stack of residual vectors, all bins at once.

    ``residuals`` is ``(m, n)``; the result is ``(m, n_bins)``. The stack is
    what lets the permutation null be computed in batches rather than one
    shuffle at a time.
    """
    products = residuals[:, left] * residuals[:, right]
    variance = np.mean(residuals**2, axis=1)
    return np.add.reduceat(products, starts, axis=1) / counts / variance[:, np.newaxis]


def _permutation_statistics(
    residual: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
    starts: np.ndarray,
    counts: np.ndarray,
    n_permutations: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """``Q`` under the exchangeability null, ``n_permutations`` times.

    Shuffling the residual **values** across fixed coordinates is the right
    null for this question: it destroys any relationship between residual and
    position while keeping both the residual marginal and the sampling
    geometry — the two things that make ``Q``'s distribution untabulated in the
    first place — exactly as observed.

    Done in batches so that a 10^5-point spectrum's pair list does not have to
    be multiplied by ``n_permutations`` all at once.
    """
    if n_permutations < 1:
        raise ResultsError(
            f"the whiteness test is calibrated by permutation, so it needs at least one "
            f"permutation; got {n_permutations}."
        )
    batch = max(1, int(4_000_000 // max(left.size, 1)))
    out = np.empty(n_permutations, dtype=float)
    done = 0
    while done < n_permutations:
        size = min(batch, n_permutations - done)
        shuffled = np.stack([rng.permutation(residual) for _ in range(size)])
        rho = _rho_of(shuffled, left, right, starts, counts)
        out[done : done + size] = np.sum(counts * rho**2, axis=1)
        done += size
    return out


def _one_dataset(available: list[str], dataset: str | None, group: str) -> str:
    chosen = select_datasets(available, None if dataset is None else [dataset], what=group)
    if len(chosen) != 1:
        raise ResultsError(
            f"this run's {group!r} group holds {sorted(available)}; name one with dataset=. "
            f"A whiteness statistic pooled across datasets with different coordinate axes would "
            f"not mean anything."
        )
    return chosen[0]


def _resolve_seed(tree: Any, seed: int | None) -> int:
    """The base seed the permutation stream derives from.

    An explicit argument wins; otherwise the run's own recorded seed, so that
    re-running the check on a stored run reproduces its p-value. An
    entropy-seeded run gets an entropy base — the honest behaviour for a run
    that did not ask to be reproducible — and the base is recorded on the
    result so the number can be reproduced after the fact anyway.

    Note what this deliberately does **not** do: derive anything from
    :func:`hash`. Python's built-in hash is salted per process, so a
    label-derived seed built on it would differ between two runs of the same
    script — the exact irreproducibility :mod:`ampere.core.rng` exists to
    prevent. The label goes to :func:`~ampere.core.rng.substream` instead.
    """
    if seed is not None:
        return int(seed)
    recorded = run_seed(tree)
    if recorded is not None:
        return int(recorded)
    return int((np.random.SeedSequence().entropy or 0) % (2**31))


# ---------------------------------------------------------------------------
# Family B: the posterior-predictive half
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class ChiSquareCheck:
    """A posterior-predictive p-value from the stored per-draw log-likelihood.

    ``diagnostics.md`` §3.3: §4.6's log-likelihood requirement "makes the
    posterior-predictive discrepancy statistics cheap already — no new
    ``InferenceData`` field, no re-running the forward model". This is that,
    and only that.

    Attributes
    ----------
    dataset
        Which dataset was checked.
    p_value
        ``mean_k Pr(chi^2_n >= T(y, theta_k))``, the Bayesian p-value for the
        chi-square discrepancy. Values near 0 mean the fit is worse than the
        model can explain; values near 1 mean it is *better*, which usually
        means the uncertainties are overstated.
    statistic
        The posterior median of ``T``.
    statistics
        ``T`` at every scored draw.
    degrees_of_freedom
        The retained sample count, which is the reference distribution's
        parameter. Deliberately **not** reduced by the number of fitted
        parameters: a posterior-predictive p-value compares the discrepancy at
        each drawn theta against the sampling distribution *at that theta*,
        where nothing has been fitted away.
    n_draws
        How many draws contributed.
    """

    dataset: str
    p_value: float
    statistic: float
    statistics: np.ndarray
    degrees_of_freedom: int
    n_draws: int

    def __repr__(self) -> str:
        return (
            f"<ChiSquareCheck {self.dataset!r}: T={self.statistic:.4g} on "
            f"{self.degrees_of_freedom} d.o.f., p={self.p_value:.3g}>"
        )


def chi_square_pvalue(tree: Any, *, dataset: str | None = None) -> ChiSquareCheck:
    """The cheap half of family B: a Bayesian p-value with no replicate draws.

    For a Gaussian family with independent noise and no noise nuisance
    parameters, the stored per-dataset log-likelihood is exactly

    .. code-block:: text

        log L(theta) = -T(theta)/2 - (n/2) log(2 pi) - sum_i log sigma_i

    so the chi-square discrepancy ``T(theta) = sum_i r_i(theta)^2`` is
    recoverable from the run alone — the uncertainties are in ``constant_data``
    and the mask with them — and its reference distribution under the model is
    ``chi^2`` on the retained count. No ``y_rep`` is needed, which is why
    ``results.md`` §7 can keep posterior-predictive replicates out of the
    default emission without making this check expensive.

    Anything outside that case is refused by name rather than approximated:
    a GP likelihood does not factorise, a fitted ``scale``/``jitter`` moves
    ``sigma`` per draw, censoring changes the density, and a non-Gaussian
    family has a different reference distribution altogether. For those,
    replicate draws are the honest route —
    :func:`~ampere.results.derived.add_posterior_predictive`.

    Raises
    ------
    ampere.core.exceptions.ResultsError
        If the run stores no per-dataset log-likelihood, if the dataset is
        ambiguous, or if its likelihood is outside the case above.
    """
    group = require_group(
        tree,
        LOG_LIKELIHOOD_GROUP,
        remedy="every run emitted by ampere.results.emit carries it, so this tree was not "
        "produced by emit().",
    )
    available = [str(name) for name in group.data_vars]
    label = _one_dataset(available, dataset, LOG_LIKELIHOOD_GROUP)
    _check_chi_square_applies(tree, label)
    sigma, retained = _stored_sigma(tree, label)
    n = int(retained)
    constant = 0.5 * n * float(np.log(2.0 * np.pi)) + float(np.sum(np.log(sigma)))
    log_likelihood = np.asarray(group[label].values, dtype=float).ravel()
    scored = log_likelihood[np.isfinite(log_likelihood)]
    if scored.size == 0:
        raise ResultsError(
            f"no draw of dataset {label!r} was scored (every stored log-likelihood is NaN or "
            f"-inf), so there is no discrepancy to check."
        )
    statistics = -2.0 * (scored + constant)
    p_values = st.chi2.sf(np.clip(statistics, 0.0, None), df=n)
    return ChiSquareCheck(
        dataset=label,
        p_value=float(np.mean(p_values)),
        statistic=float(np.median(statistics)),
        statistics=statistics,
        degrees_of_freedom=n,
        n_draws=int(scored.size),
    )


def _check_chi_square_applies(tree: Any, label: str) -> None:
    """Refuse, by name, every case where the log-likelihood is not an affine T."""
    spec = likelihood_specs(tree).get(label)
    if spec is None:
        raise ResultsError(
            f"this run records no likelihood declaration for dataset {label!r}, so there is no "
            f"way to tell whether its log-likelihood is a chi-square in disguise."
        )
    family = spec.get("family")
    noise = spec.get("noise")
    parameters = spec.get("parameters", {})
    declared = parameters.get("parameters", []) if isinstance(parameters, dict) else []
    reason = None
    if family != "gaussian":
        reason = (
            f"its family is {family!r}, whose sampling distribution for the chi-square "
            f"discrepancy is not chi-square"
        )
    elif noise != "IndependentNoise":
        reason = (
            f"its noise model is {noise!r}: a GP likelihood does not factorise over "
            f"observations, so the stored value is a joint marginal and not a sum of squares "
            f"(likelihoods.md §16). Family B is scoped away from GP fits in any case "
            f"(diagnostics.md §3.1) — use plot_gp_localisation"
        )
    elif declared:
        names = sorted(str(entry.get("name")) for entry in declared if isinstance(entry, dict))
        reason = (
            f"its noise model declares {names}, so sigma moves from draw to draw and the "
            f"log-likelihood's additive constant is not a constant"
        )
    elif spec.get("censoring") is not None:
        reason = "it declares censored samples, whose density is not the Gaussian one"
    if reason is not None:
        raise ResultsError(
            f"dataset {label!r}'s log-likelihood cannot be turned into a chi-square "
            f"discrepancy: {reason}. Draw posterior-predictive replicates instead — "
            f"ampere.results.add_posterior_predictive — which samples through the family the "
            f"likelihood actually scores with (inference.md §13)."
        )


def _stored_sigma(tree: Any, label: str) -> tuple[np.ndarray, int]:
    """The retained uncertainties of one dataset, from ``constant_data``.

    Masked samples are dropped rather than inflated to ``inf``: this is the
    normalising constant of a likelihood that never saw them, so they
    contribute neither a term nor a degree of freedom.
    """
    constant = require_group(
        tree,
        "constant_data",
        remedy="emit the run with observed=True (the default), which stores each dataset's "
        "uncertainties and mask.",
    )
    name = f"{label}_uncertainty"
    if name not in constant.variables:
        raise ResultsError(
            f"this run stores no uncertainties for dataset {label!r}, so the chi-square "
            f"discrepancy's normalisation is unknown. Attach uncertainties to the observed "
            f"container, or draw replicates with add_posterior_predictive."
        )
    sigma = np.asarray(constant[name].values, dtype=float).ravel()
    mask_name = f"{label}_mask"
    if mask_name in constant.variables:
        excluded = np.asarray(constant[mask_name].values).ravel().astype(bool)
        sigma = sigma[~excluded]
    if sigma.size == 0 or not np.all(np.isfinite(sigma)) or np.any(sigma <= 0.0):
        raise ResultsError(
            f"dataset {label!r} has no positive finite uncertainties left after masking, so the "
            f"chi-square discrepancy is undefined."
        )
    return sigma, int(sigma.size)


# ---------------------------------------------------------------------------
# Family C: the shared score
# ---------------------------------------------------------------------------


def gp_localisation_score(
    tree: Any,
    *,
    dataset: str | None = None,
    standardise: bool = True,
) -> AnomalyScore:
    """The conditioned GP mean, as the shared :class:`~ampere.core.AnomalyScore`.

    ``diagnostics.md`` §4.2 asks family C for "an aggregated anomaly score...
    coordinate-indexed, amplitude-weighted magnitude of the local deviation,
    for the same shared presentation family A uses". That is what this is: the
    magnitude of the posterior median conditioned GP mean, divided (by default)
    by the total posterior standard deviation of that mean, so the score reads
    as **how many sigma of GP the fit needed here** — higher is worse, on a
    range comparable between datasets and between fits, which is what §5's
    shared visual grammar requires of both families' scores.

    The denominator is the law of total variance — the mean of the conditional
    variances plus the variance of the conditional means across draws — and not
    the across-draw spread alone. The spread alone is the wrong quantity and
    dangerously so: a fit whose GP hyperparameters are tightly constrained has
    almost none of it, and dividing by a near-zero, structured number turns a
    clean localisation into noise. Both terms are stored, so nothing has to be
    recomputed to combine them.

    The score is a **magnitude**: the container's documented sense is
    "higher = more anomalous", and a signed quantity cannot carry that. The
    sign is not lost — it is what :func:`~ampere.results.plots.plot_gp_localisation`
    draws, and it is in the ``gp_localisation`` group this reads.

    ``interpretation_notes`` carries :data:`~ampere.results.plots.GP_LOCALISATION_CAVEAT`
    by construction. ``diagnostics.md`` §4.3 requires exactly that: the caveat
    must reach "a user extracting the score programmatically (not just viewing
    the plot)".

    Parameters
    ----------
    tree
        A run carrying a ``gp_localisation`` group.
    dataset
        Which dataset's localisation to convert; required when there is more
        than one.
    standardise
        Divide by the total posterior standard deviation. ``False`` leaves the
        score in the data's own units, which is comparable across coordinates
        but not across datasets.
    """
    from .plots import GP_LOCALISATION_CAVEAT

    group = require_group(
        tree,
        GP_LOCALISATION_GROUP,
        remedy="call ampere.results.gp_localisation(tree, problem) first, which evaluates "
        "Likelihood.conditional across the stored draws (results.md §7).",
    )
    available = sorted({str(name).rsplit("_", 1)[0] for name in group.data_vars})
    label = _one_dataset(available, dataset, GP_LOCALISATION_GROUP)
    _, coordinates = coordinate_of(group, f"{label}_mean")
    means = np.asarray(group[f"{label}_mean"].values, dtype=float)
    variances = np.asarray(group[f"{label}_variance"].values, dtype=float)
    flat_mean = means.reshape(-1, means.shape[-1])
    flat_variance = variances.reshape(-1, variances.shape[-1])
    keep = np.any(np.isfinite(flat_mean), axis=1)
    if not np.any(keep):
        raise ResultsError(
            f"no stored draw localises dataset {label!r}: every conditioned mean is NaN."
        )
    flat_mean, flat_variance = flat_mean[keep], flat_variance[keep]
    values = np.abs(np.nanmedian(flat_mean, axis=0))
    if standardise:
        total = np.nanmean(flat_variance, axis=0) + np.nanvar(flat_mean, axis=0)
        spread = np.sqrt(np.clip(total, 0.0, None))
        positive = spread > 0.0
        if not np.any(positive):
            raise ResultsError(
                f"dataset {label!r}'s conditioned GP has zero posterior variance everywhere, so "
                f"a standardised score would divide by zero. Pass standardise=False."
            )
        values = values / np.where(positive, spread, float(np.median(spread[positive])))
    mask = ~np.isfinite(values)
    return AnomalyScore(
        coordinates=coordinates,
        values=np.where(mask, 0.0, values),
        provenance=GP_LOCALISATION_PROVENANCE,
        interpretation_notes=(
            f"Dataset {label!r}: the magnitude of the posterior median conditioned GP mean"
            f"{', in units of its own total posterior standard deviation' if standardise else ''}"
            f". {GP_LOCALISATION_CAVEAT}"
        ),
        mask=mask if bool(np.any(mask)) else None,
    )


def gp_localisation_datasets(tree: Any) -> tuple[str, ...]:
    """Labels this run fitted with a GP — family C's applicable datasets."""
    return gp_datasets(tree)
