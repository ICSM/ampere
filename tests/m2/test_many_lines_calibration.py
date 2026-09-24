"""W5.8: does the warped fit's *credible interval* survive the misspecification?

``tests/m2/test_many_lines.py`` asks whether the warped Matérn recovers the
truth on one spectrum. This module asks the harder question, the one W4.2 asked
of visibilities and W3.6's :func:`~ampere.results.calibration.sbc` exists to
answer: over many spectra drawn from the prior and deviated the same way, do
the intervals the fit reports contain the truth as often as they claim?

It matters here more than anywhere else in the study, because W5.7's warp is
the one piece of the flexible likelihood with real freedom in it — six extra
sampled dimensions whose only guard is a shrinkage prior towards the identity
warp. Freedom that is free gets used to absorb signal, and a fit that absorbs
signal reports intervals that are too *narrow*. The plan's bullet names the
check by name — "SBC on injected misspecification" — and this is it.

What is injected
----------------
:class:`LineForestError`, a one-step instrument that multiplies the model flux
by ``1 + delta(lambda)`` with *delta* W5.8's own
:data:`~examples.m2_misspecification.generators.MANY_LINES` deviation: five
narrow lines in one band plus a smooth continuum error. It carries **no
parameters** — the deviation is a fact about the experiment, not something the
fit is allowed to learn — and it is present in the simulating problem and
absent from the fitted one, which is exactly what makes every fit misspecified.
Deliberately not a draw from the kernel the GP is then given: a systematic
drawn from the fitted covariance function would make the fit correctly
specified, and "a correctly specified model is calibrated" is a statement about
SBC rather than about ampere.

What is asserted, and what is reported
--------------------------------------
The assertion is the directional one ``test_visibility_calibration.py``
established, and it is a **comparison**: the standard likelihood's
central-:data:`NOMINAL` intervals contain the truth far less often than they
claim, the warped fit's do not, and the gap between them is wide. Coverage is
what is pinned and the uniformity p-value is not, for that module's reason —
over-coverage is not a calibration failure in the sense that matters for a
published error bar, since an interval that contains the truth more often than
it claims does not invite a false conclusion and one that contains it less
often does.

Why the comparison rather than "the warped fit covers at least
:data:`NOMINAL`": at :data:`SIMULATIONS` refits a coverage estimate moves in
steps of a tenth and its binomial standard error is about 0.095, so a fit that
is *exactly* calibrated returns 8 out of 10 more than a quarter of the time.
An assertion that such a fit must return 9 or 10 is an assertion about the
budget. The gap against the standard likelihood is not: it is the difference
the experiment is about, it is enormous, and it does not move with the count.
What the warped arm is additionally held to is that it does not fall more than
:data:`COVERAGE_ALLOWANCE` — two of those standard errors — below nominal,
which is the *under*-coverage this row exists to rule out.

The budget is small: :data:`SIMULATIONS` refits at the per-PR level, where
Talts et al. (2018) want a hundred or more, so these numbers are evidence of a
gross effect and nothing finer. ``sbc`` warns about exactly that and the
warning is allowed through. The ``m2_full`` row pays for
:data:`FULL_SIMULATIONS` and adds the **stationary** flexible arm beside the
other two, which is the three-way comparison this budget cannot afford.

Cost: about two minutes on the reference backend.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import pytest

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Instrument,
    Likelihood,
    Spectrum,
    Transformation,
)
from ampere.results.calibration import sbc
from examples.m2_misspecification import study
from examples.m2_misspecification.generators import MANY_LINES, generate
from examples.m2_misspecification.model import FLUX_UNIT, WAVELENGTH_UNIT

#: W5.26: every row here refits on the reference backend and asserts a claim
#: about the likelihood's calibration, not about the array library;
#: ``tests/m2/conftest.py`` skips ``study``-marked rows outside ``dev``. (The
#: already-opt-in ``m2_full`` row below carries both markers; either one
#: alone already keeps it out of the default gate.)
pytestmark = pytest.mark.study

#: Points per simulated spectrum. Small, and it can be: the claim is about
#: coverage over many refits, not about one posterior's width.
SIZE = 150

#: Refits at the per-PR level, and behind ``-m m2_full``.
SIMULATIONS = 10
FULL_SIMULATIONS = 24

#: Walkers per arm. Not one number, because an emcee ensemble needs more than
#: twice the dimension of them and the three arms have four, six and twelve
#: free parameters: one number would either starve the warped arm or spend
#: three times what the standard one needs on every refit. The *steps* are the
#: same for all three, which is the part of the budget a comparison of noise
#: models has to hold fixed.
WALKERS = {"standard": 16, "matern32": 16, "warped": 26}
STEPS = 500
BURN_IN = 250

#: Posterior draws each rank is taken against.
RANK_DRAWS = 120

#: The nominal level the coverage is read at. On
#: :data:`~ampere.results.calibration.DEFAULT_LEVELS`' grid of 21 points.
NOMINAL = 0.9

#: The parameters ranked: the continuum normalisation, which the smooth error
#: damages most, and the second line depth, which the forest sits nearest.
RANKED = ("model.A", "model.d2")

#: How much wider the warped fit's coverage must be than the standard
#: likelihood's. Measured: 0.80 against 0.00 on the worst parameter, so the gap
#: is 0.80 against a threshold of half that.
COVERAGE_MARGIN = 0.4

#: How far below :data:`NOMINAL` the warped fit's coverage may fall before it
#: counts as under-covering: two binomial standard errors at
#: :data:`SIMULATIONS` refits, rounded up. Measured: 0.80 and 0.90.
COVERAGE_ALLOWANCE = 0.2


class LineForestError(Transformation):
    """W5.8's deviation as an instrument step: ``flux -> flux * (1 + delta)``.

    A buffer rather than a parameter, and no parameters at all: the deviation
    is a fact about the experiment, and a fit that could learn it would not be
    misspecified.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)
    PRODUCES: ClassVar[type] = Spectrum

    def __init__(self, deviation: np.ndarray, *, label: str | None = None) -> None:
        super().__init__(label=label)
        self.register_buffer("deviation", np.asarray(deviation, dtype=float))

    def apply(self, samples: Any, values: Mapping[str, Any] | None) -> Spectrum:
        context = self.context(values)
        factor = 1.0 + np.asarray(context["deviation"], dtype=float)
        return samples.with_values(np.asarray(samples.values) * factor)


def _template() -> Any:
    """An empty spectrum carrying the grid and the per-point uncertainty."""
    data = generate(MANY_LINES, size=SIZE)
    return data, Spectrum(
        data.wavelength * WAVELENGTH_UNIT,
        np.zeros(data.size) * FLUX_UNIT,
        uncertainty=data.uncertainty * FLUX_UNIT,
    )


def simulating_problem() -> FittingProblem:
    """What the data are drawn from: the model **plus** the line-forest deviation."""
    data, template = _template()
    return FittingProblem(
        study.model_for("reference", data.wavelength),
        [
            Dataset(
                template,
                Instrument([LineForestError(MANY_LINES.deviation(data.wavelength))]),
                likelihood=Likelihood(GaussianFamily(), IndependentNoise()),
                label=study.DATASET_LABEL,
            )
        ],
        seed=20260918,
    )


def fitted_problem(observed: Any, *, arm: str) -> FittingProblem:
    """The model alone against *observed*, under one of the three likelihoods.

    No instrument, which is the whole point: the fit has no parameter for the
    deviation and cannot produce it at any point of its parameter space.
    """
    likelihood = (
        study.build_likelihood("standard")
        if arm == "standard"
        else study.build_likelihood("flexible", kernel=arm, solver="quasisep")
    )
    grid = np.asarray(observed.axes[0].values, dtype=float)
    return FittingProblem(
        study.model_for("reference", grid),
        [Dataset(observed, likelihood=likelihood, label=study.DATASET_LABEL)],
        seed=20260918,
    )


def run_sbc(*, arm: str, count: int) -> Any:
    """``sbc`` over *count* refits of the misspecified fit with this likelihood."""
    from ampere.inference import EmceeEngine

    def factory(replica: FittingProblem) -> Any:
        observed = replica.datasets[study.DATASET_LABEL].observed
        return EmceeEngine(fitted_problem(observed, arm=arm), walkers=WALKERS[arm])

    with pytest.warns(UserWarning, match="goodness-of-fit test"):
        return sbc(
            simulating_problem(),
            factory,
            count=count,
            draws=RANK_DRAWS,
            run_options={"steps": STEPS, "burn_in": BURN_IN},
            parameters=list(RANKED),
            seed=515151,
            label=f"many_lines, {arm} fit",
        )


def coverage_at_nominal(calibration: Any) -> np.ndarray:
    """The empirical coverage of the central :data:`NOMINAL` interval, per parameter."""
    levels = np.asarray(calibration["level"].values)
    index = int(np.argmin(np.abs(levels - NOMINAL)))
    assert abs(float(levels[index]) - NOMINAL) < 1e-9, "the nominal level is not on the grid"
    return np.asarray(calibration["coverage"].values)[index]


# ---------------------------------------------------------------------------
# The experiment is the experiment it says it is.
# ---------------------------------------------------------------------------


class TestTheInjectionIsWhatItClaims:
    def test_the_simulating_problem_carries_the_deviation_and_the_fit_does_not(self) -> None:
        simulating = simulating_problem()
        step = simulating.datasets[study.DATASET_LABEL].instrument.steps[0]
        assert isinstance(step, LineForestError)
        assert step.parameters.free_size == 0
        _, template = _template()
        fitted = fitted_problem(template, arm="warped")
        assert fitted.datasets[study.DATASET_LABEL].instrument.steps == ()

    def test_the_injected_deviation_dominates_the_noise(self) -> None:
        """If it did not, both fits would be calibrated and the row would say nothing."""
        data, _ = _template()
        injected = np.abs(data.deviated - data.truth)
        assert float(np.max(injected)) > 5.0 * float(np.median(data.uncertainty))


# ---------------------------------------------------------------------------
# The calibration claim.
# ---------------------------------------------------------------------------


class TestTheWarpedFitsIntervalsHoldUp:
    @pytest.fixture(scope="class")
    @classmethod
    def rigid(cls) -> Any:
        return run_sbc(arm="standard", count=SIMULATIONS)

    @pytest.fixture(scope="class")
    @classmethod
    def warped(cls) -> Any:
        return run_sbc(arm="warped", count=SIMULATIONS)

    def test_the_standard_likelihood_under_covers(self, rigid: Any) -> None:
        """The baseline the comparison is against, and it is not close.

        An unmodelled deviation of this size under independent noise puts the
        truth outside the 90 % interval essentially every time. Measured: 0.00
        on both parameters.
        """
        measured = coverage_at_nominal(rigid)
        print(f"\nstandard coverage at {NOMINAL}: {measured}")
        assert np.all(measured < NOMINAL)

    def test_the_warped_fit_does_not_under_cover(self, warped: Any) -> None:
        """The failure mode this row exists to rule out.

        Six sampled warp dimensions under a shrinkage prior towards the
        identity: the risk is that the warp absorbs signal and the intervals
        come back too *narrow*. Measured: 0.80 and 0.90, which at ten
        simulations is what an exactly calibrated fit returns — so what is
        pinned is that it does not fall more than
        :data:`COVERAGE_ALLOWANCE` below nominal, not that it exceeds it.
        """
        measured = coverage_at_nominal(warped)
        print(f"\nwarped coverage at {NOMINAL}: {measured}")
        assert np.all(measured >= NOMINAL - COVERAGE_ALLOWANCE)

    def test_the_gap_between_them_is_wide(self, rigid: Any, warped: Any) -> None:
        """The pinned comparison, with :data:`COVERAGE_MARGIN` of slack.

        The margin rather than two absolute thresholds, for the reason W4.5's
        fringing row gives: it is the difference the experiment is about, and
        it is the quantity that does not move when the budget does.
        """
        gap = float(np.min(coverage_at_nominal(warped)) - np.min(coverage_at_nominal(rigid)))
        print(f"\ncoverage gap (warped - standard) at {NOMINAL}: {gap:.3f}")
        assert gap >= COVERAGE_MARGIN

    def test_every_simulation_produced_a_usable_fit(self, rigid: Any, warped: Any) -> None:
        """A study that silently dropped half its refits would still pass the rows above."""
        for calibration in (rigid, warped):
            assert int(calibration.attrs["ampere_calibration_failures"]) == 0
            assert calibration["ranks"].shape == (SIMULATIONS, len(RANKED))


@pytest.mark.m2_full
def test_the_three_arms_calibration_at_the_marked_budget() -> None:
    """The comparison the per-PR budget cannot afford: three arms, more simulations.

    The stationary Matérn is the arm ``test_many_lines.py`` shows is biased on
    this scenario, so this is the question that bias raises — does the bias
    reach the interval? — asked over many spectra rather than one, at a count
    where a coverage estimate has some resolution.
    """
    measured = {
        arm: coverage_at_nominal(run_sbc(arm=arm, count=FULL_SIMULATIONS))
        for arm in ("standard", "matern32", "warped")
    }
    for arm, coverage in measured.items():
        print(f"\n{arm} coverage at {NOMINAL}: {coverage}")
    assert np.all(measured["standard"] < NOMINAL)
    assert np.min(measured["warped"]) - np.min(measured["standard"]) >= COVERAGE_MARGIN
    assert np.min(measured["warped"]) >= np.min(measured["matern32"])
    assert np.all(measured["warped"] >= NOMINAL - COVERAGE_ALLOWANCE)
