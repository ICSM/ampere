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
established: the warped fit's central-:data:`NOMINAL` intervals contain the
truth **at least** as often as they claim. Over-coverage is not a calibration
failure in the sense that matters for a published error bar — an interval that
contains the truth more often than it claims does not invite a false
conclusion, and one that contains it less often does — so coverage is what is
pinned and the uniformity p-value is not.

The budget is small: :data:`SIMULATIONS` refits at the per-PR level, where
Talts et al. (2018) want a hundred or more, so these numbers are evidence of a
gross effect and nothing finer. ``sbc`` warns about exactly that and the
warning is allowed through. The ``m2_full`` row pays for
:data:`FULL_SIMULATIONS` and adds the stationary arm beside it, which is the
comparison this budget cannot afford.

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

#: Points per simulated spectrum. Small, and it can be: the claim is about
#: coverage over many refits, not about one posterior's width.
SIZE = 150

#: Refits at the per-PR level, and behind ``-m m2_full``.
SIMULATIONS = 10
FULL_SIMULATIONS = 24

#: The emcee budget each refit gets. Twenty-six walkers because the warped arm
#: has twelve free parameters and an ensemble needs more than twice the
#: dimension of them; the stationary arm is run at the same budget so that the
#: ``m2_full`` comparison is a comparison of noise models.
WALKERS = 26
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

#: How much wider the warped fit's coverage must be than the stationary arm's,
#: in the ``m2_full`` comparison.
COVERAGE_MARGIN = 0.0


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


def fitted_problem(observed: Any, *, kernel: str) -> FittingProblem:
    """The model alone against *observed*, under the flexible likelihood with *kernel*.

    No instrument, which is the whole point: the fit has no parameter for the
    deviation and cannot produce it at any point of its parameter space.
    """
    likelihood = study.build_likelihood("flexible", kernel=kernel, solver="quasisep")
    grid = np.asarray(observed.axes[0].values, dtype=float)
    return FittingProblem(
        study.model_for("reference", grid),
        [Dataset(observed, likelihood=likelihood, label=study.DATASET_LABEL)],
        seed=20260918,
    )


def run_sbc(*, kernel: str, count: int) -> Any:
    """``sbc`` over *count* refits of the misspecified fit with this kernel."""
    from ampere.inference import EmceeEngine

    def factory(replica: FittingProblem) -> Any:
        observed = replica.datasets[study.DATASET_LABEL].observed
        return EmceeEngine(fitted_problem(observed, kernel=kernel), walkers=WALKERS)

    with pytest.warns(UserWarning, match="goodness-of-fit test"):
        return sbc(
            simulating_problem(),
            factory,
            count=count,
            draws=RANK_DRAWS,
            run_options={"steps": STEPS, "burn_in": BURN_IN},
            parameters=list(RANKED),
            seed=515151,
            label=f"many_lines, {kernel} flexible fit",
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
        fitted = fitted_problem(template, kernel="warped")
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
    def warped(cls) -> Any:
        return run_sbc(kernel="warped", count=SIMULATIONS)

    def test_the_warped_fit_does_not_under_cover(self, warped: Any) -> None:
        """The directional claim, and the only one this budget can support.

        Six sampled warp dimensions under a shrinkage prior towards the
        identity: the failure mode worth ruling out is that the warp absorbs
        signal and the intervals come back too narrow. They do not.
        """
        measured = coverage_at_nominal(warped)
        print(f"\nwarped coverage at {NOMINAL}: {measured}")
        assert np.all(measured >= NOMINAL)

    def test_every_simulation_produced_a_usable_fit(self, warped: Any) -> None:
        """A study that silently dropped half its refits would still pass the row above."""
        assert int(warped.attrs["ampere_calibration_failures"]) == 0
        assert warped["ranks"].shape == (SIMULATIONS, len(RANKED))


@pytest.mark.m2_full
def test_the_warped_fit_covers_at_least_as_well_as_the_stationary_one() -> None:
    """The comparison the per-PR budget cannot afford: both arms, more simulations.

    The stationary arm is the one ``test_many_lines.py`` shows is biased on this
    scenario, so this is the question that bias raises — does the bias reach the
    interval? — asked over many spectra rather than one.
    """
    stationary = run_sbc(kernel="matern32", count=FULL_SIMULATIONS)
    warped = run_sbc(kernel="warped", count=FULL_SIMULATIONS)
    stationary_coverage = coverage_at_nominal(stationary)
    warped_coverage = coverage_at_nominal(warped)
    print(
        f"\nat {NOMINAL}: stationary {stationary_coverage}, warped {warped_coverage}",
    )
    gap = float(np.min(warped_coverage) - np.min(stationary_coverage))
    assert gap >= COVERAGE_MARGIN
    assert np.all(warped_coverage >= NOMINAL)
