"""Milestone M2's science claim, as assertions with thresholds.

The claim, stated so that it can fail:

1. **The control works.** With no misspecification, the two likelihoods agree,
   both recover the truth, the GP's amplitude collapses towards zero, and the
   residual-whiteness test says the residuals are white. Without this row the
   rest is not evidence — a GP that always inflated the posterior would pass
   every other assertion here for the wrong reason.
2. **Under misspecification the standard likelihood is confidently wrong.** Its
   68 % interval stops containing the truth, and on at least one parameter the
   median is :data:`~examples.m2_misspecification.study.STANDARD_MIN_BIAS_WIDTHS`
   posterior widths away from it.
3. **The flexible likelihood stays honest.** Every parameter, in every
   scenario, lands within
   :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS`
   posterior widths of the truth, and its 68 % interval contains it. Note what
   is *not* claimed: not that the flexible median is closer to the truth in the
   parameter's own units — sometimes it is not — but that the uncertainty it
   reports is one a reader can believe.
4. **The diagnostics say so too.** Residual whiteness finds structure in the
   misspecified standard fits and none in the control; GP localisation puts its
   peak where a *localised* deviation was injected, and does so with a score
   far above the control's.

Budget and margins
------------------
Everything here runs at :data:`~examples.m2_misspecification.study.TEST_EMCEE`
(20 walkers, 900 steps, 450 discarded) at 200 points — about eighty seconds for
all eight runs, paid once by the session fixture. The thresholds are set well
clear of the Monte Carlo error that budget leaves: the measured run-to-run
scatter of a posterior median is up to 0.31 posterior widths there
(:data:`~examples.m2_misspecification.study.MEASURED_MEDIAN_SCATTER`), while
the assertions have margins of 0.5 widths (claim 3: threshold 1.5, worst
observed 0.98) and 1.2 widths (claim 2: threshold 2.0, worst observed 3.16).
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from examples.m2_misspecification import study
from examples.m2_misspecification.generators import scenario_named, wavelength_grid

MISSPECIFIED = ("mild", "strong_smooth", "strong_sharp")
AMPLITUDE = f"{study.DATASET_LABEL}.likelihood.amplitude"
LENGTH_SCALE = f"{study.DATASET_LABEL}.likelihood.length_scale"


def _summaries(results: Any, key: str, kind: str) -> dict[str, study.Summary]:
    return study.summarise(results[(key, kind)]["run"], names=study.PHYSICAL_NAMES)


# ---------------------------------------------------------------------------
# 1. The control
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["standard", "flexible"])
def test_without_misspecification_both_likelihoods_recover_the_truth(
    study_results: Any, kind: str
) -> None:
    summaries = _summaries(study_results, "none", kind)
    for name, summary in summaries.items():
        assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, f"{kind} {name}"
        assert summary.covers_truth, f"{kind} {name}"


def test_without_misspecification_the_gp_switches_itself_off(study_results: Any) -> None:
    """The control that makes the other three scenarios mean something.

    If the GP inflated the posterior whether or not there was anything to
    absorb, the comparison would be between a wide answer and a narrow one
    rather than between an honest answer and a wrong one. The half-normal prior
    is what lets the amplitude collapse, and this is the assertion that it does:
    the fitted amplitude in the control is more than an order of magnitude below
    the one in the 7 % fringing scenario.
    """
    control = study.summarise(study_results[("none", "flexible")]["run"], names=[AMPLITUDE])
    fringed = study.summarise(
        study_results[("strong_smooth", "flexible")]["run"], names=[AMPLITUDE]
    )
    assert control[AMPLITUDE].median < 0.1 * fringed[AMPLITUDE].median
    assert control[AMPLITUDE].median < 0.01  # Jy, against a ~1 Jy continuum


def test_the_gp_amplitude_tracks_the_injected_amplitude(study_results: Any) -> None:
    """Nothing, then 2.5 %, then 7 %: the fitted amplitude is monotone in them."""
    fitted = [
        study.summarise(study_results[(key, "flexible")]["run"], names=[AMPLITUDE])[
            AMPLITUDE
        ].median
        for key in ("none", "mild", "strong_smooth")
    ]
    assert fitted[0] < fitted[1] < fitted[2]


# ---------------------------------------------------------------------------
# 2 and 3. The headline
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key", MISSPECIFIED)
def test_the_standard_likelihood_is_confidently_wrong(study_results: Any, key: str) -> None:
    summaries = _summaries(study_results, key, "standard")
    worst = max(summary.bias_in_widths for summary in summaries.values())
    assert worst >= study.STANDARD_MIN_BIAS_WIDTHS, {
        name: round(s.bias_in_widths, 2) for name, s in summaries.items()
    }
    missed = [name for name, summary in summaries.items() if not summary.covers_truth]
    assert missed, "the standard likelihood's 68 % interval contained the truth everywhere"


@pytest.mark.parametrize("key", ["none", *MISSPECIFIED])
def test_the_flexible_likelihood_stays_honest(study_results: Any, key: str) -> None:
    summaries = _summaries(study_results, key, "flexible")
    for name, summary in summaries.items():
        assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, (
            f"{key} {name}: {summary!r}"
        )
        assert summary.covers_truth, f"{key} {name}: {summary!r}"


@pytest.mark.parametrize("key", MISSPECIFIED)
def test_the_flexible_likelihood_beats_the_standard_one_where_it_matters(
    study_results: Any, key: str
) -> None:
    """The comparison itself, parameter by parameter.

    On the worst-affected parameter of each misspecified scenario the flexible
    likelihood is closer to the truth *in posterior widths* by a wide margin.
    That is the honest form of the comparison: the standard fit's problem is
    not that its median moved but that its error bar did not.
    """
    standard = _summaries(study_results, key, "standard")
    flexible = _summaries(study_results, key, "flexible")
    worst = max(standard, key=lambda name: standard[name].bias_in_widths)
    assert flexible[worst].bias_in_widths < standard[worst].bias_in_widths


# ---------------------------------------------------------------------------
# 4. The diagnostics
# ---------------------------------------------------------------------------


def test_whiteness_finds_nothing_in_the_control(diagnoses: Any) -> None:
    result = diagnoses[("none", "standard")]
    assert result.whiteness_p_value is not None
    assert result.whiteness_p_value >= study.WHITENESS_WHITE_LEVEL


@pytest.mark.parametrize("key", MISSPECIFIED)
def test_whiteness_finds_the_misspecification(diagnoses: Any, key: str) -> None:
    result = diagnoses[(key, "standard")]
    assert result.whiteness_p_value is not None
    assert result.whiteness_p_value <= study.WHITENESS_STRUCTURE_LEVEL
    control = diagnoses[("none", "standard")]
    assert result.whiteness_statistic > control.whiteness_statistic


def test_gp_localisation_peaks_where_the_line_was_injected(diagnoses: Any) -> None:
    """The one scenario whose deviation *has* a location."""
    scenario = scenario_named("strong_sharp")
    assert scenario.localised_at is not None
    spacing = float(np.diff(wavelength_grid(200))[0])
    result = diagnoses[("strong_sharp", "flexible")]
    assert result.localisation_peak is not None
    offset = abs(result.localisation_peak - scenario.localised_at)
    assert offset <= study.LOCALISATION_TOLERANCE_POINTS * spacing, (
        f"peak at {result.localisation_peak:.6f} um, injected at {scenario.localised_at:.6f} um, "
        f"grid spacing {spacing:.2e} um"
    )


def test_gp_localisation_says_more_where_there_is_more_to_say(diagnoses: Any) -> None:
    """A peak somewhere is not a detection; a peak far above the control is."""
    control = diagnoses[("none", "flexible")]
    sharp = diagnoses[("strong_sharp", "flexible")]
    assert control.localisation_max is not None and sharp.localisation_max is not None
    assert sharp.localisation_max >= study.LOCALISATION_CONTRAST * control.localisation_max


# ---------------------------------------------------------------------------
# Provenance: the runs know what they were runs of.
# ---------------------------------------------------------------------------


def test_every_run_records_its_backend_and_engine(study_results: Any) -> None:
    for (key, kind), entry in study_results.items():
        attrs = entry["run"].attrs
        assert attrs["ampere_engine"] == "emcee", (key, kind)
        assert attrs["ampere_backend"] == "reference", (key, kind)
        assert len(attrs["ampere_spec_hash"]) == 32, (key, kind)


def test_the_two_likelihoods_are_told_apart_by_their_provenance(study_results: Any) -> None:
    """The standard fit has no GP, and the stored run says so rather than implying it."""
    from ampere.results import GP_LOCALISATION_GROUP

    for (_key, kind), entry in study_results.items():
        has_group = GP_LOCALISATION_GROUP in entry["run"].children
        assert has_group is (kind == "flexible")


# ---------------------------------------------------------------------------
# The legacy cross-check: the same conclusion with legacy's kernel.
# ---------------------------------------------------------------------------


def test_the_legacy_kernel_reaches_the_same_conclusion() -> None:
    """W2.10 changed the kernel; this is the check that the change did not change
    the answer.

    The reproduction uses Matern-3/2 because it is v2's default and because it
    is the only kernel with an O(N) path — without it the 20 000-point rung of
    the ladder does not exist. The paper study used a squared exponential. So
    the honest question is not whether the two kernels give the same
    log-likelihood (they do not, and should not), but whether the *conclusion*
    depends on the choice: does a squared-exponential flexible likelihood also
    keep the truth inside its 68 % interval where the standard one does not?

    Run densely at 200 points, because ``QuasisepGP`` refuses a
    non-quasiseparable kernel — which is itself the reason for the default.
    """
    from examples.m2_misspecification.generators import generate

    data = generate("strong_smooth", size=200)
    problem = study.build_problem(
        data, likelihood="flexible", kernel="squared_exponential", solver="dense"
    )
    summaries = study.summarise(study.run(problem, study.TEST_EMCEE), names=study.PHYSICAL_NAMES)
    for name, summary in summaries.items():
        assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, f"{name}: {summary!r}"
        assert summary.covers_truth, f"{name}: {summary!r}"
