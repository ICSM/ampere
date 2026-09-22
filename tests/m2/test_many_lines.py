"""W5.8: a deviation with two length scales, and the three kernels on it.

M2's four scenarios each inject one scale of deviation, so one stationary
length scale copes with all four and the milestone's headline — *the flexible
likelihood stays honest* — holds with a stationary Matérn-3/2 throughout. The
plan's Phase-5 "non-stationary flexible likelihood" bullet asks what happens
when the deviation has **two** scales and only one of them is in one part of
the band, and names the two answers: W5.7's warped Matérn and W4.5's ``Sum``.

``examples/m2_misspecification/many_lines.py`` is that experiment and this
module is its assertions. What is claimed, in the order the classes below make
it:

1. **The scenario is what it says it is.** The forest is confined to the band,
   the two scales are an order of magnitude apart, and both deviations are
   large compared with the noise. Without this the rest is not evidence.
2. **One length scale is not enough.** The standard likelihood is confidently
   wrong, as always. What is new is that the *stationary flexible* likelihood
   fails M2's own
   :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS`
   threshold here — the first time in this study it does — and the warped and
   summed kernels both recover it.
3. **The GP found the deviation where the deviation is.** Every flexible arm
   localises inside the line band, and its score there is far above its score
   outside: the non-stationary arms are not winning by having been handed
   somewhere else to put the residual.
4. **A sum of noise components needs a sparsity guard, and the horseshoe is
   it.** Two nearly degenerate terms on a one-component truth: under a flat
   prior the fit splits itself between them, under
   :func:`~ampere.core.shrinkage_horseshoe` it does not — and the component
   the truth *does* have survives, because a prior that shrank everything would
   pass the same assertion and be useless.

Budget and margins
------------------
The comparison runs at
:data:`~examples.m2_misspecification.study.MANY_LINES_EMCEE` (32 walkers,
1 100 steps, 550 discarded) at 200 points — about 140 seconds for the four
arms, paid once by a module-scoped fixture — and the shrinkage rows at
:data:`~examples.m2_misspecification.many_lines.SHRINKAGE_EMCEE`, which is
longer because the horseshoe's three levels give the posterior a funnel and an
ensemble sampler needs the steps to get round one. Every threshold is a
**margin** rather than a number, in the pattern
``test_fringing_kernel.py``'s ``_period_margin`` set: what is pinned is the
distance between two arms, or between an arm and a threshold the study already
had, never the value one arm happened to reach. The measured numbers at this
budget and seed are recorded beside each assertion.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from ampere.core import QuasisepGP, Sum, WarpedKernel, quantile_knots
from ampere.results import gp_localisation_score
from examples.m2_misspecification import many_lines, study
from examples.m2_misspecification.generators import (
    MANY_LINES,
    NOISE_FRACTION,
    generate,
    scenario_named,
    wavelength_grid,
)

CI_SIZE = 200

#: The three flexible arms, i.e. :data:`many_lines.ARMS` without the standard
#: likelihood.
FLEXIBLE_ARMS = many_lines.ARMS[1:]

#: The two non-stationary arms: the ones W5.7 and W4.5 added.
NON_STATIONARY = ("warped", "sum")


@pytest.fixture(scope="module")
def comparison() -> many_lines.Comparison:
    """The four fits of one ``many_lines`` spectrum, computed once for this module."""
    return many_lines.compare(size=CI_SIZE)


@pytest.fixture(scope="module")
def shrinkage() -> many_lines.Shrinkage:
    """The sparsity demonstration: two degenerate terms, two priors, one truth.

    About three and a half minutes at
    :data:`~examples.m2_misspecification.many_lines.SHRINKAGE_EMCEE`, paid once
    for the four rows that read it.
    """
    return many_lines.shrinkage(size=CI_SIZE)


def _score_in_and_out(entry: Any, band: tuple[float, float]) -> tuple[float, float]:
    """The largest GP-localisation score inside *band*, and the largest outside it."""
    label = next(iter(entry["problem"].datasets))
    score = gp_localisation_score(entry["run"], dataset=label)
    values = np.asarray(score.values, dtype=float)
    coordinates = np.asarray(score.coordinates, dtype=float).reshape(values.size, -1)[:, 0]
    inside = (coordinates >= band[0]) & (coordinates <= band[1])
    return float(np.max(values[inside])), float(np.max(values[~inside]))


# ---------------------------------------------------------------------------
# 0. The declaration: cheap, and it holds before any chain is run.
# ---------------------------------------------------------------------------


class TestTheThreeArmsCompose:
    def test_the_warped_arm_is_a_warped_matern_over_the_recorded_knots(self) -> None:
        kernel = study.build_kernel("reference", "warped")
        assert isinstance(kernel, WarpedKernel)
        assert kernel.base.FAMILY == "matern32"
        assert kernel.input_knots == study.WARP_KNOTS
        assert kernel.amplitude_knots == ()

    def test_the_recorded_knots_are_what_quantile_knots_returns(self) -> None:
        """W5.7's rule: called once, on the study's own coordinate, and recorded.

        Recorded rather than computed at fit time, because a kernel whose knots
        came from the data it is about to be fitted to would have a different
        spec hash for every dataset. This row is what stops the recorded
        numbers drifting away from the rule that produced them — and it holds
        at every rung of the ladder, because ``wavelength_grid`` spans the same
        range at every size.
        """
        for size in (200, 2_000):
            assert quantile_knots(wavelength_grid(size), len(study.WARP_KNOTS)) == pytest.approx(
                study.WARP_KNOTS
            )

    def test_the_sum_arm_is_two_materns_with_two_length_scales(self) -> None:
        kernel = study.build_kernel("reference", "sum")
        assert isinstance(kernel, Sum)
        assert [child.FAMILY for _, child in kernel.terms] == ["matern32", "matern32"]
        assert [label for label, _ in kernel.terms] == ["broad", "narrow"]

    @pytest.mark.parametrize("arm", NON_STATIONARY)
    def test_both_reach_the_o_n_solver(self, arm: str) -> None:
        """Composition is the test: ``QuasisepGP`` refuses what it cannot lower.

        This is what makes the extension affordable at the ladder's larger
        sizes, and it is the structural claim W5.7 and W4.5 each made: a warp
        of a quasiseparable kernel is quasiseparable, and so is a sum of them.
        """
        data = generate(MANY_LINES, size=CI_SIZE)
        problem = many_lines.build_arm(data, arm)
        assert isinstance(problem.datasets["default"].likelihood.noise.solver, QuasisepGP)

    def test_the_warped_arm_costs_six_extra_dimensions(self) -> None:
        """Five increments and one shrinkage scale: the degrees-of-freedom guard, counted."""
        data = generate(MANY_LINES, size=CI_SIZE)
        stationary = many_lines.build_arm(data, "matern32")
        warped = many_lines.build_arm(data, "warped")
        assert warped.free_size - stationary.free_size == len(study.WARP_KNOTS)


# ---------------------------------------------------------------------------
# 1. The scenario is what it says it is.
# ---------------------------------------------------------------------------


class TestTheScenarioHasTwoScales:
    def test_the_forest_is_confined_to_the_band(self) -> None:
        """Outside the band the only deviation is the smooth one."""
        grid = wavelength_grid(2_000)
        full = MANY_LINES.deviation(grid)
        smooth = many_lines.SMOOTH_ONLY.deviation(grid)
        lower, upper = MANY_LINES.band
        outside = (grid < lower - 6.0 * MANY_LINES.line_width) | (
            grid > upper + 6.0 * MANY_LINES.line_width
        )
        assert np.allclose(full[outside], smooth[outside], atol=1e-6)
        assert np.max(full[~outside] - smooth[~outside]) > 0.5 * MANY_LINES.amplitude

    def test_the_two_scales_are_an_order_of_magnitude_apart(self) -> None:
        """The premise of the whole comparison, as a number.

        The forest's correlation length is its line width; the smooth error's
        is a fair fraction of its period. One stationary ``length_scale``
        cannot be both, and this is the statement of by how much.
        """
        assert MANY_LINES.continuum_period / MANY_LINES.line_width > 100.0

    def test_both_deviations_dominate_the_noise(self) -> None:
        """If either were small there would be nothing for a kernel to get wrong."""
        assert MANY_LINES.amplitude > 5.0 * NOISE_FRACTION
        assert MANY_LINES.continuum_amplitude > 2.0 * NOISE_FRACTION

    def test_the_scenario_is_reachable_by_name_but_not_in_the_four(self) -> None:
        """``run_study`` still defaults to the four; this one is asked for."""
        assert scenario_named("many_lines") is MANY_LINES
        assert MANY_LINES not in study.SCENARIOS


# ---------------------------------------------------------------------------
# 2. One length scale is not enough.
# ---------------------------------------------------------------------------


class TestOneLengthScaleIsNotEnough:
    def test_the_standard_likelihood_is_confidently_wrong(
        self, comparison: many_lines.Comparison
    ) -> None:
        """As in every misspecified M2 scenario. Measured: 35.7 posterior widths."""
        worst = comparison.worst_bias["standard"]
        print(f"\nstandard worst bias: {worst:.2f} widths")
        assert worst >= study.STANDARD_MIN_BIAS_WIDTHS
        assert comparison.coverage["standard"] < len(study.PHYSICAL_NAMES)

    def test_the_stationary_flexible_likelihood_misses_too(
        self, comparison: many_lines.Comparison
    ) -> None:
        """The new result, and the reason this scenario exists.

        A stationary Matérn-3/2 has one length scale to spend. Here it is asked
        to follow a line forest *and* a smooth continuum error an order of
        magnitude broader, and whichever it chooses the other is left in the
        residual with nothing to absorb it — so the truth ends up further from
        its median than M2's own threshold allows. Measured: 2.84 widths
        against a threshold of
        :data:`~examples.m2_misspecification.study.MANY_LINES_STATIONARY_MIN_BIAS`
        (2.0), and against
        :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS`
        (1.5), which it is the first fit in this study to fail.
        """
        worst = comparison.worst_bias["matern32"]
        print(f"\nstationary worst bias: {worst:.2f} widths")
        assert worst >= study.MANY_LINES_STATIONARY_MIN_BIAS
        assert worst > study.FLEXIBLE_MAX_BIAS_WIDTHS

    @pytest.mark.parametrize("arm", NON_STATIONARY)
    def test_the_non_stationary_arms_recover_it(
        self, comparison: many_lines.Comparison, arm: str
    ) -> None:
        """M2's claim, re-made for each of the two answers.

        Every parameter inside
        :data:`~examples.m2_misspecification.study.FLEXIBLE_MAX_BIAS_WIDTHS`
        of the truth and its 68 % interval containing it — the same pair of
        assertions ``test_science.py`` makes of the flexible likelihood in the
        four original scenarios. Measured worst bias: warped 0.92, sum 0.78.
        """
        summaries = comparison.summaries[arm]
        for name, summary in summaries.items():
            assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, (
                f"{arm} {name}: {summary!r}"
            )
            assert summary.covers_truth, f"{arm} {name}: {summary!r}"

    @pytest.mark.parametrize("arm", NON_STATIONARY)
    def test_each_beats_the_stationary_arm_by_a_margin(
        self, comparison: many_lines.Comparison, arm: str
    ) -> None:
        """The comparison itself, pinned as a **ratio** rather than as two numbers.

        The ratio is what does not move when the budget does: both arms are
        fitted to the same spectrum at the same budget, so a change in the
        sampler moves both worst biases together and leaves their ratio alone.
        Measured: 3.1x for the warped arm, 3.6x for the ``Sum``.
        """
        stationary = comparison.worst_bias["matern32"]
        improved = comparison.worst_bias[arm]
        print(f"\n{arm}: {stationary:.2f} -> {improved:.2f} widths ({stationary / improved:.1f}x)")
        assert improved * study.MANY_LINES_BIAS_MARGIN <= stationary
        assert comparison.coverage[arm] >= comparison.coverage["matern32"]


# ---------------------------------------------------------------------------
# 3. The GP found the deviation where the deviation is.
# ---------------------------------------------------------------------------


class TestTheDeviationIsLocalisedWhereItIs:
    @pytest.mark.parametrize("arm", FLEXIBLE_ARMS)
    def test_every_flexible_arm_peaks_inside_the_band(
        self, comparison: many_lines.Comparison, arm: str
    ) -> None:
        """A band rather than a point, because a forest has several locations.

        ``strong_sharp``'s single line has a coordinate and
        ``test_science.py`` holds the peak to it within
        :data:`~examples.m2_misspecification.study.LOCALISATION_TOLERANCE_POINTS`
        grid spacings. A forest of five does not; what it has is a band, and
        "the peak is in the band" is the same claim with the same kind of
        margin. Measured: all three peak at 0.86295 µm, inside 0.860-0.870.
        """
        peak = comparison.diagnoses[arm].localisation_peak
        assert peak is not None
        print(f"\n{arm} localisation peak: {peak:.5f} um in band {comparison.band}")
        assert comparison.band[0] <= peak <= comparison.band[1]

    @pytest.mark.parametrize("arm", FLEXIBLE_ARMS)
    def test_the_score_in_the_band_exceeds_the_score_outside_it(
        self, comparison: many_lines.Comparison, arm: str
    ) -> None:
        """A peak in the band is not a detection unless there is less outside it.

        The contrast is within one fit rather than against the control
        scenario, which is what ``test_science.py``'s
        :data:`~examples.m2_misspecification.study.LOCALISATION_CONTRAST` row
        uses: this scenario has a deviation *everywhere* (the smooth error), so
        the question worth asking is not "did the GP find anything?" but "did
        it find more in the band than outside it?". Measured: 2.76 for the
        stationary arm, 1.90 for the warped one and 3.02 for the ``Sum``.
        """
        inside, outside = _score_in_and_out(comparison.entries[arm], comparison.band)
        print(f"\n{arm} localisation score: {inside:.2f} in band, {outside:.2f} outside")
        assert inside >= study.MANY_LINES_LOCALISATION_CONTRAST * outside


class TestTheDiagnosticsAgree:
    def test_whiteness_finds_the_misspecification_in_the_standard_fit(
        self, comparison: many_lines.Comparison
    ) -> None:
        """Measured p = 0.005, the 199-permutation floor."""
        result = comparison.diagnoses["standard"]
        assert result.whiteness_p_value is not None
        assert result.whiteness_p_value <= study.WHITENESS_STRUCTURE_LEVEL


# ---------------------------------------------------------------------------
# 4. The sparsity guard on summed noise components.
# ---------------------------------------------------------------------------


class TestTheHorseshoeShrinksASpuriousComponent:
    """The plan's "same guard generalises to every sum of noise components".

    The declaration is :func:`~ampere.core.shrinkage_horseshoe`, put on the
    kernel by :func:`~ampere.core.with_shrinkage`. What is fitted is two nearly
    degenerate Matérn terms on a truth with **one** smooth component, so the
    likelihood constrains the total and says very little about the split — and
    a prior that says nothing either will leave the fit spread across both.
    """

    def test_the_declaration_is_the_three_levels_it_claims_to_be(self) -> None:
        """Global scale, one local scale per component, one amplitude under each."""
        kernel = many_lines.horseshoe_kernel()
        assert kernel.parameters.names[:3] == (
            "shrinkage.global_scale",
            "shrinkage.first",
            "shrinkage.second",
        )
        assert kernel.parameters["first.amplitude"].references == ("shrinkage.first",)
        assert kernel.parameters["second.amplitude"].references == ("shrinkage.second",)
        for label in many_lines.DEGENERATE_LABELS:
            assert kernel.parameters[f"shrinkage.{label}"].references == ("shrinkage.global_scale",)

    def test_the_contrast_differs_in_the_amplitude_priors_and_nothing_else(self) -> None:
        """What makes it a contrast: same terms, same length scales, same solver."""
        shrunk = many_lines.horseshoe_kernel()
        plain = many_lines.degenerate_pair()
        assert [child.FAMILY for _, child in shrunk.terms] == [
            child.FAMILY for _, child in plain.terms
        ]
        for label in many_lines.DEGENERATE_LABELS:
            assert (
                shrunk.parameters[f"{label}.length_scale"]
                == plain.parameters[f"{label}.length_scale"]
            )
            assert shrunk.parameters[f"{label}.amplitude"] != plain.parameters[f"{label}.amplitude"]

    def test_putting_the_declaration_on_does_not_change_the_kernel_it_is_put_on(self) -> None:
        """:func:`~ampere.core.with_shrinkage` is functional, and this is the check."""
        plain = many_lines.degenerate_pair()
        before = plain.parameters.names
        _ = many_lines.with_shrinkage(
            plain,
            many_lines.shrinkage_horseshoe(
                tuple(f"{label}.amplitude" for label in many_lines.DEGENERATE_LABELS)
            ),
        )
        assert plain.parameters.names == before

    def test_the_horseshoe_concentrates_the_fit_on_one_component(
        self, shrinkage: many_lines.Shrinkage
    ) -> None:
        """The pinned claim, as a ratio of ratios.

        ``ratio`` is the smaller component's amplitude over the larger one's:
        1 means the fit is split evenly between two terms the truth does not
        have two of, and 0 means it chose one. What is asserted is the
        **factor** between the two priors, not either ratio on its own, for the
        reason W4.5's fringing margin gives: the factor is the effect of the
        prior, and it is what survives a change of budget. Measured here 2.93,
        and over three run seeds at a narrower prior ceiling 2.78, 2.96 and
        1.90, against a threshold of
        :data:`~examples.m2_misspecification.study.MANY_LINES_SHRINKAGE_FACTOR`.
        """
        ratios = shrinkage.ratio
        factor = ratios["flat"] / ratios["horseshoe"]
        print(
            f"\nmedian min/max: flat {ratios['flat']:.3f} -> horseshoe "
            f"{ratios['horseshoe']:.3f} ({factor:.2f}x)"
        )
        assert factor >= study.MANY_LINES_SHRINKAGE_FACTOR

    def test_the_same_claim_as_a_posterior_mass(self, shrinkage: many_lines.Shrinkage) -> None:
        """The second statistic, and the reason there are two.

        A posterior mass is far less sensitive than a quantile to how well an
        ensemble sampler explored a hierarchical geometry — and the horseshoe's
        three levels give the posterior a funnel, which is the geometry an
        ensemble sampler explores worst. Two statistics moving the same way is
        what says the effect is the prior's. Measured here 2.09, and over three
        run seeds 2.15, 2.09 and 1.68 — a relative scatter of 13 % against the
        median ratio's 22 %, which is why this is the primary statistic.
        """
        mass = shrinkage.sparse_mass
        factor = mass["horseshoe"] / mass["flat"]
        print(
            f"\nP(min/max < {many_lines.SPARSE_FRACTION}): flat {mass['flat']:.3f} -> "
            f"horseshoe {mass['horseshoe']:.3f} ({factor:.2f}x)"
        )
        assert factor >= study.MANY_LINES_SHRINKAGE_MASS_FACTOR

    def test_the_component_the_truth_has_survives(self, shrinkage: many_lines.Shrinkage) -> None:
        """A prior that shrank *everything* would pass the rows above and be useless.

        The truth has one smooth component, and the larger amplitude must still
        be there under both priors: the guard is sparsity, not silence.
        Measured: 0.0118 Jy under the flat prior and 0.0077 under the
        horseshoe, against a per-point noise sigma of about 0.0099 Jy.
        """
        for prior, largest in shrinkage.largest.items():
            print(f"\n{prior} larger amplitude: {largest:.5f} Jy")
            assert largest > study.MANY_LINES_SIGNAL_FLOOR


# ---------------------------------------------------------------------------
# The driver: the command the documentation section tells a reader to run.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("kernel", list(study.MANY_LINES_KERNELS))
def test_the_driver_runs_this_scenario_with_each_kernel(
    kernel: str, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    """``python -m examples.m2_misspecification --scenario many_lines --kernel ...``.

    A **tiny** run — the numbers are the other rows' business, and the budget is
    monkeypatched down to a few dozen steps so this costs seconds — but a real
    one through the driver's own code path, because that path has a step the
    library API does not: it prints each flexible fit's GP hyperparameters, and
    it has to find their names. The stationary kernel declares an ``amplitude``
    and a ``length_scale``; the warped arm declares six warp variables under
    ``input_warp``; the ``Sum`` declares two labelled terms. A driver that
    spells the stationary pair raises ``KeyError`` on the other two — it did,
    before this row — so what is asserted is that each arm's own hyperparameter
    names reach the table.
    """
    from examples.m2_misspecification.__main__ import main

    expected = {
        "matern32": ("amplitude", "length_scale"),
        "warped": ("base.amplitude", "input_warp.scale", "input_warp.increment0"),
        "sum": ("broad.amplitude", "narrow.length_scale"),
    }[kernel]
    tiny = study.EmceeBudget(walkers=32, steps=80, burn_in=40)
    monkeypatch.setattr(study, "TEST_EMCEE", tiny)
    monkeypatch.setattr(study, "MANY_LINES_EMCEE", tiny)
    assert main(["--scenario", "many_lines", "--kernel", kernel, "--size", "120", "--quick"]) == 0
    printed = capsys.readouterr().out
    assert "GP hyperparameters" in printed
    for name in expected:
        assert name in printed, printed


# ---------------------------------------------------------------------------
# The figure the documentation section points at.
# ---------------------------------------------------------------------------


def test_the_figure_draws_every_flexible_arm_over_the_injected_deviation(
    comparison: many_lines.Comparison, agg_backend: None, tmp_path: Any
) -> None:
    """The figure is the argument in one picture, so it must not rot silently.

    It is written at run time and never committed (ground rule 7), so what a
    test can hold it to is that it builds from the *stored* runs — one line per
    flexible arm, drawn from the ``gp_localisation`` group
    :func:`~examples.m2_misspecification.study.diagnose` derived, with nothing
    re-fitted — and that it lands on disk when asked.
    """
    from examples.m2_misspecification import figures

    figure = figures.figure_many_lines(comparison.entries, comparison.data, band=comparison.band)
    lower, upper = figure.axes
    drawn = {line.get_label() for line in upper.get_lines()}
    assert set(FLEXIBLE_ARMS) <= drawn, drawn
    assert "injected deviation" in {line.get_label() for line in lower.get_lines()}
    written = figures.save_many_lines_figure(
        comparison.entries, comparison.data, tmp_path, suffix="png"
    )
    assert [path.name for path in written] == [f"{figures.MANY_LINES_FIGURE}.png"]
    assert written[0].stat().st_size > 0


# ---------------------------------------------------------------------------
# The full budget: the same claims with ten times the data.
# ---------------------------------------------------------------------------


@pytest.mark.m2_full
def test_the_shrinkage_factors_at_the_marked_budget() -> None:
    """The same two factors where the funnel is explored properly.

    The per-PR rows pin a third below the worst of three seeds, because an
    ensemble sampler at a per-PR budget leaves that much scatter on a
    hierarchical posterior. This row pays for the budget that removes it and
    pins what the prior actually does.
    """
    result = many_lines.shrinkage(
        size=CI_SIZE, budget=study.EmceeBudget(walkers=32, steps=5_000, burn_in=2_500)
    )
    print(f"\n{result.table()}")
    ratio = result.ratio
    mass = result.sparse_mass
    assert ratio["flat"] / ratio["horseshoe"] >= study.MANY_LINES_FULL_SHRINKAGE_FACTOR
    assert mass["horseshoe"] / mass["flat"] >= study.MANY_LINES_FULL_SHRINKAGE_MASS_FACTOR


@pytest.mark.m2_full
def test_the_conclusion_holds_at_the_next_rung_of_the_ladder() -> None:
    """2 000 points, same budget: the stationary arm's failure gets worse, not better.

    The size ladder samples the *same* spectrum more finely, so the posterior
    narrows as sqrt(N) while the unmodelled deviation stays where it was —
    which is exactly the regime a misspecified fit is most damaged by. The
    claim is therefore the same claim, and the margin should widen.
    """
    comparison = many_lines.compare(size=2_000)
    print(f"\n{comparison.table()}")
    assert comparison.worst_bias["matern32"] > study.FLEXIBLE_MAX_BIAS_WIDTHS
    for arm in NON_STATIONARY:
        assert (
            comparison.worst_bias[arm] * study.MANY_LINES_BIAS_MARGIN
            <= comparison.worst_bias["matern32"]
        )
