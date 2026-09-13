"""Arm (d), the chromatic case (``phase4_placement_memo.md`` §3.6): informational.

The item's own word for this scenario's numbers is "informational" — the
claim under test is that a residual sharp in wavelength and smooth in
``(u, v)`` needs a kernel that sees both, and the three fits below *measure*
that rather than assert an ordering, because at this suite's per-PR budget
(:data:`~examples.interferometry.study.CHROMATIC_BUDGET`, three short emcee
runs) a two-parameter recovery from 168 points is noisy enough that any
single ranking is a report, not a claim. What is asserted is only that every
one of the three kernels composes and produces a usable posterior — the
structural half that has to hold before a ranking, at any budget, means
anything.

``pytest -m interferometry_full`` runs the same three fits at the
documentation budget and prints the ranking table without hiding it, which
is where a reader goes for the actual numbers (also recorded in this item's
branch report and in ``docs/source/interferometry.rst``).
"""

from __future__ import annotations

import math

import pytest

from examples.interferometry import study


class TestTheChromaticFitsCompose:
    def test_the_three_kernels_act_on_disjoint_or_matching_axes(self) -> None:
        spatial = study._chromatic_kernel("reference", "spatial")
        spectral = study._chromatic_kernel("reference", "spectral")
        product = study._chromatic_kernel("reference", "product")
        assert spatial.selected_axes(("u", "v", "spectral_axis")) == ("u", "v")
        assert spectral.selected_axes(("u", "v", "spectral_axis")) == ("spectral_axis",)
        assert product.QUASISEPARABLE is False

    def test_the_dispersed_coverage_spans_several_distinct_wavelengths(self) -> None:
        from examples.interferometry import generators as gen

        _, _, waves = gen.dispersed_visibility_coverage()
        assert len(set(waves.round(6).tolist())) == gen.N_CHANNELS


@pytest.mark.usefixtures("agg_backend")
class TestTheChromaticArmRuns:
    """Smoke-level: every kernel fits, and every number that comes out is finite."""

    @pytest.fixture(scope="class")
    @classmethod
    def results(cls) -> dict:
        return study.run_chromatic_arm()

    def test_every_kernel_produced_a_finite_posterior(self, results: dict) -> None:
        for kind, entry in results.items():
            summaries = study.summarise(entry["run"], names=list(study.gen.TRUTH))
            for name, summary in summaries.items():
                assert math.isfinite(summary.median), f"{kind}/{name} median is not finite"
                assert summary.width > 0.0, f"{kind}/{name} has a degenerate posterior"

    def test_the_ranking_is_reported(self, results: dict) -> None:
        """Prints the table ``-s`` shows; nothing here is a pinned ordering."""
        print("\narm (d): bias in posterior widths, one seed, the per-PR budget")
        for kind, entry in results.items():
            summaries = study.summarise(entry["run"], names=list(study.gen.TRUTH))
            biases = {name: round(s.bias_in_widths, 3) for name, s in summaries.items()}
            print(f"  {kind:<10s} {biases}")


@pytest.mark.interferometry_full
@pytest.mark.usefixtures("agg_backend")
class TestTheChromaticArmAtFullBudget:
    def test_the_ranking_at_the_documentation_budget(self) -> None:
        results = study.run_chromatic_arm(budget=study.DOC_BUDGET)
        print("\narm (d) at the documentation budget:")
        for kind, entry in results.items():
            summaries = study.summarise(entry["run"], names=list(study.gen.TRUTH))
            biases = {name: round(s.bias_in_widths, 3) for name, s in summaries.items()}
            print(f"  {kind:<10s} {biases}")
