"""The full size ladder, and the tight cross-backend statement. ``-m m2_full``.

Not in the PR gate. Everything here samples 2 000- or 20 000-point spectra, or
samples 200-point ones at the milestone budget, and costs tens of minutes on
three backends. What it buys is the two claims a short chain cannot make:

**The ladder.** ``DEVELOPMENT_PLAN.md`` §5's M2 asks for the study at ten to a
hundred times the paper's data size, and the interesting thing that happens
along the ladder is not that the fits get slower. It is that the *cost of being
wrong* grows: the posterior narrows as :math:`\\sqrt{N}` while the unmodelled
7 % ripple does not shrink at all, so the standard likelihood's error, measured
in its own standard deviations, grows as :math:`\\sqrt{N}`. Measured on
``strong_smooth`` at the milestone budget, the worst of the four parameters is
**11.2, 36.3 and 113.0** posterior widths from the truth at 200, 2 000 and
20 000 points — against 11.2, 35.4 and 112.0 predicted by a fixed bias and a
shrinking posterior — and the truth is outside the 68 % interval at every rung.
The flexible likelihood's worst offset does not grow at all: **0.70, 0.51 and
0.57** widths, truth covered throughout. That is the claim asserted here, and
it is the one that matters for real spectra, which have thousands of pixels
rather than two hundred.

**The tight agreement.** At the milestone budget the cross-backend disagreement
is a measurement of the two posteriors rather than of the two samplers'
remaining Monte Carlo error, so the tolerance drops from the gate's 0.5
posterior widths to
:data:`~examples.m2_misspecification.study.MILESTONE_MEDIAN_TOLERANCE`.

Run it with::

    pixi run -e dev  pytest tests/m2 -m m2_full
    pixi run -e jax  pytest tests/m2 -m m2_full
    pixi run -e torch pytest tests/m2 -m m2_full
"""

from __future__ import annotations

import pytest

from examples.m2_misspecification import study
from examples.m2_misspecification.generators import SIZES, generate

pytestmark = pytest.mark.m2_full


def _report(label: str, size: int, summaries: dict[str, study.Summary]) -> None:
    """Print the row, so ``pytest -m m2_full -s`` produces the table as well as the verdict.

    ``tests/scaling`` established the convention: a suite whose whole purpose is
    a measurement should hand the measurement to the reader rather than only
    assert something about it, because the number is what goes into
    ``docs/source/m2_misspecification.rst`` and into the work item's report.
    """
    offsets = " ".join(
        f"{name.split('.')[-1]}={summary.bias_in_widths:7.2f}"
        for name, summary in summaries.items()
    )
    covers = all(summary.covers_truth for summary in summaries.values())
    print(f"  [ladder] {label:9s} n={size:6d}  {offsets}  covers={covers}")


@pytest.mark.parametrize("size", SIZES)
def test_the_flexible_likelihood_stays_honest_at_every_rung(size: int) -> None:
    """Every parameter within the stated tolerance of the truth, at 200 to 20 000."""
    data = generate("strong_smooth", size=size)
    problem = study.build_problem(data, likelihood="flexible", solver="quasisep")
    summaries = study.summarise(
        study.run(problem, study.MILESTONE_EMCEE), names=study.PHYSICAL_NAMES
    )
    _report("flexible", size, summaries)
    for name, summary in summaries.items():
        assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, (
            f"n={size} {name}: {summary!r}"
        )
        assert summary.covers_truth, f"n={size} {name}: {summary!r}"


def test_the_standard_likelihoods_error_grows_with_the_data() -> None:
    """More data makes a misspecified fit worse, not better.

    The number that grows is ``|median - truth|`` in posterior widths: the
    ripple does not shrink, and the posterior does. Asserted as a monotone
    increase in the worst parameter's offset across the whole ladder, because
    that is the shape of the statement rather than a particular value of it.
    """
    worst = []
    for size in SIZES:
        data = generate("strong_smooth", size=size)
        problem = study.build_problem(data, likelihood="standard")
        summaries = study.summarise(
            study.run(problem, study.MILESTONE_EMCEE), names=study.PHYSICAL_NAMES
        )
        _report("standard", size, summaries)
        worst.append(max(summary.bias_in_widths for summary in summaries.values()))
    assert worst[0] < worst[1] < worst[2], dict(zip(SIZES, worst, strict=True))
    assert worst[-1] > 3.0 * worst[0]


@pytest.mark.parametrize("backend", ["torch", "jax"])
@pytest.mark.parametrize("kind", list(study.LIKELIHOODS))
def test_the_backends_agree_tightly_at_the_milestone_budget(backend: str, kind: str) -> None:
    """The number W2.10's report records, asserted rather than merely reported."""
    pytest.importorskip(backend)
    data = generate("strong_smooth", size=200)
    reference = study.summarise(
        study.run(
            study.build_problem(data, backend="reference", likelihood=kind),
            study.MILESTONE_EMCEE,
        ),
        names=study.PHYSICAL_NAMES,
    )
    native = study.summarise(
        study.run(
            study.build_problem(data, backend=backend, likelihood=kind), study.MILESTONE_NUTS
        ),
        names=study.PHYSICAL_NAMES,
    )
    for name, entry in study.agreement(reference, native).items():
        assert entry.median <= study.MILESTONE_MEDIAN_TOLERANCE, (
            f"{backend} {kind} {name}: {entry!r}"
        )
        assert entry.worst_interval <= study.MILESTONE_INTERVAL_TOLERANCE, (
            f"{backend} {kind} {name}: {entry!r}"
        )


@pytest.mark.parametrize("backend", ["torch", "jax"])
def test_nuts_samples_the_top_of_the_ladder(backend: str) -> None:
    """20 000 points, on the O(N) solver, with gradients: the rung the redesign bought.

    A short chain, because the point is that it runs at all — the dense solve
    this problem would otherwise need is a 20 000 x 20 000 Cholesky per leapfrog
    step. The recovery assertion is the loose one, since a hundred draws is not
    a posterior.
    """
    pytest.importorskip(backend)
    data = generate("strong_smooth", size=20_000)
    problem = study.build_problem(data, backend=backend, likelihood="flexible")
    run = study.run(problem, study.NutsBudget(draws=100, warmup=100, chains=1))
    summaries = study.summarise(run, names=study.PHYSICAL_NAMES)
    for name, summary in summaries.items():
        assert summary.covers_truth, f"{backend} n=20000 {name}: {summary!r}"
