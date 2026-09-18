"""The same posterior, whichever backend and whichever sampler produced it.

The comparison this module makes is **reference emcee against this
environment's NUTS**, on the same data, the same priors and the same
likelihood. It is a pairwise comparison rather than a three-way one for a
practical reason worth stating: torch and jax live in separate pixi
environments, so no single interpreter has both, and "torch agrees with jax"
is established transitively through the reference backend that every
environment does have. The reference side is re-sampled here rather than
cached, so both sides of every comparison come from one interpreter, one
numpy and one set of seeds.

What is compared, and on what scale
-----------------------------------
The four physical parameters' posterior medians and 68 % interval endpoints,
each divided by the mean of the two posteriors' own half-widths
(:func:`~examples.m2_misspecification.study.agreement`). Posterior widths are
the only scale on which "do these agree?" has an answer: two correct samplers
differ by their Monte Carlo error, and Monte Carlo error is a fraction of the
posterior, never a number of Jy.

The tolerances, and why they are these
--------------------------------------
:data:`~examples.m2_misspecification.study.CI_MEDIAN_TOLERANCE` = 0.5 posterior
widths on the median and
:data:`~examples.m2_misspecification.study.CI_INTERVAL_TOLERANCE` = 0.75 on
either 68 % endpoint, at
:data:`~examples.m2_misspecification.study.AGREEMENT_EMCEE` /
:data:`~examples.m2_misspecification.study.AGREEMENT_NUTS`. Those budgets are
longer than the science suite's, and deliberately so: at the *shorter* budget
the measured run-to-run scatter of a median reaches 0.31 posterior widths on
the worst-conditioned parameter in the study — ``B`` under the flexible
likelihood, partly degenerate with the GP — which is close enough to any
sensible tolerance that the test would be measuring the samplers' luck. At the
budgets used here the same scatter is 0.14 (emcee) and 0.07 (NUTS)
(:data:`~examples.m2_misspecification.study.MEASURED_SCATTER`), so the
difference of two runs has a standard deviation of roughly 0.09 widths and the
tolerance is about five times it. The interval tolerance is looser than the
median one for the same reason it should be: a quantile of a broad,
asymmetric posterior is the noisiest thing either sampler reports.

The tighter statement — 0.20 widths on the median — is what the *milestone*
budget achieves and what ``tests/m2/test_ladder.py`` asserts under
``-m m2_full``; the milestone runs recorded in the documentation page came in
at 0.072 widths, worst case, across both likelihoods and both backends. All
three numbers are kept because they say different things: what a per-PR gate
can afford, what a marked run asserts, and what the backends actually do.

Below the sampling, one much sharper check that costs milliseconds and would
catch almost any real disagreement first: the three backends' ``log_prob`` at
the same point, to a relative tolerance of 1e-10
(``tests/m2/test_model.py::test_the_three_backends_score_the_same_density``).
A posterior disagreement with an agreeing density is a sampler problem; a
density disagreement is a backend problem, and having both makes the two
distinguishable.
"""

from __future__ import annotations

from typing import Any

import pytest

from examples.m2_misspecification import study
from examples.m2_misspecification.generators import MANY_LINES, generate

#: The scenario the comparison is made on: the one the milestone is about, and
#: the one whose flexible posterior is the hardest to sample.
SCENARIO = "strong_smooth"

#: W5.8's arm: the scenario with two length scales, fitted with W5.7's warped
#: Matérn. It is the one declaration in the study whose parameters a gradient
#: sampler reaches *through the warp* — five increments and a shrinkage scale,
#: all ordinary ``Parameter`` objects — so it is where "NUTS over
#: the knots for free" is either true on both backends or not true at all.
WARPED_SCENARIO = MANY_LINES

pytestmark = pytest.mark.parametrize("backend", ["torch", "jax"])


@pytest.fixture(scope="module")
def comparisons() -> dict[tuple[str, str], dict[str, study.Summary]]:
    """One :func:`~examples.m2_misspecification.study.summarise` per backend and likelihood.

    Module-scoped: the reference side is sampled once per likelihood and reused
    for whichever accelerated backend this environment has.
    """
    return {}


def _summaries(
    cache: dict[tuple[str, str], dict[str, study.Summary]], backend: str, kind: str
) -> dict[str, study.Summary]:
    key = (backend, kind)
    if key not in cache:
        data = generate(SCENARIO, size=200)
        problem = study.build_problem(data, backend=backend, likelihood=kind)
        budget = study.AGREEMENT_EMCEE if backend == "reference" else study.AGREEMENT_NUTS
        cache[key] = study.summarise(study.run(problem, budget), names=study.PHYSICAL_NAMES)
    return cache[key]


@pytest.mark.parametrize("kind", list(study.LIKELIHOODS))
def test_the_posteriors_agree_within_the_stated_tolerance(
    comparisons: Any, backend: str, kind: str
) -> None:
    pytest.importorskip(backend)
    reference = _summaries(comparisons, "reference", kind)
    native = _summaries(comparisons, backend, kind)
    result = study.agreement(reference, native)
    for name, entry in result.items():
        assert entry.median <= study.CI_MEDIAN_TOLERANCE, f"{backend} {kind} {name}: {entry!r}"
        assert entry.worst_interval <= study.CI_INTERVAL_TOLERANCE, (
            f"{backend} {kind} {name}: {entry!r}"
        )


def test_the_native_run_records_its_own_backend_and_sampler(comparisons: Any, backend: str) -> None:
    """A NUTS run's provenance is a fact, not a declaration (W2.12)."""
    pytest.importorskip(backend)
    data = generate(SCENARIO, size=200)
    problem = study.build_problem(data, backend=backend, likelihood="flexible")
    assert problem.backend == backend
    run = study.run(problem, study.NutsBudget(draws=30, warmup=30, chains=1))
    assert run.attrs["ampere_backend"] == backend
    assert run.attrs["ampere_engine"] == "nuts"


def test_nuts_samples_the_warp_knots_on_this_backend(backend: str) -> None:
    """W5.8: the warp's knot variables are ordinary parameters, so NUTS gets them free.

    W5.7's claim about the warp was structural — the knot variables are
    ``Parameter`` objects like any other, so a differentiable
    backend reaches them through :func:`ampere.core.realise` with nothing added
    — and this is where it is exercised end to end rather than described. A
    short budget: what is asserted is that the six warp dimensions are in the
    posterior under the names the declaration gives them, and that the run
    records the backend and sampler that produced it.
    """
    pytest.importorskip(backend)
    data = generate(WARPED_SCENARIO, size=200)
    problem = study.build_problem(data, backend=backend, likelihood="flexible", kernel="warped")
    assert problem.backend == backend
    prefix = f"{study.DATASET_LABEL}.likelihood.input_warp"
    expected = (f"{prefix}.scale", *(f"{prefix}.increment{i}" for i in range(4)))
    run = study.run(problem, study.NutsBudget(draws=40, warmup=40, chains=1))
    assert run.attrs["ampere_backend"] == backend
    assert run.attrs["ampere_engine"] == "nuts"
    available = {str(name) for name in run["posterior"].dataset.data_vars}
    assert set(expected) <= available, sorted(available)


def test_the_warped_arm_reaches_the_same_posterior_on_this_backend(backend: str) -> None:
    """The cross-backend row, extended to W5.8's warped arm.

    Reference emcee against this environment's NUTS, on the same spectrum, the
    same priors and the same warped kernel — the comparison the rest of this
    module makes of the stationary one, now over a declaration with six more
    dimensions and a non-centred parameterisation between them. The tolerances
    are :data:`~examples.m2_misspecification.study.CI_MEDIAN_TOLERANCE` and
    :data:`~examples.m2_misspecification.study.CI_INTERVAL_TOLERANCE`, the same
    per-PR pair the stationary comparison uses.
    """
    pytest.importorskip(backend)
    data = generate(WARPED_SCENARIO, size=200)
    reference = study.summarise(
        study.run(
            study.build_problem(data, likelihood="flexible", kernel="warped"),
            study.MANY_LINES_EMCEE,
        ),
        names=study.PHYSICAL_NAMES,
    )
    native = study.summarise(
        study.run(
            study.build_problem(data, backend=backend, likelihood="flexible", kernel="warped"),
            study.MANY_LINES_NUTS,
        ),
        names=study.PHYSICAL_NAMES,
    )
    result = study.agreement(reference, native)
    for name, entry in result.items():
        print(f"\n{backend} warped {name}: {entry!r}")
    for name, entry in result.items():
        assert entry.median <= study.CI_MEDIAN_TOLERANCE, f"{backend} warped {name}: {entry!r}"
        assert entry.worst_interval <= study.CI_INTERVAL_TOLERANCE, (
            f"{backend} warped {name}: {entry!r}"
        )


def test_the_flexible_likelihood_is_honest_on_every_backend(comparisons: Any, backend: str) -> None:
    """The science claim, re-made where NUTS made the posterior.

    The claim is about the likelihood, so it must not depend on which sampler
    explored it. ``tests/m2/test_science.py`` establishes it on the reference
    backend at length; this is the one-line restatement on the other two.
    """
    pytest.importorskip(backend)
    flexible = _summaries(comparisons, backend, "flexible")
    standard = _summaries(comparisons, backend, "standard")
    for name, summary in flexible.items():
        assert summary.bias_in_widths <= study.FLEXIBLE_MAX_BIAS_WIDTHS, f"{name}: {summary!r}"
        assert summary.covers_truth, f"{name}: {summary!r}"
    worst = max(summary.bias_in_widths for summary in standard.values())
    assert worst >= study.STANDARD_MIN_BIAS_WIDTHS
