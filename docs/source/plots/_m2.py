"""Shared M2 fits for the figures of ``reading_the_diagnostics.rst``.

Not a figure script (the leading underscore keeps it out of the page): the
scripts import it, so the study is run once per build, not once per figure.

It calls the study's own public functions, as ``python -m
examples.m2_misspecification`` does (``study.run_study`` then
``study.prepare``; ``figures.save_result_figures`` is the per-figure loop these
scripts unroll), on the three scenarios the page walks through: ``none`` (the
control), ``mild`` (2.5 % fringing) and ``strong_sharp`` (a 12 % emission
line). The full-budget numbers are in the tables of ``m2_misspecification.rst``
(32 walkers x 4 000 steps); the docs budget is ``DOC_BUDGET``, 16 walkers x
400 steps with half discarded, seeded by the study's own ``SEED``.
"""

from __future__ import annotations

import functools

from examples.m2_misspecification import study

from ampere.results import add_posterior_predictive

DOC_BUDGET = study.EmceeBudget(walkers=16, steps=400, burn_in=200)
THIN = 20  # the same thinning study.prepare and figures.save_result_figures use
SCENARIOS = ("none", "mild", "strong_sharp")


@functools.cache
def runs() -> dict:
    """``run_study`` + ``prepare`` over the three scenarios; entries gain their derived groups."""
    results = study.run_study(scenarios=list(SCENARIOS), budget=DOC_BUDGET)
    study.prepare(results, thin=THIN)
    return results


@functools.cache
def with_posterior_predictive(scenario: str, kind: str):
    """The run with the posterior-predictive group added, and its dataset label."""
    entry = runs()[(scenario, kind)]
    label = next(iter(entry["problem"].datasets))
    return add_posterior_predictive(entry["run"], entry["problem"], thin=THIN), label
