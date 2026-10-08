"""The posterior-predictive check on ``strong_sharp``'s standard fit.

Reproduces ``examples.m2_misspecification.figures.save_result_figures``'s
``posterior_predictive_strong_sharp_standard`` figure at the docs budget
(``_m2.DOC_BUDGET``). Full-budget numbers: ``m2_misspecification.rst``, "The
diagnostics agree", row ``strong_sharp`` (whiteness Q = 272 for the same fit).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _m2

from ampere.results import plot_posterior_predictive

run, label = _m2.with_posterior_predictive("strong_sharp", "standard")
plot_posterior_predictive(run, datasets=[label])
