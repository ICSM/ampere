"""The control's residual plot: ``plot_residuals`` on the ``none`` scenario's standard fit.

Reproduces ``examples.m2_misspecification.figures.save_result_figures``'s
``residuals_none_standard`` figure at the docs budget (``_m2.DOC_BUDGET``).
Full-budget numbers: ``m2_misspecification.rst``, "The diagnostics agree",
row ``none`` (whiteness p = 0.185).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _m2

from ampere.results import plot_residuals

entry = _m2.runs()[("none", "standard")]
label = next(iter(entry["problem"].datasets))
plot_residuals(entry["run"], datasets=[label])
