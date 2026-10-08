"""The mild scenario's GP localisation: ``plot_gp_localisation`` on the flexible fit.

Reproduces ``examples.m2_misspecification.figures.save_result_figures``'s
``localisation_mild_flexible`` figure at the docs budget (``_m2.DOC_BUDGET``).
Full-budget numbers: ``m2_misspecification.rst``, "The diagnostics agree", row
``mild`` (peak score 3.70), and the GP table (0.0248 Jy at 0.00089 um).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _m2

from ampere.results import plot_gp_localisation

entry = _m2.runs()[("mild", "flexible")]
label = next(iter(entry["problem"].datasets))
plot_gp_localisation(entry["run"], datasets=[label])
