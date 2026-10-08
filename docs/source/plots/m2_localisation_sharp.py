"""The strong_sharp GP localisation: ``plot_gp_localisation`` on the flexible fit.

Reproduces ``examples.m2_misspecification.figures.save_result_figures``'s
``localisation_strong_sharp_flexible`` figure at the docs budget
(``_m2.DOC_BUDGET``). Full-budget numbers: ``m2_misspecification.rst``, "The
diagnostics agree", row ``strong_sharp`` (peak 12.8 at 0.86295 um), and the GP
table (0.0215 Jy at 0.00070 um).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _m2

from ampere.results import plot_gp_localisation

entry = _m2.runs()[("strong_sharp", "flexible")]
label = next(iter(entry["problem"].datasets))
plot_gp_localisation(entry["run"], datasets=[label])
