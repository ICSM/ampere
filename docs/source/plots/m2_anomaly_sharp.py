"""The anomaly score of ``strong_sharp``'s flexible fit.

``plot_anomaly_score(gp_localisation_score(...))``: reproduces
``examples.m2_misspecification.figures.save_result_figures``'s
``anomaly_strong_sharp_flexible`` figure at the docs budget
(``_m2.DOC_BUDGET``). Full-budget numbers: ``m2_misspecification.rst``, "The
diagnostics agree", row ``strong_sharp`` (peak 12.8 at 0.86295 um).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _m2

from ampere.results import gp_localisation_score, plot_anomaly_score

entry = _m2.runs()[("strong_sharp", "flexible")]
label = next(iter(entry["problem"].datasets))
plot_anomaly_score(gp_localisation_score(entry["run"], dataset=label))
