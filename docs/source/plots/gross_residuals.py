"""The gross case's residual plot: a power law fitted to the IRS spectrum of PG 1011-040.

Reproduces the ``plot_residuals(add_residuals(r, p))`` call on the
``independent/optimiser`` fit of ``docs/design/walkthroughs/persona_b2.py``, at
the docs budget (``_irs.DOC_WALKERS`` x ``_irs.DOC_STEPS``). Full-budget
numbers: the user-journeys memo, Appendix B.2 (index 1.111 +- 0.003, whiteness
Q = 1490, p = 0.005).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _irs

from ampere.results import plot_residuals

tree, _problem = _irs.independent()
plot_residuals(tree)
