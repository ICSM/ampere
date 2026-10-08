"""The gross case's GP localisation: the same power law with the permissive kernel prior.

Reproduces ``gp_localisation`` + ``plot_gp_localisation`` on the
``gp-wide/optimiser`` fit of ``docs/design/walkthroughs/persona_b2.py`` (the
memo ran ``gp-wide/prior`` through the plot; the prior start does not converge
at a docs budget), at the docs budget, with ``gp_localisation`` over draws
thinned by ``_irs.THIN``. Full-budget numbers: the user-journeys memo,
Appendix B.2 (GP wide from the optimiser: index 0.9 +- 1.0, norm on its bound,
GP 0.033 Jy at 9.4 um).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = [".", "docs/source/plots"]
import _irs

from ampere.results import plot_gp_localisation

localisation, _tree, _problem = _irs.gp_wide()
plot_gp_localisation(localisation)
