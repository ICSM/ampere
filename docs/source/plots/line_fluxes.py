"""The line-fluxes study's flexible arm: the GP localisation over rest wavelength.

Reproduces ``plot_gp_localisation(gp_localisation(run, problem))`` on the
one-temperature ladder fitted under the flexible likelihood
(``examples/line_fluxes``, arm ``"flexible"``), at the study's per-PR emcee
budget (``study.CI_BUDGET``: about twelve seconds), so the documentation
build stays well inside its minute. The conditioned GP mean is largest at the
high-excitation [Ar II] 6.99 line, where the missing hot component contributes
most.
"""

# ruff: noqa: E402

import sys
import warnings

sys.path[:0] = ["."]

from ampere.results import gp_localisation, plot_gp_localisation
from examples.line_fluxes import study

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    problem = study.build_problem("reference", "flexible")
    tree = study.run(problem, study.CI_BUDGET)
    plot_gp_localisation(gp_localisation(tree, problem))
