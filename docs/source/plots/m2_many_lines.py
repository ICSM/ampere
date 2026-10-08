"""The many-lines figure: one deviation with two length scales, three kernels on it.

Reproduces ``python -m examples.m2_misspecification.many_lines --figures DIR``
(``many_lines.compare`` then ``figures.figure_many_lines``) at the docs budget:
26 walkers x 260 steps (130 discarded), 200 points, the study's own seed,
against the 32 x 1 100 (550 discarded) of ``MANY_LINES_EMCEE``. Full-budget
numbers: ``m2_misspecification.rst``, "When one length scale is not enough"
(worst offsets 2.84, 0.92, 0.78) and "A sum of noise components needs a sparsity
guard" (flat against horseshoe, 0.0118 against 0.0077 Jy; min/max 0.254 against
0.087).
"""

# ruff: noqa: E402

import sys

sys.path[:0] = ["."]

from examples.m2_misspecification import many_lines
from examples.m2_misspecification.figures import figure_many_lines
from examples.m2_misspecification.study import EmceeBudget

comparison = many_lines.compare(budget=EmceeBudget(walkers=26, steps=260, burn_in=130))
figure_many_lines(comparison.entries, comparison.data, band=comparison.band)
