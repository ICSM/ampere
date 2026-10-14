"""The population scenario's collection: eight spectra, three deviated, one noisy.

Reproduces the left panel of ``python -m examples.m2_populations --figures DIR``
(``populations.generate_population`` then ``figures.figure_population``): the
data the three declarations of ``population.rst``'s table are fitted to, at the
driver's own seeds (spectrum *i* at ``DATA_SEED + i``, 200 points). No sampling
— the table's numbers are measured by ``tests/m2/test_populations.py`` at
``POPULATION_NUTS``, which a docs build cannot afford.
"""

# ruff: noqa: E402

import sys

sys.path[:0] = ["."]

from examples.m2_populations.figures import figure_population
from examples.m2_populations.populations import generate_population

figure_population(generate_population())
