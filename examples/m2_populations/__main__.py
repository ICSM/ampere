"""Run the comparison and print its table: ``python -m examples.m2_populations``.

The same code ``tests/m2/test_populations.py`` and the documentation figure
call, and nothing written anywhere unless ``--figures`` names a directory::

    python -m examples.m2_populations --backend jax --quick   # the per-PR budget, NUTS
    python -m examples.m2_populations --backend jax           # the marked budget
    python -m examples.m2_populations                         # emcee, ~20 min per arm
    python -m examples.m2_populations --arms independent population

The default budget on the reference path is
:data:`~examples.m2_populations.populations.POPULATION_EMCEE` (112 walkers x
1 600 steps): the problem has fifty dimensions, which an ensemble sampler pays
for. NUTS on torch or jax is the tool for it, and ``--quick`` there is
``tests/m2``'s own budget.
"""

from __future__ import annotations

import sys

from .populations import main

if __name__ == "__main__":
    sys.exit(main())
