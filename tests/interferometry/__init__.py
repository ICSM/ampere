"""W4.4's CI-runnable assertions: the interferometric-visibility study.

W5.26: this package marker is what lets ``tests/interferometry`` join
``test-all`` at all. Without it, pytest's default "prepend" import mode
names each test module by its bare basename, and three of this directory's
files collide with an identically-named file elsewhere in ``test-all``'s
one process (``test_calibration.py`` with ``tests/results``,
``test_engines.py`` with ``tests/inference``, ``test_model.py`` with
``tests/m2``) -- "import file mismatch" at collection. ``tests/conformance``
and ``tests/astrometry`` already carry this same marker for the same
reason.
"""

from __future__ import annotations
