"""Repository-wide pytest setup: the one place the root goes on ``sys.path``.

**W5.28(k), deciding and recording the question the wave-3/4 reviews asked
of the sys.path block this file replaces.** ``tests/`` is not, and stays
not, an importable package: no ``__init__.py`` anywhere under it, and none
added here. Making it one (so that pytest's own ``rootdir``-walking import
machinery would plant the repository root on ``sys.path`` for free) was
considered and rejected -- not because it would not work today, but because
it would be re-deriving, implicitly and fragile-ly, exactly the one fact
this file now states explicitly: several suites (``tests/m2``,
``tests/interferometry``, ``tests/benchmarks``, the many-lines rows of
``tests/inference/test_nuts.py``) import :mod:`examples` — a directory the
built wheel deliberately does not ship (``pyproject.toml``'s ``include =
["ampere*"]``) — so the repository root must reach ``sys.path`` by a route
that does not depend on how ``ampere`` itself was installed. Two of those
suites' own docstrings had already reasoned this through and landed on the
same answer *without* being able to point at each other's file, which is
what "the conftest sys.path pattern is used four times" was: not a defect
to design away, but one decision, taken independently, four times over,
because there was nowhere to write it down once. This file is that once.

**Why explicit, rather than pytest's own rootdir insertion.** Two concrete
ways the implicit route was seen to fail, both already lived through:
``tests/m2/conftest.py``'s original wording -- "an editable install's path
hook reaching this far would break the first time somebody ran it against a
wheel" (pytest's own insertion is keyed to *where a test file lives*, which
answers a different question from "does :mod:`examples` exist on this
machine at all") -- and ``tests/benchmarks/conftest.py``'s -- "``pixi run
bench`` invokes the bare ``pytest`` entry point and that inserts the test
file's directory and not the project's" (a package-chain walk changes
*which* directory pytest inserts, but the entry point still decides whether
it inserts one for :mod:`examples` at all). A conftest as a plain file
sidesteps both: pytest finds and runs every ``conftest.py`` from
``rootdir`` down to a collected test's own directory purely by walking the
filesystem, never by asking Python's import system where a package's root
is -- which is a guarantee ``__init__.py`` files would not add anything to,
and unconditional collection of *this* file already gives for free.

**What this replaces.** The identical block --
``_ROOT = pathlib.Path(__file__).resolve().parents[N]; if str(_ROOT) not in
sys.path: sys.path.insert(0, str(_ROOT))`` -- used to be repeated in
``tests/m2/conftest.py``, ``tests/interferometry/conftest.py``,
``tests/benchmarks/conftest.py``, ``tests/benchmarks/test_m2_sampling.py``,
``tests/benchmarks/test_m2_likelihood.py`` and
``tests/inference/test_nuts.py``, each paying to rediscover the same
reasoning. pytest imports this file before any test or conftest below it in
the tree, so one insertion here reaches every one of them; each has been
trimmed to rely on it rather than repeat it.

**What this does not replace.** ``tests/inference/test_interferometry.py``
puts a *different* directory (``tests/backends``, for a sibling fixture
module) on ``sys.path``, for a different reason -- a same-tree neighbour,
not the repository root -- and keeps doing so unchanged; the ``examples/``
single-script loaders (``tests/examples/test_cached_fit.py`` and its
siblings) insert and then remove a path around one dynamic import by name,
which is a different problem (loading one standalone script as a module)
solved a different way (temporarily, per import, because two such scripts
could collide by basename) and has nothing to gain from a permanent,
process-wide insertion.
"""

from __future__ import annotations

import pathlib
import sys

_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
