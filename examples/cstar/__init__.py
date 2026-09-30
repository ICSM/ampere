"""The carbon star on Hyperion, fitted by simulation-based inference.

The v2 twin of ``examples/cstar_model_test_sbi_v2.py`` and its ``_embedding``
variant (W6.13 (5)) -- see :mod:`.cstar` for the model, the data, the engine
and the recorded runs, :mod:`.dust` for the ``miepython`` opacities that
replace the legacy ``bhmie`` step, and :mod:`.generators` for the tracked data
and the synthetic observation.

Needs the pixi ``hyperion`` environment (``pixi install -e hyperion``):
Hyperion and miepython are example-only requirements, not an ampere extra.
Deliberately minimal: import the submodules directly --
``from examples.cstar.cstar import build_problem``.
"""

from __future__ import annotations
