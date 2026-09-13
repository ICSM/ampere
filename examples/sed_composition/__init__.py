"""One ``ModifiedBlackBody``, two datasets: a spectrum and a photometric
catalogue on the same model channel.

Memo §7.1 (``docs/design/phase4_placement_memo.md``), expanded into a runnable
example and :doc:`the tutorial page </sed_composition>` (W4.11). See
:mod:`.generators` for how the synthetic data is made and why that needs
``negotiate`` + ``compile_for`` rather than the model's own grid, and
:mod:`.sed_composition` for the problem, the fit, and the backend flag.

Deliberately minimal: unlike ``examples.m2_misspecification``, this package
does not re-export its submodules' names here. ``__main__.py`` restricts
every BLAS thread pool to one thread before numpy is imported anywhere in the
process, and ``python -m examples.sed_composition`` imports this file
*before* ``__main__.py`` runs — so if this file imported :mod:`.generators`
or :mod:`.sed_composition` (both of which import numpy) eagerly, the
restriction would already be too late. Import the submodules directly:
``from examples.sed_composition import generators`` or
``from examples.sed_composition.sed_composition import build_problem``.
"""

from __future__ import annotations
