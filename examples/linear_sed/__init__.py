"""A linear model, a photometric catalogue and two Spitzer IRS chunks.

The v2 twin of ``examples/minimal_working_example.py`` and its ``_dynesty``,
``_zeus``, ``_sbi`` and ``_sbi_embedding`` variants (W6.13 (1)) -- see
:mod:`.linear_sed` for the model, the data and the engine flag, and
:mod:`.generators` for the synthetic data and the tracked CASSIS file's
wavelength grids.

Deliberately minimal, like :mod:`examples.sed_composition`: this file does
not re-export its submodules' names. Import them directly --
``from examples.linear_sed import generators`` or
``from examples.linear_sed.linear_sed import build_problem``.
"""

from __future__ import annotations
