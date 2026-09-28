"""A four-parameter modified blackbody fitted to ten-band far-infrared photometry.

The v2 twin of ``examples/examples_paper/modifiedblackbody.py`` (W6.13 (3))
-- see :mod:`.modified_blackbody` for the model, the data and the
``--engine``/``--all`` flags, and :mod:`.generators` for the synthetic data.

Deliberately minimal, like :mod:`examples.sed_composition` and
:mod:`examples.linear_sed`: this file does not re-export its submodules'
names. Import them directly --
``from examples.modified_blackbody import generators`` or
``from examples.modified_blackbody.modified_blackbody import build_problem``.
"""

from __future__ import annotations
