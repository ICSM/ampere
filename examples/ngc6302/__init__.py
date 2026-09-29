"""The Kemper et al. (2002) NGC 6302 two-shell dust model, fitted with v2.

The v2 twin of ``examples/NGC6302.py`` and ``examples/NGC6302_zeus.py``
(W6.13 (2)) -- see :mod:`.ngc6302` for the model, the data and the
``--engine`` flag, :mod:`.generators` for the synthetic data, and
:mod:`.dust_mass` for the post-fit dust-mass table (the v2 continuation of
``examples/NGC6302-calculate-dust-mass.py``).

Deliberately minimal, like :mod:`examples.sed_composition`,
:mod:`examples.linear_sed` and :mod:`examples.modified_blackbody`: this file
does not re-export its submodules' names. Import them directly --
``from examples.ngc6302 import generators`` or
``from examples.ngc6302.ngc6302 import build_problem``.
"""

from __future__ import annotations
