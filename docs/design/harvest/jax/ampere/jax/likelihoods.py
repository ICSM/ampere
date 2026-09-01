# ============================================================================
# ASPIRATIONAL SKETCH -- harvested read-only from origin/jax (W0.7, 2026-09).
# DOES NOT IMPORT. MUST NOT BE MERGED AS-IS.
# ----------------------------------------------------------------------------
# ampere/jax/module.py (this file or one it imports) does
#     from equinox import Dataset as EqxDataset
#     from equinox import Parameter as EqxParameter
# Neither `equinox.Dataset` nor `equinox.Parameter` exists in the equinox
# public API (checked against equinox on PyPI as of this harvest) -- this
# import fails immediately. Every other module in ampere/jax/ imports
# (directly or transitively) from ampere/jax/module.py, so the whole
# ampere.jax sketch package fails to import.
# This is reference-only design intent for a future jax backend
# (see DEVELOPMENT_PLAN.md); treat it as sketch notes, not implementation.
# Full provenance: docs/design/harvest/jax/README.md
# ============================================================================

""" This module defines the base class for all `ampere.jax` likelihoods. 

The `Likelihood` class is a subclass of `Module`, and is intended to be used as
a base class for all likelihoods in the `ampere.jax` namespace. It provides a
`forward` method that should be implemented in a subclass, and a `log_prob`
method that should be implemented in a subclass.

This module also provides a few of the most important astrophyiscal likelihoods
as subclasses of `Likelihood`, including `GaussianLikelihood` and
`PoissonLikelihood`. These are intended to be used in data objects, and to be
extended to include more complex likelihoods as needed.
"""

from typing import Union, List, Dict, Any