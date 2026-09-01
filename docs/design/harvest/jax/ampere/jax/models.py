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

