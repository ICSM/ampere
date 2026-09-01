# ============================================================================
# ASPIRATIONAL SKETCH -- harvested read-only from origin/jax (W0.7, 2026-09).
# DOES NOT IMPORT AS A PACKAGE. MUST NOT BE MERGED AS-IS.
# ----------------------------------------------------------------------------
# ampere/torch/data.py (source branch filename: "ampere/torch/data,py" --
# comma typo, renamed here for clarity) does `from .likelihoods import
# Likelihood`, but ampere/torch/likelihoods.py in this same sketch is an
# empty file (0 bytes, no `Likelihood` symbol) -- that import fails.
# ampere/torch/module.py depends on `gpytorch` and `pyro`, optional
# dependencies not declared anywhere in this project's extras.
# This is reference-only design intent for a possible future torch backend
# (see DEVELOPMENT_PLAN.md); treat it as sketch notes, not implementation.
# Full provenance: docs/design/harvest/jax/README.md
# ============================================================================

""" Torch backend for Ampere. 

This module provides a torch backend for Ampere. It prarallels the `ampere.jax`
module, and provides a drop-in replacement for the current Ampere implementation
using (G)PyTorch for the underlying computations. While the JAX backend is
intended to replace the current implementation, the torch backend is intended
to complement it and provide the same functionality for torch-based models and 
inference algorithms.
"""