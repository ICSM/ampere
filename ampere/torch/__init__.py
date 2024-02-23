""" Torch backend for Ampere. 

This module provides a torch backend for Ampere. It prarallels the `ampere.jax`
module, and provides a drop-in replacement for the current Ampere implementation
using (G)PyTorch for the underlying computations. While the JAX backend is
intended to replace the current implementation, the torch backend is intended
to complement it and provide the same functionality for torch-based models and 
inference algorithms.
"""