    """ JAX backend for Ampere.

    This module provides a JAX backend for Ampere. It builds on the Equinox 
    package which provides a torch-like interface for JAX, and the GPJax package
    which provides a Gaussian Process implementation for JAX.

    This module is intended to replace the current implementation of Ampere,
    which is built in pure Python and numpy. The JAX backend will provide
    significant performance improvements, and will also allow for the use of
    GPU acceleration, a wider variety of models and more scalable inference
    algorithms.

    For now, the recoomended way to use this module is to use the `ampere.jax`
    namespace, which will provide a drop-in replacement for the current Ampere
    implementation, i.e.
    
        ```python
        import ampere.jax as ampere
        ```

    This will provide a drop-in replacement for the current Ampere implementation,
    with a modified API aiming to modernise the architecture, but using JAX for 
    the underlying computations. The current implementation will be deprecated
    in v1.0.0, being renamed to `ampere.numpy`, but will be maintained for a
    while to allow for a smooth transition.

    The JAX backend will also be complemented by a torch backend, which will
    be intended to provide the same functionality for torch-based models. 
    """