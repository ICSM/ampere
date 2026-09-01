"""
ampere.backends
===============

Optional inference backends for ampere.

Sub-packages
------------
legacy
    Thin re-exports of the numpy-based samplers (emcee, dynesty, zeus) under a
    stable import path.  No new functionality; useful when you want to write
    code that is independent of which sampler is used.
jax (future)
    JAX / numpyro backend – not yet implemented (see Issue 7).

Usage
-----
The simplest way to use the legacy backend::

    from ampere.backends.legacy import EmceeSearch, DynestySearch

Or import directly from ``ampere.infer`` as before::

    from ampere.infer import EmceeSearch
"""
