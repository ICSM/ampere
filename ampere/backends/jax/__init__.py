"""
ampere.backends.jax
===================

JAX / numpyro inference backend – **not yet implemented**.

This package is a placeholder that establishes the import path and documents
the intended integration pattern so that future work can be added without
changing the public API.

Intended integration pattern (Issue 7 / Issue 8)
-------------------------------------------------
The JAX backend will expose a ``NumpyroSearch`` class that:

1. Accepts any ampere :class:`~ampere.core.model.Model` whose
   :meth:`~ampere.core.model.Model.__call__` returns a
   :class:`~ampere.models.results.ModelResults`.
2. Builds a numpyro probabilistic program by:
   - Sampling each parameter declared in ``model.parameters`` from its
     ``prior`` (converted to a numpyro distribution).
   - Calling ``model(**params)`` with the sampled values.
   - Evaluating ``data.lnlike(params, result)`` as the likelihood.
3. Uses the ``model.population_model()`` hook (Issue 8) to optionally wrap the
   single-object model in a numpyro ``plate`` for hierarchical / population
   inference.

Dependency
----------
This backend requires ``jax``, ``numpyro``, and optionally ``blackjax``::

    pip install ampere[jax]

Raises
------
ImportError
    If JAX or numpyro are not installed when importing from this package.
"""

try:
    import jax  # noqa: F401
    import numpyro  # noqa: F401
except ImportError as _exc:
    raise ImportError(
        "The ampere JAX backend requires jax and numpyro. "
        "Install them with:  pip install ampere[jax]"
    ) from _exc

raise NotImplementedError(
    "The ampere JAX / numpyro backend is not yet implemented. "
    "It is planned for a future release (see IMPLEMENTATION_PLAN.md, Issue 7)."
)
