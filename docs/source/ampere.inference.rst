ampere.inference
================

The engines, written once against the fitting-problem surface and run on
every backend. This namespace **imports no backend**, in any module, at any
depth; its only ampere imports are :mod:`ampere.core` and
:mod:`ampere.results`. The gradient-based drivers reach a backend's
differentiable density through :func:`ampere.core.realise`, which importing
that backend registered.

There is no way to sample through these drivers and not get a stored run:
``run()`` returns the ArviZ ``DataTree``, with per-draw ``log_prior`` and
``log_likelihood``, the per-dataset decomposition, the observed data and the
full provenance attrs.

.. automodule:: ampere.inference
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Re-exported from ``ampere.inference``; shown under the module that defines it.

.. autodata:: ampere.inference.engine.DEFAULT_CACHE_SIZE
