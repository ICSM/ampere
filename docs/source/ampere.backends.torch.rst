ampere.backends.torch
=====================

Rung 2 of the capability ladder: everything the reference backend does, plus
differentiability — and therefore NUTS, stochastic variational inference and
gradient-based optimisation — plus batching through ``torch.func.vmap``, and
a per-instance ``device=`` on every piece.

Requires the ``torch`` extra — from a clone, since ampere is not on PyPI
(:doc:`install`)::

    pip install -e ".[torch]"

.. note::

   This page is built in an environment that does not have torch installed,
   with ``torch`` and ``pyro`` mocked (see ``docs/source/conf.py``). The
   signatures and the prose are the real ones; runtime-computed defaults are
   not. Build with ``pixi run -e torch docs`` to render the page against the
   real library.

.. automodule:: ampere.backends.torch
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: LoweringFallbackWarning

``LoweringFallbackWarning`` is re-exported here for the import path a torch
user reaches for, but it is one shared class across the backends and is
documented once, as :class:`ampere.core.LoweringFallbackWarning`.

Constants
---------

Re-exported from ``ampere.backends.torch``; shown under the modules that
define them. The two ``DEFAULT_*`` values are rendered from mocked ``torch``
objects here, so their displayed values are placeholders — see the note above.

.. autodata:: ampere.backends.torch.models.COORDINATE_UNIT
.. autodata:: ampere.backends.torch.models.FLUX_UNIT
.. autodata:: ampere.backends.torch.instrument.DETECTORS
.. autodata:: ampere.backends.torch._config.DEFAULT_DTYPE
.. autodata:: ampere.backends.torch._config.DEFAULT_DEVICE

.. automodule:: ampere.backends.torch.sharding
   :members:
   :undoc-members:
   :show-inheritance:
