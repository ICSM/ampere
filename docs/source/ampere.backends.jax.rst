ampere.backends.jax
===================

Rung 2 of the capability ladder, beside torch: numpyro distributions,
equinox pytrees, gradients, ``jit``, and a per-instance device.

Requires the ``jax`` extra — from a clone, since ampere is not on PyPI
(:doc:`install`)::

    pip install -e ".[jax]"

**float64 is a policy here, and you must turn it on yourself.** jax's
``jax_enable_x64`` is process-global state addressed to the application
author rather than to a library, so ampere never flips it for you; every
model, instrument step, kernel, solver and lowered parameter set in this
package raises at construction if the flag is off. Call
:func:`~ampere.backends.jax.configure_x64` before any jax work.

.. note::

   This page is built in an environment that does not have jax installed,
   with ``jax``, ``numpyro`` and ``equinox`` mocked (see
   ``docs/source/conf.py``). The signatures and the prose are the real ones;
   runtime-computed defaults are not. Build with ``pixi run -e jax docs`` to
   render the page against the real library.

.. automodule:: ampere.backends.jax
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: LoweringFallbackWarning

``LoweringFallbackWarning`` is re-exported here for the import path a jax user
reaches for, but it is one shared class across the backends and is documented
once, as :class:`ampere.core.LoweringFallbackWarning`.

Constants
---------

Re-exported from ``ampere.backends.jax``; shown under the modules that define
them. ``PRECISIONS`` is rendered from mocked ``jax`` objects here, so its
displayed value is a placeholder — see the note above.

.. autodata:: ampere.backends.jax.models.COORDINATE_UNIT
.. autodata:: ampere.backends.jax.models.FLUX_UNIT
.. autodata:: ampere.backends.jax.gp.PRECISIONS
