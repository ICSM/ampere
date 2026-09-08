ampere.backends
===============

.. automodule:: ampere.backends
   :members:

A backend is precisely three things: an array library for writing models and
transformations, a set of GP solver implementations, and the extra inference
engines its differentiability unlocks. Kernels, containers, likelihood
families, datasets and the fitting problem itself are **not** in a backend —
they are backend-neutral and live in :mod:`ampere.core`.

The three form an escalating capability ladder rather than a set of peers;
see :doc:`overview` for what that means when you compose a problem, and for
the one-backend rule every part of a problem has to satisfy.

.. toctree::
   :maxdepth: 2

   ampere.backends.reference
   ampere.backends.torch
   ampere.backends.jax
