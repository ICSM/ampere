API reference
=============

There are two APIs in this package, and they do not talk to each other.

**Ampere v2** is the redesign: a backend-neutral core of frozen contracts
(:mod:`ampere.core`), one package per array library implementing them
(:mod:`ampere.backends.reference`, :mod:`ampere.backends.torch`,
:mod:`ampere.backends.jax`), a set of engines written once against the
contracts (:mod:`ampere.inference`), and a single results format everything
emits and everything reads (:mod:`ampere.results`). :doc:`overview` explains
how the pieces fit together; these pages are the reference.

**Legacy ampere** is the v1 code that produced the published science —
``ampere.data``, ``ampere.models``, ``ampere.infer``. It is **frozen**: it
still runs, and nothing in it is being changed. It implements none of the v2
contracts and is not interoperable with them.

Ampere v2
---------

.. toctree::
   :maxdepth: 2

   ampere.core
   ampere.backends
   ampere.inference
   ampere.results

Legacy (frozen)
---------------

.. toctree::
   :maxdepth: 2

   ampere
