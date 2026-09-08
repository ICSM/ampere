Tutorials
=========

Ampere v2
---------

.. toctree::
   :maxdepth: 2

   m2_misspecification
   wstat_comparison

:doc:`m2_misspecification` is the flagship study: a deliberately misspecified
spectrum, fitted with and without the flexible likelihood, at three data
sizes and on all three backends, with the timings. :doc:`wstat_comparison`
works through registering a **user-defined likelihood family** — a profiled
Cash statistic with background — which is the extension point to reach for
when your data are not Gaussian.

Both have runnable counterparts in the repository:
``examples/m2_misspecification`` and ``examples/wstat_comparison.py``, each
covered by its own test suite so that neither can rot unnoticed.

Legacy tutorials
----------------

.. warning::

   These notebooks teach the **legacy** v1 API (``ampere.data``,
   ``ampere.models``, ``ampere.infer``), which is frozen — see :doc:`ampere`.
   They are kept because they are still the fullest worked examples of an SED
   fit and of neural posterior estimation in this repository, and because
   Phase 3 will replace them rather than delete them. They are rendered from
   their stored output; the documentation build does not execute them.

.. toctree::
   :maxdepth: 2

   notebooks/quickstart
   notebooks/Ampere_MBB_Example
   notebooks/Embedding_nets

Still to be written
-------------------

Combining different data types, conditional priors and arbitrary priors each
deserve a page of their own; :doc:`overview` and the API reference carry what
there is for now.
