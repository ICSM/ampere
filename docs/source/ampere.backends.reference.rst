ampere.backends.reference
=========================

Rung 1 of the capability ladder: pure numpy/scipy, no extra required, and
the base install's whole toolkit. It is also the conformance oracle — the
battery in ``tests/conformance/`` computes reference values here and holds
every other backend to them.

.. automodule:: ampere.backends.reference
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Re-exported from ``ampere.backends.reference``; shown under the modules that
define them.

.. autodata:: ampere.backends.reference.models.COORDINATE_UNIT
.. autodata:: ampere.backends.reference.models.FLUX_UNIT
.. autodata:: ampere.backends.reference.instrument.DETECTORS
