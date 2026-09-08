ampere.results
==============

One results format for every run, and everything that reads it. ArviZ is the
single format: every engine emits it, all plotting is written once against
it, and serialisation is netCDF. Nothing here is sampler-specific.

The namespace covers provenance (the hashes and ``ampere_*`` attributes a run
carries), emission and netCDF round-tripping, plain-data serialisation for
containers and model results, the derived groups a run does not store by
default, the post-fit diagnostics as numbers, and the plotting surface.

.. automodule:: ampere.results
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Every name below is re-exported from ``ampere.results`` — ``from
ampere.results import POSTERIOR_GROUP`` — and is shown here under the module
that defines it, which is where its documentation lives.

The group and dimension names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: ampere.results.emission.CHAIN_DIM
.. autodata:: ampere.results.emission.DRAW_DIM
.. autodata:: ampere.results.emission.POSTERIOR_GROUP
.. autodata:: ampere.results.emission.SAMPLE_STATS_GROUP
.. autodata:: ampere.results.emission.LOG_LIKELIHOOD_GROUP
.. autodata:: ampere.results.emission.LOG_LIKELIHOOD_DECOMPOSITION
.. autodata:: ampere.results.emission.OBSERVED_DATA_GROUP
.. autodata:: ampere.results.emission.CONSTANT_DATA_GROUP
.. autodata:: ampere.results.derived.POSTERIOR_PREDICTIVE_GROUP
.. autodata:: ampere.results.derived.RESIDUALS_GROUP
.. autodata:: ampere.results.derived.POINTWISE_LOG_LIKELIHOOD_GROUP
.. autodata:: ampere.results.derived.FACTORISED_DECOMPOSITION
.. autodata:: ampere.results.derived.CONDITIONAL_LOO_DECOMPOSITION
.. autodata:: ampere.results.derived.GP_LOCALISATION_GROUP

Provenance and schema versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: ampere.results.provenance.ATTR_PREFIX
.. autodata:: ampere.results.provenance.PROVENANCE_SCHEMA_VERSION
.. autodata:: ampere.results.serialisation.CONTAINER_SCHEMA_VERSION
.. autodata:: ampere.results.training.TRAINING_SET_SCHEMA_VERSION

Training sets
~~~~~~~~~~~~~

.. autodata:: ampere.results.training.SAMPLE_DIM
.. autodata:: ampere.results.training.COORDINATES_GROUP
.. autodata:: ampere.results.training.OBSERVATIONS_GROUP

Diagnostics and plotting
~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: ampere.results.diagnostics.WHITENESS_STREAM
.. autodata:: ampere.results.diagnostics.GP_LOCALISATION_PROVENANCE
.. autodata:: ampere.results.plots.GP_LOCALISATION_CAVEAT
.. autodata:: ampere.results.plots.MAX_CORNER_VARIABLES
.. autodata:: ampere.results.plots.MAX_TRACE_VARIABLES
