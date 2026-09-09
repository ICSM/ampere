ampere.results
==============

One results format for every run, and everything that reads it. ArviZ is the
single format: every engine emits it, all plotting is written once against
it, and serialisation is netCDF. Nothing here is sampler-specific.

The namespace covers provenance (the hashes and ``ampere_*`` attributes a run
carries), emission and netCDF round-tripping, plain-data serialisation for
containers and model results, the derived groups a run does not store by
default, the post-fit diagnostics as numbers, the plotting surface, and
trained-artefact caching keyed on those same provenance hashes.

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
.. autodata:: ampere.results.artefacts.ARTEFACT_CACHE_SCHEMA_VERSION

Training sets
~~~~~~~~~~~~~

.. autodata:: ampere.results.training.SAMPLE_DIM
.. autodata:: ampere.results.training.COORDINATES_GROUP
.. autodata:: ampere.results.training.OBSERVATIONS_GROUP

Trained-artefact caching
~~~~~~~~~~~~~~~~~~~~~~~~

``DEVELOPMENT_PLAN.md`` §7's trap — an SBI posterior, embedding net or
emulator reused against a problem it was not trained on is silently wrong —
closed by a content-addressed store whose key is built entirely from hashes
:mod:`ampere.results.provenance` already knows how to compute:
:func:`~ampere.results.artefact_key` builds one :class:`~ampere.results.ArtefactKey`
from a :class:`~ampere.core.dataset.FittingProblem` and a run's own settings
(the method, the estimator architecture, the budget, the rounds, the encoding
layout and the ``sbi``/torch versions), and
:class:`~ampere.results.ArtefactStore` is the ``get``/``put``/``train_or_load``
seam an inference engine (``SBIEngine(cache=...)``, W3.5) calls around
training. A stored artefact that cannot be trusted — a hand-edited sidecar, a
corrupted pickle — is a silent miss with an :class:`~ampere.results.ArtefactCacheWarning`,
never a raised error.

Diagnostics and plotting
~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: ampere.results.diagnostics.WHITENESS_STREAM
.. autodata:: ampere.results.diagnostics.GP_LOCALISATION_PROVENANCE
.. autodata:: ampere.results.plots.GP_LOCALISATION_CAVEAT
.. autodata:: ampere.results.plots.MAX_CORNER_VARIABLES
.. autodata:: ampere.results.plots.MAX_TRACE_VARIABLES
