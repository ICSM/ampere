ampere.results
==============

One results format for every run, and everything that reads it. ArviZ is the
single format: every engine emits it, all plotting is written once against
it, and serialisation is netCDF. Nothing here is sampler-specific.

The namespace covers provenance (the hashes and ``ampere_*`` attributes a run
carries), emission and netCDF round-tripping, plain-data serialisation for
containers and model results, the derived groups a run does not store by
default, the post-fit diagnostics as numbers, the plotting surface, and
trained-artefact caching keyed on those same provenance hashes. Since W3.6 it
also covers ``diagnostics.md``'s fourth family — posterior calibration, in
:mod:`ampere.results.calibration`.

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
.. autodata:: ampere.results.calibration.CALIBRATION_GROUP
.. autodata:: ampere.results.calibration.SIMULATION_DIM
.. autodata:: ampere.results.calibration.PARAMETER_DIM
.. autodata:: ampere.results.calibration.LEVEL_DIM
.. autodata:: ampere.results.calibration.TARP_LEVEL_DIM

Provenance and schema versions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autodata:: ampere.results.provenance.ATTR_PREFIX
.. autodata:: ampere.results.provenance.PROVENANCE_SCHEMA_VERSION
.. autodata:: ampere.results.serialisation.CONTAINER_SCHEMA_VERSION
.. autodata:: ampere.results.training.TRAINING_SET_SCHEMA_VERSION
.. autodata:: ampere.results.artefacts.ARTEFACT_CACHE_SCHEMA_VERSION
.. autodata:: ampere.results.calibration.CALIBRATION_SCHEMA_VERSION

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
layout, the ``sbi``/torch versions and, for TMNRE, its ``marginals``,
``truncation_epsilon`` and ``sample_with``), and :class:`~ampere.results.ArtefactStore` is the
``get``/``put``/``train_or_load`` seam an inference engine
(``SBIEngine(cache=...)``, W3.5) calls around training. A stored artefact
that cannot be trusted — a hand-edited sidecar, a corrupted pickle — is a
silent miss with an :class:`~ampere.results.ArtefactCacheWarning`, never a
raised error.

**W3.12**: the key's ``model_hash`` field is :func:`~ampere.results.model_hash`
— every model's fingerprint minus its parameter declaration, every dataset's
fingerprint minus its observed data, plus the model bindings — promoted into
:mod:`ampere.results.provenance` and recorded on every run and training set
as ``ampere_model_hash`` (``PROVENANCE_SCHEMA_VERSION`` 6), so a likelihood
family, noise model, solver or kernel swap that leaves every parameter's name
and prior unchanged is a miss here too, not only in this store:
:func:`~ampere.results.append_training_set` now refuses to grow a training
set whose ``ampere_model_hash`` disagrees with the problem given — or that
predates the attribute altogether.

Posterior calibration (family D)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``diagnostics.md`` §11's fourth diagnostic family, landed at W3.6. Where
families A to C diagnose the *model*, this one diagnoses the *inference
machinery*: does the posterior an engine produces have the coverage it claims?
Two routes, one result. :meth:`ampere.inference.SBIEngine.calibrate` is the
cheap one — an amortised posterior re-conditioned on a fresh simulation batch,
with ``sbi``'s own ``run_sbc``/``check_sbc`` and ``run_tarp``/``check_tarp``
doing the arithmetic — and :func:`~ampere.results.sbc` is the general one: the
Talts et al. loop over ``simulate`` with a full fit per simulated dataset,
expensive by design, and the only route that can validate a *likelihood*
rather than a network. Both return the same :class:`xarray.Dataset`,
:func:`~ampere.results.attach_calibration` writes it into a run as its
``calibration`` group, and :func:`~ampere.results.plot_sbc_ranks` and
:func:`~ampere.results.plot_coverage` draw it.
``examples/wstat_comparison.py``'s ``--coverage`` study is the worked example.

.. autodata:: ampere.results.calibration.CALIBRATION_STREAM
.. autodata:: ampere.results.calibration.DEFAULT_LEVELS
.. autodata:: ampere.results.calibration.SBI_ROUTE
.. autodata:: ampere.results.calibration.REFIT_ROUTE

Diagnostics and plotting
~~~~~~~~~~~~~~~~~~~~~~~~

**W3.10**: :func:`~ampere.results.plot_corner` and
:func:`~ampere.results.plot_trace` **page** above
:data:`~ampere.results.plots.MAX_CORNER_VARIABLES` /
:data:`~ampere.results.plots.MAX_TRACE_VARIABLES` rather than refuse — a list
of figures each within the cap, in merged-name order, with an array-valued
block kept whole on one page where it fits on one at all — and warn loudly
with :class:`~ampere.results.ResultsWarning`, naming the page count, the cap
and the ``var_names=`` route to a smaller figure. ``max_variables=`` still
means "how many per page"; ``paginate=False`` restores the pre-W3.10 refusal
for a caller who wants one figure or a hard failure; a call whose columns fit
the cap always returns a single :class:`~matplotlib.figure.Figure`, paginated
or not, and only a call that actually pages returns a ``list`` of them, each
carrying ``"page"`` as ``"i of n"`` in its
:func:`~ampere.results.figure_metadata`.

**W5.0**: :func:`~ampere.results.plot_trace` and the new
:func:`~ampere.results.summary` (a thin wrapper around :func:`arviz.summary`)
both warn, once and loudly, when a run's ``ampere_approximation`` root
attribute is not ``"none"`` — an R-hat, an ESS or a trace's shape describes
the optimiser or the density estimator that produced a
:class:`~ampere.inference.VIEngine` or :class:`~ampere.inference.SBIEngine`
run's draws, never sampling error against the target, and neither function
can tell a reader that on its own. :func:`~ampere.results.warn_if_approximate`
is the shared check both call.

**W5.3**: :func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals` and
:func:`~ampere.results.plot_gp_localisation` gain ``coordinate=`` (an axis
name, or a callable resolving one from a point kind's own axes — see
:func:`ampere.results._plotting.coordinate_of`) for a point kind with
several axes, such as :class:`~ampere.core.VisibilitySet` or
:class:`~ampere.core.ClosurePhases`; a kind may declare a default via
:attr:`~ampere.core.results_schema.FunctionSamples.PLOT_COORDINATE`, in which
case the argument is unnecessary. The same three, and
:func:`~ampere.results.add_posterior_predictive`,
:func:`~ampere.results.add_residuals` and
:func:`~ampere.results.gp_localisation`, gain ``component=`` — one of
:data:`~ampere.results.derived.COMPONENTS` — for a complex-valued dataset,
stored and read back as ``<label>_<component>``.
:func:`~ampere.results.plot_anomaly_score` needs neither: the
:class:`~ampere.core.AnomalyScore` it draws has already been reduced to one
real coordinate by :func:`~ampere.results.gp_localisation_score`, using the
same rule.

.. autodata:: ampere.results.diagnostics.WHITENESS_STREAM
.. autodata:: ampere.results.diagnostics.GP_LOCALISATION_PROVENANCE
.. autodata:: ampere.results.plots.GP_LOCALISATION_CAVEAT
.. autodata:: ampere.results.plots.MAX_CORNER_VARIABLES
.. autodata:: ampere.results.plots.MAX_RANK_PANELS
.. autodata:: ampere.results.plots.MAX_TRACE_VARIABLES
.. autodata:: ampere.results.derived.COMPONENTS
