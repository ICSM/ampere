"""Results: one format for every run, and everything that reads it.

``DEVELOPMENT_PLAN.md`` §4.6 makes this namespace's remit small and absolute:
**ArviZ is the single results format**, every engine emits it, all plotting is
written once against it, and serialisation is netCDF. Nothing here is
sampler-specific, because "sampler plotting parity keeps regressing" is the
problem this namespace exists to end.

The contract is ``docs/design/contracts/results.md``. In outline:

``provenance``
    The hashing recipe and the ``ampere_*`` attributes every run carries — the
    versions, the seed, the parameter-spec hash (load-bearing, not book-keeping:
    ``lowering.md`` §9.2), the data hashes, the composed bindings, the failure
    counts.
``emission``
    :func:`emit` builds the run from draws plus one
    :class:`~ampere.core.dataset.Evaluation` per draw; :class:`DrawRecorder`
    accumulates them a draw at a time; :func:`to_netcdf` / :func:`from_netcdf`
    put it on disk and get it back.
``serialisation``
    Plain-data forms for containers and ``ModelResult``\\ s — the ``(θ, result)``
    pairs design horizon (c) wants for emulator training sets, and the answer to
    ``results_schema.md`` §17 question 6.
``derived``
    The groups a run does *not* store by default (posterior-predictive
    replicates, signed residuals) and the documented rule for computing them.
``diagnostics``
    ``diagnostics.md``'s post-fit families, as numbers rather than pictures:
    family B's separation-binned, permutation-calibrated residual-whiteness
    test and its chi-square posterior-predictive p-value, and family C's
    conversion of the conditioned GP mean into the shared
    :class:`~ampere.core.AnomalyScore`.
``plots``
    The plotting surface, including ``diagnostics.md``'s families B and C.
``artefacts``
    Trained-artefact caching keyed on the problem's own provenance hashes
    (``DEVELOPMENT_PLAN.md`` §7's trap) — :class:`ArtefactStore`, its
    :class:`ArtefactKey` and :func:`artefact_key`, and the ``train_or_load``
    seam an inference engine calls around training (W3.5).

Dependencies
------------
arviz and a netCDF engine (``h5netcdf``) are **base dependencies** as of W2.2,
which is when ``ampere.inference``'s engine drivers made emission load-bearing:
Peter's 2026-09-03 ruling on ``results.md`` §15 R1 timed the promotion for
"Phase 2's engine drivers […] when a user can first emit a run", and that is
that moment. ``architecture.md`` §3's extras table records the same.

The import stays **lazy** all the same, inside the functions that need it,
still raising :class:`~ampere.core.exceptions.OptionalDependencyError` on use.
Two reasons, neither of them the original one. ``import ampere.results`` is on
the path of anything that touches this namespace at all, and arviz pulls in
xarray, pandas and its own plotting stack — perhaps a second of import time
bought for nothing by a caller who only wanted
:func:`~ampere.results.provenance.hash_container`. And a base dependency can
still be *absent*: an environment assembled by hand, or a partially installed
one, should say which package is missing and how to get it rather than fail
with a bare ``ModuleNotFoundError`` from three frames down.
"""

from __future__ import annotations

from ._plotting import figure_metadata
from .artefacts import (
    ARTEFACT_CACHE_SCHEMA_VERSION,
    ArtefactCacheWarning,
    ArtefactKey,
    ArtefactStore,
    artefact_key,
)
from .derived import (
    CONDITIONAL_LOO_DECOMPOSITION,
    FACTORISED_DECOMPOSITION,
    GP_LOCALISATION_GROUP,
    POINTWISE_LOG_LIKELIHOOD_GROUP,
    POSTERIOR_PREDICTIVE_GROUP,
    RESIDUALS_GROUP,
    add_pointwise_log_likelihood,
    add_posterior_predictive,
    add_residuals,
    gp_localisation,
    pointwise_as_log_likelihood,
)
from .diagnostics import (
    GP_LOCALISATION_PROVENANCE,
    WHITENESS_STREAM,
    ChiSquareCheck,
    WhitenessTest,
    chi_square_pvalue,
    gp_localisation_datasets,
    gp_localisation_score,
    residual_whiteness,
    separation_binned_autocorrelation,
)
from .emission import (
    CHAIN_DIM,
    CONSTANT_DATA_GROUP,
    DRAW_DIM,
    LOG_LIKELIHOOD_DECOMPOSITION,
    LOG_LIKELIHOOD_GROUP,
    OBSERVED_DATA_GROUP,
    POSTERIOR_GROUP,
    SAMPLE_STATS_GROUP,
    DrawRecorder,
    emit,
    from_netcdf,
    to_netcdf,
)
from .exceptions import ResultsError
from .plots import (
    GP_LOCALISATION_CAVEAT,
    MAX_CORNER_VARIABLES,
    MAX_TRACE_VARIABLES,
    AnomalyScoreLike,
    gp_localisation_caveat,
    plot_anomaly_score,
    plot_corner,
    plot_gp_localisation,
    plot_posterior_predictive,
    plot_residuals,
    plot_trace,
)
from .provenance import (
    ATTR_PREFIX,
    PROVENANCE_SCHEMA_VERSION,
    buffer_fingerprint,
    canonical_json,
    container_fingerprint,
    dataset_fingerprint,
    describe_likelihood,
    digest,
    hash_array,
    hash_container,
    hash_of,
    model_fingerprint,
    model_identity_hash,
    neutral_model_identity,
    package_versions,
    problem_fingerprint,
    provenance_attrs,
    spec_hashes,
)
from .training import (
    COORDINATES_GROUP,
    OBSERVATIONS_GROUP,
    SAMPLE_DIM,
    TRAINING_SET_SCHEMA_VERSION,
    TrainingSet,
    append_training_set,
    read_training_set,
    write_training_set,
)
from .serialisation import (
    CONTAINER_SCHEMA_VERSION,
    container_from_dict,
    container_to_dict,
    kind_named,
    model_result_from_dict,
    model_result_to_dict,
    register_kind,
    registered_kinds,
    training_pair_from_dict,
    training_pair_to_dict,
)

__all__ = [
    "ARTEFACT_CACHE_SCHEMA_VERSION",
    "ATTR_PREFIX",
    "CHAIN_DIM",
    "CONDITIONAL_LOO_DECOMPOSITION",
    "CONSTANT_DATA_GROUP",
    "CONTAINER_SCHEMA_VERSION",
    "COORDINATES_GROUP",
    "DRAW_DIM",
    "FACTORISED_DECOMPOSITION",
    "GP_LOCALISATION_CAVEAT",
    "GP_LOCALISATION_GROUP",
    "GP_LOCALISATION_PROVENANCE",
    "LOG_LIKELIHOOD_DECOMPOSITION",
    "LOG_LIKELIHOOD_GROUP",
    "MAX_CORNER_VARIABLES",
    "MAX_TRACE_VARIABLES",
    "OBSERVATIONS_GROUP",
    "OBSERVED_DATA_GROUP",
    "POINTWISE_LOG_LIKELIHOOD_GROUP",
    "POSTERIOR_GROUP",
    "POSTERIOR_PREDICTIVE_GROUP",
    "PROVENANCE_SCHEMA_VERSION",
    "RESIDUALS_GROUP",
    "SAMPLE_DIM",
    "SAMPLE_STATS_GROUP",
    "TRAINING_SET_SCHEMA_VERSION",
    "WHITENESS_STREAM",
    "AnomalyScoreLike",
    "ArtefactCacheWarning",
    "ArtefactKey",
    "ArtefactStore",
    "ChiSquareCheck",
    "DrawRecorder",
    "ResultsError",
    "TrainingSet",
    "WhitenessTest",
    "add_pointwise_log_likelihood",
    "add_posterior_predictive",
    "add_residuals",
    "append_training_set",
    "artefact_key",
    "buffer_fingerprint",
    "canonical_json",
    "chi_square_pvalue",
    "container_fingerprint",
    "container_from_dict",
    "container_to_dict",
    "dataset_fingerprint",
    "describe_likelihood",
    "digest",
    "emit",
    "figure_metadata",
    "from_netcdf",
    "gp_localisation",
    "gp_localisation_caveat",
    "gp_localisation_datasets",
    "gp_localisation_score",
    "hash_array",
    "hash_container",
    "hash_of",
    "kind_named",
    "model_fingerprint",
    "model_identity_hash",
    "model_result_from_dict",
    "model_result_to_dict",
    "neutral_model_identity",
    "package_versions",
    "plot_anomaly_score",
    "plot_corner",
    "plot_gp_localisation",
    "plot_posterior_predictive",
    "plot_residuals",
    "plot_trace",
    "pointwise_as_log_likelihood",
    "problem_fingerprint",
    "provenance_attrs",
    "read_training_set",
    "register_kind",
    "registered_kinds",
    "residual_whiteness",
    "separation_binned_autocorrelation",
    "spec_hashes",
    "to_netcdf",
    "training_pair_from_dict",
    "training_pair_to_dict",
    "write_training_set",
]
