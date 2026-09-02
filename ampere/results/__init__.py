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
``plots``
    The plotting surface, including ``diagnostics.md``'s families B and C.

Dependencies
------------
arviz is an **optional** dependency (``pip install "ampere[arviz]"``) and is
imported lazily, inside the functions that need it, raising
:class:`~ampere.core.exceptions.OptionalDependencyError` on use. That is
``architecture.md`` §4 rule 2 applied literally — it names ``ampere.results``
among the namespaces that must import optional dependencies lazily — and it is
what lets ``import ampere.results`` stay clean in the minimal-install job.
``architecture.md`` §3's extras table anticipates folding arviz into the base
install "when ``ampere.results`` lands"; ``docs/design/contracts/results.md``
§10 asks Peter to time that with Phase 2, when emission is actually load-bearing
and the netCDF engine (which arviz does not itself require) can be pinned with
it.
"""

from __future__ import annotations

from .derived import (
    GP_LOCALISATION_GROUP,
    POINTWISE_LOG_LIKELIHOOD_GROUP,
    POSTERIOR_PREDICTIVE_GROUP,
    RESIDUALS_GROUP,
    add_posterior_predictive,
    add_residuals,
    gp_localisation,
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
    package_versions,
    problem_fingerprint,
    provenance_attrs,
    spec_hashes,
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
    training_pair_to_dict,
)

__all__ = [
    "ATTR_PREFIX",
    "CHAIN_DIM",
    "CONSTANT_DATA_GROUP",
    "CONTAINER_SCHEMA_VERSION",
    "DRAW_DIM",
    "GP_LOCALISATION_CAVEAT",
    "GP_LOCALISATION_GROUP",
    "LOG_LIKELIHOOD_DECOMPOSITION",
    "LOG_LIKELIHOOD_GROUP",
    "OBSERVED_DATA_GROUP",
    "POINTWISE_LOG_LIKELIHOOD_GROUP",
    "POSTERIOR_GROUP",
    "POSTERIOR_PREDICTIVE_GROUP",
    "PROVENANCE_SCHEMA_VERSION",
    "RESIDUALS_GROUP",
    "SAMPLE_STATS_GROUP",
    "AnomalyScoreLike",
    "DrawRecorder",
    "ResultsError",
    "add_posterior_predictive",
    "add_residuals",
    "buffer_fingerprint",
    "canonical_json",
    "container_fingerprint",
    "container_from_dict",
    "container_to_dict",
    "dataset_fingerprint",
    "describe_likelihood",
    "digest",
    "emit",
    "from_netcdf",
    "gp_localisation",
    "gp_localisation_caveat",
    "hash_array",
    "hash_container",
    "hash_of",
    "kind_named",
    "model_fingerprint",
    "model_result_from_dict",
    "model_result_to_dict",
    "package_versions",
    "plot_anomaly_score",
    "plot_corner",
    "plot_gp_localisation",
    "plot_posterior_predictive",
    "plot_residuals",
    "plot_trace",
    "problem_fingerprint",
    "provenance_attrs",
    "register_kind",
    "registered_kinds",
    "spec_hashes",
    "to_netcdf",
    "training_pair_to_dict",
]
