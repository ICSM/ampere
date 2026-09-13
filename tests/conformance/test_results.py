"""Conformance rows for the results contract — what a run leaves behind.

The inventory is ``results.md`` §14's hand-down, verbatim: "the per-dataset
contributions sum to the scalar joint ``log_likelihood``; ``lp == log_prior +
log_likelihood`` wherever finite; a prior-rejected draw stores ``-inf``/NaN and
no failure; a netCDF round trip preserves every group, every NaN and every
provenance attribute; the spec hash changes when a prior changes and not when
the data change, and the data hash vice versa; the digest is identical across
processes with different ``PYTHONHASHSEED``; a container of each kind
round-trips by value; an array-valued parameter emits one variable with a named
dimension."

The **backend-spanning** row that section adds — two backends emitting the same
problem produce the same hashes — lives in ``test_cross_backend.py``, with the
other rows that compare two implementations rather than one against an oracle.

``ampere.results`` is arviz-gated (``architecture.md`` §4 rule 2 keeps arviz an
extra), so the whole module skips where arviz is absent — which is the case in
the ``test-py311``/``test-py312``/``test-py313`` CI environments by design.
"""

from __future__ import annotations

import dataclasses
import json
import os
import subprocess
import sys
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Cube,
    FittingProblem,
    FunctionSamples,
    Image,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    Tie,
    VisibilitySet,
)

from .composition import (
    GP_GRID,
    DatasetSpec,
    ProblemSpec,
    build_problem,
)
from .protocol import (
    ConformanceBackend,
    ModelKind,
    ModelSpec,
    TransformationKind,
    TransformationSpec,
)

arviz = pytest.importorskip("arviz", reason="the results contract needs ampere[arviz]")

from ampere.results import (  # noqa: E402  (must follow the arviz gate)
    LOG_LIKELIHOOD_GROUP,
    POSTERIOR_GROUP,
    SAMPLE_STATS_GROUP,
    DrawRecorder,
    container_from_dict,
    container_to_dict,
    emit,
    from_netcdf,
    model_result_from_dict,
    model_result_to_dict,
    provenance_attrs,
    to_netcdf,
)

CALIBRATION = TransformationSpec(TransformationKind.SCALE, label="calibration")
CALIBRATION_SITES = (
    "blue.instrument.calibration.scale",
    "red.instrument.calibration.scale",
)

JOINT = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, channels=("blue", "red"), coordinates=GP_GRID),
    datasets=(
        DatasetSpec(label="blue", channel="blue", instrument=(CALIBRATION,)),
        DatasetSpec(
            label="red", channel="red", instrument=(CALIBRATION,), data_seed=771, masked=(2,)
        ),
    ),
)

PLATED = ProblemSpec(
    model=ModelSpec(channels=("blue", "red"), coordinates=GP_GRID, plated=True),
    datasets=(
        DatasetSpec(label="blue", channel="blue"),
        DatasetSpec(label="red", channel="red", data_seed=771),
    ),
)


@pytest.fixture
def joint(backend: ConformanceBackend) -> FittingProblem:
    """Two datasets, one of them masked, so the emitted groups have content."""
    return build_problem(backend, JOINT)


def recorded(problem: FittingProblem, draws: int = 6, chains: int = 2) -> object:
    """A small run: mostly prior draws, with one the prior rejects per chain."""
    recorder = DrawRecorder(problem, chains=chains)
    rng = np.random.default_rng(20260902)
    rejected = dict(problem.reference_values)
    rejected[problem.parameters.free_names[0]] = 1.0e6
    for chain in range(chains):
        for _ in range(draws - 1):
            recorder.record(problem.prior_transform(rng.random(problem.free_size)), chain=chain)
        recorder.record(rejected, chain=chain)
    return recorder.emit(engine="conformance")


class TestEmittedDecomposition:
    """What ``sample_stats`` and the ``log_likelihood`` group must agree about."""

    def test_the_per_dataset_contributions_sum_to_the_scalar(self, joint: FittingProblem) -> None:
        tree = recorded(joint)
        stats = tree[SAMPLE_STATS_GROUP]
        per_dataset = tree[LOG_LIKELIHOOD_GROUP]
        total = sum(per_dataset[label].values for label in ("blue", "red"))
        finite = np.isfinite(stats["log_likelihood"].values)
        assert total[finite] == pytest.approx(stats["log_likelihood"].values[finite])

    def test_lp_is_the_sum_of_its_halves_wherever_finite(self, joint: FittingProblem) -> None:
        stats = recorded(joint)[SAMPLE_STATS_GROUP]
        lp = stats["lp"].values
        prior = stats["log_prior"].values
        likelihood = stats["log_likelihood"].values
        finite = np.isfinite(lp)
        assert finite.any()
        assert lp[finite] == pytest.approx((prior + likelihood)[finite])

    def test_a_prior_rejected_draw_is_minus_infinity_and_nan_and_not_a_failure(
        self, joint: FittingProblem
    ) -> None:
        """The distinction design horizon (b)'s reweighting needs, kept intact.

        ``-inf`` prior means "impossible"; NaN likelihood means "never
        evaluated". Coercing the second to ``-inf`` would lose the difference,
        and a rejected draw is not a failure — nothing went wrong.
        """
        stats = recorded(joint)[SAMPLE_STATS_GROUP]
        rejected = np.isneginf(stats["log_prior"].values)
        assert rejected.any()
        assert np.isnan(stats["log_likelihood"].values[rejected]).all()
        assert not np.isneginf(stats["log_likelihood"].values[rejected]).any()
        assert not stats["failed"].values[rejected].any()
        for label in ("blue", "red"):
            assert np.isnan(recorded(joint)[LOG_LIKELIHOOD_GROUP][label].values[rejected]).all()


class TestArrayValuedParameters:
    """A plate is one variable with a named dimension, not N scalars."""

    def test_it_emits_one_variable_carrying_the_plates_name(
        self, backend: ConformanceBackend
    ) -> None:
        problem = build_problem(backend, PLATED)
        tree = emit(
            problem,
            np.array([problem.parameters.pack(problem.reference_values)]),
            [problem.evaluate()],
            coords={"objects": ["blue", "red"]},
        )
        posterior = tree[POSTERIOR_GROUP]
        variable = posterior["model.objects.offsets"]

        assert variable.dims == ("chain", "draw", "objects")
        assert variable.shape == (1, 1, 2)
        assert "model.objects.offsets[0]" not in posterior
        assert list(posterior.coords["objects"].values) == ["blue", "red"]

    def test_the_plate_is_recorded_in_the_provenance(self, backend: ConformanceBackend) -> None:
        problem = build_problem(backend, PLATED)
        assert problem.parameters.plates == {"objects": 2}
        assert json.loads(provenance_attrs(problem)["ampere_plates"]) == {"objects": 2}


class TestNetCDFRoundTrip:
    """The stored run is the run."""

    def test_every_group_nan_and_attribute_survives(
        self, joint: FittingProblem, tmp_path: Path
    ) -> None:
        pytest.importorskip("h5netcdf", reason="netCDF serialisation needs an engine")
        tree = recorded(joint)
        path = tmp_path / "run.nc"
        assert to_netcdf(tree, path) == str(path)
        back = from_netcdf(path)

        assert set(back.children) == set(tree.children)
        for group in (POSTERIOR_GROUP, SAMPLE_STATS_GROUP, LOG_LIKELIHOOD_GROUP):
            for name, variable in tree[group].items():
                restored = back[group][name].values
                if variable.dtype.kind == "f":
                    assert np.array_equal(restored, variable.values, equal_nan=True)
                else:
                    assert np.array_equal(restored, variable.values)

        for key, value in tree.attrs.items():
            assert back.attrs[key] == value

    def test_the_nan_convention_survives_specifically(
        self, joint: FittingProblem, tmp_path: Path
    ) -> None:
        pytest.importorskip("h5netcdf", reason="netCDF serialisation needs an engine")
        tree = recorded(joint)
        path = tmp_path / "run.nc"
        to_netcdf(tree, path)
        back = from_netcdf(path)

        rejected = np.isneginf(tree[SAMPLE_STATS_GROUP]["log_prior"].values)
        assert np.isneginf(back[SAMPLE_STATS_GROUP]["log_prior"].values[rejected]).all()
        assert np.isnan(back[SAMPLE_STATS_GROUP]["log_likelihood"].values[rejected]).all()


class TestHashSensitivity:
    """The spec hash follows the declaration; the data hash follows the data."""

    def test_a_changed_prior_moves_the_spec_hash_and_not_the_data_hash(
        self, backend: ConformanceBackend
    ) -> None:
        """A tie's prior override is a declaration change and nothing else."""
        base = build_problem(
            backend,
            ProblemSpec(JOINT.model, JOINT.datasets, ties=(Tie("calibration", CALIBRATION_SITES),)),
        )
        overridden = build_problem(
            backend,
            ProblemSpec(
                JOINT.model,
                JOINT.datasets,
                ties=(Tie("calibration", CALIBRATION_SITES, prior=st.lognorm(0.5)),),
            ),
        )
        first, second = provenance_attrs(base), provenance_attrs(overridden)
        assert first["ampere_spec_hash"] != second["ampere_spec_hash"]
        assert first["ampere_data_hash"] == second["ampere_data_hash"]

    def test_changed_data_move_the_data_hash_and_not_the_spec_hash(
        self, backend: ConformanceBackend
    ) -> None:
        moved = tuple(
            dataclasses.replace(dataset, data_seed=dataset.data_seed + 1)
            for dataset in JOINT.datasets
        )
        first = provenance_attrs(build_problem(backend, JOINT))
        second = provenance_attrs(build_problem(backend, ProblemSpec(JOINT.model, moved)))
        assert first["ampere_data_hash"] != second["ampere_data_hash"]
        assert first["ampere_spec_hash"] == second["ampere_spec_hash"]

    def test_the_digest_is_identical_across_hash_seeds(self) -> None:
        """The hashes are BLAKE2b over canonical JSON, never Python's ``hash()``.

        A dict-ordering or ``hash()``-derived digest passes every in-process
        check and silently invalidates every cache key the moment
        ``PYTHONHASHSEED`` changes — which it does on every interpreter start.
        """
        script = (
            "from ampere.results import digest, hash_of;"
            "print(digest('ampere'), hash_of({'b': [1, 2], 'a': {'x': 1.5}}))"
        )
        outputs = set()
        for hash_seed in ("0", "1", "random"):
            result = subprocess.run(
                [sys.executable, "-c", script],
                capture_output=True,
                text=True,
                check=True,
                env=dict(os.environ, PYTHONHASHSEED=hash_seed),
            )
            outputs.add(result.stdout.strip())
        assert len(outputs) == 1


def containers() -> dict[str, FunctionSamples]:
    """One minimally-populated container of every registered kind."""
    grid = np.asarray(GP_GRID[:4], dtype=float)
    return {
        "spectrum": Spectrum(
            grid * u.micron,
            np.linspace(1.0, 2.0, 4) * u.Jy,
            uncertainty=np.full(4, 0.1) * u.Jy,
            mask=np.array([False, True, False, False]),
        ),
        "photometry": PhotometricPoints(
            ("W1", "W2", "W3", "W4"),
            grid * u.micron,
            np.linspace(1.0, 2.0, 4) * u.Jy,
            uncertainty=np.full(4, 0.1) * u.Jy,
            extra_coords={"limit_kind": np.array([0, 0, 1, 0])},
        ),
        "timeseries": TimeSeries(np.arange(4.0) * u.day, np.linspace(10.0, 11.0, 4) * u.mag),
        "image": Image(
            np.arange(2.0) * u.arcsec, np.arange(3.0) * u.arcsec, np.arange(6.0).reshape(2, 3)
        ),
        "cube": Cube(
            np.arange(2.0) * u.arcsec,
            np.arange(2.0) * u.arcsec,
            np.array([1.0, 2.0]) * u.micron,
            np.arange(8.0).reshape(2, 2, 2),
        ),
        "visibilities": VisibilitySet(
            np.array([1.0, 2.0]),
            np.array([3.0, 4.0]),
            np.array([2.2, 2.2]) * u.um,
            np.array([1 + 2j, 3 - 1j]),
        ),
    }


class TestContainerRoundTrips:
    """Every container kind survives the plain-data form by value."""

    @pytest.mark.parametrize("kind", sorted(containers()))
    def test_a_container_round_trips_by_value(self, kind: str) -> None:
        container = containers()[kind]
        restored = container_from_dict(container_to_dict(container))
        assert type(restored) is type(container)
        assert restored == container

    def test_a_model_result_round_trips_with_its_parameters(
        self, backend: ConformanceBackend
    ) -> None:
        problem = build_problem(backend, JOINT)
        simulation = problem.simulate(problem.reference_values)
        result = simulation.results["model"]
        assert result.parameters is not None
        restored = model_result_from_dict(model_result_to_dict(result))
        assert restored == result
