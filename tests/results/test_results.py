"""W1.8: the results contract.

The acceptance criterion is :class:`TestNetcdfRoundTrip`: an emission from a real
``FittingProblem`` evaluation — ``inference.md`` §15's two-dataset joint fit with
a tied calibration nuisance — round-tripping through netCDF with the per-sample
``log_likelihood`` and ``log_prior``, the per-dataset contributions and the
provenance attrs all intact.

The problem is rebuilt here rather than imported from ``tests/core`` because
``tests/core`` has no conftest and no shared fixture module; the construction is
the one ``docs/design/contracts/inference.md`` §15 and
``tests/core/test_dataset.py::TestJointTwoDatasetProblem`` both use.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    ComplexGaussianFamily,
    Cube,
    Dataset,
    DenseGP,
    DatasetCollection,
    FittingProblem,
    FunctionSamples,
    GaussianFamily,
    GaussianProcessNoise,
    HierarchicalPrior,
    IndependentNoise,
    Image,
    Instrument,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    Plate,
    PlateBinding,
    PoissonFamily,
    Spectrum,
    Tie,
    TimeSeries,
    Transformation,
    VisibilitySet,
)
from ampere.core.dataset import Failure, FailureReason
from ampere.core.exceptions import SchemaError
from ampere.results import (
    ATTR_PREFIX,
    CONTAINER_SCHEMA_VERSION,
    GP_LOCALISATION_CAVEAT,
    LOG_LIKELIHOOD_DECOMPOSITION,
    POINTWISE_LOG_LIKELIHOOD_GROUP,
    DrawRecorder,
    ResultsError,
    add_pointwise_log_likelihood,
    add_posterior_predictive,
    add_residuals,
    canonical_json,
    container_from_dict,
    container_to_dict,
    describe_likelihood,
    emit,
    from_netcdf,
    gp_localisation,
    gp_localisation_caveat,
    hash_array,
    hash_container,
    hash_of,
    kind_named,
    model_fingerprint,
    model_identity_hash,
    model_result_from_dict,
    model_result_to_dict,
    neutral_model_identity,
    package_versions,
    plot_anomaly_score,
    plot_corner,
    plot_gp_localisation,
    plot_posterior_predictive,
    plot_residuals,
    plot_trace,
    problem_fingerprint,
    provenance_attrs,
    register_kind,
    spec_hashes,
    to_netcdf,
    training_pair_from_dict,
    training_pair_to_dict,
)
from ampere.results.emission import _dimension_names, _index_coordinate
from ampere.results.provenance import PROVENANCE_SCHEMA_VERSION, solver_configs

arviz = pytest.importorskip("arviz", reason="ampere.results needs ampere[arviz]")


# ---------------------------------------------------------------------------
# The problem under test: inference.md §15
# ---------------------------------------------------------------------------


class Powerlaw(Model):
    """A stub model with two channels, so one model can feed two datasets."""

    def __init__(self, **grids: np.ndarray) -> None:
        self._channels = tuple(grids)
        for name, grid in grids.items():
            self.register_buffer(name, np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            {
                name: Spectrum(ctx[name] * u.micron, ctx["norm"] * ctx[name] ** ctx["index"] * u.Jy)
                for name in self._channels
            }
        )


class Calibrate(Transformation):
    """A one-parameter nuisance step: the classic per-instrument scale factor."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.register_parameter(Parameter("scale", st.lognorm(0.2)))

    def apply(self, samples: Spectrum, values: Any) -> Spectrum:
        return samples.with_values(samples.values * self.context(values)["scale"])


BLUE = np.array([1.0, 2.0, 4.0])
RED = np.array([10.0, 20.0, 40.0])
TRUTH = {"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0}
#: A point the loguniform(0.1, 10) prior on ``norm`` rejects outright.
OUTSIDE = {"model.index": -1.0, "model.norm": 1e9, "calibration": 1.0}


def blue_data(values: Any = (1.0, 0.5, 0.25)) -> Spectrum:
    return Spectrum(BLUE * u.micron, list(values) * u.Jy, uncertainty=[0.05] * 3 * u.Jy)


def joint_problem(**kwargs: Any) -> FittingProblem:
    """``inference.md`` §15's worked example, verbatim."""
    red = Spectrum(RED * u.micron, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy)
    datasets = DatasetCollection(
        {
            "blue": Dataset(blue_data(), Instrument([Calibrate()], channel="blue")),
            "red": Dataset(red, Instrument([Calibrate()], channel="red")),
        }
    )
    tie = Tie(
        "calibration",
        ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"),
    )
    return FittingProblem(
        Powerlaw(blue=BLUE, red=RED), datasets, ties=[tie], seed=20260902, **kwargs
    )


def recorded(problem: FittingProblem, *, chains: int = 2, draws: int = 4) -> Any:
    """A small but genuine run: prior draws, plus one point the prior rejects."""
    recorder = DrawRecorder(problem, chains=chains)
    rng = np.random.default_rng(20260902)
    for chain in range(chains):
        for _ in range(draws - 1):
            recorder.record(problem.prior_transform(rng.random(problem.free_size)), chain=chain)
        recorder.record(OUTSIDE, chain=chain)
    return recorder.emit(engine="emcee")


# ---------------------------------------------------------------------------
# Hashing
# ---------------------------------------------------------------------------


class TestCanonicalJson:
    def test_mapping_order_is_irrelevant(self) -> None:
        assert canonical_json({"a": 1, "b": 2}) == canonical_json({"b": 2, "a": 1})

    def test_list_order_is_significant(self) -> None:
        # Load-bearing: lowering.md §9.2 makes a lowered model's seeding depend
        # on the order parameters are emitted in, so the merged spec's order is
        # part of what identifies a run.
        assert canonical_json([1, 2]) != canonical_json([2, 1])

    def test_no_whitespace_and_sorted(self) -> None:
        assert canonical_json({"b": [1, 2], "a": "x"}) == '{"a":"x","b":[1,2]}'

    @pytest.mark.parametrize(
        ("value", "expected"),
        [
            (float("inf"), '"__inf__"'),
            (float("-inf"), '"__-inf__"'),
            (float("nan"), '"__nan__"'),
        ],
    )
    def test_non_finite_floats_become_sentinels(self, value: float, expected: str) -> None:
        assert canonical_json(value) == expected

    @pytest.mark.parametrize(
        ("spelled", "value"),
        [
            ("__inf__", float("inf")),
            ("__-inf__", float("-inf")),
            ("__nan__", float("nan")),
        ],
    )
    def test_a_string_cannot_impersonate_a_float_sentinel(self, spelled: str, value: float) -> None:
        # The value-side twin of the reserved-key rule below: without the wrap,
        # normalise(float("nan")) and normalise("__nan__") would compare equal
        # and two genuinely different records would hash alike.
        assert canonical_json(spelled) != canonical_json(value)
        assert json.loads(canonical_json(spelled)) == {"__str__": spelled}

    def test_output_is_strict_json(self) -> None:
        # allow_nan=False stays on, so anything canonical_json emits parses with
        # a strict reader -- which a netCDF attribute's consumer is.
        assert json.loads(canonical_json({"x": float("nan"), "y": [1.0, 2.0]})) == {
            "x": "__nan__",
            "y": [1.0, 2.0],
        }

    def test_numpy_scalars_and_units_normalise(self) -> None:
        assert canonical_json({"n": np.int64(3), "u": u.Jy}) == '{"n":3,"u":{"__unit__":"Jy"}}'

    def test_an_unrecordable_object_is_refused(self) -> None:
        with pytest.raises(ResultsError, match="no defined JSON form"):
            canonical_json(object())

    @pytest.mark.parametrize("key", ["__ndarray__", "__unit__", "__quantity__", "__bytes__"])
    def test_a_mapping_cannot_impersonate_an_encoding(self, key: str) -> None:
        # Without this, a hand-written {"__ndarray__": {...}} would normalise to
        # the same tree as a real array, and two different objects would hash
        # alike -- the collision the refuse-don't-drop rule exists to prevent,
        # arriving by the other door.
        with pytest.raises(ResultsError, match="reserved"):
            canonical_json({key: {"anything": 1}})


class TestArrayHashing:
    def test_content_sensitive(self) -> None:
        assert hash_array(np.arange(3.0)) != hash_array(np.array([0.0, 1.0, 3.0]))

    def test_dtype_sensitive(self) -> None:
        assert hash_array(np.arange(3.0)) != hash_array(np.arange(3.0, dtype=np.float32))

    def test_shape_sensitive(self) -> None:
        assert hash_array(np.zeros((2, 3))) != hash_array(np.zeros((3, 2)))

    def test_non_contiguous_hashes_as_its_values(self) -> None:
        base = np.arange(6.0).reshape(2, 3)
        assert hash_array(base.T) == hash_array(np.ascontiguousarray(base.T))

    def test_string_arrays_hash(self) -> None:
        assert hash_array(np.array(["W1", "W2"])) != hash_array(np.array(["W1", "W3"]))

    def test_object_arrays_hash_by_value(self) -> None:
        assert hash_array(np.array([{"a": 1}], dtype=object)) == hash_array(
            np.array([{"a": 1}], dtype=object)
        )

    def test_digest_is_stable_across_processes(self) -> None:
        # hashlib, never hash(): PYTHONHASHSEED must not reach a provenance
        # record (ampere.core.rng makes the same argument for seeds).
        code = "from ampere.results import digest; print(digest('ampere'))"
        runs = {
            subprocess.run(
                [sys.executable, "-c", code],
                capture_output=True,
                text=True,
                check=True,
                env={"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"},
            ).stdout.strip()
            for seed in ("0", "1", "random")
        }
        assert len(runs) == 1
        assert len(runs.pop()) == 32


class TestContainerAndProblemFingerprints:
    def test_data_hash_moves_with_the_values(self) -> None:
        assert hash_container(blue_data()) != hash_container(blue_data((1.0, 0.5, 0.26)))

    def test_data_hash_moves_with_the_unit(self) -> None:
        millijansky = Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.mJy)
        jansky = Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy)
        assert hash_container(millijansky) != hash_container(jansky)

    def test_data_hash_ignores_meta(self) -> None:
        # A changed comment is not a different dataset.
        annotated = Spectrum(
            BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, meta={"note": "reduced 2026-09-02"}
        )
        plain = Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy)
        assert hash_container(annotated) == hash_container(plain)

    def test_problem_hash_is_stable_for_the_same_composition(self) -> None:
        assert hash_of(problem_fingerprint(joint_problem())) == hash_of(
            problem_fingerprint(joint_problem())
        )

    @staticmethod
    def _on_grid(grid: np.ndarray) -> FittingProblem:
        observed = Spectrum(grid * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy)
        return FittingProblem(
            Powerlaw(blue=grid),
            [Dataset(observed, Instrument([Calibrate()], channel="blue"), label="blue")],
            seed=20260902,
        )

    def test_problem_hash_moves_with_a_model_buffer(self) -> None:
        # A buffer is not a parameter, so it appears nowhere in to_spec() -- but
        # a model built on a different wavelength grid predicts different
        # numbers, and DEVELOPMENT_PLAN.md §7 wants this hash to invalidate a
        # trained artefact automatically. Same data, same priors, same seed.
        first, second = self._on_grid(BLUE), self._on_grid(np.array([1.0, 2.0, 4.000001]))
        assert (
            provenance_attrs(first)["ampere_spec_hash"]
            == provenance_attrs(second)["ampere_spec_hash"]
        )
        assert hash_of(problem_fingerprint(first)) != hash_of(problem_fingerprint(second))

    def test_problem_hash_moves_with_the_model_class(self) -> None:
        # Two models declaring identical parameters can compute entirely
        # different things; only the class distinguishes them.
        class Otherwise(Powerlaw):
            def evaluate(self, **values: Any) -> ModelResult:
                ctx = self.context(values)
                return ModelResult(
                    {"blue": Spectrum(ctx["blue"] * u.micron, ctx["norm"] * ctx["blue"] * u.Jy)}
                )

        observed = blue_data()
        common: dict[str, Any] = {"seed": 20260902}
        first = FittingProblem(
            Powerlaw(blue=BLUE),
            [Dataset(observed, Instrument([Calibrate()], channel="blue"), label="blue")],
            **common,
        )
        second = FittingProblem(
            Otherwise(blue=BLUE),
            [Dataset(observed, Instrument([Calibrate()], channel="blue"), label="blue")],
            **common,
        )
        assert (
            provenance_attrs(first)["ampere_spec_hash"]
            == provenance_attrs(second)["ampere_spec_hash"]
        )
        assert hash_of(problem_fingerprint(first)) != hash_of(problem_fingerprint(second))

    def test_problem_hash_moves_with_a_transformation_buffer(self) -> None:
        # The same argument for a response matrix or a filter curve in the chain.
        class Response(Transformation):
            ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

            def __init__(self, gains: np.ndarray, **kwargs: Any) -> None:
                super().__init__(**kwargs)
                self.register_buffer("gains", np.asarray(gains, dtype=float))

            def apply(self, samples: Spectrum, values: Any) -> Spectrum:
                return samples.with_values(samples.values * self.context(values)["gains"])

        def built(gains: np.ndarray) -> FittingProblem:
            return FittingProblem(
                Powerlaw(blue=BLUE),
                [
                    Dataset(
                        blue_data(),
                        Instrument([Response(gains)], channel="blue"),
                        label="blue",
                    )
                ],
                seed=20260902,
            )

        assert hash_of(problem_fingerprint(built(np.ones(3)))) != hash_of(
            problem_fingerprint(built(np.array([1.0, 1.0, 1.1])))
        )

    def test_problem_hash_moves_with_the_ties(self) -> None:
        red = Spectrum(RED * u.micron, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy)
        untied = FittingProblem(
            Powerlaw(blue=BLUE, red=RED),
            DatasetCollection(
                {
                    "blue": Dataset(blue_data(), Instrument([Calibrate()], channel="blue")),
                    "red": Dataset(red, Instrument([Calibrate()], channel="red")),
                }
            ),
            seed=20260902,
        )
        assert hash_of(problem_fingerprint(untied)) != hash_of(problem_fingerprint(joint_problem()))

    def test_a_fixed_parameter_s_value_is_part_of_the_spec(self) -> None:
        # A fixed parameter occupies no sampler dimension but does change the
        # answer, so the spec hash has to move with it. `to_spec` records both
        # `value` and `fixed`, so this comes for free -- asserted rather than
        # assumed, because it is the one class of "constant" that is a parameter.
        one = ParameterSet([Parameter("x", value=1.0, fixed=True)])
        two = ParameterSet([Parameter("x", value=2.0, fixed=True)])
        assert one.free_size == two.free_size == 0
        assert hash_of(one.to_spec()) != hash_of(two.to_spec())

    def test_likelihood_description_separates_kernels(self) -> None:
        # The concrete reason likelihoods.md §17 Q8 matters: two likelihoods
        # with identical ParameterSets but different kernels must not share a
        # provenance record, or a cached artefact is never invalidated.
        first = Likelihood(
            PoissonFamily(),
            GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
        )
        second = Likelihood(
            PoissonFamily(),
            GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
        )
        assert describe_likelihood(first) == describe_likelihood(second)
        assert describe_likelihood(first)["parameters"] == second.parameters.to_spec()
        assert describe_likelihood(first)["solver"]
        assert describe_likelihood(first)["kernel"]["family"] == "matern32"

    def test_a_family_or_noise_buffer_reaches_the_fingerprint(self) -> None:
        # Found by W1.13's serialisation review: Likelihood forwards
        # parameters from its pieces, never buffers, so a family's background
        # template was invisible to the provenance record — the stale-cache
        # trap DEVELOPMENT_PLAN.md §7 warns about.
        class Background(GaussianFamily):
            NAME = "test_buffered_background"

            def __init__(self, template: np.ndarray) -> None:
                self.register_buffer("template", template)

        with_one = Likelihood(Background(np.array([1.0, 2.0])), IndependentNoise())
        with_other = Likelihood(Background(np.array([1.0, 3.0])), IndependentNoise())
        first, second = describe_likelihood(with_one), describe_likelihood(with_other)
        assert first["parameters"] == second["parameters"]
        assert first["buffers"]["family"] != second["buffers"]["family"]
        assert hash_of(first) != hash_of(second)


class TestProvenanceAttrs:
    def test_every_value_is_netcdf_safe(self) -> None:
        for key, value in provenance_attrs(joint_problem()).items():
            assert key.startswith(ATTR_PREFIX)
            assert isinstance(value, int | float | str), f"{key} is {type(value).__name__}"
            # `bool` is an `int` in Python but is NOT a netCDF type, so the
            # obvious isinstance check above passes one straight through to a
            # write-time failure. Excluded explicitly.
            assert not isinstance(value, bool), f"{key} is a bool"

    def test_seed_and_its_source(self) -> None:
        attrs = provenance_attrs(joint_problem())
        assert attrs["ampere_seed"] == 20260902
        assert attrs["ampere_seed_source"] == "explicit"

    def test_an_unseeded_run_says_so_rather_than_inventing_a_seed(self) -> None:
        unseeded = FittingProblem(
            Powerlaw(blue=BLUE), [Dataset(blue_data(), Instrument([], channel="blue"))]
        )
        attrs = provenance_attrs(unseeded)
        assert attrs["ampere_seed_source"] == "entropy"
        assert "ampere_seed" not in attrs

    def test_spec_hash_moves_when_a_prior_moves(self) -> None:
        class Steeper(Model):
            def __init__(self, grid: np.ndarray, location: float) -> None:
                self.register_buffer("blue", grid, unit=u.micron)
                self.register_parameter(Parameter("index", st.norm(location, 0.5)))

            def evaluate(self, **values: Any) -> ModelResult:
                ctx = self.context(values)
                return ModelResult(
                    {"blue": Spectrum(ctx["blue"] * u.micron, ctx["blue"] ** ctx["index"] * u.Jy)}
                )

        def problem(location: float) -> FittingProblem:
            return FittingProblem(
                Steeper(BLUE, location),
                [Dataset(blue_data(), Instrument([], channel="blue", input_kind=Spectrum))],
                seed=20260902,
            )

        # Same seed, same data, one prior moved: a different run.
        assert (
            provenance_attrs(problem(-1.0))["ampere_spec_hash"]
            != provenance_attrs(problem(-2.0))["ampere_spec_hash"]
        )
        assert (
            provenance_attrs(problem(-1.0))["ampere_data_hash"]
            == provenance_attrs(problem(-2.0))["ampere_data_hash"]
        )

    def test_spec_hash_moves_when_a_component_is_relabelled(self) -> None:
        # The merged names are part of the spec, so a relabelled step is a
        # different declaration even though the priors are identical.
        red = Spectrum(RED * u.micron, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy)
        relabelled = FittingProblem(
            Powerlaw(blue=BLUE, red=RED),
            DatasetCollection(
                {
                    "blue": Dataset(blue_data(), Instrument([Calibrate()], channel="blue")),
                    "red": Dataset(red, Instrument([Calibrate(label="other")], channel="red")),
                }
            ),
            seed=20260902,
        )
        assert (
            provenance_attrs(relabelled)["ampere_spec_hash"]
            != provenance_attrs(joint_problem())["ampere_spec_hash"]
        )

    def test_data_hash_moves_when_the_data_move(self) -> None:
        red = Spectrum(RED * u.micron, [0.1, 0.05, 0.025] * u.Jy, uncertainty=[0.005] * 3 * u.Jy)
        moved = FittingProblem(
            Powerlaw(blue=BLUE, red=RED),
            DatasetCollection(
                {
                    "blue": Dataset(
                        blue_data((1.0, 0.5, 0.26)), Instrument([Calibrate()], channel="blue")
                    ),
                    "red": Dataset(red, Instrument([Calibrate()], channel="red")),
                }
            ),
            ties=[
                Tie(
                    "calibration",
                    ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"),
                )
            ],
            seed=20260902,
        )
        baseline = provenance_attrs(joint_problem())
        assert provenance_attrs(moved)["ampere_data_hash"] != baseline["ampere_data_hash"]
        assert provenance_attrs(moved)["ampere_spec_hash"] == baseline["ampere_spec_hash"]

    def test_composed_bindings_are_recorded(self) -> None:
        sites = json.loads(provenance_attrs(joint_problem())["ampere_sites"])
        assert sites["calibration"] == [
            "blue.instrument.calibrate.scale",
            "red.instrument.calibrate.scale",
        ]

    def test_free_names_not_free_labels(self) -> None:
        # likelihoods.md §16(a): 10^5 scalar names is the wrong representation,
        # so the attrs carry one name per parameter and the size, never a label
        # per flat-vector element.
        attrs = provenance_attrs(joint_problem())
        assert json.loads(attrs["ampere_free_names"]) == [
            "model.index",
            "model.norm",
            "calibration",
        ]
        assert attrs["ampere_free_size"] == 3

    def test_failure_counts_and_history_travel(self) -> None:
        class Crash(RuntimeError):
            pass

        class Wrapped(Model):
            def __init__(self, grid: np.ndarray) -> None:
                self.register_buffer("grid", grid, unit=u.micron)
                self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))

            def evaluate(self, **values: Any) -> ModelResult:
                raise Crash("the RT code exited 1")

        problem = FittingProblem(
            Wrapped(BLUE), [Dataset(blue_data())], validate=False, simulator_failures=(Crash,)
        )
        for _ in range(3):
            problem.log_prob({"model.slope": 2.0})
        attrs = provenance_attrs(problem)
        assert json.loads(attrs["ampere_failure_counts"]) == {"model_failed": 3}
        history = json.loads(attrs["ampere_failures"])
        assert len(history) == 3
        assert history[0]["exception_type"] == "Crash"
        assert history[0]["where"] == "model"

    def test_capabilities_are_recorded_as_declared(self) -> None:
        assert json.loads(provenance_attrs(joint_problem())["ampere_capabilities"]) == {
            "differentiable": False,
            "batchable": False,
            "device": "cpu",
            # W2.12's fourth flag. Adding this key is what bumped
            # PROVENANCE_SCHEMA_VERSION to 4.
            "backend": "reference",
        }

    def test_component_spec_hashes_locate_the_change(self) -> None:
        hashes = spec_hashes(joint_problem())
        assert set(hashes) == {"spec", "components"}
        assert set(hashes["components"]) == {"blue", "red", "model"}

    def test_a_component_labelled_spec_cannot_shadow_the_joint_hash(self) -> None:
        # "spec" is an ordinary word for a dataset, and a flat namespace would
        # let one silently replace the joint merged-spec hash -- leaving
        # ampere_spec_hash reporting one component's declaration instead of the
        # whole run's, which is what everything downstream relies on it for.
        def labelled_spec(location: float) -> FittingProblem:
            class Line(Model):
                def __init__(self, grid: np.ndarray) -> None:
                    self.register_buffer("blue", grid, unit=u.micron)
                    self.register_parameter(Parameter("index", st.norm(location, 0.5)))

                def evaluate(self, **values: Any) -> ModelResult:
                    ctx = self.context(values)
                    return ModelResult(
                        {"blue": Spectrum(ctx["blue"] * u.micron, ctx["blue"] * u.Jy)}
                    )

            return FittingProblem(
                Line(BLUE),
                [
                    Dataset(
                        blue_data(),
                        Instrument([], channel="blue", input_kind=Spectrum),
                        label="spec",
                    )
                ],
                seed=20260902,
            )

        problem = labelled_spec(-1.0)
        hashes = spec_hashes(problem)
        assert "spec" in hashes["components"]
        assert hashes["spec"] == hash_of(problem.parameters.to_spec())
        assert provenance_attrs(problem)["ampere_spec_hash"] == hashes["spec"]
        # and the load-bearing property still holds for such a problem
        assert (
            provenance_attrs(labelled_spec(-1.0))["ampere_spec_hash"]
            != provenance_attrs(labelled_spec(-2.0))["ampere_spec_hash"]
        )

    def test_problem_hash_moves_with_the_solver_configuration(self) -> None:
        # DenseGP's jitter is added to the diagonal before factorisation, so two
        # runs of the same strategy at different jitters score the same theta
        # differently. Recording only solver.NAME would let a cached artefact
        # trained under one regularisation be served for a fit under another.
        def with_jitter(jitter: float) -> Likelihood:
            return Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(Matern32(0.3, 1.0), solver=DenseGP(jitter=jitter)),
            )

        loose, tight = with_jitter(0.0), with_jitter(0.4)
        assert describe_likelihood(loose) != describe_likelihood(tight)
        assert describe_likelihood(loose)["solver"]["config"]["jitter"] == 0.0
        assert describe_likelihood(loose)["solver"]["name"]

    def test_versions_include_the_stack_that_moves_numbers(self) -> None:
        versions = package_versions()
        assert {"python", "numpy", "scipy", "astropy"} <= set(versions)

    def test_extra_attrs_cannot_shadow_ampere_s_own(self) -> None:
        with pytest.raises(ResultsError, match="would overwrite"):
            provenance_attrs(joint_problem(), extra={"seed": 1})

    def test_extra_attrs_are_prefixed_and_serialised(self) -> None:
        attrs = provenance_attrs(joint_problem(), extra={"walkers": 32, "moves": ["stretch"]})
        assert attrs["ampere_walkers"] == 32
        assert attrs["ampere_moves"] == '["stretch"]'

    @pytest.mark.parametrize(
        ("given", "expected"),
        [
            (True, 1),
            (False, 0),
            (np.bool_(True), 1),
            (np.int64(64), 64),
            (np.float64(0.25), 0.25),
            (0.25, 0.25),
            (float("nan"), '"__nan__"'),
        ],
        ids=["true", "false", "np-bool", "np-int", "np-float", "float", "nan"],
    )
    def test_extra_attrs_are_coerced_to_netcdf_types(self, given: Any, expected: Any) -> None:
        # A bool is an int in Python but not in netCDF, and both engines refuse
        # one; a numpy scalar was being stringified for no reason.
        value = provenance_attrs(joint_problem(), extra={"x": given})["ampere_x"]
        assert value == expected
        assert type(value) is type(expected)

    def test_a_boolean_extra_attr_survives_netcdf(self, tmp_path: Any) -> None:
        # The end-to-end version of the above: `adapt=True` is exactly what
        # extra_attrs is documented for, and it used to fail at write time.
        tree = recorded(joint_problem(), chains=1, draws=2)
        tree.attrs.update(provenance_attrs(joint_problem(), extra={"adapt": True}))
        path = tmp_path / "flags.nc"
        to_netcdf(tree, path)
        assert int(from_netcdf(path).attrs["ampere_adapt"]) == 1


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------


class TestEmission:
    def test_the_groups_a_run_carries(self) -> None:
        tree = recorded(joint_problem())
        assert set(tree.children) == {
            "posterior",
            "sample_stats",
            "log_likelihood",
            "observed_data",
            "constant_data",
        }

    def test_posterior_is_keyed_by_merged_name(self) -> None:
        tree = recorded(joint_problem())
        assert set(tree["posterior"].data_vars) == {"model.index", "model.norm", "calibration"}
        assert tree["posterior"]["calibration"].dims == ("chain", "draw")
        assert tree["posterior"]["calibration"].shape == (2, 4)

    def test_sample_stats_carries_the_split(self) -> None:
        tree = recorded(joint_problem())
        stats = tree["sample_stats"]
        assert {"lp", "log_prior", "log_likelihood"} <= set(stats.data_vars)
        finite = np.isfinite(stats["lp"].values)
        assert np.allclose(
            stats["lp"].values[finite],
            (stats["log_prior"].values + stats["log_likelihood"].values)[finite],
        )

    def test_log_likelihood_is_nan_for_a_prior_rejected_draw(self) -> None:
        # inference.md §18(c): "not evaluated" and "impossible" are different
        # statements, and importance reweighting needs them apart.
        tree = recorded(joint_problem())
        stats = tree["sample_stats"]
        rejected = np.isneginf(stats["log_prior"].values)
        assert rejected.sum() == 2  # one per chain, by construction
        assert np.isnan(stats["log_likelihood"].values[rejected]).all()
        assert not np.isneginf(stats["log_likelihood"].values[rejected]).any()

    def test_contributions_are_per_dataset_and_sum_to_the_joint(self) -> None:
        tree = recorded(joint_problem())
        group = tree["log_likelihood"]
        assert set(group.data_vars) == {"blue", "red"}
        total = group["blue"].values + group["red"].values
        joint = tree["sample_stats"]["log_likelihood"].values
        finite = np.isfinite(joint)
        assert np.allclose(total[finite], joint[finite])

    def test_a_prior_rejected_draw_contributes_nothing_rather_than_zero(self) -> None:
        tree = recorded(joint_problem())
        rejected = np.isneginf(tree["sample_stats"]["log_prior"].values)
        assert np.isnan(tree["log_likelihood"]["blue"].values[rejected]).all()

    def test_the_decomposition_is_declared_on_the_group(self) -> None:
        tree = recorded(joint_problem())
        assert tree["log_likelihood"].attrs[f"{ATTR_PREFIX}decomposition"] == (
            LOG_LIKELIHOOD_DECOMPOSITION
        )
        assert (
            "does not factorise" in tree["log_likelihood"].attrs[f"{ATTR_PREFIX}decomposition_note"]
        )

    def test_failures_are_recorded_per_draw(self) -> None:
        tree = recorded(joint_problem())
        stats = tree["sample_stats"]
        assert stats["failed"].values.dtype == bool
        # A prior rejection is not a failure: zero prior mass is an answer.
        assert not stats["failed"].values.any()
        assert set(np.unique(stats["failure_reason"].values)) == {""}

    def test_observed_data_carries_the_axis_as_a_coordinate(self) -> None:
        tree = recorded(joint_problem())
        observed = tree["observed_data"]
        assert set(observed.data_vars) == {"blue", "red"}
        assert observed["blue"].dims == ("blue_spectral_axis",)
        assert np.allclose(observed.coords["blue_spectral_axis"].values, BLUE)
        assert observed["blue"].attrs["units"] == "Jy"
        assert observed.coords["blue_spectral_axis"].attrs["units"] == "micron"

    def test_uncertainties_land_in_constant_data(self) -> None:
        tree = recorded(joint_problem())
        assert np.allclose(tree["constant_data"]["blue_uncertainty"].values, 0.05)

    def test_a_mask_is_stored_with_its_convention_stated(self) -> None:
        # netCDF has no boolean type, so the sense of the ones has to be said.
        masked = Spectrum(
            BLUE * u.micron,
            [1.0, 0.5, 0.25] * u.Jy,
            uncertainty=[0.05] * 3 * u.Jy,
            mask=[False, True, False],
        )
        problem = FittingProblem(
            Powerlaw(blue=BLUE),
            [Dataset(masked, Instrument([Calibrate()], channel="blue"), label="blue")],
        )
        recorder = DrawRecorder(problem)
        recorder.record(problem.reference_values)
        tree = recorder.emit()
        assert tree["constant_data"]["blue_mask"].values.tolist() == [0, 1, 0]
        assert "excluded" in tree["constant_data"].attrs[f"{ATTR_PREFIX}mask_convention"]
        assert f"{ATTR_PREFIX}mask_convention" not in tree["observed_data"].attrs

    def test_a_data_variable_name_collision_is_refused_not_overwritten(self) -> None:
        # Data-group names join the dataset label to a role name, so two labels
        # differing only by where an underscore falls can produce the same one.
        # Silently dropping one dataset's data would be far worse than refusing.
        problem = FittingProblem(
            Powerlaw(blue=BLUE),
            DatasetCollection(
                {
                    # dataset "a"'s uncertainties are stored as "a_uncertainty"
                    "a": Dataset(
                        blue_data(),
                        Instrument([], channel="blue", input_kind=Spectrum, label="a"),
                        label="a",
                    ),
                    # and so is dataset "a_uncertainty"'s own value array
                    "a_uncertainty": Dataset(
                        blue_data(),
                        Instrument([], channel="blue", input_kind=Spectrum, label="a_uncertainty"),
                        label="a_uncertainty",
                    ),
                }
            ),
        )
        recorder = DrawRecorder(problem)
        recorder.record()
        with pytest.raises(ResultsError, match="already claimed"):
            recorder.emit()

    def test_observed_groups_can_be_left_out(self) -> None:
        recorder = DrawRecorder(joint_problem())
        recorder.record(TRUTH)
        tree = recorder.emit(observed=False)
        assert "observed_data" not in tree.children

    def test_draw_shape_is_checked_against_the_problem(self) -> None:
        problem = joint_problem()
        with pytest.raises(ResultsError, match="free dimension"):
            emit(problem, np.zeros((2, 5)), [[problem.evaluate(TRUTH)] * 2] * 1)

    @pytest.mark.parametrize("shape", [(0, 2), (1, 0)], ids=["no-chains", "no-draws"])
    def test_an_empty_run_is_refused(self, shape: tuple[int, int]) -> None:
        # Not merely tidiness: emitting it would produce zero-length sampling
        # dimensions that every consumer downstream then has to guard against.
        problem = joint_problem()
        with pytest.raises(ResultsError, match="at least one chain and one draw"):
            emit(problem, np.zeros((*shape, problem.free_size)), [])

    def test_a_single_draw_in_a_single_chain_is_a_run(self) -> None:
        problem = joint_problem()
        tree = emit(
            problem,
            np.array([problem.parameters.pack(TRUTH)]),
            [problem.evaluate(TRUTH)],
        )
        assert tree["posterior"]["calibration"].shape == (1, 1)
        assert tree.attrs["ampere_chains"] == 1 and tree.attrs["ampere_draws"] == 1

    def test_one_evaluation_per_draw_is_required(self) -> None:
        problem = joint_problem()
        with pytest.raises(ResultsError, match="One Evaluation per draw"):
            emit(problem, np.zeros((3, 3)), [problem.evaluate(TRUTH)])

    def test_log_prob_is_not_enough(self) -> None:
        problem = joint_problem()
        with pytest.raises(ResultsError, match="cannot carry the split"):
            emit(problem, np.zeros((1, 3)), [problem.log_prob(TRUTH)])

    def test_provenance_reaches_the_root(self) -> None:
        tree = recorded(joint_problem())
        assert tree.attrs["ampere_seed"] == 20260902
        assert tree.attrs["ampere_chains"] == 2
        assert tree.attrs["ampere_draws"] == 4
        assert tree.attrs["ampere_engine"] == "emcee"


class TestArrayValuedParameters:
    """likelihoods.md §16(a) and hierarchical_population.md §10.2/§10.5."""

    @staticmethod
    def latent_problem() -> FittingProblem:
        counts = Spectrum([1.0, 2.0, 3.0] * u.um, np.array([4.0, 7.0, 2.0]))

        class Rate(Model):
            def __init__(self, grid: np.ndarray) -> None:
                self.register_buffer("grid", grid, unit=u.um)
                self.register_parameter(Parameter("rate", st.loguniform(0.5, 50.0)))

            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(ctx["grid"] * u.um, np.full(ctx["grid"].shape, ctx["rate"]))

        return FittingProblem(
            Rate(np.array([1.0, 2.0, 3.0])),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(
                        PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))
                    ),
                    label="counts",
                )
            ],
        )

    def test_a_latent_block_is_one_variable_with_a_dimension(self) -> None:
        problem = self.latent_problem()
        draws = np.zeros((1, 2, problem.free_size))
        draws[..., 0] = 3.0
        tree = emit(problem, draws, [[problem.evaluate(draws[0, d]) for d in range(2)]])
        block = tree["posterior"]["counts.latent.z"]
        assert block.dims == ("chain", "draw", "counts.latent.z_dim_0")
        assert block.shape == (1, 2, 3)
        # Never 10^5 scalar names: the whole block is one variable.
        assert set(tree["posterior"].data_vars) == {"model.rate", "counts.latent.z"}

    def test_a_plate_takes_its_own_name_as_the_dimension(self) -> None:
        plate = Plate(
            "objects",
            size=2,
            hyperparameters=[Parameter("mu", st.norm(0.0, 5.0))],
            members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))],
        )
        mapping = ParameterSet.merge(
            {
                "population": ParameterSet([], plates=[plate]),
                "obj0": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
                "obj1": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
            },
            plate_bindings=[
                PlateBinding("population.objects.theta", "obj0", "theta", 0),
                PlateBinding("population.objects.theta", "obj1", "theta", 1),
            ],
        )
        member = mapping.merged["population.objects.theta"]
        assert _dimension_names(member) == ("objects",)

    def test_the_plate_coordinate_is_read_off_the_bindings(self) -> None:
        # hierarchical_population.md §10.2: "W1.8 must record the dataset labels
        # as the plate's coordinate, not an integer range."
        plate = Plate(
            "objects",
            size=2,
            hyperparameters=[Parameter("mu", st.norm(0.0, 5.0))],
            members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu"}))],
        )
        mapping = ParameterSet.merge(
            {
                "population": ParameterSet([], plates=[plate]),
                "ngc1": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
                "ngc2": ParameterSet([Parameter("cal", st.lognorm(0.1))]),
            },
            plate_bindings=[
                PlateBinding("population.objects.theta", "ngc1", "theta", 0),
                PlateBinding("population.objects.theta", "ngc2", "theta", 1),
            ],
        )
        member = mapping.merged["population.objects.theta"]
        assert _index_coordinate(mapping.bindings, member) == ["ngc1", "ngc2"]

    def test_an_unrouted_block_gets_no_invented_coordinate(self) -> None:
        problem = self.latent_problem()
        member = problem.parameters["counts.latent.z"]
        assert _index_coordinate(problem.mapping.bindings, member) is None

    def test_a_caller_may_supply_the_coordinate(self) -> None:
        problem = self.latent_problem()
        draws = np.zeros((1, 2, problem.free_size))
        draws[..., 0] = 3.0
        tree = emit(
            problem,
            draws,
            [[problem.evaluate(draws[0, d]) for d in range(2)]],
            coords={"counts.latent.z_dim_0": ["a", "b", "c"]},
        )
        assert tree["posterior"].coords["counts.latent.z_dim_0"].values.tolist() == [
            "a",
            "b",
            "c",
        ]


class TestObservedDataKinds:
    """results.md §4's dimension rule, on every container kind that has one.

    A gridded kind takes one dimension per axis; a point kind with a single axis
    takes that axis, so a spectrum plots against wavelength unaided; a point kind
    whose axes index samples *jointly* takes one sample dimension and puts its
    axes in ``constant_data``; and a complex kind is split, because netCDF has no
    complex type.
    """

    ARCSEC: ClassVar[Any] = [0.0, 1.0] * u.arcsec

    @staticmethod
    def problem_for(container: FunctionSamples, likelihood: Likelihood | None = None) -> Any:
        class Emitter(Model):
            def __init__(self, template: FunctionSamples) -> None:
                self._template = template
                self.register_parameter(Parameter("scale", st.lognorm(0.2)))

            def evaluate(self, **values: Any) -> ModelResult:
                scale = self.context(values)["scale"]
                template = self._template
                return ModelResult({"ch": template.with_values(template.values * scale)})

        return FittingProblem(
            Emitter(container),
            [
                Dataset(
                    container,
                    Instrument([], channel="ch", input_kind=type(container)),
                    likelihood,
                    label="d",
                )
            ],
        )

    def cases(self) -> dict[str, tuple[FunctionSamples, Likelihood | None]]:
        return {
            "image": (
                Image(
                    self.ARCSEC,
                    [0.0, 1.0, 2.0] * u.arcsec,
                    np.arange(6.0).reshape(2, 3) * u.Jy,
                    uncertainty=np.full((2, 3), 0.1) * u.Jy,
                ),
                None,
            ),
            "cube": (
                Cube(
                    self.ARCSEC,
                    self.ARCSEC,
                    [1.0, 2.0] * u.um,
                    np.arange(8.0).reshape(2, 2, 2) * u.Jy,
                    uncertainty=np.full((2, 2, 2), 0.1) * u.Jy,
                ),
                None,
            ),
            "visibilities": (
                VisibilitySet(
                    [1.0, 2.0], [3.0, 4.0], np.array([1 + 2j, 3 - 1j]), uncertainty=[0.1, 0.1]
                ),
                Likelihood(ComplexGaussianFamily()),
            ),
            "photometry": (
                PhotometricPoints(
                    ["W1", "W2"],
                    [3.4, 4.6] * u.um,
                    [1.0, 2.0] * u.Jy,
                    uncertainty=[0.1, 0.2] * u.Jy,
                ),
                None,
            ),
        }

    EXPECTED: ClassVar[dict[str, tuple[list[str], list[str], list[str]]]] = {
        # kind: observed variables, their dims, constant_data variables
        "image": (["d"], ["d_x", "d_y"], ["d_uncertainty"]),
        "cube": (["d"], ["d_x", "d_y", "d_spectral_axis"], ["d_uncertainty"]),
        "visibilities": (["d_imag", "d_real"], ["d_index"], ["d_u", "d_uncertainty", "d_v"]),
        "photometry": (["d"], ["d_spectral_axis"], ["d_filters", "d_uncertainty"]),
    }

    @pytest.mark.parametrize("kind", ["image", "cube", "visibilities", "photometry"])
    def test_emitted_shape_and_round_trip(self, kind: str, tmp_path: Any) -> None:
        container, likelihood = self.cases()[kind]
        problem = self.problem_for(container, likelihood)
        recorder = DrawRecorder(problem)
        recorder.record()
        tree = recorder.emit()

        variables, dims, constant = self.EXPECTED[kind]
        assert sorted(tree["observed_data"].data_vars) == variables
        for name in variables:
            assert list(tree["observed_data"][name].dims) == dims
        assert sorted(tree["constant_data"].data_vars) == constant

        path = tmp_path / f"{kind}.nc"
        to_netcdf(tree, path)
        back = from_netcdf(path)
        assert sorted(back.children) == sorted(tree.children)
        for name in variables:
            assert np.allclose(
                back["observed_data"][name].values, tree["observed_data"][name].values
            )

    def test_a_complex_container_is_split_losslessly(self) -> None:
        container, likelihood = self.cases()["visibilities"]
        problem = self.problem_for(container, likelihood)
        recorder = DrawRecorder(problem)
        recorder.record()
        observed = recorder.emit()["observed_data"]
        rebuilt = observed["d_real"].values + 1j * observed["d_imag"].values
        assert np.allclose(rebuilt, container.values)


class TestDrawRecorder:
    def test_it_evaluates_when_not_given_an_evaluation(self) -> None:
        recorder = DrawRecorder(joint_problem())
        evaluation = recorder.record(TRUTH)
        assert math.isfinite(evaluation.log_prob)
        assert len(recorder) == 1

    def test_the_stored_theta_is_the_one_that_was_scored(self) -> None:
        # The two halves of a draw must agree. `evaluate(None)` scores the
        # reference theta, so `record(None)` has to store the reference theta
        # and not the NaN that unpacking None would give.
        problem = joint_problem()
        recorder = DrawRecorder(problem)
        evaluation = recorder.record()
        tree = recorder.emit()
        assert math.isfinite(evaluation.log_prob)
        for name, value in problem.reference_values.items():
            assert float(tree["posterior"][name].values[0, 0]) == pytest.approx(float(value))
        assert float(tree["sample_stats"]["lp"].values[0, 0]) == pytest.approx(evaluation.log_prob)

    @pytest.mark.parametrize("form", ["mapping", "vector"], ids=["mapping", "vector"])
    def test_both_accepted_forms_store_the_same_draw(self, form: str) -> None:
        problem = joint_problem()
        recorder = DrawRecorder(problem)
        recorder.record(TRUTH if form == "mapping" else problem.parameters.pack(TRUTH))
        tree = recorder.emit()
        assert float(tree["posterior"]["calibration"].values[0, 0]) == pytest.approx(1.0)

    def test_ragged_chains_are_refused_rather_than_padded(self) -> None:
        recorder = DrawRecorder(joint_problem(), chains=2)
        recorder.record(TRUTH, chain=0)
        with pytest.raises(ResultsError, match="different numbers of draws"):
            recorder.emit()

    def test_an_empty_recorder_has_nothing_to_emit(self) -> None:
        with pytest.raises(ResultsError, match="nothing has been recorded"):
            DrawRecorder(joint_problem()).emit()

    def test_an_out_of_range_chain_is_refused(self) -> None:
        with pytest.raises(ResultsError, match="out of range"):
            DrawRecorder(joint_problem()).record(TRUTH, chain=3)


# ---------------------------------------------------------------------------
# The acceptance criterion
# ---------------------------------------------------------------------------


class TestNetcdfRoundTrip:
    """W1.8's acceptance criterion, on ``inference.md`` §15's joint problem."""

    @staticmethod
    def round_tripped(tmp_path: Any, tree: Any) -> Any:
        path = tmp_path / "run.nc"
        written = to_netcdf(tree, path)
        assert path.is_file() and written == str(path)
        return from_netcdf(path)

    def test_the_whole_run_survives_netcdf(self, tmp_path: Any) -> None:
        problem = joint_problem()
        tree = recorded(problem)
        back = self.round_tripped(tmp_path, tree)

        # -- the groups ---------------------------------------------------
        assert set(back.children) == set(tree.children)

        # -- per-sample log_prior, log_likelihood and their sum ------------
        for name in ("lp", "log_prior", "log_likelihood"):
            assert np.allclose(
                back["sample_stats"][name].values,
                tree["sample_stats"][name].values,
                equal_nan=True,
            )
        rejected = np.isneginf(back["sample_stats"]["log_prior"].values)
        assert rejected.any()
        assert np.isnan(back["sample_stats"]["log_likelihood"].values[rejected]).all()
        finite = np.isfinite(back["sample_stats"]["lp"].values)
        assert np.allclose(
            back["sample_stats"]["lp"].values[finite],
            (
                back["sample_stats"]["log_prior"].values
                + back["sample_stats"]["log_likelihood"].values
            )[finite],
        )

        # -- the per-dataset contributions --------------------------------
        assert set(back["log_likelihood"].data_vars) == {"blue", "red"}
        for label in ("blue", "red"):
            assert np.allclose(
                back["log_likelihood"][label].values,
                tree["log_likelihood"][label].values,
                equal_nan=True,
            )
        total = back["log_likelihood"]["blue"].values + back["log_likelihood"]["red"].values
        joint = back["sample_stats"]["log_likelihood"].values
        assert np.allclose(total[finite], joint[finite])
        assert back["log_likelihood"].attrs[f"{ATTR_PREFIX}decomposition"] == "per_dataset"

        # -- the posterior itself -----------------------------------------
        for name in ("model.index", "model.norm", "calibration"):
            assert np.allclose(back["posterior"][name].values, tree["posterior"][name].values)

        # -- the observed data --------------------------------------------
        assert np.allclose(back["observed_data"]["blue"].values, [1.0, 0.5, 0.25])
        assert np.allclose(back["observed_data"].coords["blue_spectral_axis"].values, BLUE)
        assert back["observed_data"]["blue"].attrs["units"] == "Jy"
        assert np.allclose(back["constant_data"]["blue_uncertainty"].values, 0.05)

        # -- the provenance attrs -----------------------------------------
        assert dict(back.attrs) == {
            key: (int(value) if isinstance(value, int) else value)
            for key, value in tree.attrs.items()
        }
        assert back.attrs["ampere_seed"] == 20260902
        assert back.attrs["ampere_seed_source"] == "explicit"
        assert back.attrs["ampere_engine"] == "emcee"
        assert back.attrs["ampere_backend"] == "reference"
        assert len(back.attrs["ampere_spec_hash"]) == 32
        assert len(back.attrs["ampere_data_hash"]) == 32
        assert len(back.attrs["ampere_problem_hash"]) == 32
        assert json.loads(back.attrs["ampere_sites"])["calibration"] == [
            "blue.instrument.calibrate.scale",
            "red.instrument.calibrate.scale",
        ]
        assert json.loads(back.attrs["ampere_dataset_labels"]) == ["blue", "red"]
        assert json.loads(back.attrs["ampere_tied_names"]) == ["calibration"]
        assert json.loads(back.attrs["ampere_capabilities"])["device"] == "cpu"
        assert "numpy" in json.loads(back.attrs["ampere_library_versions"])
        assert json.loads(back.attrs["ampere_likelihoods"])["blue"]["family"] == "gaussian"

        # -- and it still identifies the same composition ------------------
        assert back.attrs["ampere_spec_hash"] == provenance_attrs(problem)["ampere_spec_hash"]
        assert back.attrs["ampere_problem_hash"] == hash_of(problem_fingerprint(problem))

    @pytest.mark.parametrize("engine", ["h5netcdf", "netcdf4"])
    def test_both_netcdf_engines_write_it(self, tmp_path: Any, engine: str) -> None:
        pytest.importorskip(
            {"h5netcdf": "h5netcdf", "netcdf4": "netCDF4"}[engine],
            reason=f"the {engine} backend is not installed",
        )
        tree = recorded(joint_problem(), chains=1, draws=2)
        path = tmp_path / f"run-{engine}.nc"
        to_netcdf(tree, path, engine=engine)
        back = from_netcdf(path)
        assert back.attrs["ampere_spec_hash"] == tree.attrs["ampere_spec_hash"]

    def test_a_latent_block_survives_with_its_dimension(self, tmp_path: Any) -> None:
        problem = TestArrayValuedParameters.latent_problem()
        draws = np.zeros((1, 2, problem.free_size))
        draws[..., 0] = 3.0
        tree = emit(problem, draws, [[problem.evaluate(draws[0, d]) for d in range(2)]])
        back = self.round_tripped(tmp_path, tree)
        assert back["posterior"]["counts.latent.z"].dims == (
            "chain",
            "draw",
            "counts.latent.z_dim_0",
        )
        assert back["posterior"]["counts.latent.z"].shape == (1, 2, 3)


# ---------------------------------------------------------------------------
# Serialisation (results_schema.md §17 Q6)
# ---------------------------------------------------------------------------


class TestContainerSerialisation:
    @pytest.mark.parametrize(
        "container",
        [
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy),
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, mask=[False, True, False]),
            PhotometricPoints(
                ["W1", "W2"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy, uncertainty=[0.1, 0.2] * u.Jy
            ),
            TimeSeries([0.0, 1.0, 2.0] * u.day, [1.0, 2.0, 3.0] * u.mag),
            Image(
                [0.0, 1.0] * u.arcsec,
                [0.0, 1.0, 2.0] * u.arcsec,
                np.arange(6.0).reshape(2, 3),
            ),
            Cube(
                [0.0, 1.0] * u.arcsec,
                [0.0, 1.0] * u.arcsec,
                [1.0, 2.0] * u.um,
                np.arange(8.0).reshape(2, 2, 2),
            ),
            VisibilitySet([1.0, 2.0], [3.0, 4.0], [1 + 2j, 3 - 1j]),
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, fidelity="cheap"),
        ],
        ids=[
            "spectrum",
            "masked",
            "photometry",
            "timeseries",
            "image",
            "cube",
            "visibilities",
            "fidelity",
        ],
    )
    def test_round_trip_by_value(self, container: FunctionSamples) -> None:
        back = container_from_dict(container_to_dict(container))
        assert type(back) is type(container)
        assert back == container
        assert back.fidelity == container.fidelity

    def test_the_encoding_is_plain_data(self) -> None:
        encoded = container_to_dict(
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy)
        )
        assert json.loads(json.dumps(encoded)) == encoded
        assert encoded["kind"] == "Spectrum"
        assert encoded["unit"] == "Jy"
        assert encoded["coordinates"]["spectral_axis"]["unit"] == "micron"

    def test_a_tampered_axis_order_is_refused(self) -> None:
        encoded = container_to_dict(Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy))
        encoded["coordinates"]["spectral_axis"]["values"] = [3.0, 2.0, 1.0]
        with pytest.raises(Exception, match="increasing"):
            container_from_dict(encoded)

    def test_a_tampered_mask_is_refused_not_coerced(self) -> None:
        # A mask is strictly boolean in the base contract; coercing `[0, 2, 0]`
        # here would silently accept a record the container itself refuses.
        encoded = container_to_dict(
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, mask=[False, True, False])
        )
        encoded["mask"] = [0, 2, 0]
        with pytest.raises(Exception, match="boolean"):
            container_from_dict(encoded)

    def test_a_tampered_shape_is_refused(self) -> None:
        encoded = container_to_dict(
            Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy, uncertainty=[0.05] * 3 * u.Jy)
        )
        encoded["uncertainty"]["data"] = [0.05, 0.05]
        with pytest.raises(SchemaError, match="uncertainties"):
            container_from_dict(encoded)

    def test_a_subclass_only_invariant_is_a_known_gap(self) -> None:
        # results.md §13.12: reconstruction runs the *base* contract's checks,
        # not an invariant a subclass declares in its own __init__ -- today that
        # is exactly one thing, PhotometricPoints' unique filter names. Pinned
        # so the limitation is visible and its removal is a deliberate change.
        encoded = container_to_dict(
            PhotometricPoints(["W1", "W2"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy)
        )
        encoded["extra_coords"]["filters"]["data"] = ["W1", "W1"]
        rebuilt = container_from_dict(encoded)
        assert rebuilt.filters.tolist() == ["W1", "W1"]
        with pytest.raises(Exception, match="repeated filter names"):
            PhotometricPoints(["W1", "W1"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy)

    def test_an_unknown_kind_names_the_remedy(self) -> None:
        with pytest.raises(ResultsError, match="register_kind"):
            container_from_dict({"version": CONTAINER_SCHEMA_VERSION, "kind": "Nope"})

    def test_a_future_version_is_refused(self) -> None:
        encoded = container_to_dict(Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy))
        encoded["version"] = 99
        with pytest.raises(ResultsError, match="unsupported container record version"):
            container_from_dict(encoded)

    def test_out_of_tree_kinds_register(self) -> None:
        @register_kind
        class Polarisation(FunctionSamples):
            AXES = Spectrum.AXES

        assert kind_named("Polarisation") is Polarisation
        curve = Polarisation({"spectral_axis": BLUE * u.micron}, [0.1, 0.2, 0.3])
        assert container_from_dict(container_to_dict(curve)) == curve

    def test_unserialisable_metadata_is_refused_not_dropped(self) -> None:
        container = Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy, meta={"origin": object()})
        with pytest.raises(ResultsError, match="no plain form"):
            container_to_dict(container)

    def test_plain_metadata_travels(self) -> None:
        container = Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy, meta={"programme": "GO-1234"})
        back = container_from_dict(container_to_dict(container))
        assert back.meta["programme"] == "GO-1234"


class TestModelResultSerialisation:
    def test_channels_and_theta_round_trip(self) -> None:
        result = ModelResult(
            {
                "blue": Spectrum(BLUE * u.micron, [1.0, 0.5, 0.25] * u.Jy),
                "red": Spectrum(RED * u.micron, [0.1, 0.05, 0.025] * u.Jy),
            },
            parameters={"model.index": -1.0, "model.norm": 1.0},
        )
        back = model_result_from_dict(model_result_to_dict(result))
        assert back == result
        assert back.parameters == {"model.index": -1.0, "model.norm": 1.0}

    def test_a_training_pair_is_what_simulate_already_produces(self) -> None:
        problem = joint_problem()
        simulation = problem.simulate(TRUTH)
        pair = training_pair_to_dict(
            simulation.parameters, simulation.results["model"], failed=simulation.failed
        )
        assert pair["failed"] is False
        assert set(pair["theta"]) == {"model.index", "model.norm", "calibration"}
        rebuilt = model_result_from_dict(pair["result"])
        assert set(rebuilt) == {"blue", "red"}
        assert np.allclose(rebuilt["blue"].values, simulation.results["model"]["blue"].values)

    def test_a_training_pair_round_trips_by_value(self) -> None:
        # serialisation_review.md §4: "training_pair_from_dict completes the
        # round trip when the first consumer lands". W2.8 is that consumer.
        problem = joint_problem()
        simulation = problem.simulate(TRUTH, observe=True)
        pair = training_pair_to_dict(
            simulation.parameters,
            simulation.results["model"],
            failed=simulation.failed,
            observations=simulation.observations,
        )
        back = training_pair_from_dict(pair)
        assert back["theta"] == dict(simulation.parameters)
        assert back["result"] == simulation.results["model"]
        assert back["observations"]["blue"] == simulation.observations["blue"]

    def test_theta_keeps_its_dtype(self) -> None:
        # The first of serialisation_review.md §4's two named losses: theta
        # array values went through .tolist() and came back as Python floats.
        pair = training_pair_to_dict(
            {"population.objects.theta": np.arange(3, dtype=np.int32)},
            ModelResult(Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy)),
        )
        theta = training_pair_from_dict(pair)["theta"]["population.objects.theta"]
        assert theta.dtype == np.dtype("int32")
        assert theta.tolist() == [0, 1, 2]

    def test_a_failure_travels_with_its_detail(self) -> None:
        # The second loss: only the flag was carried, so an archived budget
        # said *that* a draw was rejected and never why.
        failure = Failure(FailureReason.MODEL_FAILED, "the RT code exited 1", where="blue")
        pair = training_pair_to_dict(
            {"model.index": -1.0},
            ModelResult(Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy)),
            failed=True,
            failure=failure,
        )
        back = training_pair_from_dict(pair)
        assert back["failed"] is True
        assert back["failure"]["reason"] == "model_failed"
        assert back["failure"]["message"] == "the RT code exited 1"
        assert back["failure"]["where"] == "blue"

    def test_a_dotted_channel_name_survives(self) -> None:
        result = ModelResult({"co.j3_2": Spectrum(BLUE * u.micron, [1.0] * 3 * u.Jy)})
        assert set(model_result_from_dict(model_result_to_dict(result))) == {"co.j3_2"}


# ---------------------------------------------------------------------------
# The declared surfaces, now that Phase 2 has filled them in
# ---------------------------------------------------------------------------


class TestPlottingSurface:
    @pytest.mark.parametrize(
        "call",
        [
            lambda: plot_corner(None),
            lambda: plot_trace(None),
            lambda: plot_posterior_predictive(None),
            lambda: plot_residuals(None),
            lambda: plot_gp_localisation(None),
            lambda: plot_anomaly_score(None),
        ],
        ids=["corner", "trace", "ppc", "residuals", "gp_localisation", "anomaly"],
    )
    def test_every_plot_is_implemented_and_refuses_a_non_run(self, call: Any) -> None:
        # W2.7 landed diagnostics.md's families B and C and W2.8 the three
        # general-purpose ones (tests/results/test_diagnostics.py and
        # test_plots.py hold them end to end); what is asserted here is that
        # results.md §8's surface is complete — every one of the six refuses
        # *as a results error* rather than as unimplemented surface, which is
        # what tells a caller the difference between "not written yet" and
        # "you passed the wrong thing".
        with pytest.raises(ResultsError):
            call()

    def test_none_of_them_raises_not_implemented_any_more(self) -> None:
        # results.md §8's closing paragraph said "every plotting function is a
        # declared signature that raises NotImplementedError"; amended at W2.8,
        # and pinned here so the sentence and the code cannot drift apart again.
        for plot in (
            plot_corner,
            plot_trace,
            plot_posterior_predictive,
            plot_residuals,
            plot_gp_localisation,
            plot_anomaly_score,
        ):
            with pytest.raises(ResultsError):
                plot(None)

    def test_the_gp_localisation_caveat_is_mandatory_and_reachable(self) -> None:
        # diagnostics.md §4.3: the caveat must reach a programmatic consumer,
        # not only someone looking at the figure.
        assert gp_localisation_caveat() == GP_LOCALISATION_CAVEAT
        assert "does not say why" in GP_LOCALISATION_CAVEAT

    def test_the_caveat_is_in_the_docstring_by_construction(self) -> None:
        # "a plotting-function requirement for W1.8, not a 'please remember to
        # mention this' note" (diagnostics.md §4.3).
        doc = " ".join((plot_gp_localisation.__doc__ or "").split())
        assert "localises where the model is deficient; it does not say why" in doc
        assert "GP_LOCALISATION_CAVEAT" in doc


class TestDerivedGroups:
    @pytest.mark.parametrize(
        "call",
        [
            lambda: add_residuals(None, joint_problem()),
            lambda: gp_localisation(None, joint_problem()),
            lambda: add_posterior_predictive(None, joint_problem()),
            lambda: add_pointwise_log_likelihood(None, joint_problem()),
        ],
        ids=["residuals", "gp_localisation", "posterior_predictive", "pointwise"],
    )
    def test_the_landed_derived_groups_refuse_a_run_they_cannot_check(self, call: Any) -> None:
        # W2.7 implemented both. Handed something that is not a run, they say
        # so — deriving a residual against the wrong problem is undetectable
        # downstream, so the provenance check comes before the arithmetic.
        with pytest.raises(ResultsError, match="provenance"):
            call()

    def test_they_are_absent_from_a_default_emission(self) -> None:
        # The cost policy, asserted: N_draws x N_obs is not paid unless asked.
        # results.md §6's "never by default" for the pointwise group is the
        # same rule and is asserted with the rest of them, because having a
        # computation available is exactly when the rule stops enforcing itself.
        tree = recorded(joint_problem(), chains=1, draws=2)
        assert "posterior_predictive" not in tree.children
        assert "residuals" not in tree.children
        assert POINTWISE_LOG_LIKELIHOOD_GROUP not in tree.children


class TestDependencyPolicy:
    def test_importing_ampere_results_does_not_import_arviz(self) -> None:
        # architecture.md §4 rule 2: ampere.results imports its optional
        # dependency lazily, so `import ampere` stays light and the
        # minimal-install CI job keeps passing.
        code = (
            "import sys, ampere.results; "
            "print(any(name == 'arviz' or name.startswith('arviz.') for name in sys.modules))"
        )
        result = subprocess.run(
            [sys.executable, "-c", code], capture_output=True, text=True, check=True
        )
        assert result.stdout.strip() == "False"


class TestModelIdentity:
    """``describe()`` in the fingerprint, and the derived neutral identity.

    Both ruled by Peter 2026-09-03 at the freeze's escalations and landed with
    W2.1 as one mechanism (``results.md`` §13.13 and §14).
    """

    class Redden(Model):
        """A model configured by a plain attribute — §13.13's worked example."""

        def __init__(self, law: str) -> None:
            self.law = law
            self.register_buffer("wavelength", BLUE, unit=u.micron)
            self.register_parameter(Parameter("av", st.halfnorm(0.0, 1.0)))

        def describe(self) -> dict[str, str]:
            return {"law": self.law}

        def evaluate(self, **values: Any) -> ModelResult:
            ctx = self.context(values)
            return ModelResult(Spectrum(ctx["wavelength"] * u.micron, ctx["wavelength"] * u.Jy))

    def test_describe_reaches_the_fingerprint(self) -> None:
        """The cache-key hole §13.13 records, closed."""
        assert model_fingerprint(self.Redden("ccm89"))["describe"] == {"law": "ccm89"}

    def test_two_configurations_no_longer_share_a_cache_key(self) -> None:
        """Before the hook these two were indistinguishable to the hash."""
        assert model_fingerprint(self.Redden("ccm89")) != model_fingerprint(self.Redden("f99"))

    def test_a_model_that_does_not_opt_in_declares_nothing(self) -> None:
        assert model_fingerprint(Powerlaw(blue=BLUE))["describe"] is None

    def test_the_neutral_identity_drops_class_and_module_and_nothing_else(self) -> None:
        model = self.Redden("ccm89")
        fingerprint = model_fingerprint(model)
        neutral = neutral_model_identity(model)
        assert set(fingerprint) - set(neutral) == {"class", "module"}
        assert all(neutral[key] == fingerprint[key] for key in neutral)

    def test_the_neutral_identity_still_sees_the_configuration(self) -> None:
        """Offering an emulator across backends must not ignore §13.13's gap."""
        assert model_identity_hash(self.Redden("ccm89")) != model_identity_hash(self.Redden("f99"))

    def test_two_implementations_of_one_declaration_share_the_neutral_identity(self) -> None:
        """ "Offer": the point of the derived identity.

        Two different classes computing the same declaration agree here, and
        disagree on the problem hash — which is what makes a cross-backend
        emulator offerable but never silently servable.
        """

        class OtherPowerlaw(Powerlaw):
            """A different class, the same declaration."""

        left, right = Powerlaw(blue=BLUE), OtherPowerlaw(blue=BLUE)
        assert model_identity_hash(left) == model_identity_hash(right)
        assert model_fingerprint(left) != model_fingerprint(right)

    def test_the_recorded_attribute_is_netcdf_safe(self) -> None:
        attrs = provenance_attrs(joint_problem())
        assert isinstance(attrs["ampere_model_identity_hashes"], str)
        assert json.loads(attrs["ampere_model_identity_hashes"]) == {
            "model": model_identity_hash(joint_problem().models["model"])
        }

    def test_the_schema_version_records_the_change(self) -> None:
        """Adding a fingerprint key changes every problem hash, so it rides a bump."""
        assert PROVENANCE_SCHEMA_VERSION >= 3


class TestTheBackendIsDerived:
    """W2.12: ``ampere_backend`` is a fact about the problem, not a claim.

    The backend joined ``differentiable``/``batchable``/``device`` as §4.5's
    fourth capability flag, so the problem's pieces declare it and this module
    reads it off them. ``backend=`` survives only as a cross-check.
    """

    def test_it_is_read_off_the_problem_by_default(self) -> None:
        attrs = provenance_attrs(joint_problem(), engine="emcee")
        assert attrs["ampere_backend"] == joint_problem().backend == "reference"

    def test_an_agreeing_explicit_value_is_accepted(self) -> None:
        attrs = provenance_attrs(joint_problem(), backend="reference")
        assert attrs["ampere_backend"] == "reference"

    def test_a_disagreeing_explicit_value_raises_rather_than_being_recorded(self) -> None:
        # Recording it would be recording a falsehood, and silently preferring
        # either side would hide a real configuration mistake.
        with pytest.raises(ResultsError) as excinfo:
            provenance_attrs(joint_problem(), backend="torch")
        message = str(excinfo.value)
        assert "'torch'" in message and "'reference'" in message

    def test_emit_refuses_the_same_disagreement(self) -> None:
        pytest.importorskip("arviz")
        problem = joint_problem()
        draws = np.zeros((1, 1, problem.free_size))
        evaluations = [[problem.evaluate(draws[0, 0])]]
        with pytest.raises(ResultsError, match="'jax'"):
            emit(problem, draws, evaluations, engine="emcee", backend="jax")

    def test_the_capabilities_payload_carries_the_fourth_flag(self) -> None:
        # This key is what bumped PROVENANCE_SCHEMA_VERSION to 4: the shape of
        # ampere_capabilities changed, so schema-3 artefacts are not
        # comparable. The constant has moved on to 5 (W2.13's three new
        # attributes) and this row asserts the payload, not the number.
        payload = json.loads(provenance_attrs(joint_problem())["ampere_capabilities"])
        assert payload["backend"] == "reference"
        assert PROVENANCE_SCHEMA_VERSION >= 4


def _gp_problem(solver: DenseGP) -> FittingProblem:
    """One GP dataset, so the solver has somewhere to report configuration from."""
    return FittingProblem(
        Powerlaw(blue=BLUE),
        DatasetCollection(
            {
                "gp": Dataset(
                    blue_data(),
                    Instrument([Calibrate()], channel="blue"),
                    Likelihood(
                        GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 1.0), solver=solver)
                    ),
                )
            }
        ),
        seed=20260907,
    )


class TestSchemaFiveAttributes:
    """W2.13: ``ampere_realised``, ``ampere_registered_lowerings``, ``ampere_solver_config``.

    All three are written on **every** run, including the gradient-free ones,
    because "this run's gradients were real" and "this run had none" are the
    two answers a reader must be able to tell apart, and silence distinguishes
    neither.
    """

    def test_the_schema_version_is_five(self) -> None:
        assert PROVENANCE_SCHEMA_VERSION == 5

    def test_a_contract_path_run_records_realised_zero(self) -> None:
        attrs = provenance_attrs(joint_problem(), engine="emcee")
        # An int, not a bool: netCDF has no boolean attribute type, and both
        # netCDF engines refuse one outright.
        assert attrs["ampere_realised"] == 0
        assert isinstance(attrs["ampere_realised"], int)

    def test_a_realised_run_says_so(self) -> None:
        attrs = provenance_attrs(joint_problem(), engine="nuts", realised=True)
        assert attrs["ampere_realised"] == 1

    def test_registered_lowerings_is_an_empty_list_by_default(self) -> None:
        attrs = provenance_attrs(joint_problem())
        assert json.loads(attrs["ampere_registered_lowerings"]) == []

    def test_registered_lowerings_records_what_it_is_given(self) -> None:
        rows = [
            {
                "kind": "prior",
                "name": "custom",
                "backend": "stub",
                "builtin": False,
                "constructor": "pkg.build",
            }
        ]
        attrs = provenance_attrs(joint_problem(), registered_lowerings=rows)
        assert json.loads(attrs["ampere_registered_lowerings"]) == rows

    def test_the_solver_config_is_empty_on_the_reference_path(self) -> None:
        problem = joint_problem()
        assert solver_configs(problem) == {}
        assert json.loads(provenance_attrs(problem)["ampere_solver_config"]) == {}

    def test_a_solver_with_configuration_is_recorded_and_not_hashed(self) -> None:
        """Fold-in 10: visible in the attrs, invisible to every hash."""

        class ConfiguredDenseGP(DenseGP):
            def provenance_config(self) -> dict[str, Any]:
                return {"dtype": "float64", "device": "cpu"}

        plain = _gp_problem(DenseGP())
        configured = _gp_problem(ConfiguredDenseGP())
        stored = json.loads(provenance_attrs(configured)["ampere_solver_config"])
        assert stored == {"gp": {"dtype": "float64", "device": "cpu"}}
        # The declaration is unchanged, so the spec hash is too -- which is
        # exactly what results.md §14 requires of two backends implementing one
        # declaration with different precision policies.
        assert (
            provenance_attrs(configured)["ampere_spec_hash"]
            == provenance_attrs(plain)["ampere_spec_hash"]
        )
        assert json.loads(
            provenance_attrs(configured)["ampere_component_spec_hashes"]
        ) == json.loads(provenance_attrs(plain)["ampere_component_spec_hashes"])
        assert "dtype" not in provenance_attrs(configured)["ampere_likelihoods"]

    def test_an_emitted_run_carries_all_three(self) -> None:
        run = recorded(joint_problem())
        for key in ("ampere_realised", "ampere_registered_lowerings", "ampere_solver_config"):
            assert key in run.attrs
