"""W2.8: ``results.md`` §11 layer 2's training-set writer, end to end.

The acceptance criterion is "a training set written, appended to, and read back
round-trips by value", and the tests below take "by value" literally: the
containers that come back must compare equal to the ones
:meth:`~ampere.core.dataset.FittingProblem.simulate` produced, coordinates,
units, masks, uncertainties and all — which is the whole point of the format
being netCDF rather than a pickle of whatever was in memory.

The other three properties the format claims are asserted too, because each is
a promise the *writer* has to keep rather than a fact about netCDF: NaN is the
native "not written" value; the coordinate arrays are stored once for the whole
set (so a set whose grid moves per sample is refused); and the spec hash sits in
the attributes, where :func:`~ampere.results.append_training_set` checks it, so
``DEVELOPMENT_PLAN.md`` §7's stale-artefact trap is a string comparison.

Everything is written into pytest's ``tmp_path``: no netCDF file belongs in this
repository (``AGENTS.md`` ground rule 7).
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    Derived,
    FittingProblem,
    Model,
    ModelResult,
    Parameter,
    PhotometricPoints,
    Spectrum,
)
from ampere.core.dataset import Failure, FailureReason, Simulation
from ampere.core.exceptions import ResultsError
from ampere.core.simulate import ScaledSigma
from ampere.results import (
    ATTR_PREFIX,
    PROVENANCE_SCHEMA_VERSION,
    TRAINING_SET_SCHEMA_VERSION,
    append_training_set,
    provenance_attrs,
    read_training_set,
    write_training_set,
)

pytest.importorskip("xarray", reason="the training-set writer needs xarray")

GRID = np.linspace(1.0, 4.0, 6)


class Powerlaw(Model):
    """``norm * x ** index``, with a per-sample fidelity tag on the output.

    The tag is there so the round trip has something beyond values to lose:
    ``fidelity`` is a container field that nothing in the arithmetic reads, and
    it is exactly the sort of thing a writer drops without anybody noticing.
    """

    def __init__(self, grid: np.ndarray) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.3)))
        self.register_parameter(Parameter("norm", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            Spectrum(
                ctx["grid"] * u.micron,
                ctx["norm"] * ctx["grid"] ** ctx["index"] * u.Jy,
                fidelity="cheap",
            )
        )


class Steeper(Powerlaw):
    """The same model with a different prior — a different declaration."""

    def __init__(self, grid: np.ndarray) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-2.0, 0.3)))
        self.register_parameter(Parameter("norm", st.loguniform(0.5, 2.0)))


def toy(model: type[Powerlaw] = Powerlaw, *, mask: np.ndarray | None = None) -> FittingProblem:
    observed = Spectrum(
        GRID * u.micron,
        (1.0 * GRID**-1.0) * u.Jy,
        uncertainty=np.full(GRID.size, 0.02) * u.Jy,
        mask=mask,
    )
    return FittingProblem(
        model(GRID),
        DatasetCollection({"sed": Dataset(observed, label="sed")}),
        seed=20260908,
        # A user's simulator raises its own exception class, and declaring it is
        # how `simulate` flags rather than raises (inference.md §13).
        simulator_failures=(RuntimeError,),
    )


def budget(problem: FittingProblem, size: int, *, observe: bool = True) -> list[Simulation]:
    return [problem.simulate(observe=observe) for _ in range(size)]


# ---------------------------------------------------------------------------
# The round trip
# ---------------------------------------------------------------------------


class TestRoundTrip:
    def test_a_written_set_reads_back_by_value(self, tmp_path: Path) -> None:
        problem = toy()
        drawn = budget(problem, 4)
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert len(stored) == 4
        for index, simulation in enumerate(drawn):
            result = stored.result(index)
            original = simulation.results["model"]
            for channel in original:
                assert result[channel] == original[channel], channel
            assert result[next(iter(original))].fidelity == "cheap"

    def test_theta_comes_back_by_merged_name(self, tmp_path: Path) -> None:
        problem = toy()
        drawn = budget(problem, 3)
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert set(stored.theta) == {"model.index", "model.norm"}
        assert np.allclose(
            stored.theta["model.index"], [s.parameters["model.index"] for s in drawn]
        )
        assert stored.parameters(0)["model.norm"] == pytest.approx(
            drawn[0].parameters["model.norm"]
        )

    def test_observations_are_carried(self, tmp_path: Path) -> None:
        # serialisation_review.md §4's second named loss: Simulation.observations
        # were not carried at all.
        problem = toy()
        drawn = budget(problem, 3, observe=True)
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.observed == ("sed",)
        assert stored.observations(1)["sed"] == drawn[1].observations["sed"]

    def test_a_budget_without_observations_carries_none(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2, observe=False), problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.observed == ()
        assert stored.observations(0) == {}

    def test_a_mask_and_its_uncertainties_survive(self, tmp_path: Path) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[1:3] = True
        problem = toy(mask=mask)
        drawn = budget(problem, 2)
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        observations = stored.observations(0)["sed"]
        assert observations.mask is not None
        assert observations.mask.tolist() == mask.tolist()
        assert observations.uncertainty is not None
        assert observations == drawn[0].observations["sed"]

    def test_the_pair_form_matches_the_in_memory_one(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        pair = read_training_set(tmp_path / "budget.nc").pair(0)
        assert sorted(pair) == ["failed", "observations", "result", "theta"]

    def test_the_provenance_recipe_is_section_nine_s(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert len(stored.spec_hash) == 32
        assert len(stored.attrs[f"{ATTR_PREFIX}problem_hash"]) == 32
        assert len(stored.attrs[f"{ATTR_PREFIX}data_hash"]) == 32
        # W3.12: ampere_model_hash joined the recipe at schema 6; schema 7
        # (W5.0), schema 8 (W5.22) and schema 9 (W6.7, ampere_start_route)
        # add no root attribute this function's own recipe touches.
        assert len(stored.attrs[f"{ATTR_PREFIX}model_hash"]) == 32
        assert stored.attrs[f"{ATTR_PREFIX}schema_version"] == PROVENANCE_SCHEMA_VERSION == 10
        assert stored.attrs[f"{ATTR_PREFIX}seed"] == 20260908
        assert stored.attrs[f"{ATTR_PREFIX}training_set_version"] == TRAINING_SET_SCHEMA_VERSION


# ---------------------------------------------------------------------------
# Extra coordinates keep their dtype (W6.14)
# ---------------------------------------------------------------------------

FILTERS = ("MCPS_B", "MCPS_V", "IRAC_1")


class Photometer(Model):
    """A model whose output is photometry, so its filter names are strings."""

    fail = False

    def __init__(self, filters: tuple[str, ...] = FILTERS) -> None:
        self.filters = filters
        self.register_parameter(Parameter("norm", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        if type(self).fail:
            raise RuntimeError("the photometer broke")
        ctx = self.context(values)
        pivots = np.linspace(1.0, 4.0, len(self.filters))
        return ModelResult(
            PhotometricPoints(
                self.filters,
                pivots * u.micron,
                ctx["norm"] * pivots**-1.0 * u.Jy,
                extra_coords={"quality": np.arange(len(self.filters), dtype=float)},
            )
        )


def photometric(filters: tuple[str, ...] = FILTERS) -> FittingProblem:
    pivots = np.linspace(1.0, 4.0, len(filters))
    observed = PhotometricPoints(
        filters,
        pivots * u.micron,
        pivots**-1.0 * u.Jy,
        uncertainty=np.full(len(filters), 0.02) * u.Jy,
        extra_coords={"quality": np.arange(len(filters), dtype=float)},
    )
    return FittingProblem(
        Photometer(filters),
        DatasetCollection({"phot": Dataset(observed, label="phot")}),
        seed=20261001,
        simulator_failures=(RuntimeError,),
    )


class TestStringCoordinates:
    def test_filter_names_round_trip_and_the_numeric_label_keeps_its_dtype(
        self, tmp_path: Path
    ) -> None:
        problem = photometric()
        drawn = budget(problem, 3)
        write_training_set(tmp_path / "bank.nc", drawn, problem)
        stored = read_training_set(tmp_path / "bank.nc")
        result = stored.result(1)["default"]
        assert result.filters.tolist() == list(FILTERS)
        assert result == drawn[1].results["model"]["default"]
        observed = stored.observations(2)["phot"]
        assert observed.filters.tolist() == list(FILTERS)
        assert observed.extra_coords["quality"].dtype == np.float64
        assert observed == drawn[2].observations["phot"]

    def test_the_attr_records_the_dtype_beside_the_unchanged_name_list(
        self, tmp_path: Path
    ) -> None:
        import xarray

        problem = photometric()
        write_training_set(tmp_path / "bank.nc", budget(problem, 2), problem)
        tree = xarray.open_datatree(tmp_path / "bank.nc").load()
        attrs = tree["observations/phot"].dataset.attrs
        assert json.loads(attrs[f"{ATTR_PREFIX}extra_coords"]) == ["filters", "quality"]
        assert json.loads(attrs[f"{ATTR_PREFIX}extra_coord_dtypes"]) == {
            "filters": "str",
            "quality": "float64",
        }

    def test_a_longer_name_in_a_later_append_still_fits(self, tmp_path: Path) -> None:
        problem = photometric()
        write_training_set(tmp_path / "bank.nc", budget(problem, 2), problem)
        longer = photometric(("MCPS_B", "MCPS_V", "A_much_longer_filter_name"))
        # The grid check compares the numeric axes, which are the same; the
        # point is that the stored width is not fixed at the first batch's.
        append_training_set(tmp_path / "bank.nc", budget(longer, 2), longer)
        stored = read_training_set(tmp_path / "bank.nc")
        assert len(stored) == 4
        assert stored.observations(3)["phot"].filters.tolist()[-1] == "A_much_longer_filter_name"

    def test_a_failed_sample_is_filled_with_the_empty_string(self, tmp_path: Path) -> None:
        import xarray

        problem = photometric()
        good = budget(problem, 1)
        Photometer.fail = True
        try:
            bad = budget(problem, 1)
        finally:
            Photometer.fail = False
        assert bad[0].failed
        write_training_set(tmp_path / "bank.nc", [*good, *bad], problem)
        tree = xarray.open_datatree(tmp_path / "bank.nc").load()
        filters = np.asarray(tree["model.default"]["extra_filters"].values)
        assert filters[0].tolist() == list(FILTERS)
        assert filters[1].tolist() == ["", "", ""]

    def test_a_string_coordinate_against_a_numeric_one_is_refused_by_name(
        self, tmp_path: Path
    ) -> None:
        problem = photometric()
        drawn = budget(problem, 2)
        original = drawn[1].results["model"]["default"]
        relabelled = PhotometricPoints(
            FILTERS,
            original.spectral_axis.values * u.micron,
            original.values * u.Jy,
            extra_coords={"quality": np.array(["x", "y", "z"])},
        )
        drawn[1] = dataclasses.replace(drawn[1], results={"model": ModelResult(relabelled)})
        with pytest.raises(ResultsError, match=r"extra coordinate 'quality'.*'float64'"):
            write_training_set(tmp_path / "bank.nc", drawn, problem)


class TestFilesWrittenBeforeTheDtypeAttr:
    def test_a_file_without_the_attr_reads_back_unchanged(self, tmp_path: Path) -> None:
        import xarray

        # A numeric-coordinate set: the only kind a pre-W6.14 writer could emit.
        problem = photometric()
        drawn = budget(problem, 3)
        numeric = [
            dataclasses.replace(
                simulation,
                results={
                    "model": ModelResult(
                        Spectrum(
                            np.linspace(1.0, 4.0, 3) * u.micron,
                            np.asarray(simulation.results["model"]["default"].values) * u.Jy,
                            extra_coords={"quality": np.arange(3, dtype=float)},
                        )
                    )
                },
                observations=None,
            )
            for simulation in drawn
        ]
        write_training_set(tmp_path / "new.nc", numeric, problem)
        tree = xarray.open_datatree(tmp_path / "new.nc").load()
        removed = 0
        for node in tree.subtree:
            if f"{ATTR_PREFIX}extra_coord_dtypes" in node.dataset.attrs:
                del node.dataset.attrs[f"{ATTR_PREFIX}extra_coord_dtypes"]
                removed += 1
        assert removed == 1
        tree.to_netcdf(tmp_path / "old.nc", engine="h5netcdf")
        new = read_training_set(tmp_path / "new.nc")
        old = read_training_set(tmp_path / "old.nc")
        assert len(old) == 3
        for index in range(3):
            assert old.result(index)["default"] == new.result(index)["default"]
            assert old.result(index)["default"].extra_coords["quality"].dtype == np.float64


# ---------------------------------------------------------------------------
# Append
# ---------------------------------------------------------------------------


class TestAppend:
    def test_the_sample_dimension_grows_and_the_old_pairs_stay(self, tmp_path: Path) -> None:
        problem = toy()
        first = budget(problem, 3)
        write_training_set(tmp_path / "budget.nc", first, problem)
        second = budget(problem, 2)
        append_training_set(tmp_path / "budget.nc", second, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert len(stored) == 5
        assert stored.attrs[f"{ATTR_PREFIX}samples"] == 5
        channel = next(iter(first[0].results["model"]))
        assert stored.result(0)[channel] == first[0].results["model"][channel]
        assert stored.result(3)[channel] == second[0].results["model"][channel]
        assert np.allclose(
            stored.theta["model.index"][:3], [s.parameters["model.index"] for s in first]
        )

    def test_a_changed_declaration_refuses(self, tmp_path: Path) -> None:
        # DEVELOPMENT_PLAN.md §7's stale-artefact trap: a budget extended after
        # the model was edited would be one file from two simulators.
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        changed = toy(Steeper)
        with pytest.raises(ResultsError, match="different declaration"):
            append_training_set(tmp_path / "budget.nc", budget(changed, 1), changed)

    def test_appending_to_a_run_rather_than_a_training_set_refuses(self, tmp_path: Path) -> None:
        from ampere.results import DrawRecorder, to_netcdf

        problem = toy()
        recorder = DrawRecorder(problem)
        recorder.record()
        to_netcdf(recorder.emit(engine="fixture"), tmp_path / "run.nc")
        with pytest.raises(ResultsError, match="not an ampere training set"):
            append_training_set(tmp_path / "run.nc", budget(problem, 1), problem)

    def test_an_empty_batch_refuses(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 1), problem)
        with pytest.raises(ResultsError, match="batch is empty"):
            append_training_set(tmp_path / "budget.nc", [], problem)

    def test_a_model_swap_with_identical_parameters_refuses_by_model_hash(
        self, tmp_path: Path
    ) -> None:
        """W3.12: the gap ``ampere_spec_hash`` alone leaves open (W3.5's finding).

        A different model class computing something else, declaring the
        *identical* parameters as ``Powerlaw`` — so the spec hash cannot
        tell the two apart, and only ``ampere_model_hash`` can.
        """

        class OtherPowerlaw(Powerlaw):
            """A different class, the same parameter declaration."""

        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        other = toy(OtherPowerlaw)
        assert (
            provenance_attrs(other)[f"{ATTR_PREFIX}spec_hash"]
            == provenance_attrs(problem)[f"{ATTR_PREFIX}spec_hash"]
        )
        assert (
            provenance_attrs(other)[f"{ATTR_PREFIX}model_hash"]
            != provenance_attrs(problem)[f"{ATTR_PREFIX}model_hash"]
        )
        with pytest.raises(ResultsError, match="model_hash"):
            append_training_set(tmp_path / "budget.nc", budget(other, 1), other)

    def test_a_pre_schema_six_file_refuses_by_name_rather_than_guessing(
        self, tmp_path: Path
    ) -> None:
        """A file written before ``ampere_model_hash`` existed has no digest to compare."""
        xarray = pytest.importorskip("xarray")

        problem = toy()
        path = tmp_path / "budget.nc"
        write_training_set(path, budget(problem, 2), problem)
        opened = xarray.open_datatree(str(path))
        tree = opened.load()
        opened.close()
        del tree.attrs[f"{ATTR_PREFIX}model_hash"]
        tree.to_netcdf(str(path))

        with pytest.raises(ResultsError, match="predates"):
            append_training_set(path, budget(problem, 1), problem)

    def test_a_pre_w511_training_set_refuses_on_append_by_name(self, tmp_path: Path) -> None:
        """**W5.11**: the axis identity moved every layout hash, once.

        A budget written before it was packed without the axis-identity
        columns, so its rows are narrower than the ones being appended. The
        refusal has to say *that*, because the alternative — reporting two
        hashes that differ — tells a reader nothing about why or what to do.
        """
        xarray = pytest.importorskip("xarray")
        import json

        from ampere.core.encoding import EncodingLayout

        problem = toy()
        path = tmp_path / "budget.nc"
        write_training_set(path, budget(problem, 2), problem)

        # A training set as an SBI run wrote it before W5.11: the layout it
        # recorded is version 1 and its datasets carry no axis codes.
        stale = EncodingLayout.from_datasets(problem.datasets).to_dict()
        stale["version"] = 1
        for record in stale["datasets"]:
            record.pop("axis_codes")
        opened = xarray.open_datatree(str(path))
        tree = opened.load()
        opened.close()
        tree.attrs[f"{ATTR_PREFIX}encoding_layout"] = json.dumps(stale)
        tree.to_netcdf(str(path))

        with pytest.raises(ResultsError, match="axis-identity") as raised:
            append_training_set(path, budget(problem, 1), problem)
        message = str(raised.value)
        assert "W5.11" in message
        assert "encoding.md" in message

    def test_a_training_set_recording_this_encoding_appends(self, tmp_path: Path) -> None:
        """The other half: a current layout is not refused."""
        xarray = pytest.importorskip("xarray")
        import json

        from ampere.core.encoding import EncodingLayout

        problem = toy()
        path = tmp_path / "budget.nc"
        write_training_set(path, budget(problem, 2), problem)
        current = EncodingLayout.from_datasets(problem.datasets).to_dict()
        opened = xarray.open_datatree(str(path))
        tree = opened.load()
        opened.close()
        tree.attrs[f"{ATTR_PREFIX}encoding_layout"] = json.dumps(current)
        tree.to_netcdf(str(path))

        append_training_set(path, budget(problem, 1), problem)
        assert len(read_training_set(path)) == 3


class TestConcatenateGroupMismatch:
    """W5.28(g): ``_concatenate`` refuses a group mismatch by name, for any group.

    Today's two optional groups (``observations/*``, keyed off the batch's
    own channels, and ``context``) each have their own guard
    (``_check_against_file``, ``_check_context``) raised before
    ``_concatenate`` is ever reached, so neither can trigger the bug this
    closes through the public API. These two tests call the private
    function directly, with hand-built trees, to pin the behaviour a *future*
    optional group would otherwise be silently exposed to: a bare
    ``KeyError`` out of ``addition[name]`` when the file has a group the
    batch does not, and a silent drop of any group the batch has that the
    file does not (the loop only ever walks the file's own group paths).
    """

    def test_a_group_the_file_has_and_the_batch_does_not_is_refused_by_name(self) -> None:
        xarray = pytest.importorskip("xarray")
        from ampere.results.training import _concatenate

        existing = xarray.DataTree.from_dict(
            {
                "theta": xarray.Dataset({"index": ("sample", [1.0, 2.0])}),
                "sample_stats": xarray.Dataset({"status": ("sample", ["ok", "ok"])}),
                "extras": xarray.Dataset({"value": ("sample", [3.0, 4.0])}),
            }
        )
        addition = xarray.DataTree.from_dict(
            {
                "theta": xarray.Dataset({"index": ("sample", [5.0])}),
                "sample_stats": xarray.Dataset({"status": ("sample", ["ok"])}),
                # "extras" is missing: a future optional group the batch was not
                # drawn with. Without this item's fix, addition["extras"] raises
                # a bare KeyError instead of this.
            }
        )
        with pytest.raises(ResultsError, match="extras"):
            _concatenate(xarray, existing, addition)

    def test_a_group_the_batch_has_and_the_file_does_not_is_refused_not_dropped(self) -> None:
        xarray = pytest.importorskip("xarray")
        from ampere.results.training import _concatenate

        existing = xarray.DataTree.from_dict(
            {
                "theta": xarray.Dataset({"index": ("sample", [1.0, 2.0])}),
                "sample_stats": xarray.Dataset({"status": ("sample", ["ok", "ok"])}),
            }
        )
        addition = xarray.DataTree.from_dict(
            {
                "theta": xarray.Dataset({"index": ("sample", [5.0])}),
                "sample_stats": xarray.Dataset({"status": ("sample", ["ok"])}),
                # "bonus" the file was never written with. Without this item's
                # fix, the loop over the file's own group paths never looks at
                # it, and it is silently absent from the merged tree.
                "bonus": xarray.Dataset({"value": ("sample", [9.0])}),
            }
        )
        with pytest.raises(ResultsError, match="bonus"):
            _concatenate(xarray, existing, addition)


class TestTheObservationContext:
    """**W5.10**: the ``context`` group, and what it deliberately does not store."""

    def test_a_contextual_budget_round_trips_its_records(self, tmp_path: Path) -> None:
        problem = toy()
        drawn = problem.simulate_many(4, observe=True, context=ScaledSigma(0.5, 2.0))
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert [dict(record or {}) for record in stored.contexts] == [
            dict(draw.context.record) for draw in drawn
        ]
        assert stored.pair(0)["context"]["kind"] == "scaled_sigma"

    def test_the_prior_is_in_the_root_attributes(self, tmp_path: Path) -> None:
        problem = toy()
        drawn = problem.simulate_many(3, observe=True, context=ScaledSigma(0.5, 2.0))
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        recorded = json.loads(stored.attrs[f"{ATTR_PREFIX}simulation_context"])
        assert recorded == ScaledSigma(0.5, 2.0).describe()

    def test_the_sigma_arrays_are_not_stored_twice(self, tmp_path: Path) -> None:
        """The drawn observation carries them; the group says which draw they came from.

        A context group that also held the sigma arrays would double a
        budget's size to say what the ``observations`` group already says.
        """
        problem = toy()
        drawn = problem.simulate_many(2, observe=True, context=ScaledSigma(2.0, 2.0))
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert set(stored.contexts[0] or {}) == {"kind", "factor"}
        carried = stored.observations(0)["sed"].uncertainty
        assert np.allclose(carried, 2.0 * 0.02)

    def test_a_budget_without_a_context_writes_no_group(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.contexts == ()
        # A plain list of draws carries no batch provenance at all, so the
        # attribute is *absent* rather than 'none' -- which is the honest
        # distinction, and why this pair checks both spellings of no context.
        assert f"{ATTR_PREFIX}simulation_context" not in stored.attrs

    def test_a_contextless_simulate_many_records_none(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", problem.simulate_many(2, observe=True), problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.contexts == ()
        assert stored.attrs[f"{ATTR_PREFIX}simulation_context"] == "none"

    def test_appending_a_contextless_batch_to_a_contextual_set_refuses(
        self, tmp_path: Path
    ) -> None:
        problem = toy()
        drawn = problem.simulate_many(2, observe=True, context=ScaledSigma(0.5, 2.0))
        write_training_set(tmp_path / "budget.nc", drawn, problem)
        with pytest.raises(ResultsError, match="observation context"):
            append_training_set(tmp_path / "budget.nc", budget(problem, 1), problem)

    def test_appending_a_contextual_batch_to_a_contextless_set_refuses(
        self, tmp_path: Path
    ) -> None:
        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        drawn = problem.simulate_many(1, observe=True, context=ScaledSigma(0.5, 2.0))
        with pytest.raises(ResultsError, match="observation context"):
            append_training_set(tmp_path / "budget.nc", drawn, problem)

    def test_two_contextual_batches_append(self, tmp_path: Path) -> None:
        problem = toy()
        first = problem.simulate_many(2, observe=True, context=ScaledSigma(0.5, 2.0))
        write_training_set(tmp_path / "budget.nc", first, problem)
        second = problem.simulate_many(2, observe=True, context=ScaledSigma(0.5, 2.0))
        append_training_set(tmp_path / "budget.nc", second, problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert len(stored) == 4
        assert len(stored.contexts) == 4
        assert all(record is not None for record in stored.contexts)


# ---------------------------------------------------------------------------
# Failures, and the properties the format claims
# ---------------------------------------------------------------------------


class Crashing(Powerlaw):
    """A model that fails on demand — the 2 % crash rate, made deterministic."""

    fail = False

    def evaluate(self, **values: Any) -> ModelResult:
        if type(self).fail:
            raise RuntimeError("the RT code exited 1")
        return super().evaluate(**values)


class TestFailures:
    def test_a_failed_draw_is_recorded_with_its_reason_not_merely_flagged(
        self, tmp_path: Path
    ) -> None:
        problem = toy(Crashing)
        good = budget(problem, 2)
        Crashing.fail = True
        try:
            bad = budget(problem, 1)
        finally:
            Crashing.fail = False
        assert bad[0].failed
        write_training_set(tmp_path / "budget.nc", [*good, *bad], problem)
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.failed.tolist() == [False, False, True]
        record = stored.failures[2]
        assert record is not None
        assert record["reason"] == str(FailureReason.MODEL_FAILED)
        assert "the RT code exited 1" in record["message"]
        assert record["exception_type"] == "RuntimeError"

    def test_a_failed_draw_has_no_result_and_says_so(self, tmp_path: Path) -> None:
        problem = toy(Crashing)
        good = budget(problem, 1)
        Crashing.fail = True
        try:
            bad = budget(problem, 1)
        finally:
            Crashing.fail = False
        write_training_set(tmp_path / "budget.nc", [*good, *bad], problem)
        stored = read_training_set(tmp_path / "budget.nc")
        with pytest.raises(ResultsError, match="failed simulation"):
            stored.result(1)

    def test_a_failed_draws_channel_is_nan_rather_than_zero(self, tmp_path: Path) -> None:
        # results.md §11's first property: "NaN is native, so a masked or
        # crashed sample needs no sentinel". Zero would be a prediction.
        import xarray

        problem = toy(Crashing)
        good = budget(problem, 1)
        Crashing.fail = True
        try:
            bad = budget(problem, 1)
        finally:
            Crashing.fail = False
        write_training_set(tmp_path / "budget.nc", [*good, *bad], problem)
        tree = xarray.open_datatree(tmp_path / "budget.nc").load()
        values = np.asarray(tree["model.default"]["values"].values)
        assert np.all(np.isfinite(values[0]))
        assert np.all(np.isnan(values[1]))

    def test_an_all_failed_batch_can_still_be_appended(self, tmp_path: Path) -> None:
        # A budget's failure rate is data, and the file already knows the
        # shapes, so a batch with nothing to take a shape from is still
        # appendable.
        problem = toy(Crashing)
        write_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        Crashing.fail = True
        try:
            append_training_set(tmp_path / "budget.nc", budget(problem, 2), problem)
        finally:
            Crashing.fail = False
        stored = read_training_set(tmp_path / "budget.nc")
        assert stored.failed.tolist() == [False, False, True, True]

    def test_a_budget_of_nothing_but_failures_refuses_to_be_written(self, tmp_path: Path) -> None:
        problem = toy(Crashing)
        Crashing.fail = True
        try:
            drawn = budget(problem, 2)
        finally:
            Crashing.fail = False
        with pytest.raises(ResultsError, match="nothing to train on"):
            write_training_set(tmp_path / "budget.nc", drawn, problem)


class TestFormatProperties:
    def test_the_coordinates_are_stored_once(self, tmp_path: Path) -> None:
        # results.md §11's second property, which is what makes a
        # coordinate-conditioned emulator's training set the size of a
        # fixed-grid one.
        import xarray

        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 5), problem)
        tree = xarray.open_datatree(tmp_path / "budget.nc").load()
        coordinates = tree["coordinates"].dataset
        assert "model.default_spectral_axis" in coordinates.variables
        assert coordinates["model.default_spectral_axis"].shape == (GRID.size,)

    def test_a_moving_grid_is_refused(self, tmp_path: Path) -> None:
        problem = toy()
        drawn = budget(problem, 2)
        moved = Spectrum(
            (GRID + 0.5) * u.micron,
            drawn[1].results["model"]["default"].values,
            fidelity="cheap",
        )
        shifted = Simulation(
            parameters=drawn[1].parameters,
            theta=drawn[1].theta,
            results={"model": ModelResult({"default": moved})},
            predicted=drawn[1].predicted,
        )
        with pytest.raises(ResultsError, match="stores each channel's coordinates once"):
            write_training_set(tmp_path / "budget.nc", [drawn[0], shifted], problem)

    def test_an_empty_budget_refuses(self, tmp_path: Path) -> None:
        with pytest.raises(ResultsError, match="at least one simulation"):
            write_training_set(tmp_path / "budget.nc", [], toy())

    def test_reading_a_missing_file_says_so(self, tmp_path: Path) -> None:
        with pytest.raises(ResultsError, match="no training set at"):
            read_training_set(tmp_path / "absent.nc")

    def test_a_future_version_is_refused(self, tmp_path: Path) -> None:
        import xarray

        problem = toy()
        write_training_set(tmp_path / "budget.nc", budget(problem, 1), problem)
        tree = xarray.open_datatree(tmp_path / "budget.nc").load()
        tree.attrs[f"{ATTR_PREFIX}training_set_version"] = 99
        tree.to_netcdf(tmp_path / "future.nc")
        with pytest.raises(ResultsError, match="unsupported training-set version"):
            read_training_set(tmp_path / "future.nc")

    def test_a_failure_record_carried_by_hand_is_refused_when_it_is_not_one(self) -> None:
        from ampere.results import training_pair_to_dict

        with pytest.raises(ResultsError, match="must be an ampere Failure"):
            training_pair_to_dict(
                {"model.index": -1.0},
                ModelResult(Spectrum(GRID * u.micron, np.ones(GRID.size) * u.Jy)),
                failure="it broke",
            )

    def test_a_failure_object_is_accepted_by_its_own_to_dict(self) -> None:
        from ampere.results import training_pair_to_dict

        pair = training_pair_to_dict(
            {"model.index": -1.0},
            ModelResult(Spectrum(GRID * u.micron, np.ones(GRID.size) * u.Jy)),
            failed=True,
            failure=Failure(FailureReason.LIKELIHOOD_FAILED, "not positive definite"),
        )
        assert pair["failure"]["reason"] == "likelihood_failed"


# ---------------------------------------------------------------------------
# Writing a budget that never fits in memory (W3.1)
# ---------------------------------------------------------------------------


class TestWritingFromChunks:
    """``simulate_many(as_chunks=True)`` reaches the file a chunk at a time.

    The property that matters is that it makes **no difference**: a set written
    from chunks is the set written from the whole batch, sample for sample.
    What differs is the peak memory (one chunk, not the budget) and the cost —
    each chunk after the first pays §13 limitation 9's ``O(existing + new)``
    append.
    """

    def test_a_batch_is_accepted_wherever_a_list_of_simulations_is(self, tmp_path: Path) -> None:
        problem = toy()
        batch = problem.simulate_many(4, observe=True)
        write_training_set(tmp_path / "batch.nc", batch, problem)
        assert len(read_training_set(tmp_path / "batch.nc")) == 4

    def test_writing_from_chunks_equals_writing_the_whole_batch(self, tmp_path: Path) -> None:
        whole = toy()
        chunked = toy()
        write_training_set(tmp_path / "whole.nc", whole.simulate_many(7, observe=True), whole)
        write_training_set(
            tmp_path / "chunks.nc",
            chunked.simulate_many(7, observe=True, as_chunks=True, chunk_size=3),
            chunked,
        )
        one = read_training_set(tmp_path / "whole.nc")
        other = read_training_set(tmp_path / "chunks.nc")
        assert len(one) == len(other) == 7
        for index in range(7):
            assert one.parameters(index) == other.parameters(index)
            left, right = one.result(index), other.result(index)
            assert sorted(left) == sorted(right)
            for channel in left:
                assert left[channel] == right[channel], (index, channel)
            assert one.observations(index)["sed"] == other.observations(index)["sed"]

    def test_the_chunks_are_pulled_lazily_so_the_budget_is_never_all_held(
        self, tmp_path: Path
    ) -> None:
        problem = toy()
        live: list[int] = []

        def counted() -> Any:
            for chunk in problem.simulate_many(6, observe=True, as_chunks=True, chunk_size=2):
                live.append(len(chunk))
                yield chunk

        write_training_set(tmp_path / "lazy.nc", counted(), problem)
        assert live == [2, 2, 2]
        assert len(read_training_set(tmp_path / "lazy.nc")) == 6

    def test_appending_takes_chunks_too(self, tmp_path: Path) -> None:
        problem = toy()
        write_training_set(tmp_path / "grow.nc", problem.simulate_many(2), problem)
        append_training_set(
            tmp_path / "grow.nc", problem.simulate_many(4, as_chunks=True, chunk_size=2), problem
        )
        assert len(read_training_set(tmp_path / "grow.nc")) == 6

    def test_a_mixed_budget_is_refused_rather_than_half_written(self, tmp_path: Path) -> None:
        problem = toy()
        batch = problem.simulate_many(2)
        with pytest.raises(ResultsError, match="mixes SimulationBatch chunks"):
            write_training_set(tmp_path / "mixed.nc", [batch, batch[0]], problem)

    def test_an_empty_chunk_iterator_is_still_an_empty_budget(self, tmp_path: Path) -> None:
        problem = toy()
        with pytest.raises(ResultsError, match="at least one simulation"):
            write_training_set(tmp_path / "empty.nc", iter(()), problem)


# ---------------------------------------------------------------------------
# W7.0: the derived node, schema 10
# ---------------------------------------------------------------------------


class DerivedNorm(Powerlaw):
    """:class:`Powerlaw` sampled in ``log_norm``, with ``norm`` derived from it."""

    def __init__(self, grid: np.ndarray) -> None:
        self.register_buffer("grid", np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.3)))
        self.register_parameter(Parameter("log_norm", st.norm(0.0, 0.3)))
        self.register_parameter(Parameter("norm", Derived("exp(log_norm)")))


class TestTheDerivedNode:
    """W7.0: a derived parameter is no training column, and is named in the attrs."""

    def test_a_derived_parameter_is_not_a_training_column(self, tmp_path: Path) -> None:
        problem = toy(DerivedNorm)
        drawn = budget(problem, 3)
        write_training_set(tmp_path / "derived.nc", drawn, problem)
        stored = read_training_set(tmp_path / "derived.nc")
        # The free vector is the encoding: the derived value is a function of it.
        assert set(stored.theta) == {"model.index", "model.log_norm"}
        assert json.loads(stored.attrs[f"{ATTR_PREFIX}derived"]) == ["model.norm"]

    def test_the_fingerprint_is_unchanged_for_a_problem_without_one(self) -> None:
        """Schema 10 moves ``ampere_problem_hash`` through the version alone.

        The digests were recorded by computing ``problem.parameters.to_spec()``
        and ``problem_fingerprint(problem)`` with ``"version"`` masked, through
        ``hash_of``, on the base commit ``5963b2f``'s code (before W7.0) for
        ``tests/conformance/test_population.py``'s population problem on the
        reference fixture, in both layouts -- so everything the fingerprint
        hashes except the schema constant is byte-identical across the bump.
        """
        from ampere.results.provenance import hash_of, problem_fingerprint
        from tests.conformance.backends.reference import ReferenceBackend
        from tests.conformance.test_population import build_population_problem

        recorded = {
            "plate": ("f943023d9574b67ef48e52d223917193", "85fcb069130f943ec684e35d746bd5fe"),
            "flat": ("630d727c2f5ece480922c2176fde5a89", "3fcb7b2b3cf54956baee2c6189923eb4"),
        }
        for layout, (spec_digest, fingerprint_digest) in recorded.items():
            problem = build_population_problem(ReferenceBackend(), layout=layout)
            assert not problem.parameters.derived_names
            fingerprint = problem_fingerprint(problem)
            assert fingerprint.pop("version") == PROVENANCE_SCHEMA_VERSION == 10
            assert hash_of(problem.parameters.to_spec()) == spec_digest, layout
            assert hash_of(fingerprint) == fingerprint_digest, layout
            assert provenance_attrs(problem)[f"{ATTR_PREFIX}derived"] == "[]"
