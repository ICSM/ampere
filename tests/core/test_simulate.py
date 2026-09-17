"""``simulate_many``: the batched forward model, the executors, chunking (W3.1).

The one contract every row here is written against is ``inference.md`` §13's
*batched form*: ``simulate_many(n, stream=s)`` returns, in order, exactly what
``[simulate(rng=child) for child in rng(s).spawn(n)]`` returns — whatever
executor ran it and however it was chunked. Everything else (failure counting,
the batch's stacked views, the reference backend's ``BATCHABLE``) is checked
against that.

Every model in this file is defined at **module scope** on purpose: the process
pool pickles the problem, and a class defined inside a test function cannot be
unpickled in a worker. That is not a limitation of the test — it is the same
constraint a user's own simulator is under, and
``FittingProblem._assert_picklable`` refuses it by name.
"""

from __future__ import annotations

import json
import math
import multiprocessing
import os
import pickle
import time
from collections.abc import Sequence
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Capabilities,
    ContainerBatch,
    ContextPrior,
    Dataset,
    DatasetError,
    FailureReason,
    FittingProblem,
    GaussianFamily,
    IndependentNoise,
    Likelihood,
    Model,
    ModelResult,
    Parameter,
    ProcessExecutor,
    ScaledSigma,
    SerialExecutor,
    SigmaArchive,
    SignalToNoise,
    Simulation,
    SimulationBatch,
    Spectrum,
    ThreadExecutor,
    chunk_bounds,
    register_realisation,
)
from ampere.core.encoding import EncodingLayout, encode
from ampere.core.exceptions import LoweringError
from ampere.core.simulate import BatchedPrediction

WAVELENGTH = np.array([1.0, 2.0, 4.0])
SEED = 20260909

# ---------------------------------------------------------------------------
# Module-scope simulators (picklable, as a worker needs)
# ---------------------------------------------------------------------------


class Flat(Model):
    """One channel, one parameter — the simplest thing a budget can simulate."""

    def __init__(self, wavelength: np.ndarray = WAVELENGTH) -> None:
        self.register_buffer("wavelength", np.asarray(wavelength, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("level", st.uniform(0.0, 10.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(
            ctx["wavelength"] * u.micron,
            np.full(ctx["wavelength"].shape, ctx["level"]) * u.Jy,
        )


class Crash(RuntimeError):
    """What a wrapped external simulator raises when the run dies."""


#: The ``level`` at which :class:`Unreliable` hangs. Deliberately **not** 5.0,
#: which is ``uniform(0, 10)``'s prior median and therefore the reference θ this
#: problem is validated at when it is composed — a hang there would be paid by
#: every test in this file that builds one.
HANGS_AT = 7.25


class Unreliable(Flat):
    """Crashes above ``level`` 9.9, and hangs at exactly :data:`HANGS_AT`.

    Two failure modes with different remedies: one the simulator reports (a
    declared ``simulator_failures`` exception, flagged as ``model_failed``) and
    one it never gets to report (a hang, which only the executor can end).
    """

    def evaluate(self, **values: Any) -> Spectrum:
        level = float(self.context(values)["level"])
        if level > 9.9:
            raise Crash("the RT code exited 1")
        if level == HANGS_AT:
            time.sleep(30.0)
        return super().evaluate(**values)


#: The ``level`` at which :class:`HardCrash` takes its worker process with it.
DIES_AT = 8.5


class HardCrash(Flat):
    """Dies the way a segfaulting compiled routine does: the process just goes.

    ``os._exit`` runs no handlers and unwinds nothing, so as far as the pool is
    concerned this is a segfault — and, crucially, it breaks the pool for
    *every* draw in flight, not only for the guilty one.
    """

    def evaluate(self, **values: Any) -> Spectrum:
        if float(self.context(values)["level"]) == DIES_AT:
            os._exit(70)
        return super().evaluate(**values)


class BatchedFlat(Flat):
    """The reference backend's ``BATCHABLE``: a table of θ in one call."""

    BATCHABLE = True
    #: Bumped every time the batch hook is used, so a row can prove it was.
    calls = 0

    def evaluate_batch(self, batch: Sequence[Any]) -> list[Spectrum]:
        type(self).calls += 1
        grid = self.context(batch[0])["wavelength"]
        levels = np.asarray([float(row["level"]) for row in batch])
        table = np.repeat(levels[:, None], grid.size, axis=1)
        return [Spectrum(grid * u.micron, row * u.Jy) for row in table]


class DeclaresOnly(Flat):
    """Declares ``BATCHABLE`` (as a torch model does for ``vmap``) and implements nothing."""

    BATCHABLE = True


class TwoChannel(Model):
    """Two channels off one model, so the batch's per-model results are exercised."""

    def __init__(self) -> None:
        self.register_buffer("grid", WAVELENGTH, unit=u.micron)
        self.register_parameter(Parameter("level", st.uniform(0.0, 10.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        grid = ctx["grid"]
        return ModelResult(
            {
                "default": Spectrum(grid * u.micron, np.full(grid.shape, ctx["level"]) * u.Jy),
                "red": Spectrum(grid * u.micron, np.full(grid.shape, 2 * ctx["level"]) * u.Jy),
            }
        )


class FakeFlat(Flat):
    """:class:`Flat`, declared on a fake native backend (W5.2's owed item 11).

    Exists to exercise :meth:`~ampere.core.dataset._NativeBatch.run` without
    torch or jax installed: a registered :func:`~ampere.core.register_realisation`
    factory is all a backend is, so a minimal one that only ever runs in this
    process is enough to stand in for one.
    """

    BACKEND = "w52-native-fake"
    BATCHABLE = True


class FakeNoise(IndependentNoise):
    """:class:`~ampere.core.IndependentNoise`, declared on the fake backend."""

    BACKEND = "w52-native-fake"
    BATCHABLE = True


class FakeGaussian(GaussianFamily):
    """:class:`~ampere.core.GaussianFamily`, declared on the fake backend."""

    BACKEND = "w52-native-fake"
    BATCHABLE = True


class _FakeRealisation:
    """A native twin minimal enough to run in ``dev``: batched prediction, native draws.

    ``simulate_batched`` reproduces :class:`FakeFlat`'s own arithmetic exactly
    (so the agreement check :meth:`~ampere.core.dataset._NativeBatch._check_agreement`
    passes) and ``sample_observations`` is the row this item is about: a
    *poisoned* level (99.0, chosen the way torch's Poisson twin's negative-rate
    check would find a bad row) fails the whole vectorised call once, exactly
    as ``torch.poisson`` checking a stacked rate array does — and the point of
    :meth:`~ampere.core.dataset._NativeBatch._sample_chunk` is that this must
    flag only the poisoned draw, not the chunk it travelled in.
    """

    def __init__(self, problem: FittingProblem) -> None:
        self._problem = problem

    @property
    def backend(self) -> str:
        return "w52-native-fake"

    @property
    def free_size(self) -> int:
        return self._problem.free_size

    def log_prob_unconstrained(self, theta: Any) -> Any:
        return self._problem.log_prob_unconstrained(np.asarray(theta))

    def simulate_batched(
        self, theta: Any, *, chunk_size: int | None = None, sharder: Any = None
    ) -> BatchedPrediction:
        levels = np.asarray(theta, dtype=float)[:, 0]
        table = np.repeat(levels[:, None], WAVELENGTH.size, axis=1)
        return BatchedPrediction(
            channels={"model": {"default": table}}, predicted={"default": table}
        )

    def sample_observations(
        self, theta: Any, predicted: Any, seeds: Sequence[int]
    ) -> dict[str, np.ndarray]:
        levels = np.asarray(theta, dtype=float)[:, 0]
        if np.any(levels == POISONED_LEVEL):
            raise LoweringError(
                "sample_observations",
                backend="w52-native-fake",
                detail="a poisoned row, by construction.",
            )
        rows = np.asarray(predicted["default"], dtype=float)
        rng = np.random.default_rng(0)
        return {"default": rows + 0.01 * rng.standard_normal(rows.shape)}


#: The level a stacked, whole-chunk check would reject -- torch's own
#: Poisson-twin bug, reproduced in miniature: the failure carries no row
#: index, and only isolating rows finds which one it was.
POISONED_LEVEL = 99.0

#: Registered once at import: every test in :class:`TestTheNativeSamplerIsolatesAFailingDraw`
#: shares it, the way ``ampere.backends.torch``/``.jax`` register once on import.
register_realisation("w52-native-fake", _FakeRealisation)


def fake_problem() -> FittingProblem:
    """A fresh problem, declared entirely on the fake backend."""
    return FittingProblem(
        FakeFlat(),
        [Dataset(observed(), likelihood=Likelihood(FakeGaussian(), FakeNoise()))],
        seed=SEED,
        capabilities=Capabilities(
            differentiable=False, batchable=True, device="cpu", backend="w52-native-fake"
        ),
    )


def observed(values: Any = (1.0, 1.0, 1.0)) -> Spectrum:
    return Spectrum(
        WAVELENGTH * u.micron,
        np.asarray(values, dtype=float) * u.Jy,
        uncertainty=[0.1, 0.1, 0.1] * u.Jy,
    )


def build(model: Model | None = None, **kwargs: Any) -> FittingProblem:
    """A fresh problem on the same seed, so two of them are the same problem."""
    kwargs.setdefault("seed", SEED)
    return FittingProblem(model if model is not None else Flat(), [Dataset(observed())], **kwargs)


def build_unreliable() -> FittingProblem:
    return FittingProblem(
        Unreliable(),
        [Dataset(observed())],
        seed=SEED,
        simulator_failures=(Crash,),
    )


def levels(count: int, **overrides: float) -> np.ndarray:
    """A ``(count, 1)`` θ table of harmless levels, with named rows overridden."""
    table = np.full((count, 1), 1.0)
    for index, value in overrides.items():
        table[int(index.lstrip("_"))] = value
    return table


def _worker_pid(_: Any) -> int:
    """Which process ran this item. Module scope, so a worker can unpickle it."""
    return os.getpid()


def values_of(batch: SimulationBatch, label: str = "default") -> np.ndarray:
    return np.asarray(batch.predicted[label].values)


def same_draws(left: SimulationBatch, right: SimulationBatch) -> bool:
    """Order and values: the whole of what "equals n calls of simulate" asserts."""
    if len(left) != len(right):
        return False
    if not np.array_equal(left.theta, right.theta):
        return False
    for one, other in zip(left, right, strict=True):
        if sorted(one.predicted) != sorted(other.predicted):
            return False
        for label in one.predicted:
            if not np.array_equal(one.predicted[label].values, other.predicted[label].values):
                return False
        if (one.observations is None) != (other.observations is None):
            return False
        if one.observations is not None and other.observations is not None:
            for label in one.observations:
                if not np.array_equal(
                    one.observations[label].values, other.observations[label].values
                ):
                    return False
    return True


def reference_loop(
    problem: FittingProblem, count: int, *, observe: bool = False
) -> SimulationBatch:
    """``n`` calls of ``simulate`` under the batch sub-stream — the contract's right-hand side."""
    children = problem.rng("simulate").spawn(count)
    return SimulationBatch(
        tuple(problem.simulate(observe=observe, rng=child) for child in children)
    )


# ---------------------------------------------------------------------------
# The reference semantics
# ---------------------------------------------------------------------------


class TestTheLoopIsTheContract:
    """§13's batched form: the same draws, in the same order, under every executor."""

    @pytest.mark.parametrize("observe", [False, True], ids=["predicted", "observed"])
    def test_the_batch_is_n_calls_of_simulate_on_the_spawned_children(self, observe: bool) -> None:
        batch = build().simulate_many(6, observe=observe)
        assert same_draws(batch, reference_loop(build(), 6, observe=observe))

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ThreadExecutor(2), ProcessExecutor(2)],
        ids=["serial", "thread-2", "process-2"],
    )
    def test_every_shipped_executor_reproduces_it(self, executor: Any) -> None:
        batch = build().simulate_many(6, observe=True, executor=executor, chunk_size=4)
        assert same_draws(batch, reference_loop(build(), 6, observe=True))

    def test_a_bare_concurrent_futures_pool_reproduces_it_too(self) -> None:
        """The protocol is ``concurrent.futures.Executor``'s, so one *is* an executor."""
        import concurrent.futures

        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            batch = build().simulate_many(6, observe=True, executor=pool)
        assert same_draws(batch, reference_loop(build(), 6, observe=True))

    def test_the_stream_advances_between_batches(self) -> None:
        problem = build()
        first = problem.simulate_many(3)
        second = problem.simulate_many(3)
        assert not np.array_equal(first.theta, second.theta)

    def test_a_named_stream_is_a_different_batch(self) -> None:
        assert not np.array_equal(
            build().simulate_many(4).theta,
            build().simulate_many(4, stream="posterior_predictive").theta,
        )

    def test_an_explicit_generator_is_spawned_from_not_drawn_from(self) -> None:
        """``rng=`` overrides the stream and keeps partition independence."""
        first = build().simulate_many(5, rng=np.random.default_rng(7), chunk_size=2)
        second = build().simulate_many(5, rng=np.random.default_rng(7), chunk_size=5)
        assert same_draws(first, second)


class TestPartitionIndependence:
    """How the work was split must not be visible in the answer."""

    @pytest.mark.parametrize("chunk_size", [1, 2, 7, 6, None], ids=str)
    def test_every_chunk_size_gives_the_same_batch(self, chunk_size: int | None) -> None:
        whole = build().simulate_many(6, observe=True)
        assert same_draws(build().simulate_many(6, observe=True, chunk_size=chunk_size), whole)

    def test_chunking_and_pooling_are_independent_of_each_other(self) -> None:
        whole = build().simulate_many(6, observe=True)
        for chunk_size in (1, 7):
            batch = build().simulate_many(
                6, observe=True, chunk_size=chunk_size, executor=ProcessExecutor(2)
            )
            assert same_draws(batch, whole)

    def test_numpy_spawning_is_cumulative_which_is_why_this_works(self) -> None:
        """The property the design rests on, asserted rather than assumed."""
        one = np.random.default_rng(11)
        many = np.random.default_rng(11)
        piecewise = [child.random() for child in one.spawn(3)] + [
            child.random() for child in one.spawn(4)
        ]
        at_once = [child.random() for child in many.spawn(7)]
        assert piecewise == at_once


class TestGivenTheta:
    """``values=`` is the SBC idiom: simulate at θ somebody else chose."""

    def test_an_array_of_theta_is_simulated_row_by_row(self) -> None:
        table = np.linspace(0.5, 9.5, 5).reshape(5, 1)
        batch = build().simulate_many(5, values=table)
        assert np.allclose(batch.theta, table)
        assert np.allclose(values_of(batch), np.repeat(table, WAVELENGTH.size, axis=1))

    def test_a_sequence_of_mappings_works_too(self) -> None:
        rows = [{"model.level": 2.0}, {"model.level": 3.0}]
        batch = build().simulate_many(2, values=rows)
        assert batch.theta.ravel().tolist() == [2.0, 3.0]

    def test_the_observations_still_come_from_the_per_draw_children(self) -> None:
        table = np.full((4, 1), 3.0)
        first = build().simulate_many(4, values=table, observe=True, chunk_size=1)
        second = build().simulate_many(4, values=table, observe=True, chunk_size=4)
        assert same_draws(first, second)
        drawn = np.asarray(first.observations["default"].values)
        assert not np.allclose(drawn[0], drawn[1])

    def test_one_mapping_is_refused_because_it_would_silently_mean_all_draws(self) -> None:
        with pytest.raises(DatasetError, match="one θ per draw"):
            build().simulate_many(3, values={"model.level": 1.0})

    def test_a_wrongly_shaped_array_names_both_shapes(self) -> None:
        with pytest.raises(DatasetError, match=r"must be \(3, 1\)"):
            build().simulate_many(3, values=np.zeros((4, 1)))

    def test_a_sequence_of_the_wrong_length_is_refused(self) -> None:
        with pytest.raises(DatasetError, match="carries one per"):
            build().simulate_many(3, values=[{"model.level": 1.0}])


# ---------------------------------------------------------------------------
# Failures
# ---------------------------------------------------------------------------


class TestFailuresAreCountedNotRaised:
    """A 2 % crash rate produces 98 usable pairs and a count, on every executor."""

    TABLE = levels(100, _13=9.99, _61=9.99)

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ProcessExecutor(2)],
        ids=["serial", "process-2"],
    )
    def test_two_percent_of_draws_crash_and_the_rest_are_usable(self, executor: Any) -> None:
        problem = build_unreliable()
        batch = problem.simulate_many(100, values=self.TABLE, executor=executor)
        assert int(batch.failed.sum()) == 2
        assert batch.failed[13] and batch.failed[61]
        assert len(batch.usable) == 98
        assert all(not draw.failed for draw in batch.usable)
        assert problem.failure_counts[FailureReason.MODEL_FAILED] == 2

    def test_the_failed_draws_keep_their_theta(self) -> None:
        batch = build_unreliable().simulate_many(100, values=self.TABLE)
        assert batch.theta[13].tolist() == [9.99]

    def test_the_failure_records_survive_and_name_the_simulator(self) -> None:
        batch = build_unreliable().simulate_many(100, values=self.TABLE)
        failure = batch.failures[13]
        assert failure is not None
        assert failure.reason is FailureReason.MODEL_FAILED
        assert failure.exception_type == "Crash"

    def test_the_parent_counts_a_pooled_budget_in_draw_order(self) -> None:
        """``inference.md`` limitation 17.7's aggregation, for this path."""
        problem = build_unreliable()
        problem.simulate_many(100, values=self.TABLE, executor=ProcessExecutor(2))
        assert sum(problem.failure_counts.values()) == 2
        assert [failure.reason for failure in problem.failures] == [
            FailureReason.MODEL_FAILED,
            FailureReason.MODEL_FAILED,
        ]
        assert all(failure.exception_type == "Crash" for failure in problem.failures)

    def test_an_undeclared_exception_still_raises_as_it_does_for_one_draw(self) -> None:
        """Equality with the loop cuts both ways: what ``simulate`` raises, this raises."""
        problem = FittingProblem(Unreliable(), [Dataset(observed())], seed=SEED)
        with pytest.raises(Crash):
            problem.simulate_many(2, values=np.full((2, 1), 9.99))


class TestExecutionFailures:
    """A draw the executor lost is its own kind of failure, with its own code."""

    def test_a_timeout_is_flagged_counted_and_keeps_its_theta(self) -> None:
        problem = build_unreliable()
        table = levels(4, _1=HANGS_AT)
        batch = problem.simulate_many(4, values=table, executor=ProcessExecutor(2, timeout=2.0))
        assert batch.failed.tolist() == [False, True, False, False]
        failure = batch.failures[1]
        assert failure is not None
        assert failure.reason is FailureReason.EXECUTION_FAILED
        assert "timeout" in failure.message
        assert batch.theta[1].tolist() == [HANGS_AT]
        assert problem.failure_counts[FailureReason.EXECUTION_FAILED] == 1
        assert len(batch.usable) == 3

    def test_a_crash_and_a_timeout_are_counted_separately(self) -> None:
        problem = build_unreliable()
        table = levels(6, _1=HANGS_AT, _4=9.99)
        batch = problem.simulate_many(6, values=table, executor=ProcessExecutor(2, timeout=2.0))
        assert int(batch.failed.sum()) == 2
        assert problem.failure_counts[FailureReason.EXECUTION_FAILED] == 1
        assert problem.failure_counts[FailureReason.MODEL_FAILED] == 1
        assert len(batch.usable) == 4

    def test_a_dead_worker_is_attributed_to_the_draw_that_killed_it(self) -> None:
        """Nothing is flagged on suspicion: the innocent neighbours are re-run."""
        problem = FittingProblem(HardCrash(), [Dataset(observed())], seed=SEED)
        batch = problem.simulate_many(5, values=levels(5, _2=DIES_AT), executor=ProcessExecutor(2))
        assert batch.failed.tolist() == [False, False, True, False, False]
        failure = batch.failures[2]
        assert failure is not None
        assert failure.reason is FailureReason.EXECUTION_FAILED
        assert "died" in failure.message
        assert len(batch.usable) == 4
        assert problem.failure_counts[FailureReason.EXECUTION_FAILED] == 1

    def test_the_reason_is_a_distinct_code_in_the_vocabulary(self) -> None:
        assert FailureReason.EXECUTION_FAILED not in (
            FailureReason.MODEL_FAILED,
            FailureReason.LIKELIHOOD_FAILED,
        )
        assert str(FailureReason.EXECUTION_FAILED) == "execution_failed"


# ---------------------------------------------------------------------------
# Picklability
# ---------------------------------------------------------------------------


class TestPicklability:
    """A process pool sends the problem to its workers, so the problem must travel."""

    def test_a_composed_problem_round_trips_through_pickle(self) -> None:
        problem = build()
        restored = pickle.loads(pickle.dumps(problem))
        assert restored.free_size == problem.free_size
        assert np.allclose(
            restored.simulate(problem.reference_values).predicted["default"].values,
            problem.simulate(problem.reference_values).predicted["default"].values,
        )

    @pytest.mark.parametrize(
        "obj",
        [observed(), Flat()(level=2.0), Dataset(observed())],
        ids=["container", "model_result", "dataset"],
    )
    def test_the_pieces_round_trip_too(self, obj: Any) -> None:
        assert pickle.loads(pickle.dumps(obj)) is not None

    def test_a_simulation_round_trips_so_a_worker_can_return_one(self) -> None:
        drawn = build().simulate(observe=True)
        restored = pickle.loads(pickle.dumps(drawn))
        assert np.array_equal(restored.theta, drawn.theta)
        assert np.array_equal(
            restored.observations["default"].values, drawn.observations["default"].values
        )

    def test_an_unpicklable_problem_is_refused_by_name_before_any_worker_starts(self) -> None:
        class Local(Flat):  # a class defined here cannot be unpickled anywhere else
            pass

        problem = FittingProblem(Local(), [Dataset(observed())], seed=SEED)
        with pytest.raises(DatasetError, match="cannot be sent to a worker process"):
            problem.simulate_many(2, executor=ProcessExecutor(1))


# ---------------------------------------------------------------------------
# BATCHABLE on the reference backend
# ---------------------------------------------------------------------------


class TestBatchEvaluation:
    """``BATCHABLE`` here means "give me the table of θ in one call"."""

    def test_the_batch_hook_is_used_under_the_serial_executor(self) -> None:
        BatchedFlat.calls = 0
        batch = build(BatchedFlat()).simulate_many(6, chunk_size=3)
        assert BatchedFlat.calls == 2  # one call per chunk, not per draw
        assert len(batch) == 6

    def test_it_gives_exactly_what_the_loop_gives(self) -> None:
        BatchedFlat.calls = 0
        batched = build(BatchedFlat()).simulate_many(6, observe=True)
        assert BatchedFlat.calls == 1
        assert same_draws(batched, reference_loop(build(Flat()), 6, observe=True))

    def test_a_pool_falls_back_to_the_loop_because_a_table_cannot_be_partitioned(self) -> None:
        BatchedFlat.calls = 0
        batch = build(BatchedFlat()).simulate_many(4, executor=ThreadExecutor(2))
        assert BatchedFlat.calls == 0
        assert len(batch) == 4

    def test_declaring_the_flag_without_implementing_the_hook_uses_the_loop(self) -> None:
        """A torch model declares ``BATCHABLE`` for ``vmap`` and implements nothing here."""
        problem = build(DeclaresOnly())
        assert DeclaresOnly.BATCHABLE
        assert not problem._batch_evaluable()
        assert same_draws(problem.simulate_many(4), reference_loop(build(Flat()), 4))

    def test_the_default_hook_refuses_by_name(self) -> None:
        with pytest.raises(NotImplementedError, match="does not implement evaluate_batch"):
            Flat().evaluate_batch([{"level": 1.0}])

    def test_a_batch_call_that_returns_the_wrong_number_of_rows_is_refused(self) -> None:
        class Short(BatchedFlat):
            def evaluate_batch(self, batch: Sequence[Any]) -> list[Spectrum]:
                return super().evaluate_batch(batch)[:-1]

        with pytest.raises(Exception, match="one result per"):
            Short().call_batch([{"level": 1.0}, {"level": 2.0}])

    def test_a_declared_failure_in_a_batch_call_flags_the_whole_chunk(self) -> None:
        class BatchCrash(BatchedFlat):
            def evaluate_batch(self, batch: Sequence[Any]) -> list[Spectrum]:
                raise Crash("the table came back empty")

        problem = FittingProblem(
            BatchCrash(), [Dataset(observed())], seed=SEED, simulator_failures=(Crash,)
        )
        batch = problem.simulate_many(4, chunk_size=2)
        assert batch.failed.all()
        assert problem.failure_counts[FailureReason.MODEL_FAILED] == 4


# ---------------------------------------------------------------------------
# Chunking, the hook, and the batch's own surface
# ---------------------------------------------------------------------------


class TestChunking:
    def test_as_chunks_yields_batches_that_know_where_they_sit(self) -> None:
        chunks = list(build().simulate_many(7, as_chunks=True, chunk_size=3))
        assert [len(chunk) for chunk in chunks] == [3, 3, 1]
        assert [chunk.offset for chunk in chunks] == [0, 3, 6]

    def test_the_chunks_rejoin_into_the_unchunked_batch(self) -> None:
        chunks = list(build().simulate_many(7, observe=True, as_chunks=True, chunk_size=3))
        assert same_draws(
            SimulationBatch.concatenate(chunks), build().simulate_many(7, observe=True)
        )

    def test_nothing_is_simulated_until_the_iterator_is_consumed(self) -> None:
        problem = build()
        chunks = problem.simulate_many(4, as_chunks=True, chunk_size=2)
        assert problem.simulate().theta.tolist() == build().simulate().theta.tolist()
        assert len(list(chunks)) == 2

    def test_the_hook_is_called_once_per_chunk_with_its_index(self) -> None:
        seen: list[int] = []
        build().simulate_many(7, chunk_size=3, on_chunk=seen.append)
        assert seen == [0, 1, 2]

    def test_the_hook_runs_before_its_chunk_which_is_what_placement_needs(self) -> None:
        order: list[str] = []

        class Recording(Flat):
            def evaluate(self, **values: Any) -> Spectrum:
                order.append("draw")
                return super().evaluate(**values)

        problem = FittingProblem(Recording(), [Dataset(observed())], seed=SEED)
        order.clear()  # composition validates the problem, which evaluates it once
        problem.simulate_many(2, chunk_size=1, on_chunk=lambda index: order.append(f"hook{index}"))
        assert order == ["hook0", "draw", "hook1", "draw"]

    def test_chunk_bounds_covers_the_budget_exactly(self) -> None:
        assert chunk_bounds(7, 3) == [(0, 3), (3, 6), (6, 7)]
        assert chunk_bounds(7, None) == [(0, 7)]
        assert chunk_bounds(0, 3) == []

    def test_a_zero_chunk_size_is_refused_by_name(self) -> None:
        with pytest.raises(DatasetError, match="at least 1"):
            build().simulate_many(3, chunk_size=0)


class TestTheBatchRecord:
    def test_an_empty_budget_is_an_empty_batch_not_an_error(self) -> None:
        batch = build().simulate_many(0)
        assert len(batch) == 0 and batch.theta.shape == (0, 0)

    def test_a_non_integer_budget_is_refused(self) -> None:
        with pytest.raises(DatasetError, match="must be an integer"):
            build().simulate_many(2.5)  # type: ignore[arg-type]

    def test_indexing_gives_the_simulation_every_section_13_consumer_wants(self) -> None:
        batch = build().simulate_many(4, observe=True)
        assert isinstance(batch[2], Simulation)
        assert np.array_equal(batch[2].theta, batch.theta[2])

    def test_a_slice_is_a_batch_that_remembers_its_offset(self) -> None:
        batch = build().simulate_many(6)
        piece = batch[2:5]
        assert isinstance(piece, SimulationBatch)
        assert len(piece) == 3 and piece.offset == 2

    def test_the_stacked_containers_agree_with_the_per_draw_ones(self) -> None:
        batch = build().simulate_many(5, observe=True)
        for index, draw in enumerate(batch):
            assert np.array_equal(
                draw.predicted["default"].values, batch.predicted["default"].values[index]
            )
            assert np.array_equal(
                draw.observations["default"].values, batch.observations["default"].values[index]
            )

    def test_observations_are_none_when_none_were_drawn(self) -> None:
        assert build().simulate_many(3).observations is None

    def test_the_raw_model_results_are_kept_per_draw(self) -> None:
        """Design horizon (c): a training set is ``(θ, ModelResult)`` pairs."""
        batch = build(TwoChannel()).simulate_many(3)
        assert sorted(batch.results) == ["model"]
        assert len(batch.results["model"]) == 3
        assert sorted(batch.results["model"][0]) == ["default", "red"]

    def test_a_failed_draw_leaves_a_nan_row_rather_than_a_gap(self) -> None:
        batch = build_unreliable().simulate_many(3, values=levels(3, _1=9.99))
        stack = batch.predicted["default"]
        assert stack.present.tolist() == [True, False, True]
        assert bool(np.all(np.isnan(stack.values[1])))
        assert stack[1] is None
        assert stack[0] is not None

    def test_iter_chunks_rechunks_an_in_memory_batch(self) -> None:
        batch = build().simulate_many(5)
        pieces = list(batch.iter_chunks(2))
        assert [len(piece) for piece in pieces] == [2, 2, 1]
        assert same_draws(SimulationBatch.concatenate(pieces), batch)

    def test_a_container_batch_rebuilds_the_draws_own_container(self) -> None:
        spectrum = observed()
        stack = ContainerBatch.stack([spectrum, spectrum.with_values([3.0, 4.0, 5.0])])
        assert stack.shape == (2, 3)
        rebuilt = stack[1]
        assert rebuilt is not None
        assert rebuilt.values.tolist() == [3.0, 4.0, 5.0]
        assert rebuilt.axes == spectrum.axes and rebuilt.unit == spectrum.unit

    def test_a_batch_in_which_every_draw_failed_has_no_stack_to_build(self) -> None:
        with pytest.raises(DatasetError, match="every draw in this batch failed"):
            ContainerBatch.stack([None, None])


class TestExecutorRefusals:
    def test_something_without_map_is_refused_by_name(self) -> None:
        with pytest.raises(DatasetError, match="map"):
            build().simulate_many(2, executor=object())  # type: ignore[arg-type]

    def test_an_executor_returning_the_wrong_number_of_results_is_refused(self) -> None:
        class Short:
            def map(self, fn: Any, items: Any) -> list[Any]:
                return [fn(item) for item in items][:-1]

        with pytest.raises(DatasetError, match="one result per item"):
            build().simulate_many(3, executor=Short())

    def test_an_executor_returning_something_else_is_refused(self) -> None:
        class Wrong:
            def map(self, fn: Any, items: Any) -> list[Any]:
                return [None for _ in items]

        with pytest.raises(DatasetError, match="ampere's returns a Simulation"):
            build().simulate_many(2, executor=Wrong())

    def test_a_pool_needs_at_least_one_worker(self) -> None:
        with pytest.raises(DatasetError, match="at least one worker"):
            ProcessExecutor(0)

    def test_a_timeout_must_be_positive(self) -> None:
        with pytest.raises(DatasetError, match="must be positive"):
            ThreadExecutor(2, timeout=0.0)

    def test_a_generator_that_cannot_spawn_is_refused_with_the_remedy(self) -> None:
        """Every generator numpy builds from a seed can spawn; something else cannot."""

        class NotAGenerator:
            def random(self) -> float:
                return 0.5

        with pytest.raises(DatasetError, match="cannot spawn independent children"):
            build().simulate_many(2, rng=NotAGenerator())  # type: ignore[arg-type]


class TestExecutorsInTheirOwnRight:
    """The executors are usable on their own; ``map`` is the whole protocol."""

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ThreadExecutor(2), ProcessExecutor(2)],
        ids=["serial", "thread", "process"],
    )
    def test_results_come_back_in_item_order(self, executor: Any) -> None:
        assert list(executor.map(abs, [-3, 1, -2, 4])) == [3, 1, 2, 4]

    @pytest.mark.parametrize(
        "executor",
        [SerialExecutor(), ThreadExecutor(2), ProcessExecutor(2)],
        ids=["serial", "thread", "process"],
    )
    def test_an_empty_item_list_is_an_empty_result_list(self, executor: Any) -> None:
        assert list(executor.map(abs, [])) == []

    def test_an_exception_in_the_mapped_callable_propagates(self) -> None:
        with pytest.raises(ZeroDivisionError):
            ThreadExecutor(2).map(lambda x: 1 / x, [1, 0, 2])

    def test_the_process_executor_is_a_context_manager(self) -> None:
        with ProcessExecutor(1) as pool:
            assert pool.map(math.sqrt, [4.0]) == [2.0]


class TestThePoolIsReusedAcrossMapCalls:
    """W3.1 slice 2's carried finding: a chunked budget must not rebuild the pool.

    ``chunk_size``'s job is bounding **memory**. Slice 1 built a pool per
    ``map`` call, so it also bounded throughput — every chunk paid for process
    creation plus one pickle of the fitting problem per worker. The rule that
    replaced it is narrow on purpose: *a pool is replaced only when a worker
    died or expired*, which is exactly the two states in which the held pool is
    not usable any more.
    """

    @staticmethod
    def counting(executor: ProcessExecutor, monkeypatch: pytest.MonkeyPatch) -> list[int]:
        """Record every pool this executor builds, without changing what it builds."""
        built: list[int] = []
        original = executor._factory

        def factory(workers: int) -> Any:
            built.append(workers)
            return original(workers)

        monkeypatch.setattr(executor, "_factory", factory)
        return built

    def test_one_pool_serves_a_multi_chunk_budget(self, monkeypatch: pytest.MonkeyPatch) -> None:
        executor = ProcessExecutor(2)
        built = self.counting(executor, monkeypatch)
        try:
            batch = build().simulate_many(6, observe=True, executor=executor, chunk_size=2)
        finally:
            executor.shutdown()
        assert len(batch) == 6
        assert built == [2], "three chunks must share one pool"

    def test_the_same_worker_processes_run_every_chunk(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        """The stronger form: not merely one pool object, the same processes."""
        executor = ProcessExecutor(2)
        try:
            first = set(executor.map(_worker_pid, [0, 1, 2, 3]))
            second = set(executor.map(_worker_pid, [0, 1, 2, 3]))
        finally:
            executor.shutdown()
        assert first == second

    def test_a_dead_worker_does_replace_the_pool(self, monkeypatch: pytest.MonkeyPatch) -> None:
        """The exception to the rule, and the only one that is not a timeout."""
        executor = ProcessExecutor(2)
        built = self.counting(executor, monkeypatch)
        problem = FittingProblem(HardCrash(), [Dataset(observed())], seed=SEED)
        try:
            batch = problem.simulate_many(5, values=levels(5, _2=DIES_AT), executor=executor)
        finally:
            executor.shutdown()
        assert batch.failed.tolist() == [False, False, True, False, False]
        assert len(built) > 1

    def test_shutdown_closes_the_held_pool(self) -> None:
        executor = ProcessExecutor(1)
        executor.map(math.sqrt, [4.0])
        assert executor._source is not None and executor._source.live
        executor.shutdown()
        assert executor._source is None

    def test_a_new_broadcast_payload_closes_the_pool(self) -> None:
        """Workers are handed the payload at start-up, so a changed one needs new workers."""
        executor = ProcessExecutor(1)
        try:
            executor.broadcast("first")
            executor.map(math.sqrt, [4.0])
            assert executor._source is not None
            executor.broadcast("second")
            assert executor._source is None
        finally:
            executor.shutdown()


class TestTheStartMethod:
    """``forkserver`` on POSIX by default (ruled, W3.1 slice 2), ``fork`` on request."""

    def test_the_default_is_forkserver_on_posix(self) -> None:
        if os.name != "posix":  # pragma: no cover - not this CI
            pytest.skip("forkserver is a POSIX start method")
        assert ProcessExecutor(2).start_method == "forkserver"

    def test_fork_is_still_available_through_mp_context(self) -> None:
        if os.name != "posix":  # pragma: no cover - not this CI
            pytest.skip("fork is a POSIX start method")
        executor = ProcessExecutor(2, mp_context=multiprocessing.get_context("fork"))
        assert executor.start_method == "fork"

    def test_a_budget_runs_under_the_default_start_method(self) -> None:
        executor = ProcessExecutor(2)
        try:
            batch = build().simulate_many(4, observe=True, executor=executor)
        finally:
            executor.shutdown()
        assert len(batch) == 4 and not batch.failed.any()

    def test_the_start_method_is_visible_in_the_repr(self) -> None:
        assert f"start_method={ProcessExecutor(2).start_method!r}" in repr(ProcessExecutor(2))


class TestTheNativePathIsAskedForByName:
    """W3.1 slice 2's core-side surface, on a backend that has no native path.

    The reference backend registers no realisation, which makes it exactly the
    right place to assert the *other* half of the contract: what
    ``simulate_many`` does when the fast path is not there. The answer is "runs
    the loop and says so", and the two rows that matter are that the answer is
    recorded rather than inferred, and that asking for the fast path by name
    gets a refusal rather than a silent slow path.
    """

    def test_the_loop_is_the_default_and_the_provenance_records_it(self) -> None:
        batch = build().simulate_many(4, observe=True)
        assert batch.provenance["simulate_batched"] is False
        assert batch.provenance["sample_backend"] == "reference"
        assert batch.provenance["simulation_context"] == "none"

    def test_a_backend_without_a_realisation_refuses_native_true(self) -> None:
        with pytest.raises(DatasetError, match=r"realisation|BATCHABLE"):
            build().simulate_many(4, native=True)

    def test_native_true_with_an_executor_is_refused_by_name(self) -> None:
        """A pool partitions the draws and a vmap evaluates them together."""
        with pytest.raises(DatasetError, match="executor"):
            build().simulate_many(4, native=True, executor=ThreadExecutor(1))

    def test_the_reserved_context_is_recorded(self) -> None:
        """The horizon note's hook: the signature exists before the machinery."""
        batch = build().simulate_many(2, context=None)
        assert batch.provenance["simulation_context"] == "none"

    def test_a_supplied_context_is_refused_rather_than_ignored(self) -> None:
        """Reserving a keyword is not the same as accepting one and dropping it."""
        with pytest.raises(DatasetError, match="reserved"):
            build().simulate_many(2, context={"sigma": 0.1})

    def test_the_provenance_of_a_mixed_budget_claims_least(self) -> None:
        """Concatenation is conservative: a mixture must not read as a native run."""
        native = SimulationBatch((), provenance={"simulate_batched": True})
        loop = SimulationBatch((), provenance={"simulate_batched": False})
        joined = SimulationBatch.concatenate([native, loop])
        assert joined.provenance["simulate_batched"] is False

    def test_the_provenance_survives_re_chunking_and_selection(self) -> None:
        batch = build().simulate_many(4, observe=True)
        assert next(iter(batch.iter_chunks(2))).provenance == batch.provenance
        assert batch.usable.provenance == batch.provenance
        assert batch[1:3].provenance == batch.provenance

    def test_an_evaluate_batch_model_is_not_a_native_batched_path(self) -> None:
        """The two are different claims and the attribute must not conflate them.

        ``Model.evaluate_batch`` is a *model* taking a table of θ;
        ``simulate_batched`` is a *lowered problem* running through a backend's
        ``vmap``. A training set that recorded the first as the second could not
        answer the question the attribute exists for.
        """
        batch = build(BatchedFlat()).simulate_many(4)
        assert batch.provenance["evaluate_batch"] is True
        assert batch.provenance["simulate_batched"] is False

    def test_a_restored_problem_still_simulates(self) -> None:
        """The round trip a forkserver worker makes, asserted where it can be read.

        ``_assert_picklable`` writes *and reads* since W3.1 slice 2, because a
        worker started by ``forkserver`` (the default now) reconstructs the
        problem rather than inheriting it. Reading it back is not enough on its
        own: the failure this row exists for was a restored problem that
        unpickled cleanly and then refused to evaluate, so the assertion is a
        simulation rather than a ``loads``.
        """
        restored = pickle.loads(pickle.dumps(build()))
        drawn = restored.simulate(observe=True)
        assert not drawn.failed
        assert drawn.observations is not None


class TestTheNativeSamplerIsolatesAFailingDraw:
    """Phase 3's owed item (W5.2): one bad draw flags one draw, natively too.

    ``_NativeBatch.run`` vectorises the native sampler over a whole chunk, so
    one θ's failure (a check on the *stacked* array — torch's Poisson twin's
    finding at W3.14) carries no row index. Before this was fixed, that
    exception propagated out of ``simulate_many`` and aborted every draw in
    the chunk, native and non-native alike, rather than flagging the single
    request that earned it — the asymmetry the numpy loop never had, since
    ``Dataset.draw_observation`` already wraps each draw in its own ``try``.
    ``_sample_chunk`` retries one row at a time after the batched call fails,
    which is what these rows hold to :class:`FakeFlat`'s minimal native twin.
    """

    def test_only_the_poisoned_draw_is_flagged(self) -> None:
        problem = fake_problem()
        batch = problem.simulate_many(
            4,
            values=np.array([[1.0], [2.0], [POISONED_LEVEL], [3.0]]),
            native=True,
            observe=True,
        )
        assert list(batch.failed) == [False, False, True, False]
        assert batch.provenance["simulate_batched"] is True
        assert batch.provenance["sample_backend"] == "w52-native-fake"

    def test_the_surviving_draws_keep_their_native_observations(self) -> None:
        """Not merely "not failed" — the native draw itself, not a fallback value."""
        problem = fake_problem()
        batch = problem.simulate_many(
            4,
            values=np.array([[1.0], [2.0], [POISONED_LEVEL], [3.0]]),
            native=True,
            observe=True,
        )
        kept = batch.usable
        assert len(kept) == 3
        observed_values = np.asarray(
            [sim.observations["default"].values for sim in kept.simulations]
        )
        predicted_values = np.asarray([sim.predicted["default"].values for sim in kept.simulations])
        # The fake sampler adds noise (0.01 * standard normal): the draw is
        # close to, but not identical to, the noise-free prediction.
        assert np.all(np.abs(observed_values - predicted_values) > 0.0)
        assert np.all(np.abs(observed_values - predicted_values) < 0.1)

    def test_a_chunk_with_no_poisoned_draw_takes_the_fast_batched_path(self) -> None:
        """The retry is paid only when the batch actually fails once."""
        problem = fake_problem()
        batch = problem.simulate_many(
            4, values=np.array([[1.0], [2.0], [3.0], [4.0]]), native=True, observe=True
        )
        assert not batch.failed.any()

    def test_the_failure_is_recorded_like_any_other(self) -> None:
        problem = fake_problem()
        problem.simulate_many(
            2, values=np.array([[1.0], [POISONED_LEVEL]]), native=True, observe=True
        )
        assert problem.failure_counts[FailureReason.LIKELIHOOD_FAILED] == 1


# ---------------------------------------------------------------------------
# W5.10: the observation context
# ---------------------------------------------------------------------------


class TestTheShippedContextPriors:
    """The three instances the item ships, each drawing what it says it does."""

    def test_they_satisfy_the_protocol(self) -> None:
        for prior in (ScaledSigma(), SigmaArchive([{"default": np.full(3, 0.2)}]), SignalToNoise()):
            assert isinstance(prior, ContextPrior), prior

    def test_scaled_sigma_multiplies_the_observed_pattern(self) -> None:
        prior = ScaledSigma(0.5, 2.0)
        containers = {"default": observed()}
        rng = np.random.default_rng(0)
        for _ in range(20):
            context = prior.draw(rng, containers)
            factor = float(context.record["factor"])
            assert 0.5 <= factor <= 2.0
            assert np.allclose(context.sigma["default"], factor * 0.1)

    def test_a_fixed_scaled_sigma_is_one_factor(self) -> None:
        context = ScaledSigma(3.0, 3.0).draw(np.random.default_rng(1), {"default": observed()})
        assert context.record == {"kind": "scaled_sigma", "factor": 3.0}
        assert np.allclose(context.sigma["default"], 0.3)

    def test_per_dataset_draws_one_factor_each(self) -> None:
        containers = {"a": observed(), "b": observed()}
        context = ScaledSigma(0.1, 10.0, per_dataset=True).draw(
            np.random.default_rng(2), containers
        )
        factors = context.record["factors"]
        assert set(factors) == {"a", "b"}
        assert factors["a"] != factors["b"]

    def test_scaled_sigma_refuses_a_dataset_with_no_uncertainties(self) -> None:
        bare = Spectrum(WAVELENGTH * u.micron, np.ones(3) * u.Jy)
        with pytest.raises(DatasetError, match="has none"):
            ScaledSigma().draw(np.random.default_rng(3), {"default": bare})

    def test_the_archive_draws_one_of_its_entries(self) -> None:
        entries = [{"default": np.full(3, 0.05)}, {"default": np.array([0.1, 0.2, 0.3])}]
        prior = SigmaArchive(entries, names=("quiet", "noisy"))
        rng = np.random.default_rng(4)
        seen = set()
        for _ in range(30):
            context = prior.draw(rng, {"default": observed()})
            index = int(context.record["index"])
            seen.add(index)
            assert context.record["name"] == ("quiet", "noisy")[index]
            assert np.allclose(context.sigma["default"], entries[index]["default"])
        assert seen == {0, 1}

    def test_an_empty_archive_is_refused(self) -> None:
        with pytest.raises(DatasetError, match="empty"):
            SigmaArchive([])

    def test_an_archive_entry_of_the_wrong_shape_is_refused_by_name(self) -> None:
        prior = SigmaArchive([{"default": np.full(5, 0.1)}])
        with pytest.raises(DatasetError, match="shape"):
            prior.draw(np.random.default_rng(5), {"default": observed()})

    def test_an_archive_entry_for_an_unknown_dataset_is_refused_by_name(self) -> None:
        prior = SigmaArchive([{"other": np.full(3, 0.1)}])
        with pytest.raises(DatasetError, match="does not have"):
            prior.draw(np.random.default_rng(6), {"default": observed()})

    def test_signal_to_noise_builds_sigma_from_the_values(self) -> None:
        containers = {"default": observed((1.0, 2.0, 4.0))}
        context = SignalToNoise(20.0, 20.0).draw(np.random.default_rng(7), containers)
        assert context.record == {"kind": "signal_to_noise", "snr": 20.0}
        assert np.allclose(context.sigma["default"], np.array([1.0, 2.0, 4.0]) / 20.0)

    def test_signal_to_noise_floors_a_near_zero_value(self) -> None:
        """The encoding refuses a retained sample whose sigma is not positive."""
        containers = {"default": observed((0.0, 2.0, 4.0))}
        context = SignalToNoise(10.0, 10.0, floor=0.5).draw(np.random.default_rng(8), containers)
        assert float(context.sigma["default"][0]) > 0.0
        # floor * median(|y|) / snr, with median(|0, 2, 4|) = 2.
        assert float(context.sigma["default"][0]) == pytest.approx(0.5 * 2.0 / 10.0)

    def test_signal_to_noise_median_reference_is_one_flat_sigma(self) -> None:
        containers = {"default": observed((1.0, 2.0, 4.0))}
        context = SignalToNoise(10.0, 10.0, reference="median").draw(
            np.random.default_rng(9), containers
        )
        assert np.allclose(context.sigma["default"], 0.2)

    def test_signal_to_noise_needs_no_observed_uncertainties(self) -> None:
        bare = Spectrum(WAVELENGTH * u.micron, np.array([1.0, 2.0, 4.0]) * u.Jy)
        context = SignalToNoise(10.0, 10.0).draw(np.random.default_rng(10), {"default": bare})
        assert np.all(context.sigma["default"] > 0.0)

    def test_each_prior_describes_itself_for_the_provenance(self) -> None:
        assert ScaledSigma(0.5, 2.0).describe()["kind"] == "scaled_sigma"
        assert SignalToNoise().describe()["kind"] == "signal_to_noise"
        record = SigmaArchive([{"default": np.full(3, 0.1)}]).describe()
        # An archive is data: its size and labels identify it, its arrays do not
        # belong in every run's attributes.
        assert record == {"kind": "sigma_archive", "entries": 1, "labels": ["default"]}


class TestABudgetDrawnUnderAContext:
    """``simulate_many(context=...)``: the reserved slot, filled (**W5.10**)."""

    def test_every_draw_carries_its_own_context(self) -> None:
        batch = build().simulate_many(6, observe=True, context=ScaledSigma(0.5, 2.0))
        factors = [float(draw.context.record["factor"]) for draw in batch]
        assert len(set(factors)) == 6
        assert all(0.5 <= factor <= 2.0 for factor in factors)

    def test_the_drawn_observation_carries_the_context_sigma(self) -> None:
        """Which is the whole reason the encoding needs no new column.

        A simulated observation's own ``uncertainty`` is what the packing's
        ``log_sigma`` column is built from, so varying it across the budget is
        what a network conditions on -- with no layout change at all.
        """
        batch = build().simulate_many(4, observe=True, context=ScaledSigma(0.5, 2.0))
        for draw in batch:
            factor = float(draw.context.record["factor"])
            sigma = np.asarray(draw.observations["default"].uncertainty)
            assert np.allclose(sigma, factor * 0.1)

    def test_a_context_is_partition_independent(self) -> None:
        """Draw *i*'s context does not depend on how the budget was chunked."""
        whole = build().simulate_many(6, observe=True, context=ScaledSigma(0.5, 2.0))
        split = build().simulate_many(6, observe=True, context=ScaledSigma(0.5, 2.0), chunk_size=2)
        assert [dict(draw.context.record) for draw in whole] == [
            dict(draw.context.record) for draw in split
        ]
        assert np.allclose(
            whole.observations["default"].values, split.observations["default"].values
        )

    def test_the_context_stream_leaves_theta_alone(self) -> None:
        """The contexts come off their own sub-stream, so a budget's theta is unchanged.

        This is the property that lets the argument exist at all: adding a
        context to a study must not silently re-draw the parameters it was
        comparing against.
        """
        with_context = build().simulate_many(6, observe=True, context=ScaledSigma(0.5, 2.0))
        without = build().simulate_many(6, observe=True)
        assert np.allclose(with_context.theta, without.theta)
        assert not np.allclose(
            with_context.observations["default"].values, without.observations["default"].values
        )

    def test_the_batch_provenance_names_the_prior(self) -> None:
        batch = build().simulate_many(3, observe=True, context=ScaledSigma(0.5, 2.0))
        recorded = json.loads(batch.provenance["simulation_context"])
        assert recorded == ScaledSigma(0.5, 2.0).describe()
        assert build().simulate_many(1, observe=True).provenance["simulation_context"] == "none"

    def test_a_budget_without_a_context_records_none_on_its_draws(self) -> None:
        batch = build().simulate_many(2, observe=True)
        assert [draw.context for draw in batch] == [None, None]

    def test_the_encoding_sees_the_context_through_log_sigma(self) -> None:
        """No layout change: the same hash, a different ``log_sigma`` per draw."""
        problem = build()
        batch = problem.simulate_many(4, observe=True, context=ScaledSigma(0.25, 4.0))
        layout = EncodingLayout.from_datasets(problem.datasets)
        encoded = encode(batch.observations, layout=layout, batched=True)
        column = np.asarray(encoded.values)[:, :, layout.group("log_sigma").offset]
        assert encoded.layout.hash == layout.hash
        assert len(set(np.round(column[:, 0], 9).tolist())) == 4


class TestTheContextRefusals:
    """Each one a claim the code could not honestly make, refused by name."""

    def test_something_that_is_not_a_context_prior(self) -> None:
        with pytest.raises(DatasetError, match="ContextPrior"):
            build().simulate_many(2, observe=True, context="loud")

    def test_a_context_without_observations(self) -> None:
        with pytest.raises(DatasetError, match="observe=False"):
            build().simulate_many(2, observe=False, context=ScaledSigma())

    def test_a_context_with_the_native_path_demanded(self) -> None:
        with pytest.raises(DatasetError, match="native=True"):
            fake_problem().simulate_many(2, observe=True, native=True, context=ScaledSigma())

    def test_the_default_falls_back_to_the_loop_rather_than_refusing(self) -> None:
        batch = fake_problem().simulate_many(2, observe=True, context=ScaledSigma(0.5, 2.0))
        assert batch.provenance["simulate_batched"] is False
        assert all(draw.context is not None for draw in batch)

    def test_a_sigma_of_the_wrong_shape_is_refused_by_the_dataset(self) -> None:
        dataset = Dataset(observed())
        with pytest.raises(DatasetError, match="has to have the observation's shape"):
            dataset.contextual_observed(np.full(5, 0.1))

    def test_a_contextual_container_keeps_everything_but_the_sigma(self) -> None:
        dataset = Dataset(observed((1.0, 2.0, 4.0)))
        swapped = dataset.contextual_observed(np.full(3, 0.5))
        assert np.allclose(swapped.values, [1.0, 2.0, 4.0])
        assert np.allclose(swapped.uncertainty, 0.5)
        assert swapped.unit == dataset.observed.unit
