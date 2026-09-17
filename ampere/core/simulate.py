"""Batched simulation: the executor protocol, the batch record, chunking.

``inference.md`` §13 gives ampere one draw of the forward model. An SBI budget
wants 10⁴ of them, and ``DEVELOPMENT_PLAN.md`` §2's *Batched simulation* row
(Peter's rulings of 2026-09-08, revised 2026-09-09) fixes how they are
produced. Three of those rulings shape this module, and each of them is a
design that was **not** taken:

1. **``vmap`` is not the design.** A vectorised map is one device; a simulation
   budget routinely wants many devices or many machines, and a single
   simulation may not fit on one device at all. So the batch is produced
   through an *executor* — an order-preserving, partition-independent mapper —
   with the serial loop as the reference semantics every executor must
   reproduce exactly. Per-chunk ``vmap`` is a throughput optimisation W3.1
   slice 2 adds *underneath* this contract, never instead of it.
2. **Slow external simulators are the canonical case, not a fallback.** A
   Fortran/C/C++/Rust routine behind a Python call, composed as a black-box
   :class:`~ampere.core.transform.Model` on the reference backend, is what most
   SBI in this field actually simulates with. :class:`ProcessExecutor` is
   therefore a shipped executor rather than an example, with a per-simulation
   ``timeout`` and crash capture, and ``examples/sbi/external_simulator.py``
   is a worked instance of the pattern.
3. **A budget must not have to fit in memory.** ``chunk_size`` bounds how many
   simulations exist at once; ``as_chunks=True`` yields
   :class:`SimulationBatch` chunks that
   :func:`~ampere.results.training.write_training_set` consumes one at a time.

What "the same as the loop" means, exactly
------------------------------------------
``FittingProblem.simulate`` advances its sub-stream, so *n* successive calls
are not partition-independent: draw 7 depends on how many draws preceded it in
this process. A batch cannot be built that way — chunk 2 would depend on
whether chunk 1 ran here or in a worker — so the batch derives **one child
generator per draw, by index**, from the batch sub-stream, using numpy's own
:meth:`numpy.random.Generator.spawn`. The contract is then exact and testable:

    ``problem.simulate_many(n, stream=s)`` produces, in order, exactly what
    ``[problem.simulate(rng=child) for child in problem.rng(s).spawn(n)]``
    produces — whatever executor ran it and however it was chunked.

Spawning is cumulative in numpy (``spawn(a)`` then ``spawn(b)`` gives the same
children as ``spawn(a + b)``), so the children are spawned a chunk at a time
without the batch ever holding *n* generators, and ``chunk_size=1`` and
``chunk_size=7`` give bit-identical batches.

θ itself is drawn **in the parent**, before the work is handed out, and the
already-advanced child generator travels with it. That is not merely an
optimisation: it is what lets a draw whose worker was killed by a timeout still
record the θ it was killed on, which is the difference between
"reject-and-record" and "lose the evidence".

Examples
--------
>>> from ampere.core.simulate import SerialExecutor
>>> SerialExecutor().map(lambda x: x * 2, [1, 2, 3])
[2, 4, 6]
"""

from __future__ import annotations

import concurrent.futures
import dataclasses
import functools
import multiprocessing
import os
import time
import types
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Protocol, TypeVar, runtime_checkable

import numpy as np

from .exceptions import DatasetError
from .results_schema import FunctionSamples, ModelResult

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .dataset import Failure, Simulation

__all__ = [
    "BatchedPrediction",
    "ChunkHook",
    "ChunkSharder",
    "ContainerBatch",
    "ContextPrior",
    "ExecutionFailure",
    "Executor",
    "ObservationContext",
    "ProcessExecutor",
    "ScaledSigma",
    "SerialExecutor",
    "SigmaArchive",
    "SignalToNoise",
    "SimulationBatch",
    "ThreadExecutor",
    "chunk_bounds",
    "worker_shared",
]

_Item = TypeVar("_Item")
_Result = TypeVar("_Result")

#: The per-chunk hook ``simulate_many(on_chunk=...)`` takes, called with the
#: chunk's index **before** that chunk is simulated.
#:
#: This is the device-placement extension point Peter's note (a) asks for. A
#: single simulation larger than one device is *not* partitioned by ampere — a
#: model-parallel simulator is a :class:`~ampere.core.transform.Model` whose
#: ``__call__`` does its own placement — but ``chunk_size=1`` plus this hook
#: guarantees such a model is never asked to hold two simulations at once, and
#: gives it somewhere to move its parameters onto the device it is about to use.
#: W3.1 slice 2 wires the torch and jax backends into it; in this slice it is a
#: documented callable and nothing more.
ChunkHook = Callable[[int], None]


@runtime_checkable
class ChunkSharder(Protocol):
    """How one chunk's stacked θ is spread over the devices a backend reports.

    **Peter's addendum of 2026-09-09**, landed at W3.1 slice 2 as an API and a
    single-device implementation, with the distributed ones designed and
    smoke-tested but not exercised — full exercise waits for the GPU item.

    The layering is the point. :class:`Executor` distributes *simulations*
    across processes and machines; a sharder distributes **one chunk's
    vectorised evaluation** across the accelerators of one host, which is a
    different axis and a different mechanism (``jax.pmap``/``shard_map``,
    ``torch.distributed``/``DTensor``). Neither subsumes the other: an SBI
    budget on eight GPUs wants both, and a budget on one CPU wants neither.

    It is deliberately not a distributed runtime. Two methods, both of which a
    single-device implementation answers trivially, so that the default costs
    nothing and the degenerate case is what CPU CI exercises:

    ``devices()``
        The devices this sharder will use, named. Length 1 is the degenerate
        case; length 0 is not allowed, because "no device" is not a thing a
        chunk can be evaluated on.
    ``shard(fn, stacked)``
        Evaluate *fn* — already vectorised over the leading axis — on
        *stacked*, returning what ``fn(stacked)`` would have returned. A
        single-device implementation *is* ``fn(stacked)``. A multi-device one
        splits the leading axis, runs the parts in parallel and reassembles, so
        the contract is a **value** contract: sharding must not change the
        answer, only where it was computed.
    """

    def devices(self) -> tuple[str, ...]:  # pragma: no cover - protocol
        """The devices this sharder spreads a chunk over, in order."""
        ...

    def shard(self, fn: Callable[[Any], Any], stacked: Any) -> Any:  # pragma: no cover - protocol
        """``fn(stacked)``, however it chooses to compute it."""
        ...


@dataclasses.dataclass(frozen=True)
class BatchedPrediction:
    """What a realisation's ``simulate_batched`` hands back: one chunk, noise-free.

    ``inference.md`` §13's *batched form*, native path. Two mappings of plain
    numpy arrays with a **leading sample axis**, and plain numpy on purpose:
    this crosses out of a backend into ``ampere.core``, which owns no array
    type but numpy's, and the consumer is going to build core containers from
    it anyway.

    Attributes
    ----------
    channels
        ``{model label: {channel: (batch, n)}}`` — every channel of every model
        the problem holds, which is what
        :attr:`SimulationBatch.results` is rebuilt from. Not merely the
        channels a dataset happens to be bound to: ``results.md`` §11 writes one
        training-set group per ``<model>.<channel>``, so a native path that
        returned fewer would silently write a smaller file than the loop.
    predicted
        ``{dataset label: (batch, n_full)}`` — the instrument-transformed
        prediction on **every** observed sample, masked ones included, because
        that is the shape ``Dataset.predict`` returns and masking is the
        consumer's (``draw_observation``'s) business.
    """

    channels: Mapping[str, Mapping[str, np.ndarray]]
    predicted: Mapping[str, np.ndarray]

    def __len__(self) -> int:
        """The batch size, read off whichever stack is present."""
        for stack in self.predicted.values():
            return int(np.shape(stack)[0])
        for produced in self.channels.values():
            for stack in produced.values():
                return int(np.shape(stack)[0])
        return 0


# ---------------------------------------------------------------------------
# The executor protocol
# ---------------------------------------------------------------------------


@runtime_checkable
class Executor(Protocol):
    """Anything that maps a callable over items and returns results **in order**.

    Deliberately the narrowest useful shape, and deliberately the shape
    :class:`concurrent.futures.Executor` already has: dask's
    ``Client.get_executor()``, ray's executor wrappers and mpi4py's
    ``MPIPoolExecutor`` all satisfy it as they stand (Peter, 2026-09-09: "dask
    and/or ray can easily provide replacements"), so scaling a budget out to a
    cluster needs no ampere-specific adapter. A queue-driven cluster array is
    the one case that wants a thin wrapper, and the wrapper is this one method.

    Two properties are contractual, and both are properties of the *result*
    rather than of the schedule:

    **Order.** ``map`` returns results positionally, result *i* for item *i*,
    regardless of completion order.

    **Purity of partition.** How the executor groups, distributes or retries
    items must not change the values it returns. ``simulate_many`` guarantees
    the randomness half of that by handing each draw its own generator, so an
    executor only has to avoid reordering.

    **What is flagged and what is raised.** A timeout or a dead worker comes
    back as an :class:`ExecutionFailure` in that item's slot, which
    ``simulate_many`` turns into a flagged
    :class:`~ampere.core.dataset.Failure`: the executor lost the work, and
    losing work is not an answer the caller should have to catch. An exception
    the *mapped callable* raises propagates instead, because that is exactly
    what the same simulator would have done on the serial path — a declared
    simulator failure is already flagged inside ``simulate``, so an exception
    that gets this far is an undeclared one, and turning it into a silently
    dropped draw would hide a bug in a plausible failure rate. It does take the
    rest of the chunk with it, since a generator-based ``map`` cannot be resumed
    after it raises.
    """

    def map(
        self, fn: Callable[[_Item], _Result], items: Iterable[_Item], /
    ) -> Iterable[_Result]:  # pragma: no cover - protocol
        """Apply *fn* to every item, returning the results in item order."""
        ...


@dataclasses.dataclass(frozen=True)
class ExecutionFailure:
    """A draw the *executor* could not deliver, as opposed to one that failed.

    The distinction matters for counting. A simulator that raises has produced
    an answer of a kind — ``inference.md`` §11's flagged failure, with the
    reason and the offending values — whereas a worker killed by a timeout or a
    segfault produced nothing at all, and merging the two would hide a broken
    cluster inside a plausible-looking prior-failure rate.

    Plain data, so it travels back from a worker process by pickle.
    """

    message: str
    exception_type: str = ""
    timed_out: bool = False


class SerialExecutor:
    """The default: one draw after another, in this process.

    The reference semantics, in the strong sense — every other executor is
    correct exactly insofar as it reproduces what this one returns. It is also
    the right choice more often than it looks: a fast analytic simulator spends
    more time being pickled than being run, and it is the only executor under
    which ``simulate_many`` can use a model's batched
    :meth:`~ampere.core.transform.Model.evaluate_batch`, since a table of θ
    evaluated in one call is by definition not partitioned across workers.
    """

    __slots__ = ()

    def map(self, fn: Callable[[_Item], _Result], items: Iterable[_Item], /) -> list[_Result]:
        """Apply *fn* to every item in order."""
        return [fn(item) for item in items]

    def __repr__(self) -> str:
        return "<SerialExecutor>"


class _PoolSource:
    """Where :func:`_windowed_map` gets a pool, and where a dead one is replaced.

    Split out at W3.1 slice 2, to make a pool outlive one ``map`` call. Slice 1
    built a pool per call, which meant a chunked budget paid the whole start-up
    cost — process creation plus one pickle of the fitting problem per worker —
    once per chunk, and ``chunk_size`` (whose job is bounding *memory*) silently
    became a throughput knob. A persistent source keeps the pool across calls
    and hands out the same one; an ephemeral source is the old behaviour, which
    is still what a thread pool wants, having nothing to amortise.

    Two ways a pool leaves: :meth:`discard`, which is the only thing that ever
    throws one away, and it is called exactly where slice 1 rebuilt — a draw
    that blew its timeout (whose worker is still running and must be killed) or
    a worker that died and broke the pool for everyone. So the rule the item
    asks for holds literally: *a pool is replaced only when a worker died or
    expired.*
    """

    __slots__ = ("_factory", "_persistent", "_pool", "_warm", "_workers")

    def __init__(
        self,
        factory: Callable[[], concurrent.futures.Executor],
        *,
        persistent: bool,
        workers: int,
        warm: bool = False,
    ) -> None:
        self._factory = factory
        self._persistent = persistent
        self._workers = workers
        self._warm = warm
        self._pool: concurrent.futures.Executor | None = None

    @property
    def workers(self) -> int:
        """How many workers the pool this source hands out was built for."""
        return self._workers

    @property
    def live(self) -> bool:
        """Whether a pool is currently held (what the reuse rows assert on)."""
        return self._pool is not None

    def acquire(self) -> concurrent.futures.Executor:
        """The pool to run on, building — and warming — one if none is held."""
        if self._pool is None:
            self._pool = self._factory()
            if self._warm:
                _warm_up(self._pool, self._workers)
        return self._pool

    def fresh(self) -> concurrent.futures.Executor:
        """An independent, throwaway pool — never the shared one.

        Crash attribution (:func:`_rerun_alone`) re-runs each survivor on a pool
        of its own precisely so that a second crash convicts a second draw; the
        shared pool must not be the thing that dies for it.
        """
        return self._factory()

    def discard(self, *, kill: bool) -> None:
        """Throw the held pool away. *kill* stops workers that are still busy."""
        pool = self._pool
        self._pool = None
        if pool is None:
            return
        if kill:
            _terminate(pool)
        pool.shutdown(wait=not kill, cancel_futures=True)

    def release(self) -> None:
        """End of a ``map`` call: an ephemeral source closes, a persistent one keeps."""
        if not self._persistent:
            self.discard(kill=False)


def _ready() -> bool:
    """The no-op a warm-up submits. Module level, so it pickles."""
    return True


def _warm_up(pool: concurrent.futures.Executor, workers: int) -> None:
    """Start every worker **before** the first deadline is taken.

    ``timeout`` is a per-simulation deadline measured from submission, which is
    only the same thing as a deadline measured from the start of the work if
    the worker is already running when the item is submitted. Under ``fork``
    that held to within milliseconds; under ``forkserver`` (the default since
    W3.1 slice 2) the first submission also pays for the fork server booting
    and the pool's initializer unpickling the fitting problem, which is seconds
    on a real problem — enough to expire an honest draw and blame the
    simulator for it.

    So the workers are started with a no-op that is *not* under a deadline, and
    with a persistent pool that cost is paid once per pool rather than once per
    chunk. Failures are swallowed deliberately: a pool that cannot start is a
    pool the real submission will report through the ordinary
    :class:`ExecutionFailure` path, and duplicating that here would only change
    which line the user sees.
    """
    try:
        for future in [pool.submit(_ready) for _ in range(max(1, workers))]:
            future.result()
    except Exception:  # pragma: no cover - reported by the real submission
        pass


def _windowed_map(
    source: _PoolSource,
    fn: Callable[[_Item], _Result],
    items: Sequence[_Item],
    *,
    timeout: float | None,
    recover: bool,
) -> list[Any]:
    """Run *items* on a pool, at most ``source.workers`` in flight, with a deadline.

    The window is what makes ``timeout`` mean "per simulation". Submitting the
    whole budget at once and timing each future from submission would expire
    every queued item at the same moment on a budget larger than the pool; with
    at most one item in flight per worker, submission time *is* start time to
    within the scheduler's own latency.

    *recover* says whether a dead pool can be replaced (processes: yes;
    threads: there is nothing to replace and nothing to kill).
    """
    try:
        return _windowed_loop(source, fn, items, timeout=timeout, recover=recover)
    finally:
        source.release()


def _windowed_loop(
    source: _PoolSource,
    fn: Callable[[_Item], _Result],
    items: Sequence[_Item],
    *,
    timeout: float | None,
    recover: bool,
) -> list[Any]:
    """:func:`_windowed_map` without the release, which its caller owns."""
    workers = source.workers
    results: list[Any] = [None] * len(items)
    pending = list(range(len(items)))
    while pending:
        pool = source.acquire()
        queue = pending
        pending = []
        inflight: list[tuple[int, concurrent.futures.Future[Any], float]] = []
        cursor = 0
        rebuild = False
        broken = False
        try:
            while cursor < len(queue) or inflight:
                while cursor < len(queue) and len(inflight) < workers:
                    index = queue[cursor]
                    cursor += 1
                    inflight.append((index, pool.submit(fn, items[index]), time.monotonic()))
                index, future, started = inflight.pop(0)
                remaining = (
                    None if timeout is None else max(0.0, timeout - (time.monotonic() - started))
                )
                try:
                    results[index] = future.result(timeout=remaining)
                except concurrent.futures.TimeoutError:
                    results[index] = ExecutionFailure(
                        message=(
                            f"the simulation did not finish within the executor's "
                            f"{timeout:g} s per-draw timeout"
                        ),
                        exception_type="TimeoutError",
                        timed_out=True,
                    )
                    if recover:
                        # A running task cannot be cancelled, and leaving it to
                        # run would keep a worker (and, at interpreter exit, the
                        # process) alive for as long as the simulator wants. The
                        # pool is therefore torn down and the *other* in-flight
                        # draws are requeued on a fresh one — they were not the
                        # ones that expired, so they are re-run rather than
                        # flagged.
                        pending = [entry[0] for entry in inflight] + queue[cursor:]
                        rebuild = True
                        break
                except concurrent.futures.BrokenExecutor:
                    # A worker died outright (a segfaulting Fortran routine is
                    # the case this exists for), and killing one worker breaks
                    # the pool for *every* draw in flight — including the one
                    # being awaited, which may well be innocent. The OS does not
                    # say which draw did it, so nothing is flagged on suspicion:
                    # every affected draw, the awaited one first, is re-run one
                    # per fresh single-worker pool, where whichever draw is
                    # guilty convicts itself and the innocent ones simply
                    # produce their results.
                    survivors = [index, *(entry[0] for entry in inflight), *queue[cursor:]]
                    broken = True
                    if recover:
                        _rerun_alone(source, fn, items, results, survivors, timeout)
                    else:
                        # A thread pool breaks only when it cannot start a
                        # worker at all — interpreter shutdown, or a hard
                        # resource limit — and there is nothing to replace it
                        # with, so the survivors are lost with it and say so.
                        for survivor in survivors:
                            results[survivor] = ExecutionFailure(
                                message=(
                                    "the executor broke while this draw was queued; the draw "
                                    "produced no result"
                                ),
                                exception_type="BrokenExecutor",
                            )
                    inflight.clear()
                    break
        finally:
            # The only two exits that throw a pool away, and they are the two
            # the reuse rule names: a draw that expired (whose worker is still
            # busy, and ``wait=False`` after ``kill`` is exactly what the
            # timeout said to do) and a worker that died and broke the pool for
            # everyone. Anything else leaves the pool for the next call.
            if rebuild or broken:
                source.discard(kill=rebuild)
    return results


def _rerun_alone(
    source: _PoolSource,
    fn: Callable[[_Item], _Result],
    items: Sequence[_Item],
    results: list[Any],
    survivors: Sequence[int],
    timeout: float | None,
) -> None:
    """Re-run *survivors* one per fresh pool, so a crash is attributed correctly."""
    for index in survivors:
        pool = source.fresh()
        try:
            future = pool.submit(fn, items[index])
            try:
                results[index] = future.result(timeout=timeout)
            except concurrent.futures.TimeoutError:
                results[index] = ExecutionFailure(
                    message=(
                        f"the simulation did not finish within the executor's "
                        f"{timeout:g} s per-draw timeout"
                    ),
                    exception_type="TimeoutError",
                    timed_out=True,
                )
                _terminate(pool)
            except concurrent.futures.BrokenExecutor:
                results[index] = ExecutionFailure(
                    message=(
                        "the worker process running this simulation died; the draw "
                        "produced no result"
                    ),
                    exception_type="BrokenProcessPool",
                )
        finally:
            pool.shutdown(wait=False, cancel_futures=True)


def _terminate(pool: concurrent.futures.Executor) -> None:
    """Kill a process pool's workers outright.

    ``shutdown(cancel_futures=True)`` cancels what has not started and waits for
    what has; there is no public way to stop a task that is already running, and
    a simulator that has blown its timeout is precisely a task that will not
    stop on its own. Reaching for ``_processes`` is therefore deliberate, and
    guarded: a pool without it (a thread pool, a foreign executor) is left
    alone.
    """
    processes = getattr(pool, "_processes", None)
    if not processes:
        return
    for process in list(processes.values()):
        try:
            process.kill()
        except (OSError, ValueError):  # pragma: no cover - already gone
            pass


class ThreadExecutor:
    """A thread pool, for simulators that spend their time waiting.

    The case is an I/O-bound wrapper — a simulator reached over a socket, a
    routine that shells out and blocks on the child, a library that releases the
    GIL for the duration of a compiled call. For anything that holds the GIL
    this is slower than :class:`SerialExecutor`, and for a genuine external
    simulator :class:`ProcessExecutor` is the right answer.

    ``timeout`` is honoured, with one honest caveat that cannot be engineered
    away: **a Python thread cannot be interrupted.** The deadline therefore
    bounds the *answer* — the draw is flagged at the deadline and its result
    discarded — but not the *work*: ``map`` still returns only once the
    abandoned thread has finished, because leaving a non-daemon worker thread
    running would hang the interpreter at exit instead. A timeout that must also
    stop the work needs :class:`ProcessExecutor`, which can kill its worker.

    Parameters
    ----------
    max_workers
        Threads in the pool. ``None`` lets :mod:`concurrent.futures` choose.
    timeout
        Seconds per simulation, or ``None`` for no deadline.
    """

    def __init__(self, max_workers: int | None = None, *, timeout: float | None = None) -> None:
        self.max_workers = None if max_workers is None else int(max_workers)
        self.timeout = None if timeout is None else float(timeout)
        if self.max_workers is not None and self.max_workers < 1:
            raise DatasetError(f"a thread pool needs at least one worker, got {max_workers!r}.")
        if self.timeout is not None and self.timeout <= 0.0:
            raise DatasetError(f"a per-simulation timeout must be positive, got {timeout!r}.")

    def map(self, fn: Callable[[_Item], _Result], items: Iterable[_Item], /) -> list[Any]:
        """Apply *fn* to every item on the pool, returning results in item order."""
        materialised = list(items)
        if not materialised:
            return []
        workers = self.max_workers or min(32, len(materialised))
        source = _PoolSource(
            functools.partial(concurrent.futures.ThreadPoolExecutor, max_workers=workers),
            persistent=False,
            workers=workers,
        )
        return _windowed_map(source, fn, materialised, timeout=self.timeout, recover=False)

    def __repr__(self) -> str:
        return f"<ThreadExecutor max_workers={self.max_workers}, timeout={self.timeout}>"


class ProcessExecutor:
    """A process pool — the route for external simulators.

    One worker process per simulation in flight, so a compiled routine that
    holds the GIL, leaks global state, or dies on a bad parameter combination
    costs one draw rather than the run. Three behaviours follow from that, and
    all three are why this is a shipped executor rather than a recipe:

    * **The problem is sent to each worker once.** ``simulate_many`` calls
      :meth:`broadcast`, which puts the
      :class:`~ampere.core.dataset.FittingProblem` in the pool's ``initializer``
      — pickled once per worker at start-up rather than once per draw. A foreign
      executor cannot be told that, so it pays the pickle per draw; for a
      simulator slow enough to want a process pool, that difference is noise,
      but it is why ampere ships its own.
    * **A timeout kills the worker.** The expired draw is flagged and its worker
      killed; the draws that were merely in flight beside it are re-run.
    * **A crash is attributed, not smeared.** When a worker dies the pool breaks
      for everyone, so the surviving in-flight draws are re-run one per fresh
      pool: whichever draw is guilty convicts itself.

    Parameters
    ----------
    max_workers
        Worker processes. ``None`` lets :mod:`concurrent.futures` choose.
    timeout
        Seconds per simulation, or ``None`` for no deadline.
    max_tasks_per_child
        Retire and replace a worker after this many draws. ``1`` gives a fresh
        process per simulation, which is the strongest isolation available and
        the right setting for a routine that cannot be run twice in one process.
    mp_context
        A :mod:`multiprocessing` context, if this executor's default start
        method is wrong for this simulator. **The default is ``forkserver`` on
        POSIX** (ruled for W3.1 slice 2), not the platform's own: forking a
        process that already holds a threaded runtime — jax says so itself, and
        torch's intra-op pools have the same shape — risks a deadlock in the
        child, and Python 3.14 changes the platform default away from ``fork``
        for exactly that reason, so ampere is not going to inherit a default
        that is about to move under it. ``forkserver`` forks children from a
        small, clean server process instead, which keeps most of ``fork``'s
        start-up saving without inheriting the parent's threads. The cost is
        that everything sent to a worker must **pickle** and be importable by
        name from a module: a simulator defined at module scope is fine, one
        defined inside a function or a test body is not. Pass
        ``get_context("fork")`` to have the old behaviour back, or
        ``get_context("spawn")`` where a library demands it — ``fork`` is still
        available, it is simply no longer what you get by saying nothing.
        ``None`` (the default) means ampere's choice, not the platform's.

    Notes
    -----
    **Not every problem can be pooled, and the refusal is the contract.** A
    problem composed on the reference backend pickles, which is the case that
    matters: a wrapped external simulator is composed there. One composed on jax
    does not — a jax array carries a process-local ``Device`` handle with no
    pickle reduction — so
    :meth:`~ampere.core.dataset.FittingProblem.simulate_many` refuses this
    executor for it by name, before a worker starts. Throughput on a device
    backend is per-chunk ``vmap`` (W3.1 slice 2), not more processes.

    Examples
    --------
    >>> with ProcessExecutor(max_workers=2) as pool:
    ...     pool.map(abs, [-1, 2, -3])
    [1, 2, 3]
    """

    def __init__(
        self,
        max_workers: int | None = None,
        *,
        timeout: float | None = None,
        max_tasks_per_child: int | None = None,
        mp_context: Any = None,
    ) -> None:
        self.max_workers = None if max_workers is None else int(max_workers)
        self.timeout = None if timeout is None else float(timeout)
        self.max_tasks_per_child = None if max_tasks_per_child is None else int(max_tasks_per_child)
        self.mp_context = _default_context() if mp_context is None else mp_context
        self._shared: Any = None
        self._source: _PoolSource | None = None
        if self.max_workers is not None and self.max_workers < 1:
            raise DatasetError(f"a process pool needs at least one worker, got {max_workers!r}.")
        if self.timeout is not None and self.timeout <= 0.0:
            raise DatasetError(f"a per-simulation timeout must be positive, got {timeout!r}.")

    def broadcast(self, value: Any) -> None:
        """Send *value* to every worker once, at start-up.

        Retrieved inside the worker with :func:`worker_shared`. ``simulate_many``
        uses this to place the fitting problem, which is why the pickling cost
        of a large problem is paid per worker rather than per draw.

        A payload is delivered by the pool's ``initializer``, so a pool already
        running was given the *previous* one. Broadcasting a different value
        therefore closes the held pool: the alternative is workers scoring a
        problem the caller has replaced.
        """
        if self._source is not None and value is not self._shared:
            self._close_pool()
        self._shared = value

    def _factory(self, workers: int) -> concurrent.futures.Executor:
        return concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=self.mp_context,
            initializer=_install_shared,
            initargs=(self._shared,),
            max_tasks_per_child=self.max_tasks_per_child,
        )

    def _pool_source(self, workers: int) -> _PoolSource:
        """The persistent source for this executor, built or reused.

        Reused when the pool it holds is wide enough for the window this call
        wants — the window bounds concurrency, so a wider pool is never wrong —
        which is what makes a chunked budget share one pool: chunks are equal
        except the last, and the last is smaller.
        """
        source = self._source
        if source is None or source.workers < workers:
            self._close_pool()
            source = _PoolSource(
                functools.partial(self._factory, workers),
                persistent=True,
                workers=workers,
                warm=True,
            )
            self._source = source
        return source

    def _close_pool(self) -> None:
        if self._source is not None:
            self._source.discard(kill=False)
            self._source = None

    def map(self, fn: Callable[[_Item], _Result], items: Iterable[_Item], /) -> list[Any]:
        """Apply *fn* to every item on the pool, returning results in item order."""
        materialised = list(items)
        if not materialised:
            return []
        workers = self.max_workers or min(len(materialised), _default_workers())
        return _windowed_map(
            self._pool_source(workers),
            fn,
            materialised,
            timeout=self.timeout,
            recover=True,
        )

    def shutdown(self) -> None:
        """Close the held pool and drop the broadcast payload.

        Idempotent, and worth calling: since W3.1 slice 2 the pool **outlives**
        a ``map`` call, so an executor kept in a long-lived object keeps worker
        processes alive until this runs (or until it is garbage-collected and
        the pool's own finaliser does it). The context-manager form calls it.
        """
        self._close_pool()
        self._shared = None

    def __enter__(self) -> ProcessExecutor:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.shutdown()

    @property
    def start_method(self) -> str:
        """The multiprocessing start method workers are created with."""
        context = self.mp_context
        if context is None:
            return multiprocessing.get_start_method()
        return str(context.get_start_method())

    def __repr__(self) -> str:
        return (
            f"<ProcessExecutor max_workers={self.max_workers}, timeout={self.timeout}, "
            f"max_tasks_per_child={self.max_tasks_per_child}, "
            f"start_method={self.start_method!r}>"
        )


def _default_workers() -> int:
    return max(1, os.cpu_count() or 1)


def _default_context() -> Any:
    """``forkserver`` on POSIX, the platform's own elsewhere (ruled, W3.1 slice 2).

    Two reasons, and neither is a preference. **Forking a threaded runtime may
    deadlock**: jax warns about it in as many words, torch's intra-op pools have
    the same shape, and a budget that hangs in a worker is the worst failure
    mode this executor has, because it looks like a slow simulator. And
    **Python 3.14 moves the platform default off ``fork`` on Linux** for that
    reason, so inheriting the platform default would mean ampere's behaviour
    changing under it at an interpreter upgrade rather than at a decision.

    ``forkserver`` over ``spawn`` because it keeps most of the start-up saving:
    the server process is paid for once and each worker is a fork of it, rather
    than a fresh interpreter re-importing the world. The price is the same as
    ``spawn``'s — everything a worker touches must pickle and be importable by
    name — which ``ProcessExecutor`` already required of the fitting problem.
    """
    if os.name != "posix":
        return None
    try:
        return multiprocessing.get_context("forkserver")
    except ValueError:  # pragma: no cover - POSIX without forkserver
        return None


#: Where :meth:`ProcessExecutor.broadcast`'s payload lands inside a worker.
_SHARED: list[Any] = [None]


def _install_shared(value: Any) -> None:
    """Worker ``initializer``: keep the broadcast payload for this process."""
    _SHARED[0] = value


def worker_shared() -> Any:
    """The value :meth:`ProcessExecutor.broadcast` sent to this worker.

    ``None`` outside a worker, or when nothing was broadcast.
    """
    return _SHARED[0]


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------


def chunk_bounds(count: int, chunk_size: int | None) -> list[tuple[int, int]]:
    """``[(start, stop), ...]`` covering ``range(count)`` in chunks of *chunk_size*.

    ``None`` means one chunk. Spelled out as a function because the chunk
    boundaries are the thing partition independence is asserted against.

    Examples
    --------
    >>> chunk_bounds(7, 3)
    [(0, 3), (3, 6), (6, 7)]
    >>> chunk_bounds(7, None)
    [(0, 7)]
    """
    if count < 0:
        raise DatasetError(f"a simulation budget cannot be negative, got {count!r}.")
    if chunk_size is None:
        return [(0, count)] if count else []
    size = int(chunk_size)
    if size < 1:
        raise DatasetError(
            f"chunk_size bounds how many simulations exist at once, so it must be at least 1, "
            f"got {chunk_size!r}. Pass chunk_size=None for a single chunk."
        )
    return [(start, min(start + size, count)) for start in range(0, count, size)]


# ---------------------------------------------------------------------------
# The observation context (W5.10)
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class ObservationContext:
    """One draw's observation context: a sigma pattern, and what drew it.

    The **reserved** ``simulate_many(context=...)`` slot, filled by W5.10.
    ``inference.md`` §13's simulator draws each observation under the observed
    container's own uncertainties; a posterior trained on that budget has seen
    exactly one error bar per sample and is amortised over *noise
    realisations* but not over *noise levels*. A context prior varies them, and
    this is one draw of it.

    Two fields, and the split is the point.

    ``sigma``
        Label to a sigma array of the observed container's shape, in the
        container's own value unit. It is substituted for the container's
        uncertainty at :meth:`~ampere.core.dataset.Dataset.draw_observation`,
        so every family's ``sample`` receives ``NoiseParams`` built from *it*
        rather than from the observation — which is what makes this a draw
        from a different noise level and not a rescaling applied afterwards.
        A label the mapping omits keeps the observed container's own sigma.
    ``record``
        A small **JSON-safe** description of what the prior drew — a factor,
        an archive index, a signal-to-noise ratio. This is what a training set
        stores per draw (:data:`~ampere.results.training.CONTEXT_GROUP`) and
        what a reader needs to answer "which context was this draw made
        under?". The sigma arrays themselves are not stored twice: a drawn
        observation *carries* its uncertainties, and they are already in the
        file.

    Plain data with no generator and no prior inside it, so it pickles to a
    worker process like every other part of a ``_DrawRequest``.
    """

    sigma: Mapping[str, np.ndarray]
    record: Mapping[str, Any]

    def __post_init__(self) -> None:
        frozen: dict[str, np.ndarray] = {}
        for label, values in dict(self.sigma).items():
            array = np.array(np.asarray(values, dtype=float), copy=True)
            array.setflags(write=False)
            frozen[str(label)] = array
        object.__setattr__(self, "sigma", types.MappingProxyType(frozen))
        object.__setattr__(self, "record", types.MappingProxyType(dict(self.record)))

    def for_label(self, label: str) -> np.ndarray | None:
        """This dataset's sigma, or ``None`` where the context does not vary it."""
        return self.sigma.get(label)

    def to_dict(self) -> dict[str, Any]:
        """The JSON-safe record a training set stores (the sigma arrays excluded)."""
        return dict(self.record)


@runtime_checkable
class ContextPrior(Protocol):
    """What ``simulate_many(context=...)`` draws each simulation's context from.

    Two methods, and no more: a prior that can draw a context and describe
    itself is everything the simulation path, the training set's provenance
    and an SBI run's attrs need.

    Implementations ship in this module — :class:`ScaledSigma`,
    :class:`SigmaArchive`, :class:`SignalToNoise` — and a user's own is any
    object with these two methods, checked structurally
    (``runtime_checkable``) rather than by inheritance, exactly as
    :class:`Executor` is.

    **The generator is given, never created.** ``inference.md`` §12: a context
    prior draws from the run's own seed derivation, so a seeded problem's
    budget is reproducible with contexts exactly as it is without them. An
    implementation that reaches for ``np.random.default_rng()`` breaks that,
    silently.
    """

    def draw(
        self, rng: np.random.Generator, observed: Mapping[str, FunctionSamples]
    ) -> ObservationContext:
        """One draw's context, given *observed*, the problem's own containers."""
        ...  # pragma: no cover - protocol

    def describe(self) -> Mapping[str, Any]:
        """A JSON-safe description of **this prior**, for the provenance."""
        ...  # pragma: no cover - protocol


def _observed_sigma(container: FunctionSamples, label: str, *, prior: str) -> np.ndarray:
    """The observed sigma a context prior scales, refused by name where there is none."""
    if container.uncertainty is None:
        raise DatasetError(
            f"the {prior} context prior varies dataset {label!r}'s uncertainties and that "
            f"dataset's observed container has none. A context prior over sigma needs an observed "
            f"sigma pattern to work from; give the dataset uncertainties, or use a context prior "
            f"that builds sigma from the values (SignalToNoise) instead."
        )
    return np.asarray(container.uncertainty, dtype=float)


def _draw_factor(rng: np.random.Generator, low: float, high: float, *, log: bool) -> float:
    """One factor in ``[low, high]``, log-uniform by default."""
    if low == high:
        return float(low)
    if log:
        return float(np.exp(rng.uniform(np.log(low), np.log(high))))
    return float(rng.uniform(low, high))


@dataclasses.dataclass(frozen=True)
class ScaledSigma:
    """Scaled copies of the **observed** sigma pattern: the simplest context prior.

    One factor per draw, applied to every dataset's uncertainties, so the
    *shape* of the error bars across a spectrum is the observation's own and
    only their overall level varies. That is the context a network most often
    needs to be amortised over — the same instrument on a brighter or a
    fainter night — and the one whose coverage is easiest to state: a
    posterior trained under ``ScaledSigma(0.5, 2.0)`` is being asked to hold
    for observations between half and twice the observed noise, and W5.10's
    acceptance row is exactly that claim and its failure outside the range.

    Parameters
    ----------
    low, high
        The factor range. ``low == high`` is a fixed rescale, which is what a
        calibration check at one noise level uses.
    log
        Draw the factor log-uniformly (the default: a factor is a scale).
    per_dataset
        Draw an independent factor for each dataset rather than one shared by
        all. Off by default, because the usual physical statement is "this
        observation is noisier", not "these instruments are independently
        noisier".
    """

    low: float = 0.5
    high: float = 2.0
    log: bool = True
    per_dataset: bool = False

    def __post_init__(self) -> None:
        if not (self.low > 0.0 and self.high >= self.low):
            raise DatasetError(
                f"ScaledSigma needs 0 < low <= high (a sigma factor is a positive scale), got "
                f"low={self.low!r}, high={self.high!r}."
            )

    def draw(
        self, rng: np.random.Generator, observed: Mapping[str, FunctionSamples]
    ) -> ObservationContext:
        """One factor (or one per dataset), applied to the observed sigma."""
        shared = _draw_factor(rng, self.low, self.high, log=self.log)
        sigma: dict[str, np.ndarray] = {}
        factors: dict[str, float] = {}
        for label, container in observed.items():
            factor = (
                _draw_factor(rng, self.low, self.high, log=self.log) if self.per_dataset else shared
            )
            factors[label] = factor
            sigma[label] = factor * _observed_sigma(container, label, prior="ScaledSigma")
        record: dict[str, Any] = {"kind": "scaled_sigma"}
        if self.per_dataset:
            record["factors"] = factors
        else:
            record["factor"] = shared
        return ObservationContext(sigma=sigma, record=record)

    def describe(self) -> Mapping[str, Any]:
        """This prior, as the provenance records it."""
        return {
            "kind": "scaled_sigma",
            "low": float(self.low),
            "high": float(self.high),
            "log": bool(self.log),
            "per_dataset": bool(self.per_dataset),
        }


@dataclasses.dataclass(frozen=True)
class SigmaArchive:
    """An archive of **real** error arrays, drawn from uniformly.

    The honest context prior where the noise has structure a factor cannot
    reproduce — a detector's read-noise floor beside a photon-limited
    continuum, a night with one bad order — and the one a real survey can
    actually supply: the error columns of the observations it already has.

    Parameters
    ----------
    entries
        A sequence of contexts. Each entry is a mapping of dataset label to a
        sigma array of that dataset's observed shape. A label an entry omits keeps
        the observed container's own sigma, so an archive over one dataset of a
        two-dataset problem is written as such rather than padded.
    names
        Optional one name per entry (a file name, a night, a programme id),
        recorded per draw so a stored budget says *which* archived error array
        each simulation used. Indices are recorded either way.
    """

    entries: Sequence[Mapping[str, np.ndarray]]
    names: Sequence[str] | None = None

    def __post_init__(self) -> None:
        if not len(self.entries):
            raise DatasetError(
                "SigmaArchive is an archive of real error arrays and this one is empty; give it "
                "at least one entry, or use ScaledSigma to vary the observed sigma instead."
            )
        if self.names is not None and len(self.names) != len(self.entries):
            raise DatasetError(
                f"SigmaArchive was given {len(self.entries)} entry(s) and "
                f"{len(self.names)} name(s); there must be one name per entry or none at all."
            )

    def draw(
        self, rng: np.random.Generator, observed: Mapping[str, FunctionSamples]
    ) -> ObservationContext:
        """One archived error array per draw, uniformly over the archive."""
        index = int(rng.integers(len(self.entries)))
        entry = self.entries[index]
        sigma: dict[str, np.ndarray] = {}
        for label, values in entry.items():
            if label not in observed:
                raise DatasetError(
                    f"SigmaArchive entry {index} carries sigma for dataset {label!r}, which this "
                    f"problem does not have; its datasets are {sorted(observed)}."
                )
            array = np.asarray(values, dtype=float)
            expected = np.asarray(observed[label].values).shape
            if array.shape != expected:
                raise DatasetError(
                    f"SigmaArchive entry {index}: dataset {label!r}'s sigma array has shape "
                    f"{array.shape} and the observed container has shape {expected}. An "
                    f"archived error array stands in for the observation's own, so it has to "
                    f"have the observation's shape."
                )
            sigma[label] = array
        record: dict[str, Any] = {"kind": "sigma_archive", "index": index}
        if self.names is not None:
            record["name"] = str(self.names[index])
        return ObservationContext(sigma=sigma, record=record)

    def describe(self) -> Mapping[str, Any]:
        """This prior, as the provenance records it.

        The arrays are **not** described: an archive is data, and a provenance
        attribute that carried it would put a survey's error columns in every
        run's attrs. Its size, its labels and its names are what identify it.
        """
        found: dict[str, Any] = {
            "kind": "sigma_archive",
            "entries": len(self.entries),
            "labels": sorted({label for entry in self.entries for label in entry}),
        }
        if self.names is not None:
            found["names"] = [str(name) for name in self.names]
        return found


@dataclasses.dataclass(frozen=True)
class SignalToNoise:
    """A parametric S/N model: sigma from the *values*, at a drawn signal-to-noise.

    The context prior for a survey that is specified rather than observed —
    "this instrument reaches S/N 20 to 100 on a source like this" — and the
    one that needs no observed uncertainties at all, which is why it is the
    answer when a dataset has none.

    sigma is ``|y| / snr`` per sample, with a **floor** of
    ``floor * median(|y|) / snr`` so that a sample whose value is near zero
    does not get a sigma of zero: the encoding refuses a retained sample whose
    sigma is not strictly positive (``encoding.md`` §6), and a likelihood that
    divides by it would not survive one either.

    Parameters
    ----------
    low, high
        The signal-to-noise range. ``low == high`` is a fixed S/N.
    log
        Draw log-uniformly (the default: S/N is a ratio).
    floor
        The floor, as a fraction of the median ``|y|``. ``0`` removes it and
        is refused where it would produce a zero sigma, by the encoding, by name.
    reference
        ``"values"`` (the default) puts sigma proportional to ``|y|``, a constant fractional
        error; ``"median"`` puts one sigma on every sample, ``median(|y|) / snr``,
        which is the flat error bar of a background-limited observation.
    """

    low: float = 10.0
    high: float = 100.0
    log: bool = True
    floor: float = 1e-3
    reference: str = "values"

    def __post_init__(self) -> None:
        if not (self.low > 0.0 and self.high >= self.low):
            raise DatasetError(
                f"SignalToNoise needs 0 < low <= high (a signal-to-noise ratio is positive), "
                f"got low={self.low!r}, high={self.high!r}."
            )
        if self.reference not in ("values", "median"):
            raise DatasetError(
                f"SignalToNoise's reference= is 'values' (sigma proportional to |y|, a "
                f"constant fractional error) or 'median' (one flat sigma, the "
                f"background-limited case), got {self.reference!r}."
            )

    def draw(
        self, rng: np.random.Generator, observed: Mapping[str, FunctionSamples]
    ) -> ObservationContext:
        """One S/N per draw, turned into a sigma array per dataset."""
        snr = _draw_factor(rng, self.low, self.high, log=self.log)
        sigma: dict[str, np.ndarray] = {}
        for label, container in observed.items():
            magnitude = np.abs(np.asarray(container.values, dtype=float))
            finite = magnitude[np.isfinite(magnitude)]
            median = float(np.median(finite)) if finite.size else 1.0
            if not (median > 0.0):
                median = 1.0
            if self.reference == "median":
                sigma[label] = np.full(magnitude.shape, median / snr)
                continue
            floor = float(self.floor) * median
            sigma[label] = np.maximum(magnitude, floor) / snr
        return ObservationContext(sigma=sigma, record={"kind": "signal_to_noise", "snr": snr})

    def describe(self) -> Mapping[str, Any]:
        """This prior, as the provenance records it."""
        return {
            "kind": "signal_to_noise",
            "low": float(self.low),
            "high": float(self.high),
            "log": bool(self.log),
            "floor": float(self.floor),
            "reference": str(self.reference),
        }


# ---------------------------------------------------------------------------
# The batch record
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class ContainerBatch:
    """One channel's or dataset's containers across a batch, stacked.

    The leading axis is the sample axis: ``values`` has shape
    ``(count, *template.shape)``. It is **not** a
    :class:`~ampere.core.results_schema.FunctionSamples`, and deliberately so —
    ``results_schema.md`` fixes a container's value shape to the one its
    coordinates imply, so a container with a sample axis would be a container
    whose values do not match its axes. What is shared across the batch (kind,
    axes, unit, extra coordinates, metadata) is held once, in ``template``,
    which is the same economy ``results.md`` §11 gives the on-disk training set;
    what varies per draw is an array.

    Indexing gives the draw's own container back, so a consumer that wants
    containers gets containers:

    >>> import astropy.units as u
    >>> from ampere.core import Spectrum
    >>> spectrum = Spectrum([1.0, 2.0] * u.um, [1.0, 2.0] * u.Jy)
    >>> batch = ContainerBatch.stack([spectrum, spectrum.with_values([3.0, 4.0])])
    >>> batch.values.shape
    (2, 2)
    >>> batch[1].values.tolist()
    [3.0, 4.0]

    Attributes
    ----------
    template
        A container of the right kind, on the right coordinates, whose values
        are the first present draw's.
    values
        ``(count, *shape)``. Rows for absent draws (a failed simulation) are
        NaN, or zero where the dtype cannot represent NaN.
    uncertainty, mask
        Stacked the same way where the containers carry them, else ``None``.
    present
        ``(count,)`` boolean: which draws contributed a container.
    """

    template: FunctionSamples
    values: np.ndarray
    uncertainty: np.ndarray | None
    mask: np.ndarray | None
    present: np.ndarray

    @classmethod
    def stack(cls, containers: Sequence[FunctionSamples | None]) -> ContainerBatch:
        """Stack a batch's containers, ``None`` marking a draw that produced none."""
        first = next((one for one in containers if one is not None), None)
        if first is None:
            raise DatasetError(
                "a container batch needs at least one draw that produced a container; every "
                "draw in this batch failed."
            )
        count = len(containers)
        shape = first.values.shape
        values = _blank(count, shape, np.asarray(first.values).dtype)
        has_uncertainty = first.uncertainty is not None
        has_mask = first.mask is not None
        uncertainty = _blank(count, shape, np.float64) if has_uncertainty else None
        mask = np.zeros((count, *shape), dtype=bool) if has_mask else None
        present = np.zeros(count, dtype=bool)
        for index, container in enumerate(containers):
            if container is None:
                continue
            values[index] = np.asarray(container.values)
            present[index] = True
            if uncertainty is not None and container.uncertainty is not None:
                uncertainty[index] = np.asarray(container.uncertainty)
            if mask is not None and container.mask is not None:
                mask[index] = np.asarray(container.mask, dtype=bool)
        for array in (values, uncertainty, mask, present):
            if array is not None:
                array.setflags(write=False)
        return cls(
            template=first, values=values, uncertainty=uncertainty, mask=mask, present=present
        )

    def __len__(self) -> int:
        return int(self.values.shape[0])

    def __getitem__(self, index: int) -> FunctionSamples | None:
        """The *index*-th draw's container, or ``None`` where the draw produced none."""
        if not self.present[index]:
            return None
        return self.template.with_values(
            self.values[index],
            uncertainty=None if self.uncertainty is None else self.uncertainty[index],
            mask=None if self.mask is None else self.mask[index],
        )

    @property
    def shape(self) -> tuple[int, ...]:
        """The stacked shape, sample axis first."""
        return tuple(int(size) for size in self.values.shape)

    def __repr__(self) -> str:
        return (
            f"<ContainerBatch {type(self.template).__name__} {self.shape}, "
            f"{int(np.count_nonzero(self.present))} present>"
        )


def _blank(count: int, shape: tuple[int, ...], dtype: np.dtype[Any]) -> np.ndarray:
    """An array whose "this draw produced nothing" value says so where it can."""
    array = np.empty((count, *shape), dtype=dtype)
    array[...] = np.nan if array.dtype.kind in "fc" else 0
    return array


@dataclasses.dataclass(frozen=True)
class SimulationBatch:
    """``count`` draws of the forward model, in order, failures included.

    What :meth:`~ampere.core.dataset.FittingProblem.simulate_many` returns.
    Every draw is here, failed ones as well: ``inference.md`` §13's
    reject-and-record is only reject-and-record if the record survives, and a
    budget's failure rate is a measurement of the prior worth keeping.

    A batch is a sequence of :class:`~ampere.core.dataset.Simulation`\\ s —
    ``batch[i]`` is the *i*-th draw and every §13 consumer keeps working
    unchanged — plus stacked views over the columns an SBI trainer actually
    wants: :attr:`theta`, :attr:`predicted`, :attr:`observations`,
    :attr:`failed`. The per-draw ``ModelResult``s are kept, not summarised,
    because design horizon (c)'s emulator training sets are ``(θ, ModelResult)``
    pairs.

    Attributes
    ----------
    simulations
        The draws, in order.
    stream
        The batch sub-stream the draws were derived from.
    offset
        Index of the first draw within the whole budget, so a chunk knows where
        it sits. ``batch.offset + i`` is draw ``batch[i]``'s global index.
    provenance
        **W3.1 slice 2**: how this batch was produced, for the training set's
        root attributes. Keys are unprefixed and
        :func:`~ampere.results.provenance_attrs` adds the ``ampere_``; two are
        written today, and both answer a question a stored budget cannot
        otherwise be asked. ``simulate_batched`` says whether the noise-free
        prediction came from a backend's vectorised path or from the loop, and
        ``sample_backend`` names the backend that drew the observations —
        because the numpy path and a native one draw from the *same*
        distribution but not from the same random stream, so a budget is only
        reproducible against the path that produced it.
    """

    simulations: tuple[Simulation, ...]
    stream: str = "simulate"
    offset: int = 0
    provenance: Mapping[str, Any] = dataclasses.field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "simulations", tuple(self.simulations))
        object.__setattr__(self, "provenance", dict(self.provenance))

    # -- sequence surface ----------------------------------------------------

    def __len__(self) -> int:
        return len(self.simulations)

    def __iter__(self) -> Iterator[Simulation]:
        return iter(self.simulations)

    def __getitem__(self, index: int | slice) -> Any:
        """``batch[i]`` is a :class:`~ampere.core.dataset.Simulation`; a slice is a batch."""
        if isinstance(index, slice):
            span = range(len(self))[index]
            return SimulationBatch(
                tuple(self.simulations[index]),
                stream=self.stream,
                offset=self.offset + (span.start if len(span) else 0),
                provenance=self.provenance,
            )
        return self.simulations[index]

    @property
    def count(self) -> int:
        """How many draws this batch holds."""
        return len(self.simulations)

    # -- stacked columns -----------------------------------------------------

    @functools.cached_property
    def theta(self) -> np.ndarray:
        """``(count, free_size)`` — the flat free-parameter vectors, in draw order."""
        if not self.simulations:
            return np.empty((0, 0), dtype=float)
        stacked = np.stack([np.asarray(draw.theta) for draw in self.simulations])
        stacked.setflags(write=False)
        return stacked

    @functools.cached_property
    def failed(self) -> np.ndarray:
        """``(count,)`` boolean: which draws must not be trained on."""
        mask = np.array([draw.failed for draw in self.simulations], dtype=bool)
        mask.setflags(write=False)
        return mask

    @property
    def failures(self) -> tuple[Failure | None, ...]:
        """The per-draw :class:`~ampere.core.dataset.Failure`, ``None`` where none."""
        return tuple(draw.failure for draw in self.simulations)

    @functools.cached_property
    def predicted(self) -> Mapping[str, ContainerBatch]:
        """Noise-free predictions per dataset label, stacked on a leading sample axis."""
        return self._stack("predicted")

    @functools.cached_property
    def observations(self) -> Mapping[str, ContainerBatch] | None:
        """Noisy draws per dataset label, or ``None`` when the batch drew none."""
        if all(draw.observations is None for draw in self.simulations):
            return None
        return self._stack("observations")

    @functools.cached_property
    def results(self) -> Mapping[str, tuple[ModelResult | None, ...]]:
        """The raw ``ModelResult`` per model label, per draw — design horizon (c)."""
        labels: list[str] = []
        for draw in self.simulations:
            for label in draw.results:
                if label not in labels:
                    labels.append(label)
        return {
            label: tuple(draw.results.get(label) for draw in self.simulations) for label in labels
        }

    def _stack(self, attribute: str) -> dict[str, ContainerBatch]:
        labels: list[str] = []
        for draw in self.simulations:
            holding = getattr(draw, attribute) or {}
            for label in holding:
                if label not in labels:
                    labels.append(label)
        stacked: dict[str, ContainerBatch] = {}
        for label in labels:
            stacked[label] = ContainerBatch.stack(
                [(getattr(draw, attribute) or {}).get(label) for draw in self.simulations]
            )
        return stacked

    # -- selection and joining -----------------------------------------------

    @property
    def usable(self) -> SimulationBatch:
        """The draws that did not fail — what an SBI trainer is given.

        ``offset`` is dropped to zero, because the selection no longer occupies
        a contiguous span of the budget and a misleading offset is worse than
        none.
        """
        return SimulationBatch(
            tuple(draw for draw in self.simulations if not draw.failed),
            stream=self.stream,
            provenance=self.provenance,
        )

    def iter_chunks(self, chunk_size: int) -> Iterator[SimulationBatch]:
        """Re-chunk an in-memory batch, for a consumer that wants it in pieces."""
        for start, stop in chunk_bounds(len(self), chunk_size):
            yield SimulationBatch(
                self.simulations[start:stop],
                stream=self.stream,
                offset=self.offset + start,
                provenance=self.provenance,
            )

    @classmethod
    def concatenate(cls, batches: Iterable[SimulationBatch]) -> SimulationBatch:
        """Join chunks back into one batch, in order.

        The provenance is the **union** over the chunks, with disagreement
        resolved conservatively: a budget in which any chunk fell back to the
        loop is a budget whose predictions did not all come from the native
        path, and saying otherwise in a training set's attributes would be a
        false record of how the file was made.
        """
        collected = list(batches)
        if not collected:
            return cls(())
        draws: list[Simulation] = []
        for batch in collected:
            draws.extend(batch.simulations)
        return cls(
            tuple(draws),
            stream=collected[0].stream,
            offset=collected[0].offset,
            provenance=_merge_provenance(collected),
        )

    def __repr__(self) -> str:
        failed = int(np.count_nonzero(self.failed)) if self.simulations else 0
        return f"<SimulationBatch {len(self)} draw(s), {failed} failed, stream={self.stream!r}>"


def _merge_provenance(batches: Sequence[SimulationBatch]) -> dict[str, Any]:
    """Join chunk provenance conservatively; see :meth:`SimulationBatch.concatenate`."""
    merged: dict[str, Any] = {}
    for batch in batches:
        for key, value in batch.provenance.items():
            if key not in merged:
                merged[key] = value
            elif merged[key] != value:
                # Two chunks disagreeing means the budget is a mixture, and the
                # only honest single value for "how was this made?" is the one
                # that claims least: False for a flag, "mixed" for a name.
                merged[key] = False if isinstance(value, bool) else "mixed"
    return merged
