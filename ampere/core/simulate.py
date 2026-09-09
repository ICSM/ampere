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
import os
import time
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from typing import TYPE_CHECKING, Any, Protocol, TypeVar, runtime_checkable

import numpy as np

from .exceptions import DatasetError
from .results_schema import FunctionSamples, ModelResult

if TYPE_CHECKING:  # pragma: no cover - typing only
    from .dataset import Failure, Simulation

__all__ = [
    "ChunkHook",
    "ContainerBatch",
    "ExecutionFailure",
    "Executor",
    "ProcessExecutor",
    "SerialExecutor",
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


def _windowed_map(
    pool_factory: Callable[[], concurrent.futures.Executor],
    fn: Callable[[_Item], _Result],
    items: Sequence[_Item],
    *,
    workers: int,
    timeout: float | None,
    recover: bool,
) -> list[Any]:
    """Run *items* on a pool, at most *workers* in flight, with a per-item deadline.

    The window is what makes ``timeout`` mean "per simulation". Submitting the
    whole budget at once and timing each future from submission would expire
    every queued item at the same moment on a budget larger than the pool; with
    at most one item in flight per worker, submission time *is* start time to
    within the scheduler's own latency.

    *recover* says whether a dead pool can be replaced (processes: yes;
    threads: there is nothing to replace and nothing to kill).
    """
    results: list[Any] = [None] * len(items)
    pending = list(range(len(items)))
    while pending:
        pool = pool_factory()
        queue = pending
        pending = []
        inflight: list[tuple[int, concurrent.futures.Future[Any], float]] = []
        cursor = 0
        rebuild = False
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
                    if recover:
                        _rerun_alone(pool_factory, fn, items, results, survivors, timeout)
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
            # ``wait=False`` after a timeout: the expired worker is still busy,
            # and waiting for it is exactly what the timeout said not to do.
            if rebuild:
                _terminate(pool)
            pool.shutdown(wait=not rebuild, cancel_futures=True)
    return results


def _rerun_alone(
    pool_factory: Callable[[], concurrent.futures.Executor],
    fn: Callable[[_Item], _Result],
    items: Sequence[_Item],
    results: list[Any],
    survivors: Sequence[int],
    timeout: float | None,
) -> None:
    """Re-run *survivors* one per fresh pool, so a crash is attributed correctly."""
    for index in survivors:
        pool = pool_factory()
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
        factory = functools.partial(concurrent.futures.ThreadPoolExecutor, max_workers=workers)
        return _windowed_map(
            factory,
            fn,
            materialised,
            workers=workers,
            timeout=self.timeout,
            recover=False,
        )

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
        A :mod:`multiprocessing` context, if the default start method is wrong
        for this simulator. Worth reaching for: the platform default on Linux is
        ``fork``, and forking a process that already holds a threaded runtime —
        jax says so itself, and torch's intra-op pools have the same shape —
        risks a deadlock in the child. ``get_context("spawn")`` is the safe
        answer where that applies, at the cost of re-importing the world in each
        worker.

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
        self.mp_context = mp_context
        self._shared: Any = None
        if self.max_workers is not None and self.max_workers < 1:
            raise DatasetError(f"a process pool needs at least one worker, got {max_workers!r}.")
        if self.timeout is not None and self.timeout <= 0.0:
            raise DatasetError(f"a per-simulation timeout must be positive, got {timeout!r}.")

    def broadcast(self, value: Any) -> None:
        """Send *value* to every worker once, at start-up.

        Retrieved inside the worker with :func:`worker_shared`. ``simulate_many``
        uses this to place the fitting problem, which is why the pickling cost
        of a large problem is paid per worker rather than per draw.
        """
        self._shared = value

    def _factory(self, workers: int) -> concurrent.futures.Executor:
        return concurrent.futures.ProcessPoolExecutor(
            max_workers=workers,
            mp_context=self.mp_context,
            initializer=_install_shared,
            initargs=(self._shared,),
            max_tasks_per_child=self.max_tasks_per_child,
        )

    def map(self, fn: Callable[[_Item], _Result], items: Iterable[_Item], /) -> list[Any]:
        """Apply *fn* to every item on the pool, returning results in item order."""
        materialised = list(items)
        if not materialised:
            return []
        workers = self.max_workers or min(len(materialised), _default_workers())
        return _windowed_map(
            functools.partial(self._factory, workers),
            fn,
            materialised,
            workers=workers,
            timeout=self.timeout,
            recover=True,
        )

    def shutdown(self) -> None:
        """Drop the broadcast payload. Pools are per ``map`` call and already closed."""
        self._shared = None

    def __enter__(self) -> ProcessExecutor:
        return self

    def __exit__(self, *exc_info: object) -> None:
        self.shutdown()

    def __repr__(self) -> str:
        return (
            f"<ProcessExecutor max_workers={self.max_workers}, timeout={self.timeout}, "
            f"max_tasks_per_child={self.max_tasks_per_child}>"
        )


def _default_workers() -> int:
    return max(1, os.cpu_count() or 1)


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
    """

    simulations: tuple[Simulation, ...]
    stream: str = "simulate"
    offset: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "simulations", tuple(self.simulations))

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
            tuple(draw for draw in self.simulations if not draw.failed), stream=self.stream
        )

    def iter_chunks(self, chunk_size: int) -> Iterator[SimulationBatch]:
        """Re-chunk an in-memory batch, for a consumer that wants it in pieces."""
        for start, stop in chunk_bounds(len(self), chunk_size):
            yield SimulationBatch(
                self.simulations[start:stop], stream=self.stream, offset=self.offset + start
            )

    @classmethod
    def concatenate(cls, batches: Iterable[SimulationBatch]) -> SimulationBatch:
        """Join chunks back into one batch, in order."""
        collected = list(batches)
        if not collected:
            return cls(())
        draws: list[Simulation] = []
        for batch in collected:
            draws.extend(batch.simulations)
        return cls(tuple(draws), stream=collected[0].stream, offset=collected[0].offset)

    def __repr__(self) -> str:
        failed = int(np.count_nonzero(self.failed)) if self.simulations else 0
        return f"<SimulationBatch {len(self)} draw(s), {failed} failed, stream={self.stream!r}>"
