"""Spreading one chunk's vectorised evaluation across this host's devices.

**Peter's addendum of 2026-09-09**, landed at W3.1 slice 2 as an API with a
single-device implementation that CPU CI exercises and a distributed one that
is designed, documented and refused-by-name without a process group; the full
exercise waits for the GPU item.

The twin of :mod:`ampere.backends.jax.sharding`, and deliberately the same
shape, but the mechanism underneath is different in a way worth stating rather
than papering over. jax's ``pmap`` is a *transformation*: one process drives
every local device and the split is a reshape. torch's distributed story is a
**process group** — one process per device, launched by ``torchrun`` or
``mp.spawn``, communicating through a collective backend — and there is no
in-process equivalent that spreads a ``torch.func.vmap`` over several GPUs.

That asymmetry has one honest consequence, and it is the reason this module is
short: **ampere does not launch a process group.** A run that wants several
GPUs is launched under ``torchrun``, and this module then evaluates each rank's
slice of the chunk and gathers the result. Launching workers on a user's behalf
is a job for a job scheduler, and a library that did it would be a library that
owns process lifetimes it cannot see the end of — the same reason
:class:`~ampere.core.simulate.Executor` takes a mapper rather than starting a
cluster.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import torch

from ampere.core.exceptions import LoweringError

from ._config import BACKEND

__all__ = ["DistributedSharder", "SingleDeviceSharder", "available_devices"]


def available_devices() -> tuple[torch.device, ...]:
    """The devices torch reports: every visible accelerator, else the CPU.

    Never used to *choose* a device — ``architecture.md`` §5's "never
    auto-detected" applies to placement as much as to precision — only to
    answer "what is there?" for a sharder asked to use everything.
    """
    if torch.cuda.is_available() and torch.cuda.device_count():
        return tuple(torch.device("cuda", index) for index in range(torch.cuda.device_count()))
    return (torch.device("cpu"),)


class SingleDeviceSharder:
    """One device, and therefore ``fn(stacked)``.

    The default, the degenerate case, and the reference implementation all at
    once. It exists as an object rather than as ``None`` so that a caller can
    *name* the single-device case, and so that the CPU-only rows in CI exercise
    the same code path a GPU run takes.

    Parameters
    ----------
    device
        Where to place the chunk before evaluating it. ``None`` leaves the
        tensors where the realisation put them, which is what
        ``simulate_batched`` does when given no sharder at all.
    """

    def __init__(self, device: Any = None) -> None:
        self._device = None if device is None else torch.device(device)

    def devices(self) -> tuple[str, ...]:
        """The one device this sharder uses, named."""
        return (str(self._device if self._device is not None else available_devices()[0]),)

    def shard(self, fn: Callable[[Any], Any], stacked: Any) -> Any:
        """``fn(stacked)``, on this sharder's device."""
        if self._device is None:
            return fn(stacked)
        return fn(stacked.to(self._device))

    def __repr__(self) -> str:
        return f"<SingleDeviceSharder devices={self.devices()}>"


class DistributedSharder:
    """One rank's slice of the chunk, gathered back into the whole.

    The multi-device implementation, written out so the API is real rather than
    sketched, and **not exercised on hardware here**: it needs an initialised
    :mod:`torch.distributed` process group, which CI has not got, so the rows in
    ``tests/gpu`` that would exercise it skip and the refusal below is what a
    CPU-only run meets. That is the state Peter's addendum asks for.

    **Every rank calls this with the same chunk.** That is the collective
    convention, not an oversight: each rank evaluates rows
    ``[rank::world_size]`` and ``all_gather`` puts the pieces back together, so
    every rank ends with the whole answer and ``simulate_many`` proceeds
    identically everywhere. The alternative — a driver rank scattering work —
    would need ampere to own which process is the driver, and a library that
    decides that has decided how a user's job script is written.

    The rows are taken **strided rather than contiguous** so that every rank
    gets a slice of the same size to within one row, without padding: a
    contiguous split of 10 rows over 4 ranks is 3/3/3/1 and the gather has to
    handle ragged pieces, while the strided one is 3/3/2/2 and each piece is
    ``ceil`` or ``floor`` of the same number. The reassembly is by index, so the
    result is in the caller's order whichever way it was cut.

    ``DTensor`` is the other route and is deliberately not taken: it shards a
    *tensor* across a mesh so that one operator's work is split, which is the
    model-parallel case ampere leaves to the model (``DEVELOPMENT_PLAN.md`` §2,
    *Batched simulation*: "a single simulation larger than one device is a
    Model whose ``__call__`` does its own placement"). What is being split here
    is a batch of independent simulations, for which a plain strided split and
    a gather is both simpler and exactly as fast.

    Parameters
    ----------
    device
        The device this rank computes on; ``None`` uses this rank's default
        CUDA device where there is one.
    group
        A :mod:`torch.distributed` process group, or ``None`` for the default.

    Raises
    ------
    LoweringError
        If :mod:`torch.distributed` is unavailable or no process group has been
        initialised — by name, with the remedy, rather than a collective that
        hangs.
    """

    def __init__(self, device: Any = None, group: Any = None) -> None:
        distributed = getattr(torch, "distributed", None)
        if distributed is None or not distributed.is_available():
            raise LoweringError(
                "sharding",
                backend=BACKEND,
                detail=(
                    "torch.distributed is not available in this build, so a chunk cannot be "
                    "spread across ranks. Use SingleDeviceSharder, which is what a "
                    "single-device run wants anyway."
                ),
            )
        if not distributed.is_initialized():
            raise LoweringError(
                "sharding",
                backend=BACKEND,
                detail=(
                    "no torch.distributed process group is initialised, so there are no ranks "
                    "to shard a chunk over. ampere does not launch one: start the run under "
                    "torchrun (or call torch.distributed.init_process_group yourself) and "
                    "build this sharder inside each rank. For one device, "
                    "SingleDeviceSharder needs none of that."
                ),
            )
        self._distributed = distributed
        self._group = group
        self._rank = int(distributed.get_rank(group=group))
        self._world = int(distributed.get_world_size(group=group))
        self._device = torch.device(device) if device is not None else _rank_device(self._rank)

    def devices(self) -> tuple[str, ...]:
        """This rank's device. One name per rank, gathered by the caller if wanted."""
        return (str(self._device),)

    @property
    def rank(self) -> int:
        """This process's rank in the group."""
        return self._rank

    @property
    def world_size(self) -> int:
        """How many ranks the chunk is spread over."""
        return self._world

    def shard(self, fn: Callable[[Any], Any], stacked: Any) -> Any:
        """Evaluate this rank's stride of *stacked* and gather the whole answer."""
        rows = int(stacked.shape[0])
        if self._world == 1 or rows == 0:
            return fn(stacked.to(self._device))
        mine = stacked[self._rank :: self._world].to(self._device)
        produced = fn(mine)
        return _gather(self._distributed, produced, rows, self._rank, self._world, self._group)

    def __repr__(self) -> str:
        return f"<DistributedSharder rank={self._rank}/{self._world} device={self._device}>"


def _rank_device(rank: int) -> torch.device:
    devices = available_devices()
    return devices[rank % len(devices)]


def _gather(
    distributed: Any,
    produced: Any,
    rows: int,
    rank: int,
    world: int,
    group: Any,
) -> Any:
    """Put every rank's stride back in the caller's order, leaf by leaf."""

    def recurse(value: Any) -> Any:
        return _gather(distributed, value, rows, rank, world, group)

    if isinstance(produced, dict):
        return {key: recurse(value) for key, value in produced.items()}
    if isinstance(produced, tuple):
        return tuple(recurse(value) for value in produced)
    if isinstance(produced, Sequence) and not isinstance(produced, (str, bytes)):
        return [recurse(value) for value in produced]
    pieces: list[Any] = [None] * world
    distributed.all_gather_object(pieces, produced.cpu(), group=group)
    whole = torch.empty(
        (rows, *tuple(produced.shape[1:])), dtype=produced.dtype, device=produced.device
    )
    for index, piece in enumerate(pieces):
        whole[index::world] = piece.to(produced.device)
    return whole
