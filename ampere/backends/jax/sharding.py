"""Spreading one chunk's vectorised evaluation across this host's devices.

**Peter's addendum of 2026-09-09**, landed at W3.1 slice 2 as an API with a
single-device implementation that CPU CI exercises and a multi-device one that
is designed, documented and smoke-tested at the API level; the full exercise
waits for the GPU item, which is where hardware to exercise it on appears.

Where this sits
---------------
There are three axes a simulation budget can be spread along, and conflating
them is how a design ends up with one knob that means three things:

* **across processes or machines** — that is
  :class:`~ampere.core.simulate.Executor`, W3.1 slice 1's protocol, and the
  right answer for a slow external simulator;
* **across the draws in one chunk** — that is ``vmap``, one device, and what
  :meth:`~ampere.backends.jax.problem.LoweredProblem.simulate_batched` does;
* **across the accelerators of one host** — that is this module.

The second and third compose: a chunk is ``vmap``ped, and a sharder decides on
how many devices that ``vmap`` runs. The first is orthogonal to both.

What a sharder may not do
-------------------------
**Change the answer.** ``shard(fn, stacked)`` must return what ``fn(stacked)``
returns; where it computed it is not part of the contract. That is why
:class:`SingleDeviceSharder` is not a placeholder but the *reference*
implementation — it is literally ``fn(stacked)``, and a multi-device sharder is
correct exactly insofar as it reproduces it.

It is also deliberately **not a distributed runtime**. There is no scheduler,
no fault tolerance and no cross-host communication here; jax already has all
three, and a wrapper that pretended to own them would be a wrapper nobody could
debug.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import Any

import jax
import jax.numpy as jnp

from ampere.core.exceptions import LoweringError

from ._config import BACKEND

__all__ = ["MeshSharder", "SingleDeviceSharder", "available_devices"]


def available_devices(platform: str | None = None) -> tuple[Any, ...]:
    """The devices jax reports, optionally for one platform only.

    Never used to *choose* a device — ``architecture.md`` §5's "never
    auto-detected" applies to placement as much as to precision — only to
    answer "what is there?" for a sharder that has been asked to use everything.
    """
    if platform is None:
        return tuple(jax.devices())
    try:
        return tuple(jax.devices(platform))
    except RuntimeError:  # pragma: no cover - platform not present
        return ()


class SingleDeviceSharder:
    """One device, and therefore ``fn(stacked)``.

    The default, the degenerate case, and the reference implementation all at
    once. It exists as an object rather than as ``None`` so that a caller can
    *name* the single-device case — a benchmark comparing one device against
    four should not have to compare "a sharder" against "no sharder" — and so
    that the CPU-only rows in CI exercise the same code path a GPU run takes.

    Parameters
    ----------
    device
        Where to place the chunk before evaluating it, or ``None`` to leave
        placement to jax (which is what every other part of this backend does
        unless told otherwise). A platform name or a ``jax.Device``.
    """

    def __init__(self, device: Any = None) -> None:
        self._device = _resolve(device)

    def devices(self) -> tuple[str, ...]:
        """The one device this sharder uses, named."""
        if self._device is None:
            return (str(jax.devices()[0]),)
        return (str(self._device),)

    def shard(self, fn: Callable[[Any], Any], stacked: Any) -> Any:
        """``fn(stacked)``, on this sharder's device."""
        if self._device is None:
            return fn(stacked)
        return fn(jax.device_put(stacked, self._device))

    def __repr__(self) -> str:
        return f"<SingleDeviceSharder devices={self.devices()}>"


class MeshSharder:
    """``jax.pmap`` over several devices, one slice of the chunk each.

    The multi-device implementation, written out so the API is real rather
    than sketched, and **not exercised on hardware here**: CI has one CPU
    device, on which this reduces to the single-device case, and the rows in
    ``tests/gpu`` that would exercise it skip without an accelerator. That is
    the state Peter's addendum asks for — designed, documented, smoke-tested —
    with the full exercise scheduled with the GPU item.

    ``pmap`` rather than ``shard_map`` because the shape here is exactly
    ``pmap``'s: one leading axis, no cross-device communication inside the
    mapped function, identical work per slice. ``shard_map`` earns its extra
    ceremony when a *single* simulation must be split across devices, which
    ampere deliberately does not do — a model too large for one device is a
    :class:`~ampere.core.transform.Model` that places itself
    (``DEVELOPMENT_PLAN.md`` §2, *Batched simulation*).

    **The chunk is padded, not rejected.** ``pmap`` needs the leading axis to
    be a multiple of the device count, and a simulation budget has no reason to
    be. The last few rows are repeated to fill the mesh and dropped from the
    result — wasted work bounded by ``len(devices) - 1`` draws per chunk, which
    is cheaper than either refusing the budget or making the caller round it.

    Parameters
    ----------
    devices
        The devices to spread over. ``None`` takes everything jax reports.

    Raises
    ------
    LoweringError
        If no device is available, which is not a state a chunk can be
        evaluated in.
    """

    def __init__(self, devices: Sequence[Any] | None = None) -> None:
        resolved = tuple(devices) if devices is not None else available_devices()
        if not resolved:
            raise LoweringError(
                "sharding",
                backend=BACKEND,
                detail=(
                    "a chunk sharder needs at least one device and jax reports none. There is "
                    "no meaningful 'no device' case: a chunk has to be evaluated somewhere."
                ),
            )
        self._devices = resolved

    def devices(self) -> tuple[str, ...]:
        """The devices this sharder spreads a chunk over, in order."""
        return tuple(str(device) for device in self._devices)

    def shard(self, fn: Callable[[Any], Any], stacked: Any) -> Any:
        """Split the leading axis over the devices, evaluate, reassemble.

        Reduces to ``fn(stacked)`` on one device, which is what makes the
        single-device case a *path through this code* rather than a branch
        around it — the degenerate case CI exercises is the same arithmetic a
        multi-device run does.
        """
        count = len(self._devices)
        rows = int(jnp.shape(stacked)[0])
        if count == 1 or rows == 0:
            return fn(jax.device_put(stacked, self._devices[0]))
        padded = (-rows) % count
        if padded:
            stacked = jnp.concatenate([stacked, jnp.repeat(stacked[-1:], padded, axis=0)])
        per_device = (rows + padded) // count
        reshaped = stacked.reshape((count, per_device, *jnp.shape(stacked)[1:]))
        # *fn* is already vectorised over the leading axis (it is the ``vmap``
        # ``simulate_batched`` built), so each device receives a slice and maps
        # over it; ``pmap`` supplies the outer axis and nothing else.
        produced = jax.pmap(fn, devices=list(self._devices))(reshaped)
        return jax.tree.map(lambda leaf: _flatten(leaf, rows), produced)

    def __repr__(self) -> str:
        return f"<MeshSharder devices={self.devices()}>"


def _flatten(leaf: Any, rows: int) -> Any:
    """Undo the ``(device, per_device, ...)`` split and drop the padding."""
    shape = jnp.shape(leaf)
    return leaf.reshape((shape[0] * shape[1], *shape[2:]))[:rows]


def _resolve(device: Any) -> Any:
    """A platform name or a ``jax.Device``, or ``None`` for "leave it to jax"."""
    if device is None or not isinstance(device, str):
        return device
    found = available_devices(device)
    if not found:
        raise LoweringError(
            "sharding",
            backend=BACKEND,
            detail=(
                f"jax reports no device on platform {device!r}. Available: "
                f"{[str(one) for one in available_devices()]}."
            ),
        )
    return found[0]
