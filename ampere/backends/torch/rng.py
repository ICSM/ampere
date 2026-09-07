"""``lowering.md`` §9 for torch: named sub-streams, and torch's generator gap.

``ampere.core.rng.substream(seed, label)`` is the shared derivation — one
integer seed per run, a pure function of ``(seed, label)`` over a stable
``blake2b`` digest — and each backend does the idiomatic thing with the integer
it returns. Here that is ``torch.Generator().manual_seed(...)``. Sharing the
*derivation* and not the mechanism is what makes "the prior-sampling stream"
mean the same thing on all three backends, and what stops adding a diagnostic
silently changing a fit's initialisation.

torch's generator gap, and the two routes round it
---------------------------------------------------
``torch.distributions.Distribution.sample`` and ``.rsample`` take
``sample_shape`` and **nothing else**: there is no ``generator=`` argument, and
the underlying ``torch.normal``/``torch.rand``/``_standard_normal`` calls are
made against the global RNG without one being threaded through. So handing a
``Generator`` to a distribution — the obvious thing to try — silently does not
seed the draw. The draws are perfectly valid; they are just not the seeded
ones, which is precisely the class of failure that is discovered months later.

``lowering.md`` §9.1 names two working routes and ranks them:

1. :func:`torch.random.fork_rng` around ``torch.manual_seed(s)`` and the
   sampling call, which restores the global state on exit so it does not leak
   into a host process. Works for any distribution.
2. Draw uniforms explicitly and push them through ``icdf``:
   ``torch.rand(shape, generator=g)`` *does* accept a generator. This is a
   genuinely ``Generator``-scoped stream with no global state touched at all.

**Route (2) where ``icdf`` exists, route (1) otherwise.** Route (2) is the only
one that is safe under threading and the only one that never mutates
process-global state — the same objection ``lowering.md`` §10.1 raises against
``torch.set_default_dtype``. Its limitation is exactly the ``icdf`` gap of
§3.6, so the two routes are complementary rather than competing and a backend
needs both. :meth:`~ampere.backends.torch.parameters.TorchParameterSpace.sample`
is where both are used.

What reproducibility means here
--------------------------------
Same seed, same backend, same library versions gives identical draws, and that
is testable. Same seed *across* backends does **not**, and ampere must not
claim otherwise: numpy's PCG64, torch's Mersenne/Philox and jax's Threefry are
different bit generators, and even where the algorithm matched, the order in
which draws are consumed would differ. Cross-backend agreement is asserted
statistically (``lowering.md`` §9.3).

Examples
--------
>>> from ampere.backends.torch.rng import generator
>>> import torch
>>> a = torch.rand(3, generator=generator(20260907, "prior"), dtype=torch.float64)
>>> b = torch.rand(3, generator=generator(20260907, "prior"), dtype=torch.float64)
>>> bool(torch.equal(a, b))
True
>>> c = torch.rand(3, generator=generator(20260907, "initialisation"), dtype=torch.float64)
>>> bool(torch.equal(a, c))
False
"""

from __future__ import annotations

import secrets

import torch

from ampere.core.rng import substream

from ._config import DEFAULT_DEVICE

__all__ = ["generator", "seed_for"]


def seed_for(seed: int | None, label: str) -> int:
    """The integer this backend seeds *label*'s stream with.

    ``seed=None`` means the run did not ask to be reproducible, which is an
    honest state rather than an error: a fresh, entropy-derived integer is
    returned and nothing about it is recorded as repeatable. Any other seed
    goes through :func:`ampere.core.rng.substream`, so the same
    ``(seed, label)`` gives the same stream in this process, the next one, and
    on another machine.
    """
    if seed is None:
        return secrets.randbelow(2**32)
    return substream(int(seed), label)


def generator(
    seed: int | None, label: str, *, device: torch.device = DEFAULT_DEVICE
) -> torch.Generator:
    """A ``torch.Generator`` for one named sub-stream of *seed*.

    Per-device, as torch's generators are: a generator made for the CPU cannot
    seed a CUDA draw, so the device is part of the request rather than
    something inferred later.
    """
    stream = torch.Generator(device=device)
    stream.manual_seed(seed_for(seed, label))
    return stream
