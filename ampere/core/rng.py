"""Seed derivation: one seed per run, named sub-streams off it.

``DEVELOPMENT_PLAN.md`` §4.5 asks the inference contracts for a "named seed
handling that lowers to each backend's model (numpy Generators, torch
Generators, jax PRNG keys), so runs are reproducible across backends", and
``docs/design/lowering.md`` §9.2 is the design that answers it. This module is
that design, implemented: a single pure function of ``(seed, label)`` that every
backend turns into its own idiomatic generator. Its home here was **ratified
at the freeze** (ruled 2026-09-03, ``lowering.md`` §12.7 = ``inference.md``
§19.5): pure stdlib+numpy, deliberately free-standing, and the alternative
was three backends agreeing by convention.

Why a shared derivation rather than a shared mechanism
------------------------------------------------------
The three RNG models are irreconcilable — numpy's :class:`numpy.random.Generator`
is stateful and accepted directly by ``scipy``'s ``rvs``; torch's is stateful,
per-device, and *not* accepted by ``torch.distributions``; jax's keys are values
that must be threaded explicitly. What can be shared is the *number*: each
backend does the idiomatic thing with the integer :func:`substream` returns
(``default_rng(n)``, ``Generator().manual_seed(n)``, ``jax.random.key(n)``), and
"the prior-sampling stream" then means the same thing on all three.

Why different concerns must not share a stream
-----------------------------------------------
Prior sampling, sampler initialisation, an SBI simulation budget and a
posterior-predictive draw are separate uses. Drawing them from one stream makes
each silently dependent on how many draws the others took, so adding a
diagnostic changes a fit's initialisation. Naming the stream removes that
coupling by construction.

Why ``blake2b`` and not :func:`hash`
------------------------------------
Python's built-in :func:`hash` is salted per process (``PYTHONHASHSEED``), so a
label-derived seed built on it would differ between two invocations of the same
script — an irreproducibility that is real, easily missed, and would be blamed
on the sampler. :mod:`hashlib` digests are stable across processes, platforms
and Python versions.

Placement
---------
``lowering.md`` §9.2 puts this "once in ``ampere.core``"; its §12.7 asks W1.13 to
ratify that placement rather than let it be assumed, since it is a small
addition to a frozen contract's surface. The module is deliberately tiny and
free-standing so that ratifying it, or moving it, is a one-line change.

Examples
--------
>>> substream(20260902, "prior") == substream(20260902, "prior")
True
>>> substream(20260902, "prior") == substream(20260902, "simulate")
False
>>> 0 <= substream(20260902, "prior") < 2**32
True
>>> generator(20260902, "prior").random() == generator(20260902, "prior").random()
True
"""

from __future__ import annotations

import hashlib

import numpy as np

__all__ = ["SEED_BYTES", "STREAM_MODULUS", "generator", "substream"]

#: Width of the derived integer. ``jax.random.fold_in`` requires a **scalar
#: 32-bit** integer (``lowering.md`` §9.2), which is the narrowest of the three
#: backends' requirements and therefore the one that fixes this.
STREAM_MODULUS = 2**32

#: Width the run seed is serialised to before hashing. Signed, so a negative
#: seed is derived from rather than rejected.
SEED_BYTES = 8

#: Personalisation string, so an ampere sub-stream digest can never collide with
#: some other blake2b use of the same label elsewhere in a pipeline.
_PERSON = b"ampere-rng"


def substream(seed: int, label: str) -> int:
    """Derive the integer seed of the named sub-stream of *seed*.

    Pure, stable across processes and platforms, and defined for every string
    label. Two different labels give unrelated streams; the same ``(seed,
    label)`` always gives the same integer.

    Parameters
    ----------
    seed
        The run's single integer seed, recorded in provenance
        (``lowering.md`` §9.2). May be negative.
    label
        Name of the concern drawing randomness — ``"prior"``,
        ``"initialisation"``, ``"simulate"``, ``"posterior_predictive"``. A
        dotted label (``"simulate.spectrum"``) is an ordinary string here; this
        function imposes no structure on it.

    Returns
    -------
    int
        In ``[0, 2**32)``.

    Raises
    ------
    TypeError
        If *seed* is not an integer or *label* is not a string. Both are
        :class:`TypeError` rather than an ampere contract error because they are
        argument-type mistakes in a pure function, not misuse of a contract.

    Examples
    --------
    >>> substream(0, "prior")
    3058495773
    >>> substream(1, "prior") != substream(0, "prior")
    True
    """
    if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
        raise TypeError(f"a run seed must be an integer, got {seed!r}")
    if not isinstance(label, str):
        raise TypeError(f"a sub-stream label must be a string, got {label!r}")
    payload = int(seed).to_bytes(SEED_BYTES, "big", signed=True) + label.encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=4, person=_PERSON).digest()
    return int.from_bytes(digest, "big")


def generator(seed: int, label: str) -> np.random.Generator:
    """The reference backend's realisation of :func:`substream`.

    ``numpy.random.default_rng(substream(seed, label))`` — spelled once so the
    reference path and the torch/jax lowerings (``lowering.md`` §9.2) are
    visibly the same policy with different mechanisms.

    Examples
    --------
    >>> rng = generator(20260902, "simulate")
    >>> float(rng.normal()) == float(generator(20260902, "simulate").normal())
    True
    """
    return np.random.default_rng(substream(seed, label))
