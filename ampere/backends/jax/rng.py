"""jax's realisation of ``ampere.core.rng``'s seed policy (``lowering.md`` §9).

One integer seed per run; named sub-streams derived from it by the *shared*
pure function :func:`ampere.core.rng.substream`; and then each backend does the
idiomatic thing with the integer. Here that is :func:`jax.random.key` — the
**typed** key, not the legacy ``PRNGKey``, which is untyped ``uint32``, carries
an extra trailing axis and no RNG-implementation information.

Two rules from §9.2 are structural rather than stylistic, and both are enforced
by this module's shape:

* **keys are threaded, never stored.** There is no key held on a module field
  anywhere in this package. A key on an :class:`equinox.Module` field makes the
  module's identity depend on RNG state, and a stale key silently reuses draws.
  Every function here *returns* a key for the caller to thread.
* **``substream`` yields a scalar 32-bit integer**, which is exactly what
  :func:`jax.random.fold_in` requires (and ``ampere.core.rng.STREAM_MODULUS``
  is 2**32 for this reason). :func:`fold` is the per-item stream.

Cross-backend note (§9.3): the same seed gives *different* draws on numpy,
torch and jax — PCG64, Mersenne/Philox and Threefry are different bit
generators. Reproducibility is per-backend, and cross-backend agreement is
asserted statistically, never element-wise.
"""

from __future__ import annotations

import jax
import jax.numpy as jnp

from ampere.core.rng import substream

__all__ = ["fold", "key", "split"]


def key(seed: int, label: str) -> jax.Array:
    """The typed PRNG key of the named sub-stream of *seed*.

    ``jax.random.key(substream(seed, label))`` — the jax half of
    ``lowering.md`` §9.2's "share the derivation, not the mechanism". The
    reference backend spells the same policy ``default_rng(substream(...))``
    and torch spells it ``Generator().manual_seed(substream(...))``.

    Parameters
    ----------
    seed
        The run's single integer seed.
    label
        The concern drawing randomness: ``"prior"``, ``"initialisation"``,
        ``"nuts.sampler"``. Different concerns must not share a stream, or
        adding a diagnostic silently changes a fit's initialisation.
    """
    return jax.random.key(substream(seed, label))


def split(rng: jax.Array, count: int = 2) -> jax.Array:
    """*count* independent keys from *rng*. A thin, typed alias for the idiom.

    Present so that no call site in this package has to remember that
    ``jax.random.split`` requires a **scalar** key and rejects a batched one.
    """
    return jax.random.split(rng, count)


def fold(rng: jax.Array, item: int) -> jax.Array:
    """The per-item stream of *rng*: ``jax.random.fold_in(rng, item)``.

    ``fold_in`` requires a scalar 32-bit integer, which is what
    :func:`~ampere.core.rng.substream` produces and what this function coerces
    *item* to, so a Python ``int``, a numpy integer and a 0-d array all behave
    the same way.
    """
    return jax.random.fold_in(rng, jnp.asarray(item, dtype=jnp.uint32))
