"""The x64 activation and the construction-time guard (``lowering.md`` §10.2).

``architecture.md`` §5 makes float64 the policy for every likelihood and GP
solve on every backend, always. On jax that policy has an awkward mechanism:
``jax_enable_x64`` is a **process-global** flag which jax's own documentation
addresses to the *application* author — "a global setting that should have one
value for your whole program, set at the top of your main file" — and never to
a library.

So ampere does not set it. ``lowering.md`` §10.2(a) (ruled 2026-09-01,
``DEVELOPMENT_PLAN.md`` §2's "jax x64 activation" row) is explicit that
flipping a process-global numerical switch as an import side effect would
change the behaviour of every other jax user in the interpreter from an
``import`` statement they may not have written. Instead:

* :func:`configure_x64` is the **explicit, idempotent** activation — call it
  before any jax work, from an application, a notebook or a test fixture. It is
  deliberately *not* named ``enable_x64``, because ``jax.enable_x64`` is a
  context manager with different (thread-local, scoped) semantics and two
  same-named functions doing different things across a library boundary is a
  support burden nobody needs (§10.2(b));
* :func:`require_x64` is the **guard**: every jax-backed model, instrument
  step, kernel, GP solver and lowered parameter set calls it at construction,
  and it **raises** when the flag is off, naming the three remedies §10.2(c)
  lists. Raising rather than warning is the load-bearing choice — a warning in
  a notebook scrolls away, and the failure it precedes is a plausible-looking
  wrong answer, not a crash.

numpyro needs no separate switch: ``numpyro.enable_x64`` is a thin wrapper over
``jax.config.update("jax_enable_x64", ...)`` with no dtype state of its own
(verified at the freeze, ``lowering.md`` §12.5), so the one flag governs the
numpyro draws too and this one guard covers both paths.
"""

from __future__ import annotations

import logging

import jax

__all__ = ["BACKEND", "configure_x64", "require_x64", "x64_enabled"]

#: This backend's one name everywhere (W2.12): the ``BACKEND`` capability flag
#: every model, step and noise model here declares, the key
#: ``ampere.core.lowering`` is consulted with, the ``backend`` a composed
#: :class:`~ampere.core.dataset.FittingProblem` reports, the ``ampere_backend``
#: a run records, and the ``name`` of the conformance fixture. One string.
BACKEND = "jax"

_LOGGER = logging.getLogger("ampere.backends.jax")

#: Whether :func:`configure_x64` has already logged. The flag itself is
#: idempotent; the log line is, too, so a library calling it defensively in
#: several places does not produce several lines.
_ANNOUNCED = False

_REMEDIES = (
    "Three ways to fix this, in the order jax's own documentation prefers:\n"
    "  1. set JAX_ENABLE_X64=1 in the environment before anything imports jax "
    "(jax reads it at its own import, so this is the only route guaranteed to "
    "precede array creation);\n"
    "  2. call ampere.backends.jax.configure_x64() before creating any jax "
    "arrays;\n"
    "  3. opt explicitly into reduced precision for this run and accept the "
    "science caveat, which is recorded in the run's provenance."
)


def x64_enabled() -> bool:
    """Whether jax is currently configured for 64-bit arrays.

    Read through ``jax.config.jax_enable_x64`` rather than remembered, because
    the flag is process-global and the host application may set (or unset) it
    at any time; a cached answer would be a guess about someone else's state.
    """
    return bool(jax.config.jax_enable_x64)


def configure_x64() -> None:
    """Turn on jax's 64-bit mode, process-wide. Idempotent; safe to call twice.

    ``lowering.md`` §10.2(b). This is the **process-global** setter
    (``jax.config.update``), not the thread-local ``jax.enable_x64`` context
    manager: the scoped form replaced a ``jax.experimental`` predecessor that
    carried an explicit "fundamentally broken … particularly when used in
    conjunction with JAX transformations" warning, and that warning's absence
    from the promoted API is a documentation gap rather than evidence the
    hazard is gone. Ampere does not stake a numerical guarantee on it.

    Call this **before any jax work** — that is the only condition under which
    it is safe. What jax does with arrays created *before* the flag flips is
    undocumented (``lowering.md`` §10.2(d)); nothing here relies on it.

    Notes
    -----
    Never called at import time by anything in ``ampere``. See this module's
    docstring for why that is a rule and not an oversight.
    """
    global _ANNOUNCED
    if x64_enabled():
        return
    jax.config.update("jax_enable_x64", True)
    if not _ANNOUNCED:
        _LOGGER.info(
            "ampere.backends.jax: enabled jax_enable_x64 (process-global). "
            "float64 is ampere's policy for likelihood and GP linear algebra "
            "(architecture.md §5)."
        )
        _ANNOUNCED = True


def require_x64(what: str) -> None:
    """Refuse to build *what* while jax is in 32-bit mode.

    ``lowering.md`` §10.2(c), and conformance §11 row 11. Called from the
    constructor of every jax-backed model, instrument step, kernel, GP solver
    and lowered parameter set.

    Parameters
    ----------
    what
        What was being constructed, named in the message — "a jax BlackBody",
        "the jax dense GP solver".

    Raises
    ------
    RuntimeError
        With a message naming the three remedies. A plain ``RuntimeError``
        rather than an ampere contract error on purpose: nothing about the
        user's *declaration* is wrong, and nothing about this backend's
        *capabilities* is missing — the interpreter is simply in a
        configuration in which ampere declines to compute, which is a runtime
        environment fault.
    """
    if x64_enabled():
        return
    raise RuntimeError(
        f"cannot build {what}: jax is running in 32-bit mode (jax_enable_x64 is off), and "
        f"ampere runs every likelihood and GP solve in float64 (architecture.md §5) because "
        f"a GP solve in float32 fails in ways that read as science bugs rather than as "
        f"numerical ones. ampere deliberately does not flip this process-global flag for you "
        f"(lowering.md §10.2(a)).\n" + _REMEDIES
    )
