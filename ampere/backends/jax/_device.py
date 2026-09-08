"""Where a piece of this backend computes — chosen, never detected.

``architecture.md`` §5's rule, in one place so that every model, instrument
step, kernel, noise model and GP solver here obeys the same one:

* **the default is the CPU**, always, on every machine. A library that noticed
  an accelerator and silently used it would make "why is this run slow?" and
  "why did this run give a different answer?" both unanswerable, and it would
  make CI (which has no accelerator, by ruling) exercise a different code path
  from the machine the science is done on;
* **a device is asked for by name** — ``device="cuda"``, ``device="tpu"`` — or
  handed over as an explicit ``jax.Device`` for the caller who has already
  chosen *which* of several. Anything else is a refusal;
* **a name this process has no device for is a refusal**, naming what is
  present. Never a fallback: a fit that was asked for a GPU and quietly ran on
  the CPU is a fit whose timings mean nothing, and whose provenance is a lie;
* **the choice is per instance, and it is a capability flag.** W2.12 made
  ``DEVICE`` a class-level flag that ``ampere.core.declared_capabilities``
  aggregates by attribute access, "all parts agree or refuse". An instance
  attribute shadows the class default, so a solver placed on a GPU beside CPU
  models is a composition-time refusal rather than a runtime surprise — which
  is the whole reason the flag is aggregated rather than merely recorded.

Why the refusal type varies
---------------------------
:func:`resolve_device` takes the exception class from its caller rather than
picking one, so the error a user sees names the contract the piece belongs to:
``TransformationError`` for a model or an instrument step
(``ampere.core.transform``), ``LikelihoodError`` for a kernel, a noise model or
a GP solver (``ampere.core.likelihood``). The message is identical; only the
class differs, and it differs the way every other refusal in that contract
does.

What placement actually means on jax
------------------------------------
jax runs an operation where its inputs live, so "placing" a piece means
:func:`jax.device_put` on the arrays it *owns* — a model's evaluation grid, an
instrument step's influence matrix, a solver's residuals and covariance. A
piece that owns no arrays (the noise-model compositions, whose numbers all
arrive from the kernel and the solver) carries the flag and nothing else, and
says so in its own docstring rather than pretending to a placement it does not
perform.
"""

from __future__ import annotations

from typing import Any

import jax

from ampere.core.exceptions import AmpereError, LikelihoodError

__all__ = ["DEVICE", "device_flag", "place_on", "resolve_device"]

#: The default device for every part of this backend. **Chosen, never
#: detected** (``architecture.md`` §5).
DEVICE = "cpu"


def resolve_device(
    device: Any,
    owner: str,
    *,
    error: type[AmpereError] = LikelihoodError,
) -> Any:
    """A jax device for *device*, by platform name. Never auto-detected.

    Parameters
    ----------
    device
        A platform name (``"cpu"``, ``"cuda"``, ``"tpu"``, …) or an explicit
        ``jax.Device``, which is passed through untouched — picking *which*
        GPU is the caller's business and ampere has no opinion about it.
    owner
        What is being placed, named in the refusal.
    error
        The exception class to raise; see this module's docstring for why the
        caller chooses it.

    Returns
    -------
    The first device jax reports for that platform.

    Raises
    ------
    AmpereError
        (of class *error*) when this process has no device of that platform,
        naming the platforms it does have.
    """
    if not isinstance(device, str):
        return device  # an explicit jax.Device, used as given
    available = jax.devices()
    for found in available:
        if found.platform == device:
            return found
    platforms = ", ".join(sorted({d.platform for d in available}))
    raise error(
        f"{owner} was asked for device {device!r}, but this jax process has no such platform. "
        f"Available: {platforms}. ampere never falls back to another device — a fit that was "
        f"asked for a GPU and quietly ran on the CPU is a fit whose timings mean nothing."
    )


def device_flag(device: Any, resolved: Any) -> str:
    """The ``DEVICE`` capability flag for a piece placed on *resolved*.

    The *asked-for* string when a name was given, so the flag reads back as the
    user wrote it; the resolved device's platform when a ``jax.Device`` was
    handed over, since that is the only name it has.
    """
    return device if isinstance(device, str) else str(resolved.platform)


def place_on(array: Any, resolved: Any) -> jax.Array:
    """*array* on *resolved*, unchanged in dtype.

    A thin wrapper so that "this backend places its buffers with
    ``device_put``" is one call spelled one way, greppable, rather than four
    modules' worth of the same idiom drifting apart.
    """
    return jax.device_put(array, resolved)
