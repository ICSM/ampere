"""dtype, device and the one backend name — ``lowering.md`` §10.1's policy, once.

Three facts every other module here imports rather than restates.

**The backend's name is ``"torch"``, and it is one string everywhere** (W2.12):
the ``BACKEND`` capability flag every model, instrument step, noise model and
solver in this package declares, the key its rows are registered under in
``ampere.core.lowering``, the ``name`` of its conformance fixture, and the
``ampere_backend`` a run records. There is no translation table.

**float64, threaded explicitly, never set globally.** ``architecture.md`` §5
makes double precision the policy for likelihood and GP linear algebra on
every backend, and ``lowering.md`` §10.1 fixes the mechanism: every tensor this
package creates takes an explicit ``dtype=`` and ``device=``, and
:func:`torch.set_default_dtype` is **never** called. It is global mutable
process state, so a library that sets it changes the numerical behaviour of
every other consumer of torch in the same interpreter — the same objection
``lowering.md`` §10.2(a) raises against ampere flipping jax's x64 flag on
import.

The consequence worth stating is the one about distributions.
``torch.distributions`` normalises its arguments through ``broadcast_all``,
which promotes bare Python floats using ``torch.get_default_dtype()`` —
float32 — but switches to the dtype *and* device of the first ``Tensor``
argument it finds. So ``Normal(0.0, 1.0)`` is float32 while
``Normal(tensor(0.0, dtype=float64), 1.0)`` is float64 throughout. Passing
every argument through :func:`as_tensor` is what makes a lowered prior's
precision a property of this package rather than of whatever the host process
last did to torch's global default.

**CPU by default; a device is chosen, never detected** (``architecture.md``
§5). :data:`DEFAULT_DEVICE` is CPU, and a machine with a GPU present takes the
same code path as CI.

W2.4 slice 3 makes that choice **per instance**. Every model, instrument step,
kernel, noise model and solver in this package takes a ``device=`` keyword
threaded exactly as ``dtype`` already was, and reports it back as its own
``DEVICE`` capability flag — an *instance* attribute shadowing the class
default, which is all ``ampere.core.declared_capabilities`` needs, since it
reads ``part.DEVICE`` by attribute access. The aggregation rule is the one
that contract already states: every part of a problem must name the same
device or composition is refused, so a whole problem is composed on one device
and nothing is moved on the user's behalf.

The string is :func:`device_name`'s, which is ``str(torch.device(...))``:
``"cpu"``, ``"cuda"``, ``"cuda:1"``. ``"cuda"`` and ``"cuda:0"`` are therefore
*different* declarations, and mixing them in one problem is refused rather than
silently reconciled — torch's own two spellings mean different things the
moment a second GPU exists, and guessing which was meant is exactly the
"silent performance collapse" ``declared_capabilities`` refuses to risk.

:func:`place` and :func:`move` are the two verbs, written once here so that
"a piece knows where it lives, and can be moved" is not re-implemented five
times.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import torch

__all__ = [
    "BACKEND",
    "DEFAULT_DEVICE",
    "DEFAULT_DTYPE",
    "as_tensor",
    "complex_dtype",
    "device_name",
    "move",
    "place",
    "resolve_device",
    "resolve_dtype",
    "to_numpy",
]

#: This backend's one name, in the sense W2.12 fixed. See the module docstring.
BACKEND = "torch"

#: The dtype every tensor this package creates is built at, unless a caller
#: asks for another. ``architecture.md`` §5's float64 policy.
DEFAULT_DTYPE: torch.dtype = torch.float64

#: Where those tensors live. CPU, chosen rather than detected.
DEFAULT_DEVICE: torch.device = torch.device("cpu")


def as_tensor(
    value: Any,
    *,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> torch.Tensor:
    """*value* as a tensor of exactly *dtype* on exactly *device*.

    The single construction point for this package, so that "float64 unless
    someone said otherwise" is enforced in one place rather than remembered at
    every call site. A tensor that already has the requested dtype and device
    is returned unchanged (``torch.as_tensor``'s own behaviour), so this is
    free on the hot path; anything else — a Python float, a numpy array, a
    tensor of another precision — is converted.

    Gradient tracking survives: ``torch.as_tensor`` on a tensor that already
    matches is the identity, and ``Tensor.to`` is a differentiable operation,
    so passing a leaf through here does not detach it.
    """
    if isinstance(value, torch.Tensor):
        return value.to(dtype=dtype, device=device)
    array = np.asarray(value)
    if array.dtype != object and not array.flags.writeable:
        # ``ampere.core``'s Buffer hands out read-only views, deliberately, so
        # that a declared constant cannot be mutated through the array a caller
        # was given. torch cannot wrap one without warning that writing to the
        # resulting tensor is undefined behaviour, so it is copied here rather
        # than the warning being silenced: a copy of a buffer is exactly what a
        # tensor that torch may move between devices ought to be.
        array = array.copy()
    return torch.as_tensor(array, dtype=dtype, device=device)


def to_numpy(value: Any) -> np.ndarray:
    """A tensor (or anything array-like) back on the numpy side of the boundary.

    ``ampere.core``'s containers, oracles and provenance all speak numpy, so
    every value this backend hands back across that boundary goes through
    here. ``detach`` is deliberate and is the reason this function exists at
    all: a tensor carrying a graph cannot be converted, and the failure
    (``RuntimeError: Can't call numpy() on Tensor that requires grad``) would
    otherwise surface far from the conversion that caused it.
    """
    if isinstance(value, torch.Tensor):
        return value.detach().cpu().numpy()
    return np.asarray(value)


def resolve_device(device: Any) -> torch.device:
    """*device* as a :class:`torch.device`, without ever detecting one.

    Accepts what ``torch.device`` accepts — a string (``"cpu"``, ``"cuda"``,
    ``"cuda:1"``), a ``torch.device``, or an integer CUDA ordinal — and nothing
    else. ``None`` is **not** accepted and does not mean "pick one":
    ``architecture.md`` §5 makes the device an explicit choice, so the absence
    of a choice is :data:`DEFAULT_DEVICE`, decided at the call site by the
    keyword's default rather than here by a probe of the machine.
    """
    if isinstance(device, torch.device):
        return device
    if isinstance(device, str):
        return torch.device(device)
    if isinstance(device, int) and not isinstance(device, bool):
        return torch.device("cuda", device)
    raise TypeError(
        f"device= takes a torch.device, a device string ('cpu', 'cuda', 'cuda:1') or a CUDA "
        f"ordinal, got {device!r}. ampere never detects a device for you (architecture.md §5): "
        f"leave it out for the CPU, or name the one you mean."
    )


def device_name(device: Any) -> str:
    """The ``DEVICE`` capability string for *device*.

    ``str(torch.device(...))`` exactly, so the flag says what torch says. See
    the module docstring for why ``"cuda"`` and ``"cuda:0"`` stay distinct.
    """
    return str(resolve_device(device))


def resolve_dtype(dtype: Any) -> torch.dtype:
    """*dtype* as a floating-point :class:`torch.dtype`, by object or by name."""
    resolved = getattr(torch, dtype, None) if isinstance(dtype, str) else dtype
    if not isinstance(resolved, torch.dtype):
        raise TypeError(f"dtype= takes a torch.dtype or the name of one, got {dtype!r}.")
    return resolved


#: The complex dtype paired with each real one. A circular complex Gaussian's
#: residual is complex and its variance is real, so the two travel together
#: through the realised path and neither may be guessed from the other at the
#: point of use (``likelihoods.md`` §4).
_COMPLEX_OF: dict[torch.dtype, torch.dtype] = {
    torch.float32: torch.complex64,
    torch.float64: torch.complex128,
}


def complex_dtype(dtype: torch.dtype) -> torch.dtype:
    """The complex dtype whose components are *dtype*.

    float64 pairs with complex128, float32 with complex64. A dtype with no
    complex partner is an error rather than a silent promotion: it would mean
    a complex family had been asked to compute in a precision torch cannot
    represent its own residuals in.
    """
    found = _COMPLEX_OF.get(dtype)
    if found is None:
        raise TypeError(
            f"there is no complex dtype whose components are {dtype}; a complex likelihood "
            f"needs float32 or float64 arithmetic."
        )
    return found


def place(owner: Any, dtype: Any, device: Any) -> None:
    """Record *owner*'s precision and device, and declare the latter.

    Three attributes, set together because they are one fact: ``dtype`` and
    ``device`` are what every tensor *owner* builds is built with, and
    ``DEVICE`` is the capability flag ``ampere.core.declared_capabilities``
    reads. Setting it here makes it an **instance** attribute shadowing the
    class-level default W2.12 gave every piece — which is the whole mechanism
    by which a device becomes per-instance without a line changing in
    ``ampere.core``.
    """
    owner.dtype = resolve_dtype(dtype)
    owner.device = resolve_device(device)
    owner.DEVICE = device_name(owner.device)


def move(owner: Any, *, dtype: Any = None, device: Any = None) -> Any:
    """Move *owner*'s registered tensors, and re-declare where they live.

    ``architecture.md`` §5's "buffers move with parameters under ``.to(...)``",
    discharged in one place: *owner* keeps its constants in a
    :class:`~ampere.backends.torch.parameters.LoweredParameters`, which is an
    ``nn.Module``, so torch's own recursion does the moving and this function
    only has to keep :func:`place`'s three attributes honest afterwards.

    In place and returning *owner*, as ``nn.Module.to`` is — a piece is a
    declaration a problem already holds by identity, so handing back a copy
    would leave the problem pointing at the unmoved one.
    """
    resolved_dtype = owner.dtype if dtype is None else resolve_dtype(dtype)
    resolved_device = owner.device if device is None else resolve_device(device)
    tensors = getattr(owner, "tensors", None)
    if tensors is not None:
        tensors.to(dtype=resolved_dtype, device=resolved_device)
    place(owner, resolved_dtype, resolved_device)
    return owner
