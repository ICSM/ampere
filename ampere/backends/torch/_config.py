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
same code path as CI. GPU execution is W2.4 slice 2.
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
