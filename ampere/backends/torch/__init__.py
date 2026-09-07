"""The torch backend: rung 2 of ``architecture.md`` §1's capability ladder.

``pip install "ampere[torch]"``. Everything here is written against the frozen
contracts of ``ampere.core`` and adds exactly what ``DEVELOPMENT_PLAN.md`` §3
says a backend is: "an array library for writing models and transformations, a
set of GP solver implementations, and the extra inference engines its
differentiability unlocks". Kernels, containers, likelihood families, datasets
and the fitting problem itself are **not** here — they are backend-neutral and
live in ``ampere.core``.

What is shipped (W2.4 slice 1)
-------------------------------
:mod:`~ampere.backends.torch.lowering`
    ``lowering.md`` §3 and §4 — the prior-family and bijection tables,
    registered as built-in rows in ``ampere.core.lowering``'s registry under
    the key ``"torch"``, with :class:`~ampere.core.exceptions.LoweringError`
    for the families torch cannot express exactly and no silent approximations
    anywhere.
:mod:`~ampere.backends.torch.parameters`
    ``lowering.md`` §5 to §8 — a merged :class:`~ampere.core.ParameterSet` lowered
    onto a nested ``torch.nn.Module``, with the prior, the bijection, the
    unconstrained density and its Jacobian computed natively and
    differentiably.
:mod:`~ampere.backends.torch.rng`
    ``lowering.md`` §9 — named sub-streams from
    :func:`ampere.core.rng.substream`, and the two routes round torch's
    missing ``generator=`` on ``Distribution.sample``.
:mod:`~ampere.backends.torch.models`
    :class:`BlackBody`, :class:`ModifiedBlackBody`, :class:`PowerLaw` — the
    ``DEVELOPMENT_PLAN.md`` §5 trio, in torch.
:mod:`~ampere.backends.torch.instrument`
    :class:`CalibrationScale`, :class:`Resample`, :class:`LSFConvolution`,
    :class:`SyntheticPhotometry` — ``transformations.md`` §10's table, in
    torch.
:mod:`~ampere.backends.torch.gp`
    :class:`DenseGP` — the exact dense solve through torch's own Cholesky, and
    therefore differentiable.

**Every piece here declares ``BACKEND = "torch"``** beside ``DIFFERENTIABLE``,
``BATCHABLE`` and ``DEVICE``, explicitly rather than by inheriting the ABCs'
defaults. That string is this backend's one name everywhere (W2.12): what
``ampere.core.declared_capabilities`` aggregates onto a problem, what a run's
``ampere_backend`` records, the key ``lowering.md`` §12.8's registry is
consulted with, and the id of this backend's conformance fixture.

Precision and device
--------------------
float64, on the CPU, threaded explicitly through every construction;
``torch.set_default_dtype`` is **never** called, because it is global mutable
process state (``lowering.md`` §10.1). GPU execution and batching are W2.4
slice 2 and are declared ``False``/``"cpu"`` honestly until they exist.

How far differentiability currently reaches — read this before relying on it
----------------------------------------------------------------------------
``DIFFERENTIABLE = True`` on the models and steps here is a claim about the
*pieces*, and it is true of them: each has a tensor-valued entry point
(``Model.evaluate_tensor``, ``Transformation.apply_tensor``,
``DenseGP.log_marginal_likelihood_tensor``,
``TorchParameterSpace.log_prior_unconstrained_tensor``) that is a
differentiable function of tensor inputs, and ``tests/backends/`` takes the
gradient of each to prove it.

What does **not** yet exist is an end-to-end ``d log_prob / dθ`` on a composed
:class:`~ampere.core.FittingProblem`, and the reason is structural rather than
an omission here. ``FittingProblem.log_prob`` runs the chain through
``ampere.core``'s containers — :class:`~ampere.core.Spectrum` and its siblings
— which hold **numpy** arrays and validate their axes in numpy; a tensor
carrying an autograd graph cannot be put in one (it raises
``RuntimeError: Can't call numpy() on Tensor that requires grad``). So the
graph is cut at the first container, and every backend's gradient stops there,
not just this one.

Two consequences, both stated rather than worked around:

* **the aggregated ``differentiable`` flag is about the pieces, not the
  composition.** That is what ``inference.md`` §10 says it is ("properties of
  the *pieces*"), and it is consistent — but a consumer reading
  ``problem.differentiable is True`` as "``torch.autograd.grad`` works on
  ``problem.log_prob``" would be wrong today. This is W2.4's principal finding
  and it is recorded as such;
* **no NUTS driver ships in slice 1.** A gradient-based sampler over a
  potential whose gradient does not exist is not a sampler, and neither
  finite-differencing it nor differentiating the prior alone would be honest.
  What the driver needs is a decision — Peter's — about where a
  tensor-valued evaluation path lives: a container that can hold a backend
  array, a backend-supplied evaluation route beside the numpy one, or a
  narrowing of what ``differentiable`` promises. The pieces above are written
  so that any of the three is a small amount of further work rather than a
  rewrite.

What slice 2 owes
-----------------
``QuasisepGP`` on the native path (GPyTorch's structured solvers against
celerite2's torch interface, *measured* against the conformance suite rather
than chosen from documentation — ``DEVELOPMENT_PLAN.md`` §6), variational
inference, batching and GPU, and the benchmark rows.

Examples
--------
>>> import numpy as np, scipy.stats as st
>>> from ampere.backends.torch import PowerLaw
>>> grid = np.array([1.0, 2.0, 4.0])
>>> model = PowerLaw(grid, norm=st.loguniform(0.1, 10.0), index=st.norm(-1.0, 0.5))
>>> model.BACKEND, model.DIFFERENTIABLE
('torch', True)
>>> result = model(norm=2.0, index=-1.0)
>>> np.round(result["default"].values, 6).tolist()
[2.0, 1.0, 0.5]

The tensor entry point is the same computation with the graph intact:

>>> import torch
>>> norm = torch.tensor(2.0, dtype=torch.float64, requires_grad=True)
>>> flux = model.evaluate_tensor(norm=norm, index=torch.tensor(-1.0, dtype=torch.float64))
>>> flux["default"].sum().backward()
>>> float(norm.grad)
1.75
"""

from __future__ import annotations

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, to_numpy
from .gp import DenseGP
from .instrument import (
    DETECTORS,
    CalibrationScale,
    LSFConvolution,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    bundled_filter_library,
)
from .lowering import (
    IcdfFallbackWarning,
    LoweredPrior,
    lower_bijection,
    lower_hierarchical,
    lower_prior,
)
from .models import (
    COORDINATE_UNIT,
    FLUX_UNIT,
    BlackBody,
    ModifiedBlackBody,
    PowerLaw,
    TorchSpectralModel,
    planck_jy,
)
from .parameters import LoweredParameters, TorchParameterSpace
from .rng import generator, seed_for

__all__ = [
    "BACKEND",
    "COORDINATE_UNIT",
    "DEFAULT_DEVICE",
    "DEFAULT_DTYPE",
    "DETECTORS",
    "FLUX_UNIT",
    "BlackBody",
    "CalibrationScale",
    "DenseGP",
    "IcdfFallbackWarning",
    "LSFConvolution",
    "LoweredParameters",
    "LoweredPrior",
    "ModifiedBlackBody",
    "PowerLaw",
    "Resample",
    "SyntheticPhotometry",
    "TorchParameterSpace",
    "TorchSpectralModel",
    "as_tensor",
    "bin_edges",
    "bundled_filter_library",
    "generator",
    "lower_bijection",
    "lower_hierarchical",
    "lower_prior",
    "planck_jy",
    "seed_for",
    "to_numpy",
]
