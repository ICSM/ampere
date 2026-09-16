"""The torch backend: rung 2 of ``architecture.md`` §1's capability ladder.

Install the ``torch`` extra -- ``pip install ".[torch]"`` from a checkout
(PyPI's ``ampere`` package is unrelated). Everything here is written against
the frozen contracts of ``ampere.core`` and adds exactly what
``DEVELOPMENT_PLAN.md`` §3 says a backend is: "an array library for writing
models and transformations, a set of GP solver implementations, and the extra
inference engines its differentiability unlocks". Kernels, containers,
likelihood families, datasets and the fitting problem itself are **not**
here — they are backend-neutral and live in ``ampere.core``.

What is shipped (W2.4 slice 1, plus W2.13's realisation)
---------------------------------------------------------
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
    therefore differentiable — with :class:`Matern32` and
    :class:`SquaredExponential`, the core kernel declarations with their
    covariances built in torch (W2.13), so a GP hyperparameter is trainable
    rather than merely declared.
:mod:`~ampere.backends.torch.noise`
    :class:`IndependentNoise` and :class:`GaussianProcessNoise` — the core
    compositions declared as this backend's, which W2.13's capability
    widening makes necessary rather than decorative.
:mod:`~ampere.backends.torch.problem`
    :func:`lower_problem` — a composed :class:`~ampere.core.FittingProblem` as
    one differentiable function of the unconstrained free vector
    (``inference.md`` §10a). Registered with ``ampere.core`` at import, which
    is how :class:`ampere.inference.NUTSEngine` reaches it without importing
    this package.

**Every piece here declares ``BACKEND = "torch"``** beside ``DIFFERENTIABLE``,
``BATCHABLE`` and ``DEVICE``, explicitly rather than by inheriting the ABCs'
defaults. That string is this backend's one name everywhere (W2.12): what
``ampere.core.declared_capabilities`` aggregates onto a problem, what a run's
``ampere_backend`` records, the key ``lowering.md`` §12.8's registry is
consulted with, and the id of this backend's conformance fixture.

Precision and device
--------------------
float64, on the CPU by default, threaded explicitly through every
construction; ``torch.set_default_dtype`` is **never** called, because it is
global mutable process state (``lowering.md`` §10.1).

**Since W2.4 slice 3 the device is a per-instance choice.** Every model,
instrument step, kernel, noise model and solver here takes a ``device=``
keyword beside its ``dtype=``, defaults to ``"cpu"``, and reports what it was
given back as its own ``DEVICE`` capability flag — an instance attribute
shadowing the class default, which is all
:func:`~ampere.core.declared_capabilities` needs, since it reads
``part.DEVICE``. So a whole problem is composed on one named device or refused
at composition with both devices named, ``.to(...)`` moves a piece's buffers
and re-declares it, and nothing is ever auto-detected: a machine with a GPU
present takes CI's path unless a caller says otherwise (``architecture.md``
§5). ``tests/gpu/`` is the API-level proof, skipped without CUDA;
``tests/backends/test_torch_device.py`` proves the declaration and composition
rules on the CPU, using ``torch.device("meta")`` as a stand-in second device.

How far differentiability reaches — and how W2.13 closed the gap
-----------------------------------------------------------------
``DIFFERENTIABLE = True`` on the pieces here is a claim about the *pieces*,
and it is true of them: each has a tensor-valued entry point
(``Model.evaluate_tensor``, ``Transformation.apply_tensor``,
``DenseGP.log_marginal_likelihood_tensor``,
``TorchParameterSpace.log_prior_unconstrained_tensor``) that is a
differentiable function of tensor inputs, and ``tests/backends/`` takes the
gradient of each to prove it.

W2.4 slice 1 shipped without an end-to-end ``d log_prob / dθ`` on a composed
:class:`~ampere.core.FittingProblem`, and recorded why: ``log_prob`` runs the
chain through ``ampere.core``'s containers, which hold **numpy** arrays, and a
tensor carrying an autograd graph cannot be put in one. So the graph is cut at
the first container, on every backend.

**That is what a realisation is for.** ``inference.md`` §10a (ruled
2026-09-07) makes the differentiable form of a problem a first-class,
registered object rather than a property the contract path was supposed to
have: :mod:`ampere.backends.torch.problem` composes the same quantity out of
the native surfaces above and hands back a differentiable function of the
unconstrained vector, and ``ampere.core.realise`` checks it against the
contract path before anyone samples with it. Two things follow:

* the aggregated ``differentiable`` flag still means what ``inference.md`` §10
  says it means — a property of the pieces — but it now covers **all** of
  them, since W2.13 widened the parts to include noise models and GP solvers.
  A problem reporting ``differentiable=True`` on this backend is one
  :func:`~ampere.core.realise` will lower;
* the NUTS driver ships. ``NUTSEngine(problem)`` on a torch problem runs
  pyro's sampler over this backend's realisation, with no density argument and
  no backend import inside ``ampere.inference``.

What slice 2 added
------------------
* :class:`QuasisepGP` — exact O(N), differentiable in the kernel
  hyperparameters, over celerite2's compiled semiseparable kernels wrapped as
  ``torch.autograd`` functions (:mod:`ampere.backends.torch._celerite`). That
  is ``DEVELOPMENT_PLAN.md`` §6's deferred library choice, settled by
  measurement against GPyTorch rather than from documentation; the
  decision-log row of 2026-09-07 carries the table. It supplies the O(N)
  ``conditional_loo`` recursion W2.3 deferred, too.
* the **prediction-aware noise models**, :class:`FractionalModelNoise` and
  :class:`FractionalModelGPNoise`, closing W2.13's carried finding that the
  ``sigma_tensor`` hook was in place and dormant.
* a **widened realisation**: censoring (the Tobit form), the non-Gaussian
  families ``ampere.core`` implements, and the quasiseparable solver inside
  the density. What is still unsupported is refused by name at construction,
  each with its own reason — see :mod:`ampere.backends.torch.problem`.
* **batching**, honestly:
  :meth:`~ampere.backends.torch.LoweredProblem.log_prob_unconstrained_batched`
  evaluates a stack of free vectors in one ``torch.func.vmap`` call, and
  ``BATCHABLE`` is ``True`` on every piece that survives it — not on
  :class:`QuasisepGP`, whose solve is a compiled extension, so a problem
  carrying it reports ``batchable=False`` and the batched call refuses.
* **the float64 opt-out**, ``GPSolver.configured(dtype=..., device=...)``:
  ``architecture.md`` §5's per-run precision choice, recorded through
  ``provenance_config()`` and deliberately absent from the spec hash.
* **variational inference** through :class:`ampere.inference.VIEngine`, which
  reaches this backend the same way ``NUTSEngine`` does.

What slice 3 added
------------------
* the **per-instance device** described above, and the GPU smoke suite;
* **``complex_gaussian``**: the circular complex Gaussian, transcribed into
  :mod:`ampere.backends.torch._families` and composing in the realised path
  with complex tensors end to end. Its *correlated* form — the circular
  complex GP, which ``likelihoods.md`` §4 declares analytic — was refused by
  name at slice 3, because ``ampere.core`` declared
  ``GP_ANALYTIC_IMPLEMENTED = False`` and refused it first: a realisation with
  no numpy oracle would be this backend inventing a likelihood. **W4.2 lifted
  that**: ``ampere.core`` has the closed form, and here it is the solver rather
  than the family body that carries it, since ``DenseGP`` now takes the
  two-column right-hand side the circular GP reduces to.
* the last **``sigma_tensor`` consumer gap**: a noise model that overrides
  ``sigma`` without supplying the native hook used to fall back silently to
  the base quadrature in the realised density — a different likelihood from
  the contract path's, arrived at without a word. It is refused by name now,
  and ``|predicted|`` is taken as a modulus rather than through a float cast,
  so a complex prediction inflates by its amplitude rather than by the
  magnitude of its real part.

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
from .astropy import TRANSLATIONS, NativeAstropyModel, from_astropy
from .gp import (
    SHO,
    DenseGP,
    HilbertSpaceGP,
    Matern12,
    Matern32,
    Matern52,
    Product,
    QuasisepGP,
    RotationTerm,
    SpectralMixture,
    SquaredExponential,
    Sum,
)
from .astrometry import (
    EpochSample,
    ReflexOrbit,
)
from .image import (
    PSFConvolution,
    TorchImageStep,
)
from .instrument import (
    DETECTORS,
    CalibrationScale,
    LSFConvolution,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    bundled_filter_library,
)
from .interferometry import (
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    ClosurePhase,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    TimeSmearing,
    TorchInterferometryStep,
    UniformDisc,
    UniformDiscVisibilities,
)
from .lowering import (
    LoweredPrior,
    LoweringFallbackWarning,
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
from .noise import (
    FractionalModelGPNoise,
    FractionalModelNoise,
    GaussianProcessNoise,
    IndependentNoise,
    # W5.9 -- joint noise over a tuple of channels.
    JointGaussianProcessNoise,
)
from .parameters import LoweredParameters, TorchParameterSpace
from .problem import LoweredProblem, lower_problem
from .rng import generator, seed_for
from .sharding import DistributedSharder, SingleDeviceSharder, available_devices

# ``inference.md`` §10a: importing this package is the user's opt-in to torch,
# and it is also the moment the torch realisation becomes reachable through
# ``ampere.core.realise`` -- which is how ``ampere.inference``'s gradient-based
# drivers get a differentiable density without importing a backend.
from ampere.core import register_realisation as _register_realisation

_register_realisation(BACKEND, lower_problem, builtin=True)

__all__ = [
    "BACKEND",
    "COORDINATE_UNIT",
    "DEFAULT_DEVICE",
    "DEFAULT_DTYPE",
    "DETECTORS",
    "FLUX_UNIT",
    "SHO",
    "TRANSLATIONS",
    "Amplitude",
    "BandwidthSmearing",
    "Binary",
    "BinaryVisibilities",
    "BlackBody",
    "CalibrationScale",
    "ClosurePhase",
    "DenseGP",
    "DistributedSharder",
    "EpochSample",
    "FourierSample",
    "FractionalModelGPNoise",
    "FractionalModelNoise",
    "GaussianProcessNoise",
    "GaussianSource",
    "GaussianSourceVisibilities",
    "HilbertSpaceGP",
    "IndependentNoise",
    "JointGaussianProcessNoise",
    "LSFConvolution",
    "LoweredParameters",
    "LoweredPrior",
    "LoweredProblem",
    "LoweringFallbackWarning",
    "Matern12",
    "Matern32",
    "Matern52",
    "ModifiedBlackBody",
    "NativeAstropyModel",
    "PSFConvolution",
    "PowerLaw",
    "Product",
    "QuasisepGP",
    "ReflexOrbit",
    "Resample",
    "RotationTerm",
    "SingleDeviceSharder",
    "SpectralMixture",
    "SquaredExponential",
    "Sum",
    "SyntheticPhotometry",
    "TimeSmearing",
    "TorchImageStep",
    "TorchInterferometryStep",
    "TorchParameterSpace",
    "TorchSpectralModel",
    "UniformDisc",
    "UniformDiscVisibilities",
    "as_tensor",
    "available_devices",
    "bin_edges",
    "bundled_filter_library",
    "from_astropy",
    "generator",
    "lower_bijection",
    "lower_hierarchical",
    "lower_prior",
    "lower_problem",
    "planck_jy",
    "seed_for",
    "to_numpy",
]
