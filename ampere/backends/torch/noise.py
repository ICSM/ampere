"""Noise models on the torch path: the core compositions and the prediction-aware pair.

Why this module exists at all
-----------------------------
Two of the four classes here add no arithmetic to ``ampere.core``'s.
:class:`IndependentNoise` declares ``scale`` and ``jitter``, and the quadrature
``sigma_eff² = (s·sigma_data)² + jitter²`` is transcribed into torch by
:mod:`ampere.backends.torch.problem` — transcribed rather than called, because
the core's own :meth:`~ampere.core.NoiseModel.sigma` coerces with ``float()``
and would cut the graph. :class:`GaussianProcessNoise` declares a kernel and a
solve strategy and delegates the whole computation to them.

The other two are the **prediction-aware** pair, and they are this backend's
share of X-1 (ruled 2026-09-03, ``likelihoods.md`` §5): a noise whose magnitude
depends on the model — ``sigma_eff² = (s·sigma_data)² + (f·predicted)²`` — is a
``NoiseModel``, not a family, and the ruling makes it *standard library* rather
than backend-neutral vocabulary, so each backend ships one. The reference and
jax backends shipped theirs at W2.1 and W2.5; torch's absence was W2.13's
carried finding ("torch ships no prediction-aware noise models; the
``sigma_tensor`` hook is dormant"), and these close it. They subclass their
``ampere.backends.reference`` counterparts for the same reason the instrument
steps and the jax noise models do: the *declaration* — what ``f`` is, how it is
resolved out of the likelihood-wide values, which composition each class is —
is not arithmetic and must not be written twice.

The ``sigma_tensor`` hook, and its shape
-----------------------------------------
W2.13 left a hook in :mod:`ampere.backends.torch.problem` for exactly this and
never used it. Slice 2 gives it its final signature, which is deliberately
**not** ``NoiseModel.sigma``'s: it takes the base sigma *already quadratured*
— ``sqrt((s·sigma_data)² + jitter²)`` for the retained samples, as a tensor —
rather than the observed container and a retain mask.

Two reasons, and the second is the load-bearing one. The container's
uncertainties are a constant, so converting them to a tensor once at
composition and reusing them is worth roughly a numpy round trip per
evaluation; and ``scale`` and ``jitter`` are **fitted parameters**, so a
``sigma_tensor`` that rebuilt the base by calling the core's own ``sigma``
would coerce them with ``float()`` and silently lose their gradient — the
exact failure mode W2.4 slice 1 found in the GP amplitude. Handing the base in
means the quadrature is written once, in torch, in the lowering, and every
parameter in it keeps its graph.

What these classes add is the **declaration**. Until W2.13 a
:class:`~ampere.core.NoiseModel` carried no capability flags and was not among
``Dataset.capability_parts``, so a problem could report ``backend="torch"``
while its noise model and GP solver ran in numpy — a gap both Phase 2 tracks
recorded. ``inference.md`` §10a (fold-in 7, ruled 2026-09-07) closed it: a
noise model and a solver carry the four flags now, and
:attr:`ampere.core.Likelihood.capability_parts` puts them on the composed
problem. So composing ``ampere.core.IndependentNoise`` into a torch problem is
a **backend disagreement**, refused at composition by name, rather than a
silent numpy island in the middle of a differentiable fit.

The class names are the core's, deliberately
--------------------------------------------
``Likelihood.to_spec`` records ``type(noise).__name__``, and
``tests/conformance/test_cross_backend.py``'s ``ampere_likelihoods`` row
compares that declaration across backends. A backend's noise model is the same
*declaration* as the core's — what differs is which library computes it — so it
must present the same name. Hence ``IndependentNoise``, not
``TorchIndependentNoise``.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.reference.noise import (
    FractionalModelGPNoise as _ReferenceFractionalModelGPNoise,
)
from ampere.backends.reference.noise import (
    FractionalModelNoise as _ReferenceFractionalModelNoise,
)
from ampere.core import GaussianProcessNoise as _CoreGaussianProcessNoise
from ampere.core import JointGaussianProcessNoise as _CoreJointGaussianProcessNoise
from ampere.core import IndependentNoise as _CoreIndependentNoise
from ampere.core import DTYPE, FunctionSamples, GPSolver, Kernel, LikelihoodError, NoiseModel

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, place, to_numpy
from .gp import DenseGP

__all__ = [
    "FractionalModelGPNoise",
    "FractionalModelNoise",
    "GaussianProcessNoise",
    "IndependentNoise",
    # W5.9 -- joint noise over a tuple of channels.
    "JointGaussianProcessNoise",
]


def _fraction(noise: NoiseModel, values: Mapping[str, Any]) -> Any:
    """Resolve ``f`` from the likelihood-wide *values*, **unconverted**.

    The mapping is filtered to this noise model's own names before ``context``:
    ``Likelihood.log_prob`` passes the flat, likelihood-wide resolved values,
    family parameters included, and completing a declaration against names it
    does not own is what the established ``IndependentNoise`` idiom avoids.

    The value is returned as it arrived — a float on the contract path, a
    tensor on the native one — because :func:`_inflate` takes either, and the
    point of this backend is that the second case keeps its gradient.
    """
    own = {key: value for key, value in values.items() if key in noise.parameters}
    return noise.context(own)["f"]


def _check_fraction(fraction: Any, owner: str) -> float:
    """The contract path's validation of ``f``. Concrete values only.

    Deliberately absent from the tensor path: a fitted ``f`` is declared with a
    prior supported on the non-negative half-line, so no sampler proposes a
    value outside it, and a Python ``if`` on a value in the middle of the
    density is what ``inference.md`` §10a rules out. The check is still worth
    having where the value is concrete, because a *fixed* ``f`` is set by hand
    and a hand-set negative one is a typo.
    """
    value = float(fraction)
    if not np.isfinite(value) or value < 0.0:
        raise LikelihoodError(
            f"{owner}: the model-uncertainty fraction f must be finite and non-negative, got "
            f"{value!r}. It multiplies the prediction, so a negative value is not merely "
            f"meaningless — it is indistinguishable from its positive twin and would make the "
            f"posterior bimodal for no reason."
        )
    return value


def _amplitude(
    predicted: Any, *, dtype: torch.dtype, device: torch.device, owner: str
) -> torch.Tensor:
    """``|predicted|`` as a real tensor of *dtype*, complex predictions included.

    ``torch.abs`` rather than the raw tensor: a complex family's prediction is
    complex, and ``likelihoods.md`` §5 leaves projecting it to an amplitude to
    the noise model rather than doing it on the model's behalf.

    The modulus is taken **before** the dtype is imposed, and that ordering is
    the whole of W2.4 slice 3's fix here. Coercing a complex prediction to
    float64 first is not a conversion at all: both
    ``Tensor.to(torch.float64)`` and ``np.asarray(z, dtype=float)`` drop the
    imaginary part with a warning rather than an error, so the amplitude a
    prediction-aware noise model inflates by would silently have been the
    *real part's* magnitude — a different, and wrong, noise model. Slice 3 is
    the first slice in which a complex prediction can reach here at all
    (``complex_gaussian`` in the realised path), so this is the fix arriving
    with the caller that needs it.
    """
    if predicted is None:
        raise LikelihoodError(
            f"{owner} needs the model prediction to compute sigma_eff, but was passed none. "
            f"Every ampere call site supplies it — Likelihood.log_prob, Likelihood.conditional "
            f"and Dataset.draw_observation all pass predicted= (likelihoods.md §5, X-1) — so a "
            f"caller invoking sigma() by hand must pass it too."
        )
    if isinstance(predicted, torch.Tensor):
        return torch.abs(predicted).to(dtype=dtype, device=device)
    array = np.asarray(predicted)
    if array.dtype.kind == "c":
        return as_tensor(np.abs(array), dtype=dtype, device=device)
    return torch.abs(as_tensor(array, dtype=dtype, device=device))


def _inflate(
    base: Any,
    fraction: Any,
    predicted: Any,
    owner: str,
    *,
    dtype: torch.dtype = DEFAULT_DTYPE,
    device: torch.device = DEFAULT_DEVICE,
) -> torch.Tensor:
    """``sqrt(base**2 + (fraction * |predicted|)**2)``, in torch.

    *dtype* and *device* are the owning noise model's own (W2.4 slice 3), so a
    noise model placed on a GPU inflates there rather than dragging the
    prediction back to the CPU once per evaluation.
    """
    scaled = as_tensor(fraction, dtype=dtype, device=device) * _amplitude(
        predicted, dtype=dtype, device=device, owner=owner
    )
    if base is None:
        # No per-sample uncertainties to combine with: the model error is the
        # whole of the noise. Well defined, and a family that cannot live
        # without observed uncertainties has already refused at composition
        # time via NoiseModel.check_compatible.
        return scaled
    return torch.sqrt(as_tensor(base, dtype=dtype, device=device) ** 2 + scaled**2)


def _check_one_device(noise: NoiseModel) -> None:
    """Refuse a GP noise model whose kernel or solver lives somewhere else.

    ``ampere.core.Likelihood.capability_parts`` puts the noise model, the
    solver and — since W3.8 — the kernel on the composed problem, so all three
    are caught by :func:`~ampere.core.declared_capabilities` once there is a
    :class:`~ampere.core.Dataset` to attach them to. This check stays, and is
    not redundant: it fires at the moment the noise model is *built*, before
    any dataset exists, and a kernel left on the CPU inside a GPU problem is
    otherwise found only by torch, per evaluation, deep inside a Cholesky — or,
    worse, not found at all, because torch will happily broadcast a CPU scalar
    against a CUDA matrix. Refusing here is the earliest honest point.
    """
    where = noise.DEVICE
    mismatched = {
        f"{type(part).__name__} on {part.DEVICE!r}"
        for part in (noise.kernel, noise.solver)
        if str(getattr(part, "DEVICE", where)) != where
    }
    if mismatched:
        raise LikelihoodError(
            f"{type(noise).__name__} is on {where!r} but {', '.join(sorted(mismatched))}. Ampere "
            f"does not move arrays between devices on your behalf (architecture.md §5), so a "
            f"kernel or solver placed elsewhere than the noise model that owns it is a mistake "
            f"rather than a plan: pass the same device= to the kernel, to the solver "
            f"(GPSolver.configured(device=...)) and to this noise model."
        )


class IndependentNoise(_CoreIndependentNoise):
    """Uncorrelated noise, declared as this backend's.

    Identical to :class:`ampere.core.IndependentNoise` in every respect but
    the capability flags. See this module's docstring for why that is worth a
    class.

    Parameters
    ----------
    scale
        Optional multiplier on the data's own uncertainties — "the catalogue
        underestimates its error bars". A prior to fit it, or a number.
    jitter
        Optional noise floor, added in quadrature.
    dtype, device
        Where this noise model's arithmetic happens, threaded exactly as on
        every other piece of this backend (W2.4 slice 3). Reported back as the
        instance's ``DEVICE`` capability flag, which is what makes
        :func:`~ampere.core.declared_capabilities` compose a whole problem on
        one device. Never auto-detected.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: Class-level default; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *,
        scale: Any = None,
        jitter: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(scale=scale, jitter=jitter)
        place(self, dtype, device)


class GaussianProcessNoise(_CoreGaussianProcessNoise):
    """The misspecification GP, declared as this backend's.

    The declaration is the core's — a kernel and a solve strategy — and all
    the arithmetic belongs to those two, so there is nothing to reimplement.
    One default matters: *solver* falls back to **this backend's**
    :class:`~ampere.backends.torch.DenseGP` rather than ``ampere.core``'s, so
    the obvious two-argument call cannot quietly put a scipy Cholesky in the
    middle of a torch problem. Passing the core solver explicitly is still
    refused, by :func:`~ampere.core.declared_capabilities`, as a backend
    disagreement.

    Parameters
    ----------
    kernel
        The covariance kernel; its hyperparameters become this noise model's.
        Use one of this backend's (:class:`~ampere.backends.torch.Matern32`,
        :class:`~ampere.backends.torch.SquaredExponential`) so the covariance
        is built in torch and the hyperparameters are trainable — with
        ``ampere.core``'s the solve is still differentiable in the *residual*
        and not in the amplitude or length scale, which is the W2.4 finding
        native kernels exist to fix.
    solver
        GP solve strategy. Defaults to this backend's ``DenseGP()``.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: Class-level default; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        scale: Any = None,
        jitter: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            kernel, DenseGP() if solver is None else solver, scale=scale, jitter=jitter
        )
        place(self, dtype, device)
        _check_one_device(self)


class FractionalModelNoise(_ReferenceFractionalModelNoise):
    """``sigma_eff**2 = (s * sigma_data)**2 + (f * predicted)**2``, in torch.

    ``likelihoods.md`` §5's X-1 noise model: ``f`` is the fractional
    uncertainty on the *model*, which is the honest way to say "I do not
    believe this forward model to better than 10 %" without inventing a family
    or inflating the data's own error bars, which would say something
    different. A :class:`~ampere.core.GaussianFamily` with one of these stays
    :attr:`~ampere.core.Marginalisation.ANALYTIC`: a prediction-dependent
    sigma is still diagonal, and the marginalisation machinery never inspects
    *how* sigma was computed.

    Parameters
    ----------
    f
        The model-uncertainty fraction. A frozen ``scipy.stats`` distribution
        to fit it — the bijection is read off the prior's own support — or a
        number to hold it fixed.
    scale
        Optional multiplier on the *data's* uncertainties. Note it multiplies
        ``sigma_data`` only, inside the quadrature, exactly as §5 writes it.
    jitter
        Optional noise floor added in quadrature to ``sigma_data``.
    dtype, device
        Where the quadrature happens (W2.4 slice 3), reported as ``DEVICE``.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: Class-level default; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        f: Any,
        *,
        scale: Any = None,
        jitter: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(f, scale=scale, jitter=jitter)
        place(self, dtype, device)

    def sigma_tensor(
        self,
        base: torch.Tensor | None,
        values: Mapping[str, Any],
        *,
        predicted: torch.Tensor,
    ) -> torch.Tensor:
        """The native surface: ``sigma_eff`` as a tensor, gradient intact.

        *base* is the quadrature of ``scale`` and ``jitter`` against the data's
        own uncertainties, computed in torch by the lowering; see this module's
        docstring for why it arrives rather than being rebuilt here.
        """
        return _inflate(
            base,
            _fraction(self, values),
            predicted,
            type(self).__name__,
            dtype=self.dtype,
            device=self.device,
        )

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray:
        """The contract surface: numpy in, numpy out, and it validates ``f``."""
        base = super(_ReferenceFractionalModelNoise, self).sigma(
            observed, retain, values, predicted=predicted
        )
        fraction = _check_fraction(_fraction(self, values), type(self).__name__)
        inflated = _inflate(
            base,
            fraction,
            predicted,
            type(self).__name__,
            dtype=self.dtype,
            device=self.device,
        )
        return to_numpy(inflated).astype(DTYPE, copy=False)


class FractionalModelGPNoise(_ReferenceFractionalModelGPNoise):
    """A misspecification GP *and* a fractional model uncertainty, in torch.

    ``likelihoods.md`` §5's composition, which is a subclass overriding only
    ``sigma``: the kernel, the solver and the ``noise_params`` assembly are all
    inherited, and the marginal likelihood is
    ``N(0, K + diag(sigma_data**2 + (f*mu)**2))`` with no new mathematics.

    The two terms say different things and are not redundant. The GP models
    *correlated* structure the forward model gets wrong — a missing feature, a
    continuum slope — while ``f`` models an uncorrelated,
    amplitude-proportional error. Fitting both lets the data decide how much of
    the misfit is smooth.

    Parameters
    ----------
    kernel
        The covariance kernel; its hyperparameters become this noise model's.
        Use one of this backend's so the covariance is built in torch.
    solver
        GP solve strategy. Defaults to this backend's
        :class:`~ampere.backends.torch.DenseGP` rather than ``ampere.core``'s,
        so the obvious two-argument call cannot quietly put a scipy Cholesky in
        the middle of a torch problem;
        :class:`~ampere.backends.torch.QuasisepGP` is the O(N) alternative.
    f
        The model-uncertainty fraction, keyword-only so it can never be
        confused with the solver.
    scale, jitter
        Applied to the data's uncertainties inside the quadrature.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: Class-level default; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        f: Any,
        scale: Any = None,
        jitter: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            kernel, DenseGP() if solver is None else solver, f=f, scale=scale, jitter=jitter
        )
        place(self, dtype, device)
        _check_one_device(self)

    def sigma_tensor(
        self,
        base: torch.Tensor | None,
        values: Mapping[str, Any],
        *,
        predicted: torch.Tensor,
    ) -> torch.Tensor:
        """The native surface: ``sigma_eff`` as a tensor, gradient intact."""
        return _inflate(
            base,
            _fraction(self, values),
            predicted,
            type(self).__name__,
            dtype=self.dtype,
            device=self.device,
        )

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray:
        base = super(_ReferenceFractionalModelGPNoise, self).sigma(
            observed, retain, values, predicted=predicted
        )
        fraction = _check_fraction(_fraction(self, values), type(self).__name__)
        inflated = _inflate(
            base,
            fraction,
            predicted,
            type(self).__name__,
            dtype=self.dtype,
            device=self.device,
        )
        return to_numpy(inflated).astype(DTYPE, copy=False)


# ---------------------------------------------------------------------------
# W5.9 -- joint noise over a tuple of channels
# ---------------------------------------------------------------------------


class JointGaussianProcessNoise(_CoreJointGaussianProcessNoise):
    """The shared-grid intrinsic coregionalisation model, declared as this backend's.

    The declaration is the core's — a kernel, a coupling ``B``, a solve
    strategy and the tuple of datasets the process spans — and all the
    arithmetic belongs to the kernel, the coupling and the solver, so there is
    nothing to reimplement. What this class adds is the flags and one default
    that matters: *solver* falls back to **this backend's**
    :class:`~ampere.backends.torch.DenseGP` rather than ``ampere.core``'s, so
    the obvious call cannot quietly put a scipy Cholesky in the middle of a
    torch problem.

    The ``T`` rotated scalar solves are composed natively by
    :mod:`ampere.backends.torch.problem`, which reads the coupling's
    eigendecomposition through :meth:`~ampere.core.ChannelCoupling.eigen` with
    ``xp=torch`` — one implementation of each parameterisation, evaluated on
    tensors rather than transcribed, so ``B``'s parameters get a gradient and
    NUTS can sample the channel coupling alongside the kernel's
    hyperparameters.

    Parameters
    ----------
    kernel, solver, datasets, coupling, scale, jitter
        As :class:`ampere.core.JointGaussianProcessNoise`. Use this backend's
        kernel and solver so the covariance is built in torch.
    dtype, device
        Where this noise model's arithmetic happens, threaded exactly as on
        every other piece of this backend. Never auto-detected.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        datasets: Any,
        coupling: Any,
        scale: Any = None,
        jitter: Any = None,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        super().__init__(
            kernel,
            DenseGP() if solver is None else solver,
            datasets=datasets,
            coupling=coupling,
            scale=scale,
            jitter=jitter,
        )
        place(self, dtype, device)
        _check_one_device(self)
