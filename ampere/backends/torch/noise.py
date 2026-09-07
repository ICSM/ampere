"""Noise models on the torch path: the two core compositions, declared as this backend's.

Why this module exists at all
-----------------------------
Nothing here adds arithmetic to ``ampere.core``'s. :class:`IndependentNoise`
declares ``scale`` and ``jitter``, and the quadrature
``sigma_eff² = (s·sigma_data)² + jitter²`` is transcribed into torch by
:mod:`ampere.backends.torch.problem` — transcribed rather than called, because
the core's own :meth:`~ampere.core.NoiseModel.sigma` coerces with ``float()``
and would cut the graph. :class:`GaussianProcessNoise` declares a kernel and a
solve strategy and delegates the whole computation to them.

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

from typing import Any, ClassVar

from ampere.core import GaussianProcessNoise as _CoreGaussianProcessNoise
from ampere.core import IndependentNoise as _CoreIndependentNoise
from ampere.core import GPSolver, Kernel

from ._config import BACKEND
from .gp import DenseGP

__all__ = ["GaussianProcessNoise", "IndependentNoise"]


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
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND


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
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        super().__init__(
            kernel, DenseGP() if solver is None else solver, scale=scale, jitter=jitter
        )
