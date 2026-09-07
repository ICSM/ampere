"""Prediction-aware noise on the jax path: ``FractionalModelNoise`` and its GP composition.

``likelihoods.md`` §5 fixes the name and the semantics —
``sigma_eff**2 = (s * sigma_data)**2 + (f * predicted)**2`` — and X-1's ruling
(2026-09-03) makes this a *standard library* noise model rather than part of
the backend-neutral vocabulary. So each backend ships one, and this is the jax
one: the same two classes, the same declaration, ``f`` an ordinary fitted
parameter, and the quadrature done in ``jax.numpy``.

Both subclass their ``ampere.backends.reference`` counterparts for the same
reason the instrument steps do (see :mod:`ampere.backends.jax.instrument`): the
declaration — what ``f`` is, how it is resolved out of the likelihood-wide
values, which composition each class is — is not arithmetic and must not be
written twice, or the two backends would eventually declare different things
and the conformance suite would call the difference numerical.

``sigma`` returns numpy
-----------------------
``NoiseModel.sigma``'s callers (``Likelihood.log_prob``,
``Likelihood.conditional``, ``Dataset.draw_observation``) are the core's, and
they work in numpy — so the contract surface converts at its boundary, exactly
as the models' ``evaluate`` does. :meth:`FractionalModelNoise.sigma_jax` is the
traced surface that :mod:`ampere.backends.jax.problem` composes, and it is
where the gradient with respect to ``f`` actually flows.

The capability flag, and a gap worth naming
-------------------------------------------
Both classes declare ``BACKEND = "jax"``. They are **not** among
``Dataset.capability_parts`` today — only the instrument steps and the compiled
models are — so nothing here contributes to the flag a composed problem
reports. W2.12 recorded that as a finding: a problem could report
``backend="jax"`` while its noise model and GP solver ran in numpy. Declaring
it anyway is what makes the answer already correct if that parts set is ever
widened, and widening it is a §4 change awaiting a ruling, not this item's to
take.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np

from ampere.backends.reference.noise import (
    FractionalModelGPNoise as _ReferenceFractionalModelGPNoise,
)
from ampere.backends.reference.noise import (
    FractionalModelNoise as _ReferenceFractionalModelNoise,
)
from ampere.core import FunctionSamples, GPSolver, Kernel, LikelihoodError, NoiseModel

from ._config import BACKEND, require_x64
from .gp import DenseGP

__all__ = ["FractionalModelGPNoise", "FractionalModelNoise"]


def _fraction(noise: NoiseModel, values: Mapping[str, Any]) -> Any:
    """Resolve ``f`` from the likelihood-wide *values*, unconverted.

    The mapping is filtered to this noise model's own names before ``context``:
    ``Likelihood.log_prob`` passes the flat, likelihood-wide resolved values,
    family parameters included, and completing a declaration against names it
    does not own is what the established ``IndependentNoise`` idiom avoids.

    The value is returned as it arrived — a float on the contract path, a jax
    array or tracer on the native one — because :func:`_inflate` is written to
    take either and the point of this backend is that the second case works.
    """
    own = {key: value for key, value in values.items() if key in noise.parameters}
    return noise.context(own)["f"]


def _check_fraction(fraction: Any, owner: str) -> float:
    """The contract path's validation of ``f``. Concrete values only.

    Deliberately absent from the traced path: a Python ``if`` on a tracer is an
    error rather than a branch, and the declaration is what keeps the condition
    from arising there — ``f`` is declared with a prior supported on the
    non-negative half-line, so no sampler proposes a value outside it. The
    check is still worth having here, where the value is concrete, because a
    *fixed* ``f`` is set by hand and a hand-set negative one is a typo.
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


def _inflate(base: Any, fraction: Any, predicted: Any, owner: str) -> jax.Array:
    """``sqrt(base**2 + (fraction * |predicted|)**2)``, in jax, float64.

    ``jnp.abs`` rather than the raw array: a complex family's prediction is
    complex128, and ``likelihoods.md`` §5 leaves projecting it to an amplitude
    to the noise model rather than doing it on the model's behalf.
    """
    if predicted is None:
        raise LikelihoodError(
            f"{owner} needs the model prediction to compute sigma_eff, but was passed none. "
            f"Every ampere call site supplies it — Likelihood.log_prob, Likelihood.conditional "
            f"and Dataset.draw_observation all pass predicted= (likelihoods.md §5, X-1) — so a "
            f"caller invoking sigma() by hand must pass it too."
        )
    scaled = jnp.asarray(fraction, dtype=jnp.float64) * jnp.abs(
        jnp.asarray(predicted, dtype=jnp.float64)
    )
    if base is None:
        # No per-sample uncertainties to combine with: the model error is the
        # whole of the noise. Well defined, and a family that cannot live
        # without observed uncertainties has already refused at composition
        # time via NoiseModel.check_compatible.
        return scaled
    return jnp.sqrt(jnp.asarray(base, dtype=jnp.float64) ** 2 + scaled**2)


class FractionalModelNoise(_ReferenceFractionalModelNoise):
    """``sigma_eff**2 = (s * sigma_data)**2 + (f * predicted)**2``, in jax.

    ``f`` is the fractional uncertainty on the *model*: the honest way to say
    "I do not believe this forward model to better than 10 %" without inventing
    a family or inflating the data's own error bars, which would say something
    different. A ``GaussianFamily`` with one of these stays
    :attr:`~ampere.core.Marginalisation.ANALYTIC` — a prediction-dependent
    sigma is still diagonal.

    Parameters
    ----------
    f
        The model-uncertainty fraction. A frozen ``scipy.stats`` distribution
        to fit it, or a number to hold it fixed.
    scale
        Optional multiplier on the *data's* uncertainties, inside the
        quadrature.
    jitter
        Optional noise floor added in quadrature to ``sigma_data``.
    """

    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, f: Any, *, scale: Any = None, jitter: Any = None) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(f, scale=scale, jitter=jitter)

    def sigma_jax(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: Any = None,
    ) -> jax.Array:
        """The native surface: ``sigma_eff`` as a jax array, gradient intact."""
        base = super(_ReferenceFractionalModelNoise, self).sigma(
            observed, retain, values, predicted=None
        )
        return _inflate(base, _fraction(self, values), predicted, type(self).__name__)

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray:
        base = super(_ReferenceFractionalModelNoise, self).sigma(
            observed, retain, values, predicted=predicted
        )
        fraction = _check_fraction(_fraction(self, values), type(self).__name__)
        return np.asarray(_inflate(base, fraction, predicted, type(self).__name__))


class FractionalModelGPNoise(_ReferenceFractionalModelGPNoise):
    """A misspecification GP *and* a fractional model uncertainty, in jax.

    ``likelihoods.md`` §5's composition: a subclass overriding only ``sigma``,
    so the marginal likelihood is
    ``N(0, K + diag(sigma_data**2 + (f*mu)**2))`` with no new mathematics. The
    two terms say different things and are not redundant — the GP models
    *correlated* structure the forward model gets wrong, ``f`` an uncorrelated,
    amplitude-proportional error — and fitting both lets the data decide how
    much of the misfit is smooth.

    Parameters
    ----------
    kernel
        The covariance kernel; its hyperparameters become this noise model's.
    solver
        GP solve strategy. Defaults to this backend's :class:`~ampere.backends.jax.DenseGP`
        rather than ``ampere.core``'s, so a jax noise model does not silently
        put a numpy Cholesky in the middle of a jax problem.
    f
        The model-uncertainty fraction, keyword-only so it can never be
        confused with the solver.
    scale, jitter
        Applied to the data's uncertainties inside the quadrature.
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
        f: Any,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        require_x64(f"a jax {type(self).__name__}")
        super().__init__(
            kernel, DenseGP() if solver is None else solver, f=f, scale=scale, jitter=jitter
        )

    def sigma_jax(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: Any = None,
    ) -> jax.Array:
        """The native surface: ``sigma_eff`` as a jax array, gradient intact."""
        base = super(_ReferenceFractionalModelGPNoise, self).sigma(
            observed, retain, values, predicted=None
        )
        return _inflate(base, _fraction(self, values), predicted, type(self).__name__)

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
        return np.asarray(_inflate(base, fraction, predicted, type(self).__name__))
