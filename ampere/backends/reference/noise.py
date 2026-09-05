"""Prediction-aware noise: ``FractionalModelNoise`` and its GP composition.

``likelihoods.md`` §5 fixes the name and the semantics —
``sigma_eff**2 = (s * sigma_data)**2 + (f * predicted)**2`` — and X-1's ruling
(2026-09-03) puts the implementation here rather than in ``ampere.core``: it is
a *standard library* noise model, in the sense ``transformations.md`` §10 names
the standard chain steps, not part of the backend-neutral vocabulary every
backend must reimplement.

This is the single most requested thing missing from legacy ampere's
likelihood. A noise whose magnitude depends on the model — a fractional model
uncertainty, a model-variance weighting — is a ``NoiseModel``, not a family:
without the prediction the only way to express one is to re-implement the
sampling distribution, which welds noise to family, cannot be reused, and
cannot reach the GP path.

Two classes, because §5 asks for two compositions:

:class:`FractionalModelNoise`
    The diagonal case. A ``GaussianFamily`` with one of these stays
    :attr:`~ampere.core.Marginalisation.ANALYTIC` — a prediction-dependent
    sigma is still diagonal, and the marginalisation machinery never inspects
    *how* sigma was computed.
:class:`FractionalModelGPNoise`
    §5's "10 % model error *and* a misspecification GP": a
    :class:`~ampere.core.GaussianProcessNoise` overriding only ``sigma``, so
    the marginal likelihood is ``N(0, K + diag(sigma_data**2 + (f*mu)**2))``
    with no new mathematics.

In both, ``f`` is an **ordinary fitted parameter** (W2.1): pass a prior to fit
it, or a number to hold it fixed.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from ampere.core import (
    DTYPE,
    FunctionSamples,
    GaussianProcessNoise,
    GPSolver,
    IndependentNoise,
    Kernel,
    LikelihoodError,
    NoiseModel,
)

from ._declare import as_parameter

__all__ = ["FractionalModelGPNoise", "FractionalModelNoise"]


def _fraction(noise: NoiseModel, values: Mapping[str, Any]) -> float:
    """Resolve ``f`` from the likelihood-wide *values*.

    The mapping is filtered to this noise model's own names before
    ``context``: ``Likelihood.log_prob`` passes the flat, likelihood-wide
    resolved values, family parameters included, and completing a declaration
    against names it does not own is exactly what the established
    ``IndependentNoise``/``GaussianProcessNoise`` idiom exists to avoid.
    """
    own = {key: value for key, value in values.items() if key in noise.parameters}
    resolved = noise.context(own)
    fraction = float(resolved["f"])
    if not np.isfinite(fraction) or fraction < 0.0:
        raise LikelihoodError(
            f"{type(noise).__name__}: the model-uncertainty fraction f must be finite and "
            f"non-negative, got {fraction!r}. It multiplies the prediction, so a negative value "
            f"is not merely meaningless — it is indistinguishable from its positive twin and "
            f"would make the posterior bimodal for no reason."
        )
    return fraction


def _inflate(
    base: np.ndarray | None,
    fraction: float,
    predicted: np.ndarray | None,
    owner: str,
) -> np.ndarray:
    """``sqrt(base**2 + (fraction * |predicted|)**2)``, in float64.

    ``np.abs`` rather than the raw array: a complex family's prediction is
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
    term = fraction * np.abs(np.asarray(predicted))
    scaled = np.asarray(term, dtype=DTYPE)
    if base is None:
        # No per-sample uncertainties to combine with: the model error is the
        # whole of the noise. Well defined, and a family that cannot live
        # without observed uncertainties has already refused at composition
        # time via NoiseModel.check_compatible.
        return scaled
    return np.sqrt(np.asarray(base, dtype=DTYPE) ** 2 + scaled**2)


class FractionalModelNoise(IndependentNoise):
    """``sigma_eff**2 = (s * sigma_data)**2 + (f * predicted)**2``.

    The noise model ``likelihoods.md`` §5 names (X-1). ``f`` is the fractional
    uncertainty on the *model*: the honest way to say "I do not believe this
    forward model to better than 10 %" without inventing a family or inflating
    the data's own error bars, which would say something different.

    Parameters
    ----------
    f
        The model-uncertainty fraction. A frozen ``scipy.stats`` distribution
        to fit it — the bijection is read off the prior's own support — or a
        number to hold it fixed.
    scale
        Optional multiplier on the *data's* uncertainties — "the catalogue
        underestimates its error bars". Inherited from
        :class:`~ampere.core.IndependentNoise`; note it multiplies
        ``sigma_data`` only, inside the quadrature, exactly as §5's formula
        writes it.
    jitter
        Optional noise floor added in quadrature to ``sigma_data``, also
        inherited.

    Notes
    -----
    ``sigma`` never assumes ``predicted`` and ``observed`` have the same
    length: ``Likelihood.pointwise_log_prob`` calls it once per retained
    sample, with a length-one prediction against the full observed container.
    Both arrays are indexed by ``retain``, which is what makes that work.
    """

    def __init__(self, f: Any, *, scale: Any = None, jitter: Any = None) -> None:
        super().__init__(scale=scale, jitter=jitter)
        self.register_parameter(as_parameter("f", f))

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray:
        base = super().sigma(observed, retain, values, predicted=predicted)
        return _inflate(base, _fraction(self, values), predicted, type(self).__name__)


class FractionalModelGPNoise(GaussianProcessNoise):
    """A misspecification GP *and* a fractional model uncertainty.

    ``likelihoods.md`` §5's composition, which is a subclass overriding only
    ``sigma``: the kernel, the solver and the ``noise_params`` assembly are all
    inherited, and the marginal likelihood is
    ``N(0, K + diag(sigma_data**2 + (f*mu)**2))`` with no new mathematics.

    The two terms say different things and are not redundant. The GP models
    *correlated* structure the forward model gets wrong — a missing feature, a
    continuum slope — while ``f`` models an uncorrelated, amplitude-proportional
    error. Fitting both lets the data decide how much of the misfit is
    smooth.

    Parameters
    ----------
    kernel
        The covariance kernel. Its hyperparameters are adopted as this noise
        model's own parameters, as for :class:`~ampere.core.GaussianProcessNoise`.
    solver
        GP solve strategy; defaults to :class:`~ampere.core.DenseGP`.
    f
        The model-uncertainty fraction, keyword-only so it can never be
        confused with the solver. A prior to fit it, a number to hold it fixed.
    scale, jitter
        As for :class:`~ampere.core.GaussianProcessNoise`, applied to the
        data's uncertainties inside the quadrature.
    """

    def __init__(
        self,
        kernel: Kernel,
        solver: GPSolver | None = None,
        *,
        f: Any,
        scale: Any = None,
        jitter: Any = None,
    ) -> None:
        super().__init__(kernel, solver, scale=scale, jitter=jitter)
        self.register_parameter(as_parameter("f", f))

    def sigma(
        self,
        observed: FunctionSamples,
        retain: np.ndarray,
        values: Mapping[str, Any],
        *,
        predicted: np.ndarray | None = None,
    ) -> np.ndarray:
        base = super().sigma(observed, retain, values, predicted=predicted)
        return _inflate(base, _fraction(self, values), predicted, type(self).__name__)
