"""Synthetic data for the SED composition example: push the truth through the
*negotiated* grid, not the model's own one.

Memo §7.1 finding 1 (``docs/design/phase4_placement_memo.md``): applying an
instrument to a model that has not adopted the negotiated grid is refused by
name. :class:`~ampere.backends.reference.SyntheticPhotometry` asks its
compiled input for the exact points it tabulated its response on
(:meth:`~ampere.core.results_schema.Axis.locate`), and a model still on its
own fallback grid does not have them — the lookup raises
:class:`~ampere.core.exceptions.SchemaError` naming a "negotiation defect
rather than a usage error", because the grid was meant to satisfy exactly the
requirement the step published. That is correct: it is what stops one
instrument's tabulation silently reading someone else's flux. But it means a
synthetic-data generator cannot just call the model on its declared grid and
push the result through an instrument — it has to do what
:class:`~ampere.core.dataset.FittingProblem` does in its own constructor:
:func:`~ampere.core.negotiate`, then :meth:`~ampere.core.transform.Model.compile_for`,
*then* evaluate.

:func:`synthetic_data` is that dance, generic over which backend built the
model and the instruments (it only uses the neutral :class:`~ampere.core.transform.Model`
and :class:`~ampere.core.transform.Instrument` surface) — see
:mod:`.sed_composition` for the backend-specific pieces themselves.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from ampere.core import Instrument, Model, PhotometricPoints, Spectrum, negotiate

__all__ = [
    "CALIBRATION_TRUTH",
    "PHOTOMETRY_FRACTIONAL_NOISE",
    "SEED",
    "SPECTRUM_FRACTIONAL_NOISE",
    "TRUTH",
    "synthetic_data",
]

#: Reproducible everywhere this example is run: the noise draw, and (through
#: ``FittingProblem(..., seed=SEED)``) the sampler's own streams.
SEED = 20260913

#: The injected physical truth — a cool, optically thin dust greybody.
#: ``scale`` is dimensionless (it is what turns :math:`B_\\nu(T)`, in Jy/sr,
#: into an observed flux density — the source's solid angle, in effect), so
#: its natural size is tiny; the flux it produces is not — see
#: :mod:`.sed_composition`'s docstring for the numbers.
TRUTH: dict[str, float] = {"temperature": 180.0, "beta": 1.6, "scale": 1.0e-15}

#: The spectrum's calibration nuisance parameter (``irs.instrument.calibration_scale.scale``
#: once bound into a problem) — a deliberate 2 % miscalibration relative to the
#: catalogue, which the fit is expected to recover.
CALIBRATION_TRUTH = 1.02

#: Fractional 1-sigma noise, applied to each instrument's own truth values.
SPECTRUM_FRACTIONAL_NOISE = 0.03
PHOTOMETRY_FRACTIONAL_NOISE = 0.05


def synthetic_data(
    model: Model,
    spectrograph: Instrument,
    camera: Instrument,
    *,
    seed: int = SEED,
) -> tuple[Spectrum, PhotometricPoints]:
    """Noisy observed containers for *spectrograph* and *camera*, from ``TRUTH``.

    Negotiates the two instruments' requirements, compiles *model* onto the
    union grid, evaluates it once at ``TRUTH``, and pushes the result through
    each instrument — exactly what :class:`~ampere.core.dataset.FittingProblem`
    does when it is built, run here by hand because there is no fitting
    problem yet: this *produces* the data one goes into. ``model`` is mutated
    in place (:meth:`~ampere.core.transform.Model.compile_for` adopts the
    negotiated grid as its template) and handed back to the caller unchanged
    in every other respect, so it can be passed on to
    :class:`~ampere.core.dataset.FittingProblem` as-is — its constructor calls
    ``compile_for`` again, on the same requirements, which is idempotent.

    Parameters
    ----------
    model
        The (uncompiled) model. Compiled onto the union grid as a side effect.
    spectrograph, camera
        The two instruments, already bound to the model's channel with
        distinct labels (:class:`~ampere.core.dataset.DatasetCollection`
        requires that once both are wrapped in datasets — see
        :mod:`.sed_composition`).
    seed
        Seeds the noise draw. The fitting problem built from the result gets
        its own seed independently — this is only about the data.

    Returns
    -------
    tuple
        ``(observed_spectrum, observed_photometry)``, ready to bind into
        :class:`~ampere.core.dataset.Dataset`\\ s.
    """
    rng = np.random.default_rng(seed)
    requirements = negotiate([spectrograph, camera])
    compiled = model.compile_for(requirements)
    truth = compiled(**TRUTH)

    # calibration_scale.scale is the spectrograph's own step's name (relative
    # to the instrument, not yet qualified by a dataset label — that
    # qualification only exists once this instrument is bound into a
    # Dataset).
    spectrum_truth = spectrograph(truth, {"calibration_scale.scale": CALIBRATION_TRUTH})
    photometry_truth = camera(truth)

    spectrum_sigma = SPECTRUM_FRACTIONAL_NOISE * spectrum_truth.values
    photometry_sigma = PHOTOMETRY_FRACTIONAL_NOISE * photometry_truth.values

    observed_spectrum = Spectrum(
        spectrum_truth.spectral_axis.values * u.um,
        (spectrum_truth.values + rng.normal(0.0, spectrum_sigma)) * u.Jy,
        uncertainty=spectrum_sigma * u.Jy,
    )
    observed_photometry = PhotometricPoints(
        photometry_truth.filters,
        photometry_truth.spectral_axis.values * u.um,
        (photometry_truth.values + rng.normal(0.0, photometry_sigma)) * u.Jy,
        uncertainty=photometry_sigma * u.Jy,
    )
    return observed_spectrum, observed_photometry
