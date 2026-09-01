"""A user-defined instrument, written as if it lived in somebody else's package.

``DEVELOPMENT_PLAN.md`` §4.3 makes user extensibility a first-class requirement
and states the acceptance criterion for it: "the test suite includes one such
extension written *out-of-tree* (imported as if third-party) to prove the
interfaces suffice". This module is that extension.

The rules it holds itself to, which is the whole point of it:

* it imports **only** ``ampere.core``'s public API — no private helper, no
  submodule path, nothing that a released ampere would not export;
* it is not part of the ``ampere`` package and is not imported by it;
* ``tests/core/test_transform.py`` loads it the way an installed third-party
  package would be loaded, and runs it through an ordinary ``Instrument``.

Between them, the classes below exercise every extension point W1.4 and W1.5
offer: a new container kind (:class:`PolarisationCurve`), a kind-preserving
transformation with a nuisance parameter and a buffer
(:class:`AtmosphericDepolarisation`), a kind-*changing*, many-to-one
transformation that must propagate a mask (:class:`PolarimeterChannels`), a
published requirement, and an instrument factory.

Nothing inside ampere knows any of this exists.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np

from ampere.core import (
    AxisRequirement,
    AxisSpec,
    FunctionSamples,
    Instrument,
    Layout,
    Order,
    Parameter,
    Spectrum,
    Transformation,
    propagate_mask,
)


class PolarisationCurve(FunctionSamples):
    """Degree of linear polarisation against strictly increasing wavelength.

    A user-defined container kind, per ``results_schema.md`` §13: three class
    attributes, and everything else — validation, masks, units, ``with_values``,
    equality — is inherited.
    """

    AXES = (
        AxisSpec(
            "spectral_axis",
            physical_types=("length", "frequency"),
            order=Order.STRICTLY_INCREASING,
        ),
    )
    LAYOUT = Layout.POINTS
    ALLOW_COMPLEX = False

    def __init__(self, spectral_axis, polarisation, **kwargs) -> None:
        super().__init__({"spectral_axis": spectral_axis}, polarisation, **kwargs)

    @property
    def spectral_axis(self):
        """The wavelength axis."""
        return self.axis("spectral_axis")


class AtmosphericDepolarisation(Transformation):
    """Multiply a spectrum by ``exp(-tau * airmass)``.

    Kind-preserving, so ``PRODUCES`` stays at its default, and grid-independent,
    so a negotiated change of the model's grid does not invalidate it. ``tau``
    is a buffer (a measured optical depth — nobody would put a prior on it);
    ``airmass`` is a nuisance parameter of the instrument, declared exactly as
    a model parameter is, and promotable to a free one without touching
    ``apply``.
    """

    ACCEPTS = (Spectrum,)

    def __init__(self, optical_depth: float, **kwargs) -> None:
        super().__init__(**kwargs)
        self.register_buffer("tau", float(optical_depth), description="grey optical depth")
        self.register_parameter(
            Parameter(
                "airmass",
                None,
                value=1.0,
                fixed=True,
                description="sec(z) at the time of observation",
            )
        )

    def apply(self, samples, values):
        context = self.context(values)
        return samples.with_values(samples.values * np.exp(-context["tau"] * context["airmass"]))


class PolarimeterChannels(Transformation):
    """Bin a spectrum into the polarimeter's broad channels.

    Kind-*changing* (``Spectrum`` in, :class:`PolarisationCurve` out) and
    many-to-one, so it is the case the mask rule exists for: an output channel
    is masked if any model sample falling in it is masked.

    Carries one free nuisance parameter, an additive instrumental polarisation.
    """

    ACCEPTS = (Spectrum,)
    PRODUCES = PolarisationCurve

    def __init__(self, edges, *, prior=None, **kwargs) -> None:
        super().__init__(**kwargs)
        edges = np.asarray(edges, dtype=float)
        if edges.ndim != 1 or edges.size < 2 or np.any(np.diff(edges) <= 0.0):
            raise ValueError("channel edges must be strictly increasing, and there must be >= 2")
        self.register_buffer("edges", edges, unit=u.micron)
        self.register_parameter(
            Parameter(
                "instrumental",
                prior,
                value=0.0,
                fixed=prior is None,
                description="additive instrumental polarisation",
            )
        )

    @property
    def centres(self) -> np.ndarray:
        """Channel centre wavelengths, in microns."""
        edges = self.buffers["edges"].array
        return 0.5 * (edges[:-1] + edges[1:])

    def weights(self, wavelength: np.ndarray) -> np.ndarray:
        """The ``(n_channels, n_samples)`` averaging matrix for *wavelength*."""
        edges = self.buffers["edges"].array
        inside = (wavelength[None, :] >= edges[:-1, None]) & (wavelength[None, :] < edges[1:, None])
        counts = inside.sum(axis=1, keepdims=True)
        return inside / np.maximum(counts, 1)

    def requirements(self):
        edges = self.buffers["edges"].array
        return (
            AxisRequirement(
                "spectral_axis",
                unit=u.um,
                intervals=(float(edges[0]), float(edges[-1])),
                max_step=float(np.diff(edges).min()) / 4.0,
                source="polarimeter channels",
            ),
        )

    def apply(self, samples, values):
        context = self.context(values)
        weights = self.weights(samples.spectral_axis.values)
        binned = weights @ samples.values + context["instrumental"]
        return PolarisationCurve(
            self.centres * u.micron,
            binned,
            mask=propagate_mask(samples, weights),
        )


def build_polarimeter(
    optical_depth: float,
    edges,
    *,
    channel: str = "polarisation",
    prior=None,
) -> Instrument:
    """A two-step polarimeter, ready to bind to *channel* of a ``ModelResult``."""
    return Instrument(
        [
            AtmosphericDepolarisation(optical_depth),
            PolarimeterChannels(edges, prior=prior),
        ],
        channel=channel,
        label="polarimeter",
    )
