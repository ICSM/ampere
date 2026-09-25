"""The synthetic interferometric source both native backends' suites are held to.

One module rather than two copies, because the *claim* the torch and jax
interferometry suites make is the same claim — "this backend's arithmetic is the
reference backend's, and its gradient is real" — and a claim written twice is a
claim that eventually differs. Nothing here imports a backend: every function
takes the backend's ``interferometry`` module as its first argument, exactly as
``tests/conformance`` takes a fixture, so the same body runs under ``-e torch``
and under ``-e jax``.

The geometry is ``tests/conformance/test_interferometry.py``'s, deliberately:
the same four telescopes, the same binary, the same field of view. A row here
that disagrees with a row there is then a disagreement about arithmetic rather
than about which source was observed, and the conformance battery's closed-form
oracles apply unchanged.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    ClosurePhases,
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    FittingProblem,
    Instrument,
    Likelihood,
    VisibilitySet,
    VonMisesFamily,
    negotiate,
)

#: The observing wavelength, micron. Monochromatic: the reference-path ``Image``
#: channel is achromatic and a chromatic sky is a ``Cube`` channel later.
WAVELENGTH = 2.2

#: Four telescopes, metres east and north of the array centre — the conformance
#: battery's own positions, chosen there so that every closure phase of
#: :data:`BINARY` is well away from zero (about -0.8 rad) and a sign error
#: cannot pass.
STATIONS = np.array([[0.0, 0.0], [-32.4, 39.8], [6.5, -30.1], [-14.7, -48.9]])

#: Hour angles the array is observed at. Four, so that the uv plane is filled
#: rather than sampled at six points: with one snapshot a separation and a flux
#: ratio are poorly determined and "NUTS recovers the binary" would be a claim
#: about the prior.
HOUR_ANGLES = np.linspace(-0.6, 0.6, 4)

#: The field of view the instrument is sensitive to, mas.
FIELD_OF_VIEW = 40.0

#: Nyquist oversampling. Four rather than the battery's ten: the rows here are
#: about gradients and agreement between two backends computing the same sum,
#: not about the quadrature's absolute accuracy, and the grid cost is quadratic.
OVERSAMPLING = 4.0

#: The synthetic binary, at its true values.
BINARY: dict[str, Any] = {
    "separation": 12.0,
    "position_angle": 0.7,
    "flux_ratio": 0.42,
    "flux": 1.7,
}

#: Both components' width, mas. A **buffer**: it is there to make the pair
#: band-limited, and one would not put a prior on a numerical device.
COMPONENT_FWHM = 2.0

#: Measurement uncertainties: Jy on a correlated flux, radians on a closure
#: phase. The von Mises concentration is ``1/sigma**2``, so 0.05 rad is a
#: well-measured triangle and the normalisation matters.
SIGMA_VISIBILITY = 0.02
SIGMA_CLOSURE = 0.05

#: The free parameters of the fitted problem, and their true values.
TRUTH: dict[str, float] = {
    "model.separation": BINARY["separation"],
    "model.flux_ratio": BINARY["flux_ratio"],
}

#: The seed every synthetic dataset here is drawn at.
SEED = 20260913


def _baselines() -> dict[tuple[int, int], np.ndarray]:
    """``b_ij = (r_j - r_i) / lambda`` in wavelengths, for every pair ``i < j``."""
    metres = WAVELENGTH * 1e-6
    return {
        (i, j): (STATIONS[j] - STATIONS[i]) / metres
        for i in range(len(STATIONS))
        for j in range(i + 1, len(STATIONS))
    }


def _rotate(u_pts: np.ndarray, v_pts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The same baselines at every hour angle, concatenated."""
    turned_u = np.concatenate([u_pts * math.cos(h) - v_pts * math.sin(h) for h in HOUR_ANGLES])
    turned_v = np.concatenate([u_pts * math.sin(h) + v_pts * math.cos(h) for h in HOUR_ANGLES])
    return turned_u, turned_v


def visibility_coverage() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """``(u, v, lambda)`` for the array's six baselines at every hour angle."""
    table = _baselines()
    pairs = sorted(table)
    stacked = np.array([table[pair] for pair in pairs])
    turned_u, turned_v = _rotate(stacked[:, 0], stacked[:, 1])
    return turned_u, turned_v, np.full(turned_u.size, WAVELENGTH)


def triangle_coverage() -> tuple[np.ndarray, ...]:
    """``(u1, v1, u2, v2, lambda)`` in the canonical ordering, at every hour angle.

    Telescopes ``i < j < k``; the stored baselines are ``ij`` and ``jk``, and
    ``ki`` is implied as their negated sum. Written out here rather than taken
    from a backend, so a row compares an ordering against the rule rather than
    against the code that implements it.
    """
    table = _baselines()
    triples = [
        (i, j, k)
        for i in range(len(STATIONS))
        for j in range(i + 1, len(STATIONS))
        for k in range(j + 1, len(STATIONS))
    ]
    first = np.array([table[(i, j)] for i, j, _ in triples])
    second = np.array([table[(j, k)] for _, j, k in triples])
    u1, v1 = _rotate(first[:, 0], first[:, 1])
    u2, v2 = _rotate(second[:, 0], second[:, 1])
    return u1, v1, u2, v2, np.full(u1.size, WAVELENGTH)


def visibilities(values: np.ndarray | None = None) -> VisibilitySet:
    """A ``VisibilitySet`` on the array's coverage, with a uniform sigma."""
    u_pts, v_pts, waves = visibility_coverage()
    filled = np.zeros(u_pts.size, dtype=np.complex128) if values is None else values
    return VisibilitySet(
        u_pts,
        v_pts,
        waves * u.micron,
        filled * u.Jy,
        uncertainty=np.full(u_pts.size, SIGMA_VISIBILITY) * u.Jy,
    )


def closure_phases(values: np.ndarray | None = None) -> ClosurePhases:
    """A ``ClosurePhases`` on the array's triangles, with a uniform sigma."""
    u1, v1, u2, v2, waves = triangle_coverage()
    filled = np.zeros(u1.size) if values is None else values
    return ClosurePhases(
        u1,
        v1,
        u2,
        v2,
        waves * u.micron,
        filled * u.rad,
        uncertainty=np.full(u1.size, SIGMA_CLOSURE) * u.rad,
    )


def chain(itf: Any, observed: Any, label: str) -> Instrument:
    """The instrument for one of the two observables, from *itf*'s own steps.

    A ``VisibilitySet`` needs the Fourier step alone; a ``ClosurePhases``
    needs it followed by the three-to-one step, on the coverage
    ``FourierSample.from_observed`` lays out in the canonical order.
    """
    sample = itf.FourierSample.from_observed(
        observed, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
    )
    steps = [sample] if isinstance(observed, VisibilitySet) else [sample, itf.ClosurePhase()]
    return Instrument(steps, channel="sky", label=label)


def truth_model(itf: Any, **overrides: Any) -> Any:
    """The binary at its true values, on its own square grid."""
    return itf.Binary.on_field(
        FIELD_OF_VIEW * u.mas,
        4,
        channels="sky",
        component_fwhm=COMPONENT_FWHM,
        **{**BINARY, **overrides},
    )


def synthetic(itf: Any, *, seed: int = SEED) -> tuple[Any, VisibilitySet, ClosurePhases]:
    """Negotiate, compile, evaluate and observe: the data a fit is given.

    Returned with the compiled truth model, because its negotiated ``(x, y)``
    grid is what the *fitted* model must be built on — a fit whose model grid
    differed from the one the data were made on would be measuring the
    discretisation rather than the source.
    """
    vis_instrument = chain(itf, visibilities(), "vis")
    t3_instrument = chain(itf, closure_phases(), "t3")
    compiled = truth_model(itf).compile_for(negotiate([vis_instrument, t3_instrument]))
    result = compiled.evaluate()
    clean_visibility = np.asarray(vis_instrument(result).values)
    clean_phase = np.asarray(t3_instrument(result).values)
    rng = np.random.default_rng(seed)
    # Circular complex noise on the visibilities, von Mises on the phases: the
    # two families the fit scores them with, so the data really are a draw from
    # the likelihood rather than from something near it.
    noisy_visibility = clean_visibility + SIGMA_VISIBILITY * (
        rng.standard_normal(clean_visibility.size) + 1j * rng.standard_normal(clean_visibility.size)
    )
    noisy_phase = np.angle(np.exp(1j * rng.vonmises(clean_phase, 1.0 / SIGMA_CLOSURE**2)))
    return compiled, visibilities(noisy_visibility), closure_phases(noisy_phase)


def fitted_model(itf: Any, compiled: Any) -> Any:
    """The binary with ``separation`` and ``flux_ratio`` free, on the negotiated grid."""
    return itf.Binary(
        compiled.buffers["x"].value * u.mas,
        compiled.buffers["y"].value * u.mas,
        channels="sky",
        component_fwhm=COMPONENT_FWHM,
        separation=st.uniform(6.0, 14.0),
        position_angle=BINARY["position_angle"],
        flux_ratio=st.uniform(0.1, 0.7),
        flux=BINARY["flux"],
    )


def two_dataset_problem(
    backend: Any, itf: Any, *, gp: bool = False, seed: int = SEED
) -> FittingProblem:
    """Visibilities and closure phases from one sky, on one backend.

    The composition Phase 4 exists to prove: two instruments bind the **same**
    ``sky`` channel by name, publish requirements on the same two axes, and
    :func:`~ampere.core.negotiate` unions them into one image grid the model
    builds once per draw.

    Parameters
    ----------
    backend
        The backend package (``ampere.backends.torch`` or ``.jax``), for its
        noise models, kernel and solver — every part of a native problem must
        declare the same backend or composition is refused.
    itf
        That backend's ``interferometry`` module.
    gp
        Whether the visibilities carry the flexible likelihood: a circular
        complex Gaussian process over ``(u, v)`` (W4.2). The closure phases
        keep independent von Mises noise either way — a GP on a *wrapped*
        observable is a latent-variable model, fitted on the native path only
        since W5.1 (``phase_problem`` in ``examples.interferometry.study`` is
        the latent arm); this fixture's arms stay the W4.4 ones.
    seed
        The noise realisation.
    """
    compiled, observed_visibility, observed_phase = synthetic(itf, seed=seed)
    vis_instrument = chain(itf, observed_visibility, "vis")
    t3_instrument = chain(itf, observed_phase, "t3")
    if gp:
        # ``axes=("u", "v")`` is the selector W4.5 put on every kernel, and
        # binding it is one of the two native-path fixes W4.2 asked this item to
        # guard: a kernel that ignored the selector here would quietly build its
        # covariance from the wrong column and still look healthy.
        kernel = backend.Matern32(
            st.loguniform(1e-3, 0.2), st.loguniform(1e6, 1e8), axes=("u", "v")
        )
        visibility_noise: Any = backend.GaussianProcessNoise(kernel, backend.DenseGP())
    else:
        visibility_noise = backend.IndependentNoise()
    datasets = DatasetCollection(
        {
            "vis": Dataset(
                observed_visibility,
                vis_instrument,
                likelihood=Likelihood(ComplexGaussianFamily(), visibility_noise),
                label="vis",
            ),
            "t3": Dataset(
                observed_phase,
                t3_instrument,
                likelihood=Likelihood(VonMisesFamily(), backend.IndependentNoise()),
                label="t3",
            ),
        }
    )
    return FittingProblem(fitted_model(itf, compiled), datasets, seed=SEED)


def reference_point(*, gp: bool = False) -> dict[str, float]:
    """A θ inside the support of every prior: the truth, plus the GP's own."""
    point = dict(TRUTH)
    if gp:
        point["vis.likelihood.amplitude"] = 0.02
        point["vis.likelihood.length_scale"] = 1.0e7
    return point
