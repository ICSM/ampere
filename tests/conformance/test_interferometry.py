"""Interferometry: the rows Phase 4's proof modality owes (W4.1).

One column per registered backend, like every other file here, and no test body
names one. A backend that has not written the interferometric vocabulary
declares ``BackendCapabilities.interferometry = False`` and every row skips
with a reason naming what is owed — the shape ``tests/conformance/README.md``
§3 asks for ("a row whose second implementation does not exist yet is a
skeleton, not an absence").

What these rows are for, in the order the work item puts them:

* the direct Fourier transform of each source model against its **closed
  form** in :mod:`tests.conformance.oracles` — not against the backend's own
  analytic model, which would be two ampere implementations agreeing;
* the closure phases of a binary against the closed form, **including the
  sign**, which is the convention that would otherwise fit a mirrored source
  just as well;
* each smearing step against brute-force fine sampling, which is what says the
  five-node quadrature is enough and that ``configure_from`` handed the
  Fourier step the right extra samples;
* mask propagation through the three-to-one step and through an average;
* two datasets on one model channel evaluating the model **once** per draw
  (``inference.md`` §8), with the pairing between the two visible in the
  composition (the plan's design horizon (h));
* ``simulate(observe=True)`` on both kinds.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    ClosurePhases,
    ComplexGaussianFamily,
    CompositionError,
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    LikelihoodError,
    RiceFamily,
    TransformationError,
    VisibilitySet,
    VonMisesFamily,
    negotiate,
)

from .oracles import (
    binary_closure_phase,
    binary_visibility,
    gaussian_visibility,
    uniform_disc_visibility,
)
from .protocol import ConformanceBackend, InterferometryPieces, Tolerances

#: The observing wavelength every row here works at, micron. Monochromatic,
#: because the spectral axis is a *coordinate* question (W4.1's amendment) and
#: the reference-path image channel is achromatic.
WAVELENGTH = 2.2

#: A four-telescope array, in metres east and north of the array centre. Four
#: telescopes are the smallest number that gives more triangles than a single
#: closure phase (six baselines, four triangles), which is what makes the
#: mask-propagation row say something: a masked baseline must take out the two
#: triangles that use it and leave the other two alone. The particular
#: positions were chosen so that every closure phase of :data:`BINARY` is well
#: away from zero (about -0.8 rad): a triangle whose closure phase is near
#: zero would agree with a sign error as readily as with the convention.
STATIONS = np.array([[0.0, 0.0], [-32.4, 39.8], [6.5, -30.1], [-14.7, -48.9]])

#: The field of view the instrument is sensitive to, mas.
FIELD_OF_VIEW = 40.0

#: How far past the array's own Nyquist rate the image is sampled. Ten is
#: enough for a band-limited source of a few mas on this array: the two
#: band-limited rows land at machine precision, ten orders inside
#: ``tolerances.cross_solver``, and the row that measures the convergence
#: shows where it comes from. At the array's own Nyquist rate the same rows
#: are wrong by tens of per cent, which is the point ``FourierSample.requirements``
#: makes about what ``oversampling = 1`` is and is not.
OVERSAMPLING = 10.0

#: The binary the phase rows are proved on. Compact enough to be resolved by
#: the array (so the closure phases are not all near zero) and broad enough in
#: its components to be band-limited on a grid this size.
BINARY = {
    "separation": 12.0,
    "position_angle": 0.7,
    "flux_ratio": 0.42,
    "flux": 1.7,
    "component_fwhm": 2.0,
}
GAUSSIAN = {"fwhm": 4.0, "flux": 1.4}
#: The two free parameters of the composed problem, at their true values.
TRUTH = {
    "model.separation": BINARY["separation"],
    "model.flux_ratio": BINARY["flux_ratio"],
}
DISC = {"diameter": 5.0, "flux": 1.1}


# ---------------------------------------------------------------------------
# Deterministic geometry
# ---------------------------------------------------------------------------


def baselines() -> dict[tuple[int, int], np.ndarray]:
    """``b_ij = (r_j - r_i) / lambda`` in wavelengths, for every pair ``i < j``.

    The sign convention the canonical ordering names: the baseline runs from
    the lower-numbered station to the higher one. Every row that cares about a
    phase depends on it.
    """
    metres_per_wavelength = WAVELENGTH * 1e-6
    return {
        (i, j): (STATIONS[j] - STATIONS[i]) / metres_per_wavelength
        for i in range(len(STATIONS))
        for j in range(i + 1, len(STATIONS))
    }


def visibility_coverage() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """``(u, v, lambda, labels)`` for every baseline of the array."""
    pairs = sorted(baselines())
    table = baselines()
    u_pts = np.array([table[pair][0] for pair in pairs])
    v_pts = np.array([table[pair][1] for pair in pairs])
    waves = np.full(len(pairs), WAVELENGTH)
    labels = np.array([f"T{i}-T{j}" for i, j in pairs])
    return u_pts, v_pts, waves, labels


def triangle_coverage() -> tuple[np.ndarray, ...]:
    """``(u1, v1, u2, v2, lambda, labels)`` in the canonical ordering.

    Telescopes ``i < j < k``; the stored baselines are ``ij`` and ``jk``, and
    ``ki`` is implied as their negated sum. Written out here rather than taken
    from the backend, so that the rows compare an ordering against the rule
    rather than against the code that implements it.
    """
    table = baselines()
    triples = [
        (i, j, k)
        for i in range(len(STATIONS))
        for j in range(i + 1, len(STATIONS))
        for k in range(j + 1, len(STATIONS))
    ]
    first = np.array([table[(i, j)] for i, j, _ in triples])
    second = np.array([table[(j, k)] for _, j, k in triples])
    labels = np.array([f"T{i}-T{j}-T{k}" for i, j, k in triples])
    return (
        first[:, 0],
        first[:, 1],
        second[:, 0],
        second[:, 1],
        np.full(len(triples), WAVELENGTH),
        labels,
    )


def observed_visibilities(values: np.ndarray | None = None) -> VisibilitySet:
    """A ``VisibilitySet`` on the array's coverage, with a uniform sigma."""
    u_pts, v_pts, waves, labels = visibility_coverage()
    filled = np.zeros(u_pts.size, dtype=np.complex128) if values is None else values
    return VisibilitySet(
        u_pts,
        v_pts,
        waves * u.micron,
        filled * u.Jy,
        uncertainty=np.full(u_pts.size, 0.02) * u.Jy,
        extra_coords={"baseline": labels},
    )


def observed_closure_phases(values: np.ndarray | None = None) -> ClosurePhases:
    """A ``ClosurePhases`` on the array's triangles, with a uniform sigma."""
    u1, v1, u2, v2, waves, labels = triangle_coverage()
    filled = np.zeros(u1.size) if values is None else values
    return ClosurePhases(
        u1,
        v1,
        u2,
        v2,
        waves * u.micron,
        filled * u.rad,
        uncertainty=np.full(u1.size, 0.05) * u.rad,
        extra_coords={"triangle": labels},
    )


def pieces_or_skip(backend: ConformanceBackend) -> InterferometryPieces:
    """This backend's interferometric classes, or a skip naming what is owed."""
    if not backend.capabilities.interferometry:
        pytest.skip(
            f"{backend.name} declares no interferometric vocabulary "
            f"(BackendCapabilities.interferometry), so the Fourier, closure-phase and smearing "
            f"steps and the three source models it would need are not there. W4.3 owes the "
            f"native twins."
        )
    return backend.interferometry()


def counting(model_class: type) -> type:
    """*model_class* with an ``evaluations`` counter, for the once-per-draw row.

    ``inference.md`` §18's reason for :class:`CountingModel` applies here for
    the same reason it applies there: there is no way to assert "the model was
    evaluated once" from outside without the model saying so. Subclassing
    whatever the fixture supplies keeps the row backend-agnostic — and keeps
    the ``BACKEND`` declaration, which a problem composed from it is checked
    against.
    """

    class Counting(model_class):  # type: ignore[valid-type, misc]
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.evaluations = 0

        def evaluate(self, **values: Any) -> Any:
            self.evaluations += 1
            return super().evaluate(**values)

    Counting.__name__ = f"Counting{model_class.__name__}"
    Counting.__qualname__ = Counting.__name__
    return Counting


def image_model(pieces: InterferometryPieces, which: str, **overrides: Any) -> Any:
    """One of the three image-emitting source models, on its default grid."""
    table = {
        "gaussian": (pieces.gaussian_source, GAUSSIAN),
        "binary": (pieces.binary, BINARY),
        "disc": (pieces.uniform_disc, DISC),
    }
    model_class, source = table[which]
    return model_class.on_field(FIELD_OF_VIEW * u.mas, 4, channels="sky", **{**source, **overrides})


def oracle_visibilities(which: str, u_pts: np.ndarray, v_pts: np.ndarray) -> np.ndarray:
    """The closed form for one of the three sources, from :mod:`.oracles`."""
    if which == "gaussian":
        return gaussian_visibility(u_pts, v_pts, **GAUSSIAN)
    if which == "binary":
        return binary_visibility(u_pts, v_pts, **BINARY)
    return uniform_disc_visibility(u_pts, v_pts, **DISC)


def transformed(
    pieces: InterferometryPieces,
    model: Any,
    *,
    steps: tuple[Any, ...] = (),
    observed: Any = None,
    oversampling: float = OVERSAMPLING,
) -> Any:
    """Negotiate, compile and run one instrument over *model*, in one call."""
    container = observed_visibilities() if observed is None else observed
    instrument = Instrument(
        [
            pieces.fourier_sample.from_observed(
                container, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=oversampling
            ),
            *steps,
        ],
        channel="sky",
        label="array",
    )
    requirements = negotiate([instrument])
    compiled = model.compile_for(requirements)
    return instrument(compiled.evaluate())


# ---------------------------------------------------------------------------
# The Fourier step
# ---------------------------------------------------------------------------


class TestFourierSample:
    """``Image -> VisibilitySet``: the kind-changing step in the middle of a chain."""

    @pytest.mark.parametrize("which", ["gaussian", "binary"])
    def test_the_transform_of_a_band_limited_image_matches_the_closed_form(
        self, backend: ConformanceBackend, tolerances: Tolerances, which: str
    ) -> None:
        """The direct DFT against van Cittert-Zernike, at ``cross_solver``.

        Band-limited sources only, and that restriction is physics rather than
        implementation: a sum over a finite grid reproduces the Fourier
        integral exactly in the limit that the image has no power beyond the
        grid's own Nyquist frequency, and a Gaussian's power falls off
        exponentially while a sharp edge's falls off as a power. See
        :meth:`test_a_sharp_edged_image_converges_rather_than_agreeing` for
        the other case, which is measured rather than asserted away.
        """
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, _, _ = visibility_coverage()
        got = transformed(pieces, image_model(pieces, which))
        expected = oracle_visibilities(which, u_pts, v_pts)
        assert got.unit == u.Jy
        error = np.max(np.abs(backend.to_numpy(got.values) - expected)) / np.max(np.abs(expected))
        assert error < tolerances.cross_solver

    @pytest.mark.parametrize("which", ["gaussian", "binary", "disc"])
    def test_the_analytic_route_matches_the_same_closed_form(
        self, backend: ConformanceBackend, tolerances: Tolerances, which: str
    ) -> None:
        """The no-steps route: a model emitting visibilities directly.

        ``interferometry.md`` §1 — "both routes are supported and neither is
        privileged". This one has no quadrature in it, so all three sources
        are held to the same tolerance, the uniform disc included.
        """
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, waves, _ = visibility_coverage()
        table = {
            "gaussian": (pieces.gaussian_source_visibilities, GAUSSIAN),
            "binary": (pieces.binary_visibilities, BINARY),
            "disc": (pieces.uniform_disc_visibilities, DISC),
        }
        model_class, source = table[which]
        model = model_class(u_pts, v_pts, waves * u.micron, channels="vis", **source)
        instrument = Instrument([], channel="vis", input_kind=VisibilitySet, label="analytic")
        got = instrument(model.evaluate())
        expected = oracle_visibilities(which, u_pts, v_pts)
        assert np.allclose(
            backend.to_numpy(got.values), expected, rtol=0.0, atol=tolerances.analytic
        )

    def test_a_sharp_edged_image_converges_rather_than_agreeing(
        self, backend: ConformanceBackend
    ) -> None:
        """A uniform disc is not band-limited, and the row says so by measuring it.

        The quadrature error of a rectangle rule over a discontinuous image
        falls roughly as the pixel scale rather than exponentially, so no grid
        this suite can afford reaches ``cross_solver``. Asserting a loose
        tolerance instead would hide the shape of the error; asserting the
        *convergence* is the honest statement, and it is also the one that
        catches a real defect — a transform with a wrong normalisation or a
        wrong sign does not converge at all.
        """
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, _, _ = visibility_coverage()
        expected = oracle_visibilities("disc", u_pts, v_pts)
        scale = np.max(np.abs(expected))
        errors = []
        for oversampling in (2.0, 32.0):
            got = transformed(pieces, image_model(pieces, "disc"), oversampling=oversampling)
            errors.append(float(np.max(np.abs(backend.to_numpy(got.values) - expected)) / scale))
        assert errors[1] < 5.0e-3
        assert errors[1] < errors[0] / 20.0

    def test_the_published_requirement_is_the_field_of_view_at_the_nyquist_step(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``intervals = +-fov/2`` and ``max_step = 1/(2 s u_max)``, per axis.

        Separable, and exactly so rather than conservatively: the DFT kernel
        factorises, so the ``x`` constraint depends on ``u`` alone and the
        ``y`` constraint on ``v`` alone (``interferometry.md`` §4).
        """
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, _, _ = visibility_coverage()
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed_visibilities(),
                    field_of_view=FIELD_OF_VIEW * u.mas,
                    oversampling=OVERSAMPLING,
                )
            ],
            channel="sky",
            label="array",
        )
        published = negotiate([instrument])["sky"]
        mas_per_rad = (1.0 * u.rad).to_value(u.mas)
        for axis, coverage in (("x", u_pts), ("y", v_pts)):
            (low, high, step, power) = published[axis].segments()[0]
            assert (low, high) == (-FIELD_OF_VIEW / 2.0, FIELD_OF_VIEW / 2.0)
            assert power is None
            expected = mas_per_rad / (2.0 * OVERSAMPLING * np.max(np.abs(coverage)))
            assert step == pytest.approx(expected, rel=tolerances.exact)

    def test_a_model_that_cannot_honour_the_pixel_scale_refuses(
        self, backend: ConformanceBackend
    ) -> None:
        """Gap I-4's loud option, on the modality it was written for.

        A model holding its own coarse grid raises rather than aliasing: power
        from beyond the Nyquist limit folded back onto the sampled baselines
        looks like real source structure, not like an error.
        """
        pieces = pieces_or_skip(backend)
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed_visibilities(),
                    field_of_view=FIELD_OF_VIEW * u.mas,
                    oversampling=OVERSAMPLING,
                )
            ],
            channel="sky",
            label="array",
        )
        requirements = negotiate([instrument])
        coarse = image_model(pieces, "gaussian", adopt_grid=False)
        with pytest.raises(CompositionError, match="aliases"):
            coarse.compile_for(requirements)

    def test_the_transform_is_the_same_object_across_draws(
        self, backend: ConformanceBackend
    ) -> None:
        """The compile-once/evaluate-many split: the axes are shared by identity."""
        pieces = pieces_or_skip(backend)
        model = image_model(pieces, "gaussian")
        container = observed_visibilities()
        instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    container, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=2.0
                )
            ],
            channel="sky",
            label="array",
        )
        compiled = model.compile_for(negotiate([instrument]))
        first = instrument(compiled.evaluate())
        second = instrument(compiled.evaluate())
        assert first.axes is second.axes
        # And the coordinates are the observed container's own, bit for bit —
        # gap I-2's rule, which check_alignment depends on.
        assert np.array_equal(first.u.values, container.u.values)
        assert np.array_equal(first.spectral_axis.values, container.spectral_axis.values)


# ---------------------------------------------------------------------------
# Closure phases
# ---------------------------------------------------------------------------


class TestClosurePhase:
    """Three visibilities to one angle, and the sign convention that goes with it."""

    def _phases(self, pieces: InterferometryPieces, **overrides: Any) -> Any:
        observed = observed_closure_phases()
        return transformed(
            pieces,
            image_model(pieces, "binary", **overrides),
            steps=(pieces.closure_phase(),),
            observed=observed,
        )

    def test_the_binary_matches_the_closed_form(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``arg(V_ij V_jk V_ki)`` against the closed form, in radians."""
        pieces = pieces_or_skip(backend)
        u1, v1, u2, v2, _, _ = triangle_coverage()
        got = self._phases(pieces)
        expected = binary_closure_phase(u1, v1, u2, v2, **BINARY)
        assert got.unit == u.rad
        # The binary must actually be resolved, or the row would pass on a
        # model that returned zeros.
        assert np.max(np.abs(expected)) > 0.1
        residual = np.angle(np.exp(1j * (backend.to_numpy(got.values) - expected)))
        assert np.max(np.abs(residual)) < tolerances.cross_solver

    def test_the_sign_convention_is_the_one_the_kind_states(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Mirroring the source negates every closure phase, and nothing else does.

        The check that a sign error cannot hide: with ``arg(V12 V23 V31)`` and
        ``b_ij = r_j - r_i``, a source reflected through the phase centre has
        closure phases of the opposite sign. A model with the *other*
        convention would agree with the closed form on a point-symmetric
        source and disagree here.
        """
        pieces = pieces_or_skip(backend)
        got = self._phases(pieces)
        mirrored = self._phases(pieces, position_angle=BINARY["position_angle"] + math.pi)
        total = backend.to_numpy(got.values) + backend.to_numpy(mirrored.values)
        assert np.max(np.abs(np.angle(np.exp(1j * total)))) < tolerances.cross_solver
        assert np.max(np.abs(backend.to_numpy(got.values))) > 0.1

    def test_a_masked_baseline_masks_every_triangle_that_uses_it(
        self, backend: ConformanceBackend
    ) -> None:
        """The many-to-one rule (``results_schema.md`` §16), on the three-to-one step.

        Masking the first baseline of the *first* triangle must take out that
        triangle and leave the others alone: an output sample is masked if any
        input that influences it is, and no more than that.
        """
        pieces = pieces_or_skip(backend)
        u1, v1, u2, v2, waves, _ = triangle_coverage()
        third = (-(u1 + u2), -(v1 + v2))
        u_pts = np.stack([u1, u2, third[0]], axis=1).reshape(-1)
        v_pts = np.stack([v1, v2, third[1]], axis=1).reshape(-1)
        mask = np.zeros(u_pts.size, dtype=bool)
        mask[0] = True
        mask[3 * 2 + 1] = True
        source = VisibilitySet(
            u_pts,
            v_pts,
            np.repeat(waves, 3) * u.micron,
            np.ones(u_pts.size, dtype=np.complex128) * u.Jy,
            mask=mask,
        )
        got = pieces.closure_phase()(source)
        assert got.mask is not None
        assert got.mask.tolist() == [True, False, True, False]

    def test_coverage_that_does_not_close_is_refused(self, backend: ConformanceBackend) -> None:
        """A triangle whose three baselines do not sum to zero is not a triangle."""
        pieces = pieces_or_skip(backend)
        u1, v1, u2, v2, waves, _ = triangle_coverage()
        third = (-(u1 + u2), -(v1 + v2))
        u_pts = np.stack([u1, u2, third[0]], axis=1).reshape(-1)
        v_pts = np.stack([v1, v2, third[1]], axis=1).reshape(-1)
        u_pts[2] *= 1.05
        broken = VisibilitySet(
            u_pts,
            v_pts,
            np.repeat(waves, 3) * u.micron,
            np.ones(u_pts.size, dtype=np.complex128) * u.Jy,
        )
        with pytest.raises(TransformationError, match="do not close"):
            pieces.closure_phase()(broken)


# ---------------------------------------------------------------------------
# The smearing steps
# ---------------------------------------------------------------------------


#: Resolving power and integration time chosen so each effect is large enough
#: to be worth modelling — a few parts in a thousand and a fifth of the
#: amplitude respectively — rather than lost in the comparison.
RESOLVING_POWER = 5.0
INTEGRATION = 5400.0


def uv_rates() -> tuple[np.ndarray, np.ndarray]:
    """A deterministic uv-track rate per baseline, wavelengths per second."""
    u_pts, v_pts, _, _ = visibility_coverage()
    omega = 7.2921e-5  # rad/s, the Earth's rotation
    return -omega * v_pts, omega * u_pts


class TestSmearing:
    """Bandwidth and time averaging, and the cross-kind ``configure_from`` behind them."""

    def _chain(self, pieces: InterferometryPieces, step: Any) -> tuple[Any, Any]:
        container = observed_visibilities()
        fourier = pieces.fourier_sample.from_observed(
            container, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
        )
        instrument = Instrument([fourier, step], channel="sky", label="array")
        model = image_model(pieces, "gaussian").compile_for(negotiate([instrument]))
        return fourier, instrument(model.evaluate())

    def test_bandwidth_smearing_matches_brute_force_fine_sampling(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Five Gauss-Legendre nodes against a thousand-point average of the closed form."""
        pieces = pieces_or_skip(backend)
        _, got = self._chain(pieces, pieces.bandwidth_smearing(resolving_power=RESOLVING_POWER))
        u_pts, v_pts, _, _ = visibility_coverage()
        offsets = np.linspace(-0.5, 0.5, 1001)
        factor = 1.0 + offsets / RESOLVING_POWER
        brute = np.array(
            [
                np.trapezoid(
                    gaussian_visibility(u_pts[i] / factor, v_pts[i] / factor, **GAUSSIAN),
                    offsets,
                )
                for i in range(u_pts.size)
            ]
        )
        unsmeared = gaussian_visibility(u_pts, v_pts, **GAUSSIAN)
        scale = np.max(np.abs(brute))
        # The effect must be bigger than the tolerance, or the row would pass
        # on a step that did nothing at all.
        assert np.max(np.abs(unsmeared - brute)) / scale > 1.0e-3
        assert (
            np.max(np.abs(backend.to_numpy(got.values) - brute)) / scale < tolerances.cross_solver
        )

    def test_time_smearing_matches_brute_force_fine_sampling(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The same, along the uv track rather than along the spoke."""
        pieces = pieces_or_skip(backend)
        du_dt, dv_dt = uv_rates()
        step = pieces.time_smearing(integration=INTEGRATION * u.s, du_dt=du_dt, dv_dt=dv_dt)
        _, got = self._chain(pieces, step)
        u_pts, v_pts, _, _ = visibility_coverage()
        offsets = np.linspace(-0.5, 0.5, 1001) * INTEGRATION
        brute = np.array(
            [
                np.trapezoid(
                    gaussian_visibility(
                        u_pts[i] + du_dt[i] * offsets, v_pts[i] + dv_dt[i] * offsets, **GAUSSIAN
                    ),
                    offsets / INTEGRATION,
                )
                for i in range(u_pts.size)
            ]
        )
        unsmeared = gaussian_visibility(u_pts, v_pts, **GAUSSIAN)
        scale = np.max(np.abs(brute))
        assert np.max(np.abs(unsmeared - brute)) / scale > 1.0e-3
        assert (
            np.max(np.abs(backend.to_numpy(got.values) - brute)) / scale < tolerances.cross_solver
        )

    def test_the_fourier_step_computes_the_extra_samples_the_smearing_asked_for(
        self, backend: ConformanceBackend
    ) -> None:
        """Gap I-3's mechanism, across a kind change.

        The smearing step cannot publish a requirement — the channel holds an
        ``Image`` — so it tells the step before it instead, and the step
        before it is the one that owns the coverage. Composing two smearing
        steps multiplies the expansion, and both reduce back to the observed
        coverage bit for bit.
        """
        pieces = pieces_or_skip(backend)
        du_dt, dv_dt = uv_rates()
        container = observed_visibilities()
        bare = pieces.fourier_sample.from_observed(container, field_of_view=FIELD_OF_VIEW * u.mas)
        Instrument([bare], channel="sky", label="bare")
        assert bare.expanded_coverage[0].size == container.n_samples

        fourier = pieces.fourier_sample.from_observed(
            container, field_of_view=FIELD_OF_VIEW * u.mas, oversampling=OVERSAMPLING
        )
        instrument = Instrument(
            [
                fourier,
                pieces.bandwidth_smearing(resolving_power=RESOLVING_POWER, nodes=3),
                pieces.time_smearing(integration=INTEGRATION * u.s, du_dt=du_dt, dv_dt=dv_dt),
            ],
            channel="sky",
            label="array",
        )
        assert fourier.expanded_coverage[0].size == container.n_samples * 3 * 5
        model = image_model(pieces, "gaussian").compile_for(negotiate([instrument]))
        got = instrument(model.evaluate())
        assert np.array_equal(got.u.values, container.u.values)
        assert np.array_equal(got.v.values, container.v.values)
        assert np.array_equal(got.spectral_axis.values, container.spectral_axis.values)

    def test_the_extra_samples_widen_the_published_pixel_scale(
        self, backend: ConformanceBackend
    ) -> None:
        """``u_max`` is taken over the *expanded* coverage, not the observed one.

        A bandwidth sub-sample at the blue edge of the channel sits further
        out along the spoke than the recorded point, and an image built for
        the recorded point alone would alias it.
        """
        pieces = pieces_or_skip(backend)
        container = observed_visibilities()
        plain = pieces.fourier_sample.from_observed(container, field_of_view=FIELD_OF_VIEW * u.mas)
        Instrument([plain], channel="sky", label="plain")
        smeared = pieces.fourier_sample.from_observed(
            container, field_of_view=FIELD_OF_VIEW * u.mas
        )
        Instrument(
            [smeared, pieces.bandwidth_smearing(resolving_power=RESOLVING_POWER)],
            channel="sky",
            label="smeared",
        )
        finer = smeared.requirements()[0].segments()[0][2]
        coarser = plain.requirements()[0].segments()[0][2]
        assert finer is not None and coarser is not None
        assert finer < coarser

    def test_a_masked_sub_sample_masks_the_average(self, backend: ConformanceBackend) -> None:
        """An average is many-to-one, and the conservative rule applies to it too."""
        pieces = pieces_or_skip(backend)
        step = pieces.bandwidth_smearing(resolving_power=RESOLVING_POWER, nodes=3)
        u_pts = np.repeat(np.array([1.0e7, 2.0e7]), 3)
        mask = np.zeros(6, dtype=bool)
        mask[4] = True
        source = VisibilitySet(
            u_pts,
            0.5 * u_pts,
            np.full(6, WAVELENGTH) * u.micron,
            np.ones(6, dtype=np.complex128) * u.Jy,
            mask=mask,
        )
        got = step(source)
        assert got.n_samples == 2
        assert got.mask is not None
        assert got.mask.tolist() == [False, True]


# ---------------------------------------------------------------------------
# The composition
# ---------------------------------------------------------------------------


class TestTwoDatasetsOnOneChannel:
    """Visibilities and closure phases from one sky, negotiated once."""

    def _problem(self, backend: ConformanceBackend, pieces: InterferometryPieces) -> Any:
        truth = image_model(pieces, "binary")
        vis_instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed_visibilities(),
                    field_of_view=FIELD_OF_VIEW * u.mas,
                    oversampling=OVERSAMPLING,
                )
            ],
            channel="sky",
            label="array_vis",
        )
        t3_instrument = Instrument(
            [
                pieces.fourier_sample.from_observed(
                    observed_closure_phases(),
                    field_of_view=FIELD_OF_VIEW * u.mas,
                    oversampling=OVERSAMPLING,
                ),
                pieces.closure_phase(),
            ],
            channel="sky",
            label="array_t3",
        )
        compiled = truth.compile_for(negotiate([vis_instrument, t3_instrument]))
        result = compiled.evaluate()
        visibilities = observed_visibilities(vis_instrument(result).values)
        phases = observed_closure_phases(t3_instrument(result).values)
        fitted = {
            **BINARY,
            "separation": st.uniform(4.0, 20.0),
            "flux_ratio": st.uniform(0.05, 0.9),
        }
        model = counting(pieces.binary)(
            compiled.buffers["x"].value * u.mas,
            compiled.buffers["y"].value * u.mas,
            channels="sky",
            **fitted,
        )
        datasets = DatasetCollection(
            {
                "vis": Dataset(
                    visibilities,
                    vis_instrument,
                    likelihood=Likelihood(ComplexGaussianFamily(), backend.independent_noise()),
                    label="vis",
                ),
                "t3": Dataset(
                    phases,
                    t3_instrument,
                    likelihood=Likelihood(VonMisesFamily(), backend.independent_noise()),
                    label="t3",
                ),
            }
        )
        return FittingProblem(model, datasets, seed=20260911), model

    def test_the_channel_records_both_instruments(self, backend: ConformanceBackend) -> None:
        """Design horizon (h)'s obligation: the pairing stays visible.

        Both datasets name the ``sky`` channel and the negotiation record
        names both instruments, so a Phase 5 noise model over a tuple of
        channels has something to bind to. Nothing here folds the two
        observables into one container.
        """
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        channel = problem.requirements["model"]["sky"]
        assert set(channel.sources) == {"array_vis", "array_t3"}
        assert sorted(channel.axes) == ["x", "y"]
        assert {type(problem.datasets[name].observed).__name__ for name in ("vis", "t3")} == {
            "VisibilitySet",
            "ClosurePhases",
        }

    def test_the_model_is_evaluated_once_per_draw(self, backend: ConformanceBackend) -> None:
        """``inference.md`` §8: once, jointly — not once per dataset."""
        pieces = pieces_or_skip(backend)
        problem, model = self._problem(backend, pieces)
        model.evaluations = 0
        problem.log_prob(TRUTH)
        assert model.evaluations == 1

    def test_the_composed_likelihood_peaks_at_the_truth(self, backend: ConformanceBackend) -> None:
        """The two likelihoods together prefer the source the data were made from."""
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        at_truth = problem.log_prob(TRUTH)
        displaced = problem.log_prob({**TRUTH, "model.separation": 16.0})
        assert np.isfinite(at_truth)
        assert at_truth > displaced

    def test_simulate_draws_observations_of_both_kinds(self, backend: ConformanceBackend) -> None:
        """``simulate(observe=True)`` on a complex kind and on a wrapped one.

        ``complex_gaussian`` has drawn since W3.14; ``von_mises`` gains its
        ``sample()`` here, with the data type that needed it, which is what
        the sampling-form principle (``likelihoods.md`` §3) asks for.
        """
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        simulation = problem.simulate(TRUTH, observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        drawn_vis = simulation.observations["vis"]
        drawn_t3 = simulation.observations["t3"]
        assert drawn_vis.values.dtype.kind == "c"
        assert drawn_t3.values.dtype.kind == "f"
        assert np.all(np.abs(drawn_t3.values) <= math.pi)
        # A draw is not the prediction: the noise actually moved the data.
        assert not np.array_equal(drawn_vis.values, simulation.predicted["vis"].values)
        assert not np.array_equal(drawn_t3.values, simulation.predicted["t3"].values)


# ---------------------------------------------------------------------------
# The amplitude route
# ---------------------------------------------------------------------------


class TestAmplitudeRoute:
    """The third route of ``interferometry.md`` §6: fit ``|V|``, not ``V``."""

    def test_the_modulus_step_produces_a_real_visibility_set(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """A real-valued ``VisibilitySet`` is legal, and this is how one is made."""
        pieces = pieces_or_skip(backend)
        u_pts, v_pts, _, _ = visibility_coverage()
        got = transformed(pieces, image_model(pieces, "binary"), steps=(pieces.amplitude(),))
        assert got.values.dtype.kind == "f"
        expected = np.abs(binary_visibility(u_pts, v_pts, **BINARY))
        error = np.max(np.abs(backend.to_numpy(got.values) - expected)) / np.max(expected)
        assert error < tolerances.cross_solver

        amplitudes = VisibilitySet(
            u_pts,
            v_pts,
            np.full(u_pts.size, WAVELENGTH) * u.micron,
            expected * u.Jy,
            uncertainty=np.full(u_pts.size, 0.02) * u.Jy,
        )
        # The step is what makes the two sides comparable: a complex
        # prediction against real amplitudes is refused (gap I-1).
        Likelihood(GaussianFamily(), backend.independent_noise()).check_alignment(got, amplitudes)

    def test_the_rice_family_is_still_what_is_owed(self, backend: ConformanceBackend) -> None:
        """The route composes; the family it wants does not exist yet, and says so.

        ``likelihoods.md`` §3's sampling-form principle is why: a family's
        likelihood and its draw arrive together with the data type that needs
        them, and Rice's data type is polarimetry, which is a later item. The
        route is the part W4.1 owes, and the refusal is by name rather than a
        silent fallback to a Gaussian on an amplitude.
        """
        pieces_or_skip(backend)
        with pytest.raises(LikelihoodError, match="declared but not implemented"):
            Likelihood(RiceFamily(), backend.independent_noise())
