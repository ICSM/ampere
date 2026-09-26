"""Gridded images: the rows the template's third modality owes (W5.5).

One column per registered backend, like every other file here, and no test body
names one. A backend that has not written the gridded vocabulary declares
``BackendCapabilities.image = False`` and every row skips with a reason naming
what is owed — the same shape ``test_interferometry.py`` and
``test_astrometry.py`` use for their own.

``docs/source/interferometry.rst`` §7 states the four rows a new modality owes.
This modality's version of them:

* **a row per step against a closed form** — a PSF convolution has no closed
  form to be checked against in general, so the oracle is the *definition*: a
  direct sum over the kernel's support, written here as three nested loops with
  nothing in common with an FFT except the answer. Both PSF forms, the
  tabulated kernel and the analytic Gaussian, at ``tolerances.cross_solver``;
* **a requirement row, and its refusal** — the padded grid the step publishes,
  the negotiated union it produces, and a model holding its own grid refusing
  by name (``CompositionError``) when that grid is coarser than the pixels
  actually observed;
* **a mask-propagation row** — one masked model pixel masks every observed
  pixel it contaminates, which on a grid is a dilation by the kernel's support
  rather than a matrix (``propagate_mask_grid``, this item's addition to
  ``transformations.md`` §13.5);
* **a draw row** — ``FittingProblem.simulate`` producing an ``Image``
  observation, the first time the ``Layout.GRID`` branch of
  ``Dataset.draw_observation`` has ever run.

Plus the cross-backend pairing the item asks for explicitly: the same
declaration handed to two backends, at ``tolerances.cross_backend``.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    CompositionError,
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Image,
    Instrument,
    Likelihood,
    LikelihoodError,
    TransformationError,
    negotiate,
)

from .conftest import cross_backend_pairs, pair_ids
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    ImagePieces,
    KernelFamily,
    SolverKind,
    Tolerances,
)

#: Pixels on a side of the observed image. Small on purpose — every row here is
#: about the transform's *rules*, and a 12x12 image exercises every one of them
#: while a direct triple loop over the kernel stays instant.
PIXELS = 12

#: The observed field, mas. A 12-pixel image across 11 mas is a pixel scale of
#: 1 mas, which is the same order as the PSF below: an image whose PSF is much
#: finer than a pixel would make the convolution a no-op and the rows vacuous.
HALF_FIELD = 5.5

#: The Gaussian PSF's full width at half maximum, mas.
PSF_FWHM = 2.0

#: The source: a Gaussian rather than a uniform disc, because a disc's sharp
#: edge is not band-limited and W4.1's own rows record that it converges as a
#: power of the pixel scale rather than to a solver tolerance. That is a fact
#: about discs, and it would confuse a row that is about convolution.
SOURCE = {"flux": 1.5, "fwhm": 3.0}

#: The surface-brightness unit W4.1's image models emit in, and therefore the
#: unit an observed image of one is measured in. ``check_alignment`` compares
#: units as well as axes, so this is not cosmetic.
BRIGHTNESS = u.Jy / u.sr

#: A two-axis Matérn-3/2 over the image's own ``(x, y)``, mas — W5.21's
#: composition rows. The length scale is of the order of the PSF, which is the
#: scale a flexible arm is meant to absorb; the value is fixed rather than
#: fitted because these rows ask only whether the composition is *accepted*.
GRID_KERNEL = CovarianceSpec(
    family=KernelFamily.MATERN32,
    amplitude=0.05,
    length_scale=3.0,
    axes=("x", "y"),
    length_scale_unit=u.mas,
)


def observed_grid() -> np.ndarray:
    """The observed pixel centres, mas. Evenly spaced, as the step requires."""
    return np.linspace(-HALF_FIELD, HALF_FIELD, PIXELS)


def observed_image(values: np.ndarray | None = None, mask: np.ndarray | None = None) -> Image:
    """An ``Image`` on :func:`observed_grid`, with a uniform uncertainty."""
    grid = observed_grid()
    filled = np.zeros((PIXELS, PIXELS)) if values is None else np.asarray(values)
    scale = float(np.max(np.abs(filled))) or 1.0
    return Image(
        grid * u.mas,
        grid * u.mas,
        filled * BRIGHTNESS,
        uncertainty=np.full((PIXELS, PIXELS), 0.02 * scale) * BRIGHTNESS,
        mask=mask,
    )


def pieces_or_skip(backend: ConformanceBackend) -> ImagePieces:
    """This backend's gridded-image classes, or a skip naming what is owed."""
    if not backend.capabilities.image:
        pytest.skip(
            f"{backend.name} declares no gridded-image vocabulary "
            f"(BackendCapabilities.image), so the PSF-convolution step these rows drive is not "
            f"there. W5.5 owes the native twins."
        )
    return backend.image()


def cross_kernel() -> np.ndarray:
    """A deliberately non-Gaussian 5x5 tabulated PSF: a plus with a bright core.

    Not separable and not symmetric under a 45-degree rotation, so a
    convolution that silently transposed its axes or applied a separable
    shortcut would show up here rather than agreeing by accident.
    """
    kernel = np.zeros((5, 5))
    kernel[2, 2] = 4.0
    kernel[1, 2] = kernel[3, 2] = 2.0
    kernel[2, 1] = kernel[2, 3] = 1.0
    kernel[0, 2] = 0.5
    return kernel


def step_for(pieces: ImagePieces, observed: Image, *, tabulated: bool = False) -> Any:
    """The PSF-convolution step, built from the observed image the supported way."""
    if tabulated:
        return pieces.psf_convolution.from_observed(observed, kernel=cross_kernel(), label="psf")
    return pieces.psf_convolution.from_observed(observed, fwhm=PSF_FWHM, label="psf")


def camera(step: Any) -> Instrument:
    """The one-step imaging instrument."""
    return Instrument([step], channel="sky", label="camera")


def source_model(pieces: ImagePieces, **overrides: Any) -> Any:
    """The Gaussian sky source on its own default grid, at the truth."""
    grid = observed_grid()
    return pieces.gaussian_source(
        grid * u.mas, grid * u.mas, channels="sky", **{**SOURCE, **overrides}
    )


def predicted_through(backend: ConformanceBackend, pieces: ImagePieces, **kwargs: Any) -> Any:
    """``(step, model image, convolved image)`` for the standard configuration."""
    observed = observed_image()
    step = step_for(pieces, observed, **kwargs)
    instrument = camera(step)
    compiled = source_model(pieces).compile_for(negotiate([instrument]))
    result = compiled.evaluate()
    return step, result["sky"], instrument(result)


def direct_convolution(image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """The definition, as a sum. **The oracle, and it shares nothing with the step.**

    ``out[i, j] = sum_ab kernel[a, b] image[i - a + ka, j - b + kb]``, with
    anything off the edge contributing zero — a zero-padded linear convolution,
    which is what the FFT route computes and the only thing it is allowed to
    compute. Written with explicit index arithmetic rather than through any
    library's ``convolve`` so that the comparison is against the definition
    rather than against a second implementation with the same conventions.
    """
    n_x, n_y = image.shape
    k_x, k_y = kernel.shape
    out = np.zeros((n_x, n_y))
    for i in range(n_x):
        for j in range(n_y):
            total = 0.0
            for a in range(k_x):
                for b in range(k_y):
                    row = i - a + k_x // 2
                    column = j - b + k_y // 2
                    if 0 <= row < n_x and 0 <= column < n_y:
                        total += kernel[a, b] * image[row, column]
            out[i, j] = total
    return out


# ---------------------------------------------------------------------------
# The step against the definition
# ---------------------------------------------------------------------------


class TestConvolutionAgainstTheDefinition:
    """The FFT route against a direct sum over the kernel's support."""

    @pytest.mark.parametrize("tabulated", [False, True], ids=["gaussian", "tabulated"])
    def test_the_fft_route_matches_a_direct_sum(
        self, backend: ConformanceBackend, tolerances: Tolerances, tabulated: bool
    ) -> None:
        pieces = pieces_or_skip(backend)
        step, model_image, convolved = predicted_through(backend, pieces, tabulated=tabulated)
        x_grid = np.asarray(model_image.x.values)
        y_grid = np.asarray(model_image.y.values)
        kernel = np.asarray(backend.to_numpy(step.psf(x_grid, y_grid)))
        expected = direct_convolution(np.asarray(backend.to_numpy(model_image.values)), kernel)
        rows, columns = step.crop_indices(x_grid, y_grid)
        got = np.asarray(backend.to_numpy(convolved.values))
        reference = expected[np.ix_(rows, columns)]
        scale = float(np.max(np.abs(reference)))
        assert scale > 0.0
        assert np.allclose(got, reference, rtol=0.0, atol=tolerances.cross_solver * scale)

    def test_the_step_conserves_flux_over_the_padded_field(
        self, backend: ConformanceBackend
    ) -> None:
        """A normalised kernel moves flux about; it does not create or destroy it.

        Checked on the *padded* grid rather than on the observed crop, because
        the crop legitimately throws away the flux the convolution pushed into
        the padding — which is exactly what the padding is there to absorb.
        """
        pieces = pieces_or_skip(backend)
        step, model_image, _ = predicted_through(backend, pieces)
        x_grid = np.asarray(model_image.x.values)
        y_grid = np.asarray(model_image.y.values)
        kernel = np.asarray(backend.to_numpy(step.psf(x_grid, y_grid)))
        assert kernel.sum() == pytest.approx(1.0, rel=1e-12)
        image = np.asarray(backend.to_numpy(model_image.values))
        blurred = direct_convolution(image, kernel)
        assert blurred.sum() == pytest.approx(image.sum(), rel=1e-6)

    def test_the_emitted_image_carries_the_observed_axes_themselves(
        self, backend: ConformanceBackend
    ) -> None:
        """``transformations.md`` §10's coordinates-from-the-container rule.

        ``check_alignment`` compares axes with ``np.array_equal``, so "equal to
        floating-point tolerance" is not enough: the emitted container must
        carry the observed coordinates, not a recomputation of them.
        """
        pieces = pieces_or_skip(backend)
        observed = observed_image()
        _, _, convolved = predicted_through(backend, pieces)
        assert isinstance(convolved, Image)
        assert np.array_equal(convolved.x.values, observed.x.values)
        assert np.array_equal(convolved.y.values, observed.y.values)
        assert convolved.values.shape == (PIXELS, PIXELS)


# ---------------------------------------------------------------------------
# What the step publishes, and what it refuses
# ---------------------------------------------------------------------------


class TestTheRequirement:
    """The padded grid, the union it produces, and the refusal it earns."""

    def test_the_requirement_pads_the_observed_field_by_the_kernel_support(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        step = step_for(pieces, observed_image())
        published = {requirement.axis: requirement for requirement in step.requirements()}
        assert set(published) == {"x", "y"}
        grid = observed_grid()
        scale = float(grid[1] - grid[0])
        pad_x, pad_y = step.support()
        for axis, pad in (("x", pad_x), ("y", pad_y)):
            requirement = published[axis]
            assert requirement.points is not None
            assert requirement.points.size == PIXELS + 2 * pad
            (low, high, max_step, power) = requirement.segments()[0]
            assert power is None
            assert max_step == pytest.approx(scale, rel=1e-3)
            assert low == pytest.approx(grid[0] - pad * scale, rel=1e-9)
            assert high == pytest.approx(grid[-1] + pad * scale, rel=1e-9)

    def test_the_negotiated_grid_is_evenly_spaced_and_holds_the_observed_pixels(
        self, backend: ConformanceBackend
    ) -> None:
        """The two published halves must agree, or the FFT route has no grid.

        ``points=`` and ``intervals``/``max_step`` state the same constraint two
        ways, and ``coordinates()`` unions them. If they disagreed by more than
        the union's own collapsing tolerance the result would be a grid with
        near-duplicate coordinates — evenly spaced nowhere, and refused by the
        step rather than by negotiation, which is a bad place to find out.
        """
        pieces = pieces_or_skip(backend)
        observed = observed_image()
        step = step_for(pieces, observed)
        needed = negotiate([camera(step)])
        for axis in ("x", "y"):
            coordinates = np.asarray(needed["sky"][axis].coordinates().to_value(u.mas))
            spacing = np.diff(coordinates)
            assert np.allclose(spacing, spacing[0], rtol=1e-9)
        # And the observed centres survive the union, which is what the crop
        # looks up.
        rows, columns = step.crop_indices(
            np.asarray(needed["sky"]["x"].coordinates().to_value(u.mas)),
            np.asarray(needed["sky"]["y"].coordinates().to_value(u.mas)),
        )
        assert rows.size == PIXELS
        assert columns.size == PIXELS

    def test_a_model_holding_a_coarser_grid_is_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        """Gap I-4's mechanism, on an image being observed rather than transformed.

        A model that holds its own grid (``adopt_grid=False``) and samples it
        more coarsely than the data were sampled cannot represent what was
        observed. ``compile_for`` says so, with ``CompositionError``, rather
        than letting the failure surface later as a missing coordinate.
        """
        pieces = pieces_or_skip(backend)
        step = step_for(pieces, observed_image())
        needed = negotiate([camera(step)])
        coarse = np.linspace(-4.0 * HALF_FIELD, 4.0 * HALF_FIELD, 5)
        fixed = pieces.gaussian_source(
            coarse * u.mas, coarse * u.mas, channels="sky", adopt_grid=False, **SOURCE
        )
        with pytest.raises(CompositionError, match="step"):
            fixed.compile_for(needed)

    def test_an_irregular_observed_grid_is_refused_at_construction(
        self, backend: ConformanceBackend
    ) -> None:
        """The FFT route's own limitation, stated as a refusal rather than a fallback."""
        pieces = pieces_or_skip(backend)
        grid = observed_grid()
        bent = grid.copy()
        bent[PIXELS // 2] += 0.3 * float(grid[1] - grid[0])
        irregular = Image(bent * u.mas, grid * u.mas, np.zeros((PIXELS, PIXELS)) * BRIGHTNESS)
        with pytest.raises(TransformationError, match="evenly spaced"):
            step_for(pieces, irregular)

    def test_a_tabulated_kernel_refuses_a_different_pixel_scale(
        self, backend: ConformanceBackend
    ) -> None:
        """A measured PSF is a buffer tied to the coordinates it was tabulated on."""
        pieces = pieces_or_skip(backend)
        step = step_for(pieces, observed_image(), tabulated=True)
        finer = np.linspace(-HALF_FIELD, HALF_FIELD, 3 * PIXELS)
        with pytest.raises(TransformationError, match="tabulated"):
            step.psf(finer, finer)


# ---------------------------------------------------------------------------
# The mask
# ---------------------------------------------------------------------------


class TestMaskPropagationOnAGrid:
    """One masked model pixel masks every observed pixel it contaminates."""

    def test_a_masked_model_pixel_masks_its_whole_neighbourhood(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        observed = observed_image()
        step = step_for(pieces, observed)
        instrument = camera(step)
        compiled = source_model(pieces).compile_for(negotiate([instrument]))
        model_image = compiled.evaluate()["sky"]
        pad_x, pad_y = step.support()
        shape = model_image.values.shape
        mask = np.zeros(shape, dtype=bool)
        centre = (shape[0] // 2, shape[1] // 2)
        mask[centre] = True
        blurred = step(model_image.with_values(model_image.values, mask=mask), None)
        assert blurred.mask is not None
        # The masked pixel sits well inside the padded grid, so its whole
        # (2k+1)x(2k+1) neighbourhood is on the grid; the crop then keeps the
        # part of it that overlaps the observed field.
        rows, columns = step.crop_indices(
            np.asarray(model_image.x.values), np.asarray(model_image.y.values)
        )
        expected = np.zeros(shape, dtype=bool)
        expected[
            centre[0] - pad_x : centre[0] + pad_x + 1,
            centre[1] - pad_y : centre[1] + pad_y + 1,
        ] = True
        assert np.array_equal(blurred.mask, expected[np.ix_(rows, columns)])
        assert bool(blurred.mask.any())

    def test_an_unmasked_image_stays_unmasked(self, backend: ConformanceBackend) -> None:
        """The usual case pays nothing: ``None`` in, ``None`` out."""
        pieces = pieces_or_skip(backend)
        _, _, convolved = predicted_through(backend, pieces)
        assert convolved.mask is None


# ---------------------------------------------------------------------------
# The composition, and a draw
# ---------------------------------------------------------------------------


def image_problem(backend: ConformanceBackend, pieces: ImagePieces) -> FittingProblem:
    """A one-dataset fitting problem whose observation is an ``Image``."""
    truth = predicted_through(backend, pieces)[2]
    observed = observed_image(np.asarray(backend.to_numpy(truth.values)))
    step = step_for(pieces, observed)
    instrument = camera(step)
    model = pieces.gaussian_source(
        observed_grid() * u.mas,
        observed_grid() * u.mas,
        channels="sky",
        flux=st.lognorm(0.4, scale=SOURCE["flux"]),
        fwhm=SOURCE["fwhm"],
    )
    datasets = DatasetCollection(
        {
            "image": Dataset(
                observed,
                instrument,
                likelihood=Likelihood(GaussianFamily(), backend.independent_noise()),
                label="image",
            )
        }
    )
    return FittingProblem(model, datasets, seed=20260916)


class TestAnImageAsAnObservation:
    """The composition the item exists for: ``Image`` on both sides of a likelihood."""

    def test_the_channel_is_recorded_and_the_likelihood_peaks_at_the_truth(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        problem = image_problem(backend, pieces)
        assert set(problem.requirements["model"]) == {"sky"}
        assert type(problem.datasets["image"].observed).__name__ == "Image"
        at_truth = problem.log_prob({"model.flux": SOURCE["flux"]})
        displaced = problem.log_prob({"model.flux": 2.5 * SOURCE["flux"]})
        assert np.isfinite(at_truth)
        assert at_truth > displaced

    def test_simulate_draws_an_image_observation(self, backend: ConformanceBackend) -> None:
        """``Dataset.draw_observation``'s ``Layout.GRID`` branch, exercised."""
        pieces = pieces_or_skip(backend)
        problem = image_problem(backend, pieces)
        simulation = problem.simulate({"model.flux": SOURCE["flux"]}, observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        drawn = simulation.observations["image"]
        assert isinstance(drawn, Image)
        assert drawn.values.shape == (PIXELS, PIXELS)
        assert drawn.values.dtype.kind == "f"
        # A draw is not the prediction: the noise actually moved the data.
        assert not np.array_equal(drawn.values, simulation.predicted["image"].values)

    @pytest.mark.parametrize(
        "kind", [SolverKind.DENSE, SolverKind.HILBERT], ids=["dense", "hilbert"]
    )
    def test_a_correlated_noise_model_is_accepted_on_a_grid_since_w5_21(
        self, backend: ConformanceBackend, kind: SolverKind
    ) -> None:
        """W5.5's blanket ``Layout.GRID`` refusal is lifted (W5.21): the shipped
        exact and reduced-rank solvers compose with a two-axis kernel over an
        ``Image`` exactly as they do over a point-set container."""
        if kind not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no {kind.value} solver")
        solver = (
            backend.gp_solver(kind, basis_size=(6, 6), boundary_factor=2.0)
            if kind is SolverKind.HILBERT
            else backend.gp_solver(kind)
        )
        noise = backend.gp_noise(backend.kernel(GRID_KERNEL), solver)
        observed = observed_image()
        likelihood = Likelihood(GaussianFamily(), noise)
        likelihood.check_alignment(observed, observed)  # no longer raises
        assert np.isfinite(likelihood.log_prob(observed, observed))

    def test_a_quasiseparable_solver_is_still_refused_by_name_on_a_grid(
        self, backend: ConformanceBackend
    ) -> None:
        """Lifting the layout gate does not touch ``REQUIRES_ORDERED_1D``: a
        kernel selecting both of an ``Image``'s axes is still refused by a
        solver that needs exactly one ordered coordinate."""
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"backend {backend.name!r} declares no quasiseparable solver")
        kernel = backend.kernel(GRID_KERNEL)
        noise = backend.gp_noise(kernel, backend.gp_solver(SolverKind.QUASISEP))
        observed = observed_image()
        with pytest.raises(LikelihoodError, match="needs one ordered coordinate axis"):
            Likelihood(GaussianFamily(), noise).check_alignment(observed, observed)


# ---------------------------------------------------------------------------
# The twins against each other
# ---------------------------------------------------------------------------


class TestCrossBackendAgreement:
    """The same declaration handed to two backends, at ``cross_backend``."""

    @pytest.mark.parametrize(("first", "second"), cross_backend_pairs(), ids=pair_ids())
    @pytest.mark.parametrize("tabulated", [False, True], ids=["gaussian", "tabulated"])
    def test_two_backends_convolve_the_same_image_the_same_way(
        self, first: ConformanceBackend, second: ConformanceBackend, tabulated: bool
    ) -> None:
        if not (first.capabilities.image and second.capabilities.image):
            pytest.skip("both backends must declare BackendCapabilities.image")
        outputs = []
        for backend in (first, second):
            _, _, convolved = predicted_through(backend, backend.image(), tabulated=tabulated)
            outputs.append(np.asarray(backend.to_numpy(convolved.values)))
        tolerance = max(
            first.capabilities.tolerances.cross_backend,
            second.capabilities.tolerances.cross_backend,
        )
        scale = float(np.max(np.abs(outputs[0])))
        assert scale > 0.0
        assert np.allclose(outputs[0], outputs[1], rtol=0.0, atol=tolerance * scale)
