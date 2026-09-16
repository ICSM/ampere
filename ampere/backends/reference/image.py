"""Gridded images on the numpy path: the PSF-convolution step.

Phase 5's gridded customer (W5.5). ``transformations.md`` §10's second table
has held an *image* slot since the freeze — "PSF convolution, ``Image`` →
``Image``; the image analogue of LSF convolution; the kernel is usually a
buffer" — and nothing had filled it: :class:`~ampere.core.Image` existed as a
kind, Phase 4's three source models emitted one, and every one of them handed
it straight to :class:`~ampere.backends.reference.interferometry.FourierSample`.
No ``Image`` had ever been a *dataset*. This module is that step, and with it
the first fit whose observed container has :attr:`~ampere.core.Layout.GRID`.

Placement follows D1's ruling as ``transformations.md`` §10 records it — "one
name per standard transformation, one module per observable, no grouping
namespace". The spectral rows live in ``instrument.py``, the Fourier row in
``interferometry.py``, epoch sampling in ``astrometry.py``; a step whose
observable is a *gridded image* accordingly lives in ``image.py``, one per
backend. (The three image-emitting source models stay where W4.1 put them, in
``interferometry.py``: they are the same objects, they are reused here
unchanged, and moving them would be a rename for its own sake.)

What is genuinely new here, beyond one more step
------------------------------------------------
**A grid is not a point set, and the mask rule had only ever been written for
a point set.** ``transformations.md`` §13.5 was explicit about it: an
influence matrix is ``(n_out, n_in)``, which fits ``Layout.POINTS``, and "a
``Layout.GRID`` container's mask must be flattened by the transformation
itself". Flattening is the wrong answer for a convolution — a 64x64 image
against a 9x9 kernel would build a sixteen-million-entry boolean matrix of
which 81 entries per row are non-zero — so W5.5 adds
:func:`~ampere.core.propagate_mask_grid` beside :func:`~ampere.core.propagate_mask`
and this step is its first caller. The *rule* is unchanged (an output pixel
touching any masked input pixel is masked); only its expression is, from a
matrix to a separable dilation by the kernel's support.

**A convolution needs the model evaluated beyond the data.** The same problem
:class:`~ampere.backends.reference.instrument.LSFConvolution` has, in two
dimensions, and with a different answer to it. An LSF publishes a padded
*interval* and lets the resampler downstream land on the observed grid; a PSF
has no resampler downstream — there is no spatial one — so it publishes a
padded grid **containing the observed pixel centres exactly**, convolves on
it, and crops. The padding is what stops the outermost observed pixels being
convolved against invented zeros, and the crop is ``transformations.md`` §10's
coordinates-from-the-container rule: the emitted ``Image`` carries the
observed axes themselves, never a recomputed copy of them.

The FFT route, and what it costs
--------------------------------
:meth:`PSFConvolution.apply` computes a zero-padded *linear* convolution with
:mod:`numpy.fft` and takes the central block, which is the same number a
direct sum over the kernel's support computes and is what
``tests/conformance/test_image.py`` holds it to. An FFT needs an evenly spaced
grid, so this step **requires one** — of its target and of whatever
negotiation hands it — rather than advertising a fast path and falling back
(:attr:`~ampere.core.Axis.regular` is how a step would advertise instead). The
refusal is loud and names the extension point: a spatial resampler, or a
direct-matrix route, is what an irregular sky grid would want, and neither is
this item's.
"""

from __future__ import annotations

import math
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
from scipy.fft import next_fast_len

from ampere.core import (
    COORDINATE_RTOL,
    DTYPE,
    Axis,
    AxisRequirement,
    Image,
    Transformation,
    TransformationError,
    propagate_mask_grid,
)

__all__ = ["COORDINATE_UNIT", "PSFConvolution"]

#: The sky-coordinate unit every image step here speaks, as
#: ``interferometry.py``'s image models do.
COORDINATE_UNIT = u.mas

#: How far beyond its centre a Gaussian PSF is tabulated, in kernel sigmas.
#: Truncating at 5 sigma loses about 1e-6 of a two-dimensional Gaussian's
#: volume, which is far below anything an image fit can see, and the truncated
#: kernel is renormalised in any case so that the step conserves flux exactly.
_PSF_PADDING_SIGMAS = 5.0

#: Relative slack on the published ``max_step``. The requirement states the
#: observed pixel scale, and the padded grid this step publishes as ``points=``
#: already samples at exactly that; the slack exists so that
#: :meth:`~ampere.core.AxisRequirement.coordinates`' own reference grid over
#: the same interval comes out with the same number of points rather than one
#: more from a last-bit ``ceil``. Anything between float noise and a part per
#: million would do; the number itself is not physics.
_STEP_SLACK = 1e-6

#: How far apart two statements of the same pixel scale may be before a
#: tabulated kernel is refused. Loose beside :data:`~ampere.core.COORDINATE_RTOL`
#: because it compares two *derived* spacings rather than two coordinates.
_SCALE_RTOL = 1e-9


def _to_mas(value: Any) -> np.ndarray:
    """*value* in mas, as a plain float array. A bare number is already mas."""
    if isinstance(value, u.Quantity):
        return np.asarray(value.to_value(COORDINATE_UNIT), dtype=DTYPE)
    return np.asarray(value, dtype=DTYPE)


def _uniform_step(grid: np.ndarray, what: str) -> float:
    """The common spacing of *grid*, or a refusal naming *what*.

    Tighter than :attr:`~ampere.core.Axis.regular`'s ``REGULARITY_RTOL``, and
    deliberately: this step matches the observed pixel centres inside a
    negotiated grid with :meth:`~ampere.core.Axis.locate`, whose tolerance is
    :data:`~ampere.core.COORDINATE_RTOL`, so a grid "regular" only to a looser
    tolerance than that is one whose pixel centres a published, evenly spaced
    requirement could not reproduce well enough to be found again. A grid built
    the way sky grids are built — ``numpy.linspace``, or a linear WCS — passes
    this comfortably.
    """
    if grid.ndim != 1 or grid.size < 2:
        raise TransformationError(
            f"PSFConvolution needs {what} as at least two coordinates along one axis, in mas; "
            f"got an array of shape {grid.shape}."
        )
    if not np.all(np.diff(grid) > 0.0):
        raise TransformationError(
            f"PSFConvolution needs {what} strictly increasing. An Image axis may legally run "
            f"either way (sky axes do), but a negotiated grid is always built ascending, so a "
            f"descending observed axis could never be matched against one. Flip the image and "
            f"its axis together before building the dataset."
        )
    even = np.linspace(grid[0], grid[-1], grid.size)
    if not np.allclose(grid, even, rtol=COORDINATE_RTOL, atol=0.0):
        raise TransformationError(
            f"PSFConvolution convolves by FFT, which needs {what} evenly spaced; the coordinates "
            f"it was given are not, to within COORDINATE_RTOL ({COORDINATE_RTOL:g}). Build the "
            f"axis with numpy.linspace (or from a linear WCS). An irregular sky grid wants a "
            f"spatial resampling step — transformations.md §10's other image slot — or a "
            f"direct-matrix convolution; neither is this step."
        )
    return float((grid[-1] - grid[0]) / (grid.size - 1))


def _fft_convolve(image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
    """Zero-padded linear convolution, central block, by FFT.

    Padded to ``n + m - 1`` on each axis (rounded up to a transform-friendly
    length) so that the periodic wrap an unpadded FFT would produce cannot fold
    the far edge of the image onto the near one. The block taken back is the
    one whose pixels line up with *image*'s own — offset by the kernel's
    half-width, which is why an odd-sized kernel with a defined centre pixel is
    required at construction.

    The zeros the padding introduces are not a boundary approximation this step
    lives with: :meth:`PSFConvolution.requirements` asks for an image that
    extends a full half-support beyond the observed field on every side, and
    :meth:`PSFConvolution.apply` crops back to the observed pixels, so no
    returned pixel ever sees one of them.
    """
    shape = tuple(n + m - 1 for n, m in zip(image.shape, kernel.shape, strict=True))
    fast = tuple(next_fast_len(length) for length in shape)
    spectrum = np.fft.rfftn(image, fast) * np.fft.rfftn(kernel, fast)
    full = np.fft.irfftn(spectrum, fast)
    centre = tuple(
        slice(m // 2, m // 2 + n) for n, m in zip(image.shape, kernel.shape, strict=True)
    )
    return np.asarray(full[centre], dtype=DTYPE)


class PSFConvolution(Transformation):
    """Convolve a model image with a point-spread function: ``Image -> Image``.

    The image analogue of
    :class:`~ampere.backends.reference.instrument.LSFConvolution`, and
    ``transformations.md`` §10's image slot, filled at W5.5. Kind-preserving
    (``PRODUCES = None``) and coordinate-*changing* in the way
    :class:`~ampere.backends.reference.instrument.Resample` is: the emitted
    image is on the **observed** grid, which is a crop of the padded grid the
    step asks negotiation for.

    Two ways to say what the PSF is, and the difference between them matters
    ------------------------------------------------------------------------
    ``kernel=`` — **a pixel kernel the caller supplies**, the form a real
    instrument uses: an odd-shaped array tabulated at the observed pixel scale,
    from a calibration observation of a point source or from an optical model.
    It is a buffer, which is ``parameters.md`` §10's answer for constant data,
    and it is ``transformations.md`` §10's *second* rule in its gridded form —
    a buffer tabulated on particular coordinates. A tabulated kernel cannot be
    rescaled without inventing information, so this step refuses a negotiated
    grid whose pixel scale is not the one the kernel was tabulated at, rather
    than interpolating quietly.

    ``fwhm=`` — **a circular Gaussian of a given full width at half maximum**,
    in mas. Also a buffer, for the same reason an LSF's width is one; and as
    with an LSF, ``promote_buffer("fwhm", prior=...)`` turns it into an
    ordinary fitted parameter with no change to this class, because
    :meth:`apply` rebuilds the kernel through
    :meth:`~ampere.core.Parameterised.context` on every evaluation. An analytic
    kernel *can* be retabulated, so this form accepts any evenly spaced
    negotiated grid at or finer than the observed pixel scale.

    A fitted width has one caveat worth stating, because it is the kind of
    thing that goes unnoticed: the **support** — how far the kernel is
    tabulated, and therefore how much padding :meth:`requirements` asks the
    model for — is fixed once, at composition time, from the declared value. A
    posterior that wanders far above it is fitting a kernel truncated tighter
    than five sigma. The truncated kernel is renormalised, so flux is still
    conserved exactly; what is lost is the faintest wings. Declare the width
    you expect, and give the padding room.

    What it publishes
    -----------------
    Per axis, one :class:`~ampere.core.AxisRequirement` carrying both halves of
    the constraint: ``points=`` at an evenly spaced grid extending the observed
    pixel centres by the kernel's half-support on each side — the coverage a
    convolution needs, *and* the exact coordinates :meth:`apply` will look up —
    together with ``intervals``/``max_step`` restating the observed pixel scale
    so that a model holding a grid of its own is refused by name when that grid
    is too coarse to represent what was observed, rather than aliasing into it.
    That refusal is the image analogue of interferometry's gap I-4, and it is
    the same mechanism: ``_ImageModel(adopt_grid=False)`` checks the published
    requirement against its own grid and raises
    :class:`~ampere.core.CompositionError`.

    Parameters
    ----------
    x, y
        The observed image's own coordinates, mas (or a
        :class:`~astropy.units.Quantity`). Evenly spaced and strictly
        increasing. Prefer :meth:`from_observed`, which takes them from the
        container.
    kernel
        The PSF as a pixel array with odd extent on both axes, tabulated at the
        observed pixel scale. Mutually exclusive with *fwhm*. Normalised to
        unit sum at construction, so the step conserves total flux.
    fwhm
        A circular Gaussian PSF's full width at half maximum, mas. Mutually
        exclusive with *kernel*.
    label
        Component label for this step within a chain.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Image,)
    #: Kind-preserving: an image goes in and an image comes out. The *grid*
    #: changes (the padding is cropped away), which ``PRODUCES`` does not speak
    #: about — ``Resample`` is the established precedent.
    PRODUCES: ClassVar[type | None] = None
    #: The fourth capability flag (W2.12), declared rather than inherited.
    BACKEND: ClassVar[str] = "reference"

    #: FWHM of a Gaussian in units of its standard deviation.
    FWHM_PER_SIGMA: ClassVar[float] = 2.0 * math.sqrt(2.0 * math.log(2.0))

    def __init__(
        self,
        x: Any,
        y: Any,
        *,
        kernel: Any = None,
        fwhm: Any = None,
        label: str | None = None,
    ) -> None:
        super().__init__(label=label)
        if (kernel is None) == (fwhm is None):
            raise TransformationError(
                "PSFConvolution takes exactly one of kernel (a pixel PSF tabulated at the "
                "observed pixel scale) or fwhm (a circular Gaussian's full width at half "
                "maximum, mas); it was given " + ("both" if fwhm is not None else "neither") + "."
            )
        x_grid = _to_mas(x).reshape(-1)
        y_grid = _to_mas(y).reshape(-1)
        steps = (
            _uniform_step(x_grid, "the observed x axis"),
            _uniform_step(y_grid, "the observed y axis"),
        )
        self.register_buffer("x_target", x_grid, unit=COORDINATE_UNIT)
        self.register_buffer("y_target", y_grid, unit=COORDINATE_UNIT)
        support: tuple[int, int]
        if kernel is not None:
            array = np.asarray(kernel, dtype=DTYPE)
            if array.ndim != 2 or array.shape[0] % 2 != 1 or array.shape[1] % 2 != 1:
                raise TransformationError(
                    f"PSFConvolution's kernel is a two-dimensional pixel array with an odd extent "
                    f"on both axes, so that one pixel is its centre; got shape {array.shape}."
                )
            if not np.all(np.isfinite(array)):
                raise TransformationError("PSFConvolution's kernel has non-finite entries.")
            total = float(array.sum())
            if not np.isfinite(total) or total <= 0.0:
                raise TransformationError(
                    f"PSFConvolution normalises its kernel to unit sum so that the step conserves "
                    f"flux, which needs a strictly positive total; this one sums to {total!r}."
                )
            self.register_buffer("psf_kernel", array / total)
            support = (array.shape[0] // 2, array.shape[1] // 2)
        else:
            width = float(_to_mas(fwhm))
            if not np.isfinite(width) or width <= 0.0:
                raise TransformationError(
                    f"PSFConvolution's fwhm must be finite and positive, in mas, got {fwhm!r}."
                )
            self.register_buffer("fwhm", width, unit=COORDINATE_UNIT)
            sigma = width / self.FWHM_PER_SIGMA
            spans = [max(math.ceil(_PSF_PADDING_SIGMAS * sigma / step), 1) for step in steps]
            support = (spans[0], spans[1])
        self._support: tuple[int, int] = support
        self._steps: tuple[float, float] = steps

    # -- construction from the data -----------------------------------------

    @classmethod
    def from_observed(
        cls,
        container: Image,
        *,
        kernel: Any = None,
        fwhm: Any = None,
        label: str | None = None,
        **placement: Any,
    ) -> Any:
        """Take the target grid from the observed image. **The supported route.**

        ``transformations.md`` §10's coordinates-from-the-container rule, in
        the form every modality's step states it: the pixel centres this step
        reproduces are the observed container's own, read off it by name, never
        recomputed from a field of view and a pixel count. ``check_alignment``
        compares axes with :func:`numpy.array_equal`, so a grid recomputed to
        the same numbers in a different order of operations fails a check whose
        message is about axes rather than about arithmetic.
        """
        if not isinstance(container, Image):
            raise TransformationError(
                f"PSFConvolution.from_observed takes the target grid from an Image, got a "
                f"{type(container).__name__}."
            )
        return cls(
            container.x.values * (container.x.unit or COORDINATE_UNIT),
            container.y.values * (container.y.unit or COORDINATE_UNIT),
            kernel=kernel,
            fwhm=fwhm,
            label=label,
            **placement,
        )

    # -- the declared surface ------------------------------------------------

    def _data(self, name: str) -> np.ndarray:
        """A declared buffer as a plain array (``instrument.py``'s ``_Step``)."""
        return np.asarray(self.buffers[name].value, dtype=DTYPE)

    def support(self) -> tuple[int, int]:
        """Half-support of the kernel in pixels, ``(x, y)``.

        A method rather than a property because it is also what
        :func:`~ampere.core.propagate_mask_grid` is handed, and naming it the
        same thing in both places is what keeps the padding this step asks for
        and the dilation it applies to a mask from drifting apart.
        """
        return self._support

    def pixel_scale(self) -> tuple[float, float]:
        """The observed pixel scale, mas, on each axis."""
        return self._steps

    def width(self) -> float | None:
        """The Gaussian FWHM in mas, or ``None`` for a tabulated kernel.

        A method rather than a property for ``parameters.md`` §10's reason: the
        declared name ``fwhm`` is a buffer's (or, after ``promote_buffer``, a
        parameter's), and neither may shadow a class attribute.
        """
        if "fwhm" in self.buffers:
            return float(self.buffers["fwhm"].value)
        if "fwhm" in self.parameters:
            return float(self.parameters["fwhm"].value)
        return None

    def target_grid(self) -> tuple[np.ndarray, np.ndarray]:
        """The observed pixel centres this step emits on, mas."""
        return self._data("x_target"), self._data("y_target")

    def padded_grid(self) -> tuple[np.ndarray, np.ndarray]:
        """The evenly spaced grid this step asks the model for, mas.

        The observed pixel centres extended by :meth:`support` pixels on each
        side, at the observed pixel scale. Built with :func:`numpy.linspace`
        over the padded span rather than by stepping outwards from the observed
        coordinates, because that is exactly how
        :meth:`~ampere.core.AxisRequirement.coordinates` builds a grid from an
        interval and a ``max_step``: publishing the two forms of the same
        constraint together only works if they agree to the last bit, and the
        cheapest way to make them agree is to compute them the same way.
        """
        built = []
        for grid, step, pad in zip(self.target_grid(), self._steps, self._support, strict=True):
            low = float(grid[0]) - pad * step
            high = float(grid[-1]) + pad * step
            built.append(np.linspace(low, high, grid.size + 2 * pad))
        return built[0], built[1]

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """The padded grid, as coordinates *and* as a coverage-plus-density pair.

        Both halves, on both axes, and each says something the other cannot.

        ``points=`` is the operative half. A convolution's output pixel is a
        weighted sum over its neighbours, so the step must be handed the
        observed pixel centres *and* their neighbourhood, at the observed pixel
        scale; a grid merely dense enough over a wide enough interval would not
        in general contain the observed centres at all, and
        :meth:`~ampere.core.Axis.locate` — which is how :meth:`apply` finds
        them in whatever union negotiation produced — would then, correctly,
        refuse. This is ``spectrum_photometry.md`` Gap 1's shape in two
        dimensions: a buffer tabulated on particular coordinates publishes them
        and looks them up.

        ``intervals``/``max_step`` restate the same pixel scale as a *density*,
        which is what makes a model holding its own fixed grid refuse by name
        when that grid is coarser than the data
        (``_ImageModel(adopt_grid=False)``, gap I-4's mechanism applied to an
        image that is being *observed* rather than Fourier-sampled). Without
        it, the same model would sail through composition and fail much later,
        inside ``locate``, with a message about a missing coordinate.

        Note what is **not** here: a ``configure_from`` hook. An LSF has one
        because it genuinely does not know its own output range — the resampler
        or the photometry step downstream sets it. A PSF convolution does know:
        its output range is the observed image, which it was built from.
        """
        published = []
        for name, grid, padded, step in zip(
            ("x", "y"), self.target_grid(), self.padded_grid(), self._steps, strict=True
        ):
            published.append(
                AxisRequirement(
                    name,
                    unit=COORDINATE_UNIT,
                    intervals=(float(padded[0]), float(padded[-1])),
                    max_step=step * (1.0 + _STEP_SLACK),
                    points=padded,
                    source=f"{self.label} (PSF convolution, {grid.size} observed pixels)",
                )
            )
        return tuple(published)

    def crop_indices(self, x_mas: Any, y_mas: Any) -> tuple[np.ndarray, np.ndarray]:
        """Where the observed pixel centres sit in the negotiated grid ``(x, y)``.

        One rule in one place, because all three backends need the same answer
        and a lookup written three times is a lookup that eventually differs.
        :meth:`~ampere.core.Axis.locate` is the rule — matching within
        :data:`~ampere.core.COORDINATE_RTOL`, which is exactly the tolerance
        negotiation's union collapses coincident coordinates at, so a centre
        this step published and a centre that survived the union are the same
        centre. It raises, by name, when they are not.
        """
        rows = Axis.build("x", np.asarray(x_mas, dtype=DTYPE), COORDINATE_UNIT)
        columns = Axis.build("y", np.asarray(y_mas, dtype=DTYPE), COORDINATE_UNIT)
        return rows.locate(self._data("x_target")), columns.locate(self._data("y_target"))

    # -- the kernel ----------------------------------------------------------

    def grid_steps(self, x_mas: Any, y_mas: Any) -> tuple[float, float]:
        """The incoming grid's pixel scale, mas, once it has been checked.

        Two checks, and both are refusals rather than fallbacks. The grid must
        be **evenly spaced**, because the FFT route needs it to be
        (:func:`_uniform_step`). And a **tabulated** kernel's own pixel scale
        must be the grid's, because a measured PSF is a buffer tied to the
        coordinates it was tabulated on — ``transformations.md`` §10's second
        rule — and rescaling one would invent information it does not carry.
        An analytic kernel has neither problem and is simply retabulated.

        Shared with the native twins, which need the same two answers in the
        same two words before their own arithmetic starts.
        """
        steps = (
            _uniform_step(np.asarray(x_mas, dtype=DTYPE).reshape(-1), "the negotiated x axis"),
            _uniform_step(np.asarray(y_mas, dtype=DTYPE).reshape(-1), "the negotiated y axis"),
        )
        if "psf_kernel" in self.buffers:
            for axis, (given, observed) in enumerate(zip(steps, self._steps, strict=True)):
                if abs(given - observed) > _SCALE_RTOL * abs(observed):
                    raise TransformationError(
                        f"PSFConvolution was built with a tabulated kernel, which is a buffer "
                        f"tied to the pixel scale it was tabulated at ({observed:g} mas on axis "
                        f"{'xy'[axis]}), and the grid it has been handed is sampled at "
                        f"{given:g} mas. Rescaling a measured PSF would invent information it "
                        f"does not carry. Either supply the kernel at the grid's own scale, or "
                        f"use the fwhm= form, whose kernel is analytic and can be retabulated."
                    )
        return steps

    def psf(self, x_mas: Any, y_mas: Any, values: Any = None) -> np.ndarray:
        """The normalised kernel, tabulated on the pixel scale of ``(x, y)``.

        Public because the conformance battery's direct-convolution row needs
        the *same* kernel the FFT route uses, so that the row measures the
        transform rather than two different Gaussians. What it must not share
        with that row is the convolution itself, which the row computes for
        itself as a plain sum over the support.
        """
        steps = self.grid_steps(x_mas, y_mas)
        if "psf_kernel" in self.buffers:
            return self._data("psf_kernel")
        width = float(self.context(values)["fwhm"])
        if not np.isfinite(width) or width <= 0.0:
            raise TransformationError(
                f"PSFConvolution's fwhm must be finite and positive at every evaluation, in mas; "
                f"this draw gave {width!r}. A fitted width wants a prior on the positive line "
                f"(promote_buffer takes a bijection for exactly this)."
            )
        sigma = width / self.FWHM_PER_SIGMA
        profiles = []
        for step, pad in zip(steps, self._support, strict=True):
            offsets = np.arange(-pad, pad + 1, dtype=DTYPE) * step
            profiles.append(np.exp(-0.5 * (offsets / sigma) ** 2))
        kernel = np.outer(profiles[0], profiles[1])
        # Renormalised over the support actually tabulated, so the step
        # conserves flux exactly whatever the truncation costs in the wings.
        return np.asarray(kernel / float(kernel.sum()), dtype=DTYPE)

    # -- the contract surface ------------------------------------------------

    def apply(self, samples: Any, values: Any) -> Image:
        """Convolve on the negotiated grid, then crop to the observed pixels.

        In that order, and not the other way round: cropping first and
        convolving after is the boundary error the padding exists to prevent.
        """
        x_in = np.asarray(samples.x.values, dtype=DTYPE)
        y_in = np.asarray(samples.y.values, dtype=DTYPE)
        kernel = self.psf(x_in, y_in, values)
        blurred = _fft_convolve(np.asarray(samples.values, dtype=DTYPE), kernel)
        selection = np.ix_(*self.crop_indices(x_in, y_in))
        # The dilation happens on the *input* grid, where the masked pixels
        # are, and only then is it cropped: a masked pixel in the padding
        # legitimately contaminates an observed pixel near the edge.
        mask = propagate_mask_grid(samples, (kernel.shape[0] // 2, kernel.shape[1] // 2))
        return Image(
            self._data("x_target") * COORDINATE_UNIT,
            self._data("y_target") * COORDINATE_UNIT,
            blurred[selection],
            unit=samples.unit,
            mask=None if mask is None else mask[selection],
        )

    def describe(self) -> dict[str, Any]:
        """The configuration provenance cannot reach through parameters alone.

        ``parameters.md``'s opt-in hook: the half-support and the pixel scale
        are derived at construction and are not buffers, yet two steps that
        differ in either compute different numbers, so a cache key that could
        not see them would be wrong.
        """
        return {"support": list(self._support), "pixel_scale": list(self._steps)}
