"""Gridded images on the jax path: the PSF-convolution twin.

The native twin of :mod:`ampere.backends.reference.image` (W5.5), as
:mod:`ampere.backends.torch.image` is on the other modern backend: the same
step, the same declaration, and the transform in ``jax.numpy`` so that a sky
model fitted to an image has a gradient under NUTS and a ``vmap``\\ ped batch
under the SBI path.

Why this subclasses the reference step
--------------------------------------
The reason this module's siblings give (:mod:`ampere.backends.jax.instrument`,
:mod:`ampere.backends.jax.interferometry`), and one that is specific to a PSF.
Here the declaration is both the expensive half and the dangerous half:
:meth:`~ampere.backends.reference.image.PSFConvolution.requirements` publishes
one constraint in two forms — a padded grid as ``points=``, and the observed
pixel scale as a ``max_step`` — and the two must agree to the last bit or the
union negotiation builds from them stops being evenly spaced, at which point
the FFT route refuses for a reason that has nothing to do with jax. The
arithmetic is one transform and one crop. So the declaration is shared and the
sums are forked, which is what makes the cross-backend rows a comparison of
two implementations rather than of one implementation with itself.

Two evaluation surfaces, as everywhere in this backend
------------------------------------------------------
``apply`` is the contract surface and returns an ``ampere.core`` container, so
a gradient stops at the container boundary; it is overridden here rather than
inherited so that a conformance row driving ``apply`` measures *this* backend's
convolution. ``apply_flux`` is the native surface — coordinates and values in,
the transformed pair out, pure and traceable — and what
:mod:`ampere.backends.jax.problem` composes into the differentiable
log-density.

The kernel is a function of the coordinates and, for the ``fwhm=`` form, of one
fitted value; the crop is a function of the coordinates alone. The coordinates
are fixed for a run, so the crop indices are cached against the grid's exact
bytes — the same rule, and the same reasoning,
:meth:`ampere.backends.jax.instrument._JaxStep._influence_jax` records for an
influence matrix.
"""

from __future__ import annotations

from typing import Any, ClassVar

import jax
import jax.numpy as jnp
import numpy as np

from ampere.backends.reference.image import COORDINATE_UNIT
from ampere.backends.reference.image import (
    PSFConvolution as _ReferencePSFConvolution,
)
from ampere.core import DTYPE, Image, TransformationError, propagate_mask_grid

from ._config import BACKEND, require_x64
from ._device import DEVICE, device_flag, place_on, resolve_device

__all__ = ["PSFConvolution"]


def _next_fast_len(length: int) -> int:
    """The next 5-smooth integer at least *length*.

    jax has no ``next_fast_len`` and scipy's is not traceable, so the rule is
    reproduced here. It has to be the **same** rule the reference backend's
    ``scipy.fft.next_fast_len`` gives, because the padded length decides how
    many implicit zeros the transform sees and therefore, at the last bit, what
    it returns; two backends padding differently would show up as a
    cross-backend disagreement with no cause visible in either.
    """
    if length <= 16:
        return max(int(length), 1)
    candidate = int(length)
    while True:
        remaining = candidate
        for factor in (2, 3, 5):
            while remaining % factor == 0:
                remaining //= factor
        if remaining == 1:
            return candidate
        candidate += 1


def _fft_convolve(image: jax.Array, kernel: jax.Array) -> jax.Array:
    """Zero-padded linear convolution over the last two axes, central block.

    ``jnp.fft.rfft2`` transforms the last two axes and broadcasts over whatever
    is in front of them, so a ``vmap``\\ ped batch costs one call. The padded
    length comes from :func:`_next_fast_len` and is a Python integer computed
    from static shapes, so it is a compile-time constant rather than a traced
    value — which is what keeps this function jittable.
    """
    n_x, n_y = image.shape[-2], image.shape[-1]
    k_x, k_y = kernel.shape[-2], kernel.shape[-1]
    fast = (_next_fast_len(n_x + k_x - 1), _next_fast_len(n_y + k_y - 1))
    spectrum = jnp.fft.rfft2(image, s=fast) * jnp.fft.rfft2(kernel, s=fast)
    full = jnp.fft.irfft2(spectrum, s=fast)
    return full[..., k_x // 2 : k_x // 2 + n_x, k_y // 2 : k_y // 2 + n_y]


class PSFConvolution(_ReferencePSFConvolution):
    """Convolve a model image with a point-spread function, in jax.

    ``Image -> Image``, kind-preserving and grid-cropping. Everything
    declarative is the reference class's, unchanged: the target-grid and kernel
    buffers, ``from_observed``, the padded grid, the published requirements,
    the pixel-scale and regularity refusals, and the crop rule. What is
    overridden is the four capability flags and the arithmetic.

    Parameters
    ----------
    x, y, kernel, fwhm, label
        As :class:`ampere.backends.reference.image.PSFConvolution`.
    device
        Where this step's constant arrays live. Never auto-detected
        (``architecture.md`` §5).
    """

    #: The gradient the whole backend exists for passes through
    #: :meth:`apply_flux`.
    DIFFERENTIABLE: ClassVar[bool] = True
    #: Whole-array ``jax.numpy`` throughout, with a leading ellipsis on the
    #: transform and the crop, so ``jax.vmap`` maps this as it maps any pure
    #: function.
    BATCHABLE: ClassVar[bool] = True
    #: Class default only; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = DEVICE
    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        require_x64(f"a jax {type(self).__name__}")
        # Popped before the reference constructor sees it: ``device=`` is this
        # backend's placement keyword and means nothing to the numpy step whose
        # declaration is inherited wholesale.
        device = kwargs.pop("device", DEVICE)
        super().__init__(*args, **kwargs)
        resolved = resolve_device(device, f"a jax {type(self).__name__}", error=TransformationError)
        self._device = resolved
        object.__setattr__(self, "DEVICE", device_flag(device, resolved))
        self._cached_crop: tuple[bytes, tuple[jax.Array, jax.Array]] | None = None

    def _place(self, value: Any) -> jax.Array:
        """*value* as a float64 jax array on this step's device."""
        return place_on(jnp.asarray(value, dtype=jnp.float64), getattr(self, "_device", None))

    def crop_arrays(self, x_mas: Any, y_mas: Any) -> tuple[jax.Array, jax.Array]:
        """:meth:`crop_indices` as device arrays, cached on the grid's bytes.

        The indices are a pure function of the negotiated grid, the grid is
        fixed for a run, and recomputing a ``locate`` against it on every
        likelihood evaluation is the hot-loop waste W2.3 measured. The key is
        exact, so a changed grid rebuilds and a repeated one does not.
        """
        x_array = np.asarray(_values_of(x_mas), dtype=DTYPE).reshape(-1)
        y_array = np.asarray(_values_of(y_mas), dtype=DTYPE).reshape(-1)
        key = x_array.tobytes() + b"|" + y_array.tobytes()
        cached = self._cached_crop
        if cached is not None and cached[0] == key:
            return cached[1]
        rows, columns = self.crop_indices(x_array, y_array)
        built = (
            place_on(jnp.asarray(rows), getattr(self, "_device", None)),
            place_on(jnp.asarray(columns), getattr(self, "_device", None)),
        )
        self._cached_crop = (key, built)
        return built

    def psf_array(self, x_mas: Any, y_mas: Any, values: Any = None) -> jax.Array:
        """The normalised kernel as a jax array, on the pixel scale of ``(x, y)``.

        The tabulated form is a constant placed once; the Gaussian form is
        rebuilt from :meth:`~ampere.core.Parameterised.context` on every call,
        in jax, so that a promoted ``fwhm`` carries a gradient. The *support*
        is the reference class's — fixed at composition from the declared
        width — which is what keeps the padding this step asked negotiation for
        and the kernel it actually builds the same size.
        """
        # The two refusals — an irregular grid, and a tabulated kernel handed a
        # different pixel scale — are the reference class's, called here so
        # that all three backends say the same words about the same mistake.
        steps = self.grid_steps(_values_of(x_mas), _values_of(y_mas))
        if "psf_kernel" in self.buffers:
            return self._place(self._data("psf_kernel"))
        width = jnp.asarray(self.context(values)["fwhm"], dtype=jnp.float64)
        sigma = width / self.FWHM_PER_SIGMA
        profiles = []
        for step, pad in zip(steps, self.support(), strict=True):
            offsets = self._place(np.arange(-pad, pad + 1, dtype=DTYPE) * step)
            profiles.append(jnp.exp(-0.5 * (offsets / sigma) ** 2))
        kernel = jnp.outer(profiles[0], profiles[1])
        return kernel / kernel.sum()

    # -- the native surface a realisation composes ---------------------------

    def apply_flux(self, flux: jax.Array, grid: Any, values: Any) -> tuple[jax.Array, Any]:
        """``(..., nx, ny) -> (..., nx_obs, ny_obs)``, and the observed grid."""
        x_mas, y_mas = grid
        kernel = self.psf_array(x_mas, y_mas, values)
        blurred = _fft_convolve(flux, kernel)
        rows, columns = self.crop_arrays(x_mas, y_mas)
        cropped = jnp.take(jnp.take(blurred, rows, axis=-2), columns, axis=-1)
        return cropped, (self._place(self._data("x_target")), self._place(self._data("y_target")))

    # -- the contract surface ------------------------------------------------

    def apply(self, samples: Any, values: Any) -> Image:
        """The same arithmetic, back on the numpy side of the container boundary."""
        x_in = np.asarray(samples.x.values, dtype=DTYPE)
        y_in = np.asarray(samples.y.values, dtype=DTYPE)
        kernel = self.psf_array(x_in, y_in, values)
        blurred = _fft_convolve(self._place(np.asarray(samples.values, dtype=DTYPE)), kernel)
        selection = np.ix_(*self.crop_indices(x_in, y_in))
        mask = propagate_mask_grid(samples, (kernel.shape[-2] // 2, kernel.shape[-1] // 2))
        return Image(
            self._data("x_target") * COORDINATE_UNIT,
            self._data("y_target") * COORDINATE_UNIT,
            np.asarray(blurred)[selection],
            unit=samples.unit,
            mask=None if mask is None else mask[selection],
        )


def _values_of(source: Any) -> Any:
    """Bare coordinates, whether *source* is an ``Axis``, an array or a tracer."""
    return np.asarray(getattr(source, "values", source), dtype=DTYPE)
