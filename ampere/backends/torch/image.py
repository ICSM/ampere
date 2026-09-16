"""Gridded images on the torch path: the PSF-convolution twin.

The native twin of :mod:`ampere.backends.reference.image` (W5.5), and the
gradient that makes an image dataset fittable under NUTS: a sky model's
parameters reach the data through a two-dimensional convolution, so a
derivative of that convolution with respect to them is what a Hamiltonian
sampler needs and what this module supplies.

Which twin pattern, and why
---------------------------
The **inheriting** one, as
:mod:`ampere.backends.torch.interferometry` uses and for the same reason,
which is worth restating because the interferometry page's rule has three
clauses and this step tests which fires. Here the *declaration* is the
expensive, dangerous half: :meth:`~ampere.backends.reference.image.PSFConvolution.requirements`
publishes a padded grid as ``points=`` **and** the observed pixel scale as a
``max_step``, and those two have to agree to the last bit or the union they go
into stops being a single evenly spaced grid — at which point the FFT route
this class exists for refuses, on a backend, for a reason that has nothing to
do with torch. The arithmetic, by contrast, is one transform and one crop.
That is clause one exactly, so these inherit and override only the four
capability flags and the sums.

What is forked is every per-evaluation number: the kernel, the transform and
the crop all happen in ``torch``, so a conformance row driving ``apply``
compares torch's convolution against a direct sum rather than numpy's against
numpy. :class:`TorchImageStep` keeps this backend's own plumbing (``dtype``,
``device``, the ``DEVICE`` flag, the ``tensors`` buffer home and ``to()``), so
an image step behaves like every other torch step under a device move.

Batching, and the leading ellipsis
----------------------------------
``BATCHABLE`` is true here and it means what it means everywhere else on this
backend: ``torch.func.vmap`` over the realised density evaluates a stack of
parameter vectors in one call. Every expression below is whole-tensor and acts
on the last two axes explicitly (``dim=(-2, -1)`` on the transforms, an
``...``-leading index on the crop), so a batch axis rides through untouched and
nothing branches on a value.
"""

from __future__ import annotations

from typing import Any, ClassVar

import numpy as np
import torch

from ampere.backends.reference.image import COORDINATE_UNIT
from ampere.backends.reference.image import (
    PSFConvolution as _ReferencePSFConvolution,
)
from ampere.core import DTYPE, Image, TransformationError, propagate_mask_grid

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, move, place, to_numpy
from .parameters import LoweredParameters

__all__ = ["PSFConvolution", "TorchImageStep"]


def _next_fast_len(length: int) -> int:
    """The next 5-smooth integer at least *length*.

    torch has no ``next_fast_len``, and ``torch.fft`` is happy with any size —
    but a transform on a length with a large prime factor is markedly slower,
    and the reference backend pads to a fast length through ``scipy.fft``. The
    two backends must pad to the **same** length or their transforms would
    differ by the rounding of the implicit zero padding, which is exactly the
    kind of cross-backend disagreement the conformance rows exist to catch, so
    this reproduces scipy's rule rather than picking its own.
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


def _fft_convolve(image: torch.Tensor, kernel: torch.Tensor) -> torch.Tensor:
    """Zero-padded linear convolution over the last two axes, central block.

    The reference module's ``_fft_convolve``, in torch, with the batch axes
    carried through: ``torch.fft.rfft2`` transforms the last two dimensions and
    broadcasts over everything in front of them, so a ``(batch, nx, ny)`` stack
    costs one call rather than a loop.
    """
    n_x, n_y = image.shape[-2], image.shape[-1]
    k_x, k_y = kernel.shape[-2], kernel.shape[-1]
    fast = (_next_fast_len(n_x + k_x - 1), _next_fast_len(n_y + k_y - 1))
    spectrum = torch.fft.rfft2(image, s=fast) * torch.fft.rfft2(kernel, s=fast)
    full = torch.fft.irfft2(spectrum, s=fast)
    return full[..., k_x // 2 : k_x // 2 + n_x, k_y // 2 : k_y // 2 + n_y]


class TorchImageStep:
    """The four capability flags, this backend's placement, and the tensors.

    A **mixin** rather than a base class, so that each concrete step still
    derives from its reference counterpart and inherits that class's whole
    declaration — the constructor validation, the padded grid, the published
    requirements, ``from_observed`` and the crop rule. See the module docstring
    for why that is the right split here.
    """

    #: The four capability flags (``inference.md`` §18), stated rather than
    #: inherited: the reference class these derive from declares
    #: ``BACKEND = "reference"`` and ``DIFFERENTIABLE = False``, and inheriting
    #: that would be this backend claiming to be another one.
    DIFFERENTIABLE: ClassVar[bool] = True
    BATCHABLE: ClassVar[bool] = True
    #: Class default only; ``__init__`` shadows it per instance.
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = BACKEND

    def __init__(
        self,
        *args: Any,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
        **kwargs: Any,
    ) -> None:
        # The reference constructor first: it validates the arguments and
        # registers the declared buffers these tensors are built from.
        super().__init__(*args, **kwargs)
        place(self, dtype, device)
        #: The torch-side home for this step's constant arrays
        #: (``lowering.md`` §7), so ``.to(...)`` is a real move.
        self.tensors = LoweredParameters()
        for name in self.buffers.names:
            self.tensors.register_buffer(
                name,
                as_tensor(
                    np.asarray(self.buffers[name].value, dtype=float),
                    dtype=self.dtype,
                    device=self.device,
                ),
                persistent=True,
            )

    def to(self, *, dtype: Any = None, device: Any = None) -> TorchImageStep:
        """Move this step's constant tensors, and re-declare where they live."""
        return move(self, dtype=dtype, device=device)

    def _real(self, value: Any) -> torch.Tensor:
        """*value* as a real tensor where this step's arithmetic happens."""
        return as_tensor(value, dtype=self.dtype, device=self.device)

    def _tensor(self, name: str) -> torch.Tensor:
        return self.tensors.get_buffer(name)


class PSFConvolution(TorchImageStep, _ReferencePSFConvolution):
    """Convolve a model image with a point-spread function, in torch.

    ``Image -> Image``, kind-preserving and grid-cropping, with the transform
    and the crop both differentiable with respect to the model's fluxes and —
    for the ``fwhm=`` form — with respect to the width itself, which is what
    makes a fitted seeing or beam size an ordinary parameter here rather than a
    special case.

    Everything declarative is the reference class's, unchanged and deliberately
    so: the target-grid and kernel buffers, ``from_observed``, the padded grid,
    the published requirements and the crop rule. See the module docstring.

    Parameters
    ----------
    x, y, kernel, fwhm, label
        As :class:`ampere.backends.reference.image.PSFConvolution`.
    dtype, device
        Where this step's constant tensors live, and in what precision. Never
        detected (``architecture.md`` §5).
    """

    def psf_tensor(self, x_mas: Any, y_mas: Any, values: Any = None) -> torch.Tensor:
        """The normalised kernel as a tensor, on the pixel scale of ``(x, y)``.

        The tabulated form is a constant, so it comes straight off
        :attr:`tensors`; the Gaussian form is rebuilt from
        :meth:`~ampere.core.Parameterised.context` on every call, in torch, so
        that a promoted ``fwhm`` carries a gradient. The *support* is the
        reference class's — fixed at composition from the declared width, which
        is what keeps the padding this step asked negotiation for and the
        kernel it actually builds the same size.
        """
        # The two refusals — an irregular grid, and a tabulated kernel handed a
        # different pixel scale — are the reference class's, called here so
        # that all three backends say the same words about the same mistake.
        steps = self.grid_steps(to_numpy(x_mas), to_numpy(y_mas))
        if "psf_kernel" in self.buffers:
            return self._tensor("psf_kernel")
        width = self._real(self.context(values)["fwhm"])
        sigma = width / self.FWHM_PER_SIGMA
        profiles = []
        for step, pad in zip(steps, self.support(), strict=True):
            offsets = torch.arange(-pad, pad + 1, dtype=self.dtype, device=self.device) * step
            profiles.append(torch.exp(-0.5 * (offsets / sigma) ** 2))
        kernel = torch.outer(profiles[0], profiles[1])
        return kernel / kernel.sum()

    # -- the native surface a realisation composes ---------------------------

    def apply_flux(self, flux: torch.Tensor, grid: Any, values: Any) -> tuple[torch.Tensor, Any]:
        """``(..., nx, ny) -> (..., nx_obs, ny_obs)``, and the observed grid.

        *grid* is the image model's ``(x, y)`` pair in mas; what comes back is
        the observed pair, because this step crops the padding away.
        """
        x_mas, y_mas = grid
        kernel = self.psf_tensor(x_mas, y_mas, values)
        blurred = _fft_convolve(flux, kernel)
        rows, columns = self.crop_indices(to_numpy(x_mas), to_numpy(y_mas))
        index_x = torch.as_tensor(np.asarray(rows), device=self.device)
        index_y = torch.as_tensor(np.asarray(columns), device=self.device)
        cropped = blurred.index_select(-2, index_x).index_select(-1, index_y)
        return cropped, (self._tensor("x_target"), self._tensor("y_target"))

    # -- the contract surface ------------------------------------------------

    def apply(self, samples: Any, values: Any) -> Image:
        """The same arithmetic, back on the numpy side of the container boundary.

        Overridden rather than inherited so that the conformance battery —
        which drives ``apply`` — measures *this* backend's convolution against
        a direct sum, not the reference backend's against itself.
        """
        x_in = np.asarray(samples.x.values, dtype=DTYPE)
        y_in = np.asarray(samples.y.values, dtype=DTYPE)
        kernel = self.psf_tensor(self._real(x_in), self._real(y_in), values)
        if not bool(torch.all(torch.isfinite(kernel))):
            raise TransformationError(
                "PSFConvolution built a non-finite kernel; a fitted fwhm must stay finite and "
                "positive at every draw."
            )
        blurred = _fft_convolve(self._real(np.asarray(samples.values, dtype=DTYPE)), kernel)
        selection = np.ix_(*self.crop_indices(x_in, y_in))
        mask = propagate_mask_grid(samples, (kernel.shape[-2] // 2, kernel.shape[-1] // 2))
        return Image(
            self._data("x_target") * COORDINATE_UNIT,
            self._data("y_target") * COORDINATE_UNIT,
            to_numpy(blurred)[selection],
            unit=samples.unit,
            mask=None if mask is None else mask[selection],
        )
