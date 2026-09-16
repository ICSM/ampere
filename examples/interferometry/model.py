"""The correct sky model: a binary, plus a fixed, fainter disc, on one grid.

``ampere`` ships no "add two image models together" operator (there is
nothing analogous to :class:`~ampere.core.kernels.Sum` for
:class:`~ampere.core.transform.Model`), and this module does not add one to
``ampere`` — the item's file ownership keeps this study out of ``ampere/``.
Instead it does what a user is meant to do: write a small
:class:`~ampere.core.transform.Model` of its own,
exactly as :mod:`examples.m2_misspecification.model` writes ``AbsorptionLines``
from the public ``ampere.core`` surface rather than reaching into a backend's
private machinery.

:class:`BinaryWithDisc` is not, however, a from-scratch reimplementation of
an image model's negotiation and Nyquist bookkeeping — that machinery
(``compile_for``'s grid adoption, the gap I-4 refusal, the native
``native_grid``/``native_flux`` surface torch and jax realise through) is
substantial and already correct in
:class:`~ampere.backends.reference.interferometry.Binary` and
:class:`~ampere.backends.reference.interferometry.GaussianSource`
(and their per-backend twins). This class **delegates** to one instance of
each — a binary with the study's two free parameters, and a disc with none
at all — and adds their brightness together on the grid both negotiate
identically, because both are given the same ``requirements`` mapping in
:meth:`compile_for`.

Why the disc has no free parameters
------------------------------------
Arm "correct" of the study asks whether a fit that is *told* the disc is
there, exactly as it is, recovers the binary honestly. That is a claim about
the binary's two parameters, not about whether the disc itself can be
fitted — giving the disc a prior would spend part of the sampler's budget on
a question the study is not asking. So every disc parameter here is
constructed from a plain number rather than a distribution, which
:func:`~ampere.backends.reference._declare.as_parameter` turns into a
**fixed** :class:`~ampere.core.Parameter` — present in the posterior's
declared vocabulary, absent from what the sampler ever moves.

One class, three backends
--------------------------
``BinaryWithDisc`` takes the backend's own ``interferometry`` module (*itf*)
as its first argument rather than being rewritten per backend the way
:mod:`examples.m2_misspecification`'s toy model is — that study's three
files exist because its physics is hand-written per array library; this
class has no arithmetic of its own to fork; it is grid negotiation and
addition, and both operations are exactly the same call whichever backend's
:class:`~ampere.core.transform.Model` it delegates to. :mod:`.model_torch`
and :mod:`.model_jax` exist to keep the item's layout (and so that
``python -m examples.interferometry`` need not import torch or jax to
build the reference variant) and each is a two-line binding of this class to
its own backend's module.
"""

from __future__ import annotations

from typing import Any

from ampere.core import Model, ModelResult

__all__ = ["BinaryWithDisc", "binary_with_disc"]


class BinaryWithDisc(Model):
    """A resolved binary plus a fixed, fainter circular disc, on one grid.

    Parameters
    ----------
    itf
        The backend's ``interferometry`` module
        (:mod:`ampere.backends.reference.interferometry` or a modern
        backend's twin) — supplies the :class:`Binary` and
        :class:`GaussianSource` classes this delegates to.
    x, y
        The model's own placeholder grid, mas — see
        :func:`examples.interferometry.generators.seed_grid`; discarded the
        moment negotiation adopts a union grid, exactly as for the two
        classes this wraps.
    disc_flux, disc_fwhm
        The disc's total flux (Jy) and full width at half maximum (mas), as
        plain numbers — always fixed (see the module docstring).
    separation, position_angle, flux_ratio, flux, component_fwhm
        Forwarded to the wrapped :class:`Binary` unchanged: a prior fits a
        parameter, a number fixes it, exactly as for :class:`Binary` itself.
    channels, adopt_grid
        Forwarded to both wrapped models; both must be bound to the same
        channel(s) for the addition below to mean anything.
    """

    def __init__(
        self,
        itf: Any,
        x: Any,
        y: Any,
        *,
        disc_flux: float,
        disc_fwhm: float,
        separation: Any = 5.0,
        position_angle: Any = 0.5,
        flux_ratio: Any = 0.4,
        flux: Any = 1.0,
        component_fwhm: float = 0.5,
        channels: str = "sky",
        adopt_grid: bool = True,
    ) -> None:
        self._binary = itf.Binary(
            x,
            y,
            separation=separation,
            position_angle=position_angle,
            flux_ratio=flux_ratio,
            flux=flux,
            component_fwhm=component_fwhm,
            channels=channels,
            adopt_grid=adopt_grid,
        )
        self._disc = itf.GaussianSource(
            x, y, fwhm=disc_fwhm, flux=disc_flux, channels=channels, adopt_grid=adopt_grid
        )
        self.channels = self._binary.channels
        # The capability flags (DEVELOPMENT_PLAN.md §4.5) are the *binary's*:
        # the disc has no free parameters, so it cannot be the reason a
        # composed problem is or is not differentiable or batchable. Set as
        # instance attributes, which is what a realisation reads
        # (``model.DIFFERENTIABLE``) regardless of what the shared
        # ``BinaryWithDisc`` class itself declares.
        self.DIFFERENTIABLE = bool(self._binary.DIFFERENTIABLE)
        self.BATCHABLE = bool(self._binary.BATCHABLE)
        self.DEVICE = self._binary.DEVICE
        self.BACKEND = self._binary.BACKEND
        for parameter in self._binary.parameters:
            self.register_parameter(parameter)

    # -- negotiation: both children see the same requirements ---------------

    def compile_for(self, requirements: Any) -> BinaryWithDisc:
        """Negotiate both children identically, so their grids cannot diverge."""
        self._binary = self._binary.compile_for(requirements)
        self._disc = self._disc.compile_for(requirements)
        return self

    # -- the contract surface: two Images, added ------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        binary_result = self._binary.evaluate(**values)
        disc_result = self._disc.evaluate()
        combined = {}
        for channel in self.channels:
            binary_image = binary_result[channel]
            disc_image = disc_result[channel]
            combined[channel] = binary_image.with_values(binary_image.values + disc_image.values)
        return ModelResult(combined)

    # -- the native surface (torch, jax): delegate, then add ------------------

    def native_grid(self, channel: str) -> Any:
        """The negotiated ``(x, y)``, identical for both children by construction."""
        return self._binary.native_grid(channel)

    def native_flux(self, channel: str, values: Any = None) -> Any:
        """The binary's native brightness plus the disc's, on *channel*.

        *values* carries this model's own (the binary's) parameter values,
        exactly as :mod:`ampere.backends.torch.problem` /
        :mod:`ampere.backends.jax.problem` route them; the disc needs none —
        every one of its parameters is fixed, so its own declared value is
        what :meth:`~ampere.core.parameter.Parameterised.context` falls back
        to when *values* is ``None``.
        """
        return self._binary.native_flux(channel, values) + self._disc.native_flux(channel, None)


def binary_with_disc(
    itf: Any, x: Any, y: Any, *, disc_flux: float, disc_fwhm: float, **binary_kwargs: Any
) -> BinaryWithDisc:
    """Convenience constructor: see :class:`BinaryWithDisc`."""
    return BinaryWithDisc(itf, x, y, disc_flux=disc_flux, disc_fwhm=disc_fwhm, **binary_kwargs)
