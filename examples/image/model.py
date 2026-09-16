"""The sky model: a compact source, plus an optional smooth background.

``ampere`` ships no "add two image models together" operator, and this module
does not add one — the item's file ownership keeps this study out of
``ampere/``. It does what a user is meant to do instead, and what
:mod:`examples.interferometry.model` already demonstrates for the same kind of
composite: write a small :class:`~ampere.core.transform.Model` that
**delegates** to two of the shipped image models and adds their brightness.

Delegation rather than reimplementation, for the reason that page records: the
negotiation machinery an image model carries (``compile_for``'s adoption of the
negotiated grid, the refusal when a model holds its own,
``native_grid``/``native_flux``) is substantial and already correct in
:class:`~ampere.backends.reference.interferometry.GaussianSource` and its two
twins. This class contributes one addition.

Why the background has no free parameters
------------------------------------------
Arm ``"correct"`` asks whether a fit that is *told* the background is there,
exactly as it is, recovers the source honestly. That is a claim about the
source's two parameters. Giving the background a prior would spend the
sampler's budget on a question the study is not asking, so every background
parameter is a plain number, which
:func:`~ampere.backends.reference._declare.as_parameter` turns into a fixed
parameter: in the posterior's declared vocabulary, never moved.

``background_flux=None`` drops the background model entirely rather than
setting its flux to zero. That is not a stylistic choice: a zero-flux Gaussian
still contributes an array of zeros to every evaluation and a fixed parameter
to every provenance record, so the misspecified arm would differ from a
genuinely background-free model in its spec hash while agreeing in its numbers.
The misspecification this study is about is a model that *does not have* the
component, so that is what the class builds.
"""

from __future__ import annotations

from typing import Any

from ampere.core import Model, ModelResult

__all__ = ["SourceWithBackground", "source_with_background"]


class SourceWithBackground(Model):
    """A compact Gaussian source, plus a fixed broader one, on one grid.

    Parameters
    ----------
    itf
        The backend's ``interferometry`` module
        (:mod:`ampere.backends.reference.interferometry` or a twin) — supplies
        the :class:`GaussianSource` class this delegates to. The image models
        live there because W4.1 put them there; W5.5 reuses them unchanged
        rather than moving them.
    x, y
        The model's own placeholder grid, mas. Discarded the moment
        negotiation adopts the PSF step's padded grid.
    flux, fwhm
        The compact source's total flux (Jy) and full width at half maximum
        (mas). A prior fits the parameter, a number fixes it.
    background_flux, background_fwhm
        The smooth component's, as plain numbers — always fixed (see the module
        docstring). ``background_flux=None`` builds no background at all, which
        is the misspecified arms' model.
    channels, adopt_grid
        Forwarded to both children; both must bind the same channel for the
        addition to mean anything.
    """

    def __init__(
        self,
        itf: Any,
        x: Any,
        y: Any,
        *,
        flux: Any = 1.0,
        fwhm: Any = 3.0,
        background_flux: float | None = None,
        background_fwhm: float = 17.0,
        channels: str = "sky",
        adopt_grid: bool = True,
    ) -> None:
        self._source = itf.GaussianSource(
            x, y, flux=flux, fwhm=fwhm, channels=channels, adopt_grid=adopt_grid
        )
        self._background = (
            None
            if background_flux is None
            else itf.GaussianSource(
                x,
                y,
                flux=float(background_flux),
                fwhm=float(background_fwhm),
                channels=channels,
                adopt_grid=adopt_grid,
            )
        )
        self.channels = self._source.channels
        # The capability flags (DEVELOPMENT_PLAN.md §4.5) are the source's: the
        # background has no free parameters, so it cannot be the reason a
        # composed problem is or is not differentiable. Instance attributes,
        # which is what a realisation reads.
        self.DIFFERENTIABLE = bool(self._source.DIFFERENTIABLE)
        self.BATCHABLE = bool(self._source.BATCHABLE)
        self.DEVICE = self._source.DEVICE
        self.BACKEND = self._source.BACKEND
        for parameter in self._source.parameters:
            self.register_parameter(parameter)

    # -- negotiation: both children see the same requirements ---------------

    def compile_for(self, requirements: Any) -> SourceWithBackground:
        """Negotiate both children identically, so their grids cannot diverge."""
        self._source = self._source.compile_for(requirements)
        if self._background is not None:
            self._background = self._background.compile_for(requirements)
        return self

    # -- the contract surface -------------------------------------------------

    def evaluate(self, **values: Any) -> ModelResult:
        source = self._source.evaluate(**values)
        if self._background is None:
            return source
        background = self._background.evaluate()
        combined = {}
        for channel in self.channels:
            image = source[channel]
            combined[channel] = image.with_values(image.values + background[channel].values)
        return ModelResult(combined)

    # -- the native surface (torch, jax) --------------------------------------

    def native_grid(self, channel: str) -> Any:
        """The negotiated ``(x, y)``, identical for both children by construction."""
        return self._source.native_grid(channel)

    def native_flux(self, channel: str, values: Any = None) -> Any:
        """The source's native brightness, plus the background's if there is one."""
        flux = self._source.native_flux(channel, values)
        if self._background is None:
            return flux
        return flux + self._background.native_flux(channel, None)


def source_with_background(itf: Any, x: Any, y: Any, **kwargs: Any) -> SourceWithBackground:
    """Convenience constructor: see :class:`SourceWithBackground`."""
    return SourceWithBackground(itf, x, y, **kwargs)
