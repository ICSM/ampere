"""The toy: a linear continuum times two Gaussian absorption lines.

Four parameters, uniform priors, an analytic forward model, and a spectral
range chosen to look like a Gaia RVS observation of the calcium triplet. It is
the paper study's model
(``examples/examples_paper/flexible_likelihood_comparison.py``), transcribed:

.. math::

    F(\\lambda) = \\bigl[A + B(\\lambda - \\lambda_0)\\bigr]
                  \\left[1 - d_1 g_1(\\lambda) - d_2 g_2(\\lambda)\\right],

with :math:`g_k` unit-height Gaussians at fixed centres and widths, and
:math:`\\lambda_0 = 0.857\\,\\mu\\mathrm{m}`. The line positions and widths are
**buffers**, not parameters: the study is about the four amplitudes, and giving
the sampler line centres to move as well would let it chase a misspecification
that the experiment is trying to make it *unable* to chase.

Why the model is written here and not in ``ampere.backends``
-------------------------------------------------------------
W2.10's scoping left this open — "if the toy model needs to ship in the
backends instead (because the realised paths only lower shipped models), say so
and put it there". It does not. Both differentiable backends' lowerings
(:mod:`ampere.backends.torch.problem`, :mod:`ampere.backends.jax.problem`)
resolve a model's native surface by **duck typing**: they walk
``model.native_grid(channel)`` / ``model.native_flux(channel, values)``, the
canonical spelling since *W5.20*, or the legacy ``grid``/``flux`` pair these
three variants deliberately keep — a model written before the ruling still
lowers unchanged, and this study is where that is exercised end to end. A
hand-written model that offers either pair lowers into a differentiable
log-density exactly as
``PowerLaw`` does, with no registration and no subclassing of a backend class.
That is a real property of the architecture — a user's own model gets NUTS —
and demonstrating it in an example is worth more than shipping a fifth model
nobody asked for. The three variants live in :mod:`.model` (numpy),
:mod:`.model_torch` and :mod:`.model_jax`.

The three are held to each other by ``tests/m2/test_model.py``, which compares
their flux arrays point by point at the same parameters; nothing but the array
library differs between them.

Units
-----
Micron on the spectral axis, Jy on the flux, converted nowhere in the hot loop
(``DEVELOPMENT_PLAN.md`` §7). The flux is a dimensionless-looking number of
order unity; ``Jy`` is the unit ampere's containers and the GP amplitude prior
carry it in.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.core import (
    DTYPE,
    ChannelRequirements,
    Model,
    ModelResult,
    Parameter,
    Spectrum,
)

__all__ = [
    "FLUX_UNIT",
    "LINE1_CENTRE",
    "LINE1_WIDTH",
    "LINE2_CENTRE",
    "LINE2_WIDTH",
    "PARAMETER_NAMES",
    "PRIOR_LIMITS",
    "REFERENCE_WAVELENGTH",
    "TRUTH",
    "WAVELENGTH_UNIT",
    "AbsorptionLines",
    "default_priors",
    "flux_at",
]

#: The spectral coordinate unit.
WAVELENGTH_UNIT = u.micron
#: The flux unit.
FLUX_UNIT = u.Jy

#: The two absorption lines: centre and Gaussian sigma, micron. Fixed, and
#: registered as buffers rather than parameters — see the module docstring.
LINE1_CENTRE = 0.8498
LINE1_WIDTH = 0.00035
LINE2_CENTRE = 0.8542
LINE2_WIDTH = 0.00030

#: Where the continuum's normalisation ``A`` is defined, micron. A buffer, for
#: the same reason ``ampere.backends.reference.PowerLaw``'s reference
#: wavelength is one: it pins what ``A`` means, and a prior on it would be a
#: prior on one's own choice of origin.
REFERENCE_WAVELENGTH = 0.857

#: The four fitted parameters, in the order every table and figure uses.
PARAMETER_NAMES: tuple[str, ...] = ("A", "B", "d1", "d2")

#: The truth the data are generated at.
TRUTH: dict[str, float] = {"A": 1.0, "B": 3.0, "d1": 0.15, "d2": 0.10}

#: Uniform prior support, ``(lower, upper)`` per parameter — the paper study's
#: ``lims``. Wide enough that the prior is not doing the inference's work, and
#: bounded, so that ``d1``/``d2`` cannot go negative and turn an absorption
#: line into an emission one the deviation would then be confused with.
PRIOR_LIMITS: dict[str, tuple[float, float]] = {
    "A": (0.5, 1.5),
    "B": (-5.0, 12.0),
    "d1": (0.0, 0.5),
    "d2": (0.0, 0.5),
}


def default_priors() -> dict[str, Any]:
    """Frozen ``scipy.stats`` uniform priors over :data:`PRIOR_LIMITS`.

    Built fresh on each call rather than held as a module constant: a frozen
    distribution is shared mutable-ish state as far as reproducibility
    arguments go, and three backends constructing three models from one object
    would make "the same prior" a claim about identity rather than about
    numbers.

    The bijection is not stated here. ``ampere.core.default_bijection_for``
    reads it off the prior's own bounded support and returns a
    :class:`~ampere.core.Logit`, which is what makes these four parameters
    safe for NUTS: no proposal can leave the support, so no draw is rejected
    for a reason a gradient sampler cannot see.
    """
    return {name: st.uniform(lower, upper - lower) for name, (lower, upper) in PRIOR_LIMITS.items()}


def flux_at(wavelength: Any, *, A: float, B: float, d1: float, d2: float) -> np.ndarray:
    """The model flux in numpy, as a free function.

    The generators need the truth before any :class:`~ampere.core.Model` exists
    (they build the data a model is later fitted to), and the cross-backend test
    needs one implementation to hold the other two to. Both want the arithmetic
    without the contract around it, so it lives here and
    :class:`AbsorptionLines` calls it.
    """
    grid = np.asarray(wavelength, dtype=DTYPE)
    continuum = A + B * (grid - REFERENCE_WAVELENGTH)
    first = np.exp(-0.5 * ((grid - LINE1_CENTRE) / LINE1_WIDTH) ** 2)
    second = np.exp(-0.5 * ((grid - LINE2_CENTRE) / LINE2_WIDTH) ** 2)
    return continuum * (1.0 - d1 * first - d2 * second)


class AbsorptionLines(Model):
    """The toy model on the reference (numpy) path.

    Follows :class:`ampere.backends.reference.models.PowerLaw` in every
    structural respect — the grid is a buffer, each physical quantity is
    accepted as a prior, a number or a ready-made
    :class:`~ampere.core.Parameter`, and :meth:`compile_for` adopts the
    negotiated grid once and refills one container per evaluation — because
    those are the contract's requirements rather than that class's habits.

    Parameters
    ----------
    wavelength
        The model's own grid, micron (bare or a
        :class:`~astropy.units.Quantity`).
    A, B, d1, d2
        Continuum normalisation at :data:`REFERENCE_WAVELENGTH`, continuum
        slope (Jy per micron), and the two lines' fractional depths. A prior to
        fit, a number to hold fixed.
    channels
        Name of the channel the emitted spectrum appears under.

    Examples
    --------
    >>> import numpy as np
    >>> grid = np.linspace(0.842, 0.872, 5)
    >>> model = AbsorptionLines(grid)
    >>> model.parameters.free_names
    ('A', 'B', 'd1', 'd2')
    >>> result = model.evaluate(A=1.0, B=3.0, d1=0.15, d2=0.10)
    >>> bool(np.allclose(result["default"].values, flux_at(grid, **TRUTH)))
    True
    """

    #: The spectral-axis requirement this model answers.
    AXIS: ClassVar[str] = "spectral_axis"
    #: The four capability flags (``inference.md`` §18, W2.12's fourth),
    #: declared explicitly. numpy: no gradients, no batching, CPU, reference.
    DIFFERENTIABLE: ClassVar[bool] = False
    BATCHABLE: ClassVar[bool] = False
    DEVICE: ClassVar[str] = "cpu"
    BACKEND: ClassVar[str] = "reference"

    def __init__(
        self,
        wavelength: Any,
        *,
        A: Any = None,
        B: Any = None,
        d1: Any = None,
        d2: Any = None,
        channels: str | Sequence[str] = "default",
    ) -> None:
        grid = _to_micron(wavelength)
        self.channels = _channels(channels)
        self.register_buffer("wavelength", grid, unit=WAVELENGTH_UNIT)
        for name, given in _declared(A=A, B=B, d1=d1, d2=d2).items():
            self.register_parameter(_as_parameter(name, given))
        self.templates: dict[str, Spectrum] = {}
        self.grids: dict[str, np.ndarray] = {channel: grid for channel in self.channels}

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated grid for each channel, once."""
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None or self.AXIS not in asked:
                continue
            grid = _to_micron(asked[self.AXIS].coordinates())
            self.templates[channel] = Spectrum(
                grid * WAVELENGTH_UNIT, np.zeros(grid.size, dtype=DTYPE), unit=FLUX_UNIT
            )
            self.grids[channel] = grid
        return self

    def grid(self, channel: str) -> np.ndarray:
        """The coordinates *channel* is evaluated on, micron."""
        return self.grids[channel]

    def flux(self, channel: str, values: Mapping[str, Any] | None = None) -> np.ndarray:
        """This model's flux on *channel*, Jy.

        Spelled the same as the two differentiable variants' native surface, so
        the three read as one model in three libraries.
        """
        context = self.context(values)
        return flux_at(
            self.grid(channel),
            A=float(context["A"]),
            B=float(context["B"]),
            d1=float(context["d1"]),
            d2=float(context["d2"]),
        )

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        emitted = {}
        for channel in self.channels:
            flux = flux_at(
                self.grid(channel),
                A=float(context["A"]),
                B=float(context["B"]),
                d1=float(context["d1"]),
                d2=float(context["d2"]),
            )
            template = self.templates.get(channel)
            emitted[channel] = (
                Spectrum(self.grid(channel) * WAVELENGTH_UNIT, flux, unit=FLUX_UNIT)
                if template is None
                else template.with_values(flux)
            )
        return ModelResult(emitted)


# ---------------------------------------------------------------------------
# Shared construction helpers. The torch and jax variants import these, so the
# three declarations cannot drift in their priors, their channel handling or
# their unit conversion — only in their arithmetic.
# ---------------------------------------------------------------------------


def _declared(**given: Any) -> dict[str, Any]:
    """Fill in :func:`default_priors` for whichever of the four were left out."""
    priors = default_priors()
    return {name: priors[name] if value is None else value for name, value in given.items()}


def _as_parameter(name: str, given: Any) -> Parameter:
    """A prior, a number or a ready-made :class:`~ampere.core.Parameter`."""
    if isinstance(given, Parameter):
        if given.name != name:
            raise ValueError(
                f"parameter {name!r} was given a Parameter named {given.name!r}; a model's "
                f"parameter names are part of its interface."
            )
        return given
    if isinstance(given, (int, float, np.floating, np.integer)) and not isinstance(given, bool):
        return Parameter(name, value=float(given), fixed=True)
    if hasattr(given, "ppf"):
        return Parameter(name, given)
    raise TypeError(
        f"parameter {name!r} must be a frozen scipy.stats distribution, a number, or an ampere "
        f"Parameter — got {type(given).__name__}."
    )


def _channels(channels: str | Sequence[str]) -> tuple[str, ...]:
    names = (channels,) if isinstance(channels, str) else tuple(str(c) for c in channels)
    if not names or len(set(names)) != len(names):
        raise ValueError(f"a spectral model needs distinct, non-empty channel names, got {names!r}")
    return names


def _to_micron(coordinates: Any) -> np.ndarray:
    """Coordinates as bare micron, whether or not they arrive as a Quantity."""
    if isinstance(coordinates, u.Quantity):
        grid = np.asarray(coordinates.to_value(WAVELENGTH_UNIT), dtype=DTYPE)
    else:
        grid = np.asarray(coordinates, dtype=DTYPE)
    if grid.ndim != 1 or grid.size == 0:
        raise ValueError(
            f"a spectral model needs a one-dimensional, non-empty wavelength grid, got shape "
            f"{grid.shape}."
        )
    return grid
