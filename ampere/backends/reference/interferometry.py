"""Interferometric visibilities and closure phases on the numpy path.

Phase 4's proof modality, as shipped code: the steps and models that turn a
model's sky brightness into what an interferometer records. ``ampere.core``
owns the two container *kinds* (:class:`~ampere.core.VisibilitySet` and
:class:`~ampere.core.ClosurePhases`, beside the kinds the freeze shipped) and
this module owns the arithmetic, which is the placement D1 ruled on
2026-09-11: a shipped observable's kind lives in core, its steps live in a
per-backend ``<observable>.py``, and there is no grouping namespace.

The composition this module exists to support, in one picture::

    Binary (Model)
      └── channel "sky" : Image           (x, y in mas; surface brightness)
             │
             ├── Instrument "gravity_vis" [FourierSample, BandwidthSmearing]
             │        → VisibilitySet, Likelihood(ComplexGaussianFamily(), IndependentNoise())
             │
             └── Instrument "gravity_t3"  [FourierSample, ClosurePhase]
                      → ClosurePhases, Likelihood(VonMisesFamily(), IndependentNoise())

Two instruments bind the **same** channel by name, publish requirements on the
same two axes, and :func:`~ampere.core.negotiate` unions them into one image
grid the model builds **once** per draw. That is the pattern
``DEVELOPMENT_PLAN.md`` §4.3 was designed for, and the pairing between the two
datasets stays visible in the composition — both name ``sky``, and
``problem.requirements["model"]["sky"].sources`` names both instruments — so
that the plan's design horizon (h), a joint noise model over a tuple of
channels, has something to bind to without a re-plumb. Nothing here folds the
two observables into one container.

Four conventions are fixed here, and each is stated where it is used as well
as here, because every one of them is a sign or a factor that would look like
real source structure if it were wrong.

**The transform.** :class:`FourierSample` computes

.. math::

    V(u, v) = \\int\\!\\!\\int I(x, y)\\,
              e^{-2\\pi i (u x + v y)}\\, \\mathrm{d}x\\,\\mathrm{d}y

with ``x``, ``y`` angular offsets **in radians** (the containers carry mas and
are converted once, at composition time) and ``u``, ``v`` baselines **in
wavelengths**, so the exponent is dimensionless. The sign is the
``exp(-2 pi i)`` convention: a component displaced to positive ``x`` acquires a
*negative* phase. The discrete form is the quadrature
``sum_jk I_jk dOmega_jk exp(-2 pi i (u x_j + v y_k))``, with ``dOmega_jk`` the
solid angle of the ``(j, k)`` cell in steradians, taken from the midpoints of
the container's own axes — so an image in Jy/sr gives visibilities in Jy, and
a grid that negotiation made uneven is integrated correctly rather than
assumed regular.

**The pixel scale.** The step publishes ``max_step = 1/(2 u_max)`` on each
image axis, which is where gap I-4 bites and why
:meth:`~ampere.core.Model.compile_for` raising matters here more than
anywhere else: a model that silently ignores it folds power from beyond the
Nyquist limit back onto the sampled baselines, where it is indistinguishable
from real structure. That is not a bad fit that announces itself.

**The closure phase.** :class:`ClosurePhase` produces
``arg(V(b1) · V(b2) · conj(V(b1 + b2)))`` for a triangle whose stored
baselines are ``b1 = ij`` and ``b2 = jk`` with telescopes ``i < j < k`` — the
canonical ordering :class:`~ampere.core.ClosurePhases` fixes. With the
opposite baseline sign convention every closure phase changes sign and a fit
to the mirrored source looks exactly as good, so the convention is asserted
against a closed form rather than against this module's own arithmetic.

**Coordinates come from the container, never from arithmetic.**
``transformations.md`` §10 (W1.11 gap I-2), and this is the modality it was
written for: ``check_alignment`` compares axes with ``numpy.array_equal``, and
a ``(u, v)`` recomputed from station coordinates, hour angle and wavelength
differs in its last bits from the one in the file. :meth:`FourierSample.from_observed`
is how the rule is kept — the coverage is a **buffer** taken from the observed
container, never recomputed, and never a parameter (one would not put a prior
on where the telescopes were).

What this module does **not** do, deliberately:

* **No FFT.** ``AxisRequirement.coordinates()`` happens to return an evenly
  spaced grid for a ``max_step`` requirement, but regularity is advertised and
  never required (``architecture.md`` §7), and ``compile_for`` may hand back
  any grid satisfying the requirement. A step that assumed an FFT would be
  incorrect, not merely slow. The gridded FFT-plus-interpolation variant is a
  later option and is refused by name rather than half-written: correctness is
  the reference backend's only goal (``architecture.md`` §2).
* **No chromatic sky.** The ``Image`` channel is achromatic: the step reads
  each sample's own wavelength (and the smearing steps use it), but the model
  produces one image for every wavelength in the dataset. A wavelength-
  dependent sky is a ``Cube`` channel and a Fourier step that publishes on
  three axes; the requirement language already covers it, and it is a later
  item.
* **No reader.** Every container here is built from arrays. OIFITS is Phase
  6's.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.special

from ampere.core import (
    DTYPE,
    AxisRequirement,
    ChannelRequirements,
    ClosurePhases,
    CompositionError,
    Image,
    Model,
    ModelResult,
    Transformation,
    TransformationError,
    VisibilitySet,
    propagate_mask,
)

from ._declare import as_parameter

__all__ = [
    "BRIGHTNESS_UNIT",
    "COORDINATE_UNIT",
    "FLUX_UNIT",
    "MAS_PER_RAD",
    "SPECTRAL_UNIT",
    "Amplitude",
    "BandwidthSmearing",
    "Binary",
    "BinaryVisibilities",
    "ClosurePhase",
    "FourierSample",
    "GaussianSource",
    "GaussianSourceVisibilities",
    "TimeSmearing",
    "UniformDisc",
    "UniformDiscVisibilities",
    "cell_solid_angle",
]

#: The angular unit every image axis and every angular size here works in.
COORDINATE_UNIT = u.mas
#: The spectral unit the wavelength axis works in.
SPECTRAL_UNIT = u.micron
#: The unit a visibility (a correlated flux) is emitted in.
FLUX_UNIT = u.Jy
#: The unit an image is emitted in: a surface brightness, so that the Fourier
#: quadrature's solid-angle factor lands on Jy.
BRIGHTNESS_UNIT = u.Jy / u.sr

#: Milliarcseconds in one radian. Converted once, here, rather than per
#: evaluation (``DEVELOPMENT_PLAN.md`` §7's units trap).
MAS_PER_RAD = float((1.0 * u.rad).to_value(u.mas))

#: Gaussian FWHM in units of its standard deviation.
_FWHM_PER_SIGMA = float(2.0 * np.sqrt(2.0 * np.log(2.0)))

#: How many quadrature nodes a smearing step averages over by default. Odd, so
#: that the centre node is the un-smeared sample itself (see
#: :func:`_quadrature`), and Gauss-Legendre, so that five nodes integrate a
#: smooth visibility across a channel to far better than the measurement.
_DEFAULT_SMEARING_NODES = 5


def _to_unit(values: Any, unit: u.UnitBase) -> np.ndarray:
    """*values* as a bare float64 array in *unit*, whether or not they arrive as a Quantity."""
    if isinstance(values, u.Quantity):
        return np.asarray(values.to_value(unit), dtype=DTYPE)
    return np.asarray(values, dtype=DTYPE)


def _quadrature(nodes: int) -> tuple[np.ndarray, np.ndarray]:
    """Gauss-Legendre nodes on ``[-1/2, 1/2]`` and their weights, **centre first**.

    Two properties are load-bearing rather than convenient.

    *Odd* order, so a node sits exactly at the centre of the interval; and the
    centre node is moved to **index 0**, so that the un-smeared sample is the
    first sub-sample of every block. That is what lets a smearing step recover
    its output coordinates by slicing its input — ``[:, 0, :]`` — rather than
    recomputing them, which would break the bit-identical coordinate rule
    (``transformations.md`` §10, gap I-2) on the very container a likelihood
    compares against the data.

    The weights are normalised to sum to one, because these steps compute a
    weighted **mean** visibility over the smearing interval, not an integral.
    """
    if nodes < 1 or nodes % 2 == 0:
        raise TransformationError(
            f"a smearing step averages over an odd number of quadrature nodes, so that the "
            f"un-smeared sample is one of them and the step can take its output coordinates "
            f"from its input rather than recomputing them; got {nodes}."
        )
    if nodes == 1:
        return np.zeros(1, dtype=DTYPE), np.ones(1, dtype=DTYPE)
    raw_nodes, raw_weights = np.polynomial.legendre.leggauss(nodes)
    # Symmetrise explicitly: the centre node must be *exactly* zero, not zero
    # to rounding, or the un-smeared coordinates come back perturbed in their
    # last bits and check_alignment fails with a message about axes.
    middle = nodes // 2
    offsets = 0.5 * np.asarray(raw_nodes, dtype=DTYPE)
    offsets[middle] = 0.0
    weights = np.asarray(raw_weights, dtype=DTYPE) / float(np.sum(raw_weights))
    order = np.concatenate([[middle], np.delete(np.arange(nodes), middle)])
    return offsets[order], weights[order]


def cell_solid_angle(x_mas: np.ndarray, y_mas: np.ndarray) -> np.ndarray:
    """Solid angle of every cell of an ``(x, y)`` image grid, in steradians.

    The quadrature weight of the Fourier sum. Cell widths are the distances
    between the midpoints of neighbouring coordinates, with the outermost two
    reflected — the usual convention for samples given as centres — taken in
    absolute value, because an ``Image`` axis is monotonic in *either*
    direction (sky axes legitimately run both ways) and a negative solid angle
    would flip the sign of the whole transform.

    Uneven spacing is handled rather than assumed away: negotiation may return
    any grid satisfying the published requirement, and the union of two
    instruments' requirements is very often not evenly spaced.
    """
    widths = [_widths(np.asarray(axis, dtype=DTYPE)) for axis in (x_mas, y_mas)]
    return np.abs(np.outer(widths[0], widths[1])) / MAS_PER_RAD**2


def _widths(centres: np.ndarray) -> np.ndarray:
    """Cell widths from sample centres: midpoints, outer two reflected."""
    if centres.size == 1:
        return np.ones(1, dtype=DTYPE)
    middle = 0.5 * (centres[1:] + centres[:-1])
    first = centres[0] - (middle[0] - centres[0])
    last = centres[-1] + (centres[-1] - middle[-1])
    edges = np.concatenate([[first], middle, [last]])
    return np.diff(edges)


# ---------------------------------------------------------------------------
# The steps
# ---------------------------------------------------------------------------


class _Step(Transformation):
    """Shared plumbing: read a declared buffer as a plain array.

    The same shape as ``instrument.py``'s private base, and separate from it on
    purpose: a buffer's name may not shadow a class attribute, so a step whose
    buffer is ``u_pts`` cannot also expose a ``u_pts`` property, and this is
    the readable spelling of ``self.buffers[name].value``.
    """

    #: The fourth capability flag (W2.12), declared rather than inherited:
    #: every piece of this backend says so for itself.
    BACKEND: ClassVar[str] = "reference"

    def _data(self, name: str) -> np.ndarray:
        return np.asarray(self.buffers[name].value, dtype=DTYPE)


@dataclasses.dataclass(frozen=True)
class _Expansion:
    """Extra ``(u, v, lambda)`` samples one averaging step needs, per output sample.

    The payload of gap I-3's mechanism. A step that averages the visibility
    over a small spread of ``(u, v)`` — bandwidth smearing along the spoke,
    time smearing along the Earth-rotation track — cannot publish an
    :class:`~ampere.core.AxisRequirement` for it, because requirements are
    statements about the *channel's* coordinates and the channel is an
    ``Image``. What it can do is tell the step before it, which is what
    :meth:`~ampere.core.Transformation.configure_from` is for.

    ``u``, ``v`` and ``wavelength`` each gain one trailing axis of length
    :attr:`nodes` relative to what they were given, and index ``0`` along that
    axis is the input sample unchanged.
    """

    u: np.ndarray
    v: np.ndarray
    wavelength: np.ndarray
    weights: np.ndarray

    @property
    def nodes(self) -> int:
        """How many sub-samples this expansion adds per input sample."""
        return int(self.weights.size)


class _AveragingStep(_Step):
    """Common half of the two smearing steps: the block reduction and the bookkeeping.

    Both steps do the same three things and differ only in *where* the extra
    samples sit: expand the ``(u, v, lambda)`` set the Fourier step evaluates
    (:meth:`expand_uv`, read by :class:`FourierSample`'s ``configure_from``),
    learn how many sub-samples the steps *after* this one added
    (:meth:`configure_from`, so that the block this one reduces can be located
    in a flat array), and average each block at evaluation time.

    The flat layout is fixed by :meth:`FourierSample.configure_from` and read
    back here: with expansions ``k1, k2, ...`` applied in chain order, sample
    ``i``'s sub-sample ``(a, b, ...)`` lives at index ``((i k1 + a) k2 + b)...``
    So the first averaging step reduces the *first* sub-axis, and what remains
    is exactly what the next one expects.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (VisibilitySet,)

    def __init__(self, *, nodes: int = _DEFAULT_SMEARING_NODES, label: str | None = None) -> None:
        super().__init__(label=label)
        offsets, weights = _quadrature(int(nodes))
        self._offsets = offsets
        self._weights = weights
        self._trailing = 1

    @property
    def nodes(self) -> int:
        """Quadrature nodes this step averages over."""
        return int(self._weights.size)

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """Count the sub-samples the steps after this one added.

        Needed to find this step's own block in a flat array: the total
        expansion is a product, and this step reduces one factor of it. Read
        from the successors' own declarations, never inferred — the same
        posture ``configure_from``'s contract states.
        """
        trailing = 1
        for step in downstream:
            expander = getattr(step, "expand_uv", None)
            if expander is not None:
                trailing *= int(getattr(step, "nodes", 1))
        self._trailing = trailing

    def expand_uv(self, u_pts: np.ndarray, v_pts: np.ndarray, wavelength: np.ndarray) -> _Expansion:
        """The extra samples this step averages over. Implemented by each step."""
        raise NotImplementedError

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        weights = self._weights
        count = weights.size
        trailing = self._trailing
        total = int(samples.n_samples)
        block = count * trailing
        if total % block:
            raise TransformationError(
                f"{type(self).__name__} was given {total} visibilities, which is not a multiple "
                f"of the {count} sub-sample(s) it averages over times the {trailing} the steps "
                f"after it added. That happens when the chain's Fourier step was not the one "
                f"configured for this chain: build the step and the instrument together, so "
                f"configure_from runs on the pair ampere will evaluate."
            )
        n_out = total // count
        shaped = np.asarray(samples.values).reshape(-1, count, trailing)
        averaged = np.einsum("a,nar->nr", weights, shaped).reshape(n_out)
        kept = (slice(None), 0, slice(None))
        coordinates = {
            name: samples.axis(name).values.reshape(-1, count, trailing)[kept].reshape(n_out)
            for name in ("u", "v", "spectral_axis")
        }
        mask = None
        if samples.mask is not None:
            mask = propagate_mask(samples, self._influence(n_out, count, trailing))
        return VisibilitySet(
            coordinates["u"],
            coordinates["v"],
            coordinates["spectral_axis"] * SPECTRAL_UNIT,
            averaged,
            unit=samples.unit,
            mask=mask,
        )

    @staticmethod
    def _influence(n_out: int, count: int, trailing: int) -> np.ndarray:
        """The ``(n_out, n_out * count)`` averaging matrix, for the mask only.

        Built only when there is a mask to carry, because the values are
        reduced by an ``einsum`` that needs no matrix and this one is
        quadratic in the sample count. ``propagate_mask`` is the rule
        (``results_schema.md`` §16): an output sample is masked if **any**
        sub-sample that influences it is masked, which for a smeared
        visibility is the right answer — a flagged channel edge contaminates
        the average it is part of.
        """
        rows = np.arange(n_out).reshape(-1, 1)
        outer, inner = np.divmod(rows, trailing)
        columns = (outer * count + np.arange(count).reshape(1, -1)) * trailing + inner
        influence = np.zeros((n_out, n_out * count), dtype=DTYPE)
        np.put_along_axis(influence, columns, 1.0, axis=1)
        return influence


class BandwidthSmearing(_AveragingStep):
    """Average the visibility across a spectral channel of finite width.

    A correlator channel is not monochromatic. The baseline **B** is fixed in
    metres, so the sampled spatial frequency ``B/lambda`` sweeps along its own
    spoke as ``lambda`` crosses the channel, and what is recorded is the mean
    visibility over that sweep. The effect is a radial smearing that reduces
    the amplitude of a resolved source at long baselines — it looks exactly
    like a larger source if it is not modelled.

    This is one of the two steps gap I-3 was written about, and the first
    *cross-kind* use of ``configure_from``: it needs more ``(u, v)`` samples,
    those samples are the Fourier step's own buffer rather than anything the
    model chooses, and a requirement on ``u`` would be refused because the
    channel holds an ``Image``. (``LSFConvolution`` already uses
    ``configure_from`` on all three backends, but within one kind.) So this
    step tells :class:`FourierSample` what to compute, and
    :class:`FourierSample` computes it.

    Parameters
    ----------
    resolving_power
        ``lambda / dlambda`` of the channel. The sweep is
        ``lambda_0 (1 +- 1/(2 R))``, so a larger ``R`` is a narrower channel
        and less smearing. A **buffer**, not a parameter: the channel width is
        the correlator's, known from the instrument.
    nodes
        Quadrature nodes, odd (see :func:`_quadrature`). Five Gauss-Legendre
        nodes integrate a visibility across an ordinary channel to far better
        than any measurement of it; raise it for a very wide channel across a
        heavily resolved source.
    label
        Component label for this step within a chain.
    """

    def __init__(
        self,
        *,
        resolving_power: float,
        nodes: int = _DEFAULT_SMEARING_NODES,
        label: str | None = None,
    ) -> None:
        super().__init__(nodes=nodes, label=label)
        power = float(resolving_power)
        if not np.isfinite(power) or power <= 0.0:
            raise TransformationError(
                f"BandwidthSmearing's resolving_power is lambda/dlambda of the channel and must "
                f"be finite and positive, got {resolving_power!r}."
            )
        self.register_buffer("resolving_power", power)

    def expand_uv(self, u_pts: np.ndarray, v_pts: np.ndarray, wavelength: np.ndarray) -> _Expansion:
        """Sub-samples along the spoke: ``lambda -> lambda(1 + s/R)``, ``u -> u lambda_0/lambda``.

        The spatial frequency is ``B/lambda`` with **B** fixed, so scaling the
        wavelength scales ``(u, v)`` by the reciprocal — the sub-samples lie on
        the *radial* line through the origin and the recorded point, which is
        why bandwidth smearing is radial and time smearing is not.
        """
        power = float(self.buffers["resolving_power"].value)
        factor = 1.0 + self._offsets / power
        scaled = 1.0 / factor
        return _Expansion(
            u=u_pts[..., None] * scaled,
            v=v_pts[..., None] * scaled,
            wavelength=wavelength[..., None] * factor,
            weights=self._weights,
        )


class TimeSmearing(_AveragingStep):
    """Average the visibility over a finite integration, along the uv track.

    While the correlator integrates, the Earth turns and the projected
    baseline moves: the recorded sample is the mean visibility along a short
    arc of the uv track rather than a value at a point. Unlike bandwidth
    smearing the direction is not radial — it is wherever the track goes —
    so the two effects are genuinely different steps rather than one with a
    parameter.

    **The uv rate must be supplied; it cannot be derived.** A
    ``VisibilitySet`` records where a sample was taken, not when or how fast
    the point was moving, and reconstructing ``du/dt`` from station
    coordinates and hour angle is precisely the recomputation gap I-2 forbids
    for the coordinates themselves. So the rates are per-sample **buffers**,
    given at construction or read from the observed container's
    ``extra_coords`` by :meth:`from_observed`, and a caller who has neither is
    refused by name rather than served a plausible guess.

    Parameters
    ----------
    integration
        Integration time per sample, seconds (or a
        :class:`~astropy.units.Quantity` of time).
    du_dt, dv_dt
        Rate of change of ``(u, v)`` along the track, in wavelengths per
        second, one per observed sample.
    nodes
        Quadrature nodes, odd (see :func:`_quadrature`).
    label
        Component label for this step within a chain.
    """

    def __init__(
        self,
        *,
        integration: Any,
        du_dt: Any,
        dv_dt: Any,
        nodes: int = _DEFAULT_SMEARING_NODES,
        label: str | None = None,
    ) -> None:
        super().__init__(nodes=nodes, label=label)
        seconds = float(_to_unit(integration, u.s))
        if not np.isfinite(seconds) or seconds <= 0.0:
            raise TransformationError(
                f"TimeSmearing's integration is the time per sample and must be finite and "
                f"positive, got {integration!r}."
            )
        rates = [np.asarray(rate, dtype=DTYPE).reshape(-1) for rate in (du_dt, dv_dt)]
        if rates[0].size != rates[1].size or rates[0].size == 0:
            raise TransformationError(
                f"TimeSmearing needs one du/dt and one dv/dt per observed sample, in wavelengths "
                f"per second; got {rates[0].size} and {rates[1].size}."
            )
        self.register_buffer("integration", seconds, unit=u.s)
        self.register_buffer("du_dt", rates[0])
        self.register_buffer("dv_dt", rates[1])

    @classmethod
    def from_observed(
        cls,
        container: VisibilitySet,
        *,
        integration: Any,
        nodes: int = _DEFAULT_SMEARING_NODES,
        label: str | None = None,
    ) -> TimeSmearing:
        """Take the uv rates from the observed container's ``extra_coords``.

        The rates belong with the data they describe, and ``extra_coords`` is
        where per-sample annotations that do not define the sampling geometry
        live (``results_schema.md`` §15.4). A container without them is
        refused with the two key names it needs, rather than with a
        ``KeyError`` from inside a step.
        """
        missing = [name for name in ("du_dt", "dv_dt") if name not in container.extra_coords]
        if missing:
            raise TransformationError(
                f"TimeSmearing.from_observed needs the uv track rates on the observed "
                f"{type(container).__name__}'s extra_coords, and {missing} are not there (it has "
                f"{sorted(container.extra_coords)}). A visibility records where a sample was "
                f"taken, not how fast the baseline was moving, and ampere will not reconstruct "
                f"that from station coordinates — pass extra_coords={{'du_dt': ..., 'dv_dt': ...}} "
                f"in wavelengths per second when you build the container, or give the rates to "
                f"TimeSmearing(...) directly."
            )
        return cls(
            integration=integration,
            du_dt=container.extra_coords["du_dt"],
            dv_dt=container.extra_coords["dv_dt"],
            nodes=nodes,
            label=label,
        )

    def expand_uv(self, u_pts: np.ndarray, v_pts: np.ndarray, wavelength: np.ndarray) -> _Expansion:
        """Sub-samples along the track: ``u -> u + (du/dt) t``, ``t`` across the integration."""
        seconds = float(self.buffers["integration"].value)
        times = self._offsets * seconds
        du_dt = self._data("du_dt")
        dv_dt = self._data("dv_dt")
        shape = u_pts.shape
        if du_dt.size != shape[0]:
            raise TransformationError(
                f"TimeSmearing was built with {du_dt.size} uv rate(s) but the chain's Fourier "
                f"step covers {shape[0]} observed sample(s). The rates are per observed sample, "
                f"in the container's own order."
            )
        broadcast: tuple[int, ...] = (shape[0],) + (1,) * (len(shape) - 1)
        return _Expansion(
            u=u_pts[..., None] + du_dt.reshape((*broadcast, 1)) * times,
            v=v_pts[..., None] + dv_dt.reshape((*broadcast, 1)) * times,
            wavelength=np.broadcast_to(wavelength[..., None], (*shape, times.size)).copy(),
            weights=self._weights,
        )


class FourierSample(_Step):
    """Sample a model image at the observed ``(u, v)`` points: ``Image -> VisibilitySet``.

    The kind-changing, coordinate-changing step in the middle of the chain
    that makes this modality a genuine test of the composition contracts. A
    direct discrete Fourier transform: no FFT, no gridding, no interpolation
    (see the module docstring for why).

    ``(u, v, lambda)`` are **buffers taken from the observed container**
    (:meth:`from_observed`), never recomputed and never parameters:
    ``architecture.md`` §6's question — "would you ever put a prior on it?" —
    answers itself, and ``check_alignment`` compares axes exactly, so a
    recomputed coordinate fails a check whose message is about axes rather
    than about arithmetic (gap I-2).

    The step carries no parameters at all. A calibration factor, a coherence
    loss or a fitted phase offset is a separate one-step transformation later
    in the chain, which is how ``transformations.md`` §5 says instrument-level
    nuisance parameters are expressed.

    What it publishes, and why it is exact rather than conservative
    (``interferometry.md`` §4): the DFT kernel factorises as
    ``exp(-2 pi i u x) exp(-2 pi i v y)``, so the sampling constraint is
    separable and a requirement language that speaks about one axis at a time
    is the *natural* shape of the constraint here, not a limitation. The field
    of view fixes ``intervals`` on both axes; the longest baseline fixes
    ``max_step`` on each.

    Parameters
    ----------
    u_pts, v_pts
        Baseline coordinates in wavelengths, one per observed sample.
    wavelength
        The wavelength of each sample (micron, or a
        :class:`~astropy.units.Quantity`). Carried through to the emitted
        container's third axis and used by the smearing steps; the image
        itself is achromatic here.
    field_of_view
        Angular extent of the region the instrument is sensitive to — the
        fibre or primary beam — in mas. The image must cover ``+-fov/2`` on
        both axes.
    oversampling
        Multiplier on the Nyquist sampling rate the published ``max_step``
        asks for. ``1`` (the default) is the array's own limit; raise it when
        the model has structure finer than the beam, which a sum at the
        Nyquist step would alias. See :meth:`requirements`.
    label
        Component label for this step within a chain.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Image,)
    PRODUCES: ClassVar[type] = VisibilitySet

    def __init__(
        self,
        u_pts: Any,
        v_pts: Any,
        wavelength: Any,
        *,
        field_of_view: Any,
        oversampling: float = 1.0,
        label: str | None = None,
    ) -> None:
        super().__init__(label=label)
        coverage = [np.asarray(axis, dtype=DTYPE).reshape(-1) for axis in (u_pts, v_pts)]
        waves = _to_unit(wavelength, SPECTRAL_UNIT).reshape(-1)
        sizes = {"u": coverage[0].size, "v": coverage[1].size, "wavelength": waves.size}
        if len(set(sizes.values())) != 1 or coverage[0].size == 0:
            raise TransformationError(
                f"FourierSample takes one u, one v and one wavelength per observed sample, got "
                f"{sizes}. These are the observed container's own coordinates; build the step "
                f"with FourierSample.from_observed(container, field_of_view=...) rather than "
                f"assembling them by hand."
            )
        if not np.all(waves > 0.0):
            raise TransformationError("FourierSample needs strictly positive wavelengths (micron).")
        fov = float(_to_unit(field_of_view, COORDINATE_UNIT))
        if not np.isfinite(fov) or fov <= 0.0:
            raise TransformationError(
                f"FourierSample's field_of_view is the angular extent of the region the "
                f"instrument sees, in mas, and must be finite and positive, got "
                f"{field_of_view!r}."
            )
        self.register_buffer("u_pts", coverage[0])
        self.register_buffer("v_pts", coverage[1])
        self.register_buffer("wavelength", waves, unit=SPECTRAL_UNIT)
        factor = float(oversampling)
        if not np.isfinite(factor) or factor < 1.0:
            raise TransformationError(
                f"FourierSample's oversampling multiplies the Nyquist sampling rate the array "
                f"itself needs, so it is at least 1, got {oversampling!r}."
            )
        self.register_buffer("field_of_view", fov, unit=COORDINATE_UNIT)
        self.register_buffer("oversampling", factor)
        self._expanded: tuple[np.ndarray, np.ndarray, np.ndarray] = (
            coverage[0],
            coverage[1],
            waves,
        )
        self._template: VisibilitySet | None = None

    # -- construction from the data -----------------------------------------

    @classmethod
    def from_observed(
        cls,
        container: VisibilitySet | ClosurePhases,
        *,
        field_of_view: Any,
        oversampling: float = 1.0,
        label: str | None = None,
    ) -> FourierSample:
        """Take the coverage from the observed container. **The supported route.**

        For a :class:`~ampere.core.VisibilitySet` the coverage is the
        container's own three axes, unchanged.

        For a :class:`~ampere.core.ClosurePhases` it is the **three baselines
        of every triangle**, in the canonical order that kind fixes: ``ij``,
        ``jk``, and the implied ``ki = -(ij + jk)``, laid out so that triangle
        ``t``'s baselines occupy positions ``3t``, ``3t + 1``, ``3t + 2``. That
        layout is what :class:`ClosurePhase` reads back, and the two are
        written here together so that they cannot drift apart.
        """
        if isinstance(container, ClosurePhases):
            u3, v3 = container.implied_baseline()
            u_pts = np.stack([container.u1.values, container.u2.values, u3], axis=1).reshape(-1)
            v_pts = np.stack([container.v1.values, container.v2.values, v3], axis=1).reshape(-1)
            waves = np.repeat(container.spectral_axis.values, 3)
            wavelength: Any = waves * (container.spectral_axis.unit or SPECTRAL_UNIT)
            return cls(
                u_pts,
                v_pts,
                wavelength,
                field_of_view=field_of_view,
                oversampling=oversampling,
                label=label,
            )
        if isinstance(container, VisibilitySet):
            axis = container.spectral_axis
            return cls(
                container.u.values,
                container.v.values,
                axis.values * (axis.unit or SPECTRAL_UNIT),
                field_of_view=field_of_view,
                oversampling=oversampling,
                label=label,
            )
        raise TransformationError(
            f"FourierSample.from_observed takes the coverage from a VisibilitySet or a "
            f"ClosurePhases, got a {type(container).__name__}. Those are the two kinds that "
            f"carry (u, v) coordinates; for anything else, pass the arrays explicitly."
        )

    # -- chain-internal negotiation ------------------------------------------

    def configure_from(self, downstream: Sequence[Transformation]) -> None:
        """Ask the steps after this one which extra ``(u, v)`` samples they need.

        Gap I-3's mechanism, used across a kind change, which is what makes it
        the first *cross-kind* instance (``LSFConvolution`` already uses
        ``configure_from`` within one kind, on all three backends). A smearing
        step downstream of this one cannot publish a requirement — the channel
        is an ``Image`` and requirements name the container's own axes — so it
        states its need to this step instead, and this step is the one that
        computes the extra samples, because the coverage is its buffer.

        The expansions compose as a product, in chain order: with factors
        ``k1, k2, ...`` the flat output is indexed ``((i k1 + a) k2 + b)...``,
        and index ``0`` along every sub-axis is the un-smeared sample, so the
        first averaging step reduces the first factor and what remains is what
        the next one expects. Nothing is inferred backwards; only the
        successors' own declarations are read, and a chain with no averaging
        step leaves the coverage exactly as the data gave it.
        """
        u_pts = self._data("u_pts")
        v_pts = self._data("v_pts")
        waves = self._data("wavelength")
        for step in downstream:
            expander = getattr(step, "expand_uv", None)
            if expander is None:
                continue
            expansion = expander(u_pts, v_pts, waves)
            u_pts, v_pts, waves = expansion.u, expansion.v, expansion.wavelength
        self._expanded = (u_pts.reshape(-1), v_pts.reshape(-1), waves.reshape(-1))
        self._template = None

    @property
    def expanded_coverage(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """The ``(u, v, lambda)`` this step actually evaluates, smearing included."""
        return self._expanded

    # -- requirements --------------------------------------------------------

    def requirements(self) -> tuple[AxisRequirement, ...]:
        """Cover the field of view, sampled at or below the Nyquist pixel scale.

        ``intervals = (-fov/2, +fov/2)`` on both axes, and
        ``max_step = 1/(2 s u_max)`` radians on each, in mas — with ``u_max``
        taken over the **expanded** coverage, so that a smearing step's extra
        samples cannot reach beyond the grid the model was asked for, and
        ``s`` the :attr:`oversampling` factor.

        A model that engages with this and cannot honour it must raise
        (``transformations.md`` §7, gap I-4's loud option): an under-sampled
        image aliases, and aliased visibilities look like real source
        structure rather than like an error.

        **What ``oversampling = 1`` is and is not.** ``1/(2 u_max)`` is the
        classical Nyquist step, and it is exactly the right answer for what
        the *array* can measure: finer buys no measurable information, coarser
        cannot represent the longest baseline. It is **not** a statement about
        the accuracy of this quadrature, and the difference matters. A sky
        with real power beyond ``u_max`` — a sharp edge, a companion far more
        compact than the beam — has that power folded back onto the sampled
        baselines by a sum at the Nyquist step, where it is indistinguishable
        from source structure: gap I-4's failure mode arriving through the
        model's own compactness rather than through a coarse grid. The array
        does not measure the folded power, but this transform invents it, so
        the choice belongs to whoever knows how compact the model is. Raise
        ``oversampling`` until the visibilities stop moving; a factor of a few
        suffices for a source comparable to the beam, and the band-limited
        conformance rows measure exactly that convergence.
        """
        half = 0.5 * float(self.buffers["field_of_view"].value)
        factor = float(self.buffers["oversampling"].value)
        u_pts, v_pts, _ = self._expanded
        return tuple(
            AxisRequirement(
                axis,
                unit=COORDINATE_UNIT,
                intervals=(-half, half),
                max_step=_nyquist_step(values, factor),
            )
            for axis, values in (("x", u_pts), ("y", v_pts))
        )

    # -- evaluation ----------------------------------------------------------

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        x_mas = samples.x.values
        y_mas = samples.y.values
        brightness = np.asarray(samples.values, dtype=DTYPE) * cell_solid_angle(x_mas, y_mas)
        u_pts, v_pts, waves = self._expanded
        # Separable, because the DFT kernel is: one (n_uv, ny) contraction and
        # one (n_uv, nx) row product, rather than an (n_uv, nx, ny) array that
        # would be the same arithmetic with the memory of a bad idea.
        x_rad = x_mas / MAS_PER_RAD
        y_rad = y_mas / MAS_PER_RAD
        along_y = np.exp(-2j * np.pi * np.outer(v_pts, y_rad)) @ brightness.T
        along_x = np.exp(-2j * np.pi * np.outer(u_pts, x_rad))
        visibility = np.einsum("nj,nj->n", along_x, along_y)
        template = self._resolve_template(samples.unit, u_pts, v_pts, waves)
        mask = propagate_mask(samples, _image_influence(samples, u_pts.size))
        return template.with_values(visibility, mask=mask)

    def _resolve_template(
        self,
        brightness_unit: u.UnitBase | None,
        u_pts: np.ndarray,
        v_pts: np.ndarray,
        waves: np.ndarray,
    ) -> VisibilitySet:
        """The output container, built once and refilled thereafter.

        The compile-once/evaluate-many split (``results_schema.md`` §10): the
        axes are validated once and shared by identity across draws. The value
        unit is the image's own multiplied by a steradian — a surface
        brightness integrated over solid angle is a flux — and it is checked
        rather than assumed, because a template built for Jy/sr and refilled
        from an image in some other unit would be the units trap in its most
        damaging form.
        """
        unit = None if brightness_unit is None else brightness_unit * u.sr
        if self._template is not None:
            if self._template.unit == unit:
                return self._template
            raise TransformationError(
                f"FourierSample was compiled for an image in "
                f"{self._template.unit if self._template.unit is not None else 'no unit'} per "
                f"steradian and is now being given one in "
                f"{brightness_unit if brightness_unit is not None else 'no unit'}. A container's "
                f"unit is fixed once at composition time; convert the model's output with "
                f".to_unit(...) rather than changing it per evaluation."
            )
        self._template = VisibilitySet(
            u_pts,
            v_pts,
            waves * SPECTRAL_UNIT,
            np.zeros(u_pts.size, dtype=np.complex128),
            unit=unit,
        )
        return self._template


def _nyquist_step(values: np.ndarray, oversampling: float) -> float | None:
    """``1/(2 s |u|_max)`` in mas, or ``None`` when the coverage is degenerate."""
    extent = float(np.max(np.abs(values))) if values.size else 0.0
    if extent <= 0.0:
        return None
    return MAS_PER_RAD / (2.0 * oversampling * extent)


def _image_influence(samples: Any, n_out: int) -> np.ndarray | None:
    """Every output visibility depends on every image pixel — or on none, unmasked.

    A Fourier transform is dense: one masked pixel touches every sample. That
    makes the honest influence matrix a block of ones, and the honest
    consequence — a single masked pixel masks the whole visibility set — is
    correct rather than convenient. Returns ``None`` for an unmasked image, so
    the usual case pays nothing.
    """
    if samples.mask is None:
        return None
    return np.ones((n_out, samples.values.size), dtype=DTYPE)


class ClosurePhase(_Step):
    """Three visibilities to one angle: ``VisibilitySet -> ClosurePhases``.

    The sum of the visibility phases around a closed triangle, which is the
    interferometric observable that survives the atmosphere: a per-telescope
    phase error enters two of the three baselines with opposite signs and
    cancels exactly.

    The input must be the three baselines of every triangle in the canonical
    order :class:`~ampere.core.ClosurePhases` fixes and
    :meth:`FourierSample.from_observed` produces — triangle ``t`` at positions
    ``3t``, ``3t + 1``, ``3t + 2``, being ``ij``, ``jk`` and the implied
    ``ki = -(ij + jk)``. The third baseline is *checked*, not assumed: a chain
    assembled by hand from the wrong coverage fails here by name rather than
    producing plausible angles for the wrong triangles.

    The output is ``arg(V_ij V_jk V_ki)`` in radians, in ``(-pi, pi]`` —
    :func:`numpy.angle`'s own range, and the range the wrapped family scores
    on. Identity labels (``triangle``, the constituent ``baseline`` names)
    are **not** synthesised here: they annotate the observed data, and a
    predicted container that invented them would be asserting something it
    does not know.

    Mask propagation is the many-to-one rule (``results_schema.md`` §16, as
    :func:`~ampere.core.propagate_mask` ratifies it): a triangle is masked if
    **any** of its three baselines is. One flagged baseline therefore masks
    every triangle it takes part in, which is what a closure phase built from
    an unusable visibility deserves.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (VisibilitySet,)
    PRODUCES: ClassVar[type] = ClosurePhases

    #: Tolerance on the implied third baseline, relative to the triangle's own
    #: longest baseline. Not exact equality: the third baseline is computed as
    #: a sum, and a caller who supplies it from a file will differ in the last
    #: bits.
    CLOSURE_RTOL: ClassVar[float] = 1e-9

    def __init__(self, *, label: str | None = None) -> None:
        super().__init__(label=label)
        self._template: ClosurePhases | None = None

    def apply(self, samples: Any, values: Any) -> ClosurePhases:
        total = int(samples.n_samples)
        if total % 3:
            raise TransformationError(
                f"ClosurePhase reads three visibilities per triangle, so it needs a multiple of "
                f"three, got {total}. Build the Fourier step before it with "
                f"FourierSample.from_observed(closure_phases, field_of_view=...), which lays the "
                f"three baselines of each triangle out in the canonical order."
            )
        n_out = total // 3
        triples = np.asarray(samples.values).reshape(n_out, 3)
        product = triples[:, 0] * triples[:, 1] * triples[:, 2]
        phase = np.angle(product)
        template = self._resolve_template(samples, n_out)
        mask = propagate_mask(samples, self._influence(n_out) if samples.mask is not None else None)
        return template.with_values(phase, mask=mask)

    def _resolve_template(self, samples: Any, n_out: int) -> ClosurePhases:
        if self._template is not None:
            return self._template
        u_pts = samples.u.values.reshape(n_out, 3)
        v_pts = samples.v.values.reshape(n_out, 3)
        waves = samples.spectral_axis.values.reshape(n_out, 3)
        self._check_canonical(u_pts, v_pts, waves)
        self._template = ClosurePhases(
            u_pts[:, 0],
            v_pts[:, 0],
            u_pts[:, 1],
            v_pts[:, 1],
            waves[:, 0] * (samples.spectral_axis.unit or SPECTRAL_UNIT),
            np.zeros(n_out, dtype=DTYPE) * u.rad,
        )
        return self._template

    def _check_canonical(self, u_pts: np.ndarray, v_pts: np.ndarray, waves: np.ndarray) -> None:
        """Refuse coverage that is not three baselines of one triangle at one wavelength."""
        scale = np.maximum(np.max(np.abs(u_pts), axis=1), np.max(np.abs(v_pts), axis=1))
        scale = np.where(scale > 0.0, scale, 1.0)
        residual = np.maximum(np.abs(u_pts.sum(axis=1)) / scale, np.abs(v_pts.sum(axis=1)) / scale)
        bad = np.flatnonzero(residual > self.CLOSURE_RTOL)
        if bad.size:
            raise TransformationError(
                f"ClosurePhase was given coverage whose three baselines do not close: triangle "
                f"{int(bad[0])} has (u1+u2+u3, v1+v2+v3) = "
                f"({float(u_pts[bad[0]].sum()):g}, {float(v_pts[bad[0]].sum()):g}), which is not "
                f"zero to {self.CLOSURE_RTOL:g} of its longest baseline. The canonical ordering "
                f"stores ij and jk and implies ki as their negated sum (ClosurePhases's "
                f"docstring); FourierSample.from_observed produces it."
            )
        if not np.all(waves == waves[:, :1]):
            raise TransformationError(
                "ClosurePhase was given a triangle whose three baselines sit at different "
                "wavelengths. A closure phase is formed from three simultaneous, co-spectral "
                "visibilities; different wavelengths are different triangles."
            )

    @staticmethod
    def _influence(n_out: int) -> np.ndarray:
        """The ``(n_out, 3 n_out)`` three-to-one adjacency, for the mask."""
        influence = np.zeros((n_out, 3 * n_out), dtype=DTYPE)
        rows = np.arange(n_out).reshape(-1, 1)
        np.put_along_axis(influence, 3 * rows + np.arange(3).reshape(1, -1), 1.0, axis=1)
        return influence


class Amplitude(_Step):
    """Take the modulus: a complex ``VisibilitySet`` to a real one.

    Kind-preserving, coordinate-preserving, parameter-free — and the reason
    :class:`~ampere.core.RiceFamily` receives real amplitudes rather than a
    complex prediction (``interferometry.md`` §6, ruled 2026-09-03). Fitting
    amplitudes rather than complex visibilities is a real scientific decision
    — you are discarding the astrometric information because the phases are
    corrupted — and making it a named step in the chain rather than a
    convention buried in a family is what ``DEVELOPMENT_PLAN.md`` §4.3 is for.

    The mask comes through :meth:`~ampere.core.FunctionSamples.with_values`,
    which inherits it: the mapping is one-to-one.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (VisibilitySet,)

    def apply(self, samples: Any, values: Any) -> VisibilitySet:
        return samples.with_values(np.abs(np.asarray(samples.values)))


# ---------------------------------------------------------------------------
# The models
# ---------------------------------------------------------------------------


class _ImageModel(Model):
    """Shared plumbing for the image-emitting trio: the grid, the channels, the template.

    Each subclass writes one method, :meth:`_brightness`, returning a surface
    brightness in Jy/sr on the ``(x, y)`` grid it is handed, in mas. The base
    handles the negotiated grid, the channel names and the hot-loop template,
    exactly as ``models.py``'s ``_SpectralModel`` does for the spectral axis.

    Sky convention, stated once: ``x`` increases towards the east and ``y``
    towards the north, both offsets from the phase centre in mas, and a
    position angle is measured **east of north**, so an offset at position
    angle ``p`` and separation ``s`` is ``(s sin p, s cos p)``. The convention
    matters because it fixes the sign of every closure phase these models
    predict.

    Parameters
    ----------
    x, y
        The model's own grid, mas. Used when nothing is negotiated, and when
        ``adopt_grid`` is ``False``.
    channels
        Name (or names) of the channel the emitted ``Image`` appears under.
    adopt_grid
        Whether :meth:`compile_for` rebuilds onto the negotiated grid (the
        default) or holds the grid it was given. ``False`` stands in for a
        model whose grid is not ours to choose — a radiative-transfer code
        with a fixed image size — and it then **refuses** a requirement it
        cannot honour rather than aliasing silently (gap I-4).
    """

    #: The fourth capability flag (W2.12), declared rather than inherited.
    BACKEND: ClassVar[str] = "reference"

    def __init__(
        self,
        x: Any,
        y: Any,
        *,
        channels: str | Sequence[str] = "sky",
        adopt_grid: bool = True,
    ) -> None:
        grids = [_to_unit(axis, COORDINATE_UNIT).reshape(-1) for axis in (x, y)]
        for name, grid in zip(("x", "y"), grids, strict=True):
            if grid.size < 2 or not np.all(np.diff(grid) > 0.0):
                raise ValueError(
                    f"an image model needs a strictly increasing {name} grid of at least two "
                    f"coordinates, in mas; got {grid.size}."
                )
        names = (channels,) if isinstance(channels, str) else tuple(str(c) for c in channels)
        if not names or len(set(names)) != len(names):
            raise ValueError(
                f"an image model needs distinct, non-empty channel names, got {names!r}."
            )
        self.channels = names
        self.adopt_grid = bool(adopt_grid)
        self.register_buffer("x", grids[0], unit=COORDINATE_UNIT)
        self.register_buffer("y", grids[1], unit=COORDINATE_UNIT)
        self.templates: dict[str, Image] = {}

    @classmethod
    def on_field(
        cls,
        field_of_view: Any,
        pixels: int,
        **kwargs: Any,
    ) -> Any:
        """A square grid of *pixels* across *field_of_view*, as a convenience."""
        half = 0.5 * float(_to_unit(field_of_view, COORDINATE_UNIT))
        grid = np.linspace(-half, half, int(pixels))
        return cls(grid, grid, **kwargs)

    # -- negotiation ---------------------------------------------------------

    def compile_for(self, requirements: Mapping[str, ChannelRequirements]) -> Model:
        """Adopt the negotiated image grid for each channel, or refuse it by name.

        The union of every bound instrument's ``x`` and ``y`` requirements is
        one grid, built once and refilled with ``with_values`` on every
        evaluation, so two instruments reading the same channel share one model
        evaluation — which is the whole argument for negotiation
        (``DEVELOPMENT_PLAN.md`` §4.3).

        With ``adopt_grid=False`` the model keeps its own grid and *checks* it
        instead. That is gap I-4's case and the reason
        :class:`~ampere.core.CompositionError` is raised rather than warned:
        an image coarser than ``1/(2 u_max)`` folds power from beyond the
        Nyquist limit back onto the sampled baselines, where it is
        indistinguishable from real source structure. The caller who wants to
        proceed anyway says so explicitly, with
        ``FittingProblem(lenient_compile=True)``.
        """
        for channel in self.channels:
            asked = requirements.get(channel)
            if asked is None:
                continue
            axes = {name: asked[name] for name in ("x", "y") if name in asked}
            if not axes:
                continue
            if self.adopt_grid:
                grids = {
                    name: _to_unit(requirement.coordinates(), COORDINATE_UNIT)
                    for name, requirement in axes.items()
                }
            else:
                grids = {name: self._data(name) for name in axes}
                for name, requirement in axes.items():
                    _refuse_unless_honoured(type(self).__name__, channel, name, requirement, grids)
            x_grid = grids.get("x", self._data("x"))
            y_grid = grids.get("y", self._data("y"))
            self.templates[channel] = Image(
                x_grid * COORDINATE_UNIT,
                y_grid * COORDINATE_UNIT,
                np.zeros((x_grid.size, y_grid.size), dtype=DTYPE),
                unit=BRIGHTNESS_UNIT,
            )
        return self

    # -- evaluation ----------------------------------------------------------

    def _data(self, name: str) -> np.ndarray:
        return np.asarray(self.buffers[name].value, dtype=DTYPE)

    def _grid(self, channel: str, context: Mapping[str, Any]) -> tuple[np.ndarray, np.ndarray]:
        template = self.templates.get(channel)
        if template is None:
            return np.asarray(context["x"], dtype=DTYPE), np.asarray(context["y"], dtype=DTYPE)
        return template.x.values, template.y.values

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        raise NotImplementedError

    def _emit(self, channel: str, grid: tuple[np.ndarray, np.ndarray], values: np.ndarray) -> Image:
        template = self.templates.get(channel)
        if template is None:
            return Image(
                grid[0] * COORDINATE_UNIT,
                grid[1] * COORDINATE_UNIT,
                values,
                unit=BRIGHTNESS_UNIT,
            )
        return template.with_values(values)

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        emitted = {
            channel: self._emit(channel, grid, self._brightness(grid[0], grid[1], context))
            for channel in self.channels
            for grid in (self._grid(channel, context),)
        }
        return ModelResult(emitted)


def _refuse_unless_honoured(
    model: str,
    channel: str,
    axis: str,
    requirement: AxisRequirement,
    grids: Mapping[str, np.ndarray],
) -> None:
    """Raise unless the model's own grid satisfies *requirement* (gap I-4)."""
    grid = grids[axis]
    steps = np.diff(grid)
    coarsest = float(np.max(np.abs(steps))) if steps.size else np.inf
    for low, high, step, _power in requirement.convert_to(COORDINATE_UNIT).segments():
        if grid[0] > low or grid[-1] < high:
            raise CompositionError(
                f"{model} holds its own image grid, and channel {channel!r} was asked to cover "
                f"{axis} in [{low:g}, {high:g}] mas while the model's grid spans "
                f"[{grid[0]:g}, {grid[-1]:g}]. A source outside the modelled field still "
                f"contributes to every visibility, so the missing coverage is not a small error. "
                f"Widen the model's grid, build it with adopt_grid=True, or say "
                f"FittingProblem(lenient_compile=True) to proceed with the unconfigured model."
            )
        if step is not None and coarsest > step:
            raise CompositionError(
                f"{model} holds its own image grid, and channel {channel!r} was asked for a {axis} "
                f"step of at most {step:g} mas while the model's coarsest is {coarsest:g} mas. "
                f"That requirement is the Nyquist limit of the longest baseline observed: an "
                f"image sampled more coarsely *aliases*, folding power from beyond the limit back "
                f"onto the sampled baselines, where it is indistinguishable from real source "
                f"structure rather than looking like an error. Refine the model's grid, build it "
                f"with adopt_grid=True, or say FittingProblem(lenient_compile=True) to proceed "
                f"with the unconfigured model."
            )


class UniformDisc(_ImageModel):
    """A uniformly bright circular disc — the standard stellar-diameter model.

    ``flux`` is the total flux density in Jy and ``diameter`` the angular
    diameter in mas; the brightness inside the disc is
    ``flux / (pi (theta/2)**2)`` in Jy/sr and zero outside.

    Its visibility has the closed form ``flux * 2 J1(pi theta rho) / (pi theta
    rho)`` (:class:`UniformDiscVisibilities`), which is an oracle for the
    Fourier step **to the accuracy of the quadrature and no further**: a
    sharp-edged disc is not band-limited, so a rectangle rule over a finite
    grid has aliases that fall only as a power of the pixel scale. That is a
    fact about discs rather than about this implementation, and the
    conformance row says so by measuring the convergence rather than by
    asserting a tolerance the shape cannot meet. The band-limited members of
    the trio — :class:`GaussianSource` and :class:`Binary` — agree with their
    closed forms to the solver tolerance.

    Parameters
    ----------
    x, y
        The model's own grid, mas.
    diameter
        Angular diameter, mas. A prior to fit it, a number to hold it fixed.
    flux
        Total flux density, Jy.
    channels, adopt_grid
        See :class:`_ImageModel`.
    """

    def __init__(
        self,
        x: Any,
        y: Any,
        *,
        diameter: Any = 3.0,
        flux: Any = 1.0,
        channels: str | Sequence[str] = "sky",
        adopt_grid: bool = True,
    ) -> None:
        super().__init__(x, y, channels=channels, adopt_grid=adopt_grid)
        self.register_parameter(as_parameter("diameter", diameter, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        radius = 0.5 * float(context["diameter"])
        if radius <= 0.0:
            raise ValueError(f"a uniform disc needs a positive diameter, got {radius * 2.0!r} mas.")
        distance = np.hypot(x_mas[:, None], y_mas[None, :])
        area = np.pi * (radius / MAS_PER_RAD) ** 2
        return np.where(distance <= radius, float(context["flux"]) / area, 0.0)


class GaussianSource(_ImageModel):
    """A circular Gaussian source — the standard resolved-but-featureless model.

    ``flux`` is the total flux density in Jy and ``fwhm`` the full width at
    half maximum in mas. Band-limited, so the direct DFT of its image agrees
    with its closed-form visibility (:class:`GaussianSourceVisibilities`) to
    the solver tolerance on any grid that satisfies the Nyquist requirement
    and covers the source.

    Parameters
    ----------
    x, y
        The model's own grid, mas.
    fwhm
        Full width at half maximum, mas.
    flux
        Total flux density, Jy.
    channels, adopt_grid
        See :class:`_ImageModel`.
    """

    def __init__(
        self,
        x: Any,
        y: Any,
        *,
        fwhm: Any = 3.0,
        flux: Any = 1.0,
        channels: str | Sequence[str] = "sky",
        adopt_grid: bool = True,
    ) -> None:
        super().__init__(x, y, channels=channels, adopt_grid=adopt_grid)
        self.register_parameter(as_parameter("fwhm", fwhm, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        sigma = float(context["fwhm"]) / _FWHM_PER_SIGMA
        if sigma <= 0.0:
            raise ValueError(f"a Gaussian source needs a positive fwhm, got {context['fwhm']!r}.")
        return _gaussian_brightness(x_mas, y_mas, 0.0, 0.0, sigma, float(context["flux"]))


class Binary(_ImageModel):
    """Two Gaussian components — the synthetic source the phase is proved on.

    The primary sits at the phase centre and the secondary at separation
    ``separation`` and position angle ``position_angle`` (radians, east of
    north); ``flux_ratio`` is the secondary's flux over the primary's and
    ``flux`` is the pair's total. Both components share one ``component_fwhm``,
    a **buffer** rather than a parameter: it is there to make the source
    band-limited — a true point pair is not representable on any grid — and
    one would not put a prior on a numerical device.

    A binary is the right source to prove a closure phase on, because its
    closure phases are non-zero, analytic, and sensitive to exactly the sign
    conventions that would otherwise go unnoticed: a point-symmetric source
    has closure phases of zero or ``pi`` and would agree with a mirrored
    model.

    Parameters
    ----------
    x, y
        The model's own grid, mas.
    separation
        Component separation, mas.
    position_angle
        Position angle of the secondary, radians east of north.
    flux_ratio
        Secondary flux over primary flux, dimensionless.
    flux
        Total flux density of the pair, Jy.
    component_fwhm
        Full width at half maximum of each component, mas. A buffer.
    channels, adopt_grid
        See :class:`_ImageModel`.
    """

    def __init__(
        self,
        x: Any,
        y: Any,
        *,
        separation: Any = 5.0,
        position_angle: Any = 0.5,
        flux_ratio: Any = 0.4,
        flux: Any = 1.0,
        component_fwhm: float = 0.5,
        channels: str | Sequence[str] = "sky",
        adopt_grid: bool = True,
    ) -> None:
        super().__init__(x, y, channels=channels, adopt_grid=adopt_grid)
        width = float(_to_unit(component_fwhm, COORDINATE_UNIT))
        if not np.isfinite(width) or width <= 0.0:
            raise ValueError(
                f"a binary's component_fwhm makes the pair band-limited and must be finite and "
                f"positive (mas), got {component_fwhm!r}."
            )
        self.register_buffer("component_fwhm", width, unit=COORDINATE_UNIT)
        self.register_parameter(as_parameter("separation", separation, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("position_angle", position_angle, unit=u.rad))
        self.register_parameter(as_parameter("flux_ratio", flux_ratio))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _brightness(
        self, x_mas: np.ndarray, y_mas: np.ndarray, context: Mapping[str, Any]
    ) -> np.ndarray:
        offset_x, offset_y = _companion_offset(context)
        ratio = float(context["flux_ratio"])
        total = float(context["flux"])
        sigma = float(context["component_fwhm"]) / _FWHM_PER_SIGMA
        primary = total / (1.0 + ratio)
        return _gaussian_brightness(x_mas, y_mas, 0.0, 0.0, sigma, primary) + _gaussian_brightness(
            x_mas, y_mas, offset_x, offset_y, sigma, primary * ratio
        )


def _companion_offset(context: Mapping[str, Any]) -> tuple[float, float]:
    """The secondary's ``(x, y)`` offset in mas, east of north."""
    separation = float(context["separation"])
    angle = float(context["position_angle"])
    return separation * np.sin(angle), separation * np.cos(angle)


def _gaussian_brightness(
    x_mas: np.ndarray,
    y_mas: np.ndarray,
    centre_x: float,
    centre_y: float,
    sigma_mas: float,
    flux: float,
) -> np.ndarray:
    """A circular Gaussian of total *flux* (Jy) as a surface brightness in Jy/sr."""
    sigma_rad = sigma_mas / MAS_PER_RAD
    squared = ((x_mas - centre_x)[:, None] ** 2 + (y_mas - centre_y)[None, :] ** 2) / sigma_mas**2
    return flux / (2.0 * np.pi * sigma_rad**2) * np.exp(-0.5 * squared)


# ---------------------------------------------------------------------------
# The analytic route: the same three sources, emitting visibilities directly
# ---------------------------------------------------------------------------


class _VisibilityModel(Model):
    """Shared plumbing for the trio that skips the Fourier step.

    ``interferometry.md`` §1: a model that computes visibilities analytically
    emits a :class:`~ampere.core.VisibilitySet` channel directly, and its
    instrument is ``Instrument([], channel=..., input_kind=VisibilitySet)`` —
    the "no steps" instrument ``transformations.md`` §5 provides. **Both
    routes are supported and neither is privileged**, and here each is the
    other's oracle: these closed forms are what the conformance battery holds
    the direct DFT to.

    The coverage is a buffer taken from the observed container
    (:meth:`from_observed`), for the reason :class:`FourierSample`'s is.
    These models publish no requirements and honour no negotiation: they
    already produce the observable, on the coordinates the data fix.
    """

    #: The fourth capability flag (W2.12), declared rather than inherited.
    BACKEND: ClassVar[str] = "reference"

    def __init__(
        self,
        u_pts: Any,
        v_pts: Any,
        wavelength: Any,
        *,
        channels: str | Sequence[str] = "vis",
    ) -> None:
        coverage = [np.asarray(axis, dtype=DTYPE).reshape(-1) for axis in (u_pts, v_pts)]
        waves = _to_unit(wavelength, SPECTRAL_UNIT).reshape(-1)
        if len({coverage[0].size, coverage[1].size, waves.size}) != 1 or waves.size == 0:
            raise ValueError(
                "an analytic visibility model needs one u, one v and one wavelength per sample."
            )
        names = (channels,) if isinstance(channels, str) else tuple(str(c) for c in channels)
        if not names or len(set(names)) != len(names):
            raise ValueError(
                f"a visibility model needs distinct, non-empty channel names, got {names!r}."
            )
        self.channels = names
        self.register_buffer("u_pts", coverage[0])
        self.register_buffer("v_pts", coverage[1])
        self.register_buffer("wavelength", waves, unit=SPECTRAL_UNIT)
        self.templates = {
            channel: VisibilitySet(
                coverage[0],
                coverage[1],
                waves * SPECTRAL_UNIT,
                np.zeros(waves.size, dtype=np.complex128),
                unit=FLUX_UNIT,
            )
            for channel in names
        }

    @classmethod
    def from_observed(
        cls,
        container: VisibilitySet | ClosurePhases,
        **kwargs: Any,
    ) -> Any:
        """Take the coverage from the observed container, exactly as the step does.

        A :class:`~ampere.core.ClosurePhases` container gives the three
        baselines of each triangle in the canonical order, so that the same
        :class:`ClosurePhase` step consumes this model's output as consumes
        :class:`FourierSample`'s.
        """
        step = FourierSample.from_observed(container, field_of_view=1.0)
        u_pts, v_pts, waves = step.expanded_coverage
        return cls(u_pts, v_pts, waves * SPECTRAL_UNIT, **kwargs)

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        raise NotImplementedError

    def evaluate(self, **values: Any) -> ModelResult:
        context = self.context(values)
        u_pts = np.asarray(context["u_pts"], dtype=DTYPE)
        v_pts = np.asarray(context["v_pts"], dtype=DTYPE)
        visibility = self._visibility(np.hypot(u_pts, v_pts), context)
        return ModelResult(
            {channel: self.templates[channel].with_values(visibility) for channel in self.channels}
        )


class UniformDiscVisibilities(_VisibilityModel):
    """``V = flux * 2 J1(pi theta rho) / (pi theta rho)`` — the disc's closed form.

    ``theta`` is the angular diameter in radians and ``rho = hypot(u, v)`` the
    baseline length in wavelengths. Real-valued, and negative beyond the first
    null, which is the sign flip an amplitude-only fit throws away.
    """

    def __init__(
        self,
        u_pts: Any,
        v_pts: Any,
        wavelength: Any,
        *,
        diameter: Any = 3.0,
        flux: Any = 1.0,
        channels: str | Sequence[str] = "vis",
    ) -> None:
        super().__init__(u_pts, v_pts, wavelength, channels=channels)
        self.register_parameter(as_parameter("diameter", diameter, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        argument = np.pi * float(context["diameter"]) / MAS_PER_RAD * rho
        # 2 J1(z)/z -> 1 as z -> 0, and the ratio is 0/0 there; jv handles the
        # numerator but not the limit, so it is written in rather than clipped.
        envelope = np.where(argument == 0.0, 1.0, 2.0 * scipy.special.j1(argument) / argument)
        return (float(context["flux"]) * envelope).astype(np.complex128)


class GaussianSourceVisibilities(_VisibilityModel):
    """``V = flux * exp(-2 pi**2 sigma**2 rho**2)`` — the Gaussian's closed form."""

    def __init__(
        self,
        u_pts: Any,
        v_pts: Any,
        wavelength: Any,
        *,
        fwhm: Any = 3.0,
        flux: Any = 1.0,
        channels: str | Sequence[str] = "vis",
    ) -> None:
        super().__init__(u_pts, v_pts, wavelength, channels=channels)
        self.register_parameter(as_parameter("fwhm", fwhm, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        sigma = float(context["fwhm"]) / _FWHM_PER_SIGMA / MAS_PER_RAD
        return (float(context["flux"]) * np.exp(-2.0 * np.pi**2 * sigma**2 * rho**2)).astype(
            np.complex128
        )


class BinaryVisibilities(_VisibilityModel):
    """``V = (flux/(1+f)) G(rho) (1 + f exp(-2 pi i (u dx + v dy)))`` — the pair's closed form.

    The primary is at the phase centre and the secondary at ``(dx, dy)``
    radians, the same geometry :class:`Binary` renders; ``G`` is the common
    Gaussian envelope of the two components. The exponent carries the
    ``exp(-2 pi i)`` sign convention this module fixes, and it is what makes
    the closure phases non-zero.
    """

    def __init__(
        self,
        u_pts: Any,
        v_pts: Any,
        wavelength: Any,
        *,
        separation: Any = 5.0,
        position_angle: Any = 0.5,
        flux_ratio: Any = 0.4,
        flux: Any = 1.0,
        component_fwhm: float = 0.5,
        channels: str | Sequence[str] = "vis",
    ) -> None:
        super().__init__(u_pts, v_pts, wavelength, channels=channels)
        width = float(_to_unit(component_fwhm, COORDINATE_UNIT))
        if not np.isfinite(width) or width <= 0.0:
            raise ValueError(
                f"a binary's component_fwhm must be finite and positive (mas), got "
                f"{component_fwhm!r}."
            )
        self.register_buffer("component_fwhm", width, unit=COORDINATE_UNIT)
        self.register_parameter(as_parameter("separation", separation, unit=COORDINATE_UNIT))
        self.register_parameter(as_parameter("position_angle", position_angle, unit=u.rad))
        self.register_parameter(as_parameter("flux_ratio", flux_ratio))
        self.register_parameter(as_parameter("flux", flux, unit=FLUX_UNIT))

    def _visibility(self, rho: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        offset_x, offset_y = _companion_offset(context)
        u_pts = np.asarray(context["u_pts"], dtype=DTYPE)
        v_pts = np.asarray(context["v_pts"], dtype=DTYPE)
        sigma = float(context["component_fwhm"]) / _FWHM_PER_SIGMA / MAS_PER_RAD
        envelope = np.exp(-2.0 * np.pi**2 * sigma**2 * rho**2)
        phase = -2j * np.pi * (u_pts * offset_x + v_pts * offset_y) / MAS_PER_RAD
        ratio = float(context["flux_ratio"])
        return (
            float(context["flux"]) / (1.0 + ratio) * envelope * (1.0 + ratio * np.exp(phase))
        ).astype(np.complex128)
