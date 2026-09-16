"""A second in-repo backend, so the battery's parametrisation is not notional.

W1.10's acceptance criterion is that "adding a backend requires only a
fixture". A suite with one fixture cannot demonstrate that: the parametrisation
would be untested machinery, and the first Phase-2 author would find out the
hard way which rows had quietly baked in an assumption about the reference
path. This module is the cheapest honest demonstration — a second backend that
runs every row with **zero test-body changes**.

It is *not* a backend in ``architecture.md`` §1's sense, and it must not be
mistaken for one. It shares the reference backend's transformations, kernels
and solvers, because it has nothing different to offer there. What it does
supply of its own is exactly the two things a real backend owns:

* **models** — the same two closed forms declared under the same parameter
  names, computed by a deliberately different route (Horner's rule; ``exp`` of
  a logarithm) in classes of their own. Two backends that agreed only because
  they shared an implementation would prove nothing. (Being honest about how
  far that goes: the power law really does land on different bits — around
  ``3e-17`` on these grids — while Horner's rule on ``offset + slope * x`` is
  the same sequence of operations and agrees exactly. The mirror's value is
  that it is a *separate class the suite never names*, not that every one of
  its numbers differs in the last place.)
* **a parameter space** — one that round-trips the declaration through
  ``to_spec``/``from_spec`` before delegating, so every parameter row is also
  a serialisation row on this fixture.

Since W2.12 it also supplies **instrument steps of its own**, which it did not
before. The reason is item 4 of that work item: the backend string is one name
everywhere, so every part a fixture composes must declare
``BACKEND == fixture.name``, and the ``test_capabilities`` row that asserts it
would be a tautology if this fixture composed reference-backend steps under the
name ``"mirror"``. The steps are one-line subclasses that change the
declaration and nothing else — the arithmetic is still the reference
backend's, honestly, because there is nothing different to offer there.

Delete this module the day a real second backend registers, or keep it: it
costs one fixture and it is the only thing standing between the battery and
silently degenerating to a single-backend suite again.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import scipy.stats as st

from ampere.backends.reference import (
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    CalibrationScale,
    ClosurePhase,
    EpochSample,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    PSFConvolution,
    ReflexOrbit,
    Resample,
    TimeSmearing,
    UniformDisc,
    UniformDiscVisibilities,
)
from ampere.core import DenseGP as _CoreDenseGP
from ampere.core import HilbertSpaceGP as _CoreHilbertSpaceGP
from ampere.core import GaussianProcessNoise as _CoreGaussianProcessNoise
from ampere.core import IndependentNoise as _CoreIndependentNoise
from ampere.core import SHO as _CoreSHO
from ampere.core import Matern12 as _CoreMatern12
from ampere.core import Matern32 as _CoreMatern32
from ampere.core import Matern52 as _CoreMatern52
from ampere.core import RotationTerm as _CoreRotationTerm
from ampere.core import QuasisepGP as _CoreQuasisepGP
from ampere.core import SquaredExponential as _CoreSquaredExponential
from ampere.core import (
    GPSolver,
    Kernel,
    Model,
    NoiseModel,
    Parameter,
    ParameterSet,
    Transformation,
)

from ._kernels import build_kernel
from ..protocol import (
    AstrometryPieces,
    ImagePieces,
    BackendCapabilities,
    CovarianceSpec,
    InterferometryPieces,
    KernelFamily,
    ModelKind,
    ModelSpec,
    SolverKind,
    TransformationKind,
    TransformationSpec,
)
from .reference import Photometry, ReferenceBackend, _CountingModel

#: This fixture's one name, in the sense W2.12 fixed: ``ConformanceBackend.name``,
#: the ``BACKEND`` every part it composes declares, ``FittingProblem.backend``,
#: ``ampere_backend`` in a run's provenance, and the key ``lowering.md`` §12.8's
#: registry would be consulted with. One string, defined once here.
BACKEND = "mirror"

__all__ = [
    "BACKEND",
    "MIRROR_ASTROMETRY",
    "MIRROR_INTERFEROMETRY",
    "DenseGP",
    "GaussianProcessNoise",
    "IndependentNoise",
    "MirrorAmplitude",
    "MirrorBackend",
    "MirrorBandwidthSmearing",
    "MirrorBinary",
    "MirrorBinaryVisibilities",
    "MirrorCalibrationScale",
    "MirrorClosurePhase",
    "MirrorEpochSample",
    "MirrorFourierSample",
    "MirrorGaussianSource",
    "MirrorGaussianSourceVisibilities",
    "MirrorLinearModel",
    "MirrorParameterSpace",
    "MirrorPhotometry",
    "MirrorPowerLawModel",
    "MirrorReflexOrbit",
    "MirrorResample",
    "MirrorTimeSmearing",
    "MirrorUniformDisc",
    "MirrorUniformDiscVisibilities",
    "QuasisepGP",
]


class MirrorLinearModel(_CountingModel):
    """``offset + slope * x`` by Horner's rule rather than as written."""

    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        self.register_parameter(Parameter("offset", st.norm(0.0, 1.0)))
        self.register_parameter(Parameter("slope", st.norm(1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        return np.polyval([context["slope"], context["offset"]], grid)


class MirrorPowerLawModel(_CountingModel):
    """``norm * (x / x_ref) ** index`` through the logarithm."""

    BACKEND: ClassVar[str] = BACKEND

    def __init__(self, spec: ModelSpec) -> None:
        super().__init__(spec)
        # The unit is part of the declaration, and W2.1's neutral model
        # identity compares buffers: a reference wavelength declared unitless
        # here and in micron there is a real disagreement, not a detail.
        self.register_buffer(
            "reference_wavelength", float(spec.reference_coordinate), unit=u.micron
        )
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))

    def _flux(self, grid: np.ndarray, context: Mapping[str, Any]) -> np.ndarray:
        ratio = np.log(grid / context["reference_wavelength"])
        return np.exp(context["index"] * ratio + np.log(context["norm"]))


class MirrorParameterSpace:
    """A parameter space reached through the serialisable declaration.

    ``ParameterSet.to_spec()`` is what ``parameters.md`` §13 calls the
    serialisable declaration and what ``results.md`` hashes into the
    provenance attrs. Realising the space from *that* rather than from the
    object it was built from means every row in ``test_parameters.py`` also
    asserts, on this fixture, that the declaration survives the trip.
    """

    def __init__(self, declaration: ParameterSet) -> None:
        self._declaration = ParameterSet.from_spec(declaration.to_spec())

    @property
    def declaration(self) -> ParameterSet:
        return self._declaration

    @property
    def free_size(self) -> int:
        return self._declaration.free_size

    def free_labels(self) -> tuple[str, ...]:
        return self._declaration.free_labels()

    def pack(self, values: Mapping[str, Any]) -> np.ndarray:
        return self._declaration.pack(values)

    def unpack(self, theta: Any) -> dict[str, Any]:
        return self._declaration.unpack(theta)

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        return self._declaration.prior_transform(unit_cube)

    def lnprior(self, values: Any) -> float:
        return self._declaration.lnprior(values)

    def constrain(self, unconstrained: Any) -> np.ndarray:
        return self._declaration.constrain(unconstrained)

    def unconstrain(self, values: Any) -> np.ndarray:
        return self._declaration.unconstrain(values)

    def lnprior_unconstrained(self, unconstrained: Any) -> float:
        return self._declaration.lnprior_unconstrained(unconstrained)


class MirrorCalibrationScale(CalibrationScale):
    """The shipped calibration step, declared as this backend's (W2.12)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorResample(Resample):
    """The shipped resampling step, declared as this backend's (W2.12)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorPhotometry(Photometry):
    """The battery's kind-changing step, declared as this backend's (W2.12)."""

    BACKEND: ClassVar[str] = BACKEND


# The interferometric vocabulary (W4.1), for the reason the three steps above
# exist: a fixture named "mirror" that composed reference-backend classes
# would declare two backends and be refused, and the interferometric rows
# would be a tautology rather than a demonstration that they are parametrised
# at all. The arithmetic is the reference path's, honestly — a discrete
# Fourier transform has nothing different to offer here — and what these
# prove is that no row names a backend.


class MirrorFourierSample(FourierSample):
    """The Fourier step, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorClosurePhase(ClosurePhase):
    """The closure-phase step, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorAmplitude(Amplitude):
    """The modulus step, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorBandwidthSmearing(BandwidthSmearing):
    """Bandwidth smearing, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorTimeSmearing(TimeSmearing):
    """Time smearing, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorUniformDisc(UniformDisc):
    """The uniform disc, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorGaussianSource(GaussianSource):
    """The Gaussian source, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorBinary(Binary):
    """The binary, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorUniformDiscVisibilities(UniformDiscVisibilities):
    """The disc's closed form, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorGaussianSourceVisibilities(GaussianSourceVisibilities):
    """The Gaussian's closed form, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorBinaryVisibilities(BinaryVisibilities):
    """The binary's closed form, declared as this backend's (W4.1)."""

    BACKEND: ClassVar[str] = BACKEND


MIRROR_INTERFEROMETRY = InterferometryPieces(
    fourier_sample=MirrorFourierSample,
    closure_phase=MirrorClosurePhase,
    amplitude=MirrorAmplitude,
    bandwidth_smearing=MirrorBandwidthSmearing,
    time_smearing=MirrorTimeSmearing,
    uniform_disc=MirrorUniformDisc,
    gaussian_source=MirrorGaussianSource,
    binary=MirrorBinary,
    uniform_disc_visibilities=MirrorUniformDiscVisibilities,
    gaussian_source_visibilities=MirrorGaussianSourceVisibilities,
    binary_visibilities=MirrorBinaryVisibilities,
)


# The astrometric vocabulary (W4.9), for the same reason as the interferometric
# one above: a fixture named "mirror" that composed reference-backend classes
# would declare two backends and be refused.


class MirrorEpochSample(EpochSample):
    """The epoch-sampling step, declared as this backend's (W4.9)."""

    BACKEND: ClassVar[str] = BACKEND


class MirrorReflexOrbit(ReflexOrbit):
    """The reflex-orbit model, declared as this backend's (W4.9)."""

    BACKEND: ClassVar[str] = BACKEND


# The gridded-image vocabulary (W5.5), for the same reason as the two above.


class MirrorPSFConvolution(PSFConvolution):
    """The PSF-convolution step, declared as this backend's (W5.5)."""

    BACKEND: ClassVar[str] = BACKEND


MIRROR_IMAGE = ImagePieces(
    psf_convolution=MirrorPSFConvolution,
    gaussian_source=MirrorGaussianSource,
    binary=MirrorBinary,
)

#: The shipped gridded-image vocabulary, as one record (W5.5).
MIRROR_ASTROMETRY = AstrometryPieces(
    epoch_sample=MirrorEpochSample,
    reflex_orbit=MirrorReflexOrbit,
)


# The noise models and solvers, declared as this fixture's. Needed since
# **W2.13**: the four capability flags widened to ``NoiseModel`` and
# ``GPSolver`` (``inference.md`` §10a, fold-in 7) and
# ``Likelihood.capability_parts`` puts them on the composed problem, so a
# fixture named ``"mirror"`` that composed the core classes — which now
# declare ``"reference"`` — would be a two-backend problem and refused. Like
# the instrument steps above, these change the declaration and nothing else:
# the arithmetic is still the reference path's, honestly, because a fixture
# that exists to prove the parametrisation has nothing different to offer in
# a Cholesky.


# These four keep the **core classes' own names**, unlike the steps above,
# and that is load-bearing rather than stylistic: ``Likelihood.to_spec``
# records ``type(noise).__name__`` and the solver's ``type(...).__name__``,
# and ``test_cross_backend``'s ``ampere_likelihoods`` row compares that
# declaration across backends. A backend's noise model and solver are the
# same *declaration* as the core's — what differs is which library computes
# it — so they must present the same name. The real backends do this
# naturally; here the core classes are imported under private aliases so
# that the shadowing is deliberate and visible.


class IndependentNoise(_CoreIndependentNoise):
    """Uncorrelated noise, declared as this backend's (W2.13)."""

    BACKEND: ClassVar[str] = BACKEND


class GaussianProcessNoise(_CoreGaussianProcessNoise):
    """The GP noise composition, declared as this backend's (W2.13)."""

    BACKEND: ClassVar[str] = BACKEND


class DenseGP(_CoreDenseGP):
    """The dense solver, declared as this backend's (W2.13)."""

    BACKEND: ClassVar[str] = BACKEND


class QuasisepGP(_CoreQuasisepGP):
    """The quasiseparable solver, declared as this backend's (W2.13)."""

    BACKEND: ClassVar[str] = BACKEND


class HilbertSpaceGP(_CoreHilbertSpaceGP):
    """The reduced-rank spectral solver, declared as this backend's (W5.4)."""

    BACKEND: ClassVar[str] = BACKEND


# The kernels, for the same reason one step further down. **W3.8** (ruled by
# Peter 2026-09-08) put the kernel on ``Likelihood.capability_parts`` too: a
# solver builds the covariance by calling ``Kernel.matrix``, so a core kernel
# inside a native solver hands back numpy and the conversion detaches the
# graph. A fixture named ``"mirror"`` composing the core kernels is therefore
# a two-backend problem now, and refused — which is this widening working,
# not a casualty of it. They keep the core classes' names for the reason the
# four above do: ``KernelSpec`` and ``Likelihood.to_spec`` record the
# declaration, and ``test_cross_backend`` compares it across backends.


class Matern12(_CoreMatern12):
    """Matérn-1/2, declared as this backend's (W3.8, family added W4.5)."""

    BACKEND: ClassVar[str] = BACKEND


class Matern32(_CoreMatern32):
    """Matérn-3/2, declared as this backend's (W3.8)."""

    BACKEND: ClassVar[str] = BACKEND


class Matern52(_CoreMatern52):
    """Matérn-5/2, declared as this backend's (W3.8, family added W4.5)."""

    BACKEND: ClassVar[str] = BACKEND


class SHO(_CoreSHO):
    """The damped oscillator, declared as this backend's (W4.5)."""

    BACKEND: ClassVar[str] = BACKEND


class RotationTerm(_CoreRotationTerm):
    """The rotation pair, declared as this backend's (W4.5)."""

    BACKEND: ClassVar[str] = BACKEND


class SquaredExponential(_CoreSquaredExponential):
    """Squared exponential, declared as this backend's (W3.8)."""

    BACKEND: ClassVar[str] = BACKEND


_MODELS = {ModelKind.LINEAR: MirrorLinearModel, ModelKind.POWER_LAW: MirrorPowerLawModel}
_KERNELS: dict[KernelFamily, type[Kernel]] = {
    KernelFamily.MATERN12: Matern12,
    KernelFamily.MATERN32: Matern32,
    KernelFamily.MATERN52: Matern52,
    KernelFamily.SHO: SHO,
    KernelFamily.ROTATION: RotationTerm,
    KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential,
}


class MirrorBackend(ReferenceBackend):
    """The second fixture. Same contracts, different arithmetic."""

    name = BACKEND
    capabilities = BackendCapabilities(
        differentiable=False,
        batchable=False,
        device="cpu",
        float64=True,
        solvers=ReferenceBackend.capabilities.solvers,
        tolerances=ReferenceBackend.capabilities.tolerances,
        interferometry=True,
        astrometry=True,
        image=True,
    )

    def model(self, spec: ModelSpec) -> Model:
        return _MODELS[spec.kind](spec)

    def transformation(self, spec: TransformationSpec) -> Transformation:
        if spec.kind is TransformationKind.SCALE:
            return MirrorCalibrationScale(st.lognorm(0.2), label=spec.label)
        if spec.kind is TransformationKind.REBIN:
            return MirrorResample(spec.target, label=spec.label)
        return MirrorPhotometry(spec.target, spec.filters, label=spec.label)

    def kernel(self, spec: CovarianceSpec) -> Kernel:
        return build_kernel(spec, _KERNELS)

    def gp_solver(
        self,
        kind: SolverKind,
        *,
        basis_size: int | Sequence[int] = 32,
        boundary_factor: float = 2.0,
    ) -> GPSolver:
        if kind is SolverKind.DENSE:
            return DenseGP()
        if kind is SolverKind.HILBERT:
            return HilbertSpaceGP(basis_size=basis_size, boundary_factor=boundary_factor)
        return QuasisepGP()

    def independent_noise(self) -> NoiseModel:
        return IndependentNoise()

    def gp_noise(self, kernel: Kernel, solver: GPSolver, *, jitter: Any = None) -> NoiseModel:
        return GaussianProcessNoise(kernel, solver, jitter=jitter)

    def interferometry(self) -> InterferometryPieces:
        return MIRROR_INTERFEROMETRY

    def image(self) -> ImagePieces:
        return MIRROR_IMAGE

    def astrometry(self) -> AstrometryPieces:
        return MIRROR_ASTROMETRY

    def parameter_space(self, declaration: ParameterSet) -> MirrorParameterSpace:
        return MirrorParameterSpace(declaration)
