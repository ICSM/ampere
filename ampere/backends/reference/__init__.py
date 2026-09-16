"""The reference backend: pure numpy/scipy, and the base install's whole toolkit.

``architecture.md`` §1 puts this at rung 1 of the capability ladder and §2
settles what it is for. Three things, and it is worth being explicit that the
third is not an afterthought:

1. **The conformance oracle.** The battery in ``tests/conformance/`` computes
   reference values here and holds every other backend to them. That is what
   stops two lockstep backends drifting (``DEVELOPMENT_PLAN.md`` §7).
2. **A useful base install.** ampere's base install -- ``pip install .`` from
   a checkout (PyPI's ``ampere`` package is unrelated) -- with no extras at
   all, is a complete numpy-only fitting environment: these models and
   instrument steps, ``ampere.core``'s flexible GP likelihood, and the
   gradient-free engines.
3. **The execution venue for adapted models.** A black-box or astropy-adapted
   model is not differentiated by this backend, but the numpy-side
   ``Instrument``/``Likelihood``/``Dataset`` machinery wrapping its raw output
   is this backend's implementation of the §4 contracts.

**Correctness is the only goal** (``architecture.md`` §2). No gradients, no
batching, no GPU, no performance target — those are what the torch and jax
backends are for, and keeping this one small is the reason a third contract
implementation is affordable at all. Everything runs in float64, which is
policy rather than a default (``architecture.md`` §5): GP solves in float32
fail in ways that read as science bugs.

Nothing here imports an optional dependency. ``pyphot`` — a base dependency —
is imported lazily inside
:meth:`~ampere.backends.reference.SyntheticPhotometry.from_library`, because a
chain with no photometry in it should not pay to open an HDF5 filter library.

What is shipped
---------------

Models (``DEVELOPMENT_PLAN.md`` §5's Phase 2 list): :class:`BlackBody`,
:class:`ModifiedBlackBody`, :class:`PowerLaw`.

Instrument steps (``transformations.md`` §10's table):
:class:`CalibrationScale`, :class:`Resample`, :class:`LSFConvolution`,
:class:`SyntheticPhotometry`.

Interferometry (Phase 4, ``interferometry.py``): the steps
:class:`FourierSample`, :class:`ClosurePhase`, :class:`BandwidthSmearing`,
:class:`TimeSmearing` and :class:`Amplitude`, the image-emitting source models
:class:`UniformDisc`, :class:`GaussianSource` and :class:`Binary`, and the same
three emitting visibilities analytically — :class:`UniformDiscVisibilities`,
:class:`GaussianSourceVisibilities`, :class:`BinaryVisibilities` — which are
the closed forms the direct transform is held to. The two container kinds they
speak in, ``VisibilitySet`` and ``ClosurePhases``, are ``ampere.core``'s.

Images (Phase 5, ``image.py``, W5.5): the PSF-convolution step
:class:`PSFConvolution`, ``Image`` -> ``Image`` — ``transformations.md`` §10's
image slot, and the first step whose observed container is a
``Layout.GRID`` one. The source models it convolves are ``interferometry.py``'s
three, reused unchanged.

Astrometry (Phase 4, ``astrometry.py``, W4.9): the epoch-sampling step
:class:`EpochSample` and the reflex-orbit model :class:`ReflexOrbit`, which
emits two :class:`~ampere.core.TimeSeries` channels (``"ra"``, ``"dec"``) —
the second modality built by the interferometry page's own template.

Noise (``likelihoods.md`` §5, X-1): :class:`FractionalModelNoise` and its GP
composition :class:`FractionalModelGPNoise`.

**Every piece here declares ``BACKEND = "reference"``** (W2.12) beside
``DIFFERENTIABLE``, ``BATCHABLE`` and ``DEVICE``, explicitly rather than by
inheriting the ABCs' default. The name is this backend's one name everywhere:
it is what ``ampere.core.declared_capabilities`` aggregates onto a problem, what
a run's ``ampere_backend`` records, the key ``lowering.md`` §12.8's registry is
consulted with, and the conformance fixture's id.

Kernels, GP solvers, families, containers and the fitting problem itself are
**not** here: they are backend-neutral and live in ``ampere.core``, which is
the claim ``inference.md`` §18 makes — "a backend supplies models and
transformations […] and nothing else".
"""

from __future__ import annotations

from .astrometry import (
    EpochSample,
    ReflexOrbit,
)
from .image import (
    PSFConvolution,
)
from .instrument import (
    DETECTORS,
    CalibrationScale,
    LSFConvolution,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    bundled_filter_library,
)
from .interferometry import (
    BRIGHTNESS_UNIT,
    MAS_PER_RAD,
    SPECTRAL_UNIT,
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    ClosurePhase,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    TimeSmearing,
    UniformDisc,
    UniformDiscVisibilities,
    cell_solid_angle,
)
from .models import (
    COORDINATE_UNIT,
    FLUX_UNIT,
    BlackBody,
    ModifiedBlackBody,
    PowerLaw,
    planck_jy,
)
from .noise import FractionalModelGPNoise, FractionalModelNoise

__all__ = [
    "BRIGHTNESS_UNIT",
    "COORDINATE_UNIT",
    "DETECTORS",
    "FLUX_UNIT",
    "MAS_PER_RAD",
    "SPECTRAL_UNIT",
    "Amplitude",
    "BandwidthSmearing",
    "Binary",
    "BinaryVisibilities",
    "BlackBody",
    "CalibrationScale",
    "ClosurePhase",
    "EpochSample",
    "FourierSample",
    "FractionalModelGPNoise",
    "FractionalModelNoise",
    "GaussianSource",
    "GaussianSourceVisibilities",
    "LSFConvolution",
    "ModifiedBlackBody",
    "PSFConvolution",
    "PowerLaw",
    "ReflexOrbit",
    "Resample",
    "SyntheticPhotometry",
    "TimeSmearing",
    "UniformDisc",
    "UniformDiscVisibilities",
    "bin_edges",
    "bundled_filter_library",
    "cell_solid_angle",
    "planck_jy",
]
