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

from .instrument import (
    DETECTORS,
    CalibrationScale,
    LSFConvolution,
    Resample,
    SyntheticPhotometry,
    bin_edges,
    bundled_filter_library,
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
    "COORDINATE_UNIT",
    "DETECTORS",
    "FLUX_UNIT",
    "BlackBody",
    "CalibrationScale",
    "FractionalModelGPNoise",
    "FractionalModelNoise",
    "LSFConvolution",
    "ModifiedBlackBody",
    "PowerLaw",
    "Resample",
    "SyntheticPhotometry",
    "bin_edges",
    "bundled_filter_library",
    "planck_jy",
]
