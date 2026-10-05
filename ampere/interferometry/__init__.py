"""Interferometry in ampere: the observable's front door (W6.12).

One import for everything a fit to interferometric data needs on the reference
(numpy) path: the two container kinds, the likelihood families they are scored
with, the reference backend's steps and source models, and the reader that
fills the containers from a file::

    from ampere.interferometry import (
        read_oifits, FourierSample, ClosurePhase, SquaredAmplitude, Binary,
    )

**What this package is.** The placement memo's option D
(``docs/design/phase4_placement_memo.md`` §2.2 to 2.4), reserved at Phase 4 and
created here with the first reader: a user-facing home *layered over* the
placement that the contracts fix, not a place code moves to. The kinds stay in
:mod:`ampere.core`, where every backend reads them; the steps and models stay
in each backend's ``interferometry`` module, under the backend rule (one name
per backend, a native problem composed entirely from one backend's pieces).
What is defined here and nowhere else is the reader, :func:`read_oifits`
(:mod:`ampere.interferometry.oifits`).

**Torch and jax users** compose from their own backend's pieces, as before:
``ampere.backends.torch.FourierSample``, ``ampere.backends.jax.ClosurePhase``
and so on. Nothing native is re-exported here, because a reference step in a
torch problem is refused at composition, and a front door that offered one
would be offering a refusal. The containers :func:`read_oifits` returns are
backend-neutral and feed any backend's ``FourierSample.from_observed``
unchanged.

**The shape the next observables follow.** ``ampere.astrometry`` and
``ampere.image`` are not created yet; when their first readers land they take
this package's shape: an ``__init__`` that re-exports the kinds and families
from :mod:`ampere.core` and the reference steps and models from
:mod:`ampere.backends.reference` by name, one module per file format holding
its reader, and the page that documents the modality gaining a "front door"
section. ``import ampere`` imports none of them; each is imported on its own.
"""

from __future__ import annotations

from ampere.backends.reference import (
    Amplitude,
    BandwidthSmearing,
    Binary,
    BinaryVisibilities,
    ClosurePhase,
    FourierSample,
    GaussianSource,
    GaussianSourceVisibilities,
    SquaredAmplitude,
    TimeSmearing,
    UniformDisc,
    UniformDiscVisibilities,
)
from ampere.core import (
    ClosurePhases,
    ComplexGaussianFamily,
    GaussianFamily,
    RiceFamily,
    VisibilitySet,
    VonMisesFamily,
)

from .oifits import OIFITSData, OIFITSError, read_oifits

__all__ = [
    "Amplitude",
    "BandwidthSmearing",
    "Binary",
    "BinaryVisibilities",
    "ClosurePhase",
    "ClosurePhases",
    "ComplexGaussianFamily",
    "FourierSample",
    "GaussianFamily",
    "GaussianSource",
    "GaussianSourceVisibilities",
    "OIFITSData",
    "OIFITSError",
    "RiceFamily",
    "SquaredAmplitude",
    "TimeSmearing",
    "UniformDisc",
    "UniformDiscVisibilities",
    "VisibilitySet",
    "VonMisesFamily",
    "read_oifits",
]
