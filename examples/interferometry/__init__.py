"""Phase 4's proof modality, as a study: a resolved binary with a fainter disc.

``DEVELOPMENT_PLAN.md`` §5 Phase 4, and ``docs/source/interferometry.rst``'s
worked case. A synthetic resolved binary with a fainter, more extended
circular-Gaussian disc is observed as both visibilities and closure phases,
and fitted three ways — with the disc in the model, with it omitted, and
with it omitted but the visibilities carrying W4.2's flexible likelihood —
asking M2's question of the modality Phase 4 built: does the flexible
likelihood keep the binary's parameters calibrated when the sky model is
wrong? A fourth, independent scenario (the item's "chromatic case",
``docs/design/phase4_placement_memo.md`` §3.6) asks whether a kernel needs to
see wavelength as its own axis at all, by measuring three candidate kernels
against a residual that is sharp in wavelength and smooth in ``(u, v)``
rather than arguing it.

The modules
-----------
:mod:`~examples.interferometry.generators`
    The truth (a binary plus a disc; a second, chromatic truth for the
    fourth scenario), reusing ``tests/backends/interferometry_fixtures.py``'s
    array geometry rather than re-deriving it.
:mod:`~examples.interferometry.model`
    :class:`~examples.interferometry.model.BinaryWithDisc`, the "correct"
    fitted model — a small, hand-written composition of two existing
    ``ampere`` models, in the shape a user's own model takes.
    :mod:`.model_torch` and :mod:`.model_jax` bind it to the two modern
    backends; imported only when asked for.
:mod:`~examples.interferometry.study`
    Composition, running (a single fit, or simulation-based calibration over
    many), and the reduction of a stored run to the numbers
    ``tests/interferometry`` asserts on.
:mod:`~examples.interferometry.figures`
    The six ``ampere.results`` plots plus the calibration figures, drawn from
    stored runs.

Running it
----------
``python -m examples.interferometry`` runs the three-arm study and the
chromatic scenario on the reference backend and prints the recovery table;
see ``--help`` for the arm, the backend and the budget. Nothing it writes is
committed: figures go to a directory the caller names.

Neither this package nor ``tests/interferometry`` adds a dependency beyond
what Phase 4 already shipped.
"""

from __future__ import annotations

from .generators import BINARY, DISC_FLUX, DISC_FWHM, SEED, TRUTH, synthetic
from .model import BinaryWithDisc, binary_with_disc
from .study import (
    ARMS,
    BACKENDS,
    CI_BUDGET,
    DOC_BUDGET,
    Summary,
    build_problem,
    coverage_at,
    model_for,
    run,
    run_calibration,
    run_study,
    summarise,
)

__all__ = [
    "ARMS",
    "BACKENDS",
    "BINARY",
    "CI_BUDGET",
    "DISC_FLUX",
    "DISC_FWHM",
    "DOC_BUDGET",
    "SEED",
    "TRUTH",
    "BinaryWithDisc",
    "Summary",
    "binary_with_disc",
    "build_problem",
    "coverage_at",
    "model_for",
    "run",
    "run_calibration",
    "run_study",
    "summarise",
    "synthetic",
]
