"""Conformance rows for the optimisers (W6.7, ``inference.md`` §10b).

The fixture registry fits this directly: one column per available backend,
each composing the same declaration from its own pieces, and no row names a
backend. Two claims:

* **The objective convention is one convention.** On every backend with a
  differentiable realisation, ``optimise(method="map")`` — the gradient route
  through ``ampere.core.realise``, the Jacobian term subtracted on the realised
  side — and ``optimise(method="scipy")`` — the gradient-free route on the
  numpy contract path — find the same constrained-space mode to ``1e-3`` in
  every free parameter's unconstrained coordinate (the coordinate both
  optimisers move, and the one whose scale does not depend on a parameter's
  units). A backend without a realisation skips, as ``TestTheRealisation``'s
  rows do.
* **The scipy route runs on every column**, and its mode is a stationary
  point of the numpy-path objective there.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import Tie, registered_realisations
from ampere.inference import optimise
from ampere.inference._optimise import constrained_objective

from .composition import GP_GRID, DatasetSpec, ProblemSpec, build_problem
from .protocol import (
    ConformanceBackend,
    ModelKind,
    ModelSpec,
    TransformationKind,
    TransformationSpec,
)

CALIBRATION = TransformationSpec(TransformationKind.SCALE, label="calibration")

#: ``test_inference.py``'s two problem shapes, restated (a test module is not
#: importable from another): one dataset with a calibration nuisance, and two
#: channels sharing a tied calibration.
SINGLE = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(DatasetSpec(instrument=(CALIBRATION,)),),
)
JOINT = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, channels=("blue", "red"), coordinates=GP_GRID),
    datasets=(
        DatasetSpec(label="blue", channel="blue", instrument=(CALIBRATION,)),
        DatasetSpec(label="red", channel="red", instrument=(CALIBRATION,), data_seed=771),
    ),
    ties=(
        Tie(
            "calibration",
            ("blue.instrument.calibration.scale", "red.instrument.calibration.scale"),
        ),
    ),
)
SHAPES = {"single": SINGLE, "joint": JOINT}

#: The agreement the item's ruling asks between the two routes' modes.
AGREEMENT = 1e-3


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_scipy_and_map_find_the_same_mode(backend: ConformanceBackend, shape: str) -> None:
    problem = build_problem(backend, SHAPES[shape])
    if problem.backend not in registered_realisations() or not problem.differentiable:
        pytest.skip(f"the {problem.backend!r} backend has no differentiable realisation")
    native = optimise(problem, method="map", starts=2)
    scipy = optimise(build_problem(backend, SHAPES[shape]), method="scipy", starts=2)
    assert native.free_labels == scipy.free_labels
    np.testing.assert_allclose(native.unconstrained, scipy.unconstrained, atol=AGREEMENT, rtol=0)
    assert native.log_prob_constrained == pytest.approx(scipy.log_prob_constrained, abs=1e-6)


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_the_scipy_mode_is_stationary(backend: ConformanceBackend, shape: str) -> None:
    problem = build_problem(backend, SHAPES[shape])
    optimum = optimise(problem, method="scipy", starts=2)
    objective = constrained_objective(problem)
    mode = optimum.unconstrained
    for i in range(mode.size):
        step = np.zeros(mode.size)
        step[i] = 1e-3
        assert objective(mode) >= objective(mode + step) - 1e-7
        assert objective(mode) >= objective(mode - step) - 1e-7
