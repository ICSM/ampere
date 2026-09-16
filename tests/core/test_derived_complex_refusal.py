"""W5.2 item 5: the three on-demand derived groups refuse a complex dataset alike.

Before this, ``ampere.results.derived.add_posterior_predictive`` refused a
complex-valued dataset explicitly (``ResultsError``, "complex-valued"), while
``add_residuals`` and ``gp_localisation`` let plain numpy casting (a complex
array assigned into a ``float`` array) silently discard the imaginary part,
with only a ``ComplexWarning`` — the kind of warning a script that redirects
stderr never sees. ``ampere.results.derived._refuse_complex`` is now the one
place all three call, so this pins that the three raised messages share one
shape (only the group name and the "what to do instead" sentence differ).

**Kept in ``tests/core/`` rather than ``tests/results/`` deliberately.**
W5.2's dispatch owns ``ampere/results/derived.py`` but not
``tests/results/`` — another wave (W5.0) touches files under
``ampere/results/`` in parallel, and ``tests/results/`` is its territory,
not this item's. ``ampere.results.DrawRecorder`` and the derived-group
functions are the ordinary public API either way, so nothing here reaches
into a private module to get around that boundary; it is simply filed next
to the item's other ``ampere/core``-adjacent changes instead of beside the
other ``tests/results`` suites this item did not touch.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianProcessNoise,
    Likelihood,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    VisibilitySet,
)
from ampere.core.exceptions import ResultsError
from ampere.results import DrawRecorder, add_posterior_predictive, add_residuals, gp_localisation

pytest.importorskip("arviz", reason="ampere.results needs arviz")

AXES = (np.array([1.0, 2.0, -3.0]), np.array([3.0, -4.0, 1.5]), np.array([2.2, 2.2, 2.2]))


class PointSource(Model):
    """A flat complex visibility -- the smallest complex-valued channel there is."""

    def __init__(self, u_axis: np.ndarray, v_axis: np.ndarray, wave: np.ndarray) -> None:
        self.register_buffer("u", np.asarray(u_axis, dtype=float))
        self.register_buffer("v", np.asarray(v_axis, dtype=float))
        self.register_buffer("wave", np.asarray(wave, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("flux", st.loguniform(0.5, 2.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            VisibilitySet(
                ctx["u"],
                ctx["v"],
                ctx["wave"] * u.micron,
                ctx["flux"] * np.ones(ctx["u"].size, dtype=complex),
            )
        )


def _vis_dataset(*, gp: bool) -> Dataset:
    observed = VisibilitySet(
        AXES[0],
        AXES[1],
        AXES[2] * u.micron,
        np.array([1 + 0j, 1 + 0j, 1 + 0j]),
        uncertainty=np.array([0.1, 0.1, 0.1]),
    )
    noise = (
        GaussianProcessNoise(Matern32(0.3, 1.0, axes=("u", "v")), solver=DenseGP()) if gp else None
    )
    return Dataset(observed, label="vis", likelihood=Likelihood(ComplexGaussianFamily(), noise))


def visibility_problem(*, gp: bool = False) -> FittingProblem:
    return FittingProblem(
        PointSource(*AXES),
        DatasetCollection({"vis": _vis_dataset(gp=gp)}),
        seed=1,
    )


def _run(problem: FittingProblem, *, draws: int = 3) -> Any:
    """A tiny genuine run, through the driver-facing recorder (as W5.2's other tests do)."""
    recorder = DrawRecorder(problem)
    for _ in range(draws):
        recorder.record({"model.flux": 1.0})
    return recorder.emit(engine="fixture")


def _message_shape(group: str) -> str:
    """The prefix common to all three refusals; only the group name varies here."""
    return (
        f"dataset 'vis' is complex-valued, and a {group} group holds one real variable "
        f"per dataset. results.md §4 splits a complex observed container into <label>_real "
        f"and <label>_imag; the derived groups do not do that yet, and which component -- "
        f"or the modulus -- is the right one to derive is the caller's choice (deferred to "
        f"W5.3)."
    )


class TestTheThreeGroupsRefuseComplexTheSameWay:
    def test_add_posterior_predictive_refuses_by_the_shared_shape(self) -> None:
        problem = visibility_problem()
        tree = _run(problem)
        with pytest.raises(ResultsError) as excinfo:
            add_posterior_predictive(tree, problem)
        assert str(excinfo.value).startswith(_message_shape("posterior-predictive"))

    def test_add_residuals_refuses_by_the_shared_shape(self) -> None:
        """Before W5.2 this cast complex to real with only a ComplexWarning."""
        problem = visibility_problem()
        tree = _run(problem)
        with pytest.raises(ResultsError) as excinfo:
            add_residuals(tree, problem)
        assert str(excinfo.value).startswith(_message_shape("residuals"))

    def test_gp_localisation_refuses_by_the_shared_shape(self) -> None:
        """Before W5.2 this cast the complex conditioned mean to real with only a warning."""
        problem = visibility_problem(gp=True)
        tree = _run(problem)
        with pytest.raises(ResultsError) as excinfo:
            gp_localisation(tree, problem)
        assert str(excinfo.value).startswith(_message_shape("gp_localisation"))
