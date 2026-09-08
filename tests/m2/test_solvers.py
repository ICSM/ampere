"""The two solvers, the two kernels, and the ladder's cost — without sampling.

Everything here is a handful of ``log_prob`` calls, so it runs in seconds and
belongs in the per-PR gate. It establishes the three facts the size ladder rests
on before any chain is started:

1. :class:`~ampere.core.QuasisepGP` and :class:`~ampere.core.DenseGP` return the
   **same number** for a Matérn-3/2 kernel. That is what makes the O(N) rung of
   the ladder a faster computation of the same likelihood rather than an
   approximation of it — W2.3's whole claim, restated on this study's problem.
2. :class:`~ampere.core.SquaredExponential` — legacy ampere's kernel — is
   **refused** by the O(N) solver, because it is not quasiseparable. That
   asymmetry is the concrete reason W2.10 reproduces the study with Matérn-3/2
   and keeps the RBF only as a 200-point, dense cross-check.
3. The dense solver becomes unaffordable where the ladder says it does.

The wall-clock table itself is ``tests/benchmarks/test_m2_likelihood.py``;
nothing here asserts a time, because a benchmark that fails on a slow runner is
a benchmark that gets deleted.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import DenseGP, QuasisepGP
from ampere.core.exceptions import LikelihoodError
from examples.m2_misspecification import study
from examples.m2_misspecification.generators import generate


@pytest.mark.parametrize("size", [200, 2_000])
def test_the_two_solvers_score_the_same_likelihood(size: int) -> None:
    """Matern-3/2 is *exactly* quasiseparable, so this is equality, not closeness."""
    data = generate("strong_smooth", size=size)
    dense = study.build_problem(data, likelihood="flexible", solver="dense")
    quasisep = study.build_problem(data, likelihood="flexible", solver="quasisep")
    assert float(quasisep.log_prob(quasisep.reference_values)) == pytest.approx(
        float(dense.log_prob(dense.reference_values)), rel=1e-10
    )


def test_the_solvers_agree_away_from_the_reference_point() -> None:
    """One point could agree by accident; a scan of the amplitude cannot."""
    data = generate("strong_sharp", size=200)
    dense = study.build_problem(data, likelihood="flexible", solver="dense")
    quasisep = study.build_problem(data, likelihood="flexible", solver="quasisep")
    assert dense.parameters.free_names == quasisep.parameters.free_names
    for amplitude in (0.01, 0.05, 0.2):
        values = {
            "model.A": 1.0,
            "model.B": 3.0,
            "model.d1": 0.15,
            "model.d2": 0.10,
            f"{study.DATASET_LABEL}.likelihood.amplitude": amplitude,
            f"{study.DATASET_LABEL}.likelihood.length_scale": 0.002,
        }
        assert float(quasisep.log_prob(values)) == pytest.approx(
            float(dense.log_prob(values)), rel=1e-10
        ), amplitude


def test_the_legacy_kernel_is_refused_by_the_on_solver() -> None:
    with pytest.raises(ValueError, match="not quasiseparable"):
        study.build_likelihood("flexible", kernel="squared_exponential", solver="quasisep")


def test_the_legacy_kernel_composes_and_scores_densely() -> None:
    """The RBF at 200 points: what legacy ampere's flexible likelihood used.

    It gives a *different* number from the Matern-3/2 at the same
    hyperparameters, and that is not a defect to be tightened away: a different
    covariance function is a different likelihood, and two kernels agreeing on
    an arbitrary point of their shared parameter space would be the surprise.
    The cross-check that matters is at the level of the *conclusion* rather than
    the number — whether the flexible likelihood's honesty survives the change
    of kernel — and that is
    ``tests/m2/test_science.py::test_the_legacy_kernel_reaches_the_same_conclusion``.
    """
    data = generate("strong_smooth", size=200)
    rbf = study.build_problem(
        data, likelihood="flexible", kernel="squared_exponential", solver="dense"
    )
    matern = study.build_problem(data, likelihood="flexible", solver="dense")
    rbf_value = float(rbf.log_prob(rbf.reference_values))
    matern_value = float(matern.log_prob(matern.reference_values))
    assert np.isfinite(rbf_value) and np.isfinite(matern_value)
    assert rbf_value != matern_value
    assert rbf.free_size == matern.free_size == 6


def test_the_quasiseparable_refusal_is_the_kernels_own() -> None:
    """Composed directly, ``ampere.core`` refuses the same pair — by name."""
    kernel = study.build_kernel(kernel="squared_exponential")
    assert kernel.spec().quasiseparable is False
    with pytest.raises(LikelihoodError):
        QuasisepGP().check_compatible(kernel, generate("none", size=32).container())
    # ... and accepts the one the study uses.
    DenseGP().check_compatible(kernel, generate("none", size=32).container())
    matern = study.build_kernel()
    assert matern.spec().quasiseparable is True
    QuasisepGP().check_compatible(matern, generate("none", size=32).container())


@pytest.mark.parametrize("size", [200, 2_000, 20_000])
def test_the_on_solver_scores_every_rung_of_the_ladder(size: int) -> None:
    """Including the top one, where the dense solver is not attempted at all."""
    data = generate("strong_smooth", size=size)
    problem = study.build_problem(data, likelihood="flexible", solver="quasisep")
    assert problem.free_size == 6
    assert np.isfinite(float(problem.log_prob(problem.reference_values)))


@pytest.mark.parametrize("size", [200, 2_000, 20_000])
def test_the_standard_likelihood_scores_every_rung(size: int) -> None:
    data = generate("strong_smooth", size=size)
    problem = study.build_problem(data, likelihood="standard")
    assert problem.free_size == 4
    assert np.isfinite(float(problem.log_prob(problem.reference_values)))
