"""W4.5's fringing demonstration: composable, exact, and **informational**.

M2's pinned assertions are elsewhere and are untouched. What is asserted here
is what the demonstration must *be able to do* for its numbers to mean
anything, which is a different and much smaller claim than what its numbers
are:

* ``Matern32 + SHO`` composes into an M2 problem and reaches the **O(N)**
  solver, at rank 4 rather than rank 2 — which is what makes the comparison
  affordable at the ladder's larger sizes;
* the composed kernel's O(N) score is the dense one's, exactly, as the
  stationary kernel's is;
* the sum's hyperparameters are namespaced and reach the posterior under
  names a reader can use.

The comparison itself — bias and coverage beside the stationary fit — is the
last test in this module, marked ``m2_full`` and asserting **nothing about
which fit wins**. It is a demonstration, and the reason it asserts nothing is
not timidity: the flexible likelihood with a stationary Matérn was already
calibrated on this scenario (that is M2's result), so a sum kernel cannot
improve the calibration, only the *shape* of the account it gives of the
residual. Pinning a bias threshold would be inventing a prediction nobody has
made. The row exists so that the comparison stays runnable and its output stays
readable; the numbers belong in a report, not in an assertion.
"""

from __future__ import annotations

import numpy as np
import pytest

from ampere.core import QuasisepGP, Sum
from examples.m2_misspecification import fringing, study
from examples.m2_misspecification.generators import generate


def _values(amplitude: float, fringe: float, period: float) -> dict[str, float]:
    """A parameter mapping for the composed problem, at the truth plus a GP state."""
    prefix = f"{study.DATASET_LABEL}.likelihood"
    return {
        "model.A": 1.0,
        "model.B": 3.0,
        "model.d1": 0.15,
        "model.d2": 0.10,
        f"{prefix}.smooth.amplitude": amplitude,
        f"{prefix}.smooth.length_scale": 0.002,
        f"{prefix}.fringe.amplitude": fringe,
        f"{prefix}.fringe.period": period,
    }


class TestTheComposedKernelComposes:
    def test_it_is_a_sum_of_a_matern_and_an_oscillator(self) -> None:
        kernel = fringing.build_fringing_kernel()
        assert isinstance(kernel, Sum)
        assert [child.FAMILY for _, child in kernel.terms] == ["matern32", "sho"]
        assert [label for label, _ in kernel.terms] == ["smooth", "fringe"]

    def test_its_hyperparameters_are_namespaced_by_term(self) -> None:
        """Two amplitudes in one kernel, and neither shadows the other."""
        kernel = fringing.build_fringing_kernel()
        assert kernel.parameters.free_names == (
            "smooth.amplitude",
            "smooth.length_scale",
            "fringe.amplitude",
            "fringe.period",
        )

    def test_it_reaches_the_quasiseparable_solver(self) -> None:
        """A sum of quasiseparable terms is quasiseparable — W4.5's structural claim.

        Composition is the test: ``QuasisepGP.check_compatible`` refuses a
        kernel it cannot lower, term by term, so a problem that builds at all
        is a problem whose GP solve is O(N).
        """
        data = generate("strong_smooth", size=200)
        problem = fringing.build_fringing_problem(data)
        noise = problem.datasets["default"].likelihood.noise
        assert isinstance(noise.solver, QuasisepGP)
        assert problem.free_size == 8

    def test_the_representation_is_rank_four(self) -> None:
        """Rank 2 for the Matérn plus rank 2 for the oscillator: the cost, stated."""
        from ampere.core.kernels import lookup_quasiseparable_term

        kernel = fringing.build_fringing_kernel()
        builder = lookup_quasiseparable_term(kernel.FAMILY)
        axis = np.linspace(0.842, 0.872, 32)
        resolved = kernel.resolve(
            {
                "smooth.amplitude": 0.05,
                "smooth.length_scale": 0.002,
                "fringe.amplitude": 0.05,
                "fringe.period": fringing.FRINGE_PERIOD,
            }
        )
        assert builder(kernel, resolved, axis).rank == 4


class TestTheTwoSolversAgreeOnIt:
    """The same equality ``test_solvers.py`` holds the stationary kernel to."""

    @pytest.mark.parametrize("size", [200, 2_000])
    def test_the_o_n_score_is_the_dense_one(self, size: int) -> None:
        from ampere.core import Dataset, DenseGP, FittingProblem, GaussianFamily, Likelihood
        from ampere.core import GaussianProcessNoise as CoreGPNoise

        data = generate("strong_smooth", size=size)
        quasisep = fringing.build_fringing_problem(data)
        dense_likelihood = Likelihood(
            GaussianFamily(),
            CoreGPNoise(fringing.build_fringing_kernel(), DenseGP()),
        )
        dense = FittingProblem(
            study.model_for("reference", data.wavelength),
            [
                Dataset(
                    data.container(),
                    likelihood=dense_likelihood,
                    label=study.DATASET_LABEL,
                )
            ],
            seed=study.SEED,
        )
        values = _values(0.05, 0.05, fringing.FRINGE_PERIOD)
        assert float(quasisep.log_prob(values)) == pytest.approx(
            float(dense.log_prob(values)), rel=1e-10
        )

    def test_they_agree_away_from_one_point(self) -> None:
        """One point could agree by accident; a scan of the fringe period cannot."""
        from ampere.core import Dataset, DenseGP, FittingProblem, GaussianFamily, Likelihood
        from ampere.core import GaussianProcessNoise as CoreGPNoise

        data = generate("strong_smooth", size=200)
        quasisep = fringing.build_fringing_problem(data)
        dense = FittingProblem(
            study.model_for("reference", data.wavelength),
            [
                Dataset(
                    data.container(),
                    likelihood=Likelihood(
                        GaussianFamily(),
                        CoreGPNoise(fringing.build_fringing_kernel(), DenseGP()),
                    ),
                    label=study.DATASET_LABEL,
                )
            ],
            seed=study.SEED,
        )
        for period in (0.0018, 0.0028, 0.0042):
            values = _values(0.03, 0.06, period)
            assert float(quasisep.log_prob(values)) == pytest.approx(
                float(dense.log_prob(values)), rel=1e-10
            )


def _period_margin(scenario: str, factor: float, size: int = 400) -> float:
    """``log_prob`` at the generating fringe period minus ``log_prob`` at *factor* times it."""
    problem = fringing.build_fringing_problem(generate(scenario, size=size))
    at_period = float(problem.log_prob(_values(0.02, 0.06, fringing.FRINGE_PERIOD)))
    off_period = float(problem.log_prob(_values(0.02, 0.06, factor * fringing.FRINGE_PERIOD)))
    return at_period - off_period


class TestTheComposedKernelSeesTheRipple:
    """The one claim about the *mathematics* that is safe to pin.

    The demonstration is worth running only if the SHO term is actually
    responding to the ripple rather than to the parameterisation. The
    statement that settles it is about the likelihood **surface**, not about a
    posterior, so it needs no sampler and has no Monte Carlo error: the
    *margin* by which the generating fringe period is preferred over a wrong
    one must grow with how much fringing is present.

    The margin rather than the raw preference, because the raw preference is
    confounded. A slow oscillator is a smoother function than a fast one, so
    at a fixed amplitude the longer period is mildly preferred even on the
    control spectrum, which has no ripple at all — measured, at 400 points:
    ``none`` prefers a 2.5x-long period by 29 nats. That is a statement about
    the oscillator's roughness and not about fringing. The *difference between
    scenarios* has that term in common and cancels it.
    """

    def test_the_likelihood_prefers_the_generating_period_when_there_is_a_fringe(
        self,
    ) -> None:
        assert _period_margin("strong_smooth", 0.4) > 0.0
        assert _period_margin("strong_smooth", 2.5) > 0.0

    @pytest.mark.parametrize("factor", [0.4, 2.5], ids=["shorter", "longer"])
    def test_the_margin_grows_with_the_amount_of_fringing(self, factor: float) -> None:
        """none (0 %) < mild (2.5 %) < strong (7 %): the ordering of the deviations."""
        margins = [_period_margin(key, factor) for key in ("none", "mild", "strong_smooth")]
        assert margins[0] < margins[1] < margins[2]


@pytest.mark.m2_full
def test_the_demonstration_runs_and_reports(capsys: pytest.CaptureFixture[str]) -> None:
    """**Informational.** Three fits of one spectrum; nothing is asserted about which wins.

    What is asserted is that the comparison runs, that every fit produced a
    posterior for all four physical parameters, and that the table names the
    three likelihoods — so the demonstration cannot rot silently. The numbers
    are printed for a reader and recorded in W4.5's report.
    """
    comparison = fringing.compare("strong_smooth", size=200)
    print(comparison.table())
    for which in ("standard", "stationary", "composed"):
        summaries = getattr(comparison, which)
        assert set(summaries) == set(study.PHYSICAL_NAMES)
        assert all(summary.n_draws > 0 for summary in summaries.values())
    rendered = capsys.readouterr().out
    assert "Matern32 + SHO" in rendered
    assert "worst bias" in rendered
