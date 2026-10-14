"""``examples/line_fluxes``: the composition, the twins, and the study's claims (W7.9).

What is here:

* the composition builds — one model channel, one ``LineIntegration`` built
  ``from_observed`` the synthetic table, the table aligned by transition label
  from the reverse order it is drawn in;
* the GP arm composes with ``DenseGP`` and agrees with ``QuasisepGP``;
* the per-backend model twins agree with the reference model on the same
  grid, and compose on torch and jax;
* **the study's claims**, at the per-PR emcee budget on the reference path
  (quoted in the report: about thirty seconds for the three arms): the
  flexible arm's ``temperature`` and ``column`` within
  :data:`~examples.line_fluxes.study.FLEXIBLE_MAX_BIAS_WIDTHS` of the cool
  component's truth, the misspecified independent arm beyond
  :data:`~examples.line_fluxes.study.MISSPECIFIED_MIN_BIAS_WIDTHS` on at least
  one of them, and the GP localisation peaking at one of the
  high-excitation (short-wavelength) lines;
* NUTS on torch and jax, marked ``heavy``;
* ``main`` prints.
"""

from __future__ import annotations

import importlib.util
import warnings
from typing import Any

import numpy as np
import pytest

from ampere.core import LineFluxes
from examples.line_fluxes import generators, study

needs_torch = pytest.mark.skipif(
    importlib.util.find_spec("torch") is None, reason="needs ampere[torch]"
)
needs_jax = pytest.mark.skipif(importlib.util.find_spec("jax") is None, reason="needs ampere[jax]")

#: The three highest-excitation transitions: where the missing hot component
#: leaves its largest residual, so where the localisation must peak.
HIGH_EXCITATION = {generators.LABELS[i] for i in np.argsort(generators.ENERGY)[-3:]}


def _prior_median(problem: Any) -> np.ndarray:
    return problem.prior_transform(np.full(problem.free_size, 0.5))


class TestTheComposition:
    def test_one_channel_one_instrument_aligned_by_transition(self) -> None:
        problem = study.build_problem("reference", "misspecified")
        assert problem.backend == "reference"
        assert set(problem.datasets) == {"lines"}
        observed = problem.datasets["lines"].observed
        assert isinstance(observed, LineFluxes)
        # Drawn in reverse label order; aligned onto the step's, which
        # from_observed took from the table itself — so the order is the table's.
        assert list(observed.transitions) == list(generators.LABELS[::-1])
        requirements = problem.requirements["model"]["lines"]
        axis = requirements.axes["spectral_axis"]
        assert axis.points is not None
        assert axis.points.size == 64 * len(generators.LABELS)

    def test_the_truth_reproduces_the_noiseless_table(self) -> None:
        """The model's rendered lines integrate back to their fluxes through the step."""
        model = study.truth_model()
        instrument = study.instrument_for("reference")
        from ampere.core import negotiate

        compiled = model.compile_for(negotiate([instrument]))
        got = np.asarray(instrument(compiled.evaluate()).values, dtype=float)
        context = model.portable_context({})
        expected = np.asarray(model.line_fluxes(context), dtype=float) * 1.0e-16
        # Six dispersions each side: the Gaussian's tails beyond the window are
        # 2e-9 of the line; the trapezoid rule on 64 points is far tighter.
        np.testing.assert_allclose(got, expected, rtol=1e-6)

    @pytest.mark.parametrize("arm", study.ARMS)
    def test_each_arm_evaluates(self, arm: str) -> None:
        problem = study.build_problem("reference", arm)
        assert np.isfinite(problem.log_prob(_prior_median(problem)))

    def test_the_dense_and_quasiseparable_flexible_arms_agree(self) -> None:
        from ampere.core import QuasisepGP

        dense = study.build_problem("reference", "flexible")
        quasisep = study.build_problem("reference", "flexible", solver=QuasisepGP())
        theta = _prior_median(dense)
        assert dense.log_prob(theta) == pytest.approx(quasisep.log_prob(theta), abs=1e-6)


class TestTheTwins:
    @needs_jax
    def test_the_jax_model_is_the_reference_model(self) -> None:
        reference = study.truth_model("reference")
        twin = study.truth_model("jax")
        grid = np.linspace(12.78, 12.85, 41)
        context = reference.portable_context({})
        np.testing.assert_allclose(
            np.asarray(twin._flux(twin.ops.asarray(grid), twin.portable_context({}))),
            reference._flux(grid, context),
            rtol=1e-12,
        )

    @needs_jax
    def test_the_flexible_arm_composes_on_jax(self) -> None:
        problem = study.build_problem("jax", "flexible")
        assert problem.backend == "jax"
        assert np.isfinite(problem.log_prob(_prior_median(problem)))

    @needs_torch
    def test_the_flexible_arm_composes_on_torch(self) -> None:
        problem = study.build_problem("torch", "flexible")
        assert problem.backend == "torch"
        assert np.isfinite(problem.log_prob(_prior_median(problem)))


@pytest.fixture(scope="module")
def fitted() -> dict[str, tuple[Any, Any]]:
    """The three arms at the per-PR budget, at the study's seed: one fit each."""
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for arm in study.ARMS:
            problem = study.build_problem("reference", arm)
            out[arm] = (problem, study.run(problem, study.CI_BUDGET))
    return out


class TestTheStudysClaims:
    def test_the_flexible_arm_recovers_the_cool_component(self, fitted: Any) -> None:
        widths = study.bias_widths(fitted["flexible"][1])
        assert max(widths.values()) < study.FLEXIBLE_MAX_BIAS_WIDTHS, widths

    def test_the_misspecified_arm_is_confidently_wrong(self, fitted: Any) -> None:
        widths = study.bias_widths(fitted["misspecified"][1])
        assert max(widths.values()) > study.MISSPECIFIED_MIN_BIAS_WIDTHS, widths
        flexible = study.bias_widths(fitted["flexible"][1])
        assert max(widths.values()) > max(flexible.values())

    def test_the_control_recovers_the_truth(self, fitted: Any) -> None:
        widths = study.bias_widths(fitted["correct"][1])
        assert max(widths.values()) < study.FLEXIBLE_MAX_BIAS_WIDTHS, widths

    def test_the_localisation_peaks_at_a_high_excitation_line(self, fitted: Any) -> None:
        problem, tree = fitted["flexible"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            peak = study.localisation_peak(tree, problem)
        assert peak in HIGH_EXCITATION, (peak, HIGH_EXCITATION)


class TestNUTS:
    @needs_jax
    @pytest.mark.heavy
    def test_the_flexible_arm_fits_under_nuts_on_jax(self) -> None:
        tree = study.run(study.build_problem("jax", "flexible"))
        assert max(study.bias_widths(tree).values()) < study.FLEXIBLE_MAX_BIAS_WIDTHS

    @needs_torch
    @pytest.mark.heavy
    def test_the_flexible_arm_fits_under_nuts_on_torch(self) -> None:
        tree = study.run(study.build_problem("torch", "flexible"))
        assert max(study.bias_widths(tree).values()) < study.FLEXIBLE_MAX_BIAS_WIDTHS


def test_main_prints_the_three_arms(capsys: pytest.CaptureFixture[str], monkeypatch: Any) -> None:
    monkeypatch.setattr(study, "CI_BUDGET", study.EmceeBudget(walkers=8, steps=30, burn_in=10))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert study.main(["--quick"]) == 0
    printed = capsys.readouterr().out
    for arm in study.ARMS:
        assert arm in printed
    assert "GP localisation peaks at" in printed
