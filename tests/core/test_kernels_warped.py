"""`WarpedKernel`: the warp itself, its refusals, its cost, and its diagnostics (W5.7).

``tests/conformance/test_likelihoods.py`` holds the warp to the §6 contract on
every registered backend — the covariance, the O(N) agreement, the bit-identical
identity, the spec hash. This file holds the things that are *about* the
reference implementation and would be noise repeated three times:

* the **parameterisation**'s properties: monotone whatever the knot variables
  are, exact at the identity, linear outside the knot range, and the warp
  composing with the kernel algebra (a warped ``Sum``, a ``Sum`` of warped
  terms);
* every **refusal**, with its message: a non-monotone knot set, a warp with no
  warp in it, a base kernel that selects axes of its own, more than one
  selected coordinate, ``value()`` on a kernel that is not stationary;
* the **cost**, on a wall clock — the conformance battery counts the
  representation's elements, which is the deterministic half of the claim, and
  this is the other half: the warped O(N) path is linear in N and is orders of
  magnitude below the dense one at a size where both still run;
* the **diagnostics**, families B and C, run on ``conditional``'s *warped*
  output with the warp recorded — the degrees-of-freedom guard's last clause,
  and the reason the warp is reported at all.
"""

from __future__ import annotations

import math
import time
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    DenseGP,
    GaussianFamily,
    GaussianProcessNoise,
    HierarchicalPrior,
    Likelihood,
    LikelihoodError,
    Matern12,
    Matern32,
    Parameter,
    QuasisepGP,
    Spectrum,
    Sum,
    VisibilitySet,
    WarpedKernel,
    quantile_knots,
)

GRID = np.linspace(0.0, 10.0, 64)
KNOTS = (0.0, 3.0, 6.5, 10.0)
AMPLITUDE_KNOTS = (0.0, 5.0, 10.0)
INCREMENTS = (0.8, -0.6, 0.35)
LEVELS = (0.3, -0.4, 0.25)


def fixed_warp(
    base: Any = None,
    *,
    input_warp: bool = True,
    amplitude_warp: bool = True,
    increments: Any = None,
    levels: Any = None,
) -> WarpedKernel:
    """The warp this module uses, with every knot variable held at a number."""
    keywords: dict[str, Any] = {}
    if input_warp:
        keywords["input_warp"] = KNOTS
        keywords["increments"] = INCREMENTS if increments is None else increments
        keywords["input_scale"] = 1.0
    if amplitude_warp:
        keywords["amplitude_warp"] = AMPLITUDE_KNOTS
        keywords["levels"] = LEVELS if levels is None else levels
        keywords["amplitude_scale"] = 1.0
    return WarpedKernel(Matern32(0.4, 2.0) if base is None else base, **keywords)


def spectrum(values: np.ndarray, sigma: float = 0.1) -> Spectrum:
    return Spectrum(
        GRID * u.micron,
        values * u.Jy,
        uncertainty=np.full(GRID.size, sigma) * u.Jy,
    )


# ---------------------------------------------------------------------------
# The parameterisation
# ---------------------------------------------------------------------------


class TestTheWarpIsMonotone:
    """Structurally, not by check: every segment slope is a softplus."""

    @pytest.mark.parametrize(
        "increments",
        [
            (0.0, 0.0, 0.0),
            (5.0, 5.0, 5.0),
            (-5.0, -5.0, -5.0),
            (12.0, -12.0, 12.0),
            (20.0, -20.0, 3.0),
            (-20.0, 20.0, -20.0),
        ],
    )
    def test_however_extreme_the_knot_variables(self, increments: tuple[float, ...]) -> None:
        """Strictly increasing across the range a sampler can actually reach."""
        kernel = fixed_warp(amplitude_warp=False, increments=increments)
        fine = np.linspace(-5.0, 15.0, 501)
        warped = np.asarray(kernel.warped_coordinate(fine, kernel.resolve(None)), dtype=float)
        assert np.all(np.isfinite(warped))
        assert np.all(np.diff(warped) > 0.0)

    @pytest.mark.parametrize("increments", [(700.0, -700.0, 0.0), (-1e3, 1e3, -1e3), (-40.0,) * 3])
    def test_absurd_knot_variables_stay_finite_and_ordered(
        self, increments: tuple[float, ...]
    ) -> None:
        """The two float64 edges, both benign, both worth pinning.

        ``+700`` is where a naive ``log(1 + exp(u))`` would overflow to ``inf``
        and the warp would stop being a warp; the stable form used here has no
        such point, which is what the finiteness assertion buys.

        Below about ``-37`` the other edge appears: the segment slope
        ``softplus(u)/softplus(0)`` drops under machine epsilon, so in the
        offset form ``x + (m - 1)(x - x_0)`` the ``m`` is lost to rounding and
        that segment becomes flat to within a ULP rather than merely very flat.
        The warp is then non-decreasing *up to rounding* rather than strictly
        increasing. That is a property of float64 and not a defect: a slope of
        1e-17 compresses the whole segment into less than one ULP of the
        coordinate, so the covariance it describes **is** the constant block a
        flat segment gives. The shrinkage prior keeps a sampler nowhere near
        here; the row records what happens if something else puts it there.
        """
        kernel = fixed_warp(amplitude_warp=False, increments=increments)
        fine = np.linspace(-5.0, 15.0, 501)
        warped = np.asarray(kernel.warped_coordinate(fine, kernel.resolve(None)), dtype=float)
        assert np.all(np.isfinite(warped))
        rounding = 8.0 * float(np.spacing(float(np.max(np.abs(warped)))))
        assert np.all(np.diff(warped) >= -rounding)

    def test_the_warp_is_linear_outside_the_knot_range(self) -> None:
        """Extension with the end segments' own slopes, so ``w`` is defined everywhere."""
        kernel = fixed_warp(amplitude_warp=False)
        below = np.linspace(-20.0, -10.0, 9)
        above = np.linspace(20.0, 30.0, 9)
        for outside in (below, above):
            warped = np.asarray(
                kernel.warped_coordinate(outside, kernel.resolve(None)), dtype=float
            )
            second = np.diff(warped, n=2)
            assert np.allclose(second, 0.0, atol=1e-12)

    def test_the_identity_warp_is_the_identity_bit_for_bit(self) -> None:
        kernel = fixed_warp(increments=0.0, levels=0.0)
        values = kernel.resolve(None)
        fine = np.linspace(-5.0, 15.0, 501)
        assert np.array_equal(np.asarray(kernel.warped_coordinate(fine, values), dtype=float), fine)
        plain = Matern32(0.4, 2.0)
        assert np.array_equal(
            kernel.matrix(fine[:, None], fine[:, None], values),
            plain.matrix(fine[:, None], fine[:, None], plain.resolve(None)),
        )

    def test_the_amplitude_warp_interpolates_its_knots_exactly(self) -> None:
        """``a`` at a knot is ``exp`` of that knot's level: the interpolant passes through."""
        kernel = fixed_warp(input_warp=False)
        record = kernel.warp_provenance(kernel.resolve(None))
        assert record is not None
        assert record["amplitude_warp_values"] == pytest.approx(
            np.exp(np.asarray(LEVELS)), abs=1e-12
        )

    def test_the_input_warp_fixes_its_first_knot(self) -> None:
        """``w(x_0) = x_0``: the anchor that leaves no global translation free."""
        kernel = fixed_warp(amplitude_warp=False)
        record = kernel.warp_provenance(kernel.resolve(None))
        assert record is not None
        assert record["input_warp_images"][0] == pytest.approx(KNOTS[0], abs=0.0)
        assert record["input_warp_knots"] == list(KNOTS)


class TestTheWarpComposesWithTheAlgebra:
    """A warp wraps anything, and anything wraps a warp — both stay quasiseparable."""

    def test_a_warped_sum_reaches_the_o_n_path(self) -> None:
        base = Sum(Matern32(0.3, 2.0), Matern12(0.15, 0.4))
        kernel = fixed_warp(base, amplitude_warp=False)
        assert kernel.QUASISEPARABLE
        observed = spectrum(1.0 + 0.3 * np.sin(GRID))
        predicted = observed.with_values(np.ones(GRID.size))
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
        quasi = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        assert quasi.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=1e-9
        )

    def test_a_sum_of_warped_terms_is_refused_on_the_o_n_path(self) -> None:
        """The other order is not representable, and says so by naming the term.

        Two terms warped differently have two coordinates, and the recursion's
        propagators come from one axis; a sum of their generators would be a
        matrix that is neither term's covariance and not obviously wrong
        either. ``DenseGP`` evaluates each term on its own coordinate and is
        unaffected, which is what the refusal offers.
        """
        kernel = Sum(
            fixed_warp(Matern32(0.3, 2.0), amplitude_warp=False),
            fixed_warp(Matern12(0.15, 0.4), input_warp=False),
        )
        observed = spectrum(1.0 + 0.3 * np.sin(GRID))
        predicted = observed.with_values(np.ones(GRID.size))
        noise = GaussianProcessNoise(kernel, QuasisepGP())
        with pytest.raises(LikelihoodError, match="term 'term0' is a WarpedKernel"):
            noise.check_compatible(GaussianFamily(), observed)
        # The dense path is genuinely fine: a sum of matrices needs no single
        # coordinate at all.
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
        assert np.isfinite(dense.log_prob(predicted, observed))

    def test_the_refusal_names_a_term_however_deep(self) -> None:
        nested = Sum(Matern32(0.3, 2.0), Sum(Matern12(0.1, 0.5), fixed_warp(Matern12(0.15, 0.4))))
        observed = spectrum(1.0 + 0.3 * np.sin(GRID))
        noise = GaussianProcessNoise(nested, QuasisepGP())
        with pytest.raises(LikelihoodError, match=r"term 'term1\.term1' is a WarpedKernel"):
            noise.check_compatible(GaussianFamily(), observed)

    def test_a_warp_of_a_warp_is_a_warp(self) -> None:
        """Composition of monotone maps is monotone, and the registry recurses."""
        inner = fixed_warp(amplitude_warp=False)
        outer = WarpedKernel(
            inner, input_warp=(0.0, 5.0, 10.0), increments=(0.4, -0.2), input_scale=1.0
        )
        values = outer.resolve(None)
        warped = np.asarray(outer.warped_coordinate(GRID, values), dtype=float)
        assert np.all(np.diff(warped) > 0.0)
        observed = spectrum(1.0 + 0.3 * np.sin(GRID))
        predicted = observed.with_values(np.ones(GRID.size))
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(outer, DenseGP()))
        quasi = Likelihood(GaussianFamily(), GaussianProcessNoise(outer, QuasisepGP()))
        assert quasi.log_prob(predicted, observed) == pytest.approx(
            dense.log_prob(predicted, observed), abs=1e-9
        )

    def test_it_warps_one_named_axis_of_a_multi_axis_container(self) -> None:
        """``axes=`` picks the coordinate; the base kernel inherits the selection."""
        rng = np.random.default_rng(5)
        size = 24
        observed = VisibilitySet(
            rng.normal(0.0, 10.0, size) * u.dimensionless_unscaled,
            rng.normal(0.0, 10.0, size) * u.dimensionless_unscaled,
            np.linspace(2.0, 2.4, size) * u.micron,
            (rng.normal(0.0, 0.2, size) + 1j * rng.normal(0.0, 0.2, size)) * u.Jy,
            uncertainty=np.full(size, 0.05) * u.Jy,
        )
        kernel = WarpedKernel(
            Matern32(0.1, 0.05),
            input_warp=(2.0, 2.2, 2.4),
            increments=(0.5, -0.3),
            input_scale=1.0,
            axes=("spectral_axis",),
        )
        DenseGP().check_compatible(kernel, observed)
        assert kernel.selected_axes([axis.name for axis in observed.axes]) == ("spectral_axis",)

        # End to end through the binding: the covariance a solver would build
        # on this container's own (n, 3) coordinate block must be the base
        # kernel on the warped *spectral* column and nothing else. The base
        # stays unbound and is handed the one column, which is why its own
        # ``axes`` is recorded (for the unit rule) but never resolved.
        bound = kernel.for_axes([axis.name for axis in observed.axes])
        points = np.column_stack([np.asarray(axis.values, dtype=float) for axis in observed.axes])
        values = bound.resolve(None)
        got = np.asarray(bound.matrix(points, points, values))
        moved = np.asarray(bound.warped_coordinate(points[:, 2], values), dtype=float)
        plain = Matern32(0.1, 0.05)
        expected = plain.matrix(moved[:, None], moved[:, None], plain.resolve(None))
        assert got == pytest.approx(expected, abs=0.0, rel=0.0)


# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------


class TestRefusals:
    """Each one names what to do instead, as every refusal in this contract does."""

    def test_a_non_monotone_knot_set(self) -> None:
        with pytest.raises(LikelihoodError, match="strictly increasing"):
            WarpedKernel(Matern32(0.4, 2.0), input_warp=(0.0, 6.0, 3.0, 10.0))

    def test_a_repeated_knot(self) -> None:
        with pytest.raises(LikelihoodError, match="strictly increasing"):
            WarpedKernel(Matern32(0.4, 2.0), amplitude_warp=(0.0, 5.0, 5.0))

    def test_a_single_knot(self) -> None:
        with pytest.raises(LikelihoodError, match="at least two knot locations"):
            WarpedKernel(Matern32(0.4, 2.0), input_warp=(1.0,))

    def test_a_warp_with_no_warp_in_it(self) -> None:
        with pytest.raises(LikelihoodError, match="neither input_warp"):
            WarpedKernel(Matern32(0.4, 2.0))

    def test_a_base_that_selects_axes_of_its_own(self) -> None:
        with pytest.raises(LikelihoodError, match="declares axes= of its own"):
            WarpedKernel(Matern32(0.4, 2.0, axes=("u",)), input_warp=KNOTS)

    def test_more_than_one_selected_coordinate(self) -> None:
        rng = np.random.default_rng(3)
        size = 12
        observed = VisibilitySet(
            rng.normal(0.0, 10.0, size) * u.dimensionless_unscaled,
            rng.normal(0.0, 10.0, size) * u.dimensionless_unscaled,
            np.linspace(2.0, 2.4, size) * u.micron,
            (rng.normal(0.0, 0.2, size) + 1j * rng.normal(0.0, 0.2, size)) * u.Jy,
            uncertainty=np.full(size, 0.05) * u.Jy,
        )
        kernel = WarpedKernel(Matern32(0.1, 5.0), input_warp=KNOTS, axes=("u", "v"))
        with pytest.raises(LikelihoodError, match="monotone map of one ordered coordinate"):
            DenseGP().check_compatible(kernel, observed)

    def test_the_closed_form_in_one_separation(self) -> None:
        kernel = fixed_warp()
        with pytest.raises(LikelihoodError, match="not stationary"):
            kernel.value(np.array([0.0, 1.0]), kernel.resolve(None))

    def test_a_wrong_number_of_knot_declarations(self) -> None:
        with pytest.raises(LikelihoodError, match="declaration\\(s\\) for 3 knot variable"):
            WarpedKernel(Matern32(0.4, 2.0), input_warp=KNOTS, increments=(0.1, 0.2))

    def test_a_single_parameter_broadcast_across_knots(self) -> None:
        with pytest.raises(LikelihoodError, match="cannot broadcast a single Parameter"):
            WarpedKernel(
                Matern32(0.4, 2.0),
                input_warp=KNOTS,
                increments=Parameter("increment0", value=0.0, fixed=True),
            )

    def test_a_non_kernel_base(self) -> None:
        with pytest.raises(LikelihoodError, match="warps a Kernel"):
            WarpedKernel("matern32", input_warp=KNOTS)  # type: ignore[arg-type]


class TestQuantileKnots:
    """The convenience that keeps knot placement a user act, recorded in the declaration."""

    def test_it_places_them_at_even_quantiles(self) -> None:
        assert quantile_knots(np.arange(11.0), 3) == (0.0, 5.0, 10.0)
        assert quantile_knots(np.arange(11.0), 2) == (0.0, 10.0)

    def test_it_follows_the_data_rather_than_the_range(self) -> None:
        """An even grid would spend its degrees of freedom on an empty tail."""
        crowded = np.concatenate([np.linspace(0.0, 1.0, 99), np.array([100.0])])
        middle = quantile_knots(crowded, 3)[1]
        assert middle < 1.0

    def test_it_refuses_fewer_than_two(self) -> None:
        with pytest.raises(LikelihoodError, match="at least two knots"):
            quantile_knots(np.arange(10.0), 1)

    def test_it_refuses_an_empty_coordinate(self) -> None:
        with pytest.raises(LikelihoodError, match="empty coordinate set"):
            quantile_knots(np.array([]), 3)


class TestTheShrinkageDeclaration:
    """The degrees-of-freedom guard, as a declaration rather than as advice."""

    def test_the_default_is_non_centred_shrinkage_about_the_identity(self) -> None:
        kernel = WarpedKernel(Matern32(0.4, 2.0), input_warp=KNOTS)
        parameters = kernel.parameters
        scale = parameters["input_warp.scale"]
        assert scale.prior is not None
        assert float(scale.prior.mean()) > 0.0
        for index in range(len(KNOTS) - 1):
            knot = parameters[f"input_warp.increment{index}"]
            assert knot.prior is not None
            # Centred on the identity warp, so leaving it costs prior mass.
            assert float(knot.prior.mean()) == pytest.approx(0.0, abs=1e-12)
            assert float(knot.prior.std()) == pytest.approx(1.0, abs=1e-12)

    def test_the_centred_form_declares_a_hierarchical_prior(self) -> None:
        kernel = WarpedKernel(
            Matern32(0.4, 2.0),
            input_warp=KNOTS,
            increments=HierarchicalPrior("norm", {"scale": "input_warp.scale"}, kwds={"loc": 0.0}),
            non_centred=False,
        )
        prior = kernel.parameters["input_warp.increment0"].prior
        assert isinstance(prior, HierarchicalPrior)
        assert prior.references == ("input_warp.scale",)

    def test_the_centred_form_defaults_to_the_hierarchical_prior(self) -> None:
        """``non_centred=False`` alone is already the shrinkage declaration."""
        kernel = WarpedKernel(Matern32(0.4, 2.0), input_warp=KNOTS, non_centred=False)
        assert isinstance(kernel.parameters["input_warp.increment1"].prior, HierarchicalPrior)

    def test_the_centred_form_refuses_a_scale_nothing_depends_on(self) -> None:
        """A flat prior on the knots with a free scale leaves a sampled dimension idle.

        In the non-centred form the kernel *uses* the scale, so it is always
        identified; in the centred form the only route from the scale to the
        knots is a ``HierarchicalPrior`` naming it. Overriding the knots'
        declaration with an ordinary prior and leaving the scale free gives a
        parameter with no posterior, which costs a sampler real work for
        nothing — so it is refused where it is declared, rather than found in
        a trace plot afterwards.
        """
        with pytest.raises(LikelihoodError, match="nothing depends on"):
            WarpedKernel(
                Matern32(0.4, 2.0),
                input_warp=KNOTS,
                increments=st.norm(0.0, 1.0),
                non_centred=False,
            )
        # Holding the scale at a number is the other way out, and is accepted.
        held = WarpedKernel(
            Matern32(0.4, 2.0),
            input_warp=KNOTS,
            increments=st.norm(0.0, 1.0),
            input_scale=0.5,
            non_centred=False,
        )
        assert "input_warp.scale" not in held.parameters.free_names

    def test_the_two_forms_are_the_same_prior(self) -> None:
        """``u = s·z`` with ``z ~ N(0, 1)`` is ``u ~ N(0, s)``: the same model, better geometry."""
        scale = 0.37
        non_centred = WarpedKernel(
            Matern32(0.4, 2.0),
            input_warp=KNOTS,
            increments=1.25,
            input_scale=scale,
            non_centred=True,
        )
        centred = WarpedKernel(
            Matern32(0.4, 2.0),
            input_warp=KNOTS,
            increments=1.25 * scale,
            input_scale=scale,
            non_centred=False,
        )
        assert np.allclose(
            np.asarray(non_centred.warped_coordinate(GRID, non_centred.resolve(None))),
            np.asarray(centred.warped_coordinate(GRID, centred.resolve(None))),
            rtol=0.0,
            atol=1e-14,
        )

    def test_the_two_forms_declare_different_specs(self) -> None:
        """Different parameters, so a chain stored under one cannot be read as the other."""
        non_centred = WarpedKernel(Matern32(0.4, 2.0), input_warp=KNOTS, non_centred=True)
        centred = WarpedKernel(
            Matern32(0.4, 2.0),
            input_warp=KNOTS,
            increments=HierarchicalPrior("norm", {"scale": "input_warp.scale"}, kwds={"loc": 0.0}),
            non_centred=False,
        )
        assert non_centred.spec().to_dict() != centred.spec().to_dict()


# ---------------------------------------------------------------------------
# Cost
# ---------------------------------------------------------------------------


def _fastest(function: Any, repeats: int = 5) -> float:
    """The best of *repeats* wall-clock timings, in seconds.

    The minimum rather than the mean: on a shared machine the noise is
    one-sided (something else took the core), so the fastest run is the closest
    thing to the cost of the work itself.
    """
    best = math.inf
    for _ in range(repeats):
        started = time.perf_counter()
        function()
        best = min(best, time.perf_counter() - started)
    return best


class TestTheWarpedSolveIsLinear:
    """The other half of the O(N) claim: a wall clock, not an assertion.

    The conformance battery counts the representation's elements, which is
    deterministic and says the *structure* is rank-J semiseparable. This says
    the structure is actually being exploited.
    """

    def _solve(self, size: int) -> Any:
        rng = np.random.default_rng(size)
        axis = np.sort(rng.uniform(0.0, 10.0, size))
        residual = rng.normal(0.0, 0.1, size)
        variance = np.full(size, 0.01)
        kernel = fixed_warp()
        values = kernel.resolve(None)
        solver = QuasisepGP()
        return lambda: solver.log_marginal_likelihood(
            kernel, axis[:, None], residual, variance, values
        )

    def test_doubling_the_data_roughly_doubles_the_cost(self) -> None:
        small = _fastest(self._solve(4_000))
        medium = _fastest(self._solve(8_000))
        large = _fastest(self._solve(16_000))
        # Linear would be 4x across the whole range; quadratic 16x, cubic 64x.
        # The bound is generous because a wall clock on a shared machine is,
        # and it still separates linear from quadratic by a factor of two.
        assert large / small < 8.0
        assert medium / small < 4.5

    def test_it_is_far_below_the_dense_solve_at_a_size_both_can_run(self) -> None:
        size = 2_000
        rng = np.random.default_rng(11)
        axis = np.sort(rng.uniform(0.0, 10.0, size))
        residual = rng.normal(0.0, 0.1, size)
        variance = np.full(size, 0.01)
        kernel = fixed_warp()
        values = kernel.resolve(None)
        quasisep = QuasisepGP()
        dense = DenseGP()
        quick = _fastest(
            lambda: quasisep.log_marginal_likelihood(
                kernel, axis[:, None], residual, variance, values
            )
        )
        slow = _fastest(
            lambda: dense.log_marginal_likelihood(
                kernel, axis[:, None], residual, variance, values
            ),
            repeats=1,
        )
        assert quick * 20.0 < slow
        # and they agree, which is what makes the speed worth anything
        assert quasisep.log_marginal_likelihood(
            kernel, axis[:, None], residual, variance, values
        ) == pytest.approx(
            dense.log_marginal_likelihood(kernel, axis[:, None], residual, variance, values),
            abs=1e-6,
        )


# ---------------------------------------------------------------------------
# Diagnostics: families B and C on the warped residuals
# ---------------------------------------------------------------------------


#: Where the "model" is deliberately wrong, in the raw coordinate.
DEFICIENT_BAND = (6.0, 8.0)


def _deficient_pair(kernel: WarpedKernel) -> tuple[Spectrum, Spectrum]:
    """A prediction that misses a bump between 6 and 8 micron, and the data."""
    rng = np.random.default_rng(20260916)
    sigma = 0.05
    inside = (GRID >= DEFICIENT_BAND[0]) & (GRID <= DEFICIENT_BAND[1])
    truth = np.ones(GRID.size)
    deficiency = np.where(inside, 0.45 * np.sin(np.pi * (GRID - 6.0) / 2.0), 0.0)
    observed = spectrum(truth + deficiency + rng.normal(0.0, sigma, GRID.size), sigma=sigma)
    return observed.with_values(truth), observed


class TestTheDiagnosticsRunOnWarpedResiduals:
    """The guard's last clause: the whiteness and localisation tests see ``w(x)``.

    A warp is flexible enough to absorb signal, so the tests that would catch
    it doing so must be run in the coordinate the residuals are *claimed* to be
    stationary in. ``Likelihood.conditional`` therefore reports there, and
    records the warp so a reader can get back to the observed axis.
    """

    def test_the_conditional_reports_in_warped_coordinates_with_the_warp_recorded(
        self,
    ) -> None:
        kernel = fixed_warp()
        predicted, observed = _deficient_pair(kernel)
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        conditional = likelihood.conditional(predicted, observed)

        assert conditional.coordinates is not None
        expected = np.asarray(kernel.warped_coordinate(GRID, kernel.resolve(None)), dtype=float)
        assert np.array_equal(conditional.coordinates, expected)
        # Not the observed axis: if it were, "run the diagnostics on the warped
        # residuals" would be a sentence with no consequence.
        assert not np.allclose(conditional.coordinates, GRID)

        assert conditional.warp is not None
        assert conditional.warp["kind"] == "warped"
        assert conditional.warp["base"] == "matern32"
        assert conditional.warp["input_warp_knots"] == list(KNOTS)
        assert conditional.warp["amplitude_warp_knots"] == list(AMPLITUDE_KNOTS)
        assert len(conditional.warp["input_warp_images"]) == len(KNOTS)

    def test_an_unwarped_kernel_reports_nothing_new(self) -> None:
        """Every declaration written before W5.7 is untouched."""
        kernel = Matern32(0.4, 2.0)
        predicted, observed = _deficient_pair(fixed_warp())
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        conditional = likelihood.conditional(predicted, observed)
        assert conditional.coordinates is None
        assert conditional.warp is None

    def test_family_b_whiteness_is_computed_in_the_warped_coordinate(self) -> None:
        """The separation bins are separations in ``w(x)``, which is the point.

        The statistic itself is ``ampere.results``'s; what W5.7 supplies is the
        axis it is handed. Running it on the raw axis and on the warped one
        gives genuinely different bin centres, so the choice is consequential
        rather than cosmetic.
        """
        from ampere.results.diagnostics import separation_binned_autocorrelation

        kernel = fixed_warp()
        predicted, observed = _deficient_pair(kernel)
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        conditional = likelihood.conditional(predicted, observed)
        assert conditional.coordinates is not None

        sigma = np.asarray(observed.uncertainty, dtype=float)
        standardised = (
            np.asarray(observed.values, dtype=float) - np.asarray(predicted.values, dtype=float)
        ) / sigma
        warped_bins, warped_rho, warped_pairs = separation_binned_autocorrelation(
            conditional.coordinates, standardised, bins=6
        )
        raw_bins, raw_rho, _ = separation_binned_autocorrelation(GRID, standardised, bins=6)

        assert np.all(np.isfinite(warped_rho))
        assert np.all(warped_pairs > 0)
        assert not np.allclose(warped_bins, raw_bins)
        # The residual is not white — that is what the injected bump is — and
        # the test says so in either coordinate; the warped one is the honest
        # one under a warped kernel.
        assert float(np.max(np.abs(warped_rho))) > 0.2
        assert float(np.max(np.abs(raw_rho))) > 0.2

    def test_family_c_localises_the_deficiency_in_the_warped_coordinate(self) -> None:
        """The conditioned mean, standardised, peaks where the model is wrong.

        ``gp_localisation_score``'s own quantity — ``|mean| / sd``, "how many
        sigma of GP the fit needed here" — computed directly from the
        conditional rather than from a stored run, because what is under test
        is the conditional's *indexing*, not ``ampere.results``'s aggregation.
        """
        kernel = fixed_warp()
        predicted, observed = _deficient_pair(kernel)
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        conditional = likelihood.conditional(predicted, observed)
        assert conditional.coordinates is not None

        score = np.abs(conditional.mean) / conditional.standard_deviation
        inside = (GRID >= DEFICIENT_BAND[0]) & (GRID <= DEFICIENT_BAND[1])
        assert float(np.max(score[inside])) > float(np.max(score[~inside]))

        # The score is indexed by the warped coordinate, which is monotone in
        # the raw one — so "where" is answerable in both, and the warp record
        # is what lets a reader map back.
        warped = conditional.coordinates
        assert np.all(np.diff(warped) > 0.0)
        peak = int(np.argmax(score))
        assert DEFICIENT_BAND[0] <= GRID[peak] <= DEFICIENT_BAND[1]


class TestWarpProvenanceIsPlainData:
    """It travels as an attribute on a diagnostics group, so it must be JSON-plain."""

    def test_every_value_is_a_plain_python_type(self) -> None:
        import json

        kernel = fixed_warp()
        record = kernel.warp_provenance(kernel.resolve(None))
        assert record is not None
        restored = json.loads(json.dumps(record))
        assert restored == record

    def test_an_unwarped_kernel_has_none(self) -> None:
        plain = Matern32(0.4, 2.0)
        assert plain.warp_provenance(plain.resolve(None)) is None
        assert np.array_equal(np.asarray(plain.warped_coordinate(GRID, plain.resolve(None))), GRID)


class TestTheDeclarationIsRecorded:
    """The knots are part of the model, so they are part of the spec."""

    def test_the_spec_carries_the_knots_and_the_parameterisation(self) -> None:
        described = fixed_warp().spec().to_dict()
        assert described["family"] == "warped"
        assert described["metadata"]["input_warp_knots"] == list(KNOTS)
        assert described["metadata"]["amplitude_warp_knots"] == list(AMPLITUDE_KNOTS)
        assert described["metadata"]["non_centred"] is True
        assert described["terms"][0]["label"] == "base"
        assert described["terms"][0]["kernel"]["family"] == "matern32"

    def test_an_unwarped_spec_has_no_metadata_key_at_all(self) -> None:
        """Which is what leaves every spec hash minted before W5.7 unchanged."""
        assert "metadata" not in Matern32(0.4, 2.0).spec().to_dict()
        assert "metadata" not in Sum(Matern32(0.4, 2.0), Matern12(0.1, 0.3)).spec().to_dict()

    def test_a_declared_prior_survives_into_the_parameters(self) -> None:
        kernel = WarpedKernel(
            Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
            input_warp=KNOTS,
            increments=st.norm(0.0, 0.5),
        )
        assert set(kernel.parameters.free_names) == {
            "base.amplitude",
            "base.length_scale",
            "input_warp.scale",
            "input_warp.increment0",
            "input_warp.increment1",
            "input_warp.increment2",
        }
