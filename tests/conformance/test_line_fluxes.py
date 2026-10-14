"""Line fluxes: the rows the modality template owes a new observable kind (W7.9).

One column per registered backend, like every other file here, and no test
body names one. A backend that has not written the line-integration step
declares ``BackendCapabilities.line_fluxes = False`` and every row skips with
a reason naming what is owed — ``test_astrometry.py``'s shape.

What these rows are for, following ruling 5 of the item's dispatch:

* (a) **the step against a closed form** — the battery's ``POWER_LAW`` model
  integrated over each window, ``∫ norm (λ/λ₀)^index c/λ² dλ``, which is
  analytic; the quadrature's error is the tabulation's, so the bound is the
  trapezoid rule's own, derived from ``samples_per_window`` in the row, and a
  denser tabulation tightens it;
* (b) **the linear continuum is exact on a linear spectrum** — ``LINEAR``
  with ``continuum="linear"`` integrates to zero;
* (c) the published requirement is ``points=`` at the tabulation, and the
  Gap-1 placement: a resampling instrument beside the step on the same
  channel widens the negotiated grid and the fluxes do not change;
* (d) a masked model sample inside a window masks that line and no other;
* (e) the composition: an observed ``LineFluxes`` in another order and in
  ``erg s⁻¹ cm⁻²`` (converted with ``.to_unit``) aligns by transition label,
  and the two refusals name the transitions;
* (f) ``simulate`` draws a ``LineFluxes``;
* (g) the flexible likelihood on the kind — ``DenseGP`` and ``QuasisepGP``
  agree on a ``LineFluxes`` container whose axis is in no particular order;
* (h) the four diagnostics' plots run on a stored run of the composition.
"""

from __future__ import annotations

import warnings
from typing import Any

import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    AxisRequirement,
    Dataset,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    LineFluxes,
    Spectrum,
    negotiate,
)
from ampere.core.exceptions import LikelihoodError, TransformationError

from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    LineFluxPieces,
    ModelKind,
    ModelSpec,
    SolverKind,
    Tolerances,
    TransformationKind,
    TransformationSpec,
)

#: The transitions every row integrates — a JWST-like trio, given out of
#: wavelength order on purpose (the kind tolerates any order, and so must the
#: step's output).
TRANSITIONS = ["[Ne II] 12.81", "[S IV] 10.51", "[Ne III] 15.56"]
#: Their rest wavelengths, micron, in :data:`TRANSITIONS`' order.
REST = np.array([12.8135, 10.5105, 15.5551])
#: The power law's reference wavelength, micron (``ModelSpec.reference_coordinate``).
REFERENCE = 10.0
#: The truth the rows are evaluated at.
TRUTH = {"model.norm": 2.0, "model.index": -1.3}
#: The energy flux unit the step emits.
W_M2 = u.W / u.m**2
#: Fractional uncertainty on an observed line flux.
FRACTION = 0.05


def pieces_or_skip(backend: ConformanceBackend) -> LineFluxPieces:
    """This backend's line-flux classes, or a skip naming what is owed."""
    if not backend.capabilities.line_fluxes:
        pytest.skip(
            f"{backend.name} declares no line-flux vocabulary "
            f"(BackendCapabilities.line_fluxes), so the line-integration step is not there. "
            f"W7.9 owes the native twins."
        )
    return backend.line_fluxes()


def model(backend: ConformanceBackend, kind: ModelKind = ModelKind.POWER_LAW) -> Any:
    """The battery's model of *kind*, on a coarse default grid the step's negotiation replaces."""
    return backend.model(
        ModelSpec(
            kind=kind,
            coordinates=(9.0, 11.0, 13.0, 15.0, 17.0),
            reference_coordinate=REFERENCE,
            compiled=True,
        )
    )


def step(pieces: LineFluxPieces, **overrides: Any) -> Any:
    """The step over :data:`TRANSITIONS`, built through the constructor."""
    arguments: dict[str, Any] = {"window": 300.0, "window_unit": "km/s", **overrides}
    return pieces.line_integration(TRANSITIONS, REST, **arguments)


def predict(backend: ConformanceBackend, instrument: Instrument, values: dict[str, float]) -> Any:
    """The chain's prediction at *values* (bare local names), negotiated for *instrument*."""
    compiled = model(backend).compile_for(negotiate([instrument]))
    return instrument(compiled.evaluate(**values))


def power_law_integral(lower: np.ndarray, upper: np.ndarray, norm: float, index: float) -> Any:
    """``∫ norm (λ/λ₀)^index · c/λ² dλ`` in W m⁻², the closed form (``index ≠ 1``)."""
    c = float(const.c.to_value(u.um / u.s))
    power = index - 1.0
    jansky = 1.0e-26  # W m^-2 Hz^-1
    return norm * c * REFERENCE ** (-index) * (upper**power - lower**power) / power * jansky


def trapezoid_bound(
    lower: np.ndarray, upper: np.ndarray, count: int, norm: float, index: float
) -> Any:
    """The trapezoid rule's error bound, ``(b - a) h² / 12 · max|g''|``.

    ``g(λ) = A λ^(index-2)``, so ``g'' = A (index-2)(index-3) λ^(index-4)``,
    largest in magnitude at the window's short end for the negative powers
    used here. ``h = (b - a) / (count - 1)``: the step's evenly spaced
    tabulation, edges included.
    """
    c = float(const.c.to_value(u.um / u.s))
    amplitude = norm * c * REFERENCE ** (-index) * 1.0e-26
    second = abs(amplitude * (index - 2.0) * (index - 3.0)) * lower ** (index - 4.0)
    width = upper - lower
    spacing = width / (count - 1)
    return width * spacing**2 / 12.0 * second


# ---------------------------------------------------------------------------
# (a) The step against the closed form
# ---------------------------------------------------------------------------


class TestTheStepAgainstTheClosedForm:
    """``POWER_LAW`` with ``continuum="none"``: the analytic integral, to the quadrature's bound."""

    HALF_WIDTH = 0.3  # micron: wide enough that the quadrature error is measurable

    def _errors(self, backend: ConformanceBackend, count: int) -> tuple[Any, Any]:
        pieces = pieces_or_skip(backend)
        line_step = step(pieces, window=self.HALF_WIDTH, window_unit="um", samples_per_window=count)
        instrument = Instrument([line_step], channel="default", label="lines")
        predicted = predict(backend, instrument, {"norm": 2.0, "index": -1.3})
        assert isinstance(predicted, LineFluxes)
        assert predicted.unit == W_M2
        assert list(predicted.transitions) == TRANSITIONS
        lower, upper = REST - self.HALF_WIDTH, REST + self.HALF_WIDTH
        exact = power_law_integral(lower, upper, 2.0, -1.3)
        got = np.asarray(backend.to_numpy(predicted.values), dtype=float)
        bound = trapezoid_bound(lower, upper, count, 2.0, -1.3)
        return np.abs(got - exact), bound

    def test_the_integral_is_the_closed_form_to_the_trapezoid_bound(
        self, backend: ConformanceBackend
    ) -> None:
        """64 samples per 0.6 µm window: the bound is ~1e-6 of the flux, and holds."""
        error, bound = self._errors(backend, 64)
        exact = power_law_integral(REST - self.HALF_WIDTH, REST + self.HALF_WIDTH, 2.0, -1.3)
        assert np.all(bound / exact < 5.0e-6)
        assert np.all(error <= bound), (error, bound)

    def test_a_denser_tabulation_tightens_the_agreement(self, backend: ConformanceBackend) -> None:
        """Four times the samples: the error falls by more than ten (the rule is second order)."""
        sparse, _ = self._errors(backend, 64)
        dense, bound = self._errors(backend, 256)
        assert np.all(dense <= bound)
        assert np.all(dense < sparse / 10.0), (sparse, dense)


# ---------------------------------------------------------------------------
# (b) The linear continuum
# ---------------------------------------------------------------------------


class TestTheLinearContinuum:
    def test_a_linear_spectrum_integrates_to_zero(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The edge-to-edge line *is* a linear spectrum, so nothing is left of it."""
        pieces = pieces_or_skip(backend)
        values = {"offset": 0.4, "slope": 0.25}
        flat = model(backend, ModelKind.LINEAR)
        subtracted = Instrument([step(pieces, continuum="linear")], channel="default", label="s")
        raw = Instrument([step(pieces)], channel="default", label="r")
        with_continuum = subtracted(
            flat.compile_for(negotiate([subtracted])).evaluate(**values)
        ).values
        without = raw(
            model(backend, ModelKind.LINEAR).compile_for(negotiate([raw])).evaluate(**values)
        )
        scale = np.abs(np.asarray(backend.to_numpy(without.values), dtype=float))
        residual = np.abs(np.asarray(backend.to_numpy(with_continuum), dtype=float))
        assert np.all(scale > 0.0)
        assert np.all(residual <= tolerances.exact * scale), residual / scale

    def test_the_declarations_are_described(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        described = step(pieces, continuum="linear", window_unit="um", window=0.1).describe()
        assert described == {"window_unit": "um", "continuum": "linear"}


# ---------------------------------------------------------------------------
# (c) The published requirement, and the Gap-1 placement
# ---------------------------------------------------------------------------


class TestTheRequirementAndThePlacement:
    def test_the_published_requirement_is_points_at_the_tabulation(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        line_step = step(pieces, samples_per_window=17)
        (requirement,) = line_step.requirements()
        assert isinstance(requirement, AxisRequirement)
        assert requirement.axis == "spectral_axis"
        assert np.array_equal(requirement.points, line_step.tabulation())
        assert line_step.tabulation().size == 3 * 17
        assert bool(np.all(np.diff(line_step.tabulation()) > 0.0))

    def test_a_second_instrument_on_the_channel_does_not_move_the_fluxes(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``spectrum_photometry.md`` Gap 1: the columns go where ``Axis.locate`` says."""
        pieces = pieces_or_skip(backend)
        lines = Instrument([step(pieces)], channel="default", label="lines")
        rebin = Instrument(
            [
                backend.transformation(
                    TransformationSpec(
                        kind=TransformationKind.REBIN,
                        target=(9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5, 16.5),
                        label="rebin",
                    )
                )
            ],
            channel="default",
            label="spectrum",
        )
        values = {"norm": 2.0, "index": -1.3}
        alone = predict(backend, lines, values)
        requirements = negotiate([lines, rebin])
        compiled = model(backend).compile_for(requirements)
        result = compiled.evaluate(**values)
        assert result["default"].n_samples > lines.steps[0].tabulation().size
        beside = lines(result)
        np.testing.assert_allclose(
            np.asarray(backend.to_numpy(beside.values), dtype=float),
            np.asarray(backend.to_numpy(alone.values), dtype=float),
            rtol=tolerances.exact,
            atol=0.0,
        )


# ---------------------------------------------------------------------------
# (d) The mask
# ---------------------------------------------------------------------------


class TestTheMask:
    def test_a_masked_sample_inside_a_window_masks_that_line_alone(
        self, backend: ConformanceBackend
    ) -> None:
        """``results_schema.md`` §16's ANY rule: a line depends on its window and nothing else."""
        pieces = pieces_or_skip(backend)
        line_step = step(pieces, samples_per_window=9)
        grid = line_step.tabulation()
        mask = np.zeros(grid.size, dtype=bool)
        # Inside [Ne III]'s window: the third block of the sorted tabulation.
        mask[2 * 9 + 4] = True
        spectrum = Spectrum(grid * u.um, np.ones(grid.size) * u.Jy, mask=mask)
        got = line_step(spectrum, None)
        assert got.mask is not None
        assert got.mask.tolist() == [False, False, True]


# ---------------------------------------------------------------------------
# (e) The composition, aligned by transition
# ---------------------------------------------------------------------------


def observed_table(
    backend: ConformanceBackend, pieces: LineFluxPieces, order: list[int]
) -> LineFluxes:
    """The truth's line fluxes in erg s⁻¹ cm⁻², in *order*, on rounded rest wavelengths.

    What a published table looks like: its own order, its own rounding of the
    rest wavelengths, the cgs unit. The user converts it once with
    ``.to_unit``; alignment by label does the rest.
    """
    instrument = Instrument([step(pieces)], channel="default", label="truth")
    truth = np.asarray(
        backend.to_numpy(predict(backend, instrument, {"norm": 2.0, "index": -1.3}).values),
        dtype=float,
    )
    erg = u.erg / u.s / u.cm**2
    flux = (truth[order] * W_M2).to(erg)
    return LineFluxes(
        [TRANSITIONS[i] for i in order],
        np.round(REST[order], 2) * u.um,
        flux,
        uncertainty=FRACTION * flux,
    )


def composition(
    backend: ConformanceBackend,
    pieces: LineFluxPieces,
    observed: LineFluxes,
    *,
    noise: Any = None,
) -> FittingProblem:
    instrument = Instrument([step(pieces)], channel="default", label="lines")
    likelihood = Likelihood(GaussianFamily(), noise or backend.independent_noise())
    dataset = Dataset(observed.to_unit(W_M2), instrument, likelihood=likelihood, label="lines")
    return FittingProblem(model(backend), {"lines": dataset}, seed=20261014)


class TestTheComposition:
    def test_a_table_in_another_order_and_in_cgs_aligns_by_transition(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [2, 0, 1])
        problem = composition(backend, pieces, table)
        aligned = problem.datasets["lines"].observed
        assert isinstance(aligned, LineFluxes)
        assert list(aligned.transitions) == TRANSITIONS
        assert np.array_equal(aligned.rest_wavelength.values, REST)
        assert aligned.unit == W_M2
        at_truth = problem.log_prob(TRUTH)
        displaced = problem.log_prob({**TRUTH, "model.index": -0.8})
        assert np.isfinite(at_truth)
        assert at_truth > displaced

    def test_from_observed_takes_the_transitions_from_the_table(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [1, 2, 0])
        built = pieces.line_integration.from_observed(table, window=300.0)
        assert list(built.transitions) == list(table.transitions)
        np.testing.assert_array_equal(built.rest_wavelengths(), table.rest_wavelength.values)

    def test_an_observed_line_the_step_does_not_integrate_is_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [0, 1, 2])
        extra = LineFluxes(
            [*TRANSITIONS, "[S III] 18.71"],
            [*REST, 18.713] * u.um,
            np.append(table.flux, table.flux[0]) * table.unit,
            uncertainty=np.append(table.uncertainty, table.uncertainty[0]) * table.unit,
        )
        with pytest.raises(LikelihoodError, match=r"\['\[S III\] 18.71'\] are not tabulated"):
            composition(backend, pieces, extra)

    def test_a_step_line_without_an_observation_is_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [0, 1])
        with pytest.raises(LikelihoodError, match=r"\['\[Ne III\] 15.56'\] have no observation"):
            composition(backend, pieces, table)


# ---------------------------------------------------------------------------
# (f) simulate
# ---------------------------------------------------------------------------


class TestSimulate:
    def test_simulate_draws_a_line_flux_table(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        problem = composition(backend, pieces, observed_table(backend, pieces, [2, 1, 0]))
        simulation = problem.simulate(TRUTH, observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        drawn = simulation.observations["lines"]
        assert isinstance(drawn, LineFluxes)
        assert list(drawn.transitions) == TRANSITIONS
        assert not np.array_equal(drawn.values, simulation.predicted["lines"].values)


# ---------------------------------------------------------------------------
# (g) The flexible likelihood over rest wavelength
# ---------------------------------------------------------------------------


class TestTheFlexibleLikelihood:
    def test_the_dense_and_quasiseparable_solvers_agree(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """One axis, in no particular order: ``QuasisepGP`` sorts, ``DenseGP`` need not."""
        pieces = pieces_or_skip(backend)
        if not {SolverKind.DENSE, SolverKind.QUASISEP} <= backend.capabilities.solvers:
            pytest.skip(f"{backend.name} does not declare both GP solvers.")
        instrument = Instrument([step(pieces)], channel="default", label="lines")
        predicted = predict(backend, instrument, {"norm": 2.0, "index": -1.3})
        rng = np.random.default_rng(7)
        flux = np.asarray(backend.to_numpy(predicted.values), dtype=float)
        observed = LineFluxes(
            TRANSITIONS,
            REST * u.um,
            flux * (1.0 + 0.05 * rng.standard_normal(flux.size)) * W_M2,
            uncertainty=FRACTION * flux * W_M2,
        )
        spec = CovarianceSpec(
            KernelFamily.MATERN32, float(0.05 * flux.mean()), 2.0, axes=("spectral_axis",)
        )
        dense = Likelihood(
            GaussianFamily(),
            backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE)),
        )
        quasisep = Likelihood(
            GaussianFamily(),
            backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)),
        )
        dense_lp = dense.log_prob(predicted, observed, values={})
        quasisep_lp = quasisep.log_prob(predicted, observed, values={})
        assert np.isfinite(dense_lp)
        assert abs(dense_lp - quasisep_lp) < tolerances.cross_solver


# ---------------------------------------------------------------------------
# (h) The four diagnostics' plots
# ---------------------------------------------------------------------------


class TestTheFourPlots:
    def test_the_four_plots_run_on_a_stored_run(self, backend: ConformanceBackend) -> None:
        """A smoke row, on the reference column: the kind needs nothing of the plots."""
        if backend.name != "reference":
            pytest.skip("the plots are backend-neutral; one column carries the smoke row.")
        pytest.importorskip("matplotlib")
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        from ampere.inference import EmceeEngine
        from ampere.results import (
            add_posterior_predictive,
            add_residuals,
            gp_localisation,
            gp_localisation_score,
            plot_anomaly_score,
            plot_gp_localisation,
            plot_posterior_predictive,
            plot_residuals,
        )

        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [2, 0, 1])
        flux = np.asarray(table.to_unit(W_M2).flux, dtype=float)
        spec = CovarianceSpec(
            KernelFamily.MATERN32, float(0.05 * flux.mean()), 2.0, axes=("spectral_axis",)
        )
        noise = backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE))
        problem = composition(backend, pieces, table, noise=noise)
        # A smoke row: a 20-step run does not converge, and a GP fit's residual
        # plot says its whiteness test is close to circular. Both warnings are
        # the library's correct answers to a deliberately tiny run.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            tree = EmceeEngine(problem, walkers=8).run(20, initial="prior")
            localised = gp_localisation(tree, problem)
            try:
                figures = [
                    plot_posterior_predictive(add_posterior_predictive(tree, problem)),
                    plot_residuals(add_residuals(tree, problem)),
                    plot_gp_localisation(localised),
                    plot_anomaly_score(gp_localisation_score(localised, dataset="lines")),
                ]
                assert all(figure is not None for figure in figures)
            finally:
                plt.close("all")


# ---------------------------------------------------------------------------
# The step's refusals
# ---------------------------------------------------------------------------


class TestTheRefusals:
    def test_overlapping_windows_are_refused_by_name(self, backend: ConformanceBackend) -> None:
        """A sample inside two windows would have two owners."""
        pieces = pieces_or_skip(backend)
        with pytest.raises(TransformationError, match=r"'\[Ne II\] 12.81'.*overlap"):
            pieces.line_integration(
                ["H I 7-6", "[Ne II] 12.81"], [12.3719, 12.8135], window=0.25, window_unit="um"
            )

    def test_a_model_not_in_jansky_is_refused_by_name(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        line_step = step(pieces, samples_per_window=5)
        grid = line_step.tabulation()
        spectrum = Spectrum(grid * u.um, np.ones(grid.size) * u.mJy)
        with pytest.raises(TransformationError, match="flux density in Jy"):
            line_step(spectrum, None)

    def test_an_unknown_continuum_or_window_unit_is_refused(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        with pytest.raises(TransformationError, match="continuum"):
            step(pieces, continuum="quadratic")
        with pytest.raises(TransformationError, match="window_unit"):
            step(pieces, window_unit="Angstrom")

    def test_from_observed_reads_the_tables_own_windows(self, backend: ConformanceBackend) -> None:
        """A table carrying ``extra_coords['window']`` needs no ``window=``."""
        pieces = pieces_or_skip(backend)
        table = LineFluxes(
            TRANSITIONS,
            REST * u.um,
            np.ones(3) * W_M2,
            extra_coords={"window": np.array([200.0, 250.0, 300.0])},
        )
        built = pieces.line_integration.from_observed(table)
        np.testing.assert_array_equal(built.windows(), [200.0, 250.0, 300.0])
        bare = LineFluxes(TRANSITIONS, REST * u.um, np.ones(3) * W_M2)
        with pytest.raises(TransformationError, match="window"):
            pieces.line_integration.from_observed(bare)


# ---------------------------------------------------------------------------
# (i) The native realisation agrees with the contract path
# ---------------------------------------------------------------------------


class TestTheRealisation:
    """The native twins' ``apply_flux``, composed by ``realise``, against the contract path."""

    @pytest.mark.parametrize("flexible", [False, True], ids=["independent", "dense_gp"])
    def test_the_realised_density_is_the_contract_density(
        self, backend: ConformanceBackend, tolerances: Tolerances, flexible: bool
    ) -> None:
        from ampere.core import realise, registered_realisations

        pieces = pieces_or_skip(backend)
        table = observed_table(backend, pieces, [1, 2, 0])
        noise = None
        if flexible:
            flux = np.asarray(table.to_unit(W_M2).flux, dtype=float)
            spec = CovarianceSpec(
                KernelFamily.MATERN32, float(0.05 * flux.mean()), 2.0, axes=("spectral_axis",)
            )
            noise = backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE))
        problem = composition(backend, pieces, table, noise=noise)
        if problem.backend not in registered_realisations():
            pytest.skip(f"the {problem.backend!r} backend registers no realisation.")
        realised = realise(problem)
        for index in (-1.6, -1.3, -0.9):
            theta = problem.parameters.pack({**TRUTH, "model.index": index})
            y = problem.unconstrain(theta)
            native = float(np.asarray(backend.to_numpy(realised.log_prob_unconstrained(y))))
            contract = float(problem.log_prob_unconstrained(y))
            assert native == pytest.approx(contract, abs=tolerances.cross_backend)


# ---------------------------------------------------------------------------
# A stored record of the kind reads back
# ---------------------------------------------------------------------------


class TestTheRecord:
    def test_a_line_flux_record_round_trips_by_kind_name(self, backend: ConformanceBackend) -> None:
        """``ampere.results`` rebuilds a container by its kind's name (training sets do)."""
        if backend.name != "reference":
            pytest.skip("the record is backend-neutral; one column carries the row.")
        from ampere.results import container_from_dict, container_to_dict, kind_named

        table = LineFluxes(
            TRANSITIONS, REST * u.um, np.ones(3) * W_M2, uncertainty=0.1 * np.ones(3) * W_M2
        )
        assert kind_named("LineFluxes") is LineFluxes
        back = container_from_dict(container_to_dict(table))
        assert isinstance(back, LineFluxes)
        assert back.transitions.tolist() == TRANSITIONS
        np.testing.assert_array_equal(back.rest_wavelength.values, REST)
