"""Astrometry: the rows the template's second modality owes (W4.9).

One column per registered backend, like every other file here, and no test
body names one. A backend that has not written the astrometric vocabulary
declares ``BackendCapabilities.astrometry = False`` and every row skips with a
reason naming what is owed — the same shape
``tests/conformance/test_interferometry.py`` uses for its own vocabulary.

What these rows are for, following the item's own list:

* the reflex-orbit model, on both its channels, against the **closed-form
  ephemeris** in :mod:`tests.conformance.oracles` — an oracle transcribed
  independently of the model, not a second call to it;
* the epoch-sampling step's published requirement (``points=`` at the exact
  observed epochs) and its identity ``apply``, mask included;
* two channels sharing one model evaluation, negotiated and compiled once per
  draw, each bound by its own dataset and instrument — the same "one model,
  several channels" shape ``test_interferometry.py`` proves for two
  *instruments* on one channel, proved here for one model's two *channels*;
* the flexible likelihood — ``GaussianProcessNoise(Matern32(...),
  QuasisepGP())`` — composing on both channels, since a ``TimeSeries``'s one
  ordered ``time`` axis is exactly what the O(N) path wants;
* ``simulate(observe=True)`` drawing both channels.
"""

from __future__ import annotations

import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    AxisRequirement,
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    TimeSeries,
    negotiate,
)

from .oracles import reflex_orbit_dec, reflex_orbit_ra
from .protocol import (
    AstrometryPieces,
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    SolverKind,
    Tolerances,
)

#: The observation epochs, days. Irregular, on purpose (real astrometric
#: scheduling is), and long enough — about two periods — that the periodic
#: term and the linear drift are both constrained.
EPOCHS = np.array([0.0, 45.0, 130.0, 210.0, 305.0, 390.0, 460.0, 540.0, 615.0, 700.0, 780.0])

#: The truth every row here is evaluated at.
ORBIT: dict[str, float] = {
    "pmra": 1.2,
    "pmdec": -0.6,
    "period": 400.0,
    "phase": 0.7,
    "amp_ra": 0.6,
    "amp_dec": 0.35,
}

#: The free parameters of a fitted problem, and their true values.
TRUTH: dict[str, float] = {
    "model.pmra": ORBIT["pmra"],
    "model.pmdec": ORBIT["pmdec"],
}

#: Per-epoch measurement uncertainty, mas.
SIGMA = 0.03


def pieces_or_skip(backend: ConformanceBackend) -> AstrometryPieces:
    """This backend's astrometric classes, or a skip naming what is owed."""
    if not backend.capabilities.astrometry:
        pytest.skip(
            f"{backend.name} declares no astrometric vocabulary "
            f"(BackendCapabilities.astrometry), so the epoch-sampling step and the reflex-orbit "
            f"model it would need are not there. W4.9 owes the native twins."
        )
    return backend.astrometry()


def observed_channel(channel: str, values: np.ndarray | None = None) -> TimeSeries:
    """A ``TimeSeries`` on :data:`EPOCHS`, with a uniform sigma."""
    filled = np.zeros(EPOCHS.size) if values is None else values
    return TimeSeries(
        EPOCHS * u.day, filled * u.mas, uncertainty=np.full(EPOCHS.size, SIGMA) * u.mas
    )


def orbit_model(pieces: AstrometryPieces, **overrides: Any) -> Any:
    """The reflex orbit on its own default epoch grid, at the truth (or overridden)."""
    return pieces.reflex_orbit(EPOCHS * u.day, **{**ORBIT, **overrides})


def chain(pieces: AstrometryPieces, observed: TimeSeries, channel: str, label: str) -> Instrument:
    """The one-step epoch-sampling instrument for one coordinate."""
    return Instrument(
        [pieces.epoch_sample.from_observed(observed)], channel=channel, label=label
    )


# ---------------------------------------------------------------------------
# The reflex-orbit model against the closed-form ephemeris
# ---------------------------------------------------------------------------


class TestReflexOrbitAgainstTheClosedForm:
    """The model's ``evaluate()``, on both channels, against an independent oracle."""

    @pytest.mark.parametrize("channel", ["ra", "dec"])
    def test_the_offset_matches_the_closed_form_ephemeris(
        self, backend: ConformanceBackend, tolerances: Tolerances, channel: str
    ) -> None:
        pieces = pieces_or_skip(backend)
        model = orbit_model(pieces)
        result = model.evaluate()
        got = backend.to_numpy(result[channel].values)
        if channel == "ra":
            expected = reflex_orbit_ra(
                EPOCHS, pmra=ORBIT["pmra"], period=ORBIT["period"], phase=ORBIT["phase"],
                amp_ra=ORBIT["amp_ra"],
            )
        else:
            expected = reflex_orbit_dec(
                EPOCHS, pmdec=ORBIT["pmdec"], period=ORBIT["period"], phase=ORBIT["phase"],
                amp_dec=ORBIT["amp_dec"],
            )
        assert np.allclose(got, expected, rtol=0.0, atol=tolerances.analytic)

    def test_period_and_phase_are_shared_by_one_evaluation(
        self, backend: ConformanceBackend
    ) -> None:
        """``period``/``phase`` tie the two channels by construction, needing no ``Tie``."""
        pieces = pieces_or_skip(backend)
        model = orbit_model(pieces, phase=1.3)
        result = model.evaluate()
        ra = backend.to_numpy(result["ra"].values)
        dec = backend.to_numpy(result["dec"].values)
        expected_ra = reflex_orbit_ra(EPOCHS, pmra=ORBIT["pmra"], period=ORBIT["period"], phase=1.3, amp_ra=ORBIT["amp_ra"])
        expected_dec = reflex_orbit_dec(EPOCHS, pmdec=ORBIT["pmdec"], period=ORBIT["period"], phase=1.3, amp_dec=ORBIT["amp_dec"])
        assert np.allclose(ra, expected_ra, rtol=0.0, atol=1e-9)
        assert np.allclose(dec, expected_dec, rtol=0.0, atol=1e-9)


# ---------------------------------------------------------------------------
# The epoch-sampling step
# ---------------------------------------------------------------------------


class TestEpochSample:
    """The instrument the design sketch's §2 asks for: no arithmetic at all."""

    def test_the_published_requirement_is_points_at_the_exact_epochs(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        step = pieces.epoch_sample.from_observed(observed_channel("ra"))
        (requirement,) = step.requirements()
        assert isinstance(requirement, AxisRequirement)
        assert requirement.axis == "time"
        assert np.array_equal(np.sort(requirement.points), np.sort(EPOCHS))

    def test_apply_is_the_identity(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        observed = observed_channel("ra")
        step = pieces.epoch_sample.from_observed(observed)
        got = step(observed, None)
        assert got is observed

    def test_a_masked_epoch_stays_masked(self, backend: ConformanceBackend) -> None:
        """The identity step still owes ``results_schema.md`` §16's rule."""
        pieces = pieces_or_skip(backend)
        mask = np.zeros(EPOCHS.size, dtype=bool)
        mask[3] = True
        observed = TimeSeries(
            EPOCHS * u.day,
            np.zeros(EPOCHS.size) * u.mas,
            uncertainty=np.full(EPOCHS.size, SIGMA) * u.mas,
            mask=mask,
        )
        step = pieces.epoch_sample.from_observed(observed)
        got = step(observed, None)
        assert got.mask is not None
        assert np.array_equal(got.mask, mask)


# ---------------------------------------------------------------------------
# Two channels, one model
# ---------------------------------------------------------------------------


def counting(model_class: type) -> type:
    """*model_class* with an ``evaluations`` counter, for the once-per-draw row."""

    class Counting(model_class):  # type: ignore[valid-type, misc]
        def __init__(self, *args: Any, **kwargs: Any) -> None:
            super().__init__(*args, **kwargs)
            self.evaluations = 0

        def evaluate(self, **values: Any) -> Any:
            self.evaluations += 1
            return super().evaluate(**values)

    Counting.__name__ = f"Counting{model_class.__name__}"
    Counting.__qualname__ = Counting.__name__
    return Counting


class TestTwoChannelsOneModel:
    """RA and Dec, negotiated and compiled once — the composition this modality proves."""

    def _problem(
        self, backend: ConformanceBackend, pieces: AstrometryPieces, *, gp: bool = False
    ) -> tuple[FittingProblem, Any]:
        ra_instrument = chain(pieces, observed_channel("ra"), "ra", "astrom_ra")
        dec_instrument = chain(pieces, observed_channel("dec"), "dec", "astrom_dec")
        truth = orbit_model(pieces)
        compiled = truth.compile_for(negotiate([ra_instrument, dec_instrument]))
        result = compiled.evaluate()
        ra_observed = observed_channel("ra", backend.to_numpy(ra_instrument(result).values))
        dec_observed = observed_channel("dec", backend.to_numpy(dec_instrument(result).values))
        model = counting(pieces.reflex_orbit)(
            EPOCHS * u.day,
            pmra=st.norm(0.0, 5.0),
            pmdec=st.norm(0.0, 5.0),
            period=ORBIT["period"],
            phase=ORBIT["phase"],
            amp_ra=ORBIT["amp_ra"],
            amp_dec=ORBIT["amp_dec"],
        )
        if gp:
            spec = CovarianceSpec(KernelFamily.MATERN32, 0.03, 120.0, axes=("time",))
            ra_noise: Any = backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP))
            dec_noise: Any = backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP))
        else:
            ra_noise = backend.independent_noise()
            dec_noise = backend.independent_noise()
        datasets = DatasetCollection(
            {
                "ra": Dataset(
                    ra_observed,
                    ra_instrument,
                    likelihood=Likelihood(GaussianFamily(), ra_noise),
                    label="ra",
                ),
                "dec": Dataset(
                    dec_observed,
                    dec_instrument,
                    likelihood=Likelihood(GaussianFamily(), dec_noise),
                    label="dec",
                ),
            }
        )
        return FittingProblem(model, datasets, seed=20260913), model

    def test_both_channels_are_recorded_on_the_model(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        assert set(problem.requirements["model"]) == {"ra", "dec"}
        assert {type(problem.datasets[name].observed).__name__ for name in ("ra", "dec")} == {
            "TimeSeries"
        }

    def test_the_model_is_evaluated_once_per_draw(self, backend: ConformanceBackend) -> None:
        """``inference.md`` §8: once, jointly — not once per channel."""
        pieces = pieces_or_skip(backend)
        problem, model = self._problem(backend, pieces)
        model.evaluations = 0
        theta = {**TRUTH, "model.period": ORBIT["period"], "model.phase": ORBIT["phase"]}
        problem.log_prob(theta)
        assert model.evaluations == 1

    def test_the_composed_likelihood_peaks_at_the_truth(self, backend: ConformanceBackend) -> None:
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        theta = {**TRUTH, "model.period": ORBIT["period"], "model.phase": ORBIT["phase"]}
        at_truth = problem.log_prob(theta)
        displaced = problem.log_prob({**theta, "model.pmra": TRUTH["model.pmra"] + 5.0})
        assert np.isfinite(at_truth)
        assert at_truth > displaced

    def test_simulate_draws_observations_of_both_channels(
        self, backend: ConformanceBackend
    ) -> None:
        pieces = pieces_or_skip(backend)
        problem, _ = self._problem(backend, pieces)
        theta = {**TRUTH, "model.period": ORBIT["period"], "model.phase": ORBIT["phase"]}
        simulation = problem.simulate(theta, observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        drawn_ra = simulation.observations["ra"]
        drawn_dec = simulation.observations["dec"]
        assert drawn_ra.values.dtype.kind == "f"
        assert drawn_dec.values.dtype.kind == "f"
        # A draw is not the prediction: the noise actually moved the data.
        assert not np.array_equal(drawn_ra.values, simulation.predicted["ra"].values)
        assert not np.array_equal(drawn_dec.values, simulation.predicted["dec"].values)

    def test_the_flexible_likelihood_composes_on_both_channels(
        self, backend: ConformanceBackend
    ) -> None:
        """``QuasisepGP`` on a ``TimeSeries``: the O(N) path this modality is chosen for."""
        pieces = pieces_or_skip(backend)
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} declares no QuasisepGP solver.")
        problem, _ = self._problem(backend, pieces, gp=True)
        theta = {**TRUTH, "model.period": ORBIT["period"], "model.phase": ORBIT["phase"]}
        assert np.isfinite(problem.log_prob(theta))

    def test_the_dense_and_quasiseparable_solvers_agree(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The O(N) recursion and the dense Cholesky score the same residual the same way."""
        pieces = pieces_or_skip(backend)
        if not {SolverKind.DENSE, SolverKind.QUASISEP} <= backend.capabilities.solvers:
            pytest.skip(f"{backend.name} does not declare both GP solvers.")
        ra_instrument = chain(pieces, observed_channel("ra"), "ra", "astrom_ra")
        model = orbit_model(pieces)
        compiled = model.compile_for(negotiate([ra_instrument]))
        result = compiled.evaluate()
        observed = observed_channel("ra", backend.to_numpy(ra_instrument(result).values))
        spec = CovarianceSpec(KernelFamily.MATERN32, 0.03, 120.0, axes=("time",))
        dense = Likelihood(
            GaussianFamily(), backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE))
        )
        quasisep = Likelihood(
            GaussianFamily(),
            backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)),
        )
        predicted = ra_instrument(compiled.evaluate())
        dense_lp = dense.log_prob(predicted, observed, values={})
        quasisep_lp = quasisep.log_prob(predicted, observed, values={})
        assert abs(dense_lp - quasisep_lp) < tolerances.cross_solver
