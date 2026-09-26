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
    RotationCoupling,
    TimeSeries,
    negotiate,
)

from ampere.core.exceptions import LikelihoodError

from .oracles import (
    coregionalised_covariance,
    coregionalised_log_density,
    kernel_matrix,
    reflex_orbit_dec,
    reflex_orbit_ra,
    rotation_coupling_matrix,
)
from .protocol import (
    AstrometryPieces,
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    SolverKind,
    Tolerances,
    approximation_envelope,
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


def observed_channel(
    channel: str, values: np.ndarray | None = None, sigma: np.ndarray | None = None
) -> TimeSeries:
    """A ``TimeSeries`` on :data:`EPOCHS`, with a uniform sigma unless *sigma* is given."""
    filled = np.zeros(EPOCHS.size) if values is None else values
    uncertainty = np.full(EPOCHS.size, SIGMA) if sigma is None else np.asarray(sigma)
    return TimeSeries(EPOCHS * u.day, filled * u.mas, uncertainty=uncertainty * u.mas)


def orbit_model(pieces: AstrometryPieces, **overrides: Any) -> Any:
    """The reflex orbit on its own default epoch grid, at the truth (or overridden)."""
    return pieces.reflex_orbit(EPOCHS * u.day, **{**ORBIT, **overrides})


def chain(pieces: AstrometryPieces, observed: TimeSeries, channel: str, label: str) -> Instrument:
    """The one-step epoch-sampling instrument for one coordinate."""
    return Instrument([pieces.epoch_sample.from_observed(observed)], channel=channel, label=label)


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
                EPOCHS,
                pmra=ORBIT["pmra"],
                period=ORBIT["period"],
                phase=ORBIT["phase"],
                amp_ra=ORBIT["amp_ra"],
            )
        else:
            expected = reflex_orbit_dec(
                EPOCHS,
                pmdec=ORBIT["pmdec"],
                period=ORBIT["period"],
                phase=ORBIT["phase"],
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
        expected_ra = reflex_orbit_ra(
            EPOCHS, pmra=ORBIT["pmra"], period=ORBIT["period"], phase=1.3, amp_ra=ORBIT["amp_ra"]
        )
        expected_dec = reflex_orbit_dec(
            EPOCHS,
            pmdec=ORBIT["pmdec"],
            period=ORBIT["period"],
            phase=1.3,
            amp_dec=ORBIT["amp_dec"],
        )
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
            ra_noise: Any = backend.gp_noise(
                backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)
            )
            dec_noise: Any = backend.gp_noise(
                backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)
            )
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
            GaussianFamily(),
            backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE)),
        )
        quasisep = Likelihood(
            GaussianFamily(),
            backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.QUASISEP)),
        )
        predicted = ra_instrument(compiled.evaluate())
        dense_lp = dense.log_prob(predicted, observed, values={})
        quasisep_lp = quasisep.log_prob(predicted, observed, values={})
        assert abs(dense_lp - quasisep_lp) < tolerances.cross_solver


# ---------------------------------------------------------------------------
# Joint noise over the two channels (W5.9)
# ---------------------------------------------------------------------------

#: The coupling every joint row below is evaluated at: an error elongated
#: along a position angle of 0.7 rad with semi-axis variances a factor of nine
#: apart, which is a systematic two independent GPs cannot express at all.
COUPLING: dict[str, float] = {
    "angle": 0.7,
    "log_variance_0": float(np.log(9.0e-4)),
    "log_variance_1": float(np.log(1.0e-4)),
}

#: The shared-grid kernel. Amplitude **fixed** at one: ``B (x) (a^2 K) =
#: (a^2 B) (x) K``, so a free amplitude beside a free ``B`` is the same degree
#: of freedom twice and ``JointGaussianProcessNoise`` refuses it by name.
JOINT_KERNEL = CovarianceSpec(KernelFamily.MATERN32, 1.0, 180.0, axes=("time",))


def joint_or_skip(backend: ConformanceBackend) -> None:
    """Skip a joint row on a backend that has not written the noise model."""
    if not backend.capabilities.joint_noise:
        pytest.skip(
            f"{backend.name} declares no joint noise model "
            f"(BackendCapabilities.joint_noise), so the shared-grid intrinsic coregionalisation "
            f"model W5.9 adds is not there."
        )


def kronecker_covariance(
    backend: ConformanceBackend, noise: Any, coordinates: np.ndarray, variance: np.ndarray
) -> np.ndarray:
    """``B (x) K_x + I_T (x) diag(variance)``, materialised — the definition of the answer.

    Channel-major: block ``(s, t)`` is ``B[s, t] K_x``. Built from the
    *declarations* — the coupling's own :meth:`~ampere.core.ChannelCoupling.
    matrix` and the kernel's own ``matrix`` — so this is not a second call to
    the thing under test but the same declaration assembled the obvious,
    ``O((NT)^3)`` way.
    """
    coupling = np.asarray(backend.to_numpy(noise.coupling_matrix(COUPLING)), dtype=float)
    kernel = np.asarray(
        backend.to_numpy(noise.kernel.matrix(coordinates, coordinates, noise.kernel.resolve({}))),
        dtype=float,
    )
    channels = coupling.shape[0]
    return np.kron(coupling, kernel) + np.eye(channels * kernel.shape[0]) * np.tile(
        variance, channels
    )


def dense_log_density(covariance: np.ndarray, residual: np.ndarray) -> float:
    """``log N(residual; 0, covariance)`` by a plain numpy Cholesky."""
    factor = np.linalg.cholesky(covariance)
    solved = np.linalg.solve(covariance, residual)
    log_determinant = 2.0 * float(np.sum(np.log(np.diag(factor))))
    return float(
        -0.5 * (float(residual @ solved) + log_determinant + residual.size * np.log(2.0 * np.pi))
    )


class TestJointChannelNoise:
    """One correlated process over ``ra`` and ``dec``: ``K = B (x) K_x`` (*W5.9*)."""

    def _pieces(
        self,
        backend: ConformanceBackend,
        pieces: AstrometryPieces,
        *,
        solver: SolverKind = SolverKind.DENSE,
        coupling: Any = None,
        kernel: CovarianceSpec | None = None,
        sigmas: np.ndarray | None = None,
        basis_size: int = 32,
        boundary_factor: float = 2.0,
    ) -> tuple[FittingProblem, Any, np.ndarray, np.ndarray]:
        """A joint problem, its noise model, the residual block and the variance.

        *sigmas* (**W5.24**), a ``(2, N)`` array, gives each channel its own
        per-epoch uncertainty: the white noise is drawn at it, the containers
        carry it, and the variance returned is the ``(N, 2)`` block rather than
        W5.9's shared vector. Left ``None``, every draw and every value is
        exactly W5.9's, so the rows written before W5.24 see the same data.
        """
        ra_instrument = chain(pieces, observed_channel("ra"), "ra", "astrom_ra")
        dec_instrument = chain(pieces, observed_channel("dec"), "dec", "astrom_dec")
        truth = orbit_model(pieces)
        compiled = truth.compile_for(negotiate([ra_instrument, dec_instrument]))
        result = compiled.evaluate()
        noiseless = {
            "ra": backend.to_numpy(ra_instrument(result).values).ravel(),
            "dec": backend.to_numpy(dec_instrument(result).values).ravel(),
        }
        rng = np.random.default_rng(20260916)
        if sigmas is None:
            observed = {
                name: observed_channel(name, values + rng.normal(0.0, SIGMA, values.shape))
                for name, values in noiseless.items()
            }
        else:
            observed = {
                name: observed_channel(
                    name,
                    values + rng.normal(0.0, sigmas[index], values.shape),
                    sigmas[index],
                )
                for index, (name, values) in enumerate(noiseless.items())
            }
        joint = backend.joint_gp_noise(
            backend.kernel(JOINT_KERNEL if kernel is None else kernel),
            backend.gp_solver(solver, basis_size=basis_size, boundary_factor=boundary_factor),
            datasets=("ra", "dec"),
            coupling=RotationCoupling(*COUPLING.values()) if coupling is None else coupling,
        )
        model = pieces.reflex_orbit(
            EPOCHS * u.day,
            pmra=st.norm(0.0, 5.0),
            pmdec=st.norm(0.0, 5.0),
            period=ORBIT["period"],
            phase=ORBIT["phase"],
            amp_ra=ORBIT["amp_ra"],
            amp_dec=ORBIT["amp_dec"],
        )
        datasets = DatasetCollection(
            {
                "ra": Dataset(
                    observed["ra"],
                    ra_instrument,
                    likelihood=Likelihood(GaussianFamily(), backend.independent_noise()),
                    label="ra",
                ),
                "dec": Dataset(
                    observed["dec"],
                    dec_instrument,
                    likelihood=Likelihood(GaussianFamily(), backend.independent_noise()),
                    label="dec",
                ),
            },
            joint={"astrom": joint},
        )
        problem = FittingProblem(model, datasets, seed=20260916)
        residual = np.concatenate(
            [
                np.asarray(observed[name].values, dtype=float).ravel() - noiseless[name]
                for name in ("ra", "dec")
            ]
        )
        variance = np.full(EPOCHS.size, SIGMA**2) if sigmas is None else np.asarray(sigmas).T ** 2
        return problem, joint, residual, variance

    def test_the_joint_density_matches_a_dense_kronecker_solve(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The rotated ``T`` solves against ``B (x) K_x`` materialised and factorised."""
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        problem, noise, residual, variance = self._pieces(backend, pieces)
        coordinates = EPOCHS.reshape(-1, 1)
        reference = dense_log_density(
            kronecker_covariance(backend, noise, coordinates, variance), residual
        )
        got = noise.log_prob(
            [residual[: EPOCHS.size], residual[EPOCHS.size :]],
            variance,
            coordinates,
            COUPLING,
        )
        assert abs(got - reference) < tolerances.cross_solver * max(1.0, abs(reference))
        # And the same number arrives through the composed problem, which is
        # the thing a fit actually evaluates.
        through = problem.log_likelihood(
            {
                "model.pmra": ORBIT["pmra"],
                "model.pmdec": ORBIT["pmdec"],
                "astrom.angle": COUPLING["angle"],
                "astrom.log_variance_0": COUPLING["log_variance_0"],
                "astrom.log_variance_1": COUPLING["log_variance_1"],
            }
        )
        assert abs(through - reference) < tolerances.cross_solver * max(1.0, abs(reference))

    def test_the_quasiseparable_solver_scores_the_same_joint_density(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``QuasisepGP`` on the rotated outputs: the O(N) path this model keeps."""
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} declares no QuasisepGP solver.")
        _, dense, residual, variance = self._pieces(backend, pieces)
        _, quasisep, _, _ = self._pieces(backend, pieces, solver=SolverKind.QUASISEP)
        coordinates = EPOCHS.reshape(-1, 1)
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        assert (
            abs(
                dense.log_prob(columns, variance, coordinates, COUPLING)
                - quasisep.log_prob(columns, variance, coordinates, COUPLING)
            )
            < tolerances.cross_solver
        )

    def test_the_rotation_recovers_t_independent_solves(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """Scored as two ordinary scalar GPs on the rotated residuals, it is the same number.

        The mechanism, checked rather than assumed: rotating by ``Q^T``
        decouples the channels, so the joint density is the sum of ``T``
        *ordinary* ``GaussianProcessNoise`` densities on the rotated residuals,
        each with the kernel scaled by its own eigenvalue. This row builds
        those ``T`` likelihoods from the shipped scalar noise model — no joint
        machinery at all — and adds them up.
        """
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        _, noise, residual, variance = self._pieces(backend, pieces)
        coordinates = EPOCHS.reshape(-1, 1)
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        eigenvalues, rotation = noise.eigen(COUPLING)
        eigenvalues = np.asarray(backend.to_numpy(eigenvalues), dtype=float).ravel()
        rotated = np.column_stack(columns) @ np.asarray(backend.to_numpy(rotation), dtype=float)
        total = 0.0
        for index, eigenvalue in enumerate(eigenvalues):
            # k(0) is amplitude**2, so an eigenvalue of the coupling is an
            # amplitude of its square root.
            scalar = backend.gp_noise(
                backend.kernel(
                    CovarianceSpec(
                        KernelFamily.MATERN32,
                        float(np.sqrt(eigenvalue)),
                        JOINT_KERNEL.length_scale,
                        axes=("time",),
                    )
                ),
                backend.gp_solver(SolverKind.DENSE),
            )
            likelihood = Likelihood(GaussianFamily(), scalar)
            observed = observed_channel("ra", rotated[:, index])
            predicted = observed_channel("ra", np.zeros(EPOCHS.size))
            total += likelihood.log_prob(predicted, observed, values={})
        got = noise.log_prob(columns, variance, coordinates, COUPLING)
        assert abs(got - total) < tolerances.cross_solver * max(1.0, abs(total))

    def test_an_identity_coupling_is_two_independent_noise_models(
        self, backend: ConformanceBackend
    ) -> None:
        """``B = I``: **bit-identical** to two independent ``GaussianProcessNoise`` datasets.

        Not "agrees to a tolerance". With ``B = I`` the rotation is the
        identity and the eigenvalues are one, so every floating-point operation
        the joint path performs is the one the scalar path performs, in the
        same order — and any difference at all would mean the joint path had
        introduced an operation that does not belong to the model.
        """
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        spec = CovarianceSpec(KernelFamily.MATERN32, 0.02, 180.0, axes=("time",))
        _, noise, residual, variance = self._pieces(
            backend,
            pieces,
            coupling=RotationCoupling(0.0, 0.0, 0.0),
            kernel=spec,
        )
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        got = noise.log_prob(
            columns,
            variance,
            EPOCHS.reshape(-1, 1),
            {"angle": 0.0, "log_variance_0": 0.0, "log_variance_1": 0.0},
        )
        separate = 0.0
        for column in columns:
            scalar = Likelihood(
                GaussianFamily(),
                backend.gp_noise(backend.kernel(spec), backend.gp_solver(SolverKind.DENSE)),
            )
            separate += scalar.log_prob(
                observed_channel("ra", np.zeros(EPOCHS.size)),
                observed_channel("ra", column),
                values={},
            )
        assert got == separate

    def test_simulate_draws_correlated_channels(self, backend: ConformanceBackend) -> None:
        """The sample covariance of ``simulate(observe=True)`` against ``B (x) K_x``.

        The row that catches the failure this whole model exists to prevent: a
        draw taken channel by channel is a draw from a covariance whose
        off-diagonal block is **zero**, and every marginal of it looks right.
        So the statistic compared here is the whole ``2N x 2N`` covariance,
        cross-channel block included, and the row also asserts that the block
        is not zero — a run that had quietly dropped the correlation would pass
        a marginal check and fail this one.
        """
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        problem, noise, _, variance = self._pieces(backend, pieces)
        theta = {
            "model.pmra": ORBIT["pmra"],
            "model.pmdec": ORBIT["pmdec"],
            "astrom.angle": COUPLING["angle"],
            "astrom.log_variance_0": COUPLING["log_variance_0"],
            "astrom.log_variance_1": COUPLING["log_variance_1"],
        }
        draws = 800
        rows = np.empty((draws, 2 * EPOCHS.size))
        simulations = problem.simulate_many(draws, values=[theta] * draws, observe=True)
        for index, simulation in enumerate(simulations):
            assert not simulation.failed
            assert simulation.observations is not None
            rows[index] = np.concatenate(
                [
                    np.asarray(simulation.observations[name].values, dtype=float).ravel()
                    - np.asarray(simulation.predicted[name].values, dtype=float).ravel()
                    for name in ("ra", "dec")
                ]
            )
        expected = kronecker_covariance(backend, noise, EPOCHS.reshape(-1, 1), variance)
        empirical = np.cov(rows, rowvar=False)
        scale = float(np.max(np.abs(expected)))
        error = float(np.max(np.abs(empirical - expected))) / scale
        assert error < 0.25, f"sample covariance is {error:.3f} away from B (x) K_x"
        # ... and the cross-channel block is genuinely there.
        cross = empirical[: EPOCHS.size, EPOCHS.size :]
        assert float(np.max(np.abs(cross))) > 0.2 * scale


# ---------------------------------------------------------------------------
# Heteroscedastic channels: the dense and reduced-rank routes (W5.24)
# ---------------------------------------------------------------------------

#: Each channel's own per-epoch sigma: a log-uniform factor of up to two
#: either side of :data:`SIGMA`, drawn once, seeded. ``(2, N)``: row ``t`` is
#: channel ``t``'s. Unequal everywhere, so no epoch lets the rotation through.
CHANNEL_SIGMAS = SIGMA * np.exp(
    np.random.default_rng(20260926).uniform(-np.log(2.0), np.log(2.0), size=(2, EPOCHS.size))
)

#: The same heteroscedasticity *along* the grid, shared by both channels: the
#: case W5.9's rotated path already scores exactly, and the one the
#: bit-identity rows hold it to.
SHARED_SIGMAS = np.vstack([CHANNEL_SIGMAS[0], CHANNEL_SIGMAS[0]])

#: The reduced-rank sweep. The box is three data half-extents wide: a
#: 180-day Matern-3/2 against a 780-day grid is poorly approximated within a
#: length scale of the boundary, and at the default factor of two the box's
#: own truncation error sits above ``approximation_final`` whatever ``m`` is
#: (measured: 3.5e-2 at m = 256). At three it falls with ``m`` down to the
#: floor.
REDUCED_RANK_SIZES = (8, 16, 32, 64)
REDUCED_RANK_BOX = 3.0

THETA_JOINT: dict[str, float] = {
    "model.pmra": ORBIT["pmra"],
    "model.pmdec": ORBIT["pmdec"],
    "astrom.angle": COUPLING["angle"],
    "astrom.log_variance_0": COUPLING["log_variance_0"],
    "astrom.log_variance_1": COUPLING["log_variance_1"],
}


def oracle_log_density(variances: np.ndarray, residual: np.ndarray) -> float:
    """The joint density at :data:`COUPLING` from :mod:`tests.conformance.oracles` alone."""
    coupling = rotation_coupling_matrix(
        COUPLING["angle"], COUPLING["log_variance_0"], COUPLING["log_variance_1"]
    )
    kernel = kernel_matrix(
        JOINT_KERNEL.family, EPOCHS, JOINT_KERNEL.amplitude, JOINT_KERNEL.length_scale
    )
    block = np.column_stack([residual[: EPOCHS.size], residual[EPOCHS.size :]])
    return coregionalised_log_density(coupling, kernel, variances, block)


class TestHeteroscedasticJointNoise:
    """Unequal per-channel sigmas: two routes behind one declaration (*W5.24*)."""

    def test_the_dense_route_matches_the_materialised_matrix(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``DenseGP`` bound: one Cholesky of ``B (x) K_x + blockdiag(diag sigma_t^2)``."""
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        problem, noise, residual, variance = self._pieces(backend, pieces, CHANNEL_SIGMAS)
        assert noise.route(variance) == "dense"
        reference = oracle_log_density(variance, residual)
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        got = noise.log_prob(columns, variance, EPOCHS.reshape(-1, 1), COUPLING)
        assert abs(got - reference) < tolerances.cross_solver * max(1.0, abs(reference))
        # The composed problem hands the group every channel's own sigma.
        through = problem.log_likelihood(THETA_JOINT)
        assert abs(through - reference) < tolerances.cross_solver * max(1.0, abs(reference))

    def test_the_reduced_rank_route_converges_to_the_dense_answer(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``HilbertSpaceGP`` bound: W5.4's convergence class, tightening with ``m``."""
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        if SolverKind.HILBERT not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} declares no reduced-rank spectral solver.")
        _, _, residual, variance = self._pieces(backend, pieces, CHANNEL_SIGMAS)
        reference = oracle_log_density(variance, residual)
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        errors = []
        for size in REDUCED_RANK_SIZES:
            _, noise, _, _ = self._pieces(
                backend,
                pieces,
                CHANNEL_SIGMAS,
                solver=SolverKind.HILBERT,
                basis_size=size,
                boundary_factor=REDUCED_RANK_BOX,
            )
            assert noise.route(variance) == "reduced_rank"
            assert noise.latent_size(EPOCHS.size, "reduced_rank") == 2 * size
            got = noise.log_prob(columns, variance, EPOCHS.reshape(-1, 1), COUPLING)
            errors.append(abs(got - reference))
        coarsest = errors[0]
        for size, error in zip(REDUCED_RANK_SIZES[1:], errors[1:], strict=True):
            allowed = approximation_envelope(coarsest, REDUCED_RANK_SIZES[0], size, tolerances)
            assert error <= allowed, f"error {error:.3e} at m={size} outside {allowed:.3e}"
        assert errors[-1] <= tolerances.approximation_final, errors

    @pytest.mark.parametrize("solver", [SolverKind.DENSE, SolverKind.QUASISEP])
    def test_equal_sigmas_are_bit_identical_to_the_rotated_path(
        self, backend: ConformanceBackend, solver: SolverKind
    ) -> None:
        """Equal channels keep W5.9's path: the ``(N, T)`` block and the vector agree exactly.

        Heteroscedastic *along* the grid, identical across the channels —
        the case the rotation is exact for. The block every W5.24 caller now
        passes must reach the very arithmetic W5.9's shared vector reached, so
        the comparison is ``==``, not a tolerance.
        """
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        if solver not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} declares no {solver} solver.")
        problem, noise, residual, variance = self._pieces(
            backend, pieces, SHARED_SIGMAS, solver=solver
        )
        assert noise.route(variance) == "rotated"
        columns = [residual[: EPOCHS.size], residual[EPOCHS.size :]]
        shared = SHARED_SIGMAS[0] ** 2
        block = noise.log_prob(columns, variance, EPOCHS.reshape(-1, 1), COUPLING)
        vector = noise.log_prob(columns, shared, EPOCHS.reshape(-1, 1), COUPLING)
        assert block == vector
        # And through the composed problem, which now builds the block itself.
        assert problem.log_likelihood(THETA_JOINT) == vector

    def test_unequal_sigmas_under_quasisep_are_refused_by_name(
        self, backend: ConformanceBackend
    ) -> None:
        """No Kronecker-free O(N) form: refused at composition, with the fix named."""
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        if SolverKind.QUASISEP not in backend.capabilities.solvers:
            pytest.skip(f"{backend.name} declares no QuasisepGP solver.")
        with pytest.raises(LikelihoodError, match="HilbertSpaceGP"):
            self._pieces(backend, pieces, CHANNEL_SIGMAS, solver=SolverKind.QUASISEP)

    def test_simulate_draws_correlated_heteroscedastic_channels(
        self, backend: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """``simulate(observe=True)`` under unequal sigmas, against the oracle's covariance.

        The whole ``2N x 2N`` sample covariance, cross-channel block included,
        is compared entry by entry against the oracle's, each entry within
        ``monte_carlo_sigmas`` of its own standard error — ``sqrt((S_ii S_jj +
        S_ij^2)/(n - 1))`` for a Gaussian sample covariance — so the margin is
        measured from the estimator rather than chosen. Two further claims
        keep the row from passing on a draw that had got one thing right by
        luck: the cross-channel block is genuinely there, and the two
        channels' white variances are the *different* ones each container
        carries rather than one shared value.
        """
        pieces = pieces_or_skip(backend)
        joint_or_skip(backend)
        problem, _, _, variance = self._pieces(backend, pieces, CHANNEL_SIGMAS)
        draws = 2000
        rows = np.empty((draws, 2 * EPOCHS.size))
        simulations = problem.simulate_many(draws, values=[THETA_JOINT] * draws, observe=True)
        for index, simulation in enumerate(simulations):
            assert not simulation.failed
            assert simulation.observations is not None
            rows[index] = np.concatenate(
                [
                    np.asarray(simulation.observations[name].values, dtype=float).ravel()
                    - np.asarray(simulation.predicted[name].values, dtype=float).ravel()
                    for name in ("ra", "dec")
                ]
            )
        coupling = rotation_coupling_matrix(
            COUPLING["angle"], COUPLING["log_variance_0"], COUPLING["log_variance_1"]
        )
        kernel = kernel_matrix(
            JOINT_KERNEL.family, EPOCHS, JOINT_KERNEL.amplitude, JOINT_KERNEL.length_scale
        )
        expected = coregionalised_covariance(coupling, kernel, variance)
        empirical = np.cov(rows, rowvar=False)
        diagonal = np.diag(expected)
        standard_error = np.sqrt((np.outer(diagonal, diagonal) + expected**2) / (draws - 1))
        z = np.abs(empirical - expected) / standard_error
        assert float(np.max(z)) < tolerances.monte_carlo_sigmas, (
            f"a sample-covariance entry sits {float(np.max(z)):.2f} standard errors from "
            f"B (x) K_x + blockdiag(diag sigma_t^2)"
        )
        cross = empirical[: EPOCHS.size, EPOCHS.size :]
        expected_cross = expected[: EPOCHS.size, EPOCHS.size :]
        assert float(np.max(np.abs(cross))) > 0.5 * float(np.max(np.abs(expected_cross)))
        # The white part differs between the channels, as the containers say.
        gp_part = np.diag(coupling[0, 0] * kernel), np.diag(coupling[1, 1] * kernel)
        white_ra = np.diag(empirical)[: EPOCHS.size] - gp_part[0]
        white_dec = np.diag(empirical)[EPOCHS.size :] - gp_part[1]
        assert (
            np.corrcoef(white_ra, variance[:, 0])[0, 1]
            > np.corrcoef(white_ra, variance[:, 1])[0, 1]
        )
        assert (
            np.corrcoef(white_dec, variance[:, 1])[0, 1]
            > np.corrcoef(white_dec, variance[:, 0])[0, 1]
        )

    def _pieces(
        self,
        backend: ConformanceBackend,
        pieces: AstrometryPieces,
        sigmas: np.ndarray,
        **kwargs: Any,
    ) -> tuple[FittingProblem, Any, np.ndarray, np.ndarray]:
        return TestJointChannelNoise()._pieces(backend, pieces, sigmas=sigmas, **kwargs)
