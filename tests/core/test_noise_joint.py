"""Joint noise over a tuple of channels: the shared-grid ICM (W5.9).

``likelihoods.md`` §7's new subsection, and §15's eleventh limitation lifted
for the case it is exact in. The rows here are the *contract* ones — the
declaration, the refusals, and the algebra against a materialised
``B ⊗ K_x``. The lockstep comparison across backends is
``tests/conformance/test_astrometry.py``'s, and the end-to-end calibration
claim is ``tests/examples/test_astrometry_example.py``'s.
"""

from __future__ import annotations

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    CholeskyCoupling,
    Dataset,
    DatasetCollection,
    DenseGP,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    JointGaussianProcessNoise,
    Likelihood,
    Matern32,
    QuasisepGP,
    RotationCoupling,
    Sum,
    TimeSeries,
)
from ampere.core.exceptions import DatasetError, LikelihoodError

#: An ordered one-dimensional grid, which is what the O(N) path wants.
GRID = np.array([0.0, 40.0, 95.0, 150.0, 210.0, 260.0, 330.0, 400.0, 460.0, 520.0])
SIGMA = 0.03

#: ``B`` far from a multiple of the identity: a factor of nine between the
#: eigenvalues at a position angle of 0.7 rad, so the cross-covariance is a
#: large part of the covariance rather than a perturbation of it.
COUPLING: dict[str, float] = {
    "angle": 0.7,
    "log_variance_0": float(np.log(9.0e-4)),
    "log_variance_1": float(np.log(1.0e-4)),
}


def kernel() -> Matern32:
    """``K_x``, amplitude fixed at one — ``B`` carries the scale."""
    return Matern32(1.0, 180.0, axes=("time",))


def coupling(**overrides: object) -> RotationCoupling:
    values = {**COUPLING, **overrides}
    return RotationCoupling(values["angle"], values["log_variance_0"], values["log_variance_1"])


def noise_model(*, solver: object = None, **kwargs: object) -> JointGaussianProcessNoise:
    return JointGaussianProcessNoise(
        kernel(),
        DenseGP() if solver is None else solver,
        datasets=("ra", "dec"),
        coupling=coupling(),
        **kwargs,
    )


def channel(values: np.ndarray | None = None, *, sigma: float = SIGMA, mask=None) -> TimeSeries:
    filled = np.zeros(GRID.size) if values is None else values
    return TimeSeries(
        GRID * u.day,
        filled * u.mas,
        uncertainty=np.full(GRID.size, sigma) * u.mas,
        mask=mask,
    )


def materialised(noise: JointGaussianProcessNoise, values=COUPLING) -> np.ndarray:
    """``B ⊗ K_x + I_T ⊗ diag(sigma^2)``, built the obvious way."""
    matrix = np.asarray(noise.coupling_matrix(values))
    covariance = noise.kernel.matrix(GRID[:, None], GRID[:, None], noise.kernel.resolve(values))
    return np.kron(matrix, covariance) + np.eye(2 * GRID.size) * SIGMA**2


def dense_log_density(covariance: np.ndarray, residual: np.ndarray) -> float:
    factor = np.linalg.cholesky(covariance)
    log_determinant = 2.0 * float(np.sum(np.log(np.diag(factor))))
    solved = np.linalg.solve(covariance, residual)
    return float(
        -0.5 * (float(residual @ solved) + log_determinant + residual.size * np.log(2.0 * np.pi))
    )


@pytest.fixture
def residual_block() -> np.ndarray:
    rng = np.random.default_rng(20260916)
    return rng.normal(0.0, 0.05, size=(GRID.size, 2))


class TestTheCouplingDeclaration:
    """``B``: a ``Parameterised`` that knows its own eigendecomposition."""

    def test_the_rotation_parameterisation_reassembles_its_matrix(self) -> None:
        at = {"angle": 0.0, "log_variance_0": 0.0, "log_variance_1": float(np.log(4.0))}
        assert np.allclose(coupling().matrix(at), np.diag([1.0, 4.0]))

    def test_the_eigendecomposition_is_exact_and_orthogonal(self) -> None:
        eigenvalues, rotation = coupling().eigen(COUPLING)
        assert np.allclose(rotation @ rotation.T, np.eye(2))
        rebuilt = (rotation * eigenvalues) @ rotation.T
        assert np.allclose(rebuilt, coupling().matrix(COUPLING))
        # Against numpy's own eigh, which is the thing the closed form replaces.
        numerical = np.linalg.eigvalsh(rebuilt)
        assert np.allclose(np.sort(eigenvalues), numerical)

    def test_the_cholesky_parameterisation_is_the_general_fallback(self) -> None:
        general = CholeskyCoupling(3, log_diagonal=[0.0, -0.5, 0.2], off_diagonal=[0.3, -0.1, 0.4])
        resolved = general.resolved({})
        matrix = general.matrix(resolved)
        assert matrix.shape == (3, 3)
        assert np.allclose(matrix, matrix.T)
        assert np.all(np.linalg.eigvalsh(matrix) > 0.0)
        lower = np.array([[1.0, 0.0, 0.0], [0.3, np.exp(-0.5), 0.0], [-0.1, 0.4, np.exp(0.2)]])
        assert np.allclose(matrix, lower @ lower.T)

    def test_a_coupling_needs_at_least_two_channels(self) -> None:
        with pytest.raises(LikelihoodError, match="two or more channels"):
            CholeskyCoupling(1, log_diagonal=[0.0])

    def test_a_couplings_parameters_take_their_priors_own_bijection(self) -> None:
        """Not ``Log``: a log-variance is negative for every ordinary value.

        The defect this row pins: :func:`ampere.core.likelihood._as_hyperparameter`
        imposes ``Log`` because a *kernel* hyperparameter is positive by
        construction, and a coupling's are not. With ``Log`` the whole
        unconstrained vector came back ``NaN`` while the density itself stayed
        finite, which is the quietest way a fit can fail.
        """
        declared = RotationCoupling(st.uniform(0.0, np.pi), st.norm(-8.0, 2.0), st.norm(-8.0, 2.0))
        unconstrained = declared.parameters.unconstrain(
            {"angle": 0.4, "log_variance_0": -8.0, "log_variance_1": -9.0}
        )
        assert np.all(np.isfinite(unconstrained))


class TestTheJointDensity:
    """The rotated ``T`` solves, against the materialised Kronecker covariance."""

    def test_the_joint_density_matches_the_dense_kronecker_solve(
        self, residual_block: np.ndarray
    ) -> None:
        noise = noise_model()
        got = noise.log_prob(residual_block, np.full(GRID.size, SIGMA**2), GRID[:, None], COUPLING)
        stacked = np.concatenate([residual_block[:, 0], residual_block[:, 1]])
        expected = dense_log_density(materialised(noise), stacked)
        assert got == pytest.approx(expected, rel=1e-10)

    def test_the_quasiseparable_solver_gives_the_same_number(
        self, residual_block: np.ndarray
    ) -> None:
        dense = noise_model()
        quasisep = noise_model(solver=QuasisepGP())
        arguments = (residual_block, np.full(GRID.size, SIGMA**2), GRID[:, None], COUPLING)
        assert dense.log_prob(*arguments) == pytest.approx(quasisep.log_prob(*arguments), abs=1e-8)

    def test_an_identity_coupling_is_bit_identical_to_two_scalar_gps(
        self, residual_block: np.ndarray
    ) -> None:
        """``B = I``: every operation is the scalar path's, in the same order."""
        unit = JointGaussianProcessNoise(
            Matern32(0.02, 180.0, axes=("time",)),
            DenseGP(),
            datasets=("ra", "dec"),
            coupling=RotationCoupling(0.0, 0.0, 0.0),
        )
        identity = {"angle": 0.0, "log_variance_0": 0.0, "log_variance_1": 0.0}
        got = unit.log_prob(residual_block, np.full(GRID.size, SIGMA**2), GRID[:, None], identity)
        scalar = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(Matern32(0.02, 180.0, axes=("time",)), DenseGP()),
        )
        separate = sum(
            scalar.log_prob(channel(), channel(residual_block[:, index]), values={})
            for index in range(2)
        )
        assert got == separate

    def test_a_non_positive_coupling_is_refused_by_name(self, residual_block: np.ndarray) -> None:
        noise = noise_model()
        with pytest.raises(LikelihoodError, match="not a covariance"):
            noise.log_prob(
                residual_block,
                np.full(GRID.size, SIGMA**2),
                GRID[:, None],
                {**COUPLING, "log_variance_0": float("nan")},
            )

    def test_the_pointwise_terms_are_one_column_per_rotated_output(
        self, residual_block: np.ndarray
    ) -> None:
        noise = noise_model()
        terms = noise.pointwise_log_prob(
            residual_block, np.full(GRID.size, SIGMA**2), GRID[:, None], COUPLING
        )
        assert terms.shape == (GRID.size, 2)
        assert np.all(np.isfinite(terms))


class TestTheDrawIsCorrelated:
    """The generative half, against the same covariance the density scores."""

    def test_the_sample_covariance_matches_b_kron_k(self) -> None:
        noise = noise_model()
        rng = np.random.default_rng(3)
        draws = 6000
        zeros = [np.zeros(GRID.size), np.zeros(GRID.size)]
        rows = np.empty((draws, 2 * GRID.size))
        for index in range(draws):
            block = noise.sample(zeros, np.full(GRID.size, SIGMA**2), GRID[:, None], COUPLING, rng)
            rows[index] = np.concatenate([block[:, 0], block[:, 1]])
        expected = materialised(noise)
        empirical = np.cov(rows, rowvar=False)
        scale = float(np.max(np.abs(expected)))
        assert float(np.max(np.abs(empirical - expected))) / scale < 0.12
        cross = empirical[: GRID.size, GRID.size :]
        assert float(np.max(np.abs(cross))) > 0.2 * scale


class TestTheRefusals:
    """Every rule that can be checked without a prediction, checked loudly."""

    def test_a_joint_noise_model_is_refused_inside_a_likelihood(self) -> None:
        with pytest.raises(LikelihoodError, match="JOINT = True"):
            Likelihood(GaussianFamily(), noise_model())

    def test_a_free_kernel_amplitude_beside_a_free_coupling_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="same degree of freedom twice"):
            JointGaussianProcessNoise(
                Matern32(st.loguniform(1e-3, 1e-1), 180.0, axes=("time",)),
                DenseGP(),
                datasets=("ra", "dec"),
                coupling=RotationCoupling(
                    st.uniform(0.0, np.pi), st.norm(-8.0, 1.0), st.norm(-8.0, 1.0)
                ),
            )

    def test_one_fixed_amplitude_in_a_sum_pins_the_scale(self) -> None:
        """A ``Sum``'s relative amplitudes are identified, so this composes."""
        composite = Sum(
            Matern32(0.02, 180.0, axes=("time",)),
            Matern32(st.loguniform(1e-3, 1e-1), 40.0, axes=("time",)),
        )
        built = JointGaussianProcessNoise(
            composite,
            DenseGP(),
            datasets=("ra", "dec"),
            coupling=RotationCoupling(
                st.uniform(0.0, np.pi), st.norm(-8.0, 1.0), st.norm(-8.0, 1.0)
            ),
        )
        assert "angle" in built.parameters

    def test_the_channel_count_must_match_the_couplings(self) -> None:
        with pytest.raises(LikelihoodError, match="2x2"):
            JointGaussianProcessNoise(
                kernel(), DenseGP(), datasets=("a", "b", "c"), coupling=coupling()
            )

    def test_repeated_dataset_labels_are_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="must be distinct"):
            JointGaussianProcessNoise(
                kernel(), DenseGP(), datasets=("ra", "ra"), coupling=coupling()
            )

    def test_noise_params_refuses_to_describe_one_channel(self) -> None:
        with pytest.raises(LikelihoodError, match="cross-covariance between the channels"):
            noise_model().noise_params(channel(), np.ones(GRID.size, dtype=bool), COUPLING)


class TestTheCollectionBinding:
    """``DatasetCollection(joint=…)``: the composition-time rules."""

    def _collection(self, **kwargs) -> DatasetCollection:
        datasets = {
            "ra": Dataset(channel(), likelihood=Likelihood(GaussianFamily()), label="ra"),
            "dec": Dataset(channel(), likelihood=Likelihood(GaussianFamily()), label="dec"),
        }
        datasets.update(kwargs.pop("datasets", {}))
        return DatasetCollection(datasets, **kwargs)

    def test_a_group_is_one_component_and_one_contribution_key(self) -> None:
        collection = self._collection(joint={"astrom": noise_model()})
        assert "astrom" in collection.components()
        assert collection.contribution_labels() == ("astrom",)
        assert collection.group_of("ra") == "astrom"
        assert collection.group_of("dec") == "astrom"

    def test_a_collection_with_no_groups_is_unchanged(self) -> None:
        collection = self._collection()
        assert collection.contribution_labels() == ("ra", "dec")
        assert collection.group_of("ra") is None

    def test_a_group_label_may_not_collide_with_a_dataset(self) -> None:
        with pytest.raises(DatasetError, match="also a dataset label"):
            self._collection(joint={"ra": noise_model()})

    def test_an_unknown_member_is_refused(self) -> None:
        stray = JointGaussianProcessNoise(
            kernel(), DenseGP(), datasets=("ra", "parallax"), coupling=coupling()
        )
        with pytest.raises(DatasetError, match="not in this collection"):
            self._collection(joint={"astrom": stray})

    def test_channels_with_different_uncertainties_are_refused(self) -> None:
        other = Dataset(channel(sigma=0.05), likelihood=Likelihood(GaussianFamily()), label="dec")
        with pytest.raises(LikelihoodError, match="different per-sample uncertainties"):
            self._collection(datasets={"dec": other}, joint={"astrom": noise_model()})

    def test_channels_on_different_grids_are_refused(self) -> None:
        shifted = TimeSeries(
            (GRID + 1e-9) * u.day,
            np.zeros(GRID.size) * u.mas,
            uncertainty=np.full(GRID.size, SIGMA) * u.mas,
        )
        other = Dataset(shifted, likelihood=Likelihood(GaussianFamily()), label="dec")
        with pytest.raises(LikelihoodError, match="not identical"):
            self._collection(datasets={"dec": other}, joint={"astrom": noise_model()})

    def test_channels_masking_different_samples_are_refused(self) -> None:
        mask = np.zeros(GRID.size, dtype=bool)
        mask[2] = True
        other = Dataset(channel(mask=mask), likelihood=Likelihood(GaussianFamily()), label="dec")
        with pytest.raises(LikelihoodError, match="mask different samples"):
            self._collection(datasets={"dec": other}, joint={"astrom": noise_model()})

    def test_a_member_carrying_its_own_noise_parameters_is_refused(self) -> None:
        other = Dataset(
            channel(),
            likelihood=Likelihood(GaussianFamily(), IndependentNoise(jitter=st.halfnorm(0, 0.1))),
            label="dec",
        )
        with pytest.raises(LikelihoodError, match="declares parameters"):
            self._collection(datasets={"dec": other}, joint={"astrom": noise_model()})

    def test_a_member_carrying_its_own_gp_is_refused(self) -> None:
        other = Dataset(
            channel(),
            likelihood=Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(Matern32(0.02, 180.0, axes=("time",)), DenseGP()),
            ),
            label="dec",
        )
        with pytest.raises(LikelihoodError, match="counted twice"):
            self._collection(datasets={"dec": other}, joint={"astrom": noise_model()})


class TestThroughAProblem:
    """The composed problem scores and simulates through the group."""

    def _problem(self, *, joint: bool = True) -> FittingProblem:
        from ampere.backends.reference import EpochSample, ReflexOrbit

        from ampere.core import Instrument

        orbit = {
            "pmra": 1.2,
            "pmdec": -0.6,
            "period": 400.0,
            "phase": 0.7,
            "amp_ra": 0.6,
            "amp_dec": 0.35,
        }
        instruments = {
            name: Instrument([EpochSample(GRID)], channel=name, label=f"astrom_{name}")
            for name in ("ra", "dec")
        }
        truth = ReflexOrbit(GRID, **orbit)
        from ampere.core import negotiate

        compiled = truth.compile_for(negotiate(list(instruments.values())))
        evaluated = compiled.evaluate()
        observed = {
            name: channel(np.asarray(step(evaluated).values).ravel())
            for name, step in instruments.items()
        }
        model = ReflexOrbit(
            GRID,
            pmra=st.norm(0.0, 5.0),
            pmdec=st.norm(0.0, 5.0),
            period=orbit["period"],
            phase=orbit["phase"],
            amp_ra=orbit["amp_ra"],
            amp_dec=orbit["amp_dec"],
        )
        datasets = DatasetCollection(
            {
                name: Dataset(
                    observed[name],
                    instruments[name],
                    likelihood=Likelihood(GaussianFamily()),
                    label=name,
                )
                for name in ("ra", "dec")
            },
            joint={"astrom": self._noise()} if joint else None,
        )
        return FittingProblem(model, datasets, seed=20260916)

    def _noise(self) -> JointGaussianProcessNoise:
        """The group with a **fitted** coupling, so its parameters join the space."""
        return JointGaussianProcessNoise(
            kernel(),
            QuasisepGP(),
            datasets=("ra", "dec"),
            coupling=RotationCoupling(
                st.uniform(0.0, np.pi), st.norm(-8.0, 1.0), st.norm(-8.0, 1.0)
            ),
        )

    def _theta(self) -> dict[str, float]:
        return {
            "model.pmra": 1.2,
            "model.pmdec": -0.6,
            "astrom.angle": COUPLING["angle"],
            "astrom.log_variance_0": COUPLING["log_variance_0"],
            "astrom.log_variance_1": COUPLING["log_variance_1"],
        }

    def test_the_groups_parameters_join_the_problems_space(self) -> None:
        problem = self._problem()
        assert set(problem.parameters.free_names) == set(self._theta())

    def test_the_contribution_is_keyed_by_the_group(self) -> None:
        problem = self._problem()
        evaluation = problem.evaluate(self._theta())
        assert set(evaluation.contributions) == {"astrom"}
        assert np.isfinite(evaluation.log_likelihood)
        assert evaluation.log_likelihood == pytest.approx(evaluation.contributions["astrom"])

    def test_simulate_draws_both_channels_together(self) -> None:
        problem = self._problem()
        simulation = problem.simulate(self._theta(), observe=True)
        assert not simulation.failed
        assert simulation.observations is not None
        assert set(simulation.observations) == {"ra", "dec"}
        for name in ("ra", "dec"):
            assert not np.array_equal(
                simulation.observations[name].values, simulation.predicted[name].values
            )

    def test_the_density_peaks_at_the_truth(self) -> None:
        problem = self._problem()
        theta = self._theta()
        assert problem.log_prob(theta) > problem.log_prob({**theta, "model.pmra": 6.0})
