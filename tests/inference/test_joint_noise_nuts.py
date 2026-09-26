"""NUTS over a joint noise group's coupling and kernel (W5.9).

``likelihoods.md`` §7's shared-grid intrinsic coregionalisation model is the
first correlated process in this contract that spans datasets, and the whole
of its value is in ``B`` — the matrix saying how the channels move together.
So the claim this file makes is the one that matters: a **gradient reaches
``B``'s parameters**, and a NUTS chain over them and the kernel's length scale
runs and recovers the injected coupling.

The failure this guards against is the one W2.14 found for the latent path and
that ``likelihoods.md`` §15's sixth limitation records at length: a lowering
that omits a term leaves the density *exactly flat* in the parameters that term
carries, so the chain runs, converges, and reports the prior. Here that would
be a coupling sampled against nothing at all.

One column per installed differentiable backend, in
``tests/inference/test_nuts.py``'s own shape; no test body names a library.
Budgets are deliberately small — this is a claim about gradients and
composition, not about posterior accuracy, which
``tests/examples/test_astrometry_example.py`` measures by simulation-based
calibration.
"""

from __future__ import annotations

import dataclasses
import importlib
from typing import Any

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    Instrument,
    Likelihood,
    RotationCoupling,
    negotiate,
)
from ampere.inference import NUTSEngine

#: Epochs. Few, because a NUTS chain is run per backend here and the claim is
#: about the gradient rather than about N.
EPOCHS = np.array([0.0, 60.0, 130.0, 200.0, 280.0, 350.0, 420.0, 500.0, 570.0, 650.0])
SIGMA = 0.03
SEED = 20260916

ORBIT = {
    "pmra": 1.2,
    "pmdec": -0.6,
    "period": 400.0,
    "phase": 0.7,
    "amp_ra": 0.6,
    "amp_dec": 0.35,
}

#: The injected coupling: an error ellipse a factor of nine from circular, at a
#: position angle of 0.7 rad.
COUPLING = {
    "angle": 0.7,
    "log_variance_0": float(np.log(9.0e-4)),
    "log_variance_1": float(np.log(1.0e-4)),
}
LENGTH_SCALE = 180.0

#: W5.24's heteroscedastic case: every channel draws its own sigma per epoch,
#: a log-uniform factor of up to two either side of :data:`SIGMA`, so the two
#: channels' diagonals differ everywhere and the group takes the reduced-rank
#: route.
SIGMA_SPREAD = float(np.log(2.0))

#: The reduced-rank solver the heteroscedastic rows bind. The box has to reach
#: well beyond the data for a length scale this long against a 650-day grid
#: (the approximation is poor within about one length scale of the boundary),
#: so the factor is three; see the parity row for what that buys.
BASIS_SIZE = 48
BOUNDARY_FACTOR = 3.0


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so a test body never names a library."""

    name: str
    module: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS,
    reason="no differentiable backend installed; NUTS needs a registered realisation",
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


def observed_channel(module: Any, values: np.ndarray, sigma: np.ndarray | None = None) -> Any:
    from ampere.core import TimeSeries

    import astropy.units as u

    del module
    return TimeSeries(
        EPOCHS * u.day,
        np.asarray(values, dtype=float) * u.mas,
        uncertainty=(np.full(EPOCHS.size, SIGMA) if sigma is None else sigma) * u.mas,
    )


def channel_sigmas() -> np.ndarray:
    """``(2, N)`` per-channel, per-epoch sigmas for the heteroscedastic rows (W5.24)."""
    rng = np.random.default_rng(SEED + 2)
    return SIGMA * np.exp(rng.uniform(-SIGMA_SPREAD, SIGMA_SPREAD, size=(2, EPOCHS.size)))


def build_problem(
    kit: Kit, *, fitted: bool = True, heteroscedastic: bool = False, dense: bool = False
) -> FittingProblem:
    """One reflex orbit, two channels, one joint noise group over both.

    The coupling's three parameters and the kernel's length scale are free; the
    kernel's **amplitude is fixed**, because ``B (x) (a^2 K) = (a^2 B) (x) K``
    makes a free amplitude and a free ``B`` the same degree of freedom twice
    and ``JointGaussianProcessNoise`` refuses the pair by name.

    ``heteroscedastic=True`` is W5.24's case: each channel carries its own
    sigma per epoch (:func:`channel_sigmas`), the white noise is drawn at
    those, and the group binds ``HilbertSpaceGP`` so it takes the
    reduced-rank route -- the NUTS-friendly one, with a ``T*m`` whitened
    block rather than the ``O((TN)^3)`` dense solve. ``dense=True`` binds
    ``DenseGP`` instead, which is the exact dense route.
    """
    module = kit.module
    instruments = {
        name: Instrument([module.EpochSample(EPOCHS)], channel=name, label=f"astrom_{name}")
        for name in ("ra", "dec")
    }
    truth = module.ReflexOrbit(EPOCHS, **ORBIT)
    compiled = truth.compile_for(negotiate(list(instruments.values())))
    evaluated = compiled.evaluate()

    # The injected data: the orbit, plus a draw from B (x) K_x, plus white
    # noise. Drawn here in numpy from the materialised Kronecker covariance, so
    # the data are generated independently of the code that scores them.
    kernel_for_truth = module.Matern32(1.0, LENGTH_SCALE, axes=("time",))
    covariance = np.asarray(
        kit.module.to_numpy(
            kernel_for_truth.matrix(EPOCHS[:, None], EPOCHS[:, None], kernel_for_truth.resolve({}))
        )
        if hasattr(kit.module, "to_numpy")
        else kernel_for_truth.matrix(
            EPOCHS[:, None], EPOCHS[:, None], kernel_for_truth.resolve({})
        ),
        dtype=float,
    )
    matrix = np.asarray(RotationCoupling(0.0, 0.0, 0.0).matrix(COUPLING), dtype=float)
    stacked = np.kron(matrix, covariance)
    stacked = stacked + np.eye(stacked.shape[0]) * 1e-10 * float(np.max(np.diag(stacked)))
    rng = np.random.default_rng(SEED)
    systematic = rng.multivariate_normal(np.zeros(stacked.shape[0]), stacked).reshape(2, -1)

    sigmas = channel_sigmas() if heteroscedastic else np.full((2, EPOCHS.size), SIGMA)
    observed = {}
    for index, name in enumerate(("ra", "dec")):
        noiseless = np.asarray(
            kit.module.to_numpy(instruments[name](evaluated).values)
            if hasattr(kit.module, "to_numpy")
            else instruments[name](evaluated).values,
            dtype=float,
        ).ravel()
        observed[name] = observed_channel(
            module,
            noiseless + systematic[index] + rng.normal(0.0, sigmas[index], noiseless.shape),
            sigmas[index] if heteroscedastic else None,
        )

    model = module.ReflexOrbit(
        EPOCHS,
        pmra=st.norm(0.0, 5.0),
        pmdec=st.norm(0.0, 5.0),
        period=ORBIT["period"],
        phase=ORBIT["phase"],
        amp_ra=ORBIT["amp_ra"],
        amp_dec=ORBIT["amp_dec"],
    )
    coupling = (
        RotationCoupling(st.uniform(0.0, np.pi), st.norm(-8.0, 1.0), st.norm(-8.0, 1.0))
        if fitted
        else RotationCoupling(
            COUPLING["angle"], COUPLING["log_variance_0"], COUPLING["log_variance_1"]
        )
    )
    noise = module.JointGaussianProcessNoise(
        module.Matern32(
            1.0, st.loguniform(50.0, 500.0) if fitted else LENGTH_SCALE, axes=("time",)
        ),
        module.DenseGP()
        if dense
        else module.HilbertSpaceGP(basis_size=BASIS_SIZE, boundary_factor=BOUNDARY_FACTOR)
        if heteroscedastic
        else module.QuasisepGP(),
        datasets=("ra", "dec"),
        coupling=coupling,
    )
    datasets = DatasetCollection(
        {
            name: Dataset(
                observed[name],
                instruments[name],
                likelihood=Likelihood(GaussianFamily(), module.IndependentNoise()),
                label=name,
            )
            for name in ("ra", "dec")
        },
        joint={"astrom": noise},
    )
    return FittingProblem(model, datasets, seed=SEED)


class TestTheGroupIsDifferentiable:
    """The gradient reaches ``B``. Everything else here depends on that."""

    def test_the_lowered_density_agrees_with_the_contract_path(self, kit: Kit) -> None:
        problem = build_problem(kit)
        lowered = kit.module.lower_problem(problem)
        theta = {
            "model.pmra": ORBIT["pmra"],
            "model.pmdec": ORBIT["pmdec"],
            "astrom.angle": COUPLING["angle"],
            "astrom.log_variance_0": COUPLING["log_variance_0"],
            "astrom.log_variance_1": COUPLING["log_variance_1"],
            "astrom.length_scale": LENGTH_SCALE,
        }
        vector = problem.parameters.pack(problem.parameters.complete(theta))
        native = float(np.asarray(lowered.log_likelihood(vector)))
        assert native == pytest.approx(problem.log_likelihood(theta), abs=1e-7)

    def test_the_gradient_is_finite_and_non_zero_in_every_coupling_parameter(
        self, kit: Kit
    ) -> None:
        problem = build_problem(kit)
        lowered = kit.module.lower_problem(problem)
        theta = {
            "model.pmra": ORBIT["pmra"],
            "model.pmdec": ORBIT["pmdec"],
            "astrom.angle": COUPLING["angle"],
            "astrom.log_variance_0": COUPLING["log_variance_0"] + 0.4,
            "astrom.log_variance_1": COUPLING["log_variance_1"] - 0.3,
            "astrom.length_scale": LENGTH_SCALE * 1.3,
        }
        unconstrained = problem.unconstrain(problem.parameters.complete(theta))
        gradient = _gradient(kit, lowered, unconstrained)
        names = problem.free_labels()
        assert np.all(np.isfinite(gradient))
        for index, name in enumerate(names):
            if name.startswith("astrom."):
                assert gradient[index] != 0.0, f"the density is flat in {name}"


def _truth() -> dict[str, float]:
    return {
        "model.pmra": ORBIT["pmra"],
        "model.pmdec": ORBIT["pmdec"],
        "astrom.angle": COUPLING["angle"],
        "astrom.log_variance_0": COUPLING["log_variance_0"],
        "astrom.log_variance_1": COUPLING["log_variance_1"],
        "astrom.length_scale": LENGTH_SCALE,
    }


class TestTheReducedRankRouteIsDifferentiable:
    """W5.24: unequal per-channel sigmas, through the reduced-rank route natively."""

    def test_the_lowered_group_takes_the_reduced_rank_route(self, kit: Kit) -> None:
        problem = build_problem(kit, heteroscedastic=True)
        lowered = kit.module.lower_problem(problem)
        (group,) = lowered._joint
        assert group.route == "reduced_rank"
        noise = problem.datasets.joint["astrom"]
        assert noise.latent_size(EPOCHS.size, "reduced_rank") == 2 * BASIS_SIZE

    @pytest.mark.parametrize("dense", [False, True], ids=["reduced_rank", "dense"])
    def test_the_lowered_density_agrees_with_the_contract_path(self, kit: Kit, dense: bool) -> None:
        problem = build_problem(kit, heteroscedastic=True, dense=dense)
        lowered = kit.module.lower_problem(problem)
        (group,) = lowered._joint
        assert group.route == ("dense" if dense else "reduced_rank")
        theta = _truth()
        vector = problem.parameters.pack(problem.parameters.complete(theta))
        native = float(np.asarray(lowered.log_likelihood(vector)))
        assert native == pytest.approx(problem.log_likelihood(theta), abs=1e-7)

    def test_the_gradient_reaches_every_coupling_parameter_and_the_kernel(self, kit: Kit) -> None:
        problem = build_problem(kit, heteroscedastic=True)
        lowered = kit.module.lower_problem(problem)
        theta = {
            **_truth(),
            "astrom.log_variance_0": COUPLING["log_variance_0"] + 0.4,
            "astrom.log_variance_1": COUPLING["log_variance_1"] - 0.3,
            "astrom.length_scale": LENGTH_SCALE * 1.3,
        }
        unconstrained = problem.unconstrain(problem.parameters.complete(theta))
        gradient = _gradient(kit, lowered, unconstrained)
        assert np.all(np.isfinite(gradient))
        for index, name in enumerate(problem.free_labels()):
            if name.startswith("astrom."):
                assert gradient[index] != 0.0, f"the density is flat in {name}"


class TestNUTSOverTheCoupling:
    """A short chain over ``B``'s three parameters and the kernel's length scale."""

    def test_nuts_runs_and_recovers_the_coupling(self, kit: Kit) -> None:
        problem = build_problem(kit)
        run = NUTSEngine(problem).run(150, warmup=150, chains=2, progress=False)
        posterior = run["posterior"].dataset
        assert set(posterior.data_vars) >= {
            "astrom.angle",
            "astrom.log_variance_0",
            "astrom.log_variance_1",
            "astrom.length_scale",
        }
        # The coupling's *scale* is what a fit of one realisation can be held
        # to: the eigenvalue labelling is exchangeable with the angle (see
        # RotationCoupling), so the total variance of B -- its trace, which is
        # invariant under that relabelling -- is the statistic.
        draws = {
            name: np.asarray(posterior[name], dtype=float).ravel()
            for name in ("astrom.log_variance_0", "astrom.log_variance_1")
        }
        trace = np.exp(draws["astrom.log_variance_0"]) + np.exp(draws["astrom.log_variance_1"])
        truth = np.exp(COUPLING["log_variance_0"]) + np.exp(COUPLING["log_variance_1"])
        lower, upper = np.percentile(trace, [2.5, 97.5])
        assert lower <= truth <= upper, (
            f"the injected trace of B, {truth:.3g}, is outside the fitted 95% interval "
            f"[{lower:.3g}, {upper:.3g}]."
        )
        assert np.all(np.isfinite(np.asarray(posterior["model.pmra"], dtype=float)))


class TestNUTSOverAHeteroscedasticCoupling:
    """W5.24: the same chain on unequal per-channel sigmas, via the reduced-rank route."""

    def test_nuts_recovers_b_and_the_length_scale(self, kit: Kit) -> None:
        problem = build_problem(kit, heteroscedastic=True)
        run = NUTSEngine(problem).run(150, warmup=150, chains=2, progress=False)
        posterior = run["posterior"].dataset
        draws = {
            name: np.asarray(posterior[name], dtype=float).ravel()
            for name in ("astrom.log_variance_0", "astrom.log_variance_1", "astrom.length_scale")
        }
        # B's scale, as the equal-sigma row holds it (the eigenvalue labelling
        # is exchangeable with the angle, so the trace is the invariant).
        trace = np.exp(draws["astrom.log_variance_0"]) + np.exp(draws["astrom.log_variance_1"])
        truth = np.exp(COUPLING["log_variance_0"]) + np.exp(COUPLING["log_variance_1"])
        lower, upper = np.percentile(trace, [2.5, 97.5])
        assert lower <= truth <= upper, (
            f"the injected trace of B, {truth:.3g}, is outside the fitted 95% interval "
            f"[{lower:.3g}, {upper:.3g}]."
        )
        lower, upper = np.percentile(draws["astrom.length_scale"], [2.5, 97.5])
        assert lower <= LENGTH_SCALE <= upper, (
            f"the injected length scale, {LENGTH_SCALE}, is outside the fitted 95% interval "
            f"[{lower:.3g}, {upper:.3g}]."
        )
        assert np.all(np.isfinite(np.asarray(posterior["model.pmra"], dtype=float)))


def _gradient(kit: Kit, lowered: Any, unconstrained: np.ndarray) -> np.ndarray:
    """``d log p / d y`` through whichever autodiff this backend uses."""
    if kit.name == "jax":
        import jax

        return np.asarray(jax.grad(lowered.log_prob_unconstrained)(unconstrained), dtype=float)
    import torch

    y = torch.as_tensor(np.asarray(unconstrained), dtype=torch.float64).requires_grad_(True)
    lowered.log_prob_unconstrained(y).backward()
    assert y.grad is not None
    return np.asarray(y.grad.detach().numpy(), dtype=float)
