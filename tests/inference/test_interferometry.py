"""The interferometric modality under every engine — Phase 4's proof (W4.3).

``tests/conformance/test_interferometry.py`` says the steps and the models are
right; ``tests/backends/test_native_interferometry.py`` says their gradients
are. This file is the claim those two exist to support: that a sky model fitted
jointly to **visibilities and closure phases** runs under every engine ampere
ships, on both modern backends, and recovers the source it was made from.

What each section is for, in the order the work item puts them:

1. **The composition.** Two instruments bind one ``sky`` channel by name,
   negotiation unions their requirements into one image grid, and the model is
   built once per draw. That is ``DEVELOPMENT_PLAN.md`` §4.3's design, and the
   modality was chosen because it exercises a kind-changing,
   coordinate-changing step, complex data, wrapped angles and two datasets on
   one channel at the same time.
2. **NUTS recovers the binary**, on both backends, at the tolerance
   ``tests/inference/test_nuts.py`` fixed for W2: the posterior mean within
   four of its own standard deviations of the injected value. A derivative of a
   discrete Fourier transform with respect to a separation is what makes that
   possible, and it is the one claim that cannot be made on the reference path.
3. **The flexible likelihood on the visibilities** — the circular complex
   Gaussian process of W4.2, which asked this item for exactly this row: a
   ``GaussianProcessNoise`` over ``(u, v)`` on a *complex* container, realised
   and sampled. Two of W4.2's native-path fixes are guarded here rather than
   asserted: the kernel must come from ``noise.kernel_for(observed)`` (so an
   ``axes=`` selector means something natively) and the coordinates must be the
   whole ``(n, d)`` block (so ``axes=("u", "v")`` is not silently the ``u``
   column alone).
4. **VI**, which is the other gradient-based engine and the cheap one.
5. **The batched forward path**, against the loop.
6. **SBI**: an NPE fit of the two-dataset problem at a smoke budget, through
   W3.3's coordinate-value-mask encoding — its first complex customer and its
   first five-axis one — passing the calibration check.
7. **The artefact cache**, keyed on a problem with two observed containers.
8. **The encoding itself**, pinned: what a ``VisibilitySet``'s three axes and a
   ``ClosurePhases``'s five actually become.

Parametrised over the backends installed here, in ``test_nuts.py``'s shape and
for its reason. Budgets are small and seeds fixed: this belongs in the per-PR
gate beside the other new-namespace suites.
"""

from __future__ import annotations

import dataclasses
import importlib
import math
import pathlib
import sys
from typing import Any, ClassVar

import numpy as np
import pytest
import scipy.stats as st

from ampere.core import LoweringError
from ampere.core.encoding import (
    EncodingError,
    EncodingLayout,
    decode,
    encode,
    encode_observations,
    unpack,
)
from ampere.core.exceptions import SchemaError
from ampere.core.simulate import SimulationBatch
from ampere.inference import EngineError, NUTSEngine, VIEngine

# ``tests/backends/interferometry_fixtures.py`` builds the synthetic source this
# file fits, and it lives there because the backend suites are its other
# consumer. ``tests/`` is not a package, so pytest puts a test file's own
# directory on ``sys.path`` and not its siblings'; the insert is explicit for the
# reason ``tests/m2/conftest.py`` gives for its own — a suite that worked because
# of somebody's editable install's path hook would break the first time it ran
# against a wheel.
_FIXTURES = pathlib.Path(__file__).resolve().parents[1] / "backends"
if str(_FIXTURES) not in sys.path:
    sys.path.insert(0, str(_FIXTURES))

import interferometry_fixtures as source  # noqa: E402


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's pieces, by name, so no test body names a library."""

    name: str
    backend: Any
    itf: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            backend = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns the flag on.
            backend.configure_x64()
        found.append(
            Kit(name, backend, importlib.import_module(f"ampere.backends.{name}.interferometry"))
        )
    return found


KITS = _installed_kits()

#: The reason a backend-dependent row skips. Not a module-level ``pytestmark``,
#: deliberately: the encoding and the artefact-cache sections are about
#: ``ampere.core`` and are composed on the **reference** backend, so they run in
#: the ``dev`` environment too — which is where a regression in them would
#: otherwise go unseen until a backend environment ran.
_NO_BACKEND = "no differentiable backend installed; this row needs a registered realisation"


@pytest.fixture(scope="module", params=[kit.name for kit in KITS] or ["none"])
def kit(request: Any) -> Kit:
    if not KITS:
        pytest.skip(_NO_BACKEND)
    return next(found for found in KITS if found.name == request.param)


@pytest.fixture(scope="module")
def reference_problem() -> Any:
    """The same two-dataset problem on the numpy path.

    ``ampere.core`` supplies the noise models, the kernel and the solver the
    reference path composes with (they declare ``BACKEND = "reference"``), so the
    fixtures module takes it as the "backend" argument unchanged. What this
    fixture is for is the rows that are about ``ampere.core`` rather than about a
    backend — the encoding's packing and the artefact key — which should not wait
    on an extra to be installed.
    """
    import ampere.core as core
    from ampere.backends.reference import interferometry as reference

    return source.two_dataset_problem(core, reference)


@pytest.fixture(scope="module")
def problem(kit: Kit) -> Any:
    """The two-dataset problem with independent noise. Module-scoped: it is read."""
    return source.two_dataset_problem(kit.backend, kit.itf)


@pytest.fixture(scope="module")
def layout(reference_problem: Any) -> EncodingLayout:
    """The problem's own packing: a complex three-axis kind beside a real five-axis one."""
    return EncodingLayout.from_datasets(reference_problem.datasets)


@pytest.fixture(scope="module")
def gp_problem(kit: Kit) -> Any:
    """The same data under the flexible likelihood (W4.2)."""
    return source.two_dataset_problem(kit.backend, kit.itf, gp=True)


# ---------------------------------------------------------------------------
# 1. The composition
# ---------------------------------------------------------------------------


class TestTheComposition:
    """Two instruments, one sky, one negotiated grid, one model evaluation."""

    def test_both_instruments_bind_the_same_channel(self, problem: Any) -> None:
        """Design horizon (h)'s obligation: the pairing stays visible."""
        channel = problem.requirements["model"]["sky"]
        assert set(channel.sources) == {"vis", "t3"}
        assert sorted(channel.axes) == ["x", "y"]
        assert {type(problem.datasets[name].observed).__name__ for name in ("vis", "t3")} == {
            "VisibilitySet",
            "ClosurePhases",
        }

    def test_the_problem_is_the_one_the_item_describes(self, problem: Any, kit: Kit) -> None:
        assert problem.backend == kit.name
        assert problem.free_size == 2
        assert problem.differentiable is True
        assert problem.batchable is True
        assert problem.parameters.free_names == ("model.separation", "model.flux_ratio")

    def test_the_two_families_are_the_two_the_observables_need(self, problem: Any) -> None:
        """A complex container scored by a circular complex Gaussian, a wrapped
        one by a von Mises — and a Gaussian on an unwrapped phase residual would
        charge a 2-degree error straddling the branch cut as a 358-degree one."""
        assert problem.datasets["vis"].likelihood.family.NAME == "complex_gaussian"
        assert problem.datasets["t3"].likelihood.family.NAME == "von_mises"

    def test_the_realised_density_agrees_with_the_contract_path(
        self, problem: Any, kit: Kit
    ) -> None:
        """``inference.md`` §10a: a realisation is checked against the numpy path."""
        lowered = kit.backend.lower_problem(problem)
        for separation, ratio in ((12.0, 0.42), (9.0, 0.2), (13.5, 0.6)):
            theta = problem.unconstrain({"model.separation": separation, "model.flux_ratio": ratio})
            assert float(lowered.log_prob_unconstrained(theta)) == pytest.approx(
                problem.log_prob_unconstrained(theta), rel=1e-11
            )

    def test_the_composed_likelihood_peaks_at_the_truth(self, problem: Any) -> None:
        at_truth = problem.log_prob(source.TRUTH)
        displaced = problem.log_prob({**source.TRUTH, "model.separation": 13.0})
        assert np.isfinite(at_truth)
        assert at_truth > displaced


# ---------------------------------------------------------------------------
# 2. NUTS recovers the binary
# ---------------------------------------------------------------------------


#: The NUTS budget. Small enough for a per-PR gate (about 25 s on torch, 12 s on
#: jax) and long enough that a two-dimensional posterior over 40 visibilities and
#: 16 closure phases is sampled rather than explored.
NUTS_DRAWS = 200
NUTS_WARMUP = 200
NUTS_CHAINS = 2

#: ``test_nuts.py``'s tolerance, unchanged: the posterior mean within four of its
#: own standard deviations of the injected value. Loose enough that a correct
#: sampler passes essentially always, and far tighter than the error a driver bug
#: produces — a mis-transposed draw array or an unconstrained-space posterior
#: recorded as a constrained one moves a summary by whole standard deviations.
RECOVERY_SIGMAS = 4.0


@pytest.fixture(scope="module")
def nuts_run(problem: Any) -> Any:
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NUTSEngine(problem).run(draws=NUTS_DRAWS, warmup=NUTS_WARMUP, chains=NUTS_CHAINS)


class TestNUTSRecoversTheBinary:
    """The acceptance row: a gradient through a DFT finds the source."""

    @pytest.mark.parametrize("name", ["model.separation", "model.flux_ratio"])
    def test_the_injected_value_is_within_the_posterior(self, nuts_run: Any, name: str) -> None:
        draws = np.asarray(nuts_run["posterior"][name])
        spread = float(draws.std())
        assert spread > 0.0
        assert abs(float(draws.mean()) - source.TRUTH[name]) < RECOVERY_SIGMAS * spread

    #: The width of each free parameter's uniform prior, for the row below.
    PRIOR_WIDTHS: ClassVar[dict[str, float]] = {
        "model.separation": 14.0,
        "model.flux_ratio": 0.7,
    }

    @pytest.mark.parametrize("name", ["model.separation", "model.flux_ratio"])
    def test_the_posterior_is_informative_rather_than_merely_wide(
        self, nuts_run: Any, name: str
    ) -> None:
        """Width, not location: a sampler that ignored the data passes the row above.

        Containment is the previous row's business, at W2's four-sigma
        tolerance; a *uniform* posterior would satisfy it for the wrong reason,
        so this one says the data actually constrained the source. It is
        deliberately not a 95 % containment check: this posterior is sharp
        enough (about 0.003 in the flux ratio) that one noise realisation in
        twenty legitimately leaves the truth outside its own central 95 %, and
        a row that failed then would be a row about the seed.
        """
        draws = np.asarray(nuts_run["posterior"][name]).reshape(-1)
        low, high = np.percentile(draws, [2.5, 97.5])
        assert (high - low) < 0.1 * self.PRIOR_WIDTHS[name]

    def test_the_run_has_the_shape_it_was_asked_for(self, nuts_run: Any) -> None:
        assert nuts_run["posterior"]["model.separation"].shape == (NUTS_CHAINS, NUTS_DRAWS)
        assert set(nuts_run["posterior"].dataset.data_vars) == {
            "model.separation",
            "model.flux_ratio",
        }

    def test_the_per_dataset_terms_are_both_there(self, nuts_run: Any) -> None:
        """Two observables, two log-likelihood groups, one fit."""
        assert set(nuts_run["log_likelihood"].dataset.data_vars) == {"vis", "t3"}

    def test_the_run_is_recorded_as_realised_on_this_backend(self, nuts_run: Any, kit: Kit) -> None:
        assert nuts_run.attrs["ampere_engine"] == "nuts"
        assert nuts_run.attrs["ampere_backend"] == kit.name
        assert nuts_run.attrs["ampere_realised"] == 1

    def test_the_sampler_found_a_healthy_geometry(self, nuts_run: Any) -> None:
        """Divergences are the diagnostic no gradient-free engine can offer, and
        a chain full of them would recover the truth for the wrong reason."""
        assert int(nuts_run.attrs["ampere_nuts_divergences"]) < 0.05 * (NUTS_DRAWS * NUTS_CHAINS)
        assert 0.0 < float(nuts_run.attrs["ampere_nuts_mean_accept_prob"]) <= 1.0

    def test_the_posterior_is_in_the_constrained_space(self, nuts_run: Any) -> None:
        separation = np.asarray(nuts_run["posterior"]["model.separation"])
        ratio = np.asarray(nuts_run["posterior"]["model.flux_ratio"])
        assert bool(np.all(separation > 6.0) and np.all(separation < 20.0))
        assert bool(np.all(ratio > 0.1) and np.all(ratio < 0.8))


# ---------------------------------------------------------------------------
# 3. The flexible likelihood on the visibilities (W4.2's row)
# ---------------------------------------------------------------------------


class TestTheCircularComplexGPOnVisibilities:
    """W4.2 asked for this: its closed form, realised and sampled, on a real modality."""

    def test_the_problem_carries_the_kernel_hyperparameters(self, gp_problem: Any) -> None:
        assert gp_problem.parameters.free_names == (
            "model.separation",
            "model.flux_ratio",
            "vis.likelihood.amplitude",
            "vis.likelihood.length_scale",
        )
        assert gp_problem.datasets["vis"].likelihood.noise.CORRELATED is True
        # And the closure phases are *not* under a GP: a GP on a wrapped
        # observable is a latent-variable model, which is Phase 5's (W5.1).
        assert bool(getattr(gp_problem.datasets["t3"].likelihood.noise, "CORRELATED", False)) is (
            False
        )

    def test_the_realised_density_agrees_with_the_contract_path(
        self, gp_problem: Any, kit: Kit
    ) -> None:
        """The GP over a *complex* residual: two real columns sharing one factorisation."""
        lowered = kit.backend.lower_problem(gp_problem)
        for amplitude, length in ((0.02, 1.0e7), (0.05, 3.0e7), (0.005, 5.0e6)):
            theta = gp_problem.unconstrain(
                {
                    **source.TRUTH,
                    "vis.likelihood.amplitude": amplitude,
                    "vis.likelihood.length_scale": length,
                }
            )
            assert float(lowered.log_prob_unconstrained(theta)) == pytest.approx(
                gp_problem.log_prob_unconstrained(theta), rel=1e-10
            )

    def test_the_axes_selector_reaches_the_native_density(self, kit: Kit) -> None:
        """W4.2's two native-path fixes, guarded rather than asserted.

        A native path that took ``axes[0]`` as the coordinates, or that read the
        kernel from the unbound declaration rather than from
        ``noise.kernel_for(observed)``, would build its covariance from the ``u``
        column alone — and would look entirely healthy. Two problems differing
        *only* in the selector must therefore give different densities, and each
        must still agree with its own contract path.
        """
        densities = []
        for axes in (("u",), ("u", "v")):
            kernel = kit.backend.Matern32(
                st.loguniform(1e-3, 0.2), st.loguniform(1e6, 1e8), axes=axes
            )
            built = _gp_problem_with(kit, kernel)
            theta = built.unconstrain(
                {
                    **source.TRUTH,
                    "vis.likelihood.amplitude": 0.03,
                    "vis.likelihood.length_scale": 2.0e7,
                }
            )
            native = float(kit.backend.lower_problem(built).log_prob_unconstrained(theta))
            assert native == pytest.approx(built.log_prob_unconstrained(theta), rel=1e-10)
            densities.append(native)
        assert abs(densities[0] - densities[1]) > 1e-3

    def test_nuts_samples_the_gp_visibility_problem(self, gp_problem: Any, kit: Kit) -> None:
        """The row W4.2 asked for, at a guard budget rather than a recovery one.

        Four dimensions, two of them a GP's hyperparameters on data with no
        injected correlation — so the *sky* parameters are identified and the
        hyperparameters are not, which is exactly the geometry a short chain
        should not be asked to characterise. What this row says is that the
        circular complex GP composes, lowers, and is sampled through a gradient
        on this backend; how well it is sampled when there is correlation to
        find is the W4.4 study's question.
        """
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = NUTSEngine(gp_problem).run(draws=30, warmup=30, chains=1)
        assert run["posterior"]["vis.likelihood.amplitude"].shape == (1, 30)
        lp = np.asarray(run["sample_stats"]["lp"])
        assert np.all(np.isfinite(lp))
        assert run.attrs["ampere_backend"] == kit.name
        # The sky parameters are still where they should be, even here.
        draws = np.asarray(run["posterior"]["model.separation"]).reshape(-1)
        assert abs(float(draws.mean()) - source.TRUTH["model.separation"]) < 1.0


def _gp_problem_with(kit: Kit, kernel: Any) -> Any:
    """The two-dataset problem with a caller-supplied kernel on the visibilities."""
    from ampere.core import (
        ComplexGaussianFamily,
        Dataset,
        DatasetCollection,
        FittingProblem,
        Likelihood,
        VonMisesFamily,
    )

    compiled, observed_visibility, observed_phase = source.synthetic(kit.itf)
    vis_instrument = source.chain(kit.itf, observed_visibility, "vis")
    t3_instrument = source.chain(kit.itf, observed_phase, "t3")
    datasets = DatasetCollection(
        {
            "vis": Dataset(
                observed_visibility,
                vis_instrument,
                likelihood=Likelihood(
                    ComplexGaussianFamily(),
                    kit.backend.GaussianProcessNoise(kernel, kit.backend.DenseGP()),
                ),
                label="vis",
            ),
            "t3": Dataset(
                observed_phase,
                t3_instrument,
                likelihood=Likelihood(VonMisesFamily(), kit.backend.IndependentNoise()),
                label="t3",
            ),
        }
    )
    return FittingProblem(source.fitted_model(kit.itf, compiled), datasets, seed=source.SEED)


# ---------------------------------------------------------------------------
# 4. VI
# ---------------------------------------------------------------------------


class TestVariationalInference:
    """The other gradient-based engine, on the same problem."""

    def test_it_finds_the_binary(self, problem: Any) -> None:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = VIEngine(problem).run(steps=400, draws=200)
        for name in ("model.separation", "model.flux_ratio"):
            draws = np.asarray(run["posterior"][name]).reshape(-1)
            spread = float(draws.std())
            assert spread > 0.0
            assert abs(float(draws.mean()) - source.TRUTH[name]) < RECOVERY_SIGMAS * spread
        assert run.attrs["ampere_engine"] == "vi"
        assert run.attrs["ampere_realised"] == 1


# ---------------------------------------------------------------------------
# 5. The batched forward path, against the loop
# ---------------------------------------------------------------------------


#: How closely the batched forward path must reproduce the loop, **scaled by the
#: quantity's own magnitude**. Scaled rather than absolute because a model image
#: is a surface brightness in Jy/sr whose peak here is about 1e16 and whose
#: Gaussian tails reach 1e-225: an absolute tolerance would be meaningless at
#: both ends. torch reproduces the loop bitwise; jax differs in the last bits,
#: which is vectorised arithmetic not being scalar arithmetic and is what
#: ``tolerances.cross_backend`` exists for.
BATCHED_RTOL = 1e-12


def _agrees(got: np.ndarray, expected: np.ndarray, what: str) -> None:
    """*got* equals *expected* to :data:`BATCHED_RTOL` of its own scale."""
    assert got.shape == expected.shape, what
    scale = float(np.max(np.abs(expected)))
    worst = float(np.max(np.abs(got - expected)))
    assert worst <= BATCHED_RTOL * max(scale, 1e-300), f"{what}: {worst:g} vs scale {scale:g}"


class TestTheBatchedForwardPath:
    """``inference.md`` §13: the loop is the semantics, and the fast path is equal to it."""

    COUNT = 4

    @staticmethod
    def loop(kit: Kit, count: int, *, observe: bool) -> SimulationBatch:
        """``n`` calls of ``simulate`` on the batch sub-stream's spawned children.

        Built on a *second* problem with the same seed, because reading
        ``rng("simulate")`` on the problem under test would advance the very
        stream the batch is about to spawn from.
        """
        other = source.two_dataset_problem(kit.backend, kit.itf)
        children = other.rng("simulate").spawn(count)
        return SimulationBatch(
            tuple(other.simulate(observe=observe, rng=child) for child in children)
        )

    def test_the_loop_is_reproduced_bitwise_with_the_fast_path_off(self, kit: Kit) -> None:
        """``native=False`` is the unamended statement, on the two-dataset problem."""
        built = source.two_dataset_problem(kit.backend, kit.itf)
        batch = built.simulate_many(self.COUNT, observe=True, native=False)
        reference = self.loop(kit, self.COUNT, observe=True)
        assert np.array_equal(np.asarray(batch.theta), np.asarray(reference.theta))
        for mine, theirs in zip(batch, reference, strict=True):
            for label in mine.predicted:
                assert np.array_equal(
                    np.asarray(mine.predicted[label].values),
                    np.asarray(theirs.predicted[label].values),
                )
            for label in mine.observations or {}:
                assert np.array_equal(
                    np.asarray(mine.observations[label].values),
                    np.asarray(theirs.observations[label].values),
                )
        assert batch.provenance["simulate_batched"] is False

    def test_the_native_batched_forward_path_equals_the_loop(self, kit: Kit) -> None:
        """The backend's own ``simulate_batched``, over a given θ table.

        This is the equality the acceptance criterion is about, measured where
        the backend hands its answer over: both datasets' instrument-transformed
        predictions and the model's whole ``sky`` channel, for every draw.
        """
        built = source.two_dataset_problem(kit.backend, kit.itf)
        reference = self.loop(kit, self.COUNT, observe=False)
        lowered = kit.backend.lower_problem(built)
        stacked = lowered.simulate_batched(np.asarray(reference.theta, dtype=float))
        assert len(stacked) == self.COUNT
        for index, simulation in enumerate(reference):
            for label, container in simulation.predicted.items():
                _agrees(
                    np.asarray(stacked.predicted[label])[index],
                    np.asarray(container.values),
                    f"dataset {label!r} draw {index}",
                )
            for model, result in simulation.results.items():
                for channel in result:
                    _agrees(
                        np.asarray(stacked.channels[model][channel])[index],
                        np.asarray(result[channel].values).ravel(),
                        f"model {model!r} channel {channel!r} draw {index}",
                    )

    def test_a_multi_axis_channel_is_flat_per_draw(self, kit: Kit) -> None:
        """``BatchedPrediction.channels`` is declared ``{model: {channel: (batch, n)}}``.

        Flat per draw, in the container's own C-order. Every channel before
        Phase 4 was one-dimensional already, so this was a no-op nobody had to
        write; an ``Image`` channel is ``(nx, ny)`` and the backend flattens it.
        """
        built = source.two_dataset_problem(kit.backend, kit.itf)
        lowered = kit.backend.lower_problem(built)
        theta = np.asarray(built.unconstrain(source.TRUTH), dtype=float).reshape(1, -1)
        stacked = lowered.simulate_batched(theta)
        image = built.models["model"].templates["sky"]
        assert np.shape(stacked.channels["model"]["sky"]) == (1, np.size(image.values))

    def test_simulate_many_native_is_blocked_by_a_core_shape_gap(self, kit: Kit) -> None:
        """**A carried defect, pinned here rather than only in a report (W4.3).**

        ``ampere.core.dataset``'s native batched path is internally inconsistent
        about a **multi-axis** channel, and this modality is the first to have
        one. ``_BatchedSampler._check_agreement`` compares the native channel
        stack against ``result[channel].values.ravel()`` — the flat form
        ``BatchedPrediction`` declares — while ``_draw`` hands the same array
        straight to ``template[channel].with_values(...)``, which needs the
        container's own ``(nx, ny)``. A one-dimensional channel satisfies both;
        an ``Image`` cannot satisfy either choice at both sites.

        The fix is one line at ``_draw``::

            channel: template[channel].with_values(
                np.asarray(native.channels[model][channel][position]).reshape(
                    np.shape(template[channel].values)
                )
            )

        verified on this branch: with it applied, ``simulate_many(native=True)``
        runs, records ``provenance["simulate_batched"] is True`` and reproduces
        the loop's predictions exactly. ``ampere/core/dataset.py`` is outside
        this item's file ownership, so the defect is carried rather than fixed,
        and this row holds the present behaviour so that the fix is visible as a
        change. **Replace it with the equality assertion when the fix lands.**
        """
        built = source.two_dataset_problem(kit.backend, kit.itf)
        with pytest.raises((SchemaError, LoweringError), match=r"shape|disagree"):
            built.simulate_many(self.COUNT, observe=True, native=True)


# ---------------------------------------------------------------------------
# 6. SBI
# ---------------------------------------------------------------------------

needs_sbi = pytest.mark.skipif(
    importlib.util.find_spec("sbi") is None, reason="SBIEngine needs ampere[sbi]"
)

#: The smoke budget. Large enough that the density estimator is *trained* rather
#: than initialised (it converges in about 200 epochs and 23 s here) and small
#: enough for a per-PR gate. The claim this fit is held to is **calibration**,
#: not sharpness: at this budget the posterior is wide, and a wide honest
#: posterior is exactly what a correct encoding of a complex container produces
#: while a scrambled one produces a narrow confident wrong answer.
SBI_BUDGET = 1200
SBI_EPOCHS = 300
SBI_DRAWS = 200

#: ``sbi``'s own floor for ``check_sbc`` and ``run_tarp``: both warn below a
#: hundred. Re-conditioning an amortised posterior is free, so the whole check
#: costs one simulation batch beside the fit's own.
CALIBRATION_COUNT = 100
CALIBRATION_DRAWS = 100


class _Narrowed:
    """A posterior wrapper whose draws are squeezed towards their own mean.

    ``tests/inference/test_sbi.py``'s device, and here for its reason: a check
    that only ever passes is not a check. Temperature-scaling a *correct*
    posterior is miscalibration of a known size and a known kind —
    over-confident, unbiased — which is exactly what TARP's negative
    area-to-curve and SBC's U-shaped rank histogram are defined to detect, where
    an under-trained network would be miscalibrated by an amount nobody controls.
    """

    def __init__(self, inner: Any, temperature: float = 0.3) -> None:
        self.inner = inner
        self.temperature = temperature

    def _squeeze(self, drawn: Any) -> Any:
        centre = drawn.mean(dim=0, keepdim=True)
        return centre + (drawn - centre) * self.temperature

    def sample_batched(self, sample_shape: Any, x: Any = None, **kwargs: Any) -> Any:
        return self._squeeze(self.inner.sample_batched(sample_shape, x=x, **kwargs))

    def sample(self, sample_shape: Any, x: Any = None, **kwargs: Any) -> Any:
        return self._squeeze(self.inner.sample(sample_shape, x=x, **kwargs))


@pytest.fixture(scope="module")
def npe_engine(problem: Any) -> Any:
    from ampere.inference import SBIEngine

    engine = SBIEngine(problem, method="npe", budget=SBI_BUDGET, layout="set", embedding="set")
    engine.run(draws=SBI_DRAWS, training={"max_num_epochs": SBI_EPOCHS})
    return engine


@pytest.fixture(scope="module")
def calibration(npe_engine: Any) -> Any:
    return npe_engine.calibrate(count=CALIBRATION_COUNT, posterior_draws=CALIBRATION_DRAWS)


@needs_sbi
class TestSBIOnVisibilitiesAndClosurePhases:
    """The encoding's first complex customer and its first five-axis one, fitted."""

    def test_the_summary_layout_refuses_a_complex_dataset_by_name(self, problem: Any) -> None:
        """W3.2's ``flat`` layout has no encoding for a complex value, and says so.

        The refusal is the right answer rather than a gap: flattening a complex
        array into a real feature vector is a *choice* (real and imaginary parts
        as two columns) that the ``set`` layout makes and this one does not. What
        matters is that the message names the remedy, because ``flat`` is the
        default and a user fitting visibilities meets this first.
        """
        from ampere.inference import SBIEngine

        engine = SBIEngine(problem, method="npe", budget=8)
        with pytest.raises((EngineError, EncodingError), match="layout='set'"):
            engine.run(draws=4, training={"max_num_epochs": 2})

    def test_the_fit_runs_under_the_set_layout(self, npe_engine: Any) -> None:
        assert npe_engine.encoding.kind == "set"
        assert npe_engine.encoding.labels == ("vis", "t3")
        assert npe_engine.encoding.complex_columns is True

    def test_the_trained_posterior_passes_the_coverage_check(self, calibration: Any) -> None:
        """The acceptance row: **calibrated** at the smoke budget.

        TARP is the statistic quoted, because it is the *joint* coverage
        diagnostic and the marginal rank histograms cannot stand in for it: a
        posterior correct in each margin and wrong in their correlation passes
        SBC and fails TARP. Measured here (torch, seed fixed): expected-coverage
        area-to-curve +0.034 against a nominal zero, its KS p-value 0.997, the
        coverage curve 0.63 and 0.69 at the nominal 0.68, and the C2ST between
        the SBC ranks and a uniform baseline 0.54 and 0.52 where 0.5 is
        "indistinguishable".

        The marginal SBC KS p-values are **not** given a 0.05 threshold here, and
        that is deliberate rather than a weakening: at a hundred simulations and
        a hundred draws they are 0.06 and 0.68 at this budget and 0.02 at a third
        of it, so a 0.05 gate would be a row about the training budget and the
        seed. The floor below catches a *catastrophically* non-uniform rank
        histogram, and the C2ST — which is stable across every budget measured —
        carries the uniformity claim.
        """
        assert calibration.sizes["simulation"] == CALIBRATION_COUNT
        ranks = np.asarray(calibration["ranks"].values)
        assert ranks.min() >= 0 and ranks.max() <= CALIBRATION_DRAWS
        assert float(np.max(calibration["c2st_ranks"].values)) < 0.7
        assert float(np.min(calibration["ks_pvalue"].values)) > 0.005
        assert abs(float(calibration.attrs["ampere_calibration_tarp_atc"])) < 0.08
        assert float(calibration.attrs["ampere_calibration_tarp_ks_pvalue"]) > 0.05
        levels = np.asarray(calibration.coords["level"].values)
        curve = np.asarray(calibration["coverage"].values)
        for index in range(curve.shape[1]):
            assert abs(float(np.interp(0.68, levels, curve[:, index])) - 0.68) < 0.12

    def test_the_calibration_records_the_layout_it_was_packed_under(
        self, calibration: Any, npe_engine: Any
    ) -> None:
        """The one thing that could go silently wrong on a complex problem.

        A calibration batch packed under a layout rebuilt from the *simulated*
        containers would standardise its columns differently from the ones the
        network trained on — and for a complex dataset the standardisation is a
        modulus median, so the difference would be real — and would then report a
        different network as calibrated.
        """
        assert calibration.attrs["ampere_calibration_route"] == "sbi"
        assert calibration.attrs["ampere_calibration_encoding_hash"] == npe_engine.encoding.hash
        assert calibration.attrs["ampere_calibration_parameterisation"] == "unconstrained"

    def test_a_temperature_scaled_posterior_fails_the_same_check(
        self, npe_engine: Any, calibration: Any
    ) -> None:
        """The arm that proves the check can fail — same engine, same batch size."""
        narrowed = npe_engine.calibrate(
            count=CALIBRATION_COUNT,
            posterior_draws=CALIBRATION_DRAWS,
            posterior=_Narrowed(npe_engine.posterior),
        )
        assert float(np.min(narrowed["ks_pvalue"].values)) < 0.01
        # Under-dispersion is a *negative* area-to-curve in TARP's convention.
        assert float(narrowed.attrs["ampere_calibration_tarp_atc"]) < -0.05


# ---------------------------------------------------------------------------
# 7. The artefact cache
# ---------------------------------------------------------------------------


class TestTheArtefactCacheOnTwoObservedContainers:
    """W3.5's key, on a problem whose data live in two containers of two kinds."""

    @pytest.fixture(scope="class")
    @staticmethod
    def pieces() -> tuple[Any, Any]:
        """The reference path's own pieces: the key is ``ampere.results``' business."""
        import ampere.core as core
        from ampere.backends.reference import interferometry as reference

        return core, reference

    @staticmethod
    def _key(built: Any) -> Any:
        from ampere.results.artefacts import artefact_key

        layout = EncodingLayout.from_datasets(built.datasets)
        return artefact_key(
            built, layout=layout.hash, method="npe", architecture="maf", budget=64, rounds=1
        )

    def test_two_identical_problems_share_one_digest(self, pieces: Any) -> None:
        first = source.two_dataset_problem(*pieces)
        second = source.two_dataset_problem(*pieces)
        assert self._key(first).digest() == self._key(second).digest()

    def test_a_different_noise_realisation_moves_the_digest(self, pieces: Any) -> None:
        first = source.two_dataset_problem(*pieces)
        other = source.two_dataset_problem(*pieces, seed=source.SEED + 1)
        assert self._key(first).digest() != self._key(other).digest()

    @pytest.mark.parametrize("moved", ["vis", "t3"])
    def test_either_observable_moving_alone_moves_the_digest(self, pieces: Any, moved: str) -> None:
        """The row a two-dataset problem needs: the key covers **both** containers.

        ``artefact_key``'s ``data_hash`` is a hash over ``{label: hash_container(...)}``,
        so a cache that keyed on the visibilities alone would serve a network
        trained against different closure phases — which is the one kind of cache
        hit that is worse than a miss.
        """
        base = source.two_dataset_problem(*pieces)
        perturbed = source.two_dataset_problem(*pieces)
        container = perturbed.datasets[moved].observed
        nudged = container.with_values(np.asarray(container.values) * 1.01)
        from ampere.core import Dataset, DatasetCollection, FittingProblem

        datasets = DatasetCollection(
            {
                label: (
                    Dataset(
                        nudged,
                        dataset.instrument,
                        likelihood=dataset.likelihood,
                        label=label,
                    )
                    if label == moved
                    else dataset
                )
                for label, dataset in perturbed.datasets.items()
            }
        )
        rebuilt = FittingProblem(
            perturbed.models["model"], datasets, seed=base.seed, lenient_compile=True
        )
        assert self._key(base).digest() != self._key(rebuilt).digest()


# ---------------------------------------------------------------------------
# 8. The encoding itself
# ---------------------------------------------------------------------------


class TestTheEncodingOnAComplexAndAFiveAxisContainer:
    """W3.3's packing, on the two kinds that are its first real test of both claims.

    Nothing here was wrong and that is the finding: the layout pads a dataset's
    coordinate columns by *axis position* in the kind's declared order, splits a
    complex value into two columns, and records ``is_complex`` as a per-set
    feature, so a ``VisibilitySet``'s three axes and a ``ClosurePhases``'s five
    coexist in one tensor with no ambiguity a network cannot resolve. These rows
    pin the layout rather than re-deriving it, because the *hash* of that layout
    is what a trained network refuses an observation on.
    """

    def test_the_widths_are_what_the_two_kinds_imply(self, layout: EncodingLayout) -> None:
        assert layout.kind == "set"
        # Five coordinate columns: the widest kind here is ``ClosurePhases``.
        assert layout.coordinates == 5
        assert layout.group("coordinate").width == 5
        assert layout.group("coordinate_features").width == 5 * 2 * layout.fourier_bands
        # Two value columns, because *some* dataset is complex...
        assert layout.complex_columns is True
        assert layout.group("value").width == 2
        assert layout.group("value_asinh").width == 2
        # ...but one log-sigma column, because sigma is the standard deviation of
        # each component of a circular complex Gaussian, not one per component.
        assert layout.group("log_sigma").width == 1

    def test_each_dataset_records_its_own_axes_and_complexness(
        self, layout: EncodingLayout
    ) -> None:
        records = {record.label: record for record in layout.datasets}
        assert records["vis"].axes == ("u", "v", "spectral_axis")
        assert records["t3"].axes == ("u1", "v1", "u2", "v2", "spectral_axis")
        assert records["vis"].is_complex is True
        assert records["t3"].is_complex is False
        assert records["vis"].has_sigma and records["t3"].has_sigma

    def test_the_is_complex_feature_distinguishes_the_two_blocks(
        self, reference_problem: Any, layout: EncodingLayout
    ) -> None:
        """The padded imaginary column of a real dataset is zero; the flag says so."""
        view = unpack(encode_observations(reference_problem.datasets, layout=layout).values, layout)
        blocks = {block.label: block for block in view.per_dataset}
        assert float(np.asarray(blocks["vis"].set_features)[0, 0, 2]) == 1.0
        assert float(np.asarray(blocks["t3"].set_features)[0, 0, 2]) == 0.0
        assert np.all(np.asarray(blocks["t3"].value)[..., 1] == 0.0)
        assert np.any(np.asarray(blocks["vis"].value)[..., 1] != 0.0)

    def test_a_constant_axis_is_not_an_absent_one(
        self, reference_problem: Any, layout: EncodingLayout
    ) -> None:
        """The subtlety a chromatic dataset will care about.

        The visibilities are monochromatic, so their ``spectral_axis`` has a zero
        range and standardises to 0.0 — the same value the two *padded* columns
        carry. The Fourier features are what keep them apart: a constant axis
        contributes ``sin(0) = 0, cos(0) = 1``, while an axis that does not exist
        contributes zeros, which is the amendment W3.3 made deliberately.
        """
        view = unpack(encode_observations(reference_problem.datasets, layout=layout).values, layout)
        block = next(one for one in view.per_dataset if one.label == "vis")
        coordinates = np.asarray(block.coordinate)[0, 0]
        assert coordinates[2] == 0.0 and coordinates[3] == 0.0 and coordinates[4] == 0.0
        bands = layout.fourier_bands
        features = np.asarray(block.coordinate_features)[0, 0]
        constant_axis = features[2 * 2 * bands : 3 * 2 * bands]
        absent_axis = features[3 * 2 * bands : 4 * 2 * bands]
        assert np.allclose(constant_axis[1::2], 1.0)
        assert np.allclose(absent_axis, 0.0)

    def test_the_round_trip_is_exact_for_both_kinds(
        self, reference_problem: Any, layout: EncodingLayout
    ) -> None:
        back = decode(encode_observations(reference_problem.datasets, layout=layout))
        for label in ("vis", "t3"):
            observed = reference_problem.datasets[label].observed
            assert back[label].values.dtype == np.asarray(observed.values).dtype
            assert np.allclose(
                back[label].values, np.asarray(observed.values), rtol=0.0, atol=1e-14
            )
            for axis in observed.axes:
                recovered = back[label].coordinates[axis.name]
                scale = max(float(np.max(np.abs(axis.values))), 1.0)
                assert np.allclose(recovered, axis.values, rtol=0.0, atol=1e-14 * scale)

    def test_a_simulated_batch_encodes_under_the_same_layout(
        self, reference_problem: Any, layout: EncodingLayout
    ) -> None:
        """Training and inference pack identically, which is the whole point."""
        batch = reference_problem.simulate_many(3, observe=True, native=False)
        assert batch.observations is not None
        encoded = encode(
            {label: batch.observations[label] for label in ("vis", "t3")},
            layout=layout,
            batched=True,
        )
        assert encoded.values.shape == (3, layout.row_cap, layout.columns_total)
        assert bool(np.all(np.isfinite(encoded.values)))
        layout.check_against(reference_problem.datasets)

    def test_the_rows_are_the_two_observables_in_order(self, layout: EncodingLayout) -> None:
        assert layout.labels == ("vis", "t3")
        first, second = layout.bounds
        assert first == (0, 24)
        assert second == (24, 24 + 16)
        assert layout.row_cap == 40


def test_the_fitted_problem_is_the_one_the_fixtures_describe() -> None:
    """The shared fixtures' own claim, so a drift is a failure rather than a surprise.

    ``tests/backends/interferometry_fixtures.py`` is read by two suites and by
    this one; the constants it exports are the geometry every row here depends
    on, and the hour-angle coverage in particular is the difference between "NUTS
    recovers the binary" and "NUTS reproduces the prior".
    """
    assert source.STATIONS.shape == (4, 2)
    assert len(source.HOUR_ANGLES) == 4
    assert source.visibility_coverage()[0].size == 24
    assert source.triangle_coverage()[0].size == 16
    assert math.isclose(source.OVERSAMPLING, 4.0)
    assert set(source.TRUTH) == {"model.separation", "model.flux_ratio"}
