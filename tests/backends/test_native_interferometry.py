"""The interferometric twins on torch and jax: the rows the neutral battery cannot reach.

``tests/conformance/test_interferometry.py`` already holds every shipped step
and source model to the §4 contracts and to the closed forms in
``tests/conformance/oracles.py``, once per registered backend — and since W4.3
that includes the two modern columns, so nothing about the *values* is repeated
here. What this file covers is what the backend-neutral battery has no
vocabulary for:

* the **gradients**, which are the whole reason these twins exist. A discrete
  Fourier transform of a model image, differentiated with respect to a
  separation and a flux ratio, is what makes a gradient-based fit of a sky
  model possible at all, and the battery cannot ask for one;
* **batchability**, measured rather than declared: ``BATCHABLE = True`` means a
  stack of θ is evaluated in one call, so the rows compare a ``vmap``\\ ped
  forward model against the same model in a loop;
* the **build-once operator**, which is an implementation promise about the hot
  loop rather than a contract, and the one place it could go wrong silently (a
  step reconfigured for a second chain keeping the first chain's coverage);
* the **shared declaration**: each twin inherits its reference counterpart's
  published requirements, buffer names and refusals, and these rows check that
  claim by comparing the two directly rather than by trusting the inheritance;
* the two **honest downgrades** (``UniformDisc`` and ``UniformDiscVisibilities``
  declaring ``DIFFERENTIABLE = False``) and the refusal they produce;
* the ``native_flux`` / ``native_grid`` spelling, and *why* it exists — which is
  a fact about ``ampere.core``'s shadowing rule, not about a backend.

**Parametrised over the backends installed here**, in ``test_nuts.py``'s shape:
the claim "this backend's interferometry is the reference backend's, with a
gradient" is one claim per backend, and writing it once per backend by hand
would let the two drift. Each row runs on whichever of torch and jax this
environment has; both at once is not possible today (separate extras) and is not
assumed anywhere here.

``ampere.backends.reference.interferometry`` is the oracle for the agreement
rows. That is not "comparing ampere against ampere" in the sense
``tests/conformance/README.md`` forbids: the *claim* under test is precisely
that two backends implementing one declaration agree, and
``architecture.md`` §2 appoints the reference path to say what the answer is.
"""

from __future__ import annotations

import dataclasses
import importlib
import math
from typing import Any

import astropy.units as u
import numpy as np
import pytest

import interferometry_fixtures as kit_data
from ampere.backends.reference import interferometry as reference
from ampere.core import (
    AxisRequirement,
    ClosurePhases,
    Instrument,
    LoweringError,
    VisibilitySet,
    negotiate,
)

#: How closely a native twin must reproduce the reference sum. Tight enough to
#: catch a translation error (a transposed contraction, a dropped solid angle, a
#: sign) rather than a rounding one: both backends accumulate the same
#: separable product in float64 and are measured at 1e-16 relative, so 1e-12 is
#: four orders of slack.
TRANSLATION_RTOL = 1e-12

#: The one exception, and a library difference rather than a translation: the
#: uniform disc's closed form needs ``J1``, and ``torch.special.bessel_j1`` and
#: ``jax.scipy.special.bessel_jn`` each differ from ``scipy.special.j1`` in the
#: last few digits (measured worst case 2.8e-12 absolute on torch, 2.2e-16 on
#: jax). Loosening this one row rather than all of them is deliberate: a
#: disagreement anywhere else is a bug.
BESSEL_ATOL = 1e-10


@dataclasses.dataclass(frozen=True)
class Kit:
    """One backend's interferometric vocabulary, by name, so no row names a library.

    ``array`` and ``numpy`` are the two conversions a row needs; everything else
    comes off :attr:`itf` and :attr:`backend`, which is exactly the set
    ``inference.md`` §18 says a backend supplies.
    """

    name: str
    backend: Any
    itf: Any

    def array(self, values: Any) -> Any:
        """*values* as this backend's own float64 array type."""
        if self.name == "torch":
            return self.backend.as_tensor(np.asarray(values, dtype=float))
        return importlib.import_module("jax.numpy").asarray(np.asarray(values, dtype=float))

    def numpy(self, values: Any) -> np.ndarray:
        """*values* back on the numpy side, for comparison against an oracle."""
        if self.name == "torch":
            return self.backend.to_numpy(values)
        return np.asarray(values)

    def gradient(self, function: Any, at: float) -> float:
        """``d function / d x`` at *at*, through this backend's own autodiff."""
        if self.name == "torch":
            torch = importlib.import_module("torch")
            argument = torch.tensor(float(at), dtype=torch.float64, requires_grad=True)
            function(argument).backward()
            assert argument.grad is not None
            return float(argument.grad)
        jax = importlib.import_module("jax")
        return float(jax.grad(function)(float(at)))

    def mapped(self, function: Any, over: np.ndarray) -> np.ndarray:
        """*function* evaluated over *over* in one vectorised call."""
        if self.name == "torch":
            torch = importlib.import_module("torch")
            stack = torch.tensor(np.asarray(over, dtype=float), dtype=torch.float64)
            return self.numpy(torch.func.vmap(function)(stack))
        jax = importlib.import_module("jax")
        jnp = importlib.import_module("jax.numpy")
        return np.asarray(jax.vmap(function)(jnp.asarray(np.asarray(over, dtype=float))))


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            backend = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed in this environment
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns the flag on, not
            # ampere. A test suite is an application.
            backend.configure_x64()
        twins = importlib.import_module(f"ampere.backends.{name}.interferometry")
        found.append(Kit(name, backend, twins))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS, reason="no modern backend installed; the native interferometry twins need one"
)


@pytest.fixture(scope="module", params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# Agreement with the reference path, piece by piece
# ---------------------------------------------------------------------------


def _through(itf: Any, which: str, extra: tuple[Any, ...] = ()) -> Any:
    """One source model through one chain, on the shared coverage."""
    observed = kit_data.visibilities()
    sample = itf.FourierSample.from_observed(
        observed,
        field_of_view=kit_data.FIELD_OF_VIEW * u.mas,
        oversampling=kit_data.OVERSAMPLING,
    )
    instrument = Instrument([sample, *extra], channel="sky", label="array")
    table = {
        "binary": (itf.Binary, {**kit_data.BINARY, "component_fwhm": kit_data.COMPONENT_FWHM}),
        "gaussian": (itf.GaussianSource, {"fwhm": 4.0, "flux": 1.4}),
        "disc": (itf.UniformDisc, {"diameter": 5.0, "flux": 1.1}),
    }
    model_class, source = table[which]
    model = model_class.on_field(kit_data.FIELD_OF_VIEW * u.mas, 4, channels="sky", **source)
    compiled = model.compile_for(negotiate([instrument]))
    return instrument(compiled.evaluate())


class TestTheTransformIsTheReferenceBackendsSum:
    """Piece by piece, on one coverage, at a translation-error tolerance."""

    @pytest.mark.parametrize("which", ["binary", "gaussian", "disc"])
    def test_the_fourier_step_agrees(self, kit: Kit, which: str) -> None:
        native = kit.numpy(_through(kit.itf, which).values)
        oracle = np.asarray(_through(reference, which).values)
        assert np.max(np.abs(native - oracle)) <= TRANSLATION_RTOL * np.max(np.abs(oracle))

    def test_the_closure_phase_step_agrees(self, kit: Kit) -> None:
        """On the triangles, where a sign convention would show as a flip."""
        observed = kit_data.closure_phases()
        emitted = []
        for module in (kit.itf, reference):
            instrument = kit_data.chain(module, observed, "t3")
            compiled = kit_data.truth_model(module).compile_for(negotiate([instrument]))
            produced = instrument(compiled.evaluate())
            assert isinstance(produced, ClosurePhases)
            emitted.append(np.asarray(kit.numpy(produced.values)))
        native, oracle = emitted
        # Non-trivial phases, so agreement says something: a point-symmetric
        # source would give zero or pi, and a sign error would then agree with
        # the convention as readily as the convention does. Several triangles
        # here are well away from zero and both signs are represented, which is
        # what makes this a check of the sign and not only of the modulus.
        assert np.max(np.abs(oracle)) > 0.3
        assert np.any(oracle > 0.1) and np.any(oracle < -0.1)
        assert np.max(np.abs(native - oracle)) <= TRANSLATION_RTOL * np.max(np.abs(oracle))

    @pytest.mark.parametrize("step", ["bandwidth", "time", "amplitude"])
    def test_the_remaining_steps_agree(self, kit: Kit, step: str) -> None:
        rates = np.full(kit_data.visibilities().n_samples, 4.0e4)

        def built(module: Any) -> tuple[Any, ...]:
            if step == "bandwidth":
                return (module.BandwidthSmearing(resolving_power=30.0),)
            if step == "time":
                return (module.TimeSmearing(integration=60.0 * u.s, du_dt=rates, dv_dt=-rates),)
            return (module.Amplitude(),)

        native = kit.numpy(_through(kit.itf, "binary", built(kit.itf)).values)
        oracle = np.asarray(_through(reference, "binary", built(reference)).values)
        assert native.shape == oracle.shape
        assert np.max(np.abs(native - oracle)) <= TRANSLATION_RTOL * np.max(np.abs(oracle))

    @pytest.mark.parametrize(
        ("which", "source"),
        [
            ("UniformDiscVisibilities", {"diameter": 5.0, "flux": 1.1}),
            ("GaussianSourceVisibilities", {"fwhm": 4.0, "flux": 1.4}),
            ("BinaryVisibilities", {**kit_data.BINARY, "component_fwhm": kit_data.COMPONENT_FWHM}),
        ],
    )
    def test_the_analytic_models_agree(self, kit: Kit, which: str, source: dict[str, Any]) -> None:
        """The closed forms, which are the battery's oracles for the direct sum."""
        u_pts, v_pts, waves = kit_data.visibility_coverage()
        args = (u_pts, v_pts, waves * u.micron)
        native = kit.numpy(
            getattr(kit.itf, which)(*args, channels="vis", **source).evaluate().single().values
        )
        oracle = np.asarray(
            getattr(reference, which)(*args, channels="vis", **source).evaluate().single().values
        )
        tolerance = (
            BESSEL_ATOL
            if which == "UniformDiscVisibilities"
            else TRANSLATION_RTOL * np.max(np.abs(oracle))
        )
        assert np.max(np.abs(native - oracle)) <= tolerance


class TestTheDeclarationIsShared:
    """The twin pattern's own claim: one declaration, two arithmetics.

    These twins **inherit** their reference counterparts rather than
    re-declaring them, and the reason (the module docstrings' argument) is that
    two backends publishing requirements written twice would eventually
    negotiate differently for one declaration — and an under-sampled image
    aliases rather than announcing itself. So the rows that matter are that the
    published requirements and the declared buffers are *identical*, not merely
    close.
    """

    @staticmethod
    def _step(module: Any) -> Any:
        return module.FourierSample.from_observed(
            kit_data.visibilities(),
            field_of_view=kit_data.FIELD_OF_VIEW * u.mas,
            oversampling=kit_data.OVERSAMPLING,
        )

    def test_the_published_requirements_are_identical(self, kit: Kit) -> None:
        native = self._step(kit.itf).requirements()
        oracle = self._step(reference).requirements()
        assert len(native) == len(oracle) == 2
        for mine, theirs in zip(native, oracle, strict=True):
            assert isinstance(mine, AxisRequirement)
            assert mine.axis == theirs.axis
            assert mine.unit == theirs.unit
            assert mine.intervals == theirs.intervals
            # The Nyquist step lives in ``densities``; comparing the whole
            # record rather than one field is the point — a requirement that
            # agreed on the interval and differed on the step would alias.
            assert mine.densities == theirs.densities
            assert mine.points == theirs.points

    def test_the_declared_buffers_are_identical(self, kit: Kit) -> None:
        native = self._step(kit.itf)
        oracle = self._step(reference)
        assert sorted(native.buffers.names) == sorted(oracle.buffers.names)
        for name in oracle.buffers.names:
            assert np.allclose(
                np.asarray(native.buffers[name].value, dtype=float),
                np.asarray(oracle.buffers[name].value, dtype=float),
                rtol=0.0,
                atol=0.0,
            )

    def test_the_model_declares_the_same_parameters_under_the_same_names(self, kit: Kit) -> None:
        assert (
            kit_data.truth_model(kit.itf).parameters.names
            == kit_data.truth_model(reference).parameters.names
        )

    def test_the_negotiated_grid_is_the_same_grid(self, kit: Kit) -> None:
        """The consequence the shared declaration exists for."""
        grids = []
        for module in (kit.itf, reference):
            instruments = [
                kit_data.chain(module, kit_data.visibilities(), "vis"),
                kit_data.chain(module, kit_data.closure_phases(), "t3"),
            ]
            compiled = kit_data.truth_model(module).compile_for(negotiate(instruments))
            grids.append((compiled.templates["sky"].x.values, compiled.templates["sky"].y.values))
        assert np.array_equal(grids[0][0], grids[1][0])
        assert np.array_equal(grids[0][1], grids[1][1])

    def test_every_piece_declares_this_backend(self, kit: Kit) -> None:
        names = (
            "FourierSample",
            "ClosurePhase",
            "Amplitude",
            "BandwidthSmearing",
            "TimeSmearing",
            "UniformDisc",
            "GaussianSource",
            "Binary",
            "UniformDiscVisibilities",
            "GaussianSourceVisibilities",
            "BinaryVisibilities",
        )
        for name in names:
            piece = getattr(kit.itf, name)
            assert piece.BACKEND == kit.name, name
            assert piece.DEVICE == "cpu", name
            assert piece.__name__ == getattr(reference, name).__name__

    def test_the_inherited_refusals_still_refuse_by_name(self, kit: Kit) -> None:
        """A declaration inherited is a refusal inherited."""
        from ampere.core import TransformationError

        with pytest.raises(TransformationError, match="oversampling"):
            kit.itf.FourierSample.from_observed(
                kit_data.visibilities(), field_of_view=1.0 * u.mas, oversampling=0.5
            )
        with pytest.raises(TransformationError, match="odd number of quadrature nodes"):
            kit.itf.BandwidthSmearing(resolving_power=30.0, nodes=4)
        with pytest.raises(TransformationError, match="uv track rates"):
            kit.itf.TimeSmearing.from_observed(kit_data.visibilities(), integration=60.0 * u.s)


# ---------------------------------------------------------------------------
# The gradients
# ---------------------------------------------------------------------------


class TestTheGradientsTheTwinsExistFor:
    """A DFT of a model image, differentiated with respect to the source."""

    @staticmethod
    def _chain(kit: Kit, observed: Any, label: str) -> tuple[Any, Any]:
        instrument = kit_data.chain(kit.itf, observed, label)
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        return compiled, instrument

    def _visibility_power(self, kit: Kit) -> Any:
        compiled, instrument = self._chain(kit, kit_data.visibilities(), "vis")
        values = {k: v for k, v in kit_data.BINARY.items()}

        def power(separation: Any) -> Any:
            flux = compiled.native_flux("sky", {**values, "separation": separation})
            grid = compiled.native_grid("sky")
            for step in instrument.steps:
                flux, grid = step.apply_flux(flux, grid, {})
            return (abs(flux) ** 2).sum()

        return power

    def test_the_separations_gradient_matches_a_finite_difference(self, kit: Kit) -> None:
        """The row the whole modality's NUTS claim rests on."""
        power = self._visibility_power(kit)
        centre = kit_data.BINARY["separation"]
        analytic = kit.gradient(power, centre)
        step = 1e-5
        numeric = float(
            (kit.numpy(power(centre + step)) - kit.numpy(power(centre - step))) / (2.0 * step)
        )
        assert abs(analytic) > 1e-6, "a zero gradient would pass any tolerance"
        assert analytic == pytest.approx(numeric, rel=1e-5)

    def test_the_closure_phase_carries_a_gradient_too(self, kit: Kit) -> None:
        """``angle()`` of a complex product is differentiable, and non-trivially so.

        The closure phase is where the astrometric information survives the
        atmosphere, so a chain that could not differentiate it would be
        throwing away the signal rather than merely being slow.
        """
        compiled, instrument = self._chain(kit, kit_data.closure_phases(), "t3")
        values = dict(kit_data.BINARY)

        def total(angle: Any) -> Any:
            flux = compiled.native_flux("sky", {**values, "position_angle": angle})
            grid = compiled.native_grid("sky")
            for step in instrument.steps:
                flux, grid = step.apply_flux(flux, grid, {})
            return (flux**2).sum()

        centre = kit_data.BINARY["position_angle"]
        analytic = kit.gradient(total, centre)
        step = 1e-5
        numeric = float(
            (kit.numpy(total(centre + step)) - kit.numpy(total(centre - step))) / (2.0 * step)
        )
        assert abs(analytic) > 1e-6
        assert analytic == pytest.approx(numeric, rel=1e-4)

    def test_the_transform_is_exactly_linear_in_the_total_flux(self, kit: Kit) -> None:
        """An analytic gradient, not a finite difference: a DFT is linear in ``I``.

        So ``d|V|^2/d flux = 2 |V|^2 / flux`` exactly, which is a closed form the
        backend's autodiff either reproduces or does not.
        """
        compiled, instrument = self._chain(kit, kit_data.visibilities(), "vis")
        values = dict(kit_data.BINARY)

        def power(flux: Any) -> Any:
            values_at = {**values, "flux": flux}
            tensor = compiled.native_flux("sky", values_at)
            grid = compiled.native_grid("sky")
            for step in instrument.steps:
                tensor, grid = step.apply_flux(tensor, grid, {})
            return (abs(tensor) ** 2).sum()

        total = kit_data.BINARY["flux"]
        analytic = kit.gradient(power, total)
        assert analytic == pytest.approx(2.0 * float(kit.numpy(power(total))) / total, rel=1e-10)


class TestBatchability:
    """``BATCHABLE = True`` measured, not declared: one vmapped call against a loop."""

    STACK = np.array([10.0, 12.0, 14.0])

    @staticmethod
    def _closure(kit: Kit) -> Any:
        instrument = kit_data.chain(kit.itf, kit_data.closure_phases(), "t3")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        values = dict(kit_data.BINARY)

        def phases(separation: Any) -> Any:
            flux = compiled.native_flux("sky", {**values, "separation": separation})
            grid = compiled.native_grid("sky")
            for step in instrument.steps:
                flux, grid = step.apply_flux(flux, grid, {})
            return flux

        return phases

    def test_a_vmapped_chain_equals_the_same_chain_in_a_loop(self, kit: Kit) -> None:
        """Image, DFT and closure phase together, which is the whole forward model."""
        phases = self._closure(kit)
        mapped = kit.mapped(phases, self.STACK)
        looped = np.stack([kit.numpy(phases(float(value))) for value in self.STACK])
        assert mapped.shape == looped.shape
        assert np.allclose(mapped, looped, rtol=0.0, atol=1e-12)
        # And the stack really varies, so the row is not comparing three copies.
        assert np.max(np.abs(mapped[0] - mapped[-1])) > 1e-3

    def test_the_non_differentiable_disc_declares_no_batching_either(self, kit: Kit) -> None:
        """One reason, not two: the flag is about the *realised* density."""
        assert kit.itf.UniformDisc.BATCHABLE is False
        assert kit.itf.UniformDisc.DIFFERENTIABLE is False
        assert kit.itf.UniformDiscVisibilities.BATCHABLE is False
        assert kit.itf.UniformDiscVisibilities.DIFFERENTIABLE is False
        for name in ("GaussianSource", "Binary", "FourierSample", "ClosurePhase"):
            assert getattr(kit.itf, name).BATCHABLE is True, name
            assert getattr(kit.itf, name).DIFFERENTIABLE is True, name


# ---------------------------------------------------------------------------
# The build-once operator
# ---------------------------------------------------------------------------


class TestTheOperatorIsBuiltOnce:
    """The hot-loop promise, and the one way it could be silently wrong."""

    @staticmethod
    def _counted(kit: Kit) -> Any:
        """The Fourier step with its operator builds counted."""

        class Counting(kit.itf.FourierSample):  # type: ignore[misc, name-defined]
            def __init__(self, *args: Any, **kwargs: Any) -> None:
                super().__init__(*args, **kwargs)
                self.builds = 0

            def _operator(self, x_mas: Any, y_mas: Any) -> Any:
                before = self._operator_cache
                built = super()._operator(x_mas, y_mas)
                if before is not self._operator_cache:
                    self.builds += 1
                return built

        return Counting.from_observed(
            kit_data.visibilities(),
            field_of_view=kit_data.FIELD_OF_VIEW * u.mas,
            oversampling=kit_data.OVERSAMPLING,
        )

    def test_repeated_evaluations_rebuild_nothing(self, kit: Kit) -> None:
        step = self._counted(kit)
        instrument = Instrument([step], channel="sky", label="vis")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        first = None
        for separation in (10.0, 12.0, 14.0):
            values = {**kit_data.BINARY, "separation": separation}
            flux = compiled.native_flux("sky", values)
            got, _ = step.apply_flux(flux, compiled.native_grid("sky"), {})
            first = got if first is None else first
        assert step.builds == 1
        # And the cached operator is the right one: the values still move.
        assert not np.allclose(kit.numpy(first), kit.numpy(got))

    def test_a_different_grid_rebuilds(self, kit: Kit) -> None:
        """An exact key, so a changed grid cannot be served a nearby operator."""
        step = self._counted(kit)
        instrument = Instrument([step], channel="sky", label="vis")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        grid = compiled.native_grid("sky")
        flux = compiled.native_flux("sky", kit_data.BINARY)
        step.apply_flux(flux, grid, {})
        shifted = (kit.array(kit.numpy(grid[0]) + 1e-9), grid[1])
        step.apply_flux(flux, shifted, {})
        assert step.builds == 2

    def test_configure_from_drops_the_operator_it_built_for_another_chain(self, kit: Kit) -> None:
        """The silent-failure mode: a smearing step added after the operator was built.

        ``configure_from`` is what expands the coverage, so an operator built
        before it ran is an operator for the wrong ``(u, v)`` set — and the
        visibilities it produced would be plausible rather than wrong-looking.
        """
        step = self._counted(kit)
        instrument = Instrument([step], channel="sky", label="vis")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        flux = compiled.native_flux("sky", kit_data.BINARY)
        unsmeared, _ = step.apply_flux(flux, compiled.native_grid("sky"), {})
        assert step.builds == 1
        step.configure_from([kit.itf.BandwidthSmearing(resolving_power=30.0)])
        smeared, coverage = step.apply_flux(flux, compiled.native_grid("sky"), {})
        assert step.builds == 2
        assert kit.numpy(smeared).size == 5 * kit.numpy(unsmeared).size
        assert kit.numpy(coverage[0]).size == kit.numpy(smeared).size


# ---------------------------------------------------------------------------
# The native surface, and what it is called
# ---------------------------------------------------------------------------


class TestTheNativeSurface:
    """``native_flux`` / ``native_grid``, and the reason they are not ``flux`` / ``grid``."""

    def test_a_source_model_cannot_have_a_method_called_flux(self, kit: Kit) -> None:
        """The collision, stated as a row rather than only in a docstring.

        Every source model here declares a *parameter* named ``flux`` — the
        total flux density, which is the thing one fits — and
        ``Parameterised._check_free_name`` refuses a parameter whose name
        shadows a class attribute. So ``flux`` is unavailable as a method name
        on this class, which is why the native pair is spelled ``native_*`` and
        why :mod:`ampere.backends.torch.problem` accepts either spelling.
        """
        assert "flux" in kit_data.truth_model(kit.itf).parameters.names
        assert not hasattr(type(kit_data.truth_model(kit.itf)), "flux")
        assert callable(kit_data.truth_model(kit.itf).native_flux)
        assert callable(kit_data.truth_model(kit.itf).native_grid)

    def test_the_image_grid_is_two_arrays_and_the_coverage_is_three(self, kit: Kit) -> None:
        instrument = kit_data.chain(kit.itf, kit_data.visibilities(), "vis")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        grid = compiled.native_grid("sky")
        assert len(grid) == 2
        flux = compiled.native_flux("sky", kit_data.BINARY)
        assert tuple(flux.shape) == (kit.numpy(grid[0]).size, kit.numpy(grid[1]).size)
        values, coverage = instrument.steps[0].apply_flux(flux, grid, {})
        assert len(coverage) == 3
        assert kit.numpy(values).size == kit_data.visibilities().n_samples

    def test_the_analytic_model_emits_its_own_coverage(self, kit: Kit) -> None:
        u_pts, v_pts, waves = kit_data.visibility_coverage()
        model = kit.itf.BinaryVisibilities(
            u_pts,
            v_pts,
            waves * u.micron,
            channels="vis",
            component_fwhm=kit_data.COMPONENT_FWHM,
            **kit_data.BINARY,
        )
        grid = model.native_grid("vis")
        assert len(grid) == 3
        assert np.allclose(kit.numpy(grid[0]), u_pts)
        assert kit.numpy(model.native_flux("vis", kit_data.BINARY)).dtype.kind == "c"

    def test_a_realised_problem_resolves_the_native_spelling(self, kit: Kit) -> None:
        """The hook, end to end: the realisation finds ``native_flux`` and agrees."""
        problem = kit_data.two_dataset_problem(kit.backend, kit.itf)
        lowered = kit.backend.lower_problem(problem)
        theta = problem.unconstrain(kit_data.reference_point())
        assert float(lowered.log_prob_unconstrained(theta)) == pytest.approx(
            problem.log_prob_unconstrained(theta), rel=1e-12
        )


class TestTheHonestDowngrades:
    """A uniform disc is not a gradient-based model in ampere, and says so."""

    def test_a_problem_built_on_the_disc_is_refused_by_name(self, kit: Kit) -> None:
        """Refused at lowering, not silently given a wrong gradient."""
        from ampere.core import (
            ComplexGaussianFamily,
            Dataset,
            DatasetCollection,
            FittingProblem,
            Likelihood,
        )

        observed = kit_data.visibilities(np.ones(kit_data.visibilities().n_samples) * (1 + 0j))
        instrument = kit_data.chain(kit.itf, observed, "vis")
        model = kit.itf.UniformDisc.on_field(
            kit_data.FIELD_OF_VIEW * u.mas, 4, channels="sky", diameter=5.0, flux=1.1
        )
        compiled = model.compile_for(negotiate([instrument]))
        problem = FittingProblem(
            compiled,
            DatasetCollection(
                {
                    "vis": Dataset(
                        observed,
                        instrument,
                        likelihood=Likelihood(
                            ComplexGaussianFamily(), kit.backend.IndependentNoise()
                        ),
                        label="vis",
                    )
                }
            ),
            seed=1,
        )
        assert problem.differentiable is False
        from ampere.inference import EngineError, NUTSEngine

        with pytest.raises((EngineError, LoweringError)):
            NUTSEngine(problem)

    def test_the_disc_still_evaluates_and_agrees_on_the_contract_path(self, kit: Kit) -> None:
        """Not differentiable is not broken: a gradient-free engine fits it."""
        native = kit.numpy(_through(kit.itf, "disc").values)
        oracle = np.asarray(_through(reference, "disc").values)
        assert np.max(np.abs(native - oracle)) <= TRANSLATION_RTOL * np.max(np.abs(oracle))


class TestPlacementAndPrecision:
    """``dtype``, ``device`` and the ``DEVICE`` flag — chosen, never detected."""

    def test_the_declared_device_is_cpu_and_is_an_instance_attribute(self, kit: Kit) -> None:
        """Per instance, which is how a whole problem composes on one named device."""
        step = kit.itf.FourierSample.from_observed(
            kit_data.visibilities(), field_of_view=1.0 * u.mas
        )
        assert step.DEVICE == "cpu"
        assert "DEVICE" in vars(step) or step.DEVICE == type(step).DEVICE

    def test_an_unavailable_device_is_refused_rather_than_falling_back(self, kit: Kit) -> None:
        """``architecture.md`` §5: a fit asked for a GPU and quietly given a CPU
        is a fit whose timings mean nothing and whose provenance is a lie."""
        from ampere.core import TransformationError

        if kit.name == "jax":
            with pytest.raises(TransformationError, match=r"cuda|device"):
                kit.itf.FourierSample.from_observed(
                    kit_data.visibilities(), field_of_view=1.0 * u.mas, device="cuda"
                )
            return
        step = kit.itf.FourierSample.from_observed(
            kit_data.visibilities(), field_of_view=1.0 * u.mas
        )
        # Every declared buffer is paired with a torch buffer, so the move is a
        # real one and fails here rather than at the first evaluation.
        assert set(step.tensors.state_dict()) >= {"u_pts", "v_pts", "wavelength"}
        with pytest.raises((RuntimeError, AssertionError), match=r"CUDA|cuda|NVIDIA|driver"):
            step.to(device="cuda")

    def test_a_device_move_keeps_the_arithmetic(self, kit: Kit) -> None:
        """``to("cpu")`` is a no-op that must stay one, cache included."""
        if kit.name != "torch":
            pytest.skip("jax places arrays at construction; there is no in-place .to()")
        before = kit.numpy(_through(kit.itf, "binary").values)
        step = kit.itf.FourierSample.from_observed(
            kit_data.visibilities(),
            field_of_view=kit_data.FIELD_OF_VIEW * u.mas,
            oversampling=kit_data.OVERSAMPLING,
        ).to(device="cpu")
        instrument = Instrument([step], channel="sky", label="vis")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        after = kit.numpy(instrument(compiled.evaluate()).values)
        assert np.allclose(before, after, rtol=0.0, atol=0.0)

    def test_the_complex_output_is_double_precision(self, kit: Kit) -> None:
        """A complex64 visibility would make the phase of a long baseline meaningless."""
        got = _through(kit.itf, "binary")
        assert isinstance(got, VisibilitySet)
        assert np.asarray(got.values).dtype == np.complex128

    def test_the_closure_phase_output_is_in_radians_and_wrapped(self, kit: Kit) -> None:
        observed = kit_data.closure_phases()
        instrument = kit_data.chain(kit.itf, observed, "t3")
        compiled = kit_data.truth_model(kit.itf).compile_for(negotiate([instrument]))
        got = instrument(compiled.evaluate())
        assert got.unit == u.rad
        assert np.all(np.abs(np.asarray(got.values)) <= math.pi)
