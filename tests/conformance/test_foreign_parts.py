"""Conformance rows for W3.8: a non-native part inside a native problem.

Ruled by Peter 2026-09-08 on W2.4 slice 3's carried finding — an
``ampere.core`` (numpy) kernel inside a torch or jax GP noise model used to
compose happily while the covariance was built in numpy, converted, and the
amplitude and length scale silently got no gradient at all. The ruling has
three parts, and this module is the per-backend proof of each:

1. **refused by default** — ``declared_capabilities`` treats it as the backend
   disagreement W2.13's fold-in 7 already gives a numpy solver, by name and
   with the remedy;
2. **accepted under an explicit opt-in when no gradient is needed** —
   ``FittingProblem(..., allow_foreign_parts=True)``, recorded in provenance,
   with the problem's ``differentiable`` becoming ``False`` so the capability
   ladder tells the truth, and sampling to the same posterior as the
   all-native twin under a gradient-free engine;
3. **always refused where a gradient is required** — ``realise``,
   ``NUTSEngine`` and ``VIEngine``, by name, whatever the flag says.

Every row is generated per registered backend, and no row names one. The rows
skip on a fixture whose *own* pieces declare the silent default backend, and
that is not a special case dressed up as one: "foreign" is defined relative to
a native backend, and a fixture that is itself the reference path has nothing
for a core kernel to be foreign *to*.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from ampere.core import Kernel, Matern32, SquaredExponential, realise
from ampere.core.exceptions import DatasetError, LoweringError
from ampere.inference import EmceeEngine, NUTSEngine, VIEngine
from ampere.inference.exceptions import EngineError
from ampere.results import provenance_attrs

from .composition import GP_GRID, DatasetSpec, NoiseKind, ProblemSpec, build_problem
from .protocol import (
    ConformanceBackend,
    CovarianceSpec,
    KernelFamily,
    ModelKind,
    ModelSpec,
    Tolerances,
)

#: A GP problem, because a kernel is the piece W3.8 was ruled on and the only
#: composed piece a solver reaches by calling back into it.
GP_SPEC = ProblemSpec(
    model=ModelSpec(kind=ModelKind.POWER_LAW, coordinates=GP_GRID),
    datasets=(
        DatasetSpec(
            label="sed",
            noise=NoiseKind.GP,
            covariance=CovarianceSpec(family=KernelFamily.MATERN32),
        ),
    ),
    seed=20260909,
)

#: ``ampere.core``'s own kernels, keyed the way a fixture's are. These are the
#: foreign pieces: the same declaration, the same arithmetic, and the silent
#: ``BACKEND = "reference"`` every piece inherits unless it says otherwise.
_CORE_KERNELS = {
    KernelFamily.MATERN32: Matern32,
    KernelFamily.SQUARED_EXPONENTIAL: SquaredExponential,
}


def core_kernel(spec: CovarianceSpec) -> Kernel:
    """``build_problem``'s ``kernel_factory``, building ``ampere.core``'s kernel."""
    return _CORE_KERNELS[spec.family](spec.amplitude, spec.length_scale)


@pytest.fixture
def native(backend: ConformanceBackend) -> ConformanceBackend:
    """*backend*, or a skip when a core kernel would not be foreign to it.

    Structural, not a name check: the fixture is asked for the kernel it would
    ordinarily compose and the row runs only if that kernel declares a backend
    of its own. A fixture whose pieces are the reference path's cannot exhibit
    a foreign part at all, and asserting a refusal that cannot happen would be
    asserting nothing.
    """
    kernel = backend.kernel(CovarianceSpec(family=KernelFamily.MATERN32))
    if str(getattr(kernel, "BACKEND", "reference")) == "reference":
        pytest.skip(
            f"{backend.name}'s own kernels declare the default backend, so an "
            f"ampere.core kernel is not foreign to them"
        )
    return backend


def foreign_problem(backend: ConformanceBackend, **kwargs: Any) -> Any:
    """:data:`GP_SPEC` on *backend*, with ``ampere.core``'s kernel substituted."""
    return build_problem(backend, GP_SPEC, kernel_factory=core_kernel, **kwargs)


class TestRefusedByDefault:
    """Part 1 of the ruling: the composition is a backend disagreement."""

    def test_the_default_refuses_the_composition(self, native: ConformanceBackend) -> None:
        with pytest.raises(DatasetError) as excinfo:
            foreign_problem(native)
        assert "different backends" in str(excinfo.value)

    def test_the_refusal_names_the_part_and_the_remedy(self, native: ConformanceBackend) -> None:
        with pytest.raises(DatasetError) as excinfo:
            foreign_problem(native)
        message = str(excinfo.value)
        # By its *qualified* name. The bare one would not discriminate: the
        # backend's own kernel is also called Matern32, deliberately, because
        # Likelihood.to_spec records the declaration and the cross-backend
        # spec-hash row compares it.
        assert "ampere.core.likelihood.Matern32" in message
        assert native.name in message
        # Both halves of the remedy: build it natively, or opt in and give up
        # the gradient.
        assert "allow_foreign_parts=True" in message
        assert "gradient-free" in message

    def test_the_all_native_twin_composes(self, native: ConformanceBackend) -> None:
        problem = build_problem(native, GP_SPEC)
        assert problem.backend == native.name
        assert problem.foreign_parts == ()


class TestTheOptIn:
    """Part 2: accepted, non-differentiable, recorded, and the same posterior."""

    def test_the_opt_in_accepts_it(self, native: ConformanceBackend) -> None:
        problem = foreign_problem(native, allow_foreign_parts=True)
        assert problem.backend == native.name
        assert problem.allow_foreign_parts

    def test_the_problem_is_no_longer_differentiable(self, native: ConformanceBackend) -> None:
        # The capability ladder telling the truth, which is the whole price of
        # the opt-in: whatever every other piece declares, a gradient cannot
        # cross the conversion this kernel forces.
        problem = foreign_problem(native, allow_foreign_parts=True)
        assert not problem.differentiable
        twin = build_problem(native, GP_SPEC)
        assert twin.differentiable == native.capabilities.differentiable

    def test_the_foreign_part_is_named_and_located(self, native: ConformanceBackend) -> None:
        problem = foreign_problem(native, allow_foreign_parts=True)
        assert problem.foreign_part_names == ("sed: ampere.core.likelihood.Matern32",)

    def test_provenance_records_the_flag_and_the_names(self, native: ConformanceBackend) -> None:
        problem = foreign_problem(native, allow_foreign_parts=True)
        attrs = provenance_attrs(problem, engine="emcee")
        assert "ampere.core.likelihood.Matern32" in attrs["ampere_foreign_parts"]
        assert "sed" in attrs["ampere_foreign_parts"]
        # And it is absent — not empty — from every other run, so a reader can
        # tell "no foreign parts" from "a run that never said".
        assert "ampere_foreign_parts" not in provenance_attrs(
            build_problem(native, GP_SPEC), engine="emcee"
        )

    def test_the_density_agrees_with_the_all_native_problem(
        self, native: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        # Before the sampler: the two declarations are one problem, so their
        # log densities are one function. If this fails the posterior row below
        # would fail for a reason that has nothing to do with W3.8.
        opted = foreign_problem(native, allow_foreign_parts=True)
        twin = build_problem(native, GP_SPEC)
        rng = np.random.default_rng(20260909)
        for _ in range(8):
            theta = opted.unconstrain(opted.parameters.sample(rng))
            assert float(opted.log_prob_unconstrained(theta)) == pytest.approx(
                float(twin.log_prob_unconstrained(theta)),
                abs=tolerances.cross_backend,
                rel=tolerances.cross_backend,
            )

    def test_emcee_samples_it_to_the_same_posterior(
        self, native: ConformanceBackend, tolerances: Tolerances
    ) -> None:
        """The gradient-free run the opt-in exists to allow.

        Both problems carry the same ``seed``, so emcee's initial ensemble and
        its whole proposal stream are the same; the draws therefore agree
        element by element rather than only in distribution, which is a far
        sharper assertion than a moment comparison and the one
        ``tolerances.cross_backend`` is for.

        ``use_realisation=False`` on both sides for the same reason: the
        all-native problem would otherwise score through its realisation and
        the opt-in one through the contract path, and this row is about the
        foreign part rather than about the fast path (which falls back, and is
        asserted to below).
        """
        opted = foreign_problem(native, allow_foreign_parts=True)
        twin = build_problem(native, GP_SPEC)
        runs = [
            EmceeEngine(problem, walkers=8, use_realisation=False).run(steps=60, burn_in=20)
            for problem in (opted, twin)
        ]
        for name in opted.parameters.free_names:
            first = np.asarray(runs[0].posterior[name].values)
            second = np.asarray(runs[1].posterior[name].values)
            assert first.shape == second.shape
            assert np.allclose(
                first, second, atol=tolerances.cross_backend, rtol=tolerances.cross_backend
            )


class TestTheGradientRoutesRefuse:
    """Part 3: no flag buys a derivative that does not exist."""

    def test_realise_refuses_by_name(self, native: ConformanceBackend) -> None:
        problem = foreign_problem(native, allow_foreign_parts=True)
        with pytest.raises(LoweringError) as excinfo:
            realise(problem)
        message = str(excinfo.value)
        assert "ampere.core.likelihood.Matern32" in message
        assert "allow_foreign_parts=True" in message

    def test_realise_refuses_a_strict_problem_too(self, native: ConformanceBackend) -> None:
        # ``strict=True`` is the other end of the same axis — it makes the
        # realisation factory refuse at construction rather than return -inf at
        # runtime — and it does not interact with this refusal at all.
        problem = foreign_problem(native, allow_foreign_parts=True, strict=True)
        with pytest.raises(LoweringError, match=r"ampere\.core\.likelihood\.Matern32"):
            realise(problem)

    @pytest.mark.parametrize("engine", [NUTSEngine, VIEngine])
    def test_the_gradient_engines_refuse_by_name(
        self, native: ConformanceBackend, engine: type
    ) -> None:
        problem = foreign_problem(native, allow_foreign_parts=True)
        with pytest.raises(EngineError) as excinfo:
            engine(problem)
        message = str(excinfo.value)
        assert "ampere.core.likelihood.Matern32" in message
        assert "sed" in message
        # The refusal has to say what the flag does *not* do, or a user reads
        # "allow_foreign_parts" and tries it again with the flag already on.
        assert "gradient-free" in message

    def test_the_fast_path_falls_back_to_the_contract_path(
        self, native: ConformanceBackend
    ) -> None:
        """``realise`` refusing is what makes the gradient-free run work.

        ``_EvaluationCache`` treats a ``LoweringError`` from ``realise`` as a
        legitimate "this problem has no usable realisation" and scores on the
        numpy contract path instead — exactly as it already does for a problem
        whose family the backend does not lower. So the refusal above is not a
        wall the gradient-free engines have to route around; it *is* the
        routing.
        """
        engine = EmceeEngine(foreign_problem(native, allow_foreign_parts=True), walkers=8)
        assert engine._cache.realisation is None
        run = engine.run(steps=20, burn_in=5)
        assert int(run.attrs["ampere_realised"]) == 0
