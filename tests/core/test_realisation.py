"""The realisation registry (W2.13 prototype): one slot per backend, checked on use.

The claim is structural: a gradient-based engine written against
``ampere.core`` alone can obtain a backend-native, differentiable density
without importing the backend, and the density it obtains is held to the numpy
contract path. These rows use a fake backend so they run in the dependency-free
``dev`` environment; ``tests/inference/test_nuts.py`` exercises the real jax
route.
"""

from __future__ import annotations

from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Capabilities,
    Dataset,
    FittingProblem,
    LoweringError,
    Model,
    Parameter,
    ParameterError,
    Realisation,
    Spectrum,
    realise,
    register_realisation,
    registered_realisations,
)


class Line(Model):
    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_parameter(Parameter("slope", st.norm(2.0, 1.0)))

    def evaluate(self, **values: Any) -> Any:
        ctx = self.context(values)
        return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)


def _problem(**kwargs: Any) -> FittingProblem:
    grid = np.array([1.0, 2.0, 3.0])
    observed = Spectrum(grid * u.um, [2.1, 3.9, 6.2] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy)
    return FittingProblem(Line(grid), [Dataset(observed)], seed=1, **kwargs)


class _Fake:
    """A realisation that computes the density on the contract path — agreeing by construction."""

    def __init__(self, problem: FittingProblem, backend: str = "fake") -> None:
        self._problem = problem
        self._backend = backend

    @property
    def backend(self) -> str:
        return self._backend

    @property
    def free_size(self) -> int:
        return self._problem.free_size

    def log_prob_unconstrained(self, theta: Any) -> Any:
        return self._problem.log_prob_unconstrained(np.asarray(theta))


class TestTheRegistry:
    def test_nothing_is_registered_for_the_reference_backend(self) -> None:
        # The reference backend has no differentiable path, and says so by name.
        with pytest.raises(
            LoweringError, match="no realisation is registered for backend 'reference'"
        ):
            realise(_problem())

    def test_the_refusal_names_the_remedies(self) -> None:
        with pytest.raises(LoweringError) as caught:
            realise(_problem())
        message = str(caught.value)
        assert "gradient-free engine" in message
        assert "ampere.backends.jax" in message

    def test_a_registered_factory_is_dispatched_on_the_problem_backend(self) -> None:
        register_realisation("fake", _Fake)
        problem = _problem(capabilities=Capabilities(backend="fake"))
        realised = realise(problem)
        assert isinstance(realised, Realisation)
        assert realised.backend == "fake"
        assert realised.free_size == problem.free_size
        assert registered_realisations() == {"fake": False}

    def test_a_second_registration_is_refused_without_override(self) -> None:
        register_realisation("fake", _Fake)
        with pytest.raises(LoweringError, match="already registered"):
            register_realisation("fake", _Fake)
        register_realisation("fake", _Fake, override=True)

    def test_builtin_rows_are_distinguished(self) -> None:
        register_realisation("fake", _Fake, builtin=True)
        assert registered_realisations() == {"fake": True}
        with pytest.raises(LoweringError, match="ampere's own"):
            register_realisation("fake", _Fake)

    def test_malformed_registrations_are_refused(self) -> None:
        with pytest.raises(ParameterError, match="backend name"):
            register_realisation("", _Fake)
        with pytest.raises(ParameterError, match="callable"):
            register_realisation("fake", None)  # type: ignore[arg-type]


class TestTheChecksOnUse:
    def test_a_result_naming_another_backend_is_refused(self) -> None:
        register_realisation("fake", lambda problem: _Fake(problem, backend="other"))
        with pytest.raises(LoweringError, match="naming backend 'other'"):
            realise(_problem(capabilities=Capabilities(backend="fake")))

    def test_a_result_with_the_wrong_free_size_is_refused(self) -> None:
        class Wrong(_Fake):
            @property
            def free_size(self) -> int:
                return 99

        register_realisation("fake", Wrong)
        with pytest.raises(LoweringError, match="free_size 99"):
            realise(_problem(capabilities=Capabilities(backend="fake")))

    def test_a_result_that_disagrees_with_the_contract_path_is_refused(self) -> None:
        class Drifted(_Fake):
            def log_prob_unconstrained(self, theta: Any) -> Any:
                return super().log_prob_unconstrained(theta) + 1.0

        register_realisation("fake", Drifted)
        with pytest.raises(LoweringError, match="disagrees with the contract path"):
            realise(_problem(capabilities=Capabilities(backend="fake")))

    def test_a_result_missing_the_surface_is_refused(self) -> None:
        register_realisation("fake", lambda problem: object())
        with pytest.raises(LoweringError, match="does not provide"):
            realise(_problem(capabilities=Capabilities(backend="fake")))

    def test_a_factory_that_cannot_evaluate_is_reported(self) -> None:
        class Broken(_Fake):
            def log_prob_unconstrained(self, theta: Any) -> Any:
                raise RuntimeError("no such thing")

        register_realisation("fake", Broken)
        with pytest.raises(LoweringError, match="could not be evaluated"):
            realise(_problem(capabilities=Capabilities(backend="fake")))
