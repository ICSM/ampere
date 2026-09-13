"""W4.7's accept criteria: conformance, NUTS, and the two new refusals, per backend.

Parametrised over whichever of ``torch``/``jax`` this environment has, exactly
as ``tests/inference/test_nuts.py`` is — the claim under test is one claim per
backend, and writing it once per backend by hand would let the two drift.
Skips entirely where neither extra is installed (``pixi run -e dev`` still
collects this file, empty).

1. **Conformance**: each of the six curated leaves, and a handful of compound
   expressions, agree with ``ampere.core.from_astropy``'s black-box adapter at
   ``tolerances.cross_backend`` (``tests/conformance/protocol.py``, 1e-9 —
   inlined below rather than imported, since ``tests/`` is not an importable
   package outside ``tests/conformance`` itself).
2. **NUTS** runs on a translated compound model and recovers an injected
   truth.
3. The **tied-parameter refusal** and the **unsupported-operator refusal**
   (:mod:`ampere.core.astropy_translations`), which
   ``tests/core/test_astropy_backend_hook.py`` does not cover (that file holds
   the leaf-decomposition refusal `translation_refusal` already tests).
"""

from __future__ import annotations

import dataclasses
import importlib
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import (
    BlackBody,
    BrokenPowerLaw1D,
    Const1D,
    Gaussian1D,
    Polynomial1D,
    PowerLaw1D,
)

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Spectrum,
)
from ampere.core import from_astropy as adapt_astropy
from ampere.core.exceptions import CapabilityError
from ampere.inference import NUTSEngine

#: ``tests/conformance/protocol.py``'s ``Tolerances.cross_backend`` — the
#: comparison class this item's Accept criterion names. Not imported: ``tests``
#: is not a package outside ``tests/conformance`` (no ``tests/__init__.py``),
#: so a cross-directory dotted import would only work by accident of how
#: pytest happens to be invoked.
CROSS_BACKEND_TOLERANCE = 1e-9

WAVELENGTH = np.geomspace(1.0, 20.0, 12)


# ---------------------------------------------------------------------------
# Backends this environment can run
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class Kit:
    name: str
    module: Any


def _installed_kits() -> list[Kit]:
    found: list[Kit] = []
    for name in ("jax", "torch"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:
            continue
        if name == "jax":
            module.configure_x64()
        found.append(Kit(name, module))
    return found


KITS = _installed_kits()

pytestmark = pytest.mark.skipif(
    not KITS, reason="needs the torch or jax extra (neither is installed here)"
)


@pytest.fixture(params=[kit.name for kit in KITS])
def kit(request: Any) -> Kit:
    return next(found for found in KITS if found.name == request.param)


# ---------------------------------------------------------------------------
# 1. Conformance: each curated leaf, standalone
# ---------------------------------------------------------------------------


def _bounded(model: Any, **bounds: tuple[float, float]) -> Any:
    for name, bound in bounds.items():
        getattr(model, name).bounds = bound
    return model


def _blackbody_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = _bounded(
        BlackBody(temperature=3000.0 * u.K, scale=1.0),
        temperature=(100.0, 9000.0),
        scale=(0.001, 100.0),
    )
    return (
        model,
        {"temperature": 3200.0, "scale": 2.5},
        {
            "output_unit": u.Jy,
            "equivalencies": [u.dimensionless_angles()],
        },
    )


def _power_law_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.5)
    model.amplitude.bounds = (0.01, 100.0)
    model.x_0.fixed = True
    model.alpha.bounds = (0.01, 5.0)
    return model, {"amplitude": 3.4, "alpha": 1.2}, {"output_unit": u.Jy}


def _broken_power_law_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = _bounded(
        BrokenPowerLaw1D(amplitude=2.0, x_break=5.0, alpha_1=1.0, alpha_2=2.0),
        amplitude=(0.01, 20.0),
        x_break=(1.0, 15.0),
        alpha_1=(0.01, 5.0),
        alpha_2=(0.01, 5.0),
    )
    return (
        model,
        {"amplitude": 1.6, "x_break": 6.0, "alpha_1": 0.8, "alpha_2": 2.3},
        {"output_unit": u.Jy},
    )


def _polynomial_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = _bounded(
        Polynomial1D(degree=2, c0=1.0, c1=0.1, c2=0.01),
        c0=(-10.0, 10.0),
        c1=(-10.0, 10.0),
        c2=(-10.0, 10.0),
    )
    return model, {"c0": 1.3, "c1": 0.2, "c2": 0.02}, {"output_unit": u.Jy}


def _gaussian_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = _bounded(
        Gaussian1D(amplitude=3.0, mean=8.0, stddev=2.0),
        amplitude=(0.01, 20.0),
        mean=(1.0, 20.0),
        stddev=(0.01, 20.0),
    )
    return model, {"amplitude": 2.7, "mean": 8.5, "stddev": 1.8}, {"output_unit": u.Jy}


def _const_case() -> tuple[Any, dict[str, float], dict[str, Any]]:
    model = _bounded(Const1D(amplitude=0.5), amplitude=(0.01, 5.0))
    return model, {"amplitude": 1.2}, {"output_unit": u.Jy}


CURATED_CASES = {
    "BlackBody": _blackbody_case,
    "PowerLaw1D": _power_law_case,
    "BrokenPowerLaw1D": _broken_power_law_case,
    "Polynomial1D": _polynomial_case,
    "Gaussian1D": _gaussian_case,
    "Const1D": _const_case,
}


class TestConformanceEachCuratedLeafAgreesWithTheBlackBoxAdapter:
    @pytest.mark.parametrize("name", sorted(CURATED_CASES))
    def test_it_agrees_at_cross_backend_tolerance(self, kit: Kit, name: str) -> None:
        model, values, options = CURATED_CASES[name]()
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, **options
        )
        adapted = adapt_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum, **options)
        got = native(**values)["default"].values
        expected = adapted(**values)["default"].values
        assert got == pytest.approx(expected, rel=CROSS_BACKEND_TOLERANCE, abs=0.0)

    @pytest.mark.parametrize("name", sorted(CURATED_CASES))
    def test_it_declares_the_capability_flags(self, kit: Kit, name: str) -> None:
        model, _, options = CURATED_CASES[name]()
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, **options
        )
        assert native.DIFFERENTIABLE is True
        assert native.BACKEND == kit.name


# ---------------------------------------------------------------------------
# Compound expressions: +, -, *, /
# ---------------------------------------------------------------------------


def _compound_sum() -> tuple[Any, dict[str, float]]:
    model = Gaussian1D(1.0, 8.0, 2.0) + Const1D(0.3)
    for name in ("amplitude_0", "mean_0", "stddev_0", "amplitude_1"):
        getattr(model, name).bounds = (0.001, 20.0)
    return model, {"amplitude_0": 2.4, "mean_0": 8.3, "stddev_0": 1.9, "amplitude_1": 0.5}


def _compound_difference() -> tuple[Any, dict[str, float]]:
    model = Gaussian1D(3.0, 8.0, 2.0) - Const1D(0.1)
    for name in ("amplitude_0", "mean_0", "stddev_0", "amplitude_1"):
        getattr(model, name).bounds = (0.001, 20.0)
    return model, {"amplitude_0": 2.8, "mean_0": 7.6, "stddev_0": 2.1, "amplitude_1": 0.15}


def _compound_product() -> tuple[Any, dict[str, float]]:
    model = Gaussian1D(3.0, 8.0, 2.0) * PowerLaw1D(1.0, 1.0, 1.0)
    for name in ("amplitude_0", "mean_0", "stddev_0", "amplitude_1", "alpha_1"):
        getattr(model, name).bounds = (0.001, 20.0)
    model.x_0_1.fixed = True
    return model, {
        "amplitude_0": 2.5,
        "mean_0": 8.1,
        "stddev_0": 1.7,
        "amplitude_1": 1.1,
        "alpha_1": 0.8,
    }


def _compound_ratio() -> tuple[Any, dict[str, float]]:
    model = Gaussian1D(3.0, 8.0, 2.0) / Const1D(2.0)
    for name in ("amplitude_0", "mean_0", "stddev_0", "amplitude_1"):
        getattr(model, name).bounds = (0.001, 20.0)
    return model, {"amplitude_0": 2.6, "mean_0": 7.9, "stddev_0": 2.2, "amplitude_1": 1.8}


COMPOUND_CASES = {
    "sum": _compound_sum,
    "difference": _compound_difference,
    "product": _compound_product,
    "ratio": _compound_ratio,
}


class TestConformanceCompoundExpressions:
    @pytest.mark.parametrize("name", sorted(COMPOUND_CASES))
    def test_it_agrees_at_cross_backend_tolerance(self, kit: Kit, name: str) -> None:
        model, values = COMPOUND_CASES[name]()
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy
        )
        adapted = adapt_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy)
        got = native(**values)["default"].values
        expected = adapted(**values)["default"].values
        assert got == pytest.approx(expected, rel=CROSS_BACKEND_TOLERANCE, abs=0.0)


# ---------------------------------------------------------------------------
# describe(): the Polynomial1D domain/window extra reaches provenance
# ---------------------------------------------------------------------------


class TestDescribe:
    def test_a_rescaled_polynomial_domain_reaches_describe(self, kit: Kit) -> None:
        model = Polynomial1D(degree=1, c0=0.0, c1=1.0, domain=(0.0, 10.0), window=(-1.0, 1.0))
        model.c0.bounds = (-5.0, 5.0)
        model.c1.bounds = (-5.0, 5.0)
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy
        )
        info = native.describe()
        assert info["backend"] == kit.name
        assert info["leaf_configuration"]["0"] == {"domain": (0.0, 10.0), "window": (-1.0, 1.0)}


# ---------------------------------------------------------------------------
# The native surface a realisation composes (W2.13): grid()/flux()
# ---------------------------------------------------------------------------


class TestTheNativeSurface:
    def test_grid_and_flux_agree_with_evaluate(self, kit: Kit) -> None:
        model, values = _compound_sum()
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy
        )
        grid = np.asarray(native.grid(native.channel))
        assert grid.shape == WAVELENGTH.shape
        flux = np.asarray(native.flux(native.channel, values))
        expected = native(**values)["default"].values
        assert flux == pytest.approx(expected)

    def test_flux_on_an_unknown_channel_is_refused_by_name(self, kit: Kit) -> None:
        model, values = _compound_sum()
        native = kit.module.from_astropy(
            model, grid=WAVELENGTH * u.micron, kind=Spectrum, output_unit=u.Jy
        )
        with pytest.raises(Exception, match="other"):
            native.flux("other", values)


# ---------------------------------------------------------------------------
# 2. NUTS on a translated compound model
# ---------------------------------------------------------------------------


TRUTH = {"amplitude_0": 3.0, "mean_0": 8.0, "stddev_0": 2.0, "amplitude_1": 0.5}
GRID = np.geomspace(1.0, 20.0, 40)
SIGMA = 0.05


def _truth_flux(grid: np.ndarray) -> np.ndarray:
    return (
        TRUTH["amplitude_0"] * np.exp(-0.5 * ((grid - TRUTH["mean_0"]) / TRUTH["stddev_0"]) ** 2)
        + TRUTH["amplitude_1"]
    )


def _noisy_data() -> Spectrum:
    rng = np.random.default_rng(20260913)
    return Spectrum(
        GRID * u.micron,
        (_truth_flux(GRID) + rng.normal(0.0, SIGMA, GRID.size)) * u.Jy,
        uncertainty=np.full(GRID.size, SIGMA) * u.Jy,
    )


DATA = _noisy_data()


def _compound_model() -> Any:
    model = Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)
    model.amplitude_0.bounds = (0.1, 10.0)
    model.mean_0.bounds = (1.0, 20.0)
    model.stddev_0.bounds = (0.1, 10.0)
    model.amplitude_1.bounds = (0.0, 5.0)
    return model


class TestNUTSOnATranslatedCompoundModel:
    def test_it_recovers_the_injected_truth(self, kit: Kit) -> None:
        native = kit.module.from_astropy(
            _compound_model(), grid=GRID * u.micron, kind=Spectrum, output_unit=u.Jy
        )
        problem = FittingProblem(
            native,
            [Dataset(DATA, likelihood=Likelihood(GaussianFamily(), kit.module.IndependentNoise()))],
            seed=20260913,
        )
        assert problem.differentiable is True
        assert problem.backend == kit.name
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            run = NUTSEngine(problem).run(draws=250, warmup=250, chains=2)
        for key, truth_key in (
            ("model.amplitude_0", "amplitude_0"),
            ("model.mean_0", "mean_0"),
            ("model.stddev_0", "stddev_0"),
            ("model.amplitude_1", "amplitude_1"),
        ):
            draws = np.asarray(run["posterior"][key])
            spread = float(draws.std())
            assert abs(float(draws.mean()) - TRUTH[truth_key]) < 5.0 * spread
        assert run.attrs["ampere_backend"] == kit.name
        assert run.attrs["ampere_engine"] == "nuts"


# ---------------------------------------------------------------------------
# 3. The two refusals this item adds
# ---------------------------------------------------------------------------


class TestTheTiedParameterRefusal:
    def test_a_tied_compound_model_is_refused_on_both_backends(self, kit: Kit) -> None:
        model = Gaussian1D(3.0, 6.0, 1.5) + Const1D(0.3)
        model.mean_0.tied = lambda m: float(m.amplitude_0.value) * 2.0
        for name in ("amplitude_0", "stddev_0", "amplitude_1"):
            getattr(model, name).bounds = (0.01, 12.0)
        with pytest.raises(CapabilityError, match="mean_0") as raised:
            kit.module.from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)
        assert isinstance(raised.value, NotImplementedError)
        assert f"ampere.backends.{kit.name}.from_astropy()" in str(raised.value)


class TestTheUnsupportedOperatorRefusal:
    def test_a_piped_compound_model_is_refused_by_name(self, kit: Kit) -> None:
        left = Gaussian1D(1.0, 2.0, 0.5)
        for name in left.param_names:
            getattr(left, name).bounds = (0.001, 20.0)
        right = Const1D(0.3)
        right.amplitude.bounds = (0.001, 20.0)
        piped = left | right
        with pytest.raises(CapabilityError, match="'\\|'") as raised:
            kit.module.from_astropy(piped, grid=WAVELENGTH * u.micron, kind=Spectrum)
        assert f"ampere.backends.{kit.name}.from_astropy()" in str(raised.value)
