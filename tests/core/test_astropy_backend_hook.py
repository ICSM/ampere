"""The opt-in translation hook's refusal rows (``astropy_compat.md`` §5), W4.6 and W4.7.

The 2026-09-01 ruling — *curated astropy→native translation is opt-in, never
silent* — has two halves. :func:`ampere.core.from_astropy` is the half that
never substitutes anything; ``ampere.backends.torch.from_astropy`` and its jax
twin are the half that substitutes **only** when asked and **only** what it has
a curated row for. W4.6 landed the hook with an **empty** table, so every call
refused, and the refusal was the whole of what there was to test. W4.7 fills
the table (``BlackBody``, ``PowerLaw1D``, ``BrokenPowerLaw1D``,
``Polynomial1D``, ``Gaussian1D``, ``Const1D`` and their compound sums,
products, differences and ratios — :mod:`ampere.core.astropy_translations`),
so what is left to test here is "everything outside the curated table still
refuses, by name" — the same refusal, now exercised against a model the table
genuinely has no row for (:class:`~astropy.modeling.functional_models.Sersic1D`)
rather than against every model.

The core half of the refusal — the leaf decomposition of a compound model, and
the message itself — is :func:`ampere.core.translation_refusal`, and is tested
here without either extra installed, because that is where the logic lives. The
two backend modules are tested in their own environments (``pixi run -e torch``
/ ``-e jax``), where they import; elsewhere those rows skip, which is the same
degradation ``tests/conformance/backends/__init__.py`` uses for the same
reason.
"""

from __future__ import annotations

import doctest
import importlib
import importlib.util
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.modeling.models import BlackBody, Const1D, Gaussian1D, PowerLaw1D, Sersic1D

import ampere.core.astropy_compat
from ampere.core import Spectrum, astropy_components, from_astropy, translation_refusal
from ampere.core.exceptions import CapabilityError

WAVELENGTH = np.geomspace(1.0, 20.0, 8)

#: The two backends that carry the hook. Named here and nowhere else in this
#: module, so adding a third is one entry.
BACKEND_MODULES = ("ampere.backends.torch", "ampere.backends.jax")


def _installed(module: str) -> bool:
    return importlib.util.find_spec(module.rsplit(".", 1)[-1]) is not None


class TestTheLeafDecomposition:
    """A compound model is translatable exactly when every one of its leaves is."""

    def test_a_simple_model_is_its_own_only_component(self) -> None:
        model = PowerLaw1D(amplitude=2.0, x_0=1.0, alpha=1.0)
        assert astropy_components(model) == (model,)

    def test_a_compound_model_decomposes_into_its_leaves(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)
        assert [type(part).__name__ for part in astropy_components(model)] == [
            "Gaussian1D",
            "Const1D",
        ]

    def test_a_nested_compound_model_decomposes_all_the_way_down(self) -> None:
        model = (Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)) * PowerLaw1D(1.0, 1.0, 1.0)
        assert [type(part).__name__ for part in astropy_components(model)] == [
            "Gaussian1D",
            "Const1D",
            "PowerLaw1D",
        ]

    def test_an_adapted_model_is_unwrapped_to_the_astropy_model_it_holds(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5)
        for name in model.param_names:
            getattr(model, name).bounds = (0.01, 10.0)
        adapted = from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)
        assert [type(part).__name__ for part in astropy_components(adapted)] == ["Gaussian1D"]


class TestTheRefusal:
    """What the hook says when it has no row — which, in W4.6, is always."""

    def test_it_names_every_untranslatable_component(self) -> None:
        model = Gaussian1D(1.0, 2.0, 0.5) + Const1D(0.3)
        message = str(translation_refusal("torch", model, {}))
        assert "Gaussian1D" in message
        assert "Const1D" in message

    def test_it_names_the_backend_and_the_black_box_route(self) -> None:
        message = str(translation_refusal("jax", PowerLaw1D(2.0, 1.0, 1.0), {}))
        assert "ampere.backends.jax.from_astropy()" in message
        assert "ampere.core.from_astropy()" in message
        assert "gradient-free engines and SBI, never NUTS or VI" in message

    def test_it_says_translation_is_opt_in_and_never_silent(self) -> None:
        """The ruling, in the message, so a reader meets it where it bites."""
        message = str(translation_refusal("torch", BlackBody(temperature=3000.0 * u.K), {}))
        assert "opt-in and never silent" in message
        assert "2026-09-01" in message

    def test_it_reports_what_the_table_currently_holds(self) -> None:
        empty = str(translation_refusal("torch", Const1D(0.3), {}))
        assert "nothing yet (W4.7 adds the curated rows)" in empty
        filled = str(translation_refusal("torch", Const1D(0.3), {Gaussian1D: object()}))
        assert "currently holds: Gaussian1D" in filled

    def test_it_is_a_capability_error_and_a_not_implemented_error(self) -> None:
        """Nothing is malformed; one path cannot serve a valid declaration (yet)."""
        refusal = translation_refusal("torch", Const1D(0.3), {})
        assert isinstance(refusal, CapabilityError)
        assert isinstance(refusal, NotImplementedError)


@pytest.mark.parametrize("module", BACKEND_MODULES)
class TestTheBackendHooks:
    """Each backend's hook, in the environment where its extra is installed."""

    def test_the_hook_is_exported_and_the_table_holds_the_curated_rows(self, module: str) -> None:
        """W4.7: the table is no longer empty, and holds exactly the curated six."""
        if not _installed(module):
            pytest.skip(f"needs the {module.rsplit('.', 1)[-1]!r} extra")
        backend = importlib.import_module(module)
        assert callable(backend.from_astropy)
        assert {cls.__name__ for cls in backend.TRANSLATIONS} == {
            "BlackBody",
            "PowerLaw1D",
            "BrokenPowerLaw1D",
            "Polynomial1D",
            "Gaussian1D",
            "Const1D",
        }

    def test_it_refuses_a_model_outside_the_curated_table(self, module: str) -> None:
        """A model the table genuinely has no row for still refuses, by name."""
        if not _installed(module):
            pytest.skip(f"needs the {module.rsplit('.', 1)[-1]!r} extra")
        backend = importlib.import_module(module)
        name = module.rsplit(".", 1)[-1]
        model = Sersic1D(amplitude=1.0, r_eff=1.0, n=4.0)
        with pytest.raises(CapabilityError) as raised:
            backend.from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)
        message = str(raised.value)
        assert f"ampere.backends.{name}.from_astropy()" in message
        assert "Sersic1D" in message

    def test_it_names_the_untranslatable_half_of_a_compound_model(self, module: str) -> None:
        """A curated leaf beside an uncurated one still refuses, naming only the latter."""
        if not _installed(module):
            pytest.skip(f"needs the {module.rsplit('.', 1)[-1]!r} extra")
        backend = importlib.import_module(module)
        model = Gaussian1D(1.0, 2.0, 0.5) + Sersic1D(amplitude=1.0, r_eff=1.0, n=4.0)
        with pytest.raises(CapabilityError, match="Sersic1D"):
            backend.from_astropy(model, grid=WAVELENGTH * u.micron, kind=Spectrum)


# ---------------------------------------------------------------------------
# The contract page and the module docstring
# ---------------------------------------------------------------------------

CONTRACT_PAGE = (
    Path(__file__).resolve().parents[2] / "docs" / "design" / "contracts" / "astropy_compat.md"
)

# IGNORE_EXCEPTION_DETAIL is deliberately *not* set, for the reason
# ``test_spec_doctests.py`` records: the page quotes ampere's refusals verbatim
# as part of its "fail loudly and specifically" argument, so those messages are
# part of what is checked.
DOCTEST_OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


class TestThisPage:
    """``docs/design/contracts/astropy_compat.md`` runs as written.

    The same discipline every other contract page is held to
    (``tests/core/test_spec_doctests.py``), applied here rather than there
    because W4.6 owns ``tests/core/test_astropy*`` and that file is shared.
    A page that drifted from the module would go red.
    """

    def test_the_contract_page_exists(self) -> None:
        assert CONTRACT_PAGE.is_file(), f"contract page missing at {CONTRACT_PAGE}"

    def test_every_worked_example_in_the_page_runs(self) -> None:
        results = doctest.testfile(
            str(CONTRACT_PAGE),
            module_relative=False,
            optionflags=DOCTEST_OPTIONS,
            verbose=False,
            report=True,
        )
        assert results.failed == 0, f"{results.failed} of {results.attempted} page examples failed"
        # A page with no runnable examples would pass vacuously; it must not.
        assert results.attempted > 15, f"only {results.attempted} examples found in {CONTRACT_PAGE}"

    def test_the_modules_own_docstring_examples_run(self) -> None:
        results = doctest.testmod(
            ampere.core.astropy_compat, optionflags=DOCTEST_OPTIONS, verbose=False
        )
        assert results.failed == 0
        assert results.attempted > 0
