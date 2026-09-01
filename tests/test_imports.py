"""Import-quarantine regression test (W0.5, issues #74-77 + W0.4 findings).

"Lazy-import quarantine" means: a missing or broken optional dependency
anywhere under ``ampere/`` must never block importing ``ampere`` or any
other, unrelated module -- the broken/heavy corner should fail, loudly and
informatively, only when something actually tries to *use* it.

This module partitions every importable file under ``ampere/`` into four
buckets and checks each one appropriately:

* **plain** (the default): must import cleanly in *any* environment,
  including the ``pip install -e .`` minimal environment with none of the
  optional extras installed.
* **optional-extra-gated** (:data:`OPTIONAL_EXTRA_MODULES`): normal
  ``pip install ampere[extra]`` gating (``ampere[sbi]``, ``ampere[zeus]``)
  -- not "broken". Skipped (not failed) when the extra isn't installed;
  must import cleanly when it is.
* **quarantined** (:data:`QUARANTINED_MODULES`): modules with a heavy or
  unfinished dependency that isn't (yet) a declared ampere extra
  (``hyperion``/``dill``, ``Starfish``) or that were never completed
  (``CCMExtinctionLaw``). These must import cleanly in *any* environment
  (their heavy imports are lazily guarded), and the affected class must
  raise an informative error -- naming the missing dependency and/or the
  ampere issue -- when actually constructed without it.
* **excluded** (:data:`EXCLUDED_MODULES`): files under ``ampere/`` that are
  not part of the package's importable surface at all (standalone scripts
  with side effects, or requiring dependencies with no ampere extra at
  all) -- out-of-scope findings, not fixed by W0.5, listed here so they
  don't silently fall through the cracks of this test's coverage.
"""

from __future__ import annotations

import importlib
import importlib.util
import pkgutil

import numpy as np
import pytest

import ampere

# ---------------------------------------------------------------------------
# Discovery
# ---------------------------------------------------------------------------


def _iter_ampere_modules():
    """Dotted name of every module/subpackage under the ``ampere`` package.

    ``onerror`` is a no-op (rather than the default, which re-raises) so
    that if a *package* (as opposed to a leaf module) ever fails to import
    while walking, this test file still collects -- and the failure shows
    up as a missing entry in ``test_every_module_is_classified`` rather
    than an opaque collection error.
    """
    return sorted(
        module_info.name
        for module_info in pkgutil.walk_packages(
            ampere.__path__, prefix="ampere.", onerror=lambda name: None
        )
    )


ALL_MODULES = _iter_ampere_modules()

# ---------------------------------------------------------------------------
# Optional-extra-gated modules: ordinary ``ampere[extra]`` gating.
# ---------------------------------------------------------------------------
OPTIONAL_EXTRA_MODULES = {
    "ampere.infer.sbi": ("torch", "sbi"),  # ampere[sbi]
    "ampere.infer.zeussearch": ("zeus",),  # ampere[zeus]
}

# ---------------------------------------------------------------------------
# Quarantined (W0.5): must import cleanly always; the dependency-gated
# class inside raises an informative error on construction. See the NOTE
# comments in each module for the detailed rationale.
# ---------------------------------------------------------------------------
QUARANTINED_MODULES = {
    "ampere.models.Hyperion",  # HyperionCStarRTModel needs hyperion+dill
    "ampere.models.QuickSED",  # QuickSEDModel needs Starfish
    "ampere.models.extinctionModels",  # F99Extinction needs dust_extinction (#75);
    # CCMExtinctionLaw is an unfinished stub (#76)
}

# ---------------------------------------------------------------------------
# Excluded: not part of ampere's importable package surface. Out-of-scope
# findings (see the W0.5 report) -- not fixed here.
# ---------------------------------------------------------------------------
EXCLUDED_MODULES = {
    # A stray debug script, not a library module: unconditionally
    # constructs a DustySpectrum and writes a Dusty input file as a *side
    # effect of being imported*, and uses a bare `from Dusty import
    # DustySpectrum` that only ever worked when run directly as `python
    # test_Dusty.py` from within ampere/models/ (Dusty.py on sys.path[0]).
    # Not reachable from any package import chain.
    "ampere.models.test_Dusty",
    # Standalone filter-generation CLI scripts requiring 'h5py', which is
    # not an ampere dependency or extra. Not reachable from any package
    # import chain (ampere/utils/__init__.py imports nothing).
    "ampere.utils.makeFilterSet",
    "ampere.utils.makeFilterSet_mod",
}

PLAIN_MODULES = sorted(
    set(ALL_MODULES) - set(OPTIONAL_EXTRA_MODULES) - QUARANTINED_MODULES - EXCLUDED_MODULES
)


def test_every_ampere_module_is_classified():
    """Guard-rail so a newly added file can't silently dodge this test's
    coverage: every module discovered under ``ampere/`` must fall into
    exactly one of the four buckets above (``PLAIN_MODULES`` is "whatever
    is left over", so this mainly checks discovery itself still works and
    that the buckets don't overlap).
    """
    assert PLAIN_MODULES, "module discovery under ampere/ found nothing -- is this broken?"
    buckets = [PLAIN_MODULES, OPTIONAL_EXTRA_MODULES, QUARANTINED_MODULES, EXCLUDED_MODULES]
    covered = set()
    for bucket in buckets:
        overlap = covered & set(bucket)
        assert not overlap, f"module(s) listed in more than one bucket: {overlap}"
        covered |= set(bucket)
    assert covered == set(ALL_MODULES)


# ---------------------------------------------------------------------------
# Top-level acceptance check (W0.5 Accept #1, verbatim)
# ---------------------------------------------------------------------------


def test_top_level_packages_import_cleanly():
    import ampere  # noqa: F401
    import ampere.models  # noqa: F401
    import ampere.infer  # noqa: F401
    import ampere.data  # noqa: F401


# ---------------------------------------------------------------------------
# Plain modules: must import cleanly, always.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("module_name", PLAIN_MODULES)
def test_plain_module_imports_cleanly(module_name):
    importlib.import_module(module_name)


# ---------------------------------------------------------------------------
# Optional-extra-gated modules.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("module_name", sorted(OPTIONAL_EXTRA_MODULES))
def test_optional_extra_module_imports_when_available(module_name):
    deps = OPTIONAL_EXTRA_MODULES[module_name]
    missing = [dep for dep in deps if importlib.util.find_spec(dep) is None]
    if missing:
        pytest.skip(f"optional extra dependency not installed: {missing}")
    importlib.import_module(module_name)


# ---------------------------------------------------------------------------
# Quarantined modules: must import cleanly regardless of environment.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("module_name", sorted(QUARANTINED_MODULES))
def test_quarantined_module_imports_cleanly(module_name):
    importlib.import_module(module_name)


def test_hyperion_raises_informative_error_without_hyperion():
    from ampere.models.Hyperion import _HYPERION_IMPORT_ERROR, HyperionCStarRTModel

    if _HYPERION_IMPORT_ERROR is None:
        pytest.skip("hyperion and dill are installed in this environment")
    with pytest.raises(ImportError, match="hyperion"):
        HyperionCStarRTModel(np.linspace(1.0, 10.0, 5))


def test_quicksed_raises_informative_error_without_starfish():
    from ampere.models.QuickSED import _STARFISH_IMPORT_ERROR, QuickSEDModel

    if _STARFISH_IMPORT_ERROR is None:
        pytest.skip("Starfish is installed in this environment")
    with pytest.raises(ImportError, match="Starfish"):
        QuickSEDModel(np.linspace(1.0, 10.0, 5))


def test_f99extinction_raises_informative_error_without_dust_extinction():
    from ampere.models.extinctionModels import F99Extinction

    if importlib.util.find_spec("dust_extinction") is not None:
        pytest.skip("dust_extinction is installed in this environment")
    with pytest.raises(ImportError, match="dust_extinction"):
        F99Extinction(np.linspace(0.3, 3.0, 10), av=1.0, rv=3.1)


def test_f99extinction_constructs_and_runs_when_dust_extinction_available():
    pytest.importorskip("dust_extinction")
    from ampere.models.extinctionModels import F99Extinction

    # F99's valid range is roughly 0.3-10 inverse-micron (~0.1-3.3 micron);
    # stay well inside it.
    wavelengths = np.linspace(0.3, 3.0, 20)
    model = F99Extinction(wavelengths)  # av, rv both free -> npars == 2
    model(1.0, 3.1)
    assert model.modelFlux.shape == wavelengths.shape
    # modelFlux is a transmission fraction (see the class docstring): in
    # (0, 1] for a positive Av.
    assert np.all(model.modelFlux > 0.0)
    assert np.all(model.modelFlux <= 1.0)


def test_ccmextinctionlaw_always_raises_informative_error():
    from ampere.models.extinctionModels import CCMExtinctionLaw

    with pytest.raises(NotImplementedError, match="#76"):
        CCMExtinctionLaw()
