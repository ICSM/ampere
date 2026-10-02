"""Every public name in the v2 namespaces carries a docstring (issue #58).

The docs build runs with warnings as errors (W6.3), which enforces that a
docstring *renders*; this module enforces that there is one to render. For
each v2 package, every name in ``__all__`` must resolve to an object whose
``__doc__`` is non-empty and whose first line is a sentence (it ends in a
full stop). The two array-library backends are checked the same way where
the library is installed, and skipped where it is not.
"""

from __future__ import annotations

import importlib
import importlib.util
import inspect
import pkgutil

import pytest

PACKAGES = [
    "ampere.core",
    "ampere.backends.reference",
    "ampere.inference",
    "ampere.results",
    "ampere.backends.torch",
    "ampere.backends.jax",
]


def _attribute_docs(package: str) -> set[str]:
    """Names that carry a ``#:`` or string-literal attribute docstring.

    A module-level constant has no ``__doc__`` of its own; what autodoc
    renders for it is the comment above the assignment, which Sphinx's own
    source analyser finds. This is the same analyser, over every module of
    the package, so the row asks exactly what the build asks.
    """
    from sphinx.pycode import ModuleAnalyzer

    documented: set[str] = set()
    root = importlib.import_module(package)
    modules = [package]
    if hasattr(root, "__path__"):
        modules += [m.name for m in pkgutil.walk_packages(root.__path__, package + ".")]
    for modname in modules:
        try:
            analyser = ModuleAnalyzer.for_module(modname)
            analyser.analyze()
        except Exception:  # pragma: no cover - an unimportable submodule
            continue
        documented |= {name for scope, name in analyser.attr_docs if scope == ""}
    return documented


@pytest.mark.parametrize("package", PACKAGES)
def test_public_names_are_documented(package: str) -> None:
    pytest.importorskip("sphinx")
    if package.endswith((".torch", ".jax")):
        lib = package.rsplit(".", 1)[1]
        if importlib.util.find_spec(lib) is None:
            pytest.skip(f"{lib} is not installed in this environment")
    module = importlib.import_module(package)
    names = getattr(module, "__all__", None)
    assert names, f"{package} declares no __all__"
    attribute_docs: set[str] | None = None
    bad: list[str] = []
    for name in names:
        obj = getattr(module, name)
        if inspect.isclass(obj) or inspect.isroutine(obj) or inspect.ismodule(obj):
            doc = (inspect.getdoc(obj) or "").strip()
            first = doc.splitlines()[0].strip() if doc else ""
            if not first.endswith("."):
                bad.append(
                    f"{package}.{name}: {first!r}" if doc else f"{package}.{name}: no docstring"
                )
        else:
            if attribute_docs is None:
                attribute_docs = _attribute_docs(package)
            if name not in attribute_docs:
                bad.append(f"{package}.{name}: constant with no #: attribute docstring")
    assert not bad, "undocumented or non-sentence docstring:\n" + "\n".join(bad)
