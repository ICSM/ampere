"""The old top-level names of the legacy packages, served as aliases.

``ampere.data``, ``ampere.models``, ``ampere.infer``, ``ampere.utils`` and
``ampere.logger`` — and every module under them — resolve to the
corresponding module under :mod:`ampere.legacy` on first import. The alias
and the target are the *same* module object, so ``__name__``, pickles and
reprs say ``ampere.legacy...``; there is no warning and no removal date (see
the legacy page of the docs).
"""

from __future__ import annotations

import importlib
import importlib.abc
import importlib.machinery
import importlib.util
import sys
from collections.abc import Sequence
from types import ModuleType

LEGACY_NAMES: tuple[str, ...] = ("data", "models", "infer", "utils", "logger")
_PREFIXES: tuple[str, ...] = tuple(f"ampere.{name}" for name in LEGACY_NAMES)


def target_of(fullname: str) -> str | None:
    """``ampere.legacy.<rest>`` for an old name, else ``None``."""
    for prefix in _PREFIXES:
        if fullname == prefix or fullname.startswith(prefix + "."):
            return "ampere.legacy" + fullname[len("ampere") :]
    return None


class _AliasLoader(importlib.abc.Loader):
    def __init__(self, module: ModuleType) -> None:
        self._module = module

    def create_module(self, spec: importlib.machinery.ModuleSpec) -> ModuleType:
        return self._module  # the import system binds this object under the old name

    def exec_module(self, module: ModuleType) -> None:
        return None  # already executed under its real name


class LegacyAliasFinder(importlib.abc.MetaPathFinder):
    def find_spec(
        self,
        fullname: str,
        path: Sequence[str] | None = None,
        target: ModuleType | None = None,
    ) -> importlib.machinery.ModuleSpec | None:
        real_name = target_of(fullname)
        if real_name is None:
            return None
        module = importlib.import_module(real_name)  # ImportError propagates, as before
        return importlib.util.spec_from_loader(
            fullname, _AliasLoader(module), is_package=hasattr(module, "__path__")
        )


def install() -> None:
    """Insert the finder ahead of :class:`~importlib.machinery.PathFinder`, once.

    An *appended* finder loses a race for any submodule import under an
    already-aliased package. Once ``ampere.infer`` has resolved to the
    *same* object as ``ampere.legacy.infer`` -- same ``__path__``, since it
    is the same module -- a fresh, unqualified request for
    ``ampere.infer.emceesearch`` hands the finders that real ``__path__``,
    and :class:`~importlib.machinery.PathFinder` (tried before anything
    appended to the end of ``sys.meta_path``) happily finds
    ``emceesearch.py`` there directly and imports a *second*, distinct
    module under the alias name -- silently breaking the "same object"
    guarantee for every submodule the parent package does not import
    eagerly (verified: ``ampere.infer.emceesearch.EmceeSearch`` and
    ``ampere.legacy.infer.emceesearch.EmceeSearch`` were two different
    classes under the append placement). Inserted ahead of ``PathFinder``
    instead, this finder gets first refusal for the five legacy prefixes at
    every depth and always answers with the one real module object; every
    other import -- including ``ampere.legacy.*`` itself -- is untouched,
    because :func:`target_of` returns ``None`` for it and this finder
    simply declines, leaving ``PathFinder`` (or an earlier finder, e.g. an
    editable install's) to answer as before.
    """
    if any(isinstance(finder, LegacyAliasFinder) for finder in sys.meta_path):
        return
    finder = LegacyAliasFinder()
    for index, existing in enumerate(sys.meta_path):
        if existing is importlib.machinery.PathFinder:
            sys.meta_path.insert(index, finder)
            break
    else:
        sys.meta_path.append(finder)
