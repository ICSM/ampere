"""The quickstart's straight line, written once for three backends (W7.13).

The worked example of the guide "Writing a model once for three backends"
(``docs/source/notebooks/portable_model.ipynb``). :mod:`.model` holds the
model, written against :class:`ampere.core.PortableModel` and its array
namespace; :mod:`.model_jax` and :mod:`.model_torch` hold its one-line native
twins. Importing this package needs neither jax nor torch -- import a twin's
module only when you build a problem on its backend.

The same source is held to a numpy oracle, to itself across backends and to
NUTS on both native backends by the conformance battery
(``tests/conformance/test_portable_model.py``), which is why the guide shows
this file rather than a copy of it.
"""

from __future__ import annotations

from .model import LinearModel

__all__ = ["LinearModel"]
