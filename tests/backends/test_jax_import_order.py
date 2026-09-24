"""W5.31: the jax backend imports cleanly after arviz (the lazy-``numpyro`` fault).

arviz 1.3 puts a lazily-loaded ``numpyro`` stub in ``sys.modules``; if this
backend's submodules then imported ``numpyro.distributions`` before anything
touched the stub, numpyro's ``__init__`` ran mid-chain and
``numpyro.distributions.distribution`` never bound, so ``numpyro.factor``
failed inside every native model. ``ampere/backends/jax/__init__.py`` now
touches ``numpyro`` first. Each row runs in a **fresh interpreter** with this
checkout first on the path, because the fault is an import-order property of
a process, not of a module.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("jax")
pytest.importorskip("arviz")

_ROOT = Path(__file__).resolve().parents[2]

_CORE_THEN_ARVIZ_THEN_BACKEND = (
    "import sys\n"
    "import ampere.core\n"
    "import arviz\n"
    "import ampere.backends.jax\n"
    "d = sys.modules['numpyro.distributions']\n"
    "assert hasattr(d, 'distribution'), type(sys.modules['numpyro']).__name__\n"
    "import numpyro\n"
    "assert type(numpyro).__name__ == 'module', type(numpyro).__name__\n"
    "numpyro.distributions.distribution.Unit  # what numpyro.factor reaches for\n"
)


def _fresh(code: str) -> subprocess.CompletedProcess[str]:
    env = {**os.environ, "PYTHONPATH": str(_ROOT)}
    return subprocess.run(
        [sys.executable, "-c", code],
        cwd=_ROOT,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )


class TestTheBackendAfterArviz:
    def test_numpyro_is_fully_loaded_when_arviz_stubbed_it_first(self) -> None:
        """The exact order that failed at W5.26's review: core, arviz, then this backend."""
        result = _fresh(_CORE_THEN_ARVIZ_THEN_BACKEND)
        assert result.returncode == 0, result.stderr[-2000:]

    def test_the_touch_is_the_reason(self) -> None:
        """The stub really is lazy after arviz alone: the row above is not vacuous."""
        result = _fresh(
            "import sys, ampere.core, arviz\n"
            "assert type(sys.modules['numpyro']).__name__ == '_LazyModule', "
            "type(sys.modules['numpyro']).__name__\n"
        )
        assert result.returncode == 0, result.stderr[-2000:]
