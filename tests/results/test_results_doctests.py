"""Execute the worked examples in the W1.8 contract spec and in ``ampere.results``.

The same arrangement ``tests/core/test_spec_doctests.py`` makes for W1.3 to W1.7:
each markdown file is fed straight to :mod:`doctest`, its examples sharing one
namespace and building on each other exactly as a reader would run them, so the
spec cannot drift from the module it documents without this suite going red.

It lives beside the other ``ampere.results`` tests rather than in
``tests/core/test_spec_doctests.py`` because it is the only spec whose examples
need an optional dependency: ``ampere.results`` is arviz-backed, and the suite
must skip rather than fail in a minimal install.
"""

from __future__ import annotations

import doctest
from pathlib import Path

import pytest

pytest.importorskip("arviz", reason="the results contract's examples need ampere[arviz]")

import ampere.results.derived
import ampere.results.emission
import ampere.results.exceptions
import ampere.results.plots
import ampere.results.provenance
import ampere.results.serialisation

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_SPEC = REPO_ROOT / "docs" / "design" / "contracts" / "results.md"

# IGNORE_EXCEPTION_DETAIL is deliberately *not* set, for the same reason the core
# suite gives: the spec quotes ampere's error messages verbatim as part of its
# "fail loudly and specifically" argument, so those messages are under test too.
OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


def test_spec_document_exists() -> None:
    assert RESULTS_SPEC.is_file(), f"contract spec missing at {RESULTS_SPEC}"


def test_results_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.8 contract spec executes as written."""
    results = doctest.testfile(
        str(RESULTS_SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    # A spec with no runnable examples would pass vacuously; it must not.
    assert results.attempted > 40, f"only {results.attempted} examples found in {RESULTS_SPEC}"


@pytest.mark.parametrize(
    "module",
    [
        ampere.results.provenance,
        ampere.results.serialisation,
        ampere.results.plots,
        ampere.results.exceptions,
    ],
    ids=["provenance", "serialisation", "plots", "exceptions"],
)
def test_module_docstring_examples_run(module: object) -> None:
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
    assert results.attempted > 0


@pytest.mark.parametrize(
    "module",
    [ampere.results.emission, ampere.results.derived],
    ids=["emission", "derived"],
)
def test_prose_only_modules_have_no_failing_examples(module: object) -> None:
    """The two modules that are declaration rather than demonstration."""
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
