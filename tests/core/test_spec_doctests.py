"""Execute the worked examples in the contract specs.

W1.3's acceptance criterion is that every worked example in
``docs/design/contracts/parameters.md`` runs and passes, so the spec cannot
drift from ``ampere.core.parameter`` without this suite going red. The
markdown file is fed straight to :mod:`doctest`; its examples share one
namespace and build on each other, exactly as a reader would run them.

The module docstrings of ``ampere.core`` are covered here too, for the same
reason.
"""

from __future__ import annotations

import doctest
from pathlib import Path

import pytest

import ampere.core.exceptions
import ampere.core.parameter

REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC = REPO_ROOT / "docs" / "design" / "contracts" / "parameters.md"

# IGNORE_EXCEPTION_DETAIL is deliberately *not* set: the spec quotes ampere's
# error messages verbatim as part of its "fail loudly and specifically"
# argument, so those messages are part of what this test checks.
OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


def test_spec_document_exists() -> None:
    assert SPEC.is_file(), f"contract spec missing at {SPEC}"


def test_parameters_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.3 contract spec executes as written."""
    results = doctest.testfile(
        str(SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    # A spec with no runnable examples would pass vacuously; it must not.
    assert results.attempted > 40, f"only {results.attempted} examples found in {SPEC}"


@pytest.mark.parametrize(
    "module", [ampere.core.parameter, ampere.core.exceptions], ids=["parameter", "exceptions"]
)
def test_module_docstring_examples_run(module: object) -> None:
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
    assert results.attempted > 0
