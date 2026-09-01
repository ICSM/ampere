"""Execute the worked examples in the contract specs.

W1.3's, W1.4's and W1.5's acceptance criteria are that every worked example in
``docs/design/contracts/parameters.md``,
``docs/design/contracts/results_schema.md`` and
``docs/design/contracts/transformations.md`` runs and passes, so no spec can
drift from the module it documents without this suite going red. Each markdown
file is fed straight to :mod:`doctest`; its examples share one namespace and
build on each other, exactly as a reader would run them.

The module docstrings of ``ampere.core`` are covered here too, for the same
reason.
"""

from __future__ import annotations

import doctest
from pathlib import Path

import pytest

import ampere.core.exceptions
import ampere.core.parameter
import ampere.core.results_schema
import ampere.core.transform

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACTS = REPO_ROOT / "docs" / "design" / "contracts"
SPEC = CONTRACTS / "parameters.md"
RESULTS_SPEC = CONTRACTS / "results_schema.md"
TRANSFORM_SPEC = CONTRACTS / "transformations.md"

# IGNORE_EXCEPTION_DETAIL is deliberately *not* set: the spec quotes ampere's
# error messages verbatim as part of its "fail loudly and specifically"
# argument, so those messages are part of what this test checks.
OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


@pytest.mark.parametrize(
    "spec",
    [SPEC, RESULTS_SPEC, TRANSFORM_SPEC],
    ids=["parameters", "results_schema", "transformations"],
)
def test_spec_document_exists(spec: Path) -> None:
    assert spec.is_file(), f"contract spec missing at {spec}"


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


def test_results_schema_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.4 contract spec executes as written."""
    results = doctest.testfile(
        str(RESULTS_SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    assert results.attempted > 40, f"only {results.attempted} examples found in {RESULTS_SPEC}"


def test_transformations_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.5 contract spec executes as written."""
    results = doctest.testfile(
        str(TRANSFORM_SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    assert results.attempted > 40, f"only {results.attempted} examples found in {TRANSFORM_SPEC}"


@pytest.mark.parametrize(
    "module",
    [
        ampere.core.parameter,
        ampere.core.results_schema,
        ampere.core.transform,
        ampere.core.exceptions,
    ],
    ids=["parameter", "results_schema", "transform", "exceptions"],
)
def test_module_docstring_examples_run(module: object) -> None:
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
    assert results.attempted > 0
