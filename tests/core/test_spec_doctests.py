"""Execute the worked examples in the contract specs.

W1.3's, W1.4's, W1.5's, W1.6's and W1.7's acceptance criteria are that every
worked example in ``docs/design/contracts/parameters.md``,
``docs/design/contracts/results_schema.md``,
``docs/design/contracts/transformations.md``,
``docs/design/contracts/likelihoods.md`` and
``docs/design/contracts/inference.md`` runs and passes, so no spec can
drift from the module it documents without this suite going red. Each markdown
file is fed straight to :mod:`doctest`; its examples share one namespace and
build on each other, exactly as a reader would run them.

The module docstrings of ``ampere.core`` are covered here too, for the same
reason, and (W6.3) so is the user documentation: every ``.. code-block::
pycon`` under ``docs/source`` is extracted and run, one namespace per page.
"""

from __future__ import annotations

import doctest
import importlib.util
import warnings
from pathlib import Path

import pytest

import ampere.core.dataset
import ampere.core.exceptions
import ampere.core.kernels
import ampere.core.likelihood
import ampere.core.lowering
import ampere.core.parameter
import ampere.core.results_schema
import ampere.core.rng
import ampere.core.transform

REPO_ROOT = Path(__file__).resolve().parents[2]
CONTRACTS = REPO_ROOT / "docs" / "design" / "contracts"
SPEC = CONTRACTS / "parameters.md"
RESULTS_SPEC = CONTRACTS / "results_schema.md"
TRANSFORM_SPEC = CONTRACTS / "transformations.md"
LIKELIHOOD_SPEC = CONTRACTS / "likelihoods.md"
INFERENCE_SPEC = CONTRACTS / "inference.md"

# IGNORE_EXCEPTION_DETAIL is deliberately *not* set: the spec quotes ampere's
# error messages verbatim as part of its "fail loudly and specifically"
# argument, so those messages are part of what this test checks.
OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


@pytest.mark.parametrize(
    "spec",
    [SPEC, RESULTS_SPEC, TRANSFORM_SPEC, LIKELIHOOD_SPEC, INFERENCE_SPEC],
    ids=["parameters", "results_schema", "transformations", "likelihoods", "inference"],
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


def test_likelihood_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.6 contract spec executes as written."""
    results = doctest.testfile(
        str(LIKELIHOOD_SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    assert results.attempted > 40, f"only {results.attempted} examples found in {LIKELIHOOD_SPEC}"


def test_inference_spec_examples_run() -> None:
    """Every ``>>>`` example in the W1.7 contract spec executes as written."""
    results = doctest.testfile(
        str(INFERENCE_SPEC),
        module_relative=False,
        optionflags=OPTIONS,
        verbose=False,
        report=True,
    )
    assert results.failed == 0, f"{results.failed} of {results.attempted} spec examples failed"
    assert results.attempted > 40, f"only {results.attempted} examples found in {INFERENCE_SPEC}"


@pytest.mark.parametrize(
    "module",
    [
        ampere.core.parameter,
        ampere.core.results_schema,
        ampere.core.transform,
        ampere.core.kernels,
        ampere.core.likelihood,
        ampere.core.dataset,
        ampere.core.rng,
        ampere.core.exceptions,
        ampere.core.lowering,
    ],
    ids=[
        "parameter",
        "results_schema",
        "transform",
        "kernels",
        "likelihood",
        "dataset",
        "rng",
        "exceptions",
        "lowering",
    ],
)
def test_module_docstring_examples_run(module: object) -> None:
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
    assert results.attempted > 0


# ---------------------------------------------------------------------------
# The user documentation (W6.3): every ``.. code-block:: pycon`` under
# ``docs/source`` is run.
# ---------------------------------------------------------------------------

DOCS_SOURCE = REPO_ROOT / "docs" / "source"

# A page whose blocks need something this environment lacks is skipped here,
# with the reason, rather than dropped from the page: the blocks stay where a
# reader finds them and the skip shows in ``pytest -rs``. Empty today -- every
# page's blocks run in the ``dev`` environment.
SKIP: dict[str, str] = {}


def _pycon_blocks(text: str) -> list[tuple[str, str | None]]:
    """A page's ``.. code-block:: pycon`` directives, in order.

    Each is ``(body, requires)``: ``requires`` is the library named by a
    ``:class: needs-<library>`` option on the directive (Sphinx renders the
    class harmlessly), or ``None``.
    """
    lines = text.splitlines()
    blocks: list[tuple[str, str | None]] = []
    i = 0
    while i < len(lines):
        if lines[i].strip() == ".. code-block:: pycon":
            indent = len(lines[i]) - len(lines[i].lstrip())
            i += 1
            requires = None
            body: list[str] = []
            while i < len(lines) and (
                not lines[i].strip() or lines[i].startswith(" " * (indent + 1))
            ):
                stripped = lines[i].strip()
                if stripped.startswith(":class: needs-") and not body:
                    requires = stripped.removeprefix(":class: needs-")
                else:
                    body.append(lines[i][indent + 4 :] if stripped else "")
                i += 1
            blocks.append(("\n".join(body).strip("\n"), requires))
        else:
            i += 1
    return blocks


def _pycon_pages() -> list[Path]:
    return sorted(
        page for page in DOCS_SOURCE.glob("*.rst") if ".. code-block:: pycon" in page.read_text()
    )


@pytest.mark.parametrize("page", _pycon_pages(), ids=lambda page: page.name)
def test_docs_page_pycon_blocks_run(
    page: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The ``pycon`` blocks of a documentation page run as written, in order.

    One namespace per page, so a later block builds on an earlier one exactly
    as the page reads.
    """
    if page.name in SKIP:
        pytest.skip(SKIP[page.name])
    # A page may write a file (an optimum to NetCDF, say) and may import the
    # repository's ``examples`` package: run from a scratch directory, with the
    # repository root importable, so nothing lands in the working tree.
    monkeypatch.chdir(tmp_path)
    monkeypatch.syspath_prepend(str(REPO_ROOT))
    found = _pycon_blocks(page.read_text())
    assert found, f"{page.name} has a pycon directive the extractor could not read"
    blocks = []
    for body, requires in found:
        # A block marked ``:class: needs-torch`` (or -jax) runs only where that
        # library is installed; elsewhere it is dropped and the drop is reported.
        if requires is not None and importlib.util.find_spec(requires) is None:
            warnings.warn(f"{page.name}: a block needing {requires} was not run here", stacklevel=1)
            continue
        blocks.append(body)
    test = doctest.DocTestParser().get_doctest("\n\n".join(blocks), {}, page.name, str(page), 0)
    runner = doctest.DocTestRunner(optionflags=OPTIONS, verbose=False)
    runner.run(test)
    results = runner.summarize(verbose=False)
    assert results.failed == 0, (
        f"{results.failed} of {results.attempted} examples failed in {page.name}"
    )
    assert results.attempted > 0, f"no examples found in {page.name}"


def test_docs_pycon_pages_are_collected() -> None:
    """The glob finds at least the pages known to carry ``pycon`` blocks."""
    names = {page.name for page in _pycon_pages()}
    assert {
        "kernels.rst",
        "solvers.rst",
        "astropy.rst",
        "optimisers.rst",
        "population.rst",
        "sed_composition.rst",
    } <= names
