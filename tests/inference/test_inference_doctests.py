"""Execute the worked examples in ``ampere.inference``.

Same arrangement as ``tests/results/test_results_doctests.py`` and
``tests/core/test_spec_doctests.py``: the examples in each module are run as
written, so documentation that has drifted from the code fails the suite rather
than misleading a reader.

These examples are unusually load-bearing for docstrings, because each one is a
complete fit: it composes a problem out of a **hand-written** model — not one
taken from ``ampere.backends`` — drives it with one of the three engines, and
asserts on the emitted run. That is ``inference.md`` §10's portability claim
demonstrated at the smallest possible scale, and it is what stops the examples
quietly becoming pseudocode.

They are real sampling runs, so this module is a few seconds rather than
milliseconds. The budgets are the smallest that still recover the truth.
"""

from __future__ import annotations

import doctest

import pytest

import ampere.inference
import ampere.inference._dynesty
import ampere.inference._emcee
import ampere.inference._zeus
import ampere.inference.engine
import ampere.inference.exceptions

OPTIONS = doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE


@pytest.mark.parametrize(
    "module",
    [
        ampere.inference,
        ampere.inference._emcee,
        ampere.inference._dynesty,
        ampere.inference._zeus,
    ],
    ids=["package", "emcee", "dynesty", "zeus"],
)
def test_module_docstring_examples_run(module: object) -> None:
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
    assert results.attempted > 0


@pytest.mark.parametrize(
    "module",
    [ampere.inference.engine, ampere.inference.exceptions],
    ids=["engine", "exceptions"],
)
def test_prose_only_modules_have_no_failing_examples(module: object) -> None:
    """The two modules that are declaration rather than demonstration."""
    results = doctest.testmod(module, optionflags=OPTIONS, verbose=False)  # type: ignore[arg-type]
    assert results.failed == 0
