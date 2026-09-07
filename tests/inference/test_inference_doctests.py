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
import ampere.inference._nuts
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


def test_the_nuts_example_runs() -> None:
    """Separate, and skipped where the ``jax`` extra is absent (W2.5).

    The other three drivers' examples run everywhere, because emcee and
    dynesty are base dependencies and zeus's import is lazy enough that the
    *docstring* still runs. NUTS is different: its example must build a jax
    model to have a gradient at all, so the example is a real fit on a real
    backend and can only run where that backend is installed. Skipping is the
    honest answer -- pretending otherwise would mean either a pseudocode
    example or a red suite in the daily-use environment.
    """
    pytest.importorskip("jax")
    pytest.importorskip("numpyro")
    results = doctest.testmod(ampere.inference._nuts, optionflags=OPTIONS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


def test_the_vi_example_runs() -> None:
    """Separate, and skipped where the ``torch`` extra is absent (W2.4 slice 2).

    The same reasoning the NUTS row states, one library along: the example is
    a real variational fit on a real backend, so it can only run where that
    backend is installed. ``VARIATIONAL_LIBRARIES`` is torch-only today; when
    the jax track adds numpyro's ``SVI`` the example stays as it is and this
    row stays as it is, because the docstring demonstrates the *driver* rather
    than a library.
    """
    pytest.importorskip("torch")
    pytest.importorskip("pyro")
    import ampere.inference._vi

    results = doctest.testmod(ampere.inference._vi, optionflags=OPTIONS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0
