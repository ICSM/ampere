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
import ampere.inference._blackjax
import ampere.inference._dynesty
import ampere.inference._emcee
import ampere.inference._nested
import ampere.inference._nuts
import ampere.inference._vi
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
    backend is installed. ``VIEngine``'s docstring demonstrates the *driver*
    rather than a library, and it has to name *some* backend to build a
    problem at all, so it names torch and this row is gated on torch — even
    though ``VARIATIONAL_LIBRARIES`` gained numpyro at W2.5 slice 2. The jax
    route is covered by ``tests/inference/test_vi.py``, which is parametrised
    over every installed backend the driver supports; a second, jax-flavoured
    copy of the same example would be a duplicate rather than coverage.
    """
    pytest.importorskip("torch")
    pytest.importorskip("pyro")
    results = doctest.testmod(ampere.inference._vi, optionflags=OPTIONS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


def test_the_nautilus_example_runs() -> None:
    """Separate, and skipped where the ``nautilus`` extra is absent (W5.14).

    The same reasoning the two rows above state. ``UltranestEngine`` carries
    no worked example of its own deliberately: the two drivers share a module
    and a run shape, so a second copy would cost another real nested-sampling
    fit in this suite to demonstrate the same thing, and
    ``tests/inference/test_nested.py`` — the engine battery — runs both
    against a closed-form evidence rather than against a docstring.
    """
    pytest.importorskip("nautilus")
    results = doctest.testmod(ampere.inference._nested, optionflags=OPTIONS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0


def test_the_blackjax_example_runs() -> None:
    """Separate, and skipped where the ``blackjax`` extra is absent (W5.14).

    The same reasoning the three rows above state, with one difference worth
    recording: ``BlackjaxEngine`` is jax-only *by nature* rather than by what
    happens to be installed, so there is no second route this example could
    have been written against. It demonstrates ``method="pathfinder"``,
    because that is the cheaper of the two — one L-BFGS path rather than a
    tuned chain — and the engine battery
    (``tests/inference/test_blackjax.py``) runs both methods against a
    closed-form posterior.
    """
    pytest.importorskip("jax")
    pytest.importorskip("blackjax")
    results = doctest.testmod(ampere.inference._blackjax, optionflags=OPTIONS, verbose=False)
    assert results.failed == 0
    assert results.attempted > 0
