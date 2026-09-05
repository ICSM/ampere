"""Exceptions and warnings raised by the engine drivers.

Deliberately small, and deliberately *here* rather than in
``ampere/core/exceptions.py``: that module is a frozen §4 contract, and
widening it is a decision-log matter (ground rule 9). This mirrors what
``ampere/results/exceptions.py`` did for :class:`~ampere.core.exceptions.
ResultsError` before Peter's 2026-09-03 ruling on ``results.md`` §15 R5 moved
it down; if the same judgement applies here, the move is two lines and no
import site changes, because :class:`EngineError` is re-exported from
``ampere.inference``.
"""

from __future__ import annotations

from ampere.core.exceptions import ContractError

__all__ = ["EngineError", "SamplingFailureWarning"]


class EngineError(ContractError):
    """An engine driver was asked for a run it cannot set up or emit.

    Raised for settings an engine genuinely cannot honour — fewer walkers than
    an ensemble move needs, a burn-in longer than the run, a problem with no
    free parameters to sample — and for a start point that could not be found
    inside the prior's support.

    Like every sibling contract error it is a :class:`ValueError`, and the
    contrast with :class:`~ampere.core.dataset.Failure` is the same one
    ``inference.md`` §11 draws: a *proposal* that cannot be scored is ``-inf``
    with a recorded reason and never an exception, whereas a *run* that cannot
    be set up should stop before it starts.
    """


class SamplingFailureWarning(UserWarning):
    """Some draws in a completed run could not be scored.

    Carries :meth:`~ampere.core.dataset.FittingProblem.failure_summary`'s
    aggregate — one line per failure class with the range of parameter values
    it happened over, never one warning per failed draw
    (``inference.md`` §11's "counting, not just recording").

    A distinct class so that a caller can silence it, promote it to an error
    (``warnings.simplefilter("error", SamplingFailureWarning)``, which is a
    reasonable thing to do in a pipeline), or assert on it in a test, without
    catching every other :class:`UserWarning` the stack emits.
    """
