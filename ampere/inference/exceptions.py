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

__all__ = [
    "BoundSaturationWarning",
    "DefaultStartWarning",
    "EngineError",
    "EnsembleSizeWarning",
    "SamplingFailureWarning",
]


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


class BoundSaturationWarning(UserWarning):
    """An optimiser converged with a coordinate at a bound of its support.

    Raised once per :func:`~ampere.inference.optimise` call by the ``"scipy"``
    route, naming every coordinate :func:`~ampere.inference.saturated_bounds`
    finds: the line search ran to the bound, so the point may be a corner of
    the box rather than a mode. The :class:`~ampere.results.Optimum` is still
    returned — an absent component's abundance at its floor is a correct
    answer — but an ensemble started around it collapses onto the bound.

    A distinct class, as :class:`SamplingFailureWarning` is, so that a caller
    can silence it, promote it to an error, or assert on it in a test without
    catching every other :class:`UserWarning` the stack emits.
    """


class DefaultStartWarning(UserWarning):
    """An ensemble engine's default start fell back to prior draws (W7.12).

    The default start of :class:`~ampere.inference.EmceeEngine` and
    :class:`~ampere.inference.ZeusEngine` is a ball at the optimiser's mode.
    Two problems cannot have one, and the run starts from the prior instead,
    saying so once: the optimiser found no start it could score, or its mode
    sits on a bound of the prior's support (a ball there is degenerate). The
    message names the reason, the bound and the remedies; ``initial="prior"``
    asks for the prior start outright and silences it. The run records the
    reason in ``ampere_start_fallback``.
    """


class EnsembleSizeWarning(UserWarning):
    """An ensemble engine was given a problem large enough to mix slowly (W7.12).

    Issued once per run at :data:`~ampere.inference.engine.ENSEMBLE_SIZE_WARNING`
    free coordinates or more, naming NUTS on a torch or jax problem as the
    tool for that size. A distinct class so that a caller who has chosen the
    ensemble deliberately can silence exactly this.
    """
