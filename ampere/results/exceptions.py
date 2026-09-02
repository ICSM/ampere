"""The results contract's own error type.

Every other contract raises its own :class:`~ampere.core.exceptions.ContractError`
subclass — ``ParameterError``, ``SchemaError``, ``TransformationError``,
``LikelihoodError``, ``DatasetError`` — and they all live together in
``ampere/core/exceptions.py``. :class:`ResultsError` belongs beside them; it is
declared here instead only because ``ampere/core/exceptions.py`` is a merged
Phase-1 contract and widening it is a decision-log matter rather than a W1.8 one
(``AGENTS.md`` ground rule 9). ``docs/design/contracts/results.md`` §13 asks
W1.13 to move it, which is a two-line change and no import-site change at all,
since ``ampere.results`` re-exports it either way.
"""

from __future__ import annotations

from ampere.core.exceptions import ContractError

__all__ = ["ResultsError"]


class ResultsError(ContractError):
    """Results emission or serialisation was asked for something it cannot do.

    A ``ContractError``, hence a ``ValueError``: like every sibling, it signals a
    composition that cannot be honoured — a container carrying metadata that will
    not serialise, a draw array whose shape disagrees with the problem's free
    size, an unknown container kind on the way back in — rather than a bug in
    ampere.

    Examples
    --------
    >>> raise ResultsError("nope")
    Traceback (most recent call last):
        ...
    ampere.results.exceptions.ResultsError: nope
    """
