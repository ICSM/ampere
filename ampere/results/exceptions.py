"""The results contract's error type — re-exported from its ruled home.

:class:`ResultsError` was declared here while ``ampere/core/exceptions.py`` was
a merged Phase-1 contract that W1.8 could not widen (``AGENTS.md`` ground
rule 9). Peter's 2026-09-03 ruling on ``docs/design/contracts/results.md`` §15
R5 moved it to ``ampere.core.exceptions``, beside every sibling contract error
(``ParameterError``, ``SchemaError``, ``TransformationError``,
``LikelihoodError``, ``DatasetError``). This module remains so that existing
import sites keep working; both paths name the same class.

Examples
--------
>>> import ampere.core.exceptions
>>> ResultsError is ampere.core.exceptions.ResultsError
True
"""

from __future__ import annotations

from ampere.core.exceptions import ResultsError

__all__ = ["ResultsError"]
