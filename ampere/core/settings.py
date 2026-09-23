"""Process-wide settings for ampere v2 (W5.30).

The package has no general configuration system, and this module is
deliberately not the start of one. It holds exactly one dataclass with one
field, read at validation time and restored by a context manager — nothing
here anticipates a second setting. A future work item that needs another
knob a user can flip and a test can restore adds a field to
:class:`Settings`; it does not build a registry.

The one setting so far controls what
:class:`ampere.core.parameter.Population` does when a ``layout="flat"``
declaration exceeds :data:`ampere.core.parameter.MAX_FLAT_MEMBERS`: refuse
(the default) or warn loudly and proceed. The limit itself — 128 — is
unaffected either way; this only changes what happens at the limit, for a
user who has read the cost and accepts it (Peter's ruling on W5.12's
"For Peter (2)", 2026-09-22).
"""

from __future__ import annotations

import contextlib
import dataclasses
from collections.abc import Iterator
from typing import Literal

__all__ = ["AmpereFlatPopulationWarning", "Settings", "override", "settings"]


class AmpereFlatPopulationWarning(UserWarning):
    """A ``layout="flat"`` :class:`~ampere.core.parameter.Population` was
    allowed to exceed :data:`~ampere.core.parameter.MAX_FLAT_MEMBERS`.

    Raised only when :attr:`Settings.flat_population_cap` is ``"warn"``; the
    default (``"refuse"``) raises a
    :class:`~ampere.core.exceptions.ParameterError` instead, and never gets
    here. Named for the package (:class:`~ampere.core.exceptions.AmpereError`
    is the exception equivalent), because this is not something a caller
    tunes into their own warning hierarchy — it fires only when a knob
    documented as costly has deliberately been turned.
    """


@dataclasses.dataclass
class Settings:
    """A small, explicit set of process-wide runtime settings.

    Instances are read at validation time by whatever they govern, never at
    import — the same value can be changed mid-session (directly, or scoped
    with :func:`override`) and the next validation sees it.

    Parameters
    ----------
    flat_population_cap
        ``"refuse"`` (the default): a ``layout="flat"``
        :class:`~ampere.core.parameter.Population` above
        :data:`~ampere.core.parameter.MAX_FLAT_MEMBERS` raises a
        :class:`~ampere.core.exceptions.ParameterError` at construction.
        ``"warn"``: the same population instead emits
        :class:`AmpereFlatPopulationWarning`, naming the member count, the
        measured cost and the ``layout="plate"`` remedy, and construction
        proceeds — the population merges and its ``lnprior`` evaluates as
        normal, just slowly (``hierarchical_population.md`` §5: about 60 ms
        per evaluation at the cap, so a hundred-thousand-evaluation ensemble
        run spends well over an hour in the prior alone).
    """

    flat_population_cap: Literal["refuse", "warn"] = "refuse"


#: The one process-wide instance. Read it directly; for a scoped change that
#: restores itself (the pattern every test here uses), use :func:`override`.
settings = Settings()


@contextlib.contextmanager
def override(**fields: object) -> Iterator[Settings]:
    """Temporarily set fields on :data:`settings`, restoring them on exit.

    An unknown field name fails atomically: :func:`getattr` reads every
    current value *before* anything is mutated, so a typo raises
    ``AttributeError`` with :data:`settings` left exactly as it was, rather
    than half-applying the change.

    >>> from ampere.core.settings import settings, override
    >>> with override(flat_population_cap="warn"):
    ...     settings.flat_population_cap
    'warn'
    >>> settings.flat_population_cap
    'refuse'
    """
    previous = {name: getattr(settings, name) for name in fields}
    for name, value in fields.items():
        setattr(settings, name, value)
    try:
        yield settings
    finally:
        for name, value in previous.items():
            setattr(settings, name, value)
