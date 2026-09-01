"""Exception types shared by every ``ampere.core`` contract.

This module exists so that the new namespaces raise *specific* errors that
tests, tooling and users can catch by type, rather than bare ``ValueError``\\ s
carrying ad-hoc messages (which is what legacy does, per-module).

Two design rules govern everything here:

1. **Fail loudly and specifically.** ``architecture.md`` §1 is explicit that a
   request ampere cannot satisfy must produce a message naming the thing that
   is missing and what to do about it — never a silent downgrade and never an
   opaque failure deep inside somebody else's library.
2. **Stay catchable by the obvious builtin.** Every contract error also
   subclasses the builtin a caller would naturally reach for
   (:class:`ValueError` for malformed declarations, :class:`ImportError` for a
   missing optional dependency), so third-party code that pre-dates ampere's
   own hierarchy keeps working.

``OptionalDependencyError``'s shape is pinned here, as ``architecture.md`` §9
asks, so that it is not reinvented once per contract spec.
"""

from __future__ import annotations

__all__ = [
    "AmpereError",
    "ContractError",
    "OptionalDependencyError",
    "ParameterError",
    "TyingError",
]


class AmpereError(Exception):
    """Base class for every exception ampere raises deliberately."""


class ContractError(AmpereError, ValueError):
    """A core contract was used in a way it does not permit.

    Subclassed per contract (see :class:`ParameterError`) so that a caller may
    catch either the specific error or the whole family. Also a
    :class:`ValueError`, because a malformed declaration *is* a bad value and
    callers should not have to know ampere's hierarchy to handle one.
    """


class ParameterError(ContractError):
    """A parameter, buffer or prior declaration is malformed or unusable.

    Raised by :mod:`ampere.core.parameter` for: names that are not usable as
    keyword arguments, duplicate names, a parameter that is both fixed and
    prior-equipped (or neither), a prior that cannot be described neutrally,
    shape/unit mismatches, and cyclic hierarchical prior references.
    """


class TyingError(ParameterError):
    """Two or more tied parameters cannot be collapsed into one.

    Raised when tied sites disagree about something that must agree — shape,
    unit, prior, or fixed value — or when a tie names a site that does not
    exist. Kept distinct from :class:`ParameterError` because tying failures
    are the ones a user is most likely to want to handle (or explain)
    separately when composing a joint fit from independently written models.
    """


class OptionalDependencyError(AmpereError, ImportError):
    """An optional dependency is required for the operation being attempted.

    Raised **on use, never on import** (``architecture.md`` §4, rule 3), so
    that ``import ampere`` never requires torch, jax, or any other heavy
    dependency.

    Parameters
    ----------
    package
        Import name of the missing package, e.g. ``"torch"``.
    extra
        The ampere extra that provides it, e.g. ``"torch"`` for
        ``pip install ampere[torch]``. ``None`` if the package is not
        available through any extra.
    context
        What was being attempted, phrased as a noun phrase, e.g.
        ``"lowering a ParameterSet to a paramax pytree"``. Used verbatim at
        the start of the message.

    Attributes
    ----------
    package, extra, context
        As above; kept as attributes so tests and tooling can assert on the
        cause without parsing the message.

    Examples
    --------
    >>> err = OptionalDependencyError(
    ...     "paramax", extra="jax", context="lowering a ParameterSet"
    ... )
    >>> print(err)
    lowering a ParameterSet requires the optional dependency 'paramax', which is not installed.
    Install it with: pip install "ampere[jax]"
    >>> err.package, err.extra
    ('paramax', 'jax')
    """

    def __init__(self, package: str, extra: str | None = None, context: str | None = None) -> None:
        self.package = package
        self.extra = extra
        self.context = context
        subject = context if context else f"this operation ({package})"
        message = f"{subject} requires the optional dependency {package!r}, which is not installed."
        if extra is not None:
            message += f'\nInstall it with: pip install "ampere[{extra}]"'
        else:
            message += f"\nInstall it with: pip install {package}"
        super().__init__(message)
