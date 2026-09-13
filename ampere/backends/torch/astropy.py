"""The opt-in astropy→native translation hook: ``contracts/astropy_compat.md`` §5.

The other half of the 2026-09-01 ruling in ``DEVELOPMENT_PLAN.md`` §2. The
default route, :func:`ampere.core.from_astropy`, evaluates the user's **actual**
``astropy.modeling`` model on the numpy path, and therefore has no gradient.
This module is how a caller asks for a *native* equivalent instead — explicitly,
by backend, in one call — and it raises for anything it has no curated row for
rather than falling back to the black box.

That refusal is the whole point. Silently returning the adapter would hand
somebody who asked for a differentiable model one that is not, and they would
discover it when ``NUTSEngine`` refused, several composition steps away from the
call that caused it. Silently substituting a lookalike for a model that *is* in
the table, without being asked, is the mirror-image failure and is what the
ruling forbids outright: a curated ``BlackBody`` is not guaranteed numerically
identical to astropy's, and a model whose implementation changed without its
author knowing is exactly the silent downgrade this architecture exists to
prevent.

**W4.6 lands the hook and an empty table; W4.7 lands the curated rows**
(``BlackBody``, ``PowerLaw1D``, ``BrokenPowerLaw1D``, ``Polynomial1D``,
``Gaussian1D``, ``Const1D`` and their compound sums and products). So today
every call here refuses, and :func:`translation_refusal` names what it could
not translate and says what the table currently holds. The jax twin of this
module is the same file with one word changed, deliberately: one name per
backend, everywhere (W2.12).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from ampere.core import (
    DEFAULT_CHANNEL,
    FunctionSamples,
    Model,
    astropy_components,
    translation_refusal,
)

from ._config import BACKEND

__all__ = ["TRANSLATIONS", "from_astropy"]

#: The curated table: astropy model class → native constructor, consulted leaf
#: by leaf so that a compound model is translatable exactly when every one of
#: its components is.
#:
#: **Empty in W4.6.** W4.7 fills it. It is a plain module-level dict rather than
#: a registry with an ``override=`` guard because it is small, closed and
#: read-only at runtime: a user who wants a native model ampere does not curate
#: writes it as an :class:`ampere.core.Model` directly, which is the honest
#: declaration and the one the capability flags then tell the truth about.
TRANSLATIONS: dict[type, Any] = {}


def from_astropy(
    model: Any,
    *,
    kind: type[FunctionSamples] | None = None,
    priors: Mapping[str, Any] | None = None,
    channel: str = DEFAULT_CHANNEL,
    grid: Any = None,
    output_unit: Any = None,
    equivalencies: Sequence[Any] = (),
) -> Model:
    """A **native** equivalent of an ``astropy.modeling`` model, or a refusal.

    The signature mirrors :func:`ampere.core.from_astropy` exactly, so that
    changing the import is the whole of the change a caller makes; what differs
    is the promise. This one returns a model built from this backend's own
    pieces — differentiable, and usable by ``NUTSEngine`` and ``VIEngine`` —
    or raises. It never returns the black-box adapter.

    Raises
    ------
    CapabilityError
        If *model*, or any component of a compound *model*, is not in
        :data:`TRANSLATIONS`. Which is every model, in W4.6. The message names
        the untranslatable components, says why there is no fallback, and
        points at :func:`ampere.core.from_astropy` for the black-box route.
    """
    components = astropy_components(model)
    if any(type(part) not in TRANSLATIONS for part in components):
        raise translation_refusal(BACKEND, model, TRANSLATIONS)
    return _compose(
        model,
        components,
        kind=kind,
        priors=priors,
        channel=channel,
        grid=grid,
        output_unit=output_unit,
        equivalencies=equivalencies,
    )


def _compose(
    model: Any,
    components: tuple[Any, ...],
    **options: Any,
) -> Model:
    """Build the native model from the curated rows. **W4.7 writes this.**

    Unreachable while :data:`TRANSLATIONS` is empty — every astropy model has
    at least one component, so :func:`from_astropy` refuses before reaching
    here — and left as an explicit refusal rather than a silent fall-through so
    that a row added without its composition rule fails loudly instead of
    returning something half-built.
    """
    raise NotImplementedError(
        f"ampere.backends.{BACKEND}.from_astropy() has table rows for "
        f"{sorted({type(part).__name__ for part in components})} but no rule for composing them "
        f"into one native model. W4.7 lands the rows and this composition together; a row "
        f"without one is a half-landed table, not a working translation."
    )
