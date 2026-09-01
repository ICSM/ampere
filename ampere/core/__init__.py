"""Backend-neutral core contracts for ampere v2.

Everything in this package is pure numpy/scipy/astropy.units/stdlib. It never
imports torch, jax, paramax or any other optional dependency, lazily or
otherwise (``architecture.md`` §3-4), because it is the shared vocabulary the
reference, torch and jax backends all implement.

Currently landed: the parameter and prior contract (W1.3,
``docs/design/contracts/parameters.md``). The remaining §4 contracts —
``results_schema``, ``transform``, ``likelihood``, ``dataset``,
``astropy_compat`` — arrive with W1.4-W1.7.
"""

from __future__ import annotations

from .exceptions import (
    AmpereError,
    ContractError,
    OptionalDependencyError,
    ParameterError,
    TyingError,
)
from .parameter import (
    SEPARATOR,
    Bijection,
    Binding,
    Buffer,
    BufferSet,
    HierarchicalPrior,
    Identity,
    Log,
    Logit,
    Parameter,
    Parameterised,
    ParameterMapping,
    ParameterSet,
    Plate,
    Prior,
    PriorSpec,
    Tie,
    default_bijection_for,
    describe_prior,
    log_density,
    prior_from_spec,
)

__all__ = [
    "SEPARATOR",
    "AmpereError",
    "Bijection",
    "Binding",
    "Buffer",
    "BufferSet",
    "ContractError",
    "HierarchicalPrior",
    "Identity",
    "Log",
    "Logit",
    "OptionalDependencyError",
    "Parameter",
    "ParameterError",
    "ParameterMapping",
    "ParameterSet",
    "Parameterised",
    "Plate",
    "Prior",
    "PriorSpec",
    "Tie",
    "TyingError",
    "default_bijection_for",
    "describe_prior",
    "log_density",
    "prior_from_spec",
]
