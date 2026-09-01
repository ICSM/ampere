"""Backend-neutral core contracts for ampere v2.

Everything in this package is pure numpy/scipy/astropy.units/stdlib. It never
imports torch, jax, paramax or any other optional dependency, lazily or
otherwise (``architecture.md`` §3-4), because it is the shared vocabulary the
reference, torch and jax backends all implement.

Currently landed: the parameter and prior contract (W1.3,
``docs/design/contracts/parameters.md``) and the ModelResult schema (W1.4,
``docs/design/contracts/results_schema.md``). The remaining §4 contracts —
``transform``, ``likelihood``, ``dataset``, ``astropy_compat`` — arrive with
W1.5-W1.7.
"""

from __future__ import annotations

from .exceptions import (
    AmpereError,
    ChannelError,
    ContractError,
    OptionalDependencyError,
    ParameterError,
    SchemaError,
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
from .results_schema import (
    DEFAULT_CHANNEL,
    REGULARITY_RTOL,
    Axis,
    AxisSpec,
    Cube,
    FunctionSamples,
    Image,
    Layout,
    ModelResult,
    Order,
    PhotometricPoints,
    Spectrum,
    TimeSeries,
    VisibilitySet,
)

__all__ = [
    "DEFAULT_CHANNEL",
    "REGULARITY_RTOL",
    "SEPARATOR",
    "AmpereError",
    "Axis",
    "AxisSpec",
    "Bijection",
    "Binding",
    "Buffer",
    "BufferSet",
    "ChannelError",
    "ContractError",
    "Cube",
    "FunctionSamples",
    "HierarchicalPrior",
    "Identity",
    "Image",
    "Layout",
    "Log",
    "Logit",
    "ModelResult",
    "OptionalDependencyError",
    "Order",
    "Parameter",
    "ParameterError",
    "ParameterMapping",
    "ParameterSet",
    "Parameterised",
    "PhotometricPoints",
    "Plate",
    "Prior",
    "PriorSpec",
    "SchemaError",
    "Spectrum",
    "Tie",
    "TimeSeries",
    "TyingError",
    "VisibilitySet",
    "default_bijection_for",
    "describe_prior",
    "log_density",
    "prior_from_spec",
]
