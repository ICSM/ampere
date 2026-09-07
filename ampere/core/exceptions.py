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
    "CapabilityError",
    "ChannelError",
    "CompositionError",
    "ContractError",
    "DatasetError",
    "LikelihoodError",
    "LoweringError",
    "LoweringFallbackWarning",
    "OptionalDependencyError",
    "ParameterError",
    "ResultsError",
    "SchemaError",
    "TransformationError",
    "TyingError",
]


class AmpereError(Exception):
    """Base class for every exception ampere raises deliberately."""


class CapabilityError(AmpereError, NotImplementedError):
    """A well-formed declaration asks for a capability this path does not have.

    The deliberate contrast is with :class:`ContractError`: nothing here is
    malformed. The declaration is valid, most routes through ampere consume
    it happily, and one specific path cannot — so the refusal must be a
    distinct type a caller can tell apart from "you wrote it wrongly", and
    the message must name the capability, not the declaration.

    The founding case (ruled 2026-09-03, ``lowering.md`` §12.1): a
    **discrete prior family**. Declaration, prior sampling, constrained-space
    ``log_prob``, ``prior_transform`` (scipy's discrete families implement
    ``ppf``) and lowering-as-distribution all work; what cannot exist is a
    continuous bijection to unconstrained space, so
    :func:`~ampere.core.parameter.default_bijection_for` raises this — and
    *only* it does, keeping the non-gradient routes that might eventually
    support discrete parameters ((variational) EM, numpyro-style enumeration,
    SBI, nested sampling, Bayesian optimisation) reachable. Discreteness is
    queryable from the canonical description
    (``describe_prior(prior).discrete``), so a future engine path branches on
    it rather than catching this.

    Also a :class:`NotImplementedError` — the builtin a caller naturally
    reaches for when an operation is unsupported rather than wrong — and
    deliberately **not** a :class:`ValueError`, because the value is fine.
    """


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


class SchemaError(ContractError):
    """A model result, channel or data container is malformed or unusable.

    Raised by :mod:`ampere.core.results_schema` for: channel names that are not
    usable, duplicate or missing channels, coordinate/value/uncertainty/mask
    shapes that disagree, coordinates that violate the ordering a container
    kind requires, units that are missing or of the wrong physical type, and
    complex values in a container kind that does not admit them.

    Kept distinct from :class:`ParameterError` because the two contracts have
    separate namespaces and separate failure modes: a channel name and a
    parameter name may collide harmlessly (see
    ``docs/design/contracts/results_schema.md`` §2).
    """


class ChannelError(SchemaError, KeyError):
    """A channel could not be bound, or was bound to the wrong kind of data.

    Raised when an instrument (or any other consumer) asks a
    :class:`~ampere.core.results_schema.ModelResult` for a channel that does
    not exist, or for one that exists but holds a different container kind than
    the consumer requires. ``DEVELOPMENT_PLAN.md`` §4.2 requires these
    mismatches to fail loudly at composition time rather than producing a
    confusing shape error inside a likelihood, and §4.3 makes name-based
    binding the instrument contract — so this is the error a consumer catches
    when it wants to offer an alternative rather than abort.

    Also a :class:`KeyError`, per this module's second rule: a
    ``ModelResult`` is a :class:`~collections.abc.Mapping`, and ``KeyError`` is
    the builtin a caller reaches for when a lookup fails. That inheritance is
    load-bearing rather than decorative — ``Mapping``'s own ``get`` and
    ``__contains__`` mixins are defined in terms of catching ``KeyError``, so
    without it ``result.get(name)`` and ``name in result`` would raise instead
    of answering.

    :meth:`__str__` is overridden because ``KeyError`` uniquely formats itself
    as ``repr(args[0])``, which would wrap every one of these messages in
    quotes and defeat the point of writing them.
    """

    def __str__(self) -> str:
        if len(self.args) == 1 and isinstance(self.args[0], str):
            return self.args[0]
        return super().__str__()


class TransformationError(ContractError):
    """A transformation or instrument chain is malformed, or misbehaved.

    Raised by :mod:`ampere.core.transform` for: a transformation handed a
    container kind it does not accept, one that returned something other than
    the kind it declares it produces, one that dropped a mask its input
    carried (``docs/design/contracts/results_schema.md`` §16 makes mask
    propagation an obligation on this contract), and malformed requirement
    declarations.

    Kept distinct from :class:`SchemaError` because a container can be
    perfectly well formed while the *chain* that produced it is not: the
    failure is in the composition, not in the data.
    """


class CompositionError(TransformationError):
    """Two pieces of a transformation chain cannot be composed.

    Raised at composition time — when an :class:`~ampere.core.transform.
    Instrument` is built, or when requirements are negotiated — for: a step
    whose input kind cannot be produced by the step before it, duplicate step
    labels within one chain, two instruments binding one channel with
    incompatible expectations, and requirements on one axis declared in
    incompatible units.

    ``DEVELOPMENT_PLAN.md`` §4.2-4.3 require these to fail loudly when the
    chain is *assembled*, not deep inside a likelihood evaluation thousands of
    samples later, which is why they are a distinct type: a caller composing a
    fit programmatically may reasonably catch this and try another chain.
    """


class LikelihoodError(ContractError):
    """A likelihood, noise model, kernel or censoring declaration is unusable.

    Raised by :mod:`ampere.core.likelihood` for: a family/noise-model
    combination whose marginalisation an engine cannot deliver, a GP solver
    strategy asked for something it does not support (an unquasiseparable
    kernel, a gridded layout, mixed coordinate units), a covariance matrix that
    will not factorise, kernel hyperparameters outside their support, censoring
    codes that are not :class:`~ampere.core.likelihood.LimitKind` values or that
    are misaligned with their container, and predicted/observed containers that
    do not describe the same samples.

    Kept distinct from :class:`SchemaError` because a container may be
    perfectly well formed and still be one this likelihood cannot consume —
    the failure is in the *composition*, not in either object.
    """


class DatasetError(ContractError):
    """A dataset, dataset collection or fitting problem is malformed.

    Raised by :mod:`ampere.core.dataset` for: an observed container whose kind
    the instrument chain cannot produce, a dataset naming a model the problem
    does not hold, colliding component labels in the joint parameter space, a
    latent declaration whose size disagrees with the number of samples that
    actually survive masking, capability declarations that contradict each
    other, and a fitting problem asked to evaluate before it is composed.

    Kept distinct from :class:`TransformationError` and :class:`LikelihoodError`
    because the failure is in the *assembly of the problem*: every piece may be
    individually well formed and still not compose into something an engine can
    run.

    Note what this is **not** for. A model that crashes, or a covariance that
    will not factorise, during an evaluation is not a malformed problem — it is
    a proposal that cannot be scored, and ``DEVELOPMENT_PLAN.md`` §4.5 requires
    those to become ``-inf`` with a recorded reason
    (:class:`~ampere.core.dataset.Failure`) rather than an exception. This error
    is for the composition-time failures that should stop a run before it
    starts.
    """


class ResultsError(ContractError):
    """Results emission or serialisation was asked for something it cannot do.

    Raised by :mod:`ampere.results` for: a container carrying metadata that
    will not serialise, a draw array whose shape disagrees with the problem's
    free size, an unknown container kind on the way back in. Like every
    sibling, it signals a composition that cannot be honoured rather than a
    bug in ampere.

    Declared in ``ampere/results/exceptions.py`` until Peter's 2026-09-03
    ruling on ``results.md`` §15 R5 moved it here beside its siblings;
    ``ampere.results`` re-exports it, so both import paths name this class.

    Examples
    --------
    >>> raise ResultsError("nope")
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.ResultsError: nope
    """


class LoweringError(AmpereError):
    """A valid declaration cannot be lowered onto the requested backend.

    Landed at the freeze (ruled 2026-09-03, ``lowering.md`` §12.4) with the
    shape §3.4 of that document specifies. Deliberately **not** a
    :class:`ContractError`: a family torch does not implement is a capability
    gap in the *backend*, not a malformed declaration by the user, and
    conflating the two would make "your prior is invalid" and "this backend
    cannot express your valid prior" indistinguishable to a caller catching
    by type. One shared class across all three backends, so tests and
    tooling assert on it uniformly.

    Raised (from Phase 2 on) when ``lowering.md`` §3.4's rule 2 applies: no
    exact construction from the target's primitives exists, so lowering
    fails early with a name attached — never a silent approximation, and
    never an automatic fallback to the reference path.

    Parameters
    ----------
    family
        The neutral name of what would not lower — a prior family
        (``"truncnorm"``), a kernel family, or a bijection class name.
    backend
        The backend that lacks it, e.g. ``"torch"``.
    parameter
        The parameter (or merged site) whose declaration triggered the
        failure, when one is known.
    detail
        Anything the caller can add — the exact missing primitive, say.

    Attributes
    ----------
    family, parameter, backend, detail
        As above, kept as fields so tooling asserts on the cause without
        parsing the message.

    Examples
    --------
    >>> err = LoweringError("truncnorm", backend="torch", parameter="temperature")
    >>> err.family, err.parameter, err.backend
    ('truncnorm', 'temperature', 'torch')
    >>> str(err).startswith("prior family 'truncnorm' (parameter 'temperature')")
    True
    """

    def __init__(
        self,
        family: str,
        *,
        backend: str,
        parameter: str | None = None,
        detail: str | None = None,
    ) -> None:
        self.family = family
        self.parameter = parameter
        self.backend = backend
        self.detail = detail
        where = f" (parameter {parameter!r})" if parameter is not None else ""
        why = detail if detail else "no exact construction exists from its primitives"
        super().__init__(
            f"prior family {family!r}{where} cannot be lowered to the {backend!r} backend: "
            f"{why}. Options: change the prior to a family the backend implements, register "
            f"your own lowering for it (register_lowering, Phase 2), or run on a backend "
            f"that has it. Ampere never substitutes an approximation silently."
        )


class LoweringFallbackWarning(UserWarning):
    """A native run is computing part of itself on the reference (numpy) path.

    **Landed at W2.13** (ruled 2026-09-07, ``DEVELOPMENT_PLAN.md`` §2's
    "Realisation surface (W2.13)" row, fold-in 8), replacing the two
    backend-local classes W2.4 and W2.5 each invented independently
    (``ampere.backends.torch.lowering.IcdfFallbackWarning`` and
    ``ampere.backends.jax.parameters.LoweringFallbackWarning``). One class, in
    the core, for the reason the backends' own docstrings gave when they asked
    for it: the fallback is a property of the *lowering contract*, not of a
    library, and a test that wants to assert "no native run silently computed
    in numpy" must be able to ``pytest.warns`` on one type across every
    backend.

    Warned when ``lowering.md`` §3.6's sanctioned fallback is taken — a prior
    family with no usable native ``icdf``, whose ``prior_transform`` therefore
    goes through ``scipy``'s ``ppf``. That is **not** the silent substitution
    §3.4 forbids: ``prior_transform`` has one mathematical definition and the
    reference path computes it exactly, so nothing about the posterior
    changes. It is warned about because a nominally torch- or jax-backed run
    that computes part of itself in numpy should never be a surprise
    discovered later.

    Warned **once per lowering** rather than once per call: the decision is
    made once, when the problem is lowered, before any sampling.
    ``FittingProblem(strict=True)`` turns it into a :class:`LoweringError`
    instead — one flag with one meaning, the run-level "I would rather fail
    than have anything smoothed over".

    A warning rather than an exception because it is not an error: this is a
    :class:`UserWarning` subclass and not an :class:`AmpereError`, and the two
    hierarchies are deliberately separate.
    """


class OptionalDependencyError(AmpereError, ImportError):
    """An optional dependency is required for the operation being attempted.

    Ratified in place at the freeze (ruled 2026-09-03 with ``lowering.md``
    §12.4's "the relevant exceptions" disposition; ``parameters.md`` §14 Q5
    asked for the ratification): the shape below is the contract.

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
