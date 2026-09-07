"""Datasets, dataset collections and the engine-facing surface (``§4.5``).

This is the contract in which ampere says *what an inference engine talks to*.
Everything below it — parameters, containers, instruments, likelihoods — is
vocabulary; this module is the sentence. ``DEVELOPMENT_PLAN.md`` §4.5 fixes the
surface, verbatim:

    A fitting problem (model + instruments + datasets) exposes:
    ``log_prob(params) -> float`` and ``log_likelihood`` / ``log_prior`` split;
    ``prior_transform(u)`` for nested sampling; ``simulate(params) -> data`` for
    SBI; capability flags: ``differentiable``, ``batchable``, ``device``,
    ``backend``.
    ... Failure signalling ... RNG policy.

Three objects, three jobs
-------------------------
:class:`Dataset`
    One observed container, the :class:`~ampere.core.transform.Instrument` that
    predicts it, and the :class:`~ampere.core.likelihood.Likelihood` that scores
    it. It owns no model: several datasets routinely observe one model through
    different instruments, so the model belongs a level up.
:class:`DatasetCollection`
    An ordered, labelled collection of datasets, plus the ``shared`` parameter
    set that is this contract's hyperprior extension point. It knows the joint
    log-likelihood is a sum and that the labels namespace the joint parameter
    space; it knows nothing about models or engines.
:class:`FittingProblem`
    Model(s) + collection + ties. It performs the **one** merge that produces
    the joint parameter space, runs the negotiate/compile/check lifecycle once
    at composition time, and exposes §4.5's surface.

The merge topology, and why it is nested
----------------------------------------
``parameters.md`` §12.4 records that :meth:`~ampere.core.parameter.ParameterSet.
merge` is not associative and that components must be merged "in one call"; its
§14 keeps a **nested** ``ParameterMapping`` expressly open as this contract's
alternative. This contract takes the nested option, under one rule:

    **The nesting rule.** Every composite object performs exactly one
    :meth:`~ampere.core.parameter.ParameterSet.merge` over its *immediate*
    children and **retains** the resulting
    :class:`~ampere.core.parameter.ParameterMapping`. Values flow down by
    calling :meth:`~ampere.core.parameter.ParameterMapping.distribute` at each
    level. No mapping is ever discarded, so no binding is ever lost.

The associativity hazard §12.4 names is *dropping* an inner merge's bindings by
re-merging its ``merged`` set and throwing the mapping away. Re-distributing
instead cannot lose a binding: it only re-derives one, one level down. Three
levels exist, and no more:

1. :class:`~ampere.core.transform.Instrument`, over its steps' labels, giving
   ``calibrate.scale``;
2. :class:`Dataset`, over :data:`INSTRUMENT_COMPONENT`,
   :data:`LIKELIHOOD_COMPONENT` and :data:`LATENT_COMPONENT`, giving
   ``instrument.calibrate.scale``;
3. :class:`FittingProblem`, over the model labels, the dataset labels and
   :data:`SHARED_COMPONENT`, giving ``sed.instrument.calibrate.scale``.

The instrument level already exists — ``transformations.md`` §14 built it — so
the choice was never "nest or not" but "how many levels, and who re-distributes".
The rationale, the alternative (a single flat merge over every leaf) and what a
reversal would cost are written out in ``docs/design/contracts/inference.md``
§4; this docstring records the rule, not the argument.

The lifecycle, stated once
--------------------------
``transformations.md`` §14 asks this contract to own *when*
:func:`~ampere.core.transform.negotiate` and
:meth:`~ampere.core.transform.Model.compile_for` are called, and
``likelihoods.md`` §16 asks it to call
:meth:`~ampere.core.likelihood.Likelihood.check_alignment` and
:meth:`~ampere.core.likelihood.Likelihood.check_engine`. All of it happens once,
in :meth:`FittingProblem.__init__`, in this order:

1. **Bind** each dataset to its model.
2. **Negotiate** — :func:`~ampere.core.transform.negotiate` over the instruments
   reading each model, giving that model's per-channel requirements.
3. **Compile** — ``model.compile_for(requirements)``; the returned model is the
   one every later evaluation uses.
4. **Merge** — one :meth:`~ampere.core.parameter.ParameterSet.merge` over the
   compiled models, the datasets and ``shared``, with the ties.
5. **Validate** — evaluate the compiled model once at a reference θ, push the
   result through every instrument and call ``check_alignment`` on each dataset,
   which is where a mis-shaped chain, a unit mismatch and an unimplemented
   latent combination are refused.

Nothing in the hot loop re-negotiates, re-compiles or re-merges at the problem
level. (An :class:`~ampere.core.transform.Instrument` still re-merges its own
steps per call, by W1.5's deliberate choice; see ``inference.md`` §11.)

Everything here is backend-neutral: numpy, ``astropy.units`` and stdlib only
(``architecture.md`` §3-4). The narrative spec, whose every example runs as a
doctest, is ``docs/design/contracts/inference.md``.

The simple path, end to end:

>>> import numpy as np, astropy.units as u, scipy.stats as st
>>> from ampere.core import Model, Parameter, Spectrum
>>> class Line(Model):
...     def __init__(self, wavelength):
...         self.register_buffer("wavelength", wavelength, unit=u.um)
...         self.register_parameter(Parameter("slope", st.norm(1.0, 1.0)))
...     def evaluate(self, **values):
...         ctx = self.context(values)
...         return Spectrum(ctx["wavelength"] * u.um, ctx["slope"] * ctx["wavelength"] * u.Jy)
>>> observed = Spectrum(
...     [1.0, 2.0, 3.0] * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.1, 0.1, 0.1] * u.Jy
... )
>>> problem = FittingProblem(Line(np.array([1.0, 2.0, 3.0])), [Dataset(observed)])
>>> problem.parameters.free_names
('model.slope',)
>>> round(problem.log_likelihood({"model.slope": 2.0}), 6)
4.15094
>>> problem.capabilities
Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
"""

from __future__ import annotations

import collections
import dataclasses
import enum
import math
import types
import warnings
from collections.abc import Iterable, Iterator, Mapping, Sequence
from typing import Any, Protocol, runtime_checkable

import numpy as np

from .exceptions import (
    CompositionError,
    DatasetError,
    LikelihoodError,
    TransformationError,
)
from .likelihood import (
    DTYPE,
    GaussianFamily,
    IndependentNoise,
    LatentDeclaration,
    Likelihood,
    Marginalisation,
)
from .parameter import (
    SEPARATOR,
    ParameterMapping,
    ParameterSet,
    Tie,
    Value,
)
from .results_schema import FunctionSamples, ModelResult
from .rng import SEED_BYTES
from .rng import generator as _generator
from .transform import Instrument, Model, negotiate

__all__ = [
    "DEFAULT_FAILURE_HISTORY",
    "INSTRUMENT_COMPONENT",
    "LATENT_COMPONENT",
    "LATENT_NAME",
    "LIKELIHOOD_COMPONENT",
    "MODEL_COMPONENT",
    "SHARED_COMPONENT",
    "Capabilities",
    "Capable",
    "Dataset",
    "DatasetCollection",
    "Evaluation",
    "Failure",
    "FailureReason",
    "FittingProblem",
    "Simulation",
    "declared_capabilities",
]

ArrayLike = Any

#: Component label under which a :class:`Dataset` merges its instrument's
#: parameters. A *role* name, not the instrument's own label: the instrument's
#: label already defaults to its channel name, which is also the commonest
#: dataset label, so using it here would produce ``sed.sed.calibrate.scale``.
INSTRUMENT_COMPONENT = "instrument"

#: Component label for the dataset's likelihood parameters (the family's and the
#: noise model's, which :class:`~ampere.core.likelihood.Likelihood` already
#: holds in one flat set).
LIKELIHOOD_COMPONENT = "likelihood"

#: Component label for a latent-GP declaration, present only when the
#: family/noise combination actually needs one for *these* data.
LATENT_COMPONENT = "latent"

#: Local name of the whitened latent block within :data:`LATENT_COMPONENT`, so
#: the merged name reads ``spectrum.latent.z``.
LATENT_NAME = "z"

#: Default component label for a :class:`FittingProblem`'s single model.
MODEL_COMPONENT = "model"

#: Default component label for a :class:`DatasetCollection`'s shared parameters.
SHARED_COMPONENT = "shared"

#: How many recent :class:`Failure` records a :class:`FittingProblem` keeps.
#: Bounded, because a long run can propose millions of unscoreable points and a
#: diagnostic must not become a memory leak; the *counts*
#: (:attr:`FittingProblem.failure_counts`) are unbounded and are what a driver
#: should report.
DEFAULT_FAILURE_HISTORY = 64


def _mask_of(container: FunctionSamples) -> np.ndarray | None:
    """A container's mask as a flat boolean array, or ``None`` if it has none."""
    if container.mask is None:
        return None
    return np.asarray(container.mask, dtype=bool).ravel()


def _masks_equal(left: np.ndarray | None, right: np.ndarray | None) -> bool:
    """Compare two masks, treating ``None`` and an all-False mask as the same."""
    if left is None and right is None:
        return True
    if left is None:
        return not bool(np.any(right))
    if right is None:
        return not bool(np.any(left))
    return left.shape == right.shape and bool(np.array_equal(left, right))


def _empty_mapping() -> Mapping[str, Any]:
    """A fresh immutable empty mapping — a dataclass default 3.11 accepts."""
    return types.MappingProxyType({})


def _check_label(name: object, kind: str) -> str:
    """Labels become merge components, so they must be bare identifiers."""
    if not isinstance(name, str) or not name.isidentifier():
        raise DatasetError(
            f"{kind} {name!r} is not usable: it must be a valid Python identifier, because it "
            f"becomes a component label in the joint parameter space and merged names are "
            f"{SEPARATOR!r}-separated identifiers."
        )
    return name


# ---------------------------------------------------------------------------
# Capability flags
# ---------------------------------------------------------------------------


@runtime_checkable
class Capable(Protocol):
    """What a model or transformation declares about how it can be run.

    ``DEVELOPMENT_PLAN.md`` §4.5 lists ``differentiable``, ``batchable`` and
    ``device`` as the fitting problem's capability flags. They are properties of
    the *pieces*, not of the problem: a problem is differentiable exactly when
    everything a gradient would have to pass through is.

    **Promoted into W1.5's ABCs at the freeze** (ruled 2026-09-03,
    ``inference.md`` §19.6): :class:`~ampere.core.transform.Model` and
    :class:`~ampere.core.transform.Transformation` carry the three as class
    attributes with the conservative defaults ``False``, ``False`` and
    ``"cpu"`` — the reference path's honest answers, reproducing the earlier
    ``getattr`` semantics exactly — so every piece a problem composes now
    declares them, silence included, and
    :func:`declared_capabilities` reads the attributes directly. Phase 2's
    torch and jax backends override them on their own subclasses; this
    Protocol remains the statement of the surface for anything duck-typed
    into :attr:`Dataset.capability_parts`.

    **A fourth flag joined them at W2.12** (decided by Fable 2026-09-07,
    ``DEVELOPMENT_PLAN.md`` §4.5): ``BACKEND``, the name of the rung of the
    capability ladder this piece runs on. It behaves exactly like ``DEVICE``
    — every part must agree, and disagreement is a configuration mistake
    rather than something ampere silently repairs — and it is what makes
    ``ampere_backend`` in a run's provenance a fact about the problem rather
    than a declaration by whoever constructed the engine.
    """

    #: Whether a gradient can be taken through this object's evaluation.
    DIFFERENTIABLE: bool
    #: Whether it evaluates a batch of parameter vectors in one call.
    BATCHABLE: bool
    #: Device its arrays live on: ``"cpu"``, ``"cuda"``, ``"cuda:0"``, ...
    DEVICE: str
    #: Backend that supplies it: ``"reference"``, ``"torch"``, ``"jax"``, ...
    #: One name per backend — the same string keys ``lowering.md`` §12.8's
    #: registry and ids the conformance fixtures.
    BACKEND: str


@dataclasses.dataclass(frozen=True)
class Capabilities:
    """What an engine may assume about a :class:`FittingProblem`.

    ``DEVELOPMENT_PLAN.md`` §4.5's flags, in one immutable record so that
    an engine driver takes them as a unit and a run's provenance records them as
    one.

    Attributes
    ----------
    differentiable
        Whether ``log_prob`` admits a gradient, i.e. whether NUTS/HMC and
        gradient-based VI or optimisation can run. Consumed by
        :meth:`~ampere.core.likelihood.Likelihood.check_engine`, which refuses a
        gradient-free engine a latent-variable likelihood.
    batchable
        Whether ``log_prob`` accepts a stack of parameter vectors and returns a
        stack of values. Ampere's reference path does not; a backend that does
        unlocks vectorised ensemble samplers and batched SBI budgets.
    device
        Where the arrays live. ``"cpu"`` on the reference path, always.
    backend
        Which rung of the capability ladder supplies the pieces:
        ``"reference"``, and in Phase 2 ``"torch"`` or ``"jax"``. **One name
        per backend, everywhere** — the same string keys ``lowering.md``
        §12.8's registry and ids the conformance fixtures, so a run's
        ``ampere_backend`` and a registered lowering's ``backend`` are
        comparable without a translation table.

    Examples
    --------
    >>> Capabilities()
    Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
    """

    differentiable: bool = False
    batchable: bool = False
    device: str = "cpu"
    backend: str = "reference"

    def __post_init__(self) -> None:
        object.__setattr__(self, "differentiable", bool(self.differentiable))
        object.__setattr__(self, "batchable", bool(self.batchable))
        if not isinstance(self.device, str) or not self.device:
            raise DatasetError(f"a device must be a non-empty string, got {self.device!r}")
        if not isinstance(self.backend, str) or not self.backend:
            raise DatasetError(f"a backend must be a non-empty string, got {self.backend!r}")

    def to_dict(self) -> dict[str, Any]:
        """A plain-data form for a run's provenance attrs (W1.8)."""
        return dataclasses.asdict(self)


def declared_capabilities(parts: Sequence[object]) -> Capabilities:
    """The capabilities of a whole, from what its parts declare.

    Conjunctive by construction, and deliberately so: a chain is differentiable
    only if *every* link is, because one black-box step is enough to stop a
    gradient. An empty collection of parts is **not** differentiable — ``all([])``
    is ``True``, and silently promising gradients for a problem with nothing in
    it is precisely the silent-capability-upgrade this architecture forbids.

    The flags are read directly (ruled 2026-09-03, ``inference.md`` §19.6):
    :class:`~ampere.core.transform.Model` and
    :class:`~ampere.core.transform.Transformation` carry ``DIFFERENTIABLE``,
    ``BATCHABLE``, ``DEVICE`` and — since W2.12 — ``BACKEND`` as class
    attributes with the conservative defaults, so every part a problem
    composes declares them — a duck-typed part must too (:class:`Capable` is
    the surface).

    ``BACKEND`` aggregates by the **device rule** rather than the conjunctive
    one, because it is an identity and not a promise: there is no
    conservative answer to "half of this problem is torch and half is jax".

    Parameters
    ----------
    parts
        The models and transformations a gradient would have to pass through.

    Raises
    ------
    DatasetError
        If the parts declare more than one device, or more than one backend.
        Ampere will not choose either for you: moving arrays between devices
        silently is how a run becomes mysteriously slow, and converting them
        between array libraries silently is how a run loses its gradients.

    Examples
    --------
    >>> class Native:
    ...     DIFFERENTIABLE = True
    ...     BATCHABLE = True
    ...     DEVICE = "cuda"
    ...     BACKEND = "torch"
    >>> declared_capabilities([Native(), Native()])
    Capabilities(differentiable=True, batchable=True, device='cuda', backend='torch')

    One conservative part withdraws the whole conjunctive claim — the ABCs'
    defaults are ``False``/``False``/``"cpu"``/``"reference"``, so a subclass
    that stays silent inherits the reference answers rather than promising
    anything — and a CPU part beside a GPU one is a device disagreement rather
    than a quiet round trip:

    >>> class NativeOnCpu:
    ...     DIFFERENTIABLE = True
    ...     BATCHABLE = True
    ...     DEVICE = "cpu"
    ...     BACKEND = "reference"
    >>> class SilentOnCpu:
    ...     DIFFERENTIABLE = False
    ...     BATCHABLE = False
    ...     DEVICE = "cpu"
    ...     BACKEND = "reference"
    >>> declared_capabilities([NativeOnCpu(), SilentOnCpu()])
    Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
    >>> declared_capabilities([])
    Capabilities(differentiable=False, batchable=False, device='cpu', backend='reference')
    >>> declared_capabilities([Native(), SilentOnCpu()])
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.DatasetError: the pieces of this problem declare different devices...

    A backend disagreement is refused the same way, and for the sharper
    reason: ampere does not convert arrays between libraries on the user's
    behalf, so a mixed problem would otherwise fail two steps later inside a
    backend, or silently drop the gradients it was assembled to provide.

    >>> class TorchOnCpu:
    ...     DIFFERENTIABLE = True
    ...     BATCHABLE = True
    ...     DEVICE = "cpu"
    ...     BACKEND = "torch"
    >>> declared_capabilities([TorchOnCpu(), SilentOnCpu()])
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.DatasetError: the pieces of this problem declare different backends...
    """
    if not parts:
        return Capabilities()
    devices = {str(part.DEVICE) for part in parts}  # type: ignore[attr-defined]
    if len(devices) > 1:
        raise DatasetError(
            f"the pieces of this problem declare different devices {sorted(devices)}. Ampere does "
            f"not move arrays between devices on your behalf — that turns a configuration mistake "
            f"into a silent performance collapse. Put every model and transformation on one "
            f"device, or pass capabilities=Capabilities(device=...) to state which one is meant."
        )
    backends = {str(part.BACKEND) for part in parts}  # type: ignore[attr-defined]
    if len(backends) > 1:
        # Name the offending pieces by class, grouped by the backend each
        # declares. W2.13 widened the parts to include noise models and GP
        # solvers, and the commonest way to reach this message is now a
        # native problem left with the core (numpy) IndependentNoise or
        # DenseGP -- which the old wording, about "models and transformations",
        # did not help anyone find.
        culprits = "; ".join(
            f"{name}: "
            + ", ".join(
                sorted({type(part).__name__ for part in parts if str(part.BACKEND) == name})  # type: ignore[attr-defined]
            )
            for name in sorted(backends)
        )
        raise DatasetError(
            f"the pieces of this problem declare different backends {sorted(backends)} "
            f"({culprits}). Ampere does not convert arrays between libraries on your behalf — a "
            f"mixed problem is a configuration mistake that would otherwise fail two steps later "
            f"inside a backend, or silently drop gradients. Build every model, transformation, "
            f"noise model and GP solver on one backend, or pass "
            f"capabilities=Capabilities(backend=...) to state which one is meant."
        )
    return Capabilities(
        differentiable=all(bool(part.DIFFERENTIABLE) for part in parts),  # type: ignore[attr-defined]
        batchable=all(bool(part.BATCHABLE) for part in parts),  # type: ignore[attr-defined]
        device=devices.pop(),
        backend=backends.pop(),
    )


# ---------------------------------------------------------------------------
# Failure signalling
# ---------------------------------------------------------------------------


class FailureReason(enum.StrEnum):
    """Why a point could not be scored.

    ``DEVELOPMENT_PLAN.md`` §4.5 requires ``log_prob`` to return ``-inf`` *with a
    recorded reason*, and the reason is only useful if it is a fixed vocabulary:
    a free-text message cannot be counted, and "37 % of your proposals failed,
    all of them ``non_positive_definite``" is the sentence a user needs. A
    :class:`enum.StrEnum`, so the code serialises into provenance attrs as
    itself.
    """

    #: The model raised while computing its ``ModelResult``. The classic
    #: external-simulator crash.
    MODEL_FAILED = "model_failed"
    #: A transformation in the instrument chain raised.
    INSTRUMENT_FAILED = "instrument_failed"
    #: The predicted container held NaN or inf on a sample the mask retains —
    #: §4.5's "external simulators ... return NaNs", caught by the likelihood
    #: and classified here.
    NON_FINITE_PREDICTION = "non_finite_prediction"
    #: The likelihood raised: a covariance that will not factorise, a kernel
    #: hyperparameter outside its support, a non-positive rate. This is the
    #: conversion ``likelihoods.md`` §16 hands to this contract.
    LIKELIHOOD_FAILED = "likelihood_failed"
    #: The likelihood returned NaN. Distinct from ``-inf``, which is a perfectly
    #: good answer meaning "impossible", and from an exception.
    NON_FINITE_LOG_LIKELIHOOD = "non_finite_log_likelihood"


@dataclasses.dataclass(frozen=True)
class Failure:
    """A recorded reason a point could not be scored.

    Deliberately a plain-data record and not the exception: it is kept in a
    history, counted, and written into a run's provenance, all of which want
    something small, picklable and comparable. The exception's *type name* is
    kept because it is what distinguishes an ampere contract failure from a
    user simulator's own error class; the traceback is not, because the same
    failure will happen thousands of times.

    Examples
    --------
    >>> failure = Failure(FailureReason.MODEL_FAILED, "the RT code exited 1", where="sed")
    >>> failure.reason
    <FailureReason.MODEL_FAILED: 'model_failed'>
    >>> failure.to_dict()["reason"]
    'model_failed'
    """

    reason: FailureReason
    message: str
    where: str = ""
    exception_type: str = ""
    values: Mapping[str, float] = dataclasses.field(default_factory=_empty_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "values", types.MappingProxyType(dict(self.values)))

    def to_dict(self) -> dict[str, Any]:
        """A JSON-compatible form for a run's provenance attrs (W1.8)."""
        return {
            "reason": str(self.reason),
            "message": self.message,
            "where": self.where,
            "exception_type": self.exception_type,
            "values": dict(self.values),
        }

    def __str__(self) -> str:
        where = f" [{self.where}]" if self.where else ""
        at = ""
        if self.values:
            at = " at " + ", ".join(f"{k}={v:.6g}" for k, v in sorted(self.values.items()))
        return f"{self.reason}{where}: {self.message}{at}"


def _scalars(values: Mapping[str, Value] | None) -> dict[str, float]:
    """The scalar entries of a routed component's values.

    Recorded on a :class:`Failure` so that "which hyperparameters broke the
    Cholesky?" is answerable from the failure record itself. Arrays are
    deliberately excluded: a latent block of 10⁵ values kept on each of a
    bounded history's entries would be a memory leak wearing a diagnostic's
    clothes, and it is the scalars — an amplitude, a length-scale, a rate — that
    localise this class of failure.
    """
    if not values:
        return {}
    found: dict[str, float] = {}
    for name, value in values.items():
        array = np.asarray(value)
        if array.ndim == 0 and array.dtype.kind in "fiub":
            found[name] = float(array)
    return found


def _failure_from(
    reason: FailureReason,
    error: BaseException,
    where: str,
    values: Mapping[str, Value] | None = None,
) -> Failure:
    return Failure(
        reason=reason,
        message=str(error).replace("\n", " "),
        where=where,
        exception_type=type(error).__name__,
        values=_scalars(values),
    )


@dataclasses.dataclass(frozen=True)
class Evaluation:
    """Everything one call to :meth:`FittingProblem.evaluate` produced.

    ``log_prob`` alone is a float, and a float cannot carry the reason it is
    ``-inf``. This record is how §4.5's "with a recorded reason" is delivered
    without changing ``log_prob``'s signature, and it is also what W1.8 stores
    per posterior draw (``DEVELOPMENT_PLAN.md`` §4.6 asks every run to keep
    per-sample ``log_likelihood`` and ``log_prior``).

    Attributes
    ----------
    log_prior
        The joint log prior. ``-inf`` outside the support.
    log_likelihood
        The joint log likelihood, or **NaN** when it was not evaluated because
        the prior had already ruled the point out. NaN rather than ``-inf``
        because "not evaluated" and "impossible" are different statements, and
        an importance-reweighting consumer (design horizon (b)) must be able to
        tell them apart.
    log_prob
        ``log_prior + log_likelihood``, or ``-inf``.
    contributions
        Per-dataset log-likelihood, by dataset label — the natural input to
        "which dataset is driving this fit?". Empty when the point was rejected
        by the prior, and **partial** when a dataset failed: it holds the
        datasets scored before the failure, in collection order, which is what
        localises the failure alongside :attr:`failure`'s own ``where``.
    failure
        Why the point could not be scored, or ``None``. A point rejected by the
        prior has **no** failure: zero prior mass is an answer, not an error.
    """

    log_prior: float
    log_likelihood: float
    log_prob: float
    contributions: Mapping[str, float] = dataclasses.field(default_factory=_empty_mapping)
    failure: Failure | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "contributions", types.MappingProxyType(dict(self.contributions)))

    @property
    def failed(self) -> bool:
        """Whether a failure was recorded (not merely that ``log_prob`` is ``-inf``)."""
        return self.failure is not None

    def __repr__(self) -> str:
        bits = f"log_prob={self.log_prob:.6g}, log_prior={self.log_prior:.6g}"
        if self.failure is not None:
            return f"<Evaluation {bits}, failed: {self.failure}>"
        return f"<Evaluation {bits}, log_likelihood={self.log_likelihood:.6g}>"


@dataclasses.dataclass(frozen=True)
class Simulation:
    """One draw of the forward model: ``(θ, x)``, or a flagged failure.

    ``DEVELOPMENT_PLAN.md`` §4.5 asks for ``simulate(params) -> data`` and that
    "``simulate`` failures are flagged so SBI can reject-and-record rather than
    train on garbage". Flagging, not raising, is the point: a simulation budget
    of 10⁴ draws with a 2 % crash rate should produce 9 800 usable pairs and a
    count, not stop at the first crash.

    Attributes
    ----------
    parameters
        The θ this simulation used, as merged names. Present even for a failure,
        because a rejected θ is itself information about the prior.
    theta
        The same θ as the flat free-parameter vector, in
        :attr:`FittingProblem.parameters` order — what an SBI package wants.
    results
        The raw :class:`~ampere.core.results_schema.ModelResult` per model
        label. The ``(θ, ModelResult)`` pair design horizon (c) wants for
        emulator training sets.
    predicted
        The noise-free prediction per dataset label, after the instrument chain.
    observations
        Noisy draws per dataset label, or ``None`` when ``observe=False``.
    failure
        Why the simulation could not be completed, or ``None``.
    """

    parameters: Mapping[str, Value]
    theta: np.ndarray
    results: Mapping[str, ModelResult] = dataclasses.field(default_factory=_empty_mapping)
    predicted: Mapping[str, FunctionSamples] = dataclasses.field(default_factory=_empty_mapping)
    observations: Mapping[str, FunctionSamples] | None = None
    failure: Failure | None = None

    def __post_init__(self) -> None:
        set_ = object.__setattr__
        set_(self, "parameters", types.MappingProxyType(dict(self.parameters)))
        set_(self, "results", types.MappingProxyType(dict(self.results)))
        set_(self, "predicted", types.MappingProxyType(dict(self.predicted)))
        if self.observations is not None:
            set_(self, "observations", types.MappingProxyType(dict(self.observations)))
        # Copy before freezing: setflags on the array we were handed would make
        # the *caller's* array read-only as a side effect of constructing this.
        theta = np.array(self.theta, dtype=DTYPE, copy=True)
        theta.setflags(write=False)
        set_(self, "theta", theta)

    @property
    def failed(self) -> bool:
        """Whether this draw should be rejected rather than trained on."""
        return self.failure is not None

    def __repr__(self) -> str:
        if self.failure is not None:
            return f"<Simulation failed: {self.failure}>"
        drawn = "" if self.observations is None else ", observed"
        return f"<Simulation {len(self.predicted)} dataset(s){drawn}>"


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


class Dataset:
    """One observation, the instrument that predicts it, the likelihood that scores it.

    The unit ``DEVELOPMENT_PLAN.md`` §4.5 composes a fitting problem out of, and
    the smallest object that can answer "how well does this prediction match
    *these* data?". It deliberately holds **no model**: the commonest joint fit
    in ampere's target scope is one physical model observed by several
    instruments, so a model inside a dataset would either be duplicated or
    shared by object identity — and shared-by-identity is exactly
    ``prior_art.md`` Tension 1's gammapy footgun. A dataset names its model with
    a *string* instead (:attr:`model`), the same way its instrument names its
    channel.

    Parameters
    ----------
    observed
        The data. Any :class:`~ampere.core.results_schema.FunctionSamples`
        subclass — a ``Spectrum``, a set of ``PhotometricPoints``, a
        ``VisibilitySet``.
    instrument
        The chain mapping a model channel to a prediction of *these* data.
        Defaults to a pure channel binding on
        :data:`~ampere.core.results_schema.DEFAULT_CHANNEL`, kind-checked
        against ``observed`` — so a model that already produces the observable
        needs no instrument at all. **The caller builds the instrument from
        the observed container's own coordinates** where a step reproduces
        them (a Fourier step's (u, v) buffer, a resampler's target grid):
        ``check_alignment`` compares axes exactly, and coordinates recomputed
        from first principles differ in their last bits (W1.11 gap I-2's
        rule, ``transformations.md`` §10).
    likelihood
        Defaults to an i.i.d. Gaussian
        (:class:`~ampere.core.likelihood.GaussianFamily` with
        :class:`~ampere.core.likelihood.IndependentNoise`), which is what the
        container's own ``uncertainty`` means.
    label
        Component label for this dataset's parameters in a joint fit. Defaults
        to the instrument's label, which itself defaults to the channel name.
    model
        Which model of a multi-model problem this dataset observes. ``None``
        means "the only one", and a problem with several models refuses a
        dataset that does not say.
    latent_name
        Local name of the whitened latent block, when the likelihood needs one.
    meta
        Free-form metadata; carried, never interpreted.

    Raises
    ------
    DatasetError
        If the instrument chain cannot produce the observed container's kind —
        caught here, at composition, rather than as a shape error inside a
        likelihood thousands of samples later.

    Examples
    --------
    >>> import numpy as np, astropy.units as u, scipy.stats as st
    >>> from ampere.core import GaussianProcessNoise, Likelihood, Matern32, Spectrum
    >>> observed = Spectrum(
    ...     [1.0, 2.0, 3.0] * u.um, [2.0, 4.0, 6.0] * u.Jy, uncertainty=[0.1] * 3 * u.Jy
    ... )
    >>> plain = Dataset(observed, label="spectrum")
    >>> plain.parameters.names
    ()
    >>> flexible = Dataset(
    ...     observed,
    ...     likelihood=Likelihood(
    ...         GaussianFamily(),
    ...         GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
    ...     ),
    ...     label="spectrum",
    ... )
    >>> flexible.parameters.free_names
    ('likelihood.amplitude', 'likelihood.length_scale')
    """

    def __init__(
        self,
        observed: FunctionSamples,
        instrument: Instrument | None = None,
        likelihood: Likelihood | None = None,
        *,
        label: str | None = None,
        model: str | None = None,
        latent_name: str = LATENT_NAME,
        meta: Mapping[str, Any] | None = None,
    ) -> None:
        if not isinstance(observed, FunctionSamples):
            raise DatasetError(
                f"a Dataset is built from an observed container, got {type(observed).__name__}. "
                f"Observed data use the same kinds a model produces (Spectrum, "
                f"PhotometricPoints, VisibilitySet, ...) — results_schema.md §16."
            )
        self.observed = observed

        if instrument is None:
            # An implicit pure-binding instrument takes the dataset's own
            # label when one is given: the instrument label is how a user
            # identifies which instrument constrained what (ruled 2026-09-03,
            # transformations.md §15 Q4), and for an instrument the dataset
            # itself conjured, the dataset's label is that identity — so two
            # labelled datasets on one channel stay composable without the
            # user naming instruments nobody wrote.
            instrument = Instrument(
                input_kind=type(observed), label=None if label is None else str(label)
            )
        if not isinstance(instrument, Instrument):
            raise DatasetError(
                f"a Dataset's instrument must be an Instrument, got {type(instrument).__name__}. "
                f"A bare Transformation becomes one with Instrument([step])."
            )
        # A dataset is a composed object, so its instrument is frozen here
        # (ruled 2026-09-03, transformations.md §15 Q5): the per-log_prob
        # re-merge disappears from the hot loop, and a step reconfigured after
        # composition is refused at the next use rather than silently ignored.
        self.instrument = instrument.freeze()

        if likelihood is None:
            likelihood = Likelihood(GaussianFamily(), IndependentNoise())
        if not isinstance(likelihood, Likelihood):
            raise DatasetError(
                f"a Dataset's likelihood must be a Likelihood, got {type(likelihood).__name__}."
            )
        self.likelihood = likelihood

        self.label = _check_label(instrument.label if label is None else label, "dataset label")
        self.model = None if model is None else _check_label(model, "model reference")
        self.meta: Mapping[str, Any] = types.MappingProxyType(dict(meta) if meta else {})

        self._check_output_kind()
        self._latent = self._declare_latent(latent_name)
        self._mapping = ParameterSet.merge(self._components())
        # The effective mask, resolved once (Peter's W1.6 §17 Q2 ruling). All
        # three are None until check_alignment has seen a predicted container.
        self._predicted_mask: np.ndarray | None = None
        self._effective_mask: np.ndarray | None = None
        self._retained_reference: int | None = None

    # -- composition-time checks ---------------------------------------------

    def _check_output_kind(self) -> None:
        produced = self.instrument.output_kind
        if not isinstance(self.observed, produced):
            raise DatasetError(
                f"dataset {self.label!r}: the instrument chain produces "
                f"{produced.__name__} but the observed data are "
                f"{type(self.observed).__name__}. A likelihood compares like with like, so the "
                f"chain must end in the observed container's own kind — add the step that "
                f"produces one (synthetic photometry turns a Spectrum into PhotometricPoints), "
                f"or bind a channel that already holds it."
            )

    def _declare_latent(self, latent_name: str) -> LatentDeclaration | None:
        """The latent block this likelihood needs *for these data*, if any.

        Asked through
        :meth:`~ampere.core.likelihood.Likelihood.marginalisation_for` rather
        than the conservative
        :attr:`~ampere.core.likelihood.Likelihood.marginalisation` property,
        because ``likelihoods.md`` §9's rule is that masking beats censoring: a
        limit on a masked sample must not conjure a latent block, and must not
        cost the fit its choice of engine.

        The size is the number of samples the *observed* mask retains. The
        predicted container may mask more (a transformation is free to), so the
        count is verified against the real thing during
        :meth:`FittingProblem` validation, where a predicted template exists.
        """
        if self.likelihood.marginalisation_for(self.observed) is not Marginalisation.LATENT:
            return None
        return self.likelihood.latent_declaration(self.observed.n_valid, latent_name)

    def _components(self) -> dict[str, ParameterSet | ParameterMapping]:
        # The instrument joins as its mapping, not its merged set (lossless
        # nesting, ruled 2026-09-02): the dataset's bindings then compose down
        # to the instrument's steps, so an inner shared_as collapse stays
        # visible at every level above. This also snapshots the instrument's
        # mapping once, here — the nesting rule's freeze point.
        components: dict[str, ParameterSet | ParameterMapping] = {
            INSTRUMENT_COMPONENT: self.instrument.mapping,
            LIKELIHOOD_COMPONENT: self.likelihood.parameters,
        }
        if self._latent is not None:
            components[LATENT_COMPONENT] = self._latent.as_parameter_set()
        return components

    # -- declarations ---------------------------------------------------------

    @property
    def channel(self) -> str:
        """The model channel this dataset's instrument binds."""
        return self.instrument.channel

    @property
    def latent(self) -> LatentDeclaration | None:
        """The latent-GP declaration these data need, or ``None``."""
        return self._latent

    @property
    def mapping(self) -> ParameterMapping:
        """This dataset's own merge: instrument, likelihood and latent.

        Computed **once**, at construction, and retained — the nesting rule.
        A dataset is a composed object: construction freezes its instrument
        (``Instrument.freeze()``, ruled 2026-09-03), so reconfiguring one of
        its steps afterwards (``promote_buffer``, say) is refused at the next
        use rather than silently ignored. Build the dataset after configuring
        its pieces.
        """
        return self._mapping

    @property
    def parameters(self) -> ParameterSet:
        """The merged set an enclosing problem takes as one component.

        Names are qualified one level: ``instrument.calibrate.scale``,
        ``likelihood.amplitude``, ``latent.z``.
        """
        return self._mapping.merged

    @property
    def capability_parts(self) -> tuple[object, ...]:
        """The objects whose capability declarations this dataset depends on.

        The instrument steps, and — **since W2.13** (ruled 2026-09-07,
        ``inference.md`` §10a, fold-in 7) — the likelihood's own parts: its
        noise model, and its GP solver when a GP is declared. Before that
        widening a problem could report ``backend="jax"`` and
        ``differentiable=True`` while its GP solve factorised in scipy, which
        both backend tracks recorded as a finding; the flags now cover the
        whole of what one evaluation passes through.
        """
        return (*self.instrument.steps, *self.likelihood.capability_parts)

    @property
    def effective_mask(self) -> np.ndarray | None:
        """Which samples this dataset excludes, resolved once at construction.

        The union of the observed and predicted containers' masks, as a boolean
        array that is ``True`` where a sample is **excluded**; ``None`` when
        nothing is excluded, which is the common case and lets a caller skip
        the indexing entirely.

        Public since **W2.13** (ruled 2026-09-07, fold-in 9). It was private,
        and both backend tracks reached past the underscore for it anyway,
        because a native path needs exactly this: the mask is resolved at
        composition and is evaluation-invariant by declaration
        (:meth:`_resolve_mask` explains why that is contract rather than
        optimisation), so which samples are retained is a constant a
        differentiable log-density can close over rather than a value it must
        recompute inside the trace.

        Returns
        -------
        numpy.ndarray or None
            Read-only. A copy is returned rather than the stored array, so a
            caller cannot make the effective mask mutable state.
        """
        if self._effective_mask is None:
            return None
        mask = np.array(self._effective_mask, dtype=bool)
        mask.flags.writeable = False
        return mask

    # -- evaluation -----------------------------------------------------------

    def route(self, values: Mapping[str, Value] | ArrayLike) -> dict[str, dict[str, Value]]:
        """Split this dataset's values into its instrument, likelihood and latent halves.

        *values* are keyed by this dataset's **own** merged names
        (``likelihood.amplitude``), which is exactly what the enclosing
        problem's ``distribute`` hands over.
        """
        return self._mapping.distribute(self._mapping.merged.complete(values))

    def predict(
        self,
        result: ModelResult,
        values: Mapping[str, Value] | ArrayLike | None = None,
        *,
        routed: Mapping[str, Mapping[str, Value]] | None = None,
    ) -> FunctionSamples:
        """Push *result* through the instrument chain to predicted-data space.

        Separated from :meth:`log_likelihood` because the predicted container is
        wanted on its own: by :meth:`FittingProblem.simulate`, by the
        posterior-predictive and residual diagnostics of ``diagnostics.md``
        family B, and by the failure classifier, which needs to know whether a
        likelihood failed because the *prediction* was NaN.

        *routed* is this dataset's values already split by :meth:`route`. A
        caller that scores as well as predicts should route once and pass it to
        both, rather than paying for ``complete`` and ``distribute`` twice per
        evaluation — which, with a latent block of 10⁵ values, is two full-array
        copies on the path this contract is sized for.
        """
        split = self._split(values, routed)
        instrument_values = split.get(INSTRUMENT_COMPONENT) or None
        predicted = self.instrument(result, instrument_values)
        if not isinstance(predicted, FunctionSamples):  # pragma: no cover - Instrument checks it
            raise TransformationError(
                f"instrument {self.instrument.label!r} produced {predicted!r}, not a container."
            )
        return predicted

    def log_likelihood_of(
        self,
        predicted: FunctionSamples,
        values: Mapping[str, Value] | ArrayLike | None = None,
        *,
        routed: Mapping[str, Mapping[str, Value]] | None = None,
    ) -> float:
        """Score an already-computed prediction against the observed data.

        The prediction is handed to the likelihood carrying the **resolved**
        effective mask (:meth:`_resolve_mask`), so the likelihood's own mask
        union is an identity rather than a recomputation, and a
        parameter-dependent mask is refused here rather than silently changing
        which data are being fitted.
        """
        resolved_predicted = self._masked_pair(predicted)
        split = self._split(values, routed)
        latent = None
        if self._latent is not None:
            latent = np.asarray(split[LATENT_COMPONENT][self._latent.parameter.name])
        return float(
            self.likelihood.log_prob(
                resolved_predicted,
                self.observed,
                split.get(LIKELIHOOD_COMPONENT),
                latent=latent,
            )
        )

    def _censored_after_masking(self, retain: np.ndarray | None = None) -> bool:
        """Whether any limit survives the mask.

        With no *retain*, the observed container's own mask answers — the
        composition-time reading :meth:`_declare_latent` needs, before any
        prediction exists. A caller holding a prediction passes the
        **effective** inclusion indicator (the observed-and-predicted union
        ``likelihoods.md`` §8 defines), so a limit a prediction-side mask
        excludes does not count — masking beats censoring for the union too.
        (The asymmetry was found at the freeze's adversarial review:
        ``log_prob`` excised such a limit while ``draw_observation`` still
        refused because of it.)
        """
        censoring = self.likelihood.censoring
        if censoring is None:
            return False
        censoring.check_against(self.observed)
        included = (
            np.asarray(self.observed.valid).ravel()
            if retain is None
            else np.asarray(retain, dtype=bool).ravel()
        )
        return bool(np.any(np.asarray(censoring.kinds)[included] != 0))

    def _split(
        self,
        values: Mapping[str, Value] | ArrayLike | None,
        routed: Mapping[str, Mapping[str, Value]] | None,
    ) -> Mapping[str, Mapping[str, Value]]:
        if routed is not None:
            return routed
        return {} if values is None else self.route(values)

    def _resolve_mask(self, predicted: FunctionSamples) -> None:
        """Take the union of the observed and predicted masks — **once**.

        Peter's ruling on ``likelihoods.md`` §17 Q2: the effective mask is the
        ``Dataset``'s to resolve at construction, not the ``Likelihood``'s to
        recompute on every call. Two things follow, and both are contract, not
        optimisation.

        **The effective mask is evaluation-invariant.** A transformation whose
        output mask depends on parameter values is *unsupported*, and
        :meth:`_masked_pair` refuses one loudly rather than letting it through.
        That is not fussiness: a mask-controlling nuisance parameter has a free
        maximum at "mask everything", because ``Likelihood.log_prob`` scores a
        fully masked pair as exactly ``0.0`` — which beats every finite
        log-likelihood the dataset could otherwise contribute. Measured on a
        three-sample toy before this check existed, the joint ``log_prob`` rose
        from -255.1 to -9.2 as the mask closed, with no failure recorded: a
        sampler would have driven the data out of its own fit and converged
        happily on a posterior informed by nothing.

        **The latent block's size is fixed here too**, which the ruling makes
        explicit rather than implicit: ``latent_declaration(n)`` takes the
        retained count, and the retained count is now a construction-time
        constant by declaration.
        """
        self._predicted_mask = _mask_of(predicted)
        excluded = ~(np.asarray(self.observed.valid).ravel() & np.asarray(predicted.valid).ravel())
        self._effective_mask = excluded if bool(np.any(excluded)) else None
        self._retained_reference = int(np.count_nonzero(~excluded))

    def _masked_pair(self, predicted: FunctionSamples) -> FunctionSamples:
        """*predicted* carrying the resolved effective mask, invariance checked.

        The union having been taken once, this stamps the answer onto the
        prediction so that ``Likelihood``'s own ``weights()`` product becomes an
        identity rather than a recomputed union. The check comes first: stamping
        without it would *hide* a parameter-dependent mask instead of refusing
        it.
        """
        if self._retained_reference is None:
            return predicted
        current = _mask_of(predicted)
        if not _masks_equal(current, self._predicted_mask):
            now = 0 if current is None else int(np.count_nonzero(current))
            before = (
                0 if self._predicted_mask is None else int(np.count_nonzero(self._predicted_mask))
            )
            raise LikelihoodError(
                f"dataset {self.label!r}: the instrument chain produced a different mask at this "
                f"parameter vector than when the problem was composed ({now} sample(s) excluded, "
                f"was {before}). The effective mask is evaluation-invariant by declaration "
                f"(inference.md §8): a log-likelihood over one subset of the data is not "
                f"comparable with one over another, and because a fully masked dataset scores "
                f"exactly 0.0, 'mask everything' would be a free maximum a sampler will find. "
                f"Mask the affected samples on the observed container instead, so the set is "
                f"fixed before the run starts."
            )
        if _masks_equal(self._effective_mask, current):
            return predicted
        return predicted.with_values(predicted.values, mask=self._effective_mask)

    def log_likelihood(
        self, result: ModelResult, values: Mapping[str, Value] | ArrayLike | None = None
    ) -> float:
        """log p(observed | model result, parameters) for this dataset alone.

        Raises rather than returning ``-inf`` on a failure: converting a
        :class:`~ampere.core.exceptions.LikelihoodError` into §4.5's
        ``-inf``-with-a-reason is :class:`FittingProblem`'s job, because only
        the problem has somewhere to record the reason.
        """
        return self.log_likelihood_of(self.predict(result, values), values)

    def observed_coordinates(self) -> dict[str, Any]:
        """The observed container's own coordinates, per axis name.

        Build a resampling or response step from *these*, not from a
        recomputed grid.
        :meth:`~ampere.core.likelihood.Likelihood.check_alignment` compares the
        predicted and observed axes for **equality**, and two grids that agree
        to nine decimal places are not equal — so a step that recreates the
        observed wavelengths with ``np.linspace`` fails composition, correctly
        but confusingly. Handing the step the container's own array makes the
        comparison trivially true and removes a whole class of near-miss.

        Returned as :class:`~astropy.units.Quantity` per axis (a bare array
        where the axis has no unit), which is what a container constructor and a
        transformation's buffer both want.
        """
        return {axis.name: axis.quantity() for axis in self.observed.axes}

    def leaf_sites(self, name: str) -> tuple[str, ...]:
        """Every leaf this dataset's merged parameter *name* actually feeds.

        A thin wrapper since lossless nesting landed (ruled 2026-09-02): the
        dataset's mapping holds the instrument's mapping as an inner
        component, so :meth:`ParameterMapping.sites_of` already composes to
        the leaves — this method just renders the paths. The composition uses
        the instrument mapping *snapshotted at construction* (the nesting
        rule's freeze point), so a step reconfigured afterwards is not picked
        up, consistently with everything else about a built dataset.
        """
        return tuple(
            f"{binding.component}{SEPARATOR}{binding.local_name}"
            for binding in self._mapping.sites_of(name)
        )

    # -- composition-time obligations -----------------------------------------

    def check_alignment(self, predicted: FunctionSamples) -> None:
        """``likelihoods.md`` §16's obligation, discharged.

        Calls :meth:`~ampere.core.likelihood.Likelihood.check_alignment` — the
        O(N) pass that belongs once, at composition — and then verifies the
        latent block against the mask the *predicted* container actually
        carries. Not optional: it is where an unimplemented latent combination
        is refused, so skipping it silently accepts a problem that cannot be
        evaluated.
        """
        self.likelihood.check_alignment(predicted, self.observed)
        self._resolve_mask(predicted)
        retained = self._retained_reference or 0
        if self._latent is None:
            return
        if retained != self._latent.size:
            raise DatasetError(
                f"dataset {self.label!r} declared {self._latent.size} latent value(s) from the "
                f"observed mask, but only {retained} sample(s) survive the union of the observed "
                f"and predicted masks. A latent-GP formulation has exactly one latent value per "
                f"retained sample, so the sampler dimension would not match the likelihood. The "
                f"usual cause is a transformation that masks samples the data do not; mask them "
                f"on the observed container too, so the count is fixed before the run starts."
            )

    def retained(self, predicted: FunctionSamples) -> int:
        """How many samples survive the union of the observed and predicted masks.

        The mask convention is ``likelihoods.md`` §8's, restated rather than
        reimplemented: a sample counts only if *both* containers call it valid.
        """
        weights = (
            np.asarray(self.observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        )
        return int(np.count_nonzero(weights > 0.0))

    def check_engine(self, *, differentiable: bool, engine: str = "this engine") -> None:
        """Refuse an engine that cannot deliver this likelihood's marginalisation.

        Passes ``observed=`` — which ``likelihoods.md`` §16 marks as
        load-bearing — so a censored sample the mask excludes is not counted
        against a gradient-free engine.
        """
        self.likelihood.check_engine(
            differentiable=differentiable, engine=engine, observed=self.observed
        )

    # -- simulation -----------------------------------------------------------

    def draw_observation(
        self,
        predicted: FunctionSamples,
        values: Mapping[str, Value] | ArrayLike | None,
        rng: np.random.Generator,
    ) -> FunctionSamples:
        """Draw one noisy realisation of *predicted* under this dataset's likelihood.

        Delegated to :meth:`LikelihoodFamily.sample` (ruled 2026-09-02,
        ``inference.md`` §19 R3), which is the generative counterpart of
        ``log_prob`` and receives the same :class:`NoiseParams` the likelihood
        scores with — so a fitted ``scale`` or ``jitter``, a GP kernel and the
        solver's own stabiliser are all in the draw exactly as they are in the
        density. :class:`GaussianFamily` implements it for both noise models; a
        family that does not implement it refuses **specifically**, naming the
        override a user should provide for an exotic observation process.

        Masked samples keep the observed container's own values: they carry zero
        information and are excluded from every likelihood, so drawing noise for
        them would be inventing data. A censoring declaration blocks a draw only
        when a limit **survives the mask** — ``likelihoods.md`` §9's rule that
        masking beats censoring, applied here as it already is in
        :meth:`_declare_latent`, so that one class does not give two answers to
        one question.
        """
        family = self.likelihood.family
        noise = self.likelihood.noise
        observed = self.observed
        routed = {} if values is None else self.route(values)
        resolved = routed.get(LIKELIHOOD_COMPONENT, {})
        weights = np.asarray(observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
        retain = weights > 0.0
        # Checked against the *effective* mask, prediction side included, so a
        # limit the union excludes does not block the draw — the same answer
        # log_prob's excision gives (masking beats censoring, likelihoods.md §9).
        if self._censored_after_masking(retain):
            raise DatasetError(
                f"dataset {self.label!r}: a censoring declaration on retained samples blocks "
                f"observation drawing — a limit is part of the observation process, and applying "
                f"the censoring operator to a draw is not implemented. Use "
                f"simulate(observe=False) and draw your own observations from the predicted "
                f"containers."
            )
        drawn = np.array(np.asarray(observed.values).ravel(), dtype=DTYPE, copy=True)
        if not np.any(retain):
            return observed.with_values(drawn.reshape(observed.shape))

        coordinates = np.column_stack(
            [np.asarray(axis.values, dtype=DTYPE) for axis in observed.axes]
        )[retain]
        realisation = np.asarray(predicted.values, dtype=DTYPE).ravel()[retain]
        # The NoiseParams are built from the noiseless prediction *before* noise
        # is added (ruled 2026-09-03, X-1 point 4): a prediction-dependent noise
        # scales with the true curve — the standard generative reading — and the
        # draw is then consistent with the density that will score it.
        params = noise.noise_params(
            observed, retain, resolved, predicted=realisation, coordinates=coordinates
        )
        try:
            realisation = family.sample(realisation, params, rng)
        except LikelihoodError as error:
            raise DatasetError(f"dataset {self.label!r}: {error}") from error
        drawn[retain] = np.asarray(realisation, dtype=DTYPE)
        return observed.with_values(drawn.reshape(observed.shape))

    def __repr__(self) -> str:
        latent = "" if self._latent is None else f", {self._latent.size} latent"
        return (
            f"<Dataset {self.label!r}: {type(self.observed).__name__}"
            f"[{self.observed.n_samples}] on channel {self.channel!r}, "
            f"{self.parameters.free_size} free parameter(s){latent}>"
        )


# ---------------------------------------------------------------------------
# DatasetCollection
# ---------------------------------------------------------------------------


class DatasetCollection(Mapping[str, Dataset]):
    """Several datasets fitted jointly, and the parameters they share.

    A :class:`collections.abc.Mapping` from label to :class:`Dataset`, so
    ``len``, iteration, ``in``, ``keys``/``values``/``items`` and ``.get`` all
    behave as expected.

    It owns two things and no more. The **joint log-likelihood is a sum** — the
    datasets are conditionally independent given the parameters, which is what
    makes a joint fit a joint fit and is the one modelling assumption this class
    makes. And the **labels namespace the joint parameter space**: each dataset
    becomes one component of the enclosing problem's merge, so two
    independently written instruments may both call their nuisance ``scale``
    without colliding.

    It deliberately does **not** own the ties. A tie may name a model site as
    easily as a dataset one, so it belongs at the only level that can see every
    site: :class:`FittingProblem`.

    Parameters
    ----------
    datasets
        A mapping of label to dataset, or an iterable of datasets (each
        contributing its own :attr:`Dataset.label`).
    shared
        The hyperprior extension point: a :class:`~ampere.core.parameter.
        ParameterSet` joining the merge as one further top-level component. Use
        it for population hyperparameters — declared here, then tied to the
        per-dataset sites that use them, or built from a
        :class:`~ampere.core.parameter.Plate`. See ``inference.md`` §9 for what
        this can and cannot express.
    shared_label
        Component label for *shared*.

    Examples
    --------
    >>> import astropy.units as u
    >>> from ampere.core import PhotometricPoints, Spectrum
    >>> sed = PhotometricPoints(
    ...     ["W1", "W2"], [3.4, 4.6] * u.um, [1.0, 2.0] * u.Jy, uncertainty=[0.1] * 2 * u.Jy
    ... )
    >>> spectrum = Spectrum(
    ...     [10.0, 11.0, 12.0] * u.um, [3.0, 3.5, 4.0] * u.Jy, uncertainty=[0.2] * 3 * u.Jy
    ... )
    >>> datasets = DatasetCollection(
    ...     {"sed": Dataset(sed, model="star"), "spectrum": Dataset(spectrum, model="star")}
    ... )
    >>> list(datasets)
    ['sed', 'spectrum']
    >>> sorted(datasets.components())
    ['sed', 'spectrum']
    """

    __slots__ = ("_datasets", "_shared", "_shared_label")

    def __init__(
        self,
        datasets: Mapping[str, Dataset] | Iterable[Dataset],
        *,
        shared: ParameterSet | None = None,
        shared_label: str = SHARED_COMPONENT,
    ) -> None:
        built: dict[str, Dataset] = {}
        pairs: Iterable[tuple[str, Dataset]]
        if isinstance(datasets, Mapping):
            pairs = datasets.items()
        elif isinstance(datasets, Dataset):
            raise DatasetError(
                "a DatasetCollection is built from several datasets; got a single Dataset. Wrap "
                "it in a list — [dataset] — or pass a mapping of label to dataset."
            )
        else:
            pairs = ((dataset.label, dataset) for dataset in datasets)
        for label, dataset in pairs:
            if not isinstance(dataset, Dataset):
                raise DatasetError(
                    f"a DatasetCollection holds Datasets, got {type(dataset).__name__} under "
                    f"label {label!r}."
                )
            checked = _check_label(label, "dataset label")
            if checked in built:
                # Loudly, and *before* the merge. Two component sets filed under
                # one key in a plain dict would silently discard the first, so
                # ParameterSet.merge's own collision detection would never fire
                # and one instrument's nuisance parameters would simply vanish
                # from the fit. This is not a hypothetical: Instrument.label
                # defaults to the channel name, so two catalogues observing one
                # SED channel collide *by default* unless one is named.
                hint = (
                    " Both labels came from an Instrument's own label, which itself defaults to "
                    "the channel name — so two instruments bound to one channel collide unless "
                    "you name them. Pass label='...' to the Dataset (or to the Instrument), or "
                    "build the collection from a mapping whose keys are the labels."
                    if dataset.label == checked
                    else ""
                )
                raise DatasetError(
                    f"two datasets are labelled {checked!r}. Labels are the component names of "
                    f"the joint parameter space, so they must be unique.{hint}"
                )
            built[checked] = dataset
        if not built:
            raise DatasetError(
                "a DatasetCollection needs at least one dataset; a fit with no data has no "
                "likelihood to evaluate."
            )
        self._datasets = built

        if shared is not None and not isinstance(shared, ParameterSet):
            raise DatasetError(
                f"a DatasetCollection's shared parameters must be a ParameterSet, got "
                f"{type(shared).__name__}. Build one from Parameters, or from Plates with "
                f"ParameterSet(plates=[...])."
            )
        self._shared = shared
        self._shared_label = _check_label(shared_label, "shared-parameter label")
        if shared is not None and self._shared_label in built:
            raise DatasetError(
                f"the shared parameters are labelled {self._shared_label!r}, which is also a "
                f"dataset label. Pass shared_label='...' to separate them."
            )

    # -- Mapping protocol -----------------------------------------------------

    def __getitem__(self, label: str) -> Dataset:
        try:
            return self._datasets[label]
        except KeyError:
            raise KeyError(
                f"no dataset labelled {label!r}; this collection holds {list(self._datasets)}."
            ) from None

    def __iter__(self) -> Iterator[str]:
        return iter(self._datasets)

    def __len__(self) -> int:
        return len(self._datasets)

    # -- declarations ---------------------------------------------------------

    @property
    def shared(self) -> ParameterSet | None:
        """The shared/hyperprior parameter set, if one was declared."""
        return self._shared

    @property
    def shared_label(self) -> str:
        """Component label of :attr:`shared`."""
        return self._shared_label

    def components(self) -> dict[str, ParameterSet | ParameterMapping]:
        """The component sets this collection contributes to the joint merge.

        Handed to :class:`FittingProblem`, which adds the models and performs
        the single merge. This class never merges on a problem's behalf: doing
        so would produce a mapping the problem then had to merge *again*, which
        is the associativity trap ``parameters.md`` §12.4 names.
        """
        components: dict[str, ParameterSet | ParameterMapping] = {
            label: dataset.mapping for label, dataset in self._datasets.items()
        }
        if self._shared is not None:
            components[self._shared_label] = self._shared
        return components

    @property
    def capability_parts(self) -> tuple[object, ...]:
        """Every object whose capability declarations the collection depends on."""
        return tuple(
            part for dataset in self._datasets.values() for part in dataset.capability_parts
        )

    def model_labels(self) -> tuple[str | None, ...]:
        """The model each dataset names, in collection order (``None`` = unstated)."""
        return tuple(dataset.model for dataset in self._datasets.values())

    # -- evaluation -----------------------------------------------------------

    def contributions(
        self,
        results: Mapping[str, ModelResult],
        routed: Mapping[str, Mapping[str, Value]],
        *,
        models: Mapping[str, str] | None = None,
    ) -> dict[str, float]:
        """Each dataset's log-likelihood, by label.

        Raises on failure; :class:`FittingProblem` is where a failure becomes
        ``-inf`` with a recorded reason, because only the problem has a place to
        record it.

        Parameters
        ----------
        results
            The evaluated :class:`~ampere.core.results_schema.ModelResult` per
            model label.
        routed
            This collection's half of the problem's ``distribute`` output:
            dataset label to that dataset's own merged values.
        models
            Dataset label to model label. Defaults to each dataset's own
            :attr:`Dataset.model`.
        """
        contributions: dict[str, float] = {}
        for label, dataset in self._datasets.items():
            model_label = dataset.model if models is None else models[label]
            if model_label is None or model_label not in results:
                raise DatasetError(
                    f"dataset {label!r} names model {model_label!r}, which was not evaluated; "
                    f"the available results are {sorted(results)}."
                )
            contributions[label] = dataset.log_likelihood(results[model_label], routed.get(label))
        return contributions

    def log_likelihood(
        self,
        results: Mapping[str, ModelResult],
        routed: Mapping[str, Mapping[str, Value]],
        *,
        models: Mapping[str, str] | None = None,
    ) -> float:
        """The joint log-likelihood: the sum over datasets."""
        return float(sum(self.contributions(results, routed, models=models).values()))

    def __repr__(self) -> str:
        shared = "" if self._shared is None else f", shared={self._shared_label!r}"
        return f"<DatasetCollection {list(self._datasets)}{shared}>"


# ---------------------------------------------------------------------------
# FittingProblem
# ---------------------------------------------------------------------------


class FittingProblem:
    """Model, data and priors composed into the object an engine runs.

    This is ``DEVELOPMENT_PLAN.md`` §4.5's contract. Any engine consuming only
    this surface — emcee, dynesty, zeus, an SBI package, a NUTS kernel on a
    modern backend — works with every backend and with legacy black-box models
    through a thin adapter, because nothing below is backend-specific.

    Parameters
    ----------
    model
        One :class:`~ampere.core.transform.Model`, or a mapping of label to
        model when several objects are fitted together (each dataset then names
        the one it observes through :attr:`Dataset.model`).
    datasets
        A :class:`DatasetCollection`, or anything one can be built from.
    ties
        Composition-time :class:`~ampere.core.parameter.Tie`\\ s. Declared here
        because this is the only level at which every site is visible; a tie
        names sites by their **full merged path**, e.g.
        ``"sed.instrument.calibrate.scale"``.
    model_label
        Component label for a single model.
    seed
        The run's one integer seed (``lowering.md`` §9.2). Named sub-streams are
        derived from it by :func:`~ampere.core.rng.substream`, so prior
        sampling, initialisation and simulation never share randomness.
        ``None`` means non-reproducible: every stream is freshly entropy-seeded.
    capabilities
        Override the flags derived from what the models and transformations
        declare. Pass this when an adapter knows something the pieces cannot say
        for themselves, or to settle a deliberate device or backend
        disagreement the pieces cannot settle between them.
    reference_values
        The θ used for composition-time validation. Defaults to the **prior
        median** — ``prior_transform`` at the centre of the unit cube — which is
        deterministic, always inside the support, and correct for hierarchical
        priors because ``prior_transform`` already evaluates them in dependency
        order.
    validate
        Run the composition-time checks eagerly (the default). ``False`` defers
        them to the first evaluation, where they still run exactly once — for a
        model whose single evaluation is expensive. They are never skipped:
        ``likelihoods.md`` §16 is explicit that ``check_alignment`` is not
        optional.
    strict
        Let every exception propagate instead of becoming ``-inf`` with a
        recorded reason. The engine-facing default is ``False``; ``True`` is for
        direct use and debugging (Peter's ruling on ``likelihoods.md`` §17 Q1).
        The intended workflow is to run non-strict, read
        :meth:`failure_summary`, then re-run with ``strict=True`` to get the
        raise at the offending draw with a full traceback.
    lenient_compile
        Downgrade a model's ``compile_for`` refusal
        (:class:`~ampere.core.exceptions.CompositionError`) to a warning and
        proceed with the unconfigured model. ``False`` by default — ruled
        2026-09-03 (``transformations.md`` §15 Q3): negotiation refuses an
        unachievable requirement by raising, and silencing that is an
        explicit, per-problem decision.
    simulator_failures
        Extra exception types that count as an unscoreable point rather than a
        bug. :class:`~ampere.core.exceptions.LikelihoodError` is always included
        — that is ``likelihoods.md`` §16's explicit hand-over of the
        non-positive-definite-covariance case — and **nothing else is**. A
        ``TypeError`` from a mis-typed model, a
        :class:`~ampere.core.exceptions.SchemaError` from a malformed container
        and a :class:`~ampere.core.exceptions.TransformationError` from a chain
        that drops a mask are all bugs, and they surface as themselves rather
        than as a mysteriously flat posterior. Declare here the exceptions a
        wrapped external code actually raises (``subprocess.CalledProcessError``,
        the wrapper's own ``SimulationFailed``), which is §4.5's
        "external simulators crash" case.
    failure_history
        How many recent :class:`Failure` records to keep.

    Examples
    --------
    See the module docstring for the two-line case, and
    ``docs/design/contracts/inference.md`` for the joint fit with a tie.
    """

    def __init__(
        self,
        model: Model | Mapping[str, Model],
        datasets: DatasetCollection | Mapping[str, Dataset] | Iterable[Dataset],
        *,
        ties: Sequence[Tie] = (),
        model_label: str = MODEL_COMPONENT,
        seed: int | None = None,
        capabilities: Capabilities | None = None,
        reference_values: Mapping[str, Value] | ArrayLike | None = None,
        validate: bool = True,
        strict: bool = False,
        lenient_compile: bool = False,
        simulator_failures: Sequence[type[BaseException]] = (),
        failure_history: int = DEFAULT_FAILURE_HISTORY,
    ) -> None:
        self._models = _normalise_models(model, model_label)
        self.datasets = (
            datasets if isinstance(datasets, DatasetCollection) else DatasetCollection(datasets)
        )
        self._bindings = self._bind_models()
        self.ties = tuple(ties)
        for tie in self.ties:
            if not isinstance(tie, Tie):
                raise DatasetError(f"ties must be Tie instances, got {type(tie).__name__}.")

        # Loud, at composition (found by the freeze's adversarial review):
        # int(1.9) silently truncating, True counting as 1, or a seed outside
        # substream's signed 64-bit derivation crashing at the first stream
        # request would all contradict the reproducibility contract this seed
        # exists for (lowering.md §9.2 — every backend derives its streams
        # from the same bytes).
        if seed is None:
            self.seed = None
        else:
            if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)):
                raise DatasetError(
                    f"a run seed must be an integer, got {seed!r}. substream "
                    f"(lowering.md §9.2) derives every named stream from it, and ampere will "
                    f"not guess what a non-integer seed means."
                )
            bound = 1 << (8 * SEED_BYTES - 1)
            if not -bound <= int(seed) < bound:
                raise DatasetError(
                    f"seed {seed!r} does not fit substream's {8 * SEED_BYTES}-byte signed "
                    f"derivation. Every backend derives its streams from the same "
                    f"{SEED_BYTES}-byte encoding, so the range is part of the "
                    f"reproducibility contract."
                )
            self.seed = int(seed)
        self._streams: dict[str, np.random.Generator] = {}
        for candidate in simulator_failures:
            if not (isinstance(candidate, type) and issubclass(candidate, BaseException)):
                raise DatasetError(
                    f"simulator_failures takes exception classes, got {candidate!r}."
                )
        # Deliberately narrow. LikelihoodError is here because likelihoods.md
        # §16 hands this contract the non-positive-definite covariance case by
        # name; everything else ampere raises is a composition bug, and turning
        # a bug into -inf produces a fit that runs, converges and is wrong.
        #
        # strict=True empties the set, so every exception propagates. That is
        # Peter's ruling on likelihoods.md §17 Q1: the engine path records and
        # returns -inf (exception control flow cannot be traced by jax, so the
        # non-raising path is forced by Phase 2 regardless), while strict
        # raising stays available for direct use and debugging. The intended
        # workflow is to run non-strict, read failure_summary(), then re-run
        # strict to get the raise at the offending draw. Emptying the tuple —
        # rather than branching at each call site — is also what lets a future
        # solver-level strict flag replace this try/except without a contract
        # change: the call sites stay as they are.
        self.strict = bool(strict)
        self._failure_types: tuple[type[BaseException], ...] = (
            () if self.strict else (LikelihoodError, *tuple(simulator_failures))
        )
        self._failures: collections.deque[Failure] = collections.deque(maxlen=int(failure_history))
        self._failure_counts: collections.Counter[FailureReason] = collections.Counter()

        # (2)-(3) Negotiate, then compile — once, before anything is merged, so
        # that a model which reconfigures itself for its instruments is the one
        # whose parameters enter the joint space.
        self._lenient_compile = bool(lenient_compile)
        self._requirements = self._negotiate()
        self._compiled = {
            label: self._compile(label, instance) for label, instance in self._models.items()
        }

        # (4) The one merge.
        self._mapping = ParameterSet.merge(self._components(), ties=self.ties)
        self._require_resolved()

        self.capabilities = (
            declared_capabilities(self._capability_parts())
            if capabilities is None
            else capabilities
        )
        if not isinstance(self.capabilities, Capabilities):
            raise DatasetError(
                f"capabilities must be a Capabilities record, got "
                f"{type(self.capabilities).__name__}."
            )

        self.reference_values: dict[str, Value] = (
            self._prior_median()
            if reference_values is None
            else dict(self._mapping.merged.complete(reference_values))
        )

        # (5) Validate.
        self._validated = False
        if validate:
            self.validate()

    # -- construction helpers -------------------------------------------------

    def _bind_models(self) -> dict[str, str]:
        """Which model each dataset observes; the default is only safe if unique."""
        bindings: dict[str, str] = {}
        sole = next(iter(self._models)) if len(self._models) == 1 else None
        for label, dataset in self.datasets.items():
            named = dataset.model
            if named is None:
                if sole is None:
                    raise DatasetError(
                        f"dataset {label!r} does not say which model it observes, and this "
                        f"problem holds {sorted(self._models)}. Pass model='...' to the Dataset; "
                        f"'the only model' is only well defined when there is one."
                    )
                bindings[label] = sole
            elif named not in self._models:
                raise DatasetError(
                    f"dataset {label!r} observes model {named!r}, which this problem does not "
                    f"hold; its models are {sorted(self._models)}."
                )
            else:
                bindings[label] = named
        orphans = sorted(set(self._models) - set(bindings.values()))
        if orphans:
            raise DatasetError(
                f"model(s) {orphans} are held by this problem but observed by no dataset, so "
                f"their parameters would be free sampler dimensions that no likelihood "
                f"constrains — the posterior would simply return their priors, and the model "
                f"would be re-evaluated on every log_prob for nothing. Give each one a dataset, "
                f"or drop it. The mirror of this check (a dataset naming no model when there are "
                f"several) is above."
            )
        return bindings

    def _negotiate(self) -> dict[str, dict[str, Any]]:
        """Step (2): the instruments' requirements, per model, per channel.

        Also the home of the instrument-label check (ruled 2026-09-03,
        ``transformations.md`` §15 Q4's residual): when more than one
        instrument reads a channel, distinct instrument labels are required —
        the label is how a user identifies which parameters are constrained
        by which data, and the requirements-provenance ``sources`` tuple is
        otherwise ambiguous. Checked here, at problem composition; the
        one-instrument default-to-channel-name case is unchanged, and
        ``negotiate`` itself merges nothing and stays silent about labels.
        """
        collected: dict[str, dict[str, Any]] = {}
        for label in self._models:
            instruments = [
                dataset.instrument
                for name, dataset in self.datasets.items()
                if self._bindings[name] == label
            ]
            by_channel: dict[str, list[str]] = {}
            for instrument in instruments:
                by_channel.setdefault(instrument.channel, []).append(instrument.label)
            for channel, labels in by_channel.items():
                duplicates = sorted({name for name in labels if labels.count(name) > 1})
                if duplicates:
                    raise DatasetError(
                        f"model {label!r}: {len(labels)} instruments read channel {channel!r} "
                        f"but share the instrument label(s) {duplicates}. Instrument.label "
                        f"defaults to the channel name, so several unnamed instruments on one "
                        f"channel are indistinguishable — in the requirements provenance "
                        f"(sources) and everywhere a user asks which instrument constrained "
                        f"what. Pass label='...' to each Instrument (the single-instrument "
                        f"default is unchanged)."
                    )
            collected[label] = dict(negotiate(instruments))
        return collected

    def _compile(self, label: str, instance: Model) -> Model:
        """Step (3): the one-off ``compile_for``, and the check that it behaved.

        A model that engages with a requirement it cannot honour raises
        ``CompositionError`` (ruled 2026-09-03, ``transformations.md`` §15
        Q3), and by default that refusal propagates — an unachievable
        requirement is a composition problem, not a warning. The explicit
        opt-out is ``lenient_compile=True``: the refusal is downgraded to a
        warning and the *unconfigured* model is used, which is the same
        declared stance as a model that ignores negotiation entirely.
        """
        try:
            compiled = instance.compile_for(self._requirements[label])
        except CompositionError as error:
            if not self._lenient_compile:
                raise
            warnings.warn(
                f"model {label!r} refused its instruments' requirements ({error}); proceeding "
                f"with the unconfigured model because lenient_compile=True. The instruments "
                f"may now be handed a grid they cannot use.",
                stacklevel=2,
            )
            return instance
        if not isinstance(compiled, Model):
            raise DatasetError(
                f"model {label!r}'s compile_for() returned {compiled!r}, not a Model. It must "
                f"return self (having stored its templates) or a new, configured instance — "
                f"transformations.md §7."
            )
        return compiled

    def _components(self) -> dict[str, ParameterSet | ParameterMapping]:
        """Step (4)'s input: every top-level component, checked for collisions."""
        components: dict[str, ParameterSet] = {}
        for label, instance in self._compiled.items():
            components[label] = instance.parameters
        for label, parameters in self.datasets.components().items():
            if label in components:
                raise DatasetError(
                    f"{label!r} labels both a model and a dataset (or the shared parameters). "
                    f"Top-level component labels share one namespace, because they qualify one "
                    f"joint parameter space; rename one, or pass model_label='...'."
                )
            components[label] = parameters
        return components

    def _require_resolved(self) -> None:
        if self._mapping.merged.is_resolved:
            return
        raise DatasetError(
            f"parameter(s) {list(self._mapping.merged.deferred_names)} are declared shared "
            f"(shared_as=...) but no site in this problem supplies their prior, so the joint "
            f"space cannot be sampled. Either give one site the prior, or collapse them with an "
            f"explicit Tie(name, sites, prior=...). Note that a shared_as label is resolved by "
            f"the merge of the level that owns both sites: two steps of one instrument are "
            f"merged by the instrument, and a label whose partner lives in another dataset must "
            f"be expressed as a Tie here instead (inference.md §6)."
        )

    def _prior_median(self) -> dict[str, Value]:
        median = self._mapping.merged.prior_transform(np.full(self.free_size, 0.5))
        return self._mapping.merged.unpack(median)

    def _capability_parts(self) -> tuple[object, ...]:
        return (*self._compiled.values(), *self.datasets.capability_parts)

    # -- declarations ---------------------------------------------------------

    @property
    def models(self) -> Mapping[str, Model]:
        """The compiled models, by label. These, not the originals, are evaluated."""
        return types.MappingProxyType(self._compiled)

    @property
    def model(self) -> Model:
        """The only model, for the single-model case."""
        if len(self._compiled) != 1:
            raise DatasetError(
                f"this problem holds {len(self._compiled)} models {sorted(self._compiled)}; ask "
                f"for the one you mean by name through .models."
            )
        return next(iter(self._compiled.values()))

    @property
    def requirements(self) -> Mapping[str, Mapping[str, Any]]:
        """What the instruments asked of each model, per channel — the negotiation record."""
        return types.MappingProxyType(
            {
                label: types.MappingProxyType(dict(value))
                for label, value in self._requirements.items()
            }
        )

    @property
    def bindings(self) -> Mapping[str, str]:
        """Dataset label to model label."""
        return types.MappingProxyType(self._bindings)

    @property
    def mapping(self) -> ParameterMapping:
        """The one merge: the joint set plus the wiring back to every component."""
        return self._mapping

    @property
    def parameters(self) -> ParameterSet:
        """The joint parameter set an engine samples."""
        return self._mapping.merged

    @property
    def free_size(self) -> int:
        """The engine's dimension."""
        return self._mapping.merged.free_size

    def free_labels(self) -> tuple[str, ...]:
        """One label per sampler dimension — corner plots, ArviZ coordinates (W1.8)."""
        return self._mapping.merged.free_labels()

    @property
    def tied_names(self) -> tuple[str, ...]:
        """Merged names with more than one binding **at the top level**.

        This reports the ties this problem itself resolved. It does *not* report
        a parameter shared further down — two steps of one instrument sharing a
        declaration-time label are collapsed by the instrument's own merge, so
        the top-level mapping sees one binding and says "not tied". Use
        :attr:`shared_names`, which descends, whenever the question is "is this
        parameter shared?" rather than "did this problem's ties fire?".
        """
        return self._mapping.tied_names

    def sites(self) -> Mapping[str, tuple[str, ...]]:
        """Every leaf each merged parameter feeds, by fully qualified path.

        A thin wrapper since lossless nesting landed (ruled 2026-09-02): the
        problem's own :attr:`~ampere.core.parameter.ParameterMapping.bindings`
        already compose through every retained inner mapping, so this groups
        them by merged name and renders the paths. ``mapping.bindings`` and
        ``mapping.tied_names`` now tell the leaf-level truth themselves; this
        method remains as the convenient rendered form.

        Returns a mapping from merged name to the paths it reaches, e.g.
        ``{"calibration": ("sed.instrument.calibrate.scale",
        "spectrum.instrument.calibrate.scale")}``. Either this or the raw
        composed bindings is right for a run's provenance — they carry the
        same information.
        """
        found: dict[str, list[str]] = {name: [] for name in self._mapping.merged.names}
        for binding in self._mapping.bindings:
            found[binding.global_name].append(f"{binding.component}{SEPARATOR}{binding.local_name}")
        return types.MappingProxyType({name: tuple(paths) for name, paths in found.items()})

    @property
    def shared_names(self) -> tuple[str, ...]:
        """Merged names driving more than one leaf, at **any** level.

        Since lossless nesting landed (ruled 2026-09-02) this is the same
        statement :attr:`tied_names` makes — the mapping's own reporting is
        leaf-level now — and it is kept as the established name for it.
        """
        return self._mapping.tied_names

    # -- the engine-facing surface (§4.5) -------------------------------------

    def log_prior(self, values: Mapping[str, Value] | ArrayLike | None = None) -> float:
        """log p(θ). ``-inf`` outside the support."""
        return self._mapping.merged.lnprior(self._resolve(values))

    def log_likelihood(self, values: Mapping[str, Value] | ArrayLike | None = None) -> float:
        """log p(data | θ), summed over datasets. ``-inf`` if the point cannot be scored."""
        resolved = self._resolve(values)
        routed = self._mapping.distribute(resolved)
        total, _contributions, failure = self._joint_log_likelihood(routed)
        self._record(failure)
        return total

    def log_prob(self, values: Mapping[str, Value] | ArrayLike | None = None) -> float:
        """log p(θ) + log p(data | θ) — the number a sampler maximises.

        Returns ``-inf`` for a point outside the prior's support **without
        evaluating the model**, which for an expensive simulator is the single
        most valuable thing this contract does. Also ``-inf``, with a recorded
        reason, for a point the forward model or the likelihood could not score;
        :meth:`evaluate` returns the reason alongside the number.
        """
        return self.evaluate(values).log_prob

    def evaluate(self, values: Mapping[str, Value] | ArrayLike | None = None) -> Evaluation:
        """:meth:`log_prob` with the split, the per-dataset terms and the reason.

        This is the method an engine driver should call: it returns everything
        W1.8 must store per posterior draw (``DEVELOPMENT_PLAN.md`` §4.6) in one
        pass, and it is where §4.5's "``-inf`` with a recorded reason" is
        actually delivered.
        """
        resolved = self._resolve(values)
        log_prior = self._mapping.merged.lnprior(resolved)
        if not math.isfinite(log_prior):
            # Zero prior mass is an answer, not a failure: no model is run and
            # nothing is recorded. log_likelihood is NaN because it was not
            # evaluated, which is a different statement from -inf.
            return Evaluation(log_prior=-math.inf, log_likelihood=math.nan, log_prob=-math.inf)
        routed = self._mapping.distribute(resolved)
        total, contributions, failure = self._joint_log_likelihood(routed)
        self._record(failure)
        log_prob = log_prior + total if math.isfinite(total) else -math.inf
        return Evaluation(
            log_prior=log_prior,
            log_likelihood=total,
            log_prob=log_prob if math.isfinite(log_prob) else -math.inf,
            contributions=contributions,
            failure=failure,
        )

    def prior_transform(self, unit_cube: ArrayLike) -> np.ndarray:
        """Map a unit-hypercube point to a free-parameter vector — nested sampling.

        Delegated to the merged set, which evaluates hierarchical priors in
        dependency order. Deterministic and backend-independent, which is why
        ``lowering.md`` §9.3 names it the right tool where determinism actually
        matters.
        """
        return self._mapping.merged.prior_transform(unit_cube)

    def unconstrain(self, values: Mapping[str, Value] | ArrayLike) -> np.ndarray:
        """Map a parameter vector into unconstrained space (gradient-based engines)."""
        return self._mapping.merged.unconstrain(values)

    def constrain(self, unconstrained: ArrayLike) -> np.ndarray:
        """The inverse of :meth:`unconstrain`."""
        return self._mapping.merged.constrain(unconstrained)

    def log_prob_unconstrained(self, unconstrained: ArrayLike) -> float:
        """:meth:`log_prob` in unconstrained space, change-of-variables term included.

        The density an HMC/NUTS kernel works with. Stated here on the reference
        path so the torch and jax lowerings have an oracle to agree with rather
        than each rediscovering the Jacobian term (``parameters.md`` §6).
        """
        log_prior = self._mapping.merged.lnprior_unconstrained(unconstrained)
        if not math.isfinite(log_prior):
            return -math.inf
        constrained = self._mapping.merged.constrain(unconstrained)
        routed = self._mapping.distribute(self._mapping.merged.unpack(constrained))
        total, _contributions, failure = self._joint_log_likelihood(routed)
        self._record(failure)
        combined = log_prior + total
        return combined if math.isfinite(combined) else -math.inf

    def sample_prior(self, rng: np.random.Generator | None = None) -> dict[str, Value]:
        """One draw from the joint prior, on the ``"prior"`` sub-stream."""
        return self._mapping.merged.sample(self.rng("prior") if rng is None else rng)

    # -- capability flags and engine checking ---------------------------------

    @property
    def differentiable(self) -> bool:
        """Whether ``log_prob`` admits a gradient (``DEVELOPMENT_PLAN.md`` §4.5)."""
        return self.capabilities.differentiable

    @property
    def batchable(self) -> bool:
        """Whether ``log_prob`` accepts a stack of parameter vectors."""
        return self.capabilities.batchable

    @property
    def device(self) -> str:
        """Where this problem's arrays live."""
        return self.capabilities.device

    @property
    def backend(self) -> str:
        """Which backend supplies this problem's pieces (W2.12).

        Derived from what the models and transformations declare, never
        asserted by whoever runs the fit: it is what
        :class:`~ampere.inference.engine.Engine` records as
        ``ampere_backend``, which is why that attribute is a fact.
        """
        return self.capabilities.backend

    def check_engine(
        self, engine: str = "this engine", *, differentiable: bool | None = None
    ) -> None:
        """Refuse an engine no dataset's likelihood can be run by.

        Discharges ``likelihoods.md`` §16's second obligation for every dataset,
        passing ``observed=`` so that a censored sample the mask excludes is not
        counted against a gradient-free engine.

        Parameters
        ----------
        engine
            Named in the error message.
        differentiable
            What the *engine* offers. Defaults to what this problem can supply;
            pass it explicitly to ask "could emcee run this?" of a problem whose
            backend happens to be differentiable.
        """
        gradient = self.differentiable if differentiable is None else bool(differentiable)
        for dataset in self.datasets.values():
            dataset.check_engine(differentiable=gradient, engine=engine)

    # -- failure signalling ---------------------------------------------------

    @property
    def failures(self) -> tuple[Failure, ...]:
        """The most recent recorded failures, oldest first. Bounded."""
        return tuple(self._failures)

    @property
    def failure_counts(self) -> Mapping[FailureReason, int]:
        """How many times each reason has occurred, unbounded.

        What an engine driver should report: "8 214 of 200 000 proposals could
        not be scored, all ``likelihood_failed``" is actionable, whereas the
        same information as ``-inf`` values is invisible.
        """
        return types.MappingProxyType(dict(self._failure_counts))

    @property
    def last_failure(self) -> Failure | None:
        """The most recent failure, or ``None``."""
        return self._failures[-1] if self._failures else None

    def failure_summary(self) -> str:
        """One aggregated sentence per failure class — what a driver should print.

        Peter's §17 Q1 ruling asks the recording half to be "aggregated so high
        failure rates do not spam". A run that could not score 8 214 of 200 000
        proposals should say so once, with the range of parameter values over
        which it happened, not emit 8 214 warnings — and the parameter range is
        the part that tells a user *which* prior is too wide.

        Returns an empty string when nothing has failed, so a driver can write
        ``if summary := problem.failure_summary(): warn(summary)``.
        """
        if not self._failure_counts:
            return ""
        lines: list[str] = []
        for reason, count in self._failure_counts.most_common():
            seen = [f for f in self._failures if f.reason is reason]
            where = sorted({f.where for f in seen if f.where})
            ranges: dict[str, tuple[float, float]] = {}
            for failure in seen:
                for name, value in failure.values.items():
                    low, high = ranges.get(name, (value, value))
                    ranges[name] = (min(low, value), max(high, value))

            detail = ""
            if ranges:
                detail = "; over " + ", ".join(
                    f"{name} in [{low:.6g}, {high:.6g}]" if low != high else f"{name}={low:.6g}"
                    for name, (low, high) in sorted(ranges.items())
                )
            located = f" in {where}" if where else ""
            sampled = "" if count <= len(seen) else f" (ranges from the last {len(seen)})"
            lines.append(f"{count} draw(s) failed: {reason}{located}{detail}{sampled}")
        lines.append("Re-run with FittingProblem(..., strict=True) to raise at the offending draw.")
        return "\n".join(lines)

    def reset_failures(self) -> None:
        """Forget the failure history and counts — call between runs."""
        self._failures.clear()
        self._failure_counts.clear()

    def _record(self, failure: Failure | None) -> None:
        if failure is None:
            return
        self._failures.append(failure)
        self._failure_counts[failure.reason] += 1

    # -- RNG policy -----------------------------------------------------------

    def rng(self, label: str) -> np.random.Generator:
        """The named sub-stream of this problem's seed.

        ``lowering.md`` §9.2's policy, realised on the reference path: one
        integer seed per run, named sub-streams derived from it by
        :func:`~ampere.core.rng.substream`, so that prior sampling, sampler
        initialisation and an SBI simulation budget never make each other
        irreproducible. The generator for a label is created once and then
        advanced, so repeated draws differ (a simulation budget must not be one
        point drawn 10⁴ times) while the whole sequence is reproducible for a
        given seed and call order.

        With ``seed=None`` every stream is entropy-seeded and nothing is
        reproducible — which is the honest behaviour for a run that did not ask
        to be.
        """
        stream = self._streams.get(label)
        if stream is None:
            stream = np.random.default_rng() if self.seed is None else _generator(self.seed, label)
            self._streams[label] = stream
        return stream

    # -- simulation (SBI) -----------------------------------------------------

    def simulate(
        self,
        values: Mapping[str, Value] | ArrayLike | None = None,
        *,
        observe: bool = False,
        rng: np.random.Generator | None = None,
        stream: str = "simulate",
    ) -> Simulation:
        """Run the forward model once: ``(θ, x)``, or a flagged failure.

        ``DEVELOPMENT_PLAN.md`` §4.5's ``simulate(params) -> data``. Failures are
        **flagged, not raised** (:attr:`Simulation.failed`), so a simulation
        budget survives a simulator that crashes on 2 % of its prior draws.

        Parameters
        ----------
        values
            The θ to simulate at. ``None`` draws one from the joint prior on the
            *stream* sub-stream — the SBI budget idiom, and the reason this
            default differs from :meth:`evaluate`'s (which uses the reference θ).
        observe
            Also draw noisy observations. See :meth:`draw_observation` for what
            is and is not implementable from the merged contracts.
        rng
            Override the sub-stream.
        stream
            Sub-stream name; see :meth:`rng`.
        """
        generator = self.rng(stream) if rng is None else rng
        resolved = (
            dict(self._mapping.merged.sample(generator))
            if values is None
            else self._resolve(values)
        )
        theta = self._mapping.merged.pack(resolved)
        routed = self._mapping.distribute(resolved)

        # Each stage is trapped separately, so a failure's recorded reason and
        # location are as specific here as they are in log_prob. Trapping the
        # model and the instrument chain together would report every chain
        # failure as MODEL_FAILED with no dataset named, which is exactly the
        # unhelpful half of "-inf with a recorded reason".
        results: dict[str, ModelResult] = {}
        predicted: dict[str, FunctionSamples] = {}
        try:
            results = self._evaluate_models(routed)
        except self._failure_types as error:
            failure = _failure_from(_reason_for(error), error, _where(error))
            self._record(failure)
            return Simulation(parameters=resolved, theta=theta, failure=failure)

        for label, dataset in self.datasets.items():
            try:
                predicted[label] = dataset.predict(
                    results[self._bindings[label]], routed.get(label)
                )
            except self._failure_types as error:
                failure = _failure_from(FailureReason.INSTRUMENT_FAILED, error, label)
                self._record(failure)
                return Simulation(
                    parameters=resolved, theta=theta, results=results, failure=failure
                )

        observations: dict[str, FunctionSamples] | None = None
        if observe:
            observations = {}
            for label, dataset in self.datasets.items():
                try:
                    observations[label] = dataset.draw_observation(
                        predicted[label], routed.get(label), generator
                    )
                except self._failure_types as error:
                    failure = _failure_from(FailureReason.LIKELIHOOD_FAILED, error, label)
                    self._record(failure)
                    return Simulation(
                        parameters=resolved,
                        theta=theta,
                        results=results,
                        predicted=predicted,
                        failure=failure,
                    )

        return Simulation(
            parameters=resolved,
            theta=theta,
            results=results,
            predicted=predicted,
            observations=observations,
        )

    # -- composition-time validation ------------------------------------------

    def validate(self) -> None:
        """Run the composition-time checks: step (5) of the lifecycle.

        Evaluates each compiled model once at :attr:`reference_values`, pushes
        the result through every instrument, and calls
        :meth:`Dataset.check_alignment`, which is where a shape mismatch, a unit
        mismatch, a censoring declaration the family cannot consume, an
        unimplemented latent combination and a mis-sized latent block are all
        refused. Idempotent; called automatically unless ``validate=False``.
        """
        if self._validated:
            return
        routed = self._mapping.distribute(dict(self.reference_values))
        results = self._evaluate_models(routed)
        self._run_alignment_checks(results, routed)

    def _run_alignment_checks(
        self,
        results: Mapping[str, ModelResult],
        routed: Mapping[str, Mapping[str, Value]],
    ) -> None:
        for label, dataset in self.datasets.items():
            predicted = dataset.predict(results[self._bindings[label]], routed.get(label))
            dataset.check_alignment(predicted)
        self._validated = True

    @property
    def validated(self) -> bool:
        """Whether the composition-time checks have run."""
        return self._validated

    # -- internals ------------------------------------------------------------

    def _resolve(self, values: Mapping[str, Value] | ArrayLike | None) -> dict[str, Value]:
        if values is None:
            return dict(self.reference_values)
        return dict(self._mapping.merged.complete(values))

    def _evaluate_models(self, routed: Mapping[str, Mapping[str, Value]]) -> dict[str, ModelResult]:
        """Evaluate each model **once** per θ, however many datasets read it.

        Load-bearing, not an optimisation. The commonest joint fit in ampere's
        target scope is one model producing several channels — a low-resolution
        SED and high-resolution windows, the RA and Dec time series of one
        reflex orbit — read by several datasets. Evaluating per dataset would
        multiply the cost of the expensive half of the loop by the number of
        datasets, and would also make a stochastic model inconsistent with
        itself inside a single ``log_prob``: two datasets would be scored
        against two *different* realisations of the same parameters. Each
        dataset takes its own channel out of the one shared result through
        :meth:`~ampere.core.results_schema.ModelResult.require`, via
        :meth:`~ampere.core.transform.Instrument.bind`.
        """
        results: dict[str, ModelResult] = {}
        for label, instance in self._compiled.items():
            try:
                results[label] = instance(dict(routed.get(label, {})))
            except Exception as error:
                # Tag and re-raise unchanged: FittingProblem decides whether an
                # exception is an unscoreable point or a bug by catching a
                # *declared* tuple of types, so wrapping a user simulator's own
                # exception class here would defeat that test — and the original
                # traceback is what a genuine bug needs.
                _tag(error, FailureReason.MODEL_FAILED, label)
                raise
        return results

    def _joint_log_likelihood(
        self, routed: Mapping[str, Mapping[str, Value]]
    ) -> tuple[float, dict[str, float], Failure | None]:
        """The sum over datasets, with §4.5's failure signalling wrapped round it."""
        try:
            results = self._evaluate_models(routed)
        except self._failure_types as error:
            where = _where(error)
            return (
                -math.inf,
                {},
                _failure_from(_reason_for(error), error, where, routed.get(where)),
            )

        # Deferred validation runs *inside* this loop, on the predicted
        # container the loop computes anyway. Running it beforehand, as an
        # earlier version did, cost two chain evaluations on the first call —
        # the very cost validate=False exists to avoid — and, worse, put
        # dataset.predict outside the try below, so a declared simulator
        # failure escaped log_prob as an exception on the first evaluation and
        # returned -inf on every later one. Whether a run died then depended on
        # whether the sampler's first proposal happened to be scoreable.
        validating = not self._validated
        contributions: dict[str, float] = {}
        total = 0.0
        for label, dataset in self.datasets.items():
            values = routed.get(label)
            split = dataset.route(values) if values is not None else None
            try:
                predicted = dataset.predict(results[self._bindings[label]], routed=split)
            except self._failure_types as error:
                return (
                    -math.inf,
                    contributions,
                    _failure_from(
                        FailureReason.INSTRUMENT_FAILED,
                        error,
                        label,
                        (split or {}).get(INSTRUMENT_COMPONENT),
                    ),
                )
            if validating:
                # Raises for a genuine composition error — a shape or unit
                # mismatch is a broken problem, not an unscoreable point — and
                # fixes the retained-sample count this dataset is held to.
                dataset.check_alignment(predicted)
            try:
                contribution = dataset.log_likelihood_of(predicted, routed=split)
            except self._failure_types as error:
                reason = _classify_likelihood_failure(dataset, predicted)
                return (
                    -math.inf,
                    contributions,
                    _failure_from(reason, error, label, (split or {}).get(LIKELIHOOD_COMPONENT)),
                )
            if math.isnan(contribution) or contribution == math.inf:
                return (
                    -math.inf,
                    contributions,
                    Failure(
                        FailureReason.NON_FINITE_LOG_LIKELIHOOD,
                        f"the likelihood returned {contribution}, which is neither a density nor "
                        f"an impossibility. NaN usually means a prediction that is NaN on a "
                        f"retained sample; +inf means a degenerate, infinitely peaked density, "
                        f"and both would otherwise be stored per draw as an Evaluation whose "
                        f"log_likelihood contradicts its log_prob.",
                        where=label,
                    ),
                )
            contributions[label] = contribution
            total += contribution
        if validating:
            self._validated = True
        return total, contributions, None

    def __repr__(self) -> str:
        models = list(self._compiled)
        return (
            f"<FittingProblem {models} -> {list(self.datasets)}, "
            f"{self.free_size} free dimension(s), {self.capabilities}>"
        )


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------


def _normalise_models(model: Model | Mapping[str, Model], model_label: str) -> dict[str, Model]:
    if isinstance(model, Model):
        return {_check_label(model_label, "model label"): model}
    if not isinstance(model, Mapping):
        raise DatasetError(
            f"a FittingProblem takes a Model, or a mapping of label to Model for a multi-object "
            f"fit, got {type(model).__name__}."
        )
    if not model:
        raise DatasetError("a FittingProblem needs at least one model.")
    normalised: dict[str, Model] = {}
    for label, instance in model.items():
        if not isinstance(instance, Model):
            raise DatasetError(
                f"model {label!r} is {type(instance).__name__}, not a Model. Subclass Model and "
                f"implement evaluate()."
            )
        normalised[_check_label(label, "model label")] = instance
    return normalised


class _TaggedFailure(Exception):
    """Carries a reason and a location from where a failure happened to where it is caught."""

    def __init__(self, reason: FailureReason, where: str) -> None:
        super().__init__(reason, where)
        self.reason = reason
        self.where = where


def _tag(error: BaseException, reason: FailureReason, where: str) -> None:
    """Attach a reason and a location to an exception without changing its type.

    Best-effort. ``BaseException`` always provides a ``__dict__``, so a
    ``__slots__`` subclass is fine, but a class overriding ``__setattr__`` — a
    frozen dataclass or ``attrs`` exception, which a wrapped external code may
    well use — refuses the assignment. Losing the tag costs a less specific
    failure reason; letting the ``AttributeError`` out of the ``except`` block
    that is trying to *handle* a failure would replace an unscoreable point with
    a crash, which is much worse.
    """
    try:
        error._ampere_failure = _TaggedFailure(reason, where)  # type: ignore[attr-defined]
    except Exception:  # an untaggable exception is still a failure
        pass


def _reason_for(error: BaseException) -> FailureReason:
    tagged = getattr(error, "_ampere_failure", None)
    return tagged.reason if isinstance(tagged, _TaggedFailure) else FailureReason.MODEL_FAILED


def _where(error: BaseException) -> str:
    tagged = getattr(error, "_ampere_failure", None)
    return tagged.where if isinstance(tagged, _TaggedFailure) else ""


def _classify_likelihood_failure(dataset: Dataset, predicted: FunctionSamples) -> FailureReason:
    """Distinguish "the simulator returned NaNs" from "the covariance would not factorise".

    ``DEVELOPMENT_PLAN.md`` §4.5 names both, and the remedies are different: one
    is the model's fault and the other the kernel's. The O(N) scan happens only
    on the cold path — a point that has already failed — so the hot loop pays
    nothing for the distinction.
    """
    weights = (
        np.asarray(dataset.observed.weights()).ravel() * np.asarray(predicted.weights()).ravel()
    )
    retained = np.asarray(predicted.values).ravel()[weights > 0.0]
    if retained.size and not bool(np.all(np.isfinite(retained))):
        return FailureReason.NON_FINITE_PREDICTION
    return FailureReason.LIKELIHOOD_FAILED
