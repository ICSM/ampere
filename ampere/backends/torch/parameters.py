"""``lowering.md`` §5 to §9 for torch: a declaration becomes an ``nn.Module`` tree.

``ampere.core.ParameterSet`` is the neutral declaration; this module is what a
torch-backed run actually evaluates. Everything here is a row of
``lowering.md``:

§1.6 the entry assertion
    Checked once, loudly, at construction: the set is resolved (no deferred
    parameters), ``merge`` has happened (no tie labels left), and every prior
    is describable. Three preconditions failing in one place beats three subtly
    different errors from three backends halfway through construction.

§5 the declaration-form table
    A **free** parameter becomes ``nn.Parameter(tensor, requires_grad=True)``
    via ``register_parameter``; a **fixed** one becomes a buffer, via
    ``register_buffer(..., persistent=True)``. An **array-valued** parameter is
    one tensor of the declared shape whose prior is i.i.d. across elements and
    summed. A **plate** lowers identically to an array-valued parameter —
    torch has no plate primitive, so the plate structure survives in the names
    and the provenance rather than in a framework object (§5.1's accepted
    asymmetry). **Units lower to nothing** (§1.4).

§6.1 the module convention
    One ``nn.Module`` per component, nested, so that ``state_dict`` keys come
    out as the merged qualified names (``spectrum.temperature``) without
    anything having to substitute a separator — torch splits ``state_dict``
    keys on ``.``, so a flat registration of a dotted name is not even legal.
    A tied parameter's merged name is deliberately unqualified
    (``parameters.md`` §14.3), so it lands on the root module, which is §12.3's
    settled answer.

§8 hierarchical priors
    Built per evaluation from the referenced parameters' current tensors, in
    the declaration's own topological order, never hoisted out of the loop.

§9 RNG
    Draws go through ``substream`` (:mod:`ampere.backends.torch.rng`), never
    torch's global generator.

§10.1 dtype
    Every tensor is created with an explicit ``dtype=``/``device=``;
    ``torch.set_default_dtype`` is never called.

The reference implementation is the oracle
------------------------------------------
``lowering.md`` §2(c) is emphatic that this class must *reproduce*
``ParameterSet.lnprior_unconstrained``, not merely produce "a valid
unconstrained density": a density differing by a non-constant function of ``y``
is a different posterior. So the arithmetic here is deliberately the same
composition — ``lnprior(constrain(y)) + Σ log|d constrain / dy|``, with the
Jacobian evaluated at the **unconstrained** point — computed with torch's
objects rather than scipy's. ``tests/conformance/test_parameters.py`` holds it
to scipy, and ``test_cross_backend.py`` holds it to the reference backend.

The two-argument Jacobian, once, in one place
----------------------------------------------
torch's ``log_abs_det_jacobian(x, y)`` takes ``x`` = the transform's *input* =
the **unconstrained** point and ``y`` = its output = the constrained one, and
returns ``log|dy/dx|``. Ampere's convention is the opposite naming (``x``
constrained, ``y`` unconstrained, one argument), so copying a call out of the
torch documentation into ampere-named variables transposes the arguments while
looking correct in both places (``lowering.md`` §2(b)). The call is made in
exactly one place below — :meth:`TorchParameterSpace.log_prior_unconstrained` —
and ``tests/backends/test_torch_lowering.py`` pins it at a point where the two
arguments differ, so a swap fails rather than producing a plausible wrong
posterior.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np
import torch
from torch import nn
from torch.distributions import transforms

from ampere.core.exceptions import LoweringError, ParameterError
from ampere.core.parameter import (
    HierarchicalPrior,
    Parameter,
    ParameterSet,
    describe_prior,
)

from ._config import BACKEND, DEFAULT_DEVICE, DEFAULT_DTYPE, as_tensor, to_numpy
from .lowering import (
    LoweredPrior,
    lower_bijection,
    lower_hierarchical,
    lower_prior,
    warn_icdf_fallback,
)

__all__ = ["LoweredParameters", "TorchParameterSpace"]


def _finite_or_minus_infinity(total: torch.Tensor) -> torch.Tensor:
    """*total*, or ``-inf`` where it is not finite, through :func:`torch.where`.

    The branch-free form of ``ParameterSet.lnprior``'s ``math.isfinite`` short
    circuit (W2.4 slice 2). Three properties matter and all three are the
    reason it is a ``where`` rather than an ``if``: it keeps the result
    attached to the autograd graph, which a fresh ``-inf`` constant would not
    and which pyro's NUTS requires of its potential; it is not data-dependent
    control flow, so ``torch.func.vmap`` can batch straight through it; and it
    collapses ``nan`` to ``-inf`` as well as ``+inf``, which is what the short
    circuit did (a prior contribution is never ``+inf``, so a ``nan`` can only
    have come from a ``-inf`` meeting something pathological, and "impossible"
    is the honest reading).
    """
    return torch.where(torch.isfinite(total), total, torch.full_like(total, -math.inf))


class LoweredParameters(nn.Module):
    """The ``nn.Module`` a :class:`~ampere.core.ParameterSet` lowers to (§6.1).

    A plain module with no ``forward``: it exists to *hold* the declaration's
    tensors under torch's own ownership rules, so that ``.to(dtype)``,
    ``.to(device)``, ``state_dict``/``load_state_dict`` and
    ``parameters()``/``buffers()`` all do the right thing without ampere
    reimplementing any of it.

    The nesting is the point. A merged name such as ``"spectrum.temperature"``
    is not a legal torch attribute name — torch splits ``state_dict`` keys on
    ``.`` — so lowering must either substitute the separator or nest real
    submodules per component. ``lowering.md`` §6.1 recommends nesting, and this
    does: ``spectrum`` is a child module and ``temperature`` a parameter on it,
    so the ``state_dict`` key *is* the merged name and torch's own device/dtype
    recursion works. The flat name and the module path staying identical is
    what keeps the ArviZ coordinate, the provenance record and the checkpoint
    key one string rather than three.

    Fixed parameters are buffers (``persistent=True``), not parameters, per
    §5.2: they and free parameters answer the same two questions — does the
    optimiser take a gradient with respect to it, and does it move with device
    and dtype — even though they are different *kinds* of thing at the contract
    level. ``persistent=True`` because a fixed parameter's value is part of the
    fit's definition and belongs in a checkpoint.
    """

    def child(self, name: str) -> LoweredParameters:
        """The named child component, created on first use."""
        existing = self._modules.get(name)
        if isinstance(existing, LoweredParameters):
            return existing
        if existing is not None or hasattr(type(self), name):
            raise LoweringError(
                name,
                backend=BACKEND,
                parameter=name,
                detail=(
                    f"component name {name!r} collides with an attribute of torch.nn.Module, so "
                    f"it cannot be a submodule; rename the component"
                ),
            )
        created = LoweredParameters()
        self.add_module(name, created)
        return created


def _place(root: LoweredParameters, name: str) -> tuple[LoweredParameters, str]:
    """The module a merged *name* belongs on, and the leaf name it takes there.

    An unqualified name (a tie label — ``parameters.md`` §14.3 keeps those
    global and unqualified) belongs to no single component and lands on the
    root, which is ``lowering.md`` §12.3's settled answer.
    """
    *path, leaf = name.split(".")
    owner = root
    for part in path:
        owner = owner.child(part)
    return owner, leaf


def _initial_value(parameter: Parameter) -> Any:
    """A sensible starting tensor value for a lowered parameter.

    The arithmetic below never reads it — every method takes θ as an argument —
    but a module whose parameters were all zero would make a ``state_dict``
    dump and any torch-side optimiser start somewhere arbitrary. So: the
    declared value where there is one, the prior's median where the family can
    give one, and zero otherwise (a hierarchical prior has no median until its
    hyperparameters do).
    """
    if parameter.value is not None:
        return parameter.value
    prior = parameter.prior
    if prior is not None and not isinstance(prior, HierarchicalPrior):
        try:
            return float(prior.median())
        except Exception:  # pragma: no cover - a described prior always has one
            return 0.0
    return 0.0


def _topological_order(declaration: ParameterSet) -> tuple[str, ...]:
    """Declaration order, with every hierarchical prior after what it references.

    ``lowering.md`` §8: "Emission order is ``ParameterSet``'s topological
    order, not declaration order. A hyperparameter is always sampled before
    anything referencing it." The core already computes exactly this order for
    its own evaluation; it is private there, so it is recomputed here from the
    public ``Parameter.references`` rather than reaching across the contract
    boundary for ``_order``.
    """
    names = list(declaration.names)
    remaining = set(names)
    ordered: list[str] = []
    while remaining:
        progressed = False
        for name in names:
            if name not in remaining:
                continue
            references = set(declaration[name].references) & set(names)
            if references & remaining:
                continue
            ordered.append(name)
            remaining.discard(name)
            progressed = True
        if not progressed:  # pragma: no cover - core rejects cycles at construction
            raise ParameterError(
                f"the declaration has a cyclic hierarchical reference among {sorted(remaining)}"
            )
    return tuple(ordered)


class TorchParameterSpace:
    """A merged :class:`~ampere.core.ParameterSet`, lowered onto torch.

    Satisfies ``tests/conformance/protocol.py``'s ``ParameterSpace``: every
    method takes and returns plain numpy on the boundary, and computes in torch
    in between. The tensor-valued twins (:meth:`log_prior_tensor`,
    :meth:`constrain_tensor`, :meth:`log_prior_unconstrained_tensor`) are the
    same arithmetic without the conversion, and they are what a gradient-based
    engine consumes: they keep the autograd graph, so
    ``torch.autograd.grad(space.log_prior_unconstrained_tensor(y), y)`` is the
    prior's score.

    Parameters
    ----------
    declaration
        The **merged** set (``lowering.md`` §1.1) — ``ParameterMapping.merged``,
        not a component's own set.
    strict
        The problem's own ``strict`` flag, threaded in from
        ``FittingProblem(strict=...)``. It governs one thing here: whether a
        prior family with no native torch ``icdf`` may fall back to the
        reference path for :meth:`prior_transform` (``lowering.md`` §3.6). One
        flag, one meaning — not a new per-lowering knob.
    dtype, device
        Threaded explicitly into every tensor (``lowering.md`` §10.1).

    Raises
    ------
    TyingError, ParameterError
        From the §1.6 entry assertion, if an unmerged or unresolved set arrives.
    LoweringError
        If any declared family or bijection has no torch row, or if ``strict``
        refuses the ``icdf`` fallback.
    """

    def __init__(
        self,
        declaration: ParameterSet,
        *,
        strict: bool = False,
        dtype: torch.dtype = DEFAULT_DTYPE,
        device: torch.device = DEFAULT_DEVICE,
    ) -> None:
        self._check_preconditions(declaration)
        self._declaration = declaration
        self._strict = bool(strict)
        self._dtype = dtype
        self._device = device
        self._order = _topological_order(declaration)

        self.module = LoweredParameters()
        self._slices: dict[str, slice] = {}
        self._priors: dict[str, LoweredPrior] = {}
        self._hierarchical: dict[str, HierarchicalPrior] = {}
        self._transforms: dict[str, transforms.Transform] = {}
        self._fixed: dict[str, torch.Tensor] = {}
        self._free: list[str] = []

        for parameter in declaration.parameters:
            owner, leaf = _place(self.module, parameter.name)
            tensor = as_tensor(_initial_value(parameter), dtype=dtype, device=device)
            if parameter.shape:
                tensor = tensor.expand(parameter.shape).clone()
            if parameter.is_fixed:
                # §5.2: a fixed parameter and a buffer answer the same two
                # runtime questions, so they take the same mechanism.
                owner.register_buffer(leaf, tensor, persistent=True)
                self._fixed[parameter.name] = tensor
                continue
            owner.register_parameter(leaf, nn.Parameter(tensor, requires_grad=True))
            self._free.append(parameter.name)
            self._slices[parameter.name] = declaration.free_slice(parameter.name)
            self._lower_prior_of(parameter)

        self._check_free_layout()
        self._warn_about_icdf()

    # -- construction helpers -------------------------------------------------

    @staticmethod
    def _check_preconditions(declaration: ParameterSet) -> None:
        """``lowering.md`` §1.6: the three standing preconditions, once, at entry.

        A lowering "may assume them and should assert them once at entry rather
        than defensively re-checking per parameter". Failing loudly here is
        worth more than three subtly different errors raised halfway through
        construction.
        """
        if not declaration.is_resolved:
            raise ParameterError(
                f"cannot lower onto the {BACKEND!r} backend: parameter(s) "
                f"{list(declaration.deferred_names)} are deferred (shared_as=... with no prior "
                f"of their own). A merged set contains only free and fixed parameters "
                f"(lowering.md §1.3); merge with ParameterSet.merge first."
            )
        unmerged = [p.name for p in declaration.parameters if p.shared_as is not None]
        if unmerged:
            raise ParameterError(
                f"cannot lower onto the {BACKEND!r} backend: parameter(s) {unmerged} still carry "
                f"a tie label. A tie is a compile-time collapse, not a runtime object "
                f"(lowering.md §1.2), so no backend should ever see one — this set has not been "
                f"through ParameterSet.merge."
            )
        declaration.to_spec()  # §1.5: refuses an opaque prior, naming the parameter.

    def _lower_prior_of(self, parameter: Parameter) -> None:
        """Lower one free parameter's prior and its unconstraining bijection."""
        prior = parameter.prior
        if isinstance(prior, HierarchicalPrior):
            # §8: nothing to build until the hyperparameters have values. The
            # bijection is still fixed at lowering time, because
            # ``unconstraining_bijection`` refuses a hierarchical family whose
            # support would move from draw to draw — so it is safe to lower it
            # once, and unsafe to rebuild it per evaluation from a distribution
            # whose support depends on the current hyperparameters.
            self._hierarchical[parameter.name] = prior
            self._transforms[parameter.name] = lower_bijection(
                parameter.unconstraining_bijection(),
                _PLACEHOLDER,
                parameter=parameter.name,
            )
            return
        lowered = lower_prior(
            describe_prior(prior),
            parameter=parameter.name,
            dtype=self._dtype,
            device=self._device,
        )
        self._priors[parameter.name] = lowered
        self._transforms[parameter.name] = lower_bijection(
            parameter.bijection, lowered, parameter=parameter.name
        )

    def _check_free_layout(self) -> None:
        """``lowering.md`` §11 row 9: ``free_size`` agrees with the declaration.

        Cheap, and it is the assertion that would catch a lowering that dropped
        a parameter, double-counted an array one, or gave a fixed parameter a
        sampler dimension.
        """
        lowered = sum(
            parameter.numel() for parameter in self.module.parameters() if parameter.requires_grad
        )
        if lowered != self._declaration.free_size:
            raise LoweringError(  # pragma: no cover - a guard against future edits
                "free_size",
                backend=BACKEND,
                detail=(
                    f"the lowered module has {lowered} trainable element(s) but the declaration "
                    f"has {self._declaration.free_size} free dimension(s)"
                ),
            )

    def _warn_about_icdf(self) -> None:
        """``lowering.md`` §3.6's contract, applied once at lowering time.

        A hierarchical prior has no lowered object until its hyperparameters
        have values, so its ``icdf`` availability is established by lowering
        the *family* once with placeholder arguments. That is sound because
        ``icdf`` is implemented (or not) by a ``torch.distributions`` class,
        not by a particular parametrisation of it — and it is necessary,
        because the alternative is deciding at the first ``prior_transform``
        call, which is after the point §3.6 says the decision is made.
        """
        missing = {self._family_of(name) for name in self._free if not self._has_native_icdf(name)}
        warn_icdf_fallback(missing, strict=self._strict, where="prior_transform")
        self._icdf_fallback = frozenset(missing)

    def _has_native_icdf(self, name: str) -> bool:
        lowered = self._priors.get(name)
        if lowered is not None:
            return lowered.has_icdf
        placeholder = {
            reference: as_tensor(1.0, dtype=self._dtype, device=self._device)
            for reference in self._hierarchical[name].references
        }
        try:
            return lower_hierarchical(
                self._hierarchical[name],
                placeholder,
                parameter=name,
                dtype=self._dtype,
                device=self._device,
            ).has_icdf
        except LoweringError:  # pragma: no cover - the family is refused at lowering
            return False

    def _family_of(self, name: str) -> str:
        prior = self._declaration[name].prior
        if isinstance(prior, HierarchicalPrior):
            return prior.family
        return describe_prior(prior).family

    # -- the ParameterSpace protocol -----------------------------------------

    @property
    def declaration(self) -> ParameterSet:
        """The backend-neutral declaration this space realises."""
        return self._declaration

    @property
    def free_size(self) -> int:
        """Number of flat free dimensions — the engine's dimension."""
        return self._declaration.free_size

    @property
    def families(self) -> frozenset[str]:
        """The neutral prior-family names this lowering consulted the registry for."""
        return frozenset(self._family_of(name) for name in self._free)

    def lowering_provenance(self) -> list[dict[str, Any]]:
        """The **user-registered** rows this lowering consulted, netCDF-safe.

        ``lowering.md`` §12.8's hardening, ready for
        ``ampere.results.provenance.provenance_attrs``'s ``extra=``. Empty for
        a lowering that used only ampere's own table, which is the point: the
        signal a reviewer wants is "did this run depend on something outside
        the conformance suite's guarantees?", and stamping ten built-in rows
        would bury it. W2.6 deferred the first-class provenance key until a
        backend drove lowering end to end; W2.13 is where that happens, and
        ``ampere.inference``'s drivers read this through the realisation.
        """
        from ampere.core.lowering import provenance_entries

        from .lowering import consulted_resolutions

        return provenance_entries(consulted_resolutions(set(self.families)))

    @property
    def icdf_fallback_families(self) -> frozenset[str]:
        """Families whose ``prior_transform`` runs on the reference path (§3.6)."""
        return self._icdf_fallback

    def free_labels(self) -> tuple[str, ...]:
        """One label per flat free dimension, array elements included.

        Taken from the declaration: a label is composition-time metadata, not
        arithmetic, and inventing a second naming scheme at the backend
        boundary would put a translation layer between the sampler's output and
        the results schema — the bug class ``lowering.md`` §8 names.
        """
        return self._declaration.free_labels()

    def pack(self, values: Mapping[str, Any]) -> np.ndarray:
        """Flatten a name-to-value mapping into the free vector."""
        return to_numpy(self.pack_tensor(values))

    def unpack(self, theta: Any) -> dict[str, Any]:
        """Expand a free vector into every parameter, fixed ones included.

        A scalar parameter comes back as a Python ``float`` and an
        array-valued one as a numpy array, which is what
        ``ParameterSet.unpack`` returns and therefore what a model's
        ``evaluate`` is written against.
        """
        values: dict[str, Any] = {}
        for name, tensor in self.unpack_tensor(theta).items():
            array = to_numpy(tensor)
            values[name] = array if self._declaration[name].shape else float(array)
        return values

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        """Map ``[0, 1]^n`` to the constrained free vector (nested sampling).

        Native ``icdf`` where torch has one; the reference (scipy) ``ppf`` for
        the families where it does not, which is ``lowering.md`` §3.6's
        sanctioned fallback: ``prior_transform`` has one mathematical
        definition, so computing it in numpy changes nothing about the
        posterior. It is warned about once, at lowering time, and refused
        outright under ``strict`` — both in :func:`~.lowering.warn_icdf_fallback`.
        """
        cube = np.asarray(unit_cube, dtype=float).reshape(-1)
        if cube.size != self.free_size:
            raise ParameterError(
                f"expected a unit-cube vector of length {self.free_size}, got {cube.size}"
            )
        flat = np.empty(self.free_size, dtype=float)
        resolved: dict[str, torch.Tensor] = dict(self._fixed)
        for name in self._order:
            parameter = self._declaration[name]
            if parameter.is_fixed:
                continue
            chunk = cube[self._slices[name]]
            # ``distribution_of`` is the same route the density takes: constant
            # for an ordinary family, rebuilt from the hyperparameters resolved
            # so far for a hierarchical one (§8's topological order is what
            # makes that well defined here).
            lowered = self.distribution_of(name, resolved)
            if lowered.has_icdf:
                drawn = to_numpy(
                    lowered.icdf(as_tensor(chunk, dtype=self._dtype, device=self._device))
                )
            elif name in self._hierarchical:
                bound = self._hierarchical[name].bind(
                    {key: to_numpy(value) for key, value in resolved.items()}
                )
                drawn = np.asarray(bound.ppf(chunk), dtype=float)
            else:
                drawn = np.asarray(parameter.prior.ppf(chunk), dtype=float)
            flat[self._slices[name]] = drawn.reshape(-1)
            resolved[name] = as_tensor(
                drawn.reshape(parameter.shape) if parameter.shape else float(drawn.reshape(-1)[0]),
                dtype=self._dtype,
                device=self._device,
            )
        return flat

    def lnprior(self, values: Any) -> float:
        """Log prior density in the constrained parameterisation."""
        return float(self.log_prior_tensor(values))

    def constrain(self, unconstrained: Any) -> np.ndarray:
        """Map the unconstrained vector into the priors' support."""
        return to_numpy(self.constrain_tensor(unconstrained))

    def unconstrain(self, values: Any) -> np.ndarray:
        """Map values in the priors' support to the unconstrained vector."""
        return to_numpy(self.unconstrain_tensor(values))

    def lnprior_unconstrained(self, unconstrained: Any) -> float:
        """Log prior density in unconstrained space, Jacobian term included."""
        return float(self.log_prior_unconstrained_tensor(unconstrained))

    # -- the tensor-valued twins ---------------------------------------------

    def _tensor(self, values: Any) -> torch.Tensor:
        vector = as_tensor(values, dtype=self._dtype, device=self._device).reshape(-1)
        if vector.numel() != self.free_size:
            raise ParameterError(
                f"expected a flat vector of length {self.free_size} "
                f"({list(self._declaration.free_names)}), got length {vector.numel()}"
            )
        return vector

    def pack_tensor(self, values: Mapping[str, Any]) -> torch.Tensor:
        """:meth:`pack`, keeping tensors (and any gradient) intact."""
        unknown = set(values) - set(self._declaration.names)
        if unknown:
            raise ParameterError(
                f"pack() got value(s) for unknown parameter(s) {sorted(unknown)}; this set "
                f"declares {list(self._declaration.names)}."
            )
        chunks: list[torch.Tensor] = []
        for name in self._free:
            try:
                raw = values[name]
            except KeyError:
                raise ParameterError(
                    f"pack() is missing a value for free parameter {name!r}"
                ) from None
            chunk = as_tensor(raw, dtype=self._dtype, device=self._device).reshape(-1)
            expected = self._declaration[name].size
            if chunk.numel() != expected:
                raise ParameterError(
                    f"parameter {name!r} expects {expected} value(s) "
                    f"(shape {self._declaration[name].shape}), got {chunk.numel()}"
                )
            chunks.append(chunk)
        if not chunks:
            return as_tensor([], dtype=self._dtype, device=self._device)
        return torch.cat(chunks)

    def unpack_tensor(self, theta: Any) -> dict[str, torch.Tensor]:
        """:meth:`unpack`, keeping tensors (and any gradient) intact.

        Fixed parameters are injected at their declared values, so a model
        always receives its full vocabulary (``parameters.md`` §7's promise
        that fixing a parameter changes which lowering row it takes and nothing
        else).
        """
        vector = self._tensor(theta)
        resolved: dict[str, torch.Tensor] = {}
        for parameter in self._declaration.parameters:
            if parameter.is_fixed:
                resolved[parameter.name] = self._fixed[parameter.name]
                continue
            chunk = vector[self._slices[parameter.name]]
            resolved[parameter.name] = (
                chunk.reshape(parameter.shape) if parameter.shape else chunk[0]
            )
        return resolved

    def _resolved(self, values: Any) -> dict[str, torch.Tensor]:
        if isinstance(values, Mapping):
            filled = dict(self._fixed)
            unknown = set(values) - set(self._declaration.names)
            if unknown:
                raise ParameterError(
                    f"got value(s) for unknown parameter(s) {sorted(unknown)}; this set "
                    f"declares {list(self._declaration.names)}."
                )
            for name, raw in values.items():
                filled[name] = as_tensor(raw, dtype=self._dtype, device=self._device)
            missing = [name for name in self._declaration.names if name not in filled]
            if missing:
                raise ParameterError(f"no value supplied for parameter(s) {missing}")
            return filled
        return self.unpack_tensor(values)

    def distribution_of(self, name: str, resolved: Mapping[str, torch.Tensor]) -> LoweredPrior:
        """The lowered prior for *name*, given the values it may reference.

        Constant for an ordinary family — lowered once, at construction — and
        rebuilt from *resolved* for a hierarchical one, which is
        ``lowering.md`` §8's rule that a hierarchical prior cannot be lowered
        once and cached.
        """
        lowered = self._priors.get(name)
        if lowered is not None:
            return lowered
        return lower_hierarchical(
            self._hierarchical[name],
            dict(resolved),
            parameter=name,
            dtype=self._dtype,
            device=self._device,
        )

    def log_prior_tensor(self, values: Any) -> torch.Tensor:
        """:meth:`lnprior` as a tensor, differentiable in *values*.

        Fixed parameters contribute nothing; an array-valued parameter
        contributes the sum over its elements (``lowering.md`` §5's row: the
        prior is i.i.d. across elements, which is what ``Independent(base,
        len(shape))`` would express structurally and what summing computes);
        hierarchical priors are bound to the current values of what they
        reference, in topological order. ``-inf`` if any contribution is
        non-finite, exactly as ``ParameterSet.lnprior`` reports.

        **Where the short circuit went** (W2.4 slice 2). This used to return
        early on the first non-finite contribution, mirroring
        ``ParameterSet.lnprior``'s own ``math.isfinite`` short circuit. The
        value is unchanged — a prior contribution is never ``+inf``, so the sum
        of a set containing a ``-inf`` is ``-inf`` or ``nan``, and both become
        ``-inf`` at the end — but the *mechanism* is now :func:`torch.where`
        rather than a Python ``if``, for the same reason
        :mod:`ampere.backends.torch.problem` gives for the density: a Python
        branch on a tensor's value is data-dependent control flow, and
        ``torch.func.vmap`` refuses it. That refusal was the only thing
        standing between this backend and an honest ``BATCHABLE = True``.
        What the short circuit bought was skipping the *remaining priors* for a
        rejected point, which costs a handful of ``log_prob`` calls on scalars
        — nothing beside the model evaluation that follows, and the density
        already declines to skip that for the same reason.
        """
        resolved = self._resolved(values)
        total = as_tensor(0.0, dtype=self._dtype, device=self._device)
        for name in self._order:
            if self._declaration[name].is_fixed:
                continue
            lowered = self.distribution_of(name, resolved)
            total = total + lowered.log_prob(resolved[name]).sum()
        return _finite_or_minus_infinity(total)

    def constrain_tensor(self, unconstrained: Any) -> torch.Tensor:
        """:meth:`constrain` as a tensor, differentiable in *unconstrained*."""
        vector = self._tensor(unconstrained)
        pieces: list[torch.Tensor] = []
        for name in self._free:
            transform = self._transforms[name]
            pieces.append(transform(vector[self._slices[name]]))
        return torch.cat(pieces) if pieces else vector

    def unconstrain_tensor(self, values: Any) -> torch.Tensor:
        """:meth:`unconstrain` as a tensor."""
        flat = self.pack_tensor(self._as_mapping(values))
        pieces: list[torch.Tensor] = []
        for name in self._free:
            transform = self._transforms[name]
            pieces.append(transform.inv(flat[self._slices[name]]))
        return torch.cat(pieces) if pieces else flat

    def _as_mapping(self, values: Any) -> dict[str, torch.Tensor]:
        return self._resolved(values)

    def log_prior_unconstrained_tensor(self, unconstrained: Any) -> torch.Tensor:
        """:meth:`lnprior_unconstrained` as a tensor, differentiable in *y*.

        ``lowering.md`` §2(c)'s definition, verbatim::

            lnprior_unconstrained(y) = lnprior(constrain(y))
                                       + Σ_i log|d constrain_i / dy_i|

        with the Jacobian evaluated **at the unconstrained point**. The
        two-argument call is ``transform.log_abs_det_jacobian(y, x)`` —
        ``y`` unconstrained (the transform's input), ``x`` constrained (its
        output) — which is torch's ``(x_input, y_output)`` order under ampere's
        opposite naming. This is the one place in the package that makes it,
        and a test pins it at a point where the two differ.
        """
        vector = self._tensor(unconstrained)
        constrained = self.constrain_tensor(vector)
        total = self.log_prior_tensor(constrained)
        for name in self._free:
            where = self._slices[name]
            transform = self._transforms[name]
            total = total + transform.log_abs_det_jacobian(vector[where], constrained[where]).sum()
        total = _finite_or_minus_infinity(total)
        return total

    # -- drawing --------------------------------------------------------------

    def sample(self, seed: int | None, label: str = "prior") -> dict[str, np.ndarray]:
        """One draw from the joint prior, on the named sub-stream (``lowering.md`` §9).

        Route (2) of §9.1 where the family has an ``icdf`` — draw uniforms with
        an explicit ``generator=`` and push them through the inverse CDF — and
        route (1), a forked global seed, where it does not. Route (2) is
        preferred because it is the only one that never touches process-global
        state and the only one that is safe under threading; its limitation is
        exactly the ``icdf`` gap, so a backend needs both.
        """
        from .rng import generator  # local: keeps the import graph acyclic

        stream = generator(seed, label)
        resolved: dict[str, np.ndarray] = {
            name: to_numpy(tensor) for name, tensor in self._fixed.items()
        }
        tensors: dict[str, torch.Tensor] = dict(self._fixed)
        for name in self._order:
            parameter = self._declaration[name]
            if parameter.is_fixed:
                continue
            shape = parameter.shape if parameter.shape else ()
            lowered = self.distribution_of(name, tensors)
            if lowered.has_icdf:
                uniforms = torch.rand(
                    shape, generator=stream, dtype=self._dtype, device=self._device
                )
                drawn = lowered.icdf(uniforms)
            else:
                with torch.random.fork_rng(devices=[]):
                    torch.manual_seed(int(torch.randint(0, 2**31 - 1, (), generator=stream)))
                    drawn = lowered.distribution.sample(torch.Size(shape))
            tensors[name] = drawn
            resolved[name] = to_numpy(drawn)
        return resolved

    def __repr__(self) -> str:
        return (
            f"<TorchParameterSpace {self.free_size} free dim(s), "
            f"dtype={self._dtype}, device={self._device}>"
        )


class _PlaceholderPrior:
    """Stands in for a lowered distribution where §4's ``biject_to`` route cannot run.

    :func:`~.lowering.lower_bijection` takes the lowered prior so that an
    *inferred* bijection can be read off its support. A hierarchical prior has
    no lowered object until its hyperparameters have values, but its bijection
    is nonetheless fixed at lowering time — ``default_bijection_for`` refuses a
    hierarchical family whose support would move from draw to draw
    (``ampere.core.parameter._default_bijection_for_hierarchical``), so what is
    passed is always a concrete ``Identity``/``Log``/``Logit`` and the support
    is never consulted. This object exists to make that explicit rather than to
    be used: touching its support is a bug, and says so.
    """

    @property
    def support(self) -> Any:  # pragma: no cover - deliberately unreachable
        raise LoweringError(
            "hierarchical",
            backend=BACKEND,
            detail=(
                "a hierarchical prior has no lowered distribution at lowering time, so its "
                "bijection must be the concrete one ampere inferred, not one read off a support"
            ),
        )


_PLACEHOLDER: Any = _PlaceholderPrior()
