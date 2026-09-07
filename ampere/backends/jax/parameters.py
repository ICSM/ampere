"""Lowering a :class:`~ampere.core.ParameterSet` onto numpyro (``lowering.md`` §5).

This module is the declaration-form table made executable: every row of §5,
with the module conventions of §6.2, the buffer rule of §7, the plate shape of
§8, the RNG policy of §9 and the ``icdf`` contract of §3.6.

What a lowering may assume (§1)
-------------------------------
:class:`LoweredParameterSet` asserts §1.6's three preconditions **once**, at
entry, rather than re-checking them per parameter:

* the set is resolved — no deferred parameter reaches lowering, because
  ``merge`` resolves one and a tie group with no prior raises ``TyingError``
  before this point (§1.3);
* every ``shared_as`` is ``None`` — a tie is a *compile-time collapse*, not a
  runtime object, so no backend should ever see, emit or name a tie label
  (§1.2). Two sites that were tied are one parameter here, with one name, one
  slice and one numpyro sample site;
* ``to_spec()`` succeeds — an opaque, duck-typed prior is refused at
  serialisation *and* at lowering, so this module works from
  :class:`~ampere.core.parameter.PriorSpec` throughout and never has to ask
  "can I describe this?" (§1.5).

Failing loudly once, here, is worth more than three subtly different errors
raised halfway through construction.

The two consumers, and why they share this object
-------------------------------------------------
``lowering.md`` §0 splits the jax column in two: a fit that runs NUTS goes
through numpyro's sample sites, while a fit that runs a gradient optimiser or
an SBI simulator goes through pytrees and never builds a numpyro model at all.
Both consume the same ``ParameterSet``, and both are here:

:meth:`LoweredParameterSet.numpyro_model`
    the sample-site view — one ``numpyro.sample`` per free parameter, in the
    core's own topological order, with real ``numpyro.plate``\\ s for plated
    parameters and ``to_event`` for bare array-valued ones (§5.1);
:meth:`LoweredParameterSet.lnprior` and friends
    the flat-vector view the engines and the conformance battery use, computed
    in jax and handed back as plain numpy at the boundary.

Order matters, and is not ours to choose
----------------------------------------
Sites are emitted in ``ParameterSet``'s **topological** order (its ``_order``),
not declaration order, so a hyperparameter is always sampled before anything
referencing it. §8 says explicitly that lowering reuses the order the core
already computes rather than recomputing it. One consequence is worth stating
because it is easy to discover late: numpyro's ``seed`` handler splits its key
once per site, sequentially, in trace order — it does not fold the site name
in — so adding a parameter changes the draws of every parameter emitted after
it. That is deterministic and reproducible, and it is why a run's provenance
must record the parameter-set spec and not merely the seed (§9.2). W1.8 already
hashes ``to_spec()`` into the attrs, so the requirement is met.

Hierarchical priors are built per evaluation
--------------------------------------------
§8: a :class:`~ampere.core.parameter.HierarchicalPrior` cannot be lowered once
and cached, because its arguments are other parameters' *values*. The
non-hierarchical rows are built once at construction; the hierarchical ones are
built inside the evaluation, from the resolved values, every time. An
implementation that hoisted them out as an optimisation would break exactly
those rows, and break them into a plausible-looking wrong answer: the prior
frozen at its initial-value hyperparameters.
"""

from __future__ import annotations

import dataclasses
import warnings
from collections.abc import Mapping, Sequence
from typing import Any

import jax
import jax.numpy as jnp
import numpy as np
import numpyro
import numpyro.distributions as npd
from numpyro.distributions.transforms import Transform

from ampere.core.exceptions import LoweringError
from ampere.core.exceptions import LoweringFallbackWarning as _LoweringFallbackWarning
from ampere.core.lowering import (
    LoweringResolution,
    lookup_bijection_lowering,
    lookup_lowering,
    provenance_entries,
)
from ampere.core.parameter import (
    HierarchicalPrior,
    Parameter,
    ParameterSet,
    PriorSpec,
    describe_prior,
)

from ._config import BACKEND, require_x64
from .bijections import log_abs_det_jacobian, lower_bijection
from .distributions import NATIVE_ICDF, has_native_icdf, lower_prior

__all__ = [
    "LoweredParameterSet",
    "LoweringFallbackWarning",
    "filter_spec",
]

_NEGATIVE_INFINITY = -jnp.inf


#: The shared warning, re-exported for the backwards-compatible import path.
#:
#: **W2.13 took the proposal slice 1 recorded here** (fold-in 8, ruled
#: 2026-09-07): the class this module used to define lives in
#: :mod:`ampere.core.exceptions` now, so ``pytest.warns`` on one type catches
#: whichever backend raised it. The name stays importable from here because
#: this is where a jax user looks for it.
LoweringFallbackWarning = _LoweringFallbackWarning


# ---------------------------------------------------------------------------
# One lowered parameter
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class _Site:
    """One free parameter, lowered: its slice, its distribution and its transform."""

    parameter: Parameter
    where: slice
    #: The numpyro distribution, built once — ``None`` for a hierarchical prior,
    #: which §8 forbids caching because its arguments are other sites' values.
    fixed_distribution: npd.Distribution | None
    #: The unconstraining transform, from ampere's own declared or inferred
    #: bijection through ``biject_to`` (see :mod:`.bijections`).
    transform: Transform
    #: Whether numpyro can invert this family's CDF **in this environment**.
    native_icdf: bool
    family: str
    #: The two registry rows this site was lowered through -- its prior family
    #: and its bijection class. Kept so a run can stamp the *user-registered*
    #: ones in provenance (``lowering.md`` §12.8's hardening: "every registered
    #: row is stamped user-registered in provenance").
    resolutions: tuple[LoweringResolution, ...]

    @property
    def name(self) -> str:
        return self.parameter.name

    @property
    def shape(self) -> tuple[int, ...]:
        return self.parameter.shape

    @property
    def size(self) -> int:
        return int(self.where.stop - self.where.start)

    def distribution(self, resolved: Mapping[str, Any]) -> npd.Distribution:
        """This site's distribution, given the values it may reference.

        Built fresh for a hierarchical prior and returned from the cache
        otherwise — the distinction §8 makes, and the reason a hoisting
        "optimisation" would be a correctness bug rather than a speed-up.
        """
        if self.fixed_distribution is not None:
            return self.fixed_distribution
        prior = self.parameter.prior
        return lower_prior(_hierarchical_spec(prior, resolved), parameter=self.name)

    def batched(self, resolved: Mapping[str, Any]) -> npd.Distribution:
        """The distribution shaped for this site's declared array shape.

        ``fn.expand(shape).to_event(len(shape))`` — one log-density scalar for
        the whole array, which is what ``parameters.md`` defines an
        array-valued parameter to be (i.i.d. across elements, summed).
        ``to_event`` is always passed its argument explicitly: bare
        ``to_event()`` consumes *every* batch dimension, which for a parameter
        that has acquired one from somewhere else silently reinterprets more
        than intended (§5.1).
        """
        base = self.distribution(resolved)
        if not self.shape:
            return base
        return base.expand(self.shape).to_event(len(self.shape))


def _hierarchical_spec(prior: HierarchicalPrior, resolved: Mapping[str, Any]) -> PriorSpec:
    """A :class:`PriorSpec` for *prior* with its references filled in.

    The referenced values may be jax tracers — inside ``numpyro_model`` they
    always are — and ``PriorSpec`` does not coerce its ``kwds``, so they travel
    through the registry to the family constructor unchanged. That is what
    makes ``dist.Normal(mu, sigma)`` with ``mu``/``sigma`` earlier sample sites
    the natural lowering §5 says it is, rather than a special case.
    """
    kwds: dict[str, Any] = dict(prior.kwds)
    for argument, reference in prior.hyperparameters.items():
        kwds[argument] = resolved[reference]
    return PriorSpec(family=prior.family, args=tuple(prior.args), kwds=kwds)


def _placeholder_spec(prior: HierarchicalPrior) -> PriorSpec:
    """*prior* with plausible constants where its references will be.

    Used only to probe ``icdf`` availability at lowering time: whether a family
    can invert its CDF is a property of the family, not of the hyperparameter
    values, but numpyro will not build a distribution without *some* values.
    ``loc``-shaped arguments get 0 and everything else 1, which is inside every
    family in the table's admissible range.
    """
    filled = {
        reference: (0.0 if argument == "loc" else 1.0)
        for argument, reference in prior.hyperparameters.items()
    }
    return _hierarchical_spec(prior, filled)


# ---------------------------------------------------------------------------
# The lowered set
# ---------------------------------------------------------------------------


class LoweredParameterSet:
    """A :class:`~ampere.core.ParameterSet` realised on numpyro and jax.

    Satisfies ``tests/conformance``'s ``ParameterSpace`` protocol — which is
    the point: ``parameters.md`` §13 hands those rows over precisely so that a
    lowered path has to reproduce the reference semantics rather than merely
    produce "a valid unconstrained density".

    Parameters
    ----------
    declaration
        The **merged** set (``mapping.merged``), never a component set: §1.1
        makes ``merge`` the resolution step, and by the time a backend sees a
        parameter, tying has already happened.
    strict
        The fitting problem's own ``strict`` flag, threaded here (§3.6, pinned
        at the freeze: one flag, one meaning). ``False`` takes the reference
        ``ppf`` fallback for families with no native ``icdf`` and warns once;
        ``True`` raises :class:`~ampere.core.exceptions.LoweringError` naming
        them instead.

    Notes
    -----
    Every method on the engine-facing surface takes and returns **plain
    numpy**, computes in jax in between, and never stores a PRNG key. The
    conversion at the boundary is deliberate and is what the protocol asks for:
    the battery asserts agreement of *values*, not of container types.
    """

    def __init__(self, declaration: ParameterSet, *, strict: bool = False) -> None:
        require_x64("a jax-lowered parameter set")
        _entry_assertion(declaration)
        self._declaration = declaration
        self._strict = bool(strict)
        self._sites = _sites_of(declaration)
        self._by_name = {site.name: site for site in self._sites}
        #: Fixed parameters, as jax arrays. §5's "fixed" row: not a sample site,
        #: no prior object constructed, and — on the pytree side — an ordinary
        #: array leaf **excluded** from the trainable partition (§6.2, §7).
        self.constants: dict[str, jax.Array] = {
            parameter.name: jnp.asarray(parameter.value, dtype=jnp.float64)
            for parameter in declaration
            if parameter.is_fixed
        }
        #: Emission order: the core's topological order, reused rather than
        #: recomputed (§8). Fixed parameters are in it and are skipped.
        # §8 says explicitly to reuse the order the core already computes
        # rather than recompute it. W2.13 made it public
        # (``ParameterSet.evaluation_order()``, fold-in 9) precisely because
        # this line used to reach past the underscore for it.
        self._order: tuple[str, ...] = tuple(declaration.evaluation_order())
        self._fallback_families = tuple(
            sorted({site.family for site in self._sites if not site.native_icdf})
        )
        self._announce_fallback()

    # -- §3.6: the icdf contract ------------------------------------------

    def _announce_fallback(self) -> None:
        """Warn once, or raise, if the prior transform will run on numpy.

        Once per *run*, not per call, because the decision is made here — at
        lowering time, before any sampling — and repeating it per proposal
        would drown the run in a message that says nothing new.
        """
        if not self._fallback_families:
            return
        families = ", ".join(repr(family) for family in self._fallback_families)
        if self._strict:
            raise LoweringError(
                self._fallback_families[0],
                backend=BACKEND,
                detail=(
                    f"prior family/families {families} have no usable native icdf on the "
                    f"{BACKEND!r} backend, so prior_transform would be computed on the reference "
                    f"(scipy) path. strict=True refuses to mix paths. Either change these priors "
                    f"to families numpyro can invert ({', '.join(sorted(NATIVE_ICDF))}), install "
                    f"what numpyro needs to invert them, or run with strict=False and accept the "
                    f"documented fallback."
                ),
            )
        warnings.warn(
            f"ampere.backends.{BACKEND}: prior family/families {families} have no usable native "
            f"icdf on the {BACKEND!r} backend, so this run's prior_transform is evaluated on the "
            f"reference (scipy) path. This is exact — prior_transform has one mathematical "
            f"definition — and cheap, since it is called once per proposal on a vector of size "
            f"{self.free_size} rather than inside the likelihood's hot loop. Pass "
            f"FittingProblem(strict=True) to refuse it instead.",
            LoweringFallbackWarning,
            stacklevel=3,
        )

    @property
    def uses_reference_prior_transform(self) -> bool:
        """Whether :meth:`prior_transform` runs on numpy (``lowering.md`` §3.6)."""
        return bool(self._fallback_families)

    @property
    def fallback_families(self) -> tuple[str, ...]:
        """The families that forced :attr:`uses_reference_prior_transform`."""
        return self._fallback_families

    # -- layout (backend-neutral facts, delegated) -------------------------

    @property
    def declaration(self) -> ParameterSet:
        """The backend-neutral declaration this space realises."""
        return self._declaration

    @property
    def free_size(self) -> int:
        return self._declaration.free_size

    def free_labels(self) -> tuple[str, ...]:
        """One label per flat free dimension.

        Delegated, and deliberately: the flat layout is a property of the
        *declaration*, and a backend that invented its own would put a
        translation layer between the sampler's output and the results schema —
        which is the bug class ``DEVELOPMENT_PLAN.md`` §4.1 exists to end. It is
        also why numpyro site names below are these names verbatim.
        """
        return self._declaration.free_labels()

    @property
    def sites(self) -> tuple[_Site, ...]:
        """The lowered free parameters, in declaration order."""
        return self._sites

    # -- §12.8's provenance hardening --------------------------------------

    @property
    def lowering_resolutions(self) -> tuple[LoweringResolution, ...]:
        """Every registry row this lowering consulted, de-duplicated.

        Both slots: the prior family of each site and the class of each site's
        bijection. Built-in rows included — filtering them is
        :func:`~ampere.core.lowering.provenance_entries`'s job, and a caller
        introspecting what was used should see everything.
        """
        seen: dict[tuple[str, str, str], LoweringResolution] = {}
        for site in self._sites:
            for resolution in site.resolutions:
                seen[(resolution.kind, resolution.name, resolution.backend)] = resolution
        return tuple(seen.values())

    def lowering_provenance(self) -> list[dict[str, Any]]:
        """The **user-registered** rows this lowering consulted, netCDF-safe.

        ``lowering.md`` §12.8's hardening, ready for
        ``ampere.results.provenance.provenance_attrs``'s ``extra=``::

            provenance_attrs(problem, extra={"registered_lowerings": lowered.lowering_provenance()})

        Empty for a run that used only ampere's own table, which is the point:
        the signal a reviewer wants is "did this run depend on something
        outside the conformance suite's guarantees?", and stamping ten built-in
        rows would bury it. Nothing in ``ampere.inference`` reads this yet —
        W2.6 deferred the *first-class* provenance key "until a real backend
        drives lowering end-to-end", and this backend does, but every row it
        uses is built-in so the key would always be empty today. See this
        branch's report.
        """
        return provenance_entries(self.lowering_resolutions)

    # -- pack / unpack ------------------------------------------------------

    def pack(self, values: Mapping[str, Any]) -> np.ndarray:
        """Flatten a name → value mapping into the free vector."""
        if not self._sites:
            return np.zeros(0, dtype=float)
        chunks = [
            jnp.reshape(jnp.asarray(values[site.name], dtype=jnp.float64), (-1,))
            for site in self._sites
        ]
        return np.asarray(jnp.concatenate(chunks))

    def unpack(self, theta: Any) -> dict[str, Any]:
        """Expand a free vector into every parameter, fixed ones included."""
        vector = jnp.asarray(theta, dtype=jnp.float64).reshape(-1)
        out: dict[str, Any] = {}
        for parameter in self._declaration:
            if parameter.is_fixed:
                out[parameter.name] = np.asarray(self.constants[parameter.name])
                if not parameter.shape:
                    out[parameter.name] = float(out[parameter.name])
                continue
            site = self._by_name[parameter.name]
            chunk = vector[site.where]
            out[parameter.name] = (
                np.asarray(chunk.reshape(site.shape)) if site.shape else float(chunk[0])
            )
        return out

    def _resolved(self, vector: jax.Array) -> dict[str, Any]:
        """Every parameter's value, keyed by name, from a free vector.

        Kept in jax throughout — this is the mapping a hierarchical prior's
        references are looked up in, so converting here would break the
        gradient path the whole backend exists for.
        """
        resolved: dict[str, Any] = dict(self.constants)
        for site in self._sites:
            chunk = vector[site.where]
            resolved[site.name] = chunk.reshape(site.shape) if site.shape else chunk[0]
        return resolved

    # -- priors -------------------------------------------------------------

    def log_prior(self, theta: Any) -> jax.Array:
        """The joint log prior at *theta*, in jax. Traceable; no Python branching.

        The reference implementation short-circuits to ``-inf`` "as soon as any
        contribution is non-finite, without evaluating the rest". That
        short-circuit is a Python ``if`` on a value, so the traced version
        cannot have it: every term is computed and the total is masked with
        :func:`jnp.where`. The *answer* is identical — a sum containing ``-inf``
        is ``-inf`` — and the cost is that a term is evaluated which the
        reference path would have skipped, which is the right trade for a
        function that must survive ``jit`` and ``grad``.
        """
        vector = jnp.asarray(theta, dtype=jnp.float64).reshape(-1)
        resolved = self._resolved(vector)
        total = jnp.asarray(0.0, dtype=jnp.float64)
        for name in self._order:
            site = self._by_name.get(name)
            if site is None:
                continue
            value = resolved[site.name]
            distribution = site.batched(resolved)
            inside = jnp.all(distribution.support(value))
            term = jnp.sum(distribution.log_prob(value))
            total = total + jnp.where(inside, term, _NEGATIVE_INFINITY)
        return jnp.where(jnp.isfinite(total), total, _NEGATIVE_INFINITY)

    def lnprior(self, values: Any) -> float:
        """``log p(θ)`` in the constrained parameterisation. ``-inf`` off support."""
        return float(self.log_prior(values))

    def prior_transform(self, unit_cube: Any) -> np.ndarray:
        """Map ``[0, 1]^n`` to the constrained free vector — nested sampling.

        An inverse-CDF composition, evaluated in the core's dependency order so
        a hyperparameter is always transformed before anything referencing it.

        When any family has no usable native ``icdf`` this delegates to the
        reference implementation wholesale, which is ``lowering.md`` §3.6's
        sanctioned fallback: the quantity has one mathematical definition, so
        computing it in numpy changes nothing about the posterior, and it is
        called once per proposal on a vector of length ``free_size`` rather than
        inside the likelihood. The decision was announced at construction.
        """
        if self.uses_reference_prior_transform:
            return self._declaration.prior_transform(unit_cube)
        cube = jnp.asarray(unit_cube, dtype=jnp.float64).reshape(-1)
        if cube.size != self.free_size:
            raise LoweringError(
                "prior_transform",
                backend=BACKEND,
                detail=f"expected a unit-cube vector of length {self.free_size}, got {cube.size}",
            )
        flat = jnp.zeros(self.free_size, dtype=jnp.float64)
        resolved: dict[str, Any] = dict(self.constants)
        for name in self._order:
            site = self._by_name.get(name)
            if site is None:
                continue
            drawn = jnp.reshape(site.distribution(resolved).icdf(cube[site.where]), (-1,))
            flat = flat.at[site.where].set(drawn)
            resolved[site.name] = drawn.reshape(site.shape) if site.shape else drawn[0]
        return np.asarray(flat)

    # -- unconstrained space -----------------------------------------------

    def constrain_jax(self, unconstrained: Any) -> jax.Array:
        """``constrain`` in jax: the traceable half of the unconstrained map."""
        y = jnp.asarray(unconstrained, dtype=jnp.float64).reshape(-1)
        out = jnp.zeros_like(y)
        for site in self._sites:
            out = out.at[site.where].set(jnp.reshape(site.transform(y[site.where]), (-1,)))
        return out

    def constrain(self, unconstrained: Any) -> np.ndarray:
        """Map the unconstrained vector into the priors' support."""
        return np.asarray(self.constrain_jax(unconstrained))

    def unconstrain(self, values: Any) -> np.ndarray:
        """Map values in the priors' support to the unconstrained vector."""
        x = jnp.asarray(self.pack(self._declaration.complete(values)), dtype=jnp.float64)
        out = jnp.zeros_like(x)
        for site in self._sites:
            out = out.at[site.where].set(jnp.reshape(site.transform.inv(x[site.where]), (-1,)))
        return np.asarray(out)

    def log_prior_unconstrained(self, unconstrained: Any) -> jax.Array:
        """``lnprior(constrain(y)) + Σ log|d constrain / dy|``, in jax.

        The density a NUTS kernel works with, and the identity ``lowering.md``
        §2(c) makes every backend reproduce — with the Jacobian evaluated **at
        the unconstrained point** and summed over every element of every free
        parameter. Not "a valid unconstrained density": a density differing by a
        non-constant function of ``y`` is a different posterior, and one
        differing by a constant is fine for MCMC and wrong for evidence.
        """
        y = jnp.asarray(unconstrained, dtype=jnp.float64).reshape(-1)
        constrained = self.constrain_jax(y)
        total = self.log_prior(constrained)
        for site in self._sites:
            chunk = y[site.where]
            total = total + jnp.sum(
                log_abs_det_jacobian(site.transform, chunk, constrained[site.where])
            )
        return jnp.where(jnp.isfinite(total), total, _NEGATIVE_INFINITY)

    def lnprior_unconstrained(self, unconstrained: Any) -> float:
        """Log prior density in unconstrained space, Jacobian term included."""
        return float(self.log_prior_unconstrained(unconstrained))

    # -- the numpyro view (§5, §8) -----------------------------------------

    def numpyro_model(self) -> Any:
        """A numpyro model function emitting one site per free parameter.

        Site names are the merged parameter names **verbatim** — already
        unique, already qualified, and already what ``free_labels`` uses as
        ArviZ coordinates. Inventing a second naming scheme at the numpyro
        boundary would put a translation layer between the sampler's output and
        the results schema (§8).

        A plated parameter lowers to a real ``numpyro.plate``, not a bare
        batched site: the two give the same joint density here, so the
        distinction is **structural** rather than numerical, and downstream
        tooling (subsampling, hierarchical diagnostics, ArviZ dimension naming)
        depends on the declaration (§5.1). ``Parameter.plate`` is the flag that
        tells them apart, and it is set precisely by ``Plate.expand()``.

        v1.3 supports one plate dimension per parameter, so exactly one plate is
        open at a time and numpyro's ``dim`` allocation — outermost ``-1``,
        inner plates leftwards, and two plates *constructed* before either is
        entered silently colliding on ``-1`` — cannot bite. It becomes live the
        moment nested plates are added, and the lowering must then pass ``dim``
        explicitly rather than rely on allocation order.
        """

        def model() -> dict[str, Any]:
            resolved: dict[str, Any] = dict(self.constants)
            for name in self._order:
                site = self._by_name.get(name)
                if site is None:
                    continue
                distribution = site.distribution(resolved)
                plate = site.parameter.plate
                if plate is not None and site.shape:
                    with numpyro.plate(plate, site.shape[0]):
                        drawn = numpyro.sample(site.name, distribution)
                elif site.shape:
                    drawn = numpyro.sample(
                        site.name, distribution.expand(site.shape).to_event(len(site.shape))
                    )
                else:
                    drawn = numpyro.sample(site.name, distribution)
                resolved[site.name] = drawn
            return resolved

        return model

    def __repr__(self) -> str:
        return (
            f"<LoweredParameterSet {self.free_size} free dimension(s), "
            f"{len(self.constants)} fixed, backend={BACKEND!r}>"
        )


# ---------------------------------------------------------------------------
# Construction helpers
# ---------------------------------------------------------------------------


def _entry_assertion(declaration: ParameterSet) -> None:
    """``lowering.md`` §1.6, run once per lowering."""
    if not declaration.is_resolved:
        raise LoweringError(
            "<deferred>",
            backend=BACKEND,
            detail=(
                f"parameter(s) {list(declaration.deferred_names)} are declared shared without a "
                f"prior of their own, so this set has not been merged. Lowering consumes "
                f"ParameterSet.merge(...).merged and nothing else (lowering.md §1.1/§1.3); a "
                f"deferred parameter is by construction not evaluable, and inventing an improper "
                f"prior for it would be a silent change of model."
            ),
        )
    unmerged = [p.name for p in declaration if p.shared_as is not None]
    if unmerged:
        raise LoweringError(
            "<tie label>",
            backend=BACKEND,
            detail=(
                f"parameter(s) {unmerged} still carry a tie label, so this set has not been "
                f"merged. A tie is a compile-time collapse, not a runtime object (lowering.md "
                f"§1.2): no backend should ever see, emit or name a tie label."
            ),
        )
    declaration.to_spec()


def _lower_site(parameter: Parameter, declaration: ParameterSet) -> _Site | None:
    """Lower one free parameter, or ``None`` for a fixed one (§5's "fixed" row)."""
    if parameter.is_fixed:
        return None
    prior = parameter.prior
    hierarchical = isinstance(prior, HierarchicalPrior)
    if hierarchical:
        family = prior.family
        probe = lower_prior(_placeholder_spec(prior), parameter=parameter.name)
        distribution = None
    else:
        spec = describe_prior(prior)
        family = spec.family
        distribution = lower_prior(spec, parameter=parameter.name)
        probe = distribution
    bijection = parameter.unconstraining_bijection()
    return _Site(
        parameter=parameter,
        where=declaration.free_slice(parameter.name),
        fixed_distribution=distribution,
        transform=lower_bijection(bijection),
        native_icdf=has_native_icdf(probe),
        family=family,
        resolutions=(
            lookup_lowering(family, BACKEND),
            lookup_bijection_lowering(type(bijection), BACKEND),
        ),
    )


def _sites_of(declaration: ParameterSet) -> tuple[_Site, ...]:
    return tuple(
        site
        for site in (_lower_site(parameter, declaration) for parameter in declaration)
        if site is not None
    )


# ---------------------------------------------------------------------------
# §6.2 and §7: the trainable/non-trainable split
# ---------------------------------------------------------------------------


def filter_spec(module: Any, trainable: Sequence[str] | None = None) -> Any:
    """A boolean pytree marking *module*'s trainable leaves, for ``eqx.partition``.

    ``lowering.md`` §6.2 and ``DEVELOPMENT_PLAN.md`` §2's "jax non-trainable
    mechanism" row (ruled 2026-09-01) fix the mechanism: **an explicit
    ``eqx.partition`` filter spec**, built as a boolean pytree mirroring the
    module, with the trainable leaves flipped to ``True``. Then
    ``trainable, frozen = eqx.partition(module, spec)`` and ``eqx.combine`` inside
    the differentiated function.

    Two things this is *not*, and both are the trap:

    * it is **never** ``eqx.field(static=True)``. A static field's value goes
      into the flattening's auxiliary data — part of the pytree *structure*,
      which jax compares and hashes to build JIT cache keys — so an array there
      is hashed on every traced call, any change to its contents triggers a full
      recompilation, and large constants are held alive in the compilation
      cache. Static is for genuinely static, small, hashable metadata: an
      integer size, a string name, a boolean flag. ``Plate.size`` is static; a
      plate member's values are not (§7);
    * it is not ``paramax.non_trainable`` or a bare ``lax.stop_gradient``.
      Both leave the leaf *in* the differentiated pytree and rest on a
      discipline someone has to remember at every call site — forgetting to
      call ``unwrap`` inside the loss silently trains the buffers, with nothing
      raising and nothing warning. A partition filter *structurally* cannot
      leak, because the leaf is not in the differentiated argument at all.

    Parameters
    ----------
    module
        Any pytree — an :class:`equinox.Module`, a dict of arrays, a
        :class:`~ampere.core.Parameterised` object's lowered leaves.
    trainable
        Names of the attributes that carry free parameters. ``None`` — the
        default — means "every inexact array leaf", which is equinox's own
        ``is_inexact_array`` filter and the right answer for a module whose
        buffers have already been separated out.

    Returns
    -------
    pytree
        Same structure as *module*, with ``True`` at each trainable leaf and
        ``False`` everywhere else.
    """
    import equinox as eqx

    if trainable is None:
        return jax.tree.map(eqx.is_inexact_array, module)
    spec = jax.tree.map(lambda _: False, module)
    wanted = tuple(trainable)
    return eqx.tree_at(
        lambda tree: [getattr(tree, name) for name in wanted],
        spec,
        replace=[True] * len(wanted),
    )
