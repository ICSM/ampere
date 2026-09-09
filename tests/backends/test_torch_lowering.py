"""``lowering.md`` §11's rows for the torch backend.

§11 lists what the conformance suite owes for lowering, and each row is "a
place a backend can be wrong without crashing, which is the criterion for
needing a mechanical check". The backend-neutral battery in
``tests/conformance`` covers the rows it can express through the seven-member
fixture protocol — the ``lnprior_unconstrained`` agreement, the flat layout,
the prior transform against scipy's quantiles. This file covers the rest,
which are all statements about *lowering* rather than about a parameter space:

* row 1, the **Jacobian argument order**, at a point where the two arguments
  differ, so that a swap fails rather than producing a plausible posterior;
* row 3, the **distribution parametrisation**, family by family, against the
  scipy original — including the ``loc``/``scale`` conversions §3.3 spells out,
  each with a non-default value, "since a case with defaults passes under the
  wrong mapping";
* row 4, **support preservation**, which is the check that catches a dropped
  ``loc`` — and here also the check that catches trusting torch's
  ``TransformedDistribution.support``, which is a codomain;
* row 5, **unsupported families raise**, naming the family and the backend,
  with no silent fallback;
* row 7, **hierarchical priors track their hyperparameters** — the check that
  catches a distribution object hoisted out of the evaluation loop;
* row 8, **buffers are not trainable**;
* row 9, **fixed parameters take no sampler dimension**;
* row 10, **seed reproducibility** within one backend.

Plus the §3.6 ``icdf``-fallback contract, which is W2.4's own acceptance
criterion, and the §1.6 entry assertion.

Oracles are ``scipy`` or closed forms. Nothing here compares ampere against
ampere except where the *claim* is an agreement between two ampere paths (§4's
``biject_to`` cross-check, and the reference oracle §2(c) makes the native path
agree with).
"""

from __future__ import annotations

import math
import warnings

import numpy as np
import pytest
import scipy.stats as st

pytest.importorskip("torch", reason="the torch backend needs ampere[torch]")

import torch
from torch.distributions import biject_to

from typing import Any

import astropy.units as u

from ampere.backends.torch import (
    BACKEND,
    DenseGP,
    FractionalModelGPNoise,
    FractionalModelNoise,
    GaussianProcessNoise,
    IndependentNoise,
    Matern32,
    PowerLaw,
    QuasisepGP,
    lower_problem,
)
from ampere.backends.torch.lowering import (
    LoweringFallbackWarning,
    lower_bijection,
    lower_hierarchical,
    lower_prior,
)
from ampere.backends.torch.parameters import TorchParameterSpace
from ampere.backends.torch.rng import generator, seed_for
from ampere.core import (
    CauchyFamily,
    Censoring,
    Dataset,
    FittingProblem,
    GaussianFamily,
    HierarchicalPrior,
    Identity,
    Likelihood,
    LimitKind,
    Log,
    Logit,
    Parameter,
    ParameterError,
    ParameterSet,
    Plate,
    PoissonFamily,
    Spectrum,
    StudentTFamily,
    describe_prior,
)
from ampere.backends.torch._families import FamilyInputs, native_log_prob
from ampere.core import Matern32 as CoreMatern32
from ampere.core import QuasisepGP as CoreQuasisepGP
from ampere.core.exceptions import LoweringError
from ampere.core.lowering import (
    lookup_bijection_lowering,
    lookup_lowering,
    register_bijection_lowering,
    run_registrant_battery,
)
from ampere.core.rng import substream


def tensor(values: object) -> torch.Tensor:
    """A float64 tensor, the only precision this backend works in."""
    return torch.as_tensor(np.asarray(values, dtype=float), dtype=torch.float64)


# ---------------------------------------------------------------------------
# §3.2's table, both tiers. Every row carries a non-default loc/scale where the
# family has one, because "a case with defaults passes under the wrong mapping"
# (§11 row 3).
# ---------------------------------------------------------------------------

FAMILIES = {
    "norm": st.norm(1.5, 2.5),
    "uniform": st.uniform(2.0, 3.0),
    "halfnorm": st.halfnorm(0.0, 2.0),
    "halfnorm-shifted": st.halfnorm(3.0, 2.0),
    "loguniform": st.loguniform(0.1, 10.0),
    "lognorm": st.lognorm(0.5, 0.0, 2.0),
    "lognorm-shifted": st.lognorm(0.5, 1.0, 2.0),
    "expon": st.expon(0.0, 2.0),
    "expon-shifted": st.expon(1.5, 2.0),
    "gamma": st.gamma(2.0, 0.0, 3.0),
    "gamma-shifted": st.gamma(2.0, 1.0, 3.0),
    "beta": st.beta(2.0, 3.0),
    "beta-stretched": st.beta(2.0, 3.0, 1.0, 4.0),
}

#: The discrete row: it lowers as a *distribution* (§3.5) and never onto a
#: gradient path, so it is checked separately from the continuous families.
POISSON = st.poisson(3.0)


@pytest.mark.parametrize("name", sorted(FAMILIES))
class TestDistributionParametrisation:
    """§11 rows 3 and 4: the density and the support, against the scipy original."""

    @staticmethod
    def points(prior: object) -> np.ndarray:
        """Nine quantiles spanning the support — inside it by construction."""
        return np.asarray(prior.ppf(np.linspace(0.05, 0.95, 9)), dtype=float)

    def test_the_log_density_matches_scipy(self, name: str) -> None:
        prior = FAMILIES[name]
        lowered = lower_prior(describe_prior(prior))
        xs = self.points(prior)
        got = lowered.log_prob(tensor(xs)).numpy()
        assert got == pytest.approx(prior.logpdf(xs), abs=1e-10)

    def test_the_support_matches_scipy(self, name: str) -> None:
        """The check that catches a dropped ``loc`` — and a trusted codomain.

        ``lowering.md`` §3.3: lowering ``halfnorm(3, 2)`` to ``HalfNormal(2)``
        would silently move the prior's support by three units. And torch's
        ``TransformedDistribution.support`` is the last transform's *codomain*,
        so the composed rows would report a support far larger than they have
        unless the lowering declares it.
        """
        prior = FAMILIES[name]
        low, high = (float(bound) for bound in prior.support())
        support = lower_prior(describe_prior(prior)).support
        lower_bound = getattr(support, "lower_bound", -math.inf)
        upper_bound = getattr(support, "upper_bound", math.inf)
        assert float(lower_bound) == pytest.approx(low, abs=1e-12)
        assert float(upper_bound) == pytest.approx(high, abs=1e-12)

    def test_no_prior_mass_outside_the_support(self, name: str) -> None:
        """``-inf`` rather than a raise or a ``nan``, as ``log_density`` gives."""
        prior = FAMILIES[name]
        lowered = lower_prior(describe_prior(prior))
        low, high = (float(bound) for bound in prior.support())
        outside = [value for value in (low - 1.0, high + 1.0) if math.isfinite(value)]
        if not outside:
            pytest.skip(f"{name} is supported on the whole real line: there is no outside")
        got = lowered.log_prob(tensor(outside)).numpy()
        assert np.all(np.isneginf(got))

    def test_the_registered_row_self_certifies(self, name: str) -> None:
        """``run_registrant_battery`` — W2.6's opt-in check — passes on every row.

        Not a duplicate of the density row above: this one goes through
        ``ampere.core.lowering``'s own comparison, so it also asserts that the
        registry row for this family really is the one the backend uses.
        """
        report = run_registrant_battery(
            FAMILIES[name],
            BACKEND,
            log_prob=lambda native, xs: native.log_prob(tensor(xs)),
            to_numpy=lambda values: values.numpy(),
        )
        assert report.passed


class TestDiscreteFamily:
    """§3.5: ``poisson`` lowers as a distribution, never onto a gradient path."""

    def test_the_probability_mass_matches_scipy(self) -> None:
        lowered = lower_prior(describe_prior(POISSON))
        counts = np.arange(0.0, 9.0)
        assert lowered.log_prob(tensor(counts)).numpy() == pytest.approx(
            POISSON.logpmf(counts), abs=1e-10
        )

    def test_it_has_no_native_icdf(self) -> None:
        """The §3.6 gap, established by asking rather than from a table."""
        assert not lower_prior(describe_prior(POISSON)).has_icdf

    def test_the_bijection_is_refused_upstream(self) -> None:
        """The refusal lives in ``default_bijection_for``, and only there.

        A discrete parameter never reaches :func:`lower_bijection` at all, so
        this backend needs no rule of its own for it — which is exactly what
        the 2026-09-03 ruling arranged.
        """
        from ampere.core.exceptions import CapabilityError

        with pytest.raises(CapabilityError, match="discrete"):
            Parameter("counts", POISSON).unconstraining_bijection()


class TestUnsupportedFamilies:
    """§11 row 5: a family with no exact construction raises, and does not fall back."""

    def test_truncnorm_raises_naming_the_family_the_parameter_and_the_backend(self) -> None:
        spec = describe_prior(st.truncnorm(-1.0, 2.0, 5.0, 3.0))
        with pytest.raises(LoweringError) as raised:
            lower_prior(spec, parameter="temperature")
        error = raised.value
        assert error.family == "truncnorm"
        assert error.parameter == "temperature"
        assert error.backend == BACKEND
        assert "TruncatedNormal" in str(error)

    def test_there_is_no_registry_row_to_fall_back_to(self) -> None:
        """§3.4 rule 3: not a silent approximation, and not a quiet reference path."""
        with pytest.raises(LoweringError, match="truncnorm"):
            lookup_lowering("truncnorm", BACKEND)

    def test_a_parameter_set_containing_one_refuses_to_lower(self) -> None:
        declaration = ParameterSet([Parameter("t", st.truncnorm(-1.0, 2.0, 5.0, 3.0))])
        with pytest.raises(LoweringError, match="truncnorm"):
            TorchParameterSpace(declaration)


# ---------------------------------------------------------------------------
# §4 and §2: bijections and the Jacobian's argument order
# ---------------------------------------------------------------------------

BIJECTIONS = {
    "identity": (Identity(), np.array([-1.3, 0.7, 2.1])),
    "log-at-zero": (Log(lower=0.0), np.array([-1.3, 0.7, 2.1])),
    "log-shifted": (Log(lower=3.0), np.array([-1.3, 0.7, 2.1])),
    "logit-unit": (Logit(lower=0.0, upper=1.0), np.array([-1.3, 0.7, 2.1])),
    "logit-interval": (Logit(lower=2.0, upper=5.0), np.array([-1.3, 0.7, 2.1])),
}


@pytest.mark.parametrize("name", sorted(BIJECTIONS))
class TestBijectionTable:
    """§4's table, against ``ampere.core``'s own bijections."""

    @staticmethod
    def lowered(name: str):
        bijection, ys = BIJECTIONS[name]
        # A prior is needed only to satisfy the signature; a *declared*
        # bijection never consults it, which the placeholder in
        # ``parameters.py`` makes explicit.
        return lower_bijection(bijection, lower_prior(describe_prior(st.norm(0.0, 1.0)))), ys

    def test_constrain_is_the_forward_direction(self, name: str) -> None:
        """§2(a): ``biject_to``'s forward maps the real line onto the support."""
        bijection, _ = BIJECTIONS[name]
        transform, ys = self.lowered(name)
        assert transform(tensor(ys)).numpy() == pytest.approx(bijection.constrain(ys), abs=1e-12)

    def test_the_inverse_is_unconstrain(self, name: str) -> None:
        bijection, _ = BIJECTIONS[name]
        transform, ys = self.lowered(name)
        xs = bijection.constrain(ys)
        assert transform.inv(tensor(xs)).numpy() == pytest.approx(ys, abs=1e-9)

    def test_the_jacobian_takes_the_unconstrained_point_first(self, name: str) -> None:
        """§11 row 1, and the reason it is written down at all.

        torch's signature is ``log_abs_det_jacobian(x, y)`` with ``x`` the
        transform's *input*. Ampere's convention names the unconstrained point
        ``y``, so the correct call is ``t.log_abs_det_jacobian(y, x)`` and the
        identifiers look transposed. The point is chosen so that the two
        arguments differ.
        """
        bijection, _ = BIJECTIONS[name]
        transform, ys = self.lowered(name)
        xs = bijection.constrain(ys)
        if name != "identity":
            # §11 insists on a point where the two arguments differ; the
            # identity is the one row where none exists, which is also why the
            # swap row below skips it.
            assert not np.allclose(xs, ys), "the row must be checked where x != y"
        got = transform.log_abs_det_jacobian(tensor(ys), tensor(xs)).numpy()
        assert got == pytest.approx(bijection.log_abs_det_jacobian(ys), abs=1e-12)

    def test_swapping_the_jacobian_arguments_is_detected(self, name: str) -> None:
        """The row above must *fail* if the arguments are swapped (§11 row 1).

        ``Identity`` is exempt: its Jacobian is zero from either direction, so
        there is nothing a swap could change — which is why §11 insists the
        check be made "at a point where the two arguments differ" and why a
        test at a symmetric point proves nothing.
        """
        if name == "identity":
            pytest.skip("the identity's Jacobian is zero whichever way round it is called")
        bijection, _ = BIJECTIONS[name]
        transform, ys = self.lowered(name)
        xs = bijection.constrain(ys)
        swapped = transform.log_abs_det_jacobian(tensor(xs), tensor(ys)).numpy()
        assert not np.allclose(swapped, bijection.log_abs_det_jacobian(ys))


@pytest.mark.parametrize("name", sorted(FAMILIES))
def test_biject_to_agrees_with_the_hand_built_table(name: str) -> None:
    """§4's own cross-check, run rather than assumed.

    §4 says to prefer ``biject_to(lowered.support)`` for an *inferred*
    bijection, with the table as documentation of what one will get. That only
    holds if the two routes agree, and they only agree if the lowered
    distribution reports its true support — which for the composed rows it does
    only because this backend declares it. So this row is simultaneously §4's
    cross-check and the regression test for the codomain trap.
    """
    prior = FAMILIES[name]
    parameter = Parameter("x", prior)
    lowered = lower_prior(describe_prior(prior))
    inferred = biject_to(lowered.support)
    declared = lower_bijection(parameter.unconstraining_bijection(), lowered)
    ys = tensor([-1.3, 0.0, 0.7, 2.1])
    assert inferred(ys).numpy() == pytest.approx(declared(ys).numpy(), abs=1e-9)
    assert inferred.log_abs_det_jacobian(ys, inferred(ys)).numpy() == pytest.approx(
        declared.log_abs_det_jacobian(ys, declared(ys)).numpy(), abs=1e-9
    )


class CustomBijection:
    """A user-supplied :class:`~ampere.core.Bijection`: numpy, and torch-hostile."""

    def constrain(self, y):
        return np.asarray(y) * 3.0

    def unconstrain(self, x):
        return np.asarray(x) / 3.0

    def log_abs_det_jacobian(self, y):
        return np.full(np.shape(y), math.log(3.0))


class TestCustomBijections:
    """§4's last row: unsupported, until someone registers one."""

    def test_a_custom_bijection_raises_naming_the_registration_hook(self) -> None:
        with pytest.raises(LoweringError) as raised:
            lower_bijection(
                CustomBijection(),
                lower_prior(describe_prior(st.norm(0.0, 1.0))),
                parameter="scale",
            )
        assert "register_bijection_lowering" in str(raised.value)
        assert raised.value.backend == BACKEND

    def test_registering_one_makes_it_work(self) -> None:
        """§12.8's hook, consumed by the backend exactly as a third party would."""
        from torch.distributions.transforms import AffineTransform

        def build(_bijection):
            return AffineTransform(loc=tensor(0.0), scale=tensor(3.0))

        try:
            lookup_bijection_lowering(CustomBijection, BACKEND)
        except LoweringError:
            register_bijection_lowering(CustomBijection, BACKEND, build)
        transform = lower_bijection(
            CustomBijection(), lower_prior(describe_prior(st.norm(0.0, 1.0)))
        )
        ys = np.array([-1.0, 0.5, 2.0])
        assert transform(tensor(ys)).numpy() == pytest.approx(
            CustomBijection().constrain(ys), abs=1e-12
        )


# ---------------------------------------------------------------------------
# §3.6: the icdf fallback contract (W2.4's own acceptance criterion)
# ---------------------------------------------------------------------------


def gamma_beta_declaration() -> ParameterSet:
    """Two families torch implements no ``icdf`` for, and one it does."""
    return ParameterSet(
        [
            Parameter("shape", st.gamma(2.0, 0.0, 3.0)),
            Parameter("fraction", st.beta(2.0, 3.0)),
            Parameter("offset", st.norm(0.0, 1.0)),
        ]
    )


class TestIcdfFallback:
    """``lowering.md`` §3.6, ruled 2026-09-03 and pinned at the freeze."""

    def test_the_warning_names_the_families_and_the_backend(self) -> None:
        with pytest.warns(LoweringFallbackWarning) as recorded:
            TorchParameterSpace(gamma_beta_declaration())
        assert len(recorded) == 1
        message = str(recorded[0].message)
        assert "gamma" in message
        assert "beta" in message
        assert BACKEND in message
        # The family torch *can* invert must not be named: the message is what
        # tells the user which priors to change.
        assert "norm" not in message.replace("lognorm", "")

    def test_the_warning_is_once_per_lowering_not_once_per_call(self) -> None:
        """ "The warning is per run […] because the fallback decision is made
        once at lowering time, before any sampling." A per-call warning would
        bury the output of a nested-sampling run under a million copies.
        """
        with pytest.warns(LoweringFallbackWarning):
            space = TorchParameterSpace(gamma_beta_declaration())
        with warnings.catch_warnings(record=True) as later:
            warnings.simplefilter("always")
            for _ in range(5):
                space.prior_transform(np.array([0.3, 0.6, 0.2]))
        assert [w for w in later if issubclass(w.category, LoweringFallbackWarning)] == []

    def test_strict_raises_instead_naming_the_same_families(self) -> None:
        with pytest.raises(LoweringError) as raised:
            TorchParameterSpace(gamma_beta_declaration(), strict=True)
        message = str(raised.value)
        assert "gamma" in message
        assert "beta" in message
        assert raised.value.backend == BACKEND
        assert "strict" in message

    def test_the_fallback_computes_the_same_quantity(self) -> None:
        """Why the fallback is legitimate at all: it is the *same* function.

        §3.6's argument is that ``prior_transform`` has one mathematical
        definition, so computing it in numpy changes nothing about the
        posterior — which is only true if the two paths agree, and that is what
        this row asserts, against ``ParameterSet.prior_transform``.
        """
        declaration = gamma_beta_declaration()
        with pytest.warns(LoweringFallbackWarning):
            space = TorchParameterSpace(declaration)
        cube = np.array([0.17, 0.43, 0.81])
        assert space.prior_transform(cube) == pytest.approx(
            declaration.prior_transform(cube), abs=1e-9
        )

    def test_a_fully_invertible_declaration_warns_about_nothing(self) -> None:
        declaration = ParameterSet(
            [
                Parameter("a", st.norm(0.0, 1.0)),
                Parameter("b", st.uniform(1.0, 2.0)),
                Parameter("c", st.halfnorm(0.0, 2.0)),
                Parameter("d", st.loguniform(0.1, 10.0)),
                Parameter("e", st.lognorm(0.5, 0.0, 2.0)),
                Parameter("f", st.expon(0.0, 2.0)),
            ]
        )
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            space = TorchParameterSpace(declaration, strict=True)
        assert [w for w in recorded if issubclass(w.category, LoweringFallbackWarning)] == []
        assert space.icdf_fallback_families == frozenset()


# ---------------------------------------------------------------------------
# §5 to §8: the declaration forms
# ---------------------------------------------------------------------------


def mixed_declaration() -> ParameterSet:
    """One free scalar of each bijection kind, one array, one fixed."""
    return ParameterSet(
        [
            Parameter("temperature", st.uniform(100.0, 9900.0)),
            Parameter("log_tau", st.norm(0.0, 1.0)),
            Parameter("width", st.halfnorm(0.0, 2.0)),
            Parameter("offset", st.norm(0.0, 0.05), shape=(3,)),
            Parameter("distance", value=1.5, fixed=True),
        ]
    )


class TestModuleConventions:
    """§6.1 and §5.2: what the lowered ``nn.Module`` looks like."""

    def test_state_dict_keys_are_the_merged_names(self) -> None:
        component = ParameterSet(
            [
                Parameter("temperature", st.uniform(100.0, 9900.0)),
                Parameter("distance", st.norm(1.5, 0.1), shared_as="distance"),
            ]
        )
        mapping = ParameterSet.merge({"sed": component, "spectrum": component})
        space = TorchParameterSpace(mapping.merged)
        assert set(space.module.state_dict()) == {
            "sed.temperature",
            "spectrum.temperature",
            "distance",
        }

    def test_a_tie_label_lands_on_the_root_module(self) -> None:
        """§12.3, settled: a tie label is unqualified, so it belongs to no component."""
        component = ParameterSet(
            [
                Parameter("temperature", st.uniform(100.0, 9900.0)),
                Parameter("distance", st.norm(1.5, 0.1), shared_as="distance"),
            ]
        )
        mapping = ParameterSet.merge({"sed": component, "spectrum": component})
        space = TorchParameterSpace(mapping.merged)
        assert "distance" in dict(space.module.named_parameters(recurse=False))

    def test_free_parameters_are_trainable_and_fixed_ones_are_buffers(self) -> None:
        """§5.2 and §11 row 8: no gradient flows to a fixed parameter."""
        space = TorchParameterSpace(mixed_declaration())
        named = dict(space.module.named_parameters())
        buffers = dict(space.module.named_buffers())
        assert set(named) == {"temperature", "log_tau", "width", "offset"}
        assert set(buffers) == {"distance"}
        assert all(parameter.requires_grad for parameter in named.values())
        assert not buffers["distance"].requires_grad

    def test_no_gradient_reaches_a_buffer(self) -> None:
        space = TorchParameterSpace(mixed_declaration())
        theta = torch.tensor(
            space.declaration.prior_transform(np.full(6, 0.5)),
            dtype=torch.float64,
            requires_grad=True,
        )
        space.log_prior_tensor(theta).backward()
        assert space.module.get_buffer("distance").grad is None

    def test_every_tensor_is_float64(self) -> None:
        """§10.1's policy, and the reason it is threaded rather than global."""
        space = TorchParameterSpace(mixed_declaration())
        for tensor_ in list(space.module.parameters()) + list(space.module.buffers()):
            assert tensor_.dtype is torch.float64

    def test_the_global_default_dtype_is_untouched(self) -> None:
        """ "Never call ``torch.set_default_dtype``" — a library that does
        changes the numerical behaviour of every other consumer of torch in the
        same interpreter.
        """
        before = torch.get_default_dtype()
        TorchParameterSpace(mixed_declaration()).prior_transform(np.full(6, 0.5))
        assert torch.get_default_dtype() is before

    def test_an_array_valued_parameter_is_one_tensor_of_the_declared_shape(self) -> None:
        space = TorchParameterSpace(mixed_declaration())
        assert tuple(space.module.get_parameter("offset").shape) == (3,)


class TestFlatLayout:
    """§11 row 9: fixed parameters take no sampler dimension."""

    def test_free_size_agrees_with_the_declaration(self) -> None:
        declaration = mixed_declaration()
        space = TorchParameterSpace(declaration)
        assert space.free_size == declaration.free_size == 6
        trainable = sum(p.numel() for p in space.module.parameters() if p.requires_grad)
        assert trainable == space.free_size

    def test_a_fixed_parameter_is_absent_from_the_free_labels(self) -> None:
        space = TorchParameterSpace(mixed_declaration())
        assert "distance" not in space.free_labels()
        assert "distance" in space.unpack(np.zeros(6))


class TestEntryAssertion:
    """§1.6: the three standing preconditions, checked once, at entry."""

    def test_a_deferred_parameter_is_refused(self) -> None:
        declaration = ParameterSet([Parameter("distance", shared_as="distance")])
        with pytest.raises(ParameterError, match="deferred"):
            TorchParameterSpace(declaration)

    def test_an_unmerged_tie_label_is_refused(self) -> None:
        declaration = ParameterSet([Parameter("distance", st.norm(1.5, 0.1), shared_as="distance")])
        with pytest.raises(ParameterError, match="tie label"):
            TorchParameterSpace(declaration)


# ---------------------------------------------------------------------------
# §8: plates and hierarchical priors
# ---------------------------------------------------------------------------


def plated_declaration() -> ParameterSet:
    plate = Plate(
        "objects",
        size=4,
        hyperparameters=[
            Parameter("mu", st.norm(0.0, 5.0)),
            Parameter("sigma", st.halfnorm(0.0, 2.0)),
        ],
        members=[Parameter("theta", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
    )
    return ParameterSet(plate.expand())


class TestPlatesAndHierarchy:
    """§8, worked through end to end on the shape the spec uses as its example."""

    def test_a_plate_lowers_to_one_tensor_of_the_plate_size(self) -> None:
        """§5.1's accepted asymmetry: torch has no plate primitive."""
        space = TorchParameterSpace(plated_declaration())
        assert tuple(space.module.get_parameter("objects.theta").shape) == (4,)
        assert space.free_size == 6

    def test_the_joint_density_matches_the_reference_path(self) -> None:
        declaration = plated_declaration()
        space = TorchParameterSpace(declaration)
        theta = declaration.prior_transform(np.linspace(0.2, 0.8, declaration.free_size))
        assert space.lnprior(theta) == pytest.approx(declaration.lnprior(theta), abs=1e-10)

    def test_the_member_density_tracks_its_hyperparameters(self) -> None:
        """§11 row 7: the check that catches a hoisted distribution object.

        A lowering that built ``Normal(mu, sigma)`` once, at lowering time,
        would freeze the prior at its initial hyperparameters and keep
        returning the same member density as ``mu`` moved — a plausible-looking
        wrong answer rather than a crash.
        """
        declaration = plated_declaration()
        space = TorchParameterSpace(declaration)
        base = {"objects.mu": 0.0, "objects.sigma": 1.0, "objects.theta": np.full(4, 0.5)}
        shifted = dict(base, **{"objects.mu": 2.0})
        assert space.lnprior(base) != pytest.approx(space.lnprior(shifted))
        assert space.lnprior(shifted) == pytest.approx(declaration.lnprior(shifted), abs=1e-10)

    def test_the_hyperparameters_take_a_gradient_through_the_member(self) -> None:
        """Why the hierarchical path builds from tensors rather than floats."""
        declaration = plated_declaration()
        space = TorchParameterSpace(declaration)
        # Deliberately away from the prior median. At the median every member
        # sits exactly at ``mu``, so d/dmu of the member density is zero by
        # symmetry and the row would pass on a lowering that had severed the
        # graph entirely.
        theta = torch.tensor(
            declaration.pack(
                {
                    "objects.mu": 0.4,
                    "objects.sigma": 1.2,
                    "objects.theta": np.array([-0.3, 0.9, 1.4, 0.1]),
                }
            ),
            dtype=torch.float64,
            requires_grad=True,
        )
        space.log_prior_tensor(theta).backward()
        assert theta.grad is not None
        mu = declaration.free_slice("objects.mu")
        # The closed form: d/dmu [log N(mu; 0, 5) + sum_i log N(theta_i; mu, sigma)]
        # = -mu / 25 + sum_i (theta_i - mu) / sigma**2.
        members = np.array([-0.3, 0.9, 1.4, 0.1])
        expected = -0.4 / 25.0 + float(np.sum(members - 0.4)) / 1.2**2
        assert float(theta.grad[mu][0]) == pytest.approx(expected, abs=1e-9)

    def test_the_prior_transform_agrees_with_the_reference_path(self) -> None:
        """A hierarchical member's quantile is a function of its hyperparameters.

        The transform resolves in topological order, so ``objects.theta``'s
        distribution is built from the ``mu`` and ``sigma`` this same cube
        produced — the same rule ``ParameterSet.prior_transform`` follows, and
        the reason it is compared against that rather than against a fixed
        family.
        """
        declaration = plated_declaration()
        space = TorchParameterSpace(declaration)
        cube = np.linspace(0.13, 0.87, declaration.free_size)
        assert space.prior_transform(cube) == pytest.approx(
            declaration.prior_transform(cube), abs=1e-9
        )

    def test_a_hierarchical_normal_needs_no_reference_fallback(self) -> None:
        """§3.6, applied to the family rather than to a parametrisation.

        ``Normal.icdf`` exists, so a hierarchical normal is inverted natively
        and nothing is warned about — which is only true if the availability
        check lowers the *family* at construction rather than waiting for the
        first call, when §3.6's decision has already had to be made.
        """
        with warnings.catch_warnings(record=True) as recorded:
            warnings.simplefilter("always")
            space = TorchParameterSpace(plated_declaration(), strict=True)
        assert [w for w in recorded if issubclass(w.category, LoweringFallbackWarning)] == []
        assert space.icdf_fallback_families == frozenset()

    def test_an_unresolvable_reference_is_refused_by_name(self) -> None:
        prior = HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"})
        with pytest.raises(LoweringError, match="not resolved yet"):
            lower_hierarchical(prior, {"mu": tensor(0.0)}, parameter="theta")

    def test_a_family_that_cannot_take_tensor_arguments_is_refused(self) -> None:
        prior = HierarchicalPrior("lognorm", {"loc": "mu"}, kwds={"s": 0.5})
        with pytest.raises(LoweringError, match="location-scale"):
            lower_hierarchical(prior, {"mu": tensor(0.0)}, parameter="theta")


# ---------------------------------------------------------------------------
# §9: RNG
# ---------------------------------------------------------------------------


class TestRandomStreams:
    """§9.2 and §11 row 10: named sub-streams, reproducible within one backend."""

    def test_the_stream_integer_is_the_shared_derivation(self) -> None:
        """The *derivation* is shared across backends; the mechanism is not."""
        assert seed_for(20260907, "prior") == substream(20260907, "prior")

    def test_the_same_seed_and_label_give_identical_draws(self) -> None:
        first = torch.rand(5, generator=generator(20260907, "prior"), dtype=torch.float64)
        second = torch.rand(5, generator=generator(20260907, "prior"), dtype=torch.float64)
        assert torch.equal(first, second)

    def test_different_labels_do_not_share_a_stream(self) -> None:
        """Adding a diagnostic must not silently change a fit's initialisation."""
        prior = torch.rand(5, generator=generator(20260907, "prior"), dtype=torch.float64)
        init = torch.rand(5, generator=generator(20260907, "initialisation"), dtype=torch.float64)
        assert not torch.equal(prior, init)

    def test_an_unseeded_run_is_honestly_not_reproducible(self) -> None:
        assert seed_for(None, "prior") != seed_for(None, "prior")

    def test_prior_draws_repeat_under_the_same_seed(self) -> None:
        space = TorchParameterSpace(mixed_declaration())
        first = space.sample(20260907)
        second = space.sample(20260907)
        for name, value in first.items():
            assert np.allclose(value, second[name])

    def test_prior_draws_land_inside_the_support(self) -> None:
        declaration = mixed_declaration()
        space = TorchParameterSpace(declaration)
        drawn = space.sample(20260907)
        assert math.isfinite(declaration.lnprior(drawn))

    def test_the_global_torch_seed_is_not_disturbed(self) -> None:
        """§9.1 route (2): a ``Generator``-scoped stream touches nothing global."""
        torch.manual_seed(1234)
        before = torch.rand(3, dtype=torch.float64)
        torch.manual_seed(1234)
        TorchParameterSpace(mixed_declaration()).sample(20260907)
        after = torch.rand(3, dtype=torch.float64)
        assert torch.equal(before, after)


# ---------------------------------------------------------------------------
# The registry itself
# ---------------------------------------------------------------------------


class TestRegistration:
    """``lowering.md`` §12.8, landed W2.6 and consumed here."""

    @pytest.mark.parametrize(
        "family",
        [
            "norm",
            "uniform",
            "halfnorm",
            "loguniform",
            "poisson",
            "lognorm",
            "expon",
            "gamma",
            "beta",
        ],
    )
    def test_every_table_row_is_registered_as_a_builtin(self, family: str) -> None:
        resolution = lookup_lowering(family, BACKEND)
        assert resolution.builtin
        assert resolution.backend == BACKEND
        assert resolution.constructor_module.startswith("ampere.backends.torch")

    @pytest.mark.parametrize("bijection", [Identity, Log, Logit])
    def test_every_bijection_row_is_registered_as_a_builtin(self, bijection: type) -> None:
        resolution = lookup_bijection_lowering(bijection, BACKEND)
        assert resolution.builtin
        assert resolution.backend == BACKEND

    def test_a_user_row_wins_over_the_builtin_one(self) -> None:
        """The registry is the authority, not a private dict in the backend."""
        from ampere.core.lowering import register_lowering
        from torch.distributions import Normal

        def wrong(spec):
            return Normal(tensor(spec.kwds["loc"] + 100.0), tensor(spec.kwds["scale"]))

        original = lookup_lowering("norm", BACKEND)
        register_lowering("norm", BACKEND, wrong, override=True)
        try:
            lowered = lower_prior(describe_prior(st.norm(0.0, 1.0)))
            assert float(lowered.distribution.mean) == pytest.approx(100.0)
        finally:
            register_lowering("norm", BACKEND, original.constructor, override=True, builtin=True)
        assert float(
            lower_prior(describe_prior(st.norm(0.0, 1.0))).distribution.mean
        ) == pytest.approx(0.0)


# ---------------------------------------------------------------------------
# The widened realisation (W2.4 slice 2)
# ---------------------------------------------------------------------------

GRID = np.geomspace(1.0, 10.0, 24)
_TRUTH = 2.0 * GRID**-1.0
_VALUES = _TRUTH + np.random.default_rng(7).normal(0.0, 0.05, GRID.size)
_UNCERTAINTY = np.full(GRID.size, 0.05)


def observed_spectrum(values: np.ndarray | None = None, uncertainty: bool = True) -> Spectrum:
    return Spectrum(
        GRID * u.um,
        (_VALUES if values is None else values) * u.Jy,
        uncertainty=(_UNCERTAINTY * u.Jy) if uncertainty else None,
    )


def power_law() -> PowerLaw:
    return PowerLaw(GRID, norm=st.lognorm(0.3, scale=2.0), index=st.norm(-1.0, 0.3))


def lowered_for(likelihood: Likelihood, observed: Spectrum | None = None) -> Any:
    problem = FittingProblem(
        power_law(),
        [Dataset(observed if observed is not None else observed_spectrum(), likelihood=likelihood)],
        seed=11,
    )
    return problem, lower_problem(problem)


def worst_disagreement(problem: FittingProblem, lowered: Any, points: int = 24) -> float:
    """The largest relative disagreement with the numpy path over *points* draws."""
    rng = np.random.default_rng(4)
    worst = 0.0
    for _ in range(points):
        y = rng.normal(0.0, 1.2, lowered.free_size)
        expected = problem.log_prob_unconstrained(y)
        got = float(lowered.log_prob_unconstrained(y))
        if not math.isfinite(expected) and not math.isfinite(got):
            continue
        worst = max(worst, abs(got - expected) / max(1.0, abs(expected)))
    return worst


def gp_case(size: int = 30) -> tuple[CoreMatern32, np.ndarray, np.ndarray, np.ndarray]:
    """A small, well-conditioned 1-D GP problem. The same shape ``test_torch_backend.py`` uses."""
    rng = np.random.default_rng(20260907)
    coordinates = np.sort(rng.uniform(0.0, 12.0, size))
    return (
        CoreMatern32(0.4, 2.0),
        coordinates,
        rng.normal(0.0, 0.3, size),
        np.full(size, 0.04),
    )


CENSORING = Censoring(
    np.where(
        np.arange(GRID.size) % 5 == 0,
        int(LimitKind.UPPER_LIMIT),
        np.where(np.arange(GRID.size) % 7 == 2, int(LimitKind.LOWER_LIMIT), 0),
    )
)

_COUNTS = np.round(np.clip(_VALUES * 10.0, 0.5, None))


class TestTheWidenedRealisation:
    """Everything slice 2 added to the density, against the numpy oracle.

    ``ampere.core`` is the definition of every quantity here, so each row is
    the same assertion — the realised density agrees with
    ``FittingProblem.log_prob_unconstrained`` at many points — applied to a
    combination W2.13 refused. The tolerance is an accumulation tolerance,
    not a modelling one: these are two transcriptions of one closed form, so
    anything above 1e-12 means they are not the same formula.
    """

    @pytest.mark.parametrize(
        ("name", "build"),
        [
            ("student_t", lambda: Likelihood(StudentTFamily(4.0), IndependentNoise())),
            (
                "student_t-fitted-nu",
                lambda: Likelihood(StudentTFamily(st.loguniform(2.0, 30.0)), IndependentNoise()),
            ),
            ("cauchy", lambda: Likelihood(CauchyFamily(), IndependentNoise())),
        ],
    )
    def test_each_non_gaussian_family_agrees_with_the_numpy_path(
        self, name: str, build: Any
    ) -> None:
        problem, lowered = lowered_for(build())
        assert worst_disagreement(problem, lowered) < 1e-12

    def test_the_poisson_family_agrees_with_the_numpy_path(self) -> None:
        """Counts, and no uncertainties at all — ``REQUIRES_UNCERTAINTY`` is False."""
        problem, lowered = lowered_for(
            Likelihood(PoissonFamily(), IndependentNoise()),
            observed=Spectrum(GRID * u.um, _COUNTS * u.Jy),
        )
        assert worst_disagreement(problem, lowered) < 1e-12

    @pytest.mark.parametrize("solver", ["dense", "quasisep"], ids=["dense", "quasisep"])
    def test_the_latent_gp_poisson_path_agrees_with_the_numpy_path(self, solver: str) -> None:
        """W2.14: the combination slice 2 refused, now held to the fixed oracle.

        ``DEVELOPMENT_PLAN.md`` §4.4's singled-out case — ``counts ~
        Poisson(rate · exp(f))``, ``f ~ GP`` — and the one no gradient-free
        engine can run, because the latent block is a sampler dimension per
        retained sample. Both solves, because the whitening transform is where
        the two representations differ most: a dense Cholesky against
        celerite's ``L √D``.
        """
        chosen = DenseGP() if solver == "dense" else QuasisepGP()
        problem, lowered = lowered_for(
            Likelihood(
                PoissonFamily(),
                GaussianProcessNoise(
                    Matern32(st.halfnorm(0.0, 1.0), st.loguniform(0.5, 5.0)), chosen
                ),
            ),
            observed=Spectrum(GRID * u.um, _COUNTS * u.Jy),
        )
        assert problem.free_size > GRID.size
        assert worst_disagreement(problem, lowered) < 1e-9

    @pytest.mark.parametrize("solver", ["dense", "quasisep"], ids=["dense", "quasisep"])
    def test_the_latent_path_gives_the_kernel_hyperparameters_a_gradient(self, solver: str) -> None:
        """The point of lowering it at all, asserted rather than assumed.

        Before W2.14 the *contract* path scored a latent likelihood at the
        whitened ``z``, so the amplitude and the length scale had no influence
        on the density and therefore no gradient to give. A realisation that
        merely reproduced that would run a NUTS chain that moved both
        parameters by prior draws alone. This row asks the realised density
        for the gradient and requires it to be non-zero in both.
        """
        chosen = DenseGP() if solver == "dense" else QuasisepGP()
        problem, lowered = lowered_for(
            Likelihood(
                PoissonFamily(),
                GaussianProcessNoise(
                    Matern32(st.halfnorm(0.0, 1.0), st.loguniform(0.5, 5.0)), chosen
                ),
            ),
            observed=Spectrum(GRID * u.um, _COUNTS * u.Jy),
        )
        names = list(problem.parameters.free_names)
        # Away from z = 0, where the derivative of f = L(theta) z in theta is
        # zero for the honest reason that f is zero whatever the kernel is.
        draw = np.random.default_rng(20260908).normal(0.0, 0.6, lowered.free_size)
        y = torch.as_tensor(draw, dtype=torch.float64).requires_grad_(True)
        value = lowered.log_prob_unconstrained(y)
        value.backward()
        assert y.grad is not None
        gradient = y.grad.detach().numpy()
        assert np.all(np.isfinite(gradient))
        for name in ("default.likelihood.amplitude", "default.likelihood.length_scale"):
            index = names.index(name)
            assert abs(float(gradient[index])) > 1e-6

    def test_a_non_positive_poisson_rate_is_minus_infinity_and_not_a_raise(self) -> None:
        """``inference.md`` §10a: the realised density never raises.

        ``ampere.core``'s ``PoissonFamily.log_prob`` raises a
        ``LikelihoodError`` for a rate ``<= 0`` — right on the contract path,
        where §11 turns it into a recorded failure. A realised density has no
        such route (nothing can be recorded inside a trace), so the same
        condition must arrive as a bare ``-inf``, still attached to the graph.
        The family body is exercised directly because the composed density
        cannot easily be driven to a non-positive rate: a lognormal norm
        underflows towards zero without reaching it.
        """
        inputs = FamilyInputs(
            predicted=torch.zeros(3, dtype=torch.float64, requires_grad=True),
            observed=torch.ones(3, dtype=torch.float64),
            sigma=None,
            values={},
            family=PoissonFamily(),
        )
        value = native_log_prob(inputs)
        assert value.grad_fn is not None
        assert float(value.detach()) == -math.inf

    @pytest.mark.parametrize("family", ["gaussian", "cauchy"])
    def test_censoring_agrees_with_the_numpy_path(self, family: str) -> None:
        """The Tobit form: ``log F(z)`` for an upper limit, ``log(1-F(z))`` below.

        Compared against ``scipy``'s ``logcdf``/``logsf`` through the contract
        path, which is the point — ``torch.special.log_ndtr`` and the
        ``atan2`` form of the Cauchy CDF are different implementations of the
        same two functions, chosen for their tails.
        """
        base = GaussianFamily() if family == "gaussian" else CauchyFamily()
        problem, lowered = lowered_for(Likelihood(base, IndependentNoise(), censoring=CENSORING))
        assert worst_disagreement(problem, lowered) < 1e-12

    def test_a_censored_gaussian_far_into_the_tail_keeps_its_precision(self) -> None:
        """Why ``log_ndtr`` rather than ``log(Phi(z))``.

        ``Phi(z)`` underflows float64 near ``z = -38``, and the log of the
        underflow is ``-inf`` rather than the ``-725`` it should be. A limit
        *lives* in that tail — a strongly discrepant upper limit is exactly the
        case a user needs the fit to handle — so the naive form would turn a
        merely improbable model into an impossible one and stop the sampler.
        """
        problem, lowered = lowered_for(
            Likelihood(GaussianFamily(), IndependentNoise(), censoring=CENSORING)
        )
        # A norm far above the data makes every upper limit violently violated.
        y = np.array([6.0, 0.0])
        got = float(lowered.log_prob_unconstrained(y))
        assert math.isfinite(got)
        assert got == pytest.approx(problem.log_prob_unconstrained(y), rel=1e-12)

    @pytest.mark.parametrize("solver", ["dense", "quasisep"], ids=["dense", "quasisep"])
    def test_the_gp_solvers_both_lower_and_agree_with_the_numpy_path(self, solver: str) -> None:
        problem, lowered = lowered_for(
            Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(
                    Matern32(st.halfnorm(0.0, 0.5), 2.0),
                    DenseGP() if solver == "dense" else QuasisepGP(),
                ),
            )
        )
        assert worst_disagreement(problem, lowered) < 1e-12

    def test_the_quasiseparable_density_is_differentiable_in_every_parameter(self) -> None:
        """The claim that makes the O(N) solver worth putting in a density.

        A gradient of exactly zero in the amplitude or the length scale is the
        signature of W2.4 slice 1's finding — a covariance built in numpy —
        and it would be invisible in any value-agreement row, because the
        *value* was always right.
        """
        _, lowered = lowered_for(
            Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(
                    Matern32(st.halfnorm(0.0, 0.5), st.loguniform(0.5, 10.0)), QuasisepGP()
                ),
            )
        )
        y = torch.zeros(lowered.free_size, dtype=torch.float64, requires_grad=True)
        lowered.log_prob_unconstrained(y).backward()
        assert y.grad is not None
        assert np.all(np.isfinite(y.grad.numpy()))
        assert np.all(np.abs(y.grad.numpy()) > 0.0)

    @pytest.mark.parametrize(
        ("name", "build"),
        [
            (
                "diagonal",
                lambda: Likelihood(GaussianFamily(), FractionalModelNoise(f=st.halfnorm(0.0, 0.2))),
            ),
            (
                "diagonal-with-scale-and-jitter",
                lambda: Likelihood(
                    GaussianFamily(),
                    FractionalModelNoise(
                        f=st.halfnorm(0.0, 0.2),
                        scale=st.lognorm(0.2),
                        jitter=st.halfnorm(0.0, 0.02),
                    ),
                ),
            ),
            (
                "gp-dense",
                lambda: Likelihood(
                    GaussianFamily(),
                    FractionalModelGPNoise(
                        Matern32(st.halfnorm(0.0, 0.5), 2.0), f=st.halfnorm(0.0, 0.2)
                    ),
                ),
            ),
            (
                "gp-quasiseparable",
                lambda: Likelihood(
                    GaussianFamily(),
                    FractionalModelGPNoise(
                        Matern32(st.halfnorm(0.0, 0.5), 2.0),
                        QuasisepGP(),
                        f=st.halfnorm(0.0, 0.2),
                    ),
                ),
            ),
        ],
    )
    def test_prediction_aware_noise_agrees_with_the_numpy_path(self, name: str, build: Any) -> None:
        """W2.13's carried finding, closed: ``sigma_tensor`` is live."""
        problem, lowered = lowered_for(build())
        assert worst_disagreement(problem, lowered) < 1e-12

    def test_the_model_uncertainty_fraction_takes_a_gradient(self) -> None:
        """``f`` multiplies the *prediction*, so its gradient runs through it.

        The check that ``sigma_tensor`` really is on the tensor path: a
        version that rebuilt the base sigma by calling ``NoiseModel.sigma``
        would coerce ``f`` with ``float()`` and this gradient would be exactly
        zero.
        """
        _, lowered = lowered_for(
            Likelihood(
                GaussianFamily(),
                FractionalModelNoise(f=st.halfnorm(0.0, 0.2), scale=st.lognorm(0.2)),
            )
        )
        y = torch.zeros(lowered.free_size, dtype=torch.float64, requires_grad=True)
        lowered.log_prob_unconstrained(y).backward()
        assert np.all(np.abs(y.grad.numpy()) > 0.0)

    def test_the_prediction_aware_classes_keep_the_core_names_and_declare_torch(self) -> None:
        """``Likelihood.to_spec`` records ``type(noise).__name__`` across backends."""
        for built in (
            FractionalModelNoise(f=0.1),
            FractionalModelGPNoise(Matern32(0.3, 2.0), f=0.1),
        ):
            assert type(built).__name__ in {"FractionalModelNoise", "FractionalModelGPNoise"}
            assert built.BACKEND == BACKEND
            assert built.DIFFERENTIABLE is True

    def test_a_fractional_noise_model_defaults_to_this_backends_solver(self) -> None:
        assert isinstance(FractionalModelGPNoise(Matern32(0.3, 2.0), f=0.1).solver, DenseGP)

    def test_a_hand_set_negative_fraction_is_refused_on_the_contract_path(self) -> None:
        from ampere.core import LikelihoodError

        noise = FractionalModelNoise(f=-0.1)
        observed = observed_spectrum()
        with pytest.raises(LikelihoodError, match="non-negative"):
            noise.sigma(
                observed,
                np.ones(GRID.size, dtype=bool),
                {"f": -0.1},
                predicted=np.ones(GRID.size),
            )


class TestBatchingAndPrecision:
    """W2.4 slice 2's ``BATCHABLE`` claim and the float64 policy's opt-out."""

    def test_the_batched_density_equals_evaluating_one_at_a_time(self) -> None:
        """The whole content of ``BATCHABLE``: same numbers, one call.

        Compared against the loop rather than against a closed form, because
        the claim under test is not "the density is right" — every other row
        here checks that — but "``vmap`` rewrote it without changing it".
        """
        _, lowered = lowered_for(Likelihood(GaussianFamily(), IndependentNoise()))
        stack = torch.as_tensor(
            np.random.default_rng(5).normal(0.0, 0.9, (7, lowered.free_size)),
            dtype=torch.float64,
        )
        one_at_a_time = torch.stack([lowered.log_prob_unconstrained(row) for row in stack])
        assert lowered.log_prob_unconstrained_batched(stack) == pytest.approx(
            one_at_a_time, abs=0.0
        )

    @pytest.mark.parametrize(
        ("name", "build"),
        [
            ("dense-gp", lambda: GaussianProcessNoise(Matern32(st.halfnorm(0.0, 0.5), 2.0))),
            ("fractional", lambda: FractionalModelNoise(f=st.halfnorm(0.0, 0.2))),
        ],
    )
    def test_the_batched_density_survives_the_richer_noise_models(
        self, name: str, build: Any
    ) -> None:
        _, lowered = lowered_for(Likelihood(GaussianFamily(), build()))
        stack = torch.as_tensor(
            np.random.default_rng(6).normal(0.0, 0.7, (4, lowered.free_size)),
            dtype=torch.float64,
        )
        one_at_a_time = torch.stack([lowered.log_prob_unconstrained(row) for row in stack])
        assert lowered.log_prob_unconstrained_batched(stack) == pytest.approx(
            one_at_a_time, abs=1e-12
        )

    def test_a_quasiseparable_problem_refuses_batching_by_name(self) -> None:
        """The measured price of the library choice, surfaced rather than hidden.

        celerite2's kernels are a compiled extension reached through a
        ``torch.autograd.Function``; ``vmap`` cannot rewrite either. Evaluating
        the stack one member at a time under a batched name would give the
        right numbers and a false cost model, so the solver declares
        ``BATCHABLE = False`` and this refuses.
        """
        problem, lowered = lowered_for(
            Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(Matern32(st.halfnorm(0.0, 0.5), 2.0), QuasisepGP()),
            )
        )
        assert problem.batchable is False
        stack = torch.zeros((3, lowered.free_size), dtype=torch.float64)
        with pytest.raises(LoweringError, match="QuasisepGP"):
            lowered.log_prob_unconstrained_batched(stack)

    def test_a_single_vector_is_refused_by_the_batched_surface(self) -> None:
        """Two methods, two shapes, decided by which was called."""
        _, lowered = lowered_for(Likelihood(GaussianFamily(), IndependentNoise()))
        with pytest.raises(LoweringError, match="stack"):
            lowered.log_prob_unconstrained_batched(torch.zeros(lowered.free_size))

    def test_the_prior_short_circuit_became_a_where_without_changing_the_answer(self) -> None:
        """The change that made batching possible, held to the contract path.

        ``log_prior_tensor`` used to return early on the first non-finite
        contribution; it now sums and converts through ``torch.where``. The
        value must be identical — including at a point outside the support,
        which is the only place the two could differ.
        """
        problem, lowered = lowered_for(Likelihood(GaussianFamily(), IndependentNoise()))
        for y in (np.zeros(lowered.free_size), np.full(lowered.free_size, 400.0)):
            assert float(lowered.log_prob_unconstrained(y)) == pytest.approx(
                problem.log_prob_unconstrained(y), abs=1e-9, nan_ok=False
            )

    def test_the_dense_solver_can_be_configured_into_float32(self) -> None:
        """``architecture.md`` §5's per-run opt-out, and its cost, both visible."""
        kernel, coordinates, residual, variance = gp_case()
        exact = DenseGP().log_marginal_likelihood(kernel, coordinates, residual, variance, {})
        reduced = DenseGP().configured(dtype=torch.float32)
        got = reduced.log_marginal_likelihood(kernel, coordinates, residual, variance, {})
        assert reduced.provenance_config()["dtype"] == "torch.float32"
        # Right to about single precision, and no better -- which is the point
        # of the trap DEVELOPMENT_PLAN.md §7 records. Asserted both ways so the
        # row fails if the opt-out silently did nothing.
        assert got == pytest.approx(exact, rel=1e-4)
        assert got != exact

    def test_configuring_a_solver_leaves_its_declaration_alone(self) -> None:
        """Why dtype and device are not dataclass fields.

        ``Likelihood.to_spec`` records a dataclass solver's fields as its
        config and ``results.md`` §14 requires the spec hash to agree across
        backends — so a configured solver must still declare exactly what the
        reference one does.
        """
        import dataclasses

        configured = DenseGP(jitter=1e-8).configured(dtype=torch.float32, device="cpu")
        assert [field.name for field in dataclasses.fields(configured)] == ["jitter"]
        assert configured.jitter == 1e-8
        assert DenseGP().provenance_config()["dtype"] == "torch.float64"

    def test_configuring_a_device_is_an_opt_in_never_a_detection(self) -> None:
        """CPU stays CPU unless a caller says otherwise (``architecture.md`` §5)."""
        assert DenseGP().provenance_config()["device"] == "cpu"
        assert DenseGP().configured(device="cpu").provenance_config()["device"] == "cpu"

    def test_a_nonsense_dtype_is_refused_by_name(self) -> None:
        from ampere.core import LikelihoodError

        with pytest.raises(LikelihoodError, match=r"torch\.dtype"):
            DenseGP().configured(dtype="not-a-dtype-at-all")
        with pytest.raises(LikelihoodError, match="floating-point"):
            DenseGP().configured(dtype=torch.int64)

    def test_the_quasiseparable_solver_refuses_to_be_reconfigured(self) -> None:
        """celerite2's kernels are float64 CPU; saying otherwise would be a lie."""
        from ampere.core import LikelihoodError

        assert QuasisepGP().configured() is not None
        with pytest.raises(LikelihoodError, match=r"float64 on the CPU"):
            QuasisepGP().configured(dtype=torch.float32)

    def test_the_solver_configuration_reaches_the_runs_attrs(self) -> None:
        """``ampere_solver_config``: recorded, never hashed (W2.13 fold-in 10)."""
        from ampere.results.provenance import provenance_attrs

        problem, _ = lowered_for(
            Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(
                    Matern32(st.halfnorm(0.0, 0.5), 2.0), DenseGP().configured(dtype=torch.float32)
                ),
            )
        )
        attrs = provenance_attrs(problem)
        assert "float32" in attrs["ampere_solver_config"]


class TestWhatTheRealisationStillRefuses:
    """Every remaining gap is refused **by name**, at construction, with a reason."""

    def test_a_censored_student_t_names_the_missing_library_function(self) -> None:
        with pytest.raises(LoweringError, match="incomplete beta"):
            lowered_for(Likelihood(StudentTFamily(4.0), IndependentNoise(), censoring=CENSORING))

    def test_a_latent_gp_solver_without_the_native_transform_is_refused_by_name(self) -> None:
        """The latent path's own native-surface requirement (W2.14).

        The blanket latent refusal this class used to carry is gone — it named
        a defect in ``ampere.core`` (a scoring path that applied no whitening
        transform, so the kernel hyperparameters entered a latent likelihood
        nowhere at all), and W2.14 fixed the defect. The surface requirement it
        leaves behind is real: a solver that declares this backend and supplies
        ``log_marginal_likelihood_native`` still cannot lower a *latent*
        dataset unless it also supplies ``latent_transform_native``, because a
        realised density that reached for the numpy contract surface would cut
        the gradient in exactly the hyperparameters the latent path exists to
        fit.
        """

        class NoWhitening(DenseGP):
            """This backend's dense solver, minus the one native method."""

            latent_transform_native = None

        counts = Spectrum(GRID * u.um, _COUNTS * u.Jy)
        likelihood = Likelihood(
            PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), NoWhitening())
        )
        with pytest.raises(LoweringError, match="latent_transform_native"):
            lowered_for(likelihood, observed=counts)

    def test_a_zero_uncertainty_gp_dataset_refuses_by_name_at_construction(self) -> None:
        """W3.0's finding, carried to this backend: refuse what the contract refuses.

        ``ampere.core.likelihood``'s ``GaussianProcessNoise.sigma`` (through
        ``_observed_sigma``) refuses a dataset with a retained uncertainty that
        is zero or negative — "an infinitely precise measurement, which no
        likelihood can normalise" — the moment the contract path evaluates it.
        W3.0 gave the jax realisation the construction-time twin of that
        refusal and recorded that torch had the identical gap; this closes it,
        with the same condition and the same wording. Reaching the refusal only
        through ``ampere.core.realise``'s one-point agreement check is not
        enough: a ``strict`` NUTS run evaluates every other point too, and
        would sample a density the contract refuses.
        """
        blank = Spectrum(GRID * u.um, _VALUES * u.Jy, uncertainty=np.zeros(GRID.size) * u.Jy)
        with pytest.raises(LoweringError, match="zero or negative uncertainties") as raised:
            lowered_for(
                Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0))),
                observed=blank,
            )
        assert "'default'" in str(raised.value)

    def test_a_gp_dataset_with_jitter_and_real_uncertainty_still_lowers(self) -> None:
        """The check is narrow: a legitimate jittered GP dataset is untouched."""
        problem, lowered = lowered_for(
            Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.3, 2.0), jitter=0.05))
        )
        theta = problem.prior_transform(np.array([0.4, 0.6]))
        assert float(lowered.log_likelihood(theta)) == pytest.approx(
            problem.log_likelihood(theta), abs=1e-8
        )

    def test_another_backends_gp_solver_is_refused_by_name(self) -> None:
        """A numpy solve in a differentiable problem is a backend disagreement."""
        from ampere.core import DatasetError

        with pytest.raises((LoweringError, DatasetError)):
            lowered_for(
                Likelihood(
                    GaussianFamily(),
                    GaussianProcessNoise(Matern32(0.3, 2.0), CoreQuasisepGP()),
                )
            )


class TestNativeBatchedSimulation:
    """W3.1 slice 2: ``simulate_batched`` and the CPU-only sharding rows.

    The cross-backend claims — agreement with the loop, partition
    independence, the distributional check on a native draw — are conformance
    rows, because they are claims about *every* backend. What is here is what
    is specific to this one: the surface's own shapes and refusals, and the
    single-device degenerate case of the sharding hook, which is the case CPU
    CI can exercise (``tests/gpu`` holds the rest).
    """

    def built(self, noise: Any = None) -> Any:
        return lowered_for(
            Likelihood(GaussianFamily(), noise if noise is not None else IndependentNoise())
        )

    def theta(self, problem: FittingProblem, draws: int) -> np.ndarray:
        rng = np.random.default_rng(20260909)
        return np.stack(
            [
                problem.prior_transform(row)
                for row in rng.uniform(0.05, 0.95, (draws, problem.free_size))
            ]
        )

    def test_it_agrees_with_the_contract_path_draw_by_draw(self) -> None:
        problem, lowered = self.built(GaussianProcessNoise(Matern32(0.3, 2.0)))
        theta = self.theta(problem, 5)
        produced = lowered.simulate_batched(theta)
        for index, row in enumerate(theta):
            expected = problem.simulate(row)
            for label, container in expected.predicted.items():
                assert np.allclose(
                    produced.predicted[label][index],
                    np.asarray(container.values),
                    rtol=0.0,
                    atol=1e-9,
                )
            for model, result in expected.results.items():
                for channel in result:
                    assert np.allclose(
                        produced.channels[model][channel][index],
                        np.asarray(result[channel].values),
                        rtol=0.0,
                        atol=1e-9,
                    )

    @pytest.mark.parametrize("chunk_size", [1, 7, None], ids=["chunk-1", "chunk-7", "whole"])
    def test_the_chunking_does_not_change_the_prediction(self, chunk_size: int | None) -> None:
        problem, lowered = self.built()
        theta = self.theta(problem, 9)
        whole = lowered.simulate_batched(theta)
        chunked = lowered.simulate_batched(theta, chunk_size=chunk_size)
        for label, values in whole.predicted.items():
            assert np.array_equal(values, chunked.predicted[label])

    def test_a_quasiseparable_problem_refuses_by_name(self) -> None:
        _, lowered = self.built(GaussianProcessNoise(Matern32(0.3, 2.0), QuasisepGP()))
        with pytest.raises(LoweringError, match="QuasisepGP"):
            lowered.simulate_batched(np.zeros((3, lowered.free_size)))

    def test_a_non_stack_is_refused(self) -> None:
        _, lowered = self.built()
        with pytest.raises(LoweringError, match="batch, "):
            lowered.simulate_batched(np.zeros(lowered.free_size))

    def test_the_single_device_sharder_is_the_degenerate_case(self) -> None:
        """CPU CI's half of Peter's addendum: the hook runs, and changes nothing."""
        from ampere.backends.torch import SingleDeviceSharder, available_devices

        problem, lowered = self.built()
        theta = self.theta(problem, 6)
        plain = lowered.simulate_batched(theta)
        for sharder in (SingleDeviceSharder(), SingleDeviceSharder("cpu")):
            assert len(sharder.devices()) == 1
            produced = lowered.simulate_batched(theta, sharder=sharder)
            for label, values in plain.predicted.items():
                assert np.array_equal(values, produced.predicted[label])
        assert available_devices()

    def test_the_distributed_sharder_refuses_without_a_process_group(self) -> None:
        """A refusal naming the remedy, never a collective that hangs."""
        from ampere.backends.torch import DistributedSharder

        if torch.distributed.is_available() and torch.distributed.is_initialized():
            pytest.skip("this process is already in a torch.distributed group")
        with pytest.raises(LoweringError, match="process group"):
            DistributedSharder()

    def test_a_native_draw_is_a_pure_function_of_its_seed(self) -> None:
        """What makes partition independence possible on a backend with no per-draw RNG."""
        problem, lowered = self.built()
        theta = self.theta(problem, 4)
        prediction = lowered.simulate_batched(theta)
        first = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 44])
        second = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 44])
        shuffled = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 45])
        for label, values in first.items():
            assert np.array_equal(values, second[label])
            assert np.array_equal(values[:3], shuffled[label][:3])
            assert not np.array_equal(values[3], shuffled[label][3])
