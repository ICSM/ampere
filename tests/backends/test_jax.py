"""The jax backend's own rows: what the neutral conformance battery cannot say.

``tests/conformance`` runs once per registered backend and is written in the
neutral vocabulary, so it can assert that this backend's numbers agree with the
reference oracle but it cannot assert anything about *numpyro*, *equinox* or
*jax* — no test body there may name a concrete backend. ``lowering.md`` §11
lists eleven rows the suite owes, and several of them are unavoidably
backend-specific:

============  ==========================================================
§11 row       where it lives
============  ==========================================================
1, 2          the Jacobian direction and ``lnprior_unconstrained``:
              :class:`TestTheJacobianDirection` here (argument order is a
              numpyro fact), and every ``test_parameters.py`` row there.
3, 4          parametrisation and support, family by family:
              :class:`TestTheDistributionTable` here — the neutral battery
              exercises three families, and the §3.2 table has ten.
5             an unsupported family raises: :class:`TestRefusals`.
6            a ``Plate`` lowers to a real ``numpyro.plate``, not a bare
              batched site: :class:`TestPlatesAndHierarchy`. This is a
              *structural* claim, so it cannot be checked numerically and
              therefore cannot be checked from the neutral suite at all.
7             hierarchical priors track their hyperparameters: same class.
8             buffers are not trainable and not static:
              :class:`TestBuffersAreLeavesNotStructure`.
9, 10         fixed parameters take no dimension; seed reproducibility:
              ``test_parameters.py`` there, and :class:`TestRandomness` here
              for the jax half of §9.
11            the x64 guard: :class:`TestTheX64Guard`.
============  ==========================================================

Plus the two things this backend exists for and the battery has no
declaration for: the §3.6 ``icdf`` contract (:class:`TestTheIcdfContract`) and
the native differentiable path (:class:`TestTheNativePath`).

x64 is turned on once, here, by an autouse session fixture — never at import of
anything in ``ampere`` (``lowering.md`` §10.2(a)).
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections.abc import Iterator
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

jax = pytest.importorskip("jax")
jnp = pytest.importorskip("jax.numpy")
eqx = pytest.importorskip("equinox")
numpyro = pytest.importorskip("numpyro")

from numpyro import handlers  # noqa: E402
from numpyro.distributions import constraints  # noqa: E402
from numpyro.distributions.transforms import biject_to  # noqa: E402

from ampere.backends import jax as ampere_jax  # noqa: E402
from ampere.backends.jax import (  # noqa: E402
    BACKEND,
    CalibrationScale,
    DenseGP,
    GaussianProcessNoise,
    IndependentNoise,
    LoweredParameterSet,
    LoweringFallbackWarning,
    Matern32,
    PowerLaw,
    QuasisepGP,
    Resample,
    SquaredExponential,
    SyntheticPhotometry,
    configure_x64,
    filter_spec,
    lower_bijection,
    lower_problem,
)
from ampere.backends.jax._declare import as_parameter  # noqa: E402
from ampere.backends.jax.distributions import has_native_icdf, lower_prior  # noqa: E402
from ampere.backends.jax.families import lower_family  # noqa: E402
from ampere.backends.jax.rng import fold, key  # noqa: E402
from ampere.core import (  # noqa: E402
    CauchyFamily,
    Censoring,
    ComplexGaussianFamily,
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    LimitKind,
    PoissonFamily,
    RiceFamily,
    StudentTFamily,
    HierarchicalPrior,
    Identity,
    Instrument,
    Likelihood,
    Log,
    Logit,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    Plate,
    Spectrum,
    Tie,
    VisibilitySet,
    negotiate,
)
from ampere.core.exceptions import LoweringError  # noqa: E402
from ampere.core.parameter import describe_prior  # noqa: E402
from ampere.core.rng import substream  # noqa: E402

REFERENCE_WAVELENGTH = 1.0


@pytest.fixture(autouse=True, scope="session")
def _x64() -> None:
    """Turn 64-bit mode on, once, from the suite rather than from the library.

    ``lowering.md`` §10.2(a): ampere never flips this process-global flag as an
    import side effect. The *application* does, at the top of its main file —
    and for a test suite the fixture is that place.
    """
    configure_x64()


# ---------------------------------------------------------------------------
# §11 row 11: the x64 guard
# ---------------------------------------------------------------------------


class TestTheX64Guard:
    """Constructing anything in 32-bit mode raises, naming the three remedies."""

    @pytest.fixture
    def without_x64(self) -> Iterator[None]:
        """Turn the flag off for the duration of one row, then put it back.

        Safe here only because nothing inside the row *computes*: each one
        constructs an object and expects the constructor to refuse. What jax
        does with arrays created before a flag flip is undocumented
        (``lowering.md`` §10.2(d)) and nothing in this file depends on it.
        """
        jax.config.update("jax_enable_x64", False)
        try:
            yield
        finally:
            jax.config.update("jax_enable_x64", True)

    def test_x64_is_reported_honestly(self, without_x64: None) -> None:
        assert ampere_jax.x64_enabled() is False

    @pytest.mark.parametrize(
        "build",
        [
            pytest.param(lambda: PowerLaw(np.array([1.0, 2.0])), id="model"),
            pytest.param(lambda: CalibrationScale(1.0), id="step"),
            pytest.param(lambda: Matern32(0.4, 2.0), id="kernel"),
            pytest.param(DenseGP, id="solver"),
            pytest.param(
                lambda: LoweredParameterSet(ParameterSet([Parameter("a", st.norm(0.0, 1.0))])),
                id="parameter-set",
            ),
        ],
    )
    def test_every_construction_path_refuses(self, without_x64: None, build: Any) -> None:
        with pytest.raises(RuntimeError, match="32-bit mode"):
            build()

    def test_the_message_names_the_three_remedies(self, without_x64: None) -> None:
        with pytest.raises(RuntimeError) as raised:
            DenseGP()
        message = str(raised.value)
        assert "JAX_ENABLE_X64=1" in message
        assert "configure_x64()" in message
        assert "reduced precision" in message

    def test_configure_x64_is_idempotent(self) -> None:
        """It is called from library code defensively; twice must be a no-op."""
        configure_x64()
        configure_x64()
        assert ampere_jax.x64_enabled() is True


# ---------------------------------------------------------------------------
# §11 rows 3 and 4: the distribution table
# ---------------------------------------------------------------------------


#: One case per row of ``lowering.md`` §3.2, both tiers. The parameter values
#: are chosen so that a *wrong* mapping gives a different answer: every
#: continuous family carries a non-default ``loc`` or ``scale`` somewhere, and
#: ``truncnorm`` carries both (§11 row 3 asks for that case by name, "since a
#: case with defaults passes under the wrong mapping").
FAMILY_CASES = {
    "norm": st.norm(1.0, 2.0),
    "uniform": st.uniform(0.5, 3.0),
    "halfnorm": st.halfnorm(0.0, 2.0),
    "halfnorm-shifted": st.halfnorm(3.0, 2.0),
    "loguniform": st.loguniform(0.1, 10.0),
    "truncnorm": st.truncnorm(-1.0, 2.0, loc=5.0, scale=3.0),
    "lognorm": st.lognorm(0.4, scale=2.0),
    "lognorm-shifted": st.lognorm(0.4, loc=1.5, scale=2.0),
    "expon": st.expon(0.0, 2.5),
    "expon-shifted": st.expon(1.5, 2.5),
    "gamma": st.gamma(2.0, 0.0, 3.0),
    "gamma-shifted": st.gamma(2.0, 1.5, 3.0),
    "beta": st.beta(2.0, 5.0),
    "beta-affine": st.beta(2.0, 5.0, loc=1.0, scale=4.0),
}


class TestTheDistributionTable:
    """Every §3.2 row, against the scipy original it was translated from."""

    @pytest.mark.parametrize("name", sorted(FAMILY_CASES))
    def test_the_log_density_matches_scipy(self, name: str) -> None:
        prior = FAMILY_CASES[name]
        lowered = lower_prior(describe_prior(prior))
        points = prior.ppf([0.05, 0.25, 0.5, 0.75, 0.95])
        assert np.asarray(lowered.log_prob(jnp.asarray(points))) == pytest.approx(
            prior.logpdf(points), abs=1e-9
        )

    @pytest.mark.parametrize("name", sorted(FAMILY_CASES))
    def test_the_support_matches_scipy(self, name: str) -> None:
        """The check that catches a dropped ``loc`` (§11 row 4).

        Asserted behaviourally rather than by reading bounds off the constraint
        object: numpyro spells the same support several ways (``positive`` and
        ``greater_than(0)`` are different types), and a test that pattern-matched
        on the type would be testing numpyro's registry rather than this
        lowering.
        """
        prior = FAMILY_CASES[name]
        lowered = lower_prior(describe_prior(prior))
        low, high = prior.support()
        inside = prior.ppf([0.01, 0.5, 0.99])
        assert bool(jnp.all(lowered.support(jnp.asarray(inside))))
        for outside in (low - 1.0, high + 1.0):
            if np.isfinite(outside):
                assert not bool(lowered.support(jnp.asarray(outside)))

    def test_truncnorm_standardisation_is_not_passed_through(self) -> None:
        """§3.3's one row where the *numbers* change, not merely their place.

        scipy's ``a``/``b`` are in units of ``scale`` about ``loc``, so
        ``truncnorm(a=-1, b=2, loc=5, scale=3)`` is supported on ``(2, 11)``.
        Passing ``a``/``b`` straight through as numpyro's ``low``/``high`` would
        truncate at ``(-1, 2)`` around a mean of 5 — a prior that excludes its
        own mode, and which reads as a badly behaved sampler rather than as a
        translation bug.
        """
        lowered = lower_prior(describe_prior(st.truncnorm(-1.0, 2.0, loc=5.0, scale=3.0)))
        assert bool(lowered.support(jnp.asarray(4.0)))
        assert not bool(lowered.support(jnp.asarray(1.0)))
        assert not bool(lowered.support(jnp.asarray(12.0)))

    def test_uniform_reads_its_second_argument_as_a_width(self) -> None:
        lowered = lower_prior(describe_prior(st.uniform(0.5, 3.0)))
        assert bool(lowered.support(jnp.asarray(3.4)))
        assert not bool(lowered.support(jnp.asarray(3.6)))

    def test_lognorm_takes_the_logarithm_of_the_scale(self) -> None:
        """Passing ``scale`` through unchanged is wrong by an exponential.

        scipy parametrises by the shape ``s`` and ``scale = exp(mu)``; numpyro's
        ``LogNormal`` takes ``loc = mu``. A lognormal's median is ``exp(loc)``,
        so the correct mapping puts it at ``scale`` (2.0) and passing ``scale``
        through as ``loc`` would put it at ``exp(2) ~ 7.39`` — the same shape of
        curve in the same place on the page, and a different prior.
        """
        lowered = lower_prior(describe_prior(st.lognorm(0.4, scale=2.0)))
        assert float(lowered.icdf(jnp.asarray(0.5))) == pytest.approx(2.0, abs=1e-9)


class TestRefusals:
    """§3.4 rule 2 and §11 row 5: no exact construction means *raise*."""

    def test_an_unknown_family_raises_naming_family_backend_and_parameter(self) -> None:
        with pytest.raises(LoweringError) as raised:
            lower_prior(describe_prior(st.cauchy(0.0, 1.0)), parameter="temperature")
        assert raised.value.family == "cauchy"
        assert raised.value.backend == BACKEND
        assert raised.value.parameter == "temperature"
        assert "norm" in str(raised.value)  # the table it does cover

    def test_a_shifted_loguniform_raises_rather_than_dropping_the_shift(self) -> None:
        """§3.3's rule 2: an inexpressible ``loc`` is never silently dropped."""
        with pytest.raises(LoweringError, match="LogUniform has no location"):
            lower_prior(describe_prior(st.loguniform(0.1, 10.0, loc=2.0)))

    def test_an_unmerged_set_is_refused_at_entry(self) -> None:
        """§1.2 and §1.6: no backend should ever see a tie label."""
        with pytest.raises(LoweringError, match="tie label"):
            LoweredParameterSet(
                ParameterSet([Parameter("d", st.norm(0.0, 1.0), shared_as="distance")])
            )

    def test_a_custom_bijection_is_refused_by_name(self) -> None:
        """§4's last row: a reference-path-only feature until someone registers it."""

        class Doubling:
            def constrain(self, y: Any) -> Any:
                return 2.0 * np.asarray(y)

            def unconstrain(self, x: Any) -> Any:
                return 0.5 * np.asarray(x)

            def log_abs_det_jacobian(self, y: Any) -> Any:
                return np.full_like(np.asarray(y, dtype=float), np.log(2.0))

        with pytest.raises(LoweringError, match="Doubling"):
            lower_bijection(Doubling())


# ---------------------------------------------------------------------------
# §11 row 1: the Jacobian direction
# ---------------------------------------------------------------------------


class TestTheJacobianDirection:
    """The silent-bug class ``lowering.md`` §2 exists to prevent.

    A reversed or transposed Jacobian does not crash, does not produce NaNs and
    does not change any shape or dtype. It produces a subtly wrong posterior
    that looks plausible. So the rows below are pinned at points where the two
    arguments *differ*, and one of them asserts that swapping them changes the
    answer — a test at a symmetric point would pass either way.
    """

    #: A point where ``x != y`` for all three bijections.
    Y = np.array([0.4, -1.2, 0.7])

    @pytest.mark.parametrize(
        "bijection",
        [Identity(), Log(lower=0.0), Log(lower=2.0), Logit(0.0, 1.0), Logit(100.0, 10000.0)],
        ids=["identity", "log", "log-shifted", "logit-unit", "logit-bounded"],
    )
    def test_constrain_is_the_transforms_forward_direction(self, bijection: Any) -> None:
        transform = lower_bijection(bijection)
        assert np.asarray(transform(jnp.asarray(self.Y))) == pytest.approx(
            np.asarray(bijection.constrain(self.Y)), rel=1e-12
        )

    @pytest.mark.parametrize(
        "bijection",
        [Identity(), Log(lower=0.0), Log(lower=2.0), Logit(0.0, 1.0), Logit(100.0, 10000.0)],
        ids=["identity", "log", "log-shifted", "logit-unit", "logit-bounded"],
    )
    def test_the_log_determinant_agrees_at_the_unconstrained_point(self, bijection: Any) -> None:
        transform = lower_bijection(bijection)
        constrained = transform(jnp.asarray(self.Y))
        got = ampere_jax.log_abs_det_jacobian(transform, self.Y, constrained)
        assert np.asarray(got) == pytest.approx(
            np.asarray(bijection.log_abs_det_jacobian(self.Y)), abs=1e-12
        )

    def test_the_row_above_fails_if_the_arguments_are_swapped(self) -> None:
        """Argument-order-sensitive by construction, as §11 row 1 requires."""
        bijection = Log(lower=0.0)
        transform = lower_bijection(bijection)
        constrained = transform(jnp.asarray(self.Y))
        swapped = transform.log_abs_det_jacobian(constrained, jnp.asarray(self.Y))
        correct = np.asarray(bijection.log_abs_det_jacobian(self.Y))
        assert not np.allclose(np.asarray(swapped), correct)

    def test_the_registry_is_asked_for_biject_to_not_transform_to(self) -> None:
        """§2(b'): the objects must be the ones ``biject_to`` returns."""
        assert type(lower_bijection(Log(lower=0.0))) is type(biject_to(constraints.positive))
        assert type(lower_bijection(Logit(0.0, 1.0))) is type(biject_to(constraints.unit_interval))


# ---------------------------------------------------------------------------
# §11 rows 6 and 7: plates and hierarchical priors
# ---------------------------------------------------------------------------


def offset_plate(size: int) -> Plate:
    return Plate(
        "objects",
        size=size,
        hyperparameters=[
            Parameter("mu", st.norm(0.0, 1.0)),
            Parameter("sigma", st.halfnorm(0.0, 1.0)),
        ],
        members=[Parameter("offsets", HierarchicalPrior("norm", {"loc": "mu", "scale": "sigma"}))],
    )


class TestPlatesAndHierarchy:
    """A ``Plate`` is a declaration, and it must survive as one."""

    @pytest.fixture
    def lowered(self) -> LoweredParameterSet:
        return LoweredParameterSet(ParameterSet(list(offset_plate(3).expand())))

    def test_a_plate_lowers_to_a_real_plate_not_a_bare_batched_site(
        self, lowered: LoweredParameterSet
    ) -> None:
        """§5.1 and §11 row 6, and the reason it is a *structural* assertion.

        A ``to_event`` site and a ``plate`` site produce the same joint density
        here, because ampere's array-valued priors are i.i.d. — so this cannot
        be checked numerically, and a backend that lowered a plate to a bare
        batched site would pass every other row. What differs is what the model
        *claims*: numpyro's subsampling, its hierarchical diagnostics and
        ArviZ's dimension naming all read the plate.
        """
        with handlers.seed(rng_seed=0):
            trace = handlers.trace(lowered.numpyro_model()).get_trace()
        assert trace["objects"]["type"] == "plate"
        member = trace["objects.offsets"]
        assert member["type"] == "sample"
        frames = member["cond_indep_stack"]
        assert [frame.name for frame in frames] == ["objects"]
        assert frames[0].size == 3

    def test_the_hyperparameters_are_bare_sites(self, lowered: LoweredParameterSet) -> None:
        with handlers.seed(rng_seed=0):
            trace = handlers.trace(lowered.numpyro_model()).get_trace()
        for name in ("objects.mu", "objects.sigma"):
            assert trace[name]["type"] == "sample"
            assert trace[name]["cond_indep_stack"] == []

    def test_site_names_are_the_merged_parameter_names_verbatim(
        self, lowered: LoweredParameterSet
    ) -> None:
        """§8: inventing a second naming scheme would put a translation layer
        between the sampler's output and the results schema."""
        with handlers.seed(rng_seed=0):
            trace = handlers.trace(lowered.numpyro_model()).get_trace()
        sites = {name for name, site in trace.items() if site["type"] == "sample"}
        assert sites == {"objects.mu", "objects.sigma", "objects.offsets"}

    def test_hierarchical_priors_track_their_hyperparameters(
        self, lowered: LoweredParameterSet
    ) -> None:
        """§11 row 7: the check that catches a distribution hoisted out of the loop.

        Two θ differing only in ``mu`` must give different member densities. An
        implementation that built ``Normal(mu, sigma)`` once at lowering time
        would return the same number twice — and would be wrong in a way no
        shape or dtype reveals.
        """
        base = np.array([0.0, 1.0, 0.1, 0.2, 0.3])
        moved = base.copy()
        moved[0] = 2.0  # objects.mu
        assert lowered.lnprior(base) != pytest.approx(lowered.lnprior(moved))

    def test_the_plate_agrees_with_the_reference_prior(self, lowered: LoweredParameterSet) -> None:
        declaration = lowered.declaration
        theta = declaration.prior_transform(np.array([0.3, 0.6, 0.2, 0.5, 0.8]))
        assert lowered.lnprior(theta) == pytest.approx(declaration.lnprior(theta), abs=1e-9)

    def test_an_array_valued_parameter_without_a_plate_is_one_event(self) -> None:
        """§5.1's other half: ``to_event(n)``, with ``n`` always explicit."""
        lowered = LoweredParameterSet(
            ParameterSet([Parameter("offset", st.norm(0.0, 0.05), shape=(3,))])
        )
        with handlers.seed(rng_seed=0):
            trace = handlers.trace(lowered.numpyro_model()).get_trace()
        site = trace["offset"]
        assert site["cond_indep_stack"] == []
        assert site["fn"].event_shape == (3,)
        assert np.shape(site["value"]) == (3,)


# ---------------------------------------------------------------------------
# §11 row 8: buffers are leaves, never structure
# ---------------------------------------------------------------------------


class _Toy(eqx.Module):
    """A module with one trainable leaf and one buffer, as §6.2 lowers them."""

    temperature: Any
    wavelength: Any


class TestBuffersAreLeavesNotStructure:
    """``lowering.md`` §7's trap, asserted rather than described."""

    def toy(self, grid: np.ndarray) -> _Toy:
        return _Toy(jnp.asarray(1000.0), jnp.asarray(grid))

    def test_changing_a_buffers_contents_does_not_change_the_pytree_structure(self) -> None:
        """The whole of why a buffer must not be an ``eqx.field(static=True)``.

        A static field's value goes into the flattening's auxiliary data — part
        of the treedef, which jax compares and hashes to build JIT cache keys —
        so a re-read opacity table or a renegotiated wavelength grid would be a
        *different structure* and trigger a full recompilation, while large
        constants stayed alive in the compilation cache.
        """
        first = self.toy(np.linspace(1.0, 10.0, 8))
        second = self.toy(np.linspace(1.0, 10.0, 8) + 0.5)
        assert jax.tree.structure(first) == jax.tree.structure(second)

    def test_a_model_buffer_is_an_ordinary_leaf(self) -> None:
        """The same claim about a real shipped model rather than a toy."""
        model = PowerLaw(np.linspace(1.0, 10.0, 8), reference_wavelength=REFERENCE_WAVELENGTH)
        other = PowerLaw(np.linspace(1.0, 20.0, 8), reference_wavelength=REFERENCE_WAVELENGTH)
        assert jax.tree.structure(model.grids) == jax.tree.structure(other.grids)
        assert jax.tree.leaves(model.grids)  # it *is* a leaf, not aux data

    def test_no_gradient_reaches_the_excluded_partition(self) -> None:
        """§6.2's deciding argument, made concrete.

        With an ``eqx.partition`` filter spec the buffer is literally absent
        (``None``) from the differentiated argument, so ``filter_grad`` returns
        ``None`` for it. That is *structural* safety: paramax's
        ``non_trainable`` and a bare ``stop_gradient`` both leave the leaf in
        the differentiated pytree and depend on someone remembering to unwrap,
        and forgetting silently trains the buffers.
        """
        toy = self.toy(np.linspace(1.0, 10.0, 8))
        spec = filter_spec(toy, trainable=("temperature",))
        trainable, frozen = eqx.partition(toy, spec)

        @eqx.filter_grad
        def loss(free: _Toy, fixed: _Toy) -> Any:
            model = eqx.combine(free, fixed)
            return jnp.sum(model.temperature * model.wavelength)

        gradient = loss(trainable, frozen)
        assert gradient.wavelength is None
        assert gradient.temperature is not None

    def test_the_default_filter_spec_marks_every_inexact_array(self) -> None:
        toy = self.toy(np.linspace(1.0, 10.0, 4))
        assert jax.tree.leaves(filter_spec(toy)) == [True, True]


# ---------------------------------------------------------------------------
# §3.6: the icdf contract
# ---------------------------------------------------------------------------


class TestTheIcdfContract:
    """Ruled 2026-09-03: sanctioned, but loud — and ``strict`` refuses it."""

    def test_availability_is_probed_rather_than_declared(self) -> None:
        """``hasattr`` is the wrong question, and this is why.

        numpyro declares ``icdf`` on the base ``Distribution``, so every family
        inherits the *name*. On 0.21.0 ``Gamma.icdf`` and ``Beta.icdf`` are
        genuinely implemented but delegate to TensorFlow Probability and raise
        ``ImportError`` where it is absent — so "the method is overridden" is
        not the question either. Whether the call works *in this environment*
        is. (``lowering.md`` §3.6's table says numpyro implements both, which is
        true of the source and not of a default install; recorded as a finding
        rather than silently worked around.)
        """
        assert has_native_icdf(lower_prior(describe_prior(st.norm(0.0, 1.0))))
        assert hasattr(lower_prior(describe_prior(st.gamma(2.0, scale=3.0))), "icdf")

    def test_a_family_with_a_native_icdf_uses_it(self) -> None:
        declaration = ParameterSet(
            [
                Parameter("t", st.uniform(100.0, 9900.0)),
                Parameter("w", st.halfnorm(0.0, 2.0)),
            ]
        )
        lowered = LoweredParameterSet(declaration)
        assert lowered.uses_reference_prior_transform is False
        cube = np.array([0.3, 0.62])
        assert lowered.prior_transform(cube) == pytest.approx(
            declaration.prior_transform(cube), abs=1e-9
        )

    def test_a_family_without_one_warns_once_naming_families_and_backend(self) -> None:
        declaration = ParameterSet([Parameter("g", st.gamma(2.0, scale=3.0))])
        with pytest.warns(LoweringFallbackWarning) as caught:
            lowered = LoweredParameterSet(declaration)
        assert len(caught) == 1, "the decision is taken once, at lowering time"
        message = str(caught[0].message)
        assert "'gamma'" in message
        assert BACKEND in message
        assert lowered.fallback_families == ("gamma",)
        assert lowered.uses_reference_prior_transform is True

    def test_the_fallback_is_exact_because_it_is_the_same_quantity(self) -> None:
        """Not §3.4's forbidden substitution: nothing about the posterior changes."""
        declaration = ParameterSet([Parameter("g", st.gamma(2.0, scale=3.0))])
        with pytest.warns(LoweringFallbackWarning):
            lowered = LoweredParameterSet(declaration)
        cube = np.array([0.37])
        assert lowered.prior_transform(cube) == pytest.approx(
            declaration.prior_transform(cube), abs=1e-12
        )

    def test_strict_refuses_it_instead(self) -> None:
        declaration = ParameterSet([Parameter("g", st.gamma(2.0, scale=3.0))])
        with pytest.raises(LoweringError) as raised:
            LoweredParameterSet(declaration, strict=True)
        assert "gamma" in str(raised.value)
        assert raised.value.backend == BACKEND

    def test_strict_comes_from_the_problem(self) -> None:
        """§3.6, pinned at the freeze: one flag, one meaning, and it is
        ``FittingProblem(strict=...)`` rather than a new per-lowering knob."""
        problem = _agreement_problem(norm=st.gamma(2.0, scale=3.0), strict=True)
        with pytest.raises(LoweringError, match="gamma"):
            lower_problem(problem)


# ---------------------------------------------------------------------------
# §9: RNG
# ---------------------------------------------------------------------------


class TestRandomness:
    """The jax half of ``lowering.md`` §9: shared derivation, native mechanism."""

    def test_the_same_seed_and_label_give_the_same_key(self) -> None:
        assert jnp.array_equal(key(20260907, "prior"), key(20260907, "prior"))

    def test_different_labels_give_different_streams(self) -> None:
        """Prior sampling and sampler initialisation must not share a stream, or
        adding a diagnostic silently changes a fit's initialisation (§9.2)."""
        assert not jnp.array_equal(key(20260907, "prior"), key(20260907, "initialisation"))

    def test_the_key_is_built_from_the_shared_derivation(self) -> None:
        assert jnp.array_equal(key(7, "prior"), jax.random.key(substream(7, "prior")))

    def test_the_key_is_typed_not_the_legacy_uint32_array(self) -> None:
        """``jax.random.key``, not ``PRNGKey``: the legacy form is untyped
        ``uint32``, carries an extra trailing axis and no RNG-implementation
        information (§9's table)."""
        assert jnp.issubdtype(key(1, "prior").dtype, jax.dtypes.prng_key)
        assert key(1, "prior").shape == ()

    def test_fold_in_takes_the_scalar_32_bit_integer_substream_produces(self) -> None:
        stream = key(1, "simulate")
        assert not jnp.array_equal(fold(stream, 0), fold(stream, 1))
        assert jnp.array_equal(fold(stream, 3), fold(stream, 3))

    def test_no_key_is_stored_on_a_lowered_object(self) -> None:
        """§9.2: "jax keys are threaded, never stored". A key on a module field
        makes the module's identity depend on RNG state, and a stale key
        silently reuses draws."""
        lowered = LoweredParameterSet(ParameterSet([Parameter("a", st.norm(0.0, 1.0))]))
        for value in vars(lowered).values():
            assert not (
                isinstance(value, jax.Array) and jnp.issubdtype(value.dtype, jax.dtypes.prng_key)
            )


# ---------------------------------------------------------------------------
# The native, differentiable path
# ---------------------------------------------------------------------------


AGREEMENT_GRID = np.geomspace(1.0, 10.0, 20)
COARSE_GRID = np.geomspace(1.5, 8.0, 7)


def _observed(grid: np.ndarray, values: np.ndarray, sigma: float, seed: int) -> Spectrum:
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.micron,
        (values + rng.normal(0.0, sigma, values.size)) * u.Jy,
        uncertainty=np.full(values.size, sigma) * u.Jy,
    )


def _power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


FINE_DATA = _observed(AGREEMENT_GRID, _power_law(AGREEMENT_GRID, 2.0, -1.2), 0.08, seed=11)
COARSE_DATA = _observed(COARSE_GRID, _power_law(COARSE_GRID, 2.0, -1.2), 0.05, seed=13)


def _agreement_problem(norm: Any = None, *, strict: bool = False) -> FittingProblem:
    """One dataset, two free parameters, on this backend's pieces."""
    return FittingProblem(
        PowerLaw(
            AGREEMENT_GRID,
            norm=st.lognorm(0.4, scale=2.0) if norm is None else norm,
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(FINE_DATA, likelihood=jax_likelihood())],
        seed=20260907,
        strict=strict,
    )


def jax_likelihood() -> Likelihood:
    """Gaussian, this backend's uncorrelated noise.

    Since W2.13 the noise model is a capability part (``inference.md`` §10a,
    fold-in 7), so ``Dataset``'s default — ``ampere.core``'s
    ``IndependentNoise``, which declares ``"reference"`` — would make every
    problem here a two-backend problem. Spelled once.
    """
    return Likelihood(GaussianFamily(), IndependentNoise())


def jax_joint_problem(seed: int | None = 20260907) -> FittingProblem:
    """``inference.md`` §15's shape, built entirely from this backend's pieces.

    The same declaration ``tests/inference/test_engines.py`` fits on the
    reference backend: one power law on two channels, observed twice through
    genuinely different instrument chains, with the two calibration factors
    tied so the problem has three free dimensions rather than four.
    """
    model = PowerLaw(
        AGREEMENT_GRID,
        norm=st.lognorm(0.4, scale=2.0),
        index=st.norm(-1.2, 0.3),
        reference_wavelength=REFERENCE_WAVELENGTH,
        channels=("blue", "red"),
    )
    return FittingProblem(
        model,
        DatasetCollection(
            {
                "blue": Dataset(
                    FINE_DATA,
                    Instrument(
                        [CalibrationScale(st.lognorm(0.05), label="calibration")], channel="blue"
                    ),
                    jax_likelihood(),
                ),
                "red": Dataset(
                    COARSE_DATA,
                    Instrument(
                        [
                            Resample(COARSE_GRID),
                            CalibrationScale(st.lognorm(0.05), label="calibration"),
                        ],
                        channel="red",
                    ),
                    jax_likelihood(),
                ),
            }
        ),
        ties=[
            Tie(
                "calibration",
                ("blue.instrument.calibration.scale", "red.instrument.calibration.scale"),
            )
        ],
        seed=seed,
    )


class TestTheNativePath:
    """``lower_problem``: the same numbers as the contract path, and a gradient."""

    def test_the_problem_reports_this_backend(self) -> None:
        problem = jax_joint_problem()
        assert problem.backend == BACKEND
        assert problem.differentiable is True

    @pytest.mark.parametrize("cube", [0.25, 0.5, 0.75])
    def test_the_log_prior_agrees_with_the_contract_path(self, cube: float) -> None:
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        theta = problem.prior_transform(np.full(problem.free_size, cube))
        assert float(lowered.log_prior(theta)) == pytest.approx(problem.log_prior(theta), abs=1e-9)

    @pytest.mark.parametrize("cube", [0.25, 0.5, 0.75])
    def test_the_log_likelihood_agrees_with_the_contract_path(self, cube: float) -> None:
        """The row that says the two surfaces compute one quantity.

        The contract path runs the models and instrument chain through
        ``ampere.core``'s containers; the native path composes ``flux`` and
        ``apply_flux`` in jax. Nothing in either is written twice on purpose —
        this is what says they did not drift.
        """
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        theta = problem.prior_transform(np.full(problem.free_size, cube))
        assert float(lowered.log_likelihood(theta)) == pytest.approx(
            problem.log_likelihood(theta), abs=1e-8
        )

    def test_the_unconstrained_density_agrees_including_the_jacobian(self) -> None:
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        y = np.array([0.3, -0.4, 0.05])
        assert float(lowered.log_prob_unconstrained(y)) == pytest.approx(
            problem.log_prob_unconstrained(y), abs=1e-8
        )

    def test_it_is_differentiable_and_the_gradient_is_finite(self) -> None:
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        y = np.array([0.3, -0.4, 0.05])
        gradient = jax.grad(lowered.log_prob_unconstrained)(jnp.asarray(y))
        assert gradient.shape == (problem.free_size,)
        assert bool(jnp.all(jnp.isfinite(gradient)))

    def test_the_gradient_matches_a_finite_difference(self) -> None:
        """An oracle that is not ampere: the density's own numerical derivative."""
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        y = np.array([0.3, -0.4, 0.05])
        analytic = np.asarray(jax.grad(lowered.log_prob_unconstrained)(jnp.asarray(y)))
        step = 1e-6
        numerical = np.empty_like(y)
        for i in range(y.size):
            plus, minus = y.copy(), y.copy()
            plus[i] += step
            minus[i] -= step
            numerical[i] = (
                problem.log_prob_unconstrained(plus) - problem.log_prob_unconstrained(minus)
            ) / (2.0 * step)
        assert analytic == pytest.approx(numerical, rel=1e-4, abs=1e-4)

    def test_it_survives_jit(self) -> None:
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        y = jnp.asarray([0.3, -0.4, 0.05])
        compiled = jax.jit(lowered.log_prob_unconstrained)
        assert float(compiled(y)) == pytest.approx(float(lowered.log_prob_unconstrained(y)))

    def test_a_gp_dataset_lowers_and_agrees(self) -> None:
        """The flexible likelihood on the native path, dense solver."""
        from ampere.core import family_named

        noise = GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP())
        problem = FittingProblem(
            PowerLaw(
                AGREEMENT_GRID,
                norm=st.lognorm(0.4, scale=2.0),
                index=st.norm(-1.2, 0.3),
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(FINE_DATA, likelihood=Likelihood(family_named("gaussian")(), noise))],
            seed=20260907,
        )
        lowered = lower_problem(problem)
        theta = problem.prior_transform(np.array([0.4, 0.6]))
        assert float(lowered.log_likelihood(theta)) == pytest.approx(
            problem.log_likelihood(theta), abs=1e-8
        )

    def test_a_zero_uncertainty_gp_dataset_refuses_by_name_at_construction(self) -> None:
        """W3.0's finding: the GP-marginal branch must refuse what the contract refuses.

        ``ampere.core.likelihood``'s ``GaussianProcessNoise.sigma`` (via
        ``_observed_sigma``) refuses a dataset with a retained uncertainty
        that is zero or negative -- "an infinitely precise measurement, which
        no likelihood can normalise" -- the moment the contract path
        evaluates it. Before this fix that refusal reached a jax problem only
        by accident, through ``ampere.core.realise``'s one-point agreement
        check (itself calling the contract path to get the reference
        log-probability); a NUTS run samples every other point too, so the
        refusal belongs on this backend's own construction path, named,
        mirroring the ``REQUIRES_UNCERTAINTY``/sigma-is-None check the non-GP
        branch already makes.
        """
        from ampere.core import family_named

        noise = GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP())
        data = Spectrum(
            AGREEMENT_GRID * u.micron,
            _power_law(AGREEMENT_GRID, 2.0, -1.2) * u.Jy,
            uncertainty=np.zeros(AGREEMENT_GRID.size) * u.Jy,
        )
        problem = FittingProblem(
            PowerLaw(
                AGREEMENT_GRID,
                norm=st.lognorm(0.4, scale=2.0),
                index=st.norm(-1.2, 0.3),
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(data, likelihood=Likelihood(family_named("gaussian")(), noise))],
            seed=20260907,
        )
        with pytest.raises(LoweringError, match="zero or negative uncertainties") as raised:
            lower_problem(problem)
        assert "'default'" in str(raised.value)

    def test_a_gp_dataset_with_jitter_and_real_uncertainty_does_not_refuse(self) -> None:
        """The check is narrow: a legitimate jittered GP dataset still lowers.

        The floor arrives through the ``jitter=`` constructor keyword, which
        this backend's ``GaussianProcessNoise`` has taken since W3.1 slice 2 —
        the parity gap W3.0 carried, since the torch class had taken it since
        W2.4 and the same three-line composition therefore raised ``TypeError``
        on one modern backend and not the other. What matters here is that
        declaring a floor alongside real, positive uncertainties does not trip
        the construction-time check.
        """
        from ampere.core import family_named

        noise = GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP(), jitter=0.1)
        problem = FittingProblem(
            PowerLaw(
                AGREEMENT_GRID,
                norm=st.lognorm(0.4, scale=2.0),
                index=st.norm(-1.2, 0.3),
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(FINE_DATA, likelihood=Likelihood(family_named("gaussian")(), noise))],
            seed=20260907,
        )
        lowered = lower_problem(problem)
        theta = problem.prior_transform(np.array([0.4, 0.6]))
        assert np.isfinite(float(lowered.log_likelihood(theta)))

    def test_a_refusal_names_what_it_cannot_lower(self) -> None:
        """Refusals happen at lowering time, never inside a trace."""
        from ampere.backends.reference import PowerLaw as ReferencePowerLaw

        problem = FittingProblem(
            ReferencePowerLaw(
                AGREEMENT_GRID,
                norm=st.lognorm(0.4, scale=2.0),
                index=-1.2,
                reference_wavelength=REFERENCE_WAVELENGTH,
            ),
            [Dataset(FINE_DATA)],
        )
        with pytest.raises(LoweringError, match="reference"):
            lower_problem(problem)

    def test_no_exception_control_flow_on_the_traced_path(self) -> None:
        """A parameter vector the model cannot be evaluated at gives ``-inf``.

        ``likelihoods.md`` §17 Q1's ruling exists because exception control flow
        cannot be traced. Under ``jit`` there is nowhere for a raise to go, so
        the failure has to arrive as a value.
        """
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        compiled = jax.jit(lowered.log_prob_unconstrained)
        assert float(compiled(jnp.asarray([1e3, 1e3, 1e3]))) == -np.inf


class TestLoweringProvenance:
    """§12.8's hardening: a user-registered row is stamped, a built-in one is not."""

    def test_a_run_on_amperes_own_table_stamps_nothing(self) -> None:
        """The signal is "did this depend on something outside the conformance
        suite's guarantees?", so ten ordinary rows must not bury it."""
        lowered = LoweredParameterSet(ParameterSet([Parameter("t", st.uniform(100.0, 9900.0))]))
        assert lowered.lowering_provenance() == []
        # ...but the rows themselves are introspectable, both slots.
        kinds = {row.kind for row in lowered.lowering_resolutions}
        assert kinds == {"prior", "bijection"}
        assert all(row.builtin for row in lowered.lowering_resolutions)

    def test_a_user_registered_row_is_stamped(self) -> None:
        from ampere.core.lowering import register_lowering

        import numpyro.distributions as npd

        def cauchy(spec: Any) -> Any:
            return npd.Cauchy(spec.kwds.get("loc", 0.0), spec.kwds.get("scale", 1.0))

        register_lowering("cauchy", BACKEND, cauchy, override=True)
        try:
            lowered = LoweredParameterSet(ParameterSet([Parameter("t", st.cauchy(0.0, 1.0))]))
            stamped = lowered.lowering_provenance()
            assert [entry["name"] for entry in stamped] == ["cauchy"]
            assert stamped[0]["backend"] == BACKEND
            assert stamped[0]["builtin"] is False
        finally:
            from ampere.core import lowering as _lowering

            _lowering._REGISTRY.pop(("prior", "cauchy", BACKEND), None)

    def test_the_lowered_problem_exposes_it_too(self) -> None:
        assert lower_problem(jax_joint_problem()).lowering_provenance() == []


class TestTheWorkedExamples:
    """The docstring examples, run as written.

    ``ampere.backends.reference`` has no doctest runner because its examples
    are prose; this package's are not — ``lower_problem``'s builds a problem,
    lowers it and differentiates it, which is the one claim in the whole
    package that is easiest to let drift into pseudocode.
    """

    @pytest.mark.parametrize(
        "module",
        ["ampere.backends.jax", "ampere.backends.jax.problem", "ampere.backends.jax.gp"],
        ids=["package", "problem", "gp"],
    )
    def test_the_examples_run(self, module: str) -> None:
        import doctest
        import importlib

        results = doctest.testmod(
            importlib.import_module(module),
            optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
            verbose=False,
        )
        assert results.failed == 0
        assert results.attempted > 0


class TestTheDenseSolver:
    """The jax Cholesky against the closed form, and its quiet-NaN guard."""

    def test_the_marginal_likelihood_matches_scipy(self) -> None:
        kernel = Matern32(0.4, 2.0)
        grid = np.linspace(1.0, 10.0, 12)
        values = kernel.resolve(None)
        covariance = np.asarray(kernel.matrix(grid, grid, values))
        sigma = np.full(grid.size, 0.1)
        residual = np.linspace(-0.3, 0.4, grid.size)
        expected = st.multivariate_normal(
            mean=np.zeros(grid.size), cov=covariance + np.diag(sigma**2)
        ).logpdf(residual)
        got = DenseGP().log_marginal_likelihood(kernel, grid, residual, sigma**2, values)
        assert got == pytest.approx(expected, abs=1e-8)

    def test_a_non_positive_definite_matrix_raises_rather_than_returning_nan(self) -> None:
        """W2.3's celerite2 finding, generalised: ``jnp.linalg.cholesky`` does not
        raise. It returns NaN, and a NaN log-likelihood in a chain is not a
        rejected proposal — it is a number that poisons every diagnostic
        downstream of it."""
        from ampere.core.exceptions import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        grid = np.array([1.0, 1.0, 1.0])  # duplicate coordinates: K is singular
        values = kernel.resolve(None)
        with pytest.raises(LikelihoodError, match="positive definite"):
            DenseGP().log_marginal_likelihood(kernel, grid, np.zeros(3), np.zeros(3), values)

    def test_the_traced_surface_returns_minus_infinity_instead(self) -> None:
        kernel = Matern32(0.4, 2.0)
        grid = np.array([1.0, 1.0, 1.0])
        values = kernel.resolve(None)
        got = jax.jit(
            lambda r: DenseGP().log_marginal_likelihood_jax(kernel, grid, r, jnp.zeros(3), values)
        )(jnp.zeros(3))
        assert float(got) == -np.inf

    def test_the_solver_is_a_frozen_dataclass_like_the_core_one(self) -> None:
        """Not decoration: ``Likelihood.describe()`` records a solver's ``config``
        only for a dataclass, and ``test_cross_backend``'s ``ampere_likelihoods``
        row compares that description across backends."""
        import dataclasses

        assert dataclasses.is_dataclass(DenseGP())
        assert [field.name for field in dataclasses.fields(DenseGP())] == ["jitter"]


# ---------------------------------------------------------------------------
# The quasiseparable solver (W2.5 slice 2)
# ---------------------------------------------------------------------------


QUASISEP_GRID = np.sort(np.random.default_rng(20260907).uniform(0.0, 40.0, 250))
QUASISEP_RESIDUAL = np.random.default_rng(5).normal(0.0, 0.3, QUASISEP_GRID.size)
QUASISEP_VARIANCE = np.full(QUASISEP_GRID.size, 0.04)


def _fitted_matern() -> Matern32:
    """A Matérn-3/2 whose hyperparameters are *free*, so a gradient has somewhere to go."""
    return Matern32(st.lognorm(0.5, scale=0.7), st.lognorm(0.5, scale=2.0))


QUASISEP_VALUES = {"amplitude": 0.7, "length_scale": 2.0}


def _quasisep_observed() -> Spectrum:
    """An ordered 1-D container for the ``check_compatible`` rows."""
    return Spectrum(
        QUASISEP_GRID * u.micron,
        QUASISEP_RESIDUAL * u.Jy,
        uncertainty=np.full(QUASISEP_GRID.size, 0.2) * u.Jy,
    )


class TestTheQuasiseparableSolver:
    """``QuasisepGP``: the same number as ``DenseGP``, in linear time, in jax.

    The neutral battery already compares the two strategies at
    ``tolerances.cross_solver`` once the fixture declares ``QUASISEP``. What it
    cannot say is anything about *celerite2* — that the import is deferred,
    that its quiet NaN is caught, that its missing batching rule is declared
    rather than discovered inside a trace, and that the representation used is
    ampere's exact rank-2 one rather than celerite2's approximate
    ``Matern32Term``. Those rows are here.
    """

    def test_the_marginal_likelihood_matches_scipy(self) -> None:
        """Against ``multivariate_normal``, not against ``DenseGP``.

        The oracle has to be outside ampere or the row proves only that two
        ampere code paths agree. This is also what says the representation is
        the *exact* Matérn-3/2: celerite2's own ``Matern32Term`` is an
        approximation and would miss by ~5e-3.
        """
        kernel = Matern32(0.7, 2.0)
        covariance = np.asarray(kernel.matrix(QUASISEP_GRID, QUASISEP_GRID, kernel.resolve(None)))
        expected = st.multivariate_normal(
            mean=np.zeros(QUASISEP_GRID.size), cov=covariance + np.diag(QUASISEP_VARIANCE)
        ).logpdf(QUASISEP_RESIDUAL)
        got = QuasisepGP().log_marginal_likelihood(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, kernel.resolve(None)
        )
        assert got == pytest.approx(expected, abs=1e-6)

    def test_it_agrees_with_the_dense_solver(self) -> None:
        kernel = _fitted_matern()
        dense = DenseGP().log_marginal_likelihood(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, QUASISEP_VALUES
        )
        quasisep = QuasisepGP().log_marginal_likelihood(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, QUASISEP_VALUES
        )
        assert quasisep == pytest.approx(dense, abs=1e-6)

    def test_the_gradient_agrees_with_the_dense_solver(self) -> None:
        """The reason this class exists at all: an O(N) solve a sampler can differentiate."""
        kernel = _fitted_matern()

        def density(solver: Any, parameters: Any) -> Any:
            return solver.log_marginal_likelihood_jax(
                kernel,
                QUASISEP_GRID,
                QUASISEP_RESIDUAL,
                QUASISEP_VARIANCE,
                {"amplitude": parameters[0], "length_scale": parameters[1]},
            )

        point = jnp.array([0.7, 2.0])
        quasisep = jax.grad(lambda p: density(QuasisepGP(), p))(point)
        dense = jax.grad(lambda p: density(DenseGP(), p))(point)
        assert np.asarray(quasisep) == pytest.approx(np.asarray(dense), abs=1e-6)
        assert np.all(np.isfinite(np.asarray(quasisep)))

    def test_unordered_coordinates_give_the_same_answer(self) -> None:
        """A Gaussian density is invariant under a simultaneous permutation.

        The solver sorts internally and undoes the permutation on the way out,
        so a caller need not know that celerite2 requires ordered coordinates.
        """
        kernel = _fitted_matern()
        order = np.random.default_rng(2).permutation(QUASISEP_GRID.size)
        solver = QuasisepGP()
        ordered = solver.log_marginal_likelihood(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, QUASISEP_VALUES
        )
        shuffled = solver.log_marginal_likelihood(
            kernel,
            QUASISEP_GRID[order],
            QUASISEP_RESIDUAL[order],
            QUASISEP_VARIANCE[order],
            QUASISEP_VALUES,
        )
        assert shuffled == pytest.approx(ordered, abs=1e-9)

    def test_conditioning_and_whitening_agree_with_the_dense_solver(self) -> None:
        kernel = _fitted_matern()
        quasisep, dense = QuasisepGP(), DenseGP()
        left = quasisep.condition(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, QUASISEP_VALUES
        )
        right = dense.condition(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, QUASISEP_VALUES
        )
        assert left.mean == pytest.approx(right.mean, abs=1e-6)
        assert left.variance == pytest.approx(right.variance, abs=1e-6)

        whitened = np.random.default_rng(9).normal(size=QUASISEP_GRID.size)
        assert quasisep.latent_transform(
            kernel, QUASISEP_GRID, whitened, QUASISEP_VALUES
        ) == pytest.approx(
            dense.latent_transform(kernel, QUASISEP_GRID, whitened, QUASISEP_VALUES), abs=1e-6
        )

    def test_a_quiet_nan_becomes_a_loud_refusal_on_the_contract_path(self) -> None:
        """celerite2 returns NaN where a Cholesky raises — measured, not assumed."""
        from ampere.core.exceptions import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        grid = np.array([1.0, 1.0, 1.0])  # duplicate coordinates: K is singular
        with pytest.raises(LikelihoodError, match="positive definite"):
            QuasisepGP().log_marginal_likelihood(
                kernel, grid, np.zeros(3), np.zeros(3), kernel.resolve(None)
            )

    def test_a_negative_variance_is_refused_before_celerite2_sees_it(self) -> None:
        from ampere.core.exceptions import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        grid = np.linspace(1.0, 4.0, 4)
        with pytest.raises(LikelihoodError, match="negative entries"):
            QuasisepGP().log_marginal_likelihood(
                kernel,
                grid,
                np.zeros(4),
                np.array([0.1, -1.0, 0.1, 0.1]),
                kernel.resolve(None),
            )

    def test_the_traced_surface_returns_minus_infinity_instead(self) -> None:
        kernel = Matern32(0.4, 2.0)
        grid = np.array([1.0, 1.0, 1.0])
        values = kernel.resolve(None)
        got = jax.jit(
            lambda r: QuasisepGP().log_marginal_likelihood_jax(
                kernel, grid, r, jnp.zeros(3), values
            )
        )(jnp.zeros(3))
        assert float(got) == -np.inf

    def test_a_negative_variance_is_minus_infinity_on_the_traced_surface(self) -> None:
        kernel = Matern32(0.4, 2.0)
        grid = np.linspace(1.0, 4.0, 4)
        values = kernel.resolve(None)
        got = QuasisepGP().log_marginal_likelihood_jax(
            kernel, grid, jnp.zeros(4), jnp.array([0.1, -1.0, 0.1, 0.1]), values
        )
        assert float(got) == -np.inf

    def test_a_non_quasiseparable_kernel_is_refused_by_the_core_check(self) -> None:
        from ampere.core.exceptions import LikelihoodError

        with pytest.raises(LikelihoodError, match="quasiseparable representation"):
            QuasisepGP().check_compatible(SquaredExponential(0.4, 2.0), _quasisep_observed())

    def test_a_kernel_with_no_builder_here_is_refused_by_name(self) -> None:
        """The table is what routes a kernel to the O(N) path, so a family with
        no entry is refused rather than silently approximated. Declaring
        ``QUASISEPARABLE`` is a claim about the *mathematics*; holding an exact
        representation is a claim about this backend, and they are different."""
        from ampere.core.exceptions import LikelihoodError

        class Undeclared(Matern32):
            FAMILY = "not_in_the_table"

        with pytest.raises(LikelihoodError, match="no exact celerite representation"):
            QuasisepGP().check_compatible(Undeclared(0.4, 2.0), _quasisep_observed())

    def test_the_leave_one_out_terms_agree_with_the_dense_closed_form(self) -> None:
        """The refusal is lifted (W2.5 slice 3), and this is what replaces it.

        ``DenseGP.conditional_loo`` forms ``(K + diag)**-1`` explicitly and
        reads its diagonal; this one accumulates the same diagonal backwards
        through the semiseparable inverse in O(N), over ``celerite2.jax.ops``'s
        own ``factor``. Nothing but agreement would show that the accumulation
        is right — a wrong ``A_ii`` is still finite, still per-sample and still
        looks exactly like a leave-one-out term.
        """
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        got = QuasisepGP().conditional_loo(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, values
        )
        expected = DenseGP().conditional_loo(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, values
        )
        assert got.shape == expected.shape
        assert got == pytest.approx(expected, abs=1e-8)

    def test_the_leave_one_out_terms_survive_unsorted_coordinates(self) -> None:
        """The permutation is undone on the way out, as it is for the density."""
        rng = np.random.default_rng(19)
        coordinates = rng.uniform(0.0, 12.0, 45)
        residual = rng.normal(0.0, 0.3, 45)
        variance = np.full(45, 0.04)
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        assert QuasisepGP().conditional_loo(
            kernel, coordinates, residual, variance, values
        ) == pytest.approx(
            DenseGP().conditional_loo(kernel, coordinates, residual, variance, values), abs=1e-8
        )

    def test_the_leave_one_out_terms_honour_the_jitter(self) -> None:
        """The jitter is part of the matrix, so both solvers must add it."""
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        assert QuasisepGP(jitter=0.05).conditional_loo(
            kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, values
        ) == pytest.approx(
            DenseGP(jitter=0.05).conditional_loo(
                kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, QUASISEP_VARIANCE, values
            ),
            abs=1e-8,
        )

    def test_a_single_sample_is_the_degenerate_case_and_is_still_right(self) -> None:
        """N = 1 has no recursion at all; the scan would have nothing to scan."""
        kernel = Matern32(0.4, 2.0)
        values = kernel.resolve(None)
        one = np.array([2.0])
        assert QuasisepGP().conditional_loo(
            kernel, one, np.array([0.3]), np.array([0.04]), values
        ) == pytest.approx(
            DenseGP().conditional_loo(kernel, one, np.array([0.3]), np.array([0.04]), values),
            abs=1e-10,
        )

    def test_a_diagonal_that_cannot_be_a_covariance_is_refused(self) -> None:
        """The contract path's preconditions apply here too."""
        from ampere.core.exceptions import LikelihoodError

        kernel = Matern32(0.4, 2.0)
        broken = QUASISEP_VARIANCE.copy()
        broken[3] = -1.0
        with pytest.raises(LikelihoodError, match="negative entries"):
            QuasisepGP().conditional_loo(
                kernel, QUASISEP_GRID, QUASISEP_RESIDUAL, broken, kernel.resolve(None)
            )

    def test_the_pointwise_group_can_now_be_emitted_under_this_solver(self) -> None:
        """What lifting the refusal actually buys (``results.md`` §6).

        ``add_pointwise_log_likelihood`` refuses by name under a solver with no
        ``conditional_loo`` and never falls back to a dense solve. Under this
        one it no longer has to.
        """
        kernel = Matern32(0.4, 2.0)
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        observed = _quasisep_observed()
        predicted = Spectrum(
            QUASISEP_GRID * u.micron,
            (np.asarray(observed.values) - QUASISEP_RESIDUAL) * u.Jy,
        )
        terms = likelihood.pointwise_log_prob(predicted, observed)
        dense = Likelihood(
            GaussianFamily(), GaussianProcessNoise(kernel, DenseGP())
        ).pointwise_log_prob(predicted, observed)
        assert terms == pytest.approx(dense, abs=1e-8)

    def test_provenance_records_which_library_produced_the_numbers(self) -> None:
        """``inference.md`` §10a fold-in 10, and the question an archive will ask.

        The numpy and jax quasiseparable paths are the same *representation*
        through different builds of the same library; nothing else in a run
        distinguishes them.
        """
        config = dict(QuasisepGP().provenance_config())
        assert config["library"] == "celerite2.jax"
        assert config["dtype"] == "float64"
        assert config["x64_policy_opt_out"] is False

    def test_it_declares_the_batching_it_cannot_do(self) -> None:
        """Not a limitation of jax: celerite2's primitives have no batching rule."""
        assert QuasisepGP().BATCHABLE is False
        assert DenseGP().BATCHABLE is True

    def test_it_is_a_frozen_dataclass_with_only_jitter_in_the_spec(self) -> None:
        """``device`` is an ``InitVar``, so it stays out of the spec hash.

        ``Likelihood.describe()`` records a dataclass solver's *fields*, and
        ``results.md`` §14 requires the spec hash to agree across backends, so
        a field this solver had and ``ampere.core.QuasisepGP`` did not would
        make two backends' declaration of one problem differ.
        """
        import dataclasses

        assert dataclasses.is_dataclass(QuasisepGP())
        assert [field.name for field in dataclasses.fields(QuasisepGP())] == ["jitter"]

    def test_celerite2_jax_is_not_imported_until_it_is_used(self) -> None:
        """``lowering.md`` §10.2(a), enforced against a dependency's side effect.

        ``celerite2.jax`` turns ``jax_enable_x64`` **on** at import. If ampere
        imported it eagerly, ``require_x64`` — the guard that is supposed to
        *refuse* when the flag is off — would be satisfied by a flag nobody in
        the user's program set. A subprocess is the only honest way to check
        it: within this session the module has long since been imported.
        """
        script = (
            "import sys\n"
            "import ampere.backends.jax\n"
            "assert 'celerite2.jax' not in sys.modules, 'imported eagerly'\n"
            "from jax import config\n"
            "assert not config.read('jax_enable_x64'), 'x64 was flipped by an import'\n"
            "print('ok')\n"
        )
        completed = subprocess.run(
            [sys.executable, "-c", script],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "JAX_ENABLE_X64": "0"},
        )
        assert completed.returncode == 0, completed.stderr
        assert "ok" in completed.stdout


class TestPrecisionAndDevice:
    """``architecture.md`` §5's float32 opt-out and the ``device=`` opt-in."""

    def test_float64_is_the_default_and_is_recorded(self) -> None:
        config = dict(DenseGP().provenance_config())
        assert config["dtype"] == "float64"
        assert config["device"] == "cpu"
        assert config["x64_policy_opt_out"] is False

    def test_the_float32_opt_out_is_recorded_as_a_departure_from_the_policy(self) -> None:
        config = dict(DenseGP(precision="float32").provenance_config())
        assert config["dtype"] == "float32"
        assert config["x64_policy_opt_out"] is True

    def test_the_float32_opt_out_actually_changes_the_arithmetic(self) -> None:
        """Not merely a label: the solve runs in single precision.

        The two answers agree to about float32's own precision and disagree by
        far more than float64's, which is exactly what an opt-out for
        throughput means and exactly why it is opt-in and recorded.
        """
        kernel = Matern32(0.4, 2.0)
        grid = np.linspace(1.0, 10.0, 40)
        residual = np.random.default_rng(4).normal(0.0, 0.2, grid.size)
        variance = np.full(grid.size, 0.05)
        values = kernel.resolve(None)
        exact = DenseGP().log_marginal_likelihood(kernel, grid, residual, variance, values)
        reduced = DenseGP(precision="float32").log_marginal_likelihood(
            kernel, grid, residual, variance, values
        )
        assert reduced == pytest.approx(exact, rel=1e-4)
        assert reduced != exact

    def test_precision_is_not_a_declaration(self) -> None:
        """It is configuration, so it must not reach the spec hash."""
        import dataclasses

        assert [f.name for f in dataclasses.fields(DenseGP(precision="float32"))] == ["jitter"]
        assert DenseGP(precision="float32") == DenseGP()

    def test_an_unknown_precision_is_refused_by_name(self) -> None:
        from ampere.core.exceptions import LikelihoodError

        with pytest.raises(LikelihoodError, match="float32"):
            DenseGP(precision="float16")

    def test_the_quasiseparable_solver_offers_no_precision_argument(self) -> None:
        """celerite2.jax is float64-only; offering the argument and ignoring it
        would be worse than not offering it."""
        with pytest.raises(TypeError):
            QuasisepGP(precision="float32")  # type: ignore[call-arg]

    def test_the_cpu_device_is_accepted_and_recorded(self) -> None:
        solver = DenseGP(device="cpu")
        assert solver.DEVICE == "cpu"
        assert dict(solver.provenance_config())["device"] == "cpu"

    def test_a_device_this_process_does_not_have_is_refused_rather_than_ignored(self) -> None:
        """Never a silent fallback: a fit asked for a GPU that quietly ran on the
        CPU is a fit whose timings mean nothing."""
        from ampere.core.exceptions import LikelihoodError

        with pytest.raises(LikelihoodError, match="no such platform"):
            DenseGP(device="definitely-not-a-platform")


class TestBatching:
    """``BATCHABLE``: ``vmap`` over the realised density, measured not asserted."""

    def test_a_vmapped_density_equals_the_same_density_in_a_loop(self) -> None:
        problem = jax_joint_problem()
        lowered = lower_problem(problem)
        assert lowered.batchable is True
        stack = np.random.default_rng(17).normal(size=(6, problem.free_size))
        batched = np.asarray(lowered.log_prob_unconstrained_batched(stack))
        loop = np.array([float(np.asarray(lowered.log_prob_unconstrained(row))) for row in stack])
        # One ulp of a log-density of order 1e4 is ~2e-12, and vmap is free to
        # accumulate in a different order; the claim is that it is the same
        # function, not that XLA reassociates identically.
        assert batched == pytest.approx(loop, rel=1e-12)

    def test_a_quasiseparable_problem_refuses_batching_by_name(self) -> None:
        """The refusal exists because the alternative is a message about a
        celerite2 primitive the user never named, raised from inside a trace."""
        problem = _gp_problem(QuasisepGP())
        lowered = lower_problem(problem)
        assert lowered.batchable is False
        with pytest.raises(LoweringError, match="QuasisepGP"):
            lowered.log_prob_unconstrained_batched(np.zeros((3, problem.free_size)))

    def test_a_non_stack_is_refused(self) -> None:
        lowered = lower_problem(jax_joint_problem())
        with pytest.raises(LoweringError, match="batch, n_dim"):
            lowered.log_prob_unconstrained_batched(np.zeros(3))


# ---------------------------------------------------------------------------
# The widened realised path (W2.5 slice 2)
# ---------------------------------------------------------------------------


def _gp_problem(solver: Any) -> FittingProblem:
    """One dataset under the flexible likelihood, on *solver*."""
    return FittingProblem(
        PowerLaw(
            AGREEMENT_GRID,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [
            Dataset(
                FINE_DATA,
                likelihood=Likelihood(
                    GaussianFamily(), GaussianProcessNoise(_fitted_matern(), solver)
                ),
            )
        ],
        seed=20260907,
    )


def _family_problem(
    family: Any,
    *,
    censoring: Any = None,
    solver: Any = None,
    counts: bool = False,
) -> FittingProblem:
    """A one-dataset problem under *family*, built entirely from this backend."""
    noise = IndependentNoise() if solver is None else GaussianProcessNoise(_fitted_matern(), solver)
    if counts:
        observed = Spectrum(
            AGREEMENT_GRID * u.micron,
            np.round(_power_law(AGREEMENT_GRID, 8.0, -0.5)) * u.Jy,
        )
        model = PowerLaw(
            AGREEMENT_GRID,
            norm=st.lognorm(0.3, scale=8.0),
            index=st.norm(-0.5, 0.2),
            reference_wavelength=REFERENCE_WAVELENGTH,
        )
    else:
        observed = FINE_DATA
        model = PowerLaw(
            AGREEMENT_GRID,
            norm=st.lognorm(0.4, scale=2.0),
            index=st.norm(-1.2, 0.3),
            reference_wavelength=REFERENCE_WAVELENGTH,
        )
    return FittingProblem(
        model,
        [Dataset(observed, likelihood=Likelihood(family, noise, censoring=censoring))],
        seed=20260907,
    )


def _limit_codes(size: int) -> Censoring:
    codes = np.zeros(size, dtype=np.int8)
    codes[3] = int(LimitKind.UPPER_LIMIT)
    codes[11] = int(LimitKind.LOWER_LIMIT)
    return Censoring(codes)


def _agrees(problem: FittingProblem, *, points: int = 12, tolerance: float = 1e-9) -> None:
    """The realised density equals the numpy contract path at many points."""
    lowered = lower_problem(problem)
    rng = np.random.default_rng(20260907)
    compared = 0
    for _ in range(points):
        y = rng.normal(0.0, 1.0, problem.free_size)
        expected = problem.log_prob_unconstrained(y)
        got = float(np.asarray(lowered.log_prob_unconstrained(y)))
        if not np.isfinite(expected):
            assert not np.isfinite(got)
            continue
        assert got == pytest.approx(expected, abs=tolerance)
        compared += 1
    assert compared > 0, "every sampled point was outside the support; the row proved nothing"


class TestTheWidenedRealisedPath:
    """Slice 1 lowered a Gaussian with dense-GP or i.i.d. noise. This is the rest.

    Every row is the same claim — the realised density agrees with
    ``FittingProblem.log_prob_unconstrained``, which is the oracle — applied to
    a declaration slice 1 refused by name. The neutral battery makes the same
    comparison for the shapes it declares; these are the shapes it does not.
    """

    def test_the_quasiseparable_solver_lowers(self) -> None:
        _agrees(_gp_problem(QuasisepGP()), tolerance=1e-6)

    def test_the_two_solvers_give_the_same_realised_density(self) -> None:
        dense, quasisep = (
            lower_problem(_gp_problem(DenseGP())),
            lower_problem(_gp_problem(QuasisepGP())),
        )
        y = np.array([0.3, -0.4, 0.2, 0.1])[: dense.free_size]
        assert float(np.asarray(quasisep.log_prob_unconstrained(y))) == pytest.approx(
            float(np.asarray(dense.log_prob_unconstrained(y))), abs=1e-6
        )

    def test_censored_gaussian_lowers(self) -> None:
        _agrees(
            _family_problem(GaussianFamily(), censoring=_limit_codes(AGREEMENT_GRID.size)),
            tolerance=1e-8,
        )

    @pytest.mark.parametrize("family", [StudentTFamily(nu=4.0), CauchyFamily()])
    def test_the_heavy_tailed_families_lower(self, family: Any) -> None:
        _agrees(_family_problem(family))

    @pytest.mark.parametrize("family", [StudentTFamily(nu=4.0), CauchyFamily()])
    def test_the_heavy_tailed_families_lower_with_censoring(self, family: Any) -> None:
        """The CDFs jax has no ``logcdf`` for, written out and checked here."""
        _agrees(_family_problem(family, censoring=_limit_codes(AGREEMENT_GRID.size)))

    def test_a_fitted_degrees_of_freedom_keeps_its_gradient(self) -> None:
        """``nu`` is an ordinary parameter, so the density must be differentiable in it."""
        problem = _family_problem(StudentTFamily(nu=st.lognorm(0.4, scale=6.0)))
        lowered = lower_problem(problem)
        gradient = np.asarray(jax.grad(lowered.log_prob_unconstrained)(np.zeros(problem.free_size)))
        assert gradient.shape == (problem.free_size,)
        assert np.all(np.isfinite(gradient))
        assert gradient[-1] != 0.0

    def test_a_censored_gradient_survives_an_exactly_zero_residual(self) -> None:
        """The NaN this backend would otherwise produce, and the reason the
        Student-t log-CDF carries a hand-written derivative.

        ``betainc``'s derivative in ``x`` diverges as ``x -> 1``, and
        ``x = nu / (nu + z**2)`` *is* 1 when a residual is exactly zero. The
        chain rule then multiplies infinity by zero, jax reports NaN, and one
        such sample poisons that term's gradient. The composite derivative is
        elementary and finite (``f(z)/F(z)``), so it is supplied rather than
        differentiated — and this row is what says so.
        """
        from ampere.backends.jax.families import lower_family as _lower

        family = StudentTFamily(nu=4.0)
        lowered = _lower(family, "sed", censored=True)
        observed = jnp.asarray([1.0, 2.0, 0.5, 3.0])
        sigma = jnp.full(4, 0.2)
        limits = jnp.asarray(
            [0, int(LimitKind.UPPER_LIMIT), int(LimitKind.LOWER_LIMIT), 0], dtype=jnp.int32
        )

        def density(predicted: Any) -> Any:
            return lowered(predicted, observed, sigma, family, {}, limits, None)

        # Every residual exactly zero: the worst case, not a random one.
        gradient = np.asarray(jax.grad(density)(observed))
        assert np.all(np.isfinite(gradient))

    def test_the_censored_student_t_derivatives_match_scipy(self) -> None:
        """The hand-written rule is checked against the definition it replaces."""
        from ampere.backends.jax.families import _student_t_logcdf, _student_t_logsf

        z = np.array([-3.0, -1.0, -0.2, 0.0, 0.5, 2.0, 4.0])
        nu = 4.0
        assert np.asarray(_student_t_logcdf(jnp.asarray(z), nu)) == pytest.approx(
            st.t.logcdf(z, nu), abs=1e-12
        )
        assert np.asarray(_student_t_logsf(jnp.asarray(z), nu)) == pytest.approx(
            st.t.logsf(z, nu), abs=1e-12
        )
        slope = np.asarray(
            jax.vmap(jax.grad(lambda value: _student_t_logcdf(value, nu)))(jnp.asarray(z))
        )
        assert slope == pytest.approx(st.t.pdf(z, nu) / st.t.cdf(z, nu), abs=1e-12)
        survival = np.asarray(
            jax.vmap(jax.grad(lambda value: _student_t_logsf(value, nu)))(jnp.asarray(z))
        )
        assert survival == pytest.approx(-st.t.pdf(z, nu) / st.t.sf(z, nu), abs=1e-12)

    def test_a_fitted_nu_under_censoring_is_refused_rather_than_faked(self) -> None:
        """jax supplies no derivative of ``betainc`` in its parameters, so there
        is no gradient in ``nu`` to be had — and a fabricated zero would be a
        fit that ran, converged and never moved ``nu``."""
        from ampere.backends.jax.families import lower_family as _lower

        with pytest.raises(LoweringError, match="no gradient in `nu`"):
            _lower(StudentTFamily(nu=st.lognorm(0.4, scale=6.0)), "sed", censored=True)

    def test_a_fitted_nu_without_censoring_is_fine(self) -> None:
        """The refusal is narrow: only the *censored* term needs the CDF."""
        from ampere.backends.jax.families import lower_family as _lower

        assert _lower(StudentTFamily(nu=st.lognorm(0.4, scale=6.0)), "sed") is not None

    def test_the_poisson_family_lowers(self) -> None:
        _agrees(_family_problem(PoissonFamily(), counts=True))

    @pytest.mark.parametrize("solver", ["dense", "quasisep"])
    def test_the_latent_gp_poisson_combination_lowers(self, solver: str) -> None:
        """``DEVELOPMENT_PLAN.md`` §4.4's singled-out case, and the one no
        gradient-free engine can run: the latent block is a sampler dimension
        per retained sample, so the problem is N + k dimensional.

        Both solves since W2.14, because the whitening transform ``f = L(θ) z``
        is now applied on this path and it is where the two representations
        differ most — a dense Cholesky against celerite's ``L √D``.
        """
        chosen = DenseGP() if solver == "dense" else QuasisepGP()
        problem = _family_problem(PoissonFamily(), solver=chosen, counts=True)
        assert problem.datasets["default"].latent is not None
        assert problem.free_size > AGREEMENT_GRID.size
        _agrees(problem, points=6)

    def test_the_latent_path_gives_the_kernel_hyperparameters_a_gradient(self) -> None:
        """The point of lowering it at all, asserted rather than assumed (W2.14).

        Until W2.14 the contract path scored a latent likelihood at the
        whitened ``z``, so the GP hyperparameters had no influence on the
        density and therefore no gradient to give, and this backend mirrored
        that deliberately. A NUTS chain on such a problem moves the amplitude
        and the length scale by prior draws alone. Evaluated away from
        ``z = 0``, where the derivative of ``f = L(θ) z`` in ``θ`` is zero for
        the honest reason that ``f`` is zero whatever ``L`` is.
        """
        problem = _family_problem(PoissonFamily(), solver=DenseGP(), counts=True)
        lowered = lower_problem(problem)
        names = list(problem.parameters.free_names)
        y = np.random.default_rng(20260908).normal(0.0, 0.6, problem.free_size)
        gradient = np.asarray(jax.grad(lowered.log_prob_unconstrained)(jnp.asarray(y)))
        assert np.all(np.isfinite(gradient))
        for name in ("default.likelihood.amplitude", "default.likelihood.length_scale"):
            assert abs(float(gradient[names.index(name)])) > 1e-6

    def test_the_terms_sum_to_the_joint_likelihood_on_a_widened_shape(self) -> None:
        problem = _family_problem(CauchyFamily(), censoring=_limit_codes(AGREEMENT_GRID.size))
        lowered = lower_problem(problem)
        y = np.zeros(problem.free_size)
        terms = lowered.log_likelihood_terms(y)
        evaluation = problem.evaluate(problem.constrain(y))
        assert float(sum(float(np.asarray(v)) for v in terms.values())) == pytest.approx(
            evaluation.log_likelihood, abs=1e-9
        )

    def test_a_family_ampere_core_does_not_implement_is_refused_by_name(self) -> None:
        """``rice`` is a declared slot whose own ``log_prob`` raises, so there is
        nothing here to agree with — and the refusal says so rather than
        pretending this backend is the limitation.

        Reached directly rather than through a problem, because ``ampere.core``
        refuses the composition *first* (``Likelihood`` will not take an
        unimplemented family at all). That is the right order and this row is
        the second line of defence: a family implemented on the reference path
        but not transcribed here must be refused by name too, and the two
        branches share the message.
        """
        with pytest.raises(LoweringError, match="rice"):
            lower_family(RiceFamily(), "sed")

    def test_another_backends_solver_is_still_refused(self) -> None:
        """W2.13's loud consequence, unchanged by the widening: a numpy solve in
        a jax problem gets no gradient at all."""
        from ampere.core import QuasisepGP as CoreQuasisepGP

        with pytest.raises(Exception, match=r"backend|QuasisepGP"):
            _gp_problem(CoreQuasisepGP())


class TestSyntheticPhotometryNatively:
    """The native photometry chain, which slice 1 shipped broken.

    ``apply_flux`` handed the inherited ``influence`` a bare array where it
    needs an :class:`~ampere.core.Axis` — the lookup is the axis's job
    (``spectrum_photometry.md`` Gap 1) — so the first evaluation of any
    realised chain ending in photometry raised ``AttributeError``. Nothing
    caught it because no conformance ``ProblemSpec`` had a photometry step;
    slice 2 adds one, and these rows pin the shipped step's own behaviour.
    """

    @staticmethod
    def _step() -> SyntheticPhotometry:
        wavelength = np.linspace(2.0, 8.0, 25)
        response = np.exp(-0.5 * ((wavelength[None, :] - np.array([[3.5], [6.0]])) / 0.6) ** 2)
        return SyntheticPhotometry(
            ["A", "B"], wavelength, response, detector="photon", label="synphot"
        )

    def test_the_native_surface_agrees_with_the_contract_surface(self) -> None:
        step = self._step()
        grid = np.linspace(2.0, 8.0, 25)
        flux = 2.0 * grid**-1.0
        spectrum = Spectrum(grid * u.micron, flux * u.Jy)
        contract = np.asarray(step.apply(spectrum, {}).values)
        native, pivots = step.apply_flux(jnp.asarray(flux), jnp.asarray(grid), {})
        assert np.asarray(native) == pytest.approx(contract, abs=1e-12)
        assert np.asarray(pivots) == pytest.approx(step.pivots(), abs=1e-12)

    def test_it_places_its_columns_by_lookup_on_a_larger_grid(self) -> None:
        """The whole point of rebuilding the axis: negotiation hands a step a
        grid larger than the one it tabulated on, and reading the values
        positionally would be silently wrong."""
        step = self._step()
        tabulated = np.linspace(2.0, 8.0, 25)
        wider = np.unique(np.concatenate([np.linspace(0.5, 12.0, 40), tabulated]))
        flux = 2.0 * wider**-1.0
        spectrum = Spectrum(wider * u.micron, flux * u.Jy)
        contract = np.asarray(step.apply(spectrum, {}).values)
        native, _ = step.apply_flux(jnp.asarray(flux), jnp.asarray(wider), {})
        assert np.asarray(native) == pytest.approx(contract, abs=1e-12)

    def test_it_is_differentiable_through_the_flux(self) -> None:
        step = self._step()
        grid = jnp.asarray(np.linspace(2.0, 8.0, 25))

        def total(scale: Any) -> Any:
            flux = scale * (2.0 * grid**-1.0)
            return jnp.sum(step.apply_flux(flux, grid, {})[0])

        gradient = float(jax.grad(total)(jnp.asarray(1.0)))
        assert np.isfinite(gradient) and gradient != 0.0


# ---------------------------------------------------------------------------
# Per-instance ``device=`` (W2.5 slice 3)
# ---------------------------------------------------------------------------


class TestPerInstanceDevice:
    """``architecture.md`` §5's placement rule, on every kind of piece.

    Slice 2 gave the two GP solvers a ``device=`` ``InitVar``. Slice 3 extends
    the same shape to the models, the instrument steps, the kernels and the
    noise models, because the flag is only useful if *every* part declares it:
    ``ampere.core.declared_capabilities`` aggregates ``DEVICE`` as "all parts
    agree or refuse", so one part that could not say where it computed would
    make that agreement unenforceable.
    """

    def _pieces(self) -> dict[str, Any]:
        grid = np.linspace(1.0, 5.0, 6)
        return {
            "PowerLaw": PowerLaw(grid),
            "Matern32": Matern32(0.4, 2.0),
            "CalibrationScale": CalibrationScale(1.0),
            "DenseGP": DenseGP(),
            "QuasisepGP": QuasisepGP(),
            "IndependentNoise": IndependentNoise(),
            "GaussianProcessNoise": GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP()),
        }

    def test_every_piece_declares_the_cpu_by_default(self) -> None:
        """Never auto-detected: the machine this runs on may well have a GPU."""
        for name, piece in self._pieces().items():
            assert piece.DEVICE == "cpu", name

    def test_the_cpu_can_be_asked_for_by_name_on_every_piece(self) -> None:
        grid = np.linspace(1.0, 5.0, 6)
        asked = [
            PowerLaw(grid, device="cpu"),
            Matern32(0.4, 2.0, device="cpu"),
            CalibrationScale(1.0, device="cpu"),
            DenseGP(device="cpu"),
            QuasisepGP(device="cpu"),
            IndependentNoise(device="cpu"),
            GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP(), device="cpu"),
        ]
        assert [piece.DEVICE for piece in asked] == ["cpu"] * len(asked)

    def test_a_device_this_process_lacks_is_refused_by_every_piece(self) -> None:
        """Never a silent fallback, and the refusal names what is present.

        The exception class follows the contract the piece belongs to — a model
        or a step is ``ampere.core.transform``'s, a kernel, a noise model or a
        solver is ``ampere.core.likelihood``'s — so the refusal reads like
        every other refusal that piece can make.
        """
        from ampere.core.exceptions import LikelihoodError, TransformationError

        grid = np.linspace(1.0, 5.0, 6)
        absent = "definitely-not-a-platform"
        for build, error in (
            (lambda: PowerLaw(grid, device=absent), TransformationError),
            (lambda: CalibrationScale(1.0, device=absent), TransformationError),
            (lambda: Matern32(0.4, 2.0, device=absent), LikelihoodError),
            (lambda: IndependentNoise(device=absent), LikelihoodError),
            (
                lambda: GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP(), device=absent),
                LikelihoodError,
            ),
            (lambda: DenseGP(device=absent), LikelihoodError),
            (lambda: QuasisepGP(device=absent), LikelihoodError),
        ):
            with pytest.raises(error, match="no such platform"):
                build()

    def test_the_model_grids_are_placed_on_the_chosen_device(self) -> None:
        """The flag is not merely a label: ``device_put`` actually ran."""
        model = PowerLaw(np.linspace(1.0, 5.0, 6), device="cpu")
        assert {d.platform for d in model.grid("default").devices()} == {"cpu"}

    def test_a_negotiated_grid_stays_on_the_chosen_device(self) -> None:
        """``compile_for`` rebuilds the jax grid, so it must place it too."""
        model = PowerLaw(np.linspace(1.0, 5.0, 6), device="cpu")
        instrument = Instrument([Resample(np.linspace(1.5, 4.5, 4))], channel="default")
        compiled = model.compile_for(negotiate([instrument]))
        assert {d.platform for d in compiled.grid("default").devices()} == {"cpu"}

    def test_a_step_places_its_influence_matrix(self) -> None:
        step = Resample(np.linspace(1.5, 4.5, 4), device="cpu")
        matrix = step._influence_jax(np.linspace(1.0, 5.0, 6))
        assert {d.platform for d in matrix.devices()} == {"cpu"}

    def test_a_kernel_places_its_covariance(self) -> None:
        grid = np.linspace(1.0, 5.0, 6)
        kernel = Matern32(0.4, 2.0, device="cpu")
        assert {d.platform for d in kernel.matrix(grid, grid, {}).devices()} == {"cpu"}
        assert {d.platform for d in kernel.diagonal(grid, {}).devices()} == {"cpu"}

    def test_a_problem_whose_pieces_disagree_about_the_device_is_refused(self) -> None:
        """The reason the flag is per instance at all.

        Faked by shadowing the flag on one solver, which is exactly what a real
        ``device="cuda"`` does on a machine that has one — the refusal is the
        same, and it happens at composition rather than inside a trace.
        """
        from ampere.core.exceptions import DatasetError

        grid = np.linspace(1.0, 5.0, 6)
        observed = Spectrum(grid * u.micron, np.ones(6) * u.Jy, uncertainty=np.full(6, 0.1) * u.Jy)
        solver = DenseGP()
        object.__setattr__(solver, "DEVICE", "cuda")
        with pytest.raises(DatasetError, match="different devices"):
            FittingProblem(
                PowerLaw(grid, norm=st.lognorm(0.4, scale=2.0)),
                [
                    Dataset(
                        observed,
                        likelihood=Likelihood(
                            GaussianFamily(),
                            GaussianProcessNoise(Matern32(0.4, 2.0), solver),
                        ),
                    )
                ],
                seed=1,
            )

    def test_an_explicit_jax_device_is_taken_as_given(self) -> None:
        """*Which* accelerator is the caller's business; ampere chooses among none."""
        device = jax.devices()[0]
        model = PowerLaw(np.linspace(1.0, 5.0, 6), device=device)
        assert model.DEVICE == device.platform

    def test_the_device_is_configuration_and_not_declaration(self) -> None:
        """It must not reach the spec hash: two runs of one declaration on two
        machines describe the same posterior. ``provenance_config`` is where it
        is recorded, and ``dataclasses.fields`` is what the hash reads."""
        import dataclasses

        assert [f.name for f in dataclasses.fields(DenseGP(device="cpu"))] == ["jitter"]
        assert DenseGP(device="cpu") == DenseGP()


# ---------------------------------------------------------------------------
# The complex Gaussian, end to end on the realised path (W2.5 slice 3)
# ---------------------------------------------------------------------------

VIS_U = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 5.5])
VIS_V = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0])


class _PointSource(Model):
    """A shifted point source, ``V(u,v) = A exp(-2 pi i s (u + v))``, in jax.

    Written here rather than shipped because ``ampere.backends.jax``'s models
    are the spectral three and the visibility modality is Phase 4's. What it
    exercises is exactly what the shipped models cannot: a **complex**
    prediction reaching the family through the realised path, with a gradient
    in both parameters.

    The parameter is ``amplitude`` and not ``flux`` because ``flux`` is the
    name of this backend's native evaluation surface, and
    ``Parameterised.register_parameter`` refuses a parameter that would shadow
    an attribute of its own class — which is the check working, not a
    limitation.
    """

    DIFFERENTIABLE = True
    BATCHABLE = True
    DEVICE = "cpu"
    BACKEND = BACKEND

    def __init__(self, u_coord: Any, v_coord: Any, *, amplitude: Any = 1.0, shift: Any = 0.0):
        self.register_buffer("u", np.asarray(u_coord, dtype=float))
        self.register_buffer("v", np.asarray(v_coord, dtype=float))
        self.register_parameter(as_parameter("amplitude", amplitude))
        self.register_parameter(as_parameter("shift", shift))

    def _visibility(self, context: Any) -> Any:
        phase = (
            -2.0
            * jnp.pi
            * jnp.asarray(context["shift"])
            * (jnp.asarray(context["u"]) + jnp.asarray(context["v"]))
        )
        return jnp.asarray(context["amplitude"]) * jnp.exp(1j * phase)

    def grid(self, channel: str) -> Any:
        return jnp.asarray(self.buffers["u"].value, dtype=jnp.float64)

    def flux(self, channel: str, values: Any = None) -> Any:
        return self._visibility(self.context(values))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(VisibilitySet(ctx["u"], ctx["v"], np.asarray(self._visibility(ctx))))


def _visibilities() -> VisibilitySet:
    exact = _PointSource(VIS_U, VIS_V, amplitude=2.0, shift=0.1).evaluate().single()
    truth = np.asarray(exact.values)
    rng = np.random.default_rng(20260908)
    noise = rng.normal(0.0, 0.05, VIS_U.size) + 1j * rng.normal(0.0, 0.05, VIS_U.size)
    return VisibilitySet(VIS_U, VIS_V, truth + noise, uncertainty=np.full(VIS_U.size, 0.05))


VISIBILITIES = _visibilities()


def _visibility_problem(noise: Any = None) -> FittingProblem:
    return FittingProblem(
        _PointSource(VIS_U, VIS_V, amplitude=st.lognorm(0.3, scale=2.0), shift=st.norm(0.1, 0.05)),
        [
            Dataset(
                VISIBILITIES,
                likelihood=Likelihood(
                    ComplexGaussianFamily(),
                    IndependentNoise() if noise is None else noise,
                ),
            )
        ],
        seed=20260908,
    )


class _GridlessPointSource(_PointSource):
    """A native model that supplies ``flux`` and forgets ``grid``.

    ``predict`` needs both, so the refusal must name both; before slice 3 it
    named only ``flux`` and a model like this died on an ``AttributeError``
    from inside the composition instead. The attribute is hidden rather than
    deleted because the inherited one is a plain method and there is no other
    way to make ``hasattr`` say no.
    """

    def __getattribute__(self, name: str) -> Any:
        if name == "grid":
            raise AttributeError(name)
        return super().__getattribute__(name)


class TestTheComplexGaussianPath:
    """``complex_gaussian`` composes into the realised density (W2.5 slice 3).

    The transcription itself landed in slice 2, in
    :mod:`ampere.backends.jax.families`. What was never shown is the thing that
    matters — that a **complex** container survives the whole realised path:
    the lowered dataset keeps ``complex128`` where every other array is
    ``float64``, the residual reaches the family as a complex difference, the
    modulus makes it real again, and a gradient comes back out. A family whose
    closed form is right and whose composition drops the imaginary part would
    pass a unit test of the closed form and fit the wrong data.
    """

    def test_the_realised_density_agrees_with_the_contract_path(self) -> None:
        _agrees(_visibility_problem(), tolerance=1e-9)

    def test_the_observed_values_keep_their_complex_dtype(self) -> None:
        """``results_schema.md``: the dtype is the declaration, so the lowering
        may not quietly make it real."""
        lowered = lower_problem(_visibility_problem())
        observed = lowered._datasets[0].observed_values
        assert observed.dtype == jnp.complex128

    def test_the_density_is_differentiable_in_both_parameters(self) -> None:
        """The whole reason the family is transcribed rather than called."""
        problem = _visibility_problem()
        lowered = lower_problem(problem)
        point = jnp.asarray(problem.unconstrain(problem.reference_values))
        gradient = np.asarray(jax.grad(lowered.log_prob_unconstrained)(point))
        assert gradient.shape == (problem.free_size,)
        assert np.all(np.isfinite(gradient))
        assert np.any(gradient != 0.0)

    def test_the_family_matches_ampere_core_at_the_same_point(self) -> None:
        """The numpy path is the oracle, here as everywhere else."""
        problem = _visibility_problem()
        theta = problem.parameters.pack(problem.reference_values)
        contract = problem.evaluate(theta)
        lowered = lower_problem(problem)
        terms = lowered.log_likelihood_terms(problem.unconstrain(theta))
        assert float(np.asarray(terms["default"])) == pytest.approx(
            contract.log_likelihood, abs=1e-9
        )

    def test_a_scaled_uncertainty_still_composes(self) -> None:
        """The noise model's ``scale`` reaches a complex family like any other."""
        _agrees(
            _visibility_problem(IndependentNoise(scale=st.lognorm(0.2, scale=1.0))),
            tolerance=1e-9,
        )

    def test_the_gp_combination_is_refused_by_name_as_the_numpy_path_refuses_it(
        self,
    ) -> None:
        """``complex_gaussian`` + GP is declared ANALYTIC and unimplemented on
        both paths (Phase 4's circular closed form), so this backend refuses
        the composition rather than inventing one."""
        from ampere.core.exceptions import LikelihoodError

        with pytest.raises(LikelihoodError):
            Likelihood(ComplexGaussianFamily(), GaussianProcessNoise(Matern32(0.4, 2.0), DenseGP()))

    def test_a_model_without_a_native_grid_is_refused_by_name(self) -> None:
        """``predict`` calls ``flux`` *and* ``grid``; the refusal now says so.

        It used to name only ``flux``, which was true of the shipped spectral
        models and misleading for anything else: a model supplying ``flux``
        alone reached ``predict`` and died on an ``AttributeError`` about
        ``grid`` from inside the composition.
        """

        problem = FittingProblem(
            _GridlessPointSource(
                VIS_U, VIS_V, amplitude=st.lognorm(0.3, scale=2.0), shift=st.norm(0.1, 0.05)
            ),
            [
                Dataset(
                    VISIBILITIES,
                    likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()),
                )
            ],
            seed=1,
        )
        with pytest.raises(LoweringError, match="grid"):
            lower_problem(problem)


class TestNativeBatchedSimulation:
    """W3.1 slice 2: ``simulate_batched`` and the CPU-only sharding rows.

    The cross-backend claims — agreement with the loop, partition
    independence, the distributional check on a native draw — are conformance
    rows, because they are claims about *every* backend. What is here is what
    is specific to this one: the surface's own shapes and refusals, and the
    single-device degenerate case of the sharding hook, which is the case CPU
    CI can exercise (``tests/gpu`` holds the multi-device rows).
    """

    def problem(self) -> FittingProblem:
        return _gp_problem(DenseGP())

    def theta(self, problem: FittingProblem, draws: int) -> np.ndarray:
        rng = np.random.default_rng(20260909)
        return np.stack(
            [
                problem.prior_transform(row)
                for row in rng.uniform(0.05, 0.95, (draws, problem.free_size))
            ]
        )

    def test_it_agrees_with_the_contract_path_draw_by_draw(self) -> None:
        problem = self.problem()
        lowered = lower_problem(problem)
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
        problem = self.problem()
        lowered = lower_problem(problem)
        theta = self.theta(problem, 9)
        whole = lowered.simulate_batched(theta)
        chunked = lowered.simulate_batched(theta, chunk_size=chunk_size)
        for label, values in whole.predicted.items():
            assert np.array_equal(values, chunked.predicted[label])

    def test_a_quasiseparable_problem_refuses_by_name(self) -> None:
        lowered = lower_problem(_gp_problem(QuasisepGP()))
        with pytest.raises(LoweringError, match="QuasisepGP"):
            lowered.simulate_batched(np.zeros((3, lowered.free_size)))

    def test_a_non_stack_is_refused(self) -> None:
        lowered = lower_problem(self.problem())
        with pytest.raises(LoweringError, match="batch, "):
            lowered.simulate_batched(np.zeros(lowered.free_size))

    def test_the_single_device_sharder_is_the_degenerate_case(self) -> None:
        """CPU CI's half of Peter's addendum: the hook runs, and changes nothing."""
        from ampere.backends.jax import MeshSharder, SingleDeviceSharder

        problem = self.problem()
        lowered = lower_problem(problem)
        theta = self.theta(problem, 6)
        plain = lowered.simulate_batched(theta)
        for sharder in (SingleDeviceSharder(), SingleDeviceSharder("cpu"), MeshSharder()):
            assert len(sharder.devices()) >= 1
            produced = lowered.simulate_batched(theta, sharder=sharder)
            for label, values in plain.predicted.items():
                assert np.allclose(values, produced.predicted[label], rtol=0.0, atol=1e-12)

    def test_a_sharder_names_a_platform_it_cannot_find(self) -> None:
        from ampere.backends.jax import SingleDeviceSharder

        with pytest.raises(LoweringError, match="no device on platform"):
            SingleDeviceSharder("no-such-platform")

    def test_a_native_draw_is_a_pure_function_of_its_seed(self) -> None:
        """What makes partition independence possible on a backend with its own RNG."""
        problem = self.problem()
        lowered = lower_problem(problem)
        theta = self.theta(problem, 4)
        prediction = lowered.simulate_batched(theta)
        first = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 44])
        second = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 44])
        shuffled = lowered.sample_observations(theta, prediction.predicted, [11, 22, 33, 45])
        for label, values in first.items():
            assert np.array_equal(values, second[label])
            assert np.array_equal(values[:3], shuffled[label][:3])
            assert not np.array_equal(values[3], shuffled[label][3])

    def test_a_family_the_core_will_not_sample_is_refused_by_name(self) -> None:
        """Peter's ruling has a ceiling: a backend samples what ``ampere.core`` samples."""
        from ampere.core import PoissonFamily

        problem = _family_problem(PoissonFamily(), counts=True)
        lowered = lower_problem(problem)
        theta = self.theta(problem, 2)
        prediction = lowered.simulate_batched(theta)
        with pytest.raises(LoweringError, match="numpy path"):
            lowered.sample_observations(theta, prediction.predicted, [1, 2])
