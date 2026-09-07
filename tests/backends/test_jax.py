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
    Resample,
    configure_x64,
    filter_spec,
    lower_bijection,
    lower_problem,
)
from ampere.backends.jax.distributions import has_native_icdf, lower_prior  # noqa: E402
from ampere.backends.jax.rng import fold, key  # noqa: E402
from ampere.core import (  # noqa: E402
    Dataset,
    DatasetCollection,
    FittingProblem,
    GaussianFamily,
    HierarchicalPrior,
    Identity,
    Instrument,
    Likelihood,
    Log,
    Logit,
    Parameter,
    ParameterSet,
    Plate,
    Spectrum,
    Tie,
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
        ["ampere.backends.jax", "ampere.backends.jax.problem"],
        ids=["package", "problem"],
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
