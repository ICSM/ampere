"""W4.5: the kernel algebra, the axis selector and the public term registry.

Organised around the item's acceptance criteria. Three groups carry the real
weight and everything else supports them:

* ``TestClosedForms`` — each new term against its textbook closed form,
  evaluated independently here rather than by calling the kernel twice;
* ``TestSemiseparableRepresentations`` — each registered representation
  rebuilt into a dense matrix from ``(c, U, V)`` and compared against the
  kernel's own ``matrix``, which is what makes "exactly quasiseparable" a
  claim about the mathematics rather than about a log-likelihood that hides it;
* ``TestAxisSelector`` — the selection, its unit rule, and the ``Product`` of
  two axis-selected kernels against the elementwise product of their matrices
  on a three-axis container.

The cross-solver and cross-backend rows live in ``tests/conformance``; this
module is the single-backend, single-implementation half.
"""

from __future__ import annotations

import math

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    SHO,
    AxisSpec,
    DenseGP,
    FunctionSamples,
    GaussianFamily,
    GaussianProcessNoise,
    Layout,
    Likelihood,
    Matern12,
    Matern32,
    Matern52,
    Order,
    Parameter,
    Product,
    QuasisepGP,
    RotationTerm,
    SpectralMixture,
    Spectrum,
    SquaredExponential,
    StationaryKernel,
    Sum,
    lookup_quasiseparable_term,
    quasiseparable_families,
    register_quasiseparable_term,
    registered_quasiseparable_terms,
    term_provenance_entries,
)
from ampere.core.exceptions import LikelihoodError
from ampere.core.kernels import (
    _forget_quasiseparable_term,
    matern12_representation,
)

SQRT3 = math.sqrt(3.0)
SQRT5 = math.sqrt(5.0)


# ---------------------------------------------------------------------------
# A three-axis point-set kind, local to these tests
# ---------------------------------------------------------------------------


class DispersedPoints(FunctionSamples):
    """Three axes, two of them dimensionless and one in micron.

    The smallest container that makes the axis selector necessary rather than
    merely available: an isotropic kernel over all three is refused by the
    single-unit rule, ``Matern32(axes=("u", "v"))`` is fine, and their
    ``Product`` is the chromatic covariance of ``phase4_placement_memo.md``
    §3.6.

    Declared here rather than taken from ``VisibilitySet`` deliberately: W4.1
    is amending that kind to three axes in parallel, and a test that waited for
    it would be a dependency the item does not have.
    """

    AXES = (
        AxisSpec("u", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("v", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
        AxisSpec("spectral_axis", physical_types=("length",), order=Order.ANY),
    )
    LAYOUT = Layout.POINTS

    def __init__(
        self,
        u_coord: object,
        v_coord: object,
        spectral_axis: object,
        values: object,
        **kwargs: object,
    ) -> None:
        super().__init__(
            {"u": u_coord, "v": v_coord, "spectral_axis": spectral_axis},
            values,
            **kwargs,  # type: ignore[arg-type]
        )


@pytest.fixture
def dispersed() -> DispersedPoints:
    rng = np.random.default_rng(20260911)
    n = 9
    return DispersedPoints(
        rng.uniform(-3.0, 3.0, n) * u.dimensionless_unscaled,
        rng.uniform(-3.0, 3.0, n) * u.dimensionless_unscaled,
        np.sort(rng.uniform(2.0, 2.4, n)) * u.um,
        rng.normal(0.0, 1.0, n) * u.Jy,
        uncertainty=np.full(n, 0.1) * u.Jy,
    )


@pytest.fixture
def coordinates() -> np.ndarray:
    rng = np.random.default_rng(20260905)
    return np.sort(rng.uniform(0.0, 12.0, 15))


# ---------------------------------------------------------------------------
# Closed forms
# ---------------------------------------------------------------------------


def _sho_closed_form(
    tau: np.ndarray, amplitude: float, period: float, quality: float
) -> np.ndarray:
    """celerite2's underdamped ``SHOTerm``, written out independently."""
    omega0 = 2.0 * math.pi / period
    split = math.sqrt(4.0 * quality**2 - 1.0)
    decay = 0.5 * omega0 / quality
    frequency = decay * split
    separation = np.abs(tau)
    return (
        amplitude**2
        * np.exp(-decay * separation)
        * (np.cos(frequency * separation) + np.sin(frequency * separation) / split)
    )


class TestClosedForms:
    """Each new term against its textbook form, and ``k(0) == amplitude**2``."""

    def test_matern12(self, coordinates: np.ndarray) -> None:
        got = Matern12(0.4, 2.0).matrix(coordinates, coordinates, {})
        separation = np.abs(coordinates[:, None] - coordinates[None, :])
        np.testing.assert_allclose(got, 0.16 * np.exp(-separation / 2.0), atol=1e-14)

    def test_matern32(self, coordinates: np.ndarray) -> None:
        got = Matern32(0.4, 2.0).matrix(coordinates, coordinates, {})
        scaled = SQRT3 * np.abs(coordinates[:, None] - coordinates[None, :]) / 2.0
        np.testing.assert_allclose(got, 0.16 * (1.0 + scaled) * np.exp(-scaled), atol=1e-14)

    def test_matern52(self, coordinates: np.ndarray) -> None:
        got = Matern52(0.4, 2.0).matrix(coordinates, coordinates, {})
        scaled = SQRT5 * np.abs(coordinates[:, None] - coordinates[None, :]) / 2.0
        expected = 0.16 * (1.0 + scaled + scaled**2 / 3.0) * np.exp(-scaled)
        np.testing.assert_allclose(got, expected, atol=1e-14)

    def test_sho(self, coordinates: np.ndarray) -> None:
        got = SHO(0.4, 1.5, 3.0).matrix(coordinates, coordinates, {})
        tau = coordinates[:, None] - coordinates[None, :]
        np.testing.assert_allclose(got, _sho_closed_form(tau, 0.4, 1.5, 3.0), atol=1e-13)

    def test_rotation_is_two_shos_at_a_period_and_its_harmonic(
        self, coordinates: np.ndarray
    ) -> None:
        """celerite2's ``RotationTerm`` parameterisation, rebuilt from its own algebra."""
        amplitude, period, q0, dq, fraction = 0.4, 1.5, 2.0, 0.5, 0.3
        got = RotationTerm(amplitude, period, q0, dq, fraction).matrix(coordinates, coordinates, {})
        tau = coordinates[:, None] - coordinates[None, :]
        power = amplitude**2 / (1.0 + fraction)
        expected = np.zeros_like(tau)
        for factor, harmonic, share in ((0.5 + q0 + dq, 1.0, 1.0), (0.5 + q0, 2.0, fraction)):
            split = math.sqrt(4.0 * factor**2 - 1.0)
            omega = harmonic * 4.0 * math.pi * factor / (period * split)
            decay = 0.5 * omega / factor
            frequency = decay * split
            separation = np.abs(tau)
            expected = expected + share * power * np.exp(-decay * separation) * (
                np.cos(frequency * separation) + np.sin(frequency * separation) / split
            )
        np.testing.assert_allclose(got, expected, atol=1e-13)

    @pytest.mark.parametrize(
        "kernel",
        [
            Matern12(0.7, 2.0),
            Matern32(0.7, 2.0),
            Matern52(0.7, 2.0),
            SquaredExponential(0.7, 2.0),
            SHO(0.7, 1.5, 3.0),
            RotationTerm(0.7, 1.5, 2.0, 0.5, 0.3),
        ],
        ids=lambda kernel: str(kernel.FAMILY),
    )
    def test_the_amplitude_is_a_marginal_standard_deviation(self, kernel: object) -> None:
        """``k(0) == amplitude**2`` for every family: one convention, no exceptions."""
        diagonal = kernel.diagonal(np.zeros(3), {})  # type: ignore[attr-defined]
        np.testing.assert_allclose(np.asarray(diagonal), 0.49, atol=1e-14)

    def test_the_rotation_terms_marginal_variance_is_the_amplitude_squared(self) -> None:
        """The power split must not change ``k(0)``; ``fraction`` only moves it about."""
        for fraction in (0.0, 0.3, 1.0, 4.0):
            kernel = RotationTerm(0.7, 1.5, 2.0, 0.5, fraction)
            assert float(kernel.diagonal(np.zeros(1), {})[0]) == pytest.approx(0.49, abs=1e-14)

    def test_an_overdamped_sho_is_refused_by_name(self) -> None:
        with pytest.raises(LikelihoodError, match=r"'quality' must be > 1/2"):
            SHO(0.4, 1.5, 0.5).matrix([[0.0]], [[1.0]], {})

    def test_the_spectral_mixture_is_the_sum_of_its_components(
        self, coordinates: np.ndarray
    ) -> None:
        mixture = SpectralMixture([0.3, 0.2], [1.5, 0.4], [3.0, 8.0])
        tau = coordinates[:, None] - coordinates[None, :]
        expected = _sho_closed_form(tau, 0.3, 1.5, 3.0) + _sho_closed_form(tau, 0.2, 0.4, 8.0)
        got = mixture.matrix(coordinates, coordinates, {})
        np.testing.assert_allclose(got, expected, atol=1e-13)


# ---------------------------------------------------------------------------
# The semiseparable representations
# ---------------------------------------------------------------------------


def _rebuild(kernel: object, coordinates: np.ndarray) -> np.ndarray:
    """``K`` from ``(c, U, V)``, in celerite2's own convention.

    ``K[n, m] = sum_j U[n, j] V[m, j] exp(-c[j] (t_n - t_m))`` for ``n > m``,
    symmetric, with ``k(0)`` on the diagonal. Written out here rather than
    taken from celerite2 so that a wrong representation fails against the
    kernel's own closed form instead of against another celerite recursion.
    """
    builder = lookup_quasiseparable_term(kernel.FAMILY)  # type: ignore[attr-defined]
    representation = builder(kernel, kernel.resolve(None), coordinates)  # type: ignore[attr-defined]
    decay = np.asarray(representation.decay)
    left = np.asarray(representation.left)
    right = np.asarray(representation.right)
    n = coordinates.size
    rebuilt = np.zeros((n, n))
    for i in range(n):
        rebuilt[i, i] = float(np.asarray(representation.marginal))
        for j in range(i):
            entry = float(
                np.sum(left[i] * right[j] * np.exp(-decay * (coordinates[i] - coordinates[j])))
            )
            rebuilt[i, j] = rebuilt[j, i] = entry
    return rebuilt


class TestSemiseparableRepresentations:
    """Every registered representation, rebuilt densely and compared with ``matrix``.

    This is the row that makes "exactly quasiseparable" mean an algebraic
    identity. A log-likelihood comparison would pass on a representation that
    was merely close; this one compares the matrices themselves, at 1e-13,
    which nothing approximate reaches.
    """

    @pytest.mark.parametrize(
        "kernel",
        [
            Matern12(0.4, 2.0),
            Matern32(0.4, 2.0),
            Matern52(0.4, 2.0),
            SHO(0.4, 1.5, 3.0),
            RotationTerm(0.4, 1.5, 2.0, 0.5, 0.3),
            Sum(Matern32(0.3, 2.0), Matern12(0.2, 0.4)),
            Sum(Matern52(0.3, 3.0), SHO(0.2, 1.5, 4.0), Matern12(0.1, 0.2)),
            SpectralMixture([0.3, 0.2], [1.5, 0.4], [3.0, 8.0]),
        ],
        ids=[
            "matern12",
            "matern32",
            "matern52",
            "sho",
            "rotation",
            "sum-of-two",
            "sum-of-three",
            "spectral-mixture",
        ],
    )
    def test_the_generators_rebuild_the_kernel_matrix(
        self, kernel: object, coordinates: np.ndarray
    ) -> None:
        expected = kernel.matrix(coordinates, coordinates, {})  # type: ignore[attr-defined]
        np.testing.assert_allclose(_rebuild(kernel, coordinates), expected, atol=1e-13)

    def test_the_ranks_are_the_documented_ones(self, coordinates: np.ndarray) -> None:
        """Rank is the cost of the O(N) solve; a silent inflation would be a regression."""
        ranks = {}
        for kernel in (
            Matern12(0.4, 2.0),
            Matern32(0.4, 2.0),
            Matern52(0.4, 2.0),
            SHO(0.4, 1.5, 3.0),
            RotationTerm(0.4, 1.5, 2.0, 0.5, 0.3),
        ):
            builder = lookup_quasiseparable_term(kernel.FAMILY)
            ranks[kernel.FAMILY] = builder(kernel, kernel.resolve(None), coordinates).rank
        assert ranks == {"matern12": 1, "matern32": 2, "matern52": 3, "sho": 2, "rotation": 4}

    def test_a_sums_rank_is_the_sum_of_its_terms(self, coordinates: np.ndarray) -> None:
        """The structural claim that makes the algebra worth having."""
        kernel = Sum(Matern52(0.3, 3.0), SHO(0.2, 1.5, 4.0), Matern12(0.1, 0.2))
        builder = lookup_quasiseparable_term(kernel.FAMILY)
        assert builder(kernel, kernel.resolve(None), coordinates).rank == 3 + 2 + 1

    @pytest.mark.parametrize(
        "kernel",
        [
            Matern12(0.4, 2.0),
            Matern32(0.4, 2.0),
            Matern52(0.4, 2.0),
            SHO(0.4, 1.5, 3.0),
            RotationTerm(0.4, 1.5, 2.0, 0.5, 0.3),
            Sum(Matern32(0.3, 2.0), SHO(0.2, 1.5, 4.0)),
        ],
        ids=["matern12", "matern32", "matern52", "sho", "rotation", "sum"],
    )
    def test_the_quasiseparable_solver_agrees_with_the_dense_one(self, kernel: object) -> None:
        """``DenseGP`` is the definition of the right answer; the O(N) path must match it."""
        rng = np.random.default_rng(20260911)
        axis = np.sort(rng.uniform(0.0, 12.0, 40))
        data = Spectrum(
            axis * u.um,
            rng.normal(0.0, 0.3, axis.size) * u.Jy,
            uncertainty=np.full(axis.size, 0.1) * u.Jy,
        )
        model = data.with_values(np.zeros(axis.size))
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
        quasisep = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        assert quasisep.log_prob(model, data) == pytest.approx(
            dense.log_prob(model, data), abs=1e-6
        )


# ---------------------------------------------------------------------------
# The algebra
# ---------------------------------------------------------------------------


class TestAlgebra:
    def test_a_sum_is_the_sum_of_the_matrices(self, coordinates: np.ndarray) -> None:
        first, second = Matern32(0.3, 2.0), Matern12(0.2, 0.4)
        np.testing.assert_allclose(
            Sum(first, second).matrix(coordinates, coordinates, {}),
            first.matrix(coordinates, coordinates, {})
            + second.matrix(coordinates, coordinates, {}),
            atol=1e-15,
        )

    def test_a_product_is_the_elementwise_product_of_the_matrices(
        self, coordinates: np.ndarray
    ) -> None:
        first, second = Matern32(1.0, 2.0), SHO(0.4, 1.5, 3.0)
        np.testing.assert_allclose(
            Product(first, second).matrix(coordinates, coordinates, {}),
            first.matrix(coordinates, coordinates, {})
            * second.matrix(coordinates, coordinates, {}),
            atol=1e-15,
        )

    def test_children_are_namespaced_so_two_of_a_family_do_not_collide(self) -> None:
        kernel = Sum(
            Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
            Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0)),
        )
        assert kernel.parameters.free_names == (
            "term0.amplitude",
            "term0.length_scale",
            "term1.amplitude",
            "term1.length_scale",
        )

    def test_labels_may_be_given_and_must_be_identifiers(self) -> None:
        named = Sum(Matern32(0.3, 2.0), Matern32(0.1, 0.2), labels=("broad", "narrow"))
        assert named.parameters.names[:2] == ("broad.amplitude", "broad.length_scale")
        with pytest.raises(LikelihoodError, match="not a Python identifier"):
            Sum(Matern32(0.3, 2.0), Matern32(0.1, 0.2), labels=("terms.0", "terms.1"))
        with pytest.raises(LikelihoodError, match="repeat a name"):
            Sum(Matern32(0.3, 2.0), Matern32(0.1, 0.2), labels=("a", "a"))

    def test_a_composites_hyperparameters_route_to_the_right_child(
        self, coordinates: np.ndarray
    ) -> None:
        """The namespacing must be a *routing*, not only a naming."""
        kernel = Sum(
            Matern32(Parameter("amplitude", value=0.3, fixed=True), 2.0),
            Matern12(Parameter("amplitude", value=0.2, fixed=True), 0.4),
        )
        values = {"term0.amplitude": 0.5, "term1.amplitude": 0.1}
        expected = Matern32(0.5, 2.0).matrix(coordinates, coordinates, {}) + Matern12(
            0.1, 0.4
        ).matrix(coordinates, coordinates, {})
        np.testing.assert_allclose(
            kernel.matrix(coordinates, coordinates, values), expected, atol=1e-15
        )

    def test_a_sum_is_quasiseparable_exactly_when_its_terms_are(self) -> None:
        assert Sum(Matern32(0.3, 2.0), Matern12(0.1, 0.2)).QUASISEPARABLE is True
        assert Sum(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)).QUASISEPARABLE is False

    def test_a_product_is_never_quasiseparable(self) -> None:
        assert Product(Matern32(0.3, 2.0), Matern12(0.1, 0.2)).QUASISEPARABLE is False

    def test_the_spec_is_a_tree_in_declaration_order(self) -> None:
        spec = Sum(Matern32(0.3, 2.0), Matern12(0.1, 0.2)).spec()
        assert spec.family == "sum"
        assert [label for label, _ in spec.terms] == ["term0", "term1"]
        assert [child.family for _, child in spec.terms] == ["matern32", "matern12"]
        assert spec.to_dict()["terms"] == [
            {
                "label": "term0",
                "kernel": {
                    "family": "matern32",
                    "hyperparameters": ["amplitude", "length_scale"],
                    "quasiseparable": True,
                },
            },
            {
                "label": "term1",
                "kernel": {
                    "family": "matern12",
                    "hyperparameters": ["amplitude", "length_scale"],
                    "quasiseparable": True,
                },
            },
        ]

    def test_reordering_a_sum_changes_its_spec(self) -> None:
        """Commutative in mathematics, not in hashing — and the reason is the names."""
        first, second = Matern32(0.3, 2.0), Matern12(0.1, 0.2)
        assert Sum(first, second).spec().to_dict() != Sum(second, first).spec().to_dict()

    def test_a_one_term_composite_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="composes two or more kernels"):
            Sum(Matern32(0.3, 2.0))

    def test_a_non_kernel_term_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="composes Kernels"):
            Sum(Matern32(0.3, 2.0), "a kernel, honest")  # type: ignore[arg-type]

    def test_a_product_is_refused_on_the_quasiseparable_path_word_for_word(self) -> None:
        data = Spectrum(
            np.linspace(1.0, 10.0, 8) * u.um,
            np.zeros(8) * u.Jy,
            uncertainty=np.full(8, 0.1) * u.Jy,
        )
        noise = GaussianProcessNoise(Product(Matern32(1.0, 2.0), Matern12(0.4, 0.5)), QuasisepGP())
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), data)
        assert str(excinfo.value) == (
            "QuasisepGP cannot lower a Product: a product of quasiseparable kernels is not "
            "quasiseparable. Where the factors act on different axes — which is what a Product "
            "is for — the result is not a function of one ordered coordinate at all, and where "
            "they act on the same one the semiseparable rank multiplies and is not recoverable "
            "from the factors' own representations. Use DenseGP, or replace the Product with a "
            "Sum, which is quasiseparable exactly when every term is."
        )

    def test_a_sum_containing_a_product_names_the_product_not_the_sum(self) -> None:
        """W5.2: the sharper refusal follows the Product wherever it nests.

        Before W5.2 this fell through to the generic ``not QUASISEPARABLE``
        refusal (`QuasisepGP needs a kernel with an exact quasiseparable
        representation, but Sum ...`), which names the *container* rather
        than the term that is actually the problem.
        """
        data = Spectrum(
            np.linspace(1.0, 10.0, 8) * u.um,
            np.zeros(8) * u.Jy,
            uncertainty=np.full(8, 0.1) * u.Jy,
        )
        product = Product(Matern32(1.0, 2.0), Matern12(0.4, 0.5), labels=("spatial", "spectral"))
        noise = GaussianProcessNoise(
            Sum(Matern32(0.3, 2.0), product, labels=("smooth", "coupled")), QuasisepGP()
        )
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), data)
        assert str(excinfo.value) == (
            "QuasisepGP cannot lower this Sum: its term 'coupled' is a Product, and a product "
            "of quasiseparable kernels is not quasiseparable. Where the factors act on "
            "different axes — which is what a Product is for — the result is not a function of "
            "one ordered coordinate at all, and where they act on the same one the "
            "semiseparable rank multiplies and is not recoverable from the factors' own "
            "representations. Use DenseGP, or replace 'coupled' with a Sum, which is "
            "quasiseparable exactly when every term is."
        )

    def test_a_sum_with_one_unlowerable_term_is_refused_by_name(self) -> None:
        data = Spectrum(
            np.linspace(1.0, 10.0, 8) * u.um,
            np.zeros(8) * u.Jy,
            uncertainty=np.full(8, 0.1) * u.Jy,
        )
        noise = GaussianProcessNoise(
            Sum(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)), QuasisepGP()
        )
        with pytest.raises(LikelihoodError, match="exact quasiseparable representation"):
            noise.check_compatible(GaussianFamily(), data)


# ---------------------------------------------------------------------------
# The axis selector
# ---------------------------------------------------------------------------


class TestAxisSelector:
    def test_the_default_is_every_axis_and_the_spec_is_unchanged(self) -> None:
        """The compatibility promise: no existing declaration or hash moves."""
        spec = Matern32(0.3, 2.0).spec()
        assert spec.axes is None
        assert spec.to_dict() == {
            "family": "matern32",
            "hyperparameters": ["amplitude", "length_scale"],
            "quasiseparable": True,
        }

    def test_a_selection_reaches_the_spec_and_so_the_hash(self) -> None:
        spec = Matern32(0.3, 2.0, axes=("u", "v")).spec()
        assert spec.axes == ("u", "v")
        assert spec.to_dict()["axes"] == ["u", "v"]
        assert spec.to_dict() != Matern32(0.3, 2.0, axes=("spectral_axis",)).spec().to_dict()

    def test_a_bare_string_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="not the single string"):
            Matern32(0.3, 2.0, axes="u")  # type: ignore[arg-type]

    def test_an_empty_or_repeating_selection_is_refused(self) -> None:
        with pytest.raises(LikelihoodError, match="empty axes"):
            Matern32(0.3, 2.0, axes=())
        with pytest.raises(LikelihoodError, match="repeats an axis name"):
            Matern32(0.3, 2.0, axes=("u", "u"))

    def test_a_selection_measures_separation_in_the_selected_axes_only(
        self, dispersed: DispersedPoints
    ) -> None:
        points = np.column_stack([np.asarray(axis.values) for axis in dispersed.axes])
        kernel = Matern32(0.3, 2.0, axes=("u", "v")).for_axes(["u", "v", "spectral_axis"])
        plane = points[:, :2]
        separation = np.sqrt(
            ((plane[:, None, :] - plane[None, :, :]) ** 2).sum(-1),
        )
        scaled = SQRT3 * separation / 2.0
        np.testing.assert_allclose(
            kernel.matrix(points, points, {}),
            0.09 * (1.0 + scaled) * np.exp(-scaled),
            atol=1e-14,
        )

    def test_binding_is_functional_so_one_kernel_serves_two_containers(self) -> None:
        kernel = Matern32(0.3, 2.0, axes=("spectral_axis",))
        first = kernel.for_axes(["u", "v", "spectral_axis"])
        second = kernel.for_axes(["spectral_axis", "u"])
        assert first is not second
        assert kernel.axes == ("spectral_axis",)
        assert first is kernel.for_axes(["u", "v", "spectral_axis"])

    def test_an_unknown_axis_name_is_refused(self, dispersed: DispersedPoints) -> None:
        noise = GaussianProcessNoise(Matern32(0.3, 2.0, axes=("w",)), DenseGP())
        with pytest.raises(LikelihoodError, match=r"it has no 'w'"):
            noise.check_compatible(GaussianFamily(), dispersed)

    def test_a_mixed_unit_container_is_still_refused_for_a_bare_kernel(
        self, dispersed: DispersedPoints
    ) -> None:
        """The pre-W4.5 rule, word for word; **the advice amended at W4.2**.

        The rule is unchanged — a Euclidean separation across mixed units is
        meaningless and is refused — and the remedy the message names is not. It
        said "use one axis, or declare a kernel that takes a length-scale per
        axis", which was the fix before the selector existed; on the modality
        this refusal actually fires for (a ``VisibilitySet``, whose axes are
        ``(u, v, spectral_axis)``) "use one axis" is not the fix and
        ``axes=("u", "v")`` is, so the message says that instead.
        """
        noise = GaussianProcessNoise(Matern32(0.3, 2.0), DenseGP())
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), dispersed)
        assert str(excinfo.value) == (
            "DenseGP measures separation as a Euclidean distance across a DispersedPoints's "
            "coordinate axes, but they carry different units ['', 'um']. A single isotropic "
            "length-scale is meaningless across mixed units; name the axes this kernel acts on "
            "with axes=(...), as in Matern32(axes=('u',)), and compose kernels on different axes "
            "with Product."
        )

    def test_a_mixed_unit_selection_is_refused_word_for_word(
        self, dispersed: DispersedPoints
    ) -> None:
        """The new rule: the single-unit check follows the selection."""
        noise = GaussianProcessNoise(Matern32(0.3, 2.0, axes=("u", "spectral_axis")), DenseGP())
        with pytest.raises(LikelihoodError) as excinfo:
            noise.check_compatible(GaussianFamily(), dispersed)
        assert str(excinfo.value) == (
            "Matern32 selects the DispersedPoints's axes ('u', 'spectral_axis'), which carry "
            "different units ['', 'um']. A single isotropic length-scale is meaningless across "
            "mixed units; select a subset whose axes share one unit, and compose the rest with "
            "Product."
        )

    def test_a_product_of_two_selections_is_accepted_and_correct(
        self, dispersed: DispersedPoints
    ) -> None:
        """The acceptance row: the chromatic covariance of the memo's §3.6.

        Two kernels on disjoint axis subsets of a three-axis container, the
        product compared against the elementwise product of the factors'
        matrices — which is the definition, evaluated through the composition
        machinery rather than beside it.
        """
        spatial = Matern32(1.0, 2.0, axes=("u", "v"))
        spectral = Matern32(0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um)
        kernel = Product(spatial, spectral)
        noise = GaussianProcessNoise(kernel, DenseGP())
        noise.check_compatible(GaussianFamily(), dispersed)

        points = np.column_stack([np.asarray(axis.values) for axis in dispersed.axes])
        bound = noise.kernel_for(dispersed)
        got = bound.matrix(points, points, {})
        expected = spatial.for_axes(["u", "v", "spectral_axis"]).matrix(
            points, points, {}
        ) * spectral.for_axes(["u", "v", "spectral_axis"]).matrix(points, points, {})
        np.testing.assert_allclose(got, expected, atol=1e-15)

    def test_the_product_reaches_the_log_likelihood(self, dispersed: DispersedPoints) -> None:
        """Not merely composable: a number comes out, and it is the dense Gaussian's."""
        kernel = Product(
            Matern32(1.0, 2.0, axes=("u", "v")),
            Matern32(0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um),
        )
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
        model = dispersed.with_values(np.zeros(dispersed.n_samples))
        points = np.column_stack([np.asarray(axis.values) for axis in dispersed.axes])
        bound = kernel.for_axes(["u", "v", "spectral_axis"])
        covariance = bound.matrix(points, points, {}) + np.diag(
            np.asarray(dispersed.uncertainty) ** 2
        )
        expected = st.multivariate_normal(cov=covariance).logpdf(np.asarray(dispersed.values))
        assert likelihood.log_prob(model, dispersed) == pytest.approx(float(expected), abs=1e-8)

    def test_a_selected_length_scale_is_checked_against_its_own_axis_unit(
        self, dispersed: DispersedPoints
    ) -> None:
        """Per-leaf units: a dimensionless (u, v) scale beside a micron spectral one."""
        good = Product(
            Matern32(1.0, 2.0, axes=("u", "v")),
            Matern32(0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um),
        )
        GaussianProcessNoise(good, DenseGP()).check_compatible(GaussianFamily(), dispersed)
        bad = Product(
            Matern32(1.0, 2.0, axes=("u", "v"), length_scale_unit=u.um),
            Matern32(0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um),
        )
        with pytest.raises(LikelihoodError, match=r"term0\.length_scale"):
            GaussianProcessNoise(bad, DenseGP()).check_compatible(GaussianFamily(), dispersed)

    def test_a_selected_axis_reaches_the_quasiseparable_path(
        self, dispersed: DispersedPoints
    ) -> None:
        """One ordered axis out of three: the O(N) solve on the spectral column."""
        kernel = Matern32(0.3, 0.05, axes=("spectral_axis",), length_scale_unit=u.um)
        model = dispersed.with_values(np.zeros(dispersed.n_samples))
        dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
        quasisep = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, QuasisepGP()))
        assert quasisep.log_prob(model, dispersed) == pytest.approx(
            dense.log_prob(model, dispersed), abs=1e-6
        )

    def test_an_unselected_multi_axis_container_is_refused_on_the_quasiseparable_path(
        self, dispersed: DispersedPoints
    ) -> None:
        noise = GaussianProcessNoise(Matern32(0.3, 2.0), QuasisepGP())
        with pytest.raises(LikelihoodError, match="different units"):
            noise.check_compatible(GaussianFamily(), dispersed)


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------


class Relabelled(Matern12):
    """A user kernel: Matérn-1/2 under another name, so its representation is known.

    Deliberately trivial. What the row proves is the *route* — that a kernel
    declared outside ampere reaches the O(N) path once its representation is
    registered — and a term whose mathematics had to be checked too would make
    a failure ambiguous.
    """

    FAMILY = "relabelled_matern12"


class TestRegistry:
    def test_the_built_in_families_are_registered(self) -> None:
        assert quasiseparable_families() == (
            "matern12",
            "matern32",
            "matern52",
            "rotation",
            "sho",
            "spectral_mixture",
            "sum",
        )

    def test_built_in_rows_say_so_and_are_not_stamped_in_provenance(self) -> None:
        rows = {row.family: row for row in registered_quasiseparable_terms()}
        assert rows["matern32"].builtin is True
        assert rows["matern32"].builder_module == "ampere.core.kernels"
        assert term_provenance_entries() == []

    def test_an_unregistered_family_is_refused_by_name(self) -> None:
        with pytest.raises(LikelihoodError, match="Families with one: matern12, matern32"):
            lookup_quasiseparable_term("no_such_family")

    def test_a_second_registration_is_refused_unless_deliberate(self) -> None:
        try:
            register_quasiseparable_term(Relabelled, matern12_representation)
            with pytest.raises(LikelihoodError, match="already registered"):
                register_quasiseparable_term(Relabelled, matern12_representation)
            row = register_quasiseparable_term(Relabelled, matern12_representation, override=True)
            assert row.builtin is False
        finally:
            _forget_quasiseparable_term(Relabelled.FAMILY)

    def test_a_built_in_row_also_requires_override(self) -> None:
        with pytest.raises(LikelihoodError, match="a built-in quasiseparable representation"):
            register_quasiseparable_term(Matern32, matern12_representation)

    def test_a_non_quasiseparable_kernel_cannot_register_one(self) -> None:
        with pytest.raises(LikelihoodError, match="declares QUASISEPARABLE = False"):
            register_quasiseparable_term(SquaredExponential, matern12_representation)

    def test_a_user_row_is_stamped_in_provenance(self) -> None:
        try:
            register_quasiseparable_term(Relabelled, matern12_representation)
            entries = term_provenance_entries()
            assert entries == [
                {
                    "kind": "quasiseparable_term",
                    "family": "relabelled_matern12",
                    "kernel": "Relabelled",
                    "builtin": False,
                    "builder": "ampere.core.kernels.matern12_representation",
                }
            ]
        finally:
            _forget_quasiseparable_term(Relabelled.FAMILY)

    def test_a_registered_user_kernel_reaches_the_quasiseparable_solver(self) -> None:
        """The acceptance row, on the reference backend.

        The cross-backend half — the same kernel, the same one registration,
        reaching ``QuasisepGP`` on torch and jax too — is a conformance row.
        """
        rng = np.random.default_rng(20260911)
        axis = np.sort(rng.uniform(0.0, 12.0, 30))
        data = Spectrum(
            axis * u.um,
            rng.normal(0.0, 0.3, axis.size) * u.Jy,
            uncertainty=np.full(axis.size, 0.1) * u.Jy,
        )
        model = data.with_values(np.zeros(axis.size))
        kernel = Relabelled(0.4, 2.0)
        noise = GaussianProcessNoise(kernel, QuasisepGP())
        with pytest.raises(LikelihoodError, match="relabelled_matern12"):
            noise.check_compatible(GaussianFamily(), data)
        try:
            register_quasiseparable_term(Relabelled, matern12_representation)
            dense = Likelihood(GaussianFamily(), GaussianProcessNoise(kernel, DenseGP()))
            quasisep = Likelihood(GaussianFamily(), noise)
            assert quasisep.log_prob(model, data) == pytest.approx(
                dense.log_prob(model, data), abs=1e-6
            )
        finally:
            _forget_quasiseparable_term(Relabelled.FAMILY)


class TestStationaryKernelIsAPublicBase:
    def test_a_user_kernel_gets_the_two_hyperparameters_for_free(self) -> None:
        class Triangular(StationaryKernel):
            """``k(r) = a**2 max(0, 1 - r/L)`` — a kernel nobody should use, declared."""

            FAMILY = "triangular"
            HYPERPARAMETERS = ("amplitude", "length_scale")

            def _covariance(self, separation: object, values: object) -> np.ndarray:
                amplitude = self._hyperparameter(values, "amplitude", allow_zero=True)  # type: ignore[arg-type]
                length_scale = self._hyperparameter(values, "length_scale")  # type: ignore[arg-type]
                scaled = np.asarray(separation) / length_scale
                return amplitude * amplitude * np.clip(1.0 - scaled, 0.0, None)

        kernel = Triangular(st.halfnorm(0.0, 1.0), 2.0)
        assert kernel.parameters.free_names == ("amplitude",)
        assert kernel.spec().quasiseparable is False
        got = kernel.matrix([[0.0], [1.0]], [[0.0], [1.0]], {"amplitude": 2.0})
        np.testing.assert_allclose(got, np.array([[4.0, 2.0], [2.0, 4.0]]), atol=1e-15)
