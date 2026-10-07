"""``ampere.results.Optimum`` — the point-estimate record (W6.7).

The dataclass on its own: its invariants, its identity hash, ``combine``,
``to_datatree`` (an ``optimum`` group and no ``posterior``), the netCDF round
trip and the ``summary`` table. The routes that build one are
``tests/inference/test_optimise.py``'s.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from ampere.core.exceptions import ResultsError
from ampere.results import OPTIMUM_GROUP, Optimum, StartSummary, from_netcdf, to_netcdf
from ampere.results.optimum import aligned


def make(**overrides: object) -> Optimum:
    fields: dict[str, object] = {
        "route": "scipy",
        "backend": "reference",
        "free_names": ("model.a", "model.b"),
        "free_labels": ("model.a", "model.b"),
        "unconstrained": np.array([0.5, -1.0]),
        "constrained": {"model.a": 0.5, "model.b": math.exp(-1.0)},
        "log_prob_constrained": -3.25,
        "log_prob_unconstrained": -4.25,
        "covariance": np.array([[0.04, 0.01], [0.01, 0.09]]),
        "covariance_refusal": None,
        "converged": True,
        "message": "Optimization terminated successfully.",
        "evaluations": 123,
        "starts": (
            StartSummary("a" * 32, -3.25, "converged"),
            StartSummary("b" * 32, -math.inf, "failed: every proposal scored -inf"),
        ),
        "provenance": {"ampere_schema_version": 9, "ampere_backend": "reference"},
    }
    fields.update(overrides)
    return Optimum(**fields)  # type: ignore[arg-type]


class TestTheRecord:
    def test_vectors_are_read_only(self) -> None:
        optimum = make()
        with pytest.raises(ValueError):
            optimum.unconstrained[0] = 1.0
        assert optimum.covariance is not None
        with pytest.raises(ValueError):
            optimum.covariance[0, 0] = 1.0

    def test_one_label_per_entry(self) -> None:
        with pytest.raises(ResultsError, match="one label per entry"):
            make(free_labels=("model.a",))

    def test_a_missing_covariance_must_say_why(self) -> None:
        with pytest.raises(ResultsError, match="covariance_refusal"):
            make(covariance=None)
        refused = make(covariance=None, covariance_refusal="Hessian not positive definite")
        assert refused.standard_deviations is None

    def test_not_both(self) -> None:
        with pytest.raises(ResultsError, match="not both"):
            make(covariance_refusal="why")

    def test_identity_moves_with_the_point_and_the_route_only(self) -> None:
        base = make()
        assert base.identity == make(message="other", evaluations=1).identity
        assert base.identity != make(route="map").identity
        assert base.identity != make(unconstrained=np.array([0.5, -1.0 + 1e-12])).identity
        assert len(base.identity) == 32

    def test_start_record(self) -> None:
        record = make().start_record()
        assert record == {
            "route": "scipy",
            "identity": make().identity,
            "log_prob_constrained": -3.25,
            "evaluations": 123,
            "converged": True,
        }

    def test_summary_names_every_label(self) -> None:
        text = make().summary()
        assert "model.a" in text and "model.b" in text and "converged" in text
        refused = make(covariance=None, covariance_refusal="saddle").summary()
        assert "covariance refused: saddle" in refused


class TestCombine:
    def test_disjoint_names_merge_block_diagonally(self) -> None:
        other = make(
            route="empirical_bayes",
            free_names=("irs.gp.scale",),
            free_labels=("irs.gp.scale",),
            unconstrained=np.array([2.0]),
            constrained={"irs.gp.scale": 2.0},
            covariance=np.array([[0.25]]),
            starts=(),
        )
        both = Optimum.combine(make(), other)
        assert both.route == "scipy+empirical_bayes"
        assert both.free_labels == ("model.a", "model.b", "irs.gp.scale")
        assert both.covariance is not None
        np.testing.assert_allclose(both.covariance[:2, :2], make().covariance)
        assert both.covariance[2, 2] == 0.25 and both.covariance[0, 2] == 0.0
        assert math.isnan(both.log_prob_constrained)
        assert both.evaluations == 246

    def test_overlap_refused_by_name(self) -> None:
        with pytest.raises(ResultsError, match=r"'model\.a' is covered by both"):
            Optimum.combine(make(), make(route="map"))

    def test_a_refused_part_refuses_the_whole(self) -> None:
        other = make(
            free_names=("c",),
            free_labels=("c",),
            unconstrained=np.array([0.0]),
            constrained={"c": 0.0},
            covariance=None,
            covariance_refusal="flat",
        )
        both = Optimum.combine(make(), other)
        assert both.covariance is None
        assert both.covariance_refusal is not None and "no covariance" in both.covariance_refusal


class TestAligned:
    def test_reorders_vector_and_covariance(self) -> None:
        vector, covariance = aligned(make(), ("model.b", "model.a"))
        np.testing.assert_array_equal(vector, [-1.0, 0.5])
        assert covariance is not None
        np.testing.assert_array_equal(covariance, [[0.09, 0.01], [0.01, 0.04]])

    def test_mismatch_refused_by_name(self) -> None:
        with pytest.raises(ResultsError, match=r"missing \['model.c'\]"):
            aligned(make(), ("model.a", "model.b", "model.c"))


class TestTheDataTree:
    def test_an_optimum_group_and_no_posterior(self) -> None:
        tree = make().to_datatree()
        assert set(tree.children) == {OPTIMUM_GROUP}
        group = tree[OPTIMUM_GROUP].dataset
        assert list(group["free_parameter"].values) == ["model.a", "model.b"]
        np.testing.assert_array_equal(group["unconstrained"].values, [0.5, -1.0])
        assert group["covariance"].shape == (2, 2)
        assert tree.attrs["ampere_optimum_route"] == "scipy"
        assert tree.attrs["ampere_schema_version"] == 9

    @pytest.mark.parametrize("refused", [False, True])
    def test_netcdf_round_trip(self, tmp_path: object, refused: bool) -> None:
        extra: dict[str, object] = {"covariance": np.eye(3) * 0.1}
        if refused:
            extra = {"covariance": None, "covariance_refusal": "not positive definite"}
        original = make(
            free_names=("model.a", "model.v"),
            free_labels=("model.a", "model.v[0]", "model.v[1]"),
            unconstrained=np.array([0.5, 1.0, 2.0]),
            constrained={"model.a": 0.5, "model.v": np.array([1.0, 2.0])},
            **extra,
        )
        path = to_netcdf(original.to_datatree(), f"{tmp_path}/optimum.nc")
        back = Optimum.from_datatree(from_netcdf(path))
        assert back.identity == original.identity
        assert back.free_names == original.free_names
        assert back.free_labels == original.free_labels
        np.testing.assert_array_equal(back.constrained["model.v"], [1.0, 2.0])
        assert back.constrained["model.a"] == 0.5
        assert back.starts == original.starts
        assert back.log_prob_constrained == original.log_prob_constrained
        assert back.covariance_refusal == original.covariance_refusal
        assert back.provenance["ampere_schema_version"] == 9
        assert back.evaluations == 123 and back.converged

    def test_a_run_is_not_an_optimum(self) -> None:
        import xarray as xr

        with pytest.raises(ResultsError, match="not an Optimum"):
            Optimum.from_datatree(xr.DataTree())


class TestTheCoordinates:
    """``Optimum.coordinates`` (W7.6): what the minimiser moved each entry in."""

    def test_omitted_means_every_entry_moved_in_u(self) -> None:
        assert make().coordinates == ("unconstrained", "unconstrained")

    def test_one_known_kind_per_label(self) -> None:
        assert make(coordinates=("constrained", "unconstrained")).coordinates == (
            "constrained",
            "unconstrained",
        )
        with pytest.raises(ResultsError, match="coordinates"):
            make(coordinates=("constrained",))
        with pytest.raises(ResultsError, match="coordinates"):
            make(coordinates=("constrained", "sideways"))

    def test_identity_does_not_see_them(self) -> None:
        assert make(coordinates=("constrained", "constrained")).identity == make().identity

    def test_combine_carries_them(self) -> None:
        other = make(
            route="empirical_bayes",
            free_names=("irs.gp.scale",),
            free_labels=("irs.gp.scale",),
            unconstrained=np.array([2.0]),
            constrained={"irs.gp.scale": 2.0},
            covariance=np.array([[0.25]]),
            starts=(),
        )
        both = Optimum.combine(make(coordinates=("constrained", "unconstrained")), other)
        assert both.coordinates == ("constrained", "unconstrained", "unconstrained")

    def test_netcdf_round_trip(self, tmp_path: object) -> None:
        original = make(coordinates=("constrained", "unconstrained"))
        tree = original.to_datatree()
        assert tree.attrs["ampere_optimum_coordinates"] == '["constrained","unconstrained"]'
        back = Optimum.from_datatree(from_netcdf(to_netcdf(tree, f"{tmp_path}/optimum.nc")))
        assert back.coordinates == original.coordinates

    def test_a_beta_tree_without_them_reads_as_unconstrained(self) -> None:
        tree = make(coordinates=("constrained", "constrained")).to_datatree()
        del tree.attrs["ampere_optimum_coordinates"]
        assert Optimum.from_datatree(tree).coordinates == ("unconstrained", "unconstrained")
