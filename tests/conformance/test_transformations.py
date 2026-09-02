"""Conformance rows for the transformation and instrument contract.

The inventory is ``transformations.md`` §14's hand-down, verbatim: "chain
composition rejects a kind mismatch; ``require`` is used for binding; a dropped
mask raises; ``propagate_mask`` matches the ANY rule on a random influence
matrix; a union of requirements satisfies each of its inputs; ``compile_for``'s
default is the identity; a compiled model's templates share axes across
evaluations."

These are composition-time rows: ``negotiate`` and ``compile_for`` run once,
before any tracing, so nothing here depends on a backend's array type. What it
*does* depend on is the backend's transformations and models declaring
themselves truthfully — which is exactly what a lowering bug breaks.
"""

from __future__ import annotations

from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest

from ampere.core import (
    AxisRequirement,
    ChannelError,
    CompositionError,
    Instrument,
    ModelResult,
    PhotometricPoints,
    Spectrum,
    TransformationError,
    negotiate,
    propagate_mask,
)
from ampere.core.transform import Transformation

from .composition import COARSE_GRID, COORDINATE_UNIT, FLUX_UNIT, GP_GRID
from .protocol import (
    ConformanceBackend,
    ModelSpec,
    TransformationKind,
    TransformationSpec,
)

SCALE = TransformationSpec(TransformationKind.SCALE, label="calibration")
REBIN = TransformationSpec(TransformationKind.REBIN, label="binner", target=COARSE_GRID)
PHOTOMETRY = TransformationSpec(
    TransformationKind.PHOTOMETRY,
    label="synphot",
    target=COARSE_GRID,
    filters=("W1", "W2", "W3", "W4"),
)


def masked_spectrum(masked: tuple[int, ...] = (2, 5)) -> Spectrum:
    """A ``Spectrum`` on :data:`GP_GRID`, with *masked* samples excluded.

    An empty *masked* gives a container with no mask at all, not an all-``False``
    one: ``Transformation.__call__``'s dropped-mask guard tests ``is None``,
    so the two are different declarations and the rows need both.
    """
    grid = np.asarray(GP_GRID, dtype=float)
    mask = None
    if masked:
        mask = np.zeros(grid.size, dtype=bool)
        mask[list(masked)] = True
    return Spectrum(
        grid * COORDINATE_UNIT,
        np.linspace(1.0, 2.0, grid.size) * FLUX_UNIT,
        uncertainty=np.full(grid.size, 0.1) * FLUX_UNIT,
        mask=mask,
    )


def at_prior_median(model: Any) -> Any:
    """Evaluate *model* at its prior median — a value every declaration has."""
    centre = np.full(model.parameters.free_size, 0.5)
    return model(model.parameters.prior_transform(centre))


class MaskDropper(Transformation):
    """A deliberately broken step: it returns a container with no mask.

    Not something a backend supplies — it exists so the row can prove the
    guard in ``Transformation.__call__`` fires, rather than only proving that
    the backend's own steps happen not to trip it.
    """

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def apply(self, samples: Any, values: Any) -> Spectrum:
        return Spectrum(
            samples.spectral_axis.quantity(), samples.values, unit=samples.unit, mask=None
        )


class TestChainComposition:
    """Kinds must compose, and a bound channel must hold what was asked for."""

    def test_a_chain_whose_kinds_do_not_meet_is_refused_at_construction(
        self, backend: ConformanceBackend
    ) -> None:
        steps = (backend.transformation(PHOTOMETRY), backend.transformation(SCALE))
        with pytest.raises(CompositionError, match="cannot be composed"):
            Instrument(steps, channel="sed")

    def test_the_same_steps_compose_the_other_way_round(self, backend: ConformanceBackend) -> None:
        instrument = Instrument(
            (backend.transformation(SCALE), backend.transformation(PHOTOMETRY)),
            channel="sed",
        )
        assert instrument.input_kind is Spectrum
        assert instrument.output_kind is PhotometricPoints

    def test_a_step_refuses_a_container_of_the_wrong_kind(
        self, backend: ConformanceBackend
    ) -> None:
        photometry = backend.transformation(PHOTOMETRY)
        points = PhotometricPoints(("A", "B"), [1.0, 2.0] * COORDINATE_UNIT, [1.0, 2.0] * FLUX_UNIT)
        with pytest.raises(CompositionError, match="accepts"):
            photometry(points, None)


class TestChannelBinding:
    """``Instrument.bind`` goes through ``ModelResult.require``."""

    def test_binding_selects_the_declared_channel(self, backend: ConformanceBackend) -> None:
        instrument = Instrument((backend.transformation(SCALE),), channel="blue", label="blue_arm")
        result = ModelResult({"blue": masked_spectrum(()), "red": masked_spectrum(())})
        assert instrument.bind(result).spectral_axis.values == pytest.approx(GP_GRID)

    def test_a_missing_channel_names_what_was_available(self, backend: ConformanceBackend) -> None:
        instrument = Instrument(
            (backend.transformation(SCALE),), channel="green", label="green_arm"
        )
        result = ModelResult({"blue": masked_spectrum(())})
        with pytest.raises(ChannelError, match="no channel named 'green'"):
            instrument.bind(result)

    def test_a_channel_of_the_wrong_kind_is_a_channel_error(
        self, backend: ConformanceBackend
    ) -> None:
        instrument = Instrument((backend.transformation(SCALE),), channel="sed", label="sed_arm")
        points = PhotometricPoints(("A", "B"), [1.0, 2.0] * COORDINATE_UNIT, [1.0, 2.0] * FLUX_UNIT)
        with pytest.raises(ChannelError):
            instrument.bind(ModelResult({"sed": points}))


class TestMaskPropagation:
    """The ANY rule, and the refusal to lose a mask on the way through."""

    def test_propagate_mask_is_the_any_rule_on_a_random_influence_matrix(self) -> None:
        """An output sample is masked iff *any* input that influences it is.

        The oracle is the rule itself, written out with an explicit loop over
        a random sparse matrix — deliberately not the vectorised form, so a
        transcription error in either would show up.
        """
        rng = np.random.default_rng(20260902)
        n_in, n_out = 17, 9
        influence = rng.normal(size=(n_out, n_in)) * (rng.random((n_out, n_in)) < 0.4)
        mask = rng.random(n_in) < 0.3
        source = Spectrum(
            np.arange(1.0, n_in + 1.0) * COORDINATE_UNIT,
            rng.normal(size=n_in) * FLUX_UNIT,
            mask=mask,
        )

        expected = np.array(
            [any(mask[j] for j in range(n_in) if influence[i, j] != 0.0) for i in range(n_out)]
        )
        assert propagate_mask(source, influence).tolist() == expected.tolist()

    def test_an_unmasked_source_propagates_nothing(self) -> None:
        assert propagate_mask(masked_spectrum(()), np.eye(len(GP_GRID))) is None
        assert propagate_mask(None) is None

    def test_a_backend_resampling_step_propagates_its_input_mask(
        self, backend: ConformanceBackend
    ) -> None:
        rebin = backend.transformation(REBIN)
        source = masked_spectrum()
        rebinned = rebin(source, None)
        expected = propagate_mask(source, rebin.influence(source.spectral_axis.values))
        assert rebinned.mask is not None
        assert rebinned.mask.tolist() == expected.tolist()
        assert rebinned.is_masked

    def test_a_step_that_drops_the_mask_is_refused(self) -> None:
        with pytest.raises(TransformationError, match="dropped the mask"):
            MaskDropper()(masked_spectrum(), None)

    def test_the_guard_does_not_fire_on_an_unmasked_input(self) -> None:
        assert MaskDropper()(masked_spectrum(()), None).mask is None


class TestRequirementNegotiation:
    """A union satisfies each of its inputs; every bound channel is published."""

    def test_a_union_covers_both_inputs(self) -> None:
        blue = AxisRequirement(
            "spectral_axis", unit=u.um, intervals=(1.0, 4.0), max_step=0.5, source="blue"
        )
        red = AxisRequirement(
            "spectral_axis", unit=u.um, intervals=(6.0, 10.0), max_step=1.0, source="red"
        )
        merged = blue.union(red)

        coordinates = merged.coordinates().to_value(u.um)
        for lower, upper in ((1.0, 4.0), (6.0, 10.0)):
            inside = coordinates[(coordinates >= lower) & (coordinates <= upper)]
            assert inside.size >= 2
            assert inside.min() == pytest.approx(lower)
            assert inside.max() == pytest.approx(upper)
        assert "blue" in merged.source and "red" in merged.source

    def test_a_union_keeps_each_intervals_own_density(self) -> None:
        """The union of a broad-coarse and a narrow-fine request is not fine everywhere."""
        coarse = AxisRequirement("spectral_axis", unit=u.um, intervals=(1.0, 10.0), max_step=3.0)
        fine = AxisRequirement("spectral_axis", unit=u.um, intervals=(4.0, 5.0), max_step=0.25)
        coordinates = coarse.union(fine).coordinates().to_value(u.um)
        window = coordinates[(coordinates >= 4.0) & (coordinates <= 5.0)]
        assert np.diff(window).max() <= 0.25 + 1e-12
        assert coordinates.size < 9.0 / 0.25

    def test_negotiate_publishes_every_bound_channel(self, backend: ConformanceBackend) -> None:
        plain = Instrument((backend.transformation(SCALE),), channel="blue", label="blue_arm")
        binned = Instrument((backend.transformation(REBIN),), channel="red", label="red_arm")
        requirements = negotiate([plain, binned])

        assert set(requirements) == {"blue", "red"}
        assert requirements["blue"].axes == {}
        asked = requirements["red"]["spectral_axis"].coordinates().to_value(u.um)
        assert asked.min() == pytest.approx(COARSE_GRID[0])
        assert asked.max() == pytest.approx(COARSE_GRID[-1])

    def test_two_instruments_on_one_channel_union_their_requirements(
        self, backend: ConformanceBackend
    ) -> None:
        first = Instrument((backend.transformation(REBIN),), channel="sed", label="one")
        second = Instrument(
            (
                backend.transformation(
                    TransformationSpec(
                        TransformationKind.REBIN, label="binner", target=(2.0, 3.0, 4.0)
                    )
                ),
            ),
            channel="sed",
            label="two",
        )
        asked = negotiate([first, second])["sed"]["spectral_axis"]
        assert "one" in asked.source and "two" in asked.source
        coordinates = asked.coordinates().to_value(u.um)
        assert coordinates.min() == pytest.approx(COARSE_GRID[0])
        assert coordinates.max() == pytest.approx(COARSE_GRID[-1])


class TestCompileFor:
    """The default is the identity; a compiled model reuses its templates."""

    def test_a_model_that_ignores_requirements_returns_itself(
        self, backend: ConformanceBackend
    ) -> None:
        model = backend.model(ModelSpec(coordinates=GP_GRID, compiled=False))
        instrument = Instrument((backend.transformation(REBIN),), channel="default")
        assert model.compile_for({}) is model
        assert model.compile_for(negotiate([instrument])) is model

    def test_an_uncompiled_model_evaluates_on_its_own_grid(
        self, backend: ConformanceBackend
    ) -> None:
        model = backend.model(ModelSpec(coordinates=GP_GRID, compiled=False))
        instrument = Instrument((backend.transformation(REBIN),), channel="default")
        compiled = model.compile_for(negotiate([instrument]))
        assert at_prior_median(compiled).single().spectral_axis.values == pytest.approx(GP_GRID)

    def test_a_compiled_model_honours_the_requested_grid(self, backend: ConformanceBackend) -> None:
        model = backend.model(ModelSpec(coordinates=GP_GRID, compiled=True))
        instrument = Instrument((backend.transformation(REBIN),), channel="default")
        requirements = negotiate([instrument])
        compiled = model.compile_for(requirements)

        asked = requirements["default"]["spectral_axis"].coordinates().to_value(u.um)
        assert at_prior_median(compiled).single().spectral_axis.values == pytest.approx(asked)

    def test_a_compiled_models_templates_share_axes_across_evaluations(
        self, backend: ConformanceBackend
    ) -> None:
        """The hot-loop contract: axes are built once and reused by identity.

        ``transformations.md`` §14 asks for this because rebuilding the
        coordinate arrays every evaluation is the cost ``compile_for`` exists
        to remove — and an implementation that quietly rebuilds them still
        passes a value comparison.
        """
        model = backend.model(ModelSpec(coordinates=GP_GRID, compiled=True))
        instrument = Instrument((backend.transformation(REBIN),), channel="default")
        compiled = model.compile_for(negotiate([instrument]))

        first = at_prior_median(compiled).single()
        second = at_prior_median(compiled).single()
        assert first.axes is second.axes
        assert first.unit == second.unit
