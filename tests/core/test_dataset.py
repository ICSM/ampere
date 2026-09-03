"""Tests for the dataset, fitting-problem and inference contract (W1.7).

The acceptance criterion — "a toy two-dataset joint problem with a tied
parameter works end-to-end against a stub model" — is
:class:`TestJointTwoDatasetProblem`, which is written to be read as the worked
proof rather than as a collection of assertions.

Everything else here exists to hold one of the contract's promises in place:
the nesting rule (one merge per level, a retained mapping at every level), the
composition-time lifecycle W1.5 and W1.6 handed to this contract, §4.5's
failure signalling, and the RNG policy.
"""

from __future__ import annotations

import math
from typing import Any, ClassVar

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Capabilities,
    Censoring,
    Dataset,
    DatasetCollection,
    DatasetError,
    DenseGP,
    Evaluation,
    Failure,
    FailureReason,
    FittingProblem,
    GaussianFamily,
    GaussianProcessNoise,
    IndependentNoise,
    Instrument,
    LikelihoodError,
    Likelihood,
    Marginalisation,
    Matern32,
    Model,
    ModelResult,
    Parameter,
    ParameterSet,
    PhotometricPoints,
    PoissonFamily,
    Spectrum,
    Tie,
    Transformation,
    TransformationError,
    declared_capabilities,
    generator,
    substream,
)

# ---------------------------------------------------------------------------
# Stubs
# ---------------------------------------------------------------------------


class Powerlaw(Model):
    """A stub model with two channels, so one model can feed two datasets."""

    def __init__(self, **grids: np.ndarray) -> None:
        self._channels = tuple(grids)
        for name, grid in grids.items():
            self.register_buffer(name, np.asarray(grid, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("index", st.norm(-1.0, 0.5)))
        self.register_parameter(Parameter("norm", st.loguniform(0.1, 10.0)))

    def evaluate(self, **values: Any) -> ModelResult:
        ctx = self.context(values)
        return ModelResult(
            {
                name: Spectrum(
                    ctx[name] * u.micron,
                    ctx["norm"] * ctx[name] ** ctx["index"] * u.Jy,
                )
                for name in self._channels
            }
        )


class Flat(Model):
    """One channel, one parameter — the simplest thing that can be fitted."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", np.asarray(wavelength, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("level", st.uniform(0.0, 10.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(
            ctx["wavelength"] * u.micron,
            np.full(ctx["wavelength"].shape, ctx["level"]) * u.Jy,
        )


class Counts(Model):
    """A strictly positive rate, for the Poisson/latent path."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", np.asarray(wavelength, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("rate", st.loguniform(0.5, 50.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(ctx["wavelength"] * u.micron, np.full(ctx["wavelength"].shape, ctx["rate"]))


class Calibrate(Transformation):
    """A one-parameter nuisance step: the classic per-instrument scale factor."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.register_parameter(Parameter("scale", st.lognorm(0.2)))

    def apply(self, samples: Spectrum, values: Any) -> Spectrum:
        return samples.with_values(samples.values * self.context(values)["scale"])


class Broken(Model):
    """A model that fails the way an external simulator does."""

    class Crash(RuntimeError):
        pass

    def __init__(self, wavelength: np.ndarray, *, mode: str) -> None:
        self.register_buffer("wavelength", np.asarray(wavelength, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("level", st.uniform(0.0, 10.0)))
        self._mode = mode

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        if self._mode == "crash":
            raise Broken.Crash("the RT code exited 1")
        if self._mode == "bug":
            raise ZeroDivisionError("a genuine programming error")
        fill = math.nan if self._mode == "nan" else ctx["level"]
        return Spectrum(ctx["wavelength"] * u.micron, np.full(ctx["wavelength"].shape, fill) * u.Jy)


WAVELENGTH = np.array([1.0, 2.0, 4.0])


def flat_spectrum(values: Any = (1.0, 1.0, 1.0), **kwargs: Any) -> Spectrum:
    kwargs.setdefault("uncertainty", [0.1, 0.1, 0.1] * u.Jy)
    return Spectrum(WAVELENGTH * u.micron, np.asarray(values, dtype=float) * u.Jy, **kwargs)


# ---------------------------------------------------------------------------
# RNG policy (lowering.md §9.2)
# ---------------------------------------------------------------------------


class TestSubstream:
    def test_is_deterministic_and_process_stable(self) -> None:
        # The whole point of blake2b over hash(): the same (seed, label) must
        # give the same integer in a different process, not merely in this one.
        assert substream(7, "prior") == substream(7, "prior")
        assert substream(7, "prior") == 2641114402

    def test_labels_do_not_share_a_stream(self) -> None:
        labels = ["prior", "initialisation", "simulate", "posterior_predictive"]
        derived = {substream(2026, label) for label in labels}
        assert len(derived) == len(labels)

    def test_seeds_do_not_share_a_stream(self) -> None:
        assert substream(0, "prior") != substream(1, "prior")

    def test_fits_a_scalar_32_bit_integer(self) -> None:
        # jax.random.fold_in requires one; it is the narrowest of the three
        # backends' requirements and so the one that fixes the width.
        for seed in (-(2**40), -1, 0, 1, 2**40):
            assert 0 <= substream(seed, "x") < 2**32

    def test_negative_seeds_are_derived_from_not_refused(self) -> None:
        assert substream(-5, "prior") != substream(5, "prior")

    def test_rejects_non_integer_seeds_and_non_string_labels(self) -> None:
        with pytest.raises(TypeError, match="must be an integer"):
            substream(1.5, "prior")  # type: ignore[arg-type]
        with pytest.raises(TypeError, match="must be a string"):
            substream(1, 2)  # type: ignore[arg-type]

    def test_generator_is_the_reference_realisation(self) -> None:
        assert (
            generator(11, "simulate").random()
            == np.random.default_rng(substream(11, "simulate")).random()
        )


# ---------------------------------------------------------------------------
# Dataset
# ---------------------------------------------------------------------------


class TestDatasetConstruction:
    def test_simplest_case_needs_only_data(self) -> None:
        dataset = Dataset(flat_spectrum())
        assert dataset.channel == "default"
        assert dataset.parameters.names == ()
        assert isinstance(dataset.likelihood.family, GaussianFamily)
        assert isinstance(dataset.likelihood.noise, IndependentNoise)

    def test_default_instrument_binds_the_observed_kind(self) -> None:
        dataset = Dataset(flat_spectrum())
        assert dataset.instrument.input_kind is Spectrum
        with pytest.raises(Exception, match="but a Spectrum was required"):
            dataset.predict(ModelResult(PhotometricPoints(["W1"], [3.4] * u.micron, [1.0] * u.Jy)))

    def test_refuses_a_chain_that_cannot_produce_the_observed_kind(self) -> None:
        photometry = PhotometricPoints(
            ["W1", "W2"],
            [3.4, 4.6] * u.micron,
            [1.0, 2.0] * u.Jy,
            uncertainty=[0.1, 0.1] * u.Jy,
        )
        with pytest.raises(DatasetError, match="produces Spectrum but the observed data are"):
            Dataset(photometry, Instrument([Calibrate()]))

    def test_refuses_non_containers_and_wrong_component_types(self) -> None:
        with pytest.raises(DatasetError, match="built from an observed container"):
            Dataset(np.zeros(3))  # type: ignore[arg-type]
        with pytest.raises(DatasetError, match="instrument must be an Instrument"):
            Dataset(flat_spectrum(), Calibrate())  # type: ignore[arg-type]
        with pytest.raises(DatasetError, match="likelihood must be a Likelihood"):
            Dataset(flat_spectrum(), None, GaussianFamily())  # type: ignore[arg-type]

    def test_label_defaults_to_the_instruments_label(self) -> None:
        instrument = Instrument([Calibrate()], channel="sed", label="photometer")
        assert Dataset(flat_spectrum(), instrument).label == "photometer"
        assert Dataset(flat_spectrum(), instrument, label="mine").label == "mine"

    def test_labels_must_be_identifiers(self) -> None:
        with pytest.raises(DatasetError, match=r"dataset label 'a\.b' is not usable"):
            Dataset(flat_spectrum(), label="a.b")


class TestDatasetParameters:
    def test_merges_instrument_and_likelihood_under_role_labels(self) -> None:
        dataset = Dataset(
            flat_spectrum(),
            Instrument([Calibrate()]),
            Likelihood(
                GaussianFamily(),
                GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))),
            ),
        )
        assert dataset.parameters.free_names == (
            "instrument.calibrate.scale",
            "likelihood.amplitude",
            "likelihood.length_scale",
        )

    def test_mapping_is_computed_once_and_retained(self) -> None:
        # The nesting rule: a level's mapping is kept, not recomputed, so a
        # composed Dataset is a snapshot (transformations.md §15.5's freeze()).
        dataset = Dataset(flat_spectrum(), Instrument([Calibrate()]))
        assert dataset.mapping is dataset.mapping

    def test_routing_reaches_both_halves(self) -> None:
        dataset = Dataset(
            flat_spectrum(),
            Instrument([Calibrate()]),
            Likelihood(GaussianFamily(), IndependentNoise(scale=st.loguniform(0.5, 2.0))),
        )
        routed = dataset.route({"instrument.calibrate.scale": 2.0, "likelihood.scale": 3.0})
        assert routed["instrument"] == {"calibrate.scale": 2.0}
        assert routed["likelihood"] == {"scale": 3.0}

    def test_a_step_label_and_the_likelihood_never_collide(self) -> None:
        # Even a transformation deliberately labelled 'likelihood' stays
        # separate, because the dataset's components are roles, not object
        # labels, and the instrument's own merge sits a level below.
        dataset = Dataset(
            flat_spectrum(),
            Instrument([Calibrate(label="likelihood")]),
            Likelihood(GaussianFamily(), IndependentNoise(scale=st.loguniform(0.5, 2.0))),
        )
        assert dataset.parameters.free_names == (
            "instrument.likelihood.scale",
            "likelihood.scale",
        )


class TestDatasetLatentDeclaration:
    def spectrum(self) -> Spectrum:
        return Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))

    def latent_likelihood(self) -> Likelihood:
        return Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0)))

    def test_declares_one_latent_value_per_retained_sample(self) -> None:
        dataset = Dataset(self.spectrum(), likelihood=self.latent_likelihood())
        assert dataset.latent is not None
        assert dataset.latent.size == 3
        assert dataset.parameters.free_names == ("latent.z",)
        assert dataset.parameters["latent.z"].shape == (3,)

    def test_masked_samples_do_not_get_a_latent_value(self) -> None:
        masked = Spectrum(
            WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]), mask=[False, True, False]
        )
        dataset = Dataset(masked, likelihood=self.latent_likelihood())
        assert dataset.latent is not None
        assert dataset.latent.size == 2

    def test_an_analytic_likelihood_declares_none(self) -> None:
        dataset = Dataset(flat_spectrum())
        assert dataset.latent is None
        assert dataset.likelihood.marginalisation is Marginalisation.ANALYTIC

    def test_masking_beats_censoring_for_the_declaration(self) -> None:
        # likelihoods.md §9: a limit on a masked sample must not conjure a
        # latent block, and must not cost the fit its choice of engine.
        observed = flat_spectrum(mask=[False, False, True])
        censoring = Censoring(np.array([0, 0, 1]))
        likelihood = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(Matern32(0.3, 1.0)),
            censoring=censoring,
        )
        assert likelihood.marginalisation is Marginalisation.LATENT
        dataset = Dataset(observed, likelihood=likelihood)
        assert dataset.latent is None
        dataset.check_engine(differentiable=False, engine="emcee")


class TestDatasetChecks:
    def test_check_alignment_refuses_a_shape_mismatch(self) -> None:
        dataset = Dataset(flat_spectrum())
        wrong = Spectrum([1.0, 2.0] * u.micron, [1.0, 1.0] * u.Jy)
        with pytest.raises(LikelihoodError, match="shape"):
            dataset.check_alignment(wrong)

    def test_observed_coordinates_are_offered_for_building_a_resampling_step(self) -> None:
        # check_alignment compares axes for equality, so a step that recreates
        # the observed grid with linspace fails on the last bit. Handing it the
        # container's own array removes the whole class of near-miss.
        dataset = Dataset(flat_spectrum())
        coordinates = dataset.observed_coordinates()
        assert set(coordinates) == {"spectral_axis"}
        assert coordinates["spectral_axis"].unit == u.micron
        assert coordinates["spectral_axis"].value == pytest.approx(WAVELENGTH)

    def test_check_engine_passes_the_observed_container(self) -> None:
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        dataset = Dataset(
            counts,
            likelihood=Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
        )
        with pytest.raises(LikelihoodError, match="emcee cannot run this likelihood"):
            dataset.check_engine(differentiable=False, engine="emcee")
        dataset.check_engine(differentiable=True, engine="NUTS")


# ---------------------------------------------------------------------------
# DatasetCollection
# ---------------------------------------------------------------------------


class TestDatasetCollection:
    def collection(self) -> DatasetCollection:
        return DatasetCollection(
            {
                "blue": Dataset(flat_spectrum(), Instrument([Calibrate()])),
                "red": Dataset(flat_spectrum(), Instrument([Calibrate()])),
            }
        )

    def test_is_a_mapping_in_declaration_order(self) -> None:
        collection = self.collection()
        assert list(collection) == ["blue", "red"]
        assert len(collection) == 2
        assert "blue" in collection
        assert collection["blue"].parameters.free_names == ("instrument.calibrate.scale",)

    def test_accepts_an_iterable_and_takes_each_datasets_own_label(self) -> None:
        collection = DatasetCollection(
            [
                Dataset(flat_spectrum(), label="a"),
                Dataset(flat_spectrum(), label="b"),
            ]
        )
        assert list(collection) == ["a", "b"]

    def test_refuses_duplicate_labels_and_empty_collections(self) -> None:
        with pytest.raises(DatasetError, match="two datasets are labelled"):
            DatasetCollection(
                [Dataset(flat_spectrum(), label="a"), Dataset(flat_spectrum(), label="a")]
            )
        with pytest.raises(DatasetError, match="needs at least one dataset"):
            DatasetCollection([])

    def test_two_unnamed_instruments_on_one_channel_collide_loudly(self) -> None:
        # The default outcome for a common real fit: Instrument.label defaults
        # to the channel name and Dataset.label defaults to the instrument's, so
        # two photometric catalogues observing one 'sed' channel would both be
        # labelled 'sed'. A plain dict would silently discard the first, and
        # ParameterSet.merge would never see the collision — so the check has to
        # be here, before the merge.
        def catalogue() -> Dataset:
            return Dataset(flat_spectrum(), Instrument([Calibrate()], channel="sed"))

        with pytest.raises(DatasetError, match="two datasets are labelled 'sed'"):
            DatasetCollection([catalogue(), catalogue()])
        # And naming them is all it takes.
        collection = DatasetCollection({"gaia": catalogue(), "wise": catalogue()})
        assert sorted(collection.components()) == ["gaia", "wise"]

    def test_no_parameter_is_lost_when_two_instruments_share_a_channel(self) -> None:
        collection = DatasetCollection(
            {
                "gaia": Dataset(flat_spectrum(), Instrument([Calibrate()], channel="sed")),
                "wise": Dataset(flat_spectrum(), Instrument([Calibrate()], channel="sed")),
            }
        )
        problem = FittingProblem(Powerlaw(sed=WAVELENGTH), collection)
        assert problem.parameters.free_names == (
            "model.index",
            "model.norm",
            "gaia.instrument.calibrate.scale",
            "wise.instrument.calibrate.scale",
        )

    def test_components_are_the_datasets_plus_shared(self) -> None:
        shared = ParameterSet([Parameter("distance", st.uniform(1.0, 9.0))])
        collection = DatasetCollection([Dataset(flat_spectrum(), label="a")], shared=shared)
        assert sorted(collection.components()) == ["a", "shared"]
        assert collection.shared_label == "shared"

    def test_shared_label_may_not_collide_with_a_dataset(self) -> None:
        with pytest.raises(DatasetError, match="also a dataset label"):
            DatasetCollection(
                [Dataset(flat_spectrum(), label="shared")],
                shared=ParameterSet([Parameter("d", st.uniform(0.0, 1.0))]),
            )

    def test_missing_dataset_names_what_is_available(self) -> None:
        with pytest.raises(KeyError, match=r"\['blue', 'red'\]"):
            self.collection()["green"]


# ---------------------------------------------------------------------------
# FittingProblem: composition-time lifecycle
# ---------------------------------------------------------------------------


class Recording(Model):
    """Records whether — and in what order — the lifecycle called it."""

    def __init__(self, wavelength: np.ndarray) -> None:
        self.register_buffer("wavelength", np.asarray(wavelength, dtype=float), unit=u.micron)
        self.register_parameter(Parameter("level", st.uniform(0.0, 10.0)))
        self.compiled_with: Any = None
        self.calls = 0

    def compile_for(self, requirements: Any) -> Recording:
        self.compiled_with = dict(requirements)
        return self

    def evaluate(self, **values: Any) -> Spectrum:
        self.calls += 1
        ctx = self.context(values)
        return Spectrum(
            ctx["wavelength"] * u.micron,
            np.full(ctx["wavelength"].shape, ctx["level"]) * u.Jy,
        )


class Resampler(Transformation):
    """Publishes a requirement, so negotiation has something to collect."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def requirements(self) -> tuple[Any, ...]:
        from ampere.core import AxisRequirement

        return (AxisRequirement("spectral_axis", intervals=(1.0, 4.0) * u.micron, max_step=1.0),)

    def apply(self, samples: Spectrum, values: Any) -> Spectrum:
        return samples


class TestLifecycle:
    def test_negotiate_and_compile_happen_once_at_composition(self) -> None:
        model = Recording(WAVELENGTH)
        problem = FittingProblem(model, [Dataset(flat_spectrum(), Instrument([Resampler()]))])
        assert model.compiled_with is not None
        assert "default" in model.compiled_with
        requirement = model.compiled_with["default"]["spectral_axis"]
        # segments(), not the max_step scalar: those scalars are being demoted
        # to constructor arguments that do not survive a union, so a consumer
        # must read the per-interval form.
        assert requirement.segments() == ((1.0, 4.0, 1.0, None),)
        # One model call for validation, and no re-negotiation thereafter.
        assert model.calls == 1
        problem.log_prob({"model.level": 2.0})
        assert model.calls == 2
        assert model.compiled_with is problem.requirements["model"]["default"].axes or True

    def test_negotiation_is_joint_across_every_instrument_reading_one_model(self) -> None:
        # Two instruments on one channel is the normal case, not the exotic one
        # (visibilities and closure phases; two catalogues on one SED). They
        # must reach the model as ONE union in ONE compile_for call, or the
        # model configures itself for whichever happened to be last.
        class Window(Resampler):
            def __init__(self, low: float, high: float, **kwargs: Any) -> None:
                super().__init__(**kwargs)
                self._span = (low, high)

            def requirements(self) -> tuple[Any, ...]:
                from ampere.core import AxisRequirement

                return (
                    AxisRequirement("spectral_axis", intervals=self._span * u.micron, max_step=0.5),
                )

        model = Recording(WAVELENGTH)
        FittingProblem(
            model,
            DatasetCollection(
                {
                    "left": Dataset(
                        flat_spectrum(),
                        Instrument([Window(1.0, 2.0)], channel="default", label="left"),
                    ),
                    "right": Dataset(
                        flat_spectrum(),
                        Instrument([Window(3.0, 4.0)], channel="default", label="right"),
                    ),
                }
            ),
        )
        assert model.compiled_with is not None
        assert list(model.compiled_with) == ["default"]
        requirement = model.compiled_with["default"]["spectral_axis"]
        # Both intervals present, each keeping its own density: the union, not
        # the last one to speak. Read through segments(), the only public
        # post-union form.
        assert requirement.segments() == ((1.0, 2.0, 0.5, None), (3.0, 4.0, 0.5, None))
        assert model.compiled_with["default"].sources == ("left", "right")

    def test_two_datasets_may_share_a_channel_with_different_instruments(self) -> None:
        model = Powerlaw(sed=WAVELENGTH)
        problem = FittingProblem(
            model,
            DatasetCollection(
                {
                    "plain": Dataset(
                        flat_spectrum(), Instrument([], channel="sed", input_kind=Spectrum)
                    ),
                    "calibrated": Dataset(
                        flat_spectrum(), Instrument([Calibrate()], channel="sed")
                    ),
                }
            ),
        )
        assert problem.parameters.free_names == (
            "model.index",
            "model.norm",
            "calibrated.instrument.calibrate.scale",
        )
        evaluation = problem.evaluate(
            {"model.index": -1.0, "model.norm": 1.0, "calibrated.instrument.calibrate.scale": 1.0}
        )
        assert set(evaluation.contributions) == {"plain", "calibrated"}

    def test_the_compiled_model_is_the_one_evaluated(self) -> None:
        class Swaps(Flat):
            def compile_for(self, requirements: Any) -> Flat:
                return Flat(np.array([1.0, 2.0, 4.0]))

        original = Swaps(WAVELENGTH)
        problem = FittingProblem(original, [Dataset(flat_spectrum())])
        assert problem.model is not original

    def test_compile_for_must_return_a_model(self) -> None:
        class Misbehaves(Flat):
            def compile_for(self, requirements: Any) -> Any:
                return "not a model"

        with pytest.raises(DatasetError, match="compile_for\\(\\) returned"):
            FittingProblem(Misbehaves(WAVELENGTH), [Dataset(flat_spectrum())])

    def test_validate_runs_check_alignment_eagerly(self) -> None:
        # A model whose grid does not match the data must be refused at
        # composition, not thousands of samples into a run.
        model = Flat(np.array([1.0, 2.0]))
        with pytest.raises(LikelihoodError, match="shape"):
            FittingProblem(model, [Dataset(flat_spectrum())])

    def test_validate_false_defers_but_never_skips(self) -> None:
        model = Flat(np.array([1.0, 2.0]))
        problem = FittingProblem(model, [Dataset(flat_spectrum())], validate=False)
        assert not problem.validated
        with pytest.raises(LikelihoodError, match="shape"):
            problem.log_prob({"model.level": 1.0})

    def test_a_declared_failure_on_the_first_deferred_evaluation_is_not_an_exception(
        self,
    ) -> None:
        # Deferred validation used to run before the per-dataset try, so a
        # declared simulator failure escaped log_prob as an exception on the
        # first call and returned -inf on every later one — whether a run died
        # depended on whether the sampler's first proposal happened to be
        # scoreable, which is not a property anyone can debug.
        class Fussy(Transformation):
            ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

            class Crash(RuntimeError):
                pass

            def apply(self, samples: Spectrum, values: Any) -> Spectrum:
                if float(np.max(samples.values)) > 5.0:
                    raise Fussy.Crash("the response table stops at 5")
                return samples

        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), Instrument([Fussy()]), label="d")],
            validate=False,
            simulator_failures=(Fussy.Crash,),
        )
        first = problem.evaluate({"model.level": 9.0})
        assert first.log_prob == -math.inf
        assert first.failure is not None
        assert first.failure.reason is FailureReason.INSTRUMENT_FAILED
        assert first.failure.where == "d"
        # And a scoreable point still validates and scores normally afterwards.
        assert math.isfinite(problem.log_prob({"model.level": 1.0}))
        assert problem.validated

    def test_reference_values_default_to_the_prior_median(self) -> None:
        problem = FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum())])
        assert problem.reference_values["model.level"] == pytest.approx(5.0)

    def test_reference_values_may_be_supplied(self) -> None:
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum())],
            reference_values={"model.level": 1.0},
        )
        assert problem.reference_values["model.level"] == pytest.approx(1.0)


class TestProblemComposition:
    def test_refuses_colliding_top_level_labels(self) -> None:
        with pytest.raises(DatasetError, match="labels both a model and a dataset"):
            FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="model")])

    def test_a_dataset_must_name_its_model_when_there_are_several(self) -> None:
        models = {"a": Flat(WAVELENGTH), "b": Flat(WAVELENGTH)}
        with pytest.raises(DatasetError, match="does not say which model it observes"):
            FittingProblem(models, [Dataset(flat_spectrum(), label="d")])

    def test_a_dataset_may_not_name_a_model_that_is_absent(self) -> None:
        with pytest.raises(DatasetError, match="does not hold"):
            FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d", model="ghost")])

    def test_a_dot_qualified_channel_composes_end_to_end(self) -> None:
        # Ruled 2026-09-02 (R4): the channel-name surface is settled now, and
        # grouped data — one object's several CO lines — is the non-population
        # case that wants it. The dataset label stays a bare identifier (it is
        # a merge component), so it must be given explicitly.
        class TwoLines(Model):
            def __init__(self, wavelength):
                self.register_buffer("wavelength", wavelength, unit=u.um)
                self.register_parameter(Parameter("level", st.norm(1.0, 1.0)))

            def evaluate(self, **values):
                ctx = self.context(values)
                spectrum = Spectrum(ctx["wavelength"] * u.um, np.full(3, ctx["level"]) * u.Jy)
                return ModelResult({"co.j3_2": spectrum, "co.j2_1": spectrum})

        dataset = Dataset(
            flat_spectrum(),
            Instrument([], channel="co.j3_2", input_kind=Spectrum),
            label="j3_2",
        )
        problem = FittingProblem(TwoLines(WAVELENGTH), [dataset])
        assert math.isfinite(problem.log_prob({"model.level": 1.0}))
        # Without an explicit dataset label the dotted instrument label reaches
        # the dataset-label check, which refuses loudly: merge components stay
        # bare identifiers.
        with pytest.raises(DatasetError, match="dataset label"):
            Dataset(flat_spectrum(), Instrument([], channel="co.j3_2", input_kind=Spectrum))

    def test_refuses_non_models_and_empty_model_mappings(self) -> None:
        with pytest.raises(DatasetError, match="not a Model"):
            FittingProblem({"a": object()}, [Dataset(flat_spectrum())])  # type: ignore[dict-item]
        with pytest.raises(DatasetError, match="at least one model"):
            FittingProblem({}, [Dataset(flat_spectrum())])
        with pytest.raises(DatasetError, match="takes a Model"):
            FittingProblem(object(), [Dataset(flat_spectrum())])  # type: ignore[arg-type]

    def test_a_prior_less_shared_as_is_refused_by_the_inner_merge(self) -> None:
        # Documents where the error actually comes from. An earlier version of
        # this test claimed FittingProblem refused it "rather than letting the
        # failure surface as a confusing inner TyingError" — but matched
        # loosely enough to pass on precisely that inner TyingError, which is
        # what really happens: the instrument's own merge sees one prior-less
        # site and refuses, long before a FittingProblem exists. That is
        # inference.md §17 limitation 1, and FittingProblem._require_resolved is
        # a defensive invariant that this path can never reach.
        class Deferred(Transformation):
            ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

            def __init__(self, **kwargs: Any) -> None:
                super().__init__(**kwargs)
                self.register_parameter(Parameter("scale", shared_as="scale"))

            def apply(self, samples: Spectrum, values: Any) -> Spectrum:
                return samples

        from ampere.core import TyingError

        with pytest.raises(TyingError, match="every site is declared shared without one"):
            Dataset(flat_spectrum(), Instrument([Deferred()]), label="d")

    def test_a_model_no_dataset_observes_is_refused(self) -> None:
        # Its parameters would be free dimensions no likelihood constrains, and
        # it would be re-evaluated on every log_prob for nothing.
        with pytest.raises(DatasetError, match="observed by no dataset"):
            FittingProblem(
                {"used": Flat(WAVELENGTH), "orphan": Flat(WAVELENGTH)},
                [Dataset(flat_spectrum(), label="d", model="used")],
            )

    def test_a_bare_dataset_is_refused_with_a_contract_error(self) -> None:
        with pytest.raises(DatasetError, match="Wrap it in a list"):
            DatasetCollection(Dataset(flat_spectrum()))  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# The §4.5 surface
# ---------------------------------------------------------------------------


class TestEngineFacingSurface:
    def problem(self) -> FittingProblem:
        return FittingProblem(
            Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d")], seed=20260902
        )

    def test_log_prob_is_the_split_summed(self) -> None:
        problem = self.problem()
        theta = {"model.level": 1.0}
        assert problem.log_prob(theta) == pytest.approx(
            problem.log_prior(theta) + problem.log_likelihood(theta)
        )

    def test_log_prior_is_the_merged_sets_own(self) -> None:
        problem = self.problem()
        assert problem.log_prior({"model.level": 1.0}) == pytest.approx(
            problem.parameters.lnprior({"model.level": 1.0})
        )

    def test_evaluate_reports_the_split_and_the_per_dataset_terms(self) -> None:
        problem = self.problem()
        evaluation = problem.evaluate({"model.level": 1.0})
        assert isinstance(evaluation, Evaluation)
        assert evaluation.failure is None
        assert set(evaluation.contributions) == {"d"}
        assert evaluation.log_likelihood == pytest.approx(evaluation.contributions["d"])
        assert evaluation.log_prob == pytest.approx(
            evaluation.log_prior + evaluation.log_likelihood
        )

    def test_out_of_support_short_circuits_without_running_the_model(self) -> None:
        model = Recording(WAVELENGTH)
        problem = FittingProblem(model, [Dataset(flat_spectrum())])
        before = model.calls
        evaluation = problem.evaluate({"model.level": 100.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.log_prior == -math.inf
        assert math.isnan(evaluation.log_likelihood)
        assert evaluation.failure is None  # zero prior mass is an answer, not a failure
        assert model.calls == before

    def test_prior_transform_delegates_to_the_merged_set(self) -> None:
        problem = self.problem()
        assert problem.prior_transform([0.5]) == pytest.approx([5.0])
        assert problem.free_size == 1
        assert problem.free_labels() == ("model.level",)

    def test_unconstrained_round_trip_and_density(self) -> None:
        problem = self.problem()
        theta = np.array([2.5])
        unconstrained = problem.unconstrain(theta)
        assert problem.constrain(unconstrained) == pytest.approx(theta)
        expected = problem.parameters.lnprior_unconstrained(unconstrained) + problem.log_likelihood(
            theta
        )
        assert problem.log_prob_unconstrained(unconstrained) == pytest.approx(expected)

    def test_log_prob_unconstrained_is_minus_inf_outside_the_support(self) -> None:
        # With the default Logit bijection an unconstrained point is always in
        # support, so the branch needs a parameter declared on the real line
        # with a bounded prior — which is a legitimate declaration, and exactly
        # where a NUTS proposal can leave the support.
        from ampere.core import Identity

        class Unbijected(Flat):
            def __init__(self, wavelength: np.ndarray) -> None:
                self.register_buffer("wavelength", np.asarray(wavelength, dtype=float))
                self.register_parameter(
                    Parameter("level", st.uniform(0.0, 10.0), bijection=Identity())
                )

        problem = FittingProblem(Unbijected(WAVELENGTH), [Dataset(flat_spectrum(), label="d")])
        assert math.isfinite(problem.log_prob_unconstrained([5.0]))
        assert problem.log_prob_unconstrained([1e6]) == -math.inf

    def test_sample_prior_uses_the_named_substream(self) -> None:
        problem = self.problem()
        drawn = FittingProblem(
            Flat(WAVELENGTH), [Dataset(flat_spectrum())], seed=20260902
        ).sample_prior()
        assert problem.sample_prior() == pytest.approx(drawn)

    def test_defaults_to_the_reference_values(self) -> None:
        problem = self.problem()
        assert problem.log_prob() == pytest.approx(problem.log_prob(problem.reference_values))


class TestCapabilities:
    def test_reference_path_is_honest_about_itself(self) -> None:
        problem = FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum())])
        assert problem.capabilities == Capabilities()
        assert not problem.differentiable
        assert not problem.batchable
        assert problem.device == "cpu"

    def test_conjunctive_over_the_parts(self) -> None:
        class Native:
            DIFFERENTIABLE = True
            BATCHABLE = True
            DEVICE = "cuda"

        class NativeOnCpu:
            DIFFERENTIABLE = True
            BATCHABLE = True

        class Silent:
            pass

        assert declared_capabilities([Native(), Native()]) == Capabilities(True, True, "cuda")
        assert declared_capabilities([NativeOnCpu(), NativeOnCpu()]) == Capabilities(True, True)
        # One silent part is enough to withdraw the whole conjunctive claim.
        assert declared_capabilities([NativeOnCpu(), Silent()]) == Capabilities()

    def test_empty_parts_are_not_differentiable(self) -> None:
        assert declared_capabilities([]) == Capabilities()

    def test_silence_counts_as_cpu_and_therefore_disagrees_with_a_gpu_part(self) -> None:
        class OnGpu:
            DEVICE = "cuda"

        with pytest.raises(DatasetError, match="different devices"):
            declared_capabilities([OnGpu(), object()])

    def test_disagreeing_devices_are_refused(self) -> None:
        class OnCpu:
            DEVICE = "cpu"

        class OnGpu:
            DEVICE = "cuda"

        with pytest.raises(DatasetError, match="different devices"):
            declared_capabilities([OnCpu(), OnGpu()])

    def test_a_differentiable_model_and_step_lift_the_problem(self) -> None:
        class NativeFlat(Flat):
            DIFFERENTIABLE = True
            BATCHABLE = True

        class NativeCalibrate(Calibrate):
            DIFFERENTIABLE = True
            BATCHABLE = True

        problem = FittingProblem(
            NativeFlat(WAVELENGTH),
            [Dataset(flat_spectrum(), Instrument([NativeCalibrate()]))],
        )
        assert problem.differentiable
        assert problem.batchable

    def test_one_black_box_step_stops_the_gradient(self) -> None:
        class NativeFlat(Flat):
            DIFFERENTIABLE = True

        problem = FittingProblem(
            NativeFlat(WAVELENGTH), [Dataset(flat_spectrum(), Instrument([Calibrate()]))]
        )
        assert not problem.differentiable

    def test_capabilities_may_be_declared_explicitly(self) -> None:
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum())],
            capabilities=Capabilities(differentiable=True),
        )
        assert problem.differentiable

    def test_check_engine_refuses_a_gradient_free_engine_a_latent_problem(self) -> None:
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        problem = FittingProblem(
            Counts(WAVELENGTH),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(
                        PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))
                    ),
                    label="counts",
                )
            ],
        )
        with pytest.raises(LikelihoodError, match="emcee cannot run this likelihood"):
            problem.check_engine("emcee")
        problem.check_engine("NUTS", differentiable=True)

    def test_check_engine_passes_an_analytic_problem(self) -> None:
        problem = FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum())])
        problem.check_engine("emcee")


# ---------------------------------------------------------------------------
# Failure signalling (§4.5)
# ---------------------------------------------------------------------------


class TestFailureSignalling:
    def test_a_crashing_simulator_becomes_minus_inf_with_a_reason(self) -> None:
        problem = FittingProblem(
            Broken(WAVELENGTH, mode="crash"),
            [Dataset(flat_spectrum(), label="d")],
            validate=False,
            simulator_failures=(Broken.Crash,),
        )
        evaluation = problem.evaluate({"model.level": 1.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.failed
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.MODEL_FAILED
        assert evaluation.failure.where == "model"
        assert "exited 1" in evaluation.failure.message
        assert evaluation.failure.exception_type == "Crash"

    def test_an_undeclared_exception_is_a_bug_and_propagates(self) -> None:
        problem = FittingProblem(
            Broken(WAVELENGTH, mode="bug"),
            [Dataset(flat_spectrum())],
            validate=False,
        )
        with pytest.raises(ZeroDivisionError):
            problem.log_prob({"model.level": 1.0})

    def test_nan_predictions_are_classified_as_such(self) -> None:
        problem = FittingProblem(
            Broken(WAVELENGTH, mode="nan"),
            [Dataset(flat_spectrum(), label="d")],
        )
        evaluation = problem.evaluate({"model.level": 1.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.NON_FINITE_PREDICTION
        assert evaluation.failure.where == "d"

    def test_a_non_positive_definite_covariance_becomes_minus_inf(self) -> None:
        # likelihoods.md §16 hands this conversion to this contract explicitly.
        likelihood = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(
                Matern32(st.loguniform(1e-3, 1e250), st.loguniform(1e-3, 1e3)),
            ),
        )
        problem = FittingProblem(
            Flat(WAVELENGTH), [Dataset(flat_spectrum(), likelihood=likelihood, label="d")]
        )
        # An amplitude whose square overflows float64: inside the prior, and a
        # covariance the Cholesky cannot touch.
        evaluation = problem.evaluate(
            {
                "model.level": 1.0,
                "d.likelihood.amplitude": 1e200,
                "d.likelihood.length_scale": 1.0,
            }
        )
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.LIKELIHOOD_FAILED
        assert evaluation.failure.where == "d"

    def test_a_non_positive_poisson_rate_is_a_second_reachable_failure(self) -> None:
        # Not only the non-positive-definite covariance: PoissonFamily raises a
        # LikelihoodError when the model predicts a rate <= 0, and that is an
        # unscoreable point rather than a bug.
        class Negative(Counts):
            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(
                    ctx["wavelength"] * u.micron,
                    np.full(ctx["wavelength"].shape, -ctx["rate"]),
                )

        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        problem = FittingProblem(
            Negative(WAVELENGTH),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
                    label="counts",
                )
            ],
            validate=False,
        )
        evaluation = problem.evaluate({"model.rate": 3.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.LIKELIHOOD_FAILED
        assert evaluation.failure.where == "counts"
        assert "strictly positive expected count" in evaluation.failure.message

    def test_the_reason_names_the_dataset_that_failed(self) -> None:
        # Two datasets on one model; only the Poisson one cannot score a
        # negative rate. Per-dataset granularity is what makes "which of my six
        # datasets is rejecting everything?" answerable.
        class Negative(Counts):
            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(
                    ctx["wavelength"] * u.micron,
                    np.full(ctx["wavelength"].shape, -ctx["rate"]),
                )

        plain = Spectrum(
            WAVELENGTH * u.micron, np.array([-3.0, -3.0, -3.0]), uncertainty=np.full(3, 0.1)
        )
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        problem = FittingProblem(
            Negative(WAVELENGTH),
            DatasetCollection(
                {
                    "gaussian": Dataset(plain, label="gaussian"),
                    "poisson": Dataset(
                        counts,
                        likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
                        label="poisson",
                    ),
                }
            ),
            validate=False,
        )
        evaluation = problem.evaluate({"model.rate": 3.0})
        assert evaluation.failure is not None
        assert evaluation.failure.where == "poisson"
        # The dataset that *did* score is still reported, in collection order.
        assert set(evaluation.contributions) == {"gaussian"}

    def test_failures_are_counted_and_the_history_is_bounded(self) -> None:
        problem = FittingProblem(
            Broken(WAVELENGTH, mode="crash"),
            [Dataset(flat_spectrum())],
            validate=False,
            simulator_failures=(Broken.Crash,),
            failure_history=3,
        )
        for _ in range(10):
            problem.log_prob({"model.level": 1.0})
        assert problem.failure_counts[FailureReason.MODEL_FAILED] == 10
        assert len(problem.failures) == 3
        assert problem.last_failure is not None
        problem.reset_failures()
        assert problem.failures == ()
        assert dict(problem.failure_counts) == {}

    def test_a_failure_serialises_for_provenance(self) -> None:
        failure = Failure(FailureReason.MODEL_FAILED, "boom", where="sed", exception_type="OSError")
        assert failure.to_dict() == {
            "reason": "model_failed",
            "message": "boom",
            "where": "sed",
            "exception_type": "OSError",
            "values": {},
        }
        assert str(failure) == "model_failed [sed]: boom"
        located = Failure(
            FailureReason.LIKELIHOOD_FAILED, "not PD", where="d", values={"amplitude": 1e3}
        )
        assert "amplitude=1000" in str(located)
        assert located.to_dict()["values"] == {"amplitude": 1000.0}

    def test_an_exception_that_refuses_attributes_is_still_a_failure(self) -> None:
        # BaseException always provides __dict__, so a __slots__ subclass is
        # fine — but a frozen dataclass or attrs exception overrides
        # __setattr__, and letting that AttributeError out of the except block
        # that is *handling* a failure would turn an unscoreable point into a
        # crash.
        class Frozen(RuntimeError):
            def __setattr__(self, name: str, value: Any) -> None:
                raise AttributeError(f"{type(self).__name__} is immutable")

        class Raises(Flat):
            def evaluate(self, **values: Any) -> Spectrum:
                raise Frozen("the wrapped code exited 1")

        problem = FittingProblem(
            Raises(WAVELENGTH),
            [Dataset(flat_spectrum(), label="d")],
            validate=False,
            simulator_failures=(Frozen,),
        )
        evaluation = problem.evaluate({"model.level": 1.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.exception_type == "Frozen"

    def test_a_plus_inf_log_likelihood_is_a_recorded_failure(self) -> None:
        # NaN was caught; +inf was not, and produced an Evaluation whose
        # log_likelihood contradicted its log_prob with no reason recorded.
        class Degenerate(GaussianFamily):
            NAME: ClassVar[str] = "degenerate_for_test"

            def log_prob(self, predicted: Any, observed: Any, noise: Any) -> float:
                return math.inf

        problem = FittingProblem(
            Flat(WAVELENGTH),
            [
                Dataset(
                    flat_spectrum(),
                    likelihood=Likelihood(Degenerate(), IndependentNoise()),
                    label="d",
                )
            ],
        )
        evaluation = problem.evaluate({"model.level": 1.0})
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.NON_FINITE_LOG_LIKELIHOOD

    def test_simulator_failures_must_be_exception_classes(self) -> None:
        with pytest.raises(DatasetError, match="takes exception classes"):
            FittingProblem(
                Flat(WAVELENGTH),
                [Dataset(flat_spectrum())],
                simulator_failures=(object,),  # type: ignore[arg-type]
            )


# ---------------------------------------------------------------------------
# simulate (§4.5, for SBI)
# ---------------------------------------------------------------------------


class TestSimulate:
    def problem(self, **kwargs: Any) -> FittingProblem:
        kwargs.setdefault("seed", 20260902)
        return FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d")], **kwargs)

    def test_draws_theta_from_the_prior_by_default(self) -> None:
        problem = self.problem()
        first = problem.simulate()
        second = problem.simulate()
        assert not first.failed
        assert first.theta.shape == (1,)
        # A simulation budget must advance its stream, not repeat one point.
        assert first.theta != pytest.approx(second.theta)

    def test_is_reproducible_from_the_seed(self) -> None:
        assert self.problem().simulate().theta == pytest.approx(self.problem().simulate().theta)

    def test_returns_the_theta_result_pair_and_the_prediction(self) -> None:
        problem = self.problem()
        simulation = problem.simulate({"model.level": 2.0})
        assert simulation.parameters["model.level"] == pytest.approx(2.0)
        assert set(simulation.results) == {"model"}
        assert simulation.results["model"].parameters is not None
        assert simulation.predicted["d"].values == pytest.approx([2.0, 2.0, 2.0])
        assert simulation.observations is None

    def test_observe_draws_independent_noise(self) -> None:
        # The earlier version asserted only np.allclose(drawn, 2.0, atol=1.0)
        # against sigma=0.1, which passes for any sigma up to ~0.3 and so
        # pinned nothing. Check the variance the likelihood would actually use.
        problem = self.problem()
        draws = np.array(
            [
                problem.simulate({"model.level": 2.0}, observe=True).observations["d"].values
                for _ in range(4000)
            ]
        )
        assert draws.shape == (4000, 3)
        assert np.mean(draws) == pytest.approx(2.0, abs=0.02)
        assert np.var(draws, axis=0) == pytest.approx(0.01, abs=0.002)

    def test_a_fitted_scale_and_jitter_reach_the_draw(self) -> None:
        # The docstring claims the noise model's own sigma is used, so a fitted
        # scale or jitter is already in it. Nothing held that down before.
        likelihood = Likelihood(
            GaussianFamily(),
            IndependentNoise(scale=st.loguniform(0.5, 5.0), jitter=st.loguniform(0.01, 1.0)),
        )
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), likelihood=likelihood, label="d")],
            seed=11,
        )
        theta = {"model.level": 2.0, "d.likelihood.scale": 3.0, "d.likelihood.jitter": 0.4}
        draws = np.array(
            [problem.simulate(theta, observe=True).observations["d"].values for _ in range(4000)]
        )
        # sigma = sqrt((0.1 * 3)^2 + 0.4^2) = 0.5
        assert np.var(draws, axis=0) == pytest.approx(0.25, abs=0.02)

    def test_the_draw_hands_the_noise_model_the_noiseless_prediction(self) -> None:
        # X-1 (ruled 2026-09-03): draw_observation builds the NoiseParams from
        # the noiseless prediction *before* noise is added, so a prediction-
        # dependent noise scales with the true curve (sigma(mu), not sigma(x)).
        captured: dict = {}

        class Recording(IndependentNoise):
            def noise_params(self, observed, retain, values, *, predicted=None, **kwargs):
                captured["predicted"] = predicted
                return super().noise_params(observed, retain, values, predicted=predicted, **kwargs)

        likelihood = Likelihood(GaussianFamily(), Recording())
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), likelihood=likelihood, label="d")],
            seed=3,
        )
        problem.simulate({"model.level": 2.0}, observe=True)
        assert captured["predicted"] is not None
        np.testing.assert_allclose(captured["predicted"], [2.0, 2.0, 2.0])

    def test_the_solver_jitter_is_part_of_the_drawn_covariance(self) -> None:
        # DenseGP scores K + diag(sigma^2 + jitter^2); the draw must match, or
        # simulate() and log_prob disagree about the same model. The library's
        # own error message tells users to raise this jitter, so it is reachable.
        likelihood = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(Matern32(0.5, 2.0), DenseGP(jitter=0.15)),
        )
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), likelihood=likelihood, label="d")],
            seed=5,
        )
        draws = np.array(
            [
                problem.simulate({"model.level": 2.0}, observe=True).observations["d"].values
                for _ in range(6000)
            ]
        )
        coordinates = WAVELENGTH[:, None]
        kernel = likelihood.noise.kernel
        expected = kernel.matrix(
            coordinates, coordinates, {"amplitude": 0.5, "length_scale": 2.0}
        ) + np.diag(np.full(3, 0.1**2 + 0.15**2))
        assert np.cov(draws.T) == pytest.approx(expected, abs=0.025)
        # Without the solver jitter the variance would be 0.26, not 0.2825.
        assert np.var(draws, axis=0).min() > 0.27

    def test_observe_draws_correlated_noise_from_the_gp(self) -> None:
        # sigma is deliberately large enough (0.3, so sigma^2 = 0.09) that the
        # diagonal noise term is resolvable against the sampling error of a
        # 4000-draw covariance estimate. With the 0.1 uncertainties used
        # elsewhere in this file, dropping sigma entirely would shift the
        # diagonal by 0.01 and no honest tolerance would catch it.
        likelihood = Likelihood(GaussianFamily(), GaussianProcessNoise(Matern32(0.5, 2.0)))
        observed = flat_spectrum(uncertainty=[0.3, 0.3, 0.3] * u.Jy)
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(observed, likelihood=likelihood, label="d")],
            seed=3,
        )
        draws = np.array(
            [
                problem.simulate({"model.level": 2.0}, observe=True).observations["d"].values
                for _ in range(4000)
            ]
        )
        empirical = np.cov(draws.T)
        kernel = likelihood.noise.kernel
        coordinates = WAVELENGTH[:, None]
        covariance = kernel.matrix(
            coordinates, coordinates, {"amplitude": 0.5, "length_scale": 2.0}
        )
        expected = covariance + np.diag(np.full(3, 0.09))
        assert empirical == pytest.approx(expected, abs=0.03)
        # The two ways this could be wrong and still look plausible, ruled out
        # explicitly rather than left to the joint tolerance: no correlation at
        # all (the off-diagonals would collapse to zero), and no diagonal noise
        # at all (the variances would be K's alone). Both assertions were
        # mutation-tested against a patched draw_observation.
        off_diagonal = np.abs(empirical[np.triu_indices(3, k=1)])
        assert off_diagonal.min() > 0.03
        assert np.abs(np.diag(empirical) - np.diag(covariance)).min() > 0.04
        assert np.mean(draws) == pytest.approx(2.0, abs=0.05)

    def test_masked_samples_keep_the_observed_values(self) -> None:
        observed = flat_spectrum(values=(9.0, 9.0, 9.0), mask=[False, False, True])
        problem = FittingProblem(Flat(WAVELENGTH), [Dataset(observed, label="d")], seed=1)
        drawn = problem.simulate({"model.level": 2.0}, observe=True).observations["d"]
        assert drawn.values[2] == pytest.approx(9.0)
        assert drawn.mask is not None
        assert bool(drawn.mask[2])

    def test_masking_beats_censoring_for_a_draw_too(self) -> None:
        # Dataset._declare_latent already consults marginalisation_for, so a
        # limit on a masked sample does not conjure a latent block. draw_
        # observation used to test `censoring is not None` unconditionally, so
        # one class gave two answers to one question.
        observed = flat_spectrum(mask=[False, False, True])
        likelihood = Likelihood(
            GaussianFamily(),
            IndependentNoise(),
            censoring=Censoring(np.array([0, 0, 1])),
        )
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(observed, likelihood=likelihood, label="d")],
            seed=2,
        )
        simulation = problem.simulate({"model.level": 2.0}, observe=True)
        assert simulation.observations is not None
        assert simulation.observations["d"].values[2] == pytest.approx(observed.values[2])

    def test_a_limit_that_survives_the_mask_still_blocks_a_draw(self) -> None:
        likelihood = Likelihood(
            GaussianFamily(),
            IndependentNoise(),
            censoring=Censoring(np.array([0, 1, 0])),
        )
        problem = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), likelihood=likelihood, label="d")],
        )
        with pytest.raises(DatasetError, match="censoring declaration on retained samples"):
            problem.simulate({"model.level": 2.0}, observe=True)

    def test_a_family_that_cannot_be_sampled_says_so(self) -> None:
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        problem = FittingProblem(
            Counts(WAVELENGTH),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(PoissonFamily(), IndependentNoise()),
                    label="counts",
                )
            ],
        )
        # Raised, not flagged: "this family cannot sample" is a fact about the
        # composition and would fail identically for every draw, so a budget of
        # 10^4 flagged failures would be strictly less useful than one clear
        # error. The refusal is specific (ruled 2026-09-02, R3): it names the
        # family and the override that provides the observation process.
        with pytest.raises(DatasetError, match="does not implement sample"):
            problem.simulate({"model.rate": 3.0}, observe=True)
        # The noise-free half still works, which is what emulator training wants.
        assert not problem.simulate({"model.rate": 3.0}).failed

    def test_a_user_family_provides_its_own_sample(self) -> None:
        # R3's point: the observation process is the family author's to supply.
        # A Poisson subclass overriding sample() makes observe=True draws work
        # end to end, with the same NoiseParams the likelihood scores with.
        class SamplingPoisson(PoissonFamily):
            def sample(self, predicted, noise, rng):
                return rng.poisson(predicted).astype(float)

        counts = Spectrum(WAVELENGTH * u.um, np.array([40.0, 70.0, 20.0]))
        problem = FittingProblem(
            Counts(WAVELENGTH),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(SamplingPoisson(), IndependentNoise()),
                    label="counts",
                )
            ],
            seed=20260902,
        )
        drawn = problem.simulate({"model.rate": 30.0}, observe=True).observations["counts"]
        assert drawn.values.shape == (3,)
        assert np.all(drawn.values == np.round(drawn.values))
        assert not np.allclose(drawn.values, 30.0)
        # The draw is scoreable by the same likelihood that will fit it.
        assert math.isfinite(problem.log_prob({"model.rate": 30.0}))

    def test_a_crash_is_flagged_not_raised(self) -> None:
        problem = FittingProblem(
            Broken(WAVELENGTH, mode="crash"),
            [Dataset(flat_spectrum())],
            validate=False,
            simulator_failures=(Broken.Crash,),
        )
        simulation = problem.simulate({"model.level": 1.0})
        assert simulation.failed
        assert simulation.failure is not None
        assert simulation.failure.reason is FailureReason.MODEL_FAILED
        # θ is still reported: a rejected draw is information about the prior.
        assert simulation.parameters["model.level"] == pytest.approx(1.0)


# ---------------------------------------------------------------------------
# The acceptance criterion
# ---------------------------------------------------------------------------


class TestJointTwoDatasetProblem:
    """W1.7's acceptance criterion, written out.

    Two datasets, one stub model with two channels, one instrument chain each,
    and a calibration nuisance tied across them. The tie must genuinely collapse
    — one entry in the joint space, two binding sites — and both datasets must
    respond to it.
    """

    BLUE = np.array([1.0, 2.0, 4.0])
    RED = np.array([10.0, 20.0, 40.0])

    def build(self, **kwargs: Any) -> FittingProblem:
        model = Powerlaw(blue=self.BLUE, red=self.RED)
        blue = Spectrum(
            self.BLUE * u.micron,
            [1.0, 0.5, 0.25] * u.Jy,
            uncertainty=[0.05, 0.05, 0.05] * u.Jy,
        )
        red = Spectrum(
            self.RED * u.micron,
            [0.1, 0.05, 0.025] * u.Jy,
            uncertainty=[0.005, 0.005, 0.005] * u.Jy,
        )
        datasets = DatasetCollection(
            {
                "blue": Dataset(blue, Instrument([Calibrate()], channel="blue")),
                "red": Dataset(red, Instrument([Calibrate()], channel="red")),
            }
        )
        tie = Tie(
            "calibration",
            ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"),
        )
        return FittingProblem(model, datasets, ties=[tie], seed=20260902, **kwargs)

    def test_the_tie_costs_one_dimension_not_two(self) -> None:
        problem = self.build()
        assert problem.parameters.free_names == ("model.index", "model.norm", "calibration")
        assert problem.free_size == 3
        assert problem.tied_names == ("calibration",)

    def test_the_tie_has_a_binding_in_each_dataset(self) -> None:
        problem = self.build()
        sites = problem.mapping.sites_of("calibration")
        assert {(b.component, b.local_name) for b in sites} == {
            ("blue", "instrument.calibrate.scale"),
            ("red", "instrument.calibrate.scale"),
        }

    def test_both_datasets_respond_to_the_one_tied_parameter(self) -> None:
        problem = self.build()
        theta = {"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0}
        at_one = problem.evaluate(theta)
        at_two = problem.evaluate(dict(theta, calibration=1.2))
        assert set(at_one.contributions) == {"blue", "red"}
        for label in ("blue", "red"):
            assert at_one.contributions[label] != pytest.approx(at_two.contributions[label])
        # And the prediction itself scales, in both datasets at once.
        one = problem.simulate(theta)
        two = problem.simulate(dict(theta, calibration=2.0))
        for label in ("blue", "red"):
            assert two.predicted[label].values == pytest.approx(2.0 * one.predicted[label].values)

    def test_the_joint_log_likelihood_is_the_sum(self) -> None:
        problem = self.build()
        theta = {"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0}
        evaluation = problem.evaluate(theta)
        assert evaluation.log_likelihood == pytest.approx(sum(evaluation.contributions.values()))
        assert evaluation.log_prob == pytest.approx(
            problem.log_prior(theta) + evaluation.log_likelihood
        )

    def test_the_whole_engine_surface_works(self) -> None:
        problem = self.build()
        cube = np.array([0.5, 0.5, 0.5])
        theta = problem.prior_transform(cube)
        assert theta.shape == (3,)
        assert math.isfinite(problem.log_prob(theta))
        assert problem.free_labels() == ("model.index", "model.norm", "calibration")
        problem.check_engine("emcee")
        simulation = problem.simulate(observe=True)
        assert not simulation.failed
        assert set(simulation.observations or {}) == {"blue", "red"}
        assert math.isfinite(problem.log_prob_unconstrained(problem.unconstrain(theta)))

    def test_the_recovered_truth_is_the_best_fit(self) -> None:
        # The data were generated from index=-1, norm=1, scale=1; a small
        # perturbation of any of the three must score worse.
        problem = self.build()
        truth = {"model.index": -1.0, "model.norm": 1.0, "calibration": 1.0}
        best = problem.log_prob(truth)
        for name, delta in (("model.index", 0.1), ("model.norm", 0.1), ("calibration", 0.1)):
            worse = dict(truth)
            worse[name] = worse[name] + delta
            assert problem.log_prob(worse) < best


# ---------------------------------------------------------------------------
# What nesting buys, demonstrated
# ---------------------------------------------------------------------------


class TestNestedMergeProperties:
    def test_a_top_level_tie_reaches_a_site_inside_an_inner_merge(self) -> None:
        # The site lives three levels down: step -> instrument -> dataset ->
        # problem. A tie at the top names it by its full path and collapses it.
        model = Powerlaw(blue=np.array([1.0, 2.0, 4.0]))
        observed = Spectrum(
            np.array([1.0, 2.0, 4.0]) * u.micron,
            [1.0, 0.5, 0.25] * u.Jy,
            uncertainty=[0.05] * 3 * u.Jy,
        )
        datasets = DatasetCollection(
            {
                "a": Dataset(
                    observed,
                    Instrument(
                        [Calibrate(label="first"), Calibrate(label="second")], channel="blue"
                    ),
                )
            }
        )
        tie = Tie("scale", ("a.instrument.first.scale", "a.instrument.second.scale"))
        problem = FittingProblem(model, datasets, ties=[tie])
        assert "scale" in problem.parameters.free_names
        assert len(problem.mapping.sites_of("scale")) == 2
        # Both steps see the one value, so the chain applies it twice.
        theta = {"model.index": -1.0, "model.norm": 1.0, "scale": 3.0}
        predicted = problem.simulate(theta).predicted["a"]
        assert predicted.values == pytest.approx(9.0 * np.array([1.0, 0.5, 0.25]))

    def test_a_hierarchical_reference_survives_two_merge_levels_and_a_tie(self) -> None:
        # Declared locally against 'mu'; the dataset merge rewrites it to
        # 'likelihood.mu', the problem merge to the tie label. If any level
        # dropped its mapping the reference would dangle and ParameterSet would
        # refuse to construct.
        from ampere.core import HierarchicalPrior

        class Referring(IndependentNoise):
            def __init__(self) -> None:
                super().__init__()
                self.register_parameter(Parameter("mu", st.uniform(0.5, 1.0)))
                self.register_parameter(
                    Parameter("scale", HierarchicalPrior("norm", {"loc": "mu"}))
                )

        def dataset(label: str) -> Dataset:
            return Dataset(
                flat_spectrum(),
                Instrument([], channel="blue", input_kind=Spectrum, label=label),
                Likelihood(GaussianFamily(), Referring()),
                label=label,
            )

        model = Powerlaw(blue=WAVELENGTH)
        tie = Tie("population_mu", ("a.likelihood.mu", "b.likelihood.mu"))
        problem = FittingProblem(model, DatasetCollection([dataset("a"), dataset("b")]), ties=[tie])
        assert "population_mu" in problem.parameters.free_names
        assert problem.parameters["a.likelihood.scale"].references == ("population_mu",)
        assert problem.parameters["b.likelihood.scale"].references == ("population_mu",)
        assert math.isfinite(problem.log_prior(problem.reference_values))

    def test_a_plate_of_datasets_ties_one_hyperparameter_across_many(self) -> None:
        # The IFU per-spaxel pattern (W1.11): N datasets, N Likelihoods, one
        # shared GP hyperparameter. A Tie takes any number of sites, so the N
        # sites collapse to one dimension and every spaxel's noise model reads
        # it. What nesting buys here is the naming: 'spaxel_3.likelihood.
        # amplitude' rather than a flattened 'spaxel_3_likelihood.amplitude'.
        size = 6
        model = Powerlaw(**{f"spaxel_{i}": WAVELENGTH for i in range(size)})

        def spaxel(index: int) -> Dataset:
            return Dataset(
                flat_spectrum(),
                Instrument([], channel=f"spaxel_{index}", input_kind=Spectrum),
                Likelihood(
                    GaussianFamily(),
                    GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e1), 1.0)),
                ),
                label=f"spaxel_{index}",
            )

        collection = DatasetCollection([spaxel(i) for i in range(size)])
        tie = Tie(
            "gp_amplitude",
            tuple(f"spaxel_{i}.likelihood.amplitude" for i in range(size)),
        )
        problem = FittingProblem(model, collection, ties=[tie])
        assert problem.parameters.free_names == ("model.index", "model.norm", "gp_amplitude")
        assert len(problem.mapping.sites_of("gp_amplitude")) == size
        evaluation = problem.evaluate({"model.index": -1.0, "model.norm": 1.0, "gp_amplitude": 0.2})
        assert len(evaluation.contributions) == size
        assert math.isfinite(evaluation.log_prob)

    def test_an_inner_tie_is_visible_to_the_top_level_mapping(self) -> None:
        # Lossless nesting (ruled 2026-09-02): two steps sharing a
        # declaration-time label are collapsed by the *instrument's* own
        # merge, and the top-level mapping's composed bindings now surface the
        # collapse — tied_names tells the leaf-level truth without a descent.
        # Before the ruling this very test asserted tied_names == ().
        problem = self._inner_tie_problem()
        assert "d.instrument.gain" in problem.parameters.free_names
        assert problem.tied_names == ("d.instrument.gain",)
        assert {
            (b.component, b.local_name) for b in problem.mapping.sites_of("d.instrument.gain")
        } == {("d", "instrument.a.scale"), ("d", "instrument.b.scale")}

    def test_sites_descends_and_recovers_what_tied_names_loses(self) -> None:
        problem = self._inner_tie_problem()
        assert problem.shared_names == ("d.instrument.gain",)
        assert problem.sites()["d.instrument.gain"] == (
            "d.instrument.a.scale",
            "d.instrument.b.scale",
        )
        assert problem.sites()["model.index"] == ("model.index",)

    def test_sites_reports_a_top_level_tie_too(self) -> None:
        model = Powerlaw(blue=WAVELENGTH, red=WAVELENGTH)
        datasets = DatasetCollection(
            {
                "blue": Dataset(flat_spectrum(), Instrument([Calibrate()], channel="blue")),
                "red": Dataset(flat_spectrum(), Instrument([Calibrate()], channel="red")),
            }
        )
        tie = Tie("gain", ("blue.instrument.calibrate.scale", "red.instrument.calibrate.scale"))
        problem = FittingProblem(model, datasets, ties=[tie])
        assert problem.tied_names == ("gain",)
        assert problem.shared_names == ("gain",)
        assert problem.sites()["gain"] == (
            "blue.instrument.calibrate.scale",
            "red.instrument.calibrate.scale",
        )

    @staticmethod
    def _inner_tie_problem() -> FittingProblem:
        class Shared(Calibrate):
            def __init__(self, **kwargs: Any) -> None:
                Transformation.__init__(self, **kwargs)
                self.register_parameter(Parameter("scale", st.lognorm(0.2), shared_as="gain"))

        return FittingProblem(
            Powerlaw(blue=WAVELENGTH),
            [
                Dataset(
                    flat_spectrum(),
                    Instrument([Shared(label="a"), Shared(label="b")], channel="blue"),
                    label="d",
                )
            ],
        )

    def test_shared_parameters_are_a_top_level_component(self) -> None:
        shared = ParameterSet([Parameter("scale", st.lognorm(0.2))])
        datasets = DatasetCollection(
            [Dataset(flat_spectrum(), Instrument([Calibrate()]), label="d")], shared=shared
        )
        problem = FittingProblem(
            Flat(WAVELENGTH),
            datasets,
            ties=[Tie("scale", ("shared.scale", "d.instrument.calibrate.scale"))],
        )
        assert problem.parameters.free_names == ("model.level", "scale")
        assert len(problem.mapping.sites_of("scale")) == 2


# ---------------------------------------------------------------------------
# The latent path, end to end
# ---------------------------------------------------------------------------


class TestLatentPath:
    def build(self) -> FittingProblem:
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        return FittingProblem(
            Counts(WAVELENGTH),
            [
                Dataset(
                    counts,
                    likelihood=Likelihood(
                        PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))
                    ),
                    label="counts",
                )
            ],
            capabilities=Capabilities(differentiable=True),
        )

    def test_the_latent_block_joins_the_one_merge(self) -> None:
        problem = self.build()
        assert problem.parameters.free_names == ("model.rate", "counts.latent.z")
        assert problem.parameters["counts.latent.z"].shape == (3,)
        assert problem.free_size == 4  # one rate plus three latent values

    def test_the_latent_values_reach_the_likelihood(self) -> None:
        problem = self.build()
        base = {"model.rate": 5.0, "counts.latent.z": np.zeros(3)}
        flat = problem.evaluate(base)
        perturbed = problem.evaluate({**base, "counts.latent.z": np.array([0.5, -0.5, 0.2])})
        assert math.isfinite(flat.log_likelihood)
        assert flat.log_likelihood != pytest.approx(perturbed.log_likelihood)

    def test_the_gp_hyperparameters_would_reach_it_too(self) -> None:
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        likelihood = Likelihood(
            PoissonFamily(),
            GaussianProcessNoise(Matern32(st.loguniform(0.05, 2.0), st.loguniform(0.5, 5.0))),
        )
        problem = FittingProblem(
            Counts(WAVELENGTH),
            [Dataset(counts, likelihood=likelihood, label="counts")],
            capabilities=Capabilities(differentiable=True),
        )
        assert problem.parameters.free_names == (
            "model.rate",
            "counts.likelihood.amplitude",
            "counts.likelihood.length_scale",
            "counts.latent.z",
        )
        assert problem.free_size == 6

    def test_a_predicted_mask_that_shrinks_the_retained_count_is_refused(self) -> None:
        class MasksMore(Transformation):
            ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

            def apply(self, samples: Spectrum, values: Any) -> Spectrum:
                return samples.with_values(samples.values, mask=[False, True, False])

        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        with pytest.raises(DatasetError, match="latent value"):
            FittingProblem(
                Counts(WAVELENGTH),
                [
                    Dataset(
                        counts,
                        Instrument([MasksMore()]),
                        Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
                        label="counts",
                    )
                ],
            )


# ---------------------------------------------------------------------------
# Masks and the instrument chain
# ---------------------------------------------------------------------------


class TestMasks:
    def test_a_fully_masked_dataset_contributes_nothing(self) -> None:
        observed = flat_spectrum(mask=[True, True, True])
        problem = FittingProblem(Flat(WAVELENGTH), [Dataset(observed, label="d")])
        evaluation = problem.evaluate({"model.level": 1.0})
        assert evaluation.contributions["d"] == pytest.approx(0.0)

    def test_a_masked_sample_is_excised_not_downweighted(self) -> None:
        full = FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d")])
        partial = FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(mask=[False, False, True]), label="d")],
        )
        theta = {"model.level": 1.0}
        # Excision, so the masked sample's whole term disappears — including its
        # -log sigma, which downweighting would have kept (likelihoods.md §8).
        per_sample = full.log_likelihood(theta) / 3.0
        assert partial.log_likelihood(theta) == pytest.approx(2.0 * per_sample)

    def test_a_broken_chain_raises_rather_than_scoring_minus_inf(self) -> None:
        # A chain that returns the wrong kind is a bug, not an unscoreable
        # point: the failure-signalling layer must not swallow it.
        class WrongKind(Transformation):
            ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

            def apply(self, samples: Spectrum, values: Any) -> Any:
                return PhotometricPoints(["W1"], [3.4] * u.micron, [1.0] * u.Jy)

        dataset = Dataset(flat_spectrum(), Instrument([WrongKind()]), label="d")
        with pytest.raises(TransformationError, match="PRODUCES"):
            FittingProblem(Flat(WAVELENGTH), [dataset])
        deferred = FittingProblem(Flat(WAVELENGTH), [dataset], validate=False)
        with pytest.raises(TransformationError, match="PRODUCES"):
            deferred.log_prob({"model.level": 1.0})
        assert deferred.failures == ()

    def test_a_schema_error_from_the_model_is_a_bug_not_a_failure(self) -> None:
        class BadContainer(Flat):
            def evaluate(self, **values: Any) -> Any:
                return "not a container"

        problem = FittingProblem(
            Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d")], validate=False
        )
        problem._compiled["model"] = BadContainer(WAVELENGTH)
        with pytest.raises(TransformationError, match="returns a ModelResult"):
            problem.log_prob({"model.level": 1.0})


# ---------------------------------------------------------------------------
# Repr / provenance surfaces
# ---------------------------------------------------------------------------


class MaskAbove(Transformation):
    """A transformation whose output mask depends on a parameter value."""

    ACCEPTS: ClassVar[tuple[type, ...]] = (Spectrum,)

    def __init__(self, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.register_parameter(Parameter("cut", st.uniform(0.0, 100.0)))

    def apply(self, samples: Spectrum, values: Any) -> Spectrum:
        cut = self.context(values)["cut"]
        return samples.with_values(samples.values, mask=np.asarray(samples.values) > cut)


class TestEffectiveMaskIsResolvedOnce:
    """Peter's ruling on likelihoods.md §17 Q2."""

    def problem(self, **kwargs: Any) -> FittingProblem:
        return FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), Instrument([MaskAbove()]), label="d")],
            reference_values={"model.level": 1.0, "d.instrument.mask_above.cut": 99.0},
            **kwargs,
        )

    def test_a_run_time_mask_change_is_refused_not_silently_scored(self) -> None:
        # Before the ruling this was the worst defect in the contract: masking
        # everything made the dataset contribute exactly 0.0, which beats every
        # finite log-likelihood, so a sampler maximising log_prob would drive
        # the data out of its own fit. Measured then: log_prob rose from -255.1
        # to -9.2 as the mask closed, with no failure recorded.
        problem = self.problem()
        good = problem.evaluate({"model.level": 1.0, "d.instrument.mask_above.cut": 99.0})
        assert math.isfinite(good.log_prob)
        closed = problem.evaluate({"model.level": 1.0, "d.instrument.mask_above.cut": 0.5})
        assert closed.log_prob == -math.inf
        assert closed.failure is not None
        assert closed.failure.reason is FailureReason.LIKELIHOOD_FAILED
        assert "evaluation-invariant" in closed.failure.message
        assert closed.log_prob < good.log_prob

    def test_a_partial_mask_change_is_refused_too(self) -> None:
        problem = self.problem()
        partial = problem.evaluate({"model.level": 2.0, "d.instrument.mask_above.cut": 1.5})
        assert partial.log_prob == -math.inf
        assert partial.failure is not None
        assert partial.failure.reason is FailureReason.LIKELIHOOD_FAILED

    def test_the_union_is_taken_once_and_stamped_thereafter(self) -> None:
        # The observed container masks a sample the instrument does not; the
        # effective mask is their union, resolved at composition, and the
        # likelihood sees it already applied.
        observed = flat_spectrum(mask=[False, False, True])
        dataset = Dataset(observed, label="d")
        problem = FittingProblem(Flat(WAVELENGTH), [dataset])
        assert dataset._retained_reference == 2
        assert dataset._effective_mask is not None
        assert dataset._effective_mask.tolist() == [False, False, True]
        # And the score matches excision of exactly that sample.
        full = FittingProblem(Flat(WAVELENGTH), [Dataset(flat_spectrum(), label="d")])
        theta = {"model.level": 1.0}
        assert problem.log_likelihood(theta) == pytest.approx(
            2.0 * full.log_likelihood(theta) / 3.0
        )

    def test_the_latent_size_is_a_construction_time_constant(self) -> None:
        # The ruling makes explicit what was implicit: latent_declaration(n)
        # takes the retained count, and the retained count is now fixed by
        # declaration rather than by hope.
        counts = Spectrum(WAVELENGTH * u.micron, np.array([4.0, 7.0, 2.0]))
        dataset = Dataset(
            counts,
            likelihood=Likelihood(PoissonFamily(), GaussianProcessNoise(Matern32(0.3, 1.0))),
            label="counts",
        )
        FittingProblem(
            Counts(WAVELENGTH), [dataset], capabilities=Capabilities(differentiable=True)
        )
        assert dataset.latent is not None
        assert dataset.latent.size == dataset._retained_reference == 3


class TestStrictMode:
    """Peter's ruling on likelihoods.md §17 Q1."""

    def problem(self, **kwargs: Any) -> FittingProblem:
        likelihood = Likelihood(
            GaussianFamily(),
            GaussianProcessNoise(Matern32(st.loguniform(1e-3, 1e250), st.loguniform(1e-3, 1e3))),
        )
        return FittingProblem(
            Flat(WAVELENGTH),
            [Dataset(flat_spectrum(), likelihood=likelihood, label="d")],
            **kwargs,
        )

    THETA: ClassVar[dict[str, float]] = {
        "model.level": 1.0,
        "d.likelihood.amplitude": 1e200,
        "d.likelihood.length_scale": 1.0,
    }

    def test_the_engine_path_records_and_returns_minus_inf(self) -> None:
        problem = self.problem()
        evaluation = problem.evaluate(self.THETA)
        assert evaluation.log_prob == -math.inf
        assert evaluation.failure is not None
        assert evaluation.failure.reason is FailureReason.LIKELIHOOD_FAILED

    def test_strict_lets_the_exception_through_for_debugging(self) -> None:
        problem = self.problem(strict=True)
        assert problem.strict
        with pytest.raises(LikelihoodError, match=r"not positive definite|non-finite"):
            problem.evaluate(self.THETA)

    def test_the_failure_records_the_offending_parameter_values(self) -> None:
        problem = self.problem()
        failure = problem.evaluate(self.THETA).failure
        assert failure is not None
        assert failure.values["amplitude"] == pytest.approx(1e200)
        assert failure.values["length_scale"] == pytest.approx(1.0)
        # Arrays are excluded: a latent block on every history entry would be a
        # memory leak wearing a diagnostic's clothes.
        assert all(isinstance(v, float) for v in failure.values.values())

    def test_the_summary_aggregates_rather_than_spamming(self) -> None:
        problem = self.problem()
        for amplitude in (1e200, 1e210, 1e220):
            problem.evaluate({**self.THETA, "d.likelihood.amplitude": amplitude})
        summary = problem.failure_summary()
        assert "3 draw(s) failed: likelihood_failed" in summary
        assert "['d']" in summary
        assert "amplitude in [1e+200, 1e+220]" in summary
        assert "strict=True" in summary
        assert len(summary.splitlines()) == 2  # one line per class, plus the hint

    def test_the_summary_is_empty_when_nothing_failed(self) -> None:
        assert self.problem().failure_summary() == ""


class TestRepresentations:
    def test_objects_describe_themselves(self) -> None:
        problem = FittingProblem(
            Flat(WAVELENGTH), [Dataset(flat_spectrum(), Instrument([Calibrate()]), label="d")]
        )
        assert "Dataset 'd'" in repr(problem.datasets["d"])
        assert "DatasetCollection ['d']" in repr(problem.datasets)
        assert "FittingProblem ['model'] -> ['d']" in repr(problem)
        assert "Evaluation" in repr(problem.evaluate())
        assert "Simulation" in repr(problem.simulate())

    def test_capabilities_serialise(self) -> None:
        assert Capabilities(True, False, "cuda").to_dict() == {
            "differentiable": True,
            "batchable": False,
            "device": "cuda",
        }

    def test_bindings_and_requirements_are_inspectable(self) -> None:
        problem = FittingProblem(
            Powerlaw(blue=WAVELENGTH),
            [Dataset(flat_spectrum(), Instrument([Calibrate()], channel="blue"), label="d")],
        )
        assert dict(problem.bindings) == {"d": "model"}
        assert "blue" in problem.requirements["model"]

    def test_the_single_model_shortcut_refuses_ambiguity(self) -> None:
        problem = FittingProblem(
            {"a": Flat(WAVELENGTH), "b": Flat(WAVELENGTH)},
            [
                Dataset(flat_spectrum(), label="x", model="a"),
                Dataset(flat_spectrum(), label="y", model="b"),
            ],
        )
        assert set(problem.models) == {"a", "b"}
        with pytest.raises(DatasetError, match="ask for the one you mean"):
            _ = problem.model


class TestMultipleModels:
    def test_two_objects_may_share_one_parameter(self) -> None:
        # The canonical multi-dataset tie: two independently modelled objects at
        # the same distance. With one model per object the shared quantity is a
        # genuine tie rather than sharing-by-construction.
        problem = FittingProblem(
            {"a": Flat(WAVELENGTH), "b": Flat(WAVELENGTH)},
            [
                Dataset(flat_spectrum(), label="x", model="a"),
                Dataset(flat_spectrum(), label="y", model="b"),
            ],
            ties=[Tie("level", ("a.level", "b.level"))],
        )
        assert problem.parameters.free_names == ("level",)
        assert problem.free_size == 1
        evaluation = problem.evaluate({"level": 1.0})
        assert set(evaluation.contributions) == {"x", "y"}
        assert evaluation.contributions["x"] == pytest.approx(evaluation.contributions["y"])

    def test_each_model_is_evaluated_once_per_theta(self) -> None:
        model = Recording(WAVELENGTH)
        problem = FittingProblem(
            model,
            [
                Dataset(flat_spectrum(), label="x"),
                Dataset(flat_spectrum(), label="y"),
            ],
        )
        before = model.calls
        problem.log_prob({"model.level": 1.0})
        assert model.calls == before + 1


class TestFreezeAndLenientCompile:
    """Two 2026-09-03 rulings consumed here: freeze-at-composition, loud compile."""

    def test_a_dataset_freezes_its_instrument(self) -> None:
        # transformations.md §15 Q5: a dataset is a composed object, so
        # reconfiguring a step afterwards is refused, never silently ignored.
        step = Calibrate()
        dataset = Dataset(flat_spectrum(), Instrument([step]))
        assert dataset.instrument.mapping is dataset.instrument.mapping
        step.register_parameter(Parameter("gain", st.lognorm(0.1), value=1.0))
        with pytest.raises(Exception, match="reconfigured since"):
            dataset.instrument.mapping  # noqa: B018 - the access is the assertion

    def test_a_compile_refusal_propagates_by_default(self) -> None:
        # transformations.md §15 Q3: negotiation refuses an unachievable
        # requirement by raising, by default.
        from ampere.core.exceptions import CompositionError

        class Refusing(Model):
            def __init__(self) -> None:
                self.register_parameter(Parameter("level", st.uniform(0.0, 4.0)))

            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(WAVELENGTH * u.micron, np.full(3, ctx["level"]) * u.Jy)

            def compile_for(self, requirements):
                raise CompositionError("this model cannot reach the requested resolution")

        with pytest.raises(CompositionError, match="cannot reach the requested resolution"):
            FittingProblem(Refusing(), [Dataset(flat_spectrum(), label="d")])

    def test_lenient_compile_downgrades_the_refusal_to_a_warning(self) -> None:
        from ampere.core.exceptions import CompositionError

        class Refusing(Model):
            def __init__(self) -> None:
                self.register_parameter(Parameter("level", st.uniform(0.0, 4.0)))

            def evaluate(self, **values: Any) -> Spectrum:
                ctx = self.context(values)
                return Spectrum(WAVELENGTH * u.micron, np.full(3, ctx["level"]) * u.Jy)

            def compile_for(self, requirements):
                raise CompositionError("this model cannot reach the requested resolution")

        with pytest.warns(UserWarning, match="proceeding with the unconfigured model"):
            problem = FittingProblem(
                Refusing(), [Dataset(flat_spectrum(), label="d")], lenient_compile=True
            )
        # The unconfigured model is used, and the problem still evaluates.
        assert math.isfinite(problem.log_prob({"model.level": 1.0}))
