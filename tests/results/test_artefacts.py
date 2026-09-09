"""W3.5: trained-artefact caching keyed on the problem's own hashes, end to end.

``DEVELOPMENT_PLAN.md`` §7's trap, restated by the item: an SBI posterior,
embedding net or emulator reused against a problem it was not trained on is
silently wrong. This is the test that a cache built from
``ampere.results.provenance``'s hashes actually closes that trap rather than
leaving the exact hole the module docstring found in the item's literal text
(``spec_hashes``' ``"spec"`` entry alone does not move when only the
likelihood family, noise model, solver or kernel change) — so
:class:`~ampere.results.ArtefactKey`'s three hashes are tested **separately**:
a change to the prior, the model or the data must move exactly the ingredient
that changed and no other, which is what lets a miss "say which ingredient
changed" (the item text) rather than only that something did.

``TestRealSBIPosteriorRoundTrip`` is the one class that needs the ``sbi``
extra: a real, trained ``sbi`` 0.27 ``NeuralPosterior`` is pickled through
:class:`~ampere.results.ArtefactStore` and unpickled back, and its
``potential()`` and ``sample()`` are checked to agree with the original,
seeded, exactly — the round trip sbi's own package was verified (live,
against the installed 0.27.0) to support before this module was written to
rely on it. Everything else here needs neither ``sbi`` nor torch and runs in
``dev``, which is also what
``test_importing_artefacts_pulls_in_neither_torch_nor_sbi`` checks directly.

W3.3 (the coordinate-value-mask encoding) had not merged when this module was
written, so ``SBIEngine(cache=...)`` does not exist yet; every test below
exercises :func:`~ampere.results.artefact_key` and
:class:`~ampere.results.ArtefactStore` directly, which is exactly the seam
``SBIEngine`` will call once it does.

Everything is written into pytest's ``tmp_path``: no cache file belongs in
this repository (``AGENTS.md`` ground rule 7).
"""

from __future__ import annotations

import importlib.util
import inspect
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import Dataset, FittingProblem, Model, Parameter, Spectrum
from ampere.core.exceptions import ResultsError
from ampere.results import ArtefactCacheWarning, ArtefactKey, ArtefactStore, artefact_key

HAS_SBI = importlib.util.find_spec("sbi") is not None

#: One decorator rather than a skipif on every row, mirroring
#: ``tests/inference/test_sbi.py``'s ``needs_sbi`` for the same reason: the
#: module still *collects*, and everything not marked with it, in ``dev``.
needs_sbi = pytest.mark.skipif(
    not HAS_SBI,
    reason="needs the 'sbi' extra (pixi run -e sbi ...)",
)

SEED = 20260909
GRID = np.array([1.0, 2.0, 3.0])


class Line(Model):
    """A one-parameter linear model with an adjustable, non-parameter buffer.

    ``scale`` is a plain buffer, not a :class:`~ampere.core.Parameter` — the
    same kind of thing ``model_fingerprint``'s buffer hashing and ``describe``
    hook exist for (``results.md`` §13.13). Two ``Line`` instances built with
    the same ``slope_mean`` but a different ``scale`` therefore have an
    *identical* parameter declaration (so :func:`~ampere.results.spec_hashes`'
    ``"spec"`` entry agrees) but compute something different — precisely the
    "model changed, parameter spec did not" case
    ``ampere.results.artefacts``'s module docstring names as the reason
    ``ArtefactKey.model_hash`` exists beside ``prior_hash``.
    """

    def __init__(
        self, wavelength: np.ndarray, *, slope_mean: float = 2.0, scale: float = 1.0
    ) -> None:
        self.register_buffer("wavelength", wavelength, unit=u.um)
        self.register_buffer("scale", np.asarray(scale, dtype=float))
        self.register_parameter(Parameter("slope", st.norm(slope_mean, 1.0)))

    def evaluate(self, **values: Any) -> Spectrum:
        ctx = self.context(values)
        return Spectrum(
            ctx["wavelength"] * u.um, ctx["scale"] * ctx["slope"] * ctx["wavelength"] * u.Jy
        )


def _problem(
    *,
    slope_mean: float = 2.0,
    scale: float = 1.0,
    values: tuple[float, ...] = (2.0, 4.0, 6.0),
    seed: int | None = SEED,
) -> FittingProblem:
    """A tiny, cheap-to-hash toy problem with three independently movable knobs."""
    observed = Spectrum(GRID * u.um, list(values) * u.Jy, uncertainty=[0.3, 0.3, 0.3] * u.Jy)
    model = Line(GRID, slope_mean=slope_mean, scale=scale)
    return FittingProblem(model, [Dataset(observed)], seed=seed)


def _key(problem: FittingProblem, **overrides: Any) -> ArtefactKey:
    """:func:`artefact_key` with sensible defaults, overridable one at a time."""
    settings: dict[str, Any] = {
        "layout": "flat",
        "method": "npe",
        "architecture": "maf",
        "budget": 200,
        "rounds": 1,
    }
    settings.update(overrides)
    return artefact_key(problem, **settings)


# ---------------------------------------------------------------------------
# The key: every ingredient moves the digest, and only the ingredient that
# actually changed moves its own named hash.
# ---------------------------------------------------------------------------


class TestArtefactKey:
    def test_identical_problem_and_settings_give_an_identical_key(self) -> None:
        first = _key(_problem())
        second = _key(_problem())
        assert first.digest() == second.digest()
        assert first.ingredients() == second.ingredients()

    def test_prior_change_moves_prior_hash_and_only_prior_hash(self) -> None:
        base = _key(_problem())
        changed = _key(_problem(slope_mean=3.0))
        assert base.digest() != changed.digest()
        assert base.prior_hash != changed.prior_hash
        assert base.model_hash == changed.model_hash
        assert base.data_hash == changed.data_hash

    def test_model_change_moves_model_hash_and_only_model_hash(self) -> None:
        """The gap ``spec_hashes`` alone would have left open (module docstring)."""
        base = _key(_problem())
        changed = _key(_problem(scale=2.0))
        assert base.digest() != changed.digest()
        assert base.prior_hash == changed.prior_hash
        assert base.model_hash != changed.model_hash
        assert base.data_hash == changed.data_hash

    def test_data_change_moves_data_hash_and_only_data_hash(self) -> None:
        base = _key(_problem())
        changed = _key(_problem(values=(2.0, 4.0, 6.1)))
        assert base.digest() != changed.digest()
        assert base.prior_hash == changed.prior_hash
        assert base.model_hash == changed.model_hash
        assert base.data_hash != changed.data_hash

    @pytest.mark.parametrize(
        "overrides",
        [
            {"layout": "coordinate_value_mask"},
            {"method": "nle"},
            {"architecture": "nsf"},
            {"budget": 500},
            {"rounds": 2},
        ],
        ids=["layout", "method", "architecture", "budget", "rounds"],
    )
    def test_each_estimator_setting_moves_the_digest(self, overrides: dict[str, Any]) -> None:
        base = _key(_problem())
        changed = _key(_problem(), **overrides)
        assert base.digest() != changed.digest()

    def test_seed_change_moves_the_digest(self) -> None:
        base = _key(_problem(seed=1))
        changed = _key(_problem(seed=2))
        assert base.digest() != changed.digest()
        assert base.seed == 1
        assert changed.seed == 2

    def test_unseeded_problem_records_the_entropy_sentinel(self) -> None:
        key = _key(_problem(seed=None))
        assert key.seed == "entropy"

    def test_versions_ingredient_reflects_which_packages_are_installed(self) -> None:
        key = _key(_problem())
        assert "python" in key.versions
        if HAS_SBI:
            assert "sbi" in key.versions
            assert "torch" in key.versions
        else:
            assert "sbi" not in key.versions
            assert "torch" not in key.versions

    @needs_sbi
    def test_recording_fewer_packages_moves_the_digest(self) -> None:
        """The ``sbi``/``torch`` versions are as load-bearing as anything else.

        Needs the extra installed: in ``dev``, neither package is present, so
        ``package_versions`` records the same (empty) thing whichever names
        are asked for and the two keys would collide — a fact about ``dev``
        having nothing installed to record, not about this ingredient.
        """
        both = _key(_problem())
        one = artefact_key(
            _problem(),
            layout="flat",
            method="npe",
            architecture="maf",
            budget=200,
            rounds=1,
            packages=("sbi",),
        )
        assert both.digest() != one.digest()

    @pytest.mark.parametrize(
        "overrides",
        [
            {"layout": ""},
            {"method": ""},
            {"architecture": ""},
            {"budget": 0},
            {"rounds": 0},
        ],
        ids=["layout", "method", "architecture", "budget", "rounds"],
    )
    def test_a_partial_key_is_refused(self, overrides: dict[str, Any]) -> None:
        with pytest.raises(ResultsError):
            _key(_problem(), **overrides)

    def test_ingredients_are_json_safe(self, tmp_path: Path) -> None:
        """What :meth:`ArtefactStore.put` writes as the sidecar must round-trip through JSON."""
        key = _key(_problem())
        path = tmp_path / "ingredients.json"
        path.write_text(json.dumps(key.ingredients()))
        assert json.loads(path.read_text()) == key.ingredients()


# ---------------------------------------------------------------------------
# The store: hit, miss, a named diff on a miss, and a corrupted entry that
# never raises.
# ---------------------------------------------------------------------------


class TestArtefactStore:
    def test_get_is_none_before_anything_is_stored(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        assert store.get(_key(_problem())) is None

    def test_put_then_get_round_trips_a_plain_object(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        store.put(key, {"posterior": "stub", "loss": [1.0, 0.5]})
        assert store.get(key) == {"posterior": "stub", "loss": [1.0, 0.5]}

    def test_a_different_key_is_a_miss_even_in_the_same_store(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        store.put(_key(_problem()), "trained-for-budget-200")
        assert store.get(_key(_problem(), budget=500)) is None

    def test_train_or_load_trains_exactly_once_and_hits_thereafter(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        calls: list[int] = []

        def train() -> dict[str, int]:
            calls.append(1)
            return {"trained_at_call": len(calls)}

        first, hit_1 = store.train_or_load(key, train)
        second, hit_2 = store.train_or_load(key, train)

        assert hit_1 is False
        assert hit_2 is True
        assert len(calls) == 1
        assert first == second == {"trained_at_call": 1}

    def test_second_run_with_identical_problem_and_settings_trains_nothing(
        self, tmp_path: Path
    ) -> None:
        """The item's own Accept line, at the ``train_or_load`` seam directly."""
        store = ArtefactStore(tmp_path)
        first_problem = _problem()
        second_problem = _problem()  # a *second*, independently-built, identical problem
        calls: list[int] = []

        def train() -> str:
            calls.append(1)
            return "the one true posterior"

        posterior_1, hit_1 = store.train_or_load(_key(first_problem), train)
        posterior_2, hit_2 = store.train_or_load(_key(second_problem), train)

        assert hit_1 is False
        assert hit_2 is True
        assert len(calls) == 1
        assert posterior_1 == posterior_2

    def test_changed_ingredient_is_a_miss_whose_diff_names_it(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        base_key = _key(_problem())
        store.put(base_key, "artefact-v1")

        changed_key = _key(_problem(), budget=500)
        assert store.get(changed_key) is None

        diff = store.diff(changed_key)
        assert diff["budget"] == (200, 500)
        assert "method" not in diff  # unchanged ingredients are not reported

    def test_changed_model_diff_names_model_hash_not_prior_or_data(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        store.put(_key(_problem()), "artefact-v1")

        changed_key = _key(_problem(scale=2.0))
        diff = store.diff(changed_key)

        assert "model_hash" in diff
        assert "prior_hash" not in diff
        assert "data_hash" not in diff

    def test_diff_is_empty_with_nothing_stored_yet(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        assert store.diff(_key(_problem())) == {}

    def test_corrupted_artefact_is_a_silent_miss_with_a_warning(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        store.put(key, "a perfectly good artefact")
        (tmp_path / f"{key.digest()}.pkl").write_bytes(b"not a pickle at all")

        with pytest.warns(ArtefactCacheWarning):
            result = store.get(key)
        assert result is None

    def test_corrupted_artefact_makes_train_or_load_retrain(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        calls: list[int] = []

        def train() -> str:
            calls.append(1)
            return "fresh artefact"

        store.train_or_load(key, train)
        (tmp_path / f"{key.digest()}.pkl").write_bytes(b"garbage")

        with pytest.warns(ArtefactCacheWarning):
            artefact, hit = store.train_or_load(key, train)

        assert hit is False
        assert artefact == "fresh artefact"
        assert len(calls) == 2

    def test_hand_edited_sidecar_is_a_silent_miss_with_a_warning(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        store.put(key, "an artefact")
        sidecar_path = tmp_path / f"{key.digest()}.json"
        payload = json.loads(sidecar_path.read_text())
        payload["ingredients"]["budget"] = 999999
        sidecar_path.write_text(json.dumps(payload))

        with pytest.warns(ArtefactCacheWarning):
            result = store.get(key)
        assert result is None

    def test_unreadable_sidecar_json_is_a_silent_miss_with_a_warning(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        store.put(key, "an artefact")
        (tmp_path / f"{key.digest()}.json").write_text("{not valid json")

        with pytest.warns(ArtefactCacheWarning):
            result = store.get(key)
        assert result is None

    def test_no_force_parameter_exists_anywhere(self) -> None:
        """ "There is no 'force' that skips the spec-hash comparison" (the item text)."""
        for member in (ArtefactStore.get, ArtefactStore.put, ArtefactStore.train_or_load):
            assert "force" not in inspect.signature(member).parameters

    def test_put_writes_only_inside_the_given_root(self, tmp_path: Path) -> None:
        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        store.put(key, "an artefact")
        written = {path.name for path in tmp_path.iterdir()}
        assert f"{key.digest()}.pkl" in written
        assert f"{key.digest()}.json" in written


# ---------------------------------------------------------------------------
# A real, trained sbi posterior: the round trip the item's Accept line asks for.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def trained_posterior() -> Any:
    """A real, small, seeded sbi 0.27 NPE posterior — trained once for this module."""
    import torch
    from sbi.inference import NPE
    from sbi.utils import BoxUniform

    torch.manual_seed(SEED)
    prior = BoxUniform(low=-2.0 * torch.ones(2), high=2.0 * torch.ones(2))
    theta = prior.sample((200,))
    x = theta + 0.1 * torch.randn_like(theta)

    trainer = NPE(prior=prior)
    trainer.append_simulations(theta, x)
    estimator = trainer.train(max_num_epochs=5, show_train_summary=False)
    posterior = trainer.build_posterior(estimator)
    posterior.set_default_x(torch.zeros(2))
    return posterior


@needs_sbi
class TestRealSBIPosteriorRoundTrip:
    """Verified live against the installed sbi 0.27.0 before relying on it:

    a trained ``NeuralPosterior`` pickles and unpickles with plain
    :mod:`pickle` and its ``potential()``/``sample()`` agree with the
    original exactly under a fixed seed — no ``dill``, no state-dict
    reconstruction. That is the whole of why
    :class:`~ampere.results.ArtefactStore` never imports ``sbi`` or
    ``torch``: pickling a posterior is generic, and the posterior's own
    classes carry what they need to reconstruct themselves.
    """

    def test_cached_posterior_round_trips_sample_and_potential(
        self, tmp_path: Path, trained_posterior: Any
    ) -> None:
        import torch

        store = ArtefactStore(tmp_path)
        key = _key(_problem())
        calls: list[int] = []

        def train() -> Any:
            calls.append(1)
            return trained_posterior

        first, hit_1 = store.train_or_load(key, train)
        assert hit_1 is False
        assert first is trained_posterior  # a miss returns train()'s own object, unpickled or not

        second, hit_2 = store.train_or_load(key, train)
        assert hit_2 is True
        assert len(calls) == 1  # the second call never retrained
        assert second is not trained_posterior  # it came back through pickle, not by reference

        torch.manual_seed(123)
        original_draws = trained_posterior.sample((25,), show_progress_bars=False)
        torch.manual_seed(123)
        cached_draws = second.sample((25,), show_progress_bars=False)
        assert torch.equal(original_draws, cached_draws)

        original_potential = trained_posterior.potential(original_draws, track_gradients=False)
        cached_potential = second.potential(original_draws, track_gradients=False)
        assert torch.equal(original_potential, cached_potential)

    def test_a_different_method_setting_is_a_miss_that_names_method(
        self, tmp_path: Path, trained_posterior: Any
    ) -> None:
        store = ArtefactStore(tmp_path)
        store.put(_key(_problem(), method="npe"), trained_posterior)
        miss_key = _key(_problem(), method="nle")
        assert store.get(miss_key) is None
        diff = store.diff(miss_key)
        assert diff["method"] == ("npe", "nle")


# ---------------------------------------------------------------------------
# The namespace stays clean: this module never pulls in torch or sbi.
# ---------------------------------------------------------------------------


def test_importing_artefacts_pulls_in_neither_torch_nor_sbi() -> None:
    """Run as a subprocess: this interpreter may already have imported both.

    Mirrors ``tests/inference/test_sbi.py``'s identical check for
    ``ampere.inference`` — the cache module makes the same promise for the
    same reason: it stores opaque bytes, and never needs the packages that
    produced them to do so.
    """
    probe = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "import sys, ampere.results.artefacts; "
                "print(sorted(n for n in ('torch', 'sbi') if n in sys.modules))"
            ),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert probe.stdout.strip() == "[]"
