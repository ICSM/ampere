"""W5.20: the parameter namespace is core's, and identical on every backend.

``parameters.md`` §10 says which names a parameter or buffer may not take.
Until W5.20 the *enforcement* was ``hasattr(type(self), name)`` over the whole
MRO, so the answer depended on which backend's base class a model inherited —
``grid`` was reserved on a torch spectral model and free on a hand-written one,
``grid_tensor`` was reserved on torch and jax astrometry and nowhere else, and
an interferometric source model could declare a parameter called ``flux`` on
the reference class but not on its twin (W4.3's finding). A core contract's
namespace rule must not be spelled differently per backend, so the rule is now
a stated, finite set.

What these rows hold:

* the **same answer everywhere** — a model declaring parameters called
  ``flux``, ``grid``, ``grid_tensor`` and ``AXIS`` composes *and lowers* on
  every backend, and ``to`` / ``type`` / ``apply`` are refused on every backend
  with the same core message;
* the **pinned torch list** still covers ``torch.nn.Module``'s real namespace,
  so a torch upgrade that adds an attribute fails here rather than in a user's
  lowering;
* torch's lowering **refuses in ampere's vocabulary**, never by letting torch's
  own ``KeyError`` escape;
* **no spec hash moves**: parameter names are unchanged by W5.20, only their
  legality, so a declaration's hash is the one it had before the change.

**Parametrised over the backends installed here**, in
``test_native_interferometry.py``'s shape: the reference backend always runs,
and whichever of torch and jax this environment has joins it.
"""

from __future__ import annotations

import importlib
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

from ampere.core import (
    Dataset,
    FittingProblem,
    GaussianFamily,
    Likelihood,
    Parameter,
    Spectrum,
)
from ampere.core import IndependentNoise as CoreIndependentNoise
from ampere.core.exceptions import LoweringError, ParameterError
from ampere.core.parameter import TORCH_MODULE_NAMES, reserved_names
from ampere.results.provenance import hash_of

#: The names W5.20 frees: each is an attribute of *some* backend's model class
#: (``flux``/``grid`` on the spectral and astrometry models, ``grid_tensor`` on
#: torch's and jax's, ``AXIS`` on both spectral base classes) and none is in the
#: core reserved set, so all four are legal parameter names everywhere.
FREED = ("flux", "grid", "grid_tensor", "AXIS")

#: The names W5.20 reserves *everywhere*: ``torch.nn.Module``'s namespace, which
#: torch's lowering could never have carried. Refused in core, so the answer is
#: the same on the reference and jax paths too.
REFUSED = ("to", "type", "apply")

#: Wavelengths every row here works on, and the observed spectrum built from
#: them. Small and fixed: nothing below is about the numbers.
GRID = np.linspace(1.0, 5.0, 6)
OBSERVED = Spectrum(
    GRID * u.micron,
    np.full(GRID.size, 2.0) * u.Jy,
    uncertainty=np.full(GRID.size, 0.1) * u.Jy,
)

#: ``hash_of(PowerLaw(...).parameters.to_spec())`` for the declaration
#: :meth:`Kit.model` builds, measured on this branch *before* the namespace
#: change landed and unchanged after it. W5.20 moves no parameter name, no
#: prior, no shape and no unit; it moves only which names are *legal*, and
#: method names (``flux`` -> ``native_flux``) are not hashed at all. This row is
#: what says so out loud, since a §4 contract changed
#: (``DEVELOPMENT_PLAN.md`` ground rule 9).
PINNED_SPEC_HASH = "59d52abd29714a756bdf088f10a53a8b"


class Kit:
    """One backend's ``PowerLaw`` and its own uncorrelated noise model.

    A problem must be composed entirely from one backend's pieces
    (``inference.md`` §10a), and since W2.13 the noise model is one of them —
    ``ampere.core``'s ``IndependentNoise`` declares ``"reference"``, so a torch
    problem built with it would be a two-backend problem and refuse to realise.
    The reference backend has no module-level ``IndependentNoise`` of its own;
    core's *is* its one.
    """

    def __init__(self, name: str) -> None:
        self.name = name
        self.module = importlib.import_module(f"ampere.backends.{name}")
        if name == "reference":
            self.noise: Any = CoreIndependentNoise
        else:
            self.noise = self.module.IndependentNoise

    def model(self, extra: tuple[str, ...] = ()) -> Any:
        """This backend's ``PowerLaw``, with *extra* parameters declared on top.

        The shipped declaration is used rather than a bespoke class so that the
        rows below are about the namespace and nothing else: the physics, the
        buffer, the channels and the native surface are whatever this backend
        already ships.
        """
        model = self.module.PowerLaw(GRID, norm=st.lognorm(0.4, scale=2.0), index=st.norm(-1.2, 0.3))
        for name in extra:
            model.register_parameter(Parameter(name, st.norm(0.0, 1.0)))
        return model

    def problem(self, model: Any) -> FittingProblem:
        """*model* against one Gaussian spectrum, in this backend's pieces."""
        return FittingProblem(
            model,
            [Dataset(OBSERVED, likelihood=Likelihood(GaussianFamily(), self.noise()))],
            seed=20260915,
        )


def _installed() -> list[str]:
    """``"reference"``, plus whichever modern backends this environment has."""
    found = ["reference"]
    for name in ("torch", "jax"):
        try:
            module = importlib.import_module(f"ampere.backends.{name}")
        except ImportError:  # the extra is not installed here
            continue
        if name == "jax":
            # ``lowering.md`` §10.2(a): the application turns x64 on, not
            # ampere. A test suite is an application.
            module.configure_x64()
        found.append(name)
    return found


BACKENDS = _installed()


@pytest.fixture(params=BACKENDS)
def backend(request: Any) -> Kit:
    return Kit(request.param)


# ---------------------------------------------------------------------------
# One answer, every backend
# ---------------------------------------------------------------------------


class TestTheSameAnswerOnEveryBackend:
    @pytest.mark.parametrize("name", FREED)
    def test_a_backend_class_attribute_no_longer_blocks_a_parameter(
        self, backend: Kit, name: str
    ) -> None:
        """The heart of the item, one row per freed name per backend.

        Each of these *is* an attribute of at least one shipped model class.
        Under the old rule a model inheriting that class could not declare a
        parameter of the same name while a model inheriting a different one
        could — which made ``parameters.md`` §10 mean three different things.
        """
        model = backend.model((name,))
        assert name in model.parameters
        assert name in model.parameters.free_names

    def test_the_freed_names_all_lower_natively_together(self, backend: Kit) -> None:
        """Not merely declarable: composable, and lowerable where that applies.

        ``flux``, ``grid``, ``grid_tensor`` and ``AXIS`` are declared at once,
        so this also covers the torch case the old rule hid — every parameter
        becomes an attribute of an ``nn.Module`` under lowering, and none of
        these four collides with one.
        """
        problem = backend.problem(backend.model(FREED))
        # Qualified by the merge, and on torch that means a ``model`` child
        # module with a ``flux`` parameter on it — the nesting ``lowering.md``
        # §6.1 specifies, which is where a reserved leaf would have bitten.
        free = set(problem.parameters.free_names)
        assert {f"model.{name}" for name in FREED} <= free
        theta = problem.prior_transform(np.full(problem.free_size, 0.5))
        contract = problem.log_prob(theta)
        assert np.isfinite(contract)
        if backend.name == "reference":
            return
        lowered = importlib.import_module(f"ampere.backends.{backend.name}.problem")
        realised = lowered.lower_problem(problem)
        assert float(realised.log_prob(theta)) == pytest.approx(contract, abs=1e-8)

    @pytest.mark.parametrize("name", REFUSED)
    def test_a_reserved_name_is_refused_with_the_core_message(
        self, backend: Kit, name: str
    ) -> None:
        """Refused by core, so the wording is one wording rather than three.

        ``to`` is a real example rather than a contrived one: it is a method on
        torch's model classes *and* on ``nn.Module``, so before W5.20 it was
        refused on torch at declaration, refused on jax nowhere, and refused on
        torch's lowering by a ``KeyError`` in torch's own vocabulary.
        """
        with pytest.raises(ParameterError, match="is a reserved name"):
            backend.model((name,))

    def test_the_reserved_set_does_not_depend_on_the_backend(self, backend: Kit) -> None:
        """A backend may add public attributes to its models freely.

        Every public attribute this backend's ``PowerLaw`` has that is *not* in
        the core reserved set is a legal parameter name on it — which is the
        amended §10 sentence, checked rather than asserted.
        """
        public = {name for name in dir(type(backend.model())) if not name.startswith("_")}
        assert public - reserved_names(), "this model would make a degenerate check"
        for name in sorted(public - reserved_names()):
            assert name in backend.model((name,)).parameters


# ---------------------------------------------------------------------------
# No spec hash moves
# ---------------------------------------------------------------------------


def test_the_declaration_hash_is_the_one_it_had_before_w5_20() -> None:
    """Ground rule 9's evidence: the contract changed, the hashes did not.

    W5.20 renames *methods* (``flux`` -> ``native_flux``) and re-states which
    *names* are legal. Neither is hashed: ``ParameterSet.to_spec`` records
    names, priors, shapes and units, and a method name appears in none of them.
    """
    assert hash_of(Kit("reference").model().parameters.to_spec()) == PINNED_SPEC_HASH


# ---------------------------------------------------------------------------
# torch: the pinned list, and the lowering guard
# ---------------------------------------------------------------------------


torch_only = pytest.mark.skipif(
    "torch" not in BACKENDS, reason="the pinned nn.Module list is checked where torch is installed"
)


@torch_only
class TestThePinnedTorchNamespace:
    def test_the_literal_still_covers_torch_nn_module(self) -> None:
        """The loud failure a torch upgrade must produce.

        ``ampere.core`` cannot import torch (``architecture.md`` §4 rule 2), so
        the list is a literal. A literal that nobody checks is a literal that
        goes stale, and the failure mode it would produce — a parameter core
        accepted and torch's lowering could not carry — is exactly the
        backend-dependent behaviour this item exists to end.
        """
        nn = importlib.import_module("torch.nn")
        public = {name for name in dir(nn.Module) if not name.startswith("_")}
        assert public <= set(TORCH_MODULE_NAMES), sorted(public - set(TORCH_MODULE_NAMES))

    def test_it_covers_the_instance_attributes_too(self) -> None:
        """``training`` is set in ``__init__``, so ``dir(Module)`` misses it."""
        nn = importlib.import_module("torch.nn")
        instance = {name for name in vars(nn.Module()) if not name.startswith("_")}
        assert instance <= set(TORCH_MODULE_NAMES), sorted(instance - set(TORCH_MODULE_NAMES))


@torch_only
class TestTorchLoweringRefusesByTheSameList:
    """Defence in depth: core refuses these names, and so does the lowering.

    The component half of this has been refused with a ``LoweringError`` since
    W2.4; the *leaf* half fell through to ``nn.Module.register_parameter``,
    whose ``KeyError("attribute 'type' already exists")`` names torch's
    attribute and not ampere's parameter, and is not a ``LoweringError`` the
    conformance battery can match on.
    """

    def test_a_reserved_leaf_is_refused_by_name(self) -> None:
        parameters = importlib.import_module("ampere.backends.torch.parameters")
        root = parameters.LoweredParameters()
        with pytest.raises(LoweringError, match="type"):
            root.check_leaf("type")

    def test_the_refusal_is_amperes_and_not_torchs(self) -> None:
        parameters = importlib.import_module("ampere.backends.torch.parameters")
        root = parameters.LoweredParameters()
        with pytest.raises(LoweringError) as raised:
            root.check_leaf("training")
        assert "core reserved set" in str(raised.value)

    def test_a_free_leaf_is_accepted(self) -> None:
        parameters = importlib.import_module("ampere.backends.torch.parameters")
        root = parameters.LoweredParameters()
        for name in FREED:
            assert root.check_leaf(name) == name

    def test_a_component_name_is_still_refused(self) -> None:
        """The half that already worked, kept honest by the same list."""
        parameters = importlib.import_module("ampere.backends.torch.parameters")
        root = parameters.LoweredParameters()
        with pytest.raises(LoweringError, match="collides with an attribute"):
            root.child("state_dict")
