"""The jax backend: numpyro distributions, equinox pytrees, gradients.

Rung 2 of ``architecture.md`` §1's capability ladder, beside torch: everything
the reference backend does, plus **differentiability** — which is what unlocks
NUTS, gradient-based VI and gradient-based optimisation — and, in later slices,
batching and GPU.

Import discipline
-----------------
``architecture.md`` §4 rule 2: ``import ampere``, ``import ampere.core``,
``import ampere.inference`` and ``import ampere.results`` must never require
jax, and ``ampere/backends/__init__.py`` imports nothing at all. **This
subpackage is the only place in ampere that imports jax at module top level**,
and importing it is the user's explicit opt-in::

    from ampere.backends import jax as ampere_jax   # needs the `jax` extra

The extra is ``pip install "ampere[jax]"``: jax, numpyro (the distributions,
the ``biject_to`` registry and the NUTS kernel) and equinox (the
``partition``/``combine`` filter-spec mechanism). CPU jaxlib is what that
resolves to and is what CI wants.

float64 is not a default here, it is a policy — and you must turn it on
--------------------------------------------------------------------------
``architecture.md`` §5 requires float64 for all likelihood and GP linear
algebra on every backend, always; a GP solve in float32 fails in ways that read
as science bugs. jax's ``jax_enable_x64`` is process-global and its own
documentation addresses it to the *application* author, never to a library, so
**ampere never flips it for you** (``lowering.md`` §10.2(a), ruled 2026-09-01).
Instead:

>>> from ampere.backends.jax import configure_x64
>>> configure_x64()          # idempotent; call before any jax work

and every model, instrument step, kernel, solver and lowered parameter set here
**raises** at construction if the flag is off, naming the three remedies
(§10.2(c)). Raising rather than warning is deliberate: a warning in a notebook
scrolls away, and the failure it precedes is a plausible-looking wrong answer
rather than a crash.

What is shipped
---------------
Models (``DEVELOPMENT_PLAN.md`` §5's Phase 2 list): :class:`BlackBody`,
:class:`ModifiedBlackBody`, :class:`PowerLaw` — the same declarations as the
reference backend's, computed in ``jax.numpy``.

Instrument steps (``transformations.md`` §10's table):
:class:`CalibrationScale`, :class:`Resample`, :class:`LSFConvolution`,
:class:`SyntheticPhotometry`.

Kernels and solvers: :class:`Matern32`, :class:`SquaredExponential` and
:class:`DenseGP` — the neutral ``ampere.core`` declarations with their linear
algebra done in jax. ``QuasisepGP`` is **not** here yet: slice 2 chooses
between tinygp's ``QuasisepSolver`` and celerite2.jax after evaluating both
against the conformance suite (``DEVELOPMENT_PLAN.md`` §6), and shipping a
solver slot that quietly delegated to numpy would defeat the point.

Noise (``likelihoods.md`` §5): the prediction-aware
:class:`FractionalModelNoise` and its GP composition
:class:`FractionalModelGPNoise` (X-1), plus :class:`IndependentNoise` and
:class:`GaussianProcessNoise` — this backend's declarations of the two core
compositions, same names, same declarations, ``BACKEND = "jax"``. The last
two exist because **W2.13** widened the capability flags to ``NoiseModel``
and ``GPSolver`` (``inference.md`` §10a, fold-in 7), so composing
``ampere.core.IndependentNoise`` into a jax problem is now a backend
disagreement rather than a silent numpy island.

Lowering: :func:`~ampere.backends.jax.parameters.LoweredParameterSet` is
``lowering.md`` §5's declaration-form table made executable — numpyro sample
sites, ``biject_to`` bijections, real ``numpyro.plate``\\ s, and the §3.6
``icdf`` contract.

**Every piece here declares ``BACKEND = "jax"``** beside ``DIFFERENTIABLE``,
``BATCHABLE`` and ``DEVICE``, explicitly rather than by inheriting the ABCs'
defaults. That string is this backend's one name everywhere (W2.12): what
``ampere.core.declared_capabilities`` aggregates onto a problem, what a run's
``ampere_backend`` records, the key ``ampere.core.lowering``'s registry is
consulted with, and the conformance fixture's id.
"""

from __future__ import annotations

from ._config import BACKEND, configure_x64, require_x64, x64_enabled
from .bijections import log_abs_det_jacobian, lower_bijection
from .distributions import lower_prior
from .gp import DenseGP, Matern32, SquaredExponential
from .instrument import CalibrationScale, LSFConvolution, Resample, SyntheticPhotometry
from .models import COORDINATE_UNIT, FLUX_UNIT, BlackBody, ModifiedBlackBody, PowerLaw, planck_jy
from .noise import (
    FractionalModelGPNoise,
    FractionalModelNoise,
    GaussianProcessNoise,
    IndependentNoise,
)
from .parameters import LoweredParameterSet, LoweringFallbackWarning, filter_spec
from .problem import LoweredProblem, lower_problem

# ``inference.md`` §10a: importing this package is the user's opt-in to jax,
# and it is also the moment the jax realisation becomes reachable through
# ``ampere.core.realise`` -- which is how ``ampere.inference``'s gradient-based
# drivers get a differentiable density without importing a backend.
from ampere.core import register_realisation as _register_realisation

_register_realisation(BACKEND, lower_problem, builtin=True)

__all__ = [
    "BACKEND",
    "COORDINATE_UNIT",
    "FLUX_UNIT",
    "BlackBody",
    "CalibrationScale",
    "DenseGP",
    "FractionalModelGPNoise",
    "FractionalModelNoise",
    "GaussianProcessNoise",
    "IndependentNoise",
    "LSFConvolution",
    "LoweredParameterSet",
    "LoweredProblem",
    "LoweringFallbackWarning",
    "Matern32",
    "ModifiedBlackBody",
    "PowerLaw",
    "Resample",
    "SquaredExponential",
    "SyntheticPhotometry",
    "configure_x64",
    "filter_spec",
    "log_abs_det_jacobian",
    "lower_bijection",
    "lower_prior",
    "lower_problem",
    "planck_jy",
    "require_x64",
    "x64_enabled",
]
