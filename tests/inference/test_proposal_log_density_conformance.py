"""W5.32 (k): one conformance row per producer of ``sample_stats.proposal_log_density``.

``results.md`` §9 (W5.0, amended by W5.14 and W5.29) makes a project-wide
promise: **every engine whose stored draws are not draws from the target**
records the proposal's own log-density per draw, in
``sample_stats.proposal_log_density``, *in the same coordinates the stored*
``log_prior``/``log_likelihood`` *already are* — the constrained
free-parameter vector, not an internal unconstrained one a fitted guide or a
trained density estimator may have worked in. That one requirement is what
makes::

    weight = exp(log_prior + log_likelihood - proposal_log_density)

a valid, self-normalising importance weight computed from the stored groups
alone, on any engine, without knowing which one produced the run.

Three producers exist in the tree today (confirmed by ``grep -rn
'"proposal_log_density"' ampere/``, which finds exactly three call sites that
*write* the key -- a fourth hit, in ``ampere/results/population.py``, only
*reads* it back for reweighting):

* :class:`~ampere.inference.VIEngine` (``_vi.py`` ~570), every guide family;
* :class:`~ampere.inference.BlackjaxEngine`'s pathfinder route (``_blackjax.py``
  ~442);
* :class:`~ampere.inference.SBIEngine`'s NPE draws (``_sbi.py`` ~2055).

W5.14's carried note counted **four** producers and left the fourth
unnamed; this module's own grep (above) finds three writers and one reader,
so either the population reweighter's read was the "fourth", or W5.14's
count was one too many -- there is no fourth *producer* to find in the
current tree.

Two claims per producer, not one, because "finite and shaped right" is not
"correct": all three convert their proposal's density from the
**unconstrained** coordinates it was fitted in into the constrained ones
``results.md`` §9 asks for, by subtracting
:func:`~ampere.inference.engine.unconstrained_jacobian_correction` -- one
shared formula, so a producer skipping the subtraction would silently bias
every reweighting by exactly the change of variables it forgot (and, on an
unconstrained-support parameter, would not: the identity bijection's
jacobian is zero, so the fixture below deliberately gives ``norm`` a
positive-constrained, ``st.lognorm`` prior instead, so this module can
actually catch that bug rather than pass regardless of it). The Laplace
guide's own fitted density is a multivariate Gaussian *by construction*, on
any problem, so it is known in closed form independently of whether it
happens to be a good approximation to the true posterior -- checked
directly against an independently built ``scipy`` density from the
engine's own recorded ``guide_loc``/``guide_scale_tril``, in the spirit of
``tests/inference/test_vi.py``'s own row for the same guide (there checked
on a conjugate problem where the fit is *also* the exact posterior; this
module does not need that stronger property). Blackjax's pathfinder and
SBI's NPE are checked by the importance-weight identity instead, since
neither's proposal is a density this module can write down independently:
the weight it implies must still be finite and positive on every draw.

Each producer needs a different extra (jax for VI's laplace here, `jax` +
`blackjax` for the pathfinder route, `sbi` for NPE), so this is
``tests/inference`` rather than ``tests/conformance``: the three producers
share no fixture-parametrised backend the way the conformance battery's
column-per-backend rows do. A producer absent from this environment skips
by name, naming the extra it needs.

**A note on the item text that dispatched this module**: it describes the
convention as "the log-density of the *unconstrained* draw under the
approximation". ``results.md`` §9's own text (quoted above) and all three
producers' own comments say the opposite -- *constrained* coordinates, with
the unconstrained-space density explicitly moved there by subtracting the
jacobian term. This module tests the contract as ``results.md`` §9 and the
tree's own code actually state it; there is no producer/spec disagreement
to record.
"""

from __future__ import annotations

import importlib.util
import warnings
from typing import Any

import astropy.units as u
import numpy as np
import pytest
import scipy.stats as st

import ampere.backends.reference as ReferenceModule
from ampere.core import Dataset, FittingProblem, GaussianFamily, Likelihood, Spectrum
from ampere.core import IndependentNoise as CoreIndependentNoise
from ampere.inference import VIEngine
from ampere.inference.engine import unconstrained_jacobian_correction

HAS_JAX = importlib.util.find_spec("jax") is not None
HAS_BLACKJAX = HAS_JAX and importlib.util.find_spec("blackjax") is not None
HAS_TORCH = importlib.util.find_spec("torch") is not None
HAS_SBI = HAS_TORCH and importlib.util.find_spec("sbi") is not None

needs_jax = pytest.mark.skipif(not HAS_JAX, reason="needs the 'jax' extra (pixi run -e jax ...)")
needs_blackjax = pytest.mark.skipif(
    not HAS_BLACKJAX,
    reason="needs the 'blackjax' extra, folded into the 'jax' pixi environment",
)
needs_sbi = pytest.mark.skipif(not HAS_SBI, reason="needs the 'sbi' extra (pixi run -e sbi ...)")

#: A small power-law fit, repeated here rather than imported from
#: ``test_vi.py``/``test_blackjax.py``'s own near-identical fixtures:
#: ``tests/`` carries no ``__init__.py`` (``tests/conftest.py``'s own
#: docstring), so a sibling test module is not import-safe across files
#: without the dynamic-loader dance ``tests/inference/test_sbi.py``'s
#: ``_sibling`` uses for ``examples/`` -- overkill for one small,
#: self-contained fixture.
#:
#: ``norm``'s prior is deliberately **not** ``st.norm`` (unconstrained
#: support, whose constraining bijection is the identity and whose jacobian
#: correction is therefore always exactly zero -- a producer that forgot
#: the subtraction entirely would still pass every row below). ``st.lognorm``
#: constrains ``norm`` to the positive reals through a real log bijection, so
#: the jacobian term below is genuinely non-zero and a dropped or
#: sign-flipped correction is something this module can actually catch.
SEED = 20260927
REFERENCE_WAVELENGTH = 1.0
GRID = np.geomspace(1.0, 10.0, 20)
SIGMA = 0.2
INDEX = -1.0
PRIOR_SHAPE, PRIOR_SCALE = 0.3, 2.0


def _power_law(grid: np.ndarray, norm: float, index: float) -> np.ndarray:
    return norm * (grid / REFERENCE_WAVELENGTH) ** index


def _noisy(grid: np.ndarray, truth: np.ndarray, sigma: float, seed: int) -> Spectrum:
    rng = np.random.default_rng(seed)
    return Spectrum(
        grid * u.um,
        (truth + rng.normal(0.0, sigma, grid.size)) * u.Jy,
        uncertainty=np.full(grid.size, sigma) * u.Jy,
    )


DATA = _noisy(GRID, _power_law(GRID, 2.0, INDEX), SIGMA, seed=7)


def _problem(module: Any, *, seed: int = SEED) -> FittingProblem:
    """One free parameter (``norm``, positive-constrained); ``index`` fixed.

    ``getattr(module, "IndependentNoise", None)``: the reference backend
    (used for the SBI row) does not re-export its own -- ``BACKEND="reference"``
    is already :class:`ampere.core.IndependentNoise`'s own default, so there
    is nothing to re-export -- unlike torch/jax, which each declare a
    same-shaped subclass for the capability flags alone.
    """
    noise_cls = getattr(module, "IndependentNoise", None) or CoreIndependentNoise
    return FittingProblem(
        module.PowerLaw(
            GRID,
            norm=st.lognorm(PRIOR_SHAPE, scale=PRIOR_SCALE),
            index=INDEX,
            reference_wavelength=REFERENCE_WAVELENGTH,
        ),
        [Dataset(DATA, likelihood=Likelihood(GaussianFamily(), noise_cls()))],
        seed=seed,
    )


def _assert_shape_and_finiteness(run: Any) -> np.ndarray:
    """One value per draw (``chain`` x ``draw``), finite -- the basics every producer owes."""
    stats = run["sample_stats"].dataset
    proposal = np.asarray(stats["proposal_log_density"])
    lp = np.asarray(stats["lp"])
    assert proposal.shape == lp.shape
    assert np.all(np.isfinite(proposal))
    return proposal.ravel()


def _importance_weights(run: Any) -> np.ndarray:
    """``exp(log_prior + log_likelihood - proposal_log_density)``, one per draw."""
    stats = run["sample_stats"].dataset
    log_prior = np.asarray(stats["log_prior"]).ravel()
    log_likelihood = np.asarray(stats["log_likelihood"]).ravel()
    proposal = np.asarray(stats["proposal_log_density"]).ravel()
    return np.exp(log_prior + log_likelihood - proposal)


@needs_jax
def test_vi_laplace_matches_the_closed_form_density_of_its_own_fit() -> None:
    """The one approximation whose density is known in closed form.

    A Laplace guide's fit is a multivariate Gaussian by construction, and
    its own recorded ``guide_loc``/``guide_scale_tril`` describe it
    completely regardless of how good the fit is -- ``proposal_log_density``
    must equal that Gaussian's density, moved into constrained coordinates,
    on every draw. ``norm``'s positive-constrained (``st.lognorm``) prior
    makes the jacobian term genuinely non-zero, so this is a real check of
    the coordinate move and not a coincidence of an identity bijection.
    """
    import ampere.backends.jax as backend

    backend.configure_x64()
    problem = _problem(backend)
    engine = VIEngine(problem)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = engine.run(draws=500, steps=1500, guide="laplace")
    proposal = _assert_shape_and_finiteness(run)

    constrained = np.asarray(run["posterior"]["model.norm"]).ravel()
    unconstrained_2d = np.array(
        [problem.parameters.unconstrain({"model.norm": value}) for value in constrained]
    )  # (draws, 1) -- one free parameter, and unconstrained_jacobian_correction wants rows
    unconstrained = unconstrained_2d.ravel()
    jacobian = unconstrained_jacobian_correction(problem, unconstrained_2d)
    loc = float(np.asarray(engine.guide_loc).ravel()[0])
    scale = float(np.asarray(engine.guide_scale_tril).reshape(1, 1)[0, 0])
    expected_unconstrained_density = st.norm(loc, scale).logpdf(unconstrained)
    # proposal is already constrained-space; adding the jacobian back
    # recovers the guide's own unconstrained-space density (engine.py's
    # docstring: "Subtracting this term ... moves it into the same
    # coordinates" -- this is that subtraction, inverted, to compare against
    # scipy in the space the guide was actually fitted in).
    np.testing.assert_allclose(proposal + jacobian, expected_unconstrained_density, atol=1e-6)


@needs_blackjax
def test_blackjax_pathfinder_satisfies_the_importance_weight_identity() -> None:
    """No closed form here (Pathfinder's Gaussian is an L-BFGS by-product,
    not exact), so checked by the identity ``results.md`` §9 exists for."""
    import ampere.backends.jax as backend

    from ampere.inference import BlackjaxEngine

    backend.configure_x64()
    problem = _problem(backend)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = BlackjaxEngine(problem, method="pathfinder").run(draws=4000)
    _assert_shape_and_finiteness(run)
    weights = _importance_weights(run)
    assert np.all(np.isfinite(weights))
    assert np.all(weights > 0.0)


@needs_sbi
def test_sbi_npe_satisfies_the_importance_weight_identity() -> None:
    """SBI's trained density estimator, same identity -- on the reference
    backend, since NPE fits any problem regardless of array library."""
    from ampere.inference import SBIEngine

    problem = _problem(ReferenceModule)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        run = SBIEngine(problem, method="npe", budget=300).run(
            draws=500, training={"max_num_epochs": 30}
        )
    _assert_shape_and_finiteness(run)
    weights = _importance_weights(run)
    assert np.all(np.isfinite(weights))
    assert np.all(weights > 0.0)
