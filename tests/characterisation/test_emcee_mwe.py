"""Characterisation test for the emcee flow (``examples/minimal_working_example.py``).

Wraps ``ampere.infer.emceesearch.EmceeSearch`` end-to-end: build the same
linear-model + photometry + spectrum problem as the example, run a short
seeded MCMC, and check the run completes, produces the expected output
shapes, diagnostic plots and a picklable result, and recovers posterior
summary statistics close to previously-recorded golden values.

Deliberate deviation from the legacy example's defaults: we pass
``preopt=False`` and construct the walker guess directly (tightly around
the values used to generate the synthetic data) instead of relying on
``EmceeSearch.optimise``'s default ``preopt=True`` path. That path selects
its scipy.optimize starting point via an unseeded
``numpy.random.default_rng()`` call
(``ampere/infer/emceesearch.py:234-258``), which we found makes the
resulting MAP estimate -- and hence the whole chain -- non-reproducible
between runs (see the W0.4 report's out-of-scope findings). Skipping it
also removes an expensive numerical-gradient ``scipy.optimize.minimize``
call, keeping this test fast.
"""

from __future__ import annotations

import pickle

import numpy as np

from ampere.infer.emceesearch import EmceeSearch
from emcee import moves

from conftest import (
    SEED,
    assert_means_within_tolerance,
    build_linear_sed_problem,
    seed_default_rng,
    xfail_if_pyphot_incompatible,
)

NWALKERS = 32
NSAMPLES = 300
BURNIN = 150

# Golden posterior means, captured from a real deterministic run of this
# exact test (fixed seed, preopt disabled -- see module docstring) against
# a pyphot version compatible with legacy ampere.data.photometry. Order:
# [slope, intercept, calVar, cov scale factor, cov scale length].
#
# Tolerances are deliberately loose: slope/intercept are only weakly (and
# anti-correlated-ly) constrained by a ~300-sample chain, and the two noise
# nuisance parameters are only weakly constrained by construction (the
# synthetic data has i.i.d., not correlated, noise). The tolerances are
# chosen to comfortably survive numpy/scipy/emcee point-version drift while
# still catching a badly broken likelihood (e.g. slope/calVar off by an
# order of magnitude, or NaNs).
GOLDEN_MEANS = [0.97514, 1.37913, 0.99975, 0.06655, 0.00735]
TOLERANCES = [0.2, 1.5, 0.05, 0.2, 0.05]
PARAM_LABELS = ["slope", "intercept", "calVar", "cov scale factor", "cov scale length"]


@xfail_if_pyphot_incompatible
def test_emcee_minimal_working_example(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    seed_default_rng(monkeypatch, SEED)

    problem = build_linear_sed_problem(SEED)
    model = problem["model"]
    dataset = problem["dataset"]

    move_schedule = [(moves.DEMove(), 0.8), (moves.DESnookerMove(), 0.2)]
    optimizer = EmceeSearch(
        model=model, data=dataset, nwalkers=NWALKERS, moves=move_schedule, vectorize=False
    )
    # emcee keeps its own isolated RandomState, independent of numpy's
    # global legacy RNG (see the ``random_state`` property on
    # emcee.EnsembleSampler); seed it explicitly for reproducibility.
    optimizer.sampler.random_state = np.random.RandomState(SEED + 1).get_state()

    rng = np.random.default_rng(SEED + 2)
    base_guess = np.array(
        [problem["true_slope"], problem["true_intercept"], 1.0, 0.05, 0.05]
    )
    guess = [
        base_guess + 1e-3 * rng.standard_normal(optimizer.npars)
        for _ in range(optimizer.nwalkers)
    ]

    optimizer.optimise(
        nsamples=NSAMPLES, burnin=BURNIN, guess=guess, preopt=False, progress=False
    )

    # --- output shapes ---
    assert optimizer.allSamples.shape == (NWALKERS, NSAMPLES, optimizer.npars)
    assert optimizer.samples.shape == (NWALKERS * (NSAMPLES - BURNIN), optimizer.npars)
    assert np.all(np.isfinite(optimizer.samples))

    # --- posterior summary statistics vs golden values ---
    means = optimizer.samples.mean(axis=0)
    assert_means_within_tolerance(means, GOLDEN_MEANS, TOLERANCES, PARAM_LABELS)

    # --- postProcess: diagnostic plots ---
    optimizer.postProcess()
    produced = {p.name for p in tmp_path.glob("*.png")}
    for suffix in ("_walkers.png", "_corner.png", "_lnprob.png", "_posteriorpredictive.png"):
        assert any(name.endswith(suffix) for name in produced), (
            f"expected a plot ending in {suffix!r}, got {produced}"
        )
    # plot_covmats only plots Spectrum entries (dataset[1] here, index 1)
    assert any(name.endswith("_covMat_1.png") for name in produced), produced

    # --- pickle round-trip, written to tmp_path only ---
    pickle_path = tmp_path / "test_pickle_emcee.pkl"
    with open(pickle_path, "wb") as f:
        pickle.dump(optimizer, f)
    assert pickle_path.exists()
    with open(pickle_path, "rb") as f:
        reloaded = pickle.load(f)
    assert isinstance(reloaded, EmceeSearch)
    np.testing.assert_array_equal(reloaded.samples, optimizer.samples)
