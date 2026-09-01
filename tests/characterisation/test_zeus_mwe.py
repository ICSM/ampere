"""Characterisation test for the zeus flow
(``examples/minimal_working_example_zeus.py``).

Wraps ``ampere.infer.zeussearch.ZeusSearch`` end-to-end: build the same
linear-model + photometry + spectrum problem as the example, run a short
seeded ensemble-slice-sampling run, and check the run completes, produces
the expected output shapes, diagnostic plots and a picklable result, and
recovers posterior summary statistics close to previously-recorded golden
values.

Skipped cleanly (via ``skipif``) when ``zeus-mcmc`` is not installed, e.g.
in an environment built without the ``ampere[zeus]`` extra.

Deliberate deviations from the legacy example's defaults, and a known
reproducibility limitation:

* Like the emcee test, we pass ``preopt=False`` and construct the walker
  guess directly -- see ``test_emcee_mwe.py``'s module docstring for why
  (an unseeded ``numpy.random.default_rng()`` call in
  ``ampere/infer/zeussearch.py:99-116``). This is additionally required
  for zeus specifically because zeus's ``EnsembleSampler`` *hard-fails*
  (``ValueError: Invalid walker initial positions!``) if that unseeded
  preopt path produces a poorly-converged MAP estimate that leaves any
  walker with non-finite log-probability -- we hit this in practice while
  developing this test.
* Unlike emcee (which keeps its own seedable ``RandomState``) and dynesty
  (which takes an explicit ``rstate``), zeus's ``EnsembleSampler`` draws
  directly from numpy's global legacy RNG (``np.random.uniform`` etc. in
  ``zeus/ensemble.py`` and ``zeus/moves.py``). We seed that global RNG and
  also pin ``sampler.tune``/``sampler.mu`` (zeus's adaptive step-size
  tuning otherwise makes even the *number* of likelihood evaluations
  non-deterministic). Despite this, repeated runs still produce visibly
  different posterior summaries (confirmed empirically: same seeds, same
  wall-clock-independent code path, yet different results across
  invocations) -- the root cause was not identified within W0.4's scope
  (flagged as an out-of-scope finding in the W0.4 report). Tolerances
  below are therefore wider than the other three flows' to absorb this
  residual run-to-run jitter, while still catching a badly broken
  likelihood.
"""

from __future__ import annotations

import pickle

import numpy as np

from conftest import (
    SEED,
    assert_means_within_tolerance,
    build_linear_sed_problem,
    requires_zeus,
    xfail_if_pyphot_incompatible,
)

NWALKERS = 16
NSAMPLES = 60
BURNIN = 30

# Golden posterior means, captured from real runs of this test (see module
# docstring: zeus's own exploration is not fully reproducible even with
# global-RNG seeding, so these are representative central values with
# generous tolerances rather than an exact target).
# Order: [slope, intercept, calVar, cov scale factor, cov scale length].
GOLDEN_MEANS = [0.978, 1.34, 0.9995, 0.070, 0.0072]
TOLERANCES = [0.25, 2.5, 0.05, 0.3, 0.05]
PARAM_LABELS = ["slope", "intercept", "calVar", "cov scale factor", "cov scale length"]


@requires_zeus
@xfail_if_pyphot_incompatible
def test_zeus_minimal_working_example(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    from ampere.infer.zeussearch import ZeusSearch

    problem = build_linear_sed_problem(SEED)
    model = problem["model"]
    dataset = problem["dataset"]

    np.random.seed(SEED + 4)  # zeus draws directly from the global legacy RNG
    optimizer = ZeusSearch(model=model, data=dataset, nwalkers=NWALKERS, vectorize=False)
    optimizer.sampler.tune = False
    optimizer.sampler.mu = 1.0

    rng = np.random.default_rng(SEED + 5)
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
    assert optimizer.allSamples.shape == (NSAMPLES, NWALKERS, optimizer.npars)
    assert optimizer.samples.shape == (
        (NSAMPLES - BURNIN) * NWALKERS,
        optimizer.npars,
    )
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
    assert any(name.endswith("_covMat_1.png") for name in produced), produced

    # --- pickle round-trip, written to tmp_path only ---
    pickle_path = tmp_path / "test_pickle_zeus.pkl"
    with open(pickle_path, "wb") as f:
        pickle.dump(optimizer, f)
    assert pickle_path.exists()
    with open(pickle_path, "rb") as f:
        reloaded = pickle.load(f)
    assert isinstance(reloaded, ZeusSearch)
    np.testing.assert_array_equal(reloaded.samples, optimizer.samples)
