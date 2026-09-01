"""Characterisation test for the dynesty flow
(``examples/minimal_working_example_dynesty.py``).

Wraps ``ampere.infer.dynestysearch.DynestyNestedSampler`` end-to-end: build
the same linear-model + photometry + spectrum problem as the example, run a
short nested-sampling run with a fixed ``rstate``, and check the run
completes, produces the expected output shapes, diagnostic plots and a
picklable result, and recovers posterior summary statistics close to
previously-recorded golden values.

Unlike emcee/zeus, ``DynestyNestedSampler`` forwards an explicit ``rstate``
straight through to ``dynesty.NestedSampler`` (it isn't consumed by
anything else in ``ampere``), so this flow is exactly, deterministically
reproducible given a fixed ``rstate`` -- confirmed empirically by running it
twice and diffing the posterior mean/std bit-for-bit.
"""

from __future__ import annotations

import pickle

import numpy as np

from ampere.infer.dynestysearch import DynestyNestedSampler

from conftest import (
    SEED,
    assert_means_within_tolerance,
    build_linear_sed_problem,
    xfail_if_pyphot_incompatible,
)

NLIVE = 50
DLOGZ = 5.0

# Golden posterior means, captured from a real deterministic run of this
# exact test. Order: [slope, intercept, calVar, cov scale factor, cov scale
# length]. Tolerances are deliberately loose -- see test_emcee_mwe.py's
# module docstring for the rationale, which applies equally here (same
# synthetic dataset, comparable chain size).
GOLDEN_MEANS = [0.97550, 1.33972, 0.99828, 0.06908, 0.00632]
TOLERANCES = [0.2, 1.5, 0.05, 0.2, 0.05]
PARAM_LABELS = ["slope", "intercept", "calVar", "cov scale factor", "cov scale length"]


@xfail_if_pyphot_incompatible
def test_dynesty_minimal_working_example(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    problem = build_linear_sed_problem(SEED)
    model = problem["model"]
    dataset = problem["dataset"]

    rstate = np.random.default_rng(SEED + 3)
    optimizer = DynestyNestedSampler(model=model, data=dataset, nlive=NLIVE, rstate=rstate)
    optimizer.optimise(dlogz=DLOGZ, print_progress=False)

    # --- output shapes ---
    samples = optimizer.results.samples
    assert samples.ndim == 2
    assert samples.shape[0] > NLIVE
    assert samples.shape[1] == optimizer.npars
    assert np.all(np.isfinite(samples))
    assert np.isfinite(optimizer.results.logz[-1])

    # --- posterior summary statistics vs golden values ---
    from dynesty import utils as dyfunc

    weights = np.exp(optimizer.results.logwt - optimizer.results.logz[-1])
    mean, _cov = dyfunc.mean_and_cov(samples, weights)
    assert_means_within_tolerance(mean, GOLDEN_MEANS, TOLERANCES, PARAM_LABELS)

    # --- postProcess: diagnostic plots ---
    optimizer.postProcess()
    produced = {p.name for p in tmp_path.glob("*.png")}
    for suffix in ("_corner.png", "_trace.png", "_posteriorpredictive.png"):
        assert any(name.endswith(suffix) for name in produced), (
            f"expected a plot ending in {suffix!r}, got {produced}"
        )
    assert any(name.endswith("_covMat_1.png") for name in produced), produced

    # --- pickle round-trip, written to tmp_path only ---
    pickle_path = tmp_path / "test_pickle_dynesty.pkl"
    with open(pickle_path, "wb") as f:
        pickle.dump(optimizer, f)
    assert pickle_path.exists()
    with open(pickle_path, "rb") as f:
        reloaded = pickle.load(f)
    assert isinstance(reloaded, DynestyNestedSampler)
    np.testing.assert_array_equal(reloaded.results.samples, optimizer.results.samples)
