"""Characterisation test for the SBI (SNPE) flow
(``examples/minimal_working_example_sbi.py``).

Wraps ``ampere.infer.sbi.SBI_SNPE`` end-to-end: build the same linear-model
+ photometry + spectrum problem as the example, run a short seeded
simulation-based-inference round, and check the run completes, produces the
expected output shapes and a picklable result, and recovers posterior
summary statistics close to previously-recorded golden values.
``postProcess()`` (diagnostic plots) is checked in a separate, ``xfail``-ed
test -- see below.

Skipped cleanly (via ``skipif``) when ``torch``/``sbi`` are not installed,
e.g. in an environment built without the ``ampere[sbi]`` extra.

Notes on determinism and budget, and two legacy bugs found while writing
this test (documented in full in the W0.4 report):

* ``ampere/infer/sbi.py`` calls ``numpy.random.default_rng()`` with no seed
  in several places (``SBI_SNPE.sample``,
  ``SBI_SNPE._check_prior_normalisation``); we neutralise this with the
  same ``seed_default_rng`` monkeypatch used by the other characterisation
  tests. Neural-network training itself (``sbi``'s ``SNPE.train()``) can
  still introduce a little residual run-to-run floating-point jitter beyond
  what seeding ``torch``/``numpy`` controls (observed empirically); the
  tolerances below absorb that.
* ``SBI_SNPE(..., check_prior_normalisation=False)`` -- offered specifically
  to *skip* the (expensive) prior-normalisation Monte-Carlo integration --
  is actually broken in this version of ampere: skipping it also skips
  setting ``self._prior_is_normalised``, which ``SBI_SNPE.log_prob`` (called
  by ``sbi`` while validating the prior, immediately inside ``__init__``)
  unconditionally reads, raising ``AttributeError`` before the object is
  even constructed. We therefore leave ``check_prior_normalisation=True``
  (the example's default) and instead pass a much smaller
  ``n_prior_norm_samples`` to keep that step's Monte-Carlo integration
  cheap, which is a supported, documented parameter of ``SBI_SNPE``.
* ``optimizer.postProcess()`` -- called by the example script exactly as
  we call it below -- reliably raises (``IndexError`` or ``ValueError``,
  depending on the run) with the currently-resolved ``sbi`` package
  (0.27.0). See ``test_sbi_postprocess_plots`` for the root cause and
  reference.
"""

from __future__ import annotations

import pickle

import numpy as np
import pytest

from conftest import (
    SEED,
    assert_means_within_tolerance,
    build_linear_sed_problem,
    requires_sbi,
    seed_default_rng,
    xfail_if_pyphot_incompatible,
)

N_PRIOR_NORM_SAMPLES = 2000
NSAMPLES = 400
NSAMPLES_POST = 1000
N_ROUNDS = 1

# Golden posterior means, captured from real runs of this test. Order:
# [slope, intercept, calVar, cov scale factor, cov scale length].
# With only 400 simulations and a single SNPE round, slope/intercept
# (strongly anti-correlated, wide flat prior) are only weakly constrained,
# and NN-training jitter (see module docstring) adds further run-to-run
# spread on top of that -- hence the very wide intercept tolerance, close
# to its full prior width of 20. calVar is tightly constrained by
# construction (calUnc=0.0025) and is checked much more strictly.
GOLDEN_MEANS = [1.031, 0.0, 1.000, 0.192, 0.0087]
TOLERANCES = [0.5, 7.0, 0.05, 0.3, 0.03]
PARAM_LABELS = ["slope", "intercept", "calVar", "cov scale factor", "cov scale length"]


def _build_and_run_sbi(monkeypatch):
    seed_default_rng(monkeypatch, SEED)

    import torch

    from ampere.infer.sbi import SBI_SNPE

    torch.manual_seed(SEED + 6)

    problem = build_linear_sed_problem(SEED)
    optimizer = SBI_SNPE(
        model=problem["model"],
        data=problem["dataset"],
        n_prior_norm_samples=N_PRIOR_NORM_SAMPLES,
    )
    optimizer.optimise(nsamples=NSAMPLES, nsamples_post=NSAMPLES_POST, n_rounds=N_ROUNDS)
    return optimizer


@requires_sbi
@xfail_if_pyphot_incompatible
def test_sbi_minimal_working_example(tmp_path, monkeypatch):
    from ampere.infer.sbi import SBI_SNPE

    monkeypatch.chdir(tmp_path)
    optimizer = _build_and_run_sbi(monkeypatch)

    # --- output shapes ---
    samples = optimizer.samples
    samples_np = (
        samples.detach().numpy() if hasattr(samples, "detach") else np.asarray(samples)
    )
    assert samples_np.shape == (NSAMPLES_POST, optimizer.npars)
    assert np.all(np.isfinite(samples_np))

    # --- posterior summary statistics vs golden values ---
    means = samples_np.mean(axis=0)
    assert_means_within_tolerance(means, GOLDEN_MEANS, TOLERANCES, PARAM_LABELS)

    # --- pickle round-trip, written to tmp_path only ---
    pickle_path = tmp_path / "test_pickle_sbi.pkl"
    with open(pickle_path, "wb") as f:
        pickle.dump(optimizer, f)
    assert pickle_path.exists()
    with open(pickle_path, "rb") as f:
        reloaded = pickle.load(f)
    assert isinstance(reloaded, SBI_SNPE)


@requires_sbi
@xfail_if_pyphot_incompatible
@pytest.mark.xfail(
    reason=(
        "ampere/infer/mixins.py:542 (SBIPostProcessor.print_summary, called "
        "from SBI_SNPE.postProcess) indexes self.bestPars as if it were a "
        "flat (npars,) vector: 'self.bestPars[i]' for i in range(self.npars). "
        "self.bestPars comes from SBIPostProcessor.get_map() "
        "(ampere/infer/mixins.py:292), i.e. 'self.posterior.map()'. With the "
        "currently-resolved sbi package (0.27.0), DirectPosterior.map() "
        "returns a batched tensor whose shape varies from call to call "
        "(observed both (1,) and (1, npars) across otherwise-identical "
        "runs -- an sbi-side quirk, not something we control), so this "
        "raises either 'IndexError: index 1 is out of bounds for dimension "
        "0 with size 1' or 'ValueError: only one element tensors can be "
        "converted to Python scalars' (from the '%.5f' formatting), "
        "depending on the run. This reproduces with "
        "examples/minimal_working_example_sbi.py's own call to "
        "optimizer.postProcess() -- not something introduced by this test. "
        "Not fixed here (frozen legacy code) -- see the W0.4 report. "
        "(No 'raises=' pinned above, deliberately, since the exact "
        "exception type is not stable.)"
    ),
    strict=False,
)
def test_sbi_postprocess_plots(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    optimizer = _build_and_run_sbi(monkeypatch)

    optimizer.postProcess()
    produced = {p.name for p in tmp_path.glob("*.png")}
    for suffix in ("_corner.png", "_posteriorpredictive.png"):
        assert any(name.endswith(suffix) for name in produced), (
            f"expected a plot ending in {suffix!r}, got {produced}"
        )
