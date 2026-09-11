Simulation-based inference
===========================

Some models have no likelihood you can write down: a compiled
radiative-transfer code, a Monte Carlo photon-transport simulator, anything
where you can generate synthetic data at a given :math:`\theta` but cannot
evaluate :math:`p(\text{data} \mid \theta)`. :class:`~ampere.inference.SBIEngine`
is ampere's answer — neural simulation-based inference over `sbi
<https://sbi-dev.github.io/sbi/>`_ 0.27, fitting the same
:class:`~ampere.core.FittingProblem` every other engine fits, on the same
reference backend, by training a neural network on simulated
:math:`(\theta, x)` pairs instead of consuming ``log_prob``.

Six scripts in ``examples/sbi/`` cover it end to end, and this page walks
through them in the order a user would actually reach for them: fitting a
black-box simulator (plus the bare native-problem case beside it), caching
the trained posterior, truncated marginal ratio estimation for tighter
posteriors, checking that a posterior is calibrated, and the embedding and
reproducibility questions that cut across all of them. Every one of them is
also exercised by the test suite, named at the point it is introduced below,
so none of this page can rot unnoticed.

Requires the ``sbi`` extra (``pixi install -e sbi``; ``pip install -e
".[sbi]"``), which resolves to **sbi 0.27.0 with a CPU torch 2.13 and no
pyro-ppl** — so it unlocks :class:`~ampere.inference.SBIEngine` without
unlocking :class:`~ampere.inference.NUTSEngine`/:class:`~ampere.inference.VIEngine`
on a torch problem. ``sbi`` and ``torch`` are imported lazily inside
``SBIEngine``'s methods, so ``import ampere.inference`` in the base install
pulls in neither.

A black-box simulator, end to end
----------------------------------

The case the whole layer exists for is a compiled simulator behind a Python
call, and ``examples/sbi/`` builds one honestly rather than faking it with a
Python function. :download:`toy_powerlaw.py
<../../examples/sbi/toy_powerlaw.py>` is not ampere code at all — it imports
nothing from ampere — and behaves the way the Fortran, C++ and Rust
radiative-transfer codes ampere is asked to wrap actually behave: parameters
arrive as command-line arguments and an input file, output goes to a named
file, progress goes to stdout and diagnostics to stderr, and a parameter
combination it cannot handle (here, ``norm`` above 9.5) is a non-zero exit
status with a message on stderr, not an exception. Its physics is a power
law, :math:`f(x) = \text{norm} \cdot x^{\text{index}}`, chosen only so the
wrapper around it can be checked against a closed form.

:download:`external_simulator.py <../../examples/sbi/external_simulator.py>`
is that wrapper: ``ExternalPowerlaw``, a black-box
:class:`~ampere.core.Model` whose ``evaluate`` marshals the grid to a file,
builds a command line, runs the program in a scratch directory of its own,
and reads the answer back off disk. A non-zero exit becomes the declared
``SimulatorFailed`` exception, which
:class:`~ampere.core.dataset.FittingProblem`'s ``simulator_failures=`` turns
into a flagged, counted failure rather than a fatal one. Run it directly to
see the mechanics without a fit::

    $ python examples/sbi/external_simulator.py 24 2
    24 draw(s) simulated, 22 usable (2 failed)
         2 x model_failed
      first failure: model_failed [model]: toy_powerlaw.py exited 2 at index=-1.40318, norm=9.53144; stderr: toy_powerlaw: FATAL: failed to converge at norm=9.53144 (the opacity table does not extend this far); stdout: toy_powerlaw: index=-1.40318 norm=9.53144
      observations stacked as (24, 12)
      theta stacked as (24, 2)

That budget ran under a pool of subprocesses,
:meth:`~ampere.core.dataset.FittingProblem.simulate_many`'s
:class:`~ampere.core.simulate.ProcessExecutor`, which **defaults to the
``forkserver`` start method** on POSIX (`fork` remains available through
``mp_context=`` for a caller who wants it; a plain `fork` can deadlock a
threaded runtime such as torch's or jax's, which is why it is not the
default). Failure accounting is the parent's — the counts above are right
whichever worker produced the failing draw, because ``simulate_many`` owns
both ends of the pool and replays the records in draw order — and a
per-simulation ``timeout=`` on the executor kills a draw that overruns and
flags it ``execution_failed``, a different count from a simulator that
reported its own failure. This script stops at the ``(θ, x)`` pairs, on
purpose: it is not an SBI fit, only what one is built on.
``tests/examples/test_external_simulator.py`` runs it under the pool,
including a simulated crash and a simulated timeout.

:download:`fit_external_simulator.py
<../../examples/sbi/fit_external_simulator.py>` is the fit.  It composes the
*same* black-box model and hands the problem straight to
:class:`~ampere.inference.SBIEngine` — ``executor=`` is the only line that
mentions parallelism, and the engine never sees a worker::

    $ pixi run -e sbi python examples/sbi/fit_external_simulator.py 40 2
    npe on reference: 40 simulation(s), 40 usable (0 failed)
      network: maf, embedding none, 229 epoch(s), summary layout flat (12 feature(s))
      prior in unconstrained coordinates; context none
      posterior:
        model.index    -1.1817 +- 0.6794   [-1.9856, -0.4061]
        model.norm     +2.1317 +- 0.5568   [+1.5584, +2.7035]
      true log p(theta, x) over the draws: -207.361 +- 235.414
      estimator log-density (normalised): +0.071 +- 1.142
      54.8 s wall clock

(a real run at that budget, quoted verbatim; a production fit wants a larger
one — the module's own default is 400 simulations over 4 workers).

What :meth:`~ampere.inference.SBIEngine.run` returns is an ordinary ArviZ
run, like every other engine's, and that is the point worth dwelling on:

* the ``posterior`` group is drawn back in the **constrained**, user-facing
  coordinates, even though the density estimator was trained in
  **unconstrained** ones (``ampere_sbi_parameterisation``) — a bounded prior
  needs no restriction and the estimator sees :math:`\mathbb{R}^n`;
* ``sample_stats`` carries **two log-densities per draw**: ``lp``,
  ``log_prior`` and ``log_likelihood`` are the *true* ones, scored on the
  numpy contract path after the fit (which is why they can be compared with
  the estimator's own opinion, above — a comparison only possible because
  this toy's likelihood happens to be writable, which is exactly why it is a
  good demonstration and a bad advertisement); ``ampere_sbi_log_prob`` beside
  them is the estimator's own, and ``ampere_sbi_log_prob_kind`` says whether
  it is normalised (NPE) or known only up to the evidence (NLE, NRE);
* a rejected simulation never reaches the network: failed draws are dropped
  and counted (``ampere_sbi_failures``, ``ampere_sbi_usable_simulations``),
  which is ``inference.md`` §13's reject-and-record signal arriving where a
  user can read it;
* the run's full provenance is an ``ampere_sbi_*`` attribute for every
  setting that could otherwise be silently wrong: the method
  (``ampere_sbi_method``), the trainer (``ampere_sbi_trainer``), the density
  estimator architecture and its training loss curve
  (``ampere_sbi_density_estimator``, ``ampere_sbi_training_loss``,
  ``ampere_sbi_final_training_loss``), the summary layout and its width
  (``ampere_sbi_summary_layout``, ``ampere_sbi_summary_features``), the
  executor and chunk size (``ampere_sbi_executor``,
  ``ampere_sbi_chunk_size``), and the ``sbi``/torch versions
  (``ampere_sbi_version``).

Pass ``--training-set pairs.nc`` to write every simulated pair to a netCDF
training set (``results.md`` §11), failures included, a chunk at a time, so a
budget larger than memory reaches disk without being held in it.
``tests/inference/test_sbi.py``'s ``TestTheShippedExample`` runs this script
at a small budget on every gate.

**The black-box case is not the only one.** When the model is built from a
*native* backend instead — :mod:`ampere.backends.torch` or
:mod:`ampere.backends.jax` — nothing about the fit changes, and that is worth
seeing on its own: :download:`npe_native.py
<../../examples/sbi/npe_native.py>` is the bare case, an ordinary native
``FittingProblem`` handed straight to ``SBIEngine`` — no subprocess, no
``cache=``, no truncation, one call. ``simulate_many``'s native batched path
(W3.1 slice 2) runs the whole training budget automatically — ``native=None``,
"use it if it works" — without ``SBIEngine`` ever asking for it by name.
``tests/examples/test_npe_native.py`` runs it at a small budget on every
gate.

Caching a trained posterior
-----------------------------

Training a network is the expensive part of an SBI fit, and
``DEVELOPMENT_PLAN.md`` §7 names the trap directly: an SBI posterior reused
against a problem it was not trained on is silently wrong. ``SBIEngine``
closes it with a ``cache=`` argument —
:class:`~ampere.results.ArtefactStore`, keyed by
:func:`~ampere.results.artefact_key` on hashes the engine already
computes — so the *same* problem, method, architecture, budget and encoding
train once and every identical rerun loads instead.

:download:`cached_fit.py <../../examples/sbi/cached_fit.py>` demonstrates the
store directly, one level below where ``SBIEngine(cache=...)`` calls it —
useful in its own right for a caller who wants the caching seam without a
full engine, and the reason this script is written the way it is rather than
around ``SBIEngine`` (a historical accident of dispatch order, noted in its
own docstring: it predates W3.3's encoding layout, which is why it uses
``layout="flat"`` explicitly and sbi's own prior rather than the engine's
unconstrained bridge). A first run trains and writes; a second, identical
run is a hit; changing a setting is a miss that *names what moved* rather
than leaving a caller to notice on their own::

    $ pixi run -e sbi python examples/sbi/cached_fit.py --budget 150 --epochs 15
    cache miss: trained a new posterior (6.55s)
      nothing was cached here before
      cache directory: ~/.cache/ampere-examples/sbi-cached-fit
      five posterior draws of slope: [2.517, 3.48, 3.03, 1.542, 1.603]

    $ pixi run -e sbi python examples/sbi/cached_fit.py --budget 150 --epochs 15
    cache hit: trained nothing (2.71s to load and verify)
      ...

    $ pixi run -e sbi python examples/sbi/cached_fit.py --budget 400 --epochs 15
    cache miss: trained a new posterior (4.16s)
      ingredient(s) that changed since the last write here: budget (150 -> 400)
      ...

(three real, consecutive runs; the cache directory sits under the platform
cache directory, never inside the repository — ground rule 7: no binary
artefacts in git). Nothing about the store looks inside what it caches — it
treats the trained ``sbi`` posterior as opaque bytes, verified on load by
actually sampling from it, not merely trusted to have deserialised.

An engine-level cache hit or miss is visible the same way on an ordinary
``SBIEngine(problem, cache=store).run(...)`` fit, as ``ampere_sbi_cache_hit``
and ``ampere_sbi_cache_key`` in the run's attrs — ``cached_fit.py`` is worth
reading for the mechanism, not because it is how a caller normally reaches
for it.

Truncated marginal ratio estimation (TMNRE)
----------------------------------------------

A single amortised fit trains a network that can answer *any* observation
from the same instrument, which is expensive when only one observation is on
the table. ``method="tmnre"`` (Miller et al. 2021) trades that amortisation
for concentration: each round trains a likelihood-to-evidence ratio
classifier — one per 1-D parameter (``marginals=1``) or per unordered pair
(``marginals=2``) rather than one for the whole joint, so no estimator has to
represent correlations it does not need to — restricts the *prior* to the
box where those marginals exceed ``truncation_epsilon`` (default ``1e-4``,
roughly ±4σ on the worked example, and boxes stop shrinking noticeably after
round 2) times their own maximum, and repeats. Expressed through ``sbi``
0.27's own ``NRE`` trainers and a ``RestrictedPrior`` subclass — **not** a
revived swyft implementation: swyft pins ``pytorch-lightning<=1.9.5``, which
cannot install beside the torch the ``sbi`` extra resolves to, and Peter
ruled against reviving it once that was known (2026-09-09).

:download:`tmnre_fit.py <../../examples/sbi/tmnre_fit.py>` prints every claim
the method makes rather than asserting them in a docstring::

    $ pixi run -e sbi python examples/sbi/tmnre_fit.py --budget 150 --rounds 2 --draws 30 --sample-with mcmc
    tmnre on reference: 2 round(s) x 150 simulations, 2 marginal estimator(s), 21.5s
      amortised: no -- a truncated estimator is chosen at this observation and is not amortised across others
      truncation box per round (constrained coordinates; truth {'model.norm': 3.0, 'model.index': -0.8}):
        round 1: model.norm [0.818, 10.525], model.index [-2.254, 0.745]  (log-volume +2.036, prior mass accepted 100.0%, contains the truth: yes)
        round 2: model.norm [2.048, 3.957], model.index [-1.293, -0.543]  (log-volume -0.705, prior mass accepted 98.8%, contains the truth: yes)
      marginals (estimated 1-D posteriors over the final box):
        model.norm   90% [   2.802,    3.237]  ...sparkline...
        model.index  90% [  -0.938,   -0.758]  ...sparkline...
      posterior (joint estimator, sampled by mcmc):
        model.norm      2.870 +/- 0.500  (truth 3.000, -0.26 sigma)
        model.index    -0.968 +/- 0.177  (truth -0.800, -0.95 sigma)

(a real run at a reduced budget; the sparklines above are elided, the script
draws them as one-line Unicode bar charts). Three things a run states about
itself, because they are easy to get wrong silently:

* **the boxes nest, and the script says whether each still holds the
  truth** — a box that lost it would be the method failing in the one way
  that matters, since no later round can put back mass an earlier round threw
  away;
* **the ``posterior`` group is an ordinary posterior**, not the marginals
  read back as one — a joint ratio estimator trained across every round
  (``discard_prior_samples=True`` would save none of the budget, since it was
  measured to cost about a third of it for no gain) multiplied by the *final*
  truncated prior, sampled by MCMC over the same trained estimator by default
  or, with ``--sample-with rejection``, by rejection (i.i.d., what a single
  "chain" promises) — the default was switched on 2026-09-10 because
  rejection's cost rises sharply as truncation succeeds: 353.6 s against
  43.8 s for the same three rounds at the module's default budget, and *both*
  numbers get worse as the method works, not better;
* **``ampere_sbi_amortised`` is ``0`` for every TMNRE run**, whatever its
  round count — the box is chosen at this observation, so the trained
  estimator must not be re-conditioned on another one, which is the price of
  the concentration the boxes show. ``method="nre"`` at one round stays
  amortised and is the right choice when the same network must serve many
  observations.

The run also carries a fifth, non-default group beside the usual four:
``marginals``, one estimated 1-D (or pair) posterior per parameter over the
final box, and ``ampere_sbi_truncation`` — one JSON record per round with the
box, its log-volume and the proposal's acceptance rate, which is what the
script's table above is printing.
``tests/inference/test_sbi.py``'s ``TestTheShippedTMNREExample`` runs this
script at a small budget on every gate.

Calibration: is the posterior honest?
----------------------------------------

None of the above says whether the *trained network* is telling the truth —
a badly calibrated posterior can look reasonable on one observation and be
systematically too narrow, too wide, or biased on the population of
observations the instrument actually sees. That is
``diagnostics.md`` §11's fourth diagnostic family, and for an ``SBIEngine``
run it is nearly free:
:meth:`~ampere.inference.SBIEngine.calibrate` re-conditions the run's own
trained posterior on a fresh simulation batch, which for an amortised
estimator costs one batch and **no retraining at all**. The arithmetic is
``sbi``'s own ``run_sbc``/``check_sbc`` and ``run_tarp``/``check_tarp``; what
the method owns is handing them the right tensors, encoded with **the run's
own** :class:`~ampere.core.encoding.EncodingLayout` rather than one rebuilt
from the calibration batch — otherwise a check could pass by testing a
different, re-standardised network from the one the run actually produced.

A real, small run (the toy linear model above, ``method="npe"`` at budget
300)::

    >>> run = engine.run(draws=50, training={"max_num_epochs": 40})
    >>> calibration = engine.calibrate(count=80, posterior_draws=200, tarp=True)
    >>> sorted(calibration.data_vars)
    ['c2st_ranks', 'coverage', 'ks_pvalue', 'ranks', 'tarp_coverage']
    >>> float(calibration["ks_pvalue"])
    0.136

Under a calibrated posterior the rank of the truth among the posterior draws
is uniform, so :func:`~ampere.results.plot_sbc_ranks` draws a flat histogram
with the 99 % binomial null band and names the failure shape on the panel
(**U-shaped**: too narrow; **inverted-U**: too wide; a **slope**: biased) and
:func:`~ampere.results.plot_coverage` plots the empirical-against-nominal
credible-level curve, TARP's own curve alongside it in black where it was
requested. Both take the run or the bare ``calibration`` group.
:func:`~ampere.results.attach_calibration` writes the group into a run so it
round-trips through netCDF like every other group.

**Two qualifications the status rows record, stated plainly rather than
buried.** First, a TMNRE run's posterior is defined *on its truncation box*,
so ``calibrate()`` checks it against the *truncated* prior
(``ampere_calibration_reference``) — checking against the original prior
would report miscalibration for every truth the box legitimately excludes,
which is an artefact of the method rather than a property of the estimator.
Second, an SBI run emits **one chain of i.i.d. draws**, not the
multi-chain structure a NUTS or VI run's ``sample_stats`` carries — the
usual multi-chain diagnostics (:math:`\hat R`, ESS) do not apply, and
calibration is the tool that stands in their place for this engine.
``examples/wstat_comparison.py --coverage`` is the general (non-SBI) route's
worked example, over :func:`~ampere.results.sbc`, at real scale — it is
where the profiled Cash-with-background statistic's spectral index is shown
failing simulation-based calibration (Kolmogorov–Smirnov :math:`p =
4\times10^{-4}`) where the equivalent two-dataset Bayesian fit, on the same
simulations, does not (:math:`p = 0.18`); see :doc:`wstat_comparison`.

Embedding choices
--------------------

Every ``SBIEngine`` fit needs a fixed-size summary of the observed data to
condition the network on, and the *default* is a fixed-size vector: each
dataset's observed values, masked samples dropped, concatenated in
``datasets`` order (``layout="flat"``, ``ampere_sbi_summary_layout``). That
is enough for a single fitting problem, where the data layout never changes
between simulation and inference — which is why it stays the default rather
than something a user must opt out of.

For a problem whose datasets vary in size or sampling between simulations —
irregular sampling, missing data, amortising across differently-configured
instruments — the coordinate–value–mask **encoding**
(:class:`~ampere.core.encoding.EncodingLayout`, ``docs/design/contracts/
encoding.md``) packs every dataset into one tensor: standardised coordinates,
Fourier coordinate features, whitened values, log σ and one mask column per
sample, frozen and hashed from the **observed data alone**, before any
simulation. Two embeddings read it (``layout="set"``, ``embedding="set"`` or
``"transformer"``): a permutation-invariant pooling network and an attention
one, both `sbi` 0.27 nets wrapped so the one mask column becomes whatever the
underlying net expects (a NaN-filled row for the set embedding, an
``attention_mask`` for the transformer — worked around here because sbi 0.27
drops it unless the transformer is built causal). Both pool every retained
token — a mask-weighted mean, not `sbi`'s own last-token read, which W3.11
found silently wrong under padding — into a vector of width
``max(2 · free_size, 32)`` by default (a floor against an under-parameterised
first layer emitting all zeros on a one- or two-parameter problem;
``output_dim=`` always overrides it).

That width and readout are **defaults, not findings** (Peter's rider,
2026-09-09): how the right embedding width and pooling strategy depend on
the data, the model and the problem's structure is deferred to a later
embedding study, alongside a matching study of TMNRE's truncation ``ε``
(``DEVELOPMENT_PLAN.md`` §6).

Reproducibility
-------------------

A seeded :class:`~ampere.core.FittingProblem` makes every other engine's run
repeat bitwise. Until W3.15, an ``SBIEngine`` run did not: nothing seeded
torch's own global generator, so a network's initial weights, its trainer's
batch order, and — for an MCMC-sampled NLE, NRE or TMNRE posterior — the
sampler's own momentum draws varied between two runs of the *same* seeded
problem. ``run()`` now seeds torch (and, where an MCMC posterior draws
through it, numpy's legacy global generator) from the problem's own
sub-stream before every network it builds, before each round's training, and
before the final posterior draw — cache hit or not — restoring both
generators once the seeded step finishes, so nothing leaks into unrelated
code sharing the process. The seed actually used is recorded as
``ampere_sbi_torch_seed`` (absent when ``problem.seed is None``, which asks
for fresh randomness every time)::

    >>> run.attrs["ampere_sbi_torch_seed"]
    1282168454

``tests/inference/test_sbi.py``'s ``TestTorchIsSeededFromTheProblem`` asserts
the bitwise repeat, for both an NPE fit and a TMNRE one sampled by MCMC.

.. seealso::

   :doc:`ampere.inference` covers the full API surface, including
   :class:`~ampere.inference.SBIEngine`'s constructor arguments and every
   ``ampere_sbi_*`` attribute; :doc:`ampere.results` covers
   :class:`~ampere.results.ArtefactStore` and the calibration group in
   detail; :doc:`advanced` links back here from "very slow models" and
   embedding networks; :doc:`migrating` maps the legacy ``SBI_SNPE`` class
   onto :class:`~ampere.inference.SBIEngine`.
