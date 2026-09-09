WStat, and why ampere does not ship it
=======================================

X-ray astronomers fitting a source spectrum together with its background
almost always reach for the profiled Cash-with-background statistic — the
"W statistic" in XSPEC's terminology. It needs no background *model* at
all: the background rate in each channel is replaced by its
maximum-likelihood value, profiled out analytically, so the fit only ever
varies the source parameters.

**Ampere does not ship a WStat likelihood family**, and this is a
deliberate, documented decision (ruled 2026-09-03; see
``DEVELOPMENT_PLAN.md`` §2's "WStat / profile likelihoods" row and
``docs/design/modalities/awkward_instrument.md`` §5 and §9 question 3). A
profiled statistic is not a marginal likelihood — the background's own
uncertainty is replaced by a point estimate rather than integrated out — so
it does not produce a Bayesian posterior in the sense the rest of ampere
does. Ampere's documentation takes the opinionated line that the correct
approach is the **two-dataset Bayesian formulation**: source and background
as two datasets, sharing one background model, fit jointly.

That said, WStat is expressible in ampere, because the likelihood contract's
extension surface is exactly "a user adds a family without touching ampere"
(``docs/design/contracts/likelihoods.md`` §3). The worked example below
builds it as a user family, builds the two-dataset joint fit alongside it,
runs both, and compares them — pros, cons, and results, not false balance.

The full, runnable code is :download:`examples/wstat_comparison.py
<../../examples/wstat_comparison.py>`. It is included here as a literal
script rather than an executed notebook deliberately: this repository's
policy is that no run outputs, figures or other binary artefacts are
committed to git, and Sphinx's ``nbsphinx_allow_errors = True`` setting
means a notebook that failed to execute would not fail this documentation
build either — so a script the reader (and
``tests/examples/test_wstat_comparison.py``) can actually run, rather than a
page whose apparent success proves nothing, is the honest choice here. Run
it yourself with::

    python examples/wstat_comparison.py

which prints the comparison below, deterministically (every random draw —
the synthetic data and both engines' own streams — is seeded), in around
ten to twenty seconds on the reference backend. ``--coverage`` adds the
repeated-trial coverage study described at the end of this page.

The two routes
---------------

**Route 1 — the profiled statistic, as a user family**
(``ProfiledCashWithBackground``). One dataset (the source-region counts);
the background-region counts are data the family carries as its own
buffer, aligned per channel, together with the exposure/area ratio between
the two regions. For channel :math:`i`, the profiled background estimate is

.. math::

    \hat b_i = \frac{C_i + \sqrt{C_i^2 + 4\rho(\rho+1) B_i m_i}}{2\rho(\rho+1)},
    \qquad C_i = \rho (S_i + B_i) - (\rho + 1) m_i,

with :math:`\rho` the ratio, :math:`m_i` the model's predicted source
counts, and :math:`S_i`, :math:`B_i` the observed source- and
background-region counts. ``log_prob`` sums
``log Poisson(S_i; m_i + ratio * b_hat_i) + log Poisson(B_i; b_hat_i)``
over channels.

Three things this family gets right, on purpose, because they are exactly
what the likelihood contract asks of a user family carrying its own aligned
data:

* **Safe under masking.** The background buffer is excised with
  ``NoiseParams.retain`` — the boolean indicator over the *full*,
  unmasked containers every ampere noise model sets — so a masked source
  channel cannot silently pair against the wrong background channel.
* **Its ``sample()`` refuses**, and for a specific, structural reason
  rather than only the contract's general one: the profiled background
  estimate is a function of the very counts a draw would need to produce,
  so there is no forward generative model to sample from without already
  having the data it profiles over.
* **Its per-sample log-likelihood terms are flagged as degraded.** Each
  channel's term already used that channel's own counts to fix its own
  nuisance parameter, so — even though the terms sum exactly to
  ``log_prob``, unlike a GP's — they are not valid predictive densities.
  Handing them to ``arviz.loo``/``waic`` would be circular.

**Route 2 — the two-dataset Bayesian joint fit** (``XraySourceAndBackground``,
``build_joint_problem``). Source and background regions as two
``Dataset``\ s in one ``DatasetCollection``, observing two channels of one
model that shares the background parameter between them (`inference.md` §8's
pattern, specialised to the case where the shared physics needs no channel
of its own). Both channels are ordinary ``PoissonFamily`` likelihoods — no
profiling, no user family, nothing bespoke.

The comparison
---------------

Both routes are run with ``emcee`` on the same deterministic, seeded
synthetic data — 32 channels, a faint power-law source over a flat
background, chosen so several channels have very few or zero counts, which
is precisely the regime where a profiled nuisance and a marginalised one can
part company. A representative run's output::

    WStat (profiled, user family) vs the two-dataset Bayesian joint fit
    ======================================================================

    parameter          truth          wstat median [68%]          joint median [68%]
    src_norm           6.000        6.534 [5.727, 7.442]        6.592 [5.666, 7.488]
    src_index         -1.500     -1.701 [-1.911, -1.514]     -1.803 [-2.024, -1.614]

    On the source index, this run's two posteriors agree: each median falls inside
    the other approach's 68% interval; the two intervals are close to the same width.

**Recommendation: prefer the two-dataset Bayesian joint fit.** This holds
regardless of how closely the two routes' numbers happen to agree in any one
run. The joint fit is a proper marginal likelihood: its posterior is a
posterior, its per-sample terms are valid predictive densities for
``arviz.loo``/``waic``, and both regions can be simulated forward for
posterior-predictive checks. WStat can do none of that structurally, not as
a matter of a particular run's luck.

**The trade-off, stated plainly rather than as false balance.** The joint
fit pays for those guarantees with a background *model* — here, one
flat-shape parameter. If that shape assumption is wrong, the joint fit is
biased in a way WStat cannot be, because WStat's per-channel nuisance needs
no shape assumption at all; that is its one genuine advantage, and the
reason it remains a reasonable quick-look tool. The literature on profile
likelihoods in counting experiments (Cash 1979 and its X-ray descendants)
documents that profiling a background out, rather than marginalising it,
tends to bias and overstate the precision of the parameters of interest
specifically in the low-count regime this example uses. A single seeded run
— including the one printed above — cannot itself demonstrate that bias
reliably; only a repeated-trial coverage study can. The next section is that
study.

The repeated-trial coverage study
----------------------------------

Ampere's fourth diagnostic family — posterior calibration
(``docs/design/contracts/diagnostics.md`` §11) — *is* a repeated-trial
coverage study, so running one here costs the example a function call rather
than a bespoke experiment. :func:`ampere.results.sbc` is the Talts et al.
(2018) loop: draw :math:`\theta` from the prior, simulate a dataset there,
fit it from scratch, and record where the truth fell among the posterior
draws. Under a correctly calibrated posterior that **rank** is uniform, so
its histogram is flat and the empirical coverage of every central credible
interval sits on the diagonal::

    python examples/wstat_comparison.py --coverage --full

Both routes are run against the same generative model — the two-dataset
formulation, with source- and background-region counts drawn from
independent Poisson processes at the drawn parameters — and, crucially,
against the **same simulated datasets**. That makes the joint fit a
*control* rather than a second experiment: it fits the very model the data
came from, so its ranks are uniform by construction and its row is what
Monte Carlo noise looks like at this budget.

Two details of how the study is set up matter more than they look.

* **The joint problem needs a ``sample()``.** ``PoissonFamily`` refuses to
  provide one by default — a ``log_prob`` says how a datum is scored, not
  how one is generated, and ampere will not guess an observation process.
  For a counting experiment there is nothing to guess, so the example's
  ``CountingPoisson`` subclass supplies it in three lines. That is the
  supported route the refusal itself names, and it is why the *WStat*
  problem still cannot be simulated from: its profiled background estimate
  is a function of the very counts a draw would have to produce.
* **Simulation-based calibration averages over the prior it draws from.**
  Run under the example's broad ``NORM_PRIOR``
  (:math:`\log\mathcal{U}(0.5, 40)`), most prior draws are *bright*
  sources, where profiling a background out is harmless — and the study
  duly reports both routes as calibrated, with every uniformity p-value
  above 0.19 and coverage curves that agree to within Monte Carlo noise
  (reproduce it with ``--coverage --full --broad-prior``). That is a true
  answer to a question nobody asked. The study therefore runs under
  ``COVERAGE_NORM_PRIOR`` (:math:`\log\mathcal{U}(0.5, 3)`), which is the
  faint regime the concern is actually about. The lesson generalises: a
  calibration study is only as informative as the prior it integrates over,
  and "well calibrated" always carries an implicit "over this prior".

In the faint regime, over 128 simulated datasets with each rank taken
against 200 posterior draws::

    route   parameter              KS p  mean rank  rank var  68% cov.  95% cov.
    wstat   model.norm            0.457      0.518    0.0759     0.736     0.945
    wstat   model.index        0.000423      0.558    0.0722     0.716     0.922
    joint   model.src_norm        0.797      0.496    0.0761     0.719     0.938
    joint   model.src_index       0.183      0.539    0.0731     0.759     0.938

Under uniformity the mean rank fraction is 0.5 with a standard error of
0.026 at 128 trials, and its variance is :math:`1/12 = 0.0833`. The rank
histograms themselves, ten equal bins each, are the picture behind those
numbers (the script prints them; :func:`ampere.results.plot_sbc_ranks` draws
them properly, with the binomial null band, for anyone who wants a figure)::

    wstat   model.norm         8  13  11  13  17  14  13  12  15  12   (uniform expects 12.8)
    wstat   model.index        8  11   6   5  19  17  20  12  18  12   (uniform expects 12.8)
    joint   model.src_norm    11   9  15  12  17  14  17  10  11  12   (uniform expects 12.8)
    joint   model.src_index    9   8  11  15  10  16  17  17  12  13   (uniform expects 12.8)

**The profiled route's spectral index fails the uniformity test; the
control's, on the same simulations, does not.** WStat's ``model.index``
ranks give a Kolmogorov–Smirnov p-value of :math:`4\times10^{-4}` against
0.18 for the joint fit's ``model.src_index``, with a mean rank fraction of
0.558 (2.3 standard errors above 0.5) and a variance below the uniform's.
The histogram says the same thing in one line: 8, 11, 6, 5 in the lower half
against 20, 12, 18, 12 in the upper. Read as a picture, the truth lands in
the upper part of the WStat posterior more often than it should, and that
posterior is a little too tight. Read as a statement about the fit,
profiling the background out has both biased the source index and overstated
its precision — exactly the two failures the literature attributes to a
profiled nuisance at low counts, and neither of them visible in the single
run printed above, whose two posteriors agreed.

Two honest qualifications, because a coverage study that overstates itself
is worse than none. First, the control's own index ranks are shifted in the
same direction by about two thirds as much (0.539, 1.5 standard errors) and
are not significant; the two routes are ranked on the same simulations, so
part of WStat's shift is noise the control shares. What the study
establishes is that the profiled route's departure is *detectable* at this
budget while the correctly specified fit's is not — it does not pin down the
size of the effect to better than 128 trials allow. Second, the source
*normalisation* shows no such failure on either route (p = 0.46 and 0.80):
the bias here is in the spectral shape, not the flux scale.

None of this changes the recommendation, and that is the point of having
made it structurally: the joint fit was the right choice before the study,
for reasons that do not depend on how the ranks fall. What the study adds is
that the cost of the profiled route is now measured rather than asserted.

The reduced version of the study — six simulations, short chains — runs in
``tests/examples/test_wstat_comparison.py``, so the machinery the numbers
above came from is exercised on every commit; the full budget above takes
about twenty minutes on the reference backend, which is why it is behind a
flag rather than in ``main()``.

.. note::

   Rank histograms and coverage curves are drawn by
   :func:`ampere.results.plot_sbc_ranks` and
   :func:`ampere.results.plot_coverage`, and
   ``python examples/wstat_comparison.py --coverage --figures PREFIX``
   writes them. They are not committed to this repository: no run outputs or
   other binary artefacts are (see above).
