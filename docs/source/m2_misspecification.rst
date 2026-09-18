Model misspecification, and what the flexible likelihood does about it
=======================================================================

Ampere's distinguishing feature is a likelihood that does not assume your model
is right. This page is the evidence that it works, and the measurement of what
it costs.

The experiment is deliberately small enough to be checked by eye and controlled
enough to be checked by arithmetic. A four-parameter toy — a linear continuum
times two Gaussian absorption lines, on a Gaia-RVS-like band from 0.842 to
0.872 µm — is fitted to data generated from *itself*, and then to the same data
after a deviation the model has no parameter for has been multiplied into it.
Each spectrum is fitted twice: once with the **standard** likelihood
(independent Gaussian noise, the ordinary chi-square) and once with the
**flexible** one (a Matérn-3/2 Gaussian process over the residuals).

The runnable study is ``examples/m2_misspecification``, a small package in the
repository rather than a single script — the generators, the toy model in its
three array libraries, the study driver and the figure code are separate
modules because the tests and the benchmarks import them one at a time::

    python -m examples.m2_misspecification                # 200 points, reference backend
    python -m examples.m2_misspecification --size 20000   # the top of the size ladder
    python -m examples.m2_misspecification --backend jax  # NUTS, if the extra is installed
    python -m examples.m2_misspecification --figures /tmp/m2

Every assertion below is made by ``tests/m2`` (``pixi run test-m2``), at a
short budget, on every pull request; the numbers quoted are from the longer
milestone budget, which ``pytest tests/m2 -m m2_full`` re-asserts. This page
carries no figures, because this repository commits no binary artefacts; the
``--figures`` flag above writes all of them.

The four scenarios
------------------

Each is a multiplicative deviation applied to the noise-free truth before 1 %
Gaussian noise is added. The noise realisation is shared across the four, so
the only thing that changes between them is the deviation.

================= ====================================================
``none``          nothing — the control
``mild``          2.5 % sinusoidal fringing, 0.0028 µm period
``strong_smooth`` 7 % sinusoidal fringing, same period
``strong_sharp``  12 % Gaussian emission line at 0.8630 µm, σ = 0.35 nm
================= ====================================================

The control is not a formality. A GP that widened the posterior whether or not
there was anything to absorb would pass every other test on this page for the
wrong reason, so the first thing to check is that it does not.

Two choices worth stating
-------------------------

**The kernel is Matérn-3/2, not the squared exponential.** Two structural
reasons. A Matérn-3/2 sample path is once differentiable rather than analytic,
which is a better description of a real model deficiency than infinite
smoothness. And it is *exactly* quasiseparable, so
:class:`~ampere.core.QuasisepGP` computes the same likelihood as
:class:`~ampere.core.DenseGP` in O(N) rather than O(N³) — which is the only
reason the 20 000-point rung of the ladder below exists at all. The squared
exponential is kept, dense-only, at 200 points, as the cross-check against what
ampere's own earlier study used; it reaches the same conclusion.

**The "standard" likelihood has no GP in it.** An earlier version of this study
produced its standard case by pinning a GP's hyperparameter priors to
:math:`10^{-10}`. That is not the same object: it keeps two extra dimensions in
the sampler, keeps the GP solve in the hot loop, and its "zero" is a prior whose
mass sits at :math:`10^{-10}` rather than at nothing. The honest statement of
"no flexible likelihood" is :class:`~ampere.core.IndependentNoise`, which is
also what a reader would actually write.

The result
----------

The statistic is ``|median − truth|`` divided by the posterior's **own** 68 %
half-width — how many of its own standard deviations the fit is away from the
right answer. It is the right statistic because it distinguishes the two ways a
misspecified fit can fail. A likelihood that is *wrong* puts the truth many
widths from its median: confident and mistaken. A likelihood that has *absorbed*
the misspecification may still have a displaced median, but its posterior is
wide enough to contain the truth, and its error bar is then something a reader
can act on.

Reference backend, emcee, 200 points, 32 walkers x 4 000 steps (2 000
discarded):

============= ========== ====== ====== ====== ====== ==============
scenario      likelihood *A*    *B*    *d1*   *d2*   68 % covers?
============= ========== ====== ====== ====== ====== ==============
none          standard     0.23   0.12   0.48   0.54 yes
none          flexible     0.17   0.04   0.60   0.46 yes
mild          standard     0.07   0.54   2.83   3.51 **no**
mild          flexible     0.12   0.15   0.83   0.74 yes
strong_smooth standard     0.54   1.62   7.39  11.19 **no**
strong_smooth flexible     0.26   0.15   0.70   0.61 yes
strong_sharp  standard     4.71   3.36   1.00   1.33 **no**
strong_sharp  flexible     0.63   0.40   0.34   0.14 yes
============= ========== ====== ====== ====== ====== ==============

Read the ``strong_smooth`` row carefully, because it is the one that makes the
point. The standard likelihood puts ``d2`` **eleven of its own standard
deviations** away from the truth. That is not a slightly wrong answer with an
honest error bar; it is a precise answer that is wrong, and nothing in the fit
itself says so. The flexible likelihood's ``d2`` is 0.61 widths away — its
median has moved too, but its interval contains the truth, because the GP has
taken the ripple and the parameters' uncertainty has grown to match.

Note what is *not* claimed. The flexible likelihood does not recover a better
point estimate in the parameter's own units; sometimes it does not. What it
recovers is calibrated uncertainty, and ``tests/m2`` asserts exactly that:
every parameter within 1.5 posterior widths of the truth and inside the 68 %
interval, in every scenario.

The GP switches itself off when there is nothing to do
-------------------------------------------------------

The control, quantitatively. The fitted GP amplitude tracks the injected one
and collapses when there is none:

============= ====================== ==========================
scenario      GP amplitude (Jy)      GP length scale (µm)
============= ====================== ==========================
none          0.0012 ± 0.0011        0.0025 ± 0.0022
mild          0.0248 ± 0.0050        0.00089 ± 0.00020
strong_smooth 0.0942 ± 0.0288        0.00154 ± 0.00038
strong_sharp  0.0215 ± 0.0034        0.00070 ± 0.00013
============= ====================== ==========================

Against a continuum of about 1 Jy, the control's amplitude is a tenth of a per
cent — twenty times below the mild scenario's and eighty below the strong
one's, and consistent with zero at one standard deviation. The half-normal
prior is what makes that possible: its mass piles up at zero, so "no
deficiency" is a conclusion the fit can reach rather than a corner it is pushed
away from. The length scale in the control is correspondingly unconstrained —
there is no structure to set a scale — while in the fringing scenarios it
settles at a fraction of the 0.0028 µm ripple period, which is what a
stationary Matérn kernel needs in order to follow a sinusoid.

The diagnostics agree
---------------------

Two post-fit families, each answering the question it is scoped to
(``diagnostics.md`` §3.1 and §4). :func:`ampere.results.residual_whiteness` asks
whether a **standard**-likelihood fit left structure behind — a GP fit's
residuals are whitened by construction, so asking it there would be close to
circular. :func:`ampere.results.gp_localisation` asks *where* a flexible fit
needed its GP.

============= ============================== =================================
scenario      whiteness *p* (standard fit)   GP localisation (flexible fit)
============= ============================== =================================
none          0.185 — consistent with white  peak score 0.28, no structure
mild          0.005 — the floor, *Q* = 376   peak score 3.70
strong_smooth 0.005 — the floor, *Q* = 520   peak score 2.08
strong_sharp  0.005 — the floor, *Q* = 272   **peak at 0.86295 µm**, score 12.8
============= ============================== =================================

The 0.005 values are the resolution floor of a 199-permutation test, not a
coincidence: the smallest *p*-value that many permutations can report is
:math:`1/200`. The last row is the sharpest statement on this page — the
emission line was injected at 0.86300 µm, and GP localisation, which was told
nothing about it, put its peak 0.05 nm away. That is a third of one grid
spacing at 200 points, and its peak score is forty-five times the control's, so
"it peaked somewhere" and "it found something" are distinguishable.

Three backends, one posterior
-----------------------------

The same problem is composed from the reference backend's parts and sampled
with emcee, and from the torch and jax backends' parts and sampled with NUTS
over their realised (differentiable) log-densities. Two comparisons are made.

At the level of the **density**, the three agree to a relative tolerance of
:math:`10^{-10}` at a fixed parameter vector — in practice to the last few
bits — for the standard likelihood, the dense GP and the O(N) GP alike.

At the level of the **posterior**, disagreement is measured in posterior
widths, which is the only scale on which the question has an answer: two correct
samplers differ by their Monte Carlo error, and Monte Carlo error is a fraction
of the posterior rather than a number of Jy. At the milestone budget (emcee
32 × 4 000; NUTS 4 chains × 1 000 draws after 1 000 warm-up), on
``strong_smooth`` at 200 points, the **worst** disagreement over both
likelihoods and all four parameters is:

============ ================= ============================
comparison   median            68 % interval endpoint
============ ================= ============================
torch vs ref 0.072 widths      0.109 widths
jax vs ref   0.053 widths      0.080 widths
============ ================= ============================

torch and jax are compared with each other through the reference backend rather
than directly, because the two live in separate environments and no interpreter
has both. ``tests/m2`` asserts 0.5 widths on the median at its own shorter
budget and 0.20 under ``-m m2_full``; both are several times the Monte Carlo
error at the budget used, which is what makes them assertions about the backends
rather than about the samplers' luck.

The size ladder, and why it matters
-----------------------------------

The study runs at 200, 2 000 and 20 000 points over the same band, so the ladder
samples one spectrum more finely rather than observing more of it. The physics,
the deviation and the noise level are held fixed; what changes is the amount of
data.

The consequence is the reason misspecification matters more, not less, for real
spectra. The posterior narrows as :math:`\sqrt{N}`; the unmodelled 7 % ripple
does not shrink at all. So the standard likelihood's error, measured in its own
standard deviations, *grows* along the ladder. Measured, on ``strong_smooth``,
at the milestone budget — the worst of the four parameters at each rung:

========== ================================= ==============================
N          standard: worst offset (widths)   flexible: worst offset (widths)
========== ================================= ==============================
200        11.2 (``d2``), truth excluded     0.70, truth covered
2 000      36.3 (``d2``), truth excluded     0.51, truth covered
20 000     113.0 (``d2``), truth excluded    0.57, truth covered
========== ================================= ==============================

The standard likelihood's error grows almost exactly as :math:`\sqrt{N}`
(11.2 → 36.3 → 113.0 against a predicted 11.2 → 35.4 → 112.0), which is what it
must do if the bias is fixed and the posterior is shrinking. At 20 000 points
the fit reports ``d2`` with an error bar a hundred times too small to contain
the right answer. The flexible likelihood's offset does **not** grow: it is
0.7, 0.5 and 0.6 widths, and its 68 % interval contains the truth at every
rung. That is the claim that matters for real spectra, which have thousands of
pixels rather than two hundred, and it is asserted by
``tests/m2/test_ladder.py`` under ``-m m2_full`` — which prints exactly the
table above.

The ladder is only reachable because of the O(N) solver. A dense Gaussian
process at 20 000 points is a 20 000 × 20 000 Cholesky factorisation per
likelihood evaluation, and NUTS needs one per leapfrog step.

When one length scale is not enough
-----------------------------------

Every one of the four scenarios injects **one** scale of deviation, which is
why one stationary length scale copes with all four and everything above holds
with a plain Matérn-3/2. W5.8 adds the scenario where it does not:
``many_lines``, a forest of five narrow lines confined to 0.860–0.870 µm — whose
correlation length is their width, 0.00035 µm — plus a smooth continuum error
across the whole band, whose correlation length is a hundred times larger. A
kernel with one ``length_scale`` must choose, and whichever it chooses the
other feature is left in the residual with nothing to absorb it.

Three kernels are compared on it, all three quasiseparable and so all three
O(N): the stationary Matérn-3/2 the study has used throughout, a **warped**
Matérn-3/2 (one monotone map of the coordinate, so one length scale covers two)
and a **sum** of two Matérn-3/2 terms (two length scales, added). Reference
backend, 200 points, 32 walkers x 1 100 steps (550 discarded):

============ ======================================= ================== ========
likelihood   kernel                                  worst offset       covered
============ ======================================= ================== ========
standard     —                                       35.70              1 / 4
flexible     ``Matern32``                            2.84               3 / 4
flexible     ``WarpedKernel(Matern32, input_warp=)`` 0.92               4 / 4
flexible     ``Matern32 + Matern32``                 0.78               4 / 4
============ ======================================= ================== ========

The second row is the one to read carefully. It is the **only** place in this
study where a flexible likelihood misses the 1.5-posterior-width threshold the
rest of the page is built on, and it misses it for a structural reason rather
than a numerical one: the kernel is the wrong *shape* for the deviation, not
too small or too slow. Both of the other two recover it, by factors of 3.1 and
3.7 on the worst offset, and both keep the truth inside all four 68 % intervals.

All three flexible fits localise the deviation **inside the line band** — all
three peak at 0.86295 µm, and each reports a higher anomaly score inside the
band than outside it — so the two non-stationary kernels are not winning by
having been handed somewhere else to put the residual. The standard fit's
residual-whiteness p-value is 0.005, the 199-permutation floor, as in the other
misspecified scenarios.

What the warp actually learned is worth reading off its posterior, because it
is the mechanism rather than the outcome. The six knots are the plain quantiles
of the band, 0.842 to 0.872 in steps of 0.006, so they were **not** placed
where the answer is; the fitted increments are

.. code-block:: text

    input_warp.scale         1.559 +- 0.25
    input_warp.increment0   -1.649 +- 0.54     0.842 - 0.848
    input_warp.increment1   -1.675 +- 0.50     0.848 - 0.854
    input_warp.increment2   -1.613 +- 0.52     0.854 - 0.860
    input_warp.increment3    1.345 +- 0.50     0.860 - 0.866
    input_warp.increment4    1.117 +- 0.50     0.866 - 0.872

and a segment's slope is :math:`\zeta(s\,u_k)/\zeta(0)` with
:math:`\zeta(u) = \log(1+e^u)`. The first three segments come out at a slope
of about 0.10 — the coordinate is compressed, so distances shrink and the
effective length scale there is ten times the base — and the last two at about
3.2, so the effective length scale in the line band is three times *shorter*.
The base length scale is 0.0029 µm, which puts the fit at roughly 0.028 µm
where the smooth arch is and 0.0009 µm where the forest is: a ratio of thirty
between two halves of one band, found by the data, from one length scale and
five increments under a prior centred on the identity warp. The three knots
that switch sign do so at 0.860 µm, which is where the forest starts.

The figure, written by ``--figures``, is the argument in one picture: the
injected deviation over each arm's conditioned GP mean, with the line band
shaded.

.. code-block:: text

    python -m examples.m2_misspecification.many_lines
    python -m examples.m2_misspecification.many_lines --figures /tmp/m2
    python -m examples.m2_misspecification --scenario many_lines --kernel warped

A sum of noise components needs a sparsity guard
-------------------------------------------------

The kernel algebra that makes the third arm possible is also the freedom to add
a component the data do not need, and the guard against it is a prior:
``ampere.core.regularised_horseshoe``, put on a kernel by
``ampere.core.with_shrinkage``. It is the recommended prior for any ``Sum`` of
noise terms — one global scale shared by every component, one local scale per
component under it, and each component's amplitude under its local scale.

The demonstration fits **two nearly degenerate** Matérn terms to a spectrum
whose deviation has exactly one smooth component, so the likelihood pins their
total and says almost nothing about how it is divided. What divides it is the
prior:

============ ==================== ==================== ====================
prior        larger amplitude     median min/max       P(min/max < 0.1)
============ ==================== ==================== ====================
flat         0.0118 Jy            0.254                0.257
horseshoe    0.0077 Jy            0.087                0.537
============ ==================== ==================== ====================

The redundant component's share falls by a factor of 2.9 and the posterior mass
at "the fit chose one component" rises by 2.1 — while the component the truth
*does* have survives, which matters: a prior that shrank everything would move
both numbers the same way and be useless. Both factors are asserted well below
what they measure, because the horseshoe's three levels give the posterior a
funnel and an ensemble sampler explores one unevenly; the numbers above were
reproduced across three run seeds and the assertions sit a third below the
worst of them.

``python -m examples.m2_misspecification.many_lines --shrinkage`` prints that
table. ``tests/m2/test_many_lines_calibration.py`` asks the harder question of
the warped arm — over many spectra drawn from the prior and deviated the same
way, do its credible intervals contain the truth as often as they claim? — by
simulation-based calibration, and pins that they do not come back too narrow.

What it costs
-------------

Measured by ``tests/benchmarks/test_m2_likelihood.py`` under
pytest-benchmark (``pixi run -e <env> bench``), which writes ``benchmark.json``
as a CI artefact per environment; the table below is a transcription of one such
artefact rather than a hand-timed number. **One log-density evaluation of this
study's own problem** at a fixed parameter vector, in milliseconds, median of
the timed rounds on one x86-64 development machine. Absolute values are
machine-dependent and the ratios much less so, which is why the paragraphs
below quote ratios; ``benchmark.json`` from your own ``pixi run bench`` carries
the machine and the interpreter alongside the numbers.

The contract path — :meth:`ampere.core.FittingProblem.log_prob`, the only
comparison legacy and all three backends can all run:

============================== ========== ========== ==========
row                            N = 200    N = 2 000  N = 20 000
============================== ========== ========== ==========
legacy ampere, dense RBF             1.90     352.2       --
v2 reference, no GP                  0.41       0.48       1.06
v2 reference, ``DenseGP``            1.17     143.8       --
v2 reference, ``QuasisepGP``         0.62       0.77       2.24
v2 torch, ``DenseGP``                1.60      97.2       --
v2 torch, ``QuasisepGP``             1.04       1.20       4.04
v2 jax, ``DenseGP``                  3.89     187.6       --
v2 jax, ``QuasisepGP``              28.2       27.8      32.8
============================== ========== ========== ==========

Four things to read off it. Legacy and v2's dense solver are the same order at
200 points, and both become impossible by 20 000 — legacy's kernel has no
quasiseparable form, so it has no third column and cannot have one. The O(N)
solver turns a 187× penalty at 2 000 points into none at all, and its cost at
20 000 points is roughly twice that of having no GP whatsoever, on a hundred
times the paper study's data. The dense rows on every backend grow by two to
three orders of magnitude between the first two columns, which is what
:math:`O(N^3)` looks like when you plot it. And the jax rows are flat rather
than fast, because on this path they are measuring XLA dispatch through the
numpy container boundary rather than arithmetic — 28 ms whatever N is. That is
not the path a jax user runs, and the next table is.

The realised path is the differentiable log-density
:func:`ampere.core.realise` builds and NUTS consumes, timed as **value and
gradient**, which is what one leapfrog step costs. There is no reference or
legacy row because neither has a gradient — the point of the redesign rather
than a gap in the table.

================================= ========== ========== ==========
row                               N = 200    N = 2 000  N = 20 000
================================= ========== ========== ==========
v2 torch, ``DenseGP``                   4.19     559.1       --
v2 torch, ``QuasisepGP``                2.33       2.75       8.21
v2 jax, ``DenseGP`` (jitted)            4.82     566.7       --
v2 jax, ``QuasisepGP`` (jitted)         0.28       0.61       4.70
================================= ========== ========== ==========

The jax O(N) row is the headline number of the whole redesign: **a value and a
gradient of a 20 000-point Gaussian-process likelihood in 4.7 ms**, where legacy
ampere needed 350 ms for a value alone at 2 000 points and could not reach
20 000 at all. It is also seven times faster than the same backend's contract
path at that size, and a hundred times faster at 200 points, which is the
measurement that says why the realised path exists at all: the dispatch cost
that dominates the contract row is paid once at trace time and never again.

The dense rows in this table are the same :math:`O(N^3)` wall the contract path
hit, now with a gradient attached — 560 ms per leapfrog step at 2 000 points, on
both backends, against 3 ms for the quasiseparable solve. A NUTS run takes
hundreds of those steps per draw.

A short sampling run per backend is tracked alongside
(``tests/benchmarks/test_m2_sampling.py``), because a sampler's cost is the
per-evaluation cost times a step count that depends on the geometry, and NUTS
and an ensemble sampler differ in both factors at once. At 200 points, on the
flexible likelihood: emcee 16 × 150 takes 1.5 s, jax NUTS 50 + 50 takes 4.3 s
and torch NUTS 50 + 50 takes 6.8 s — for 2 400 ensemble evaluations against a
hundred *gradient-guided* draws, which is a comparison of cost that says
nothing on its own about how many effective samples each bought.

Reproducing it
--------------

.. code-block:: text

    pixi run test-m2                       # the assertions, ~10 min
    pixi run -e jax test-m2                # and the jax cross-backend rows
    pixi run -e torch test-m2              # and the torch ones
    pixi run bench                         # the benchmark table -> benchmark.json
    pytest tests/m2 -m m2_full             # the full ladder, tens of minutes

Every number on this page comes from one of those commands. The thresholds the
tests assert live in one place —
``examples/m2_misspecification/study.py``'s ``FLEXIBLE_MAX_BIAS_WIDTHS``,
``STANDARD_MIN_BIAS_WIDTHS``, ``WHITENESS_STRUCTURE_LEVEL`` and their
neighbours — so this page, the tests and the study cite one set of numbers
rather than three transcriptions of it.
