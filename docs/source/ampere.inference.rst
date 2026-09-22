ampere.inference
================

The engines, written once against the fitting-problem surface and run on
every backend. This namespace **imports no backend**, in any module, at any
depth; its only ampere imports are :mod:`ampere.core` and
:mod:`ampere.results`. The gradient-based drivers reach a backend's
differentiable density through :func:`ampere.core.realise`, which importing
that backend registered.

There is no way to sample through these drivers and not get a stored run:
``run()`` returns the ArviZ ``DataTree``, with per-draw ``log_prior`` and
``log_likelihood``, the per-dataset decomposition, the observed data and the
full provenance attrs.

Nested sampling: three engines, one convention
-----------------------------------------------

Three drivers here do nested sampling, and a run archived from any of them is
comparable with a run archived from the others without a reader knowing which
produced it. They share the live-point default (``max(100, 25 (n_dim + 1))``),
the equal-weight resampling rule, and the engine-neutral evidence triple
``ampere_log_evidence`` / ``ampere_log_evidence_err`` /
``ampere_evidence_method`` (the last is ``"nested_sampling"`` for all three —
the *family*, since ``ampere_engine`` already names the engine).

* :class:`~ampere.inference.DynestyEngine` is in the **base install** and is
  the one to reach for first. Static or dynamic, ellipsoidal or slice
  sampling, and no extra to install.
* :class:`~ampere.inference.NautilusEngine` (the ``nautilus`` extra) is
  importance nested sampling with a neural boundary: it typically reaches a
  given effective sample size in fewer likelihood calls than an ellipsoidal
  decomposition, which is the budget that matters when the forward model is a
  radiative-transfer code. It needs **at least two free parameters** — the
  library refuses fewer, and the driver says so at construction — and it
  reports no evidence uncertainty of its own, so ampere records the standard
  importance-sampling estimate ``1/sqrt(n_eff)`` and notes in
  ``ampere_nautilus_log_z_err_source`` that it did. ``options={"n_networks":
  0}`` turns the neural boundary off, which is the right setting for a cheap,
  low-dimensional problem.
* :class:`~ampere.inference.UltranestEngine` (the ``ultranest`` extra) is
  MLFriends region sampling with a bootstrapped termination criterion. It is
  the conservative one: designed not to miss a mode, with an evidence error
  whose bootstrap and tail halves are recorded separately
  (``ampere_ultranest_logzerr_bs`` and ``_tail``) and an insertion-order test
  (``ampere_ultranest_insertion_order_converged``) that says when they should
  not be believed. Like zeus it draws from numpy's process-global generator,
  so its run is not thread-safe against other code drawing from ``np.random``
  at the same time; ampere seeds and restores that state around the run, so a
  seeded problem still repeats exactly.

All three resample their weighted dead points to equal weight for the
``posterior`` group — ArviZ's ``posterior`` has no weight axis, so storing
weighted draws there would silently misreport every summary — and record the
original count (``ampere_<engine>_dead_points``) beside the resampled one. The
raw weighted output stays on ``engine.sampler``.

``tests/inference/test_nested.py`` is the **engine battery** these two landed
with: each engine's evidence checked against a two-parameter conjugate
linear-Gaussian problem whose marginal likelihood is written out in closed
form, each engine SBC-ranked through
:func:`~ampere.results.calibration.sbc`, and one cost record
(``ampere_engine_evaluations``) asserted per run.

Which variational guide to ask for
-----------------------------------

:class:`~ampere.inference.VIEngine`'s ``guide=`` argument **is** the
approximation being made, and since W5.14 there are four of them. All four are
fitted in the *unconstrained* space, so the constrained posterior each implies
is already warped by ampere's own bijections; all four record what they were
in ``ampere_vi_guide`` and in the engine-neutral ``ampere_approximation``.

* ``"normal"`` — a diagonal Gaussian. O(d) parameters, fast, and wrong in
  exactly the way a correlated posterior is correlated: it underestimates the
  marginal variances. ``ampere_approximation = "mean_field"``.
* ``"multivariate"`` — a full-covariance Gaussian fitted by maximising the
  ELBO, O(d²) parameters. The first thing to reach for when a corner plot
  shows the tilted ellipse every SED fit produces.
  ``ampere_approximation = "multivariate"``.
* ``"laplace"`` — a full-covariance Gaussian too, but taken from the
  *curvature*: SVI optimises to the MAP point and the covariance is then the
  inverse Hessian there, computed once. Cheap where the ELBO fit is dear,
  **exact when the posterior really is Gaussian**, and arbitrarily wrong when
  it is not — it describes one point, not the mass. It needs a second
  derivative of the realised density, which is a real demand on a backend's
  lowering. ``ampere_approximation = "laplace"``.
* ``"flow"`` — an inverse-autoregressive flow over a standard-normal base, and
  the only family here that can represent a skewed or heavy-tailed posterior.
  It costs a small neural network per transform and wants a longer
  optimisation than the Gaussian families (its ELBO trace is how you tell
  whether it got one). It needs **at least two free parameters** — there is
  nothing for an autoregressive flow to be autoregressive over below that, and
  the driver refuses a one-dimensional problem by name.
  ``ampere_approximation = "normalising_flow"``.

One practical difference is worth knowing before reading a run back. A
Gaussian guide's fit is kept on the engine as plain arrays
(``engine.guide_loc`` with ``engine.guide_scale`` or
``engine.guide_scale_tril``), so its density can be rebuilt with
``scipy.stats`` independently of which library fitted it. **A flow's cannot**:
its fit is a neural network's weights, all three attributes stay ``None``, and
the run says so in ``ampere_vi_guide_parameters`` (``"loc, scale"``,
``"loc, scale_tril"`` or ``"none"``). What every family does carry is the
per-draw ``sample_stats.proposal_log_density``, so
``exp(log_prior + log_likelihood - proposal_log_density)`` is an importance
weight buildable from the stored groups whatever the guide was.

The blackjax route, on jax
---------------------------

:class:`~ampere.inference.BlackjaxEngine` (the ``blackjax`` extra) is **jax
only**, by nature rather than by what happens to be installed: blackjax is a
jax library and there is no maintained torch counterpart to borrow, so a
problem on another backend is refused at construction by name. One dependency
buys two methods ampere does not otherwise have.

* ``method="mclmc"`` — microcanonical Langevin Monte Carlo (Robnik & Seljak):
  a deterministic isokinetic trajectory with periodic momentum decoherence and
  **no accept/reject step**, reported at several times NUTS's efficiency per
  gradient on smooth high-dimensional posteriors, which is what a latent GP
  over a spectrum or a hierarchical plate is. It has no default step size, so
  tuning is not optional: ``warmup`` is the tuner's budget in integrator
  steps, and it chooses the step size, the momentum decoherence length ``L``
  and a diagonal preconditioner. Because there is no Metropolis correction the
  chain carries a discretisation bias, traded against the tuner's energy
  variance target; the run records both facts
  (``ampere_blackjax_adjusted``, ``ampere_blackjax_desired_energy_var``), and
  ``ampere_approximation`` is still ``"none"`` because that key answers "do
  chain diagnostics apply?" — and for a Markov chain they do.
* ``method="pathfinder"`` — Zhang et al.'s quasi-Newton approximation: L-BFGS
  on the log density, a Gaussian built from the inverse-Hessian factors at
  every iterate, and the ELBO-best one kept. It emits one "chain" of i.i.d.
  draws like a variational fit, writes ``ampere_approximation =
  "pathfinder"``, and carries ``proposal_log_density`` so the approximation is
  correctable from the stored groups.

The second use for Pathfinder is as a **start-point search**:
``run(..., initial="pathfinder")`` fits it first and starts each MCLMC chain
at one of its draws, which is a far better start than a prior draw on a
posterior the prior barely covers. The whole thing is still reproducible from
the problem's seed — jax has no global RNG, so every stream is an explicit key
and two runs of a seeded problem agree bitwise.

Neither method estimates a marginal likelihood, and neither writes the
engine-neutral evidence triple. Pathfinder's ELBO is a single path's lower
bound of unknown tightness and is kept under
``ampere_blackjax_pathfinder_elbo``, where nothing will mistake it for a
nested sampler's estimate.

Simulation-based inference
--------------------------

:class:`~ampere.inference.SBIEngine` is the one driver here that does not
consume ``log_prob`` while it fits. Its training surface is
:meth:`~ampere.core.dataset.FittingProblem.simulate_many` and nothing else,
which is why it is the engine for a model whose likelihood cannot be written
down — a compiled radiative-transfer code behind a Python call, composed as a
black-box :class:`~ampere.core.Model` on the reference backend.
``executor=`` and ``chunk_size=`` pass straight through to ``simulate_many``,
so that budget runs under a process pool, a thread pool, dask, ray or
``MPIPoolExecutor`` without the driver knowing which;
``examples/sbi/external_simulator.py`` is that case end to end and
``examples/sbi/fit_external_simulator.py`` fits it.

Three things about an SBI run are worth knowing before reading one:

* **the prior is bridged into unconstrained coordinates.** The density
  estimator sees ℝⁿ, so a bounded prior needs no ``RestrictedPrior`` and an
  NPE posterior has no leakage to correct; the drawn posterior is mapped back
  through ``constrain`` before anything is stored, so the ``posterior`` group
  is in the user's own coordinates like every other run's. The run records
  ``ampere_sbi_parameterisation``;
* **the run carries three log-densities per draw.** ``lp``, ``log_prior`` and
  ``log_likelihood`` in ``sample_stats`` are the *true* ones, scored on the
  numpy contract path after the fit; ``ampere_sbi_log_prob`` beside them is
  the *estimator's* own, in the **unconstrained** coordinates it was trained
  on, which is what makes calibration and (with care about that coordinate
  choice) importance reweighting possible. ``ampere_sbi_log_prob_kind`` says
  whether it is normalised (NPE) or known only up to the evidence (NLE, NRE).
  **W5.0** adds ``sample_stats.proposal_log_density``, the results contract's
  engine-neutral name for the same idea, moved into the same **constrained**
  coordinates ``log_prior``/``log_likelihood`` are already in — so
  ``exp(log_prior + log_likelihood - proposal_log_density)`` is a valid
  importance weight from the stored groups alone, on any engine, and
  ``ampere_approximation`` is ``"density_estimator"`` for every method here
  (npe, nle, nre and tmnre alike);
* **the summary the network saw has a name.** ``ampere_sbi_summary_layout``
  is ``"flat"`` today — each dataset's observed values, masked samples
  dropped, concatenated in ``datasets`` order — and a network trained on one
  layout must never be believed about another;
* **a seeded problem's run is reproducible, network and all.** ``run()``
  seeds torch's own global generator (and, for an MCMC-sampled NLE, NRE or
  TMNRE posterior, ``numpy``'s legacy one, which is what ``sbi``'s default
  MCMC method actually draws through) from the problem's seed before every
  network it builds, before each round's training and before the final
  posterior draw — cache hit or not — so two runs of the same seeded problem
  give bitwise-identical draws. Both are restored once the seeded step
  finishes, so nothing about this leaks into unrelated code sharing the
  process. The seed used is recorded as ``ampere_sbi_torch_seed``, absent
  when ``problem.seed is None``, which asks for fresh randomness every time
  instead.

Requires the ``sbi`` extra; both it and torch are imported inside ``run``, so
``import ampere.inference`` in the base install pulls in neither.

Truncated marginal ratio estimation
------------------------------------

``method="tmnre"`` is Miller et al. (2021)'s TMNRE, expressed through ``sbi``
0.27's own ``NRE`` trainers and ``RestrictedPrior`` rather than through the
archived swyft implementation (which pins ``pytorch-lightning <= 1.9.5`` and
cannot be installed beside the torch the ``sbi`` extra resolves to). Two ideas,
and a run shows both:

* **marginal ratio estimation** — one classifier per 1-D ``θ_i``, and per
  unordered pair at ``marginals=2``, rather than one for the joint ratio. Each
  is trained on the same ``x`` against the corresponding *columns* of ``θ``, so
  none of them has to represent the joint's correlations, which is what makes
  the method work as the parameter count grows;
* **truncation** — after each round the *prior* is restricted to the
  hyperrectangle where those 1-D marginals exceed ``truncation_epsilon`` times
  their own maximum, the next round simulates inside it, and every estimator is
  retrained. The restricted prior is the original prior renormalised on a
  subset rather than a learned proposal, so the ratio is unchanged inside the
  box and no importance correction appears anywhere.

A TMNRE run emits ordinary joint draws: a joint estimator is trained alongside
the marginals in the final round, on the truncated prior, and its posterior is
sampled by MCMC over the same estimator by default (``sample_with="rejection"``
gives i.i.d. draws instead, at a cost that rises as truncation succeeds — ruled
2026-09-10 on the measured 353.6 s against 43.8 s).
Beside them the run carries a ``marginals`` group — each estimator's log-ratio
and estimated marginal posterior on a grid over the final box, in both the
unconstrained and the constrained parameterisation — and its truncation history
in ``ampere_sbi_truncation``, one record per round with the box, its
log-volume, the round's counts and the proposal's acceptance rate. The boxes
are nested by construction.

The cost is stated in the run: **a truncated estimator is not amortised.** The
box is chosen at the observed data, so the estimator must not be re-conditioned
on another observation, and ``ampere_sbi_amortised`` is ``0`` for every TMNRE
run whatever its round count. ``examples/sbi/tmnre_fit.py`` is the method end
to end on the toy joint problem.

.. automodule:: ampere.inference
   :members:
   :imported-members:
   :undoc-members:
   :show-inheritance:

Constants
---------

Re-exported from ``ampere.inference``; shown under the module that defines it.

.. autodata:: ampere.inference.engine.DEFAULT_CACHE_SIZE

.. autodata:: ampere.inference._sbi.METHODS

.. autodata:: ampere.inference._sbi.EMBEDDINGS

.. autodata:: ampere.inference._sbi.SET_EMBEDDINGS

.. autodata:: ampere.inference._sbi.LAYOUTS

.. autodata:: ampere.inference._sbi.SUMMARY_LAYOUT

.. autodata:: ampere.inference._sbi.TMNRE_SAMPLERS

.. autodata:: ampere.inference._tmnre.MARGINAL_ORDERS

.. autodata:: ampere.inference._tmnre.DEFAULT_TRUNCATION_EPSILON
