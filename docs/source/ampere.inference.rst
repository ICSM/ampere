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
* **the run carries two log-densities per draw.** ``lp``, ``log_prior`` and
  ``log_likelihood`` in ``sample_stats`` are the *true* ones, scored on the
  numpy contract path after the fit; ``ampere_sbi_log_prob`` beside them is
  the *estimator's* own, which is what makes calibration and importance
  reweighting possible. ``ampere_sbi_log_prob_kind`` says whether it is
  normalised (NPE) or known only up to the evidence (NLE, NRE);
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
sampled by rejection (``sample_with="mcmc"`` is the fallback for a narrow box).
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
