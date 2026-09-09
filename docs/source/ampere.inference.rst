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
  layout must never be believed about another.

Requires the ``sbi`` extra; both it and torch are imported inside ``run``, so
``import ampere.inference`` in the base install pulls in neither.

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

.. autodata:: ampere.inference._sbi.SUMMARY_LAYOUT
