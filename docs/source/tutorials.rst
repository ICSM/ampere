Tutorials
=========

Ampere v2
---------

.. toctree::
   :maxdepth: 2

   sed_composition
   photometry_spectra
   m2_misspecification
   interferometry
   astrometry
   image
   wstat_comparison
   sbi
   population
   notebooks/quickstart
   notebooks/Ampere_MBB_Example

:doc:`notebooks/quickstart` is the shortest route from nothing to a v2 fit —
a straight line fitted to two photometric bands and a *Spitzer*/IRS spectrum,
first with independent noise and then with the flexible likelihood, ending in
the :mod:`ampere.results` plots — and :doc:`notebooks/Ampere_MBB_Example` fits
a modified blackbody to ten bands of photometry with emcee and then zeus
behind one interface. Both are executed when the documentation is built (a
cell that raises fails the build), state their sampling budgets and wall
times in their last cell, and are the notebook forms of ``examples/linear_sed``
and ``examples/modified_blackbody``.
:doc:`sed_composition` is the simplest composition there is — one model, a
spectrum and a photometric catalogue, two instruments on one channel — and
the page to start with if you have not built a multi-dataset fit before.
:doc:`photometry_spectra` extends that composition to three observations of
one source — two spectrographs of different resolving power, each with its
own uncertain calibration factor, and a catalogue — with the factors tied or
left free and the flexible likelihood as the complement for the residual
calibration alone cannot explain.
:doc:`m2_misspecification` is the flagship study: a deliberately misspecified
spectrum, fitted with and without the flexible likelihood, at three data
sizes and on all three backends, with the timings. :doc:`interferometry` is
the template for adding a **new modality** — a kind-changing,
coordinate-changing step, a complex container, a wrapped angular family, and
the M2 question asked of a resolved binary observed as visibilities and
closure phases, plus the chromatic case a flexible likelihood over ``(u, v)``
alone cannot express. :doc:`astrometry` is the second worked modality —
built by following that template, to test the template's own claim — a
reflex orbit on two ``TimeSeries`` channels, an epoch-sampling instrument,
and the time-domain flexible likelihood on the O(N) ``QuasisepGP`` path a
single ordered coordinate axis is exactly the right shape for; its closing
section reports what the template did not say. :doc:`image` is the third,
and the first whose observed container is **gridded** — a PSF convolution on
an ``Image``, the mask rule rewritten for a grid, and the
``DenseGP``-against-Hilbert-space-GP benchmark Phase 5's "chosen by
measurement" rule asks for; its closing section is about what a ``Layout``
other than ``POINTS`` changes. :doc:`wstat_comparison`
works through registering a **user-defined likelihood family** — a profiled
Cash statistic with background — which is the extension point to reach for
when your data are not Gaussian. :doc:`sbi` is for the opposite case — a
model with no likelihood to write down at all — fitting a black-box
simulator with :class:`~ampere.inference.SBIEngine`, caching the trained
posterior, truncated marginal ratio estimation, and checking that the
result is calibrated. :doc:`population` is Phase 5's addition — one
population of objects declared with :class:`~ampere.core.Population`,
fitted jointly on the native path and by reweighting archived single-object
fits, the two routes ``hierarchical_population.md`` left open.

The first seven have runnable counterparts in the repository:
``examples/sed_composition``, ``examples/m2_misspecification``,
``examples/interferometry``, ``examples/astrometry``, ``examples/image``,
``examples/wstat_comparison.py`` and ``examples/sbi/``, each covered by its
own test suite so that none can rot unnoticed. The v2 twins of six classic
legacy examples are listed in ``examples/README.md`` and set beside their
legacy originals, with snippets, in :ref:`the migration guide
<migrating-side-by-side>`;
among them ``examples/phoenix_star`` and ``examples/star_disc`` fit stars with
ampere's own PHOENIX emulator, a PCA plus one Gaussian process per weight that
evaluates identically on the reference, torch and jax backends, and
``examples/cstar`` fits a dusty carbon star with the Hyperion radiative-transfer
code by simulation-based inference, pooling Hyperion runs across worker
processes (it needs the pixi ``hyperion`` environment; see :doc:`install`).
:doc:`population` has none —
its code is small enough to walk through inline, and it is exercised
instead by ``tests/inference/test_population_nuts.py`` and
``tests/results/test_population.py``.

Legacy tutorials
----------------

.. warning::

   This notebook teaches the **legacy** v1 API (``ampere.infer.sbi``), which is
   frozen — see :doc:`legacy`. It is kept as the reference for the legacy
   embedding-network dictionary (:doc:`advanced` cites it), and because
   Phase 3 replaced it rather than deleted it: :doc:`sbi` is its v2
   counterpart. It is rendered from what it contains and is **never executed**
   by the documentation build, because it needs ``torch`` and ``sbi``, which
   the docs environment deliberately lacks. The other two notebooks that used
   to sit here, ``quickstart`` and ``Ampere_MBB_Example``, are now written on v2
   and are listed above with the v2 tutorials; the legacy code they used to
   show is covered name by name in :doc:`migrating`.

.. toctree::
   :maxdepth: 2

   notebooks/Embedding_nets

Still to be written
-------------------

Conditional priors and arbitrary priors each deserve a page of their own;
:doc:`overview` and the API reference carry what there is for now.
Combining different data types is no longer on this list —
:doc:`sed_composition` is that page.
