"""The RHMF exploratory trial (W6.9): does pre-fit robust factorisation localise misspecification?

``diagnostics.md`` §2 describes a family of **pre-fit** screens: factorise a
*collection* of comparable spectra into a low-rank model with
`Robusta-HMF <https://github.com/TomHilder/robusta-hmf>`_ (robust
heteroskedastic matrix factorisation, in JAX) and read the iteratively
reweighted least-squares weights as "how much did this low-rank model trust
this data point". Where the weights fall, a smooth model of the collection is
going to struggle, and so, plausibly, will whatever physical model is fitted
next. W2.7 deferred the family because ``robusta-hmf`` did not clear the
maturity bar of ``DEVELOPMENT_PLAN.md`` §2; W6.9 is the trial that item asked
for, run against a pinned commit, with the adoptability re-check repeated.

**This is a trial, not a feature.** ``ampere.diagnostics`` does not exist and
nothing under ``ampere/`` imports ``robusta_hmf``. The extra is non-default and
exists for this directory alone; the adapter is three functions in
:mod:`examples.rhmf_trial.trial`, not a public API, and no hyperparameter
default is promoted by anything it prints.

What it does
------------
:func:`~examples.rhmf_trial.trial.to_matrix`
    Aligned :class:`~ampere.core.Spectrum` containers to ``(Y, W, coordinates)``:
    values, inverse-variance weights, and an ampere mask as zero weight
    (Lesson R1).
:func:`~examples.rhmf_trial.trial.fit_rhmf`
    ``Robusta(rank, robust_scale).fit(Y, W)``.
:func:`~examples.rhmf_trial.trial.anomaly_score`
    The robust weights as an :class:`ampere.core.AnomalyScore` with
    ``provenance="rhmf_prefit"`` — the renderer
    (:func:`ampere.results.plot_anomaly_score`) already accepts it — either per
    feature for one object or per object by a low quantile across features
    (Lesson R2).

Two trials run from ``python -m examples.rhmf_trial``:

1. **Spectra.** A collection of M2 spectra (controls with no deviation, plus
   every scenario of
   :data:`examples.m2_misspecification.generators.EXTENDED_SCENARIOS`),
   factorised over a grid of ``rank`` and ``robust_scale``. It reports, per
   scenario, the score inside the injected deviation against outside it, with
   the controls' ratio beside it.
2. **Image.** W5.5's synthetic image, with its omitted smooth background as
   the injected deviation, under two flattenings: each image of a collection
   as one row, and the rows of a single image as the objects.

Running it
----------
::

    pixi run -e rhmf python -m examples.rhmf_trial --quick --out /tmp/rhmf
    pixi run -e rhmf python -m examples.rhmf_trial --out /tmp/rhmf

``--out`` names the run directory for the tables (CSV) and figures (PNG); it is
required and is never inside git. ``--quick`` finishes in under two minutes.
Without the extra the script raises
:class:`ampere.core.exceptions.OptionalDependencyError` naming it.

The findings, the hyperparameter sensitivity and the adoptability re-check are
in ``docs/design/contracts/diagnostics.md`` §7 ("Amended W6.9") and the
decision-log row of the same name in ``DEVELOPMENT_PLAN.md`` §2.
"""
