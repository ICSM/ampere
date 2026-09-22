Adding a modality: astrometric time series
===========================================

Phase 4's second worked modality (W4.9), and its purpose is different from
the first: :doc:`interferometry` proved the redesign end to end by being
hard — a kind-changing, coordinate-changing step, complex data, a wrapped
angular family. This page's job is to be *easy*, and to say honestly whether
it was — Peter's ruling on 2026-09-11 was "this will be an essential test of
the code", meaning the test is of :doc:`interferometry` itself: can a second
modality nobody wrote the contracts for be built by following that page,
with no new contract surface? Section order below follows that page's,
section for section, so a gap is visible as a gap rather than as a
rewritten structure. If you have not read :doc:`interferometry` yet, read it
first — this page cites it throughout as "the template".

The runnable example is ``examples/astrometry``
(:download:`generators.py <../../examples/astrometry/generators.py>`,
:download:`astrometry.py <../../examples/astrometry/astrometry.py>`,
:download:`__main__.py <../../examples/astrometry/__main__.py>`): a synthetic
source drifting under proper motion and wobbling under an unseen companion's
reflex motion, observed at a dozen or so irregular epochs in both right
ascension and declination, fitted with emcee on the reference backend or
NUTS on torch/jax. Run it yourself::

    python -m examples.astrometry                      # reference, emcee
    python -m examples.astrometry --backend torch       # NUTS
    python -m examples.astrometry --gp                  # the flexible likelihood

``tests/astrometry``, ``tests/backends/test_reference_astrometry.py``,
``tests/backends/test_native_astrometry.py`` and
``tests/examples/test_astrometry_example.py`` are this page's own coverage.

1. The kind: nothing to add
-------------------------------

The template's §1 says a new modality usually starts with a
:class:`~ampere.core.results_schema.FunctionSamples` subclass. This one does
not: :class:`~ampere.core.TimeSeries` already existed, unamended, since
before Phase 4 — one ``time`` axis,
``order=Order.STRICTLY_INCREASING``, ``Layout.POINTS``, real-valued. A 2D sky
position is two channels of it, not a new kind: ``results_schema.md`` §15.2's
rule for a vector-valued observable ("one value array per container",
written for Stokes I/Q/U/V) applies equally to right ascension and
declination, so :class:`~ampere.backends.reference.astrometry.ReflexOrbit`
emits ``{"ra": TimeSeries, "dec": TimeSeries}`` from one model rather than
inventing a two-column container. **The first finding, and a genuinely
positive one**: a template written for one modality's kind correctly
anticipated that the next modality might need no kind at all — the three
class attributes the template's §1 asks a reader to weigh (``AXES``,
``LAYOUT``, ``ALLOW_COMPLEX``) is exactly the checklist that says "this
already exists" as readily as it says "declare a new one".

2. The step: what it takes from the observed data, and what it publishes
--------------------------------------------------------------------------

:class:`~ampere.backends.reference.astrometry.EpochSample` is
``TimeSeries -> TimeSeries``, kind-preserving (``PRODUCES`` left at its
default of ``None``) — the "epoch sampling" row
``transformations.md`` §10's standard-library table already named before
this item built it. Like :class:`~ampere.backends.reference.interferometry.FourierSample`,
its coordinates are a **buffer taken from the observed container at
construction** (``EpochSample.from_observed(container)``), never recomputed:
the same rule, applied to a time coordinate rather than a baseline. Unlike
``FourierSample``, it has no arithmetic at all — ``apply`` is the identity,
because once negotiation has adopted the published ``points=`` requirement
the container it is handed is already on exactly the right epochs. **The
step publishes a requirement and does nothing else**, which is the
template's §2 lesson taken to its logical end: a step's job is sometimes
entirely the declaration, with an empty ``apply``.

3. ``configure_from``: not needed here
------------------------------------------

The template's §3 is about a step reading ahead in its own chain
(``configure_from``, for a smearing step that needs to know what the
Fourier step ahead of it computed). Every instrument chain in this modality
is one step long — ``EpochSample`` alone — so there is nothing for
``configure_from`` to do, and this module does not override it. **This is
the first place the template does not tell a reader what to do**: nothing
in :doc:`interferometry` says whether a one-step chain needs to say so
explicitly, or whether "the hook exists and this step does not use it" is
itself worth documenting. This page says so, here, because a reader
following the template section by section would otherwise wonder whether
they had missed something.

4. The two-channel composition
-----------------------------------

Where :doc:`interferometry`'s §4 composes **two instruments on one model
channel** (visibilities and closure phases, both reading the ``sky``
channel), this modality composes the mirror image: **one model, two
channels**, each with its own instrument and its own dataset:

.. code-block:: python

    from ampere.core import (Dataset, DatasetCollection, FittingProblem, Instrument,
                             Likelihood, GaussianFamily, IndependentNoise, negotiate)
    from ampere.backends.reference import ReflexOrbit, EpochSample

    # 1. One model, two channels "ra" and "dec".
    model = ReflexOrbit(epochs, pmra=..., pmdec=..., period=..., phase=..., amp_ra=..., amp_dec=...)

    # 2. One instrument per channel, each with one step.
    ra_instrument = Instrument([EpochSample(epochs)], channel="ra", label="astrom_ra")
    dec_instrument = Instrument([EpochSample(epochs)], channel="dec", label="astrom_dec")

    # 3. Two datasets, independent Gaussian noise (the rigid comparison).
    datasets = DatasetCollection({
        "ra": Dataset(observed_ra, ra_instrument,
                      likelihood=Likelihood(GaussianFamily(), IndependentNoise()), label="ra"),
        "dec": Dataset(observed_dec, dec_instrument,
                       likelihood=Likelihood(GaussianFamily(), IndependentNoise()), label="dec"),
    })

    # 4. One negotiation, one compile, the model built once per draw.
    problem = FittingProblem(model, datasets, seed=20260913)
    problem.requirements["model"]["ra"].sources    # ('astrom_ra',)
    problem.requirements["model"]["dec"].sources   # ('astrom_dec',)

Because the two instruments bind **different** channels, the label
collision :doc:`sed_composition` and :doc:`interferometry` both have to walk
through — two instruments defaulting to the same label on one shared
channel — cannot happen here at all; explicit labels are given anyway, for
symmetry with the other pages, not because anything would collide without
them. ``period`` and ``phase`` are shared between the two channels *at the
language level*: both are read from the same ``values`` mapping inside one
:meth:`~ampere.backends.reference.astrometry.ReflexOrbit.evaluate` call, so
the orbital parameters are tied by construction and need no explicit
``Tie`` — the same "one model, several channels" pattern the design sketch's
§1 describes for the reflex orbit's shared period and phase.

``ReflexOrbit.compile_for`` adopts each channel's own negotiated epoch grid
**independently** — a strict generalisation of the design sketch's own
worked example, which reads the epochs off the ``"ra"`` channel alone and
assumes ``"dec"`` shares them. A caller who genuinely observes RA and Dec at
different epochs (a real possibility — a scanning astrometric mission may
not measure both coordinates simultaneously) gets the physically correct
answer; one who observes them together, as ``examples/astrometry`` does,
loses nothing.

5. The likelihood: the O(N) path this modality is chosen to exercise
-------------------------------------------------------------------------

:class:`~ampere.core.TimeSeries` is ``Layout.POINTS`` with one ordered
one-dimensional coordinate axis — exactly what
:attr:`~ampere.core.GPSolver.REQUIRES_ORDERED_1D` asks for, so
:class:`~ampere.core.QuasisepGP` composes on both channels with no refusal
at all:

.. code-block:: python

    kernel = Matern32(st.halfnorm(scale=0.05), st.loguniform(10.0, 1000.0), axes=("time",))
    noise = GaussianProcessNoise(kernel, QuasisepGP())

This is the template's §5 lesson working in the *other* direction from
:doc:`interferometry`'s: there, a ``(u, v)`` point has no ordering under any
axis choice and ``QuasisepGP`` refuses by name, so the flagship GP runs on
``DenseGP`` alone; here, the O(N) solver is not merely available but is the
natural choice, since a light curve or an astrometric time series is
exactly the shape ``DEVELOPMENT_PLAN.md`` §2 built the quasiseparable path
for. ``tests/conformance/test_astrometry.py`` holds ``DenseGP`` and
``QuasisepGP`` to agreement at ``cross_solver`` tolerance on both channels.

The family is :class:`~ampere.core.GaussianFamily` — real-valued, no
wrapping, no complex container — so unlike :doc:`interferometry`'s wrapped
closure phases, a GP composes here as an **analytic** marginalisation on
every backend and every engine, emcee included. There is no "blind alley"
of a latent composition to avoid, because nothing about this modality's
values needs one.

6. A measured finding: period aliasing, not a template gap
-----------------------------------------------------------------

``examples/astrometry``'s own study found something the template's models
never could: a reflex orbit is **periodic**, and every model
:doc:`interferometry` fits (a Gaussian, a uniform disc, a binary) is not. A
period prior wide enough to reach past the observed epochs' own baseline
(a ``loguniform(50, 2000)`` was tried first, against a twelve-epoch,
860-day baseline) is a genuinely multi-modal problem: an emcee ensemble and
a multi-chain NUTS run both, measured directly, locked onto spurious
periods and stayed there — not a slow mixing problem a longer chain fixes,
but a real second mode a chain that started near it has no reason to leave.
This is not a defect in the composition (every conformance row above holds
regardless), and it is not new astrophysics either — period aliasing under
sparse sampling is a known problem in the literature this modality
represents. It is a fact about *this kind of model* that a template built
from non-periodic sources had no occasion to teach: **a modality whose
model is periodic needs an informed period prior** (a `norm(400, 30)`
around a period a periodogram or a previous epoch has already suggested,
which is what ``examples/astrometry`` fits), not a search over decades, or
its own recovery test becomes a report on prior-driven aliasing rather than
on the composition. See §10 below (What the template did not say) for this
as the item's own carried finding.

7. What the conformance suite owes a new modality
-----------------------------------------------------

``tests/conformance/test_astrometry.py`` pays the same four rows
:doc:`interferometry`'s §7 lists, adapted to this modality's shape:

- **a row per model against a closed form** — the reflex orbit's ``ra`` and
  ``dec`` offsets against the closed-form ephemeris in
  :mod:`tests.conformance.oracles`, and a row confirming ``period``/``phase``
  are read from one shared evaluation;
- **the step's published requirement**, checked directly — ``points=`` at
  the exact observed epochs — and its identity ``apply``, mask included;
- **a two-channel row** — one model, two channels, both recorded in
  ``problem.requirements``, evaluated once per draw (not once per channel);
- **a draw row** — ``FittingProblem.simulate`` producing both channels from
  one model evaluation.

A fifth row this modality adds, beyond the template's four:
``DenseGP``/``QuasisepGP`` agreement on a ``TimeSeries``, which
:doc:`interferometry` could not write (its own O(N) path is refused by
name) and which is the whole reason W4.9 exists.

8. The native twins: the inheriting pattern, unchanged
-----------------------------------------------------------

``ampere.backends.{torch,jax}.astrometry`` follow the **inheriting**
pattern the template's §8 chose for interferometry — each class derives
from its reference counterpart and overrides only the capability flags,
this backend's device placement, and the arithmetic:

.. code-block:: python

    class ReflexOrbit(_ReferenceReflexOrbit):
        ...  # only _offset is new; compile_for and the two-channel template are inherited

The template's §8 states a *rule* for choosing between inheriting and
re-declaring: inherit when the declaration is the dangerous half and the
arithmetic is cheap; re-declare when the declaration is thin and the
arithmetic is the whole of the step. By that rule, this modality is
**ambiguous in the opposite direction** from interferometry's own case:
``EpochSample``'s declaration is a few lines and its arithmetic is *none*
at all (the identity), and ``ReflexOrbit``'s declaration (the two-channel
template, ``compile_for``) is not obviously more dangerous to duplicate
than its four-line ``sin``/``cos`` arithmetic is cheap to fork. This item
inherited anyway, for consistency with the one pattern the codebase has
used so far and because there is genuinely nothing to gain from
duplicating a two-channel template that must stay bit-identical across
backends — but the template's own stated rule does not settle this case
either way, which is this page's second concrete finding for W4.8 (see §10).

9. The study: recovering an injected orbit
------------------------------------------------

``examples/astrometry`` injects a proper motion, a 400-day period, a phase
and a reflex semi-amplitude in each coordinate, observes it at a dozen or
so irregular epochs (0.03 mas uncertainty per epoch), and fits all six
parameters with an informed period prior (§6). Measured at the pinned seed
(``examples.astrometry.generators.SEED``), every parameter lands inside its
central 95 % credible interval on every backend:

.. code-block:: text

    emcee, reference (24 walkers x 1000 draws, ~16 s):
        pmra       1.200 +/- 0.003   (truth 1.2)     95%[1.194, 1.207]
        pmdec     -0.594 +/- 0.004   (truth -0.6)    95%[-0.602, -0.583]
        period    395.0  +/- 12.4    (truth 400)     95%[361.6, 401.5]
        phase       1.34 +/- 1.74    (truth 0.7)     95%[0.64, 5.96]
        amp_ra     0.672 +/- 0.198   (truth 0.6)     95%[0.583, 1.216]
        amp_dec    0.399 +/- 0.129   (truth 0.35)    95%[0.336, 0.762]

    NUTS, torch (2 chains x 200 draws, ~34 s):
        pmra       1.2004 +/- 0.0030 95%[1.1946, 1.2062]
        pmdec     -0.5985 +/- 0.0049 95%[-0.6068, -0.5896]
        period    388.7   +/- 11.1   95%[376.9, 401.4]
        amp_ra     0.5914 +/- 0.0104 95%[0.5710, 0.6089]
        amp_dec    0.3344 +/- 0.0184 95%[0.3048, 0.3651]

    NUTS, jax (2 chains x 150 draws, ~7 s):
        (agrees with torch to three significant figures throughout)

``tests/astrometry/test_recovery.py`` pins this claim at a per-PR budget
(the same shape, smaller): every one of the six parameters inside its
central 95 % interval, on the reference backend under emcee and on torch
and jax under NUTS. NUTS needs an order of magnitude fewer draws for
comparable or tighter intervals than emcee — unsurprising given the
gradient, and the same relationship :doc:`sed_composition` and
:doc:`interferometry` both report.

10. One correlated process over both channels (W5.9)
-----------------------------------------------------------

The design sketch this modality came from recorded a gap and dispositioned
it to Phase 5: a **single** correlated noise process over RA and Dec
together, rather than the two independent ones W4.9 shipped. W5.9 closes
it, and astrometry is the first customer because this page's composition
already has the shape the model needs — one evaluation producing two
channels on one epoch grid.

The model is the **shared-grid intrinsic coregionalisation model**::

    K = B ⊗ K_x

with ``K_x`` an ordinary :class:`~ampere.core.Matern32` along the epochs and
``B`` a 2×2 positive-definite matrix saying how the two sky axes move
together. ``B`` is declared in the coordinates the physics is stated in —
:class:`~ampere.core.RotationCoupling` takes a position angle and two
log-variances, so ``B`` is an error ellipse on the sky with a major axis, a
minor axis and an orientation.

Why this is not two GPs with tied hyperparameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two independent GPs can reproduce each axis's *marginal* scatter exactly.
What they cannot express is the **cross-covariance**: that the RA and Dec
residuals at one epoch move together. A centroiding systematic with a
preferred direction — the ordinary case for a ground-based astrometric
solution — is almost entirely cross-covariance, so "two GPs with tied
hyperparameters" is not an approximation of the joint model. It is a
different model, which happens to agree about every marginal.

Where it is declared, and why not on the likelihood
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :class:`~ampere.core.Likelihood` scores one dataset and this scores two,
so the joint noise model goes on the
:class:`~ampere.core.DatasetCollection`::

    DatasetCollection(
        {"ra": Dataset(...), "dec": Dataset(...)},
        joint={"astrom": JointGaussianProcessNoise(
            Matern32(1.0, 150.0, axes=("time",)), QuasisepGP(),
            datasets=("ra", "dec"),
            coupling=RotationCoupling(angle_prior, variance_prior, variance_prior),
        )},
    )

Each member dataset keeps a bare ``Likelihood(GaussianFamily())``. The group
joins the joint parameter space as one further component, so its parameters
are ``astrom.angle``, ``astrom.log_variance_0`` and
``astrom.log_variance_1``, and it contributes **one** log-likelihood term
under the label ``"astrom"`` in place of ``"ra"`` and ``"dec"``'s separate
ones — the first use of ``DatasetCollection.contributions`` as something
other than one term per dataset.

It is still O(N)
~~~~~~~~~~~~~~~~~~~~~~~~

Diagonalise ``B = Q Λ Qᵀ`` and rotate the two residual vectors by ``Qᵀ``.
The rotated outputs are independent, each an ordinary scalar GP with
covariance ``λ_s K_x + diag(σ²)``, so the joint density is two
``QuasisepGP`` solves rather than one dense factorisation of a 2N×2N
matrix. Exact, not approximate: ``tests/conformance/test_astrometry.py``
scores the same residual both ways on every backend and holds them to the
solver tolerance.

One restriction follows from the rotation itself and is checked at
composition: the two channels must carry the **same per-epoch
uncertainties**. ``Qᵀ ⊗ I`` leaves ``diag(σ²)`` diagonal only where every
channel's ``σ`` is the same vector; heteroscedasticity *along* the epoch
grid is unaffected. Unequal per-channel errors, mismatched grids and the
general LMC belong to the dense/reduced-rank follow-on
(``likelihoods.md`` §15).

What it buys: the calibration study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``python -m examples.astrometry --sbc joint`` (and ``--sbc independent``,
``--sbc rigid``) runs simulation-based calibration of each arm against data
carrying an injected correlated centroiding systematic. Both arms are given
the noise process they are entitled to know — the joint arm the whole of
``B``, the comparison arm each channel's *correct marginal* amplitude — so
that the one thing which differs between them is the cross-covariance.

What is ranked is not a parameter. It is
``(pmra + pmdec)/sqrt(2)``, the **diagonal of the proper-motion plane**, and
the reason is the mechanism:

* A cross-channel systematic does not bias either channel's own parameter.
  It **correlates their errors**. ``pmra`` is measured from ``ra`` alone and
  ``pmdec`` from ``dec`` alone, both with the same weight along the epoch
  grid — a proper motion is a linear trend — so an error shared by the two
  sky axes makes the two parameter errors move together, at about 0.86 here.
* Each *marginal* posterior is therefore still about the right width, and
  the study reports that: both arms' marginal coverage on ``pmra`` and
  ``pmdec`` comes out nominal. That is a real and easily missed finding —
  checking marginals alone would have said the independent model was fine.
* The **joint** posterior is where the difference lives, and the projection
  that sees it is the one along which the two errors add. Its true variance
  is ``v(1 + rho)``; a model that believes the errors independent reports
  ``v``, understating the interval by ``sqrt(1 + rho)``, about 1.36. That is
  undercoverage of the *direction* of a measured proper motion, which is a
  quantity astronomers publish.

``tests/examples/test_astrometry_example.py`` pins the comparison, behind
the ``astrometry_full`` marker because two arms of forty-eight refits each
is well over ten minutes; a reduced-budget sibling runs on every PR and
checks the machinery rather than the claim.

A note on the parameterisation. Adding ``pi/2`` to the angle and exchanging
the two log-variances gives the same ``B``, so the *matrix* is identified
while the three parameters are identified only up to that relabelling —
the label switching a mixture model has. The density is unaffected; the
marginal on the angle is bimodal. Summarise ``B`` itself, or fix the angle
where the instrument's own is known, which is what the calibration study
does. **W5.28(c)** gives the first route a function:
:func:`~ampere.results.diagnostics.coupling_matrix_summary` reassembles
``B`` per posterior draw and hands it to :func:`~ampere.results.summary`'s
own machinery, rather than reporting the angle's own mean and HDI across
the relabelling jump.

11. What the template now says, and where
-------------------------------------------------

This section originally listed every gap, ambiguity or missing instruction
this item met while following :doc:`interferometry` section by section, for
W4.8 to fix in the template itself. **W4.8 did**, so this section now says
where each fix landed rather than that it was missing — the item's real
deliverable turned into the template's own text, and this page keeps the
record of what changed and why.

- **§3 (``configure_from``)**: the template showed a step *using* the hook
  but never said what a modality with a one-step chain should do about it.
  :doc:`interferometry` §3 now states it directly, citing this page's
  ``EpochSample`` as the example: nothing is the right answer for a
  one-step chain, and a modality that says so explicitly (in its docstring
  or its page, as this one now does above) has made a decision rather than
  skipped a section.
- **§8 (inherit vs. re-declare)**: the stated rule ("declaration dangerous
  and expensive => inherit; arithmetic dominant => re-declare") did not
  resolve this modality's case, where the declaration is cheap *and* the
  arithmetic is cheap. :doc:`interferometry` §8 now carries the third
  clause this page asked for: when both halves are cheap, inherit anyway,
  for the one-pattern-per-codebase argument, not because duplicating either
  half would be dangerous.
- **A periodic model is a hazard class the template had no case for.**
  :doc:`interferometry` now closes with "A hazard the models here never
  raise: periodic models and prior width", naming this page's §6 finding
  directly: a modality whose model is periodic (or otherwise genuinely
  multi-modal in its own right, independent of the likelihood/noise
  machinery) needs its own discussion of prior width, separate from the
  composition questions the rest of the template is about.
- **The "two instruments, one channel" vs. "one model, two channels"
  distinction is now named.** :doc:`interferometry` §4 is retitled "The
  two-dataset composition: two instruments on one channel" and opens by
  naming both shapes and pointing here for the other one, so a reader
  arriving from that page alone no longer has to discover by surprise that
  a multi-dataset composition can look like this page's instead.
- **What the template got right, stated as a finding too, and left
  unchanged**: §1's "three class attributes are the whole of what a kind
  declares" correctly predicted that a kind might already exist and need no
  declaration at all — the strongest possible form of "no new contract
  surface" is not needing a new class, and the template's own emphasis on
  those three attributes (rather than on the mechanics of writing a new
  subclass) is what made that visible immediately rather than after a false
  start. Nothing about this needed fixing, so :doc:`interferometry` §1 is
  unchanged.

See also
------------

* :doc:`interferometry` — the template this page follows and tests, section
  for section.
* :doc:`sed_composition` — the simple composition case both pages cite.
* ``docs/design/modalities/astrometric_timeseries.md`` — the design sketch
  this item's ``ReflexOrbit`` and ``EpochSample`` are drawn from, including
  the joint 2-vector GP gap dispositioned to Phase 5. **W5.9 closed that
  gap**; section 10 above is the result.
