Adding a modality: interferometric visibilities
=================================================

Phase 4 asked one modality to prove the redesign end to end: a
kind-changing, coordinate-changing step in the middle of the chain, complex
data, a wrapped angular family, and two datasets sharing one sky model. This
page is written as the template for the *next* modality (W4.9's astrometric
time series, and whatever comes after it): each section below is a piece a
new modality supplies, in the order a reader adds one. If you have not built
a combined fit before, read :doc:`sed_composition` first — it is the same
five nouns with no kind change and no complex container, and this page cites
it throughout as "the case you already know".

The runnable example is ``examples/interferometry``
(:download:`generators.py <../../examples/interferometry/generators.py>`,
:download:`model.py <../../examples/interferometry/model.py>`,
:download:`study.py <../../examples/interferometry/study.py>`,
:download:`figures.py <../../examples/interferometry/figures.py>`,
:download:`__main__.py <../../examples/interferometry/__main__.py>`): a
synthetic resolved binary with a fainter, more extended disc, observed as
both visibilities and closure phases, fitted three ways to ask M2's question
of this modality, plus a fourth, chromatic scenario. Run it yourself::

    python -m examples.interferometry                    # 3 arms, reference, emcee
    python -m examples.interferometry --chromatic         # arm (d)
    python -m examples.interferometry --calibration       # the SBC coverage row
    python -m examples.interferometry --backend torch     # NUTS

``tests/interferometry`` and ``tests/examples/test_interferometry_study.py``
are this page's own coverage.

1. The kind: three class attributes, and a convention worth stating
---------------------------------------------------------------------

A modality's data starts as a new :class:`~ampere.core.results_schema.FunctionSamples`
subclass — :class:`~ampere.core.VisibilitySet` and
:class:`~ampere.core.ClosurePhases` here. Three class attributes are the
whole of what a kind declares about its shape:

.. code-block:: python

    class VisibilitySet(FunctionSamples):
        AXES = (
            AxisSpec("u", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
            AxisSpec("v", physical_types=("dimensionless",), equivalent_units=(u.rad**-1,)),
            AxisSpec("spectral_axis", physical_types=("length", "frequency", "energy"), order=Order.ANY),
        )
        LAYOUT = Layout.POINTS
        ALLOW_COMPLEX = True

``AXES`` names the coordinates and their units (a photometric point has one
axis, a spectrum one, an image two, a visibility three — one wavelength *per
sample*, ``Order.ANY`` because a ``(u, v)`` point set has no natural order);
``LAYOUT`` says whether the axes are a shared grid or a scattered point set;
``ALLOW_COMPLEX`` says whether a value may be complex. Nothing else about a
kind is negotiable, and nothing else needs to be — the whole of
``ampere.core``'s generic machinery (masking, unit checking, the likelihood
contract) reads only these three things.

:class:`~ampere.core.ClosurePhases` is the same three attributes, and it is
also the example of a kind carrying a **convention** rather than only a
shape. A triangle of telescopes has three baselines but only two are
independent (the third is minus their sum), and there are three ways to pick
which two to store. Two files that picked differently would look unrelated
to a kernel and would fail the bit-identical coordinate pairing a likelihood
needs (``transformations.md`` §10), so the kind fixes one: telescopes labelled
``i < j < k`` by the array's own order, ``(u1, v1)`` the baseline ``ij`` and
``(u2, v2)`` the baseline ``jk``, the third baseline ``ki`` implied and never
stored. **A kind that has more than one way to represent the same
observation must pick one and say so in its docstring** — that is the
transferable lesson, not the specific ordering.

The chromatic amendment (§6 below) is also a lesson about a kind's shape:
Peter's question of 2026-09-11 found that a *chromatic* model error is sharp
in wavelength and smooth in ``(u, v)``, which the two-axis form
:class:`VisibilitySet` shipped with could not express at all — not
approximately, not by widening a kernel, but literally not, because a kernel
sees a container's axes and nothing else. The fix was a third axis, not a
new kind, and it needed a decision-log entry and a conformance update in the
same PR (ground rule 9). **When a new modality's residual structure needs a
coordinate the kind does not carry, add the axis; do not reach for the
kernel to make up the difference.**

2. The step: what it takes from the observed data, and what it publishes
--------------------------------------------------------------------------

The kind-changing step here is
:class:`~ampere.backends.reference.interferometry.FourierSample`,
``Image -> VisibilitySet``: a direct discrete Fourier transform, no FFT, no
gridding, chosen so its output is exact rather than approximate. Every step
declares ``ACCEPTS``/``PRODUCES`` — the two kinds it converts between, which
is what lets :func:`~ampere.core.negotiate` know an instrument chain type-checks
before anything runs:

.. code-block:: python

    class FourierSample(_Step):
        ACCEPTS: ClassVar[tuple[type, ...]] = (Image,)
        PRODUCES: ClassVar[type] = VisibilitySet

The ``(u, v, wavelength)`` triple a ``FourierSample`` samples at is a
**buffer taken from the observed container** at construction
(``FourierSample.from_observed(container, field_of_view=...)``), never
recomputed and never a parameter. ``architecture.md`` §6's question — "would
you ever put a prior on it?" — answers itself, and taking the coordinates
from the data rather than from the model is what makes
``check_alignment``'s comparison a comparison of two independent things
rather than a tautology (gap I-2). This is the second transferable lesson: a
new modality's kind-changing step reads its own coordinates off the observed
data it is given, by name, at construction — it does not invent them and it
does not ask the model for them.

What the step *publishes* is the other half. A ``FourierSample`` asks for an
image covering ``+-field_of_view/2`` on both axes, sampled at
``max_step = 1/(2 s u_max)`` where ``s`` is an oversampling factor and
``u_max`` the longest baseline observed — the Nyquist limit of the array,
computed once from the data and handed to whichever model is bound to the
channel. An image sampled more coarsely **aliases**: power from beyond the
limit folds back onto the sampled baselines, indistinguishable from real
source structure (gap I-4), so the requirement is enforced as a refusal
(:class:`~ampere.core.exceptions.CompositionError`) rather than a warning,
with an explicit opt-out (``FittingProblem(lenient_compile=True)``) for a
caller who wants the unconfigured model anyway. **A step that changes what a
model must provide says so by publishing a requirement, not by silently
tolerating whatever it is given.**

3. ``configure_from``: one step asking another for its samples
-----------------------------------------------------------------

:class:`~ampere.backends.reference.interferometry.BandwidthSmearing` and
:class:`~ampere.backends.reference.interferometry.TimeSmearing` average a
visibility over several sub-samples of ``(u, v)`` — the array's finite
bandwidth or integration time smears a point into a small patch — and they
need the *preceding* ``FourierSample`` to have computed those sub-samples
before they can average them. ``configure_from`` is the hook: after a chain
is assembled, each step is asked, in order, whether it wants to look at what
comes after it. ``FourierSample`` uses it to learn how many extra samples
per point the smearing step needs and widens its own published pixel scale
accordingly; the smearing step uses it to read the extra samples back.

This is the first **cross-kind** use of ``configure_from`` in ``ampere``
(``LSFConvolution`` already uses it *within* a kind, to learn its output
range from a ``Resample`` after it) — worth citing because it is the
generalisation a new modality is likely to need if its steps come in a
chain longer than one: a step publishing a requirement that only the step
after it can compute is exactly this shape, whatever the kinds involved.

**If your instrument chain is only one step long, there is nothing for
``configure_from`` to do, and that is the right answer, not a gap to fill.**
:doc:`astrometry` §3 found this the hard way: its ``EpochSample`` step is the
whole of its chain, so it does not override ``configure_from`` at all, and
the template as written gave no way to tell whether that silence meant
"correct" or "forgotten". It means correct — a step with no successor has
nothing to read ahead of, and leaving the hook at its no-op default is the
entire implementation a one-step chain needs. State this explicitly when you
build a new modality: a chain that never overrides ``configure_from``, and
says so in its own docstring or page, has made a decision, not skipped a
section of this template.

4. The two-dataset composition: two instruments on one channel
-------------------------------------------------------------------

**Name the shape, because it is one of at least two that ship.** This
section builds **two instruments on one model channel** — visibilities and
closure phases, both reading the ``sky`` channel — which is one composition
shape, not *the* composition shape. :doc:`astrometry` §4 builds the mirror
image, **one model, two channels** (an orbit's ``ra`` and ``dec``, each with
its own instrument and its own dataset), and the two are genuinely
different: this shape's label collision (below) cannot happen in the other
one, because two instruments bound to *different* channels never compete
for a default label. A reader arriving here first should not assume every
multi-dataset composition looks like this one.

This is where the page stops being new and becomes :doc:`sed_composition`
again — literally the same code, with the kind swapped:

.. code-block:: python

    from ampere.core import (Dataset, DatasetCollection, FittingProblem, Instrument,
                             Likelihood, ComplexGaussianFamily, VonMisesFamily,
                             IndependentNoise, VisibilitySet, ClosurePhases, negotiate)
    from ampere.backends.reference import Binary, FourierSample, ClosurePhase

    # 1. One model, one channel "sky".
    model = Binary(x, y, separation=st.uniform(6.0, 20.0), flux_ratio=st.uniform(0.1, 0.8), channels="sky")

    # 2. Two instruments on that channel, distinct labels, coordinates from the data.
    vis_instrument = Instrument(
        [FourierSample.from_observed(observed_vis, field_of_view=40 * u.mas)], channel="sky", label="vis")
    t3_instrument = Instrument(
        [FourierSample.from_observed(observed_t3, field_of_view=40 * u.mas), ClosurePhase()],
        channel="sky", label="t3")

    # 3. Two datasets, two likelihoods -- a complex family on the visibilities, a
    #    wrapped one on the closure phases.
    datasets = DatasetCollection({
        "vis": Dataset(observed_vis, vis_instrument,
                       likelihood=Likelihood(ComplexGaussianFamily(), IndependentNoise()), label="vis"),
        "t3":  Dataset(observed_t3, t3_instrument,
                       likelihood=Likelihood(VonMisesFamily(), IndependentNoise()), label="t3"),
    })

    # 4. One negotiation, one compile, the model built once per draw.
    problem = FittingProblem(model, datasets, seed=20260913)
    problem.requirements["model"]["sky"].sources   # ('vis', 't3')

Both instruments ask ``negotiate`` for an image covering their own field of
view at their own Nyquist step; the union of the two requirements is one
grid, and ``Binary`` is compiled onto it once. Every draw of the sampler
therefore evaluates the binary's image **once**, not twice, and each
instrument reads its own baselines out of that one evaluation — the whole
argument for negotiation (``DEVELOPMENT_PLAN.md`` §4.3), unaffected by which
kind the two instruments happen to emit.

The one new thing an unwary reader meets here that :doc:`sed_composition`
does not have is the label collision that page's own section walks through:
two instruments bound to one channel default their label to the channel
name and collide unless named — ``label="vis"``/``label="t3"`` above is the
fix, not a stylistic choice.

``examples/interferometry/generators.py`` builds this study's own truth — a
binary plus a fainter, partially-resolved circular disc — the same way
:mod:`examples.sed_composition.generators` does: negotiate, ``compile_for``,
evaluate at the truth, push through each instrument, add noise. See that
module's page section "Why synthetic data needs the negotiated grid" for the
failure this dance avoids; it is unchanged here.

5. The likelihoods: a complex Gaussian, a wrapped angle, a circular GP
--------------------------------------------------------------------------

Visibilities are complex, so the family scoring them is
:class:`~ampere.core.ComplexGaussianFamily` — circular complex noise, the
standard interferometric assumption, one real ``sigma`` applying
independently to the real and imaginary parts. Closure phases are angles in
``(-pi, pi]``, so theirs is :class:`~ampere.core.VonMisesFamily`, the wrapped
family for a value with no ordinary residual.

The flexible likelihood — the point of the whole exercise — is
:class:`~ampere.core.GaussianProcessNoise` with an ``axes`` selector, W4.2's
circular complex Gaussian process:

.. code-block:: python

    kernel = Matern32(st.halfnorm(scale=0.1), st.loguniform(2e7, 2e8), axes=("u", "v"))
    noise = GaussianProcessNoise(kernel, DenseGP())

``axes=("u", "v")`` matters because a ``VisibilitySet``'s three axes carry
**mixed units** — dimensionless ``(u, v)`` and a spectral wavelength — and an
isotropic kernel over all three is refused by ``GPSolver.check_compatible``'s
single-unit rule; the selector is what lets a kernel act on a named subset.
The GP is **always on the visibilities, never on the closure phases**, in
every arm of this study: a GP composed with :class:`~ampere.core.ComplexGaussianFamily`
marginalises analytically (one Cholesky of the covariance, no latent block),
while a GP composed with a *wrapped* family like
:class:`~ampere.core.VonMisesFamily` is a **latent** composition — reachable
only on a modern backend under NUTS or VI, never on the reference path with
emcee (``phase4_placement_memo.md`` §7.2 calls this "the blind alley" its own
walk-through found; the item this page documents runs on emcee, so the
choice is structural rather than a preference). Note also the O(N) refusal
this modality forces: a ``(u, v)`` point has no ordered one-dimensional
coordinate whatever a kernel selects, so :class:`~ampere.core.QuasisepGP`
refuses a visibility kernel by name — ``DenseGP`` only, here.

**A GP-flexible likelihood answers one question**: does it keep the
*physical* parameters calibrated when the sky model is wrong, by absorbing
the misspecification's own smooth correlation in the residual, rather than
letting it bias the parameters that are actually of interest? That is M2's
question, asked of this modality — see §7 below for the numbers.

6. The chromatic argument (memo §3.6), measured rather than argued
-----------------------------------------------------------------------

Peter's question: if a model is missing a spectral feature confined to a
patch of sky, is the residual correlated across wavelength, across
``(u, v)``, or does one absorb the other? Write the missing component as
``dI(x, y, lambda) = P(x, y) S(lambda)`` — a patch with a band profile — and
its effect on baseline **B** at wavelength lambda is
``dV(B, lambda) = S(lambda) F(B/lambda)``: a product of two functions on two
*different* coordinates, sharp in wavelength, smooth in spatial frequency.
Along one baseline's spoke, wavelength and ``|u|`` move together, so
projecting onto ``(u, v)`` alone half-works; on a *different* baseline the
same band sits at a different radius, so two points close in ``(u, v)`` can
be at wavelengths a factor of two apart — one in the band, one out of it. An
isotropic ``(u, v)`` kernel is asked to do two incompatible things at once
and cannot; the fix is a kernel that sees wavelength as its own axis, and a
``Product`` of a spatial block and a spectral block:
``k(u, v, lambda) = k_uv(u, v) . k_lambda(lambda)``.

``examples/interferometry``'s chromatic scenario measures the claim rather
than arguing it: an achromatic binary (the fitted model throughout — the
disc of arms (a)-(c) plays no part here) plus a compact, unresolved patch
whose flux follows a narrow Gaussian line in wavelength, observed at several
distinct wavelengths of the *same* baselines (dispersed coverage:
``(u, v) = B/lambda`` differs per channel even though the physical baseline
does not). Three noise models are fitted to the same data: ``Matern32(axes=("u",
"v"))``, ``Matern32(axes=("spectral_axis",))``, and their ``Product`` —
fixing one kernel's amplitude to ``1.0`` in the product, since a product's
marginal variance is the product of its terms' and two free amplitudes
over-parameterise it by one degree of freedom.

**Measured, one seed, the per-PR budget** (``|median - truth|`` in posterior
widths on ``(separation, flux_ratio)``; see the branch report for the
documentation-budget numbers)::

    spatial-only   0.50   2.16
    spectral-only  1.17   5.21
    product        1.97   1.51

The item's own word for this scenario is *informational*, and the honest
report is that the product does not straightforwardly dominate at this
budget — the spatial-only kernel's flux-ratio bias is smaller than the
product's in this particular draw, though the spectral-only kernel is
clearly the worst of the three throughout. Two things are worth separating
here: whether a kernel finds the *right structure* (does its own
length-scale posterior land near the injected line width?) and whether that
translates into the smallest bias on two parameters from a 168-point, single
draw — the second is a noisier question than the first at this budget, and
this page reports the numbers rather than smoothing them into a cleaner
story than the data support. ``tests/interferometry/test_chromatic.py``
asserts only that all three kernels compose and produce a finite posterior;
``pytest -m interferometry_full`` re-runs the three at the documentation
budget for a less noisy comparison.

7. What the conformance suite owes a new modality
-----------------------------------------------------

``tests/conformance/test_interferometry.py`` is the checklist this modality
paid, and it is the transferable list for the next one:

- **a row per step against a closed form** — the direct DFT of a band-limited
  image against :class:`~ampere.backends.reference.interferometry.GaussianSourceVisibilities`'s
  analytic visibility, the closure phase of a binary against its closed
  form and its sign convention;
- **a mask-propagation row** — a masked input baseline masks every closure
  triangle it takes part in (the many-to-one rule,
  ``results_schema.md`` §16);
- **a two-dataset row** — two instruments on one channel, both recorded in
  ``problem.requirements``, the model evaluated once per draw;
- **a draw row** — ``FittingProblem.simulate`` producing both kinds from one
  model, and ``simulate_many(native=True)`` agreeing with the loop on the
  modern backends.

A new modality that cannot state its own version of these four rows has not
yet demonstrated that its steps, its kind and its composition behave the way
every other modality's do.

8. The native twins: the inheriting pattern, and why
--------------------------------------------------------

Two patterns exist in this codebase for associating a backend's step with
its reference counterpart — there is no step registry
(``phase4_placement_memo.md`` §1.2), so association is by **the same class
name in the backend's own module**, plus the native
``apply_flux``/``native_flux`` surface and the ``BACKEND`` declaration.
``ampere.backends.torch.instrument`` **re-declares** every step and is held
to the reference backend by the conformance battery; the interferometry
steps and models instead **inherit** their reference counterpart on both
modern backends and override only the four capability flags and the
arithmetic:

.. code-block:: python

    class GaussianSource(_JaxImageModel, _ReferenceGaussianSource):
        ...  # only _brightness_native is new; compile_for, negotiation, from_observed are inherited

The reason is specific to this modality: here the *declaration* is
substantial (the Nyquist requirement, the expanded-coverage bookkeeping
:class:`FourierSample` does for the smearing steps, ``compile_for``'s grid
adoption) and the arithmetic is comparatively small, the reverse of an
instrument step like an LSF convolution. Writing the declaration twice would
risk two backends negotiating *slightly* different pixel scales for the same
requirement — invisible until an image sampled just below Nyquist aliases,
which looks like real source structure rather than a bug (gap I-4). **When a
new modality's per-backend declaration is the expensive, dangerous half and
the arithmetic the cheap half, inherit; when the declaration is thin and the
arithmetic is the whole of the step, re-declare** (``ampere.backends.torch.instrument``'s
own reasoning, for the opposite case).

That rule has two clauses and :doc:`astrometry` needed a third. Its
``EpochSample`` has no arithmetic at all (``apply`` is the identity) *and*
its declaration is a few lines, so neither clause fires — both halves are
cheap, and the rule as written does not say what to do. **When both the
declaration and the arithmetic are cheap, inherit anyway, for the
one-pattern-per-codebase argument, not because duplicating either half would
be dangerous**: there is nothing to gain from a second, hand-maintained copy
of a two-channel template that must stay bit-identical across backends, and
a codebase with one twinning convention is easier to read than one with two
conventions and a rule for choosing between them that only sometimes
applies.

``examples/interferometry/model.py``'s ``BinaryWithDisc`` — the "correct"
sky model of §9 below — takes the same lesson one step further: rather than
adding a fourth subclass to ``ampere.backends.*``, it is a small,
hand-written :class:`~ampere.core.transform.Model` that **delegates** to one
``Binary`` and one ``GaussianSource`` instance and adds their brightness,
following ``native_grid``/``native_flux`` through unchanged. It composes
under NUTS on every backend that declares those two classes, and it needed
no change to ``ampere`` at all — the transferable claim being that a user's
own model can be built this way, out of existing pieces, without touching
the library.

9. The study: does the flexible likelihood survive an unmodelled disc?
---------------------------------------------------------------------------

``examples/interferometry/generators.py`` builds one truth: a resolved
binary (the same synthetic source ``tests/backends/interferometry_fixtures.py``
builds, reused rather than re-derived — see that module for the array
geometry) plus a fainter, partially-resolved circular
:class:`~ampere.backends.reference.interferometry.GaussianSource` disc,
sized against the array's own resolution so it neither sits invisibly under
the binary's own signal nor resolves out to nothing on every baseline
observed. It is observed once, as both visibilities and closure phases, and
fitted three ways:

============  ==============================  ================
Arm           Model                           Visibility noise
============  ==============================  ================
correct       binary + disc (disc fixed)       independent
incomplete    binary alone (disc omitted)      independent
flexible      binary alone (disc omitted)      Matern-3/2 GP over (u, v)
============  ==============================  ================

The closure phases stay under independent von Mises noise throughout — the
flagship GP is always on the visibilities (§5).

**The claim is about coverage, not about one draw's luck**, so
``examples/interferometry/study.run_calibration`` asks it through simulation-based
calibration (:func:`~ampere.results.calibration.sbc`, in
``tests/m2/test_visibility_calibration.py``'s shape): draw a fresh binary
from the fitted models' own prior, simulate a fresh two-dataset observation
from the **correct** truth, refit under one arm's formulation, repeat twelve
times. Measured at the pinned seed, central-90 % coverage on
``(separation, flux_ratio)``:

.. code-block:: text

    correct      [1.00, 1.00]   -- calibrated, as a well-specified fit should be
    incomplete   [0.25, 0.00]   -- confidently wrong: the omitted disc biases both parameters
    flexible     [0.92, 0.92]   -- recovered: the GP absorbed the disc's smooth (u, v) excess

``tests/interferometry/test_calibration.py`` pins the direction and a margin
(twelve simulations is a smoke budget in Talts et al.'s sense — evidence of a
gross effect, not a fine one), not these exact figures. The whole run — three
arms, twelve simulations each — takes about two and a half minutes on the
reference backend.

Three engines, as the item asked
-----------------------------------

``examples/interferometry/study.run`` dispatches on the problem's own
backend, exactly as :mod:`examples.m2_misspecification.study` does: emcee on
the reference backend (every number above), NUTS on torch through
:class:`~ampere.inference.NUTSEngine` (``tests/interferometry/test_engines.py``
fits both the "correct" and "flexible" arms at a smoke budget — the composite
model's ``native_grid``/``native_flux`` surface, exercised), and
:class:`~ampere.inference.SBIEngine` for neural posterior estimation on the
"flexible" arm, with ``layout="set"`` (a complex, two-dataset container needs
it — W3.3's coordinate-value-mask encoding is the first thing that reaches a
visibility). Its calibration row quotes **TARP and the coverage curve**
rather than the marginal rank histograms' KS p-value, for
``tests/inference/test_interferometry.py``'s reason: a posterior correct in
each parameter's own margin and wrong about their correlation passes a
marginal test and fails TARP, which is the failure mode that matters for a
two-parameter fit whose parameters are exactly as correlated as a binary's
separation and flux ratio are.

The plots: a found limitation, not assumed
-----------------------------------------------

Four of the "six shipped plots" —
:func:`~ampere.results.plot_posterior_predictive`,
:func:`~ampere.results.plot_residuals`,
:func:`~ampere.results.plot_gp_localisation` and
:func:`~ampere.results.plot_anomaly_score` — read a stored group through
:func:`ampere.results._plotting.coordinate_of`, which needs **exactly one
ordered coordinate axis** and refuses a "point kind with several axes" by
name, citing ``results.md`` §4 and ``DEVELOPMENT_PLAN.md`` §4.4/§4.8's own
staging: *"Gridded and multi-axis kinds are the Phase 5 staging."* Both of
this study's kinds are exactly that — a :class:`~ampere.core.VisibilitySet`
sample is a joint ``(u, v, spectral_axis)`` point with no natural order, and
a :class:`~ampere.core.ClosurePhases` sample five such coordinates — so all
four of those renderers refuse both datasets, on every arm, today. This is
not the complex-valued gap W4.2's own carried note anticipated ("W4.4's
figures must pick a component or modulus of the complex conditional mean")
— that note undersold it: even a real, single-component view of the same
data still has no ordered axis to plot against, because the refusal is about
the *coordinate*, not the *value type*.

What this leaves, and what ``examples/interferometry/figures.py`` renders,
is :func:`~ampere.results.plot_corner` and :func:`~ampere.results.plot_trace`
(posterior-only — neither needs a data coordinate) for every arm, plus
:func:`~ampere.results.plot_sbc_ranks` and :func:`~ampere.results.plot_coverage`
for the calibration study of the previous section, which live in *parameter*
space and do not hit the same wall. ``python -m examples.interferometry
--figures DIR`` (add ``--calibration`` for the coverage figures) writes
them; nothing is committed (ground rule 7).

**The transferable lesson for the next modality**: if your kind is a
multi-axis point set (an interferometric visibility, a set of astrometric
positions), budget for exactly this gap before promising "the six plots" —
``ampere.results``'s family B/C diagnostics are staged for Phase 5 and are
not yet a modality-agnostic surface. Checking ``ampere/results/`` for this
restriction before relying on any of the four above is worth the two minutes
it costs.

A hazard the models here never raise: periodic models and prior width
-----------------------------------------------------------------------------

Every source model this page fits — a Gaussian, a uniform disc, a binary —
has a smooth, unimodal likelihood in its own parameters over any prior wide
enough to be "uninformative". Nothing here says what to do about a model
that is not, and :doc:`astrometry` §6 found one: a reflex orbit is
**periodic**, and a period prior wide enough to reach past the observed
epochs' own baseline is genuinely multi-modal — an emcee ensemble and a
multi-chain NUTS run both locked onto a spurious period and stayed there,
not because of slow mixing but because a chain that started near a real
second mode has no reason to leave it. **A modality whose model is periodic
(or otherwise genuinely multi-modal in its own right, independent of the
likelihood or noise machinery) needs its own discussion of prior width,
separate from the composition questions the rest of this page is about**: an
informed prior centred on a value a periodogram or a previous epoch has
already suggested, not a search over decades hoping the sampler finds the
right one. Check whether your new modality's model has this shape before
trusting a wide "uninformative" prior to behave the way it does for every
model on this page.

See also
------------

* :doc:`sed_composition` — the simple case this page cites throughout: one
  model, two instruments, one channel, no kind change and no complex
  container.
* :doc:`astrometry` — the second worked modality, built by following this
  page section by section; its closing section records where this page
  needed fixing, which is what W4.8 did.
* :doc:`m2_misspecification` — the flagship study this page's arms
  reproduce the *shape* of on a different observable: a correctly and an
  incorrectly specified model, with and without the flexible likelihood.
* ``docs/design/phase4_placement_memo.md`` §3.6-3.7 and §7.2 — the chromatic
  argument in full, and the worked composition this page's §4 is drawn from.
