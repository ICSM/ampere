Adding a modality: images, and what a grid changes
====================================================

Phase 5's gridded customer (W5.5), and the third modality built by following
:doc:`interferometry`'s template section by section. Its purpose is narrower
than either of the two before it and can be stated in one sentence: **no image
had ever been a dataset**. :class:`~ampere.core.Image` has been a kind since
W1.4, Phase 4's three source models all emit one, and every one of them handed
it straight to a Fourier step — so the kind had only ever been an
*intermediate*, never an observation. This page is the first fit whose observed
container has ``Layout.GRID``, and its closing section reports what the
template did not say about that.

Read :doc:`interferometry` first; the section order below is its, so that a gap
shows up as a gap rather than as a rewritten structure. :doc:`astrometry` is
the other page written this way, and it is worth reading both: astrometry's
finding was that a modality can be *easy*, and this one's is that a grid is a
different shape of hard from either.

The runnable example is ``examples/image``
(:download:`generators.py <../../examples/image/generators.py>`,
:download:`model.py <../../examples/image/model.py>`,
:download:`grid_gp.py <../../examples/image/grid_gp.py>`,
:download:`study.py <../../examples/image/study.py>`,
:download:`figures.py <../../examples/image/figures.py>`,
:download:`__main__.py <../../examples/image/__main__.py>`): a compact Gaussian
source on a smooth, much broader background, observed through a Gaussian PSF,
fitted three ways to ask M2's question of a gridded observable. Run it
yourself::

    python -m examples.image                    # 3 arms, reference, emcee
    python -m examples.image --calibration      # the SBC coverage row
    python -m examples.image --benchmark        # DenseGP vs HSGP at three N
    python -m examples.image --backend torch    # NUTS

``tests/conformance/test_image.py`` and
``tests/examples/test_image_study.py`` are this page's own coverage.

1. The kind: nothing to add, and one attribute that changes everything
------------------------------------------------------------------------

Like :doc:`astrometry`'s, this modality needed no new kind:
:class:`~ampere.core.Image` already declared everything it needed —

.. code-block:: python

    class Image(FunctionSamples):
        AXES = (
            AxisSpec("x", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
            AxisSpec("y", physical_types=("angle", "length"), order=Order.STRICTLY_MONOTONIC),
        )
        LAYOUT = Layout.GRID

— and it did not need amending. But one of those three attributes is different
from anything the template's earlier modalities used, and it is the whole
subject of this page: ``LAYOUT = Layout.GRID``.

A ``Layout.POINTS`` container's axes **are** its samples: a spectrum with 400
coordinates has 400 values, and stacking its axes column-wise gives the
``(N, n_axes)`` matrix everything downstream wants. A ``Layout.GRID``
container's axes are the separable grids its samples are the *product* of: a
64x64 image has two axes of 64 coordinates and 4096 values. Everywhere in the
contracts that "the axes" stood in for "the samples" — and there turned out to
be three such places, two of which nobody had noticed — a grid is a different
computation.

**The transferable lesson**, and it is the one thing to take from this page if
you take nothing else: ``LAYOUT`` is not documentation. Adding a modality whose
``LAYOUT`` differs from every shipped modality's is a bigger change than adding
one with a new ``AXES`` tuple, because ``AXES`` is data the generic machinery
reads and ``LAYOUT`` is a branch it has to have.

2. The step: what it takes from the observed data, and what it publishes
--------------------------------------------------------------------------

:class:`~ampere.backends.reference.image.PSFConvolution` is ``Image -> Image``,
kind-preserving (``PRODUCES`` left at ``None``) — the "PSF convolution" row
``transformations.md`` §10's second table named at the freeze and nothing had
filled. It lives in a new ``image.py`` per backend rather than beside the
spectral steps in ``instrument.py``, which is D1's "one module per observable"
applied as the §10 amendment states it. (The image *models* it convolves stay
in ``interferometry.py``, where W4.1 put them: they are the same objects,
reused unchanged.)

Like every other step in the table, it takes its coordinates **from the
observed container**::

    step = PSFConvolution.from_observed(observed_image, fwhm=2.5)

The PSF itself comes two ways, and the difference between them is the second
rule of §10 in its gridded form. ``kernel=`` is a pixel array tabulated at the
observed pixel scale — a measured PSF, a buffer tied to the coordinates it was
measured on — and the step **refuses** a grid sampled at a different scale
rather than interpolating a calibration product into one it is not. ``fwhm=``
is a circular Gaussian, analytic, so it is simply retabulated on whatever grid
arrives; and because ``apply`` builds it through
:meth:`~ampere.core.Parameterised.context` every evaluation,
``promote_buffer("fwhm", prior=...)`` turns the seeing into an ordinary fitted
parameter with no change to the class, exactly as it does for an LSF's width.

What it publishes is where a convolution differs from everything before it, and
it is worth being precise about why. A convolution's output pixel is a weighted
sum over its neighbours, so the step needs the model evaluated **beyond** the
observed field — the same problem
:class:`~ampere.backends.reference.instrument.LSFConvolution` has, in two
dimensions. The LSF's answer is to publish a padded *interval* and leave a
resampler downstream to land on the observed grid. A PSF has no resampler
downstream; there is no spatial one. So it publishes both halves of one
constraint:

* ``points=`` at an evenly spaced grid extending the observed pixel centres by
  the kernel's half-support on each side — the coverage the convolution needs,
  **and** the exact coordinates :meth:`apply` will look up with
  :meth:`~ampere.core.Axis.locate`;
* ``intervals``/``max_step`` restating the observed pixel scale as a *density*,
  which is what makes a model holding its own grid refuse by name when that
  grid is too coarse to represent what was observed.

Then ``apply`` convolves on the padded grid and **crops** to the observed
pixels, in that order — cropping first and convolving after is the boundary
error the padding exists to prevent.

Publishing one constraint two ways has a cost that is easy to miss: the two
halves have to agree to the last bit, or the union negotiation builds from them
contains near-duplicate coordinates and is evenly spaced nowhere — at which
point the FFT route refuses, on whichever backend happened to run, for a reason
that has nothing to do with that backend. The step builds its padded grid the
way :meth:`~ampere.core.AxisRequirement.coordinates` builds one, from a
``linspace`` over the same interval, so they agree by construction rather than
by luck, and a conformance row asserts the negotiated grid is still regular.
**If your new modality's step publishes the same constraint in two forms, make
one of them compute the other.**

3. ``configure_from``: not needed, for a reason worth stating
---------------------------------------------------------------

An LSF convolution has a ``configure_from`` hook because it genuinely does not
know its own output range — a resampler or a photometry step downstream sets
it, and gap I-3's mechanism is how it finds out. A PSF convolution does know:
its output range is the observed image, which it was built from. So this step
has no hook, and the absence is the point. :doc:`astrometry` §3 reached the
same conclusion for the same reason, and between them the rule is now clear:
**a step needs ``configure_from`` when something downstream of it decides its
output range, and not otherwise.**

4. The composition: one channel, one dataset, and a model built by addition
------------------------------------------------------------------------------

The simplest composition of the three worked modalities: one sky channel, one
instrument, one dataset. What is not simple is the *model*, and the way the
example builds it is the transferable part.

``ampere`` ships no "add two image models together" operator — there is nothing
analogous to :class:`~ampere.core.kernels.Sum` for a
:class:`~ampere.core.transform.Model` — and ``examples/image`` does not add one.
It writes a small model of its own that **delegates** to two shipped
:class:`GaussianSource` instances and adds their brightness, following
``examples/interferometry/model.py``'s ``BinaryWithDisc`` exactly. Neither the
negotiation machinery nor the native ``native_grid``/``native_flux`` surface is
reimplemented; both are the children's, threaded through. It needed no change
to ``ampere`` at all.

One detail in that class is a finding rather than a style choice. The
misspecified arms are built with ``background_flux=None``, which **omits the
background model entirely** rather than setting its flux to zero. A zero-flux
Gaussian would still contribute an array of zeros to every evaluation and a
fixed parameter to every provenance record, so the misspecified arm would differ
from a genuinely background-free model in its spec hash while agreeing in its
numbers — two runs that are the same fit and hash differently. **A
misspecification is a model that does not have the component, not a model that
has it at zero.**

5. The likelihood: where a grid is currently blocked
--------------------------------------------------------

This is the section that does not yet read like the two pages before it, and
the honest thing is to say so plainly.

``ampere.core`` refuses a correlated noise model on **any** ``Layout.GRID``
container, in two places:
:meth:`~ampere.core.GPSolver.check_compatible` ("the v1 GP solvers work on
point-set containers […]; gridded 2D+ data are the subject of the SVGP / SKI /
Vecchia strategy slots") and ``Likelihood._coordinates`` ("a correlated noise
model needs point-set coordinates"). Both messages name the Phase 5 slot they
are waiting for, and both predate an image ever being an observation.

Neither refusal is about mathematics. A stationary kernel over ``(x, y)`` is a
function of coordinates and a grid has coordinates — it just keeps them
separably, so the ``(N, 2)`` matrix a kernel wants has to be broadcast out of
the axes rather than stacked from them. That broadcast has existed and been
correct since W3.3, in ``ampere.core.encoding``'s coordinate matrix for the SBI
encoder; W5.5 lifted the same rule into :func:`~ampere.core.sample_coordinates`
so that :meth:`Dataset.draw_observation` could draw an image at all.

So ``examples/image/grid_gp.py`` lifts the two gates **out of tree**, from the
public API, in the way ``transformations.md`` §11 exists to prove is possible:
three small subclasses that override only the layout check and leave every
other rule — the kernel's axis selection, the quasiseparable refusals, the
ordered-1D rule — to the base classes. Read that module before the study; it is
a demonstration that nothing but the gate is missing, not a design, and it
carries the two-edit library change it stands in for.

With the gates lifted, a two-axis ``Matern32(axes=("x", "y"))`` composes over
an image and scores exactly as it does over a spectrum.

6. The study: does the flexible likelihood survive an unmodelled background?
-------------------------------------------------------------------------------

One truth — a compact Gaussian source on a smooth, much broader background —
observed once through the PSF, fitted three ways: ``correct`` (the background
in the model at its true value), ``incomplete`` (the background absent), and
``flexible`` (the background absent, and a two-axis Matérn-3/2 GP over
``(x, y)`` on the residuals).

The background is chosen to be smooth **in the observed coordinates**
deliberately, and that is the chromatic lesson of :doc:`interferometry` §1
restated for two spatial axes: a kernel sees a container's axes and nothing
else, so a component the flexible likelihood is meant to absorb has to be
smooth in the coordinates the container actually carries. A periodic detector
fringe would not be, and no widening of the kernel would make it so.

``python -m examples.image --calibration`` runs the SBC row;
``tests/examples/test_image_study.py`` pins it. Twelve simulations is a smoke
budget in Talts et al.'s sense — evidence of a gross effect, not a fine one —
and what is pinned is the *direction*, not a figure.

7. The benchmark: the phase's own "chosen by measurement"
-------------------------------------------------------------

``DEVELOPMENT_PLAN.md`` §5's rule for Phase 5 is that a solver is chosen by
measurement, and this modality is where the measurement happens, because a
gridded dataset is the first one where N is large enough for the choice to
matter. ``python -m examples.image --benchmark`` scores the *same* flexible
likelihood on the *same* image under the exact :class:`~ampere.core.DenseGP`
and under W5.4's approximate
:class:`~ampere.core.HilbertSpaceGP` with a tensor-product basis, at three
image sizes.

Note what neither solver is doing: **neither exploits the grid.** Both treat an
image as N scattered points that happen to lie on a lattice. A Kronecker or SKI
solver would use the lattice, and that is W5.6's bake-off; these two are the
honest baseline its candidates have to beat.

8. What the conformance suite owes a gridded modality
---------------------------------------------------------

``tests/conformance/test_image.py`` is the checklist this modality paid, in the
template's §7 shape, with one substitution worth explaining:

- **a row per step against a closed form** — a PSF convolution has none. There
  is no Bessel function to compare against, because the answer depends on the
  source, so the oracle is the **definition**: a direct sum over the kernel's
  support, written as explicit index arithmetic that shares nothing with an FFT
  except the answer. **When a new modality's step has no closed form, the
  definition is the oracle, and it must be written so that it could not agree
  by construction** — no calling the same helper, no borrowing the same
  convolution routine;
- **a mask-propagation row** — one masked model pixel masks exactly its
  ``(2k+1)²`` neighbourhood, cropped to the observed field;
- **a composition row** — an ``Image`` on both sides of a likelihood, the
  channel recorded, the density peaking at the truth;
- **a draw row** — ``FittingProblem.simulate`` producing an ``Image``
  observation.

Plus this modality's own: the requirement's padding equals the kernel's
half-support, the negotiated union is still evenly spaced, and three refusals
(a model holding a coarser grid, an irregular observed axis, a tabulated kernel
handed a different pixel scale).

9. The native twins: the inheriting pattern, first clause
--------------------------------------------------------------

:doc:`interferometry` §8's rule has three clauses, and this step fires the
first cleanly, which is worth recording because :doc:`astrometry` needed a new
one. Here the declaration is the expensive **and** the dangerous half — one
constraint published in two forms that have to agree to the last bit — and the
arithmetic is one transform and one crop. So ``ampere.backends.torch.image``
and ``ampere.backends.jax.image`` derive from the reference class and override
only the four capability flags and the sums.

Two things are shared rather than rewritten, and the reason generalises:
``grid_steps()`` (the regularity and pixel-scale refusals) and
``crop_indices()`` (the ``Axis.locate`` lookup) were factored **onto the
reference class** so that all three backends say the same words about the same
mistake and find the same pixels. **Anything that is a rule rather than
arithmetic belongs on the shared class, even when it lives inside a method the
twins override.**

One thing had to be reproduced rather than shared, and it is a small trap:
``scipy.fft.next_fast_len`` is neither available in torch nor traceable in jax,
and the padded transform length decides how many implicit zeros the transform
sees. Two backends padding differently would show up as a cross-backend
disagreement with no cause visible in either, so each twin reimplements the same
5-smooth rule and says why. **Check whether your transform has a
performance-only parameter that is silently also a numerical one.**

10. What the template did not say about a gridded kind
---------------------------------------------------------

Four findings, in the order a reader building the next gridded modality would
hit them.

**(a) ``Layout`` is a branch, not a label — and three places had no branch.**
The template's §1 asks a reader to weigh three class attributes and treats all
three as declarations the generic machinery reads. ``LAYOUT`` is not like that:
``Layout.GRID`` means "the axes are not the samples", and every piece of
machinery that flattened a container's axes into an ``(N, n_axes)`` matrix had
been written as though they were. Three places did it, and only one of them —
the SBI encoder, W3.3 — had the grid branch. The other two were
``Dataset.draw_observation`` (which raised an ``IndexError`` about a boolean
mask, saying nothing about layouts) and ``Likelihood._coordinates`` (which
refuses by name, correctly, but as a deferred slot rather than as a gap). **A
new modality with an unfamiliar ``LAYOUT`` should grep for every place the
container's axes are stacked before writing anything else.**

**(b) A mask rule written for a point set does not survive a grid, and the fix
is not to generalise the matrix.** ``transformations.md`` §13.5 said a
``Layout.GRID`` container's mask "must be flattened by the transformation
itself", which is honest but wrong-shaped for a convolution: the influence
matrix for a 64x64 image against a 9x9 kernel has sixteen million entries, 81
of them non-zero per row. The rule that matters ("an output sample touching any
masked input is masked") is unchanged; its *expression* is a separable
dilation, ``O(k N)``, which is what :func:`~ampere.core.propagate_mask_grid`
computes. **A generic rule stated in terms of one data structure may need a
second expression, not a second rule.**

**(c) The plotting wall is a different wall on a grid.**
:doc:`interferometry` §9's closing lesson is that four of the six shipped plots
refuse a multi-axis point kind, and it tells the next modality to budget for
that. An image hits the same refusal for a different reason, and the difference
changes what fixing it would mean. Interferometry's kinds have no natural order,
so there is genuinely nothing to plot a residual *against*. An image's two axes
**are** ordered; what is wrong is the renderer's output *shape* — a
one-dimensional figure where the natural picture is a two-dimensional panel. So
widening ``coordinate_of`` to accept two axes would unblock neither until
someone decides what a residual picture is. ``examples/image/figures.py`` draws
those panels by hand to show what is missing. **Two different gaps can hide
behind one refusal message; check which one you have before assuming a shared
fix.**

**(d) A gridded dataset is where N stops being small, and the contracts had
noticed only half of it.** Every modality before this one has hundreds to
thousands of samples; a modest image has tens of thousands. W5.4's
:class:`~ampere.core.HilbertSpaceGP` was written for exactly that and says so
("the scaling answer in two and three axes"), and yet the gate in front of it
was still shut, because the two were specified in different items and nothing
ran an image through a GP until this one did. **When two items are each written
to meet the other halfway, something has to actually run the join.**

See also
------------

* :doc:`interferometry` — the template this page follows, section by section.
* :doc:`astrometry` — the other page written that way; its closing section
  reports a different set of template gaps.
* :doc:`m2_misspecification` — the flagship study whose shape this page's three
  arms reproduce on a gridded observable.
* :doc:`kernels` — what a two-axis kernel's ``axes=`` selection means, and why
  a kernel sees only the container's own coordinates.
