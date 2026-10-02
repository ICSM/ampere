The kernel algebra
====================

:doc:`concept` introduces the flexible likelihood and its default kernel in
a paragraph; this page is the full reference for :mod:`ampere.core.kernels`
(landed at W4.5, extended at W5.7): the quasiseparable families and their ranks,
``Sum``/``Product``/``SpectralMixture``, the ``axes=`` selector and the
per-leaf unit rule, and ``register_quasiseparable_term`` for reaching the
O(N) path from outside ``ampere``. ``tests/core/test_kernels.py`` is this
page's own coverage, and every code block below is drawn from a test or a
docstring that already runs there.

1. The families, and their ranks
------------------------------------------

A kernel is a covariance function plus, optionally, an exact representation
on the O(N) path. :func:`~ampere.core.quasiseparable_families` lists those
that have one:

.. code-block:: pycon

    >>> from ampere.core import quasiseparable_families
    >>> quasiseparable_families()
    ('matern12', 'matern32', 'matern52', 'rotation', 'sho', 'spectral_mixture',
     'sum', 'warped')

Each is a semiseparable matrix of a fixed **rank** — the size of the
generator vectors celerite2's solver factorises the covariance into, which
is what fixes the constant in front of the O(N) solve:

.. list-table::
   :header-rows: 1
   :widths: 22 10 68

   * - Kernel
     - Rank
     - Notes
   * - :class:`~ampere.core.Matern12`
     - 1
     - The Ornstein-Uhlenbeck / exponential kernel — exact, and the cheapest
       quasiseparable term there is.
   * - :class:`~ampere.core.Matern32`
     - 2
     - The default kernel (:doc:`concept`): once differentiable, a better
       description of a real model deficiency than infinite smoothness, and
       exactly quasiseparable.
   * - :class:`~ampere.core.Matern52`
     - 3
     - **Exact**, not approximate — there is no exact Matérn-5/2 in the
       celerite *basis* (``e^{-c*tau}(a cos d*tau + b sin d*tau)``), but the
       solver underneath factorises any rank-*J* **semiseparable** matrix,
       and a degree-2 polynomial in :math:`t_n - t_m` is a rank-3 bilinear
       form in :math:`(1, t, t^2)`. ``likelihoods.md`` §15.3's older claim
       that Matérn-5/2 has no exact quasiseparable form is withdrawn.
   * - :class:`~ampere.core.SHO`
     - 2
     - The damped simple-harmonic-oscillator term (celerite's ``SHOTerm``,
       ``Q > 1/2``) — the right shape for a residual with a *period*:
       interference fringing, an instrumental ripple.
   * - :class:`~ampere.core.RotationTerm`
     - 4
     - A pair of ``SHO`` terms, celerite2's own construction for a
       non-sinusoidal periodic residual.
   * - :class:`~ampere.core.SpectralMixture`
     - :math:`2k`
     - A sum of :math:`k` ``SHO`` terms with free frequencies — the
       flexible-likelihood analogue of a Gaussian-mixture spectral density.
   * - ``Sum``
     - :math:`\sum_i \text{rank}_i`
     - Not a kernel family in its own right; registered so that a **sum of
       quasiseparable terms** reaches the O(N) path by concatenating its
       children's generators (§2 below).
   * - :class:`~ampere.core.WarpedKernel`
     - as its base
     - **W5.7**, and likewise not a family in its own right: a monotone input
       warp and a diagonal amplitude warp around any base kernel, so the
       flexible likelihood stops being stationary without leaving the O(N)
       path. ``likelihoods.md`` §6 is the full account.

A ``Sum`` of noise terms should carry the **sparsity prior** that goes with it,
so that a component the data do not need is switched off by the prior rather
than fitted to whatever is left over. That prior is
:func:`~ampere.core.shrinkage_horseshoe` — one global scale shared by every
component, one local scale per component under it, each amplitude under its own
local scale — and :func:`~ampere.core.with_shrinkage` puts it on a kernel
without changing anything else about it::

    kernel = with_shrinkage(
        Sum(Matern32(...), Matern32(...), labels=("broad", "narrow")),
        shrinkage_horseshoe(("broad.amplitude", "narrow.amplitude")),
    )

It is the recommended prior for any sum of noise components (**W5.8**);
``likelihoods.md`` §6 gives the account, including which part of Piironen &
Vehtari's regularised horseshoe is declarable today and which is not, and
``examples/m2_misspecification/many_lines.py`` is the demonstration.

:class:`~ampere.core.SquaredExponential` ships too, retained deliberately
(``likelihoods.md`` §14) as the point of comparison M2 needs, but it
declares ``QUASISEPARABLE = False``: infinitely smooth, so no polynomial
generator represents it exactly, and it composes on :class:`~ampere.core.DenseGP`
only — or, since **W5.4**, on :class:`~ampere.core.HilbertSpaceGP`.

Every stationary family above has a second representation too: its **power
spectral density**, :meth:`~ampere.core.Kernel.spectral_density`, which is
what :class:`~ampere.core.HilbertSpaceGP` builds a reduced-rank
approximation from. The three Matérns, the squared exponential, the ``SHO``
and any ``Sum`` of those have one in closed form — and so, for free, does
``SpectralMixture``, which *is* a ``Sum`` of ``SHO`` terms here; a ``Product``
and a ``RotationTerm`` do not, and say so by name at composition. It is the one route that makes a squared exponential scale —
its spectral density is a Gaussian, so the approximation converges
exponentially in the basis size, which is the opposite of the trade the O(N)
path offers.

.. code-block:: pycon

    >>> from ampere.core import HilbertSpaceGP, Matern32, SquaredExponential
    >>> float(Matern32(0.4, 2.0).spectral_density(0.0, {"amplitude": 0.4,
    ...                                                 "length_scale": 2.0}))
    0.7390...
    >>> HilbertSpaceGP(basis_size=64).EXACT
    False


A default is not a restriction. Every family above is an ordinary
:class:`~ampere.core.Kernel` with priors on its own hyperparameters, so
picking one (or several, see §2) is exactly as easy as picking the default:

.. code-block:: python

    from ampere.core import (
        GaussianProcessNoise, HilbertSpaceGP, Matern32, QuasisepGP, SHO,
    )

    kernel = Matern32(amplitude=st.halfnorm(scale=0.1),
                       length_scale=st.loguniform(0.5, 10.0))
    noise = GaussianProcessNoise(kernel, QuasisepGP())      # exact, O(N)

    # or, in two or three axes, or for a kernel with no quasiseparable form:
    reduced = GaussianProcessNoise(kernel, HilbertSpaceGP(basis_size=64))

2. ``Sum``, ``Product`` and ``SpectralMixture``
---------------------------------------------------

Two composite kernels, and one composite that is really a sum in disguise:

.. code-block:: pycon

    >>> from ampere.core import Matern12, Matern32, SquaredExponential, Sum, Product
    >>> Sum(Matern32(0.3, 2.0), Matern12(0.1, 0.2)).QUASISEPARABLE
    True
    >>> Sum(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)).QUASISEPARABLE
    False
    >>> Product(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)).QUASISEPARABLE
    False

**A sum of quasiseparable terms is quasiseparable.** ``Sum``'s own
``QUASISEPARABLE`` flag is computed from its children rather than declared —
``all(term.QUASISEPARABLE for term in terms)`` — so a broad Matérn plus a
narrow oscillator (M2's ``fringing`` scenario, refitted at W4.5 with
``Matern32 + SHO``) still costs O(N): the generators concatenate and the
ranks add (§1's table), which is what :func:`~ampere.core.kernels.sum_representation`
does mechanically.

**A product of quasiseparable kernels is not, in general, quasiseparable**,
and ``Product`` says so by declaring ``QUASISEPARABLE = False`` outright
rather than trying and failing per instance — a kernel product is a
Hadamard product of two semiseparable matrices, which is semiseparable only
in special cases celerite2's solver does not exploit. ``Product`` therefore
composes on :class:`~ampere.core.DenseGP` only, and the refusal on the O(N)
path names the reason rather than merely "not registered":

.. code-block:: pycon

    >>> import numpy as np, astropy.units as u
    >>> from ampere.core import GaussianFamily, GaussianProcessNoise, QuasisepGP, Spectrum
    >>> data = Spectrum(np.linspace(1.0, 10.0, 8) * u.um, np.zeros(8) * u.Jy,
    ...                 uncertainty=np.full(8, 0.1) * u.Jy)
    >>> noise = GaussianProcessNoise(Product(Matern32(0.3, 2.0), SquaredExponential(0.1, 0.2)), QuasisepGP())
    >>> noise.check_compatible(GaussianFamily(), data)
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: QuasisepGP cannot lower a Product: a product
    of quasiseparable kernels is not quasiseparable. ... Use DenseGP, or replace the
    Product with a Sum, which is quasiseparable exactly when every term is.

Where the two factors act on **disjoint axes** — the chromatic case §4
below is built from — the product still has no O(N) form, but it composes
cleanly on the dense path, which is the case that matters: :doc:`interferometry`
§6 is a ``Product`` of a spatial and a spectral kernel fitted on the dense
solver throughout.

``SpectralMixture`` is not a third composition rule; it *is* ``Sum``, of
``SHO`` terms with free frequencies — one amplitude, one period and one
quality per component, the component count the only real decision:

.. code-block:: pycon

    >>> from ampere.core import SpectralMixture
    >>> mixture = SpectralMixture([0.2, 0.1], [1.0, 0.3], [4.0, 8.0])
    >>> mixture.spec().family
    'spectral_mixture'
    >>> mixture.parameters.names[:3]
    ('component0.amplitude', 'component0.period', 'component0.quality')
    >>> mixture.QUASISEPARABLE
    True

3. Non-stationarity: ``WarpedKernel``
-----------------------------------------------------------

Every family above is **stationary**: the covariance depends only on the
separation between two points, not on where they sit. :class:`~ampere.core.WarpedKernel`
(**W5.7**) is the flexible likelihood's answer to a residual that is not —
one or two wrappers around any base kernel, each preserving the exact O(N)
solve rather than trading it away:

* **Input warping** maps the coordinate through a monotone piecewise-linear
  spline over a small number of fixed knots, ``k(x, x') -> k(w(x), w(x'))``.
  A length scale that should be short in one band and long in another is
  this, not a new family: :math:`w` compresses the coordinate where the
  process varies quickly and stretches it where it does not. The knot
  *locations* are fixed and explicit in the declaration — :func:`~ampere.core.quantile_knots`
  is the convenience for choosing them from the data — and only the segment
  slopes between them are fitted parameters.
* **Amplitude warping** scales the marginal standard deviation by a second,
  independent piecewise-linear spline in :math:`\log a(x)`, giving
  :math:`D K D` with :math:`D = \mathrm{diag}(a(x))` for a residual whose
  *size*, not its correlation structure, varies across the band.

Both act on the same one coordinate a kernel already runs its recursion on,
and both preserve quasiseparability **exactly**, for reasons that are
structural rather than numerical (``likelihoods.md`` §6 has the full
account): a monotone warp of the coordinate does not change the sorting
permutation ``QuasisepGP``'s ordering precondition needs, so the base
kernel's generators evaluated at :math:`w(x)` are still a valid
factorisation; and a diagonal congruence :math:`D K D` of a rank-:math:`J`
semiseparable matrix is rank-:math:`J` semiseparable, with the diagonal
absorbed straight into the generators (:math:`U \to \mathrm{diag}(a)\,U`,
:math:`V \to \mathrm{diag}(a)\,V`). Nothing about the solver, the
approximation or the representation changes — a warped Matérn-3/2 is still
exactly quasiseparable and still runs on :class:`~ampere.core.QuasisepGP`:

.. code-block:: pycon

    >>> import scipy.stats as st
    >>> from ampere.core import Matern32, WarpedKernel, quantile_knots
    >>> base = Matern32(st.loguniform(1e-3, 1e1), st.loguniform(0.1, 10.0))
    >>> warped = WarpedKernel(base, input_warp=(0.0, 5.0, 10.0))
    >>> warped.QUASISEPARABLE
    True
    >>> sorted(name for name in warped.parameters.names if "warp" in name)
    ['input_warp.increment0', 'input_warp.increment1', 'input_warp.scale']

The knot variables (``input_warp.incrementN``, ``amplitude_warp.levelN``)
and each warp's shared shrinkage scale (``input_warp.scale``,
``amplitude_warp.scale``) are **ordinary** :class:`~ampere.core.Parameter`
objects — nothing about fitting, lowering or provenance treats them
differently from a kernel's own hyperparameters, so NUTS reaches all of
them on the torch and jax backends for free. Their default priors are not
arbitrary: every knot variable is hierarchical under a single half-normal
scale per warp (non-centred by default, the better posterior geometry for
NUTS), which puts the *identity* warp — the base kernel, exactly, bit for
bit — at the prior's centre of mass, so the data must pay to move away from
it. This is the same shape :func:`~ampere.core.shrinkage_horseshoe` gives a
``Sum``'s component amplitudes (§1 above): few knots, and a hierarchical
prior that shrinks toward "no effect" rather than a hard limit on how many
knots are allowed.

**The worked case.** ``examples/m2_misspecification/many_lines.py``
(W5.8) is the M2 scenario a single length scale cannot cover — a forest of
narrow lines beside a smooth, much broader continuum error — and it is
where a warped Matérn-3/2 is compared against the plain kernel and against
a ``Sum`` of two Matérn-3/2 terms. :doc:`m2_misspecification` has the full
account and the numbers; in short, the plain stationary kernel is the one
place in the whole M2 study that *fails* the page's own 1.5-posterior-width
threshold (2.84), while the warped kernel (0.92) and the two-length-scale
sum (0.78) both recover it, all three localising the deviation inside the
line band.

**Caution.** Both warps are piecewise-linear **within** their knot range
and extended **linearly** beyond the end knots — cheap and monotone, but
not flat, so a warp whose knots do not cover the data's own range keeps
extrapolating past the last one instead of levelling off. The reviewer's
note carried from W5.7's review is the rule of thumb: choose knots (by eye,
or with :func:`~ampere.core.quantile_knots`) that **cover the data range**,
rather than relying on the extrapolation to do the right thing outside it —
a linearly extrapolated :math:`\log a` is the amplitude-warp case this
bites hardest, since it is exponentiated.

4. The ``axes=`` selector, and the per-leaf unit rule
-----------------------------------------------------------

Ruled 2026-09-11 (``phase4_placement_memo.md`` §3.6 item 3): a kernel acts
on a **named subset** of a container's axes, not on all of them by default
alone. Every :class:`~ampere.core.Kernel` and :class:`~ampere.core.KernelSpec`
takes ``axes=``:

.. code-block:: pycon

    >>> Matern32(0.3, 2.0).axes is None            # every axis (the pre-W4.5 default)
    True
    >>> Matern32(0.3, 2.0, axes=("u", "v")).axes
    ('u', 'v')

This exists because a container may carry axes in **mixed units** — a
:class:`~ampere.core.VisibilitySet`'s dimensionless ``(u, v)`` and its
spectral ``spectral_axis`` — and an isotropic kernel over all three would
average a baseline length against a wavelength, which is meaningless and
which ``ampere`` refuses rather than silently computes:

.. code-block:: pycon

    >>> from ampere.core import DenseGP, VisibilitySet
    >>> a_visibility_set = VisibilitySet(
    ...     [120.0, -35.0, 88.0, -210.0], [45.0, 190.0, -66.0, 12.0],
    ...     [1.3, 1.3, 0.87, 0.87] * u.mm, [1.0 + 0.2j, 0.6 - 0.3j, 0.4 + 0.0j, 0.1 - 0.05j],
    ...     uncertainty=[0.02, 0.02, 0.03, 0.05],
    ... )
    >>> GaussianProcessNoise(Matern32(0.3, 2.0), DenseGP()).check_compatible(
    ...     GaussianFamily(), a_visibility_set)
    Traceback (most recent call last):
        ...
    ampere.core.exceptions.LikelihoodError: DenseGP measures separation as a Euclidean
    distance across a VisibilitySet's coordinate axes, but they carry different units
    [...]. A single isotropic length-scale is meaningless across mixed units; name the
    axes this kernel acts on with axes=(...), as in Matern32(axes=('u',)), and compose kernels on different axes
    with Product.

The remedy is the selector: ``Matern32(axes=("u", "v"))`` is the isotropic
``(u, v)`` kernel, restricted to the two dimensionless axes, so the
single-unit rule applies to the *subset* rather than to the whole
container. The binding is **functional**, through :meth:`~ampere.core.Kernel.for_axes`,
so one kernel instance may be shared between datasets whose axes differ —
nothing is mutated in place — and the selection is part of the kernel's
own hash (``KernelSpec.to_dict()`` includes ``axes`` only when a selection
was made), so a pre-W4.5 spec hash is unaffected.

5. The chromatic case: a ``Product`` on disjoint axes
-----------------------------------------------------------

The axis selector's own reason for existing (memo §3.6, measured rather
than argued in :doc:`interferometry` §6): a missing sky component with a
band profile produces a visibility residual that is **sharp in wavelength
and smooth in spatial frequency** — two different correlation structures on
two different coordinates, which one isotropic kernel cannot express at
all, approximately or otherwise, because a kernel sees a container's axes
and nothing else. The fix composes two axis-selected kernels with
``Product``:

.. code-block:: python

    from ampere.core import Matern32, Product

    spatial = Matern32(amplitude=st.halfnorm(scale=0.1),
                        length_scale=st.loguniform(2e7, 2e8), axes=("u", "v"))
    spectral = Matern32(amplitude=1.0,                      # fixed: see below
                         length_scale=st.loguniform(0.01, 1.0), axes=("spectral_axis",))
    kernel = Product(spatial, spectral)     # k(u, v, lambda) = k_uv(u, v) . k_lambda(lambda)

This is W4.2's circular complex GP's chromatic form — the same closed-form
marginalisation as the isotropic case (:class:`~ampere.core.ComplexGaussianFamily`
plus :class:`~ampere.core.GaussianProcessNoise` is ``ANALYTIC``, one
Cholesky of ``K(theta) + diag(sigma**2)``), with ``K`` now built from the
product kernel's dense covariance rather than an isotropic one. One of the
two amplitudes must be fixed rather than fitted: a product's marginal
variance is the product of its terms', so two free amplitudes
over-parameterise the covariance by one degree of freedom.
:doc:`interferometry` §6 measures the claim rather than arguing it, on a
synthetic chromatic patch observed at dispersed ``(u, v)`` coverage, and
reports the honest result — the product does not straightforwardly
dominate a spatial-only kernel at the per-PR budget, though the
spectral-only kernel is clearly the worst of the three throughout.

6. Registering your own term: ``register_quasiseparable_term``
---------------------------------------------------------------------

A user-defined :class:`~ampere.core.Kernel` subclass works on
:class:`~ampere.core.DenseGP` the moment it implements ``matrix``/
``diagonal`` — the ABC is public and that is all the dense solver needs.
Reaching :class:`~ampere.core.QuasisepGP` needs one more thing: a
**celerite representation**, registered once, for every backend, through
:func:`~ampere.core.register_quasiseparable_term`. Before W4.5 this table
was private; it is now the public route the design horizon notes asked
for, in the lowering/realisation registries' own shape — one slot per
kernel family, no silent overwrite, a built-in row distinguished from a
user one.

The out-of-tree example ``tests/core/test_kernels.py`` uses is deliberately
trivial — a kernel that is Matérn-1/2 under another name, so its
mathematics need no checking and the test proves the *route* rather than
the arithmetic:

.. code-block:: pycon

    >>> from ampere.core import Matern12, register_quasiseparable_term
    >>> from ampere.core.kernels import matern12_representation
    >>> class Relabelled(Matern12):
    ...     """A user kernel: Matérn-1/2 under another name."""
    ...     FAMILY = "relabelled_matern12"
    >>> row = register_quasiseparable_term(Relabelled, matern12_representation)
    >>> row.family, row.builtin
    ('relabelled_matern12', False)

A builder is ``(kernel, values, axis) -> CeleriteRepresentation`` — ``values``
the kernel's resolved hyperparameters, ``axis`` the sorted coordinate in the
backend's own array type, and ``kernel.ops`` the array namespace to build in
(:class:`~ampere.core.kernels.ArrayOps`, so one registration serves all
three backends: numpy through celerite2's own ``GaussianProcess``, torch
through ``celerite2.backprop``, jax through ``celerite2.jax.ops``). A second
registration for one family is refused unless ``override=True`` — a
silently replaced covariance is a changed model, the same discipline the
lowering registry applies — and a kernel that declares
``QUASISEPARABLE = False`` cannot register a representation at all: the
contradiction is caught here rather than left for `QuasisepGP` to discover
later. `registered_quasiseparable_terms()` and `term_provenance_entries()`
report what is registered, the second only counting user rows, so a run's
provenance can be checked against the terms it actually used.

See also
------------

* :doc:`concept` — the flexible likelihood in one page, with the kernel
  algebra summarised for a first read.
* :doc:`interferometry` — the ``axes=`` selector's own reason for existing,
  and the chromatic case measured end to end.
* :doc:`advanced` — the noise-model section, for where this page fits
  among ampere's other extension points.
* ``docs/design/contracts/likelihoods.md`` §6-§8 — the frozen kernel
  contract this page documents the shipped implementation of.
