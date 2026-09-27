GP solver strategies
=====================

:doc:`kernels` is the covariance function; this page is *how the algebra is
done* — the swappable :class:`~ampere.core.GPSolver` strategy behind
:class:`~ampere.core.GaussianProcessNoise` and (since W5.9) behind
:class:`~ampere.core.JointGaussianProcessNoise`. Every strategy computes the
same quantity, ``log N(residual; 0, K(theta) + diag(sigma^2))``, and
:class:`~ampere.core.DenseGP` is the definition of the right answer against
which every other one is measured (``likelihoods.md`` §7 is the frozen
contract this page documents the shipped implementation of).

.. code-block:: pycon

    >>> from ampere.core import (
    ...     DenseGP, QuasisepGP, HilbertSpaceGP, EquispacedFourierGP, VecchiaResponseGP,
    ... )
    >>> (DenseGP.EXACT, QuasisepGP.EXACT, HilbertSpaceGP.EXACT)
    (True, True, False)
    >>> (EquispacedFourierGP.EXACT, VecchiaResponseGP.EXACT)
    (False, False)

1. Exact: ``DenseGP`` and ``QuasisepGP``
------------------------------------------

:class:`~ampere.core.DenseGP` requires nothing about the kernel or the
container: any kernel, any layout, a dense ``N×N`` Cholesky, ``O(N³)``. It
is what every other strategy is checked against, and the fallback the error
message of every refusal below names.

:class:`~ampere.core.QuasisepGP` is the celerite2-backed O(N) path, and it
requires two things at composition, checked before anything about
implementation status: **one ordered axis** (``REQUIRES_ORDERED_1D`` counts
the axes the kernel *selects*, not how many the container has), and an
**exactly quasiseparable kernel** — one of :doc:`kernels`' families and
``Sum``\ s of them, since a ``Product`` of quasiseparable terms is not itself
quasiseparable and is refused by name rather than silently costing more:

.. code-block:: pycon

    >>> from ampere.core import GaussianProcessNoise, SquaredExponential
    >>> noise = GaussianProcessNoise(SquaredExponential(0.3, 2.0), QuasisepGP())
    >>> noise.solver.REQUIRES_QUASISEPARABLE
    True

Both are implemented on all three backends, and both are checked
bit-tightly against each other in the conformance battery (the
``cross_solver`` tolerance, ``1e-6`` by default) rather than against a
looser one, because they compute the same number by genuinely different
recursions.

**``Layout.GRID`` since W5.21.** Until W5.21, both solvers' composition-time
check refused *any* correlated noise model over a ``Layout.GRID`` container
outright — a gate that predated an image ever being observed, not a
mathematical restriction. A stationary kernel over a grid's axes is a
function of coordinates exactly as it is over a point set, and
:func:`~ampere.core.sample_coordinates` already broadcasts a grid's separable
axes into the ``(N, d)`` matrix a kernel wants (W5.5). ``GPSolver.check_compatible``
now accepts ``Layout.POINTS`` **or** ``Layout.GRID``, so
:class:`~ampere.core.DenseGP` and :class:`~ampere.core.HilbertSpaceGP` (below)
compose over an :class:`~ampere.core.Image` exactly as they do over a
:class:`~ampere.core.Spectrum` — see :doc:`image` §5. The ordered-1D rule is
unchanged: :class:`~ampere.core.QuasisepGP` still refuses a kernel that
selects both of a grid's axes, because it counts the kernel's *selected*
axes rather than the container's layout, so a 2-D grid needs ``DenseGP`` or
one of the approximate strategies below.

``conditional_loo`` — the leave-one-out decomposition ``results.md`` §6
names — is where the two exact strategies diverge: ``DenseGP`` computes it
in closed form from the same Cholesky (Sundararajan & Keerthi 2001); the
**reference** ``QuasisepGP`` defers it (celerite2's public numpy interface
exposes no O(N) route to the diagonal a leave-one-out term needs) and
refuses by name rather than quietly costing O(N²) under an O(N) name. Both
differentiable backends' own ``QuasisepGP`` supply it anyway, by the same
O(N) backward accumulation their autodiff already does over the
factorisation — the deferral was a coupling to celerite2's numpy
factorisation convention, not a fact about the mathematics.

2. Approximate: reduced rank, ``EXACT = False``
--------------------------------------------------

:class:`~ampere.core.HilbertSpaceGP` (**W5.4**) is the phase's one shipped
reduced-rank solver, on all three backends. On a box containing the data, a
stationary kernel with a closed-form spectral density (the three Matérns,
:class:`~ampere.core.SquaredExponential`, :class:`~ampere.core.SHO`, and any
``Sum`` of those) is diagonalised by the Dirichlet Laplacian's eigenfunctions,
``K ≈ Φ diag(S) Φᵀ`` with ``m`` basis members, and every quantity
follows from Woodbury at ``O(N m + m³)``. Two things the plan's Phase 5
bullet got wrong in passing, both corrected at W5.4: the approximation's
parameters are **dataclass fields**, not merely recorded metadata, so they
enter the spec hash — two runs at different ``m`` are not the same model —
and the latent block a NUTS draw samples is sized ``m``, not ``N``:

.. code-block:: pycon

    >>> solver = HilbertSpaceGP(basis_size=32, boundary_factor=2.0)
    >>> solver.provenance_config()
    {'basis_size': [32], 'boundary_factor': 2.0}
    >>> from ampere.core import Matern32
    >>> solver.latent_size(Matern32(0.3, 2.0), 500)
    32

``conditional_loo`` is **exact in the approximation** here — the same
Sundararajan & Keerthi identity ``DenseGP`` uses, applied to the Woodbury
factorisation ``HilbertSpaceGP`` actually holds — rather than deferred the
way the reference ``QuasisepGP``'s is. That is the general rule this
contract holds every ``EXACT = False`` strategy to: ``conditional_loo`` is
either exact in the approximation, or refused by name; never silently
approximate on top of an approximation.

**Two prototypes that stay in the tree, and are not promoted.**
:class:`~ampere.core.EquispacedFourierGP` and
:class:`~ampere.core.VecchiaResponseGP` were built to measure against
``HilbertSpaceGP`` at W5.6's bake-off, on the reference backend only — no
torch or jax twin, so a problem built with either is a numpy-only problem.
Measured at N = 4096 (Matérn-3/2, a 8 mas length scale), the disagreement
with ``DenseGP`` was 2.1e-1 nats for HSGP at m = 256 against 1.1e-2 for EFGP
at m = 289 and 30.2 for Vecchia at k = 30; at N = 65 536, HSGP took 0.41 s
against EFGP's 0.64 s and Vecchia's 14.3 s. **EFGP** is the one worth a
follow-on: at N = 16 384 it overtakes HSGP's cost past m ~ 576 (0.59 s
against 0.92 s at m ~ 1024, at 4.7x less memory), because its Toeplitz
normal equations solve by FFT-accelerated CG rather than a dense Cholesky
(0.16 s against 10.95 s at m = 6561, on memory alone 2.6 MiB against 1.31
GiB) — but nothing in the shipped modalities is at that scale yet, so it
stays a prototype rather than a third backend twin. **Vecchia** stays a
**slot** (``VecchiaGP``) rather than being promoted from its
``VecchiaResponseGP`` prototype: it is the right method for a short-range,
*rough* process in 2–3D that a spectral method needs a large ``m`` for, and
none of Phase 5's own modalities are that yet — the measurement is written
down for the day one is.

**The tolerance class.** An ``EXACT = False`` solver cannot be held to a
fixed number against ``DenseGP`` — how close it comes depends on the
kernel, the data and the approximation's own parameters — so the
conformance suite asserts a **convergence** instead
(``tests/conformance/protocol.py``'s ``Tolerances.approximation_order`` /
``approximation_floor`` / ``approximation_final`` and its
``approximation_envelope`` function, exercised by
``tests/conformance/test_likelihoods.py``'s
``TestApproximateSolverConvergence``, around line 1798). Refining the
approximation must reduce the disagreement with ``DenseGP`` at a rate the
method's own spectral analysis predicts for the kernel under test (an
isotropic Matérn-:math:`\nu`'s error falls as :math:`m^{-2\nu}`), down to a
floor the method's finite box or grid imposes rather than to zero, and the
*finest* setting in a sweep must actually reach a stated final tolerance. A
wrong implementation — the wrong box, a spectral density with the wrong
dimension, a Woodbury solve missing a factor — fails this in a way a fixed
tolerance would not catch, because more basis members converge to the wrong
process just as happily as to the right one.

3. Joint: noise over a tuple of channels
--------------------------------------------

:class:`~ampere.core.JointGaussianProcessNoise` (**W5.9**) is a solver
question too, not only a noise-model one: it scores ``T`` channels of one
model on **one shared grid** together, ``K = B ⊗ K_x``, where ``B`` is
a ``T×T`` :class:`~ampere.core.ChannelCoupling` and ``K_x`` is an ordinary
kernel on the grid. Diagonalising ``B = Q Λ Qᵀ`` and rotating the
residuals by ``Qᵀ`` decouples the joint problem into ``T`` scalar GPs
sharing ``K_x``, so the **bound solver** — ``DenseGP`` by default,
:class:`~ampere.core.QuasisepGP` on an ordered one-dimensional grid for the
O(N) path — decides once for every rotated output, exactly as it would for
one dataset:

.. code-block:: pycon

    >>> import numpy as np, scipy.stats as st
    >>> from ampere.core import JointGaussianProcessNoise, RotationCoupling
    >>> noise = JointGaussianProcessNoise(
    ...     Matern32(1.0, st.loguniform(1.0, 1e3), axes=("time",)),
    ...     QuasisepGP(),
    ...     datasets=("ra", "dec"),
    ...     coupling=RotationCoupling(
    ...         st.uniform(0.0, np.pi), st.norm(-7.0, 2.0), st.norm(-7.0, 2.0)
    ...     ),
    ... )
    >>> noise.datasets
    ('ra', 'dec')
    >>> noise.JOINT, noise.CORRELATED
    (True, True)

:class:`~ampere.core.RotationCoupling` is the ``T = 2`` case (a rotation
angle and two log-variances — the astrometric and polarimetric
parameterisation, in the physical coordinates rather than as three raw
matrix entries) and :class:`~ampere.core.CholeskyCoupling` the general
``T`` (``B = L Lᵀ``). ``K_x``'s own amplitude is refused as a free
parameter at construction, because a free amplitude and a free ``B`` are
exactly degenerate.

**Heteroscedastic channels (W5.24).** W5.9 shipped with one restriction:
every channel had to share the same per-sample uncertainty vector, because
``Qᵀ ⊗ I`` leaves ``diag(σ_t²)`` diagonal only when every channel's ``σ_t``
is the same vector. W5.24 lifts it with two routes, chosen by
``JointGaussianProcessNoise.route(variance)`` from the channels' variances
and the bound solver — never by an explicit argument, so equal channels
always keep the cheaper rotated path bit-identically:

* **Dense**, bound to :class:`~ampere.core.DenseGP`: ``B ⊗ (K_x + j²I)
  + blockdiag(diag(σ_t²))`` factorised directly, exact for any
  ``σ``, ``O((TN)³)``.
* **Reduced-rank**, bound to :class:`~ampere.core.HilbertSpaceGP` (or
  ``EquispacedFourierGP`` on the reference path only): Woodbury against the
  per-channel diagonal with ``K_x``'s own Φ̃ feature map, exact
  in the approximation at ``O(TN·(Tm)²)``, with a ``T·m`` whitened
  block.

Unequal variances under any other solver — :class:`~ampere.core.QuasisepGP`
above all, which has no Kronecker-free O(N) form once the diagonal couples
the rotated outputs — are refused at composition, by name, with the fix
named. **A jitter-convention note, carried from W5.24's review**: the two
unequal-variance routes add the solver's ``jitter`` differently — ``B ⊗
j²I`` on the dense route, a per-channel diagonal term on the reduced-rank
one — both defaulting to zero, so the difference is invisible until a
non-default jitter and unequal channels are combined; the general
(non-Kronecker) linear model of coregionalisation and mismatched grids
remain a follow-on beyond both routes.

See also
------------

* :doc:`kernels` — the covariance functions these strategies factorise;
  §3's ``WarpedKernel`` composes with every exact strategy above unchanged.
* :doc:`image` §5–§7 — ``Layout.GRID`` composing with ``DenseGP`` and
  ``HilbertSpaceGP``, measured.
* :doc:`astrometry` §10 — ``JointGaussianProcessNoise`` and heteroscedastic
  channels, the worked modality.
* ``docs/design/contracts/likelihoods.md`` §7 — the frozen solver-strategy
  contract this page documents the shipped implementation of, including the
  strategies that remain declared slots (``WindowedSparseGP``,
  ``InducingPointGP``, ``StructuredGridGP``, ``VecchiaGP``).
