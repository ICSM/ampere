"""**A found blocker, worked around in this example rather than in the library.**

Read this module before the rest of the study; it is the one place where
``examples/image`` does something a user should not have to do.

What is blocked
---------------
``ampere.core`` refuses a correlated noise model on **any** ``Layout.GRID``
container, in two places, both of which predate an ``Image`` ever being an
observation and both of which say so:

``GPSolver.check_compatible``
    "the v1 GP solvers work on point-set containers (Spectrum, TimeSeries,
    PhotometricPoints, VisibilitySet); gridded 2D+ data are the subject of the
    SVGP / SKI / Vecchia strategy slots (DEVELOPMENT_PLAN.md §4.4, Phase 5)."

``Likelihood._coordinates``
    "a correlated noise model needs point-set coordinates, but a {kind} has a
    grid layout. Gridded 2D+ data are the SVGP / SKI / Vecchia strategy slots
    (Phase 5)."

Neither refusal is about mathematics. A stationary kernel over ``(x, y)`` is a
function of coordinates, and a grid has coordinates — it just keeps them
separably, as one array per axis, so the ``(N, 2)`` matrix a kernel wants has
to be broadcast out of them rather than column-stacked. That broadcast already
exists and has been correct since W3.3: ``ampere.core.encoding``'s
``_coordinate_matrix`` does it for the SBI encoder, and W5.5 lifted the same
rule into ``ampere.core.sample_coordinates`` so that
``Dataset.draw_observation`` could draw an image at all.

So the two refusals are **the Phase 5 slot they name**, still closed. W5.5's
own item asks for the flexible arm on ``DenseGP`` at small N and on
``HilbertSpaceGP(basis_size=(m, m))`` at realistic N — and W5.4's
``HilbertSpaceGP`` docstring says in as many words that the method "is the
scaling answer in **two and three axes**, where ``QuasisepGP`` does not apply
at all". The pieces are all present; the gate is shut.

What this module does about it, and what it must not be mistaken for
--------------------------------------------------------------------
W5.5's file ownership puts ``ampere/core/likelihood.py`` out of scope, so the
gate is not opened here. Instead this module subclasses the two classes that
hold it — from the public API, out of tree, in exactly the way
``transformations.md`` §11 exists to prove is possible — and lifts **only** the
layout check, leaving every other rule (the kernel's axis selection, the
quasiseparable refusals, the ordered-1D rule) to the base classes. That is
three small overrides, and each is written so that a point-set container takes
the base class's path unchanged.

It is a *demonstration that nothing but the gate is missing*, not a design.
The library change these stand in for is two edits in
``ampere/core/likelihood.py``: allow ``Layout.GRID`` in
``GPSolver.check_compatible`` when the kernel's selected axes are the
container's, and build ``Likelihood._coordinates`` with
``ampere.core.sample_coordinates`` instead of a bare ``column_stack``. It needs
a decision-log entry (ground rule 9, ``likelihoods.md`` §7's limitation), which
is why it is proposed in W5.5's report rather than taken here. **Delete this
module the day that lands**, and the study will import ``DenseGP`` and
``HilbertSpaceGP`` directly.

What it deliberately does not claim
-----------------------------------
Nothing here makes a grid *cheap*. ``DenseGP`` on a 128x128 image is a
16,384-square Cholesky, which is why the study measures rather than asserts,
and why the HSGP arm exists at all. The SKI and Vecchia slots the refusal names
are about exploiting a grid's structure, which neither of these solvers does;
that is W5.6's bake-off, and this module has nothing to say about it.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from ampere.core import (
    DenseGP,
    FunctionSamples,
    HilbertSpaceGP,
    Kernel,
    Layout,
    Likelihood,
    sample_coordinates,
)

__all__ = ["GridDenseGP", "GridHilbertSpaceGP", "GridLikelihood"]


class _GriddedSolver:
    """Lift *only* the ``Layout.POINTS`` gate in :meth:`GPSolver.check_compatible`.

    A mixin, so each concrete solver still inherits every other check its base
    class makes — and so that a point-set container still takes the base
    class's path, byte for byte, which is what keeps this module from changing
    the answer to any question anyone was already asking.
    """

    def check_compatible(self, kernel: Kernel, observed: FunctionSamples) -> None:
        if observed.LAYOUT is Layout.GRID:
            # The one check that still applies on a grid: the kernel must
            # select axes this container actually has, in compatible units.
            kernel.check_axes(observed, owner=self.NAME)  # type: ignore[attr-defined]
            return
        super().check_compatible(kernel, observed)  # type: ignore[misc]


class GridDenseGP(_GriddedSolver, DenseGP):
    """:class:`~ampere.core.DenseGP`, allowed to see an ``Image``.

    Exact, ``O(N³)``, and the study's small-N arm. At 24x24 that is a
    576-square Cholesky per evaluation, which is affordable; at 128x128 it is
    16,384-square, which is the point the benchmark is making.
    """


class GridHilbertSpaceGP(_GriddedSolver, HilbertSpaceGP):
    """:class:`~ampere.core.HilbertSpaceGP`, allowed to see an ``Image``.

    W5.4's reduced-rank solver, whose tensor-product basis is exactly what two
    axes want: ``basis_size=(m, m)`` gives ``m_total = m²`` basis functions and
    a cost of ``O(N m_total + m_total³)``, which is why ``m`` **per axis** stays
    small — eight per axis is sixty-four in total, and the cube of sixty-four
    is nothing beside the cube of sixteen thousand. See that class's docstring
    for the boundary factor, and for why neither it nor ``m`` is chosen for you.
    """


class GridLikelihood(Likelihood):
    """:class:`~ampere.core.Likelihood`, able to build a grid's coordinate matrix.

    One override, and it delegates: :func:`~ampere.core.sample_coordinates` is
    ``ampere.core``'s own layout-aware stacking, so the coordinates a kernel
    sees here are the same ones the SBI encoder and
    ``Dataset.draw_observation`` see. A point-set container takes the base
    class's path.
    """

    def _coordinates(self, observed: FunctionSamples, retain: np.ndarray) -> np.ndarray:
        if observed.LAYOUT is Layout.GRID:
            return np.ascontiguousarray(sample_coordinates(observed)[retain])
        return super()._coordinates(observed, retain)


def gridded(likelihood_family: Any, noise: Any, **kwargs: Any) -> GridLikelihood:
    """A :class:`GridLikelihood` over *likelihood_family* and *noise*."""
    return GridLikelihood(likelihood_family, noise, **kwargs)
