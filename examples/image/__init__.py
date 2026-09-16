"""An image as a dataset: the PSF-convolution modality, end to end (W5.5).

The third worked modality in this repository, after
:mod:`examples.interferometry` and the astrometric time series, and the first
whose observed container is **gridded** — a :class:`~ampere.core.Image`, two
separable axes, ``Layout.GRID``. ``docs/source/image.rst`` is this package's
prose; run it with::

    python -m examples.image                      # three arms, reference, emcee
    python -m examples.image --calibration        # the SBC coverage row
    python -m examples.image --benchmark          # DenseGP vs HSGP at three N
    python -m examples.image --backend torch      # NUTS

What is here
------------
:mod:`.generators`
    The truth — a compact Gaussian source on a smooth, much broader background
    — and the noisy image of it a PSF-convolving camera records.
:mod:`.model`
    ``SourceWithBackground``: a small user-written
    :class:`~ampere.core.transform.Model` that delegates to two of W4.1's
    shipped image models and adds them, rather than a new class in ``ampere``.
:mod:`.model_torch`, :mod:`.model_jax`
    Two-line bindings of that class to each modern backend's own models.
:mod:`.grid_gp`
    **Read this one first.** The two ``Layout.POINTS`` gates that refuse a
    correlated noise model on an ``Image`` today, lifted out of tree so that
    the flexible arm can run, with the library change they stand in for written
    out in full.
:mod:`.study`
    The three arms, the SBC calibration route, and the ``DenseGP``-against-HSGP
    benchmark that is Phase 5's own "chosen by measurement" rule made concrete.
:mod:`.figures`
    What can be drawn: a corner plot, a trace, and hand-rolled image panels.
    Four of the six shipped plots refuse a two-axis kind today — see the
    module, and ``docs/source/image.rst``'s closing section.

Nothing here is committed as output (ground rule 7): ``--figures DIR`` writes
to a directory you name.
"""

from __future__ import annotations

__all__ = ["figures", "generators", "grid_gp", "model", "study"]
