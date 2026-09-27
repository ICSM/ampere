"""Run the astrometry example: ``python -m examples.astrometry``.

    python -m examples.astrometry                      # reference, emcee
    python -m examples.astrometry --backend torch       # NUTS
    python -m examples.astrometry --gp                  # the flexible likelihood
    python -m examples.astrometry --joint               # the joint channel noise (W5.9)
    python -m examples.astrometry --sbc joint           # the calibration study
    python -m examples.astrometry --walkers 32 --steps 2000
    python -m examples.astrometry --wide-prior          # the aliasing hazard, its remedy (W5.15)
    python -m examples.astrometry --wide-prior --engine nautilus  # or ultranest, in `-e nested`
    python -m examples.astrometry --wide-prior --engine dynesty --live-points 1000 --dlogz 0.05

``--wide-prior``, ``--engine``, ``--live-points`` and ``--dlogz`` are the
period-aliasing arm's own flags -- :func:`astrometry.build_problem`'s and
:func:`astrometry.fit`'s docstrings carry the full account (W5.15).
``--engine`` chooses among the three nested samplers (``dynesty``,
``nautilus``, ``ultranest``) and defaults to dynesty once ``--wide-prior`` is
given; ``--live-points`` sets every nested sampler's own live-set size
(``None`` defaults to :func:`~ampere.inference.engine.default_live_points`);
``--dlogz`` is dynesty's and ultranest's stopping criterion (nautilus takes
its own ``f_live``/``n_eff`` instead).

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, exactly as :mod:`examples.sed_composition`'s own
``__main__.py`` does and for the same reason: a dozen-epoch, six-parameter
problem is far too small for multi-threaded BLAS to pay for its own
thread-pool overhead.
"""

from __future__ import annotations

import os

for _threads in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_threads, "1")

from .astrometry import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
