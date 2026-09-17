"""Run the astrometry example: ``python -m examples.astrometry``.

    python -m examples.astrometry                      # reference, emcee
    python -m examples.astrometry --backend torch       # NUTS
    python -m examples.astrometry --gp                  # the flexible likelihood
    python -m examples.astrometry --joint               # the joint channel noise (W5.9)
    python -m examples.astrometry --sbc joint           # the calibration study
    python -m examples.astrometry --walkers 32 --steps 2000

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
