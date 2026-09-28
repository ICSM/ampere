"""Run the photometry-and-spectra example: ``python -m examples.photometry_spectra``.

    python -m examples.photometry_spectra                          # untied, independent noise
    python -m examples.photometry_spectra --tie                    # one shared calibration factor
    python -m examples.photometry_spectra --gp                     # the flexible likelihood
    python -m examples.photometry_spectra --tie --gp
    python -m examples.photometry_spectra --figures /tmp/figures

Restricts every BLAS thread pool to one thread **before numpy is imported
anywhere in this process**, exactly as
:mod:`examples.sed_composition.__main__` does and for the same reason (see
its docstring): the composition's grids are a couple of thousand points at
most, too small for multi-threaded BLAS to repay its own thread-pool
overhead.
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

from .photometry_spectra import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
