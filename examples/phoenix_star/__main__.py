"""Run the PHOENIX star example: ``python -m examples.phoenix_star``.

    python -m examples.phoenix_star --engine sbi                   # pixi run -e sbi
    python -m examples.phoenix_star --engine nuts --backend torch  # pixi run -e torch
    python -m examples.phoenix_star --engine emcee

Caps every BLAS thread pool at four threads before numpy is imported (a
caller's own setting wins), as :mod:`examples.linear_sed`'s ``__main__`` does
at one: the grids here are a few thousand points.
"""

from __future__ import annotations

import os

for _threads in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_threads, "4")

from .phoenix_star import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
