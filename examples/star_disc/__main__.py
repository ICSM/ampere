"""Run the star-plus-disc example: ``python -m examples.star_disc``.

    python -m examples.star_disc                # the votable + the RVS spectrum
    python -m examples.star_disc --synthetic    # the coverage run's data
    python -m examples.star_disc --irs PATH     # plus a fetched CASSIS file

Caps every BLAS thread pool at four threads before numpy is imported (a
caller's own setting wins), as :mod:`examples.phoenix_star`'s ``__main__``.
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

from .star_disc import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    raise SystemExit(main())
