"""Run the carbon-star example: ``python -m examples.cstar``.

    pixi run -e hyperion python -m examples.cstar                  # the data, legacy budget
    pixi run -e hyperion python -m examples.cstar --embedding
    pixi run -e hyperion python -m examples.cstar --synthetic --photons quick --rounds 1

Caps every BLAS/OpenMP thread pool at four threads **before numpy is imported
anywhere in this process**: Hyperion runs serially, one model per worker
process, and the parallelism is across workers (``--workers``), so a pool per
worker would only oversubscribe the machine. Torch's own intra-op pool, used
by the network's training, is capped the same way. Scoped to a subprocess
invocation of this module, never imposed on a caller who imports
:mod:`.cstar` as a library.
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

from .cstar import main  # noqa: E402 -- must follow the thread-pool env vars above

if __name__ == "__main__":
    try:
        import torch

        torch.set_num_threads(4)
    except ImportError:
        pass
    raise SystemExit(main())
