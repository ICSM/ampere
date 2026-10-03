# scripts/cluster

The GPU rows (`tests/gpu`) on a cluster with A100s. The procedure, what the log
must show and why CI never runs these rows are in `docs/development.md`, "The GPU
rows on the cluster"; this directory is only the tools.

| File | Runs where | Does |
| --- | --- | --- |
| `config.example.yaml` | committed template | copy to `config.yaml` (gitignored) and fill in |
| `run.sh` | your machine | `bootstrap`, `submit <tag>`, `fetch <tag>`; only `ssh`/`rsync` |
| `bootstrap.sh` | login node (fed by `run.sh bootstrap`) | directories, pixi, deploy key, clone, `pixi install -e gpu --locked` |
| `gpu_rows.sbatch` | compute node (via `sbatch`) | `nvidia-smi`, `module spider cuda`, the driver-floor check, `pixi run -e gpu gpu` |

Nothing that identifies an account, a project or a key is committed; every such
value is read from `config.yaml`. The scripts refuse to run on the example values.
