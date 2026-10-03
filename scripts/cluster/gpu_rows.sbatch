#!/usr/bin/env bash
# W6.16: run ampere's GPU rows (tests/gpu) on one A100. Submitted by
# scripts/cluster/run.sh, which exports the variables below and names the log
# with --output. No partition, account or QoS is named, by site policy: the
# cluster that runs this is the cluster of the login host it was submitted from.
#
# Required environment (exported by run.sh): AMPERE_REPO, AMPERE_PIXI_HOME,
# AMPERE_TAG, AMPERE_CLUSTER_NAME. Optional: AMPERE_MODULES (space-separated).
#SBATCH --job-name=ampere-gpu-rows
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=ampere-gpu-rows-%j.log
set -euo pipefail

: "${AMPERE_REPO:?not set -- submit through scripts/cluster/run.sh}"
: "${AMPERE_PIXI_HOME:?not set -- submit through scripts/cluster/run.sh}"
: "${AMPERE_TAG:?not set -- submit through scripts/cluster/run.sh}"

# The CUDA 12 wheels (torch cu126, jax[cuda12]) bundle their libraries and need
# only a driver at or above this floor (NVIDIA's CUDA 12 minor-version
# compatibility). Below it the rows would skip or fail for a reason that has
# nothing to do with ampere, so the job refuses to continue.
DRIVER_FLOOR=525

echo "== ampere GPU rows: tag ${AMPERE_TAG}, cluster ${AMPERE_CLUSTER_NAME:-unknown}, job ${SLURM_JOB_ID:-none} on $(hostname -s), $(date -u +%Y-%m-%dT%H:%M:%SZ)"

# 1. What this node actually has -- first, so the log records it whatever follows.
echo "== nvidia-smi"
nvidia-smi
echo "== module spider cuda"
# Lmod writes its listing to stderr; the log wants it. Not fatal if absent.
if command -v module >/dev/null 2>&1; then
  module spider cuda 2>&1 || echo "(module spider found no cuda module)"
else
  echo "(no module command on this node)"
fi

driver="$(nvidia-smi --query-gpu=driver_version --format=csv,noheader | head -n 1)"
driver_major="${driver%%.*}"
echo "== driver ${driver}, floor ${DRIVER_FLOOR}"
if [ "${driver_major:-0}" -lt "$DRIVER_FLOOR" ]; then
  echo "REFUSED: driver ${driver} is below the CUDA 12 wheels' floor (${DRIVER_FLOOR}). Re-solve the gpu environment against what nvidia-smi reports (docs/development.md, 'The GPU rows on the cluster')." >&2
  exit 3
fi

# 2. Modules, if the local config named any (module lines belong here, never in dotfiles).
if [ -n "${AMPERE_MODULES:-}" ]; then
  # shellcheck disable=SC2086  # deliberate word-splitting: a list of module names
  module load ${AMPERE_MODULES}
  echo "== loaded: ${AMPERE_MODULES}"
fi

# 3. The rows.
export PIXI_HOME="$AMPERE_PIXI_HOME"
export PIXI_CACHE_DIR="$AMPERE_PIXI_HOME/cache"
export OMP_NUM_THREADS=4
cd "$AMPERE_REPO"
echo "== checkout: $(git describe --tags --always) ($(git rev-parse HEAD))"
out="$(mktemp)"
trap 'rm -f "$out"' EXIT
set +e
"$AMPERE_PIXI_HOME/bin/pixi" run -e gpu gpu -q --tb=short -p no:cacheprovider 2>&1 | tee "$out"
status="${PIPESTATUS[0]}"
set -e
echo "== pytest exit status ${status} (tag ${AMPERE_TAG}, job ${SLURM_JOB_ID:-none})"
# The summary line last, so `tail -n 1` of the log is the line the status row quotes.
grep -E '[0-9]+ (passed|failed|skipped|error)' "$out" | tail -n 1 || echo "(no pytest summary line found)"
exit "$status"
