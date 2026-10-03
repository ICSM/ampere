#!/usr/bin/env bash
# W6.16: one-time (and re-runnable) set-up of ampere's GPU environment on the
# cluster. Run it ON THE LOGIN NODE, normally through
#   scripts/cluster/run.sh bootstrap
# which ssh-es the login node and feeds this file to bash with the variables
# below set from your local config.yaml.
#
# It is idempotent and has two phases, because the deploy key must be added on
# GitHub by a person between them:
#   1. lay out the directories, install the static pixi binary, generate a
#      read-only deploy key and print its public half. If the clone fails
#      (the key is not yet registered) it stops there, exit status 0.
#   2. re-run it once the key is added: clone ICSM/ampere, then
#      `pixi install -e gpu --locked` (a few gigabytes; the install is the
#      reason everything lives under the project space and not under /home,
#      which is a hard 20 GB).
# It does NOT smoke-test the GPU: that is the job's first act (gpu_rows.sbatch).
#
# Required environment: AMPERE_FRED (project space root, e.g. /fred/<project>),
# AMPERE_REPO, AMPERE_PIXI_HOME. Optional: PIXI_VERSION (default below).
set -euo pipefail

: "${AMPERE_FRED:?set AMPERE_FRED (the project space root)}"
: "${AMPERE_REPO:?set AMPERE_REPO (the clone's directory)}"
: "${AMPERE_PIXI_HOME:?set AMPERE_PIXI_HOME (pixi's home and cache)}"
PIXI_VERSION="${PIXI_VERSION:-0.68.1}"
REPO_URL="git@github.com:ICSM/ampere.git"

case "$AMPERE_FRED" in
  CHANGE_ME* | "") echo "bootstrap: AMPERE_FRED is unset or the example value" >&2; exit 2 ;;
esac

root="$AMPERE_FRED/ampere"
mkdir -p "$root" "$root/logs" "$AMPERE_PIXI_HOME/bin" "$AMPERE_PIXI_HOME/cache" "$HOME/.ssh"
echo "bootstrap: layout under $root"

# --- pixi: a static binary under the project space -------------------------
pixi_bin="$AMPERE_PIXI_HOME/bin/pixi"
if [ ! -x "$pixi_bin" ]; then
  url="https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-x86_64-unknown-linux-musl.tar.gz"
  echo "bootstrap: installing pixi ${PIXI_VERSION} from $url"
  curl -fsSL "$url" | tar -xz -C "$AMPERE_PIXI_HOME/bin" pixi
fi
export PIXI_HOME="$AMPERE_PIXI_HOME"
export PIXI_CACHE_DIR="$AMPERE_PIXI_HOME/cache"
echo "bootstrap: $("$pixi_bin" --version) at $pixi_bin"

# --- the read-only deploy key ----------------------------------------------
key="$HOME/.ssh/ampere_deploy_ed25519"
if [ ! -f "$key" ]; then
  ssh-keygen -t ed25519 -N "" -C "ampere-deploy-$(hostname -s)" -f "$key"
  chmod 600 "$key"
fi
echo
echo "bootstrap: the PUBLIC half of the deploy key (add it at"
echo "  https://github.com/ICSM/ampere/settings/keys, title e.g. 'cluster', and"
echo "  leave 'Allow write access' UNCHECKED):"
echo
cat "$key.pub"
echo

# --- the clone ---------------------------------------------------------------
export GIT_SSH_COMMAND="ssh -i $key -o IdentitiesOnly=yes -o StrictHostKeyChecking=accept-new"
if [ -d "$AMPERE_REPO/.git" ]; then
  echo "bootstrap: clone already present at $AMPERE_REPO"
elif git clone "$REPO_URL" "$AMPERE_REPO"; then
  echo "bootstrap: cloned into $AMPERE_REPO"
else
  echo "bootstrap: the clone failed -- most likely the deploy key above is not yet"
  echo "  registered on the repository. Add it, then run this script again."
  exit 0
fi
git -C "$AMPERE_REPO" config core.sshCommand "$GIT_SSH_COMMAND"

# --- the environment ---------------------------------------------------------
cd "$AMPERE_REPO"
git fetch --tags --quiet
echo "bootstrap: pixi install -e gpu --locked (at $(git rev-parse --short HEAD))"
"$pixi_bin" install -e gpu --locked
echo "bootstrap: done. Next: scripts/cluster/run.sh submit <tag>"
