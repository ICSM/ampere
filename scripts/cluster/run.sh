#!/usr/bin/env bash
# W6.16: the laptop-side driver for ampere's GPU rows on the cluster.
#
#   scripts/cluster/run.sh bootstrap        one-time set-up on the login node
#   scripts/cluster/run.sh submit <tag>     check out <tag> there and sbatch the rows
#   scripts/cluster/run.sh fetch  <tag>     copy the log back, through the data mover
#
# It only ever ssh-es the login node (and rsyncs through the data-mover host);
# it never runs the tests itself, and `submit` and `bootstrap` ask before they
# act. Every account-specific value is read from scripts/cluster/config.yaml
# (gitignored; the template is config.example.yaml) and the script refuses to
# run while any value is still CHANGE_ME.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
config="${AMPERE_CLUSTER_CONFIG:-$here/config.yaml}"
logdir="$HOME/.cache/ampere-gates/gpu"

usage() {
  sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 2
}

[ -f "$config" ] || { echo "run.sh: no $config -- copy config.example.yaml to config.yaml and fill it in" >&2; exit 2; }
case "$config" in
  *config.example.yaml) echo "run.sh: refusing to run with the example configuration" >&2; exit 2 ;;
esac

# cfg KEY: the value of a flat "key: value" line (comments and blanks skipped).
cfg() {
  sed -n "s/^$1:[[:space:]]*\\(.*[^[:space:]]\\)[[:space:]]*\$/\\1/p" "$config" | head -n 1
}
for k in cluster_name ssh_alias data_mover_alias user project_id fred_dir repo_dir pixi_home; do
  v="$(cfg "$k")"
  if [ -z "$v" ] || [ "$v" = "CHANGE_ME" ]; then
    echo "run.sh: '$k' in $config is empty or still CHANGE_ME" >&2
    exit 2
  fi
done
cluster_name="$(cfg cluster_name)"
ssh_alias="$(cfg ssh_alias)"
mover="$(cfg data_mover_alias)"
fred="$(cfg fred_dir)"
repo="$(cfg repo_dir)"
pixi_home="$(cfg pixi_home)"
# modules: a YAML flow list, "[]" or "[a/1, b/2]"; flattened to a space list.
modules="$(cfg modules | tr -d '[],"'"'" | xargs)"

confirm() {
  printf '%s [y/N] ' "$1" >&2
  read -r reply
  [ "$reply" = "y" ] || [ "$reply" = "Y" ] || { echo "run.sh: not confirmed; nothing done" >&2; exit 1; }
}

check_tag() {
  case "${1:-}" in
    "" | *[!A-Za-z0-9._+-]*) echo "run.sh: a tag (letters, digits, . _ + -) is required" >&2; exit 2 ;;
  esac
}

cmd="${1:-}"
case "$cmd" in
  bootstrap)
    echo "run.sh: this feeds scripts/cluster/bootstrap.sh to bash on '$ssh_alias'." >&2
    echo "  It lays out directories under $fred/ampere, installs pixi, generates a deploy" >&2
    echo "  key and prints its public half; run it again after the key is registered." >&2
    confirm "Run the bootstrap on $cluster_name?"
    # Config values are quoted for the remote shell with printf %q.
    # shellcheck disable=SC2029  # the remote command is meant to expand here
    ssh "$ssh_alias" "AMPERE_FRED=$(printf %q "$fred") AMPERE_REPO=$(printf %q "$repo") AMPERE_PIXI_HOME=$(printf %q "$pixi_home") bash -s" \
      < "$here/bootstrap.sh"
    ;;

  submit)
    tag="${2:-}"
    check_tag "$tag"
    remote_script="$repo/scripts/cluster/gpu_rows.sbatch"
    echo "run.sh: checking out $tag in $repo on $cluster_name ..." >&2
    # shellcheck disable=SC2029
    ssh "$ssh_alias" "cd $(printf %q "$repo") && git fetch --tags --quiet && git checkout --quiet $(printf %q "$tag") && echo \"at \$(git describe --tags --always) (\$(git rev-parse --short HEAD))\" && cat $(printf %q "$remote_script")"
    echo >&2
    echo "run.sh: the above is the script that will be submitted (modules: '${modules:-none}')." >&2
    confirm "Submit it for $tag on $cluster_name?"
    log_pattern="$fred/ampere/logs/gpu-$tag-%j.log"
    # shellcheck disable=SC2029
    out="$(ssh "$ssh_alias" "cd $(printf %q "$repo") && sbatch --parsable --output=$(printf %q "$log_pattern") --export=ALL,AMPERE_REPO=$(printf %q "$repo"),AMPERE_PIXI_HOME=$(printf %q "$pixi_home"),AMPERE_TAG=$(printf %q "$tag"),AMPERE_CLUSTER_NAME=$(printf %q "$cluster_name"),AMPERE_MODULES=$(printf %q "$modules") $(printf %q "$remote_script")")"
    jobid="${out%%;*}"
    mkdir -p "$logdir"
    printf '%s %s tag=%s job=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$cluster_name" "$tag" "$jobid" >> "$logdir/submissions.txt"
    echo "run.sh: submitted job $jobid for $tag (logged in $logdir/submissions.txt)"
    echo "  watch it:  ssh $ssh_alias sacct -j $jobid --format=JobID,State,Elapsed,ExitCode"
    echo "  then:      scripts/cluster/run.sh fetch $tag"
    echo "  (do not watch squeue in a loop -- site policy)"
    ;;

  fetch)
    tag="${2:-}"
    check_tag "$tag"
    mkdir -p "$logdir/raw"
    echo "run.sh: fetching gpu-$tag-*.log through $mover ..." >&2
    rsync -t "$mover:$fred/ampere/logs/gpu-$tag-*.log" "$logdir/raw/"
    # shellcheck disable=SC2012  # file names are ours: gpu-<tag>-<jobid>.log
    newest="$(ls -t "$logdir"/raw/gpu-"$tag"-*.log | head -n 1)"
    cp "$newest" "$logdir/$tag.log"
    echo "run.sh: $logdir/$tag.log  (from $(basename "$newest"))"
    echo "--- last lines (the summary line is the status row's quote):"
    tail -n 3 "$logdir/$tag.log"
    ;;

  *) usage ;;
esac
