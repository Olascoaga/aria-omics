#!/usr/bin/env bash
#
# Snapshot the conda envs that ARIA actually depends on, as explicit
# lockfiles, so a peer reviewer can reproduce the run with:
#
#   conda create --name <env> --file envs/<env>.linux-64.lock
#   pip install -r envs/<env>.pip.lock     # if a pip lockfile exists
#
# Strategy:
#   - For every env in the whitelist that is INSTALLED on this machine,
#     run `conda list --explicit` against that env and write
#     envs/<env>.linux-64.lock. This is what `conda-lock --kind explicit`
#     would write but without invoking the solver, which hangs on
#     bioconda+pip mixes in our setup.
#   - Capture pip side as envs/<env>.pip.lock via `pip freeze` from the
#     same env. Skipped if `pip` reports nothing.
#   - Envs that are not installed are reported as deferred — the v4.4
#     report will surface that explicitly. scATAC env is deferred to v4.5
#     by design.
#
# Usage:
#   scripts/generate_locks.sh                 # snapshot the v4.4 whitelist
#   scripts/generate_locks.sh aria-rna-env    # snapshot one env
#   scripts/generate_locks.sh --all           # snapshot every env in envs/
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENVS_DIR="$ROOT/envs"

# Whitelist for v4.4. v4.4 Stage C uses aria-rna-env only; hic and
# integration are listed because they are part of the v4.4 contract for
# any future bulk-ATAC / integration rerun. If they are not installed
# here they will be skipped with a clear message.
WHITELIST=(
  "aria-rna-env"
  "aria-hic-env"
  "aria-integration-env"
)
DEFERRED=(
  "aria-chromatin-env"   # deferred to v4.5 scATAC sprint
)

mode="whitelist"
explicit_target=""

if [[ $# -ge 1 ]]; then
  case "$1" in
    --all)
      mode="all"
      ;;
    --help|-h)
      sed -n '2,28p' "$0"
      exit 0
      ;;
    *)
      mode="explicit"
      explicit_target="$1"
      ;;
  esac
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found." >&2
  exit 1
fi

env_installed() {
  conda env list 2>/dev/null \
    | awk 'NF && $1 !~ /^#/ {print $1}' \
    | grep -Fxq "$1"
}

snapshot_one() {
  local env_name="$1"
  local lock_file="$ENVS_DIR/${env_name}.linux-64.lock"
  local pip_file="$ENVS_DIR/${env_name}.pip.lock"

  if ! env_installed "$env_name"; then
    echo "skip $env_name: env not installed locally (snapshot deferred)" >&2
    return 0
  fi

  echo "==> snapshotting $env_name -> $lock_file"
  # `conda list --explicit` produces the exact same format conda-lock
  # writes with --kind explicit, plus an @EXPLICIT marker. It does not
  # invoke the solver.
  conda list --name "$env_name" --explicit > "$lock_file"

  # pip side. If the env has no pip packages, skip the .pip.lock file
  # rather than emit an empty one — keeps the envs/ dir clean.
  local pip_path
  pip_path="$(conda run --name "$env_name" --no-capture-output which pip 2>/dev/null || true)"
  if [[ -n "$pip_path" ]]; then
    local pip_dump
    pip_dump="$(conda run --name "$env_name" --no-capture-output pip freeze 2>/dev/null || true)"
    if [[ -n "$pip_dump" ]]; then
      printf "%s\n" "$pip_dump" > "$pip_file"
      echo "    + pip lock: $pip_file"
    fi
  fi
}

case "$mode" in
  whitelist)
    for env_name in "${WHITELIST[@]}"; do
      snapshot_one "$env_name"
    done
    for env_name in "${DEFERRED[@]}"; do
      echo "skip $env_name (deferred to v4.5; rerun with: scripts/generate_locks.sh $env_name)" >&2
    done
    ;;
  explicit)
    snapshot_one "$explicit_target"
    ;;
  all)
    for env_file in "$ENVS_DIR"/aria-*-env.yml; do
      env_name="$(basename "${env_file%.yml}")"
      snapshot_one "$env_name"
    done
    ;;
esac

echo "done."
