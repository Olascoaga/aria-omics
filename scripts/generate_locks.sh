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
#   scripts/generate_locks.sh                 # snapshot the published-env whitelist
#   scripts/generate_locks.sh aria-rna-env    # snapshot one env
#   scripts/generate_locks.sh --all           # snapshot every env in envs/
#   scripts/generate_locks.sh --requirements  # only refresh the top-level
#                                             #   requirements.lock (pip core)
#
# P2-1: the published runtime envs (RNA / ATAC=chromatin / integration) are
# conda-managed; their linux-64 explicit locks are SNAPSHOTS of the installed
# env (conda-lock's solver hangs on the bioconda+pip mix here). Multi-platform
# locks and locks for envs that are not installed on this machine are produced
# by the hermetic Docker lane (P2-2), not fabricated here. The pip core fallback
# is requirements.lock (orchestrator/aria-env).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENVS_DIR="$ROOT/envs"

# Published runtime envs (P2-1). aria-rna-env is the validated RNA baseline;
# ingestion owns raw scRNA FASTQ quantification; chromatin (scATAC) and
# integration are published as part of the 4.6/4.7 contract. Envs not installed
# locally are skipped with a clear message (their locks come from the Docker
# lane), never fabricated.
WHITELIST=(
  "aria-rna-env"
  "aria-ingestion-env"
  "aria-chromatin-env"
  "aria-integration-env"
)
DEFERRED=(
  "aria-hic-env"   # Hi-C dispatch is OFF (P0-3); env lock is not a release gate
)

# The orchestrator env whose pip freeze backs the top-level requirements.lock.
CORE_ENV="aria-env"

mode="whitelist"
explicit_target=""

if [[ $# -ge 1 ]]; then
  case "$1" in
    --all)
      mode="all"
      ;;
    --requirements)
      mode="requirements"
      ;;
    --help|-h)
      sed -n '2,33p' "$0"
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

# Snapshot the orchestrator env's pip closure into the top-level
# requirements.lock (pip-only fallback). Honest: skipped if CORE_ENV is absent.
snapshot_requirements() {
  local req_file="$ROOT/requirements.lock"
  if ! env_installed "$CORE_ENV"; then
    echo "skip requirements.lock: $CORE_ENV not installed locally" >&2
    return 0
  fi
  echo "==> snapshotting $CORE_ENV pip closure -> requirements.lock"
  local py_ver freeze
  py_ver="$(conda run --name "$CORE_ENV" --no-capture-output python --version 2>/dev/null | awk '{print $2}')"
  freeze="$(conda run --name "$CORE_ENV" --no-capture-output pip freeze 2>/dev/null \
            | grep -ivE '^-e |^aria-omics| @ file://' || true)"
  {
    echo "# ARIA core/orchestrator pip lock (P2-1) — pip fallback for the aria-env."
    echo "# Fully-pinned snapshot of the validated orchestrator environment so a"
    echo "# pip-only install is reproducible. Provenance: pip freeze of conda env"
    echo "# '$CORE_ENV' on linux-64, Python ${py_ver:-unknown}, $(date +%Y-%m-%d)."
    echo "#"
    echo "# This covers the CORE (LLM/orchestration/IO) stack. The heavy scientific"
    echo "# runtime envs (RNA/ATAC/integration) are conda-managed and locked"
    echo "# separately under envs/<env>.linux-64.lock (+ .pip.lock); multi-platform"
    echo "# and non-installed-env locks are produced by the hermetic Docker lane."
    echo "# Regenerate with: scripts/generate_locks.sh --requirements"
    echo "#"
    printf "%s\n" "$freeze"
  } > "$req_file"
}

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
      echo "skip $env_name (deferred; rerun with: scripts/generate_locks.sh $env_name)" >&2
    done
    snapshot_requirements
    ;;
  requirements)
    snapshot_requirements
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
