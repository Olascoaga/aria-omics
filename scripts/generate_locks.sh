#!/usr/bin/env bash
#
# Snapshot the conda envs that ARIA actually depends on, as explicit
# lockfiles, so a peer reviewer can reproduce the run with:
#
#   conda create --name <env> --file envs/<env>.linux-64.lock
#   pip install -r envs/<env>.pip.lock     # if a pip lockfile exists
#
# Strategy:
#   - Read the active lock targets from aria.utils.environment_specs.
#   - For every required env that is INSTALLED on this machine,
#     run `conda list --explicit` against that env and write
#     envs/<env>.linux-64.lock. This is what `conda-lock --kind explicit`
#     would write but without invoking the solver, which hangs on
#     bioconda+pip mixes in our setup.
#   - Capture only packages whose conda record channel is `pypi`, normalized as
#     portable exact version pins or immutable VCS/archive references. Never
#     persist build-host file URLs.
#   - Envs that are not installed are reported as deferred, never fabricated.
#
# Usage:
#   scripts/generate_locks.sh                 # snapshot active registry targets
#   scripts/generate_locks.sh aria-rna-env    # snapshot one env
#   scripts/generate_locks.sh --all           # snapshot every env in envs/
#   scripts/generate_locks.sh --requirements  # only refresh the top-level
#                                             #   requirements.lock (pip core)
#
# P2-1/A5: active runtime and benchmark envs are conda-managed; their linux-64
# explicit locks are SNAPSHOTS of the installed
# env (conda-lock's solver hangs on the bioconda+pip mix here). Multi-platform
# locks and locks for envs that are not installed on this machine are produced
# by the hermetic Docker lane (P2-2), not fabricated here. The pip core fallback
# is requirements.lock (orchestrator/aria-env).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENVS_DIR="$ROOT/envs"

PYTHON_BIN="${PYTHON_BIN:-python}"
mapfile -t LOCK_ENVS < <(
  PYTHONPATH="$ROOT" "$PYTHON_BIN" -m aria.utils.environment_specs lock-envs
)

# The orchestrator env whose portable package snapshot backs requirements.lock.
CORE_ENV="aria-env"

mode="registry"
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
  freeze="$(conda run --name "$CORE_ENV" --no-capture-output \
            python -m pip list --format=freeze 2>/dev/null \
            | grep -ivE '^-e |^aria-omics' || true)"
  {
    echo "# ARIA core/orchestrator pip lock (P2-1) — pip fallback for the aria-env."
    echo "# Fully-pinned snapshot of the validated orchestrator environment so a"
    echo "# pip-only install is reproducible. Provenance: package list of conda env"
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
  local lock_tmp="${lock_file}.tmp"
  local records_tmp="${pip_file}.records.tmp"
  local pip_tmp="${pip_file}.tmp"

  if ! env_installed "$env_name"; then
    echo "skip $env_name: env not installed locally (snapshot deferred)" >&2
    return 0
  fi

  echo "==> snapshotting $env_name -> $lock_file"
  # `conda list --explicit` produces the exact same format conda-lock
  # writes with --kind explicit, plus an @EXPLICIT marker. It does not
  # invoke the solver.
  conda list --name "$env_name" --explicit > "$lock_tmp"
  mv "$lock_tmp" "$lock_file"

  # Pip's traditional environment snapshot can preserve build-host references such as
  # file:///tmp/... . Conda's JSON marks true PyPI installs explicitly, so use
  # that authoritative channel and normalize it through the shared registry.
  conda list --name "$env_name" --json > "$records_tmp"
  PYTHONPATH="$ROOT" conda run --name "$env_name" --no-capture-output \
    python -m aria.utils.environment_specs pip-lock < "$records_tmp" > "$pip_tmp"
  rm -f "$records_tmp"
  if [[ -s "$pip_tmp" ]]; then
    mv "$pip_tmp" "$pip_file"
    echo "    + portable pip lock: $pip_file"
  else
    rm -f "$pip_tmp" "$pip_file"
  fi
}

case "$mode" in
  registry)
    for env_name in "${LOCK_ENVS[@]}"; do
      snapshot_one "$env_name"
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
