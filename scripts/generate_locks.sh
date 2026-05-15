#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

for env_file in "$ROOT"/envs/aria-*-env.yml; do
  lock_file="${env_file%.yml}.linux-64.lock"
  conda-lock lock \
    --file "$env_file" \
    --kind explicit \
    --platform linux-64 \
    --lockfile "$lock_file"
done

