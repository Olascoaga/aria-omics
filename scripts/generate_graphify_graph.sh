#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
OUT_DIR="${OUT_DIR:-$ROOT/docs/architecture/graphify}"
GRAPHIFY_BIN="${GRAPHIFY_BIN:-graphify}"

if ! command -v "$GRAPHIFY_BIN" >/dev/null 2>&1; then
  if [[ -x /home/medusa/anaconda3/bin/graphify ]]; then
    GRAPHIFY_BIN=/home/medusa/anaconda3/bin/graphify
  else
    echo "graphify not found. Install graphifyy or set GRAPHIFY_BIN=/path/to/graphify." >&2
    exit 1
  fi
fi

TMP_ROOT="$(mktemp -d /tmp/aria_graphify.XXXXXX)"
trap 'rm -rf "$TMP_ROOT"' EXIT

CORPUS="$TMP_ROOT/corpus"
RUN_OUT="$TMP_ROOT/out"
mkdir -p "$CORPUS" "$RUN_OUT" "$OUT_DIR"

git -C "$ROOT" archive --format=tar --output="$TMP_ROOT/repo.tar" HEAD -- \
  . \
  ':(exclude)docs/architecture/graphify/GRAPH_REPORT.md' \
  ':(exclude)docs/architecture/graphify/GRAPH_TREE.html' \
  ':(exclude)docs/architecture/graphify/graph.html' \
  ':(exclude)docs/architecture/graphify/graph.json' \
  ':(exclude)docs/architecture/graphify/manifest.json'
tar -xf "$TMP_ROOT/repo.tar" -C "$CORPUS"

"$GRAPHIFY_BIN" extract "$CORPUS" --out "$RUN_OUT" --no-cluster
"$GRAPHIFY_BIN" cluster-only "$RUN_OUT" --no-label

cp "$RUN_OUT/graphify-out/graph.json" "$OUT_DIR/graph.json"
cp "$RUN_OUT/graphify-out/graph.html" "$OUT_DIR/graph.html"
cp "$RUN_OUT/graphify-out/GRAPH_REPORT.md" "$OUT_DIR/GRAPH_REPORT.md"
cp "$RUN_OUT/graphify-out/manifest.json" "$OUT_DIR/manifest.json"

"$GRAPHIFY_BIN" tree \
  --graph "$OUT_DIR/graph.json" \
  --output "$OUT_DIR/GRAPH_TREE.html" \
  --root "$CORPUS" \
  --label ARIA

python - "$OUT_DIR" "$CORPUS" "$RUN_OUT" "$ROOT" <<'PY'
from pathlib import Path
import sys

out_dir = Path(sys.argv[1])
corpus = sys.argv[2]
run_out = sys.argv[3]
root = sys.argv[4]

for path in out_dir.iterdir():
    if path.is_file() and path.suffix in {".json", ".md", ".html"}:
        text = path.read_text(encoding="utf-8")
        text = text.replace(corpus, root).replace(run_out, root)
        path.write_text(text, encoding="utf-8")
PY

echo "Graphify graph written to $OUT_DIR"
