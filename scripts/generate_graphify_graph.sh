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

if ! "$GRAPHIFY_BIN" extract "$CORPUS" --out "$RUN_OUT" --no-cluster; then
  echo "[graphify extract] failed; falling back to code-only graphify update" >&2
  "$GRAPHIFY_BIN" update "$CORPUS" --no-cluster
  rm -rf "$RUN_OUT/graphify-out"
  mkdir -p "$RUN_OUT"
  cp -R "$CORPUS/graphify-out" "$RUN_OUT/graphify-out"
fi

# Structure-only: drop the inferred/semantic (confidence=INFERRED) and the
# rationale/concept/document layers so the map reflects ARIA's REAL code
# structure and relationships (deterministic, no LLM). Runs on the raw
# extraction BEFORE clustering so the report/tree/html are built from the
# filtered graph.
python "$ROOT/scripts/graphify_structure_filter.py" \
  "$RUN_OUT/graphify-out/graph.json"

if [[ ! -f "$RUN_OUT/graphify-out/manifest.json" ]]; then
  python - "$RUN_OUT/graphify-out/graph.json" \
    "$RUN_OUT/graphify-out/manifest.json" "$CORPUS" <<'PY'
from pathlib import Path
import hashlib
import json
import sys

graph_path = Path(sys.argv[1])
manifest_path = Path(sys.argv[2])
corpus = Path(sys.argv[3])
graph = json.loads(graph_path.read_text(encoding="utf-8"))
files = sorted({
    n.get("source_file")
    for n in graph.get("nodes", [])
    if n.get("file_type") == "code" and n.get("source_file")
})
manifest = {}
for rel in files:
    path = corpus / rel
    if not path.is_file():
        continue
    data = path.read_bytes()
    digest = hashlib.md5(data).hexdigest()
    manifest[str(path)] = {
        "mtime": path.stat().st_mtime,
        "ast_hash": digest,
        "semantic_hash": digest,
    }
manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
PY
fi

# cluster-only regenerates the human-readable report + html (god-files,
# communities) FROM the structure-only graph.json. graphify computes communities
# at report time and does not persist them per-node, so graph.json stays the
# complete, clean structural graph (its "fewer nodes than before" safety refusal
# is expected here and harmless — the filter intentionally shrank the graph).
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
