#!/usr/bin/env bash
set -euo pipefail

ROOT="$(git rev-parse --show-toplevel)"
OUT_DIR="${OUT_DIR:-$ROOT/docs/architecture/graphify}"
GRAPHIFY_BIN="${GRAPHIFY_BIN:-graphify}"

if ! command -v "$GRAPHIFY_BIN" >/dev/null 2>&1; then
  if [[ -x "$HOME/anaconda3/bin/graphify" ]]; then
    GRAPHIFY_BIN="$HOME/anaconda3/bin/graphify"
  else
    echo "graphify not found. Install graphify or set GRAPHIFY_BIN=/path/to/graphify." >&2
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

# Copy each artifact only if it was produced. `cluster-only` skips graph.html when the
# graph exceeds the viz node limit; under `set -e` a hard `cp graph.html` would abort
# the script BEFORE the path-scrub + README stamp below (this is why the absolute-path
# leak survived earlier regens). A skipped graph.html keeps its previous committed copy,
# which the scrub step still cleans.
for _f in graph.json graph.html GRAPH_REPORT.md manifest.json; do
  if [[ -f "$RUN_OUT/graphify-out/$_f" ]]; then
    cp "$RUN_OUT/graphify-out/$_f" "$OUT_DIR/$_f"
  fi
done

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

# F11: committed Graphify artifacts must be relocatable and must NEVER embed a
# machine-absolute path. The temp corpus/run_out roots (and the absolute repo root)
# previously leaked into manifest/report/html as `/home/<user>/...` on every regen.
# Map them to a repo-relative form: "<root>/aria/x.py" -> "aria/x.py"; a bare
# "<root>" (e.g. the report title) -> the stable repo label "ARIA".
abs_roots = [r.rstrip("/") for r in (corpus, run_out, root)]
for path in out_dir.iterdir():
    if path.is_file() and path.suffix in {".json", ".md", ".html"}:
        text = path.read_text(encoding="utf-8")
        for ar in abs_roots:
            text = text.replace(ar + "/", "")
        for ar in abs_roots:
            text = text.replace(ar, "ARIA")
        path.write_text(text, encoding="utf-8")
PY

# F11: keep the README snapshot pointer fresh (it had drifted to a stale commit).
# Stamp the HEAD the graph was built from + the date so it never goes stale on regen.
HEAD_SHORT="$(git -C "$ROOT" rev-parse --short HEAD)"
TODAY="$(date +%Y-%m-%d)"
README="$OUT_DIR/README.md"
if [[ -f "$README" ]]; then
  python - "$README" "$HEAD_SHORT" "$TODAY" <<'PY'
import re
import sys
from pathlib import Path

readme, head, today = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
text = readme.read_text(encoding="utf-8")
text = re.sub(r"(?m)^- Commit: `[0-9a-f]+`.*$",
              f"- Commit: `{head}` (HEAD the graph was built from)", text)
text = re.sub(r"(?m)^- Generated: .*$", f"- Generated: {today}", text)
readme.write_text(text, encoding="utf-8")
PY
fi

echo "Graphify graph written to $OUT_DIR"
