# ARIA Graphify Map (structure-only)

This directory contains the generated Graphify map for the tracked ARIA
repository snapshot. It is a **structure-only** graph: it reflects ARIA's REAL
code structure and relationships — no inferred or LLM-derived nodes/edges.

## What "structure-only" means

`scripts/graphify_structure_filter.py` runs on the raw extraction and keeps ONLY:

- **nodes** with `file_type == "code"` (files, classes, functions, methods);
- **edges** with `confidence == "EXTRACTED"` and a structural relation:
  `contains`, `imports`, `imports_from`, `calls`, `method`, `inherits`,
  `defines`, `re_exports`, `implements`, `references`.

It **drops**: every `confidence == "INFERRED"` edge (`uses` and a few inferred
`calls`/`references`), the `rationale_for` edges, and all `rationale` / `concept`
/ `document` nodes; it zeroes the LLM token counters and the hyperedge layer.
The graph carries a `"structure_only": true` provenance stamp.

One deterministic AST completeness guard restores the real members of
`aria/memory/memory.py:ARIAMemory` when the upstream extractor emits only an
aggregate class node. The added nodes and `method` edges come directly from the
tracked Python AST; no inferred/semantic layer is introduced. In particular,
`.export_experiment_snapshot()` must be queryable in every generated graph.

## Snapshot

- Commit: `c19c2e5` (HEAD the graph was built from)
- Generated: 2026-07-18
- Corpus: tracked repository files only, generated from `git archive HEAD`
- Private operational memory (`memory/`) and local agent settings are excluded
- Structure-only graph: **9062 code nodes / 17641 EXTRACTED structural edges** in `graph.json`.
- Clustered report: **9062 nodes / 17641 edges / 593 communities** in `GRAPH_REPORT.md` (581 shown, 12 thin omitted).
- No LLM extraction or labeling: generation uses the local AST-only update path,
  token counters remain at zero, and `--no-label` keeps "Community N"
  placeholders.

## Files

- `graph.json`: GraphRAG/queryable graph.
- `GRAPH_TREE.html`: interactive collapsible file/symbol tree.
- `GRAPH_REPORT.md`: generated hub/community report.
- `manifest.json`: source-file hashes for the clean tracked snapshot.

The force-directed `graph.html` is intentionally absent: the current graph is
larger than Graphify's safe visualization limit. The generator removes an older
copy instead of presenting a stale visualization.

## Use

Consult this generated map when a change spans multiple agents, scripts,
narrators, tests, or documentation surfaces. Use it together with the curated
impact map in `../code_graph.md`: Graphify helps discover broad relationships;
the curated map remains the release-risk checklist.

```bash
graphify explain "BulkRNAAgent" \
  --graph docs/architecture/graphify/graph.json

graphify query "What touches pseudobulk DE?" \
  --graph docs/architecture/graphify/graph.json --budget 2000

graphify path "scRNAAgent" "rna_pseudobulk_de.py" \
  --graph docs/architecture/graphify/graph.json
```

Open `GRAPH_TREE.html` for an interactive, file-oriented navigation surface.

## Regenerate

Use the repository script so the graph is built from a clean tracked snapshot
instead of local private files. The script runs `graphify update --force
--no-cluster` (local AST extraction) → `scripts/graphify_structure_filter.py` →
`graphify cluster-only` (report and clustering) → `graphify tree`:

```bash
scripts/generate_graphify_graph.sh
```

The structure-only filter and post-processing normalize Graphify's code graph
to the stable `edges` schema used by the committed artifact.
