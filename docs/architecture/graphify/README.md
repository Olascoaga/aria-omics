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

## Snapshot

- Commit: `a14559f` (HEAD the graph was built from)
- Generated: 2026-06-21
- Corpus: tracked repository files only, generated from `git archive HEAD`
- Private operational memory (`memory/`) and local agent settings are excluded
- Structure-only graph: **5027 code nodes / 10441 EXTRACTED structural edges** in
  `graph.json` (the structure-only filter drops non-code nodes and the
  inferred/rationale layers, then keeps only structural edges)
- Clustered report: 5027 nodes / 10441 edges / 303 communities in
  `GRAPH_REPORT.md` (289 shown, 14 thin omitted; community detection ignores a
  few isolated nodes — normal)
- No LLM layer in the final graph: 0 input/output token counters; community naming
  is skipped (`--no-label` → "Community N" placeholders). The upstream extractor
  may still print token logs before the structure-only filter removes rationale
  and inferred layers.

## Files

- `graph.json`: GraphRAG/queryable graph.
- `graph.html`: interactive force graph.
- `GRAPH_TREE.html`: collapsible file/symbol tree.
- `GRAPH_REPORT.md`: generated hub/community report.
- `manifest.json`: source-file hashes for the clean tracked snapshot.

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

Open `graph.html` for the dense relationship map and `GRAPH_TREE.html` for a
more file-oriented navigation surface.

## Regenerate

Use the repository script so the graph is built from a clean tracked snapshot
rather than local private files. The script runs `graphify extract` →
`scripts/graphify_structure_filter.py` (the structure-only reduction) →
`graphify cluster-only` (report/html) → `graphify tree`:

```bash
scripts/generate_graphify_graph.sh
```

Graphify currently emits NumPy 2 compatibility warnings from optional
`pandas`/`pyarrow`/`numexpr` imports during report generation in the base Conda
environment. The warning is noisy but the graph, report, and HTML are produced.
If the semantic extractor hits an API quota, the script falls back to
`graphify update` and the structure-only filter normalizes that code-only
`links` schema into the same `edges` schema used by the committed graph.
