# ARIA Graphify Map

This directory contains the generated Graphify map for the tracked ARIA
repository snapshot.

## Snapshot

- Commit: `d7b97144`
- Generated: 2026-05-31
- Corpus: tracked repository files only, generated from `git archive HEAD`
- Private operational memory (`memory/`) and local agent settings are excluded
- Graph: 2299 extracted nodes / 6732 extracted edges in `graph.json`
- Clustered report: 2296 nodes / 5428 edges / 101 communities in `GRAPH_REPORT.md`
- Semantic extraction cost: 28,388 input tokens / 2,116 output tokens, about `$0.0205`

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
/home/medusa/anaconda3/bin/graphify explain "BulkRNAAgent" \
  --graph docs/architecture/graphify/graph.json

/home/medusa/anaconda3/bin/graphify query "What touches pseudobulk DE?" \
  --graph docs/architecture/graphify/graph.json --budget 2000

/home/medusa/anaconda3/bin/graphify path "scRNAAgent" "rna_pseudobulk_de.py" \
  --graph docs/architecture/graphify/graph.json
```

Open `graph.html` for the dense relationship map and `GRAPH_TREE.html` for a
more file-oriented navigation surface.

## Regenerate

Use the repository script so the graph is built from a clean tracked snapshot
rather than local private files:

```bash
scripts/generate_graphify_graph.sh
```

Graphify currently emits NumPy 2 compatibility warnings from optional
`pandas`/`pyarrow`/`numexpr` imports during report generation in the base Conda
environment. The warning is noisy but the graph, report, and HTML are produced.
