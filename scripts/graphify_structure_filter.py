#!/usr/bin/env python
"""Reduce a Graphify graph.json to ARIA's REAL code structure only.

Graphify's raw extraction mixes three layers:
  1. real code structure  — file_type=="code" nodes + EXTRACTED structural edges
     (imports/imports_from/calls/contains/method/inherits/defines/re_exports/
     implements/references);
  2. an inferred/semantic layer — confidence=="INFERRED" edges (`uses`, some
     `calls`/`references`);
  3. a rationale/doc layer — `rationale`/`concept`/`document` nodes + the
     `rationale_for` edges synthesised from comments/docstrings.

This filter keeps ONLY layer 1 so the map reflects what the code actually is and
how it actually connects — no inferred or LLM-derived nodes/edges. It is purely
deterministic (no model calls). Edges whose endpoints were dropped are removed;
LLM token counters and hyperedges are zeroed.

Usage: graphify_structure_filter.py <graph.json> [--out <path>]
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter

# EXTRACTED structural relations that reflect real code structure (Option A:
# structure + references).
STRUCTURAL_RELATIONS = {
    "contains", "imports", "imports_from", "calls", "method",
    "inherits", "defines", "re_exports", "implements", "references",
}


def filter_graph(graph: dict) -> tuple[dict, dict]:
    nodes = graph.get("nodes", [])
    edges = graph.get("edges", [])

    keep_ids = {n["id"] for n in nodes if n.get("file_type") == "code"}
    new_nodes = [n for n in nodes if n["id"] in keep_ids]
    new_edges = [
        e for e in edges
        if e.get("confidence") == "EXTRACTED"
        and e.get("relation") in STRUCTURAL_RELATIONS
        and e.get("source") in keep_ids
        and e.get("target") in keep_ids
    ]

    stats = {
        "nodes_before": len(nodes),
        "nodes_after": len(new_nodes),
        "nodes_dropped_by_type": dict(Counter(
            n.get("file_type") for n in nodes if n["id"] not in keep_ids
        )),
        "edges_before": len(edges),
        "edges_after": len(new_edges),
        "edges_dropped_inferred": sum(
            1 for e in edges if e.get("confidence") == "INFERRED"
        ),
        "edges_kept_by_relation": dict(Counter(e["relation"] for e in new_edges)),
    }

    out = dict(graph)
    out["nodes"] = new_nodes
    out["edges"] = new_edges
    out["hyperedges"] = []          # drop the inferred hyperedge layer
    out["input_tokens"] = 0         # no LLM layer in a structure-only graph
    out["output_tokens"] = 0
    out["structure_only"] = True    # provenance stamp
    return out, stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("graph", help="path to graph.json")
    ap.add_argument("--out", default=None, help="output path (default: in place)")
    args = ap.parse_args(argv)

    graph = json.loads(open(args.graph, encoding="utf-8").read())
    filtered, stats = filter_graph(graph)

    out_path = args.out or args.graph
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(filtered, fh, indent=2)

    print("Graphify structure-only filter:")
    print(f"  nodes: {stats['nodes_before']} -> {stats['nodes_after']} "
          f"(dropped non-code: {stats['nodes_dropped_by_type']})")
    print(f"  edges: {stats['edges_before']} -> {stats['edges_after']} "
          f"(dropped INFERRED: {stats['edges_dropped_inferred']})")
    print(f"  kept edge relations: {stats['edges_kept_by_relation']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
