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
import ast
import json
import re
import sys
from collections import Counter
from pathlib import Path

# EXTRACTED structural relations that reflect real code structure (Option A:
# structure + references).
STRUCTURAL_RELATIONS = {
    "contains", "imports", "imports_from", "calls", "method",
    "inherits", "defines", "re_exports", "implements", "references",
}


def _slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def ensure_python_class_members(
    graph: dict,
    source_root: str | Path,
    source_file: str,
    class_name: str,
) -> dict[str, int]:
    """Deterministically add AST-extracted members omitted by Graphify.

    Graphify occasionally emits only a repository-level class node for shared
    infrastructure modules. This narrow fallback preserves real Python structure:
    it parses a named class, reuses existing nodes when present, and adds only
    EXTRACTED ``contains``/``method`` edges. No semantic or inferred layer enters
    the committed map.
    """
    root = Path(source_root)
    path = root / source_file
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=source_file)
    class_ast = next(
        (node for node in ast.walk(tree)
         if isinstance(node, ast.ClassDef) and node.name == class_name),
        None,
    )
    if class_ast is None:
        raise ValueError(f"class {class_name!r} not found in {source_file}")

    nodes = graph.setdefault("nodes", [])
    edges = graph.setdefault("edges", [])
    used_ids = {str(node.get("id")) for node in nodes}

    def unique_id(base: str) -> str:
        candidate = base
        suffix = 2
        while candidate in used_ids:
            candidate = f"{base}_{suffix}"
            suffix += 1
        used_ids.add(candidate)
        return candidate

    module = next(
        (node for node in nodes
         if node.get("source_file") == source_file
         and node.get("label") == Path(source_file).name),
        None,
    )
    module_added = 0
    if module is None:
        module = {
            "id": unique_id(f"python_{_slug(source_file)}"),
            "label": Path(source_file).name,
            "file_type": "code",
            "source_file": source_file,
        }
        nodes.append(module)
        module_added = 1

    class_node = next(
        (node for node in nodes
         if node.get("source_file") == source_file
         and node.get("label") == class_name),
        None,
    )
    class_added = 0
    if class_node is None:
        class_node = {
            "id": unique_id(f"python_{_slug(source_file)}_{_slug(class_name)}"),
            "label": class_name,
            "file_type": "code",
            "source_file": source_file,
        }
        nodes.append(class_node)
        class_added = 1

    def add_edge(source: str, target: str, relation: str, line: int) -> None:
        if any(
            edge.get("source") == source
            and edge.get("target") == target
            and edge.get("relation") == relation
            for edge in edges
        ):
            return
        edges.append({
            "source": source,
            "target": target,
            "relation": relation,
            "confidence": "EXTRACTED",
            "source_file": source_file,
            "source_location": f"L{line}",
            "weight": 1.0,
        })

    add_edge(module["id"], class_node["id"], "contains", class_ast.lineno)
    methods_added = 0
    for member in class_ast.body:
        if not isinstance(member, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        label = f".{member.name}()"
        method = next(
            (node for node in nodes
             if node.get("source_file") == source_file
             and node.get("label") in {label, f"{class_name}.{member.name}()"}),
            None,
        )
        if method is None:
            method = {
                "id": unique_id(f"{class_node['id']}_{_slug(member.name)}"),
                "label": label,
                "file_type": "code",
                "source_file": source_file,
            }
            nodes.append(method)
            methods_added += 1
        add_edge(class_node["id"], method["id"], "method", member.lineno)

    return {
        "module_added": module_added,
        "class_added": class_added,
        "methods_added": methods_added,
    }


def filter_graph(graph: dict) -> tuple[dict, dict]:
    nodes = graph.get("nodes", [])
    # `graphify extract` writes `edges`; `graphify update` (code-only, no LLM)
    # writes the same structural records under `links`. Normalize both so the
    # structure-only reducer can run deterministically even when the semantic
    # extractor is unavailable.
    edges = graph.get("edges")
    if edges is None:
        edges = graph.get("links", [])

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
    out.pop("links", None)
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
    ap.add_argument(
        "--source-root",
        default=None,
        help="tracked source root used by deterministic AST enrichments",
    )
    ap.add_argument(
        "--ensure-class-members",
        action="append",
        default=[],
        metavar="PATH:CLASS",
        help="ensure real AST members for a class Graphify omitted",
    )
    args = ap.parse_args(argv)

    graph = json.loads(open(args.graph, encoding="utf-8").read())
    filtered, stats = filter_graph(graph)
    enrichments = []
    if args.ensure_class_members and not args.source_root:
        ap.error("--source-root is required with --ensure-class-members")
    for spec in args.ensure_class_members:
        try:
            source_file, class_name = spec.rsplit(":", 1)
        except ValueError:
            ap.error(f"invalid --ensure-class-members value: {spec!r}")
        enrichment = ensure_python_class_members(
            filtered, args.source_root, source_file, class_name
        )
        enrichments.append((spec, enrichment))

    out_path = args.out or args.graph
    with open(out_path, "w", encoding="utf-8") as fh:
        json.dump(filtered, fh, indent=2)

    print("Graphify structure-only filter:")
    print(f"  nodes: {stats['nodes_before']} -> {stats['nodes_after']} "
          f"(dropped non-code: {stats['nodes_dropped_by_type']})")
    print(f"  edges: {stats['edges_before']} -> {stats['edges_after']} "
          f"(dropped INFERRED: {stats['edges_dropped_inferred']})")
    print(f"  kept edge relations: {stats['edges_kept_by_relation']}")
    for spec, enrichment in enrichments:
        print(f"  AST class enrichment {spec}: {enrichment}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
