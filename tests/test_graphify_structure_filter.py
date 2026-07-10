"""Deterministic Graphify structure enrichment guards."""
from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "graphify_structure_filter.py"


def _module():
    spec = importlib.util.spec_from_file_location("graphify_structure_filter", SCRIPT)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_ensure_class_members_adds_ariamemory_export_method():
    module = _module()
    graph = {
        "nodes": [{
            "id": "aria_memory",
            "label": "ARIAMemory",
            "file_type": "code",
            "source_file": "aria/memory/memory.py",
        }],
        "edges": [],
    }

    stats = module.ensure_python_class_members(
        graph, ROOT, "aria/memory/memory.py", "ARIAMemory"
    )

    methods = {
        node["label"]: node["id"]
        for node in graph["nodes"]
        if node.get("source_file") == "aria/memory/memory.py"
    }
    method_id = methods[".export_experiment_snapshot()"]
    assert any(
        edge.get("source") == "aria_memory"
        and edge.get("target") == method_id
        and edge.get("relation") == "method"
        and edge.get("confidence") == "EXTRACTED"
        for edge in graph["edges"]
    )
    assert stats["methods_added"] >= 1
