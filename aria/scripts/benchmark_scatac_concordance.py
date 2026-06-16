"""scATAC P3a — external concordance benchmark runner (ARIA vs SnapATAC2/...).

Dispatchable IPC runner for the external-tool concordance benchmark
(master-plan W3.1). It loads ARIA's scATAC outputs and a reference tool's outputs
(both as plain JSON — ``clusters`` / ``da_peaks`` / ``motifs``), scores them with
``aria.benchmarks.concordance_atac``, and writes a manifest.

No-fabrication contract (ADR-002): if the external tool's outputs are absent the
runner returns ``status="not_run"`` / ``error_type="MissingDependency"`` instead
of inventing a perfect concordance. P3a is buildable now because it CONSUMES the
comparator's outputs; the infra-gated P3b drop-in is the SnapATAC2/ArchR/Signac
driver that PRODUCES them inside ``aria-bench-env`` (installing those tools +
public datasets is Samael's call and is not done here).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from aria.benchmarks.concordance_atac import score_atac_concordance
from aria.scripts._base import run_script


def _load_outputs(path_str: str | None) -> dict[str, Any] | None:
    if not path_str:
        return None
    path = Path(path_str)
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def benchmark_scatac_concordance(params: dict) -> dict:
    tool = str(params.get("tool", "external"))
    aria = _load_outputs(params.get("aria_outputs_json"))
    external = _load_outputs(params.get("external_outputs_json"))

    if aria is None:
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": "ARIA scATAC outputs JSON missing or unreadable "
                      "(expected clusters/da_peaks/motifs).",
        }
    if external is None:
        return {
            "status": "not_run",
            "error_type": "MissingDependency",
            "tool": tool,
            "reason": (
                f"{tool} outputs JSON missing or unreadable. The external "
                "comparator must run in aria-bench-env (P3b, infra-gated) before "
                "concordance can be scored; ARIA does not fabricate a score."),
        }

    manifest = score_atac_concordance(
        aria, external, tool=tool,
        top_k_motifs=params.get("top_k_motifs"),
        slack=int(params.get("slack", 0)),
    )
    manifest["benchmark"] = "P3a_scatac_external_concordance"
    manifest["benchmark_version"] = "v1"
    manifest["scope"] = "external_tool_concordance"
    manifest.setdefault("messages", []).append(
        "P3a scores ARIA scATAC outputs against an external tool's outputs "
        "(cluster ARI/NMI, DA-peak genomic overlap, motif top-k/RBO). It is "
        "execution/concordance evidence, NOT a de-alpha of the scATAC modality."
    )

    output_dir = params.get("output_dir")
    if output_dir:
        from aria.version import __version__, collect_version_metadata
        outdir = Path(output_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        manifest["aria_version"] = __version__
        manifest["provenance"] = collect_version_metadata()
        manifest_name = str(params.get(
            "manifest_name", f"p3a_concordance_{tool.lower()}.json"))
        manifest_path = outdir / manifest_name
        manifest.setdefault("artifacts", {})["manifest_json"] = str(manifest_path)
        # T10 (tri-audit 2026-06-14): keep machine-absolute paths out of public
        # benchmark JSONs. Reuse the A1 relativizer (single implementation).
        from aria.scripts.benchmark_a1_external_comparators import _relativize_paths
        manifest = _relativize_paths(manifest)
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return manifest


if __name__ == "__main__":
    run_script(benchmark_scatac_concordance)
