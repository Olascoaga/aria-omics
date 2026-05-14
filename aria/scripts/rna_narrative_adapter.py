"""
ARIA scRNA Narrative Adapter
-----------------------------
Pure-Python helper that converts the JSON report produced by the scRNA E2E
harness (tests/test_scrna_e2e.py) into the agent_results shape that
NarrativeAgent expects.

This is the bridge that lets the harness path drive NarrativeAgent without
running the full Orchestrator + MessageBus stack — useful for offline
report generation on existing pseudobulk / multi-sample workspaces.

Supported modes:
  - "pseudobulk"   (between-condition DE per cell-type via pyDESeq2)
  - "standard"     (per-cluster markers + pathway enrichment)

Outputs a dict with the shape:
    {
      "agent_results": {"scrna_agent": {...}},
      "exp_context":   {organism, genome, modalities, design, ...},
      "intent":        {summary, question, ...},
      "decisions":     [],
      "findings_list": [],
    }

This module intentionally has no LLM, no bus, and no scanpy dependency.
It is safe to import in any environment (only stdlib + pathlib).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional, Union

log = logging.getLogger("aria.narrative_adapter")


def load_e2e_report(path: Union[str, Path]) -> dict:
    """Load and validate an e2e_report.json from the harness workspace."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"e2e report not found: {p}")
    with open(p) as f:
        report = json.load(f)
    if "stages" not in report:
        raise ValueError(
            f"Invalid e2e report at {p}: missing 'stages' key. "
            f"Top-level keys: {list(report.keys())}"
        )
    return report


def adapt(report: Union[dict, str, Path],
          workspace: Optional[Union[str, Path]] = None) -> dict:
    """
    Convert an e2e_report.json (or already-loaded dict) into the
    NarrativeAgent input bundle.

    Args:
        report:     Either a dict (already loaded JSON) or a path to the
                    e2e_report.json file.
        workspace:  Override the workspace path (used as fallback for
                    locating h5ad / CSVs). Defaults to the workspace
                    recorded in the report.

    Returns:
        Dict with keys: agent_results, exp_context, intent, decisions,
        findings_list. Suitable to feed directly into the same arguments
        NarrativeAgent.run() accepts.
    """
    if isinstance(report, (str, Path)):
        report = load_e2e_report(report)

    stages = report.get("stages", {}) or {}
    workspace = Path(workspace or report.get("workspace", "."))
    mode      = report.get("mode", "standard")
    organism  = report.get("organism", "Unknown")
    tissue    = report.get("tissue_hint", "")
    files_in  = report.get("data", []) or []

    findings: dict = {}
    output_h5ad: Optional[str] = None

    # QC ──────────────────────────────────────────────────────────────────
    qc = stages.get("qc")
    if qc and qc.get("status") in ("success", "done"):
        per_sample = qc.get("per_sample", []) or []
        # Aggregate per-sample stats into the shape NarrativeAgent reads.
        n_after = sum((s.get("n_cells_after") or 0) for s in per_sample)
        # `n_cells_before` is not always emitted by rna_qc, but we can
        # reconstruct it from pct_removed when present:
        #     before = after / (1 - pct_removed/100)
        n_before = 0
        for s in per_sample:
            nb = s.get("n_cells_before")
            if nb:
                n_before += nb
            else:
                na = s.get("n_cells_after") or 0
                pr = s.get("pct_removed")
                if pr is not None and pr < 100:
                    n_before += int(round(na / (1.0 - pr / 100.0)))
                else:
                    n_before += na
        pct_rm = (
            100.0 * (n_before - n_after) / n_before
            if n_before else 0.0
        )
        findings["qc"] = {
            "n_cells_before": n_before,
            "n_cells_after":  n_after,
            "pct_removed":    round(pct_rm, 1),
            "n_samples":      qc.get("n_samples", len(per_sample)),
            "per_sample":     per_sample,
            "mt_threshold":   _first_mt_threshold(per_sample),
        }

    # Concat (multi-sample) ───────────────────────────────────────────────
    concat = stages.get("concat")
    if concat and concat.get("status") in ("success", "done"):
        findings["concat"] = {
            "n_samples":      concat.get("n_samples"),
            "n_cells_total":  concat.get("n_cells_total"),
            "n_genes_shared": concat.get("n_genes_shared"),
            "batch_col":      concat.get("batch_col"),
            "output_path":    concat.get("output_path"),
        }
        output_h5ad = concat.get("output_path") or output_h5ad

    # Integration ─────────────────────────────────────────────────────────
    integ = stages.get("integration")
    if integ and integ.get("status") in ("success", "done"):
        findings["integration"] = {
            "status":                  "done",
            "method":                  integ.get("method", "harmony"),
            "n_batches":               integ.get("n_batches"),
            "batch_col":               integ.get("batch_col"),
            "silhouette_before":       integ.get("silhouette_before"),
            "silhouette_after":        integ.get("silhouette_after"),
            "batch_correction_delta":  integ.get("batch_correction_delta"),
            "rep_used":                integ.get("rep_used"),
        }
        output_h5ad = integ.get("output_path") or output_h5ad

    # Clustering ──────────────────────────────────────────────────────────
    clu = stages.get("clustering")
    if clu and clu.get("status") in ("success", "done"):
        findings["clustering"] = {
            "n_clusters":    clu.get("n_clusters"),
            "resolution":    clu.get("resolution"),
            "cluster_sizes": clu.get("cluster_sizes", {}),
            "top_markers":   clu.get("top_markers", {}),
            "rep_used":      clu.get("rep_used"),
            "groupby":       clu.get("groupby"),
            "predef_clusters": clu.get("predef_clusters", False),
        }
        findings["clustering_decision"] = {
            "n_clusters":  clu.get("n_clusters"),
            "recommended": clu.get("resolution"),
            "groupby":     clu.get("groupby"),
            "predef_clusters": clu.get("predef_clusters", False),
            "n_candidates": len(
                (stages.get("advise", {}) or {}).get("candidates", []) or []
            ) or 1,
            "justification":
                "Selected by silhouette score across candidate resolutions.",
        }
        output_h5ad = clu.get("output_path") or output_h5ad

    # CellTypist annotation ───────────────────────────────────────────────
    ct = stages.get("celltypist")
    if ct and ct.get("status") in ("success", "done"):
        per_cl = ct.get("per_cluster", {}) or {}
        # Flatten per-cluster majority labels into the {cluster_id: label} map.
        ctype_map = {
            str(cid): (info.get("majority_label") or info.get("label") or "")
            for cid, info in per_cl.items()
        }
        findings["cell_types"] = {
            "cell_types":      ctype_map,
            "model_used":      ct.get("model_used"),
            "label_col":       ct.get("label_col"),
            "n_unique_labels": ct.get("n_unique_labels"),
            "n_cells":         ct.get("n_cells"),
            "predictions_path": ct.get("predictions_path"),
        }
        output_h5ad = ct.get("output_path") or output_h5ad

    # Per-cluster DE (standard mode) ──────────────────────────────────────
    de = stages.get("de")
    if de and de.get("status") in ("success", "done"):
        findings["differential_expression"] = {
            "n_significant":       de.get("n_significant_total"),
            "n_significant_genes": de.get("n_significant_total"),
            "n_clusters":          de.get("n_clusters"),
            "n_sig_by_cluster":    de.get("n_sig_by_cluster", {}),
            "de_genes_by_cluster": de.get("de_genes_by_cluster", {}),
            "padj_max":            de.get("padj_max"),
            "lfc_min":             de.get("lfc_min"),
            "groupby":             de.get("groupby"),
        }

    # Per-cluster pathways (standard mode) ────────────────────────────────
    pw = stages.get("pathways")
    pw_target_key = "pathways"
    if pw and pw.get("status") in ("success", "done"):
        # When in pseudobulk mode, the pathways stage describes per
        # (group::comparison) enrichment instead of per-cluster markers.
        if mode == "pseudobulk":
            pw_target_key = "pseudobulk_pathways"
        findings[pw_target_key] = {
            "organism":    pw.get("organism", organism),
            "databases":   pw.get("databases", {}),
            "per_cluster": pw.get("per_cluster", {}),
        }

    # Cell-cell communication (v4.3.7) ────────────────────────────────────
    ccc = stages.get("cell_communication") or stages.get("cellcomm")
    if ccc and ccc.get("status") in ("success", "done"):
        findings["cell_communication"] = {
            "status":             "done",
            "method":             ccc.get("method"),
            "n_cell_types":       ccc.get("n_cell_types"),
            "n_interactions":     ccc.get("n_interactions"),
            "n_autocrine_dropped": ccc.get("n_autocrine_dropped", 0),
            "top_interactions":   ccc.get("top_interactions", []) or [],
            "top_pairs":          ccc.get("top_pairs", []) or [],
            "output_path":        ccc.get("output_path"),
        }
        # cellcomm mode's input h5ad is the annotated/clustered one; use
        # it as the figure source so UMAP fallback kicks in.
        cc_inputs = report.get("cellcomm_inputs", {}) or {}
        output_h5ad = output_h5ad or cc_inputs.get("h5ad")

    # Trajectory (v4.3.6) ─────────────────────────────────────────────────
    traj = stages.get("trajectory")
    if traj and traj.get("status") in ("success", "done"):
        findings["trajectory"] = {
            "status":     "done",
            "groupby":    traj.get("groupby"),
            "paga":       traj.get("paga", {}) or {},
            "pseudotime": traj.get("pseudotime", {}) or {},
            "velocity":   traj.get("velocity", {}) or {},
            "output_path": traj.get("output_path"),
        }
        # The trajectory script writes a new h5ad with PAGA + DPT
        # populated. Prefer it for downstream figure generation.
        output_h5ad = traj.get("output_path") or output_h5ad

    # Pseudobulk DE (v4.3.4+) ─────────────────────────────────────────────
    pb = stages.get("pseudobulk")
    if pb and pb.get("status") in ("success", "done"):
        findings["pseudobulk_de"] = {
            "groupby":       pb.get("groupby"),
            "condition_col": pb.get("condition_col"),
            "replicate_col": pb.get("replicate_col"),
            "covariates":    pb.get("covariates", []) or [],
            "thresholds":    pb.get("thresholds", {}) or {},
            "n_groups":      pb.get("n_groups"),
            "per_group":     pb.get("per_group", {}) or {},
        }
        # For pseudobulk runs, the source h5ad is the input — fall back
        # to that since there's no clustering output here.
        pb_inputs = report.get("pseudobulk_inputs", {}) or {}
        output_h5ad = output_h5ad or pb_inputs.get("h5ad")

    # Assemble the agent_results envelope ─────────────────────────────────
    scrna_envelope = {
        "status":      "done",
        "findings":    {"scRNA": {"findings": findings}},
        "output_h5ad": output_h5ad,
        "_workspace":  str(workspace),
        "_mode":       mode,
    }
    agent_results = {"scrna_agent": scrna_envelope}

    # Build exp_context + intent stubs so NarrativeAgent has something
    # to print in the meta strip and methods section.
    exp_context = {
        "organism":    organism,
        "genome":      _genome_for_organism(organism),
        "tissue":      tissue,
        "modalities":  {"scRNA": [str(f) for f in files_in]},
        "experiment_type": "scRNA_pseudobulk" if mode == "pseudobulk"
                          else "scRNA",
    }

    pb_inputs = report.get("pseudobulk_inputs", {}) or {}
    cc_inputs = report.get("cellcomm_inputs", {}) or {}
    if mode == "cellcomm":
        cb = cc_inputs.get("groupby", "?")
        user_question = (
            f"Which cell types signal to which in this dataset "
            f"(grouped by {cb}, autocrine excluded)?"
        )
        exp_context["design"] = {
            "groupby":  cb,
            "n_perms":  cc_inputs.get("n_perms"),
            "method":   findings.get("cell_communication", {}).get("method"),
        }
    elif mode == "pseudobulk":
        cond = pb_inputs.get("condition", "?")
        cmps = pb_inputs.get("comparisons", []) or []
        cmp_str = ", ".join(f"{a} vs {b}" for a, b in cmps) if cmps else "?"
        user_question = (
            f"What changes between {cond} levels ({cmp_str}) "
            f"in each cell type?"
        )
        exp_context["design"] = {
            "condition":      cond,
            "replicate":      pb_inputs.get("replicate"),
            "groupby":        pb_inputs.get("groupby"),
            "comparisons":    cmps,
            "covariates":     pb_inputs.get("covariates", []),
        }
    else:
        user_question = (
            f"Characterize cell types and their markers in this "
            f"{tissue or 'tissue'} scRNA-seq dataset."
        )
    exp_context["user_question"] = user_question

    intent = {
        "summary":              user_question,
        "biological_entities":  _entities_from_pb(pb_inputs),
        "experiment_type":      exp_context["experiment_type"],
    }

    return {
        "agent_results": agent_results,
        "exp_context":   exp_context,
        "intent":        intent,
        "decisions":     [],
        "findings_list": [],
    }


def _first_mt_threshold(per_sample: list) -> Optional[float]:
    """Pull the first sample's mt_threshold, if available."""
    for s in per_sample or []:
        v = s.get("mt_threshold") or s.get("mt_high")
        if v is not None:
            return v
    return None


def _genome_for_organism(organism: str) -> str:
    """Map a Latin name to the default genome build ARIA assumes."""
    o = (organism or "").lower()
    if "homo sapiens" in o or "human" in o:
        return "GRCh38"
    if "mus musculus" in o or "mouse" in o:
        return "GRCm39"
    return "unknown"


def _entities_from_pb(pb_inputs: dict) -> list:
    """Derive a short list of biological entities for slug/title."""
    ents = []
    cond = pb_inputs.get("condition")
    if cond:
        ents.append(str(cond))
    for a, b in pb_inputs.get("comparisons", []) or []:
        ents.extend([str(a), str(b)])
    return ents[:4]


if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("usage: rna_narrative_adapter.py <e2e_report.json>")
        sys.exit(1)
    bundle = adapt(sys.argv[1])
    print(json.dumps({k: (list(v.keys()) if isinstance(v, dict) else v)
                       for k, v in bundle.items()}, indent=2, default=str))
