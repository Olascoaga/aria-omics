#!/usr/bin/env python3
"""V48 Slice I7 — End-to-end real-data validation of Cross-Modal Evidence Graph,
IntegrationNarrator, Methods, Run-Ledger Reconciliation, and Dual PNG/SVG Figures.

Executes on paired 10x Multiome (`GSE278576` `hc11`):
  1. `aria/scripts/integration_paired_qc.py` (`aria-chromatin-env`) — I1 joint QC
     and bijective `/obs` cell-type label join.
  2. `aria/scripts/integration_wnn.py` (`aria-chromatin-env`) — I2 joint WNN +
     inline I3 5-axis anti-degradation guardrail.
  3. `aria/scripts/integration_mofa.py` (`aria-integration-env`) — I4 MOFA+
     decomposition with depth/cell-cycle confounder quarantine.
  4. `aria/scripts/integration_peak2gene.py` (`aria-integration-env`) — I5
     Signac-faithful GC + accessibility matched trans-chromosome cis linking.
  5. `aria/scripts/integration_eregulons.py` (`aria-integration-env`) — I6
     SCENIC+-style candidate eRegulon assembly (`pos_tf_pos_p2g` / `neg_tf_pos_p2g`).
  6. Two independent executions of `NarrativeAgent.run()` (`reproducible_mode=True`)
     verifying:
       - Bit-exact `cross_modal_evidence_graph` determinism,
       - Zero C1-C4 / W-CLAIM / W-LEDGER violations across all 6 `multiome.*` blocks,
       - Claim 6 monotone LLM-variation invariance (paraphrase invariance +
         adversarial causal injection withholding),
       - Dual `.png` + vector `.svg` (`<text>` preserved) publication figures,
       - Staged CSV tables (`mofa_factor_loadings.csv`, `peak_gene_links.csv`,
         `eregulons.csv`, `eregulon_triplets.csv`),
       - HTML-safe copy-pasteable Methods section and RO-Crate provenance bundle.

Emits a sanitized public JSON artifact with zero machine-absolute paths
(enforced by `tests/test_public_artifacts_no_absolute_paths.py`).
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from aria.agents.narrative.claim_compiler import compile_public_claims  # noqa: E402
from aria.agents.narrative.compose_prose import compose_block_prose  # noqa: E402
from aria.agents.narrative.narrators.integration import (  # noqa: E402
    IntegrationNarrator,
)
from aria.agents.narrative.run_ledger import (  # noqa: E402
    build_run_ledger,
    ensure_report_ledger_nodes,
)
from aria.agents.narrative.validators import find_causal_language  # noqa: E402
from aria.agents.narrative_agent import NarrativeAgent  # noqa: E402
from aria.memory.memory import ARIAMemory  # noqa: E402
from aria.version import __version__, collect_version_metadata  # noqa: E402

_FORBIDDEN_PATH_PREFIXES = (
    "/home/",
    "/Users/",
    "/root/",
    "/tmp/",
    "/mnt/",
    "/private/",
    "/data/",
    "/var/",
    "/opt/",
)


def _resolve_conda_bin() -> str:
    found = shutil.which("conda")
    if found:
        return found
    cand = Path(sys.executable).resolve().parents[2] / "bin" / "conda"
    if cand.is_file():
        return str(cand)
    return "conda"


def _run_ipc_script(
    env: str,
    script_rel: str,
    params: dict[str, Any],
    work_dir: Path,
    tag: str,
) -> tuple[dict[str, Any], float]:
    in_json = work_dir / f"{tag}_in.json"
    out_json = work_dir / f"{tag}_out.json"
    in_json.write_text(json.dumps(params, indent=2), encoding="utf-8")
    t0 = time.time()
    proc = subprocess.run(
        [
            _resolve_conda_bin(),
            "run",
            "--no-capture-output",
            "-n",
            env,
            "python",
            str(ROOT / script_rel),
            str(in_json),
            str(out_json),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    elapsed = round(time.time() - t0, 2)
    if proc.returncode != 0 or not out_json.is_file():
        raise RuntimeError(
            f"{script_rel} in {env} exited with {proc.returncode}:\n"
            f"STDOUT:\n{proc.stdout[-2500:]}\nSTDERR:\n{proc.stderr[-2500:]}"
        )
    return json.loads(out_json.read_text(encoding="utf-8")), elapsed


def _sanitize_public_obj(val: Any) -> Any:
    """Recursively strip machine-absolute filesystem paths from public manifests."""
    if isinstance(val, dict):
        return {str(k): _sanitize_public_obj(v) for k, v in val.items()}
    if isinstance(val, list):
        return [_sanitize_public_obj(v) for v in val]
    if isinstance(val, str):
        if val.startswith(_FORBIDDEN_PATH_PREFIXES):
            p = Path(val)
            try:
                return str(p.resolve().relative_to(ROOT.resolve()))
            except Exception:
                return p.name
        return val
    return val


def _assert_no_forbidden_paths(obj: Any, trail: str = "root") -> None:
    if isinstance(obj, dict):
        for k, v in obj.items():
            _assert_no_forbidden_paths(v, f"{trail}.{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            _assert_no_forbidden_paths(v, f"{trail}[{i}]")
    elif isinstance(obj, str):
        if obj.startswith(_FORBIDDEN_PATH_PREFIXES):
            raise AssertionError(f"Forbidden absolute path at {trail}: {obj}")


def _run_narrative_once(
    integration_envelope: dict[str, Any],
    reports_root: Path,
    experiment_id: str,
) -> tuple[dict[str, Any], Path, dict[str, Any], str]:
    mem = ARIAMemory(db_path=":memory:")
    agent = NarrativeAgent(memory=mem, llm=None)
    agent.reports_dir = reports_root
    agent.reports_dir.mkdir(parents=True, exist_ok=True)

    env_copy = copy.deepcopy(integration_envelope)
    t0 = time.time()
    res = agent.run(
        experiment_id=experiment_id,
        context={
            "exp_context": {
                "organism": "Homo sapiens",
                "genome": "GRCh38",
                "user_question": (
                    "Characterize paired 10x Multiome joint manifold structure, "
                    "latent variance programs, cis peak-to-gene associations, "
                    "and candidate eRegulons in human hippocampus."
                ),
                "reproducible_mode": True,
            },
            "biological_intent": {
                "summary": "Paired human hippocampus 10x Multiome integrative regulatory audit"
            },
            "agent_results": {"integration_agent": env_copy},
            "findings": [],
        },
    )
    res["_narrative_runtime_sec"] = round(time.time() - t0, 2)
    if res.get("status") != "done":
        raise RuntimeError(f"NarrativeAgent failed: {res}")
    report_html = Path(str(res["report_path"]))
    report_dir = report_html.parent
    meth_path = report_dir / "methodology.json"
    rocrate_path = report_dir / "ro-crate-metadata.json"
    if not (report_html.is_file() and meth_path.is_file() and rocrate_path.is_file()):
        raise RuntimeError(f"Missing report bundle files in {report_dir}")

    meth = json.loads(meth_path.read_text(encoding="utf-8"))
    html_text = report_html.read_text(encoding="utf-8")
    return res, report_dir, meth, html_text


def _verify_claim6_monotone_invariance(
    integration_envelope: dict[str, Any],
) -> dict[str, Any]:
    """Verify paraphrase invariance + causal injection withholding on real hc11 blocks."""
    narrator = IntegrationNarrator()
    agent_results = {"integration_agent": copy.deepcopy(integration_envelope)}
    ctx = {
        "exp_context": {"organism": "Homo sapiens", "genome": "GRCh38"},
        "agent_results": agent_results,
    }
    base_blocks = narrator.collect(
        "integration_agent", copy.deepcopy(integration_envelope), ctx
    )
    ledger = build_run_ledger({}, agent_results)
    ensure_report_ledger_nodes(ledger, base_blocks)
    base_comp = compile_public_claims(base_blocks, exp_ctx={}, run_ledger=ledger)

    # 1. Benign associative claim paraphrase does not alter deterministic rendered_prose.
    ereg_orig = next(b for b in base_blocks if b.id == "multiome.eregulons")
    prose_orig = compose_block_prose(ereg_orig)
    ereg_para = copy.deepcopy(ereg_orig)
    m = ereg_para.metrics
    ereg_para.claim = (
        f"Across {m['n_shared_cells']} shared cells, {m['n_candidate_modules_tested']} "
        f"modules were tested and {m['n_eregulons']} BH-FDR candidate eRegulons "
        f"({m['n_pos_tf_pos_p2g']} pos_tf_pos_p2g, {m['n_neg_tf_pos_p2g']} neg_tf_pos_p2g) "
        f"spanning {m['n_unique_tfs']} unique TFs and {m['n_triplets']} triplets "
        f"were retained ({m['n_paralog_ambiguous']} paralog-ambiguous, "
        f"{m['n_depth_confounded']} depth-confounded)."
    )
    para_comp = compile_public_claims([ereg_para], exp_ctx={}, run_ledger=ledger)
    paraphrase_invariant = (
        compose_block_prose(ereg_para) == prose_orig
        and len(para_comp.claims) == 1
        and len(para_comp.withheld) == 0
    )

    # 2. Adversarial causal verb injection on multiome.eregulons must be withheld.
    ereg_adv = copy.deepcopy(ereg_orig)
    top_tf = (ereg_adv.metrics.get("top_tfs") or ["TF"])[0]
    ereg_adv.claim = (
        f"{top_tf} directly regulates downstream transcription across "
        f"{m['n_shared_cells']} shared cells."
    )
    adv_comp = compile_public_claims([ereg_adv], exp_ctx={}, run_ledger=ledger)
    causal_injection_withheld = (
        len(adv_comp.claims) == 0 and len(adv_comp.withheld) == 1
    )

    return {
        "n_baseline_published": len(base_comp.claims),
        "n_baseline_withheld": len(base_comp.withheld),
        "paraphrase_invariant": bool(paraphrase_invariant),
        "causal_injection_target": "multiome.eregulons",
        "causal_injection_withheld": bool(causal_injection_withheld),
        "n_published_after_causal_injection": len(adv_comp.claims),
        "n_withheld_after_causal_injection": len(adv_comp.withheld),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5mu", required=True, help="Input paired 10x Multiome .h5mu.")
    ap.add_argument(
        "--reference-h5ad",
        default="",
        help="Optional reference .h5ad for bijective /obs label join in I1.",
    )
    ap.add_argument("--sample-id", default="hc11")
    ap.add_argument("--audit-group-col", default="subclass")
    ap.add_argument("--reference-sample-col", default="orig.ident")
    ap.add_argument("--gtf", required=True, help="Reference GTF annotation.")
    ap.add_argument("--genome-fasta", required=True, help="Indexed reference FASTA.")
    ap.add_argument("--motif-collection", default="JASPAR2024_CORE_vertebrates")
    ap.add_argument("--chromatin-env", default="aria-chromatin-env")
    ap.add_argument("--integration-env", default="aria-integration-env")
    ap.add_argument("--p2g-min-corr", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--output-json",
        default="docs/benchmark_results/multiome_concordance/v48_i7_cross_modal_narrative_hc11.json",
    )
    args = ap.parse_args(argv)

    with tempfile.TemporaryDirectory(prefix="v48_i7_narrative_") as td:
        work = Path(td)
        timings: dict[str, float] = {}

        print("[1/6] Running I1 integration_paired_qc...", flush=True)
        qc_params: dict[str, Any] = {
            "h5mu_path": str(args.h5mu),
            "output_dir": str(work / "i1_qc"),
            "sample_id": str(args.sample_id),
            "organism": "Homo sapiens",
            "genome": "hg38",
        }
        if args.reference_h5ad and Path(args.reference_h5ad).is_file():
            qc_params.update(
                {
                    "reference_h5ad": str(args.reference_h5ad),
                    "reference_sample_col": str(args.reference_sample_col),
                    "reference_label_col": str(args.audit_group_col),
                    "target_label_col": str(args.audit_group_col),
                }
            )
        qc_res, timings["i1_paired_qc_sec"] = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_paired_qc.py",
            qc_params,
            work,
            "i1_qc",
        )
        if qc_res.get("status") != "success":
            raise RuntimeError(f"I1 paired_qc failed: {qc_res}")
        post_qc_h5mu = Path(str(qc_res["output_path"]))
        print(
            f"      Retained {qc_res['n_cells_after']}/{qc_res['n_cells_before']} "
            f"shared cells ({timings['i1_paired_qc_sec']}s)",
            flush=True,
        )

        print("[2/6] Running I2/I3 integration_wnn + 5-axis guardrail...", flush=True)
        wnn_params: dict[str, Any] = {
            "h5mu_path": str(post_qc_h5mu),
            "output_dir": str(work / "i2_wnn"),
            "k": 20,
            "resolution": 0.5,
            "seed": int(args.seed),
        }
        if args.audit_group_col:
            wnn_params["audit_group_col"] = str(args.audit_group_col)
        wnn_res, timings["i2_i3_wnn_sec"] = _run_ipc_script(
            args.chromatin_env,
            "aria/scripts/integration_wnn.py",
            wnn_params,
            work,
            "i2_wnn",
        )
        if wnn_res.get("status") != "success":
            raise RuntimeError(f"I2/I3 WNN failed: {wnn_res}")
        print(
            f"      WNN clusters={wnn_res['n_joint_clusters']}, "
            f"guardrail={wnn_res['guardrail_status']} ({timings['i2_i3_wnn_sec']}s)",
            flush=True,
        )

        print("[3/6] Running I4 integration_mofa...", flush=True)
        mofa_res, timings["i4_mofa_sec"] = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_mofa.py",
            {
                "h5mu_path": str(post_qc_h5mu),
                "output_dir": str(work / "i4_mofa"),
                "n_factors": 10,
                "n_top_genes": 3000,
                "n_top_peaks": 5000,
                "seed": int(args.seed),
                "write_h5mu": False,
            },
            work,
            "i4_mofa",
        )
        if mofa_res.get("status") != "success":
            raise RuntimeError(f"I4 MOFA+ failed: {mofa_res}")
        print(
            f"      MOFA+ active={mofa_res['n_active_factors']}/{mofa_res['n_factors']}, "
            f"excluded={len(mofa_res.get('excluded_factors') or [])} ({timings['i4_mofa_sec']}s)",
            flush=True,
        )

        print("[4/6] Running I5 integration_peak2gene...", flush=True)
        p2g_res, timings["i5_peak2gene_sec"] = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_peak2gene.py",
            {
                "h5mu_path": str(post_qc_h5mu),
                "output_dir": str(work / "i5_p2g"),
                "gtf_path": str(args.gtf),
                "genome_fasta": str(args.genome_fasta),
                "distance_kb": 500,
                "min_corr": float(args.p2g_min_corr),
                "max_padj": 0.05,
                "seed": int(args.seed),
            },
            work,
            "i5_p2g",
        )
        if p2g_res.get("status") != "success":
            raise RuntimeError(f"I5 peak2gene failed: {p2g_res}")
        print(
            f"      Peak2Gene links={p2g_res['n_links']} "
            f"(tested={p2g_res['n_tested']}, {timings['i5_peak2gene_sec']}s)",
            flush=True,
        )

        print("[5/6] Running I6 integration_eregulons...", flush=True)
        ereg_res, timings["i6_eregulons_sec"] = _run_ipc_script(
            args.integration_env,
            "aria/scripts/integration_eregulons.py",
            {
                "h5mu_path": str(post_qc_h5mu),
                "output_dir": str(work / "i6_ereg"),
                "gtf_path": str(args.gtf),
                "genome_fasta": str(args.genome_fasta),
                "motif_collection": str(args.motif_collection),
                "distance_kb": 500,
                "motif_pvalue": 5e-5,
                "p2g_min_corr": 0.05,
                "p2g_max_padj": 0.05,
                "min_tf_gene_corr": 0.05,
                "min_target_genes": 5,
                "max_regulon_padj": 0.05,
                "standard_chroms_only": True,
                "seed": int(args.seed),
            },
            work,
            "i6_ereg",
        )
        if ereg_res.get("status") != "success":
            raise RuntimeError(f"I6 eRegulons failed: {ereg_res}")
        print(
            f"      eRegulons={ereg_res['n_eregulons']}, "
            f"triplets={ereg_res['n_triplets']} ({timings['i6_eregulons_sec']}s)",
            flush=True,
        )

        integration_envelope = {
            "status": "done",
            "strategy": "wnn_then_mofa",
            "findings": {
                "multiome_gate": {
                    "allowed": True,
                    "passed": True,
                    "object_type": "paired_mudata",
                    "same_cell": True,
                },
                "paired_qc": qc_res,
                "wnn": wnn_res,
                "mofa": mofa_res,
                "peak_to_gene": p2g_res,
                "eregulons": ereg_res,
            },
        }

        print("[6/6] Running I7 NarrativeAgent E2E x2 + governance & figure audits...", flush=True)
        res_1, report_dir_1, meth_1, html_1 = _run_narrative_once(
            integration_envelope, work / "reports_1", "v48_i7_hc11_run1"
        )
        res_2, _, meth_2, _ = _run_narrative_once(
            integration_envelope, work / "reports_2", "v48_i7_hc11_run2"
        )
        timings["i7_narrative_run1_sec"] = res_1["_narrative_runtime_sec"]
        timings["i7_narrative_run2_sec"] = res_2["_narrative_runtime_sec"]

        graph_1 = meth_1["cross_modal_evidence_graph"]
        graph_2 = meth_2["cross_modal_evidence_graph"]
        graph_bit_exact = (
            json.dumps(graph_1, sort_keys=True) == json.dumps(graph_2, sort_keys=True)
        )
        if not graph_bit_exact:
            raise AssertionError("cross_modal_evidence_graph differed across repeat runs")
        if not (
            graph_1["ran"] is True
            and graph_1["non_causal_contract"] is True
            and graph_1["circular_convergence_prevented"] is True
        ):
            raise AssertionError(f"Evidence graph contract flags failed: {graph_1}")

        # Verify executive_summary + all 6 multiome.* blocks published with 0 withheld & 0 ledger violations.
        nblocks = meth_1.get("narrative_blocks") or []
        claims = meth_1.get("claims") or []
        block_ids = [b.get("id") for b in nblocks]
        multiome_block_ids = [bid for bid in block_ids if str(bid).startswith("multiome.")]
        expected_ids = [
            "multiome.paired_qc",
            "multiome.wnn",
            "multiome.mofa",
            "multiome.peak_to_gene",
            "multiome.eregulons",
            "multiome.triangulation",
        ]
        if multiome_block_ids != expected_ids:
            raise AssertionError(f"Expected {expected_ids}, got {multiome_block_ids} (all={block_ids})")
        if meth_1["claim_compilation"]["n_withheld"] != 0:
            raise AssertionError(f"Unexpected withheld claims: {meth_1['claim_compilation']}")
        if meth_1["run_ledger"]["claim_ledger_verification"]["n_violations"] != 0:
            raise AssertionError(
                f"Ledger verification violations: {meth_1['run_ledger']['claim_ledger_verification']}"
            )
        for c in claims:
            if c.get("tier") not in {"descriptive", "associative"}:
                raise AssertionError(f"Non-associative tier on {c}: {c.get('tier')}")
            if c.get("language_violation") is not None:
                raise AssertionError(f"Language violation on {c}")
            lnode = str(c.get("ledger_node_id") or "")
            bid = str(c.get("block_id") or "")
            if bid.startswith("multiome.") and not lnode.startswith("ledger://integration/"):
                raise AssertionError(f"Unexpected integration ledger_node_id on {c}")
            if bid == "executive_summary" and lnode != "ledger://report/executive_summary":
                raise AssertionError(f"Unexpected exec summary ledger_node_id on {c}")
        for nb in nblocks:
            if find_causal_language(str(nb.get("claim") or "")) is not None:
                raise AssertionError(f"Causal claim in {nb.get('id')}")
            if find_causal_language(str(nb.get("rendered_prose") or "")) is not None:
                raise AssertionError(f"Causal prose in {nb.get('id')}")

        # Verify all 4 dual PNG + SVG figures exist, SVGs preserve editable <text>, and HTML embeds base64 PNGs.
        expected_figures = (
            "wnn_joint_vs_modality_umap",
            "wnn_modality_weights_by_cluster",
            "mofa_variance_explained",
            "regulatory_triangulation_summary",
        )
        n_embedded_base64_pngs = html_1.count("data:image/png;base64,")
        if n_embedded_base64_pngs < 4:
            raise AssertionError(
                f"Expected >=4 self-contained base64 PNG embeds in report.html, got {n_embedded_base64_pngs}"
            )
        figure_audit: dict[str, dict[str, Any]] = {}
        for fkey in expected_figures:
            png_f = report_dir_1 / "figures" / f"{fkey}.png"
            svg_f = report_dir_1 / "figures" / f"{fkey}.svg"
            if not (png_f.is_file() and svg_f.is_file()):
                raise AssertionError(f"Missing dual figure {fkey} in {report_dir_1 / 'figures'}")
            png_bytes = png_f.stat().st_size
            svg_bytes = svg_f.stat().st_size
            svg_txt = svg_f.read_text(encoding="utf-8")
            has_text = "<svg" in svg_txt and "<text" in svg_txt
            if png_bytes <= 1000 or svg_bytes <= 1000 or not has_text:
                raise AssertionError(f"Invalid dual figure {fkey}: png={png_bytes}, svg={svg_bytes}")
            figure_audit[fkey] = {
                "png_rel": f"figures/{fkey}.png",
                "svg_rel": f"figures/{fkey}.svg",
                "png_bytes": int(png_bytes),
                "svg_bytes": int(svg_bytes),
                "svg_has_editable_text": bool(has_text),
                "embedded_base64_in_report_html": True,
            }

        # Verify all 4 supplementary CSV tables are staged in report_dir/tables/.
        expected_tables = (
            "mofa_factor_loadings.csv",
            "peak_gene_links.csv",
            "eregulons.csv",
            "eregulon_triplets.csv",
        )
        table_audit: dict[str, dict[str, Any]] = {}
        for tname in expected_tables:
            tp = report_dir_1 / "tables" / tname
            if not tp.is_file():
                raise AssertionError(f"Missing staged CSV table: {tp}")
            lines = [ln for ln in tp.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if len(lines) < 2:
                raise AssertionError(f"Staged table {tname} has <2 rows")
            table_audit[tname] = {
                "rel_path": f"tables/{tname}",
                "n_rows_including_header": len(lines),
                "bytes": int(tp.stat().st_size),
            }

        # Verify HTML report safety & Methods completeness.
        if "&amp;lt;" in html_1 or "<script>alert" in html_1:
            raise AssertionError("HTML escaping or tag sanitization check failed in report.html")
        methods_paragraphs = IntegrationNarrator().methods(
            "integration_agent", integration_envelope
        )
        if len(methods_paragraphs) != 5:
            raise AssertionError(f"Expected 5 Methods paragraphs, got {len(methods_paragraphs)}")

        claim6_audit = _verify_claim6_monotone_invariance(integration_envelope)
        if not (
            claim6_audit["paraphrase_invariant"]
            and claim6_audit["causal_injection_withheld"]
        ):
            raise AssertionError(f"Claim 6 invariance audit failed: {claim6_audit}")

        manifest = _sanitize_public_obj(
            {
                "benchmark": "v48_i7_cross_modal_narrative_governance",
                "slice": "V48_I7_cross_modal_evidence_graph_and_narrator",
                "aria_version": __version__,
                "provenance": collect_version_metadata(),
                "dataset": {
                    "accession": "GSE278576",
                    "sample_id": str(args.sample_id),
                    "tissue": "human_hippocampus_multiome_10x",
                    "n_cells_raw": int(qc_res["n_cells_before"]),
                    "n_cells_post_qc": int(qc_res["n_cells_after"]),
                    "n_rna_genes": int(qc_res["n_genes"]),
                    "n_atac_peaks": int(qc_res["n_peaks"]),
                },
                "stage_timings_sec": timings,
                "upstream_stage_summary": {
                    "paired_qc": {
                        "pct_removed": qc_res["pct_removed"],
                        "metrics_not_computed": qc_res.get("metrics_not_computed", []),
                    },
                    "wnn": {
                        "n_joint_clusters": wnn_res["n_joint_clusters"],
                        "n_rna_only_clusters": wnn_res["n_rna_only_clusters"],
                        "n_atac_only_clusters": wnn_res["n_atac_only_clusters"],
                        "mean_rna_weight": wnn_res["mean_rna_weight"],
                        "mean_atac_weight": wnn_res["mean_atac_weight"],
                        "wnn_vs_rna_ari": wnn_res["wnn_vs_rna_ari"],
                        "wnn_vs_atac_ari": wnn_res["wnn_vs_atac_ari"],
                        "n_discordant_clusters": wnn_res["n_discordant_clusters"],
                        "guardrail_status": wnn_res["guardrail_status"],
                    },
                    "mofa": {
                        "n_factors": mofa_res["n_factors"],
                        "n_active_factors": mofa_res["n_active_factors"],
                        "n_shared_factors": len(mofa_res.get("shared_factors") or []),
                        "n_excluded_factors": len(mofa_res.get("excluded_factors") or []),
                        "variance_explained": mofa_res["variance_explained"],
                        "variance_explained_additive": mofa_res["variance_explained_additive"],
                    },
                    "peak_to_gene": {
                        "distance_kb": p2g_res["distance_kb"],
                        "min_corr": p2g_res["min_corr"],
                        "max_padj": p2g_res["max_padj"],
                        "n_tested": p2g_res["n_tested"],
                        "n_links": p2g_res["n_links"],
                        "n_positive_correlations": p2g_res["n_positive_correlations"],
                        "n_negative_correlations": p2g_res["n_negative_correlations"],
                        "n_promoter_proximal": p2g_res["n_promoter_proximal"],
                        "n_distal_cis": p2g_res["n_distal_cis"],
                    },
                    "eregulons": {
                        "motif_pvalue": ereg_res["motif_pvalue"],
                        "n_candidate_modules_tested": ereg_res["n_candidate_modules_tested"],
                        "n_eregulons": ereg_res["n_eregulons"],
                        "n_triplets": ereg_res["n_triplets"],
                        "n_unique_tfs": ereg_res["n_unique_tfs"],
                        "n_pos_tf_pos_p2g": ereg_res["n_pos_tf_pos_p2g"],
                        "n_neg_tf_pos_p2g": ereg_res["n_neg_tf_pos_p2g"],
                        "n_paralog_ambiguous": ereg_res["n_paralog_ambiguous"],
                        "n_depth_confounded": ereg_res["n_depth_confounded"],
                    },
                },
                "narrative_governance_audit": {
                    "published_block_ids": block_ids,
                    "n_published_blocks": len(claims),
                    "n_withheld_blocks": int(meth_1["claim_compilation"]["n_withheld"]),
                    "ledger_n_violations": int(
                        meth_1["run_ledger"]["claim_ledger_verification"]["n_violations"]
                    ),
                    "all_estimands_associative": True,
                    "zero_causal_language_matches": True,
                    "methods_paragraphs_count": len(methods_paragraphs),
                    "html_double_escape_free": True,
                    "ro_crate_emitted": True,
                    "claim6_monotone_invariance": claim6_audit,
                },
                "cross_modal_evidence_graph_audit": {
                    "schema_version": graph_1["schema_version"],
                    "ran": graph_1["ran"],
                    "repeat_run_bit_exact": graph_bit_exact,
                    "non_causal_contract": graph_1["non_causal_contract"],
                    "circular_convergence_prevented": graph_1["circular_convergence_prevented"],
                    "executed_integration_stages": graph_1["executed_integration_stages"],
                    "absent_not_tested": graph_1["absent_not_tested"],
                    "wnn_guardrail": graph_1["wnn_guardrail"],
                    "mofa_summary": graph_1["mofa_summary"],
                    "summary_counts": graph_1["summary_counts"],
                    "top_regulatory_axes": graph_1["regulatory_axes"][:15],
                },
                "dual_png_svg_figures": figure_audit,
                "staged_csv_tables": table_audit,
            }
        )

        _assert_no_forbidden_paths(manifest)
        out_path = (
            Path(args.output_json)
            if Path(args.output_json).is_absolute()
            else (ROOT / args.output_json)
        )
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        print(
            f"Wrote sanitized V48 I7 validation manifest: "
            f"{out_path.relative_to(ROOT) if out_path.is_relative_to(ROOT) else out_path}",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
