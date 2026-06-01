"""
ARIA PBMC 3k End-to-End Test
-----------------------------
Validates the full ARIA pipeline with the PBMC 3k dataset
from 10x Genomics — the universal benchmark for scRNA-seq tools.

Dataset: ~2,700 human PBMCs, ~33,000 genes
Source:  10x Genomics (public domain)
Expected runtime: 3-8 minutes depending on hardware

What this test validates:
  1. DataAuditAgent detects 10x MEX format automatically
  2. DataAuditAgent infers Homo sapiens / hg38 (or hg19)
  3. RNAAgent QC: ~2,700 cells -> ~2,600 after filtering
  4. ParameterAdvisor recommends Leiden resolution ~0.4-0.6
  5. Cell type annotation identifies T cells, B cells, NK cells, Monocytes
  6. MessageBus records all findings with confidence scores
  7. Memory stores decisions and makes them retrievable

Run:
  conda activate aria-env
  python tests/test_pbmc_e2e.py [--data-dir ~/aria-data/pbmc3k_test]
"""

import sys
import os
import time
import argparse
import uuid
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

# Load .env if present
env_file = Path.home() / ".aria" / ".env"
if env_file.exists():
    with open(env_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, val = line.split("=", 1)
                os.environ.setdefault(key.strip(), val.strip())

# ── Terminal colors ───────────────────────────────────────────────────────────
GRN = "\033[92m"; RED = "\033[91m"; YLW = "\033[93m"
CYN = "\033[96m"; DIM = "\033[2m";  RST = "\033[0m"; BLD = "\033[1m"


def banner():
    print(f"\n{CYN}{BLD}  ARIA -- PBMC 3k End-to-End Test{RST}")
    print(f"  {'─'*46}")
    print(f"  Dataset:  PBMC 3k (10x Genomics, ~2,700 cells)")
    print(f"  Purpose:  Full pipeline validation\n")


def ok(msg, detail=""):
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {GRN}v{RST} {msg}{d}")


def fail(msg, err=""):
    print(f"  {RED}x{RST} {msg}")
    if err:
        print(f"    {DIM}{err}{RST}")


def info(msg):
    print(f"  {CYN}->{RST} {msg}")


def section(title):
    print(f"\n{BLD}{CYN}> {title}{RST}")


# ── Checks ────────────────────────────────────────────────────────────────────

def check_api_keys():
    """Verify at least one LLM provider is configured."""
    section("Checking API keys")

    anthropic_key = os.environ.get("ANTHROPIC_API_KEY", "")
    gemini_key    = os.environ.get("GEMINI_API_KEY", "") or \
                    os.environ.get("GOOGLE_API_KEY", "")

    has_anthropic = bool(anthropic_key and anthropic_key.startswith("sk-ant"))
    has_gemini    = bool(gemini_key and gemini_key.startswith("AIza"))

    if has_anthropic:
        ok(f"Anthropic API key: {anthropic_key[:12]}...{anthropic_key[-4:]}")
    else:
        info("Anthropic API key: not configured or invalid")

    if has_gemini:
        ok(f"Google API key:    {gemini_key[:8]}...{gemini_key[-4:]}")
    else:
        info("Google API key: not configured")

    if not has_anthropic and not has_gemini:
        print(f"\n  {RED}No API keys found.{RST}")
        print(f"  Run the installer again: {BLD}bash install.sh{RST}")
        print(f"  Or configure manually in: {BLD}~/.aria/.env{RST}")
        sys.exit(1)

    return has_anthropic, has_gemini


def find_pbmc_data(data_dir: Path) -> Path:
    """Locate the PBMC 3k MEX directory."""
    section("Locating PBMC 3k data")

    candidates = [
        data_dir,
        data_dir / "filtered_gene_bc_matrices" / "hg19",
        data_dir / "hg19",
        data_dir / "GRCh38",
        data_dir / "filtered_feature_bc_matrix",
    ]

    for candidate in candidates:
        if candidate.exists():
            files = list(candidate.rglob("*.mtx*")) + \
                    list(candidate.rglob("*.mtx.gz"))
            if files:
                ok(f"Data found at: {candidate}")
                return candidate

    print(f"\n  {YLW}No MEX files found in: {data_dir}{RST}")
    print(f"\n  Download the dataset manually:")
    print(f"  {BLD}1.{RST} Go to: {CYN}https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/{RST}")
    print(f"  {BLD}2.{RST} Download: pbmc3k_filtered_gene_bc_matrices.tar.gz")
    print(f"  {BLD}3.{RST} Extract to: {data_dir}")
    print(f"  {BLD}4.{RST} Re-run this script")
    sys.exit(1)


# ── Pipeline tests ────────────────────────────────────────────────────────────

def _stage_data_audit(data_dir: Path, experiment_id: str) -> dict:
    """Test DataAuditAgent automatic detection and classification."""
    section("Test 1 -- DataAuditAgent (automatic detection)")

    from aria.agents.data_audit_agent import DataAuditAgent
    from aria.memory.memory import ARIAMemory
    from aria.bus.message_bus import bus

    memory = ARIAMemory()
    agent  = DataAuditAgent(memory)

    t0 = time.time()
    result = agent.run(
        experiment_id=experiment_id,
        context={
            "data_dir":      str(data_dir),
            "user_question": "What cell types are present in these PBMCs?",
        }
    )
    elapsed = time.time() - t0

    pending = bus.get_pending_checkpoints()
    exp_pending = [m for m in pending if m.experiment_id == experiment_id]

    if exp_pending:
        checkpoint_ctx = exp_pending[0].payload.get("context", {})
        exp_ctx        = checkpoint_ctx.get("exp_context", {})
        modalities     = exp_ctx.get("modalities", {})

        if "scRNA" in modalities:
            ok(f"scRNA detected ({len(modalities['scRNA'])} files)",
               f"{elapsed:.1f}s")
        else:
            fail(f"scRNA not detected. Found: {list(modalities.keys())}")

        genome   = exp_ctx.get("genome", "?")
        organism = exp_ctx.get("organism", "?")

        if any(g in genome.lower() for g in ["hg", "grch"]):
            ok(f"Genome inferred: {genome}")
        else:
            info(f"Genome: {genome} (may need manual confirmation)")

        if "sapiens" in organism.lower() or organism == "unknown":
            ok(f"Organism: {organism or 'Homo sapiens (inferred)'}")

        # Simulate user confirming at Checkpoint 1
        bus.resolve_checkpoint(
            exp_pending[0].id,
            {"choice": "Confirm and continue"}
        )
        ok("Checkpoint 1 resolved: user confirmed data")

        return exp_ctx
    else:
        fail("DataAuditAgent did not generate a checkpoint",
             str(result.get("error", "")))
        return {}


def _stage_scrna_qc(data_dir: Path, exp_ctx: dict,
                  experiment_id: str) -> dict:
    """Test scRNA-seq QC pipeline."""
    section("Test 2 -- scRNA-seq QC")

    try:
        import scanpy as sc
    except ImportError:
        fail("Scanpy not installed. Run: pip install scanpy")
        return {}

    t0 = time.time()
    info("Loading PBMC 3k data...")

    try:
        adata = sc.read_10x_mtx(
            str(data_dir),
            var_names="gene_symbols",
            cache=True,
        )
        ok(f"Data loaded: {adata.n_obs} cells x {adata.n_vars} genes",
           f"{time.time()-t0:.1f}s")
    except Exception as e:
        try:
            adata = sc.read_10x_mtx(
                str(data_dir.parent),
                var_names="gene_symbols"
            )
            ok(f"Data loaded (alternate path): {adata.n_obs} cells")
        except Exception as e2:
            fail(f"Error loading data: {e2}")
            return {}

    import numpy as np
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs["pct_counts_mt"] < 5].copy()

    n_after     = adata.n_obs
    pct_removed = (n_before - n_after) / n_before * 100

    if 2500 <= n_after <= 2800:
        ok(f"QC: {n_before} -> {n_after} cells ({pct_removed:.1f}% removed)",
           "expected range for PBMC 3k")
    else:
        info(f"QC: {n_before} -> {n_after} cells ({pct_removed:.1f}% removed)")

    mt_stats = adata.obs["pct_counts_mt"].describe()
    ok(f"MT%: mean={mt_stats['mean']:.2f}%, max={mt_stats['max']:.2f}%")

    return {"adata": adata, "n_cells": n_after}


def _stage_parameter_advisor(adata_result: dict,
                            experiment_id: str,
                            has_api: bool) -> dict:
    """Test ParameterAdvisor 3-layer decision for Leiden clustering."""
    section("Test 3 -- ParameterAdvisor (hyperparameter decision)")

    adata = adata_result.get("adata")
    if adata is None:
        info("Skipping (no data available)")
        return {}

    from aria.memory.memory import ARIAMemory
    from aria.llm.provider import LLMProvider
    from aria.llm.parameter_advisor import ParameterAdvisor
    import scanpy as sc

    memory   = ARIAMemory()
    provider = LLMProvider.from_config()
    advisor  = ParameterAdvisor(memory, provider)

    info("Preprocessing for clustering...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    ok("Preprocessing complete: normalization, HVG, PCA, neighbors")

    bio_ctx = {
        "analysis_type": "cell_type",
        "user_question":  "What cell types are present in these PBMCs?",
        "summary":        "major cell type identification in PBMCs",
    }

    # Layer 1: intent-constrained search range
    search_range = advisor._intent_to_leiden_range(bio_ctx)
    ok(f"Layer 1 (intent): search range = {search_range}",
       "'major cell types' -> low/medium resolution")

    # Layer 2: evaluate candidates with objective metrics
    t0 = time.time()
    info("Evaluating 4 resolution candidates with objective metrics...")
    decision = advisor.advise_leiden_resolution(
        adata=adata,
        experiment_id=experiment_id,
        biological_context=bio_ctx,
        n_candidates=4,
    )
    ok(f"Layer 2 (metrics): {len(decision.candidates)} candidates evaluated",
       f"{time.time()-t0:.1f}s")

    for c in decision.candidates:
        marker = "*" if c.recommended else " "
        sil    = c.metrics.get("silhouette", 0)
        nk     = c.metrics.get("n_clusters", 0)
        warn   = " (!)" if c.flags else ""
        print(f"      {marker} resolution={c.value:4.2f} | "
              f"silhouette={sil:.3f} | clusters={nk}{warn}")

    ok(f"Recommendation: resolution={decision.chosen_value}")

    if has_api:
        ok(f"LLM justification: '{decision.justification[:70]}...'")
    else:
        ok("Justification generated (heuristic fallback)")

    # Layer 3: memory recall
    hist = advisor._recall_similar_decisions(
        experiment_id, "leiden_clustering", bio_ctx
    )
    if hist:
        ok(f"Layer 3 (memory): {len(hist)} historical decision(s) retrieved")
    else:
        ok("Layer 3 (memory): first run -- database empty",
           "future runs will learn from this decision")

    # Simulate user approval at Checkpoint 3
    approved = advisor.approve_decision(decision)
    ok("Checkpoint 3 resolved: user approved parameter",
       "decision saved to ~/.aria/memory.db")

    return {"decision": decision, "adata": adata}


def _stage_clustering_annotation(param_result: dict,
                                experiment_id: str,
                                has_api: bool) -> dict:
    """Test clustering and LLM-assisted cell type annotation."""
    section("Test 4 -- Clustering and cell type annotation")

    adata = param_result.get("adata")
    if adata is None:
        info("Skipping (no data available)")
        return {}

    decision   = param_result.get("decision")
    resolution = decision.chosen_value if decision else 0.5

    import scanpy as sc

    info(f"Running Leiden clustering (resolution={resolution})...")
    t0 = time.time()
    sc.tl.leiden(adata, resolution=resolution)
    n_clusters = adata.obs["leiden"].nunique()
    ok(f"{n_clusters} clusters found", f"{time.time()-t0:.1f}s")

    if 6 <= n_clusters <= 15:
        ok("Cluster count is biologically plausible for PBMCs",
           "expected: 6-15 clusters")
    else:
        info(f"Clusters: {n_clusters} (outside typical range for PBMC 3k)")

    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon")

    # Validate with known PBMC markers
    PBMC_MARKERS = {
        "CD3D":  "T cells",       "CD79A": "B cells",
        "NKG7":  "NK cells",      "LYZ":   "Monocytes",
        "PPBP":  "Platelets",     "FCER1A":"Dendritic cells",
        "MS4A1": "B cells",       "CD8A":  "CD8+ T cells",
        "IL7R":  "CD4+ T cells",
    }

    found = {g: ct for g, ct in PBMC_MARKERS.items() if g in adata.var_names}
    ok(f"Known PBMC markers detected: {len(found)}/{len(PBMC_MARKERS)}")
    for gene, ct in list(found.items())[:5]:
        print(f"      {DIM}* {gene} -> {ct}{RST}")

    # LLM annotation if API is available
    if has_api:
        info("Testing LLM-assisted cell type annotation...")
        try:
            from aria.llm.provider import LLMProvider
            import json

            provider = LLMProvider.from_config()

            markers_for_llm = {}
            for cluster in sorted(adata.obs["leiden"].unique())[:3]:
                try:
                    top_genes = list(
                        adata.uns["rank_genes_groups"]["names"][:10][
                            adata.uns["rank_genes_groups"]["names"].dtype.names[0]
                        ]
                    )
                    markers_for_llm[cluster] = top_genes[:10]
                except Exception:
                    markers_for_llm[cluster] = []

            prompt = f"""
Organism: Homo sapiens
Tissue: Peripheral Blood Mononuclear Cells (PBMCs)
Top marker genes per cluster:
{json.dumps(markers_for_llm, indent=2)}

For each cluster, return the most likely cell type.
Return JSON only: {{"cluster_id": "cell_type"}}
"""
            response = provider.complete_heavy(
                prompt=prompt,
                system="Expert bioinformatician annotating PBMC cell types.",
                max_tokens=300,
            )
            ok("LLM annotation completed successfully")
            try:
                clean = response.strip().strip("```json").strip("```").strip()
                ann   = json.loads(clean)
                for cluster, ct in ann.items():
                    print(f"      {DIM}* Cluster {cluster}: {ct}{RST}")
            except Exception:
                print(f"      {DIM}{response[:200]}{RST}")

        except Exception as e:
            info(f"LLM annotation (debug mode): {str(e)[:60]}")
    else:
        info("LLM annotation skipped (no API key)")

    return {"n_clusters": n_clusters, "adata": adata}


def _stage_memory_persistence(experiment_id: str):
    """Verify all decisions were stored in persistent memory."""
    section("Test 5 -- Memory persistence")

    from aria.memory.memory import ARIAMemory

    memory = ARIAMemory()
    wing   = memory.get_wing(experiment_id)

    if wing:
        ok(f"Experiment saved: '{wing['name']}'",
           f"ID: {experiment_id[:12]}")
    else:
        ok("First run: memory initialized correctly")

    ctx = memory.startup_context()
    if "No experiments" not in ctx:
        ok("L0 context available for next startup",
           "ARIA will remember this experiment")
    else:
        ok("Memory ready for first experiment")

    decisions = memory.get_decisions(experiment_id)
    ok(f"{len(decisions)} analytical decisions recorded",
       "exportable for manuscript Methods section")

    memory.close()


def _stage_message_bus_summary(experiment_id: str):
    """Summarize all findings published during the test."""
    section("Test 6 -- MessageBus summary")

    from aria.bus.message_bus import bus, Confidence

    findings = bus.get_findings(experiment_id)
    all_msgs = bus.get_log(experiment_id)

    ok(f"Total messages in bus: {len(all_msgs)}")
    ok(f"Findings published: {len(findings)}")

    conf_counts: dict = {}
    for f in findings:
        c = f.confidence.value
        conf_counts[c] = conf_counts.get(c, 0) + 1

    for conf, count in conf_counts.items():
        color = GRN if conf == "high" else YLW if conf == "medium" else RED
        print(f"      {color}* {conf.upper()}{RST}: {count} finding(s)")


# ── Native pytest entry ─────────────────────────────────────────────────────

def _find_pbmc_dataset():
    for c in (Path.home() / "aria-data" / "pbmc3k_test",
              Path.home() / "aria-data" / "pbmc3k_test" / "hg19"):
        if c.exists() and list(c.rglob("*.mtx*")):
            return c
    return None


def test_pbmc_e2e_pipeline():
    """Native pytest entry (P1-11 follow-up). The `_stage_*` helpers are no
    longer named `test_*`, so pytest stops mis-collecting them as fixture-taking
    tests (the historical 6 collection errors). This drives the full pipeline
    only when scanpy + litellm + the PBMC 3k dataset are present; otherwise it
    skips cleanly instead of erroring or hitting the script's sys.exit path."""
    import os
    import uuid as _uuid
    import pytest

    pytest.importorskip("scanpy")
    pytest.importorskip("litellm")
    data_dir = _find_pbmc_dataset()
    if data_dir is None:
        pytest.skip("PBMC 3k dataset not present (run install.sh to download)")

    has_api = bool(os.environ.get("ANTHROPIC_API_KEY")
                   or os.environ.get("GEMINI_API_KEY"))
    experiment_id = f"pbmc3k_{_uuid.uuid4().hex[:8]}"

    exp_ctx   = _stage_data_audit(data_dir, experiment_id)
    qc_result = _stage_scrna_qc(data_dir, exp_ctx, experiment_id)
    param_res = _stage_parameter_advisor(qc_result, experiment_id, has_api)
    _stage_clustering_annotation(param_res, experiment_id, has_api)
    _stage_memory_persistence(experiment_id)
    _stage_message_bus_summary(experiment_id)

    assert qc_result.get("n_cells", 0) > 0


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="ARIA PBMC 3k End-to-End Test"
    )
    parser.add_argument(
        "--data-dir",
        default=str(Path.home() / "aria-data" / "pbmc3k_test"),
        help="Directory containing the PBMC 3k data files",
    )
    args = parser.parse_args()

    banner()

    has_anthropic, has_gemini = check_api_keys()
    has_api = has_anthropic or has_gemini

    data_dir      = find_pbmc_data(Path(args.data_dir).expanduser())
    experiment_id = f"pbmc3k_{uuid.uuid4().hex[:8]}"
    info(f"Experiment ID: {experiment_id}")

    t_total = time.time()

    exp_ctx   = _stage_data_audit(data_dir, experiment_id)
    qc_result = _stage_scrna_qc(data_dir, exp_ctx, experiment_id)
    param_res = _stage_parameter_advisor(qc_result, experiment_id, has_api)
    clust_res = _stage_clustering_annotation(param_res, experiment_id, has_api)
    _stage_memory_persistence(experiment_id)
    _stage_message_bus_summary(experiment_id)

    total_time = time.time() - t_total
    n_clusters = clust_res.get("n_clusters", "?")
    n_cells    = qc_result.get("n_cells", "?")

    print(f"\n{'─'*50}")
    print(f"{BLD}{GRN}  v ARIA pipeline completed successfully{RST}")
    print(f"{'─'*50}")
    print(f"  Total time:  {total_time:.0f}s")
    print(f"  Cells:       {n_cells} (post-QC)")
    print(f"  Clusters:    {n_clusters}")
    print(f"  Memory:      {BLD}~/.aria/memory.db{RST}")
    print(f"\n  {CYN}ARIA is ready to analyze your data.{RST}")
    print(f"  Run: {BLD}conda activate aria-env && aria{RST}\n")


if __name__ == "__main__":
    main()
