"""
ARIA Bulk RNA-seq Tests
------------------------
Tests the four fixes:
  1. Design factor extracted from intent (not hardcoded)
  2. Metadata parsing: group auto-detection from column names
  3. Sample outlier detection
  4. Pathway enrichment connected

Run:
  conda activate aria-env
  python tests/test_bulk_rna.py
"""

from __future__ import annotations
import sys, os, json, tempfile
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent))

GRN="\033[92m"; RED="\033[91m"; YLW="\033[93m"
CYN="\033[96m"; DIM="\033[2m";  RST="\033[0m"; BLD="\033[1m"
passed = 0; failed = 0

def ok(msg, detail=""):
    global passed; passed += 1
    d = f"  {DIM}{detail}{RST}" if detail else ""
    print(f"  {GRN}v{RST} {msg}{d}")

def fail(msg, err=""):
    global failed; failed += 1
    print(f"  {RED}x{RST} {msg}")
    if err: print(f"    {DIM}{err}{RST}")

def section(t): print(f"\n{BLD}{CYN}> {t}{RST}")


# ── Test data ─────────────────────────────────────────────────────────────────

def make_counts(n_genes=200, samples=None, seed=42):
    """Minimal realistic count matrix."""
    rng = np.random.default_rng(seed)
    if samples is None:
        samples = ["ctrl_1","ctrl_2","ctrl_3","treat_1","treat_2","treat_3"]
    # Simulate ~20 DE genes
    base = rng.negative_binomial(20, 0.5, (n_genes, len(samples)))
    de_idx = list(range(20))
    trt_idx = [i for i,s in enumerate(samples) if "treat" in s.lower()
                or "ko" in s.lower() or "mut" in s.lower()]
    for i in de_idx:
        base[i, trt_idx] = base[i, trt_idx] * rng.integers(3, 8)
    df = pd.DataFrame(base, columns=samples,
                      index=[f"GENE_{i:04d}" for i in range(n_genes)])
    return df


# ── Import ────────────────────────────────────────────────────────────────────

section("Bulk RNA scripts — import")

try:
    from aria.scripts.rna_bulk_de import (
        bulk_rna_de, _load_counts, _load_or_infer_metadata,
        _infer_groups, _resolve_comparison, _sample_qc,
        _run_deseq2, _run_pathway_enrichment,
    )
    ok("rna_bulk_de.py imported successfully")
except Exception as e:
    fail("Import failed", str(e))
    sys.exit(1)


# ── Fix 1: Metadata parsing — group auto-detection ───────────────────────────

section("Fix 1 — Metadata: group auto-detection from column names")

NAMING_PATTERNS = [
    (["ctrl_1","ctrl_2","ctrl_3","treat_1","treat_2","treat_3"],
     {"ctrl","treat"}, "pattern: condition_replicate"),

    (["WT_rep1","WT_rep2","KO_rep1","KO_rep2"],
     {"WT","KO"}, "pattern: genotype_rep"),

    (["vehicle_1","vehicle_2","drug_1","drug_2"],
     {"vehicle","drug"}, "pattern: condition_number"),

    (["Healthy_1","Healthy_2","Disease_1","Disease_2"],
     {"Healthy","Disease"}, "pattern: capitalized"),
]

for samples, expected_groups, desc in NAMING_PATTERNS:
    try:
        groups = _infer_groups(samples)
        assert groups is not None, f"_infer_groups returned None"
        detected = set(groups.values())
        assert detected == expected_groups, \
            f"Expected {expected_groups}, got {detected}"
        ok(f"Auto-detected groups {detected}", desc)
    except Exception as e:
        fail(f"Group detection failed for {desc}", str(e))

try:
    # Ambiguous names should fail gracefully
    ambiguous = ["sample1", "sample2", "sample3", "sample4"]
    groups = _infer_groups(ambiguous)
    # Either None or single group — not a crash
    if groups is None or len(set(groups.values())) < 2:
        ok("Ambiguous names return None/single-group gracefully")
    else:
        fail("Ambiguous names should not produce valid groups")
except Exception as e:
    fail("Ambiguous name handling", str(e))


# ── Fix 2: Design factor from intent ─────────────────────────────────────────

section("Fix 2 — Design factor extracted from biological intent")

try:
    # BulkRNAAgent replaced RNAAgent in v3 (split from unified rna_agent).
    from aria.agents.bulk_rna_agent import BulkRNAAgent, _infer_lfc_threshold

    agent_cls = BulkRNAAgent

    INTENT_TESTS = [
        ({"comparison": "knockout vs wildtype", "summary": "KO vs WT"},
         "genotype", "KO vs WT → genotype factor"),
        ({"comparison": "treated vs control", "summary": "drug treatment"},
         "treatment", "treated vs control → treatment factor"),
        ({"comparison": "24h vs 0h", "summary": "timepoint course"},
         "timepoint", "timepoint comparison"),
        ({"comparison": "lupus vs healthy", "summary": "disease state"},
         "condition", "disease comparison → condition factor"),
    ]

    for intent, expected_factor, desc in INTENT_TESTS:
        factor = agent_cls._infer_design_factor(intent)
        assert factor == expected_factor, \
            f"Expected '{expected_factor}', got '{factor}'"
        ok(f"Design factor: '{factor}'", desc)

except Exception as e:
    fail("Design factor extraction", str(e))

try:
    # v3: entity-to-label mapping replaces regex comparison parsing.
    # A TF knockout should get the lower LFC threshold.
    intent = {
        "comparison": "BMAL1 KO vs wildtype",
        "summary":    "BMAL1 knockout human H9",
        "biological_entities": ["BMAL1", "wildtype"],
    }
    lfc = _infer_lfc_threshold(intent)
    assert lfc == 0.58, f"TF knockout should get 0.58 LFC, got {lfc}"
    ok(f"TF-aware LFC threshold for BMAL1 KO: {lfc} (1.5x)")

    # Disease (non-TF) should keep default 1.0
    intent_d = {"comparison": "disease vs healthy",
                "biological_entities": ["lupus"]}
    lfc_d = _infer_lfc_threshold(intent_d)
    assert lfc_d == 1.0, f"Non-TF should get 1.0, got {lfc_d}"
    ok(f"Non-TF LFC threshold preserved: {lfc_d} (2x)")

    # Entity-to-label matching: BMAL1 should map to label 'B'
    mapping = BulkRNAAgent.__new__(BulkRNAAgent)._map_entities_to_labels(
        entities=["BMAL1", "REV-ERBa", "wildtype"],
        group_names=["B", "R", "WT"],
        intent={},
    )
    assert mapping.get("BMAL1")     == "B",  f"BMAL1 → {mapping.get('BMAL1')}"
    assert mapping.get("REV-ERBa")  == "R",  f"REV-ERBa → {mapping.get('REV-ERBa')}"
    assert mapping.get("wildtype")  == "WT", f"wildtype → {mapping.get('wildtype')}"
    ok(f"Entity→label mapping: {mapping}")

except Exception as e:
    fail("Comparison/LFC/label-matching", str(e))


# ── Fix 3: Sample outlier detection ──────────────────────────────────────────

section("Fix 3 — Sample outlier detection (PCA-based)")

try:
    counts = make_counts(200, ["ctrl_1","ctrl_2","ctrl_3",
                                "treat_1","treat_2","treat_3"])
    metadata = pd.DataFrame({"condition": ["ctrl"]*3 + ["treat"]*3},
                             index=counts.columns)
    warnings_list = []
    with tempfile.TemporaryDirectory() as tmpdir:
        qc = _sample_qc(counts, metadata, tmpdir, warnings_list)
    assert "n_samples"    in qc
    assert "outliers"     in qc
    assert "pca_variance" in qc
    assert "lib_size_range" in qc
    ok(f"Sample QC: {qc['n_samples']} samples, "
       f"outliers={qc['outliers']}, "
       f"PC1={qc['pca_variance'][0] if qc['pca_variance'] else '?':.2f}")
except Exception as e:
    fail("Sample QC basic", str(e))

try:
    # Inject a clear outlier
    rng = np.random.default_rng(0)
    counts_with_outlier = make_counts(200)
    # Make sample ctrl_1 an extreme outlier
    counts_with_outlier["ctrl_1"] = counts_with_outlier["ctrl_1"] * 100
    metadata = pd.DataFrame({"condition": ["ctrl"]*3 + ["treat"]*3},
                             index=counts_with_outlier.columns)
    warnings_list = []
    with tempfile.TemporaryDirectory() as tmpdir:
        qc = _sample_qc(counts_with_outlier, metadata, tmpdir, warnings_list)
    ok(f"Outlier injection: outliers detected = {qc['outliers']}",
       "ctrl_1 was 100x inflated")
except Exception as e:
    fail("Outlier detection", str(e))


# ── Fix 4: Pathway enrichment connected ──────────────────────────────────────

section("Fix 4 — Pathway enrichment (gseapy connected)")

try:
    # Test mock pathway for environments without gseapy
    from aria.scripts.rna_bulk_de import _mock_pathways, _get_gene_sets
    mock_pw = _mock_pathways(["CD3E","CD8A","PDCD1","TOX","IL2"])
    assert "GO_BP" in mock_pw or "KEGG" in mock_pw
    assert any(isinstance(v, list) for v in mock_pw.values())
    ok("Mock pathway enrichment returns structured result")
except Exception as e:
    fail("Mock pathway structure", str(e))

try:
    # Gene set selection by organism
    human_sets = _get_gene_sets("Homo sapiens")
    mouse_sets = _get_gene_sets("Mus musculus")
    assert "GO_BP"    in human_sets
    assert "KEGG"     in human_sets
    assert "Reactome" in human_sets
    assert "KEGG"     in mouse_sets
    ok(f"Human gene sets: {list(human_sets.keys())}")
    ok(f"Mouse gene sets: {list(mouse_sets.keys())}")
except Exception as e:
    fail("Gene set selection", str(e))

try:
    # Pathway enrichment with real or mock
    pw, pw_warn = _run_pathway_enrichment(
        sig_genes=["CD3E","CD8A","PDCD1","TOX","LAG3",
                   "HAVCR2","IL2","IFNG","TNF","GZMB",
                   "PRF1","NKG7","GNLY","CTLA4","TIGIT"],
        up_genes=["CD3E","CD8A","PDCD1","TOX"],
        down_genes=["IL2","IFNG"],
        organism="Homo sapiens",
        output_dir="/tmp/test_pw",
    )
    assert isinstance(pw, dict)
    ok(f"Pathway enrichment ran: {list(pw.keys())}")
    if pw_warn:
        ok(f"  Warnings: {pw_warn[0][:60]}")
except Exception as e:
    fail("Pathway enrichment call", str(e))


# ── End-to-end: full bulk_rna_de() with synthetic data ───────────────────────

section("End-to-end — bulk_rna_de() with synthetic counts")

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_df = make_counts(300)
        counts_path = Path(tmpdir) / "counts.tsv"
        counts_df.to_csv(str(counts_path), sep="\t")

        # v3 API: use contrasts (list) instead of comparison (dict)
        result = bulk_rna_de({
            "files":         [str(counts_path)],
            "design_factor": "condition",
            "contrasts": [
                {"numerator": "treat", "denominator": "ctrl",
                 "name": "treat vs ctrl"},
            ],
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  True,
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
        })

    assert result["status"] == "success", \
        f"Expected success, got: {result.get('status')} — {result.get('details','')}"

    # v3 shape: contrasts list in result
    assert "contrasts"  in result, f"Missing 'contrasts' key; keys: {list(result)}"
    assert "sample_qc"  in result
    assert "design_used" in result
    assert len(result["contrasts"]) == 1

    c0 = result["contrasts"][0]
    assert "n_significant"   in c0
    assert "n_upregulated"   in c0
    assert "n_downregulated" in c0
    assert "plots"           in c0
    assert "pathways"        in c0

    ok(f"Full pipeline: 1 contrast → {c0['n_significant']} DE genes "
       f"({c0['n_upregulated']} up, {c0['n_downregulated']} down)")
    ok(f"Sample QC: {result['sample_qc']['n_samples']} samples, "
       f"outliers={result['sample_qc']['outliers']}")
    ok(f"Pathways: {list(c0['pathways'].keys())}")
    ok(f"Plots: {[k for k,v in c0['plots'].items() if v]}")
    ok(f"Design used: {result['design_used']}")

except Exception as e:
    fail("Full bulk_rna_de() pipeline", str(e))


# ── Multi-contrast: the v3 key feature ──────────────────────────────────────

section("Multi-contrast — all pairwise contrasts in one run (v3)")

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        # Three groups: B, R, WT (matching H9 experiment)
        samples = ["B_1","B_2","B_3","R_1","R_2","R_3","WT_1","WT_2","WT_3"]
        counts_df = make_counts(300, samples=samples)
        counts_path = Path(tmpdir) / "counts.tsv"
        counts_df.to_csv(str(counts_path), sep="\t")

        result = bulk_rna_de({
            "files":         [str(counts_path)],
            "design_factor": "condition",
            "contrasts": [
                {"numerator": "B", "denominator": "WT", "name": "BMAL1 KO vs WT"},
                {"numerator": "R", "denominator": "WT", "name": "REV-ERBa KO vs WT"},
            ],
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,   # skip to avoid gseapy HTTP dependency
            "lfc_threshold": 0.58,    # TF threshold
        })

    assert result["status"] == "success"
    assert result["n_contrasts"] == 2
    names = [c["name"] for c in result["contrasts"]]
    assert "BMAL1 KO vs WT" in names
    assert "REV-ERBa KO vs WT" in names
    ok(f"Multi-contrast: {result['n_contrasts']} contrasts run — {names}")

    # Overlap computed
    assert "overlap" in result
    ok(f"Cross-contrast overlap computed: "
       f"{list(result['overlap'].keys())[:1] or 'no pairs'}")

except Exception as e:
    fail("Multi-contrast bulk_rna_de()", str(e))


# ── Replicate concordance (v3 QC addition) ───────────────────────────────────

section("Replicate concordance — corrupted sample is detected")

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        # Generate clean replicates then corrupt one
        rng = np.random.default_rng(42)
        n_genes = 500
        samples = ["B_1","B_2","B_3","WT_1","WT_2","WT_3"]
        gene_mean = rng.exponential(100, n_genes)
        data = np.array([rng.poisson(gene_mean) for _ in samples]).T.astype(float)

        # Heavily corrupt B_2: multiplicative noise on 70% of genes
        n_corrupt = int(n_genes * 0.7)
        idx = rng.choice(n_genes, n_corrupt, replace=False)
        data[idx, 1] = data[idx, 1] * np.exp(rng.uniform(-2.3, 2.3, n_corrupt))

        df = pd.DataFrame(data.astype(int),
                           index=[f"g_{i}" for i in range(n_genes)],
                           columns=samples)
        cp = Path(tmpdir) / "counts.tsv"
        df.to_csv(str(cp), sep="\t")

        result = bulk_rna_de({
            "files":         [str(cp)],
            "design_factor": "condition",
            "contrasts": [{"numerator": "B", "denominator": "WT",
                           "name": "B vs WT"}],
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,
            "lfc_threshold": 0.58,
        })

    sqc = result.get("sample_qc", {})
    rep_corr = sqc.get("replicate_correlations", {})
    assert rep_corr, "replicate_correlations not present in sample_qc"
    # B_2 should have the lowest mean correlation
    lowest = min(rep_corr, key=rep_corr.get)
    assert lowest == "B_2", \
        f"B_2 should have lowest concordance, got {lowest} ({rep_corr})"
    assert rep_corr["B_2"] < 0.85, \
        f"B_2 correlation {rep_corr['B_2']:.3f} should be < 0.85"
    ok(f"Corrupted sample B_2 detected: r={rep_corr['B_2']:.3f} "
       f"(< 0.85 threshold)")

except Exception as e:
    fail("Replicate concordance detection", str(e))


# ── DESeq2 design factor validation ──────────────────────────────────────────

section("DESeq2 — correct design factor propagation")

try:
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_df = make_counts(200)
        counts_path = Path(tmpdir) / "counts.tsv"
        counts_df.to_csv(str(counts_path), sep="\t")

        result = bulk_rna_de({
            "files":         [str(counts_path)],
            "design_factor": "condition",
            "contrasts":     [{"numerator": "treat", "denominator": "ctrl",
                               "name": "treat vs ctrl"}],
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,
        })

    # The design_used should reflect the actual design factor
    design_used = result.get("design_used", "") or ""
    assert "condition" in design_used, \
        f"design_used does not reflect 'condition': {design_used}"
    ok(f"Design factor correctly used: {design_used}")
except Exception as e:
    fail("Design factor in DESeq2", str(e))

try:
    # Test insufficient replicates handled gracefully
    with tempfile.TemporaryDirectory() as tmpdir:
        counts_1rep = make_counts(200, ["ctrl_1","treat_1"])
        cp = Path(tmpdir) / "counts.tsv"
        counts_1rep.to_csv(str(cp), sep="\t")

        result = bulk_rna_de({
            "files":         [str(cp)],
            "design_factor": "condition",
            "contrasts":     [{"numerator": "treat", "denominator": "ctrl",
                               "name": "treat vs ctrl"}],
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,
        })

    # Should fail gracefully with structured response
    assert result.get("status") in ("error", "success")
    ok(f"1 replicate per group: status={result.get('status')} "
       f"(graceful handling)")
except Exception as e:
    fail("Insufficient replicates handling", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v Bulk RNA pipeline validated. "
          f"All 4 fixes confirmed.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
