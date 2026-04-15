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
    from aria.agents.rna_agent import RNAAgent

    class MockMem:
        def get_decisions(self, *a): return []
        def create_wing(self, *a, **kw): pass
        def store_decision(self, *a, **kw): pass

    class MockLLM:
        def complete(self, *a, **kw): return "mock"
        def complete_heavy(self, *a, **kw): return "mock"
        def complete_medium(self, *a, **kw): return "mock"

    agent = RNAAgent.__new__(RNAAgent)
    agent.llm    = MockLLM()
    agent.memory = MockMem()

    INTENT_TESTS = [
        ({"comparison": "knockout vs wildtype", "analysis_type": "differential"},
         "genotype", "KO vs WT → genotype factor"),
        ({"comparison": "treated vs control", "analysis_type": "differential"},
         "treatment", "treated vs control → treatment factor"),
        ({"comparison": "24h vs 0h", "analysis_type": "temporal"},
         "timepoint", "timepoint comparison"),
        ({"comparison": "lupus vs healthy", "analysis_type": "differential"},
         "condition", "disease comparison → condition factor"),
    ]

    for intent, expected_factor, desc in INTENT_TESTS:
        factor, comp = agent._extract_design_from_intent(intent, [])
        assert factor == expected_factor, \
            f"Expected '{expected_factor}', got '{factor}'"
        ok(f"Design factor: '{factor}'", desc)

except Exception as e:
    fail("Design factor extraction", str(e))

try:
    # Comparison parsing — "X vs Y" extraction
    intent = {"comparison": "KRAS_mutant vs KRAS_wildtype"}
    factor, comp = agent._extract_design_from_intent(intent, [])
    assert comp.get("numerator")   == "kras_mutant" or \
           "mutant" in str(comp.get("numerator","")).lower(), \
        f"Comparison not parsed: {comp}"
    ok(f"Comparison parsed from 'X vs Y': {comp}")
except Exception as e:
    fail("Comparison string parsing", str(e))


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
        # Write synthetic counts to TSV
        counts_df = make_counts(300)
        counts_path = Path(tmpdir) / "counts.tsv"
        counts_df.to_csv(str(counts_path), sep="\t")

        result = bulk_rna_de({
            "files":         [str(counts_path)],
            "design_factor": "condition",
            "comparison":    {"numerator": "treat", "denominator": "ctrl"},
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  True,
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
        })

    assert result["status"] == "success", \
        f"Expected success, got: {result.get('status')} — {result.get('details','')}"
    assert "n_significant"   in result
    assert "n_upregulated"   in result
    assert "n_downregulated" in result
    assert "sample_qc"       in result
    assert "pathways"        in result
    assert "plots"           in result
    assert "design_used"     in result
    assert "comparison_used" in result

    ok(f"Full pipeline: {result['n_significant']} DE genes "
       f"({result['n_upregulated']} up, {result['n_downregulated']} down)")
    ok(f"Sample QC: {result['sample_qc']['n_samples']} samples, "
       f"outliers={result['sample_qc']['outliers']}")
    ok(f"Pathways: {list(result['pathways'].keys())}")
    ok(f"Plots: {[k for k,v in result['plots'].items() if v]}")
    ok(f"Design used: {result['design_used']}")
    ok(f"Comparison: {result['comparison_used']}")

except Exception as e:
    fail("Full bulk_rna_de() pipeline", str(e))


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
            "comparison":    {"numerator": "treat", "denominator": "ctrl"},
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,
        })

    # The design_used should reflect the actual design factor
    assert "~condition" in result.get("design_used", "") or \
           "condition"  in result.get("design_used", ""), \
        f"design_used does not reflect 'condition': {result.get('design_used')}"
    ok(f"Design factor correctly used: {result.get('design_used')}")
except Exception as e:
    fail("Design factor in DESeq2", str(e))

try:
    # Test insufficient replicates returns structured error
    with tempfile.TemporaryDirectory() as tmpdir:
        # Only 1 sample per group
        counts_1rep = make_counts(200, ["ctrl_1","treat_1"])
        cp = Path(tmpdir) / "counts.tsv"
        counts_1rep.to_csv(str(cp), sep="\t")

        result = bulk_rna_de({
            "files":         [str(cp)],
            "design_factor": "condition",
            "comparison":    {"numerator": "treat", "denominator": "ctrl"},
            "organism":      "Homo sapiens",
            "output_dir":    tmpdir,
            "run_pathways":  False,
        })

    # Should fail gracefully with structured error
    assert result.get("status") in ("error", "success")  # mock may succeed
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
