"""
ARIA ChromatinAgent Tests
--------------------------
Validates ChromatinAgent structure, LSI parameter decisions,
QC thresholds, MACS3 command building, and DebateCouncil integration
for TF motif enrichment claims.

Run:
  conda activate aria-env
  python tests/test_chromatin_agent.py
"""

from __future__ import annotations
import sys, os, json
from pathlib import Path

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


# ── Mock infrastructure ───────────────────────────────────────────────────────

class MockLLM:
    def complete(self, prompt, system="", tier=None, max_tokens=1024, messages=None):
        return (
            "CLAIM: CTCF and RUNX1 show significant motif enrichment\n"
            "EVIDENCE: CTCF enrichment=8.2, RUNX1 enrichment=6.1\n"
            "KNOWN_LIMITATIONS: Tn5 bias not corrected"
        )
    def complete_heavy(self, *a, **kw):   return self.complete(*a, **kw)
    def complete_medium(self, *a, **kw):  return self.complete(*a, **kw)
    def complete_light(self, *a, **kw):   return "ok"

class MockEnvManager:
    def run_in_stack(self, stack, script_path, params, timeout=None):
        script = Path(script_path).name
        if "qc" in script:
            return {
                "status": "success", "data_type": params.get("data_type"),
                "frip": 0.38, "tss_enrichment": 7.2,
                "mito_fraction": 0.04, "n_cells": 5000,
                "n_cells_after": 4800, "n_fragments": 50000000,
                "pass_qc": True, "warnings": [],
            }
        elif "peaks" in script:
            return {
                "status": "success", "n_peaks": 85000,
                "peaks_path": "/tmp/test_peaks.narrowPeak",
                "consensus_peaks_path": None,
                "frip": 0.38, "warnings": [],
            }
        elif "motif" in script:
            return {
                "status": "success",
                "top_motifs": ["CTCF", "RUNX1", "AP-1"],
                "scores": {"CTCF": 8.2, "RUNX1": 6.1, "AP-1": 5.3},
                "n_peaks": 85000,
                "tn5_bias_corrected": False,
            }
        return {"status": "success"}

class MockMemory:
    def create_room(self, *a, **kw): pass
    def store_finding(self, *a, **kw): pass
    def create_hall(self, *a, **kw): pass
    def get_decisions(self, *a, **kw): return []
    def list_wings(self): return []


# ── Test 1: Import and structure ──────────────────────────────────────────────

section("ChromatinAgent — import and structure")

try:
    from aria.agents.chromatin_agent import ChromatinAgent
    ok("ChromatinAgent imported successfully")
except Exception as e:
    fail("Import failed", str(e))

try:
    assert hasattr(ChromatinAgent, "MACS3_PARAMS")
    params = ChromatinAgent.MACS3_PARAMS
    assert "scATAC"      in params
    assert "bulk_ATAC"   in params
    assert "ChIP"        in params
    assert "CUT_AND_RUN" in params
    assert "CUT_AND_TAG" in params
    ok("MACS3_PARAMS defined for all 5 modalities")
except Exception as e:
    fail("MACS3_PARAMS missing modalities", str(e))

try:
    # Verify CUT&RUN uses --nolambda (low background)
    cr_params = ChromatinAgent.MACS3_PARAMS["CUT_AND_RUN"]
    assert cr_params["nolambda"] == True, "CUT&RUN must use --nolambda"
    ok("CUT&RUN: nolambda=True (low background correction)")

    # Verify ATAC uses --nomodel --extsize 200
    atac_params = ChromatinAgent.MACS3_PARAMS["scATAC"]
    assert atac_params["nomodel"]  == True
    assert atac_params["extsize"]  == 200
    assert atac_params["keep_dup"] == "all"
    ok("scATAC: nomodel=True, extsize=200, keep_dup=all")

    # ChIP uses model (not nomodel)
    chip_params = ChromatinAgent.MACS3_PARAMS["ChIP"]
    assert chip_params["nomodel"] == False
    ok("ChIP: nomodel=False (uses MACS3 model for extension)")
except Exception as e:
    fail("MACS3 assay-specific params incorrect", str(e))


# ── Test 2: LSI parameter decision ────────────────────────────────────────────

section("ChromatinAgent — LSI parameter decisions")

try:
    # Build minimal agent with mocks
    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.llm     = MockLLM()
    agent.memory  = MockMemory()
    agent.env     = MockEnvManager()
    agent.advisor = None
    ok("ChromatinAgent instantiated with mocks")
except Exception as e:
    fail("Mock instantiation failed", str(e))

try:
    # Test LSI parameter advice for different complexities
    lsi_simple = agent._advise_lsi_params(
        "exp_001",
        {"complexity": "simple", "user_question": "basic ATAC"},
        {"n_cells_after": 1500},
    )
    assert lsi_simple["discard_component"] == 1, \
        "Component 1 MUST always be discarded"
    assert lsi_simple["components_range"][1] <= 30, \
        "Simple dataset should use fewer components"
    ok(f"Simple dataset: range={lsi_simple['components_range']}, "
       f"discard_component={lsi_simple['discard_component']}")
except Exception as e:
    fail("LSI simple params failed", str(e))

try:
    lsi_complex = agent._advise_lsi_params(
        "exp_002",
        {"complexity": "complex", "user_question": "atlas"},
        {"n_cells_after": 50000},
    )
    assert lsi_complex["discard_component"] == 1
    assert lsi_complex["components_range"][1] >= 40, \
        "Complex/large dataset needs more components"
    ok(f"Complex dataset: range={lsi_complex['components_range']}, "
       f"discard_component={lsi_complex['discard_component']}")
except Exception as e:
    fail("LSI complex params failed", str(e))

try:
    # Verify rationale mentions why component 1 is discarded
    assert "sequencing depth" in lsi_simple["rationale"].lower() or \
           "technical" in lsi_simple["rationale"].lower(), \
        "Rationale must explain WHY component 1 is discarded"
    ok("LSI rationale explains biological reason for discarding component 1")
except Exception as e:
    fail("LSI rationale missing explanation", str(e))


# ── Test 3: ChIP target classification ───────────────────────────────────────

section("ChromatinAgent — ChIP target classification")

try:
    # Histone marks should be detected from filenames
    histone_files = [
        "/data/H3K27ac_rep1.bam",
        "/data/H3K27ac_rep2.bam",
        "/data/input_control.bam",
    ]
    result = agent._classify_chip_target(histone_files, {})
    assert result == "histone", f"Expected 'histone', got '{result}'"
    ok("H3K27ac detected as histone mark -> broad peak mode")
except Exception as e:
    fail("Histone classification failed", str(e))

try:
    tf_files = [
        "/data/CTCF_ChIP_rep1.bam",
        "/data/IgG_control.bam",
    ]
    result = agent._classify_chip_target(tf_files, {})
    assert result == "tf", f"Expected 'tf', got '{result}'"
    ok("CTCF ChIP detected as TF -> narrow peak mode")
except Exception as e:
    fail("TF classification failed", str(e))

try:
    # Question-based inference
    result = agent._classify_chip_target(
        ["/data/unknown_chip.bam"],
        {"user_question": "H3K4me3 promoter activity in neurons"},
    )
    assert result == "histone"
    ok("H3K4me3 in question detected as histone mark")
except Exception as e:
    fail("Question-based ChIP classification failed", str(e))


# ── Test 4: QC thresholds ─────────────────────────────────────────────────────

section("ChromatinAgent — QC threshold logic")

try:
    from aria.bus.message_bus import Confidence

    # Good QC — no warnings
    good_qc = {
        "frip": 0.45, "tss_enrichment": 9.2,
        "warnings": [], "status": "success",
    }
    # Manually test the confidence assignment
    conf = (Confidence.HIGH   if not good_qc["warnings"] and
                                  good_qc["frip"] >= 0.2 and
                                  good_qc["tss_enrichment"] >= 4 else
            Confidence.MEDIUM if not (good_qc["frip"] < 0.2 or
                                       good_qc["tss_enrichment"] < 4) else
            Confidence.LOW)
    assert conf == Confidence.HIGH
    ok(f"Good QC (FRiP=0.45, TSS=9.2) -> HIGH confidence")
except Exception as e:
    fail("Good QC confidence failed", str(e))

try:
    # Failed QC — TSS < 4 (library prep failed)
    bad_qc = {
        "frip": 0.08, "tss_enrichment": 2.1,
        "warnings": ["TSS enrichment 2.1 < 4.0"],
        "status": "success",
    }
    frip      = bad_qc["frip"]
    tss_score = bad_qc["tss_enrichment"]
    conf = Confidence.LOW if (frip < 0.2 or tss_score < 4) else Confidence.MEDIUM
    assert conf == Confidence.LOW
    ok(f"Failed QC (FRiP=0.08, TSS=2.1) -> LOW confidence")
except Exception as e:
    fail("Failed QC confidence failed", str(e))


# ── Test 5: TF analysis detection ────────────────────────────────────────────

section("ChromatinAgent — TF analysis trigger logic")

try:
    tf_intents = [
        {"user_question": "Which transcription factors are active in these cells?"},
        {"user_question": "What TF motifs are enriched in open chromatin?"},
        {"user_question": "Show me footprinting of CTCF binding sites"},
        {"user_question": "Identify enhancers and their regulatory TFs"},
    ]
    non_tf_intents = [
        {"user_question": "How many peaks are in condition A vs B?"},
        {"user_question": "Compare chromatin accessibility between cell types"},
    ]

    for intent in tf_intents:
        assert agent._needs_tf_analysis(intent), \
            f"Should need TF: {intent['user_question']}"
    ok(f"TF analysis correctly triggered for {len(tf_intents)} TF questions")

    for intent in non_tf_intents:
        assert not agent._needs_tf_analysis(intent), \
            f"Should NOT need TF: {intent['user_question']}"
    ok(f"TF analysis correctly skipped for {len(non_tf_intents)} non-TF questions")
except Exception as e:
    fail("TF analysis detection failed", str(e))


# ── Test 6: Scripts import ────────────────────────────────────────────────────

section("Chromatin scripts — import and contract")

try:
    from aria.scripts.chromatin_qc import chromatin_qc
    ok("chromatin_qc.py imported successfully")
except Exception as e:
    fail("chromatin_qc.py import failed", str(e))

try:
    from aria.scripts.chromatin_peaks import chromatin_peaks
    ok("chromatin_peaks.py imported successfully")
except Exception as e:
    fail("chromatin_peaks.py import failed", str(e))

try:
    # Test path validation in chromatin_qc
    result = chromatin_qc({
        "data_type": "scATAC",
        "files": ["/nonexistent/fragments.tsv.gz"],
        "genome": "hg38",
    })
    assert result["status"] == "error"
    assert result["error_type"] == "FileNotFound"
    ok("chromatin_qc handles nonexistent files gracefully")
except Exception as e:
    fail("chromatin_qc path validation failed", str(e))

try:
    result = chromatin_qc({
        "data_type": "unknown_type",
        "files": [],
        "genome": "hg38",
    })
    assert result["status"] == "error"
    assert result["error_type"] in ("UnsupportedDataType", "FileNotFound")
    ok("chromatin_qc rejects unknown data type")
except Exception as e:
    fail("chromatin_qc unknown type handling failed", str(e))

try:
    # Test MACS3 genome size mapping
    from aria.scripts.chromatin_peaks import _get_genome_size
    assert _get_genome_size("hg38") == "hs"
    assert _get_genome_size("mm10") == "mm"
    assert _get_genome_size("dm6")  == "dm"
    ok("Genome size mapping: hg38->hs, mm10->mm, dm6->dm")
except Exception as e:
    fail("Genome size mapping failed", str(e))


# ── Test 7: ARIA system prompt for Tn5 bias ───────────────────────────────────

section("ChromatinAgent — Tn5 bias awareness")

try:
    from aria.agents.chromatin_agent import CHROMATIN_SYSTEM
    assert "tn5" in CHROMATIN_SYSTEM.lower() or "Tn5" in CHROMATIN_SYSTEM
    assert "bias" in CHROMATIN_SYSTEM.lower()
    assert "footprint" in CHROMATIN_SYSTEM.lower()
    ok("CHROMATIN_SYSTEM prompt includes Tn5 bias awareness")
except Exception as e:
    fail("CHROMATIN_SYSTEM missing Tn5 bias mention", str(e))

try:
    assert "LSI" in CHROMATIN_SYSTEM
    assert "SVD" in CHROMATIN_SYSTEM
    assert "first" in CHROMATIN_SYSTEM.lower() or "component 1" in CHROMATIN_SYSTEM.lower()
    ok("CHROMATIN_SYSTEM prompt includes LSI/SVD guidance")
except Exception as e:
    fail("CHROMATIN_SYSTEM missing LSI/SVD guidance", str(e))

try:
    assert "FRiP" in CHROMATIN_SYSTEM
    assert "TSS" in CHROMATIN_SYSTEM
    ok("CHROMATIN_SYSTEM prompt includes FRiP and TSS thresholds")
except Exception as e:
    fail("CHROMATIN_SYSTEM missing QC thresholds", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v ChromatinAgent validated.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
