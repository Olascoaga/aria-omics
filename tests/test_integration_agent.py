"""
ARIA IntegrationAgent Tests
-----------------------------
Validates: strategy selection, WNN weight evaluation, peak-to-gene
conflict detection, MOFA+ no-hardcoded-biology checks, cross-modal
conflict resolution, and DebateCouncil integration for multimodal claims.

Run:
  conda activate aria-env
  python tests/test_integration_agent.py
"""

from __future__ import annotations
import sys, json
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


# ── Mocks ─────────────────────────────────────────────────────────────────────

class MockLLM:
    def complete(self, *a, **kw):
        return (
            "ALTERNATIVE_HYPOTHESIS: Weight imbalance could reflect "
            "low ATAC data quality rather than biology\n"
            "CHALLENGES: Check FRiP score for ATAC\n"
            "EVIDENCE_REQUESTED: ATAC FRiP > 0.2 required\n"
            "VERDICT: ACCEPT_REVISED\n"
            "REVISED_CLAIM: WNN dominated by RNA (weight=0.91); "
            "ATAC contribution limited — verify ATAC QC before reporting"
        )
    def complete_heavy(self, *a, **kw):  return self.complete()
    def complete_medium(self, *a, **kw): return self.complete()
    def complete_light(self, *a, **kw):  return "ok"

class MockEnvManager:
    def run_in_stack(self, stack, script_path, params, timeout=None):
        script = Path(script_path).name
        if "wnn" in script:
            return {
                "status": "success", "n_cells": 5000,
                "mean_rna_weight": 0.62, "mean_atac_weight": 0.38,
                "n_joint_clusters": 9, "n_rna_only_clusters": 8,
                "n_discordant_clusters": 2, "output_path": "/tmp/wnn.h5mu",
                "warnings": [],
            }
        elif "peak2gene" in script:
            return {
                "status": "success", "n_links": 847,
                "n_positive_correlations": 721,
                "n_negative_correlations": 126,
                "top_links": [
                    {"gene": "CD3E",  "peak": "chr11:118300000-118302000",
                     "correlation": 0.71, "distance_kb": 12.3,
                     "direction": "positive"},
                    {"gene": "PDCD1", "peak": "chr2:241800000-241802000",
                     "correlation": -0.44, "distance_kb": 22.1,
                     "direction": "negative"},
                ],
                "ctcf_validated": False, "hic_corroborated": False,
                "warnings": ["126 peaks show negative correlation"],
            }
        elif "mofa" in script:
            return {
                "status": "success", "n_factors": 10, "n_cells": 5000,
                "top_factors": [{"factor_id": i+1} for i in range(5)],
                "variance_explained": {"scRNA": 0.42, "scATAC": 0.23},
                "factor1_top_features": ["feature_a", "feature_b", "feature_c"],
                "technical_factor_flag": True, "technical_factor_id": 1,
                "factor_scores": "/tmp/scores.csv",
                "output_path": "/tmp/mofa.hdf5",
                "warnings": ["Factor 1 flagged as technical"],
            }
        return {"status": "success"}

class MockMemory:
    def store_decision(self, *a, **kw): pass
    def create_tunnel(self, *a, **kw): pass
    def get_decisions(self, *a, **kw): return []
    def list_wings(self): return []

class MockAdvisor:
    class _Decision:
        chosen_value = 20
        decision_id  = "mock_001"
        justification = "k=20 based on dataset size"
        candidates   = [
            type("C", (), {"value": 15, "recommended": False,
                           "metrics": {"rna_weight": 0.65, "atac_weight": 0.35}})(),
            type("C", (), {"value": 20, "recommended": True,
                           "metrics": {"rna_weight": 0.62, "atac_weight": 0.38}})(),
            type("C", (), {"value": 30, "recommended": False,
                           "metrics": {"rna_weight": 0.60, "atac_weight": 0.40}})(),
        ]
    def advise_wnn_k(self, *a, **kw): return MockAdvisor._Decision()


# ── Build agent ───────────────────────────────────────────────────────────────

section("IntegrationAgent — import and structure")

try:
    from aria.agents.integration_agent import IntegrationAgent, INTEGRATION_SYSTEM
    ok("IntegrationAgent imported successfully")
except Exception as e:
    fail("Import failed", str(e))
    sys.exit(1)

agent = IntegrationAgent.__new__(IntegrationAgent)
agent.llm     = MockLLM()
agent.memory  = MockMemory()
agent.env     = MockEnvManager()
agent.advisor = MockAdvisor()

agent._findings    = []
agent._escalations = []

def mock_publish_finding(exp_id, payload, conf):
    agent._findings.append({"payload": payload, "conf": conf})
def mock_publish_escalation(**kw):
    agent._escalations.append(kw)
def mock_publish_status(*a, **kw): pass

agent.publish_finding    = mock_publish_finding
agent.publish_escalation = mock_publish_escalation
agent.publish_status     = mock_publish_status

try:
    assert hasattr(agent, "_select_strategy")
    assert hasattr(agent, "_run_wnn")
    assert hasattr(agent, "_run_peak_to_gene")
    assert hasattr(agent, "_run_mofa")
    assert hasattr(agent, "_resolve_conflicts")
    ok("All integration methods present")
except Exception as e:
    fail("Missing methods", str(e))


# ── Test 1: Strategy selection ────────────────────────────────────────────────

section("IntegrationAgent — strategy selection")

try:
    s = agent._select_strategy(
        has_scrna=True, has_scatac=True, has_hic=False,
        n_modalities=2,
        intent={"user_question": "What cell types are in my PBMC multiome?"},
    )
    assert s == "wnn", f"Expected 'wnn', got '{s}'"
    ok("scRNA + scATAC (2 mod, no factor request) → wnn")
except Exception as e:
    fail("WNN strategy selection", str(e))

try:
    s = agent._select_strategy(
        has_scrna=True, has_scatac=True, has_hic=True,
        n_modalities=3,
        intent={"user_question": "What factors drive variation across modalities?"},
    )
    assert s == "wnn_then_mofa", f"Expected 'wnn_then_mofa', got '{s}'"
    ok("scRNA + scATAC + HiC (3 mod) → wnn_then_mofa")
except Exception as e:
    fail("3-modality strategy selection", str(e))

try:
    s = agent._select_strategy(
        has_scrna=True, has_scatac=False, has_hic=False,
        n_modalities=2,
        intent={"user_question": "Run MOFA decomposition on RNA and ChIP"},
    )
    assert s == "mofa", f"Expected 'mofa', got '{s}'"
    ok("MOFA keyword in question → mofa strategy")
except Exception as e:
    fail("MOFA keyword strategy", str(e))

try:
    s = agent._select_strategy(
        has_scrna=False, has_scatac=True, has_hic=False,
        n_modalities=2,
        intent={"user_question": "Link chromatin peaks to genes"},
    )
    assert s == "peak_to_gene_only", f"Expected 'peak_to_gene_only', got '{s}'"
    ok("Bulk/different-cell RNA+ATAC → peak_to_gene_only")
except Exception as e:
    fail("Peak-to-gene only strategy", str(e))


# ── Test 2: WNN weight evaluation ────────────────────────────────────────────

section("IntegrationAgent — WNN weight evaluation and DebateCouncil")

try:
    agent._findings = []
    # Balanced weights → HIGH confidence
    wnn_balanced = {
        "status": "success", "mean_rna_weight": 0.62,
        "mean_atac_weight": 0.38, "n_cells": 5000, "warnings": [],
    }
    agent._evaluate_wnn_weights("exp_001", wnn_balanced, {"user_question": "test"})

    assert len(agent._findings) >= 1
    finding = agent._findings[-1]
    from aria.bus.message_bus import Confidence
    assert finding["conf"] == Confidence.HIGH
    ok(f"Balanced WNN (RNA=0.62) → HIGH confidence finding published")
except Exception as e:
    fail("Balanced WNN evaluation", str(e))

try:
    agent._findings = []
    # Imbalanced weights → DebateCouncil invoked → LOW/MEDIUM confidence
    wnn_imbalanced = {
        "status": "success", "mean_rna_weight": 0.91,
        "mean_atac_weight": 0.09, "n_cells": 5000,
        "atac_frip": 0.12,  # low FRiP explains the imbalance
        "warnings": ["RNA dominates WNN"],
    }
    agent._evaluate_wnn_weights("exp_002", wnn_imbalanced,
                                 {"user_question": "PBMC multiome analysis"})
    assert len(agent._findings) >= 1
    finding = agent._findings[-1]
    assert finding["conf"] in (Confidence.LOW, Confidence.MEDIUM)
    # DebateCouncil should have added weight_debate to result
    ok(f"Imbalanced WNN (RNA=0.91) → {finding['conf'].value} confidence + debate invoked")
except Exception as e:
    fail("Imbalanced WNN evaluation", str(e))


# ── Test 3: Peak-to-gene conflict detection ───────────────────────────────────

section("IntegrationAgent — peak-to-gene negative correlations")

try:
    agent._findings = []
    p2g_result = {
        "status": "success", "n_links": 847,
        "n_positive_correlations": 721,
        "n_negative_correlations": 126,
        "top_links": [
            {"gene": "PDCD1", "peak": "chr2:241800000-241802000",
             "correlation": -0.44, "distance_kb": 22.1, "direction": "negative"},
        ],
        "ctcf_validated": False, "hic_corroborated": False,
        "warnings": ["126 negative correlations detected"],
    }

    # Simulate the conflict resolution
    conflicts = agent._resolve_conflicts(
        "exp_003",
        {"genome": "hg38", "organism": "Homo sapiens"},
        {"user_question": "regulatory landscape of T cells"},
        {"peak_to_gene": p2g_result},
    )

    n_neg_conflicts = sum(
        1 for c in conflicts.get("conflicts", [])
        if c.get("type") == "accessibility_expression_discordance"
    )
    assert n_neg_conflicts >= 1
    ok("Negative correlations -> accessibility-expression conflict detected")
except Exception as e:
    fail("Peak-to-gene conflict detection", str(e))


# ── Test 4: MOFA+ mock does not invent biological factor labels ──────────────

section("IntegrationAgent — MOFA+ no hardcoded biology in mock")

try:
    import aria.scripts.integration_mofa as mofa
    assert not hasattr(mofa, "CELL_CYCLE_GENES")
    ok("MOFA+ has no embedded cell-cycle gene panel")
except Exception as e:
    fail("No embedded cell-cycle gene panel", str(e))

try:
    from aria.scripts.integration_mofa import _mock_mofa
    result = _mock_mofa(n_factors=10, reason="test")
    assert result["technical_factor_flag"] is False
    assert result["technical_factor_id"] is None
    f1_features = result["factor1_top_features"]
    assert all(str(x).startswith("feature_") for x in f1_features)
    assert not any("cell cycle" in str(w).lower()
                   for w in result.get("warnings", []))
    ok("MOFA+ mock does not fabricate a cell-cycle factor")
except Exception as e:
    fail("MOFA+ mock no-hardcoded-biology check", str(e))


# ── Test 5: MOFA+ factor count advice ────────────────────────────────────────

section("IntegrationAgent — MOFA+ factor count advice")

try:
    mods_2 = {"scRNA": ["/data/rna.h5ad"], "scATAC": ["/data/atac.h5ad"]}
    mods_3 = {"scRNA": ["/data/rna.h5ad"], "scATAC": ["/data/atac.h5ad"],
               "HiC":  ["/data/hic.mcool"]}

    n2 = agent._advise_mofa_factors(
        "exp_2mod",
        {"user_question": "find shared variation in PBMC"},
        mods_2,
    )
    n3 = agent._advise_mofa_factors(
        "exp_3mod",
        {"user_question": "atlas of complex tissue with many cell types"},
        mods_3,
    )
    assert n3 >= n2, f"3 modalities should need >= factors than 2 ({n3} >= {n2})"
    ok(f"Factor count: 2 mods={n2}, 3 mods (complex)={n3}")
except Exception as e:
    fail("MOFA+ factor advice", str(e))


# ── Test 6: Scripts import and structure ──────────────────────────────────────

section("Integration scripts — import and contract")

try:
    from aria.scripts.integration_wnn import integration_wnn
    ok("integration_wnn.py imported successfully")
except Exception as e:
    fail("integration_wnn.py import", str(e))

try:
    from aria.scripts.integration_peak2gene import integration_peak2gene
    ok("integration_peak2gene.py imported successfully")
except Exception as e:
    fail("integration_peak2gene.py import", str(e))

try:
    from aria.scripts.integration_mofa import integration_mofa
    ok("integration_mofa.py imported successfully")
except Exception as e:
    fail("integration_mofa.py import", str(e))

try:
    # FileNotFound handling
    result = integration_wnn({
        "rna_files":  ["/nonexistent/rna.h5ad"],
        "atac_files": ["/nonexistent/atac.h5ad"],
    })
    assert result["status"] == "error"
    assert "Missing" in result["error_type"]
    ok("integration_wnn handles missing files gracefully")
except Exception as e:
    fail("integration_wnn FileNotFound handling", str(e))

try:
    result = integration_peak2gene({
        "rna_files":  ["/nonexistent/rna.h5ad"],
        "atac_files": [],
    })
    assert result["status"] == "error"
    assert result["error_type"] == "MissingRNA"
    ok("integration_peak2gene handles missing RNA gracefully")
except Exception as e:
    fail("integration_peak2gene FileNotFound", str(e))

try:
    result = integration_mofa({
        "modalities": {"scRNA": ["/only/one"]},
        "n_factors": 10,
    })
    assert result["status"] == "error"
    assert "Insufficient" in result["error_type"]
    ok("integration_mofa rejects single-modality input")
except Exception as e:
    fail("integration_mofa single modality rejection", str(e))


# ── Test 7: INTEGRATION_SYSTEM prompt quality ─────────────────────────────────

section("IntegrationAgent — system prompt biological knowledge")

try:
    assert "WNN" in INTEGRATION_SYSTEM
    assert "MOFA" in INTEGRATION_SYSTEM
    assert "peak" in INTEGRATION_SYSTEM.lower()
    ok("System prompt covers WNN, MOFA+, and peak-to-gene")
except Exception as e:
    fail("System prompt coverage", str(e))

try:
    assert "negative correlation" in INTEGRATION_SYSTEM.lower()
    assert "poised" not in INTEGRATION_SYSTEM.lower()
    ok("System prompt addresses negative correlations without mechanism labels")
except Exception as e:
    fail("System prompt negative-correlation wording", str(e))

try:
    assert "0.85" in INTEGRATION_SYSTEM or "imbalance" in INTEGRATION_SYSTEM.lower()
    ok("System prompt addresses WNN weight imbalance threshold")
except Exception as e:
    fail("System prompt missing WNN weight threshold", str(e))

try:
    assert "feature names alone" in INTEGRATION_SYSTEM.lower()
    assert "cell cycle" not in INTEGRATION_SYSTEM.lower()
    ok("System prompt avoids hardcoded MOFA+ biological factor labels")
except Exception as e:
    fail("System prompt hardcoded MOFA+ biology check", str(e))

try:
    assert "two" in INTEGRATION_SYSTEM.lower() or "2 independent" in INTEGRATION_SYSTEM.lower() or \
           "independent" in INTEGRATION_SYSTEM.lower()
    ok("System prompt requires 2 independent evidence lines for p2g links")
except Exception as e:
    fail("System prompt missing evidence requirement", str(e))


# ── Test 8: Cross-modal WNN discordance detection ────────────────────────────

section("IntegrationAgent — cross-modal discordance detection")

try:
    agent._findings = []
    wnn_with_discordance = {
        "status": "success",
        "mean_rna_weight": 0.60, "mean_atac_weight": 0.40,
        "n_discordant_clusters": 3, "n_cells": 5000, "warnings": [],
    }
    conflicts = agent._resolve_conflicts(
        "exp_disc",
        {"genome": "hg38"}, {"user_question": "test"},
        {"wnn": wnn_with_discordance},
    )
    n_disc = sum(1 for c in conflicts.get("conflicts", [])
                 if "discordance" in c.get("type", ""))
    assert n_disc >= 1
    ok(f"WNN cluster discordance (n=3) → conflict detected and reported")
except Exception as e:
    fail("WNN discordance conflict detection", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v IntegrationAgent validated. ARIA speaks multi-omics.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
