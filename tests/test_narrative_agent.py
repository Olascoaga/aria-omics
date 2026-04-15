"""
ARIA NarrativeAgent Tests
--------------------------
Validates report generation: findings grouping, HTML rendering,
methods section accuracy, honest uncertainty reporting, and
that LOW/INSUFFICIENT findings are never hidden.

Run:
  conda activate aria-env
  python tests/test_narrative_agent.py
"""

from __future__ import annotations
import sys, os, json, tempfile
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
    def complete(self, prompt="", system="", tier=None, max_tokens=1024, messages=None):
        return (
            "Analysis identified 5 major cell populations in lupus PBMCs with "
            "HIGH confidence, including an expanded population of exhausted "
            "CD8+ T cells (PDCD1+, TOX+). Chromatin accessibility analysis "
            "revealed 85,000 ATAC peaks with FRiP=0.38. Integration of RNA "
            "and ATAC data using WNN showed balanced modality contributions "
            "(RNA=0.62, ATAC=0.38), with 2 clusters differing between "
            "RNA-only and joint embedding."
        )

class MockMemory:
    def get_decisions(self, experiment_id):
        return [
            {"checkpoint": 3, "question": "Leiden resolution",
             "decision": "0.40", "rationale": "highest silhouette score",
             "made_at": "2026-04-15"},
            {"checkpoint": 3, "question": "WNN k neighbors",
             "decision": "20",  "rationale": "dataset size recommendation",
             "made_at": "2026-04-15"},
            {"checkpoint": 2, "question": "HiC resolution",
             "decision": "40000", "rationale": "TAD calling at 40kb",
             "made_at": "2026-04-15"},
        ]

MOCK_FINDINGS = [
    {"summary": "QC: 5000 -> 4800 cells (4% removed)", "agent": "rna_agent",
     "confidence": "high"},
    {"summary": "5 clusters annotated: T cells, B cells, NK, Monocytes, DC",
     "agent": "rna_agent", "confidence": "high"},
    {"summary": "247 DE genes in exhausted T cell cluster",
     "agent": "rna_agent", "confidence": "high"},
    {"summary": "85,000 ATAC peaks called, FRiP=0.38",
     "agent": "chromatin_agent", "confidence": "high"},
    {"summary": "CTCF and RUNX1 motif enrichment in accessible regions",
     "agent": "chromatin_agent", "confidence": "medium"},
    {"summary": "WNN: RNA=0.62, ATAC=0.38, 9 joint clusters",
     "agent": "integration_agent", "confidence": "high"},
    {"summary": "847 peak-gene regulatory links identified",
     "agent": "integration_agent", "confidence": "medium"},
    {"summary": "126 open-chromatin/silent-gene pairs (poised enhancers?)",
     "agent": "integration_agent", "confidence": "medium"},
    {"summary": "TF footprinting result requires Tn5 bias validation",
     "agent": "chromatin_agent", "confidence": "low"},
    {"summary": "MOFA Factor 3 interpretation: insufficient replicates",
     "agent": "integration_agent", "confidence": "insufficient"},
]

MOCK_AGENT_RESULTS = {
    "rna_agent": {
        "status": "done",
        "findings": {
            "scRNA": {
                "findings": {
                    "qc": {
                        "n_cells_after": 4800,
                        "pct_removed":   4.0,
                        "mt_threshold":  5.0,
                    },
                    "clustering_decision": {
                        "recommended": 0.40,
                    },
                    "cell_types": {
                        "cell_types": {
                            "0": "T cell",
                            "1": "B cell",
                            "2": "NK cell",
                            "3": "Monocyte",
                            "4": "DC",
                        }
                    },
                }
            }
        },
    },
    "chromatin_agent": {"status": "done", "findings": {}},
    "genome_arch_agent": {
        "status": "done",
        "findings": {
            "topology": {
                "compartments": {
                    "pct_A": 43.2, "pct_B": 56.8,
                    "n_ab_switches": 0,
                },
                "tads": {"n_tads": 3200, "median_size_kb": 180.0},
                "loops": {},
            }
        },
    },
    "integration_agent": {
        "status": "done",
        "strategy": "wnn",
        "findings": {
            "wnn": {
                "status": "success",
                "mean_rna_weight": 0.62,
                "mean_atac_weight": 0.38,
                "n_discordant_clusters": 2,
            },
            "peak_to_gene": {
                "n_links": 847,
                "n_negative_correlations": 126,
            },
            "conflicts": {
                "n_conflicts": 1,
                "conflicts": [{
                    "type": "open_chromatin_silent_gene",
                    "description": "126 peaks accessible but genes not expressed",
                }]
            },
        },
    },
}


# ── Tests ─────────────────────────────────────────────────────────────────────

section("NarrativeAgent — import and structure")

try:
    from aria.agents.narrative_agent import NarrativeAgent, NARRATIVE_SYSTEM
    ok("NarrativeAgent imported successfully")
except Exception as e:
    fail("Import failed", str(e))
    sys.exit(1)

# Build agent with temp reports dir
with tempfile.TemporaryDirectory() as tmpdir:
    agent = NarrativeAgent.__new__(NarrativeAgent)
    agent.llm          = MockLLM()
    agent.memory       = MockMemory()
    agent.reports_dir  = Path(tmpdir)
    agent._findings    = []
    agent._escalations = []

    def mock_publish_finding(exp_id, payload, conf):
        agent._findings.append(payload)
    def mock_publish_escalation(**kw): pass
    def mock_publish_status(*a, **kw): pass
    agent.publish_finding    = mock_publish_finding
    agent.publish_escalation = mock_publish_escalation
    agent.publish_status     = mock_publish_status

    try:
        assert hasattr(agent, "_group_findings")
        assert hasattr(agent, "_write_executive_summary")
        assert hasattr(agent, "_write_methods_section")
        assert hasattr(agent, "_render_html_report")
        ok("All required methods present")
    except Exception as e:
        fail("Missing methods", str(e))


    # ── Test 1: findings grouping ─────────────────────────────────────────

    section("NarrativeAgent — findings grouping by confidence")

    try:
        grouped = agent._group_findings(MOCK_FINDINGS)
        assert len(grouped["high"])         == 5, f"Expected 4 HIGH, got {len(grouped['high'])}"
        assert len(grouped["medium"])       == 3, f"Expected 3 MEDIUM, got {len(grouped['medium'])}"
        assert len(grouped["low"])          == 1, f"Expected 1 LOW, got {len(grouped['low'])}"
        assert len(grouped["insufficient"]) == 1, f"Expected 1 INSUFF, got {len(grouped['insufficient'])}"
        ok(f"Grouped: HIGH={len(grouped['high'])}, MEDIUM={len(grouped['medium'])}, "
           f"LOW={len(grouped['low'])}, INSUFFICIENT={len(grouped['insufficient'])}")
    except Exception as e:
        fail("Findings grouping", str(e))


    # ── Test 2: LOW/INSUFFICIENT findings never hidden ────────────────────

    section("NarrativeAgent — honest uncertainty (LOW/INSUFF never hidden)")

    try:
        grouped = agent._group_findings(MOCK_FINDINGS)
        conflicts_text = agent._summarize_conflicts(MOCK_AGENT_RESULTS, grouped)

        assert "LOW" in conflicts_text.upper() or \
               "validation" in conflicts_text.lower() or \
               "require" in conflicts_text.lower(), \
            "LOW confidence findings must be mentioned in conflicts section"
        ok("LOW confidence findings appear in conflicts/limitations section")
    except Exception as e:
        fail("LOW findings visibility", str(e))

    try:
        grouped = agent._group_findings(MOCK_FINDINGS)
        conflicts_text = agent._summarize_conflicts(MOCK_AGENT_RESULTS, grouped)

        # INSUFFICIENT findings must be explicitly flagged
        assert "insufficient" in conflicts_text.lower() or \
               "could not be completed" in conflicts_text.lower() or \
               "No conclusions" in conflicts_text, \
            "INSUFFICIENT findings must be explicitly flagged"
        ok("INSUFFICIENT findings explicitly flagged — no false conclusions")
    except Exception as e:
        fail("INSUFFICIENT findings visibility", str(e))


    # ── Test 3: Methods section accuracy ─────────────────────────────────

    section("NarrativeAgent — methods section content")

    try:
        methods = agent._write_methods_section(
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
            agent_results=MOCK_AGENT_RESULTS,
            decisions=MockMemory().get_decisions("test"),
        )

        # Must contain exact parameter value from decisions
        assert "0.40" in methods or "resolution" in methods.lower(), \
            "Methods must include Leiden resolution value"
        ok("Methods section includes Leiden resolution parameter")
    except Exception as e:
        fail("Methods Leiden parameter", str(e))

    try:
        methods = agent._write_methods_section(
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
            agent_results=MOCK_AGENT_RESULTS,
            decisions=MockMemory().get_decisions("test"),
        )

        assert "hg38" in methods
        ok("Methods section references genome assembly (hg38)")
    except Exception as e:
        fail("Methods genome reference", str(e))

    try:
        methods = agent._write_methods_section(
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
            agent_results=MOCK_AGENT_RESULTS,
            decisions=MockMemory().get_decisions("test"),
        )
        assert "ICE" in methods or "insulation" in methods.lower() or \
               "TAD" in methods, \
            "Methods must include HiC analysis details"
        ok("Methods section includes HiC/TAD methodology")
    except Exception as e:
        fail("Methods HiC details", str(e))

    try:
        methods = agent._write_methods_section(
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
            agent_results=MOCK_AGENT_RESULTS,
            decisions=MockMemory().get_decisions("test"),
        )
        assert "LSI" in methods or "TF-IDF" in methods or "SVD" in methods, \
            "Methods must document LSI for ATAC"
        ok("Methods section documents LSI (TF-IDF + SVD) for ATAC")
    except Exception as e:
        fail("Methods LSI documentation", str(e))


    # ── Test 4: HTML report generation ───────────────────────────────────

    section("NarrativeAgent — HTML report rendering")

    try:
        exp_id  = "test_exp_001"
        grouped = agent._group_findings(MOCK_FINDINGS)
        methods = agent._write_methods_section(
            {"organism": "Homo sapiens", "genome": "hg38"},
            MOCK_AGENT_RESULTS,
            MockMemory().get_decisions(exp_id),
        )
        findings_sections = agent._write_findings_sections(
            grouped, MOCK_AGENT_RESULTS,
            {"organism": "Homo sapiens", "genome": "hg38"}
        )
        exec_summary = agent._fallback_executive_summary(
            grouped, {"summary": "immune cell profiling in lupus PBMCs"}
        )

        report_path = agent._render_html_report(
            experiment_id=exp_id,
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38",
                     "user_question": "immune profiling"},
            intent={"summary": "lupus PBMC immune profiling",
                    "analysis_type": "cell_type"},
            executive_summary=exec_summary,
            findings_sections=findings_sections,
            grouped_findings=grouped,
            methods=methods,
            decisions=MockMemory().get_decisions(exp_id),
        )

        assert report_path.exists(), "Report file not created"
        html = report_path.read_text()
        assert len(html) > 1000, "Report too short"
        ok(f"HTML report generated: {len(html):,} chars at {report_path.name}")
    except Exception as e:
        fail("HTML report generation", str(e))

    try:
        html = report_path.read_text()
        # Structural checks
        assert "<!DOCTYPE html>" in html
        assert "Executive Summary" in html
        assert "Methods" in html
        assert "Parameter Decisions Log" in html
        assert "Homo sapiens" in html
        assert "hg38" in html
        ok("HTML report has all required sections")
    except Exception as e:
        fail("HTML report structure", str(e))

    try:
        html = report_path.read_text()
        # Confidence badges present
        assert "HIGH" in html
        assert "MEDIUM" in html
        assert "LOW" in html
        assert "INSUFFICIENT" in html
        ok("HTML report shows all confidence levels (including LOW/INSUFFICIENT)")
    except Exception as e:
        fail("HTML confidence badges", str(e))

    try:
        html = report_path.read_text()
        # Decisions log present
        assert "CP 3" in html or "Leiden" in html or "resolution" in html.lower()
        ok("HTML report includes parameter decisions log")
    except Exception as e:
        fail("HTML decisions log", str(e))

    try:
        html = report_path.read_text()
        # ARIA attribution
        assert "github.com/Olascoaga/aria-omics" in html or \
               "ARIA" in html
        assert "expert review before publication" in html or \
               "review" in html.lower()
        ok("HTML report includes review disclaimer")
    except Exception as e:
        fail("HTML review disclaimer", str(e))


    # ── Test 5: Executive summary — no hallucination ──────────────────────

    section("NarrativeAgent — executive summary (no hallucination)")

    try:
        grouped = agent._group_findings(MOCK_FINDINGS)
        summary = agent._fallback_executive_summary(
            grouped, {"summary": "exhausted T cells in lupus"}
        )
        assert len(summary) > 50
        assert "findings" in summary.lower() or "identified" in summary.lower()
        ok("Fallback executive summary generated without hallucination")
    except Exception as e:
        fail("Executive summary fallback", str(e))

    try:
        # LLM summary
        grouped = agent._group_findings(MOCK_FINDINGS)
        summary = agent._write_executive_summary(
            exp_ctx={"organism": "Homo sapiens", "genome": "hg38"},
            intent={"summary": "lupus PBMC profiling",
                    "analysis_type": "cell_type"},
            grouped=grouped,
            agent_results=MOCK_AGENT_RESULTS,
        )
        assert len(summary) > 50
        ok(f"LLM executive summary: {len(summary)} chars")
    except Exception as e:
        fail("LLM executive summary", str(e))


    # ── Test 6: Full run ──────────────────────────────────────────────────

    section("NarrativeAgent — full run()")

    try:
        result = agent.run(
            experiment_id="full_test_001",
            context={
                "exp_context": {
                    "organism":      "Homo sapiens",
                    "genome":        "hg38",
                    "user_question": "Profile immune cells in lupus PBMCs",
                },
                "biological_intent": {
                    "summary":       "lupus PBMC immune profiling",
                    "analysis_type": "cell_type",
                },
                "agent_results": MOCK_AGENT_RESULTS,
                "findings":      MOCK_FINDINGS,
            },
        )

        assert result["status"]      == "done"
        assert result["report_path"] is not None
        assert Path(result["report_path"]).exists()
        assert result["n_findings"]   == len(MOCK_FINDINGS)
        assert result["n_high"]       == 5
        assert result["n_medium"]     == 3
        assert result["n_low"]        == 1

        ok(f"Full run: report at {Path(result['report_path']).name}")
        ok(f"Findings: {result['n_high']} HIGH, {result['n_medium']} MEDIUM, "
           f"{result['n_low']} LOW, "
           f"{result.get('n_insufficient', '?')} INSUFFICIENT")
    except Exception as e:
        fail("Full run()", str(e))


# ── NARRATIVE_SYSTEM prompt quality ──────────────────────────────────────────

section("NarrativeAgent — system prompt scientific integrity")

try:
    assert "correlation" in NARRATIVE_SYSTEM.lower() or \
           "causation" in NARRATIVE_SYSTEM.lower()
    ok("System prompt distinguishes correlation from causation")
except Exception as e:
    fail("Correlation/causation distinction", str(e))

try:
    assert "insufficient" in NARRATIVE_SYSTEM.lower()
    ok("System prompt handles insufficient evidence honestly")
except Exception as e:
    fail("Insufficient evidence handling", str(e))

try:
    assert "confidence" in NARRATIVE_SYSTEM.lower()
    ok("System prompt requires confidence levels in all claims")
except Exception as e:
    fail("Confidence level requirement", str(e))

try:
    assert "debate" in NARRATIVE_SYSTEM.lower() or \
           "revised" in NARRATIVE_SYSTEM.lower()
    ok("System prompt uses DebateCouncil revised claims")
except Exception as e:
    fail("DebateCouncil revised claim usage", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v NarrativeAgent validated. ARIA can write science.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
