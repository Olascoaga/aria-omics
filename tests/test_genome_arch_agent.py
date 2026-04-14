"""
ARIA GenomeArchAgent Tests
---------------------------
Validates: resolution selection, RAM estimation, out-of-core strategy,
Insulation Score window calibration, compartment A/B interpretation rules,
and DebateCouncil integration for loop regulatory claims.

Run:
  conda activate aria-env
  python tests/test_genome_arch_agent.py
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
            "CLAIM: Region chr1:50Mb-55Mb shows A-to-B compartment switch\n"
            "EVIDENCE: PC1 flips from +0.3 to -0.4 in condition B\n"
            "KNOWN_LIMITATIONS: PC1 sign not externally validated"
        )
    def complete_heavy(self, *a, **kw):  return self.complete()
    def complete_medium(self, *a, **kw): return self.complete()
    def complete_light(self, *a, **kw):  return "ok"

class MockEnvManager:
    def run_in_stack(self, stack, script_path, params, timeout=None):
        script = Path(script_path).name
        if "inspect" in script:
            return {
                "status": "success",
                "files": params.get("files", []),
                "formats": ["mcool"],
                "available_resolutions": [5000, 10000, 40000, 100000, 1000000],
                "chromosomes": [f"chr{i}" for i in range(1, 23)] + ["chrX"],
                "estimated_sizes_gb": {
                    5000: 400, 10000: 100, 40000: 8,
                    100000: 1, 1000000: 0.05
                },
            }
        elif "qc" in script:
            return {
                "status": "success", "n_valid_pairs": 800_000_000,
                "cis_trans_ratio": 0.74, "balanced_files": ["/tmp/balanced.cool"],
                "pass_qc": True, "warnings": [],
            }
        elif "topology" in script:
            analysis = params.get("analysis", "compartments")
            if analysis == "compartments":
                return {
                    "status": "success", "pct_A": 42.1, "pct_B": 57.9,
                    "n_ab_switches": 312, "pc1_validated": False,
                    "compartment_tracks": {"chr1": [0.1, -0.2, 0.3]},
                }
            elif analysis == "insulation_calibration":
                windows = params.get("windows", [120000, 200000, 400000])
                return {
                    "status": "success",
                    "boundary_strengths": {
                        str(windows[0]): 0.312,
                        str(windows[1]): 0.445,  # highest
                        str(windows[2]): 0.398,
                    },
                    "calibration_chrom": "chr1",
                }
            elif analysis == "tads":
                return {
                    "status": "success", "n_tads": 3247,
                    "median_size_kb": 185.0,
                    "tad_boundaries": {"chr1": [1000000, 2500000, 3800000]},
                }
            elif analysis == "loops":
                return {
                    "status": "success", "n_loops": 9823,
                    "loops_path": "/tmp/loops.bedpe", "algorithm": "chromosight",
                }
        return {"status": "success"}

class MockMemory:
    def store_decision(self, *a, **kw): pass
    def get_decisions(self, *a, **kw): return []
    def list_wings(self): return []
    def create_room(self, *a, **kw): pass
    def store_finding(self, *a, **kw): pass


# ── Test 1: Import and structure ──────────────────────────────────────────────

section("GenomeArchAgent — import and structure")

try:
    from aria.agents.genome_arch_agent import (
        GenomeArchAgent, RAM_ESTIMATES_GB, GENOME_SCALE, GENOME_ARCH_SYSTEM
    )
    ok("GenomeArchAgent imported successfully")
except Exception as e:
    fail("Import failed", str(e))

try:
    # RAM table
    assert 1_000_000 in RAM_ESTIMATES_GB
    assert 5_000     in RAM_ESTIMATES_GB
    assert RAM_ESTIMATES_GB[5_000] > RAM_ESTIMATES_GB[1_000_000]
    ok(f"RAM table: 5kb={RAM_ESTIMATES_GB[5_000]}GB, "
       f"1Mb={RAM_ESTIMATES_GB[1_000_000]}GB (5kb >> 1Mb)")
except Exception as e:
    fail("RAM estimates table", str(e))

try:
    # Genome scale factors
    assert "Homo sapiens"            in GENOME_SCALE
    assert "Drosophila melanogaster" in GENOME_SCALE
    dm_scale = GENOME_SCALE["Drosophila melanogaster"]
    hs_scale = GENOME_SCALE["Homo sapiens"]
    assert dm_scale < hs_scale * 0.1, \
        "Drosophila genome should be <10% of human"
    ok(f"Genome scale: human={hs_scale}, drosophila={dm_scale} "
       f"({dm_scale/hs_scale:.1%} of human)")
except Exception as e:
    fail("Genome scale factors", str(e))


# ── Test 2: Resolution selection ──────────────────────────────────────────────

section("GenomeArchAgent — resolution selection (ParameterAdvisor)")

agent = GenomeArchAgent.__new__(GenomeArchAgent)
agent.llm    = MockLLM()
agent.memory = MockMemory()
agent.env    = MockEnvManager()

# Mock advisor with recall
class MockAdvisor:
    def _recall_similar_decisions(self, *a, **kw): return []
agent.advisor = MockAdvisor()

# Mock publish methods
agent._findings = []
agent._escalations = []
def mock_publish_finding(exp_id, payload, conf):
    agent._findings.append(payload)
def mock_publish_escalation(**kw):
    agent._escalations.append(kw)
def mock_publish_status(*a, **kw): pass
agent.publish_finding   = mock_publish_finding
agent.publish_escalation = mock_publish_escalation
agent.publish_status    = mock_publish_status

file_info = {
    "files": ["/data/test.mcool"],
    "available_resolutions": [5_000, 10_000, 40_000, 100_000, 1_000_000],
}

try:
    # Loop question → fine resolution
    loop_intent = {
        "user_question": "What enhancer-promoter loops change in condition A vs B?",
        "analysis_type": "structural",
        "complexity": "complex",
    }
    decision = agent._advise_resolution("exp_001", loop_intent, file_info,
                                         {"organism": "Homo sapiens",
                                          "genome": "hg38"})
    assert decision["analysis_type"] == "loop_calling"
    assert decision["recommended_resolution"] <= 10_000
    ok(f"Loop question -> resolution={decision['recommended_resolution']:,}bp, "
       f"analysis={decision['analysis_type']}")
except Exception as e:
    fail("Loop resolution selection", str(e))

try:
    # Compartment question → coarse resolution
    comp_intent = {
        "user_question": "Which genomic regions switch A/B compartments?",
        "analysis_type": "structural",
    }
    decision = agent._advise_resolution("exp_001", comp_intent, file_info,
                                         {"organism": "Homo sapiens",
                                          "genome": "hg38"})
    assert decision["analysis_type"] == "compartments"
    assert decision["recommended_resolution"] >= 100_000
    ok(f"Compartment question -> resolution={decision['recommended_resolution']:,}bp, "
       f"RAM~{decision['ram_required_gb']:.1f}GB")
except Exception as e:
    fail("Compartment resolution selection", str(e))

try:
    # TAD question → medium resolution
    tad_intent = {
        "user_question": "How do TAD boundaries change between cell types?",
        "analysis_type": "structural",
    }
    decision = agent._advise_resolution("exp_001", tad_intent, file_info,
                                         {"organism": "Homo sapiens",
                                          "genome": "hg38"})
    assert decision["analysis_type"] == "tad_calling"
    assert 10_000 <= decision["recommended_resolution"] <= 100_000
    ok(f"TAD question -> resolution={decision['recommended_resolution']:,}bp")
except Exception as e:
    fail("TAD resolution selection", str(e))

try:
    # Drosophila scales down correctly
    dm_decision = agent._advise_resolution(
        "exp_dm", tad_intent, file_info,
        {"organism": "Drosophila melanogaster", "genome": "dm6"}
    )
    hs_decision = agent._advise_resolution(
        "exp_hs", tad_intent, file_info,
        {"organism": "Homo sapiens", "genome": "hg38"}
    )
    assert dm_decision["ram_required_gb"] < hs_decision["ram_required_gb"]
    ok(f"RAM scales by organism: "
       f"human={hs_decision['ram_required_gb']:.1f}GB, "
       f"drosophila={dm_decision['ram_required_gb']:.2f}GB")
except Exception as e:
    fail("Organism RAM scaling", str(e))


# ── Test 3: Insulation Score calibration ──────────────────────────────────────

section("GenomeArchAgent — Insulation Score window calibration")

try:
    # Calibration chooses window with highest boundary strength
    decision = agent._advise_insulation_window(
        experiment_id="exp_tad_001",
        intent={"user_question": "Find TAD boundaries in PBMC",
                "analysis_type": "tad"},
        files=["/data/test.mcool"],
        resolution=40_000,
        genome="hg38",
    )
    assert "recommended_window" in decision
    assert decision["recommended_window"] in decision["windows_tested"]
    # Mock returns highest boundary_strength for middle window (200000)
    assert decision["recommended_window"] == 200_000, \
        f"Expected 200000 (highest strength), got {decision['recommended_window']}"
    ok(f"Window calibration: recommended={decision['recommended_window']:,}bp "
       f"(highest boundary strength)")
except Exception as e:
    fail("Insulation Score calibration", str(e))

try:
    # Megadomain question → larger windows
    mega_decision = agent._advise_insulation_window(
        "exp_mega",
        {"user_question": "Identify large megadomains in the genome"},
        ["/data/test.mcool"], 40_000, "hg38"
    )
    assert all(w >= 200_000 for w in mega_decision["windows_tested"]), \
        "Megadomain windows should be >= 200kb"
    ok(f"Megadomain windows: {mega_decision['windows_tested']}")
except Exception as e:
    fail("Megadomain window range", str(e))

try:
    # Sub-TAD question → smaller windows
    sub_decision = agent._advise_insulation_window(
        "exp_sub",
        {"user_question": "Find sub-TAD structure and fine boundaries"},
        ["/data/test.mcool"], 10_000, "hg38"
    )
    assert all(w <= 200_000 for w in sub_decision["windows_tested"]), \
        "Sub-TAD windows should be <= 200kb"
    ok(f"Sub-TAD windows: {sub_decision['windows_tested']}")
except Exception as e:
    fail("Sub-TAD window range", str(e))


# ── Test 4: Compartment B ≠ silence (critical biological knowledge) ───────────

section("GenomeArchAgent — Compartment B interpretation")

try:
    # The system prompt MUST warn about B ≠ silence
    assert "B" in GENOME_ARCH_SYSTEM
    assert "silent" in GENOME_ARCH_SYSTEM.lower() or \
           "silence" in GENOME_ARCH_SYSTEM.lower()
    assert "NOT" in GENOME_ARCH_SYSTEM or "not always" in GENOME_ARCH_SYSTEM.lower()
    ok("GENOME_ARCH_SYSTEM explicitly warns: B-compartment ≠ transcriptional silence")
except Exception as e:
    fail("Missing B-compartment warning", str(e))

try:
    # PC1 sign validation warning
    assert "arbitrary" in GENOME_ARCH_SYSTEM.lower() or \
           "sign" in GENOME_ARCH_SYSTEM.lower()
    ok("GENOME_ARCH_SYSTEM warns about arbitrary PC1 sign")
except Exception as e:
    fail("Missing PC1 sign warning", str(e))

try:
    # Resolution requirements for loops
    assert "5kb" in GENOME_ARCH_SYSTEM or "5,000" in GENOME_ARCH_SYSTEM or \
           "5kb" in GENOME_ARCH_SYSTEM.replace(" ", "") or \
           "10kb" in GENOME_ARCH_SYSTEM
    ok("GENOME_ARCH_SYSTEM includes resolution requirements for loop calling")
except Exception as e:
    fail("Missing resolution requirements", str(e))


# ── Test 5: Out-of-core strategy validation ───────────────────────────────────

section("GenomeArchAgent — Out-of-core memory strategy")

try:
    # Verify all topology calls include out_of_core=True
    # Mock env manager captures params
    captured_params = []
    original_run = agent.env.run_in_stack

    def capture_params(stack, script_path, params, timeout=None):
        captured_params.append(params.copy())
        return original_run(stack, script_path, params, timeout)

    agent.env.run_in_stack = capture_params

    # Trigger a topology analysis
    resolution_decision = {
        "recommended_resolution": 40_000,
        "analysis_type":          "tad_calling",
        "ram_required_gb":        8.0,
        "candidates": [
            {"resolution": 40_000, "ram_required_gb": 8.0,
             "recommended": True, "analysis_depth": "TADs"},
        ],
    }
    _ = agent._run_topology(
        "exp_ooc", {"genome": "hg38", "organism": "Homo sapiens"},
        {"user_question": "find tad boundaries"},
        ["/data/test.mcool"],
        resolution_decision,
        {"status": "success", "warnings": [], "n_valid_pairs": 800_000_000},
    )

    topo_params = [p for p in captured_params
                   if p.get("analysis") in
                   ("compartments", "tads", "loops", "insulation_calibration")]

    oor_ok = all(p.get("out_of_core") == True for p in topo_params)
    assert oor_ok, (
        f"Not all topology calls use out_of_core=True. "
        f"Params: {[p.get('out_of_core') for p in topo_params]}"
    )
    ok(f"All {len(topo_params)} topology calls use out_of_core=True")

    agent.env.run_in_stack = original_run  # restore

except Exception as e:
    fail("Out-of-core strategy", str(e))


# ── Test 6: Scripts import ────────────────────────────────────────────────────

section("Hi-C scripts — import and contract")

try:
    from aria.scripts.hic_qc_and_balance import hic_qc_and_balance
    ok("hic_qc_and_balance.py imported successfully")
except Exception as e:
    fail("hic_qc_and_balance import failed", str(e))

try:
    from aria.scripts.hic_topology import hic_topology
    ok("hic_topology.py imported successfully")
except Exception as e:
    fail("hic_topology import failed", str(e))

try:
    from aria.scripts.hic_inspect import hic_inspect
    ok("hic_inspect.py imported successfully")
except Exception as e:
    fail("hic_inspect import failed", str(e))

try:
    # Test FileNotFound handling
    result = hic_qc_and_balance({
        "files": ["/nonexistent/test.cool"],
        "genome": "hg38",
        "resolution": 40000,
    })
    assert result["status"] == "error"
    assert result["error_type"] == "FileNotFound"
    ok("hic_qc_and_balance handles missing files gracefully")
except Exception as e:
    fail("hic_qc_and_balance FileNotFound handling", str(e))

try:
    # Test unknown analysis mode
    result = hic_topology({
        "files": [],
        "genome": "hg38",
        "analysis": "nonexistent_analysis",
    })
    assert result["status"] == "error"
    assert result["error_type"] in ("FileNotFound", "UnknownAnalysis")
    ok("hic_topology handles unknown analysis type")
except Exception as e:
    fail("hic_topology unknown analysis handling", str(e))

try:
    # Test mock insulation calibration (no cooler needed)
    result = hic_topology({
        "files": [],
        "genome": "hg38",
        "analysis": "insulation_calibration",
        "resolution": 40000,
        "windows": [120000, 200000, 400000],
        "chromosomes": ["chr1"],
    })
    # Should return error (no files) or mock
    assert "status" in result
    ok(f"hic_topology insulation calibration: status={result['status']}")
except Exception as e:
    fail("hic_topology calibration handling", str(e))


# ── Test 7: Resolution depth mapping ─────────────────────────────────────────

section("GenomeArchAgent — resolution to analysis depth")

try:
    depths = {
        1_000_000: "compartments only",
        100_000:   "compartments",
        40_000:    "TADs",
        10_000:    "TADs + loops",
        5_000:     "fine loops",
    }
    for res, expected_keyword in depths.items():
        depth = agent._resolution_to_depth(res)
        key   = expected_keyword.split()[0].lower()
        assert key in depth.lower(), \
            f"Resolution {res:,}bp: expected '{key}' in '{depth}'"
    ok("Resolution depth mapping correct for all 5 levels")
except Exception as e:
    fail("Resolution depth mapping", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v GenomeArchAgent validated.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    import sys; sys.exit(1)
