"""
ARIA DebateCouncil Tests
-------------------------
Validates the peer review system without requiring live API calls.
Tests structure, parsing, verdict logic, and mock debate flow.

Run:
  conda activate aria-env
  python tests/test_debate_council.py
"""

from __future__ import annotations
import sys, os
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


# ── Test 1: Import and structure ──────────────────────────────────────────────

section("DebateCouncil — import and structure")

try:
    from aria.agents.debate_council import (
        DebateCouncil, DebateVerdict, DebateResult, DebateRound,
        PROPOSER_SYSTEM, CRITIC_SYSTEM, CONSENSUS_SYSTEM
    )
    ok("DebateCouncil imported successfully")
except Exception as e:
    fail("Import failed", str(e))

try:
    assert "ALTERNATIVE HYPOTHESIS FIRST" in CRITIC_SYSTEM
    assert "cross-validation" in CRITIC_SYSTEM.lower() or "CROSS-VALIDATION" in CRITIC_SYSTEM
    assert "ACCEPT" in CRITIC_SYSTEM
    assert "REJECT" in CRITIC_SYSTEM
    assert "TOPOLOGICAL" in CRITIC_SYSTEM
    ok("Critic system prompt contains all required adversarial directives")
except Exception as e:
    fail("Critic system prompt missing directives", str(e))

try:
    assert "alternative hypothesis" in PROPOSER_SYSTEM.lower() or \
           "KNOWN_LIMITATIONS" in PROPOSER_SYSTEM
    ok("Proposer system prompt includes limitations requirement")
except Exception as e:
    fail("Proposer system prompt issue", str(e))


# ── Test 2: Verdict parsing ───────────────────────────────────────────────────

section("DebateCouncil — critic response parsing")

class MockLLM:
    """Mock LLM that returns deterministic responses for testing."""
    def complete(self, prompt, system="", tier=None, max_tokens=1024, messages=None):
        if "CRITIC_SYSTEM" in system or "alternative hypothesis" in system.lower():
            return self._critic_response()
        elif "NarrativeAgent" in system or "consensus" in system.lower():
            return self._consensus_response()
        else:
            return self._proposer_response()

    def _proposer_response(self):
        return (
            "CLAIM: Cluster 3 represents terminally exhausted CD8+ T cells\n"
            "EVIDENCE: PDCD1 log2FC=2.1, TOX log2FC=1.8, TCF7 not significant\n"
            "KNOWN_LIMITATIONS: Small cluster size (n=87); no protein validation"
        )

    def _critic_response(self):
        return (
            "ALTERNATIVE_HYPOTHESIS: Cluster 3 could be precursor-exhausted (Tpex) "
            "not terminal Tex — TCF7 must be explicitly negative, not just non-significant\n"
            "CHALLENGES: Absence of GZMB expression not confirmed; "
            "co-expression of TOX and TCF7 not evaluated\n"
            "EVIDENCE_REQUESTED: Explicit TCF7 log2FC value; GZMB expression level\n"
            "VERDICT: ACCEPT_REVISED\n"
            "REVISED_CLAIM: Cluster 3 shows markers consistent with exhausted CD8+ T cells "
            "(PDCD1+, TOX+); terminal vs precursor status requires TCF7 confirmation"
        )

    def _consensus_response(self):
        return (
            "Cluster 3 displays a transcriptional signature consistent with exhausted "
            "CD8+ T cells, marked by elevated PDCD1 and TOX expression. "
            "Terminal exhaustion (Tex) versus precursor-exhausted (Tpex) status "
            "could not be definitively established without explicit TCF7 quantification.\n\n"
            "Limitations:\n"
            "- TCF7 status requires explicit quantification for Tex/Tpex distinction\n"
            "- Small cluster size (n=87) limits statistical power\n"
            "- Protein-level validation (flow cytometry) recommended before publishing"
        )

mock_llm = MockLLM()

try:
    council = DebateCouncil(llm=mock_llm, max_rounds=3)
    ok("DebateCouncil instantiated with mock LLM")
except Exception as e:
    fail("Instantiation failed", str(e))

try:
    # Test critic parsing
    critic_resp = mock_llm._critic_response()
    verdict, alternative, challenges, revised = council._parse_critic(critic_resp)

    assert verdict == DebateVerdict.ACCEPT_REVISED, f"Expected ACCEPT_REVISED, got {verdict}"
    assert "Tpex" in alternative or "precursor" in alternative.lower()
    assert len(challenges) > 0
    assert "TCF7" in revised or "exhausted" in revised.lower()

    ok(f"Critic parsing: verdict={verdict.value}",
       f"alternative='{alternative[:60]}...'")
except Exception as e:
    fail("Critic parsing failed", str(e))


# ── Test 3: Full debate flow ──────────────────────────────────────────────────

section("DebateCouncil — full debate simulation")

try:
    result = council.resolve(
        topic="cell_type_annotation",
        initial_claim="Cluster 3 represents terminally exhausted CD8+ T cells",
        evidence={
            "top_markers":  ["PDCD1", "TOX", "HAVCR2", "LAG3"],
            "log2fc":       {"PDCD1": 2.1, "TOX": 1.8, "TCF7": 0.3},
            "n_cells":      87,
            "cluster_id":   "3",
        },
        biological_context={
            "organism":      "Homo sapiens",
            "user_question": "What immune cell subtypes are present in lupus PBMCs?",
            "analysis_type": "cell_type",
        },
    )

    assert isinstance(result, DebateResult)
    assert result.verdict in list(DebateVerdict)
    assert result.consensus
    assert result.confidence in ("high", "medium", "low", "insufficient")
    assert isinstance(result.rounds, list)
    assert len(result.rounds) >= 1

    ok(f"Full debate completed: verdict={result.verdict.value}, "
       f"confidence={result.confidence}, rounds={len(result.rounds)}")
    ok(f"Consensus generated ({len(result.consensus)} chars)")
    ok(f"Was revised by Critic: {result.was_revised}")

    if result.limitations:
        ok(f"Limitations extracted: {len(result.limitations)}")
        for lim in result.limitations[:2]:
            print(f"      {DIM}* {lim[:80]}{RST}")

except Exception as e:
    fail("Full debate simulation failed", str(e))


# ── Test 4: Verdict confidence mapping ───────────────────────────────────────

section("DebateCouncil — verdict to confidence mapping")

try:
    assert council._verdict_to_confidence(DebateVerdict.ACCEPT, 1)         == "high"
    assert council._verdict_to_confidence(DebateVerdict.ACCEPT, 3)         == "medium"
    assert council._verdict_to_confidence(DebateVerdict.ACCEPT_REVISED, 2) == "medium"
    assert council._verdict_to_confidence(DebateVerdict.REJECT, 2)         == "low"
    assert council._verdict_to_confidence(DebateVerdict.INSUFFICIENT, 3)   == "insufficient"
    ok("Confidence mapping: ACCEPT(1)=high, ACCEPT(3)=medium, "
       "REJECT=low, INSUFFICIENT=insufficient")
except Exception as e:
    fail("Confidence mapping failed", str(e))


# ── Test 5: Integration with Confidence enum from MessageBus ─────────────────

section("DebateCouncil — integration with MessageBus Confidence")

try:
    from aria.bus.message_bus import Confidence

    mapping = {
        "high":         Confidence.HIGH,
        "medium":       Confidence.MEDIUM,
        "low":          Confidence.LOW,
        "insufficient": Confidence.INSUFFICIENT,
    }
    for str_conf, enum_conf in mapping.items():
        assert enum_conf.value == str_conf
    ok("DebateCouncil confidence strings map correctly to MessageBus Confidence enum")
except Exception as e:
    fail("MessageBus Confidence mapping failed", str(e))


# ── Test 6: Adversarial directives validation ─────────────────────────────────

section("DebateCouncil — adversarial directives from Gemini's review")

try:
    # Verify the three key directives Gemini specified are present
    assert "alternative hypothesis" in CRITIC_SYSTEM.lower(), \
        "Missing: force alternative hypothesis"
    assert "cross-validation" in CRITIC_SYSTEM.lower() or \
           "absence" in CRITIC_SYSTEM.lower(), \
        "Missing: demand cross-validation / negative evidence"
    assert "topolog" in CRITIC_SYSTEM.lower(), \
        "Missing: topological awareness for algorithmic fit"
    ok("All 3 adversarial directives from Gemini's review present in Critic prompt")
    ok("Directive 1: Alternative hypothesis formulation required")
    ok("Directive 2: Cross-validation / negative evidence required")
    ok("Directive 3: Topological awareness for algorithmic fit")
except AssertionError as e:
    fail("Missing adversarial directive", str(e))
except Exception as e:
    fail("Directive validation failed", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}v DebateCouncil validated. ARIA has peer review.{RST}\n")
else:
    print(f"\n{YLW}Some tests need attention.{RST}\n")
    sys.exit(1)
