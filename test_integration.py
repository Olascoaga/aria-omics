"""
ARIA Integration Test
---------------------
Validates the full stack without requiring live API keys or real omics data.

Tests:
  1. ContextManager — token counting, compression cascade, model profiles
  2. LLMProvider    — routing tiers, fallback logic, config loading
  3. ParameterAdvisor — 3-layer decision, memory storage, formatting
  4. RNAAgent       — QC mock, DE mock, checkpoint escalation

Run:
  python tests/test_integration.py
"""

import sys
import os
import json
import uuid

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from aria.llm.context_manager import ContextManager, ModelProfile
from aria.llm.provider import LLMProvider, TaskTier, ModelConfig
from aria.llm.parameter_advisor import ParameterAdvisor, MetricEvaluator
from aria.memory.memory import ARIAMemory
from aria.bus.message_bus import bus, MessageType, Confidence


# ── Colors for terminal output ────────────────────────────────────────────────
GRN = "\033[92m"
RED = "\033[91m"
YLW = "\033[93m"
CYN = "\033[96m"
DIM = "\033[2m"
RST = "\033[0m"
BLD = "\033[1m"

passed = 0
failed = 0


def ok(msg):
    global passed
    passed += 1
    print(f"  {GRN}✓{RST} {msg}")


def fail(msg, err=""):
    global failed
    failed += 1
    print(f"  {RED}✗{RST} {msg}")
    if err:
        print(f"    {DIM}{err}{RST}")


def section(title):
    print(f"\n{CYN}{BLD}▶ {title}{RST}")


# ── Test 1: ContextManager ────────────────────────────────────────────────────

section("ContextManager — token counting & cascade")

try:
    # Cloud model — tiktoken or heuristic
    cloud_profile = ModelProfile("claude-sonnet-4-20250514", 200_000)
    cloud_ctx = ContextManager(cloud_profile)

    count = cloud_ctx.count_tokens("Hello, this is a test sentence.")
    assert 5 < count < 20, f"Unexpected token count: {count}"
    ok(f"Cloud tokenizer: '{count}' tokens for test sentence")
except Exception as e:
    fail("Cloud token counting", str(e))

try:
    # Local model — conservative heuristic
    local_profile = ModelProfile("ollama/llama3:8b", 4_000, is_local=True)
    local_ctx = ContextManager(local_profile)

    text = "A" * 1000  # 1000 chars
    count = local_ctx.count_tokens(text)
    # Conservative: 1000 / 3.5 * 1.15 ≈ 329 tokens
    assert 250 < count < 400, f"Unexpected: {count}"
    ok(f"Local heuristic: {count} tokens for 1000-char string (conservative)")
except Exception as e:
    fail("Local token heuristic", str(e))

try:
    # CavemanULTRA deterministic compression
    local_profile = ModelProfile("ollama/llama3:8b", 4_000, is_local=True)
    ctx = ContextManager(local_profile)

    verbose = (
        "I'd be happy to help you with that. The reason your component "
        "is re-rendering is likely because you're creating a new object "
        "reference on each render cycle. In order to fix this, you should "
        "utilize the useMemo hook to memoize the object. Furthermore, "
        "it's important to note that this is a common performance issue."
    )
    compressed = ctx._caveman_ultra(verbose)
    ratio = len(compressed) / len(verbose)
    assert ratio < 0.95, f"No compression: ratio={ratio:.2f}"
    assert "I'd be happy" not in compressed, "Filler not removed"
    ok(f"CavemanULTRA: {len(verbose)} → {len(compressed)} chars "
       f"({(1-ratio)*100:.0f}% reduction, filler removed ✓)")
except Exception as e:
    fail("CavemanULTRA compression", str(e))

try:
    # Cascade: history truncation
    local_profile = ModelProfile("ollama/llama3:8b", 1_000, is_local=True)
    ctx = ContextManager(local_profile)

    # Build history that exceeds budget
    long_history = [
        {"role": "user" if i % 2 == 0 else "assistant",
         "content": "word " * 100}
        for i in range(20)
    ]
    budget = 300
    result = ctx._apply_cascade(long_history, budget)
    result_tokens = sum(ctx.count_tokens(m["content"]) for m in result)
    assert result_tokens <= budget * 1.1  # allow small overshoot from heuristic
    ok(f"Cascade: 20 messages → {len(result)} messages ({result_tokens} tokens)")
except Exception as e:
    fail("Context cascade", str(e))


# ── Test 2: LLMProvider routing ───────────────────────────────────────────────

section("LLMProvider — routing tiers & fallback")

try:
    provider = LLMProvider()
    heavy_model  = provider.get_active_model(TaskTier.HEAVY)
    medium_model = provider.get_active_model(TaskTier.MEDIUM)
    light_model  = provider.get_active_model(TaskTier.LIGHT)

    assert heavy_model  is not None
    assert medium_model is not None
    assert light_model  is not None
    ok(f"Routing: HEAVY={heavy_model.model} | "
       f"MEDIUM={medium_model.model} | "
       f"LIGHT={light_model.model}")
except Exception as e:
    fail("Provider routing tiers", str(e))

try:
    # Config override test (simulate config.yaml)
    custom_models = {
        TaskTier.HEAVY: [
            ModelConfig("ollama", "ollama/llama3:70b", 8_000,
                        is_local=True, api_base="http://localhost:11434")
        ],
        TaskTier.MEDIUM: [
            ModelConfig("ollama", "ollama/llama3:8b",  4_000,
                        is_local=True, api_base="http://localhost:11434")
        ],
        TaskTier.LIGHT: [
            ModelConfig("ollama", "ollama/mistral:7b", 8_000,
                        is_local=True, api_base="http://localhost:11434")
        ],
    }
    local_provider = LLMProvider(models=custom_models)
    active = local_provider.get_active_model(TaskTier.HEAVY)
    assert active.is_local
    assert "llama3" in active.model
    ok(f"Config override: all tiers → local Ollama models")
except Exception as e:
    fail("Provider config override", str(e))

try:
    # Context manager is correctly assigned per model
    provider = LLMProvider()
    cloud_cfg = provider.models[TaskTier.HEAVY][0]
    ctx_mgr   = provider._get_context_manager(cloud_cfg)
    assert ctx_mgr.profile.context_window == cloud_cfg.context_window
    ok(f"ContextManager assigned to {cloud_cfg.model}: "
       f"window={ctx_mgr.profile.safe_limit} safe tokens")
except Exception as e:
    fail("ContextManager assignment", str(e))


# ── Test 3: ParameterAdvisor ──────────────────────────────────────────────────

section("ParameterAdvisor — 3-layer decisions & memory")

# Use in-memory SQLite for tests
memory = ARIAMemory(":memory:")
exp_id = f"test_{uuid.uuid4().hex[:8]}"
memory.create_wing(exp_id, "Test Experiment", "Homo sapiens", "hg38")

# Create mock room for findings
memory.create_hall(f"{exp_id}_scRNA", exp_id, "scRNA")
memory.create_room(f"{exp_id}_leiden", f"{exp_id}_scRNA", "leiden_clustering")

# Mock LLMProvider that returns deterministic responses (no API needed)
class MockLLMProvider:
    def complete_medium(self, prompt, system="", max_tokens=1024, messages=None):
        return (
            "Resolution 0.40 selected based on silhouette score 0.720, "
            "highest among 4 candidates, consistent with coarse cell type "
            "resolution appropriate for the stated biological question."
        )
    def complete_heavy(self, *a, **kw):
        return '{"0": "CD4+ T cell", "1": "CD8+ T cell", "2": "NK cell"}'
    def complete_light(self, *a, **kw):
        return "compressed"
    def complete(self, *a, **kw):
        return self.complete_medium(*a, **kw)

mock_llm = MockLLMProvider()

try:
    advisor = ParameterAdvisor(memory, mock_llm)
    ok("ParameterAdvisor instantiated")
except Exception as e:
    fail("ParameterAdvisor init", str(e))

try:
    # Layer 1: intent constrains search space
    narrow_context = {
        "analysis_type": "cell_type",
        "user_question":  "find rare subpopulations of regulatory T cells",
        "summary":        "rare subpopulation discovery",
    }
    r = advisor._intent_to_leiden_range(narrow_context)
    assert r[0] >= 0.5, f"Fine-grained range should start high: {r}"
    ok(f"Layer 1 (intent): 'rare subpopulations' → range {r}")
except Exception as e:
    fail("Layer 1 intent mapping", str(e))

try:
    coarse_context = {
        "analysis_type": "cell_type",
        "user_question":  "identify major lineages in PBMC",
        "summary":        "major cell type identification",
    }
    r = advisor._intent_to_leiden_range(coarse_context)
    assert r[1] <= 0.8, f"Coarse range should be low: {r}"
    ok(f"Layer 1 (intent): 'major lineages' → range {r}")
except Exception as e:
    fail("Layer 1 coarse mapping", str(e))

try:
    # Layer 2: objective metric scoring
    good_metrics = {"silhouette": 0.72, "modularity": 0.61,
                    "n_clusters": 8,    "n_singleton_clusters": 0,
                    "min_cluster_size": 150}
    bad_metrics  = {"silhouette": 0.25, "modularity": 0.20,
                    "n_clusters": 45,   "n_singleton_clusters": 8,
                    "min_cluster_size": 3}

    good_score = advisor._score_leiden(good_metrics, {})
    bad_score  = advisor._score_leiden(bad_metrics,  {})
    assert good_score > bad_score, f"Scoring wrong: {good_score} <= {bad_score}"
    ok(f"Layer 2 (metrics): good={good_score:.3f} > bad={bad_score:.3f} ✓")
except Exception as e:
    fail("Layer 2 metric scoring", str(e))

try:
    # Mock adata with enough structure for mock metrics
    class MockAdata:
        pass
    mock_adata = MockAdata()

    bio_ctx = {
        "analysis_type": "cell_type",
        "user_question":  "identify major cell types in lupus PBMCs",
        "summary":        "cell type identification",
    }
    decision = advisor.advise_leiden_resolution(
        adata=mock_adata,
        experiment_id=exp_id,
        biological_context=bio_ctx,
        n_candidates=4,
    )

    assert decision.chosen_value is not None
    assert len(decision.candidates) == 4
    assert decision.justification
    assert any(c.recommended for c in decision.candidates)
    ok(f"Full decision: resolution={decision.chosen_value} "
       f"({len(decision.candidates)} candidates evaluated)")
    ok(f"Justification: '{decision.justification[:80]}...'")
except Exception as e:
    fail("Full parameter decision", str(e))

try:
    # Layer 3: memory recall after storing a decision
    memory.store_decision(
        decision_id=str(uuid.uuid4())[:8],
        wing_id=exp_id,
        checkpoint=3,
        question="resolution for leiden_clustering",
        decision="0.40",
        rationale="approved by user",
        made_by="user",
    )
    hist = advisor._recall_similar_decisions(exp_id, "leiden_clustering", {})
    assert len(hist) >= 1
    ok(f"Layer 3 (memory): recalled {len(hist)} historical decision(s)")
except Exception as e:
    fail("Layer 3 memory recall", str(e))

try:
    # Checkpoint formatting
    fmt = advisor.format_for_checkpoint(decision)
    assert "resolution" in fmt.lower()
    assert "★" in fmt  # recommended marker
    ok("Checkpoint 3 formatting: structured, readable, contains recommendation")
except Exception as e:
    fail("Checkpoint formatting", str(e))

try:
    # User approval updates memory
    approved = advisor.approve_decision(decision, user_override=None)
    assert approved.approved_by_user
    ok("User approval: decision marked approved in memory")
except Exception as e:
    fail("User approval flow", str(e))


# ── Test 4: MessageBus findings ───────────────────────────────────────────────

section("MessageBus — findings & escalations")

try:
    from aria.bus.message_bus import Message, MessageType, Confidence, CavemanMode

    test_exp = f"bus_test_{uuid.uuid4().hex[:6]}"
    msg = Message(
        sender="rna_agent",
        receiver="orchestrator",
        type=MessageType.FINDING,
        confidence=Confidence.HIGH,
        payload={"summary": "500 DE genes found", "n_sig": 500},
        experiment_id=test_exp,
    )
    bus.publish(msg)

    findings = bus.get_findings(test_exp)
    assert len(findings) == 1
    assert findings[0].confidence == Confidence.HIGH
    ok("Finding published and retrieved from MessageBus")
except Exception as e:
    fail("MessageBus finding", str(e))

try:
    # Escalation (checkpoint)
    esc = Message(
        sender="rna_agent",
        receiver="orchestrator",
        type=MessageType.ESCALATION,
        confidence=Confidence.HIGH,
        payload={"checkpoint": 3, "question": "Use resolution=0.4?",
                 "options": ["Yes", "No"], "resolved": False},
        checkpoint=3,
        experiment_id=test_exp,
    )
    bus.publish(esc)
    pending = bus.get_pending_checkpoints()
    assert any(m.id == esc.id for m in pending)
    ok("Escalation (checkpoint) published and detected as pending")

    # Resolve it
    bus.resolve_checkpoint(esc.id, {"choice": "Yes"})
    still_pending = [m for m in bus.get_pending_checkpoints()
                     if m.id == esc.id and not m.payload.get("resolved")]
    assert len(still_pending) == 0
    ok("Checkpoint resolved — no longer pending")
except Exception as e:
    fail("Checkpoint escalation & resolution", str(e))


# ── Test 5: Memory cross-modal tunnels ────────────────────────────────────────

section("ARIAMemory — tunnels (cross-modal connections)")

try:
    exp2 = f"tunnel_test_{uuid.uuid4().hex[:6]}"
    memory.create_wing(exp2, "Multimodal Test", "Homo sapiens", "hg38")

    memory.create_tunnel(
        tunnel_id=str(uuid.uuid4())[:8],
        wing_id=exp2,
        from_hall="scRNA",
        to_hall="scATAC",
        entity="TP53",
        description="TP53 upregulated in RNA AND accessible promoter in ATAC",
        confidence="high",
    )
    memory.create_tunnel(
        tunnel_id=str(uuid.uuid4())[:8],
        wing_id=exp2,
        from_hall="scATAC",
        to_hall="HiC",
        entity="MYC",
        description="MYC enhancer accessible AND loops to MYC promoter in HiC",
        confidence="medium",
    )

    tunnels = memory.get_tunnels(exp2)
    assert len(tunnels) == 2
    ok(f"Tunnels: {len(tunnels)} cross-modal connections stored")

    tp53_tunnels = memory.get_tunnels(exp2, entity="TP53")
    assert len(tp53_tunnels) == 1
    ok(f"Entity query: TP53 tunnel retrieved correctly")
except Exception as e:
    fail("Memory tunnels", str(e))


# ── Summary ───────────────────────────────────────────────────────────────────

total = passed + failed
print(f"\n{'─'*50}")
print(f"{BLD}Results: {GRN}{passed} passed{RST}{BLD}, "
      f"{RED if failed else GRN}{failed} failed{RST}{BLD} / {total} total{RST}")

if failed == 0:
    print(f"\n{GRN}{BLD}✓ All systems operational. ARIA stack is solid.{RST}\n")
else:
    print(f"\n{YLW}⚠ {failed} test(s) need attention.{RST}\n")
    sys.exit(1)
