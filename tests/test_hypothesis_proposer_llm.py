"""S5 guards for the LLM proposer + the full bulk RNA -> hypotheses chain.

The proposer is the single place the LLM is given freedom. It is tested with a
FAKE complete callable (no network): the real model only enters at runtime. The
end-to-end test proves audited bulk RNA findings flow through adapter -> proposer
-> the four gates -> quarantine, and that a model that invents a fact or drops a
confound is rejected, never published.
"""

from __future__ import annotations

import json

from aria.agents.hypothesis_agent import HypothesisAgent
from aria.agents.narrative.hypothesis import (
    LLMProposer,
    bulk_rna_evidence,
    build_proposer_prompt,
    classify_parse,
    parse_hypotheses,
)
from aria.llm.provider import TaskTier


def _agent_results() -> dict:
    return {
        "bulk_rna_agent": {
            "findings": {
                "low_power_warning": True,
                "contrasts": [
                    {
                        "name": "old_vs_young",
                        "top_genes": [
                            {"symbol": "GATA1", "log2fc": 2.3},
                            {"symbol": "KLF1", "log2fc": 1.9},
                        ],
                    }
                ],
            }
        }
    }


def _ledger() -> dict:
    return {
        "entries": [
            {"node_id": "ledger://bulk/differential_expression", "status": "ran"},
            {"node_id": "ledger://bulk/pathway_enrichment", "status": "ran"},
        ]
    }


def _good_hypothesis_json() -> str:
    return json.dumps(
        [
            {
                "id": "g1",
                "mechanism": "co-upregulation of GATA1 and KLF1 may reflect a "
                "shared erythroid program",
                "entities": ["GATA1", "KLF1"],
                "observation_refs": ["ledger://bulk/differential_expression"],
                "experiment": {
                    "perturbation": "GATA1 knockdown",
                    "readout": "KLF1 expression by qPCR",
                    "predicted_direction": "decrease",
                    "refuting_outcome": "KLF1 unchanged",
                },
                "devils_advocate": {
                    "simpler_explanation": "both respond to the same stimulus",
                    "confounds": ["low replication"],
                },
            }
        ]
    )


# ── parse + provenance ──────────────────────────────────────────────────────

def test_parse_hypotheses_handles_good_json():
    hyps = parse_hypotheses(_good_hypothesis_json())
    assert len(hyps) == 1
    assert hyps[0].entities == ["GATA1", "KLF1"]


def test_parse_hypotheses_is_lenient_on_garbage():
    assert parse_hypotheses("not json at all") == []
    assert parse_hypotheses("") == []
    assert parse_hypotheses("```json\n[]\n```") == []


def test_parse_hypotheses_strips_code_fences():
    fenced = "```json\n" + _good_hypothesis_json() + "\n```"
    assert len(parse_hypotheses(fenced)) == 1


def test_proposer_returns_empty_for_no_signals():
    proposer = LLMProposer(lambda p, s: _good_hypothesis_json())
    assert proposer([], {}) == []


def test_proposer_stamps_model_provenance():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    proposer = LLMProposer(
        lambda p, s: _good_hypothesis_json(), model_label="test-model"
    )
    hyps = proposer(signals, {})
    prov = hyps[0].provenance
    assert prov["model_label"] == "test-model"
    assert "prompt_sha256" in prov
    assert "input_evidence_sha256" in prov


def test_prompt_lists_only_audited_evidence():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    prompt = build_proposer_prompt(signals, {"biological_question": "aging?"}, 3)
    assert "GATA1" in prompt and "KLF1" in prompt
    assert "ledger://bulk/differential_expression" in prompt
    assert "aging?" in prompt


# ── full chain through the agent ────────────────────────────────────────────

def test_end_to_end_accepts_grounded_gated_hypothesis():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    agent = HypothesisAgent(
        proposer=LLMProposer(lambda p, s: _good_hypothesis_json())
    )
    out = agent.generate(signals, _ledger())
    assert [h["id"] for h in out["hypotheses"]] == ["g1"]
    assert out["honest_null"] is False
    assert out["quarantine"][0]["hypothesis_node"] == "hypothesis://g1"
    assert out["hypotheses"][0]["provenance"]["model_label"] == "aria-llm"


def test_end_to_end_rejects_model_inventing_a_gene():
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    invented = json.dumps(
        [
            {
                "id": "bad",
                "mechanism": "FOXP3 may regulate the program",
                "entities": ["FOXP3"],
                "observation_refs": ["ledger://bulk/differential_expression"],
                "experiment": {
                    "perturbation": "FOXP3 KO",
                    "readout": "qPCR",
                    "predicted_direction": "decrease",
                    "refuting_outcome": "no change",
                },
                "devils_advocate": {
                    "simpler_explanation": "noise",
                    "confounds": ["low replication"],
                },
            }
        ]
    )
    agent = HypothesisAgent(proposer=LLMProposer(lambda p, s: invented))
    out = agent.generate(signals, _ledger())
    assert out["hypotheses"] == []
    assert out["honest_null"] is True
    gates = {f["gate"] for f in out["rejected"][0]["failures"]}
    assert "grounding" in gates


# ── A+C: token budget + non-silent truncation diagnostics ───────────────────

def _verbose_hypotheses_json(n: int = 4) -> str:
    """A realistic multi-hypothesis response (mechanism + 4-field experiment +
    devils_advocate), the shape that overran the old 2048-token budget."""
    items = []
    for i in range(n):
        items.append(
            {
                "id": f"v{i}",
                "mechanism": (
                    "co-upregulation of the two flagged markers may reflect a "
                    "shared upstream regulatory program that could be probed by "
                    "targeted perturbation rather than direct interaction " * 3
                ),
                "entities": ["GATA1", "KLF1"],
                "observation_refs": ["ledger://bulk/differential_expression"],
                "experiment": {
                    "perturbation": "knock down the candidate regulator and "
                    "profile by RNA-seq " * 3,
                    "readout": "expression of the two markers by qPCR " * 3,
                    "predicted_direction": "decrease in both markers " * 3,
                    "refuting_outcome": "markers remain unchanged " * 3,
                },
                "devils_advocate": {
                    "simpler_explanation": "both respond to a shared stimulus " * 3,
                    "confounds": ["low replication"],
                },
            }
        )
    return "```json\n" + json.dumps(items, indent=2) + "\n```"


def _truncated_response() -> str:
    """A response cut off mid-object, exactly the max_tokens failure mode."""
    full = _verbose_hypotheses_json(4)
    return full[: int(len(full) * 0.55)]


def test_from_provider_uses_generous_token_budget():
    # A (root cause): the proposer must request enough tokens for a full set of
    # verbose hypotheses, not the old 2048 that truncated the JSON.
    seen: list[int] = []

    class _RecordingLLM:
        def complete(self, *, prompt, system, tier, max_tokens):
            assert tier is TaskTier.HEAVY
            seen.append(max_tokens)
            return _verbose_hypotheses_json(4)

    signals = bulk_rna_evidence(_agent_results(), _ledger())
    proposer = LLMProposer.from_provider(_RecordingLLM())
    hyps = proposer(signals, {})
    assert seen == [8192]
    assert len(hyps) == 4


def test_classify_parse_distinguishes_failure_modes():
    assert classify_parse("", 0)["status"] == "empty_response"
    assert classify_parse("```json\n[]\n```", 0)["status"] == "no_candidates"
    assert classify_parse(_truncated_response(), 0)["status"] == "parse_error"
    assert classify_parse(_verbose_hypotheses_json(4), 4)["status"] == "ok"


def test_proposer_records_truncation_diagnostic():
    # C: a truncated response yields zero candidates AND a non-silent diagnostic.
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    proposer = LLMProposer(lambda p, s: _truncated_response())
    assert proposer(signals, {}) == []
    assert proposer.last_diagnostics["status"] == "parse_error"
    assert proposer.last_diagnostics["raw_chars"] > 0


def test_agent_honest_null_names_a_parse_failure():
    # C: when the model answered but parsing failed, the honest-null is labelled
    # a generation failure (null_reason) instead of staying mute.
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    agent = HypothesisAgent(proposer=LLMProposer(lambda p, s: _truncated_response()))
    out = agent.generate(signals, _ledger())
    assert out["honest_null"] is True
    assert out["n_candidates"] == 0
    assert out["null_reason"] == "parse_error"
    assert out["proposer_diagnostics"]["status"] == "parse_error"
    # A genuine empty-list decline must NOT be labelled a failure.
    agent_ok = HypothesisAgent(proposer=LLMProposer(lambda p, s: "[]"))
    out_ok = agent_ok.generate(signals, _ledger())
    assert out_ok["honest_null"] is True
    assert "null_reason" not in out_ok


def test_end_to_end_rejects_model_dropping_a_confound():
    # The evidence carries a low-replication confound; a hypothesis that does not
    # acknowledge it is rejected by the devils_advocate gate.
    signals = bulk_rna_evidence(_agent_results(), _ledger())
    no_confound = json.dumps(
        [
            {
                "id": "drop",
                "mechanism": "GATA1 and KLF1 may share a regulator",
                "entities": ["GATA1", "KLF1"],
                "observation_refs": ["ledger://bulk/differential_expression"],
                "experiment": {
                    "perturbation": "GATA1 knockdown",
                    "readout": "KLF1 qPCR",
                    "predicted_direction": "decrease",
                    "refuting_outcome": "no change",
                },
                "devils_advocate": {
                    "simpler_explanation": "shared stimulus",
                    "confounds": [],
                },
            }
        ]
    )
    agent = HypothesisAgent(proposer=LLMProposer(lambda p, s: no_confound))
    out = agent.generate(signals, _ledger())
    assert out["hypotheses"] == []
    gates = {f["gate"] for f in out["rejected"][0]["failures"]}
    assert "devils_advocate" in gates
