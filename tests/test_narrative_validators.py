from pathlib import Path

import pytest


def _success_block(**overrides):
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    data = {
        "id": "scrna.pseudobulk.Monocytes.STIM_vs_CTRL",
        "modality": "scRNA-seq",
        "analysis": "pseudobulk_de",
        "block_type": "result",
        "title": "Monocytes STIM_vs_CTRL",
        "status": "success",
        "confidence": "high",
        "claim": "Monocytes had 140 DE genes.",
        "evidence": [EvidenceItem("global-FDR DE genes", 140)],
    }
    data.update(overrides)
    return NarrativeBlock(**data)


def test_success_claim_and_evidence_are_fail_hard():
    from aria.agents.narrative.types import NarrativeBlock
    from aria.agents.narrative.validators import (
        NarrativeValidationError,
        validate_block,
    )

    block = object.__new__(NarrativeBlock)
    block.id = "bad"
    block.status = "success"
    block.claim = ""
    block.evidence = []
    block.confidence = "high"
    block.caveats = []
    block.error = None
    block.warnings = []
    block.metadata = {}
    block.figures = []
    block.tables = []
    block.analysis = "qc"

    with pytest.raises(NarrativeValidationError, match="requires claim"):
        validate_block(block)


def test_non_success_must_explain_failure():
    from aria.agents.narrative.types import NarrativeBlock
    from aria.agents.narrative.validators import (
        NarrativeValidationError,
        validate_block,
    )

    block = NarrativeBlock(
        id="scrna.marker_discovery",
        modality="scRNA-seq",
        analysis="marker_discovery",
        block_type="error",
        title="Marker discovery",
        status="error",
        confidence="insufficient",
        claim="",
    )

    with pytest.raises(NarrativeValidationError, match="requires claim"):
        validate_block(block)

    block.error = "Wilcoxon timed out."
    validate_block(block)


def test_low_confidence_gets_limitation_warning():
    from aria.agents.narrative.validators import validate_block

    block = _success_block(confidence="low")
    validate_block(block)
    assert any("limitations section" in w for w in block.warnings)


def test_causal_language_degrades_and_adds_caveat_without_causal_evidence():
    from aria.agents.narrative.validators import validate_block

    block = _success_block(
        claim="ISG15 drives monocyte response.",
        confidence="high",
    )
    validate_block(block)
    assert block.confidence == "medium"
    assert any("drives" in caveat.text for caveat in block.caveats)

    causal = _success_block(
        claim="Perturbation directly regulates target expression.",
        confidence="high",
        metadata={"causal_evidence": True},
    )
    validate_block(causal)
    assert causal.confidence == "high"
    assert not causal.caveats


def test_trajectory_without_velocity_gets_caveat_not_failure():
    from aria.agents.narrative.validators import validate_block

    block = _success_block(
        id="scrna.trajectory",
        analysis="trajectory",
        block_type="exploratory",
        title="Trajectory",
        claim="PAGA/DPT ordered groups.",
        metadata={"velocity_computed": False},
    )
    validate_block(block)
    assert any("PAGA/DPT is exploratory" in c.text for c in block.caveats)


def test_referenced_tables_and_figures_must_exist_at_render_time(tmp_path):
    from aria.agents.narrative.validators import (
        NarrativeValidationError,
        validate_block,
    )

    fig = tmp_path / "figures" / "plot.png"
    tbl = tmp_path / "tables" / "table.tsv"
    fig.parent.mkdir()
    tbl.parent.mkdir()
    fig.write_bytes(b"png")
    tbl.write_text("x\n", encoding="utf-8")

    block = _success_block(
        figures=[{"path": "figures/plot.png"}],
        tables=[{"path": "tables/table.tsv"}],
    )
    validate_block(block, base_dir=tmp_path, check_files=True)

    missing = _success_block(figures=[{"path": "figures/missing.png"}])
    with pytest.raises(NarrativeValidationError, match="referenced figure"):
        validate_block(missing, base_dir=tmp_path, check_files=True)


def test_registry_collects_first_accepting_narrator():
    from aria.agents.narrative.registry import NarrativeRegistry

    class Narrator:
        name = "dummy"

        def accepts(self, agent_name, agent_result):
            return agent_name == "dummy_agent"

        def collect(self, agent_name, agent_result, context=None):
            return [_success_block(id="dummy.result", modality="dummy")]

        def methods(self, agent_name, agent_result, context=None):
            return []

        def figures(self, agent_name, agent_result, report_dir=None):
            return []

        def tables(self, agent_name, agent_result, report_dir=None):
            return []

    registry = NarrativeRegistry()
    registry.register(Narrator())
    blocks = registry.collect_blocks({"dummy_agent": {"status": "done"}})
    assert [b.id for b in blocks] == ["dummy.result"]
