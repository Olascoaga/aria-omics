import pytest


def test_success_block_requires_claim_and_evidence():
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock

    with pytest.raises(ValueError, match="requires claim"):
        NarrativeBlock(
            id="scrna.qc",
            modality="scRNA-seq",
            analysis="qc",
            block_type="qc",
            title="QC",
            status="success",
            confidence="high",
            claim="",
            evidence=[EvidenceItem(label="cells", value=100)],
        )

    with pytest.raises(ValueError, match="requires evidence"):
        NarrativeBlock(
            id="scrna.qc",
            modality="scRNA-seq",
            analysis="qc",
            block_type="qc",
            title="QC",
            status="success",
            confidence="high",
            claim="QC passed",
        )


def test_block_round_trips_to_dict_for_methodology_json():
    from aria.agents.narrative.types import (
        Caveat,
        EvidenceItem,
        NarrativeBlock,
    )

    block = NarrativeBlock(
        id="scrna.pseudobulk.Monocytes.STIM_vs_CTRL",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="Monocytes STIM_vs_CTRL",
        status="success",
        confidence="medium",
        claim="Monocytes had 140 global-FDR DE genes.",
        evidence=[
            EvidenceItem(
                label="global-FDR DE genes",
                value=140,
                source="pseudobulk_de",
                path="tables/scrna_pseudobulk_de_summary.tsv",
            )
        ],
        caveats=[Caveat("Interpret as differential expression.", "info")],
        metrics={"n_up": 80, "n_down": 60},
        figures=[{"path": "figures/de.png", "caption": "DE bar"}],
        tables=[{"path": "tables/de.tsv", "label": "DE genes"}],
        methods=["pyDESeq2 pseudobulk"],
        warnings=["low-ish power"],
        metadata={"block_type_convention": "result"},
    )

    data = block.to_dict()
    rebuilt = NarrativeBlock.from_dict(data)

    assert rebuilt == block
    assert rebuilt.to_dict()["evidence"][0]["label"] == "global-FDR DE genes"
    assert rebuilt.to_dict()["caveats"][0]["severity"] == "info"
