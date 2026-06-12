def _synthetic_scrna_findings():
    return {
        "qc": {
            "n_cells_before": 29173,
            "n_cells_after": 29173,
            "pct_removed": 0.0,
        },
        "differential_expression": {
            "status": "error",
            "error_type": "Timeout",
            "details": "Execution exceeded 3600s limit.",
        },
        "differential_abundance": {
            "per_comparison": {
                "condition_a_vs_condition_b": {
                    "status": "success",
                    "n_significant": 1,
                    "per_cell_type": [{
                        "name": "GroupA",
                        "significant": True,
                        "direction": "up",
                    }],
                }
            }
        },
        "pseudobulk_de": {
            "groupby": "cluster",
            "condition_col": "condition",
            "replicate_col": "sample_id",
            "n_groups": 3,
            "per_group": {
                "GroupA": {
                    "n_pseudosamples": 16,
                    "per_comparison": {
                        "condition_a_vs_condition_b": {
                            "status": "success",
                            "n_significant": 180,
                            "n_significant_local": 220,
                            "n_significant_global": 180,
                            "n_up": 100,
                            "n_up_global": 95,
                            "n_down": 80,
                            "n_down_global": 85,
                            "corrected_for_composition": True,
                            "power_estimate_at_lfc_min": 0.78,
                            "top_genes": [
                                {"gene": "GENE_UP_1", "log2fc": 2.1},
                                {"gene": "GENE_UP_2", "log2fc": 1.8},
                                {"gene": "GENE_DOWN_1", "log2fc": -1.2},
                            ],
                        }
                    },
                },
                "GroupB": {
                    "n_pseudosamples": 16,
                    "per_comparison": {
                        "condition_a_vs_condition_b": {
                            "status": "success",
                            "n_significant_global": 12,
                            "n_significant_local": 20,
                            "n_up_global": 8,
                            "n_down_global": 4,
                            "corrected_for_composition": False,
                            "top_genes": [{"gene": "GENE_UP_3", "log2fc": 1.4}],
                        }
                    },
                },
            },
        },
        "pseudobulk_pathways": {
            "per_cluster": {
                "GroupA::condition_a_vs_condition_b": {
                    "n_significant": 3,
                    "results": {
                        "GO_BP": [{
                            "term": "pathway_alpha_response",
                            "adjusted_p": 1e-6,
                        }]
                    },
                }
            }
        },
        "cell_communication": {
            "status": "done",
            "method": "liana_rank_aggregate (specificity_rank)",
            "n_cell_types": 3,
            "n_interactions": 10,
            "n_autocrine_dropped": 4,
            "top_interactions": [{
                "source": "GroupA",
                "target": "GroupB",
                "ligand": "LIGAND_A",
                "receptor": "RECEPTOR_B",
            }],
        },
        "trajectory": {
            "status": "done",
            "paga": {"n_connections": 6, "n_strong": 1},
            "pseudotime": {
                "computed": True,
                "pseudotime_by_group": {"GroupA": 0.2, "GroupB": 0.8},
            },
            "velocity": {"computed": False},
        },
    }


def test_scrna_narrator_generates_blocks_for_all_synthetic_results():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator
    from aria.agents.narrative.validators import validate_blocks

    agent_result = {
        "status": "done",
        "findings": {"scRNA": {"findings": _synthetic_scrna_findings()}},
    }
    blocks = validate_blocks(
        ScrnaNarrator().collect("scrna_agent", agent_result)
    )
    ids = {block.id for block in blocks}

    assert "scrna.qc" in ids
    assert "scrna.marker_discovery" in ids
    assert "scrna.composition.condition_a_vs_condition_b" in ids
    assert "scrna.pseudobulk.GroupA.condition_a_vs_condition_b" in ids
    assert "scrna.pseudobulk.GroupB.condition_a_vs_condition_b" in ids
    assert "scrna.pathway.GroupA_condition_a_vs_condition_b" in ids
    assert "scrna.cellcomm" in ids
    assert "scrna.trajectory" in ids

    mono = next(b for b in blocks
                if b.id == "scrna.pseudobulk.GroupA.condition_a_vs_condition_b")
    assert mono.claim == "GroupA condition_a_vs_condition_b had 180 global-FDR DE genes."
    assert any(ev.value == "GENE_UP_1" or "GENE_UP_1" in ev.label for ev in mono.evidence)
    assert mono.metrics["corrected_for_composition"] is True
    assert any("composition covariate" in caveat.text for caveat in mono.caveats)

    marker = next(b for b in blocks if b.id == "scrna.marker_discovery")
    assert marker.status == "error"
    assert marker.error == "Execution exceeded 3600s limit."

    trajectory = next(b for b in blocks if b.id == "scrna.trajectory")
    assert trajectory.confidence == "low"
    assert any("PAGA/DPT is exploratory" in c.text for c in trajectory.caveats)


def test_scrna_narrator_surfaces_degraded_abundance_fdr_caveat():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    findings = _synthetic_scrna_findings()
    comp = findings["differential_abundance"]["per_comparison"][
        "condition_a_vs_condition_b"
    ]
    comp.update({
        "degraded": True,
        "confidence": "degraded",
        "fdr_family": "donor_level_clr_only",
        "n_fdr_tests": 2,
        "n_fisher_diagnostic": 1,
        "degradation_reason": (
            "1 cell type used Fisher exact as a cell-level diagnostic and "
            "was excluded from donor-level FDR."
        ),
        "caveats": [
            "1 cell type used Fisher exact as a cell-level diagnostic and "
            "was excluded from donor-level FDR."
        ],
    })

    agent_result = {"status": "done",
                    "findings": {"scRNA": {"findings": findings}}}
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    da = next(
        b for b in blocks
        if b.id == "scrna.composition.condition_a_vs_condition_b"
    )

    assert da.confidence == "low"
    assert any("Fisher exact" in caveat.text for caveat in da.caveats)


def test_scrna_narrator_flags_lognorm_recovered_counts():
    """F-SCI-LOGNORM (audit 2026-05-28): when DESeq2 counts were
    reverse-engineered from log-normalized values, every pseudobulk block must
    carry a visible recovery caveat."""
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    findings = _synthetic_scrna_findings()
    findings["pseudobulk_de"]["lognorm_recovered"] = True
    findings["pseudobulk_de"]["count_source"] = "recovered_from_lognorm"

    agent_result = {"status": "done",
                    "findings": {"scRNA": {"findings": findings}}}
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    pb = next(b for b in blocks
              if b.id == "scrna.pseudobulk.GroupA.condition_a_vs_condition_b")
    assert any("reverse-engineered" in c.text.lower() for c in pb.caveats), \
        [c.text for c in pb.caveats]
    # F-SCI-LOGNORM: recovered-count DE must not be presented at the same trust
    # level as raw-count DE.
    assert pb.confidence == "low"


def test_scrna_narrator_uses_per_cluster_fdr_label_when_strategy_set():
    """F-SCI-FDR (audit 2026-05-28): when the pseudobulk script declared the
    per-cluster FDR strategy, the narrative claim must name that family (not the
    legacy global-FDR wording) and report the primary count."""
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    findings = _synthetic_scrna_findings()
    comp = findings["pseudobulk_de"]["per_group"]["GroupA"]["per_comparison"][
        "condition_a_vs_condition_b"
    ]
    comp["fdr_strategy"] = "per_cluster"
    comp["n_significant"] = 220          # primary = per-cluster (local) count
    comp["n_significant_local"] = 220

    agent_result = {"status": "done",
                    "findings": {"scRNA": {"findings": findings}}}
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    pb = next(b for b in blocks
              if b.id == "scrna.pseudobulk.GroupA.condition_a_vs_condition_b")
    assert pb.claim == "GroupA condition_a_vs_condition_b had 220 per-cluster FDR DE genes."
    # Both families must remain visible as audit evidence.
    labels = {ev.label: ev.value for ev in pb.evidence}
    assert labels.get("global-FDR DE genes") == 180
    assert labels.get("local-FDR DE genes") == 220


def test_scrna_narrator_surfaces_design_matrix_warnings():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    findings = _synthetic_scrna_findings()
    comp = findings["pseudobulk_de"]["per_group"]["GroupA"]["per_comparison"][
        "condition_a_vs_condition_b"
    ]
    comp["design_check"] = {
        "status": "warnings",
        "issues": [{
            "severity": "warning",
            "check": "n1_design_cells",
            "message": "Some condition x covariate design cells have one sample.",
        }],
    }

    agent_result = {"status": "done",
                    "findings": {"scRNA": {"findings": findings}}}
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    pb = next(b for b in blocks
              if b.id == "scrna.pseudobulk.GroupA.condition_a_vs_condition_b")

    assert any("Design-matrix warning" in c.text for c in pb.caveats)


def test_scrna_narrator_emits_data_quality_block_for_qc_flags():
    """X8/X9: integration-overcorrection and annotation-coherence flags must
    surface as a visible data-quality limitation block."""
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    findings = _synthetic_scrna_findings()
    findings["integration_qc"] = {
        "status": "warnings",
        "issues": [{
            "severity": "warning",
            "check": "possible_overcorrection",
            "message": "Cluster silhouette is negative; overcorrection.",
            "recommendation": "Lower integration strength.",
        }],
    }
    findings["annotation_qc"] = {
        "status": "unverified",
        "issues": [{
            "severity": "warning",
            "check": "annotation_unverified",
            "message": "Reused obs annotations were not marker-verified.",
            "recommendation": "Validate canonical markers manually.",
        }],
    }
    agent_result = {"status": "done",
                    "findings": {"scRNA": {"findings": findings}}}
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    dq = next((b for b in blocks if b.id == "scrna.data_quality"), None)
    assert dq is not None
    assert dq.block_type == "limitation"
    texts = " ".join(c.text for c in dq.caveats)
    assert "overcorrection" in texts.lower()
    assert "marker-verified" in texts.lower()


def test_scrna_narrator_methods_reuse_legacy_methods():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    agent_result = {
        "status": "done",
        "findings": {"scRNA": {"findings": _synthetic_scrna_findings()}},
    }
    methods = ScrnaNarrator().methods("scrna_agent", agent_result)
    assert len(methods) == 1
    assert "pseudobulk aggregation" in methods[0]


# ── B-DD1 (scRNA-lane audit 2026-06-11): per-cluster marker discovery is a
# double-dipped (selection-then-test) within-data ranking. It must be surfaced
# as a DESCRIPTIVE ranking with an explicit circularity caveat — never as
# inferential significance — and the ClaimCompiler must cap it to descriptive.

def _synthetic_marker_success_findings():
    findings = _synthetic_scrna_findings()
    findings["differential_expression"] = {
        "status": "success",
        "groupby": "leiden",
        "n_clusters": 4,
        "n_significant": 320,
        "n_significant_genes": 320,
        "n_sig_by_cluster": {"0": 120, "1": 90, "2": 70, "3": 40},
        "de_genes_by_cluster": {
            "0": [{"gene": "MARKER_A", "log2fc": 3.1, "padj": 1e-9}],
        },
        "padj_max": 0.05,
        "lfc_min": 0.5,
    }
    return findings


def test_scrna_marker_block_is_descriptive_with_double_dip_caveat():
    from aria.agents.narrative.narrators.scrna import ScrnaNarrator

    agent_result = {
        "status": "done",
        "findings": {"scRNA": {"findings": _synthetic_marker_success_findings()}},
    }
    blocks = ScrnaNarrator().collect("scrna_agent", agent_result)
    marker = next((b for b in blocks if b.id == "scrna.marker_discovery"), None)
    assert marker is not None
    assert marker.status == "success"
    # Markers are an exploratory ranking, not a stats-gated DE result.
    assert marker.block_type == "exploratory"
    assert marker.analysis == "marker_discovery"
    assert marker.confidence == "low"

    # The claim must NOT present markers as significance.
    assert "significant" not in marker.claim.lower()

    # An explicit double-dipping / circularity caveat must be present.
    caveat_text = " ".join(c.text for c in marker.caveats).lower()
    assert ("double-dip" in caveat_text or "double dip" in caveat_text
            or "selection" in caveat_text or "circular" in caveat_text)


def test_claim_compiler_caps_marker_discovery_to_descriptive():
    from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
    from aria.agents.narrative.claim_compiler import classify_claim

    block = NarrativeBlock(
        id="scrna.marker_discovery",
        modality="scRNA-seq",
        analysis="marker_discovery",
        block_type="exploratory",
        title="Per-cluster marker discovery",
        status="success",
        confidence="low",
        claim="Per-cluster marker discovery ranked candidate markers across 4 clusters.",
        evidence=[EvidenceItem(label="clusters", value=4, source="marker_discovery")],
    )
    c = classify_claim(block)
    assert c.tier == "descriptive"
    assert c.licensed_language == "descriptive"
