"""v4.6 scATAC steps 5-6 — ChromatinNarrator (LSI / DA / motifs) + agent wiring.

The narrator turns the step 2-4 script outputs into validated NarrativeBlocks
(W-CLAIM strict render must pass), routes them into the Chromatin report section,
and the agent's `.h5mu` matrix flow dispatches the scripts and stores findings
under the keys the narrator + run-ledger read. Honest skips (no replicates, no
genome) become limitation blocks, never fabricated results (ADR-002/ADR-011).
"""

import pytest

from aria.agents.narrative.narrators.chromatin import ChromatinNarrator


# Findings shaped exactly like the real HC11 step 2-4 outputs.
def _findings(motifs_ran=True, pseudobulk_ran=False):
    da_by_cluster = {str(i): (1 + i * 3) for i in range(8)}  # all nonzero
    f = {
        "qc": {
            "status": "success", "data_type": "scATAC",
            "n_cells": 3143, "n_peaks": 60990, "mito_fraction": 0.04,
            "frip": None, "tss_enrichment": None, "qc_complete": False,
            "pass_qc": None,
            "metrics_not_computed": ["frip", "tss_enrichment"],
        },
        "lsi": {
            "status": "success", "input_kind": "h5mu", "atac_modality": "atac",
            "n_cells_total": 3143, "n_cells_used": 3143, "n_peaks": 60990,
            "n_components_computed": 50, "dropped_components": [0],
            "n_components_used": 49, "rep_used": "X_lsi", "n_clusters": 8,
            "cluster_sizes": {str(i): 300 for i in range(8)},
            "resolution": 1.0, "sketch_used": False,
            "output_path": "/tmp/x/lsi_clustered.h5ad",
        },
        "differential_accessibility": {
            "status": "success", "groupby": "leiden", "n_clusters": 8,
            "padj_max": 0.05, "lfc_min": 0.5,
            "per_cluster": {
                "ran": True, "n_da_total": 13294,
                "n_da_by_cluster": da_by_cluster,
                "output_csv": "/tmp/x/chromatin_da_per_cluster.csv",
            },
            "pseudobulk": (
                {"ran": True, "comparisons": [
                    {"test": "B", "reference": "A", "status": "success",
                     "n_sig": 120, "n_up": 80, "n_down": 40}],
                 "output_csv": "/tmp/x/pb.csv"}
                if pseudobulk_ran else
                {"ran": False,
                 "reason": "no usable condition column (None not in obs)"}
            ),
        },
        "motifs": (
            {"status": "success", "ran": True, "method": "hypergeometric",
             "genome_fasta": "/g/GRCh38.fa",
             "motif_source": {"collection": "JASPAR2024_CORE_vertebrates",
                              "release": "2024", "n_motifs": 879,
                              "sha256": "dd494278"},
             "n_groups": 8,
             "per_group": {str(i): {"n_enriched": (518 if i == 2 else 0),
                                    "top_motifs": []} for i in range(8)},
             "output_csv": "/tmp/x/chromatin_motif_enrichment.csv"}
            if motifs_ran else
            {"status": "success", "ran": False,
             "reason": "motif collection not staged; run scripts/fetch_motifs.py"}
        ),
    }
    return f


def _collect(findings):
    n = ChromatinNarrator()
    return n.collect("chromatin_agent", {"status": "done", "findings": findings},
                     {})


def test_narrator_emits_clustering_da_and_motif_blocks():
    blocks = _collect(_findings())
    ids = {b.id for b in blocks}
    assert "chromatin.clustering" in ids
    assert "chromatin.differential_accessibility.per_cluster" in ids
    assert "chromatin.motifs" in ids
    # pseudobulk honestly skipped -> a limitation block, not a fabricated result
    assert "chromatin.differential_accessibility.pseudobulk_skipped" in ids


def test_clustering_claim_numbers_are_all_in_evidence():
    block = next(b for b in _collect(_findings())
                 if b.id == "chromatin.clustering")
    vals = {str(e.value) for e in block.evidence}
    assert {"8", "3143", "60990", "49", "1"} <= vals
    assert block.confidence == "medium"


def test_motif_block_is_association_only_with_provenance():
    block = next(b for b in _collect(_findings()) if b.id == "chromatin.motifs")
    text = " ".join(c.text.lower() for c in block.caveats)
    assert "association-only" in text
    assert "chromvar" in text          # per-cell activity disclosed as out-of-scope
    labels = {e.label for e in block.evidence}
    assert "Motif collection" in labels and "Motifs tested" in labels


def test_blocks_pass_strict_wclaim_render():
    # The strict render is the W-CLAIM gate: every claim must be supported by
    # its evidence card or it raises NarrativeValidationError.
    from aria.agents.narrative.render_blocks import render_blocks
    html = render_blocks(_collect(_findings()), strict=True)
    assert "Chromatin" in html or "chromatin" in html.lower()


def test_blocks_group_into_chromatin_section():
    from aria.agents.narrative.render_blocks import group_blocks_by_prefix
    groups = group_blocks_by_prefix(_collect(_findings()))
    assert "chromatin" in groups and groups["chromatin"]


def test_motif_skip_is_a_limitation_block():
    blocks = _collect(_findings(motifs_ran=False))
    mb = next(b for b in blocks if b.id == "chromatin.motifs.skipped")
    assert mb.status == "limitation" and mb.claim == ""
    assert "not run" in mb.caveats[0].text.lower()


def test_registry_collects_and_validates_chromatin_blocks():
    from aria.agents.narrative.registry import registry_with
    reg = registry_with((ChromatinNarrator(),))
    blocks = reg.collect_blocks(
        {"chromatin_agent": {"status": "done", "findings": _findings()}}, {})
    assert any(b.id == "chromatin.motifs" for b in blocks)


def test_methods_cover_each_ran_stage():
    n = ChromatinNarrator()
    methods = n.methods("chromatin_agent",
                        {"status": "done", "findings": _findings()}, {})
    joined = " ".join(methods).lower()
    assert "tf-idf" in joined and "differential accessibility" in joined
    assert "motif enrichment" in joined


# ── Agent .h5mu matrix flow (step 6 wiring) ───────────────────────────────────

class _FakeEnv:
    def __init__(self):
        self.calls = []

    def run_in_stack(self, stack, script_path, params):
        self.calls.append((script_path, params))
        f = _findings()
        if script_path.endswith("chromatin_qc.py"):
            return f["qc"]
        if script_path.endswith("chromatin_lsi_clustering.py"):
            return f["lsi"]
        if script_path.endswith("chromatin_diffacc.py"):
            return f["differential_accessibility"]
        if script_path.endswith("chromatin_motifs.py"):
            return f["motifs"]
        return {"status": "error", "error_type": "Unexpected"}


def _bare_agent(env):
    from aria.agents.chromatin_agent import ChromatinAgent
    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    return agent


def test_agent_h5mu_flow_dispatches_and_stores_findings():
    env = _FakeEnv()
    agent = _bare_agent(env)
    exp_ctx = {"genome": "hg38", "organism": "Homo sapiens",
               "genome_fasta": "/g/GRCh38.fa",
               "motif_collection": "JASPAR2024_CORE_vertebrates"}
    res = agent._run_scatac(  # routes to _run_scatac_matrix on a .h5mu
        "exp1", exp_ctx, {}, ["/data/hc11_paired.h5mu"])
    assert res["status"] == "done"
    findings = res["findings"]
    assert set(findings) >= {"qc", "lsi", "differential_accessibility", "motifs"}

    # dispatched scripts in order, on the chromatin stack
    scripts = [Path_basename(c[0]) for c in env.calls]
    assert scripts == ["chromatin_qc.py", "chromatin_lsi_clustering.py",
                       "chromatin_diffacc.py", "chromatin_motifs.py"]
    # genome + motif collection forwarded to the motif script
    motif_params = next(p for s, p in env.calls if s.endswith("chromatin_motifs.py"))
    assert motif_params["genome_fasta"] == "/g/GRCh38.fa"
    assert motif_params["motif_collection"] == "JASPAR2024_CORE_vertebrates"

    # and the narrator turns those findings into the chromatin blocks
    blocks = ChromatinNarrator().collect("chromatin_agent", res, {})
    assert {"chromatin.clustering", "chromatin.motifs"} <= {b.id for b in blocks}


def Path_basename(p):
    from pathlib import Path
    return Path(p).name
