"""scATAC P4.3 governance tie-in: TOBIAS differential TF binding becomes a GOVERNED
associative narrative block with footprint<->RNA cross-evidence. The claim is capped
associative (no causal language, counts-only so W-CLAIM is satisfied), and the
cross-modal concordance is honest (unavailable RNA is not assumed concordant).
"""

from __future__ import annotations

import pytest


def _summary():
    return {
        "parsed": True, "n_motifs_tested": 879, "n_ranked_candidates": 873,
        "ranking_basis": {"fdr_controlled": False, "replicate_inference": False,
                          "pseudobulk": True},
        "top_toward_Monocyte": [
            {"tf": "CEBPB", "change": 1.01, "pvalue": 1e-195, "n_sites": 8987},
            {"tf": "FOSL2::JUND", "change": 0.48, "pvalue": 1e-188, "n_sites": 25363},
        ],
        "top_toward_T_cell": [
            {"tf": "EGR1", "change": -0.27, "pvalue": 1e-180, "n_sites": 56877},
            {"tf": "GHOST_TF", "change": -0.10, "pvalue": 1e-9, "n_sites": 30},
        ],
    }


def test_footprint_rna_concordance_scores_cross_modal():
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {
        "CEBPB": {"Monocyte": 3.0, "T_cell": 0.5},   # concordant (footprint+RNA Mono)
        "JUND":  {"Monocyte": 2.0, "T_cell": 1.0},   # dimer component, concordant Mono
        "EGR1":  {"Monocyte": 0.2, "T_cell": 1.5},   # concordant (footprint+RNA T_cell)
        # GHOST_TF has no RNA -> not checkable
    }
    res = footprint_rna_concordance(_summary(), rna, "Monocyte", "T_cell")
    assert res["n_checked"] == 3          # CEBPB, FOSL2::JUND(JUND), EGR1
    assert res["n_concordant"] == 3
    assert res["concordance_rate"] == pytest.approx(1.0)
    ghost = res["toward_T_cell"][1]
    assert ghost["rna_checked"] is False and ghost["rna_concordant"] is None


def test_footprint_rna_concordance_marks_discordant():
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {"CEBPB": {"Monocyte": 0.1, "T_cell": 2.0}}  # footprint Mono but RNA T_cell
    summary = {"parsed": True, "n_significant": 1, "n_motifs_tested": 1,
               "top_toward_Monocyte": [{"tf": "CEBPB", "change": 1.0,
                                        "pvalue": 1e-9, "n_sites": 50}],
               "top_toward_T_cell": []}
    res = footprint_rna_concordance(summary, rna, "Monocyte", "T_cell")
    assert res["n_checked"] == 1 and res["n_concordant"] == 0
    assert res["toward_Monocyte"][0]["rna_concordant"] is False


def test_footprint_block_is_governed_associative():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.validators import find_causal_language
    from aria.agents.narrative.synthesis.footprint_rna import footprint_rna_concordance

    rna = {"CEBPB": {"Monocyte": 3.0, "T_cell": 0.5},
           "JUND": {"Monocyte": 2.0, "T_cell": 1.0},
           "EGR1": {"Monocyte": 0.2, "T_cell": 1.5}}
    cross = footprint_rna_concordance(_summary(), rna, "Monocyte", "T_cell")
    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary(), "rna_cross_evidence": cross}

    block = ChromatinNarrator()._footprint_block(footprinting)
    assert block.status == "success" and block.confidence == "medium"
    assert "873" in block.claim and "879" in block.claim
    # cross-modal counts surfaced honestly
    assert "3" in block.claim  # n_concordant
    # associative only: no causal language in the claim; mandatory caveat present
    assert not find_causal_language(block.claim)
    assert any("associative" in c.text.lower() for c in block.caveats)
    # no specific TF name in the claim (W-CLAIM named-entity safety)
    assert "CEBPB" not in block.claim and "EGR1" not in block.claim


def test_footprint_block_surfaces_aggregate_figures():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary(),
                    "aggregate_plots": {
                        "CEBPB": {"png": "/tmp/x/CEBPB_aggregate.png",
                                  "svg": "/tmp/x/CEBPB_aggregate.svg"},
                        "EGR1": {"png": "/tmp/x/EGR1_aggregate.png"}}}
    block = ChromatinNarrator()._footprint_block(footprinting)
    paths = {f["path"] for f in block.figures}
    assert "/tmp/x/CEBPB_aggregate.png" in paths
    assert len(block.figures) == 2
    assert all("footprint" in f["caption"].lower() for f in block.figures)


def test_collect_emits_footprint_block_when_present():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    footprinting = {"ran": True, "group_a": "Monocyte", "group_b": "T_cell",
                    "differential_summary": _summary()}
    blocks = ChromatinNarrator().collect(
        "ChromatinAgent", {"findings": {"footprinting": footprinting}})
    ids = [b.id for b in blocks]
    assert "chromatin.differential_tf_footprinting" in ids


# ════════════════════════════════════════════════════════════════════════════
# ADR-078 Track 2 — production orchestration, narration and hypothesis grounding
# ════════════════════════════════════════════════════════════════════════════

import csv
import sys
from pathlib import Path


def _agent(env=None):
    from aria.agents.chromatin_agent import ChromatinAgent

    agent = ChromatinAgent.__new__(ChromatinAgent)
    agent.env = env
    agent.publish_status = lambda *a, **k: None
    agent.publish_finding = lambda *a, **k: None
    agent.publish_escalation = lambda *a, **k: None
    agent._publish_qc_finding = lambda *a, **k: None
    agent._publish_peaks_finding = lambda *a, **k: None
    return agent


class _RecordingEnv:
    """Fake EnvironmentManager: records every dispatch, answers per script."""

    def __init__(self, answers=None):
        self.calls = []
        self.answers = answers or {}

    def run_in_stack(self, stack, script_path, params, timeout=None):
        name = Path(script_path).name
        self.calls.append((stack, name, params))
        answer = self.answers.get(name)
        if callable(answer):
            return answer(params)
        if answer is not None:
            return dict(answer)
        return {"status": "success", "ran": False, "reason": "fake"}

    def names(self):
        return [name for _, name, _ in self.calls]

    def params_for(self, name):
        return next(p for _, n, p in self.calls if n == name)


def _no_local_assets(monkeypatch):
    from aria.utils import genomes, motifs

    monkeypatch.setattr(genomes, "resolve_local_genome_fasta",
                        lambda assembly: (None, None))
    monkeypatch.setattr(motifs, "load_local_motif_collection",
                        lambda collection: None)


def _inferential_footprinting(tmp_path=None, **overrides):
    fp = {
        "ran": True, "group_a": "Treated", "group_b": "Control",
        "group_label": "Conditions", "group_kind": "conditions",
        "validation_level": "production",
        "differential_summary": {
            "parsed": True, "n_motifs_tested": 50, "n_significant": 3,
            "alpha": 0.05,
            "inference": {
                "status": "success",
                "test": "welch_t_test_on_replicate_mean_scores",
                "multiple_testing": "benjamini_hochberg_across_motifs",
                "replicates_per_condition": {"Treated": 4, "Control": 4},
                "low_power_warning": False,
                "null_label_controls": {"n_permutations": 68},
            },
            "top_toward_Treated": [
                {"tf": "CEBPB", "motif_id": "MA0466.3", "change": 0.9,
                 "pvalue": 1e-6, "padj": 1e-4},
                {"tf": "FOSL2::JUND(var.2)", "motif_id": "MA1145.1",
                 "change": 0.5, "pvalue": 1e-4, "padj": 0.01},
            ],
            "top_toward_Control": [
                {"tf": "EGR1", "motif_id": "MA0162.4", "change": -0.4,
                 "pvalue": 1e-3, "padj": 0.03},
                {"tf": "NOT_SIG", "motif_id": "MA9999.1", "change": -0.2,
                 "pvalue": 0.2, "padj": 0.4},
            ],
        },
        "motif_source": {"collection": "JASPAR2024_CORE_vertebrates",
                         "release": "2024", "n_motifs": 879},
    }
    fp.update(overrides)
    return fp


# ── Part A: orchestration ────────────────────────────────────────────────────

def test_footprint_common_params_auto_resolves_genome_fai_and_motifs(
        tmp_path, monkeypatch):
    from aria.utils import genomes, motifs

    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n")
    (tmp_path / "genome.fa.fai").write_text("chr1\t4\t6\t4\t5\n")
    meme = tmp_path / "jaspar.meme"
    meme.write_text("MEME version 4\n")
    seen = {}

    def _resolve(assembly):
        seen["assembly"] = assembly
        return str(fasta), "managed:ARIA_GENOME_DIR"

    def _load(collection):
        seen["collection"] = collection
        return str(meme), {"collection": collection, "release": "2024",
                           "n_motifs": 879}

    monkeypatch.setattr(genomes, "resolve_local_genome_fasta", _resolve)
    monkeypatch.setattr(motifs, "load_local_motif_collection", _load)

    params = _agent()._footprint_common_params({"genome": "hg38"}, "out")
    assert params["genome_fasta"] == str(fasta)
    assert params["genome_fai"] == f"{fasta}.fai"
    assert params["motif_meme"] == str(meme)
    assert params["_motif_source"]["collection"] == motifs.DEFAULT_COLLECTION
    assert seen == {"assembly": "hg38", "collection": motifs.DEFAULT_COLLECTION}
    assert params["make_volcano_plot"] is True and params["resume"] is True

    # Explicit user paths always win and skip resolution entirely.
    monkeypatch.setattr(genomes, "resolve_local_genome_fasta",
                        lambda a: pytest.fail("must not resolve"))
    monkeypatch.setattr(motifs, "load_local_motif_collection",
                        lambda c: pytest.fail("must not resolve"))
    explicit = _agent()._footprint_common_params(
        {"genome_fasta": "/x/g.fa", "genome_fai": "/x/g.fa.fai",
         "motif_meme": "/x/m.meme", "footprint_make_volcano_plot": False,
         "footprint_resume": False}, "out")
    assert explicit["genome_fasta"] == "/x/g.fa"
    assert explicit["motif_meme"] == "/x/m.meme"
    assert explicit["_motif_source"] is None
    assert explicit["make_volcano_plot"] is False and explicit["resume"] is False


def _scatac_matrix_env(clustered):
    return _RecordingEnv({
        "chromatin_lsi_clustering.py": {"status": "success",
                                        "output_path": str(clustered)},
        "chromatin_diffacc.py": {"status": "success"},
        "chromatin_motifs.py": {"status": "success", "ran": False},
        "chromatin_regulatory.py": {"status": "success", "ran": False},
        "chromatin_footprint_tobias.py": {"status": "success", "ran": False,
                                          "reason": "TOBIAS not installed",
                                          "validation_level": "production"},
    })


@pytest.mark.parametrize("enabled", [None, False])
def test_scatac_enable_footprinting_opt_out(tmp_path, monkeypatch, enabled):
    _no_local_assets(monkeypatch)
    clustered = tmp_path / "clustered.h5ad"
    env = _scatac_matrix_env(clustered)
    ctx = {"genome": "hg38", "genome_fasta": str(tmp_path / "g.fa")}
    if enabled is not None:
        ctx["enable_footprinting"] = enabled
    res = _agent(env)._run_scatac_matrix(
        "exp", ctx, {}, str(tmp_path / "m.h5mu"), precomputed_qc={"status": "ok"})
    dispatched = "chromatin_footprint_tobias.py" in env.names()
    assert dispatched is (enabled is None)
    assert ("footprinting" in res["findings"]) is (enabled is None)
    if enabled is None:
        assert res["findings"]["footprinting"]["validation_level"] == "production"


@pytest.mark.parametrize("enabled", [True, False])
def test_bulk_enable_footprinting_opt_out(monkeypatch, enabled):
    agent = _agent(_RecordingEnv({
        "chromatin_bulk_diffacc.py": {"status": "success", "ran": True,
                                      "comparisons": []}}))
    for name in ("_run_bulk_atac_motifs", "_run_bulk_atac_annotation",
                 "_run_bulk_atac_ora"):
        monkeypatch.setattr(agent, name, lambda *a, **k: None)
    called = []
    monkeypatch.setattr(agent, "_run_bulk_atac_footprinting",
                        lambda *a, **k: called.append(1) or {"ran": False})
    findings = {}
    agent._run_bulk_atac_da_branch(
        "exp", {"enable_footprinting": enabled}, {}, {}, findings,
        {"status": "success", "counts_matrix_path": "c.tsv",
         "sample_metadata_path": "m.tsv"},
        peaks_result={"status": "success"}, files=["/d/a.bam"])
    assert bool(called) is enabled
    assert ("footprinting" in findings) is enabled


def _synthetic_h5ad(path, *, with_condition=True):
    import anndata
    import numpy as np
    import pandas as pd

    n = 8
    obs = {
        "donor": ["d1", "d2", "d3", "d4", "d1", "d2", "d3", "d4"],
        "batch": ["b1", "b1", "b2", "b2", "b1", "b1", "b2", "b2"],
        "predicted_celltype": ["Mono"] * 4 + ["T"] * 4,
        "leiden": ["0", "1"] * 4,
    }
    if with_condition:
        obs["condition"] = pd.Categorical(
            ["stim"] * 4 + ["ctrl"] * 3 + [None])
    var = pd.DataFrame(index=["chr1:100-200", "chr1_300_400", "chr2-500-600",
                              "not_a_peak", "chr3:900-800"])
    adata = anndata.AnnData(
        X=np.ones((n, 5), dtype="float32"),
        obs=pd.DataFrame(obs, index=[f"BC{i}-1" for i in range(n)]), var=var)
    adata.write_h5ad(path)
    return path


def test_derives_clean_peaks_bed_from_h5ad_var_names(tmp_path):
    h5ad = _synthetic_h5ad(tmp_path / "clustered.h5ad")
    bed = _agent()._derive_peaks_bed_from_h5ad(h5ad, tmp_path / "fp" / "peaks.bed")
    rows = [line.split("\t") for line in Path(bed).read_text().splitlines()]
    assert rows == [["chr1", "100", "200"], ["chr1", "300", "400"],
                    ["chr2", "500", "600"]]


def _read_tsv(path):
    with open(path) as fh:
        return list(csv.reader(fh, delimiter="\t"))


def test_derives_barcode_design_from_confirmed_condition_not_transferred_labels(
        tmp_path):
    h5ad = _synthetic_h5ad(tmp_path / "clustered.h5ad")
    ctx = {"design": {"main_factor": "condition"},
           "comparisons": [["stim", "ctrl"]], "replicate_col": "donor"}
    tsv, ga, gb = _agent()._derive_scatac_barcode_design(
        ctx, h5ad, tmp_path / "fp" / "derived_barcode_design.tsv")
    rows = _read_tsv(tsv)
    assert rows[0] == ["barcode", "group", "replicate", "batch"]
    assert (ga, gb) == ("stim", "ctrl")
    assert rows[1] == ["BC0-1", "stim", "d1", "b1"]
    assert len(rows) == 1 + 7  # the null-condition barcode is dropped
    groups = {r[1] for r in rows[1:]}
    assert groups == {"stim", "ctrl"}  # never predicted_celltype / leiden

    # The derived TSV round-trips through the backend's design loader.
    from aria.scripts.chromatin_footprint_tobias import _load_barcode_design

    design = _load_barcode_design(tsv)
    assert design["BC5-1"] == {"group": "ctrl", "replicate": "d2", "batch": "b1"}


def test_barcode_design_never_guesses_contrast_or_trusts_transferred_labels(
        tmp_path):
    h5ad = _synthetic_h5ad(tmp_path / "with_cond.h5ad")
    # A two-level condition column without an explicit contrast: no guessed pair.
    tsv, ga, gb = _agent()._derive_scatac_barcode_design(
        {"design": {"main_factor": "condition"}}, h5ad, tmp_path / "a.tsv")
    assert Path(tsv).is_file() and ga is None and gb is None

    # Only transferred labels / Leiden available -> no inferential grouping.
    bare = _synthetic_h5ad(tmp_path / "bare.h5ad", with_condition=False)
    assert _agent()._derive_scatac_barcode_design(
        {"groupby": "predicted_celltype"}, bare, tmp_path / "b.tsv") is None
    # ... unless the user explicitly asks for that column.
    tsv, ga, gb = _agent()._derive_scatac_barcode_design(
        {"footprint_groupby": "predicted_celltype", "group_a": "Mono",
         "group_b": "T"}, bare, tmp_path / "c.tsv")
    assert (ga, gb) == ("Mono", "T")
    assert {r[1] for r in _read_tsv(tsv)[1:]} == {"Mono", "T"}


def test_h5ad_obs_reader_works_without_anndata(tmp_path, monkeypatch):
    from aria.agents.chromatin_agent import ChromatinAgent

    h5ad = _synthetic_h5ad(tmp_path / "clustered.h5ad")
    expected = ChromatinAgent._read_h5ad_obs(h5ad)
    monkeypatch.setitem(sys.modules, "anndata", None)  # import -> ImportError
    names, columns = ChromatinAgent._read_h5ad_obs(h5ad)
    assert names == expected[0]
    assert columns["condition"] == expected[1]["condition"]
    assert columns["condition"][-1] is None
    assert columns["donor"] == expected[1]["donor"]


def test_scatac_footprinting_dispatch_derives_inputs_and_stamps_production(
        tmp_path, monkeypatch):
    from aria.utils import genomes, motifs

    h5ad = _synthetic_h5ad(tmp_path / "clustered.h5ad")
    monkeypatch.setattr(genomes, "resolve_local_genome_fasta",
                        lambda a: ("/g/genome.fa", "managed"))
    monkeypatch.setattr(motifs, "load_local_motif_collection",
                        lambda c: ("/m/jaspar.meme", {"collection": c}))
    ran = {"status": "success", "ran": True, "group_a": "stim", "group_b": "ctrl",
           "differential_summary": {"parsed": True,
                                    "top_toward_stim": [{"tf": "CEBPB"}],
                                    "top_toward_ctrl": [{"tf": "EGR1"}]}}
    env = _RecordingEnv({"chromatin_footprint_tobias.py": ran})
    ctx = {"design": {"main_factor": "condition"}, "replicate_col": "donor",
           "comparisons": [["stim", "ctrl"]],
           "fragments_files": ["/f/s1.tsv.gz", "/f/s2.tsv.gz"],
           "sample_ids": ["s1", "s2"],
           "rna_group_means": {"CEBPB": {"stim": 3.0, "ctrl": 1.0},
                               "EGR1": {"stim": 0.1, "ctrl": 2.0}}}
    result = _agent(env)._run_scatac_footprinting(
        "exp", ctx, str(h5ad), str(tmp_path / "aria_chromatin"))
    params = env.params_for("chromatin_footprint_tobias.py")
    assert params["peaks_bed"].endswith("footprinting/derived_peaks.bed")
    assert params["barcode_groups"].endswith(
        "footprinting/derived_barcode_design.tsv")
    assert (params["group_a"], params["group_b"]) == ("stim", "ctrl")
    assert params["fragments_files"] == ["/f/s1.tsv.gz", "/f/s2.tsv.gz"]
    assert params["fragments_file"] == "/f/s1.tsv.gz"
    assert params["sample_ids"] == ["s1", "s2"]
    assert params["genome_fasta"] == "/g/genome.fa"
    assert "_motif_source" not in params
    assert result["validation_level"] == "production"
    assert result["motif_source"]["collection"] == motifs.DEFAULT_COLLECTION
    assert result["rna_cross_evidence"]["n_concordant"] == 2


def test_rna_cross_evidence_failure_never_breaks_footprinting(tmp_path):
    result = {"ran": True, "group_a": "A", "group_b": "B",
              "differential_summary": {"parsed": True,
                                       "top_toward_A": [{"tf": "CEBPB"}]}}
    bad = tmp_path / "rna.h5ad"
    bad.write_text("not an h5ad")
    _agent()._attach_footprint_rna_cross_evidence(result, {"rna_data_path": str(bad)})
    assert "rna_cross_evidence" not in result
    descriptive = {"ran": False, "differential_summary": {"parsed": True}}
    _agent()._attach_footprint_rna_cross_evidence(
        descriptive, {"rna_group_means": {"CEBPB": {"A": 1, "B": 0}}})
    assert "rna_cross_evidence" not in descriptive


def test_bulk_raw_fastq_propagates_aligned_sample_ids_into_replicate_design(
        tmp_path, monkeypatch):
    _no_local_assets(monkeypatch)
    bams = []
    for sid in ("S1", "S2", "S3", "S4", "S5", "S6"):
        bam = tmp_path / f"{sid}.filtered.bam"
        bam.write_bytes(b"")
        bams.append(str(bam))
    preprocessing = {
        "status": "success", "ran": True, "analysis": "bulk_atac_alignment",
        "bam_files": [{"name": Path(b).name.split(".")[0],
                       "sample_id": Path(b).name.split(".")[0],
                       "status": "success", "bam": b} for b in bams],
    }
    meta = {f"S{i}": {"condition": "treated" if i <= 3 else "control",
                      "replicate": f"r{(i - 1) % 3 + 1}",
                      "batch": "b1" if i in (1, 4) else "b2"}
            for i in range(1, 7)}
    env = _RecordingEnv({
        "chromatin_qc.py": {"status": "success", "data_type": "bulk_ATAC"},
        "chromatin_peaks.py": {"status": "success",
                               "peaks_path": str(tmp_path / "peaks.narrowPeak")},
        "chromatin_peak_counts.py": {"status": "success",
                                     "counts_matrix_path": "c.tsv",
                                     "sample_metadata_path": "m.tsv"},
        "chromatin_bulk_diffacc.py": {
            "status": "success", "ran": True, "data_type": "bulk_ATAC",
            "comparisons": [{"test": "treated", "reference": "control",
                             "status": "success"}]},
    })
    agent = _agent(env)
    monkeypatch.setattr(agent, "_run_bulk_atac_preprocessing",
                        lambda *a, **k: (list(bams), dict(preprocessing)))
    for name in ("_run_bulk_atac_motifs", "_run_bulk_atac_annotation",
                 "_run_bulk_atac_ora"):
        monkeypatch.setattr(agent, name, lambda *a, **k: None)
    fastqs = [str(tmp_path / f"S{i}_R1.fastq.gz") for i in range(1, 7)]
    res = agent._run_bulk_atac(
        "exp", {"sample_metadata": meta, "comparisons": [["treated", "control"]],
                "design": {"main_factor": "condition", "batch_covariate": "batch"}},
        {}, fastqs)
    assert res["status"] == "done"
    assert env.params_for("chromatin_peak_counts.py")["sample_ids"] == [
        "S1", "S2", "S3", "S4", "S5", "S6"]
    params = env.params_for("chromatin_footprint_tobias.py")
    assert params["replicate_bams"] == {
        "treated": {"r1": [bams[0]], "r2": [bams[1]], "r3": [bams[2]]},
        "control": {"r1": [bams[3]], "r2": [bams[4]], "r3": [bams[5]]},
    }
    assert params["replicate_batches"] == {
        "treated": {"r1": "b1", "r2": "b2", "r3": "b2"},
        "control": {"r1": "b1", "r2": "b2", "r3": "b2"},
    }
    assert "reference_fasta" not in params  # BAM input: no CRAM reference
    assert res["findings"]["footprinting"]["validation_level"] == "production"


def test_bulk_replicate_design_still_refuses_filename_identity():
    agent = _agent()
    ctx = {"sample_metadata": {"lane1": {"condition": "A", "replicate": "r1",
                                         "batch": "b1"}}}
    assert agent._bulk_replicate_bams(ctx, ["lane1.bam"]) == {}
    assert agent._bulk_replicate_batches(ctx, ["lane1.bam"]) == {}


def test_bulk_replicate_batches_leaves_disagreeing_technical_libraries_unadjusted():
    agent = _agent()
    ctx = {"sample_ids": ["a1", "a2", "b1"],
           "sample_metadata": {
               "a1": {"condition": "A", "replicate": "d1", "seq_batch": "x"},
               "a2": {"condition": "A", "replicate": "d1", "seq_batch": "y"},
               "b1": {"condition": "B", "replicate": "d2", "seq_batch": "x"}}}
    batches = agent._bulk_replicate_batches(
        ctx, ["a1.bam", "a2.bam", "b1.bam"], batch_col="seq_batch")
    assert batches == {"B": {"d2": "x"}}
    assert agent._bulk_replicate_batches(
        {**ctx, "replicate_batches": {"A": {"d1": "z"}}}, []) == {"A": {"d1": "z"}}


def test_bulk_cram_footprinting_forwards_reference_fasta(monkeypatch):
    _no_local_assets(monkeypatch)
    env = _RecordingEnv()
    _agent(env)._run_bulk_atac_footprinting(
        "exp", {"genome_fasta": "/ref/genome.fa"}, {"comparisons": []},
        {"peaks_path": "/p/peaks.bed"}, ["/d/s1.cram"], contract={})
    assert env.params_for("chromatin_footprint_tobias.py")[
        "reference_fasta"] == "/ref/genome.fa"


# ── Part B: narration ────────────────────────────────────────────────────────

def test_footprint_block_surfaces_volcano_and_production_metadata():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator
    from aria.agents.narrative.validators import find_causal_language

    fp = _inferential_footprinting(
        volcano_plot={"png": "/r/footprint_volcano.png",
                      "svg": "/r/footprint_volcano.svg"},
        aggregate_plots={"CEBPB": {"png": "/r/CEBPB_aggregate.png"}})
    block = ChromatinNarrator()._footprint_block(fp)
    ids = [f["id"] for f in block.figures]
    assert ids == ["footprint_volcano", "footprint_CEBPB"]
    assert all("footprint" in f["caption"].lower() for f in block.figures)
    assert block.metadata["validation_level"] == "production"
    assert not find_causal_language(block.claim)
    assert "CEBPB" not in block.claim and "EGR1" not in block.claim
    # Without an explicit level, the production default applies.
    fp.pop("validation_level")
    assert ChromatinNarrator()._footprint_block(fp).metadata[
        "validation_level"] == "production"


@pytest.mark.parametrize("test_name,expected", [
    ("paired_t_test_on_within_donor_mean_score_differences", "paired t-test"),
    ("ols_t_test_on_replicate_mean_scores_batch_adjusted",
     "batch-adjusted linear model"),
    ("welch_t_test_on_replicate_mean_scores", "Welch t-test"),
])
def test_footprint_block_names_the_fitted_replicate_test(test_name, expected):
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    fp = _inferential_footprinting()
    fp["differential_summary"]["inference"]["test"] = test_name
    block = ChromatinNarrator()._footprint_block(fp)
    method = next(e.value for e in block.evidence if e.label == "Inferential method")
    assert method.startswith(expected)


def test_narrator_footprint_tables_and_figures_from_generate_figures(tmp_path):
    from aria.agents._narrative_chromatin import generate_figures
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    volcano = tmp_path / "footprint_volcano.png"
    volcano.write_bytes(b"png")
    agg = tmp_path / "FOSL2JUND_aggregate.png"
    agg.write_bytes(b"png")
    fp = _inferential_footprinting(
        volcano_plot={"png": str(volcano)},
        aggregate_plots={"FOSL2::JUND": {"png": str(agg)},
                         "MISSING": {"png": str(tmp_path / "nope.png")}},
        bindetect_results=str(tmp_path / "bindetect_results.txt"))
    fp["differential_summary"]["results_table"] = str(
        tmp_path / "replicate_footprint_inference.tsv")
    findings = {"footprinting": fp}
    generate_figures(findings, None, tmp_path / "figs")
    figs = findings["figures"]
    assert figs["footprint_volcano"] == str(volcano)
    assert figs["footprint_FOSL2__JUND"] == str(agg)
    assert not any("MISSING" in k for k in figs)

    tables = {t["id"]: t["path"] for t in ChromatinNarrator().tables(
        "chromatin_agent", {"findings": findings})}
    assert tables["chromatin.footprint_inference"].endswith(
        "replicate_footprint_inference.tsv")
    assert tables["chromatin.footprint_bindetect"].endswith("bindetect_results.txt")
    # Not-run footprinting surfaces nothing.
    assert not any(k.startswith("chromatin.footprint") for k in (
        t["id"] for t in ChromatinNarrator().tables(
            "chromatin_agent", {"findings": {"footprinting": {"ran": False}}})))


def test_methods_paragraph_bulk_inferential_and_scatac_descriptive():
    from aria.agents.narrative.narrators.chromatin import ChromatinNarrator

    narrator = ChromatinNarrator()
    bulk = narrator.methods("chromatin_agent", {"findings": {
        "qc": {"status": "success", "data_type": "bulk_ATAC"},
        "footprinting": _inferential_footprinting()}})
    text = next(p for p in bulk if "TOBIAS ATACorrect" in p)
    for token in ("ScoreBigwig", "BINDetect", "JASPAR2024_CORE_vertebrates",
                  "Welch t-test", "Benjamini–Hochberg", "<= 0.05",
                  "68 design-respecting", "Treated: 4"):
        assert token in text
    assert "not FDR-controlled" not in text

    descriptive = _inferential_footprinting(motif_source=None)
    descriptive["differential_summary"]["inference"] = {
        "status": "descriptive_only",
        "reason": "missing_explicit_replicate_or_donor_identity"}
    sc = narrator.methods("chromatin_agent", {"findings": {
        "lsi": {"status": "success"}, "footprinting": descriptive}})
    text = next(p for p in sc if "TOBIAS ATACorrect" in p)
    assert "descriptive and not FDR-controlled" in text
    assert "missing_explicit_replicate_or_donor_identity" in text
    assert "supplied local MEME motif collection" in text
    assert not any("TOBIAS" in p for p in narrator.methods(
        "chromatin_agent", {"findings": {"footprinting": {"ran": False}}}))


# ── Part C: hypothesis grounding ─────────────────────────────────────────────

def test_footprint_signals_emit_only_bh_significant_tfs_with_caveats():
    from aria.agents.narrative.hypothesis.adapters._chromatin_common import (
        FOOTPRINT_NODE, footprint_signals,
    )
    from aria.agents.narrative.hypothesis.caveats import (
        FOOTPRINT_ASSOCIATIVE, PARALOG_MOTIF_AMBIGUITY,
    )

    sigs = footprint_signals(_inferential_footprinting(), modality="bulk_ATAC")
    by = {s.entity: s for s in sigs}
    assert list(by) == ["CEBPB", "FOSL2", "JUND", "EGR1"]  # NOT_SIG dropped
    assert all(s.entity_kind == "tf_motif" and s.audited_node_ref == FOOTPRINT_NODE
               for s in sigs)
    assert by["CEBPB"].direction == "up" and by["EGR1"].direction == "down"
    assert by["CEBPB"].measure == "footprint_score(Treated_vs_Control)"
    assert by["CEBPB"].context == "Treated_vs_Control"
    assert by["CEBPB"].caveats_inherited == [FOOTPRINT_ASSOCIATIVE]
    for sub in ("FOSL2", "JUND"):
        assert by[sub].caveats_inherited == [FOOTPRINT_ASSOCIATIVE,
                                             PARALOG_MOTIF_AMBIGUITY]


def test_footprint_signals_refuse_descriptive_and_flag_power_and_batch(tmp_path):
    from aria.agents.narrative.hypothesis.adapters._chromatin_common import (
        footprint_signals,
    )
    from aria.agents.narrative.hypothesis.caveats import BATCH, LOW_REPLICATION

    for status in ("descriptive_only", "not_run", None):
        fp = _inferential_footprinting()
        fp["differential_summary"]["inference"]["status"] = status
        assert footprint_signals(fp, modality="scATAC") == []
    assert footprint_signals({"ran": False}, modality="scATAC") == []
    assert footprint_signals(None, modality="scATAC") == []

    fp = _inferential_footprinting()
    inf = fp["differential_summary"]["inference"]
    inf.update(low_power_warning=True, batch_adjusted=True,
               test="ols_t_test_on_replicate_mean_scores_batch_adjusted")
    sig = footprint_signals(fp, modality="scATAC")[0]
    assert LOW_REPLICATION in sig.caveats_inherited
    assert BATCH in sig.caveats_inherited

    # The persisted BH table (all motifs) is preferred over the top-N lists.
    table = tmp_path / "replicate_footprint_inference.tsv"
    table.write_text(
        "motif_id\ttf\tmean_score_a\tmean_score_b\tchange\tpvalue\tpadj\tn_sites\n"
        "M1\tRUNX1\t1\t0\t1.2\t1e-5\t0.002\t50\n"
        "M2\tSPI1\t0\t1\t-0.8\t1e-4\t0.002\t50\n"
        "M3\tNULLTF\t0\t0\t0.01\t0.9\t0.95\t50\n")
    fp = _inferential_footprinting()
    fp["differential_summary"]["results_table"] = str(table)
    assert [s.entity for s in footprint_signals(fp, modality="scATAC")] == [
        "RUNX1", "SPI1"]


def test_bulk_atac_evidence_grounds_footprints_even_without_motif_enrichment():
    from aria.agents.narrative.hypothesis import bulk_atac_evidence
    from aria.agents.narrative.hypothesis.caveats import LOW_REPLICATION
    from aria.agents.narrative.run_ledger import build_run_ledger

    agent_results = {"chromatin_agent": {"findings": {
        "differential_accessibility": {
            "data_type": "bulk_ATAC", "ran": True, "status": "success",
            "comparisons": [{"test": "Treated", "reference": "Control",
                             "low_power_warning": True}]},
        "footprinting": _inferential_footprinting(),
    }}}
    sigs = bulk_atac_evidence(agent_results)
    assert {s.entity for s in sigs} == {"CEBPB", "FOSL2", "JUND", "EGR1"}
    assert all(s.modality == "bulk_ATAC" for s in sigs)
    assert all(LOW_REPLICATION in s.caveats_inherited for s in sigs)

    ledger = build_run_ledger({}, agent_results)
    assert {s.entity for s in bulk_atac_evidence(agent_results, ledger)} == {
        "CEBPB", "FOSL2", "JUND", "EGR1"}
    # A ledger in which the footprint node did not run grounds nothing.
    assert bulk_atac_evidence(agent_results, {"nodes": []}) == []

    descriptive = _inferential_footprinting()
    descriptive["differential_summary"]["inference"]["status"] = "descriptive_only"
    assert bulk_atac_evidence({"chromatin_agent": {"findings": {
        "footprinting": descriptive}}}) == []


def test_scatac_evidence_adds_footprint_layer():
    from aria.agents.narrative.hypothesis import scatac_evidence

    fp = _inferential_footprinting(group_a="Monocyte", group_b="T_cell")
    summary = fp["differential_summary"]
    summary["top_toward_Monocyte"] = summary.pop("top_toward_Treated")
    summary["top_toward_T_cell"] = summary.pop("top_toward_Control")
    sigs = scatac_evidence({"chromatin_agent": {"findings": {"footprinting": fp}}})
    assert [s.entity for s in sigs] == ["CEBPB", "FOSL2", "JUND", "EGR1"]
    assert all(s.modality == "scATAC" and s.context == "Monocyte_vs_T_cell"
               for s in sigs)


def test_devils_advocate_requires_owning_footprint_caveats():
    from aria.agents.narrative.hypothesis import (
        Hypothesis, build_signals_by_entity, check_devils_advocate,
    )
    from aria.agents.narrative.hypothesis.adapters._chromatin_common import (
        footprint_signals,
    )
    from aria.agents.narrative.hypothesis.caveats import (
        CAVEAT_CODES, FOOTPRINT_ASSOCIATIVE, caveat_gloss,
    )

    assert FOOTPRINT_ASSOCIATIVE in CAVEAT_CODES
    assert "Tn5-protected motif occupancy" in caveat_gloss(FOOTPRINT_ASSOCIATIVE)
    index = build_signals_by_entity(
        footprint_signals(_inferential_footprinting(), modality="bulk_ATAC"))

    def _hyp(entities, confounds):
        return Hypothesis(id="h", mechanism="m", entities=entities,
                          devils_advocate={"simpler_explanation": "a paralog TF",
                                           "confounds": confounds})

    hidden = check_devils_advocate(_hyp(["CEBPB"], []), index)
    assert not hidden.passed and "footprint_associative" in hidden.reason
    assert check_devils_advocate(
        _hyp(["CEBPB"], ["Tn5 footprint occupancy is associative, not causal "
                         "regulation"]), index).passed
    assert check_devils_advocate(
        _hyp(["CEBPB"], ["footprint_associative"]), index).passed
    # A dimer subunit also carries paralog ambiguity, which must be owned too.
    dimer = check_devils_advocate(_hyp(["JUND"], ["tn5 occupancy"]), index)
    assert not dimer.passed and "paralog_motif_ambiguity" in dimer.reason
    assert check_devils_advocate(
        _hyp(["JUND"], ["tn5 occupancy", "paralogous dimer subunits share the "
                                         "motif"]), index).passed
