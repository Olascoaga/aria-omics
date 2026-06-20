"""B1 (GEO connector for ATAC): the connector must recognize and bucket the two
ATAC modalities, not only RNA.

Before B1 the connector classified only bulk_RNA/scRNA/unknown and bucketed only
counts/h5ad/h5/mtx — ATAC studies were invisible (mis-typed as RNA, fragments
mis-bucketed as count tables). These guards pin the ATAC detection + file
classification so all four ARIA-processable modalities are reachable from GEO.
"""

from aria.connectors.geo_connector import (
    _classify_suppl_file,
    _infer_data_type,
)


# ── data-type inference ───────────────────────────────────────────────────

def test_bulk_atac_from_library_strategy():
    md = {"title": "K562 ATAC-seq", "samples": [],
          "library_strategy": "ATAC-seq"}
    assert _infer_data_type(md, {}) == "bulk_ATAC"


def test_scatac_from_single_cell_keywords():
    md = {"title": "Single-cell ATAC-seq of PBMC (10x)", "samples": [],
          "library_strategy": "ATAC-seq"}
    assert _infer_data_type(md, {}) == "scATAC"


def test_scatac_from_snatac_keyword():
    md = {"title": "snATAC-seq hippocampus", "samples": [],
          "library_strategy": ""}
    assert _infer_data_type(md, {}) == "scATAC"


def test_scatac_from_fragments_file():
    md = {"title": "chromatin accessibility", "samples": [],
          "library_strategy": ""}
    files = {"fragments": ["/x/GSM_fragments.tsv.gz"]}
    assert _infer_data_type(md, files) == "scATAC"


def test_atac_beats_rna_branch_for_10x_atac():
    # A 10x scATAC study carries the same "10x"/"single cell" tokens the RNA
    # branch keys on; ATAC must win.
    md = {"title": "10x single cell ATAC of brain", "samples": [],
          "library_strategy": "ATAC-seq"}
    assert _infer_data_type(md, {}) == "scATAC"


def test_rna_still_classified():
    assert _infer_data_type(
        {"title": "bulk RNA-seq", "samples": [], "library_strategy": "RNA-Seq"},
        {}) == "bulk_RNA"
    assert _infer_data_type(
        {"title": "scRNA-seq 10x", "samples": [], "library_strategy": ""},
        {}) == "scRNA"


def test_unknown_when_no_signal():
    assert _infer_data_type({"title": "mystery", "samples": [],
                             "library_strategy": ""}, {}) == "unknown"


def test_original_atac_paper_descriptive_title():
    # Real GSE47753 (Buenrostro 2013): library_strategy="OTHER" (pre-dates the
    # GEO ATAC-seq term) and a descriptive title without the literal "atac".
    md = {"title": ("Transposition of native chromatin for fast and sensitive "
                    "epigenomic profiling of open chromatin"),
          "samples": [], "library_strategy": "OTHER"}
    assert _infer_data_type(md, {}) == "bulk_ATAC"


def test_atacseq_no_separator():
    assert _infer_data_type({"title": "ATACseq of K562", "samples": [],
                             "library_strategy": ""}, {}) == "bulk_ATAC"


def test_dnase_open_chromatin_not_atac():
    # Negative control: DNase-seq also profiles "open chromatin" but must NOT be
    # mis-typed as ATAC (we deliberately do not key on "open chromatin").
    md = {"title": "DNase-seq open chromatin of GM12878", "samples": [],
          "library_strategy": "DNase-Hypersensitivity"}
    assert _infer_data_type(md, {}) == "unknown"


# ── supplementary-file bucketing ──────────────────────────────────────────

def test_fragments_bucketed_not_counts():
    # The decisive bug: fragments.tsv.gz matches the generic .tsv.gz counts
    # pattern and must be caught first.
    assert _classify_suppl_file("GSM123_fragments.tsv.gz") == "fragments"


def test_peaks_buckets():
    assert _classify_suppl_file("GSM_peaks.narrowPeak.gz") == "peaks"
    assert _classify_suppl_file("sample.broadPeak") == "peaks"
    assert _classify_suppl_file("GSM_peaks.bed.gz") == "peaks"


def test_bam_bucket():
    assert _classify_suppl_file("aligned.bam") == "bam"
    assert _classify_suppl_file("aligned.cram") == "bam"


def test_rna_buckets_unchanged():
    assert _classify_suppl_file("counts.h5ad") == "h5ad"
    assert _classify_suppl_file("filtered_feature_bc_matrix.h5") == "h5"
    assert _classify_suppl_file("matrix.mtx.gz") == "mtx"
    assert _classify_suppl_file("GSE_counts.tsv.gz") == "counts"
    assert _classify_suppl_file("readme.pdf") is None


# ── data-audit routing: GEO data_type -> ARIA modality ────────────────────

def test_geo_bucket_map_routes_atac_to_chromatin():
    from aria.agents.data_audit_agent import _geo_bucket_map
    sc = _geo_bucket_map("scATAC")
    assert sc["fragments"] == "scATAC" and sc["bam"] == "scATAC"
    bulk = _geo_bucket_map("bulk_ATAC")
    assert bulk["bam"] == "bulk_ATAC" and bulk["peaks"] == "bulk_ATAC"
    # ATAC must NOT be collapsed into the RNA lane.
    assert "bulk_RNA" not in sc.values() and "scRNA" not in bulk.values()


def test_geo_bucket_map_rna_unchanged():
    from aria.agents.data_audit_agent import _geo_bucket_map
    assert _geo_bucket_map("scRNA")["counts"] == "scRNA"
    assert _geo_bucket_map("bulk_RNA")["counts"] == "bulk_RNA"
    assert _geo_bucket_map("unknown")["counts"] == "bulk_RNA"
