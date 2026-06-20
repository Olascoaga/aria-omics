"""F10 (preprint audit 2026-06-19): the bulk ATAC marker peak-to-gene mapping is
reproducible from declared parameters (it was hand-curated in the V47 artifact).

Also guards that the V47 DA artifact now scopes its claim honestly (plumbing/direction,
not calibration) and labels its `markers` block as illustrative + reproducible-via-the
script.
"""
import json
from pathlib import Path

from aria.benchmarks.bulk_atac_marker_mapping import (
    map_peaks_to_gene_markers,
    parse_peak_interval,
)

ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (ROOT / "docs" / "benchmark_results" / "bulk_atac"
            / "v47_k562_gm12878_da_validation.json")


def test_parse_peak_interval():
    assert parse_peak_interval("chr1:100-200") == ("chr1", 100, 200)
    assert parse_peak_interval("chrX:5-5") == ("chrX", 5, 5)
    assert parse_peak_interval("not-a-peak") is None
    assert parse_peak_interval("chr1:200-100") is None  # end < start
    assert parse_peak_interval(123) is None


def _da_rows():
    return [
        {"peak": "chr1:990-1010", "log2FoldChange": -4.0},   # 30 bp from TSS 1000
        {"peak": "chr1:1100-1200", "log2FoldChange": -5.0},  # 150 bp from TSS 1000
        {"peak": "chr1:9000-9100", "log2FoldChange": 2.0},   # far from TSS 1000
        {"peak": "chr2:1000-1100", "log2FoldChange": 9.0},   # wrong chromosome
        {"peak": "bad", "log2FoldChange": 1.0},              # unparseable
    ]


def test_mapping_is_deterministic_and_window_scoped():
    genes = [{"gene": "GENEA", "chrom": "chr1", "tss": 1000, "expected": "negative"}]
    # narrow window -> only the two nearby peaks (mean of -4 and -5)
    narrow = map_peaks_to_gene_markers(_da_rows(), genes, window_bp=500)
    assert narrow[0]["n_peaks_in_window"] == 2
    assert narrow[0]["mean_log2fc"] == -4.5
    assert narrow[0]["called_direction"] == "negative"
    assert narrow[0]["expected"] == "negative"   # pass-through preserved
    # same inputs -> identical output (reproducible)
    assert map_peaks_to_gene_markers(_da_rows(), genes, window_bp=500) == narrow
    # wider window pulls in the far same-chrom peak, changing the mean
    wide = map_peaks_to_gene_markers(_da_rows(), genes, window_bp=20000)
    assert wide[0]["n_peaks_in_window"] == 3
    assert wide[0]["mean_log2fc"] != narrow[0]["mean_log2fc"]


def test_mapping_no_peaks_returns_none_never_imputed():
    genes = [{"gene": "GENEZ", "chrom": "chr9", "tss": 42}]
    out = map_peaks_to_gene_markers(_da_rows(), genes, window_bp=100)
    assert out[0]["n_peaks_in_window"] == 0
    assert out[0]["mean_log2fc"] is None
    assert "called_direction" not in out[0]


def test_v47_artifact_scopes_claim_and_marks_markers_illustrative():
    data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    scope = data.get("validation_scope")
    assert isinstance(scope, dict)
    assert scope.get("validates")
    # the synthetic FDR/recovery gate is referenced as the calibration evidence
    assert "synthetic" in json.dumps(scope).lower()
    bs = data["biological_sanity"]
    mapping = bs.get("marker_peak_to_gene_mapping")
    assert isinstance(mapping, dict)
    assert "illustrative" in json.dumps(mapping).lower()
    assert "bulk_atac_marker_mapping" in mapping.get("reproducible_method", "")
