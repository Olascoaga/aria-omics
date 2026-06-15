"""scATAC P1 preprint hardening: doublets, batch QC, peak provenance."""

import pytest


def test_atac_doublet_detector_flags_synthetic_doublets_and_keeps_clean_cells():
    from aria.utils.chromatin_robustness import detect_atac_doublets

    clean = []
    for i in range(56):
        row = [0.0] * 120
        for j in range(5):
            row[(i + j * 11) % 120] = 1.0
        clean.append(row)

    doublets = []
    for i in range(4):
        row = [0.0] * 120
        for j in range(90):
            row[(i + j) % 120] = 4.0
        doublets.append(row)

    res = detect_atac_doublets(clean + doublets, min_cells=50, max_rate=0.12)
    assert res["ran"] is True
    assert res["n_doublets"] == 4
    assert res["doublet_rate"] > 0
    assert sum(res["_doublet_mask"][:56]) == 0
    assert sum(res["_doublet_mask"][56:]) == 4


def test_atac_doublet_detector_skips_honestly_when_too_few_cells():
    from aria.utils.chromatin_robustness import detect_atac_doublets

    res = detect_atac_doublets([[1, 0], [0, 1]], min_cells=50)
    assert res["ran"] is False
    assert "needs >=" in res["reason"]
    assert res["n_doublets"] == 0
    assert res["_doublet_mask"] == [False, False]


def test_atac_batch_qc_warns_on_hidden_batch_and_batch_dominated_lsi():
    pytest.importorskip("sklearn")
    from aria.utils.chromatin_robustness import assess_atac_batch_embedding

    embedding = [[0.0, 0.0], [0.1, 0.0], [0.0, 0.1],
                 [5.0, 5.0], [5.1, 5.0], [5.0, 5.1]]
    batch = ["lane1", "lane1", "lane1", "lane2", "lane2", "lane2"]
    clusters = ["0", "0", "0", "1", "1", "1"]

    res = assess_atac_batch_embedding(
        obs_columns=["condition", "donor", "sequencing_lane"],
        embedding=embedding,
        batch_labels=batch,
        cluster_labels=clusters,
        condition_col="condition",
        replicate_col="donor",
        declared_batch=None,
    )
    checks = {i["check"] for i in res["issues"]}
    assert res["status"] == "warnings"
    assert "unmodeled_batch" in checks
    assert "residual_batch_effect" in checks
    assert res["metrics"]["batch_silhouette"] > 0.10


def test_consensus_peak_provenance_distinguishes_verified_from_unverified():
    from aria.utils.chromatin_robustness import assess_consensus_peak_provenance

    peaks = ["chr1:10-50", "chr1:70-100", "chr2:5-30"]
    verified = assess_consensus_peak_provenance(
        peaks,
        metadata={
            "method": "overlap_unified",
            "n_samples": 4,
            "overlap_fraction": 0.72,
            "rare_peak_policy": "preserve peaks present in >=2 samples",
        },
        input_kind="h5ad",
    )
    assert verified["status"] == "verified"
    assert verified["issues"] == []

    unverified = assess_consensus_peak_provenance(
        ["peakA", "peakA", "chr1:10-20"], input_kind="h5ad")
    checks = {i["check"] for i in unverified["issues"]}
    assert unverified["status"] == "unverified"
    assert "duplicate_peak_names" in checks
    assert "non_coordinate_peak_names" in checks
    assert "consensus_peak_provenance_unverified" in checks
