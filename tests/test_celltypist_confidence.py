"""N-ANNO2 (scRNA annotation audit 2026-06-12): the per-cluster annotation
confidence must be GENUINE, not structurally 1.0.

With majority_voting=True and over_clustering=leiden, CellTypist assigns one
label per Leiden cluster, so a per-cluster `frequency` computed over the
MV-collapsed labels is 1.0 by construction. The confidence proxy must instead
reflect the RAW per-cell predicted labels and the model's probability matrix.
"""

import math

import pytest

from aria.scripts.rna_celltypist import _summarize_annotation_confidence


def test_frequency_reflects_raw_label_agreement_not_collapsed_majority():
    # Cluster "0": majority-voted (assigned) label collapsed every cell to
    # "T cell", but the RAW per-cell predictions disagree (6 T, 4 NK) and the
    # model was not confident (prob of assigned label ~0.45).
    cluster_ids = ["0"] * 10
    assigned = ["T cell"] * 10                      # MV-collapsed
    raw = ["T cell"] * 6 + ["NK cell"] * 4          # genuine per-cell calls
    conf = [0.45] * 10

    per_cluster = _summarize_annotation_confidence(
        cluster_ids, raw, assigned, conf,
    )

    info = per_cluster["0"]
    assert info["label"] == "T cell"
    # Structural-1.0 bug would report 1.0; genuine agreement is 6/10.
    assert info["frequency"] == pytest.approx(0.6)
    assert info["mean_confidence"] == pytest.approx(0.45)
    # The runner-up RAW label must be visible, not silently collapsed away.
    assert info["alt_labels"], "alt_labels must surface the raw disagreement"
    assert any(a["label"] == "NK cell" for a in info["alt_labels"])


def test_confident_homogeneous_cluster_stays_high():
    cluster_ids = ["1"] * 8
    assigned = ["B cell"] * 8
    raw = ["B cell"] * 8
    conf = [0.97] * 8

    info = _summarize_annotation_confidence(cluster_ids, raw, assigned, conf)["1"]
    assert info["frequency"] == pytest.approx(1.0)
    assert info["mean_confidence"] == pytest.approx(0.97)
    assert info["alt_labels"] == []


def test_missing_probabilities_degrade_to_none_confidence():
    cluster_ids = ["0"] * 4
    assigned = ["T cell"] * 4
    raw = ["T cell"] * 3 + ["NK cell"]
    conf = [float("nan")] * 4   # probability matrix unavailable

    info = _summarize_annotation_confidence(cluster_ids, raw, assigned, conf)["0"]
    assert info["frequency"] == pytest.approx(0.75)
    assert info["mean_confidence"] is None   # not faked as 0 or 1
