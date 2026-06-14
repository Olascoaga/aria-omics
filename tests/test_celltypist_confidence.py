"""N-ANNO2 (scRNA annotation audit 2026-06-12): the per-cluster annotation
confidence must be GENUINE, not structurally 1.0.

With majority_voting=True and over_clustering=leiden, CellTypist assigns one
label per Leiden cluster, so a per-cluster `frequency` computed over the
MV-collapsed labels is 1.0 by construction. The confidence proxy must instead
reflect the RAW per-cell predicted labels and the model's probability matrix.
"""

import math
import sys
import types

import pytest

from aria.scripts.rna_celltypist import (
    _resolve_celltypist_model,
    rna_celltypist,
    _summarize_annotation_confidence,
)


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


# ── N-ANNO3 (scRNA annotation audit 2026-06-12): annotating any tissue with the
# immune default model when no hint was given is a silent wrong-model risk. The
# fallback must be disclosed, not hidden.

def test_explicit_model_is_not_flagged():
    model, source, warning = _resolve_celltypist_model(
        "Adult_Human_Kidney.pkl", "Homo sapiens", "")
    assert model == "Adult_Human_Kidney.pkl"
    assert source == "explicit"
    assert warning is None


def test_tissue_hint_match_is_not_flagged():
    model, source, warning = _resolve_celltypist_model(None, "Homo sapiens", "pbmc")
    assert model == "Immune_All_Low.pkl"
    assert source == "tissue_hint"
    assert warning is None


def test_no_hint_falls_back_to_immune_default_with_warning():
    model, source, warning = _resolve_celltypist_model(None, "Homo sapiens", "")
    assert model == "Immune_All_Low.pkl"
    assert source == "default_immune_fallback"
    assert warning and "immune" in warning.lower()


def test_unknown_tissue_hint_also_flagged():
    model, source, warning = _resolve_celltypist_model(None, "Homo sapiens", "tumor")
    assert source == "default_immune_fallback"
    assert warning is not None


def test_mouse_no_hint_flags_human_immune_model_mismatch():
    model, source, warning = _resolve_celltypist_model(None, "Mus musculus", "")
    assert source == "default_immune_fallback"
    # The default is a human immune model on a mouse dataset — must be disclosed.
    assert warning and ("mouse" in warning.lower() or "human" in warning.lower())


def test_celltypist_blocks_unknown_tissue_without_override(tmp_path):
    result = rna_celltypist({
        "data_path": str(tmp_path / "input.h5ad"),
        "organism": "Homo sapiens",
        "tissue_hint": "tumor",
        "output_dir": str(tmp_path),
    })

    assert result["status"] == "error"
    assert result["error_type"] == "DefaultImmuneModelRequiresAck"
    assert result["requires_ack"] is True
    assert result["ack_param"] == "allow_default_immune_model"
    assert result["model_source"] == "default_immune_fallback"
    assert "Immune_All_Low.pkl" in result["details"]
    assert result["tissue_hint"] == "tumor"


def test_celltypist_model_download_respects_air_gapped(monkeypatch, tmp_path):
    calls = []

    class FakeModels:
        @staticmethod
        def download_models(**kwargs):
            calls.append(kwargs)
            raise AssertionError("download_models must not be called")

    monkeypatch.setitem(
        sys.modules,
        "celltypist",
        types.SimpleNamespace(models=FakeModels),
    )
    monkeypatch.setitem(
        sys.modules,
        "scanpy",
        types.SimpleNamespace(pp=types.SimpleNamespace()),
    )
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")

    result = rna_celltypist({
        "data_path": str(tmp_path / "input.h5ad"),
        "organism": "Homo sapiens",
        "model": "ARIA_T3_Not_Cached_Model.pkl",
        "output_dir": str(tmp_path),
    })

    assert calls == []
    assert result["status"] == "error"
    assert result["error_type"] == "EgressBlocked"
    assert "celltypist_model_hub" in result["details"]
