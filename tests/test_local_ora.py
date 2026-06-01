"""P1-7 — local, offline hypergeometric ORA (W-PRIV).

These tests run in the light PR lane (no gseapy, no network). They prove:
  * the hypergeometric p-value matches a hand-computed reference;
  * GMT parsing + versioned-manifest loading work;
  * `run_ora` enriches a planted term and BH-corrects;
  * Enrichr is opt-in (default off) and a missing local library degrades
    honestly (no fabrication) rather than silently calling the network.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from aria.utils import ora


# ── hypergeometric core ──────────────────────────────────────────────────────

def _ref_sf(k, M, n, N):
    """Reference P(X>=k) via math.comb (exact, slow) for cross-checking."""
    from math import comb
    denom = comb(M, N)
    return sum(
        comb(n, i) * comb(M - n, N - i) / denom
        for i in range(k, min(n, N) + 1)
    )


def test_hypergeom_sf_matches_reference():
    for (k, M, n, N) in [(3, 100, 10, 20), (1, 50, 5, 5), (5, 200, 30, 25)]:
        assert ora.hypergeom_sf(k, M, n, N) == pytest.approx(
            _ref_sf(k, M, n, N), rel=1e-9, abs=1e-12
        )


def test_hypergeom_sf_edges():
    assert ora.hypergeom_sf(0, 100, 10, 20) == 1.0          # k<=0 -> certain
    assert ora.hypergeom_sf(21, 100, 10, 20) == 0.0         # k > min(n,N)


# ── GMT parsing + versioned library loading ──────────────────────────────────

def _write_library(tmp_path, name, terms, manifest=None):
    base = tmp_path / name
    base.mkdir(parents=True)
    lines = []
    for term, genes in terms.items():
        lines.append("\t".join([term, "", *genes]))
    (base / f"{name}.gmt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    if manifest is not None:
        (base / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    return base


def test_parse_gmt(tmp_path):
    base = _write_library(tmp_path, "L", {"T1": ["A", "b", "C"], "T2": ["D"]})
    sets = ora.parse_gmt(base / "L.gmt")
    assert sets == {"T1": ["A", "B", "C"], "T2": ["D"]}   # upper-cased


def test_load_local_library_with_manifest(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    _write_library(
        tmp_path, "GO_BP_2021", {"path_x": ["A", "B", "C"]},
        manifest={"library": "GO_BP_2021", "source": "Enrichr",
                  "release": "2021", "date": "2026-06-01"},
    )
    loaded = ora.load_local_library("GO_BP_2021")
    assert loaded is not None
    gene_sets, version = loaded
    assert gene_sets == {"path_x": ["A", "B", "C"]}
    assert version["release"] == "2021"
    assert version["source"] == "Enrichr"
    assert version["n_terms"] == 1


def test_load_local_library_absent(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    assert ora.load_local_library("NotThere") is None


# ── enrichment behavior ──────────────────────────────────────────────────────

def test_run_ora_enriches_planted_term():
    # A universe of 1000 genes; one term of 20 genes; the query is 18 of those
    # 20 plus 2 random — strongly over-represented in that term, not in a decoy.
    bg = [f"G{i}" for i in range(1000)]
    hit_term = [f"G{i}" for i in range(20)]
    decoy = [f"G{i}" for i in range(500, 520)]
    gene_sets = {"HIT": hit_term, "DECOY": decoy}
    query = [f"G{i}" for i in range(18)] + ["G900", "G901"]

    res = ora.run_ora(query, gene_sets, bg, padj_max=0.05, top=20)
    terms = {r["term"] for r in res}
    assert "HIT" in terms
    assert "DECOY" not in terms
    hit = next(r for r in res if r["term"] == "HIT")
    assert hit["overlap"] == "18/20"
    assert hit["odds_ratio"] > 1
    assert 0.0 <= hit["padj"] < 0.05


def test_run_ora_no_significant_returns_empty():
    bg = [f"G{i}" for i in range(1000)]
    gene_sets = {"T": [f"G{i}" for i in range(20)]}
    # query disjoint from the term -> no over-representation
    res = ora.run_ora(["G500", "G501", "G502"], gene_sets, bg, padj_max=0.05)
    assert res == []


def test_local_ora_for_databases_reports_missing(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    _write_library(tmp_path, "LIB_A", {"T": ["A", "B", "C"]},
                   manifest={"release": "v1"})
    results, versions, missing = ora.local_ora_for_databases(
        ["A", "B"], {"A": "LIB_A", "B": "LIB_B_absent"},
        ["A", "B", "C", "D", "E"], padj_max=1.0,
    )
    assert "A" in results
    assert versions["A"]["release"] == "v1"
    assert missing == ["B"]


# ── Enrichr opt-in (default off) ─────────────────────────────────────────────

def test_enrichr_is_opt_in(monkeypatch):
    monkeypatch.delenv(ora.ENRICHR_OPT_IN_ENV, raising=False)
    assert ora.enrichr_opt_in() is False
    monkeypatch.setenv(ora.ENRICHR_OPT_IN_ENV, "1")
    assert ora.enrichr_opt_in() is True
    monkeypatch.setenv(ora.ENRICHR_OPT_IN_ENV, "0")
    assert ora.enrichr_opt_in() is False


# ── bulk DE script: local ORA is the default, no network/gseapy needed ───────

def _import_bulk():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "rbd_under_test", "aria/scripts/rna_bulk_de.py"
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write_human_libraries(tmp_path, hit_term, hit_genes):
    """Write the three human libraries _get_gene_sets() asks for, each holding
    the planted hit term so local ORA can recover it."""
    for name in ("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022"):
        _write_library(tmp_path, name, {hit_term: hit_genes,
                                        "DECOY": ["Z1", "Z2", "Z3"]})


def test_bulk_run_pathway_enrichment_defaults_to_local(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    monkeypatch.delenv(ora.ENRICHR_OPT_IN_ENV, raising=False)
    rbd = _import_bulk()

    universe = [f"G{i}" for i in range(1000)]
    hit_genes = [f"G{i}" for i in range(20)]
    _write_human_libraries(tmp_path, "PLANTED_TERM", hit_genes)
    sig = [f"G{i}" for i in range(18)]   # 18/20 of the planted term

    pathways, warnings, meta = rbd._run_pathway_enrichment(
        sig_genes=sig, up_genes=sig, down_genes=[],
        organism="Homo sapiens", output_dir=str(tmp_path),
        symbol_map=None, background_genes=universe,
    )

    assert meta["method"] == "local_hypergeometric"
    assert meta["databases_skipped"] == []
    assert meta["databases_enrichr"] == []
    assert set(meta["databases_local"]) == {"GO_BP", "KEGG", "Reactome"}
    # the exact gene-set release is recorded for methodology.json
    assert meta["gene_set_versions"]["GO_BP"]["library"] == "GO_Biological_Process_2021"
    assert any(r["term"] == "PLANTED_TERM" for r in pathways["GO_BP"])


def test_bulk_run_pathway_enrichment_skips_without_gmt_or_optin(tmp_path, monkeypatch):
    # empty GMT dir + Enrichr not opted in -> honest skip, no fabrication, no net
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path / "empty"))
    monkeypatch.delenv(ora.ENRICHR_OPT_IN_ENV, raising=False)
    rbd = _import_bulk()

    pathways, warnings, meta = rbd._run_pathway_enrichment(
        sig_genes=["A", "B", "C"], up_genes=["A"], down_genes=["B"],
        organism="Homo sapiens", output_dir=str(tmp_path),
        symbol_map=None, background_genes=["A", "B", "C", "D", "E"],
    )

    assert pathways == {}
    assert meta["method"] == "none"
    assert set(meta["databases_skipped"]) == {"GO_BP", "KEGG", "Reactome"}
    assert meta["databases_enrichr"] == []
    assert any("opt-in" in w or "fetch_genesets" in w for w in warnings)


def test_per_cluster_local_ora_success(tmp_path, monkeypatch):
    monkeypatch.setenv(ora.GMT_DIR_ENV, str(tmp_path))
    monkeypatch.delenv(ora.ENRICHR_OPT_IN_ENV, raising=False)
    hit_genes = [f"G{i}" for i in range(20)]
    for name in ("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022"):
        _write_library(tmp_path, name, {"PLANTED": hit_genes, "DECOY": ["Z1", "Z2"]},
                       manifest={"library": name, "release": "2021"})
    from aria.scripts.rna_pathway_per_cluster import rna_pathway_per_cluster

    universe = [f"G{i}" for i in range(1000)]
    records = [{"gene": f"G{i}", "log2fc": 3.0, "padj": 0.001} for i in range(18)]
    out = rna_pathway_per_cluster({
        "de_genes_by_cluster": {"c1": records},
        "organism": "Homo sapiens",
        "background_genes": universe,
        "output_dir": str(tmp_path),
    })

    assert out["status"] == "success"
    assert out["ora_method"] == "local_hypergeometric"
    assert out["gene_set_versions"]["GO_BP"]["release"] == "2021"
    assert out["databases_skipped"] == []
    hits = out["per_cluster"]["c1"]["results"]["GO_BP"]
    assert any(r["term"] == "PLANTED" for r in hits)
