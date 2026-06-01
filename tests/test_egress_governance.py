"""W-PRIV (P1-7 / P1-8): air-gapped mode must govern ALL network egress.

`ARIA_AIR_GAPPED` previously only blocked the LLM layer. The ORA path
(`gseapy.enrichr`) still shipped DE gene lists to an external server, and the
GEO/SRA connector still fetched from NCBI -- so a "sensitive / air-gapped" run
leaked. A shared egress gate now refuses Enrichr ORA and connector fetches under
air-gapped mode.

P1-7 then made local hypergeometric ORA the default, so air-gapped no longer
means "no ORA": ORA runs locally and offline against versioned GMTs, and only
the *Enrichr* (network) path is refused -- even when explicitly opted in. These
tests assert that stronger invariant: under air-gapped no gene list reaches the
network, while local ORA still works.
"""

import json

import pytest

from aria.utils import ora as _ora
from aria.utils.privacy import (
    EgressBlocked,
    air_gapped_enabled,
    assert_egress_allowed,
    egress_allowed,
)


def _write_human_gmts(gmt_dir, term, genes):
    """Write the three human libraries the ORA engine asks for, locally."""
    for name in ("GO_Biological_Process_2021", "KEGG_2021_Human", "Reactome_2022"):
        base = gmt_dir / name
        base.mkdir(parents=True)
        (base / f"{name}.gmt").write_text(
            "\t".join([term, "", *genes]) + "\n", encoding="utf-8"
        )
        (base / "manifest.json").write_text(
            json.dumps({"library": name, "release": "2021"}), encoding="utf-8"
        )


def test_egress_gate_tracks_air_gapped(monkeypatch):
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    assert egress_allowed() is True
    assert_egress_allowed("enrichr")  # does not raise

    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    assert air_gapped_enabled() is True
    assert egress_allowed() is False
    with pytest.raises(EgressBlocked):
        assert_egress_allowed("enrichr")


def test_bulk_ora_refuses_enrichr_egress_when_air_gapped(monkeypatch, tmp_path):
    # Even with Enrichr EXPLICITLY opted in, air-gapped must refuse the network
    # ORA. With no local GMTs the databases are skipped (no fabrication, no
    # egress) and the caveat names air-gapped as the reason.
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    monkeypatch.setenv("ARIA_ALLOW_ENRICHR", "1")
    monkeypatch.setenv("ARIA_GMT_DIR", str(tmp_path / "empty"))
    from aria.scripts.rna_bulk_de import _run_pathway_enrichment

    pathways, warnings, meta = _run_pathway_enrichment(
        sig_genes=["TP53", "STAT1"], up_genes=["TP53"], down_genes=["STAT1"],
        organism="Homo sapiens", output_dir=str(tmp_path),
        background_genes=["TP53", "STAT1", "EGFR"],
    )
    assert pathways == {}
    assert meta["databases_enrichr"] == []         # no network ORA happened
    assert any("air-gapped" in w.lower() or "ARIA_AIR_GAPPED" in w
               for w in warnings)


def test_bulk_local_ora_still_runs_when_air_gapped(monkeypatch, tmp_path):
    # P1-7: local hypergeometric ORA is offline, so it runs even air-gapped.
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    monkeypatch.setenv("ARIA_GMT_DIR", str(tmp_path))
    _write_human_gmts(tmp_path, "PLANTED", [f"G{i}" for i in range(20)])
    from aria.scripts.rna_bulk_de import _run_pathway_enrichment

    universe = [f"G{i}" for i in range(1000)]
    sig = [f"G{i}" for i in range(18)]
    pathways, warnings, meta = _run_pathway_enrichment(
        sig_genes=sig, up_genes=sig, down_genes=[],
        organism="Homo sapiens", output_dir=str(tmp_path),
        background_genes=universe,
    )
    assert meta["method"] == "local_hypergeometric"
    assert any(r["term"] == "PLANTED" for r in pathways["GO_BP"])


def test_per_cluster_ora_refuses_enrichr_when_air_gapped(monkeypatch, tmp_path):
    # Opted into Enrichr + air-gapped + no local GMT -> honest skip, no egress.
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    monkeypatch.setenv("ARIA_ALLOW_ENRICHR", "1")
    monkeypatch.setenv("ARIA_GMT_DIR", str(tmp_path / "empty"))
    from aria.scripts.rna_pathway_per_cluster import rna_pathway_per_cluster

    out = rna_pathway_per_cluster({
        "de_genes_by_cluster": {"c1": [{"gene": "TP53", "log2fc": 2.0, "padj": 0.01}]},
        "organism": "Homo sapiens",
        "output_dir": str(tmp_path),
    })
    assert out["status"] == "skipped"
    assert "air_gapped" in out["reason"]


def test_geo_fetch_refuses_when_air_gapped(monkeypatch, tmp_path):
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    from aria.connectors.geo_connector import GEOConnector

    conn = GEOConnector(cache_dir=str(tmp_path))
    with pytest.raises(EgressBlocked):
        conn.fetch("GSE12345")
