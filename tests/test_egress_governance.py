"""W-PRIV (P1-7 / P1-8): air-gapped mode must govern ALL network egress.

`ARIA_AIR_GAPPED` previously only blocked the LLM layer. The ORA path
(`gseapy.enrichr`) still shipped DE gene lists to an external server, and the
GEO/SRA connector still fetched from NCBI -- so a "sensitive / air-gapped" run
leaked. A shared egress gate now refuses Enrichr ORA and connector fetches under
air-gapped mode, degrading honestly (ORA skipped with a caveat; connector
refuses) instead of leaking.
"""

import pytest

from aria.utils.privacy import (
    EgressBlocked,
    air_gapped_enabled,
    assert_egress_allowed,
    egress_allowed,
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


def test_bulk_ora_skips_enrichr_when_air_gapped(monkeypatch, tmp_path):
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    from aria.scripts.rna_bulk_de import _run_pathway_enrichment

    pathways, warnings = _run_pathway_enrichment(
        sig_genes=["TP53", "STAT1"], up_genes=["TP53"], down_genes=["STAT1"],
        organism="Homo sapiens", output_dir=str(tmp_path),
    )
    assert pathways == {}
    assert any("air-gapped" in w.lower() or "ARIA_AIR_GAPPED" in w
               for w in warnings)


def test_per_cluster_ora_skips_when_air_gapped(monkeypatch, tmp_path):
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    from aria.scripts.rna_pathway_per_cluster import rna_pathway_per_cluster

    out = rna_pathway_per_cluster({
        "de_genes_by_cluster": {"c1": ["TP53", "STAT1"]},
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
