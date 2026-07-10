"""Preprint-readiness audit — Q2 (interim mitigation of blocker B7).

TOBIAS pools every cell/BAM of a condition into a single pseudobulk, so the
per-site `<A>_<B>_pvalue` cannot support an FDR-controlled significance claim.
Until the full B7 fix (per-replicate model + BH multiplicity + null controls),
`summarize_bindetect` must present a DESCRIPTIVE candidate ranking with an
explicit `ranking_basis` disclosure — never `n_significant` — and the narrator
must not call these motifs "significant".

Tracker: memory/audit/ARIA_PLAN_AUDITORIA_preprint_journal_2026-07-09.md
"""
from pathlib import Path

from aria.scripts.chromatin_footprint_tobias import summarize_bindetect
from aria.agents.narrative.narrators.chromatin import ChromatinNarrator


def _write_bindetect(tmp_path: Path, a: str, b: str) -> str:
    """Minimal TOBIAS bindetect_results.txt with the columns summarize reads."""
    header = ["name", f"{a}_{b}_change", f"{a}_{b}_pvalue", "total_tfbs"]
    rows = [
        ["KLF1", "1.20", "1e-100", "9000"],   # ranked toward a
        ["GATA1", "0.80", "1e-50", "4000"],   # ranked toward a
        ["IRF8", "-0.90", "1e-40", "6000"],   # ranked toward b
        ["NOISE", "0.05", "1e-9", "5"],       # dropped: too few sites
        ["FLAT", "0.10", "0.9", "8000"],      # dropped: p above threshold
    ]
    p = tmp_path / "bindetect_results.txt"
    p.write_text("\t".join(header) + "\n"
                 + "\n".join("\t".join(r) for r in rows) + "\n")
    return str(p)


def test_summarize_reports_descriptive_ranking_not_significance(tmp_path):
    summary = summarize_bindetect(_write_bindetect(tmp_path, "K562", "GM12878"),
                                  "K562", "GM12878")
    assert summary["parsed"] is True
    # No significance count leaks; the count is a descriptive candidate ranking.
    assert "n_significant" not in summary
    assert summary["n_ranked_candidates"] == 3      # 2 dropped (sites, p)
    assert summary["n_motifs_tested"] == 5
    # Explicit honest disclosure of what the ranking is (and is not).
    basis = summary["ranking_basis"]
    assert basis["fdr_controlled"] is False
    assert basis["replicate_inference"] is False
    assert basis["pseudobulk"] is True
    assert "not FDR-controlled significance".lower() in basis["note"].lower()


def test_narrator_block_does_not_claim_significance(tmp_path):
    footprinting = {
        "ran": True, "data_type": "scATAC",
        "group_a": "Monocyte", "group_b": "T_cell",
        "differential_summary": summarize_bindetect(
            _write_bindetect(tmp_path, "Monocyte", "T_cell"),
            "Monocyte", "T_cell"),
    }
    blocks = ChromatinNarrator().collect(
        "chromatin_agent", {"findings": {"footprinting": footprinting}})
    block = next(b for b in blocks
                 if b.id == "chromatin.differential_tf_footprinting")
    low = block.claim.lower()
    assert "descriptively ranked" in low
    assert "not fdr-controlled significance" in low
    # The evidence label no longer sells a significance count.
    assert not any(e.label.startswith("Significant") for e in block.evidence)
    assert any(e.label.startswith("Top-ranked") for e in block.evidence)
    # The caveat must warn about pseudoreplication + no FDR.
    caveat_txt = " ".join(c.text.lower() for c in block.caveats)
    assert "pseudobulk" in caveat_txt
    assert "replication" in caveat_txt
    assert "not fdr-controlled" in caveat_txt
