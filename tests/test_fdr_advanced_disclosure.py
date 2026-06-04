"""P1-2 residual honest closure (ADR-027): IHW + s-values are NOT implemented,
and that is disclosed transparently rather than faked.

There is no validated pure-Python IHW (Ignatiadis-Huber) estimator and pydeseq2
0.5.4 exposes no s-values; a hand-rolled covariate-weighted BH could silently
break FDR control (ADR-002), so ARIA ships neither and says so. The primary FDR
stays the pre-registered Benjamini-Hochberg family.
"""

from aria.utils.stats import fdr_advanced_methods_disclosure


def test_disclosure_marks_ihw_and_svalues_unavailable_with_reasons():
    d = fdr_advanced_methods_disclosure()
    assert d["ihw"]["available"] is False
    assert d["ihw"]["status"] == "not_implemented"
    assert d["ihw"]["reason"]            # a non-empty scientific reason
    assert d["s_values"]["available"] is False
    assert d["s_values"]["status"] == "not_available"
    assert "pydeseq2" in d["s_values"]["reason"].lower()


def test_disclosure_keeps_bh_as_primary_and_never_claims_ihw():
    d = fdr_advanced_methods_disclosure()
    assert "benjamini" in d["primary_method"].lower()
    # No path may report IHW/s-values as if they ran.
    assert d["ihw"]["status"] != "applied"
    assert d["s_values"]["status"] != "applied"


def test_disclosure_is_deterministic():
    assert fdr_advanced_methods_disclosure() == fdr_advanced_methods_disclosure()
