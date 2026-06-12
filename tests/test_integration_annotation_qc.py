"""X8 integration QC red-flags + X9 annotation-coherence checks."""

from aria.utils.integration_qc import assess_integration_quality
from aria.utils.annotation_qc import assess_annotation_coherence


# ── X8 integration QC ──────────────────────────────────────────────────────

def test_integration_clean_when_mixed_and_structure_preserved():
    res = assess_integration_quality(
        silhouette_before=0.30, silhouette_after=0.02, cluster_silhouette=0.55)
    assert res["status"] == "clean"
    assert res["issues"] == []


def test_integration_flags_residual_batch_effect():
    res = assess_integration_quality(
        silhouette_before=0.40, silhouette_after=0.35, cluster_silhouette=0.5)
    checks = {i["check"] for i in res["issues"]}
    assert "residual_batch_effect" in checks


def test_integration_flags_worsened_mixing():
    res = assess_integration_quality(
        silhouette_before=0.05, silhouette_after=0.20, cluster_silhouette=0.4)
    checks = {i["check"] for i in res["issues"]}
    assert "integration_worsened_mixing" in checks


def test_integration_flags_overcorrection_on_negative_cluster_silhouette():
    res = assess_integration_quality(
        silhouette_before=0.40, silhouette_after=0.01, cluster_silhouette=-0.05)
    checks = {i["check"] for i in res["issues"]}
    assert "possible_overcorrection" in checks
    msg = next(i["message"] for i in res["issues"]
               if i["check"] == "possible_overcorrection")
    assert "overcorrection" in msg.lower()


def test_integration_classic_overcorrection_signature_is_blocking():
    """N-INT1: strong batch mixing AND collapsed cluster structure is the classic
    overcorrection signature and must escalate from a soft warning to blocking."""
    res = assess_integration_quality(
        silhouette_before=0.40, silhouette_after=0.01, cluster_silhouette=-0.05)
    oc = next(i for i in res["issues"] if i["check"] == "possible_overcorrection")
    assert oc["severity"] == "blocking"


def test_integration_overcorrection_without_strong_mixing_stays_warning():
    # Structure collapsed but mixing barely changed -> not the classic signature.
    res = assess_integration_quality(
        silhouette_before=0.05, silhouette_after=0.04, cluster_silhouette=-0.05)
    oc = next(i for i in res["issues"] if i["check"] == "possible_overcorrection")
    assert oc["severity"] == "warning"


def test_integration_handles_missing_metrics():
    res = assess_integration_quality(None, None, None)
    assert res["status"] == "clean"


# ── X9 annotation coherence ─────────────────────────────────────────────────

def test_annotation_unverified_when_reused_without_markers():
    res = assess_annotation_coherence(
        top_markers={"A": [], "B": []}, cluster_sizes={"A": 100, "B": 80},
        reused=True, markers_verified=False)
    assert res["status"] == "unverified"
    assert any(i["check"] == "annotation_unverified" for i in res["issues"])


def test_annotation_flags_label_without_distinct_markers():
    res = assess_annotation_coherence(
        top_markers={
            "TypeA": ["G1", "G2", "G3", "G4"],
            "TypeB": ["G5", "G6", "G7"],
            "Noise": [],          # no distinct signature
        },
        cluster_sizes={"TypeA": 500, "TypeB": 300, "Noise": 40},
        reused=True, markers_verified=True, min_markers=3)
    assert res["status"] == "warnings"
    issue = next(i for i in res["issues"]
                 if i["check"] == "labels_without_distinct_markers")
    assert "Noise" in issue["message"]
    assert res["per_label"]["Noise"]["distinct"] is False
    assert res["per_label"]["TypeA"]["distinct"] is True


def test_annotation_clean_when_all_labels_distinct():
    res = assess_annotation_coherence(
        top_markers={"A": ["G1", "G2", "G3"], "B": ["G4", "G5", "G6", "G7"]},
        cluster_sizes={"A": 100, "B": 120},
        reused=True, markers_verified=True, min_markers=3)
    assert res["status"] == "clean"
    assert res["issues"] == []
