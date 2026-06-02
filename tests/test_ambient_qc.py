"""P1-4: ambient-RNA contamination detector — WARN only, never corrects.

`assess_ambient_contamination` is data-driven: it measures how often each
cluster's top-marker genes appear as a top marker across OTHER clusters. A gene
that is "top" in a large fraction of clusters is non-specific (ambient/soup-like
or housekeeping). When many top-marker slots are filled by such ubiquitous
genes, ARIA warns that ambient contamination may be inflating cross-cluster
signal and recommends the optional decontamination step. No hardcoded gene
lists (ADR-011); no correction is applied.
"""

from aria.utils.ambient_qc import assess_ambient_contamination


def test_clean_when_markers_are_cluster_specific():
    res = assess_ambient_contamination(
        top_markers={
            "A": ["G1", "G2", "G3", "G4"],
            "B": ["G5", "G6", "G7", "G8"],
            "C": ["G9", "G10", "G11", "G12"],
        },
    )
    assert res["status"] == "clean"
    assert res["issues"] == []
    assert res["metrics"]["ubiquitous_genes"] == []


def test_flags_ubiquitous_markers_across_clusters():
    # SOUP1 / SOUP2 are "top" in every cluster -> ambient-like leakage.
    res = assess_ambient_contamination(
        top_markers={
            "A": ["SOUP1", "SOUP2", "A3", "A4"],
            "B": ["SOUP1", "SOUP2", "B3", "B4"],
            "C": ["SOUP1", "SOUP2", "C3", "C4"],
            "D": ["SOUP1", "SOUP2", "D3", "D4"],
        },
        ubiquity_fraction=0.5,
    )
    checks = {i["check"] for i in res["issues"]}
    assert "possible_ambient_contamination" in checks
    assert set(res["metrics"]["ubiquitous_genes"]) >= {"SOUP1", "SOUP2"}
    # Detector only: it recommends decontamination, it does not claim to correct.
    msg = next(i for i in res["issues"]
               if i["check"] == "possible_ambient_contamination")
    assert "decont" in msg["recommendation"].lower() or \
           "ambient" in msg["recommendation"].lower()


def test_unverified_when_fewer_than_two_clusters():
    res = assess_ambient_contamination(top_markers={"A": ["G1", "G2", "G3"]})
    assert res["status"] in ("unverified", "clean")
    assert res["issues"] == []


def test_graceful_on_empty_markers():
    res = assess_ambient_contamination(top_markers={})
    assert res["status"] in ("unverified", "clean")
    assert res["issues"] == []

    res2 = assess_ambient_contamination(top_markers=None)
    assert res2["status"] in ("unverified", "clean")
    assert res2["issues"] == []


def test_empty_marker_lists_do_not_count_as_ubiquitous():
    res = assess_ambient_contamination(
        top_markers={"A": [], "B": [], "C": []},
    )
    assert res["metrics"]["ubiquitous_genes"] == []
    assert res["status"] in ("unverified", "clean")
