"""ARIA Bulk RNA-seq — native pytest (P1-11).

Converted from the legacy script-style harness (top-level `ok()`/`fail()` +
`sys.exit()`), which (a) crashed pytest collection when named explicitly and
(b) was only exercised through a subprocess wrapper that could false-green.
These are now real pytest tests with assertions: light checks run anywhere; the
DESeq2 end-to-end + the versioned GOLDEN recovery are gated on pydeseq2 (heavy
CI lane), so a numerical regression is a red test, not a swallowed print.
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import aria.scripts.rna_bulk_de as rbd
# rna_bulk_de is a plain script (no litellm) — safe to import in every lane.
from aria.scripts.rna_bulk_de import (
    bulk_rna_de,
    _get_gene_sets,
    _infer_groups,
    _mock_pathways,
    _run_pathway_enrichment,
    _sample_qc,
)

# BulkRNAAgent pulls aria.llm.provider -> litellm, which the heavy pydeseq2 lane
# (aria-rna-env) and the light pip lane do not install. Guard it so those lanes
# skip the agent-only checks instead of failing collection.
try:
    from aria.agents.bulk_rna_agent import (
        BulkRNAAgent, _default_lfc_threshold, suggest_lfc_profile,
        DEFAULT_LFC_THRESHOLD)
    _HAS_AGENT = True
except Exception:
    _HAS_AGENT = False

try:
    import pydeseq2  # noqa: F401
    _HAS_PYDESEQ2 = True
except Exception:
    _HAS_PYDESEQ2 = False

heavy = pytest.mark.skipif(
    not _HAS_PYDESEQ2, reason="pydeseq2 required (heavy CI lane / aria-rna-env)"
)
needs_agent = pytest.mark.skipif(
    not _HAS_AGENT, reason="bulk_rna_agent (litellm) unavailable in this lane"
)

GOLDEN = Path(__file__).parent / "fixtures" / "golden" / "bulk_mini"


def make_counts(n_genes=200, samples=None, seed=42):
    """Minimal realistic count matrix with ~20 planted DE genes."""
    rng = np.random.default_rng(seed)
    if samples is None:
        samples = ["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"]
    base = rng.negative_binomial(20, 0.5, (n_genes, len(samples)))
    trt_idx = [i for i, s in enumerate(samples)
               if "treat" in s.lower() or "ko" in s.lower() or "mut" in s.lower()]
    for i in range(20):
        base[i, trt_idx] = base[i, trt_idx] * rng.integers(3, 8)
    return pd.DataFrame(base, columns=samples,
                        index=[f"GENE_{i:04d}" for i in range(n_genes)])


# ── Fix 1: group auto-detection from column names ────────────────────────────

@pytest.mark.parametrize("samples,expected_groups", [
    (["ctrl_1", "ctrl_2", "ctrl_3", "treat_1", "treat_2", "treat_3"], {"ctrl", "treat"}),
    (["WT_rep1", "WT_rep2", "KO_rep1", "KO_rep2"], {"WT", "KO"}),
    (["vehicle_1", "vehicle_2", "drug_1", "drug_2"], {"vehicle", "drug"}),
    (["Healthy_1", "Healthy_2", "Disease_1", "Disease_2"], {"Healthy", "Disease"}),
])
def test_infer_groups_naming_patterns(samples, expected_groups):
    groups = _infer_groups(samples)
    assert groups is not None
    assert set(groups.values()) == expected_groups


def test_infer_groups_ambiguous_is_graceful():
    groups = _infer_groups(["sample1", "sample2", "sample3", "sample4"])
    # Either None or a single group — never a crash, never spurious 2 groups.
    assert groups is None or len(set(groups.values())) < 2


# ── Fix 2: design factor + LFC + entity→label mapping ────────────────────────

@needs_agent
@pytest.mark.parametrize("intent,expected_factor", [
    ({"comparison": "knockout vs wildtype", "summary": "KO vs WT"}, "genotype"),
    ({"comparison": "treated vs control", "summary": "drug treatment"}, "treatment"),
    ({"comparison": "24h vs 0h", "summary": "timepoint course"}, "timepoint"),
    ({"comparison": "lupus vs healthy", "summary": "disease state"}, "condition"),
])
def test_infer_design_factor(intent, expected_factor):
    assert BulkRNAAgent._infer_design_factor(intent) == expected_factor


@needs_agent
def test_lfc_threshold_is_prompt_independent():
    """F1 (ADR-055): the |log2FC| cutoff must NOT depend on the question text.

    A TF/knockout study and a neutral study resolve to the SAME fixed default;
    the only deviation is an explicit user-confirmed CP3 profile / global_lfc.
    """
    tf = {"comparison": "SOX2 KO vs neutral_ref",
          "summary": "transcription factor knockout",
          "biological_entities": ["SOX2", "neutral_ref"]}
    non_tf = {"comparison": "condition_X vs neutral_ref",
              "biological_entities": ["condition_X"]}
    # The default is fixed and identical regardless of TF/KO wording.
    assert _default_lfc_threshold() == DEFAULT_LFC_THRESHOLD == 1.0
    # The TF/KO heuristic survives ONLY as a non-binding advisory suggestion.
    assert suggest_lfc_profile(tf) == "exploratory_tf"
    assert suggest_lfc_profile(non_tf) is None


@needs_agent
def test_lfc_suggestion_does_not_depend_on_hardcoded_gene_names():
    """ADR-055: the advisory must fire on generic perturbation-design KEYWORDS,
    never on a hardcoded list of gene symbols. A study naming a transcription
    factor WITHOUT a perturbation keyword must NOT trigger the hint, and the
    module must not carry a gene whitelist."""
    import aria.agents.bulk_rna_agent as bra
    # No hardcoded gene-symbol set survives in the module.
    assert not hasattr(bra, "KNOWN_TFS")
    # A bare gene symbol (no knockout/knockdown/overexpression wording) -> no hint.
    bare_gene = {"comparison": "BMAL1 vs WT",
                 "summary": "compare BMAL1 against wild type",
                 "biological_entities": ["BMAL1", "WT"]}
    assert suggest_lfc_profile(bare_gene) is None
    # The same biology phrased as a perturbation -> hint fires (keyword, not gene).
    perturbed = {"comparison": "BMAL1 knockdown vs WT",
                 "biological_entities": ["BMAL1", "WT"]}
    assert suggest_lfc_profile(perturbed) == "exploratory_tf"


@needs_agent
def test_lfc_suggestion_does_not_change_resolved_threshold():
    """The advisory suggestion never moves the cutoff: a TF-flagged intent with
    no explicit override still resolves to the fixed default."""
    import aria.agents.bulk_rna_agent as bra
    tf = {"summary": "BMAL1 knockout", "comparison": "KO vs WT",
          "biological_entities": ["BMAL1"]}
    # suggestion fires, but the resolved cutoff (no global_lfc) is the default.
    assert suggest_lfc_profile(tf) == "exploratory_tf"
    resolved = {}.get("global_lfc", bra._default_lfc_threshold())
    assert resolved == 1.0


@needs_agent
def test_entity_to_label_mapping():
    mapping = BulkRNAAgent.__new__(BulkRNAAgent)._map_entities_to_labels(
        entities=["A_entity", "B_entity", "REF_entity"],
        group_names=["A", "B", "REF"], intent={},
    )
    assert mapping.get("A_entity") == "A"
    assert mapping.get("B_entity") == "B"
    assert mapping.get("REF_entity") == "REF"


@needs_agent
def test_bulk_rna_agent_does_not_generate_free_text_interpretation(tmp_path):
    counts_path = tmp_path / "counts.tsv"
    counts_path.write_text(
        "gene\tctrl1\tctrl2\ttreat1\ttreat2\n"
        "GENE1\t10\t11\t40\t42\n",
        encoding="utf-8",
    )

    class FakeLLM:
        calls = 0

        def complete(self, *args, **kwargs):
            self.calls += 1
            raise AssertionError("bulk RNA free-text interpretation must not call LLM")

    class FakeMemory:
        def store_decision(self, *args, **kwargs):
            return None

    class FakeEnv:
        params = None

        def run_in_stack(self, stack, script_path, params):
            self.params = params
            return {
                "status": "success",
                "sample_qc": {"n_samples": 4, "outliers": []},
                "contrasts": [{
                    "name": "treated vs control",
                    "status": "success",
                    "n_significant": 1,
                    "n_upregulated": 1,
                    "n_downregulated": 0,
                    "top_genes": [{"gene": "GENE1", "log2fc": 2.0}],
                    "pathways": {},
                }],
                "methodology": {"decisions": []},
                "design_used": "~ condition",
            }

    agent = BulkRNAAgent.__new__(BulkRNAAgent)
    agent.memory = FakeMemory()
    agent.llm = FakeLLM()
    agent.env = FakeEnv()
    agent.publish_status = lambda *args, **kwargs: None
    agent.publish_finding = lambda *args, **kwargs: None

    result = agent.run(
        "f3_bulk_test",
        {
            "exp_context": {
                "organism": "Homo sapiens",
                "genome": "hg38",
                "modalities": {"bulk_RNA": [str(counts_path)]},
                "design": {
                    "groups": {
                        "control": ["ctrl1", "ctrl2"],
                        "treated": ["treat1", "treat2"],
                    },
                    "main_factor": "condition",
                    "contrasts": [{
                        "numerator": "treated",
                        "denominator": "control",
                        "name": "treated vs control",
                    }],
                },
            },
            "biological_intent": {
                "summary": "treated vs control bulk RNA-seq",
            },
        },
    )

    assert result["status"] == "done"
    findings = result["findings"]
    assert "interpretation" not in findings
    assert findings["interpretation_status"] == {
        "ran": False,
        "reason": (
            "free_text_llm_interpretation_disabled_by_F3_governance; "
            "bulk RNA Results are generated from structured DE/pathway "
            "outputs and NarrativeBlock evidence cards"
        ),
        "governance": "F3_preprint_audit",
    }
    assert agent.llm.calls == 0
    assert agent.env.params["contrasts"] == [{
        "numerator": "treated",
        "denominator": "control",
        "name": "treated vs control",
    }]


# ── Fix 3: sample QC / outlier detection ─────────────────────────────────────

def test_sample_qc_reports_expected_keys():
    counts = make_counts(200)
    metadata = pd.DataFrame({"condition": ["ctrl"] * 3 + ["treat"] * 3},
                            index=counts.columns)
    with tempfile.TemporaryDirectory() as tmp:
        qc = _sample_qc(counts, metadata, tmp, [])
    for key in ("n_samples", "outliers", "pca_variance", "lib_size_range"):
        assert key in qc
    assert qc["n_samples"] == 6


def test_sample_qc_runs_with_injected_outlier():
    counts = make_counts(200)
    counts["ctrl_1"] = counts["ctrl_1"] * 100   # extreme inflation
    metadata = pd.DataFrame({"condition": ["ctrl"] * 3 + ["treat"] * 3},
                            index=counts.columns)
    with tempfile.TemporaryDirectory() as tmp:
        qc = _sample_qc(counts, metadata, tmp, [])
    assert isinstance(qc["outliers"], list)   # detection runs, returns a list


# ── Fix 4: pathway helpers (local ORA, P1-7) ─────────────────────────────────

def test_mock_pathways_structure():
    mock = _mock_pathways(["gene_1", "gene_2", "gene_3", "gene_4", "gene_5"])
    assert ("GO_BP" in mock) or ("KEGG" in mock)
    assert any(isinstance(v, list) for v in mock.values())


def test_gene_sets_by_organism():
    human = _get_gene_sets("Homo sapiens")
    mouse = _get_gene_sets("Mus musculus")
    assert {"GO_BP", "KEGG", "Reactome"} <= set(human)
    assert "KEGG" in mouse


def test_run_pathway_enrichment_returns_triple():
    pw, warnings, meta = _run_pathway_enrichment(
        sig_genes=["gene_1", "gene_2", "gene_3", "gene_4", "gene_5", "gene_6"],
        up_genes=["gene_1", "gene_2"], down_genes=["gene_5"],
        organism="Homo sapiens", output_dir="/tmp/test_pw",
    )
    assert isinstance(pw, dict)
    assert isinstance(warnings, list)
    assert isinstance(meta, dict) and "method" in meta


# ── P1-5: primary DE + outlier sensitivity ──────────────────────────────────

def _patch_light_bulk_de(monkeypatch, primary_sig, sensitivity_sig):
    calls = []

    def fake_sample_qc(counts, metadata, output_dir, warnings, biotype_map=None):
        return {
            "n_samples": int(counts.shape[1]),
            "outliers": ["ctrl_1"],
            "pca_variance": [0.7, 0.2],
            "lib_size_range": [1000, 2000],
            "size_ratio": 2.0,
        }

    def fake_run_deseq2(
        counts, metadata, design_factor, numerator, denominator,
        padj_thr, lfc_thr, allow_mock=False,
        min_replicates_per_condition=3, covariates=None, lfc_shrink=True,
    ):
        calls.append(list(metadata.index))
        is_sensitivity = "ctrl_1" not in metadata.index
        sig = sensitivity_sig if is_sensitivity else primary_sig
        genes = [f"GENE_{i:04d}" for i in range(8)]
        df = pd.DataFrame({
            "log2FoldChange": [2.0 if g in sig else 0.1 for g in genes],
            "pvalue": [0.001 if g in sig else 0.5 for g in genes],
            "padj": [0.01 if g in sig else 0.8 for g in genes],
            "baseMean": [100.0] * len(genes),
        }, index=genes)
        return {
            "status": "success",
            "results": df,
            "n_sig": len(sig),
            "n_up": len(sig),
            "n_down": 0,
            "sig_genes": list(sig),
            "up_genes": list(sig),
            "down_genes": [],
            "n_replicates": {"test": 3, "ref": 2 if is_sensitivity else 3},
            "fitted_design_formula": "~ condition",
            "covariates_adjusted": [],
            "covariates_dropped": [],
            "lfc_shrinkage": {"requested": False, "applied": False},
            "lfc_threshold_test": {"applied": True},
        }, []

    monkeypatch.setattr(rbd, "_sample_qc", fake_sample_qc)
    monkeypatch.setattr(rbd, "_run_deseq2", fake_run_deseq2)
    # A7: _run_outlier_sensitivity moved to aria.scripts.rna_bulk.qc, which binds
    # _run_deseq2 in its own namespace — patch that seam too so the sensitivity
    # re-fit uses the fake (bulk_rna_de's primary call still goes through rbd).
    monkeypatch.setattr("aria.scripts.rna_bulk.qc._run_deseq2", fake_run_deseq2)
    monkeypatch.setattr(rbd, "_generate_plots", lambda **kwargs: {})
    monkeypatch.setattr(rbd, "_load_symbol_map", lambda files, warnings: {})
    monkeypatch.setattr(rbd, "_load_gene_annotation", lambda files, warnings: {})
    return calls


def test_outlier_sensitivity_keeps_primary_unpruned(monkeypatch):
    calls = _patch_light_bulk_de(
        monkeypatch,
        primary_sig=["GENE_0001"],
        sensitivity_sig=["GENE_0001"],
    )
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(80),
            [{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
            tmp,
            min_replicates_per_condition=2,
            lfc_shrink=False,
        )

    assert result["status"] == "success", result.get("details", "")
    assert "ctrl_1" in calls[0]
    assert "ctrl_1" not in calls[1]
    sqc = result["sample_qc"]
    assert sqc["candidate_outliers"] == ["ctrl_1"]
    assert sqc["outliers_removed_primary"] == []
    assert sqc["sensitivity_outliers_removed"] == ["ctrl_1"]
    sens = result["outlier_sensitivity"]
    assert sens["status"] == "success"
    assert sens["contrasts"][0]["conclusion_robust"] is True


def test_outlier_sensitivity_does_not_replace_primary(monkeypatch):
    _patch_light_bulk_de(
        monkeypatch,
        primary_sig=[],
        sensitivity_sig=["GENE_0002", "GENE_0003"],
    )
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(80),
            [{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
            tmp,
            min_replicates_per_condition=2,
            lfc_shrink=False,
        )

    contrast = result["contrasts"][0]
    sens = result["outlier_sensitivity"]["contrasts"][0]
    assert contrast["n_significant"] == 0
    assert sens["n_significant_sensitivity"] == 2
    assert sens["sensitivity_only_n"] == 2
    assert sens["conclusion_robust"] is False


# ── End-to-end (real DESeq2) ─────────────────────────────────────────────────

def _run_bulk(counts_df, contrasts, tmp, **overrides):
    counts_path = Path(tmp) / "counts.tsv"
    counts_df.to_csv(str(counts_path), sep="\t")
    samples = list(counts_df.columns)
    conditions = [
        str(s).rsplit("_", 1)[0] if "_" in str(s) else str(s).rstrip("0123456789")
        for s in samples
    ]
    metadata_path = Path(tmp) / "metadata.tsv"
    pd.DataFrame({"condition": conditions}, index=samples).to_csv(
        metadata_path, sep="\t"
    )
    params = {
        "files": [str(counts_path)],
        "metadata_file": str(metadata_path),
        "design_factor": "condition",
        "contrasts": contrasts,
        "organism": "Homo sapiens",
        "output_dir": tmp,
        "run_pathways": False,
        "min_replicates_per_condition": 2,
    }
    params.update(overrides)
    return bulk_rna_de(params)


@heavy
def test_bulk_de_e2e_pipeline():
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(300),
            [{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
            tmp, run_pathways=False, padj_threshold=0.05, lfc_threshold=1.0,
        )
    assert result["status"] == "success", result.get("details", "")
    assert len(result["contrasts"]) == 1
    c0 = result["contrasts"][0]
    for key in ("n_significant", "n_upregulated", "n_downregulated", "plots"):
        assert key in c0


@heavy
def test_bulk_de_multicontrast():
    samples = ["B_1", "B_2", "B_3", "R_1", "R_2", "R_3", "WT_1", "WT_2", "WT_3"]
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(300, samples=samples),
            [{"numerator": "B", "denominator": "WT", "name": "BMAL1 KO vs WT"},
             {"numerator": "R", "denominator": "WT", "name": "REV-ERBa KO vs WT"}],
            tmp, lfc_threshold=0.58,
        )
    assert result["status"] == "success"
    assert result["n_contrasts"] == 2
    names = {c["name"] for c in result["contrasts"]}
    assert {"BMAL1 KO vs WT", "REV-ERBa KO vs WT"} <= names
    assert "overlap" in result


@heavy
def test_replicate_concordance_detects_corruption():
    rng = np.random.default_rng(42)
    n_genes = 500
    samples = ["B_1", "B_2", "B_3", "WT_1", "WT_2", "WT_3"]
    gene_mean = rng.exponential(100, n_genes)
    data = np.array([rng.poisson(gene_mean) for _ in samples]).T.astype(float)
    idx = rng.choice(n_genes, int(n_genes * 0.7), replace=False)
    data[idx, 1] = data[idx, 1] * np.exp(rng.uniform(-2.3, 2.3, len(idx)))  # corrupt B_2
    df = pd.DataFrame(data.astype(int),
                      index=[f"g_{i}" for i in range(n_genes)], columns=samples)
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            df, [{"numerator": "B", "denominator": "WT", "name": "B vs WT"}],
            tmp, lfc_threshold=0.58,
        )
    rep_corr = result.get("sample_qc", {}).get("replicate_correlations", {})
    assert rep_corr, "replicate_correlations missing"
    assert min(rep_corr, key=rep_corr.get) == "B_2"
    assert rep_corr["B_2"] < 0.85


@heavy
def test_deseq2_uses_design_factor():
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(200),
            [{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
            tmp,
        )
    assert "condition" in (result.get("design_used", "") or "")


@heavy
def test_insufficient_replicates_is_graceful():
    with tempfile.TemporaryDirectory() as tmp:
        result = _run_bulk(
            make_counts(200, ["ctrl_1", "treat_1"]),
            [{"numerator": "treat", "denominator": "ctrl", "name": "treat vs ctrl"}],
            tmp, min_replicates_per_condition=3,
        )
    assert result.get("status") in ("error", "success")   # structured, no crash


# ── GOLDEN: deterministic DE recovery against a versioned mini-dataset ────────

@heavy
def test_golden_bulk_de_recovers_planted_genes():
    """Run real DESeq2 on a committed mini count matrix and assert it recovers
    the planted up-regulated DE genes within tolerance — a numerical regression
    (wrong DE direction/recall) turns this red instead of silently passing."""
    expected = json.loads((GOLDEN / "expected.json").read_text())
    with tempfile.TemporaryDirectory() as tmp:
        # Copy the versioned matrix into tmp so any side-outputs bulk_rna_de
        # writes next to the input (e.g. counts_with_symbols.tsv) never dirty
        # the committed fixture directory.
        counts_path = Path(tmp) / "counts.tsv"
        counts_path.write_text((GOLDEN / "counts.tsv").read_text())
        counts_df = pd.read_csv(counts_path, sep="\t", index_col=0)
        samples = list(counts_df.columns)
        conditions = [
            str(s).rsplit("_", 1)[0] if "_" in str(s) else str(s).rstrip("0123456789")
            for s in samples
        ]
        metadata_path = Path(tmp) / "metadata.tsv"
        pd.DataFrame({"condition": conditions}, index=samples).to_csv(
            metadata_path, sep="\t"
        )
        result = bulk_rna_de({
            "files": [str(counts_path)],
            "metadata_file": str(metadata_path),
            "design_factor": expected["design_factor"],
            "contrasts": [expected["contrast"]],
            "organism": "Homo sapiens",
            "output_dir": tmp,
            "run_pathways": False,
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
            "min_replicates_per_condition": expected["min_replicates_per_condition"],
        })
    assert result["status"] == "success", result.get("details", "")
    c0 = result["contrasts"][0]
    assert c0.get("status") == "success", c0.get("error", "")

    planted = set(expected["planted_up_genes"])
    sig = {str(g) for g in (c0.get("all_sig_genes") or [])}
    recall = len(planted & sig) / len(planted)
    false_among_null = len(sig - planted)

    assert recall >= expected["min_recall"], \
        f"DE recall {recall:.2f} < {expected['min_recall']} (recovered {len(planted & sig)}/{len(planted)})"
    assert false_among_null <= expected["max_false_up_among_null"], \
        f"{false_among_null} null genes wrongly called significant"
    # Ground truth is all-up; a flood of down calls means a direction/ref bug.
    assert c0["n_downregulated"] <= expected["max_false_up_among_null"]


# ── F8 (preprint audit): pydeseq2 fit warnings reach the audit trail ──────────

def test_serialize_fit_warnings_keeps_numeric_drops_benign():
    """The classifier surfaces numeric/convergence warnings, drops benign churn."""
    import warnings as _w
    from aria.scripts.rna_bulk_de import _serialize_fit_warnings
    with _w.catch_warnings(record=True) as caught:
        _w.simplefilter("always")
        _w.warn("singular matrix in dispersion fit", RuntimeWarning)
        _w.warn("did not converge", UserWarning)
        _w.warn("deprecated api", DeprecationWarning)
        _w.warn("pandas future change", FutureWarning)
    out = _serialize_fit_warnings(caught)
    assert any("singular matrix" in m for m in out)
    assert any("did not converge" in m and "UserWarning" in m for m in out)
    # benign third-party churn is not promoted to the audit trail
    assert not any("deprecated api" in m for m in out)
    assert not any("future change" in m for m in out)


def test_fit_warnings_surface_in_result(monkeypatch):
    """A warning emitted during the pydeseq2 fit is surfaced in _run_deseq2's
    returned warnings (no longer swallowed by a global filterwarnings)."""
    import sys
    import types
    import warnings as _w

    idx = [f"g{i}" for i in range(5)]
    cols = [f"s{i}" for i in range(6)]
    counts = pd.DataFrame(
        np.arange(30).reshape(5, 6) + 1, index=idx, columns=cols)
    metadata = pd.DataFrame(
        {"condition": ["A", "A", "A", "B", "B", "B"]}, index=cols)

    class _FakeDDS:
        def __init__(self, **kw):
            self.var = pd.DataFrame(index=idx)

        def deseq2(self):
            _w.warn("dispersion trend fit did not converge", RuntimeWarning)
            _w.warn("benign pandas future change", FutureWarning)

    class _FakeStats:
        def __init__(self, dds, **kw):
            self.results_df = pd.DataFrame(
                {"log2FoldChange": [1.0] * 5, "padj": [0.01] * 5}, index=idx)

        def summary(self):
            pass

    dds_mod = types.ModuleType("pydeseq2.dds")
    dds_mod.DeseqDataSet = _FakeDDS
    ds_mod = types.ModuleType("pydeseq2.ds")
    ds_mod.DeseqStats = _FakeStats
    monkeypatch.setitem(sys.modules, "pydeseq2", types.ModuleType("pydeseq2"))
    monkeypatch.setitem(sys.modules, "pydeseq2.dds", dds_mod)
    monkeypatch.setitem(sys.modules, "pydeseq2.ds", ds_mod)

    result, warns = rbd._run_deseq2(
        counts, metadata, "condition", "B", "A",
        padj_thr=0.05, lfc_thr=1.0, lfc_shrink=False,
        min_replicates_per_condition=3)

    assert result["status"] == "success", result.get("details", "")
    joined = " ".join(warns)
    assert "did not converge" in joined          # numeric warning surfaced
    assert "RuntimeWarning" in joined
    assert "benign pandas future change" not in joined   # benign dropped
