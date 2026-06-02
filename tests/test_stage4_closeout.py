"""Stage 4 closeout guards: privacy, stats gates, and robustness manifests."""

from __future__ import annotations

import json
import time

from aria.agents.narrative.claim_compiler import classify_claim
from aria.agents.narrative.robustness import build_robustness_multiverse
from aria.agents.narrative.types import EvidenceItem, NarrativeBlock
from aria.llm.provider import LLMProvider, ModelConfig, TaskTier
from aria.utils.environment_manager import EnvironmentManager
from aria.utils.privacy import redact_sensitive_params


def _de_block(**metrics):
    return NarrativeBlock(
        id="scrna.pseudobulk.GroupA.condA_vs_condB",
        modality="scRNA-seq",
        analysis="pseudobulk_de",
        block_type="result",
        title="GroupA condA_vs_condB",
        status="success",
        confidence="medium",
        claim="GroupA condA_vs_condB had DE genes.",
        evidence=[EvidenceItem("DE genes", metrics.get("n_significant", 1), "test")],
        metrics=metrics,
    )


def test_privacy_redacts_paths_and_secrets():
    redacted = redact_sensitive_params({
        "data_path": "/private/study/sample.h5ad",
        "api_token": "secret-value",
        "nested": {"output_dir": "/tmp/aria/report"},
    })
    assert redacted["data_path"] == "<path:sample.h5ad>"
    assert redacted["api_token"] == "<redacted>"
    assert redacted["nested"]["output_dir"] == "<path:report>"


def test_failed_run_archive_redacts_input_json(tmp_path, monkeypatch):
    monkeypatch.delenv("ARIA_PRESERVE_FAILED_INPUTS", raising=False)
    mgr = EnvironmentManager(workspace_dir=str(tmp_path))
    input_file = tmp_path / "input_abcd.json"
    input_file.write_text(
        json.dumps({"data_path": "/private/study/sample.h5ad", "alpha": 0.05}),
        encoding="utf-8",
    )
    output_file = tmp_path / "output_abcd.json"

    mgr._archive_failed_run("abcd", "rna", input_file, output_file, {"status": "error"})

    run_dir = tmp_path / "failed" / "rna_abcd"
    assert not (run_dir / "input.json").exists()
    payload = json.loads((run_dir / "input.redacted.json").read_text())
    assert payload["data_path"] == "<path:sample.h5ad>"
    assert "/private/study" not in (run_dir / "input.redacted.json").read_text()


def test_air_gapped_provider_filters_cloud_models(monkeypatch, tmp_path):
    monkeypatch.setenv("ARIA_AIR_GAPPED", "1")
    local = ModelConfig("ollama", "ollama/test", 8000, is_local=True)
    cloud = ModelConfig("openai", "gpt-test", 8000, is_local=False)
    provider = LLMProvider(
        models={
            TaskTier.HEAVY: [cloud, local],
            TaskTier.MEDIUM: [local],
            TaskTier.LIGHT: [local],
        },
        cache_dir=str(tmp_path),
    )
    assert [m.model for m in provider.models[TaskTier.HEAVY]] == ["ollama/test"]


def test_llm_cache_ttl_expires_stale_entries(monkeypatch, tmp_path):
    monkeypatch.delenv("ARIA_AIR_GAPPED", raising=False)
    monkeypatch.setenv("ARIA_LLM_CACHE_TTL_DAYS", "1")
    local = ModelConfig("ollama", "ollama/test", 8000, is_local=True)
    provider = LLMProvider(
        models={TaskTier.HEAVY: [local], TaskTier.MEDIUM: [local], TaskTier.LIGHT: [local]},
        cache_dir=str(tmp_path),
    )
    key = provider._cache_key("m", "s", "p", 1, 0.0, 0, provider._cache_version_salt)
    path = provider._cache_path(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({
        "text": "stale",
        "created_unix": time.time() - 3 * 86400,
        "version_salt": provider._cache_version_salt,
    }))
    assert provider._cache_get(key) is None


def test_claim_compiler_stats_gate_downgrades_weak_de_claim():
    weak = classify_claim(_de_block(
        n_significant=10,
        power_estimate_at_effective_alpha=0.2,
    ))
    assert weak.tier == "descriptive"
    assert any("Stats-evidence gate" in lim for lim in weak.limitations)

    supported = classify_claim(_de_block(
        n_significant=10,
        power_estimate_at_effective_alpha=0.8,
    ))
    assert supported.tier == "associative"


def test_robustness_multiverse_manifest_from_pseudobulk_results():
    manifest = build_robustness_multiverse({
        "scrna_agent": {
            "findings": {
                "pseudobulk_de": {
                    "per_group": {
                        "GroupA": {
                            "per_comparison": {
                                "condA_vs_condB": {
                                    "status": "success",
                                    "corrected_for_composition": True,
                                    "robustness_multiverse": {
                                        "stable_significant_genes": 3,
                                        "stability_basis": "gene_id_intersection",
                                        "stable_gene_ids": ["GeneA", "GeneB", "GeneC"],
                                        "fdr_axis_evaluated": True,
                                        "n_local": 5,
                                        "n_global": 3,
                                        "fdr_family_variants": {
                                            "per_cluster": {"n_significant": 5},
                                            "global": {"n_significant": 3},
                                        },
                                    },
                                }
                            }
                        }
                    }
                }
            }
        }
    })
    assert manifest["status"] == "available"
    assert manifest["entries"][0]["stable_significant_genes"] == 3
    assert manifest["entries"][0]["stability_status"] == "computed"
    assert manifest["entries"][0]["stability_basis"] == "gene_id_intersection"
    assert manifest["entries"][0]["stable_gene_ids"] == ["GeneA", "GeneB", "GeneC"]
    assert manifest["entries"][0]["fdr_axis_evaluated"] is True
    assert manifest["entries"][0]["composition_covariate"] == "included"


def test_robustness_multiverse_does_not_invent_intersection_from_counts():
    manifest = build_robustness_multiverse({
        "scrna_agent": {
            "findings": {
                "pseudobulk_de": {
                    "per_group": {
                        "GroupA": {
                            "per_comparison": {
                                "condA_vs_condB": {
                                    "status": "success",
                                    "corrected_for_composition": False,
                                    "n_significant_local": 8,
                                    "n_significant_global": 5,
                                }
                            }
                        }
                    }
                }
            }
        }
    })

    entry = manifest["entries"][0]
    assert entry["n_local_fdr"] == 8
    assert entry["n_global_fdr"] == 5
    assert entry["stable_significant_genes"] is None
    assert entry["stability_status"] == "not_computed"
    assert entry["stability_basis"] is None
