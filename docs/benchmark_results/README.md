# ARIA Benchmark & Concordance Artifacts

This directory stores ARIA's committed, machine-relativized benchmark manifests, figures, rubrics, and cryptographic preprint-freeze receipts. Every JSON artifact in this tree is guarded by `tests/test_public_artifacts_no_absolute_paths.py` so zero machine-absolute paths (`/home/...`, `/tmp/...`) can ever be committed.

---

## Directory Layout

| Subdirectory | Role & Authority |
| :--- | :--- |
| **`preprint_v1/`** | **Authoritative Preprint Freeze Root (`aria/benchmarks/preprint_freeze.py`).** Contains `inventory.json`, `receipts/*.json` (SHA-256 + `source_snapshot_hash` + `workflow_hash` bound to clean git state), and verified artifacts for **Claims 1–7** (`claim_1/` through `claim_5/`). |
| **`multiome_concordance/`** | **V48 Paired 10x Multiome External Concordance & Orchestrated Validation (`GSE278576` `hc11_paired.h5mu`, 2,517 post-QC cells, 13 hippocampal subclasses):** WNN vs R `Seurat 5.3.0` (`v48_i2`), 5-axis WNN anti-degradation guardrail (`v48_i3`), MOFA+ vs `muon.tl.mofa` (`v48_i4`), paired `peak2gene` vs R `Signac::LinkPeaks` (`v48_i5`), candidate `eRegulons` vs `scenicplus==1.0.0` (`v48_i6`), non-circular `IntegrationNarrator` E2E (`v48_i7`), live orchestrated CP3.5 dispatch (`v48_i8`), and cross-modal `HypothesisAgent S10` Arms A–D (`v48_i9`). |
| **`scatac_concordance/`** | **scATAC External Concordance (`ADR-048`, `ADR-050`):** `SnapATAC2` clustering/DA concordance (`p3a`), 5v5 multi-donor `GSE278576` pseudobulk DA vs R `DESeq2`/`edgeR`/`limma` (`p3`), `Signac::LinkPeaks` (`p4_peak2gene`), and `Signac::GeneActivity` (`p4_geneactivity`). |
| **`bulk_atac/`** | **V47 / ADR-077 Bulk ATAC End-to-End & External R Concordance Validation:** Multi-replicate (`n=4v4` + technical replicates + batch covariate + 20 null permutations) concordance against R `DESeq2`, `edgeR`, and `limma-voom` (`adr077_multireplicate_r_concordance.json`), plus ENCODE `K562` vs `GM12878` raw alignment (`b2a`), replicate-gated DESeq2 DA (`v47_*_da_validation.json`), and motif enrichment (`v47_*_motif_enrichment.json`), guarded by `tests/test_bulk_atac_v47_provenance.py`. |
| **`bulk_atac_footprint/` & `scatac_footprint/`** | **TF Footprinting (`TOBIAS` BINDetect, `B7`):** Replicate-aware inference (`ENCODE` `chr21+22` + label permutations) and descriptive single-pool candidate ranking, guarded by `tests/test_preprint_audit_b7_footprint_inference.py`. |
| **`b2_claim/`** | **Track B2 Claim Faithfulness Rubric:** `B2_RUBRIC.md` defining the frozen 5-class failure taxonomy (`unsupported`, `overclaim`, `fabricated`, `causal_inflation`, `missing_caveat`) used by both `governance_b2.py` and the multi-annotator kit (`b2_annotation.py`). |
| **`b3_frozen_reagent/`** | **Track B3 Frozen Open-Weight Reagent Ablation Supplement:** SHA-256-pinned open-weight reagent deterministic replay across ARIA's 6 core LLM callsites (`B3_FROZEN_REAGENT_SUPPLEMENT.md`, JSON manifest, and supplemental SVG figure). |
| **`geo_connector/`** | **Track E6 Atomic GEO/SRA Retrieval Validation:** `e6_atomic_retrieval.json` across all 4 primary modalities (`bulk_RNA`, `scRNA`, `scATAC`, `bulk_ATAC`), guarded by `tests/test_preprint_audit_q1_doc_claims.py`. |
| **`hypothesis/`** | **Single-Modality `HypothesisAgent` Red-Team Evaluation:** `hypothesis_redteam.json` (`S0`–`S9` adversarial gate evaluation). |

---

## Regenerating Preprint Freeze Evidence

```bash
conda activate aria-env

# Inspect current receipt status across all 21 required lanes
python scripts/run_preprint_freeze_inventory.py

# Re-execute a registered lane from a clean checkout
python scripts/run_preprint_freeze_lane.py <lane_id>
```

See [Preprint-v1 Evidence Freeze](../architecture/preprint_freeze_v1.md) and [`scripts/README.md`](../../scripts/README.md) for the complete lane mapping.
