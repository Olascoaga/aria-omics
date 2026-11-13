# ARIA Benchmark & Concordance Artifacts

This directory stores ARIA's committed, machine-relativized benchmark manifests, figures, rubrics, and cryptographic preprint-freeze receipts. Every JSON artifact in this tree is guarded by `tests/test_public_artifacts_no_absolute_paths.py` so zero machine-absolute paths (`/home/...`, `/tmp/...`) can ever be committed.

---

## Directory Layout

| Subdirectory | Role & Authority |
| :--- | :--- |
| **`preprint_v1/`** | **Authoritative Preprint Freeze Root (`aria/benchmarks/preprint_freeze.py`).** Contains `inventory.json`, `receipts/*.json` (SHA-256 + `source_snapshot_hash` + `workflow_hash` bound to clean git state), and verified artifacts for **Claims 1–7** (`claim_1/` through `claim_5/`). |
| **`multiome_concordance/`** | **V48, ADR-079 Paired Same-Cell Multiome & ADR-080 Cross-Modal Contrast/Donor Concordance Validation (`GSE278576` human hippocampus):** 10-donor non-same-cell RNA $\times$ ATAC contrast + donor-matched cis concordance & donor-deranged null receipt (`adr080_cross_modal_concordance_validation.json`), 4-donor balanced (`2x2`) cohort union + donor-aware Harmony WNN + ANOVA $\eta^2$ MOFA+ + barcode-permuted null `peak2gene` receipt (`adr079_multidonor_multiome_validation.json`), WNN vs R `Seurat 5.3.0` (`v48_i2`), 5-axis WNN anti-degradation guardrail (`v48_i3`), MOFA+ vs `muon.tl.mofa` (`v48_i4`), paired `peak2gene` vs R `Signac::LinkPeaks` (`v48_i5`), candidate `eRegulons` vs `scenicplus==1.0.0` (`v48_i6`), non-circular `IntegrationNarrator` E2E (`v48_i7`), live orchestrated dispatch (`v48_i8`), and cross-modal `HypothesisAgent S10` Arms A–D (`v48_i9`). |
| **`scatac_concordance/`** | **scATAC External Concordance (`ADR-048`, `ADR-050`):** `SnapATAC2` clustering/DA concordance (`p3a`), 5v5 multi-donor `GSE278576` pseudobulk DA vs R `DESeq2`/`edgeR`/`limma` (`p3`), `Signac::LinkPeaks` (`p4_peak2gene`), and `Signac::GeneActivity` (`p4_geneactivity`). |
| **`bulk_atac/`** | **V47 / ADR-077 Bulk ATAC End-to-End & External R Concordance Validation:** Multi-replicate (`n=4v4` + technical replicates + batch covariate + 20 null permutations) concordance against R `DESeq2`, `edgeR`, and `limma-voom` (`adr077_multireplicate_r_concordance.json`), plus ENCODE `K562` vs `GM12878` raw alignment (`b2a`), replicate-gated DESeq2 DA (`v47_*_da_validation.json`), and motif enrichment (`v47_*_motif_enrichment.json`), guarded by `tests/test_bulk_atac_v47_provenance.py`. |
| **`bulk_atac_footprint/` & `scatac_footprint/`** | **TF Footprinting (`TOBIAS` BINDetect, `B7` / `ADR-078`):** External R (`t.test` Welch & paired + `lm` batch OLS + BH) + live `TOBIAS 0.17.3` concordance & 20-null calibration (`adr078_tobias_replicate_concordance.json`), `ENCODE` `chr21+22` replicate inference (`b7_*`), and descriptive single-pool candidate ranking, guarded by `tests/test_preprint_audit_b7_footprint_inference.py`. |
| **`epigenomics/`** | **Epigenomic Protein–DNA Interaction Validation (`ChIP`, `CUT_AND_RUN`, `CUT_AND_TAG`, `ADR-081`):** `adr081_chip_cutrun_cuttag_validation.json` certifying wet-lab mouse `mm39` + `sacCer3` spike-in `CUT&RUN` (`SEACR` stringent/relaxed, Gaussian-copula `IDR`, BigWig export, spike-in `pyDESeq2` vs R Bioconductor `DESeq2` `r = 1.0000`), unidirectional global-loss spike-in recovery (`240/240` recalled, `0` false gains vs `60` naive median-of-ratios false gains) + `0 / 3,000` null false positives, wet-lab `mm10` `ChIP-seq` broad megadomains vs external `EDD` (`95.92%` bp recall), `CUT&Tag` `SEACR` calling, and `C1–C4` claim governance. |
| **`hic/`** | **3D Genome Architecture Validation (`HiC` / `Micro-C`, `ADR-082`):** `adr082_hic_3d_genome_validation.json` certifying `.pairs -> .mcool` multi-resolution preprocessing + `A4` resume, `ICE` balancing & `HiCRep` `SCC`, standalone `cooltools 0.7.1` concordance (`E1` Pearson `r = 1.0000`, Crane insulation Pearson `r = 1.0000`, `100%` TAD boundary recall, `dotfinder.dots` loop + GTF anchor gene recovery), `4v4` empirical-Bayes `differential_3d` `A_to_B` / `B_to_A` switching + `1v1` descriptive guard + `0` null false discoveries, dual PNG+SVG figures, and `C1–C4` + `HypothesisAgent` governance. |
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
