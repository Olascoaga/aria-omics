# Maintenance, Resource & Benchmark Runner Scripts

`scripts/` contains developer, release-maintenance, resource-bootstrap, and reproducible benchmark/validation CLI runners. Runtime subprocesses dispatched by ARIA's agents during user analyses live in `aria/scripts/` (not here) and are governed by `aria/utils/script_contracts.py`.

---

## 1. Release, Lockfile & Resource Utilities

| Script | Purpose |
| :--- | :--- |
| `generate_locks.sh` | Maintainer-only generator of exact `@EXPLICIT` Conda (`envs/*.linux-64.lock`) and PyPI (`envs/*.pip.lock`) lockfiles for all active stacks in `ENVIRONMENT_SPECS`. |
| `update_version_metadata.py` | Synchronizes version badges and validates `docs/release_notes_v{__version__}.md` against `aria/version.py`. |
| `fetch_genesets.py` | Downloads and hashes local GMT pathway collections for offline hypergeometric ORA (`~/.aria/genesets/`). |
| `fetch_motifs.py` | Stages offline MEME v4 vertebrate TF motif collections (`JASPAR2024`) with SHA-256 provenance. |
| `fetch_seqc_maqc_reference.py` | Bootstraps the Bioconductor `seqc` RefSeq + TaqMan qPCR + Thermo ERCC reference bundles into `~/.aria/benchmarks/`. |

---

## 2. Preprint Evidence Freeze Engine (`FASE 8`)

| Script | Purpose |
| :--- | :--- |
| `run_preprint_freeze_inventory.py` | Recomputes SHA-256 artifact hashes and `source_snapshot_hash` / `workflow_hash` provenance across all 21 Claims 1–7 lanes -> `docs/benchmark_results/preprint_v1/inventory.json`. |
| `run_preprint_freeze_lane.py` | Executes a single registered freeze lane from a clean source tree and writes its cryptographic receipt under `docs/benchmark_results/preprint_v1/receipts/`. |

---

## 3. Benchmark A — Statistical Calibration & External Concordance

| Script | Modality / Comparator |
| :--- | :--- |
| `run_a1_bulk_de_benchmark.py` | Synthetic-truth bulk RNA-seq DE calibration + `lfcThreshold` recall/precision frontier (`Claim 1`). |
| `run_a1_external_comparators.py` | Bulk DE vs R `DESeq2`, `edgeR-QLF`, and `limma-voom` in `aria-bench-env` (`Claim 1`). |
| `run_a1_seqc_maqc_benchmark.py` | Single-site (`BGI`) SEQC/MAQC bulk DE vs TaqMan qPCR truth (`Claim 1`). |
| `run_a1_seqc_multisite_benchmark.py` | Five-site (`BGI`, `CNL`, `MAY`, `AGR`, `NVS`) cross-lab SEQC reproducibility (`Claim 1`). |
| `run_a1_ercc_dose_response.py` | ERCC 92 spike-in fold-change recovery and 6.3-order dynamic-range linearity (`Claim 1`). |
| `run_c1_h9_fastq_e2e.py` | Real 9-sample paired-FASTQ -> STAR -> featureCounts -> DESeq2 -> governed report capsule E2E (`Claim 1`). |
| `run_a2_pseudobulk_benchmark.py` | Donor-aware scRNA pseudobulk DE vs naive cell-level pseudoreplication anti-pattern (`Claim 2`). |
| `run_a2_external_comparators.py` | Real Kang et al. 2018 (8v8 PBMC IFN-$\beta$) ARIA pseudobulk DE vs `muscat` (`edgeR-QLF`) (`Claim 2`). |
| `aria_pseudobulk_da_from_tsv.py` | Shared TSV-to-`_run_deseq2` bridge used by multi-donor chromatin DA comparators (`Claim 2` / `Claim 5`). |
| `aria_atac_pseudobulk_matrix.py` | Builds donor-by-peak pseudobulk count matrices for chromatin null and DA benchmarks (`Claim 5`). |

---

## 4. Benchmark B — Design, Narrative & Null Governance

| Script | Claim / Governance Axis |
| :--- | :--- |
| `run_b1_design_governance.py` | 23-case adversarial experimental design inference/escalation/refusal matrix (`Claim 3`). |
| `run_c3_design_blind_gold.py` | Exports blind sheets and scores independent human design-gold CSVs (`Claim 3`). |
| `run_b2_claim_narrative.py` | 4-arm claim faithfulness & causal-overreach ablation (`governed`, `guards_off`, `naive`, `template`) (`Claim 4`). |
| `run_b2_annotation_kit.py` | Exports blind multi-annotator CSVs and computes Cohen's $\kappa$, Krippendorff's $\alpha$, and adjudicated gold (`Claim 4`). |
| `run_c4_report_e2e_false_narrative.py` | End-to-end bus -> claim compiler -> HTML report faithfulness harness + human scoring CLI (`Claim 4`). |
| `run_b4_null_narrative.py` | Null-evidence narrative fabrication and false-flag calibration (`Claim 5`). |
| `run_c5_multimodal_null_permutations.py` | 20 RNA + 20 ATAC label-permutation DESeq2 -> `compile_public_claims` null calibration (`Claim 5`). |
| `run_hypothesis_blind_eval.py` | Adversarial red-team evaluation harness for single-modality `HypothesisAgent` (`S0`–`S9`). |
| `run_e6_geo_sra_validation.py` | Multi-modality closed-world GEO/SRA staging & manifest validation runner (`E6`). |

---

## 5. Chromatin & Paired 10x Multiome Concordance Suites (V46–V48)

| Script | Validation Scope |
| :--- | :--- |
| `run_scatac_da_benchmark.py` | Synthetic replicated scATAC & bulk ATAC pseudobulk DA recovery/FDR (`W-CALIB`). |
| `run_scatac_multisample_validation.py` | Multi-donor consensus-peak scATAC pseudobulk DA live execution harness (`ADR-042` / `ADR-043`). |
| `run_scatac_multisample_da_concordance.py` | 5 young vs 5 old `GSE278576` donor `.h5mu` scATAC DA vs R `DESeq2`/`edgeR`/`limma` (`Claim 2`, `ADR-048`). |
| `run_scatac_concordance_benchmark.py` | scATAC TF-IDF/LSI clustering & marker concordance vs `SnapATAC2`. |
| `run_scatac_peak2gene_concordance.py` | Chromatin `peak2gene` concordance vs R `Signac::LinkPeaks` (`ADR-050`). |
| `run_scatac_geneactivity_concordance.py` | Peak-window gene-activity concordance vs R `Signac::GeneActivity`. |
| `run_scatac_footprint_tobias.py` & `run_bulk_atac_footprint_tobias.py` | `TOBIAS` Tn5-corrected footprinting validation runners (`scATAC` and `ENCODE` bulk ATAC `B7`). |
| `run_tobias_footprint_concordance.py` | Multi-replicate (`n=4v4` Welch, paired-donor, batch-adjusted OLS + 20 null cohorts) + live `TOBIAS 0.17.3` footprinting concordance vs R `stats` (`ADR-078`). |
| `run_bulk_atac_concordance.py` | Multi-replicate (`n=4v4` + tech-rep + batch + 20 null permutations) bulk ATAC DA concordance vs R `DESeq2`, `edgeR`, and `limma-voom` (`ADR-077`). |
| `run_bulk_atac_marker_mapping.py` & `run_bulk_atac_v47_provenance.py` | Bulk ATAC `K562` vs `GM12878` regulatory marker concordance and artifact provenance verifier. |
| `run_wnn_concordance.py` | Paired 10x Multiome WNN concordance vs R `Seurat 5.3.0` `FindMultiModalNeighbors` (`V48 I2`, `ADR-065`). |
| `run_wnn_guardrail_validation.py` | 5-axis WNN anti-degradation guardrail validation on `GSE278576` `hc11` (`V48 I3`, `ADR-066`). |
| `run_mofa_concordance.py` | Paired MOFA+ latent factor & loading concordance vs `muon.tl.mofa` (`V48 I4`, `ADR-067`). |
| `run_peak2gene_paired_concordance.py` | Unified paired `peak2gene_core` concordance vs R `Signac::LinkPeaks` (`V48 I5`, `ADR-068`). |
| `run_eregulon_concordance.py` | Candidate `eRegulons` 4-tier concordance vs upstream `scenicplus==1.0.0` (`V48 I6`, `ADR-069`). |
| `run_v48_i7_narrative_validation.py` | Non-circular Cross-Modal Evidence Graph, `IntegrationNarrator`, and dual PNG+SVG E2E runner (`V48 I7`, `ADR-070`). |
| `run_v48_i8_orchestrated_validation.py` | Live `OrchestratorAgent` CP3.5 acknowledgement + `IntegrationAgent` E2E validation (`V48 I8`, `ADR-071`). |
| `run_v48_i9_cross_modal_hypotheses.py` | Cross-modal `HypothesisAgent S10` Arms A–D + live LLM E2E validation (`V48 I9`, `ADR-072`). |
| `run_adr079_multidonor_multiome_validation.py` | 4-donor `GSE278576` cohort union, donor-aware Harmony WNN, MOFA+ donor $\eta^2$ screening, barcode-permuted null `peak2gene`, and `C3` ledger production validation (`ADR-079`). |
